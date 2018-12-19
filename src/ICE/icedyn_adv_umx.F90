MODULE icedyn_adv_umx
   !!==============================================================================
   !!                       ***  MODULE  icedyn_adv_umx  ***
   !! sea-ice : advection using the ULTIMATE-MACHO scheme
   !!==============================================================================
   !! History :  3.6  !  2014-11  (C. Rousset, G. Madec)  Original code
   !!            4.0  !  2018     (many people)           SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv_umx   : update the tracer trend with the 3D advection trends using a TVD scheme
   !!   ultimate_x(_y)    : compute a tracer value at velocity points using ULTIMATE scheme at various orders
   !!   macho             : ???
   !!   nonosc            : compute monotonic tracer fluxes by a non-oscillatory algorithm 
   !!----------------------------------------------------------------------
   USE phycst         ! physical constant
   USE dom_oce        ! ocean domain
   USE sbc_oce , ONLY : nn_fsbc   ! update frequency of surface boundary condition
   USE ice            ! sea-ice variables
   USE icevar         ! sea-ice: operations
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_adv_umx   ! called by icedyn_adv.F90
      
   REAL(wp) ::   z1_6   = 1._wp /   6._wp   ! =1/6
   REAL(wp) ::   z1_120 = 1._wp / 120._wp   ! =1/120

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: amaxu, amaxv
   
   ! advect H all the way (and get V=H*A at the end)
   LOGICAL :: ll_thickness = .FALSE.
   
   ! look for 9 points around in nonosc limiter  
   LOGICAL :: ll_9points = .FALSE.  ! false=better h?

   ! use HgradU in nonosc limiter  
   LOGICAL :: ll_HgradU = .TRUE.   ! no effect?

   ! if T interpolated at u/v points is negative, then interpolate T at u/v points using the upstream scheme
   LOGICAL :: ll_neg = .TRUE.      ! keep TRUE
   
   ! limit the fluxes
   LOGICAL :: ll_zeroup1 = .FALSE. ! false ok if Hbig otherwise needed for 2D sinon on a des valeurs de H trop fortes !!
   LOGICAL :: ll_zeroup2 = .FALSE. ! false ok for 1D, 2D, 3D
   LOGICAL :: ll_zeroup4 = .FALSE. ! false ok for 1D, 2D, 3D
   LOGICAL :: ll_zeroup5 = .FALSE. ! false ok for 1D, 2D

   ! fluxes that are limited are u*H, then (u*H)*(ua/u) is used for V (only for nonosc)
   LOGICAL :: ll_clem   = .TRUE.   ! simpler than rachid and works

   ! First advect H as H*=Hdiv(u), then use H* for H(n+1)=H(n)-div(uH*)
   LOGICAL :: ll_gurvan = .FALSE.  ! must be false for 1D case !!

   ! First guess as div(uH) (-true-) or Hdiv(u)+ugradH (-false-)
   LOGICAL :: ll_1stguess_clem = .FALSE. ! better negative values but less good h

   ! advect (or not) open water. If not, retrieve it from advection of A
   LOGICAL :: ll_ADVopw = .FALSE.  !
   
   ! alternate directions for upstream
   LOGICAL :: ll_upsxy = .TRUE.

   ! alternate directions for high order
   LOGICAL :: ll_hoxy = .TRUE.
   
   ! prelimiter: use it to avoid overshoot in H
   LOGICAL :: ll_prelimiter_zalesak = .TRUE.  ! from: Zalesak(1979) eq. 14 => true is better for 1D but false is better in 3D (for h and negative values) => pb in x-y?
   LOGICAL :: ll_prelimiter_devore  = .FALSE.  ! from: Devore eq. 11 (worth than zalesak)

   ! iterate on the limiter (only nonosc)
   LOGICAL :: ll_limiter_it2 = .FALSE.
   

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_adv_umx( kn_umx, kt, pu_ice, pv_ice,  &
      &                        pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pe_s, pe_i )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_adv_umx  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of 
      !!                 tracers and add it to the general trend of tracer equations
      !!                 using an "Ultimate-Macho" scheme
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kn_umx     ! order of the scheme (1-5=UM or 20=CEN2)
      INTEGER                     , INTENT(in   ) ::   kt         ! time step
      REAL(wp), DIMENSION(:,:)    , INTENT(in   ) ::   pu_ice     ! ice i-velocity
      REAL(wp), DIMENSION(:,:)    , INTENT(in   ) ::   pv_ice     ! ice j-velocity
      REAL(wp), DIMENSION(:,:)    , INTENT(inout) ::   pato_i     ! open water area
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_i       ! ice volume
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_s       ! snw volume
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   psv_i      ! salt content
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   poa_i      ! age content
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pa_i       ! ice concentration
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pa_ip      ! melt pond fraction
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_ip      ! melt pond volume
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pe_s       ! snw heat content
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pe_i       ! ice heat content
      !
      INTEGER  ::   ji, jj, jk, jl, jt      ! dummy loop indices
      INTEGER  ::   icycle                  ! number of sub-timestep for the advection
      REAL(wp) ::   zamsk                   ! 1 if advection of concentration, 0 if advection of other tracers
      REAL(wp) ::   zdt
      REAL(wp), DIMENSION(1)       ::   zcflprv, zcflnow   ! send zcflnow and receive zcflprv
      REAL(wp), DIMENSION(jpi,jpj) ::   zudy, zvdx
      REAL(wp), DIMENSION(jpi,jpj) ::   zati1, zati2



      REAL(wp), DIMENSION(jpi,jpj)     ::   zcu_box, zcv_box
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zua_ho, zva_ho
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z1_ai, z1_aip
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zhvar

      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_dyn_adv_umx: Ultimate-Macho advection scheme'
      !
      !
      ! --- If ice drift field is too fast, use an appropriate time step for advection (CFL test for stability) --- !
      !     When needed, the advection split is applied at the next time-step in order to avoid blocking global comm.
      !     ...this should not affect too much the stability... Was ok on the tests we did...
      zcflnow(1) =                  MAXVAL( ABS( pu_ice(:,:) ) * rdt_ice * r1_e1u(:,:) )
      zcflnow(1) = MAX( zcflnow(1), MAXVAL( ABS( pv_ice(:,:) ) * rdt_ice * r1_e2v(:,:) ) )
      
      ! non-blocking global communication send zcflnow and receive zcflprv
      CALL mpp_delay_max( 'icedyn_adv_umx', 'cflice', zcflnow(:), zcflprv(:), kt == nitend - nn_fsbc + 1 )

      IF( zcflprv(1) > .5 ) THEN   ;   icycle = 2
      ELSE                         ;   icycle = 1
      ENDIF
      
      zdt = rdt_ice / REAL(icycle)

      ! --- transport --- !
      zudy(:,:) = pu_ice(:,:) * e2u(:,:)
      zvdx(:,:) = pv_ice(:,:) * e1v(:,:)

      ! --- define velocity for advection: u*grad(H) --- !
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            IF    ( pu_ice(ji,jj) * pu_ice(ji-1,jj) <= 0._wp ) THEN   ;   zcu_box(ji,jj) = 0._wp
            ELSEIF( pu_ice(ji,jj)                   >  0._wp ) THEN   ;   zcu_box(ji,jj) = pu_ice(ji-1,jj)
            ELSE                                                      ;   zcu_box(ji,jj) = pu_ice(ji  ,jj)
            ENDIF

            IF    ( pv_ice(ji,jj) * pv_ice(ji,jj-1) <= 0._wp ) THEN   ;   zcv_box(ji,jj) = 0._wp
            ELSEIF( pv_ice(ji,jj)                   >  0._wp ) THEN   ;   zcv_box(ji,jj) = pv_ice(ji,jj-1)
            ELSE                                                      ;   zcv_box(ji,jj) = pv_ice(ji,jj  )
            ENDIF
         END DO
      END DO

      IF( ll_zeroup2 ) THEN
         IF(.NOT. ALLOCATED(amaxu))       ALLOCATE(amaxu (jpi,jpj,jpl))
         IF(.NOT. ALLOCATED(amaxv))       ALLOCATE(amaxv (jpi,jpj,jpl))
      ENDIF
      !---------------!
      !== advection ==!
      !---------------!
      DO jt = 1, icycle

!!$         IF( ll_ADVopw ) THEN
!!$            zamsk = 1._wp
!!$            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zudy, zvdx, zcu_box, zcv_box, pato_i(:,:), pato_i(:,:) )   ! Open water area 
!!$            zamsk = 0._wp
!!$         ELSE
            zati1(:,:) = SUM( pa_i(:,:,:), dim=3 )
!!$         ENDIF
         
         WHERE( pa_i(:,:,:) >= epsi20 )   ;   z1_ai(:,:,:) = 1._wp / pa_i(:,:,:)
         ELSEWHERE                        ;   z1_ai(:,:,:) = 0.
         END WHERE
            !
         WHERE( pa_ip(:,:,:) >= epsi20 )  ;   z1_aip(:,:,:) = 1._wp / pa_ip(:,:,:)
         ELSEWHERE                        ;   z1_aip(:,:,:) = 0.
         END WHERE
         !
         IF( ll_zeroup2 ) THEN
            DO jl = 1, jpl
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     amaxu(ji,jj,jl)=MAX( pa_i(ji,jj,jl), pa_i(ji,jj-1,jl), pa_i(ji,jj+1,jl), &
                        &                                 pa_i(ji+1,jj,jl), pa_i(ji+1,jj-1,jl), pa_i(ji+1,jj+1,jl) )
                     amaxv(ji,jj,jl)=MAX( pa_i(ji,jj,jl), pa_i(ji-1,jj,jl), pa_i(ji+1,jj,jl), &
                        &                                 pa_i(ji,jj+1,jl), pa_i(ji-1,jj+1,jl), pa_i(ji+1,jj+1,jl) )
                  END DO
               END DO
            END DO
            CALL lbc_lnk_multi('icedyn_adv_umx', amaxu, 'T', 1., amaxv, 'T', 1.)
         ENDIF
         !
         DO jl = 1, jpl
            zua_ho(:,:,jl) = zudy(:,:)
            zva_ho(:,:,jl) = zvdx(:,:)
         END DO
         
         zamsk = 1._wp
         CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, pa_i, pa_i, zua_ho, zva_ho ) ! Ice area
         zamsk = 0._wp
         !
         IF( ll_thickness ) THEN
            DO jl = 1, jpl
               zua_ho(:,:,jl) = zudy(:,:)
               zva_ho(:,:,jl) = zvdx(:,:)
            END DO
         ENDIF
            !
         zhvar(:,:,:) = pv_i(:,:,:) * z1_ai(:,:,:)
         CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho , zva_ho , zcu_box, zcv_box, zhvar, pv_i )    ! Ice volume
         IF( ll_thickness )   pv_i(:,:,:) = zhvar(:,:,:) * pa_i(:,:,:)
         !
         zhvar(:,:,:) = pv_s(:,:,:) * z1_ai(:,:,:)
         CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho , zva_ho , zcu_box, zcv_box, zhvar, pv_s )    ! Snw volume
         IF( ll_thickness )   pv_s(:,:,:) = zhvar(:,:,:) * pa_i(:,:,:)
         !
         zhvar(:,:,:) = psv_i(:,:,:) * z1_ai(:,:,:)
         CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho , zva_ho , zcu_box, zcv_box, zhvar, psv_i )    ! Salt content
         !
         zhvar(:,:,:) = poa_i(:,:,:) * z1_ai(:,:,:)
         CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho , zva_ho , zcu_box, zcv_box, zhvar, poa_i )    ! Age content
         !
         DO jk = 1, nlay_i
            zhvar(:,:,:) = pe_i(:,:,jk,:) * z1_ai(:,:,:)
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pe_i(:,:,jk,:) ) ! Ice heat content
         END DO
         !
         DO jk = 1, nlay_s
            zhvar(:,:,:) = pe_s(:,:,jk,:) * z1_ai(:,:,:)
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pe_s(:,:,jk,:) ) ! Snw heat content
         END DO
            !
         IF ( ln_pnd_H12 ) THEN
            !
            DO jl = 1, jpl
               zua_ho(:,:,jl) = zudy(:,:)
               zva_ho(:,:,jl) = zvdx(:,:)
            END DO
            
            zamsk = 1._wp
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, pa_ip, pa_ip, zua_ho, zva_ho ) ! mp fraction
            zamsk = 0._wp
            !
            zhvar(:,:,:) = pv_ip(:,:,:) * z1_ai(:,:,:)
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho , zva_ho , zcu_box, zcv_box, zhvar, pv_ip ) ! mp volume
         ENDIF
         !
         !
!!$         IF( .NOT. ll_ADVopw ) THEN
            zati2(:,:) = SUM( pa_i(:,:,:), dim=3 )
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  pato_i(ji,jj) = pato_i(ji,jj) - ( zati2(ji,jj) - zati1(ji,jj) ) &                                                  ! Open water area
                     &                          - ( zudy(ji,jj) - zudy(ji-1,jj) + zvdx(ji,jj) - zvdx(ji,jj-1) )*r1_e1e2t(ji,jj)*zdt
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', pato_i(:,:), 'T',  1. )
!!$         ENDIF
         !
      END DO
      !
   END SUBROUTINE ice_dyn_adv_umx

   
   SUBROUTINE adv_umx( pamsk, kn_umx, jt, kt, pdt, pu, pv, puc, pvc, pubox, pvbox, pt, ptc, pua_ho, pva_ho )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE adv_umx  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of 
      !!       tracers and add it to the general trend of tracer equations
      !!
      !! **  Method  :   TVD scheme, i.e. 2nd order centered scheme with
      !!       corrected flux (monotonic correction)
      !!       note: - this advection scheme needs a leap-frog time scheme
      !!
      !! ** Action : - pt  the after advective tracer
      !!----------------------------------------------------------------------
      REAL(wp)                    , INTENT(in   )           ::   pamsk          ! advection of concentration (1) or other tracers (0)
      INTEGER                     , INTENT(in   )           ::   kn_umx         ! order of the scheme (1-5=UM or 20=CEN2)
      INTEGER                     , INTENT(in   )           ::   jt             ! number of sub-iteration
      INTEGER                     , INTENT(in   )           ::   kt             ! number of iteration
      REAL(wp)                    , INTENT(in   )           ::   pdt            ! tracer time-step
      REAL(wp), DIMENSION(:,:  ), INTENT(in   )           ::   pu   , pv      ! 2 ice velocity components => u*e2
      REAL(wp), DIMENSION(:,:,:), INTENT(in   )           ::   puc  , pvc     ! 2 ice velocity components => u*e2 or u*a*e2u
      REAL(wp), DIMENSION(:,:  ), INTENT(in   )           ::   pubox, pvbox   ! upstream velocity
      REAL(wp), DIMENSION(:,:,:), INTENT(inout)           ::   pt             ! tracer field
      REAL(wp), DIMENSION(:,:,:), INTENT(inout)           ::   ptc            ! tracer content field
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out), OPTIONAL ::   pua_ho, pva_ho ! high order u*a fluxes
      !
      INTEGER  ::   ji, jj, jl       ! dummy loop indices  
      REAL(wp) ::   ztra             ! local scalar
      INTEGER  ::   kn_limiter = 1   ! 1=nonosc ; 2=superbee ; 3=h3(rachid)
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zfu_ho , zfv_ho , zt_u, zt_v, zpt
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zfu_ups, zfv_ups, zt_ups   ! only for nonosc 
      !!----------------------------------------------------------------------
      !
      !  upstream (_ups) advection with initial mass fluxes
      ! ---------------------------------------------------

      IF( ll_gurvan .AND. pamsk==0. ) THEN
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  pt(ji,jj,jl) = ( pt (ji,jj,jl) + pt(ji,jj,jl) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj)     &
                     &                           + pt(ji,jj,jl) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) )   &
                     &           * tmask(ji,jj,1)
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'icedyn_adv_umx', pt, 'T', 1. )
      ENDIF

      
      IF( .NOT. ll_upsxy ) THEN

         ! fluxes in both x-y directions
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ll_clem ) THEN
                     zfu_ups(ji,jj,jl) = MAX( pu(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pu(ji,jj), 0._wp ) * pt(ji+1,jj,jl)
                     zfv_ups(ji,jj,jl) = MAX( pv(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pv(ji,jj), 0._wp ) * pt(ji,jj+1,jl)
                  ELSE
                     zfu_ups(ji,jj,jl) = MAX( puc(ji,jj,jl), 0._wp ) * pt(ji,jj,jl) + MIN( puc(ji,jj,jl), 0._wp ) * pt(ji+1,jj,jl)
                     zfv_ups(ji,jj,jl) = MAX( pvc(ji,jj,jl), 0._wp ) * pt(ji,jj,jl) + MIN( pvc(ji,jj,jl), 0._wp ) * pt(ji,jj+1,jl)
                  ENDIF
               END DO
            END DO
         END DO

      ELSE
         !
         IF( MOD( (kt - 1) / nn_fsbc , 2 ) ==  MOD( (jt - 1) , 2 ) ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
            ! flux in x-direction
            DO jl = 1, jpl
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     IF( ll_clem ) THEN
                        zfu_ups(ji,jj,jl) = MAX( pu(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pu(ji,jj), 0._wp ) * pt(ji+1,jj,jl)
                     ELSE
                        zfu_ups(ji,jj,jl) = MAX( puc(ji,jj,jl), 0._wp ) * pt(ji,jj,jl) + MIN( puc(ji,jj,jl), 0._wp ) * pt(ji+1,jj,jl)
                     ENDIF
                  END DO
               END DO
            END DO
            
            ! first guess of tracer content from u-flux
            DO jl = 1, jpl
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     IF( ll_clem ) THEN
                        IF( ll_gurvan ) THEN
                           zpt(ji,jj,jl) = ( pt(ji,jj,jl) - ( zfu_ups(ji,jj,jl) - zfu_ups(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                        ELSE
                           zpt(ji,jj,jl) = ( pt(ji,jj,jl) - ( zfu_ups(ji,jj,jl) - zfu_ups(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                              &            + pt(ji,jj,jl) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                              &            ) * tmask(ji,jj,1)
                        ENDIF
                     ELSE
                        zpt(ji,jj,jl) = ( ptc(ji,jj,jl) - ( zfu_ups(ji,jj,jl) - zfu_ups(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj) ) &
                           &         * tmask(ji,jj,1)
                     ENDIF
                     !!                  IF( ji==26 .AND. jj==86) THEN
                     !!                     WRITE(numout,*) '************************'
                     !!                     WRITE(numout,*) 'zpt upstream',zpt(ji,jj)
                     !!                  ENDIF
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zpt, 'T', 1. )
            !
            ! flux in y-direction
            DO jl = 1, jpl
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     IF( ll_clem ) THEN
                        zfv_ups(ji,jj,jl) = MAX( pv(ji,jj), 0._wp ) * zpt(ji,jj,jl) + MIN( pv(ji,jj), 0._wp ) * zpt(ji,jj+1,jl)
                     ELSE
                        zfv_ups(ji,jj,jl) = MAX( pvc(ji,jj,jl), 0._wp ) * zpt(ji,jj,jl) + MIN( pvc(ji,jj,jl), 0._wp ) * zpt(ji,jj+1,jl)
                     ENDIF
                  END DO
               END DO
            END DO
         !
         ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
            ! flux in y-direction
            DO jl = 1, jpl
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     IF( ll_clem ) THEN
                        zfv_ups(ji,jj,jl) = MAX( pv(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pv(ji,jj), 0._wp ) * pt(ji,jj+1,jl)
                     ELSE
                        zfv_ups(ji,jj,jl) = MAX( pvc(ji,jj,jl), 0._wp ) * pt(ji,jj,jl) + MIN( pvc(ji,jj,jl), 0._wp ) * pt(ji,jj+1,jl)
                     ENDIF
                  END DO
               END DO
            END DO

            ! first guess of tracer content from v-flux
            DO jl = 1, jpl
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     IF( ll_clem ) THEN
                        IF( ll_gurvan ) THEN
                           zpt(ji,jj,jl) = ( pt(ji,jj,jl) - ( zfv_ups(ji,jj,jl) - zfv_ups(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                        ELSE
                           zpt(ji,jj,jl) = ( pt(ji,jj,jl) - ( zfv_ups(ji,jj,jl) - zfv_ups(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) &
                              &            + pt(ji,jj,jl) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) ) &
                              &            * tmask(ji,jj,1)
                        ENDIF
                     ELSE
                        zpt(ji,jj,jl) = ( ptc(ji,jj,jl) - ( zfv_ups(ji,jj,jl) - zfv_ups(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) ) &
                           &            * tmask(ji,jj,1)
                     ENDIF
                     !!                  IF( ji==26 .AND. jj==86) THEN
                     !!                     WRITE(numout,*) '************************'
                     !!                     WRITE(numout,*) 'zpt upstream',zpt(ji,jj)
                     !!                  ENDIF
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zpt, 'T', 1. )
            !
            ! flux in x-direction
            DO jl = 1, jpl
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     IF( ll_clem ) THEN
                        zfu_ups(ji,jj,jl) = MAX( pu(ji,jj), 0._wp ) * zpt(ji,jj,jl) + MIN( pu(ji,jj), 0._wp ) * zpt(ji+1,jj,jl)
                     ELSE
                        zfu_ups(ji,jj,jl) = MAX( puc(ji,jj,jl), 0._wp ) * zpt(ji,jj,jl) + MIN( puc(ji,jj,jl), 0._wp ) * zpt(ji+1,jj,jl)
                     ENDIF
                  END DO
               END DO
            END DO
            !
         ENDIF
         
      ENDIF

      IF( ll_clem .AND. kn_limiter /= 1 )  &
         & CALL ctl_stop( 'STOP', 'icedyn_adv_umx: ll_clem incompatible with limiters other than nonosc' )

      IF( ll_zeroup2 ) THEN
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  IF( amaxu(ji,jj,jl) == 0._wp )   zfu_ups(ji,jj,jl) = 0._wp
                  IF( amaxv(ji,jj,jl) == 0._wp )   zfv_ups(ji,jj,jl) = 0._wp
               END DO
            END DO
         END DO
      ENDIF

      ! guess after content field with upstream scheme
      DO jl = 1, jpl
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               ztra          = - (   zfu_ups(ji,jj,jl) - zfu_ups(ji-1,jj  ,jl)   &
                  &                + zfv_ups(ji,jj,jl) - zfv_ups(ji  ,jj-1,jl) ) * r1_e1e2t(ji,jj)
               IF( ll_clem ) THEN
                  IF( ll_gurvan ) THEN
                     zt_ups(ji,jj,jl) = ( pt (ji,jj,jl) + pdt * ztra ) * tmask(ji,jj,1)
                  ELSE
                     zt_ups(ji,jj,jl) = ( pt (ji,jj,jl) + pdt * ztra + ( pt(ji,jj,jl) * pdt * ( pu(ji,jj) - pu(ji-1,jj) )   &
                        &                                            +   pt(ji,jj,jl) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) ) &
                        &                                              * r1_e1e2t(ji,jj) * (1.-pamsk) ) * tmask(ji,jj,1)
                  ENDIF
               ELSE
                  zt_ups(ji,jj,jl) = ( ptc(ji,jj,jl) + pdt * ztra ) * tmask(ji,jj,1)
               ENDIF
               !!            IF( ji==26 .AND. jj==86) THEN
               !!               WRITE(numout,*) '**************************'
               !!               WRITE(numout,*) 'zt upstream',zt_ups(ji,jj)
               !!            ENDIF
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', zt_ups, 'T', 1. )

      ! High order (_ho) fluxes 
      ! -----------------------
      SELECT CASE( kn_umx )
         !
      CASE ( 20 )                          !== centered second order ==!
         !
         CALL cen2( pamsk, kn_limiter, jt, kt, pdt, pt, pu, pv, puc, pvc, ptc, zfu_ho, zfv_ho,  &
            &       zt_ups, zfu_ups, zfv_ups )
         !
      CASE ( 1:5 )                         !== 1st to 5th order ULTIMATE-MACHO scheme ==!
         !
         CALL macho( pamsk, kn_limiter, kn_umx, jt, kt, pdt, pt, pu, pv, puc, pvc, pubox, pvbox, ptc, zt_u, zt_v, zfu_ho, zfv_ho,  &
            &        zt_ups, zfu_ups, zfv_ups )
         !
      END SELECT

      IF( ll_thickness ) THEN
         ! final trend with corrected fluxes
         ! ------------------------------------
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  IF( ll_gurvan ) THEN
                     ztra       = - ( zfu_ho(ji,jj,jl) - zfu_ho(ji-1,jj,jl) + zfv_ho(ji,jj,jl) - zfv_ho(ji,jj-1,jl) ) * r1_e1e2t(ji,jj) 
                  ELSE
                     ztra       = ( - ( zfu_ho(ji,jj,jl) - zfu_ho(ji-1,jj,jl) + zfv_ho(ji,jj,jl) - zfv_ho(ji,jj-1,jl) )  & 
                        &           + ( pt(ji,jj,jl) * ( pu(ji,jj) - pu(ji-1,jj) ) * (1.-pamsk) ) &
                        &           + ( pt(ji,jj,jl) * ( pv(ji,jj) - pv(ji,jj-1) ) * (1.-pamsk) ) ) * r1_e1e2t(ji,jj)
                  ENDIF
                  pt(ji,jj,jl) = ( pt(ji,jj,jl) + pdt * ztra ) * tmask(ji,jj,1)
                  
                  IF( pt(ji,jj,jl) < -epsi20 ) THEN
                     WRITE(numout,*) 'T<0 ',pt(ji,jj,jl)
                  ENDIF
                  
                  IF( pt(ji,jj,jl) < 0._wp .AND. pt(ji,jj,jl) >= -epsi20 )   pt(ji,jj,jl) = 0._wp
                  
                  !!               IF( ji==26 .AND. jj==86) THEN
                  !!                  WRITE(numout,*) 'zt high order',pt(ji,jj)
                  !!               ENDIF
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'icedyn_adv_umx', pt, 'T',  1. )
      ENDIF
   
      ! Rachid trick
      ! ------------
      IF( ll_clem ) THEN
         IF( pamsk == 0. ) THEN
            DO jl = 1, jpl
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     IF( ABS( puc(ji,jj,jl) ) > 0._wp .AND. ABS( pu(ji,jj) ) > 0._wp ) THEN
                        zfu_ho (ji,jj,jl) = zfu_ho (ji,jj,jl) * puc(ji,jj,jl) / pu(ji,jj)
                        zfu_ups(ji,jj,jl) = zfu_ups(ji,jj,jl) * puc(ji,jj,jl) / pu(ji,jj)
                     ELSE
                        zfu_ho (ji,jj,jl) = 0._wp
                        zfu_ups(ji,jj,jl) = 0._wp
                     ENDIF
                     !
                     IF( ABS( pvc(ji,jj,jl) ) > 0._wp .AND. ABS( pv(ji,jj) ) > 0._wp ) THEN
                        zfv_ho (ji,jj,jl) = zfv_ho (ji,jj,jl) * pvc(ji,jj,jl) / pv(ji,jj)
                        zfv_ups(ji,jj,jl) = zfv_ups(ji,jj,jl) * pvc(ji,jj,jl) / pv(ji,jj)
                     ELSE
                        zfv_ho (ji,jj,jl) = 0._wp  
                        zfv_ups(ji,jj,jl) = 0._wp  
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF

      IF( ll_zeroup5 ) THEN
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = 2, fs_jpim1   ! vector opt.
                  zpt(ji,jj,jl) = ( ptc(ji,jj,jl) - ( zfu_ho(ji,jj,jl) - zfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                            - ( zfv_ho(ji,jj,jl) - zfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj)  ) * tmask(ji,jj,1)
                  IF( zpt(ji,jj,jl) < 0. ) THEN
                     zfu_ho(ji  ,jj,jl) = zfu_ups(ji  ,jj,jl)
                     zfu_ho(ji-1,jj,jl) = zfu_ups(ji-1,jj,jl)
                     zfv_ho(ji  ,jj,jl) = zfv_ups(ji  ,jj,jl)
                     zfv_ho(ji,jj-1,jl) = zfv_ups(ji,jj-1,jl)
                  ENDIF
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'icedyn_adv_umx', zfu_ho, 'U',  -1., zfv_ho, 'V',  -1. )
      ENDIF

      ! output high order fluxes u*a
      ! ----------------------------
      IF( PRESENT( pua_ho ) ) THEN
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  pua_ho(ji,jj,jl) = zfu_ho(ji,jj,jl)
                  pva_ho(ji,jj,jl) = zfv_ho(ji,jj,jl)
               END DO
            END DO
         END DO
      ENDIF


      IF( .NOT.ll_thickness ) THEN
         ! final trend with corrected fluxes
         ! ------------------------------------
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1 
                  ztra = - ( zfu_ho(ji,jj,jl) - zfu_ho(ji-1,jj,jl) + zfv_ho(ji,jj,jl) - zfv_ho(ji,jj-1,jl) ) * r1_e1e2t(ji,jj) * pdt  
                  
                  ptc(ji,jj,jl) = ( ptc(ji,jj,jl) + ztra ) * tmask(ji,jj,1)
                  
                  !!               IF( ji==26 .AND. jj==86) THEN
                  !!                  WRITE(numout,*) 'ztc high order',ptc(ji,jj)
                  !!               ENDIF
                  
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'icedyn_adv_umx', ptc, 'T',  1. )
      ENDIF
      
      !
   END SUBROUTINE adv_umx

   SUBROUTINE cen2( pamsk, kn_limiter, jt, kt, pdt, pt, pu, pv, puc, pvc, ptc, pfu_ho, pfv_ho, &
      &             pt_ups, pfu_ups, pfv_ups )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE macho  ***
      !!     
      !! **  Purpose :   compute  
      !!
      !! **  Method  :   ... ???
      !!                 TIM = transient interpolation Modeling 
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      REAL(wp)                    , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      INTEGER                     , INTENT(in   ) ::   kn_limiter       ! limiter
      INTEGER                     , INTENT(in   ) ::   jt               ! number of sub-iteration
      INTEGER                     , INTENT(in   ) ::   kt               ! number of iteration
      REAL(wp)                    , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pt               ! tracer fields
      REAL(wp), DIMENSION(:,:  ), INTENT(in   ) ::   pu, pv           ! 2 ice velocity components
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   puc, pvc         ! 2 ice velocity * A components
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ptc              ! tracer content at before time step 
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfu_ho, pfv_ho   ! high order fluxes 
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pt_ups           ! upstream guess of tracer content 
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pfu_ups, pfv_ups ! upstream fluxes 
      !
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      LOGICAL  ::   ll_xy = .TRUE.
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zzt
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ll_xy ) THEN   !-- no alternate directions --!
         !
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ll_clem ) THEN
                     pfu_ho(ji,jj,jl) = 0.5 * pu(ji,jj) * ( pt(ji,jj,jl) + pt(ji+1,jj,jl) )
                     pfv_ho(ji,jj,jl) = 0.5 * pv(ji,jj) * ( pt(ji,jj,jl) + pt(ji,jj+1,jl) )
                  ELSE
                     pfu_ho(ji,jj,jl) = 0.5 * puc(ji,jj,jl) * ( pt(ji,jj,jl) + pt(ji+1,jj,jl) )
                     pfv_ho(ji,jj,jl) = 0.5 * pvc(ji,jj,jl) * ( pt(ji,jj,jl) + pt(ji,jj+1,jl) )
                  ENDIF
               END DO
            END DO
         END DO
         IF    ( kn_limiter == 1 ) THEN
            IF( ll_clem ) THEN
               CALL nonosc_2d( pamsk, pdt, pu, puc, pv, pvc, ptc, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
            ELSE
               CALL nonosc_2d( pamsk, pdt, pu, puc, pv, pvc, ptc, ptc, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
            ENDIF
         ELSEIF( kn_limiter == 2 ) THEN
            CALL limiter_x( pdt, pu, puc, pt, pfu_ho )
            CALL limiter_y( pdt, pv, pvc, pt, pfv_ho )
         ELSEIF( kn_limiter == 3 ) THEN
            CALL limiter_x( pdt, pu, puc, pt, pfu_ho, pfu_ups )
            CALL limiter_y( pdt, pv, pvc, pt, pfv_ho, pfv_ups )
         ENDIF
         !
      ELSE                    !-- alternate directions --!
         !
         IF( MOD( (kt - 1) / nn_fsbc , 2 ) ==  MOD( (jt - 1) , 2 ) ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
            !
            ! flux in x-direction
            DO jl = 1, jpl
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     IF( ll_clem ) THEN
                        pfu_ho(ji,jj,jl) = 0.5 * pu(ji,jj) * ( pt(ji,jj,jl) + pt(ji+1,jj,jl) )
                     ELSE
                        pfu_ho(ji,jj,jl) = 0.5 * puc(ji,jj,jl) * ( pt(ji,jj,jl) + pt(ji+1,jj,jl) )
                     ENDIF
                  END DO
               END DO
            END DO
            IF( kn_limiter == 2 )   CALL limiter_x( pdt, pu, puc, pt, pfu_ho )
            IF( kn_limiter == 3 )   CALL limiter_x( pdt, pu, puc, pt, pfu_ho, pfu_ups )

            ! first guess of tracer content from u-flux
            DO jl = 1, jpl
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     IF( ll_clem ) THEN
                        zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)        &
                           &            + pt(ji,jj,jl) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) ) &
                           &         * tmask(ji,jj,1)
                     ELSE                     
                        zzt(ji,jj,jl) = ( ptc(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                     ENDIF
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zzt, 'T', 1. )

            ! flux in y-direction
            DO jl = 1, jpl
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     IF( ll_clem ) THEN
                        pfv_ho(ji,jj,jl) = 0.5 * pv(ji,jj) * ( zzt(ji,jj,jl) + zzt(ji,jj+1,jl) )
                     ELSE                     
                        pfv_ho(ji,jj,jl) = 0.5 * pvc(ji,jj,jl) * ( zzt(ji,jj,jl) + zzt(ji,jj+1,jl) )
                     ENDIF
                  END DO
               END DO
            END DO
            IF( kn_limiter == 2 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho )
            IF( kn_limiter == 3 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho, pfv_ups )

         ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
            !
            ! flux in y-direction
            DO jl = 1, jpl
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     IF( ll_clem ) THEN
                        pfv_ho(ji,jj,jl) = 0.5 * pv(ji,jj) * ( pt(ji,jj,jl) + pt(ji,jj+1,jl) )
                     ELSE                     
                        pfv_ho(ji,jj,jl) = 0.5 * pvc(ji,jj,jl) * ( pt(ji,jj,jl) + pt(ji,jj+1,jl) )
                     ENDIF
                  END DO
               END DO
            END DO
            IF( kn_limiter == 2 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho )
            IF( kn_limiter == 3 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho, pfv_ups )
            !
            ! first guess of tracer content from v-flux
            DO jl = 1, jpl
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     IF( ll_clem ) THEN
                        zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) &
                           &                     + pt(ji,jj,jl) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) ) &
                           &         * tmask(ji,jj,1)
                     ELSE
                        zzt(ji,jj,jl) = ( ptc(ji,jj,jl) - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                     ENDIF
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zzt, 'T', 1. )
            !
            ! flux in x-direction
            DO jl = 1, jpl
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     IF( ll_clem ) THEN
                        pfu_ho(ji,jj,jl) = 0.5 * pu(ji,jj) * ( zzt(ji,jj,jl) + zzt(ji+1,jj,jl) )
                     ELSE
                        pfu_ho(ji,jj,jl) = 0.5 * puc(ji,jj,jl) * ( zzt(ji,jj,jl) + zzt(ji+1,jj,jl) )
                     ENDIF
                  END DO
               END DO
            END DO
            IF( kn_limiter == 2 )   CALL limiter_x( pdt, pu, puc, pt, pfu_ho )
            IF( kn_limiter == 3 )   CALL limiter_x( pdt, pu, puc, pt, pfu_ho, pfu_ups )

         ENDIF
         IF( ll_clem ) THEN
            IF( kn_limiter == 1 )   CALL nonosc_2d( pamsk, pdt, pu, puc, pv, pvc, ptc, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
         ELSE
            IF( kn_limiter == 1 )   CALL nonosc_2d( pamsk, pdt, pu, puc, pv, pvc, ptc, ptc, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
         ENDIF
         
      ENDIF
   
   END SUBROUTINE cen2

   
   SUBROUTINE macho( pamsk, kn_limiter, kn_umx, jt, kt, pdt, pt, pu, pv, puc, pvc, pubox, pvbox, ptc, pt_u, pt_v, pfu_ho, pfv_ho, &
      &              pt_ups, pfu_ups, pfv_ups )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE macho  ***
      !!     
      !! **  Purpose :   compute  
      !!
      !! **  Method  :   ... ???
      !!                 TIM = transient interpolation Modeling 
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      REAL(wp)                    , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      INTEGER                     , INTENT(in   ) ::   kn_limiter       ! limiter
      INTEGER                     , INTENT(in   ) ::   kn_umx           ! order of the scheme (1-5=UM or 20=CEN2)
      INTEGER                     , INTENT(in   ) ::   jt               ! number of sub-iteration
      INTEGER                     , INTENT(in   ) ::   kt               ! number of iteration
      REAL(wp)                    , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pt               ! tracer fields
      REAL(wp), DIMENSION(:,:  ), INTENT(in   ) ::   pu, pv           ! 2 ice velocity components
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   puc, pvc         ! 2 ice velocity * A components
      REAL(wp), DIMENSION(:,:  ), INTENT(in   ) ::   pubox, pvbox     ! upstream velocity
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ptc              ! tracer content at before time step 
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pt_u, pt_v       ! tracer at u- and v-points 
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfu_ho, pfv_ho   ! high order fluxes 
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pt_ups           ! upstream guess of tracer content 
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pfu_ups, pfv_ups ! upstream fluxes 
      !
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp) ::   ztra
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zzt, zzfu_ho, zzfv_ho
      !!----------------------------------------------------------------------
      !
      IF( MOD( (kt - 1) / nn_fsbc , 2 ) ==  MOD( (jt - 1) , 2 ) ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
         !
         !                                                        !--  ultimate interpolation of pt at u-point  --!
         CALL ultimate_x( kn_umx, pdt, pt, pu, puc, pt_u, pfu_ho )
         !                                                        !--  limiter in x --!
         IF( kn_limiter == 2 )   CALL limiter_x( pdt, pu, puc, pt, pfu_ho )
         IF( kn_limiter == 3 )   CALL limiter_x( pdt, pu, puc, pt, pfu_ho, pfu_ups )
         !                                                        !--  advective form update in zzt  --!

         IF( ll_1stguess_clem ) THEN

            ! first guess of tracer content from u-flux
            DO jl = 1, jpl
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     IF( ll_clem ) THEN
                        IF( ll_gurvan ) THEN
                           zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                        ELSE
                           zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                              &            + pt(ji,jj,jl) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                              &         ) * tmask(ji,jj,1)
                        ENDIF
                     ELSE
                        zzt(ji,jj,jl) = ( ptc(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                     ENDIF
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zzt, 'T', 1. )

         ELSE

            DO jl = 1, jpl
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     IF( ll_gurvan ) THEN
                        zzt(ji,jj,jl) = pt(ji,jj,jl) - pubox(ji,jj   ) * pdt * ( pt_u(ji,jj,jl) - pt_u(ji-1,jj,jl) ) * r1_e1t(ji,jj)  &
                           &                         - pt   (ji,jj,jl) * pdt * ( pu  (ji,jj) - pu  (ji-1,jj) ) * r1_e1e2t(ji,jj)
                     ELSE
                        zzt(ji,jj,jl) = pt(ji,jj,jl) - pubox(ji,jj   ) * pdt * ( pt_u(ji,jj,jl) - pt_u(ji-1,jj,jl) ) * r1_e1t(ji,jj)  &
                           &                         - pt   (ji,jj,jl) * pdt * ( pu  (ji,jj) - pu  (ji-1,jj) ) * r1_e1e2t(ji,jj) * pamsk
                     ENDIF
                     zzt(ji,jj,jl) = zzt(ji,jj,jl) * tmask(ji,jj,1)
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zzt, 'T', 1. )
         ENDIF
         !
         !                                                        !--  ultimate interpolation of pt at v-point  --!
         IF( ll_hoxy ) THEN
            CALL ultimate_y( kn_umx, pdt, zzt, pv, pvc, pt_v, pfv_ho )
         ELSE
            CALL ultimate_y( kn_umx, pdt, pt, pv, pvc, pt_v, pfv_ho )
         ENDIF
         !                                                        !--  limiter in y --!
         IF( kn_limiter == 2 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho )
         IF( kn_limiter == 3 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho, pfv_ups )
         !         
         !
      ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
         !
         !                                                        !--  ultimate interpolation of pt at v-point  --!
         CALL ultimate_y( kn_umx, pdt, pt, pv, pvc, pt_v, pfv_ho )
         !                                                        !--  limiter in y --!
         IF( kn_limiter == 2 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho )
         IF( kn_limiter == 3 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho, pfv_ups )
         !                                                        !--  advective form update in zzt  --!
         IF( ll_1stguess_clem ) THEN
            
            ! first guess of tracer content from v-flux 
            DO jl = 1, jpl
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     IF( ll_clem ) THEN
                        IF( ll_gurvan ) THEN
                           zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                        ELSE
                           zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) &
                              &            + pt(ji,jj,jl) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                              &         ) * tmask(ji,jj,1)
                        ENDIF
                     ELSE
                        zzt(ji,jj,jl) = ( ptc(ji,jj,jl) - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) ) &
                           &         * tmask(ji,jj,1)
                     ENDIF
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zzt, 'T', 1. )
            
         ELSE
            
            DO jl = 1, jpl
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     IF( ll_gurvan ) THEN
                        zzt(ji,jj,jl) = pt(ji,jj,jl) - pvbox(ji,jj   ) * pdt * ( pt_v(ji,jj,jl) - pt_v(ji,jj-1,jl) ) * r1_e2t(ji,jj)  &
                           &                         - pt   (ji,jj,jl) * pdt * ( pv  (ji,jj) - pv  (ji,jj-1) ) * r1_e1e2t(ji,jj)
                     ELSE
                        zzt(ji,jj,jl) = pt(ji,jj,jl) - pvbox(ji,jj   ) * pdt * ( pt_v(ji,jj,jl) - pt_v(ji,jj-1,jl) ) * r1_e2t(ji,jj)  &
                           &                         - pt   (ji,jj,jl) * pdt * ( pv  (ji,jj) - pv  (ji,jj-1) ) * r1_e1e2t(ji,jj) * pamsk
                     ENDIF
                     zzt(ji,jj,jl) = zzt(ji,jj,jl) * tmask(ji,jj,1)
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zzt, 'T', 1. )
         ENDIF
         !
         !                                                        !--  ultimate interpolation of pt at u-point  --!
         IF( ll_hoxy ) THEN
            CALL ultimate_x( kn_umx, pdt, zzt, pu, puc, pt_u, pfu_ho )
         ELSE
            CALL ultimate_x( kn_umx, pdt, pt, pu, puc, pt_u, pfu_ho )
         ENDIF
         !                                                        !--  limiter in x --!
         IF( kn_limiter == 2 )   CALL limiter_x( pdt, pu, puc, pt, pfu_ho )
         IF( kn_limiter == 3 )   CALL limiter_x( pdt, pu, puc, pt, pfu_ho, pfu_ups )
         !
         !
      ENDIF

     
      IF( kn_limiter == 1 ) THEN
         IF( .NOT. ll_limiter_it2 ) THEN
            IF( ll_clem ) THEN
               CALL nonosc_2d ( pamsk, pdt, pu, puc, pv, pvc, ptc, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
            ELSE
               CALL nonosc_2d ( pamsk, pdt, pu, puc, pv, pvc, ptc, ptc, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
            ENDIF
         ELSE
            zzfu_ho(:,:,:) = pfu_ho(:,:,:)
            zzfv_ho(:,:,:) = pfv_ho(:,:,:)
            ! 1st iteration of nonosc (limit the flux with the upstream solution)
            IF( ll_clem ) THEN
               CALL nonosc_2d ( pamsk, pdt, pu, puc, pv, pvc, ptc, pt, pt_ups, pfu_ups, pfv_ups, zzfu_ho, zzfv_ho )
            ELSE
               CALL nonosc_2d ( pamsk, pdt, pu, puc, pv, pvc, ptc, ptc, pt_ups, pfu_ups, pfv_ups, zzfu_ho, zzfv_ho )
            ENDIF
            ! guess after content field with high order
            DO jl = 1, jpl
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ztra          = - ( zzfu_ho(ji,jj,jl) - zzfu_ho(ji-1,jj,jl) + zzfv_ho(ji,jj,jl) - zzfv_ho(ji,jj-1,jl) ) * r1_e1e2t(ji,jj)
                     zzt(ji,jj,jl) =   ( ptc(ji,jj,jl) + pdt * ztra ) * tmask(ji,jj,1)
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zzt, 'T', 1. )
            ! 2nd iteration of nonosc (limit the flux with the limited high order solution)
            IF( ll_clem ) THEN
               CALL nonosc_2d ( pamsk, pdt, pu, puc, pv, pvc, ptc, pt, zzt, zzfu_ho, zzfv_ho, pfu_ho, pfv_ho )
            ELSE
               CALL nonosc_2d ( pamsk, pdt, pu, puc, pv, pvc, ptc, ptc, zzt, zzfu_ho, zzfv_ho, pfu_ho, pfv_ho )
            ENDIF
         ENDIF
      ENDIF
      !
   END SUBROUTINE macho


   SUBROUTINE ultimate_x( kn_umx, pdt, pt, pu, puc, pt_u, pfu_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ultimate_x  ***
      !!     
      !! **  Purpose :   compute  
      !!
      !! **  Method  :   ... ???
      !!                 TIM = transient interpolation Modeling 
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kn_umx    ! order of the scheme (1-5=UM or 20=CEN2)
      REAL(wp)                    , INTENT(in   ) ::   pdt       ! tracer time-step
      REAL(wp), DIMENSION(:,:  ), INTENT(in   ) ::   pu        ! ice i-velocity component
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   puc       ! ice i-velocity * A component
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pt_u      ! tracer at u-point 
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfu_ho    ! high order flux 
      !
      INTEGER  ::   ji, jj, jl             ! dummy loop indices
      REAL(wp) ::   zcu, zdx2, zdx4    !   -      -
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   ztu1, ztu2, ztu3, ztu4
      !!----------------------------------------------------------------------
      !
      !                                                     !--  Laplacian in i-direction  --!
      DO jl = 1, jpl
         DO jj = 2, jpjm1         ! First derivative (gradient)
            DO ji = 1, fs_jpim1
               ztu1(ji,jj,jl) = ( pt(ji+1,jj,jl) - pt(ji,jj,jl) ) * r1_e1u(ji,jj) * umask(ji,jj,1)
            END DO
            !                     ! Second derivative (Laplacian)
            DO ji = fs_2, fs_jpim1
               ztu2(ji,jj,jl) = ( ztu1(ji,jj,jl) - ztu1(ji-1,jj,jl) ) * r1_e1t(ji,jj)
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', ztu2, 'T', 1. )
      !
      !                                                     !--  BiLaplacian in i-direction  --!
      DO jl = 1, jpl
         DO jj = 2, jpjm1         ! Third derivative
            DO ji = 1, fs_jpim1
               ztu3(ji,jj,jl) = ( ztu2(ji+1,jj,jl) - ztu2(ji,jj,jl) ) * r1_e1u(ji,jj) * umask(ji,jj,1)
            END DO
            !                     ! Fourth derivative
            DO ji = fs_2, fs_jpim1
               ztu4(ji,jj,jl) = ( ztu3(ji,jj,jl) - ztu3(ji-1,jj,jl) ) * r1_e1t(ji,jj)
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', ztu4, 'T', 1. )
      !
      !
      SELECT CASE (kn_umx )
      !
      CASE( 1 )                                                   !==  1st order central TIM  ==! (Eq. 21)
         !        
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (                           pt(ji+1,jj,jl) + pt(ji,jj,jl)   &
                     &                                    - SIGN( 1._wp, pu(ji,jj) ) * ( pt(ji+1,jj,jl) - pt(ji,jj,jl) ) )
               END DO
            END DO
         END DO
         !
      CASE( 2 )                                                   !==  2nd order central TIM  ==! (Eq. 23)
         !
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
                  pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (                                pt(ji+1,jj,jl) + pt(ji,jj,jl)   &
                     &                                               -              zcu   * ( pt(ji+1,jj,jl) - pt(ji,jj,jl) ) ) 
               END DO
            END DO
         END DO
         !  
      CASE( 3 )                                                   !==  3rd order central TIM  ==! (Eq. 24)
         !
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
                  zdx2 = e1u(ji,jj) * e1u(ji,jj)
!!rachid       zdx2 = e1u(ji,jj) * e1t(ji,jj)
                  pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (         (                      pt  (ji+1,jj,jl) + pt  (ji,jj,jl)        &
                     &                                               -              zcu   * ( pt  (ji+1,jj,jl) - pt  (ji,jj,jl) )  )   &
                     &        + z1_6 * zdx2 * ( zcu*zcu - 1._wp ) * (                         ztu2(ji+1,jj,jl) + ztu2(ji,jj,jl)        &
                     &                                               - SIGN( 1._wp, zcu ) * ( ztu2(ji+1,jj,jl) - ztu2(ji,jj,jl) )  )   )
               END DO
            END DO
         END DO
         !
      CASE( 4 )                                                   !==  4th order central TIM  ==! (Eq. 27)
         !
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
                  zdx2 = e1u(ji,jj) * e1u(ji,jj)
!!rachid       zdx2 = e1u(ji,jj) * e1t(ji,jj)
                  pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (         (                pt  (ji+1,jj,jl) + pt  (ji,jj,jl)        &
                     &                                               -          zcu * ( pt  (ji+1,jj,jl) - pt  (ji,jj,jl) )  )   &
                     &        + z1_6 * zdx2 * ( zcu*zcu - 1._wp ) * (                   ztu2(ji+1,jj,jl) + ztu2(ji,jj,jl)        &
                     &                                               - 0.5_wp * zcu * ( ztu2(ji+1,jj,jl) - ztu2(ji,jj,jl) )  )   )
               END DO
            END DO
         END DO
         !
      CASE( 5 )                                                   !==  5th order central TIM  ==! (Eq. 29)
         !
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
                  zdx2 = e1u(ji,jj) * e1u(ji,jj)
                  !!rachid       zdx2 = e1u(ji,jj) * e1t(ji,jj)
                  zdx4 = zdx2 * zdx2
                  pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (            (                   pt  (ji+1,jj,jl) + pt  (ji,jj,jl)       &
                     &                                                     -          zcu * ( pt  (ji+1,jj,jl) - pt  (ji,jj,jl) ) )   &
                     &        + z1_6   * zdx2 * ( zcu*zcu - 1._wp ) *     (                   ztu2(ji+1,jj,jl) + ztu2(ji,jj,jl)       &
                     &                                                     - 0.5_wp * zcu * ( ztu2(ji+1,jj,jl) - ztu2(ji,jj,jl) ) )   &
                     &        + z1_120 * zdx4 * ( zcu*zcu - 1._wp ) * ( zcu*zcu - 4._wp ) * ( ztu4(ji+1,jj,jl) + ztu4(ji,jj,jl)       &
                     &                                               - SIGN( 1._wp, zcu ) * ( ztu4(ji+1,jj,jl) - ztu4(ji,jj,jl) ) ) )
               END DO
            END DO
         END DO
         !
      END SELECT
      !                                                     !-- High order flux in i-direction  --!
      IF( ll_neg ) THEN
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( pt_u(ji,jj,jl) < 0._wp ) THEN
                     pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (                           pt(ji+1,jj,jl) + pt(ji,jj,jl)   &
                        &                                    - SIGN( 1._wp, pu(ji,jj) ) * ( pt(ji+1,jj,jl) - pt(ji,jj,jl) ) )
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      DO jl = 1, jpl
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               IF( ll_clem ) THEN
                  pfu_ho(ji,jj,jl) = pu(ji,jj) * pt_u(ji,jj,jl)
               ELSE
                  pfu_ho(ji,jj,jl) = puc(ji,jj,jl) * pt_u(ji,jj,jl)
               ENDIF
            END DO
         END DO
      END DO
      !
   END SUBROUTINE ultimate_x
   
 
   SUBROUTINE ultimate_y( kn_umx, pdt, pt, pv, pvc, pt_v, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ultimate_y  ***
      !!     
      !! **  Purpose :   compute  
      !!
      !! **  Method  :   ... ???
      !!                 TIM = transient interpolation Modeling 
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kn_umx    ! order of the scheme (1-5=UM or 20=CEN2)
      REAL(wp)                    , INTENT(in   ) ::   pdt       ! tracer time-step
      REAL(wp), DIMENSION(:,:  ), INTENT(in   ) ::   pv        ! ice j-velocity component
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pvc       ! ice j-velocity*A component
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pt_v      ! tracer at v-point 
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfv_ho    ! high order flux 
      !
      INTEGER  ::   ji, jj, jl       ! dummy loop indices
      REAL(wp) ::   zcv, zdy2, zdy4    !   -      -
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   ztv1, ztv2, ztv3, ztv4
      !!----------------------------------------------------------------------
      !
      !                                                     !--  Laplacian in j-direction  --!
      DO jl = 1, jpl
         DO jj = 1, jpjm1         ! First derivative (gradient)
            DO ji = fs_2, fs_jpim1
               ztv1(ji,jj,jl) = ( pt(ji,jj+1,jl) - pt(ji,jj,jl) ) * r1_e2v(ji,jj) * vmask(ji,jj,1)
            END DO
         END DO
         DO jj = 2, jpjm1         ! Second derivative (Laplacian)
            DO ji = fs_2, fs_jpim1
               ztv2(ji,jj,jl) = ( ztv1(ji,jj,jl) - ztv1(ji,jj-1,jl) ) * r1_e2t(ji,jj)
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', ztv2, 'T', 1. )
      !
      !                                                     !--  BiLaplacian in j-direction  --!
      DO jl = 1, jpl
         DO jj = 1, jpjm1         ! First derivative
            DO ji = fs_2, fs_jpim1
               ztv3(ji,jj,jl) = ( ztv2(ji,jj+1,jl) - ztv2(ji,jj,jl) ) * r1_e2v(ji,jj) * vmask(ji,jj,1)
            END DO
         END DO
         DO jj = 2, jpjm1         ! Second derivative
            DO ji = fs_2, fs_jpim1
               ztv4(ji,jj,jl) = ( ztv3(ji,jj,jl) - ztv3(ji,jj-1,jl) ) * r1_e2t(ji,jj)
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', ztv4, 'T', 1. )
      !
      !
      SELECT CASE (kn_umx )
         !
      CASE( 1 )                                                !==  1st order central TIM  ==! (Eq. 21)
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (                          ( pt(ji,jj+1,jl) + pt(ji,jj,jl) )  &
                     &                                     - SIGN( 1._wp, pv(ji,jj) ) * ( pt(ji,jj+1,jl) - pt(ji,jj,jl) ) )
               END DO
            END DO
         END DO
         !
      CASE( 2 )                                                !==  2nd order central TIM  ==! (Eq. 23)
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
                  pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (     ( pt(ji,jj+1,jl) + pt(ji,jj,jl) )  &
                     &                                     - zcv * ( pt(ji,jj+1,jl) - pt(ji,jj,jl) ) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'icedyn_adv_umx', pt_v, 'V',  1. )
         !
      CASE( 3 )                                                !==  3rd order central TIM  ==! (Eq. 24)
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
                  zdy2 = e2v(ji,jj) * e2v(ji,jj)
!!rachid       zdy2 = e2v(ji,jj) * e2t(ji,jj)
                  pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (                              ( pt  (ji,jj+1,jl) + pt  (ji,jj,jl)       &
                     &                                     -                        zcv   * ( pt  (ji,jj+1,jl) - pt  (ji,jj,jl) ) )   &
                     &        + z1_6 * zdy2 * ( zcv*zcv - 1._wp ) * (                         ztv2(ji,jj+1,jl) + ztv2(ji,jj,jl)       &
                     &                                               - SIGN( 1._wp, zcv ) * ( ztv2(ji,jj+1,jl) - ztv2(ji,jj,jl) ) ) )
               END DO
            END DO
         END DO
         !
      CASE( 4 )                                                !==  4th order central TIM  ==! (Eq. 27)
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
                  zdy2 = e2v(ji,jj) * e2v(ji,jj)
!!rachid       zdy2 = e2v(ji,jj) * e2t(ji,jj)
                  pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (                        ( pt  (ji,jj+1,jl) + pt  (ji,jj,jl)       &
                     &                                               -          zcv * ( pt  (ji,jj+1,jl) - pt  (ji,jj,jl) ) )   &
                     &        + z1_6 * zdy2 * ( zcv*zcv - 1._wp ) * (                   ztv2(ji,jj+1,jl) + ztv2(ji,jj,jl)       &
                     &                                               - 0.5_wp * zcv * ( ztv2(ji,jj+1,jl) - ztv2(ji,jj,jl) ) ) )
               END DO
            END DO
         END DO
         !
      CASE( 5 )                                                !==  5th order central TIM  ==! (Eq. 29)
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
                  zdy2 = e2v(ji,jj) * e2v(ji,jj)
!!rachid       zdy2 = e2v(ji,jj) * e2t(ji,jj)
                  zdy4 = zdy2 * zdy2
                  pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (                              ( pt  (ji,jj+1,jl) + pt  (ji,jj,jl)      &
                     &                                                     -          zcv * ( pt  (ji,jj+1,jl) - pt  (ji,jj,jl) ) )  &
                     &        + z1_6   * zdy2 * ( zcv*zcv - 1._wp ) *     (                   ztv2(ji,jj+1,jl) + ztv2(ji,jj,jl)      &
                     &                                                     - 0.5_wp * zcv * ( ztv2(ji,jj+1,jl) - ztv2(ji,jj,jl) ) )  &
                     &        + z1_120 * zdy4 * ( zcv*zcv - 1._wp ) * ( zcv*zcv - 4._wp ) * ( ztv4(ji,jj+1,jl) + ztv4(ji,jj,jl)      &
                     &                                               - SIGN( 1._wp, zcv ) * ( ztv4(ji,jj+1,jl) - ztv4(ji,jj,jl) ) ) )
               END DO
            END DO
         END DO
         !
      END SELECT
      !                                                     !-- High order flux in j-direction  --!
      IF( ll_neg ) THEN
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( pt_v(ji,jj,jl) < 0._wp ) THEN
                     pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (                          ( pt(ji,jj+1,jl) + pt(ji,jj,jl) )  &
                        &                                     - SIGN( 1._wp, pv(ji,jj) ) * ( pt(ji,jj+1,jl) - pt(ji,jj,jl) ) )
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      DO jl = 1, jpl
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               IF( ll_clem ) THEN
                  pfv_ho(ji,jj,jl) = pv(ji,jj) * pt_v(ji,jj,jl)
               ELSE
                  pfv_ho(ji,jj,jl) = pvc(ji,jj,jl) * pt_v(ji,jj,jl)
               ENDIF
            END DO
         END DO
      END DO
      !
   END SUBROUTINE ultimate_y
     

   SUBROUTINE nonosc_2d( pamsk, pdt, pu, puc, pv, pvc, ptc, pt, pt_low, pfu_low, pfv_low, pfu_ho, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc  ***
      !!     
       !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
      !!       scheme and the before field by a nonoscillatory algorithm 
      !!
      !! **  Method  :   ... ???
      !!       warning : pt and pt_low must be masked, but the boundaries
      !!       conditions on the fluxes are not necessary zalezak (1979)
      !!       drange (1995) multi-dimensional forward-in-time and upstream-
      !!       in-space based differencing for fluid
      !!----------------------------------------------------------------------
      REAL(wp)                     , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      REAL(wp)                     , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION (:,:  ), INTENT(in   ) ::   pu               ! ice i-velocity => u*e2
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   puc              ! ice i-velocity *A => u*e2*a
      REAL(wp), DIMENSION (:,:  ), INTENT(in   ) ::   pv               ! ice j-velocity => v*e1
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   pvc              ! ice j-velocity *A => v*e1*a
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   ptc, pt, pt_low  ! before field & upstream guess of after field
      REAL(wp), DIMENSION (:,:,:), INTENT(inout) ::   pfv_low, pfu_low ! upstream flux
      REAL(wp), DIMENSION (:,:,:), INTENT(inout) ::   pfv_ho, pfu_ho   ! monotonic flux
      !
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp) ::   zpos, zneg, zbig, zsml, z1_dt, zpos2, zneg2   ! local scalars
      REAL(wp) ::   zau, zbu, zcu, zav, zbv, zcv, zup, zdo, zsign, zcoef        !   -      -
      REAL(wp), DIMENSION(jpi,jpj    ) :: zbup, zbdo
      REAL(wp), DIMENSION(jpi,jpj,jpl) :: zbetup, zbetdo, zti_low, ztj_low, zzt
      !!----------------------------------------------------------------------
      zbig = 1.e+40_wp
      zsml = epsi20
      
      IF( ll_zeroup2 ) THEN
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  IF( amaxu(ji,jj,jl) == 0._wp )   pfu_ho(ji,jj,jl) = 0._wp
                  IF( amaxv(ji,jj,jl) == 0._wp )   pfv_ho(ji,jj,jl) = 0._wp
               END DO
            END DO
         END DO
      ENDIF

      IF( ll_zeroup4 ) THEN
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  IF( pfu_low(ji,jj,jl) == 0._wp )   pfu_ho(ji,jj,jl) = 0._wp
                  IF( pfv_low(ji,jj,jl) == 0._wp )   pfv_ho(ji,jj,jl) = 0._wp
               END DO
            END DO
         END DO
      ENDIF


      IF( ll_zeroup1 ) THEN
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  IF( ll_gurvan ) THEN
                     zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                        &                           - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                  ELSE
                     zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                        &                           - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                        &                     + pt(ji,jj,jl) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                        &                     + pt(ji,jj,jl) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                        &         ) * tmask(ji,jj,1)
                  ENDIF
                  IF( zzt(ji,jj,jl) < 0._wp ) THEN
                     pfu_ho(ji,jj,jl)   = pfu_low(ji,jj,jl)
                     pfv_ho(ji,jj,jl)   = pfv_low(ji,jj,jl)
                     WRITE(numout,*) '*** 1 negative high order zzt ***',ji,jj,zzt(ji,jj,jl)
                  ENDIF
                  !!               IF( ji==26 .AND. jj==86) THEN
                  !!                  WRITE(numout,*) 'zzt high order',zzt(ji,jj)
                  !!                  WRITE(numout,*) 'pfu_ho',(pfu_ho(ji,jj,jl)) * r1_e1e2t(ji,jj) * pdt
                  !!                  WRITE(numout,*) 'pfv_ho',(pfv_ho(ji,jj,jl)) * r1_e1e2t(ji,jj) * pdt
                  !!                  WRITE(numout,*) 'pfu_hom1',(pfu_ho(ji-1,jj,jl)) * r1_e1e2t(ji,jj) * pdt
                  !!                  WRITE(numout,*) 'pfv_hom1',(pfv_ho(ji,jj-1,jl)) * r1_e1e2t(ji,jj) * pdt
                  !!               ENDIF
                  IF( ll_gurvan ) THEN
                     zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                        &                           - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                  ELSE
                     zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                        &                           - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                        &                     + pt(ji,jj,jl) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                        &                     + pt(ji,jj,jl) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                        &         ) * tmask(ji,jj,1)
                  ENDIF
                  IF( zzt(ji,jj,jl) < 0._wp ) THEN
                     pfu_ho(ji-1,jj,jl) = pfu_low(ji-1,jj,jl)
                     pfv_ho(ji,jj-1,jl) = pfv_low(ji,jj-1,jl)
                     WRITE(numout,*) '*** 2 negative high order zzt ***',ji,jj,zzt(ji,jj,jl)
                  ENDIF
                  IF( ll_gurvan ) THEN
                     zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                        &                           - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                  ELSE
                     zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                        &                           - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                        &                     + pt(ji,jj,jl) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                        &                     + pt(ji,jj,jl) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                        &         ) * tmask(ji,jj,1)
                  ENDIF
                  IF( zzt(ji,jj,jl) < 0._wp ) THEN
                     WRITE(numout,*) '*** 3 negative high order zzt ***',ji,jj,zzt(ji,jj,jl)
                  ENDIF
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'icedyn_adv_umx', pfu_ho, 'U', -1., pfv_ho, 'V', -1. )
      ENDIF

      
      ! antidiffusive flux : high order minus low order
      ! --------------------------------------------------
      DO jl = 1, jpl
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               pfu_ho(ji,jj,jl) = pfu_ho(ji,jj,jl) - pfu_low(ji,jj,jl)
               pfv_ho(ji,jj,jl) = pfv_ho(ji,jj,jl) - pfv_low(ji,jj,jl)
            END DO
         END DO
      END DO

      ! extreme case where pfu_ho has to be zero
      ! ----------------------------------------
      !                                    pfu_ho
      !                           *         --->
      !                        |      |  *   |        | 
      !                        |      |      |    *   |    
      !                        |      |      |        |    *
      !            t_low :       i-1     i       i+1       i+2   
      IF( ll_prelimiter_zalesak ) THEN
         
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1 
                  zti_low(ji,jj,jl)= pt_low(ji+1,jj  ,jl)
                  ztj_low(ji,jj,jl)= pt_low(ji  ,jj+1,jl)
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'icedyn_adv_umx', zti_low, 'T', 1., ztj_low, 'T', 1. )

         !! this does not work ??
         !!            DO jj = 2, jpjm1
         !!               DO ji = fs_2, fs_jpim1
         !!                  IF( SIGN( 1., pfu_ho(ji,jj) ) /= SIGN( 1., pt_low (ji+1,jj  ) - pt_low (ji  ,jj) ) .AND.     &
         !!               &      SIGN( 1., pfv_ho(ji,jj) ) /= SIGN( 1., pt_low (ji  ,jj+1) - pt_low (ji  ,jj) )           &
         !!               &    ) THEN
         !!                     IF( SIGN( 1., pfu_ho(ji,jj) ) /= SIGN( 1., zti_low(ji+1,jj ) - zti_low(ji  ,jj) ) .AND.   &
         !!               &         SIGN( 1., pfv_ho(ji,jj) ) /= SIGN( 1., ztj_low(ji,jj+1 ) - ztj_low(ji  ,jj) )         &
         !!               &       ) THEN
         !!                        pfu_ho(ji,jj) = 0.  ;   pfv_ho(ji,jj) = 0.
         !!                     ENDIF
         !!                     IF( SIGN( 1., pfu_ho(ji,jj) ) /= SIGN( 1., pt_low (ji  ,jj) - pt_low (ji-1,jj  ) ) .AND.  &
         !!               &         SIGN( 1., pfv_ho(ji,jj) ) /= SIGN( 1., pt_low (ji  ,jj) - pt_low (ji  ,jj-1) )        &
         !!               &       )  THEN
         !!                        pfu_ho(ji,jj) = 0.  ;   pfv_ho(ji,jj) = 0.
         !!                     ENDIF
         !!                  ENDIF
         !!                END DO
         !!            END DO            

         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  IF ( pfu_ho(ji,jj,jl) * ( pt_low(ji+1,jj,jl) - pt_low(ji,jj,jl) ) <= 0. .AND.  &
                     & pfv_ho(ji,jj,jl) * ( pt_low(ji,jj+1,jl) - pt_low(ji,jj,jl) ) <= 0. ) THEN
                     !
                     IF(  pfu_ho(ji,jj,jl) * ( zti_low(ji+1,jj,jl) - zti_low(ji,jj,jl) ) <= 0 .AND.  &
                        & pfv_ho(ji,jj,jl) * ( ztj_low(ji,jj+1,jl) - ztj_low(ji,jj,jl) ) <= 0)  pfu_ho(ji,jj,jl)=0. ; pfv_ho(ji,jj,jl)=0.
                     !
                     IF(  pfu_ho(ji,jj,jl) * ( pt_low(ji  ,jj,jl) - pt_low(ji-1,jj,jl) ) <= 0 .AND.  &
                        & pfv_ho(ji,jj,jl) * ( pt_low(ji  ,jj,jl) - pt_low(ji,jj-1,jl) ) <= 0)  pfu_ho(ji,jj,jl)=0. ; pfv_ho(ji,jj,jl)=0.
                     !
                  ENDIF
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'icedyn_adv_umx', pfu_ho, 'U', -1., pfv_ho, 'V', -1. )   ! lateral boundary cond.

      ELSEIF( ll_prelimiter_devore ) THEN
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1 
                  zti_low(ji,jj,jl)= pt_low(ji+1,jj  ,jl)
                  ztj_low(ji,jj,jl)= pt_low(ji  ,jj+1,jl)
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'icedyn_adv_umx', zti_low, 'T', 1., ztj_low, 'T', 1. )

         z1_dt = 1._wp / pdt
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  zsign = SIGN( 1., pt_low(ji+1,jj,jl) - pt_low(ji,jj,jl) )
                  pfu_ho(ji,jj,jl) =  zsign * MAX( 0. , MIN( ABS(pfu_ho(ji,jj,jl)) , &
                     &                          zsign * ( pt_low (ji  ,jj,jl) - pt_low (ji-1,jj,jl) ) * e1e2t(ji  ,jj) * z1_dt , &
                     &                          zsign * ( zti_low(ji+1,jj,jl) - zti_low(ji  ,jj,jl) ) * e1e2t(ji+1,jj) * z1_dt ) )

                  zsign = SIGN( 1., pt_low(ji,jj+1,jl) - pt_low(ji,jj,jl) )
                  pfv_ho(ji,jj,jl) =  zsign * MAX( 0. , MIN( ABS(pfv_ho(ji,jj,jl)) , &
                     &                          zsign * ( pt_low (ji,jj  ,jl) - pt_low (ji,jj-1,jl) ) * e1e2t(ji,jj  ) * z1_dt , &
                     &                          zsign * ( ztj_low(ji,jj+1,jl) - ztj_low(ji,jj  ,jl) ) * e1e2t(ji,jj+1) * z1_dt ) )
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'icedyn_adv_umx', pfu_ho, 'U', -1., pfv_ho, 'V', -1. )   ! lateral boundary cond.

      ENDIF


      ! Search local extrema
      ! --------------------
      ! max/min of pt & pt_low with large negative/positive value (-/+zbig) outside ice cover
      z1_dt = 1._wp / pdt
      DO jl = 1, jpl
         
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF    ( pt(ji,jj,jl) <= 0._wp .AND. pt_low(ji,jj,jl) <= 0._wp ) THEN
                  zbup(ji,jj) = -zbig
                  zbdo(ji,jj) =  zbig
               ELSEIF( pt(ji,jj,jl) <= 0._wp .AND. pt_low(ji,jj,jl) > 0._wp ) THEN
                  zbup(ji,jj) = pt_low(ji,jj,jl)
                  zbdo(ji,jj) = pt_low(ji,jj,jl)
               ELSEIF( pt(ji,jj,jl) > 0._wp .AND. pt_low(ji,jj,jl) <= 0._wp ) THEN
                  zbup(ji,jj) = pt(ji,jj,jl)
                  zbdo(ji,jj) = pt(ji,jj,jl)
               ELSE
                  zbup(ji,jj) = MAX( pt(ji,jj,jl) , pt_low(ji,jj,jl) )
                  zbdo(ji,jj) = MIN( pt(ji,jj,jl) , pt_low(ji,jj,jl) )
               ENDIF
            END DO
         END DO

         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               !
               IF( .NOT. ll_9points ) THEN
                  zup  = MAX( zbup(ji,jj), zbup(ji-1,jj  ), zbup(ji+1,jj  ), zbup(ji  ,jj-1), zbup(ji  ,jj+1) )  ! search max/min in neighbourhood
                  zdo  = MIN( zbdo(ji,jj), zbdo(ji-1,jj  ), zbdo(ji+1,jj  ), zbdo(ji  ,jj-1), zbdo(ji  ,jj+1) )
                  !
               ELSE
                  zup  = MAX( zbup(ji,jj), zbup(ji-1,jj  ), zbup(ji+1,jj  ), zbup(ji  ,jj-1), zbup(ji  ,jj+1), &  ! search max/min in neighbourhood
                     &                     zbup(ji-1,jj-1), zbup(ji+1,jj+1), zbup(ji+1,jj-1), zbup(ji-1,jj+1)  )
                  zdo  = MIN( zbdo(ji,jj), zbdo(ji-1,jj  ), zbdo(ji+1,jj  ), zbdo(ji  ,jj-1), zbdo(ji  ,jj+1), &
                     &                     zbdo(ji-1,jj-1), zbdo(ji+1,jj+1), zbdo(ji+1,jj-1), zbdo(ji-1,jj+1)  )
               ENDIF
               !
               zpos = MAX( 0., pfu_ho(ji-1,jj,jl) ) - MIN( 0., pfu_ho(ji  ,jj,jl) ) &  ! positive/negative part of the flux
                  & + MAX( 0., pfv_ho(ji,jj-1,jl) ) - MIN( 0., pfv_ho(ji,jj  ,jl) )
               zneg = MAX( 0., pfu_ho(ji  ,jj,jl) ) - MIN( 0., pfu_ho(ji-1,jj,jl) ) &
                  & + MAX( 0., pfv_ho(ji,jj  ,jl) ) - MIN( 0., pfv_ho(ji,jj-1,jl) )
               !
               IF( ll_HgradU .AND. .NOT.ll_gurvan ) THEN
                  zneg2 = (   pt(ji,jj,jl) * MAX( 0., pu(ji,jj) - pu(ji-1,jj) ) + pt(ji,jj,jl) * MAX( 0., pv(ji,jj) - pv(ji,jj-1) ) &
                     &    ) * ( 1. - pamsk )
                  zpos2 = ( - pt(ji,jj,jl) * MIN( 0., pu(ji,jj) - pu(ji-1,jj) ) - pt(ji,jj,jl) * MIN( 0., pv(ji,jj) - pv(ji,jj-1) ) &
                     &    ) * ( 1. - pamsk )
               ELSE
                  zneg2 = 0. ; zpos2 = 0.
               ENDIF
               !
               !                                  ! up & down beta terms
               IF( (zpos+zpos2) > 0. ) THEN ; zbetup(ji,jj,jl) = MAX( 0._wp, zup - pt_low(ji,jj,jl) ) / (zpos+zpos2) * e1e2t(ji,jj) * z1_dt
               ELSE                         ; zbetup(ji,jj,jl) = 0. ! zbig
               ENDIF
               !
               IF( (zneg+zneg2) > 0. ) THEN ; zbetdo(ji,jj,jl) = MAX( 0._wp, pt_low(ji,jj,jl) - zdo ) / (zneg+zneg2) * e1e2t(ji,jj) * z1_dt
               ELSE                         ; zbetdo(ji,jj,jl) = 0. ! zbig
               ENDIF
               !
               ! if all the points are outside ice cover
               IF( zup == -zbig )   zbetup(ji,jj,jl) = 0. ! zbig
               IF( zdo ==  zbig )   zbetdo(ji,jj,jl) = 0. ! zbig            
               !

!!            IF( ji==26 .AND. jj==86) THEN
!               WRITE(numout,*) '-----------------'
!               WRITE(numout,*) 'zpos',zpos,zpos2
!               WRITE(numout,*) 'zneg',zneg,zneg2
!               WRITE(numout,*) 'puc/pu',ABS(puc(ji,jj))/MAX(epsi20, ABS(pu(ji,jj)))
!               WRITE(numout,*) 'pvc/pv',ABS(pvc(ji,jj))/MAX(epsi20, ABS(pv(ji,jj)))
!               WRITE(numout,*) 'pucm1/pu',ABS(puc(ji-1,jj))/MAX(epsi20, ABS(pu(ji-1,jj)))
!               WRITE(numout,*) 'pvcm1/pv',ABS(pvc(ji,jj-1))/MAX(epsi20, ABS(pv(ji,jj-1)))
!               WRITE(numout,*) 'pfu_ho',(pfu_ho(ji,jj)+pfu_low(ji,jj)) * r1_e1e2t(ji,jj) * pdt
!               WRITE(numout,*) 'pfv_ho',(pfv_ho(ji,jj)+pfv_low(ji,jj)) * r1_e1e2t(ji,jj) * pdt
!               WRITE(numout,*) 'pfu_hom1',(pfu_ho(ji-1,jj)+pfu_low(ji-1,jj)) * r1_e1e2t(ji,jj) * pdt
!               WRITE(numout,*) 'pfv_hom1',(pfv_ho(ji,jj-1)+pfv_low(ji,jj-1)) * r1_e1e2t(ji,jj) * pdt
!               WRITE(numout,*) 'pfu_low',pfu_low(ji,jj) * r1_e1e2t(ji,jj) * pdt
!               WRITE(numout,*) 'pfv_low',pfv_low(ji,jj) * r1_e1e2t(ji,jj) * pdt
!               WRITE(numout,*) 'pfu_lowm1',pfu_low(ji-1,jj) * r1_e1e2t(ji,jj) * pdt
!               WRITE(numout,*) 'pfv_lowm1',pfv_low(ji,jj-1) * r1_e1e2t(ji,jj) * pdt
!               
!               WRITE(numout,*) 'pt',pt(ji,jj)
!               WRITE(numout,*) 'ptim1',pt(ji-1,jj)
!               WRITE(numout,*) 'ptjm1',pt(ji,jj-1)
!               WRITE(numout,*) 'pt_low',pt_low(ji,jj)
!               WRITE(numout,*) 'zbetup',zbetup(ji,jj)
!               WRITE(numout,*) 'zbetdo',zbetdo(ji,jj)
!               WRITE(numout,*) 'zup',zup
!               WRITE(numout,*) 'zdo',zdo
!            ENDIF
            !
            END DO
         END DO
      END DO
      CALL lbc_lnk_multi( 'icedyn_adv_umx', zbetup, 'T', 1., zbetdo, 'T', 1. )   ! lateral boundary cond. (unchanged sign)

      
      ! monotonic flux in the y direction
      ! ---------------------------------
      DO jl = 1, jpl
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zau = MIN( 1._wp , zbetdo(ji,jj,jl) , zbetup(ji+1,jj,jl) )
               zbu = MIN( 1._wp , zbetup(ji,jj,jl) , zbetdo(ji+1,jj,jl) )
               zcu = 0.5  + SIGN( 0.5 , pfu_ho(ji,jj,jl) )
               !
               zcoef = ( zcu * zau + ( 1._wp - zcu ) * zbu )

               pfu_ho(ji,jj,jl) = pfu_ho(ji,jj,jl) * zcoef + pfu_low(ji,jj,jl)

!!            IF( ji==26 .AND. jj==86) THEN
!!               WRITE(numout,*) 'coefU',zcoef
!!               WRITE(numout,*) 'pfu_ho',(pfu_ho(ji,jj,jl)) * r1_e1e2t(ji,jj) * pdt
!!               WRITE(numout,*) 'pfu_hom1',(pfu_ho(ji-1,jj,jl)) * r1_e1e2t(ji,jj) * pdt
!!            ENDIF

            END DO
         END DO

         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zav = MIN( 1._wp , zbetdo(ji,jj,jl) , zbetup(ji,jj+1,jl) )
               zbv = MIN( 1._wp , zbetup(ji,jj,jl) , zbetdo(ji,jj+1,jl) )
               zcv = 0.5  + SIGN( 0.5 , pfv_ho(ji,jj,jl) )
               !
               zcoef = ( zcv * zav + ( 1._wp - zcv ) * zbv )

               pfv_ho(ji,jj,jl) = pfv_ho(ji,jj,jl) * zcoef + pfv_low(ji,jj,jl)

!!            IF( ji==26 .AND. jj==86) THEN
!!               WRITE(numout,*) 'coefV',zcoef
!!               WRITE(numout,*) 'pfv_ho',(pfv_ho(ji,jj,jl)) * r1_e1e2t(ji,jj) * pdt
!!               WRITE(numout,*) 'pfv_hom1',(pfv_ho(ji,jj-1,jl)) * r1_e1e2t(ji,jj) * pdt
!!            ENDIF
            END DO
         END DO

         ! clem test
         DO jj = 2, jpjm1
            DO ji = 2, fs_jpim1   ! vector opt.
               IF( ll_gurvan ) THEN
                  zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                           - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
               ELSE
                  zzt(ji,jj,jl) = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                           - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                     + pt(ji,jj,jl) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                     &                     + pt(ji,jj,jl) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                     &         ) * tmask(ji,jj,1)
               ENDIF
               IF( zzt(ji,jj,jl) < -epsi20 ) THEN
                  WRITE(numout,*) 'T<0 nonosc',zzt(ji,jj,jl)
               ENDIF
            END DO
         END DO

      END DO

      !
      !
   END SUBROUTINE nonosc_2d

   SUBROUTINE limiter_x( pdt, pu, puc, pt, pfu_ho, pfu_ups )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE limiter_x  ***
      !!     
      !! **  Purpose :   compute flux limiter 
      !!----------------------------------------------------------------------
      REAL(wp)                   , INTENT(in   )           ::   pdt          ! tracer time-step
      REAL(wp), DIMENSION (:,:  ), INTENT(in   )           ::   pu           ! ice i-velocity => u*e2
      REAL(wp), DIMENSION (:,:,:), INTENT(in   )           ::   puc          ! ice i-velocity *A => u*e2*a
      REAL(wp), DIMENSION (:,:,:), INTENT(in   )           ::   pt           ! ice tracer
      REAL(wp), DIMENSION (:,:,:), INTENT(inout)           ::   pfu_ho       ! high order flux
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ), OPTIONAL ::   pfu_ups      ! upstream flux
      !
      REAL(wp) ::   Cr, Rjm, Rj, Rjp, uCFL, zpsi, zh3, zlimiter, Rr
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp), DIMENSION (jpi,jpj,jpl) ::   zslpx       ! tracer slopes 
      !!----------------------------------------------------------------------
      !
      DO jl = 1, jpl
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zslpx(ji,jj,jl) = ( pt(ji+1,jj,jl) - pt(ji,jj,jl) ) * umask(ji,jj,1)
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', zslpx, 'U', -1.)   ! lateral boundary cond.
      
      DO jl = 1, jpl
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               uCFL = pdt * ABS( pu(ji,jj) ) * r1_e1e2t(ji,jj)
               
               Rjm = zslpx(ji-1,jj,jl)
               Rj  = zslpx(ji  ,jj,jl)
               Rjp = zslpx(ji+1,jj,jl)

               IF( PRESENT(pfu_ups) ) THEN

                  IF( pu(ji,jj) > 0. ) THEN   ;   Rr = Rjm
                  ELSE                        ;   Rr = Rjp
                  ENDIF

                  zh3 = pfu_ho(ji,jj,jl) - pfu_ups(ji,jj,jl)     
                  IF( Rj > 0. ) THEN
                     zlimiter =  MAX( 0., MIN( zh3, MAX(-Rr * 0.5 * ABS(pu(ji,jj)),  &
                        &        MIN( 2. * Rr * 0.5 * ABS(pu(ji,jj)),  zh3,  1.5 * Rj * 0.5 * ABS(pu(ji,jj)) ) ) ) )
                  ELSE
                     zlimiter = -MAX( 0., MIN(-zh3, MAX( Rr * 0.5 * ABS(pu(ji,jj)),  &
                        &        MIN(-2. * Rr * 0.5 * ABS(pu(ji,jj)), -zh3, -1.5 * Rj * 0.5 * ABS(pu(ji,jj)) ) ) ) )
                  ENDIF
                  pfu_ho(ji,jj,jl) = pfu_ups(ji,jj,jl) + zlimiter

               ELSE
                  IF( Rj /= 0. ) THEN
                     IF( pu(ji,jj) > 0. ) THEN   ;   Cr = Rjm / Rj
                     ELSE                        ;   Cr = Rjp / Rj
                     ENDIF
                  ELSE
                     Cr = 0.
                     !IF( pu(ji,jj) > 0. ) THEN   ;   Cr = Rjm * 1.e20
                     !ELSE                        ;   Cr = Rjp * 1.e20
                     !ENDIF
                  ENDIF

                  ! -- superbee --
                  zpsi = MAX( 0., MAX( MIN(1.,2.*Cr), MIN(2.,Cr) ) )
                  ! -- van albada 2 --
                  !!zpsi = 2.*Cr / (Cr*Cr+1.)

                  ! -- sweby (with beta=1) --
                  !!zpsi = MAX( 0., MAX( MIN(1.,1.*Cr), MIN(1.,Cr) ) )
                  ! -- van Leer --
                  !!zpsi = ( Cr + ABS(Cr) ) / ( 1. + ABS(Cr) )
                  ! -- ospre --
                  !!zpsi = 1.5 * ( Cr*Cr + Cr ) / ( Cr*Cr + Cr + 1. )
                  ! -- koren --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( (1.+2*Cr)/3., 2. ) ) )
                  ! -- charm --
                  !IF( Cr > 0. ) THEN   ;   zpsi = Cr * (3.*Cr + 1.) / ( (Cr + 1.) * (Cr + 1.) )
                  !ELSE                 ;   zpsi = 0.
                  !ENDIF
                  ! -- van albada 1 --
                  !!zpsi = (Cr*Cr + Cr) / (Cr*Cr +1)
                  ! -- smart --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( 0.25+0.75*Cr, 4. ) ) )
                  ! -- umist --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( 0.25+0.75*Cr, MIN(0.75+0.25*Cr, 2. ) ) ) )

                  ! high order flux corrected by the limiter
                  pfu_ho(ji,jj,jl) = pfu_ho(ji,jj,jl) - ABS( pu(ji,jj) ) * ( (1.-zpsi) + uCFL*zpsi ) * Rj * 0.5

               ENDIF
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', pfu_ho, 'U', -1.)   ! lateral boundary cond.
      !
   END SUBROUTINE limiter_x

   SUBROUTINE limiter_y( pdt, pv, pvc, pt, pfv_ho, pfv_ups )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE limiter_y  ***
      !!     
      !! **  Purpose :   compute flux limiter 
      !!----------------------------------------------------------------------
      REAL(wp)                   , INTENT(in   )           ::   pdt          ! tracer time-step
      REAL(wp), DIMENSION (:,:  ), INTENT(in   )           ::   pv           ! ice i-velocity => u*e2
      REAL(wp), DIMENSION (:,:,:), INTENT(in   )           ::   pvc          ! ice i-velocity *A => u*e2*a
      REAL(wp), DIMENSION (:,:,:), INTENT(in   )           ::   pt           ! ice tracer
      REAL(wp), DIMENSION (:,:,:), INTENT(inout)           ::   pfv_ho       ! high order flux
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ), OPTIONAL ::   pfv_ups      ! upstream flux
      !
      REAL(wp) ::   Cr, Rjm, Rj, Rjp, vCFL, zpsi, zh3, zlimiter, Rr
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp), DIMENSION (jpi,jpj,jpl) ::   zslpy       ! tracer slopes 
      !!----------------------------------------------------------------------
      !
      DO jl = 1, jpl
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zslpy(ji,jj,jl) = ( pt(ji,jj+1,jl) - pt(ji,jj,jl) ) * vmask(ji,jj,1)
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', zslpy, 'V', -1.)   ! lateral boundary cond.

      DO jl = 1, jpl
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               vCFL = pdt * ABS( pv(ji,jj) ) * r1_e1e2t(ji,jj)

               Rjm = zslpy(ji,jj-1,jl)
               Rj  = zslpy(ji,jj  ,jl)
               Rjp = zslpy(ji,jj+1,jl)

               IF( PRESENT(pfv_ups) ) THEN

                  IF( pv(ji,jj) > 0. ) THEN   ;   Rr = Rjm
                  ELSE                        ;   Rr = Rjp
                  ENDIF

                  zh3 = pfv_ho(ji,jj,jl) - pfv_ups(ji,jj,jl)     
                  IF( Rj > 0. ) THEN
                     zlimiter =  MAX( 0., MIN( zh3, MAX(-Rr * 0.5 * ABS(pv(ji,jj)),  &
                        &        MIN( 2. * Rr * 0.5 * ABS(pv(ji,jj)),  zh3,  1.5 * Rj * 0.5 * ABS(pv(ji,jj)) ) ) ) )
                  ELSE
                     zlimiter = -MAX( 0., MIN(-zh3, MAX( Rr * 0.5 * ABS(pv(ji,jj)),  &
                        &        MIN(-2. * Rr * 0.5 * ABS(pv(ji,jj)), -zh3, -1.5 * Rj * 0.5 * ABS(pv(ji,jj)) ) ) ) )
                  ENDIF
                  pfv_ho(ji,jj,jl) = pfv_ups(ji,jj,jl) + zlimiter

               ELSE

                  IF( Rj /= 0. ) THEN
                     IF( pv(ji,jj) > 0. ) THEN   ;   Cr = Rjm / Rj
                     ELSE                        ;   Cr = Rjp / Rj
                     ENDIF
                  ELSE
                     Cr = 0.
                     !IF( pv(ji,jj) > 0. ) THEN   ;   Cr = Rjm * 1.e20
                     !ELSE                        ;   Cr = Rjp * 1.e20
                     !ENDIF
                  ENDIF

                  ! -- superbee --
                  zpsi = MAX( 0., MAX( MIN(1.,2.*Cr), MIN(2.,Cr) ) )
                  ! -- van albada 2 --
                  !!zpsi = 2.*Cr / (Cr*Cr+1.)

                  ! -- sweby (with beta=1) --
                  !!zpsi = MAX( 0., MAX( MIN(1.,1.*Cr), MIN(1.,Cr) ) )
                  ! -- van Leer --
                  !!zpsi = ( Cr + ABS(Cr) ) / ( 1. + ABS(Cr) )
                  ! -- ospre --
                  !!zpsi = 1.5 * ( Cr*Cr + Cr ) / ( Cr*Cr + Cr + 1. )
                  ! -- koren --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( (1.+2*Cr)/3., 2. ) ) )
                  ! -- charm --
                  !IF( Cr > 0. ) THEN   ;   zpsi = Cr * (3.*Cr + 1.) / ( (Cr + 1.) * (Cr + 1.) )
                  !ELSE                 ;   zpsi = 0.
                  !ENDIF
                  ! -- van albada 1 --
                  !!zpsi = (Cr*Cr + Cr) / (Cr*Cr +1)
                  ! -- smart --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( 0.25+0.75*Cr, 4. ) ) )
                  ! -- umist --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( 0.25+0.75*Cr, MIN(0.75+0.25*Cr, 2. ) ) ) )

                  ! high order flux corrected by the limiter
                  pfv_ho(ji,jj,jl) = pfv_ho(ji,jj,jl) - ABS( pv(ji,jj) ) * ( (1.-zpsi) + vCFL*zpsi ) * Rj * 0.5

               ENDIF
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', pfv_ho, 'V', -1.)   ! lateral boundary cond.
      !
   END SUBROUTINE limiter_y

#else
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_adv_umx
