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

   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: amaxu, amaxv
   
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
      REAL(wp) ::   zcfl , zdt
      REAL(wp), DIMENSION(jpi,jpj) ::   zudy, zvdx, zcu_box, zcv_box, zua_ho, zva_ho
      REAL(wp), DIMENSION(jpi,jpj) ::   zhvar
      REAL(wp), DIMENSION(jpi,jpj) ::   zati1, zati2, z1_ai, z1_aip
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_dyn_adv_umx: Ultimate-Macho advection scheme'
      !
      !
      ! --- If ice drift field is too fast, use an appropriate time step for advection (CFL test for stability) --- !        
      zcfl  =            MAXVAL( ABS( pu_ice(:,:) ) * rdt_ice * r1_e1u(:,:) )
      zcfl  = MAX( zcfl, MAXVAL( ABS( pv_ice(:,:) ) * rdt_ice * r1_e2v(:,:) ) )
      IF( lk_mpp )   CALL mpp_max( zcfl )

      IF( zcfl > 0.5 ) THEN   ;   icycle = 2 
      ELSE                    ;   icycle = 1 
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
         IF(.NOT. ALLOCATED(amaxu))       ALLOCATE(amaxu (jpi,jpj))
         IF(.NOT. ALLOCATED(amaxv))       ALLOCATE(amaxv (jpi,jpj))
      ENDIF
      !---------------!
      !== advection ==!
      !---------------!
      DO jt = 1, icycle

         IF( ll_ADVopw ) THEN
            zamsk = 1._wp
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zudy, zvdx, zcu_box, zcv_box, pato_i(:,:), pato_i(:,:) )   ! Open water area 
            zamsk = 0._wp
         ELSE
            zati1(:,:) = SUM( pa_i(:,:,:), dim=3 )
         ENDIF
         
         DO jl = 1, jpl
            !
            WHERE( pa_i(:,:,jl) >= epsi20 )   ;   z1_ai(:,:) = 1._wp / pa_i(:,:,jl)
            ELSEWHERE                         ;   z1_ai(:,:) = 0.
            END WHERE
            !
            WHERE( pa_ip(:,:,jl) >= epsi20 )  ;   z1_aip(:,:) = 1._wp / pa_ip(:,:,jl)
            ELSEWHERE                         ;   z1_aip(:,:) = 0.
            END WHERE
            !
            IF( ll_zeroup2 ) THEN
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     amaxu(ji,jj)=MAX( pa_i(ji,jj,jl), pa_i(ji,jj-1,jl), pa_i(ji,jj+1,jl), &
                        &                              pa_i(ji+1,jj,jl), pa_i(ji+1,jj-1,jl), pa_i(ji+1,jj+1,jl) )
                     amaxv(ji,jj)=MAX( pa_i(ji,jj,jl), pa_i(ji-1,jj,jl), pa_i(ji+1,jj,jl), &
                        &                              pa_i(ji,jj+1,jl), pa_i(ji-1,jj+1,jl), pa_i(ji+1,jj+1,jl) )
                 END DO
               END DO
               CALL lbc_lnk_multi(amaxu, 'T', 1., amaxv, 'T', 1.)
            ENDIF
            !
            zamsk = 1._wp
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zudy, zvdx, zcu_box, zcv_box, pa_i(:,:,jl), pa_i(:,:,jl), &        ! Ice area
               &          zua_ho, zva_ho )
            zamsk = 0._wp
            !
            IF( ll_thickness ) THEN
               zua_ho(:,:) = zudy(:,:)
               zva_ho(:,:) = zvdx(:,:)
            ENDIF
            !
            zhvar(:,:) = pv_i(:,:,jl) * z1_ai(:,:)
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho , zva_ho , zcu_box, zcv_box, zhvar(:,:), pv_i (:,:,jl) )    ! Ice volume
            IF( ll_thickness )   pv_i(:,:,jl) = zhvar(:,:) * pa_i(:,:,jl)
            !
            zhvar(:,:) = pv_s(:,:,jl) * z1_ai(:,:)
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho , zva_ho , zcu_box, zcv_box, zhvar(:,:), pv_s (:,:,jl) )    ! Snw volume
            IF( ll_thickness )   pv_s(:,:,jl) = zhvar(:,:) * pa_i(:,:,jl)
            !
            zhvar(:,:) = psv_i(:,:,jl) * z1_ai(:,:)
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho , zva_ho , zcu_box, zcv_box, zhvar(:,:), psv_i(:,:,jl) )    ! Salt content
            !
            zhvar(:,:) = poa_i(:,:,jl) * z1_ai(:,:)
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho , zva_ho , zcu_box, zcv_box, zhvar(:,:), poa_i(:,:,jl) )    ! Age content
            !
            DO jk = 1, nlay_i
               zhvar(:,:) = pe_i(:,:,jk,jl) * z1_ai(:,:)
               CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar(:,:), pe_i(:,:,jk,jl) ) ! Ice heat content
            END DO
            !
            DO jk = 1, nlay_s
               zhvar(:,:) = pe_s(:,:,jk,jl) * z1_ai(:,:)
               CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar(:,:), pe_s(:,:,jk,jl) ) ! Snw heat content
            END DO
            !
            IF ( ln_pnd_H12 ) THEN
               !
               zamsk = 1._wp
               CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zudy, zvdx, zcu_box, zcv_box, pa_ip(:,:,jl), pa_ip(:,:,jl), &   ! mp fraction
                  &          zua_ho, zva_ho )
               zamsk = 0._wp

               zhvar(:,:) = pv_ip(:,:,jl) * z1_ai(:,:)
               CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho , zva_ho , zcu_box, zcv_box, zhvar(:,:), pv_ip(:,:,jl) ) ! mp volume
            ENDIF
            !
            !
         END DO
         !
         IF( .NOT. ll_ADVopw ) THEN
            zati2(:,:) = SUM( pa_i(:,:,:), dim=3 )
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  pato_i(ji,jj) = pato_i(ji,jj) - ( zati2(ji,jj) - zati1(ji,jj) ) &                                                  ! Open water area
                     &                          - ( zudy(ji,jj) - zudy(ji-1,jj) + zvdx(ji,jj) - zvdx(ji,jj-1) )*r1_e1e2t(ji,jj)*zdt
               END DO
            END DO
            CALL lbc_lnk( pato_i(:,:), 'T',  1. )
         ENDIF
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
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   )           ::   pu   , pv      ! 2 ice velocity components => u*e2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   )           ::   puc  , pvc     ! 2 ice velocity components => u*e2 or u*a*e2u
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   )           ::   pubox, pvbox   ! upstream velocity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)           ::   pt             ! tracer field
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)           ::   ptc            ! tracer content field
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out), OPTIONAL ::   pua_ho, pva_ho ! high order u*a fluxes
      !
      INTEGER  ::   ji, jj           ! dummy loop indices  
      REAL(wp) ::   ztra             ! local scalar
      INTEGER  ::   kn_limiter = 1   ! 1=nonosc ; 2=superbee ; 3=h3(rachid)
      REAL(wp), DIMENSION(jpi,jpj) ::   zfu_ho , zfv_ho , zt_u, zt_v, zpt
      REAL(wp), DIMENSION(jpi,jpj) ::   zfu_ups, zfv_ups, zt_ups   ! only for nonosc 
      !!----------------------------------------------------------------------
      !
      !  upstream (_ups) advection with initial mass fluxes
      ! ---------------------------------------------------

      IF( ll_gurvan .AND. pamsk==0. ) THEN
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               pt(ji,jj) = ( pt (ji,jj) + pt(ji,jj) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) &
                  &                     + pt(ji,jj) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
            END DO
         END DO
         CALL lbc_lnk( pt, 'T', 1. )
      ENDIF

      
      IF( .NOT. ll_upsxy ) THEN

         ! fluxes in both x-y directions
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               IF( ll_clem ) THEN
                  zfu_ups(ji,jj) = MAX( pu(ji,jj), 0._wp ) * pt(ji,jj) + MIN( pu(ji,jj), 0._wp ) * pt(ji+1,jj)
                  zfv_ups(ji,jj) = MAX( pv(ji,jj), 0._wp ) * pt(ji,jj) + MIN( pv(ji,jj), 0._wp ) * pt(ji,jj+1)
               ELSE
                  zfu_ups(ji,jj) = MAX( puc(ji,jj), 0._wp ) * pt(ji,jj) + MIN( puc(ji,jj), 0._wp ) * pt(ji+1,jj)
                  zfv_ups(ji,jj) = MAX( pvc(ji,jj), 0._wp ) * pt(ji,jj) + MIN( pvc(ji,jj), 0._wp ) * pt(ji,jj+1)
               ENDIF
            END DO
         END DO

      ELSE
         !
         IF( MOD( (kt - 1) / nn_fsbc , 2 ) ==  MOD( (jt - 1) , 2 ) ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
            ! flux in x-direction
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ll_clem ) THEN
                     zfu_ups(ji,jj) = MAX( pu(ji,jj), 0._wp ) * pt(ji,jj) + MIN( pu(ji,jj), 0._wp ) * pt(ji+1,jj)
                  ELSE
                     zfu_ups(ji,jj) = MAX( puc(ji,jj), 0._wp ) * pt(ji,jj) + MIN( puc(ji,jj), 0._wp ) * pt(ji+1,jj)
                  ENDIF
               END DO
            END DO
            
            ! first guess of tracer content from u-flux
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( ll_clem ) THEN
                     IF( ll_gurvan ) THEN
                        zpt(ji,jj) = ( pt(ji,jj) - ( zfu_ups(ji,jj) - zfu_ups(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                     ELSE
                        zpt(ji,jj) = ( pt(ji,jj) - ( zfu_ups(ji,jj) - zfu_ups(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                           &                     + pt(ji,jj) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                           &         ) * tmask(ji,jj,1)
                     ENDIF
                  ELSE
                     zpt(ji,jj) = ( ptc(ji,jj) - ( zfu_ups(ji,jj) - zfu_ups(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj) ) &
                        &         * tmask(ji,jj,1)
                  ENDIF
!!                  IF( ji==26 .AND. jj==86) THEN
!!                     WRITE(numout,*) '************************'
!!                     WRITE(numout,*) 'zpt upstream',zpt(ji,jj)
!!                  ENDIF
               END DO
            END DO
            CALL lbc_lnk( zpt, 'T', 1. )
            !
            ! flux in y-direction
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ll_clem ) THEN
                     zfv_ups(ji,jj) = MAX( pv(ji,jj), 0._wp ) * zpt(ji,jj) + MIN( pv(ji,jj), 0._wp ) * zpt(ji,jj+1)
                  ELSE
                     zfv_ups(ji,jj) = MAX( pvc(ji,jj), 0._wp ) * zpt(ji,jj) + MIN( pvc(ji,jj), 0._wp ) * zpt(ji,jj+1)
                  ENDIF
               END DO
            END DO

         !
         ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
            ! flux in y-direction
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ll_clem ) THEN
                     zfv_ups(ji,jj) = MAX( pv(ji,jj), 0._wp ) * pt(ji,jj) + MIN( pv(ji,jj), 0._wp ) * pt(ji,jj+1)
                  ELSE
                     zfv_ups(ji,jj) = MAX( pvc(ji,jj), 0._wp ) * pt(ji,jj) + MIN( pvc(ji,jj), 0._wp ) * pt(ji,jj+1)
                  ENDIF
               END DO
            END DO

            ! first guess of tracer content from v-flux
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( ll_clem ) THEN
                     IF( ll_gurvan ) THEN
                        zpt(ji,jj) = ( pt(ji,jj) - ( zfv_ups(ji,jj) - zfv_ups(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                     ELSE
                        zpt(ji,jj) = ( pt(ji,jj) - ( zfv_ups(ji,jj) - zfv_ups(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) &
                        &                        + pt(ji,jj) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) ) &
                        &         * tmask(ji,jj,1)
                     ENDIF
                  ELSE
                     zpt(ji,jj) = ( ptc(ji,jj) - ( zfv_ups(ji,jj) - zfv_ups(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) ) &
                        &         * tmask(ji,jj,1)
                  ENDIF
!!                  IF( ji==26 .AND. jj==86) THEN
!!                     WRITE(numout,*) '************************'
!!                     WRITE(numout,*) 'zpt upstream',zpt(ji,jj)
!!                  ENDIF
                END DO
            END DO
            CALL lbc_lnk( zpt, 'T', 1. )
            !
            ! flux in x-direction
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ll_clem ) THEN
                     zfu_ups(ji,jj) = MAX( pu(ji,jj), 0._wp ) * zpt(ji,jj) + MIN( pu(ji,jj), 0._wp ) * zpt(ji+1,jj)
                  ELSE
                     zfu_ups(ji,jj) = MAX( puc(ji,jj), 0._wp ) * zpt(ji,jj) + MIN( puc(ji,jj), 0._wp ) * zpt(ji+1,jj)
                  ENDIF
               END DO
            END DO
            !
         ENDIF
         
      ENDIF

      IF( ll_clem .AND. kn_limiter /= 1 )  &
         & CALL ctl_stop( 'STOP', 'icedyn_adv_umx: ll_clem incompatible with limiters other than nonosc' )

      IF( ll_zeroup2 ) THEN
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               IF( amaxu(ji,jj) == 0._wp )   zfu_ups(ji,jj) = 0._wp
               IF( amaxv(ji,jj) == 0._wp )   zfv_ups(ji,jj) = 0._wp
            END DO
         END DO
      ENDIF

      ! guess after content field with upstream scheme
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            ztra          = - (   zfu_ups(ji,jj) - zfu_ups(ji-1,jj  )   &
               &                + zfv_ups(ji,jj) - zfv_ups(ji  ,jj-1) ) * r1_e1e2t(ji,jj)
            IF( ll_clem ) THEN
               IF( ll_gurvan ) THEN
                  zt_ups(ji,jj) = ( pt (ji,jj) + pdt * ztra ) * tmask(ji,jj,1)
               ELSE
                  zt_ups(ji,jj) = ( pt (ji,jj) + pdt * ztra + ( pt(ji,jj) * pdt * ( pu(ji,jj) - pu(ji-1,jj) )   &
                     &                                      +   pt(ji,jj) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) ) &
                     &                                        * r1_e1e2t(ji,jj) * (1.-pamsk) ) * tmask(ji,jj,1)
               ENDIF
            ELSE
               zt_ups(ji,jj) = ( ptc(ji,jj) + pdt * ztra ) * tmask(ji,jj,1)
            ENDIF
!!            IF( ji==26 .AND. jj==86) THEN
!!               WRITE(numout,*) '**************************'
!!               WRITE(numout,*) 'zt upstream',zt_ups(ji,jj)
!!            ENDIF
         END DO
      END DO
      CALL lbc_lnk( zt_ups, 'T', 1. )

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
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               IF( ll_gurvan ) THEN
                  ztra       = - ( zfu_ho(ji,jj) - zfu_ho(ji-1,jj) + zfv_ho(ji,jj) - zfv_ho(ji,jj-1) ) * r1_e1e2t(ji,jj) 
               ELSE
                  ztra       = ( - ( zfu_ho(ji,jj) - zfu_ho(ji-1,jj) + zfv_ho(ji,jj) - zfv_ho(ji,jj-1) )  & 
                     &           + ( pt(ji,jj) * ( pu(ji,jj) - pu(ji-1,jj) ) * (1.-pamsk) ) &
                     &           + ( pt(ji,jj) * ( pv(ji,jj) - pv(ji,jj-1) ) * (1.-pamsk) ) ) * r1_e1e2t(ji,jj)
               ENDIF
               pt(ji,jj) = ( pt(ji,jj) + pdt * ztra ) * tmask(ji,jj,1)

               IF( pt(ji,jj) < -epsi20 ) THEN
                  WRITE(numout,*) 'T<0 ',pt(ji,jj)
               ENDIF
               
               IF( pt(ji,jj) < 0._wp .AND. pt(ji,jj) >= -epsi20 )   pt(ji,jj) = 0._wp
               
!!               IF( ji==26 .AND. jj==86) THEN
!!                  WRITE(numout,*) 'zt high order',pt(ji,jj)
!!               ENDIF
            END DO
         END DO
         CALL lbc_lnk( pt, 'T',  1. )
      ENDIF
   
      ! Rachid trick
      ! ------------
      IF( ll_clem ) THEN
         IF( pamsk == 0. ) THEN
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ABS( puc(ji,jj) ) > 0._wp .AND. ABS( pu(ji,jj) ) > 0._wp ) THEN
                     zfu_ho (ji,jj) = zfu_ho (ji,jj) * puc(ji,jj) / pu(ji,jj)
                     zfu_ups(ji,jj) = zfu_ups(ji,jj) * puc(ji,jj) / pu(ji,jj)
                  ELSE
                     zfu_ho (ji,jj) = 0._wp
                     zfu_ups(ji,jj) = 0._wp
                  ENDIF
                  !
                  IF( ABS( pvc(ji,jj) ) > 0._wp .AND. ABS( pv(ji,jj) ) > 0._wp ) THEN
                     zfv_ho (ji,jj) = zfv_ho (ji,jj) * pvc(ji,jj) / pv(ji,jj)
                     zfv_ups(ji,jj) = zfv_ups(ji,jj) * pvc(ji,jj) / pv(ji,jj)
                  ELSE
                     zfv_ho (ji,jj) = 0._wp  
                     zfv_ups(ji,jj) = 0._wp  
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF( ll_zeroup5 ) THEN
         DO jj = 2, jpjm1
            DO ji = 2, fs_jpim1   ! vector opt.
               zpt(ji,jj) = ( ptc(ji,jj) - ( zfu_ho(ji,jj) - zfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                  &                      - ( zfv_ho(ji,jj) - zfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj)  ) * tmask(ji,jj,1)
               IF( zpt(ji,jj) < 0. ) THEN
                  zfu_ho(ji,jj)   = zfu_ups(ji,jj)
                  zfu_ho(ji-1,jj) = zfu_ups(ji-1,jj)
                  zfv_ho(ji,jj)   = zfv_ups(ji,jj)
                  zfv_ho(ji,jj-1) = zfv_ups(ji,jj-1)
               ENDIF
            END DO
         END DO
         CALL lbc_lnk_multi( zfu_ho, 'U',  -1., zfv_ho, 'V',  -1. )
      ENDIF

      ! output high order fluxes u*a
      ! ----------------------------
      IF( PRESENT( pua_ho ) ) THEN
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               pua_ho(ji,jj) = zfu_ho(ji,jj)
               pva_ho(ji,jj) = zfv_ho(ji,jj)
            END DO
         END DO
      ENDIF


      IF( .NOT.ll_thickness ) THEN
         ! final trend with corrected fluxes
         ! ------------------------------------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1 
               ztra = - ( zfu_ho(ji,jj) - zfu_ho(ji-1,jj) + zfv_ho(ji,jj) - zfv_ho(ji,jj-1) ) * r1_e1e2t(ji,jj) * pdt  

               ptc(ji,jj) = ( ptc(ji,jj) + ztra ) * tmask(ji,jj,1)
                             
!!               IF( ji==26 .AND. jj==86) THEN
!!                  WRITE(numout,*) 'ztc high order',ptc(ji,jj)
!!               ENDIF
               
            END DO
         END DO
         CALL lbc_lnk( ptc, 'T',  1. )
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
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pt               ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pu, pv           ! 2 ice velocity components
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   puc, pvc         ! 2 ice velocity * A components
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   ptc              ! tracer content at before time step 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pfu_ho, pfv_ho   ! high order fluxes 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pt_ups           ! upstream guess of tracer content 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pfu_ups, pfv_ups ! upstream fluxes 
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      LOGICAL  ::   ll_xy = .TRUE.
      REAL(wp), DIMENSION(jpi,jpj) ::   zzt
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ll_xy ) THEN   !-- no alternate directions --!
         !
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               IF( ll_clem ) THEN
                  pfu_ho(ji,jj) = 0.5 * pu(ji,jj) * ( pt(ji,jj) + pt(ji+1,jj) )
                  pfv_ho(ji,jj) = 0.5 * pv(ji,jj) * ( pt(ji,jj) + pt(ji,jj+1) )
               ELSE
                  pfu_ho(ji,jj) = 0.5 * puc(ji,jj) * ( pt(ji,jj) + pt(ji+1,jj) )
                  pfv_ho(ji,jj) = 0.5 * pvc(ji,jj) * ( pt(ji,jj) + pt(ji,jj+1) )
               ENDIF
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
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ll_clem ) THEN
                     pfu_ho(ji,jj) = 0.5 * pu(ji,jj) * ( pt(ji,jj) + pt(ji+1,jj) )
                  ELSE
                     pfu_ho(ji,jj) = 0.5 * puc(ji,jj) * ( pt(ji,jj) + pt(ji+1,jj) )
                  ENDIF
               END DO
            END DO
            IF( kn_limiter == 2 )   CALL limiter_x( pdt, pu, puc, pt, pfu_ho )
            IF( kn_limiter == 3 )   CALL limiter_x( pdt, pu, puc, pt, pfu_ho, pfu_ups )

            ! first guess of tracer content from u-flux
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( ll_clem ) THEN
                     zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)        &
                           &                  + pt(ji,jj) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) ) &
                           &         * tmask(ji,jj,1)
                  ELSE                     
                     zzt(ji,jj) = ( ptc(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                  ENDIF
               END DO
            END DO
            CALL lbc_lnk( zzt, 'T', 1. )

            ! flux in y-direction
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ll_clem ) THEN
                     pfv_ho(ji,jj) = 0.5 * pv(ji,jj) * ( zzt(ji,jj) + zzt(ji,jj+1) )
                  ELSE                     
                     pfv_ho(ji,jj) = 0.5 * pvc(ji,jj) * ( zzt(ji,jj) + zzt(ji,jj+1) )
                  ENDIF
               END DO
            END DO
            IF( kn_limiter == 2 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho )
            IF( kn_limiter == 3 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho, pfv_ups )

         ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
            !
            ! flux in y-direction
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ll_clem ) THEN
                     pfv_ho(ji,jj) = 0.5 * pv(ji,jj) * ( pt(ji,jj) + pt(ji,jj+1) )
                  ELSE                     
                     pfv_ho(ji,jj) = 0.5 * pvc(ji,jj) * ( pt(ji,jj) + pt(ji,jj+1) )
                  ENDIF
               END DO
            END DO
            IF( kn_limiter == 2 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho )
            IF( kn_limiter == 3 )   CALL limiter_y( pdt, pv, pvc, pt, pfv_ho, pfv_ups )
            !
            ! first guess of tracer content from v-flux
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( ll_clem ) THEN
                     zzt(ji,jj) = ( pt(ji,jj) - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) &
                        &                     + pt(ji,jj) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) ) &
                        &         * tmask(ji,jj,1)
                  ELSE
                     zzt(ji,jj) = ( ptc(ji,jj) - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                  ENDIF
               END DO
            END DO
            CALL lbc_lnk( zzt, 'T', 1. )
            !
            ! flux in x-direction
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( ll_clem ) THEN
                     pfu_ho(ji,jj) = 0.5 * pu(ji,jj) * ( zzt(ji,jj) + zzt(ji+1,jj) )
                  ELSE
                     pfu_ho(ji,jj) = 0.5 * puc(ji,jj) * ( zzt(ji,jj) + zzt(ji+1,jj) )
                  ENDIF
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
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pt               ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pu, pv           ! 2 ice velocity components
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   puc, pvc         ! 2 ice velocity * A components
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pubox, pvbox     ! upstream velocity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   ptc              ! tracer content at before time step 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pt_u, pt_v       ! tracer at u- and v-points 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pfu_ho, pfv_ho   ! high order fluxes 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pt_ups           ! upstream guess of tracer content 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pfu_ups, pfv_ups ! upstream fluxes 
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) ::   ztra
      REAL(wp), DIMENSION(jpi,jpj) ::   zzt, zzfu_ho, zzfv_ho
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
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( ll_clem ) THEN
                     IF( ll_gurvan ) THEN
                        zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                     ELSE
                        zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                           &                     + pt(ji,jj) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                           &         ) * tmask(ji,jj,1)
                     ENDIF
                  ELSE
                     zzt(ji,jj) = ( ptc(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                  ENDIF
               END DO
            END DO
            CALL lbc_lnk( zzt, 'T', 1. )

         ELSE

            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( ll_gurvan ) THEN
                     zzt(ji,jj) = pt(ji,jj) - pubox(ji,jj) * pdt * ( pt_u(ji,jj) - pt_u(ji-1,jj) ) * r1_e1t(ji,jj)  &
                        &                   - pt   (ji,jj) * pdt * ( pu  (ji,jj) - pu  (ji-1,jj) ) * r1_e1e2t(ji,jj)
                  ELSE
                     zzt(ji,jj) = pt(ji,jj) - pubox(ji,jj) * pdt * ( pt_u(ji,jj) - pt_u(ji-1,jj) ) * r1_e1t(ji,jj)  &
                        &                   - pt   (ji,jj) * pdt * ( pu  (ji,jj) - pu  (ji-1,jj) ) * r1_e1e2t(ji,jj) * pamsk
                  ENDIF
                  zzt(ji,jj) = zzt(ji,jj) * tmask(ji,jj,1)
               END DO
            END DO
            CALL lbc_lnk( zzt, 'T', 1. )
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
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( ll_clem ) THEN
                     IF( ll_gurvan ) THEN
                        zzt(ji,jj) = ( pt(ji,jj) - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                     ELSE
                        zzt(ji,jj) = ( pt(ji,jj) - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) &
                           &                     + pt(ji,jj) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                           &         ) * tmask(ji,jj,1)
                     ENDIF
                  ELSE
                     zzt(ji,jj) = ( ptc(ji,jj) - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) ) &
                        &         * tmask(ji,jj,1)
                  ENDIF
                END DO
            END DO
            CALL lbc_lnk( zzt, 'T', 1. )
            
         ELSE
            
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  IF( ll_gurvan ) THEN
                     zzt(ji,jj) = pt(ji,jj) - pvbox(ji,jj) * pdt * ( pt_v(ji,jj) - pt_v(ji,jj-1) ) * r1_e2t(ji,jj)  &
                        &                   - pt   (ji,jj) * pdt * ( pv  (ji,jj) - pv  (ji,jj-1) ) * r1_e1e2t(ji,jj)
                  ELSE
                     zzt(ji,jj) = pt(ji,jj) - pvbox(ji,jj) * pdt * ( pt_v(ji,jj) - pt_v(ji,jj-1) ) * r1_e2t(ji,jj)  &
                        &                   - pt   (ji,jj) * pdt * ( pv  (ji,jj) - pv  (ji,jj-1) ) * r1_e1e2t(ji,jj) * pamsk
                  ENDIF
                  zzt(ji,jj) = zzt(ji,jj) * tmask(ji,jj,1)
               END DO
            END DO
            CALL lbc_lnk( zzt, 'T', 1. )
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
            zzfu_ho(:,:) = pfu_ho(:,:)
            zzfv_ho(:,:) = pfv_ho(:,:)
            ! 1st iteration of nonosc (limit the flux with the upstream solution)
            IF( ll_clem ) THEN
               CALL nonosc_2d ( pamsk, pdt, pu, puc, pv, pvc, ptc, pt, pt_ups, pfu_ups, pfv_ups, zzfu_ho, zzfv_ho )
            ELSE
               CALL nonosc_2d ( pamsk, pdt, pu, puc, pv, pvc, ptc, ptc, pt_ups, pfu_ups, pfv_ups, zzfu_ho, zzfv_ho )
            ENDIF
            ! guess after content field with high order
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  ztra       = - ( zzfu_ho(ji,jj) - zzfu_ho(ji-1,jj) + zzfv_ho(ji,jj) - zzfv_ho(ji,jj-1) ) * r1_e1e2t(ji,jj)
                  zzt(ji,jj) =   ( ptc(ji,jj) + pdt * ztra ) * tmask(ji,jj,1)
               END DO
            END DO
            CALL lbc_lnk( zzt, 'T', 1. )
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
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pu        ! ice i-velocity component
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   puc       ! ice i-velocity * A component
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pt_u      ! tracer at u-point 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pfu_ho    ! high order flux 
      !
      INTEGER  ::   ji, jj             ! dummy loop indices
      REAL(wp) ::   zcu, zdx2, zdx4    !   -      -
      REAL(wp), DIMENSION(jpi,jpj) ::   ztu1, ztu2, ztu3, ztu4
      !!----------------------------------------------------------------------
      !
      !                                                     !--  Laplacian in i-direction  --!
      DO jj = 2, jpjm1         ! First derivative (gradient)
         DO ji = 1, fs_jpim1
            ztu1(ji,jj) = ( pt(ji+1,jj) - pt(ji,jj) ) * r1_e1u(ji,jj) * umask(ji,jj,1)
         END DO
         !                     ! Second derivative (Laplacian)
         DO ji = fs_2, fs_jpim1
            ztu2(ji,jj) = ( ztu1(ji,jj) - ztu1(ji-1,jj) ) * r1_e1t(ji,jj)
         END DO
      END DO
      CALL lbc_lnk( ztu2, 'T', 1. )
      !
      !                                                     !--  BiLaplacian in i-direction  --!
      DO jj = 2, jpjm1         ! Third derivative
         DO ji = 1, fs_jpim1
            ztu3(ji,jj) = ( ztu2(ji+1,jj) - ztu2(ji,jj) ) * r1_e1u(ji,jj) * umask(ji,jj,1)
         END DO
         !                     ! Fourth derivative
         DO ji = fs_2, fs_jpim1
            ztu4(ji,jj) = ( ztu3(ji,jj) - ztu3(ji-1,jj) ) * r1_e1t(ji,jj)
         END DO
      END DO
      CALL lbc_lnk( ztu4, 'T', 1. )
      !
      !
      SELECT CASE (kn_umx )
      !
      CASE( 1 )                                                   !==  1st order central TIM  ==! (Eq. 21)
         !        
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (                              pt(ji+1,jj) + pt(ji,jj)   &
                  &                                    - SIGN( 1._wp, pu(ji,jj) ) * ( pt(ji+1,jj) - pt(ji,jj) ) )
            END DO
         END DO
         !
      CASE( 2 )                                                   !==  2nd order central TIM  ==! (Eq. 23)
         !
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (                                   pt(ji+1,jj) + pt(ji,jj)   &
                  &                                               -              zcu   * ( pt(ji+1,jj) - pt(ji,jj) ) ) 
            END DO
         END DO
         !  
      CASE( 3 )                                                   !==  3rd order central TIM  ==! (Eq. 24)
         !
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
               zdx2 = e1u(ji,jj) * e1u(ji,jj)
!!rachid       zdx2 = e1u(ji,jj) * e1t(ji,jj)
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (         (                         pt  (ji+1,jj) + pt  (ji,jj)        &
                  &                                               -              zcu   * ( pt  (ji+1,jj) - pt  (ji,jj) )  )   &
                  &        + z1_6 * zdx2 * ( zcu*zcu - 1._wp ) * (                         ztu2(ji+1,jj) + ztu2(ji,jj)        &
                  &                                               - SIGN( 1._wp, zcu ) * ( ztu2(ji+1,jj) - ztu2(ji,jj) )  )   )
            END DO
         END DO
         !
      CASE( 4 )                                                   !==  4th order central TIM  ==! (Eq. 27)
         !
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
               zdx2 = e1u(ji,jj) * e1u(ji,jj)
!!rachid       zdx2 = e1u(ji,jj) * e1t(ji,jj)
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (         (                   pt  (ji+1,jj) + pt  (ji,jj)        &
                  &                                               -          zcu * ( pt  (ji+1,jj) - pt  (ji,jj) )  )   &
                  &        + z1_6 * zdx2 * ( zcu*zcu - 1._wp ) * (                   ztu2(ji+1,jj) + ztu2(ji,jj)        &
                  &                                               - 0.5_wp * zcu * ( ztu2(ji+1,jj) - ztu2(ji,jj) )  )   )
            END DO
         END DO
         !
      CASE( 5 )                                                   !==  5th order central TIM  ==! (Eq. 29)
         !
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
               zdx2 = e1u(ji,jj) * e1u(ji,jj)
!!rachid       zdx2 = e1u(ji,jj) * e1t(ji,jj)
               zdx4 = zdx2 * zdx2
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (               (                   pt  (ji+1,jj) + pt  (ji,jj)       &
                  &                                                     -          zcu * ( pt  (ji+1,jj) - pt  (ji,jj) ) )   &
                  &        + z1_6   * zdx2 * ( zcu*zcu - 1._wp ) *     (                   ztu2(ji+1,jj) + ztu2(ji,jj)       &
                  &                                                     - 0.5_wp * zcu * ( ztu2(ji+1,jj) - ztu2(ji,jj) ) )   &
                  &        + z1_120 * zdx4 * ( zcu*zcu - 1._wp ) * ( zcu*zcu - 4._wp ) * ( ztu4(ji+1,jj) + ztu4(ji,jj)       &
                  &                                               - SIGN( 1._wp, zcu ) * ( ztu4(ji+1,jj) - ztu4(ji,jj) ) ) )
            END DO
         END DO
         !
      END SELECT
      !                                                     !-- High order flux in i-direction  --!
      IF( ll_neg ) THEN
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               IF( pt_u(ji,jj) < 0._wp ) THEN
                  pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (                              pt(ji+1,jj) + pt(ji,jj)   &
                     &                                    - SIGN( 1._wp, pu(ji,jj) ) * ( pt(ji+1,jj) - pt(ji,jj) ) )
               ENDIF
            END DO
         END DO
      ENDIF

      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            IF( ll_clem ) THEN
               pfu_ho(ji,jj) = pu(ji,jj) * pt_u(ji,jj)
            ELSE
               pfu_ho(ji,jj) = puc(ji,jj) * pt_u(ji,jj)
            ENDIF
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
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pv        ! ice j-velocity component
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pvc       ! ice j-velocity*A component
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pt_v      ! tracer at v-point 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pfv_ho    ! high order flux 
      !
      INTEGER  ::   ji, jj       ! dummy loop indices
      REAL(wp) ::   zcv, zdy2, zdy4    !   -      -
      REAL(wp), DIMENSION(jpi,jpj) ::   ztv1, ztv2, ztv3, ztv4
      !!----------------------------------------------------------------------
      !
      !                                                     !--  Laplacian in j-direction  --!
      DO jj = 1, jpjm1         ! First derivative (gradient)
         DO ji = fs_2, fs_jpim1
            ztv1(ji,jj) = ( pt(ji,jj+1) - pt(ji,jj) ) * r1_e2v(ji,jj) * vmask(ji,jj,1)
         END DO
      END DO
      DO jj = 2, jpjm1         ! Second derivative (Laplacian)
         DO ji = fs_2, fs_jpim1
            ztv2(ji,jj) = ( ztv1(ji,jj) - ztv1(ji,jj-1) ) * r1_e2t(ji,jj)
         END DO
      END DO
      CALL lbc_lnk( ztv2, 'T', 1. )
      !
      !                                                     !--  BiLaplacian in j-direction  --!
      DO jj = 1, jpjm1         ! First derivative
         DO ji = fs_2, fs_jpim1
            ztv3(ji,jj) = ( ztv2(ji,jj+1) - ztv2(ji,jj) ) * r1_e2v(ji,jj) * vmask(ji,jj,1)
         END DO
      END DO
      DO jj = 2, jpjm1         ! Second derivative
         DO ji = fs_2, fs_jpim1
            ztv4(ji,jj) = ( ztv3(ji,jj) - ztv3(ji,jj-1) ) * r1_e2t(ji,jj)
         END DO
      END DO
      CALL lbc_lnk( ztv4, 'T', 1. )
      !
      !
      SELECT CASE (kn_umx )
      !
      CASE( 1 )                                                !==  1st order central TIM  ==! (Eq. 21)
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (                             ( pt(ji,jj+1) + pt(ji,jj) )  &
                  &                                     - SIGN( 1._wp, pv(ji,jj) ) * ( pt(ji,jj+1) - pt(ji,jj) ) )
            END DO
         END DO
         !
      CASE( 2 )                                                !==  2nd order central TIM  ==! (Eq. 23)
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (        ( pt(ji,jj+1) + pt(ji,jj) )  &
                  &                                     - zcv * ( pt(ji,jj+1) - pt(ji,jj) ) )
            END DO
         END DO
         CALL lbc_lnk( pt_v, 'V',  1. )
         !
      CASE( 3 )                                                !==  3rd order central TIM  ==! (Eq. 24)
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
               zdy2 = e2v(ji,jj) * e2v(ji,jj)
!!rachid       zdy2 = e2v(ji,jj) * e2t(ji,jj)
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (                                 ( pt  (ji,jj+1) + pt  (ji,jj)       &
                  &                                     -                        zcv   * ( pt  (ji,jj+1) - pt  (ji,jj) ) )   &
                  &        + z1_6 * zdy2 * ( zcv*zcv - 1._wp ) * (                         ztv2(ji,jj+1) + ztv2(ji,jj)       &
                  &                                               - SIGN( 1._wp, zcv ) * ( ztv2(ji,jj+1) - ztv2(ji,jj) ) ) )
            END DO
         END DO
         !
      CASE( 4 )                                                !==  4th order central TIM  ==! (Eq. 27)
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
               zdy2 = e2v(ji,jj) * e2v(ji,jj)
!!rachid       zdy2 = e2v(ji,jj) * e2t(ji,jj)
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (                           ( pt  (ji,jj+1) + pt  (ji,jj)       &
                  &                                               -          zcv * ( pt  (ji,jj+1) - pt  (ji,jj) ) )   &
                  &        + z1_6 * zdy2 * ( zcv*zcv - 1._wp ) * (                   ztv2(ji,jj+1) + ztv2(ji,jj)       &
                  &                                               - 0.5_wp * zcv * ( ztv2(ji,jj+1) - ztv2(ji,jj) ) ) )
            END DO
         END DO
         !
      CASE( 5 )                                                !==  5th order central TIM  ==! (Eq. 29)
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
               zdy2 = e2v(ji,jj) * e2v(ji,jj)
!!rachid       zdy2 = e2v(ji,jj) * e2t(ji,jj)
               zdy4 = zdy2 * zdy2
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (                                 ( pt  (ji,jj+1) + pt  (ji,jj)      &
                  &                                                     -          zcv * ( pt  (ji,jj+1) - pt  (ji,jj) ) )  &
                  &        + z1_6   * zdy2 * ( zcv*zcv - 1._wp ) *     (                   ztv2(ji,jj+1) + ztv2(ji,jj)      &
                  &                                                     - 0.5_wp * zcv * ( ztv2(ji,jj+1) - ztv2(ji,jj) ) )  &
                  &        + z1_120 * zdy4 * ( zcv*zcv - 1._wp ) * ( zcv*zcv - 4._wp ) * ( ztv4(ji,jj+1) + ztv4(ji,jj)      &
                  &                                               - SIGN( 1._wp, zcv ) * ( ztv4(ji,jj+1) - ztv4(ji,jj) ) ) )
            END DO
         END DO
         !
      END SELECT
      !                                                     !-- High order flux in j-direction  --!
      IF( ll_neg ) THEN
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1
               IF( pt_v(ji,jj) < 0._wp ) THEN
                  pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (                             ( pt(ji,jj+1) + pt(ji,jj) )  &
                     &                                     - SIGN( 1._wp, pv(ji,jj) ) * ( pt(ji,jj+1) - pt(ji,jj) ) )
               ENDIF
            END DO
         END DO
      ENDIF

      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            IF( ll_clem ) THEN
               pfv_ho(ji,jj) = pv(ji,jj) * pt_v(ji,jj)
            ELSE
               pfv_ho(ji,jj) = pvc(ji,jj) * pt_v(ji,jj)
            ENDIF
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
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   pu               ! ice i-velocity => u*e2
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   puc              ! ice i-velocity *A => u*e2*a
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   pv               ! ice j-velocity => v*e1
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   pvc              ! ice j-velocity *A => v*e1*a
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   ptc, pt, pt_low  ! before field & upstream guess of after field
      REAL(wp), DIMENSION (jpi,jpj), INTENT(inout) ::   pfv_low, pfu_low ! upstream flux
      REAL(wp), DIMENSION (jpi,jpj), INTENT(inout) ::   pfv_ho, pfu_ho   ! monotonic flux
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) ::   zpos, zneg, zbig, zsml, z1_dt, zpos2, zneg2   ! local scalars
      REAL(wp) ::   zau, zbu, zcu, zav, zbv, zcv, zup, zdo, zsign, zcoef        !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: zbetup, zbetdo, zbup, zbdo, zti_low, ztj_low, zzt
      !!----------------------------------------------------------------------
      zbig = 1.e+40_wp
      zsml = epsi20

      IF( ll_zeroup2 ) THEN
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               IF( amaxu(ji,jj) == 0._wp )   pfu_ho(ji,jj) = 0._wp
               IF( amaxv(ji,jj) == 0._wp )   pfv_ho(ji,jj) = 0._wp
            END DO
         END DO
      ENDIF
      
      IF( ll_zeroup4 ) THEN
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               IF( pfu_low(ji,jj) == 0._wp )   pfu_ho(ji,jj) = 0._wp
               IF( pfv_low(ji,jj) == 0._wp )   pfv_ho(ji,jj) = 0._wp
            END DO
         END DO
      ENDIF


      IF( ll_zeroup1 ) THEN
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               IF( ll_gurvan ) THEN
                  zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                     - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
               ELSE
                  zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                     - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                     + pt(ji,jj) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                     &                     + pt(ji,jj) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                     &         ) * tmask(ji,jj,1)
               ENDIF
               IF( zzt(ji,jj) < 0._wp ) THEN
                  pfu_ho(ji,jj)   = pfu_low(ji,jj)
                  pfv_ho(ji,jj)   = pfv_low(ji,jj)
                  WRITE(numout,*) '*** 1 negative high order zzt ***',ji,jj,zzt(ji,jj)
               ENDIF
!!               IF( ji==26 .AND. jj==86) THEN
!!                  WRITE(numout,*) 'zzt high order',zzt(ji,jj)
!!                  WRITE(numout,*) 'pfu_ho',(pfu_ho(ji,jj)) * r1_e1e2t(ji,jj) * pdt
!!                  WRITE(numout,*) 'pfv_ho',(pfv_ho(ji,jj)) * r1_e1e2t(ji,jj) * pdt
!!                  WRITE(numout,*) 'pfu_hom1',(pfu_ho(ji-1,jj)) * r1_e1e2t(ji,jj) * pdt
!!                  WRITE(numout,*) 'pfv_hom1',(pfv_ho(ji,jj-1)) * r1_e1e2t(ji,jj) * pdt
!!               ENDIF
               IF( ll_gurvan ) THEN
                  zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                     - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
               ELSE
                  zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                     - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                     + pt(ji,jj) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                     &                     + pt(ji,jj) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                     &         ) * tmask(ji,jj,1)
               ENDIF
               IF( zzt(ji,jj) < 0._wp ) THEN
                  pfu_ho(ji-1,jj) = pfu_low(ji-1,jj)
                  pfv_ho(ji,jj-1) = pfv_low(ji,jj-1)
                  WRITE(numout,*) '*** 2 negative high order zzt ***',ji,jj,zzt(ji,jj)
               ENDIF
               IF( ll_gurvan ) THEN
                  zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                     - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
               ELSE
                  zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                     - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj)  &
                     &                     + pt(ji,jj) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                     &                     + pt(ji,jj) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                     &         ) * tmask(ji,jj,1)
               ENDIF
               IF( zzt(ji,jj) < 0._wp ) THEN
                  WRITE(numout,*) '*** 3 negative high order zzt ***',ji,jj,zzt(ji,jj)
               ENDIF
            END DO
         END DO
         CALL lbc_lnk_multi( pfu_ho, 'U', -1., pfv_ho, 'V', -1. )
      ENDIF

      
      ! antidiffusive flux : high order minus low order
      ! --------------------------------------------------
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            pfu_ho(ji,jj) = pfu_ho(ji,jj) - pfu_low(ji,jj)
            pfv_ho(ji,jj) = pfv_ho(ji,jj) - pfv_low(ji,jj)
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
         
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1 
               zti_low(ji,jj)= pt_low(ji+1,jj  )
               ztj_low(ji,jj)= pt_low(ji  ,jj+1)
            END DO
         END DO
         CALL lbc_lnk_multi( zti_low, 'T', 1., ztj_low, 'T', 1. )

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

         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               IF ( pfu_ho(ji,jj) * ( pt_low(ji+1,jj) - pt_low(ji,jj) ) <= 0. .AND.  &
                  & pfv_ho(ji,jj) * ( pt_low(ji,jj+1) - pt_low(ji,jj) ) <= 0. ) THEN
                  !
                  IF(  pfu_ho(ji,jj) * ( zti_low(ji+1,jj) - zti_low(ji,jj) ) <= 0 .AND.  &
                     & pfv_ho(ji,jj) * ( ztj_low(ji,jj+1) - ztj_low(ji,jj) ) <= 0)  pfu_ho(ji,jj)=0. ; pfv_ho(ji,jj)=0.
                  !
                  IF(  pfu_ho(ji,jj) * ( pt_low(ji  ,jj) - pt_low(ji-1,jj) ) <= 0 .AND.  &
                     & pfv_ho(ji,jj) * ( pt_low(ji  ,jj) - pt_low(ji,jj-1) ) <= 0)  pfu_ho(ji,jj)=0. ; pfv_ho(ji,jj)=0.
                  !
               ENDIF
            END DO
         END DO
         CALL lbc_lnk_multi( pfu_ho, 'U', -1., pfv_ho, 'V', -1. )   ! lateral boundary cond.

      ELSEIF( ll_prelimiter_devore ) THEN
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1 
               zti_low(ji,jj)= pt_low(ji+1,jj  )
               ztj_low(ji,jj)= pt_low(ji  ,jj+1)
            END DO
         END DO
         CALL lbc_lnk_multi( zti_low, 'T', 1., ztj_low, 'T', 1. )

         z1_dt = 1._wp / pdt
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zsign = SIGN( 1., pt_low(ji+1,jj) - pt_low(ji,jj) )
               pfu_ho(ji,jj) =  zsign * MAX( 0. , MIN( ABS(pfu_ho(ji,jj)) , &
                  &                          zsign * ( pt_low (ji  ,jj) - pt_low (ji-1,jj) ) * e1e2t(ji  ,jj) * z1_dt , &
                  &                          zsign * ( zti_low(ji+1,jj) - zti_low(ji  ,jj) ) * e1e2t(ji+1,jj) * z1_dt ) )

               zsign = SIGN( 1., pt_low(ji,jj+1) - pt_low(ji,jj) )
               pfv_ho(ji,jj) =  zsign * MAX( 0. , MIN( ABS(pfv_ho(ji,jj)) , &
                  &                          zsign * ( pt_low (ji,jj  ) - pt_low (ji,jj-1) ) * e1e2t(ji,jj  ) * z1_dt , &
                  &                          zsign * ( ztj_low(ji,jj+1) - ztj_low(ji,jj  ) ) * e1e2t(ji,jj+1) * z1_dt ) )
            END DO
         END DO
         CALL lbc_lnk_multi( pfu_ho, 'U', -1., pfv_ho, 'V', -1. )   ! lateral boundary cond.
            
      ENDIF
         
      
      ! Search local extrema
      ! --------------------
      ! max/min of pt & pt_low with large negative/positive value (-/+zbig) outside ice cover
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF    ( pt(ji,jj) <= 0._wp .AND. pt_low(ji,jj) <= 0._wp ) THEN
               zbup(ji,jj) = -zbig
               zbdo(ji,jj) =  zbig
            ELSEIF( pt(ji,jj) <= 0._wp .AND. pt_low(ji,jj) > 0._wp ) THEN
               zbup(ji,jj) = pt_low(ji,jj)
               zbdo(ji,jj) = pt_low(ji,jj)
            ELSEIF( pt(ji,jj) > 0._wp .AND. pt_low(ji,jj) <= 0._wp ) THEN
               zbup(ji,jj) = pt(ji,jj)
               zbdo(ji,jj) = pt(ji,jj)
            ELSE
               zbup(ji,jj) = MAX( pt(ji,jj) , pt_low(ji,jj) )
               zbdo(ji,jj) = MIN( pt(ji,jj) , pt_low(ji,jj) )
            ENDIF
         END DO
      END DO

 
      z1_dt = 1._wp / pdt
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
            zpos = MAX( 0., pfu_ho(ji-1,jj) ) - MIN( 0., pfu_ho(ji  ,jj) ) &  ! positive/negative part of the flux
               & + MAX( 0., pfv_ho(ji,jj-1) ) - MIN( 0., pfv_ho(ji,jj  ) )
            zneg = MAX( 0., pfu_ho(ji  ,jj) ) - MIN( 0., pfu_ho(ji-1,jj) ) &
               & + MAX( 0., pfv_ho(ji,jj  ) ) - MIN( 0., pfv_ho(ji,jj-1) )
            !
            IF( ll_HgradU .AND. .NOT.ll_gurvan ) THEN
               zneg2 = (   pt(ji,jj) * MAX( 0., pu(ji,jj) - pu(ji-1,jj) ) + pt(ji,jj) * MAX( 0., pv(ji,jj) - pv(ji,jj-1) ) &
                  &    ) * ( 1. - pamsk )
               zpos2 = ( - pt(ji,jj) * MIN( 0., pu(ji,jj) - pu(ji-1,jj) ) - pt(ji,jj) * MIN( 0., pv(ji,jj) - pv(ji,jj-1) ) &
                  &    ) * ( 1. - pamsk )
            ELSE
               zneg2 = 0. ; zpos2 = 0.
            ENDIF
            !
            !                                  ! up & down beta terms
            IF( (zpos+zpos2) > 0. ) THEN ; zbetup(ji,jj) = MAX( 0._wp, zup - pt_low(ji,jj) ) / (zpos+zpos2) * e1e2t(ji,jj) * z1_dt
            ELSE                         ; zbetup(ji,jj) = 0. ! zbig
            ENDIF
            !
            IF( (zneg+zneg2) > 0. ) THEN ; zbetdo(ji,jj) = MAX( 0._wp, pt_low(ji,jj) - zdo ) / (zneg+zneg2) * e1e2t(ji,jj) * z1_dt
            ELSE                         ; zbetdo(ji,jj) = 0. ! zbig
            ENDIF
            !
            ! if all the points are outside ice cover
            IF( zup == -zbig )   zbetup(ji,jj) = 0. ! zbig
            IF( zdo ==  zbig )   zbetdo(ji,jj) = 0. ! zbig            
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
      CALL lbc_lnk_multi( zbetup, 'T', 1., zbetdo, 'T', 1. )   ! lateral boundary cond. (unchanged sign)

      
      ! monotonic flux in the y direction
      ! ---------------------------------
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            zau = MIN( 1._wp , zbetdo(ji,jj) , zbetup(ji+1,jj) )
            zbu = MIN( 1._wp , zbetup(ji,jj) , zbetdo(ji+1,jj) )
            zcu = 0.5  + SIGN( 0.5 , pfu_ho(ji,jj) )
            !
            zcoef = ( zcu * zau + ( 1._wp - zcu ) * zbu )
            
            pfu_ho(ji,jj) = pfu_ho(ji,jj) * zcoef + pfu_low(ji,jj)

!!            IF( ji==26 .AND. jj==86) THEN
!!               WRITE(numout,*) 'coefU',zcoef
!!               WRITE(numout,*) 'pfu_ho',(pfu_ho(ji,jj)) * r1_e1e2t(ji,jj) * pdt
!!               WRITE(numout,*) 'pfu_hom1',(pfu_ho(ji-1,jj)) * r1_e1e2t(ji,jj) * pdt
!!            ENDIF

         END DO
      END DO

      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            zav = MIN( 1._wp , zbetdo(ji,jj) , zbetup(ji,jj+1) )
            zbv = MIN( 1._wp , zbetup(ji,jj) , zbetdo(ji,jj+1) )
            zcv = 0.5  + SIGN( 0.5 , pfv_ho(ji,jj) )
            !
            zcoef = ( zcv * zav + ( 1._wp - zcv ) * zbv )
            
            pfv_ho(ji,jj) = pfv_ho(ji,jj) * zcoef + pfv_low(ji,jj)

!!            IF( ji==26 .AND. jj==86) THEN
!!               WRITE(numout,*) 'coefV',zcoef
!!               WRITE(numout,*) 'pfv_ho',(pfv_ho(ji,jj)) * r1_e1e2t(ji,jj) * pdt
!!               WRITE(numout,*) 'pfv_hom1',(pfv_ho(ji,jj-1)) * r1_e1e2t(ji,jj) * pdt
!!            ENDIF
         END DO
      END DO

      ! clem test
      DO jj = 2, jpjm1
         DO ji = 2, fs_jpim1   ! vector opt.
            IF( ll_gurvan ) THEN
               zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                  &                     - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
            ELSE
               zzt(ji,jj) = ( pt(ji,jj) - ( pfu_ho(ji,jj) - pfu_ho(ji-1,jj) ) * pdt * r1_e1e2t(ji,jj)  &
                  &                     - ( pfv_ho(ji,jj) - pfv_ho(ji,jj-1) ) * pdt * r1_e1e2t(ji,jj)  &
                  &                     + pt(ji,jj) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                  &                     + pt(ji,jj) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
                  &         ) * tmask(ji,jj,1)
            ENDIF
            IF( zzt(ji,jj) < -epsi20 ) THEN
               WRITE(numout,*) 'T<0 nonosc',zzt(ji,jj)
            ENDIF
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
      REAL(wp)                     , INTENT(in   )           ::   pdt          ! tracer time-step
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   )           ::   pu           ! ice i-velocity => u*e2
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   )           ::   puc          ! ice i-velocity *A => u*e2*a
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   )           ::   pt           ! ice tracer
      REAL(wp), DIMENSION (jpi,jpj), INTENT(inout)           ::   pfu_ho       ! high order flux
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ), OPTIONAL ::   pfu_ups      ! upstream flux
      !
      REAL(wp) ::   Cr, Rjm, Rj, Rjp, uCFL, zpsi, zh3, zlimiter, Rr
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp), DIMENSION (jpi,jpj) ::   zslpx       ! tracer slopes 
      !!----------------------------------------------------------------------
      !
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zslpx(ji,jj) = ( pt(ji+1,jj) - pt(ji,jj) ) * umask(ji,jj,1)
         END DO
      END DO
      CALL lbc_lnk( zslpx, 'U', -1.)   ! lateral boundary cond.
      
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            uCFL = pdt * ABS( pu(ji,jj) ) * r1_e1e2t(ji,jj)

            Rjm = zslpx(ji-1,jj)
            Rj  = zslpx(ji  ,jj)
            Rjp = zslpx(ji+1,jj)

            IF( PRESENT(pfu_ups) ) THEN

               IF( pu(ji,jj) > 0. ) THEN   ;   Rr = Rjm
               ELSE                        ;   Rr = Rjp
               ENDIF

               zh3 = pfu_ho(ji,jj) - pfu_ups(ji,jj)     
               IF( Rj > 0. ) THEN
                  zlimiter =  MAX( 0., MIN( zh3, MAX(-Rr * 0.5 * ABS(puc(ji,jj)),  &
                     &        MIN( 2. * Rr * 0.5 * ABS(puc(ji,jj)),  zh3,  1.5 * Rj * 0.5 * ABS(puc(ji,jj)) ) ) ) )
               ELSE
                  zlimiter = -MAX( 0., MIN(-zh3, MAX( Rr * 0.5 * ABS(puc(ji,jj)),  &
                     &        MIN(-2. * Rr * 0.5 * ABS(puc(ji,jj)), -zh3, -1.5 * Rj * 0.5 * ABS(puc(ji,jj)) ) ) ) )
               ENDIF
               pfu_ho(ji,jj) = pfu_ups(ji,jj) + zlimiter
               
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
               pfu_ho(ji,jj) = pfu_ho(ji,jj) - ABS( puc(ji,jj) ) * ( (1.-zpsi) + uCFL*zpsi ) * Rj * 0.5

            ENDIF
         END DO
      END DO
      CALL lbc_lnk( pfu_ho, 'U', -1.)   ! lateral boundary cond.
      !
   END SUBROUTINE limiter_x

   SUBROUTINE limiter_y( pdt, pv, pvc, pt, pfv_ho, pfv_ups )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE limiter_y  ***
      !!     
      !! **  Purpose :   compute flux limiter 
      !!----------------------------------------------------------------------
      REAL(wp)                     , INTENT(in   )           ::   pdt          ! tracer time-step
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   )           ::   pv           ! ice i-velocity => u*e2
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   )           ::   pvc          ! ice i-velocity *A => u*e2*a
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   )           ::   pt           ! ice tracer
      REAL(wp), DIMENSION (jpi,jpj), INTENT(inout)           ::   pfv_ho       ! high order flux
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ), OPTIONAL ::   pfv_ups      ! upstream flux
      !
      REAL(wp) ::   Cr, Rjm, Rj, Rjp, vCFL, zpsi, zh3, zlimiter, Rr
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp), DIMENSION (jpi,jpj) ::   zslpy       ! tracer slopes 
      !!----------------------------------------------------------------------
      !
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zslpy(ji,jj) = ( pt(ji,jj+1) - pt(ji,jj) ) * vmask(ji,jj,1)
         END DO
      END DO
      CALL lbc_lnk( zslpy, 'V', -1.)   ! lateral boundary cond.
      
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            vCFL = pdt * ABS( pv(ji,jj) ) * r1_e1e2t(ji,jj)

            Rjm = zslpy(ji,jj-1)
            Rj  = zslpy(ji,jj  )
            Rjp = zslpy(ji,jj+1)

            IF( PRESENT(pfv_ups) ) THEN

               IF( pv(ji,jj) > 0. ) THEN   ;   Rr = Rjm
               ELSE                        ;   Rr = Rjp
               ENDIF

               zh3 = pfv_ho(ji,jj) - pfv_ups(ji,jj)     
               IF( Rj > 0. ) THEN
                  zlimiter =  MAX( 0., MIN( zh3, MAX(-Rr * 0.5 * ABS(pvc(ji,jj)),  &
                     &        MIN( 2. * Rr * 0.5 * ABS(pvc(ji,jj)),  zh3,  1.5 * Rj * 0.5 * ABS(pvc(ji,jj)) ) ) ) )
               ELSE
                  zlimiter = -MAX( 0., MIN(-zh3, MAX( Rr * 0.5 * ABS(pvc(ji,jj)),  &
                     &        MIN(-2. * Rr * 0.5 * ABS(pvc(ji,jj)), -zh3, -1.5 * Rj * 0.5 * ABS(pvc(ji,jj)) ) ) ) )
               ENDIF
               pfv_ho(ji,jj) = pfv_ups(ji,jj) + zlimiter

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
               pfv_ho(ji,jj) = pfv_ho(ji,jj) - ABS( pvc(ji,jj) ) * ( (1.-zpsi) + vCFL*zpsi ) * Rj * 0.5

            ENDIF
         END DO
      END DO
      CALL lbc_lnk( pfv_ho, 'V', -1.)   ! lateral boundary cond.
      !
   END SUBROUTINE limiter_y

#else
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_adv_umx
