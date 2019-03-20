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
   !!   nonosc_ice        : compute monotonic tracer fluxes by a non-oscillatory algorithm 
   !!----------------------------------------------------------------------
   USE phycst         ! physical constant
   USE dom_oce        ! ocean domain
   USE sbc_oce , ONLY : nn_fsbc   ! update frequency of surface boundary condition
   USE ice            ! sea-ice variables
   USE icevar         ! sea-ice: operations
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_adv_umx   ! called by icedyn_adv.F90
      
   REAL(wp) ::   z1_6   = 1._wp /   6._wp   ! =1/6
   REAL(wp) ::   z1_120 = 1._wp / 120._wp   ! =1/120
   
   ! limiter: 1=nonosc_ice, 2=superbee, 3=h3(rachid)
   INTEGER ::   kn_limiter = 1

   ! if T interpolated at u/v points is negative, then interpolate T at u/v points using the upstream scheme
   !   clem: if set to true, the 2D test case "diagonal advection" does not work (I do not understand why)
   !         but in realistic cases, it avoids having very negative ice temperature (-50) at low ice concentration 
   LOGICAL ::   ll_neg = .TRUE.
   
   ! alternate directions for upstream
   LOGICAL ::   ll_upsxy = .TRUE.

   ! alternate directions for high order
   LOGICAL ::   ll_hoxy = .TRUE.
   
   ! prelimiter: use it to avoid overshoot in H
   !   clem: if set to true, the 2D test case "diagnoal advection" does not work (I do not understand why)
   LOGICAL ::   ll_prelimiter_zalesak = .FALSE.  ! from: Zalesak(1979) eq. 14 => better for 1D. Not well defined in 2D


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
      REAL(wp), DIMENSION(1)           ::   zcflprv, zcflnow   ! send zcflnow and receive zcflprv
      REAL(wp), DIMENSION(jpi,jpj)     ::   zudy, zvdx, zcu_box, zcv_box 
      REAL(wp), DIMENSION(jpi,jpj)     ::   zati1, zati2
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zua_ho, zva_ho
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z1_ai, z1_aip
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zhvar
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_dyn_adv_umx: Ultimate-Macho advection scheme'
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

      !---------------!
      !== advection ==!
      !---------------!
      DO jt = 1, icycle

         ! record at_i before advection (for open water)
         zati1(:,:) = SUM( pa_i(:,:,:), dim=3 )
         
         ! inverse of A and Ap
         WHERE( pa_i(:,:,:) >= epsi20 )   ;   z1_ai(:,:,:) = 1._wp / pa_i(:,:,:)
         ELSEWHERE                        ;   z1_ai(:,:,:) = 0.
         END WHERE
         WHERE( pa_ip(:,:,:) >= epsi20 )  ;   z1_aip(:,:,:) = 1._wp / pa_ip(:,:,:)
         ELSEWHERE                        ;   z1_aip(:,:,:) = 0.
         END WHERE
         !
         ! set u*a=u for advection of A only 
         DO jl = 1, jpl
            zua_ho(:,:,jl) = zudy(:,:)
            zva_ho(:,:,jl) = zvdx(:,:)
         END DO
         
         zamsk = 1._wp
         CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, pa_i, pa_i, zua_ho, zva_ho ) !-- Ice area
         zamsk = 0._wp
         !
         zhvar(:,:,:) = pv_i(:,:,:) * z1_ai(:,:,:)
         CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pv_i )                !-- Ice volume
         !
         zhvar(:,:,:) = pv_s(:,:,:) * z1_ai(:,:,:)
         CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pv_s )                !-- Snw volume
         !
         zhvar(:,:,:) = psv_i(:,:,:) * z1_ai(:,:,:)
         CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, psv_i )               !-- Salt content
         !
         DO jk = 1, nlay_i
            zhvar(:,:,:) = pe_i(:,:,jk,:) * z1_ai(:,:,:)
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pe_i(:,:,jk,:) )   !-- Ice heat content
         END DO
         !
         DO jk = 1, nlay_s
            zhvar(:,:,:) = pe_s(:,:,jk,:) * z1_ai(:,:,:)
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pe_s(:,:,jk,:) )   !-- Snw heat content
         END DO
         !
         IF( iom_use('iceage') .OR. iom_use('iceage_cat') ) THEN                                                              !-- Ice Age
            ! clem: in theory we should use the formulation below to advect the ice age, but the code is unable to deal with
            !       fields that do not depend on volume (here oa_i depends on concentration). It creates abnormal ages that
            !       spread into the domain. Therefore we cheat and consider that ice age should be advected as ice concentration
            !!zhvar(:,:,:) = poa_i(:,:,:) * z1_ai(:,:,:)
            !!CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, poa_i )
            ! set u*a=u for advection of ice age
            DO jl = 1, jpl
               zua_ho(:,:,jl) = zudy(:,:)
               zva_ho(:,:,jl) = zvdx(:,:)
            END DO
            zamsk = 1._wp
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, poa_i, poa_i )
            zamsk = 0._wp
         ENDIF
         !
         IF ( ln_pnd_H12 ) THEN                                                                                               !-- melt ponds
            ! set u*a=u for advection of Ap only 
            DO jl = 1, jpl
               zua_ho(:,:,jl) = zudy(:,:)
               zva_ho(:,:,jl) = zvdx(:,:)
            END DO
            !
            zamsk = 1._wp
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, pa_ip, pa_ip, zua_ho, zva_ho ) ! fraction
            zamsk = 0._wp
            !
            zhvar(:,:,:) = pv_ip(:,:,:) * z1_aip(:,:,:)
            CALL adv_umx( zamsk, kn_umx, jt, kt, zdt, zudy, zvdx, zua_ho, zva_ho, zcu_box, zcv_box, zhvar, pv_ip )                 ! volume
         ENDIF
         !
         zati2(:,:) = SUM( pa_i(:,:,:), dim=3 )
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               pato_i(ji,jj) = pato_i(ji,jj) - ( zati2(ji,jj) - zati1(ji,jj) ) &                                              !-- Open water area
                  &                          - ( zudy(ji,jj) - zudy(ji-1,jj) + zvdx(ji,jj) - zvdx(ji,jj-1) ) * r1_e1e2t(ji,jj) * zdt
            END DO
         END DO
         CALL lbc_lnk( 'icedyn_adv_umx', pato_i(:,:), 'T',  1. )
         !
      END DO
      !
   END SUBROUTINE ice_dyn_adv_umx

   
   SUBROUTINE adv_umx( pamsk, kn_umx, jt, kt, pdt, pu, pv, puc, pvc, pubox, pvbox, pt, ptc, pua_ho, pva_ho )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE adv_umx  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of 
      !!                 tracers and add it to the general trend of tracer equations
      !!
      !! **  Method  :   - calculate upstream fluxes and upstream solution for tracer H
      !!                 - calculate tracer H at u and v points (Ultimate)
      !!                 - calculate the high order fluxes using alterning directions (Macho?)
      !!                 - apply a limiter on the fluxes (nonosc_ice)
      !!                 - convert this tracer flux to a tracer content flux (uH -> uV)
      !!                 - calculate the high order solution for tracer content V
      !!
      !! ** Action : solve 2 equations => a) da/dt = -div(ua)
      !!                                  b) dV/dt = -div(uV) using dH/dt = -u.grad(H)
      !!             in eq. b), - fluxes uH are evaluated (with UMx) and limited (with nonosc_ice). This step is necessary to get a good H.
      !!                        - then we convert this flux to a "volume" flux this way => uH*ua/u
      !!                             where ua is the flux from eq. a)
      !!                        - at last we estimate dV/dt = -div(uH*ua/u)
      !!
      !! ** Note : - this method can lead to small negative V (since we only limit H) => corrected in icedyn_adv.F90 conserving mass etc.
      !!           - negative tracers at u-v points can also occur from the Ultimate scheme (usually at the ice edge) and the solution for now
      !!             is to apply an upstream scheme when it occurs. A better solution would be to degrade the order of
      !!             the scheme automatically by applying a mask of the ice cover inside Ultimate (not done).
      !!           - Eventhough 1D tests give very good results (typically the one from Schar & Smolarkiewiecz), the 2D is less good.
      !!             Large values of H can appear for very small ice concentration, and when it does it messes the things up since we
      !!             work on H (and not V). It probably comes from the prelimiter of zalesak which is coded for 1D and not 2D.
      !!             Therefore, after advection we limit the thickness to the largest value of the 9-points around (only if ice
      !!             concentration is small).
      !! To-do: expand the prelimiter from zalesak to make it work in 2D
      !!----------------------------------------------------------------------
      REAL(wp)                        , INTENT(in   )           ::   pamsk          ! advection of concentration (1) or other tracers (0)
      INTEGER                         , INTENT(in   )           ::   kn_umx         ! order of the scheme (1-5=UM or 20=CEN2)
      INTEGER                         , INTENT(in   )           ::   jt             ! number of sub-iteration
      INTEGER                         , INTENT(in   )           ::   kt             ! number of iteration
      REAL(wp)                        , INTENT(in   )           ::   pdt            ! tracer time-step
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   )           ::   pu   , pv      ! 2 ice velocity components => u*e2
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   )           ::   puc  , pvc     ! 2 ice velocity components => u*e2 or u*a*e2u
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   )           ::   pubox, pvbox   ! upstream velocity
      REAL(wp), DIMENSION(:,:,:)      , INTENT(inout)           ::   pt             ! tracer field
      REAL(wp), DIMENSION(:,:,:)      , INTENT(inout)           ::   ptc            ! tracer content field
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out), OPTIONAL ::   pua_ho, pva_ho ! high order u*a fluxes
      !
      INTEGER  ::   ji, jj, jl       ! dummy loop indices  
      REAL(wp) ::   ztra             ! local scalar
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zfu_ho , zfv_ho , zpt
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zfu_ups, zfv_ups, zt_ups
      !!----------------------------------------------------------------------
      !
      ! Upstream (_ups) fluxes 
      ! -----------------------
      CALL upstream( pamsk, jt, kt, pdt, pt, pu, pv, zt_ups, zfu_ups, zfv_ups )
      
      ! High order (_ho) fluxes 
      ! -----------------------
      SELECT CASE( kn_umx )
         !
      CASE ( 20 )                          !== centered second order ==!
         !
         CALL cen2( pamsk, jt, kt, pdt, pt, pu, pv, zt_ups, zfu_ups, zfv_ups, zfu_ho, zfv_ho )
         !
      CASE ( 1:5 )                         !== 1st to 5th order ULTIMATE-MACHO scheme ==!
         !
         CALL macho( pamsk, kn_umx, jt, kt, pdt, pt, pu, pv, pubox, pvbox, zt_ups, zfu_ups, zfv_ups, zfu_ho, zfv_ho )
         !
      END SELECT
      !
      !              --ho    --ho
      ! new fluxes = u*H  *  u*a / u
      ! ----------------------------
      IF( pamsk == 0._wp ) THEN
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
      !                                   --ho
      ! in case of advection of A: output u*a
      ! -------------------------------------
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
      !
      ! final trend with corrected fluxes
      ! ---------------------------------
      DO jl = 1, jpl
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1 
               ztra = - ( zfu_ho(ji,jj,jl) - zfu_ho(ji-1,jj,jl) + zfv_ho(ji,jj,jl) - zfv_ho(ji,jj-1,jl) )  
               !
               ptc(ji,jj,jl) = ( ptc(ji,jj,jl) + ztra * r1_e1e2t(ji,jj) * pdt ) * tmask(ji,jj,1)               
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', ptc, 'T',  1. )
      !
   END SUBROUTINE adv_umx


   SUBROUTINE upstream( pamsk, jt, kt, pdt, pt, pu, pv, pt_ups, pfu_ups, pfv_ups )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE upstream  ***
      !!     
      !! **  Purpose :   compute the upstream fluxes and upstream guess of tracer
      !!----------------------------------------------------------------------
      REAL(wp)                        , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      INTEGER                         , INTENT(in   ) ::   jt               ! number of sub-iteration
      INTEGER                         , INTENT(in   ) ::   kt               ! number of iteration
      REAL(wp)                        , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pt               ! tracer fields
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   ) ::   pu, pv           ! 2 ice velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pt_ups           ! upstream guess of tracer 
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfu_ups, pfv_ups ! upstream fluxes 
      !
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp) ::   ztra          ! local scalar
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zpt
      !!----------------------------------------------------------------------

      IF( .NOT. ll_upsxy ) THEN         !** no alternate directions **!
         !
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  pfu_ups(ji,jj,jl) = MAX( pu(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pu(ji,jj), 0._wp ) * pt(ji+1,jj,jl)
                  pfv_ups(ji,jj,jl) = MAX( pv(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pv(ji,jj), 0._wp ) * pt(ji,jj+1,jl)
               END DO
            END DO
         END DO
         !
      ELSE                              !** alternate directions **!
         !
         IF( MOD( (kt - 1) / nn_fsbc , 2 ) ==  MOD( (jt - 1) , 2 ) ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
            !
            DO jl = 1, jpl              !-- flux in x-direction
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     pfu_ups(ji,jj,jl) = MAX( pu(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pu(ji,jj), 0._wp ) * pt(ji+1,jj,jl)
                  END DO
               END DO
            END DO
            !
            DO jl = 1, jpl              !-- first guess of tracer from u-flux
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ztra = - ( pfu_ups(ji,jj,jl) - pfu_ups(ji-1,jj,jl) )              &
                        &   + ( pu     (ji,jj   ) - pu     (ji-1,jj   ) ) * pt(ji,jj,jl) * (1.-pamsk)
                     !
                     zpt(ji,jj,jl) = ( pt(ji,jj,jl) + ztra * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zpt, 'T', 1. )
            !
            DO jl = 1, jpl              !-- flux in y-direction
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     pfv_ups(ji,jj,jl) = MAX( pv(ji,jj), 0._wp ) * zpt(ji,jj,jl) + MIN( pv(ji,jj), 0._wp ) * zpt(ji,jj+1,jl)
                  END DO
               END DO
            END DO
            !
         ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
            !
            DO jl = 1, jpl              !-- flux in y-direction
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     pfv_ups(ji,jj,jl) = MAX( pv(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pv(ji,jj), 0._wp ) * pt(ji,jj+1,jl)
                  END DO
               END DO
            END DO
            !
            DO jl = 1, jpl              !-- first guess of tracer from v-flux
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ztra = - ( pfv_ups(ji,jj,jl) - pfv_ups(ji,jj-1,jl) )  &
                        &   + ( pv     (ji,jj   ) - pv     (ji,jj-1   ) ) * pt(ji,jj,jl) * (1.-pamsk)
                     !
                     zpt(ji,jj,jl) = ( pt(ji,jj,jl) + ztra * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zpt, 'T', 1. )
            !
            DO jl = 1, jpl              !-- flux in x-direction
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     pfu_ups(ji,jj,jl) = MAX( pu(ji,jj), 0._wp ) * zpt(ji,jj,jl) + MIN( pu(ji,jj), 0._wp ) * zpt(ji+1,jj,jl)
                  END DO
               END DO
            END DO
            !
         ENDIF
         
      ENDIF
      !
      DO jl = 1, jpl                    !-- after tracer with upstream scheme
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               ztra = - (   pfu_ups(ji,jj,jl) - pfu_ups(ji-1,jj  ,jl)   &
                  &       + pfv_ups(ji,jj,jl) - pfv_ups(ji  ,jj-1,jl) ) &
                  &   + (   pu     (ji,jj   ) - pu     (ji-1,jj     )   &
                  &       + pv     (ji,jj   ) - pv     (ji  ,jj-1   ) ) * pt(ji,jj,jl) * (1.-pamsk)
               !
               pt_ups(ji,jj,jl) = ( pt(ji,jj,jl) + ztra * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'icedyn_adv_umx', pt_ups, 'T', 1. )

   END SUBROUTINE upstream

   
   SUBROUTINE cen2( pamsk, jt, kt, pdt, pt, pu, pv, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE cen2  ***
      !!     
      !! **  Purpose :   compute the high order fluxes using a centered
      !!                 second order scheme 
      !!----------------------------------------------------------------------
      REAL(wp)                        , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      INTEGER                         , INTENT(in   ) ::   jt               ! number of sub-iteration
      INTEGER                         , INTENT(in   ) ::   kt               ! number of iteration
      REAL(wp)                        , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pt               ! tracer fields
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   ) ::   pu, pv           ! 2 ice velocity components
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pt_ups           ! upstream guess of tracer 
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pfu_ups, pfv_ups ! upstream fluxes 
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfu_ho, pfv_ho   ! high order fluxes 
      !
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp) ::   ztra          ! local scalar
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zpt
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ll_hoxy ) THEN           !** no alternate directions **!
         !
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  pfu_ho(ji,jj,jl) = 0.5_wp * pu(ji,jj) * ( pt(ji,jj,jl) + pt(ji+1,jj  ,jl) )
                  pfv_ho(ji,jj,jl) = 0.5_wp * pv(ji,jj) * ( pt(ji,jj,jl) + pt(ji  ,jj+1,jl) )
               END DO
            END DO
         END DO
         !
         IF    ( kn_limiter == 1 ) THEN
            CALL nonosc_ice( pamsk, pdt, pu, pv, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
         ELSEIF( kn_limiter == 2 .OR. kn_limiter == 3 ) THEN
            CALL limiter_x( pdt, pu, pt, pfu_ups, pfu_ho )
            CALL limiter_y( pdt, pv, pt, pfv_ups, pfv_ho )
         ENDIF
         !
      ELSE                              !** alternate directions **!
         !
         IF( MOD( (kt - 1) / nn_fsbc , 2 ) ==  MOD( (jt - 1) , 2 ) ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
            !
            DO jl = 1, jpl              !-- flux in x-direction
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     pfu_ho(ji,jj,jl) = 0.5_wp * pu(ji,jj) * ( pt(ji,jj,jl) + pt(ji+1,jj,jl) )
                  END DO
               END DO
            END DO
            IF( kn_limiter == 2 .OR. kn_limiter == 3 )   CALL limiter_x( pdt, pu, pt, pfu_ups, pfu_ho )

            DO jl = 1, jpl              !-- first guess of tracer from u-flux
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ztra = - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) )              &
                        &   + ( pu    (ji,jj   ) - pu    (ji-1,jj   ) ) * pt(ji,jj,jl) * (1.-pamsk)
                     !
                     zpt(ji,jj,jl) = ( pt(ji,jj,jl) + ztra * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zpt, 'T', 1. )

            DO jl = 1, jpl              !-- flux in y-direction
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     pfv_ho(ji,jj,jl) = 0.5_wp * pv(ji,jj) * ( zpt(ji,jj,jl) + zpt(ji,jj+1,jl) )
                  END DO
               END DO
            END DO
            IF( kn_limiter == 2 .OR. kn_limiter == 3 )   CALL limiter_y( pdt, pv, pt, pfv_ups, pfv_ho )

         ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
            !
            DO jl = 1, jpl              !-- flux in y-direction
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     pfv_ho(ji,jj,jl) = 0.5_wp * pv(ji,jj) * ( pt(ji,jj,jl) + pt(ji,jj+1,jl) )
                  END DO
               END DO
            END DO
            IF( kn_limiter == 2 .OR. kn_limiter == 3 )   CALL limiter_y( pdt, pv, pt, pfv_ups, pfv_ho )
            !
            DO jl = 1, jpl              !-- first guess of tracer from v-flux
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ztra = - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) )  &
                        &   + ( pv    (ji,jj   ) - pv    (ji,jj-1   ) ) * pt(ji,jj,jl) * (1.-pamsk)
                     !
                     zpt(ji,jj,jl) = ( pt(ji,jj,jl) + ztra * pdt * r1_e1e2t(ji,jj) ) * tmask(ji,jj,1)
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'icedyn_adv_umx', zpt, 'T', 1. )
            !
            DO jl = 1, jpl              !-- flux in x-direction
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     pfu_ho(ji,jj,jl) = 0.5_wp * pu(ji,jj) * ( zpt(ji,jj,jl) + zpt(ji+1,jj,jl) )
                  END DO
               END DO
            END DO
            IF( kn_limiter == 2 .OR. kn_limiter == 3 )   CALL limiter_x( pdt, pu, pt, pfu_ups, pfu_ho )

         ENDIF
         IF( kn_limiter == 1 )   CALL nonosc_ice( pamsk, pdt, pu, pv, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
         
      ENDIF
   
   END SUBROUTINE cen2

   
   SUBROUTINE macho( pamsk, kn_umx, jt, kt, pdt, pt, pu, pv, pubox, pvbox, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE macho  ***
      !!     
      !! **  Purpose :   compute the high order fluxes using Ultimate-Macho scheme  
      !!
      !! **  Method  :   ...
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      REAL(wp)                        , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      INTEGER                         , INTENT(in   ) ::   kn_umx           ! order of the scheme (1-5=UM or 20=CEN2)
      INTEGER                         , INTENT(in   ) ::   jt               ! number of sub-iteration
      INTEGER                         , INTENT(in   ) ::   kt               ! number of iteration
      REAL(wp)                        , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pt               ! tracer fields
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   ) ::   pu, pv           ! 2 ice velocity components
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   ) ::   pubox, pvbox     ! upstream velocity
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pt_ups           ! upstream guess of tracer 
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pfu_ups, pfv_ups ! upstream fluxes 
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfu_ho, pfv_ho   ! high order fluxes 
      !
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zt_u, zt_v, zpt
      !!----------------------------------------------------------------------
      !
      IF( MOD( (kt - 1) / nn_fsbc , 2 ) ==  MOD( (jt - 1) , 2 ) ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
         !
         !                                                        !--  ultimate interpolation of pt at u-point  --!
         CALL ultimate_x( kn_umx, pdt, pt, pu, zt_u, pfu_ho )
         !                                                        !--  limiter in x --!
         IF( kn_limiter == 2 .OR. kn_limiter == 3 )   CALL limiter_x( pdt, pu, pt, pfu_ups, pfu_ho )
         !                                                        !--  advective form update in zpt  --!
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  zpt(ji,jj,jl) = ( pt(ji,jj,jl) - (  pubox(ji,jj   ) * ( zt_u(ji,jj,jl) - zt_u(ji-1,jj,jl) ) * r1_e1t  (ji,jj) &
                     &                              + pt   (ji,jj,jl) * ( pu  (ji,jj   ) - pu  (ji-1,jj   ) ) * r1_e1e2t(ji,jj) &
                     &                                                                                        * pamsk           &
                     &                             ) * pdt ) * tmask(ji,jj,1) 
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'icedyn_adv_umx', zpt, 'T', 1. )
         !
         !                                                        !--  ultimate interpolation of pt at v-point  --!
         IF( ll_hoxy ) THEN
            CALL ultimate_y( kn_umx, pdt, zpt, pv, zt_v, pfv_ho )
         ELSE
            CALL ultimate_y( kn_umx, pdt, pt , pv, zt_v, pfv_ho )
         ENDIF
         !                                                        !--  limiter in y --!
         IF( kn_limiter == 2 .OR. kn_limiter == 3 )   CALL limiter_y( pdt, pv, pt, pfv_ups, pfv_ho )
         !         
         !
      ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
         !
         !                                                        !--  ultimate interpolation of pt at v-point  --!
         CALL ultimate_y( kn_umx, pdt, pt, pv, zt_v, pfv_ho )
         !                                                        !--  limiter in y --!
         IF( kn_limiter == 2 .OR. kn_limiter == 3 )   CALL limiter_y( pdt, pv, pt, pfv_ups, pfv_ho )
         !                                                        !--  advective form update in zpt  --!
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  zpt(ji,jj,jl) = ( pt(ji,jj,jl) - (  pvbox(ji,jj   ) * ( zt_v(ji,jj,jl) - zt_v(ji,jj-1,jl) ) * r1_e2t  (ji,jj) &
                     &                              + pt   (ji,jj,jl) * ( pv  (ji,jj   ) - pv  (ji,jj-1   ) ) * r1_e1e2t(ji,jj) &
                     &                                                                                        * pamsk           &
                     &                             ) * pdt ) * tmask(ji,jj,1) 
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'icedyn_adv_umx', zpt, 'T', 1. )
         !
         !                                                        !--  ultimate interpolation of pt at u-point  --!
         IF( ll_hoxy ) THEN
            CALL ultimate_x( kn_umx, pdt, zpt, pu, zt_u, pfu_ho )
         ELSE
            CALL ultimate_x( kn_umx, pdt, pt , pu, zt_u, pfu_ho )
         ENDIF
         !                                                        !--  limiter in x --!
         IF( kn_limiter == 2 .OR. kn_limiter == 3 )   CALL limiter_x( pdt, pu, pt, pfu_ups, pfu_ho )
         !
      ENDIF

      IF( kn_limiter == 1 )   CALL nonosc_ice( pamsk, pdt, pu, pv, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
      !
   END SUBROUTINE macho


   SUBROUTINE ultimate_x( kn_umx, pdt, pt, pu, pt_u, pfu_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ultimate_x  ***
      !!     
      !! **  Purpose :   compute tracer at u-points 
      !!
      !! **  Method  :   ...
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kn_umx    ! order of the scheme (1-5=UM or 20=CEN2)
      REAL(wp)                        , INTENT(in   ) ::   pdt       ! tracer time-step
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   ) ::   pu        ! ice i-velocity component
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pt_u      ! tracer at u-point 
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfu_ho    ! high order flux 
      !
      INTEGER  ::   ji, jj, jl             ! dummy loop indices
      REAL(wp) ::   zcu, zdx2, zdx4        !   -      -
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
                  pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (                                pt(ji+1,jj,jl) + pt(ji,jj,jl)   &
                     &                                         - SIGN( 1._wp, pu(ji,jj) ) * ( pt(ji+1,jj,jl) - pt(ji,jj,jl) ) )
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
                     &                                                            - zcu   * ( pt(ji+1,jj,jl) - pt(ji,jj,jl) ) ) 
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
!!rachid          zdx2 = e1u(ji,jj) * e1t(ji,jj)
                  pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (         (                      pt  (ji+1,jj,jl) + pt  (ji,jj,jl)     &
                     &                                                            - zcu   * ( pt  (ji+1,jj,jl) - pt  (ji,jj,jl) ) ) &
                     &        + z1_6 * zdx2 * ( zcu*zcu - 1._wp ) *    (                      ztu2(ji+1,jj,jl) + ztu2(ji,jj,jl)     &
                     &                                               - SIGN( 1._wp, zcu ) * ( ztu2(ji+1,jj,jl) - ztu2(ji,jj,jl) ) ) )
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
!!rachid          zdx2 = e1u(ji,jj) * e1t(ji,jj)
                  pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (         (                      pt  (ji+1,jj,jl) + pt  (ji,jj,jl)     &
                     &                                                            - zcu   * ( pt  (ji+1,jj,jl) - pt  (ji,jj,jl) ) ) &
                     &        + z1_6 * zdx2 * ( zcu*zcu - 1._wp ) *    (                      ztu2(ji+1,jj,jl) + ztu2(ji,jj,jl)     &
                     &                                                   - 0.5_wp * zcu   * ( ztu2(ji+1,jj,jl) - ztu2(ji,jj,jl) ) ) )
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
!!rachid          zdx2 = e1u(ji,jj) * e1t(ji,jj)
                  zdx4 = zdx2 * zdx2
                  pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (        (                       pt  (ji+1,jj,jl) + pt  (ji,jj,jl)     &
                     &                                                            - zcu   * ( pt  (ji+1,jj,jl) - pt  (ji,jj,jl) ) ) &
                     &        + z1_6   * zdx2 * ( zcu*zcu - 1._wp ) * (                       ztu2(ji+1,jj,jl) + ztu2(ji,jj,jl)     &
                     &                                                   - 0.5_wp * zcu   * ( ztu2(ji+1,jj,jl) - ztu2(ji,jj,jl) ) ) &
                     &        + z1_120 * zdx4 * ( zcu*zcu - 1._wp ) * ( zcu*zcu - 4._wp ) * ( ztu4(ji+1,jj,jl) + ztu4(ji,jj,jl)     &
                     &                                               - SIGN( 1._wp, zcu ) * ( ztu4(ji+1,jj,jl) - ztu4(ji,jj,jl) ) ) )
               END DO
            END DO
         END DO
         !
      END SELECT
      !
      ! if pt at u-point is negative then use the upstream value
      !    this should not be necessary if a proper sea-ice mask is set in Ultimate
      !    to degrade the order of the scheme when necessary (for ex. at the ice edge)
      IF( ll_neg ) THEN
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( pt_u(ji,jj,jl) < 0._wp ) THEN
                     pt_u(ji,jj,jl) = 0.5_wp * umask(ji,jj,1) * (                                pt(ji+1,jj,jl) + pt(ji,jj,jl)   &
                        &                                         - SIGN( 1._wp, pu(ji,jj) ) * ( pt(ji+1,jj,jl) - pt(ji,jj,jl) ) )
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
      !                                                     !-- High order flux in i-direction  --!
      DO jl = 1, jpl
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               pfu_ho(ji,jj,jl) = pu(ji,jj) * pt_u(ji,jj,jl)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE ultimate_x
   
 
   SUBROUTINE ultimate_y( kn_umx, pdt, pt, pv, pt_v, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ultimate_y  ***
      !!     
      !! **  Purpose :   compute tracer at v-points 
      !!
      !! **  Method  :   ...
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kn_umx    ! order of the scheme (1-5=UM or 20=CEN2)
      REAL(wp)                        , INTENT(in   ) ::   pdt       ! tracer time-step
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   ) ::   pv        ! ice j-velocity component
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pt_v      ! tracer at v-point 
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfv_ho    ! high order flux 
      !
      INTEGER  ::   ji, jj, jl         ! dummy loop indices
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
                  pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (                                pt(ji,jj+1,jl) + pt(ji,jj,jl)   &
                     &                                         - SIGN( 1._wp, pv(ji,jj) ) * ( pt(ji,jj+1,jl) - pt(ji,jj,jl) ) )
               END DO
            END DO
         END DO
         !
      CASE( 2 )                                                !==  2nd order central TIM  ==! (Eq. 23)
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
                  pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (                                pt(ji,jj+1,jl) + pt(ji,jj,jl)   &
                     &                                                            - zcv *   ( pt(ji,jj+1,jl) - pt(ji,jj,jl) ) )
               END DO
            END DO
         END DO
         !
      CASE( 3 )                                                !==  3rd order central TIM  ==! (Eq. 24)
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
                  zdy2 = e2v(ji,jj) * e2v(ji,jj)
!!rachid          zdy2 = e2v(ji,jj) * e2t(ji,jj)
                  pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (      (                         pt  (ji,jj+1,jl) + pt  (ji,jj,jl)     &
                     &                                                            - zcv   * ( pt  (ji,jj+1,jl) - pt  (ji,jj,jl) ) ) &
                     &        + z1_6 * zdy2 * ( zcv*zcv - 1._wp ) * (                         ztv2(ji,jj+1,jl) + ztv2(ji,jj,jl)     &
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
!!rachid          zdy2 = e2v(ji,jj) * e2t(ji,jj)
                  pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (      (                         pt  (ji,jj+1,jl) + pt  (ji,jj,jl)     &
                     &                                                            - zcv   * ( pt  (ji,jj+1,jl) - pt  (ji,jj,jl) ) ) &
                     &        + z1_6 * zdy2 * ( zcv*zcv - 1._wp ) * (                         ztv2(ji,jj+1,jl) + ztv2(ji,jj,jl)     &
                     &                                                   - 0.5_wp * zcv   * ( ztv2(ji,jj+1,jl) - ztv2(ji,jj,jl) ) ) )
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
!!rachid          zdy2 = e2v(ji,jj) * e2t(ji,jj)
                  zdy4 = zdy2 * zdy2
                  pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (                              ( pt  (ji,jj+1,jl) + pt  (ji,jj,jl)     &
                     &                                                            - zcv   * ( pt  (ji,jj+1,jl) - pt  (ji,jj,jl) ) ) &
                     &        + z1_6   * zdy2 * ( zcv*zcv - 1._wp ) * (                       ztv2(ji,jj+1,jl) + ztv2(ji,jj,jl)     &
                     &                                                   - 0.5_wp * zcv   * ( ztv2(ji,jj+1,jl) - ztv2(ji,jj,jl) ) ) &
                     &        + z1_120 * zdy4 * ( zcv*zcv - 1._wp ) * ( zcv*zcv - 4._wp ) * ( ztv4(ji,jj+1,jl) + ztv4(ji,jj,jl)     &
                     &                                               - SIGN( 1._wp, zcv ) * ( ztv4(ji,jj+1,jl) - ztv4(ji,jj,jl) ) ) )
               END DO
            END DO
         END DO
         !
      END SELECT
      !
      ! if pt at v-point is negative then use the upstream value
      !    this should not be necessary if a proper sea-ice mask is set in Ultimate
      !    to degrade the order of the scheme when necessary (for ex. at the ice edge)
      IF( ll_neg ) THEN
         DO jl = 1, jpl
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  IF( pt_v(ji,jj,jl) < 0._wp ) THEN
                     pt_v(ji,jj,jl) = 0.5_wp * vmask(ji,jj,1) * (                              ( pt(ji,jj+1,jl) + pt(ji,jj,jl) )  &
                        &                                         - SIGN( 1._wp, pv(ji,jj) ) * ( pt(ji,jj+1,jl) - pt(ji,jj,jl) ) )
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
      !                                                     !-- High order flux in j-direction  --!
      DO jl = 1, jpl
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               pfv_ho(ji,jj,jl) = pv(ji,jj) * pt_v(ji,jj,jl)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE ultimate_y
     

   SUBROUTINE nonosc_ice( pamsk, pdt, pu, pv, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc_ice  ***
      !!     
      !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
      !!       scheme and the before field by a non-oscillatory algorithm 
      !!
      !! **  Method  :   ...
      !!----------------------------------------------------------------------
      REAL(wp)                   , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      REAL(wp)                   , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION (:,:  ), INTENT(in   ) ::   pu               ! ice i-velocity => u*e2
      REAL(wp), DIMENSION (:,:  ), INTENT(in   ) ::   pv               ! ice j-velocity => v*e1
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   pt, pt_ups       ! before field & upstream guess of after field
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   pfv_ups, pfu_ups ! upstream flux
      REAL(wp), DIMENSION (:,:,:), INTENT(inout) ::   pfv_ho, pfu_ho   ! monotonic flux
      !
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp) ::   zpos, zneg, zbig, zup, zdo, z1_dt              ! local scalars
      REAL(wp) ::   zau, zbu, zcu, zav, zbv, zcv, zcoef, zzt       !   -      -
      REAL(wp), DIMENSION(jpi,jpj    ) :: zbup, zbdo
      REAL(wp), DIMENSION(jpi,jpj,jpl) :: zbetup, zbetdo, zti_ups, ztj_ups
      !!----------------------------------------------------------------------
      zbig = 1.e+40_wp
      
      ! antidiffusive flux : high order minus low order
      ! --------------------------------------------------
      DO jl = 1, jpl
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               pfu_ho(ji,jj,jl) = pfu_ho(ji,jj,jl) - pfu_ups(ji,jj,jl)
               pfv_ho(ji,jj,jl) = pfv_ho(ji,jj,jl) - pfv_ups(ji,jj,jl)
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
      !            t_ups :       i-1     i       i+1       i+2   
      IF( ll_prelimiter_zalesak ) THEN
         
         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1 
                  zti_ups(ji,jj,jl)= pt_ups(ji+1,jj  ,jl)
                  ztj_ups(ji,jj,jl)= pt_ups(ji  ,jj+1,jl)
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'icedyn_adv_umx', zti_ups, 'T', 1., ztj_ups, 'T', 1. )

         DO jl = 1, jpl
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  IF ( pfu_ho(ji,jj,jl) * ( pt_ups(ji+1,jj  ,jl) - pt_ups(ji,jj,jl) ) <= 0._wp .AND.  &
                     & pfv_ho(ji,jj,jl) * ( pt_ups(ji  ,jj+1,jl) - pt_ups(ji,jj,jl) ) <= 0._wp ) THEN
                     !
                     IF(  pfu_ho(ji,jj,jl) * ( zti_ups(ji+1,jj  ,jl) - zti_ups(ji,jj,jl) ) <= 0._wp .AND.  &
                        & pfv_ho(ji,jj,jl) * ( ztj_ups(ji  ,jj+1,jl) - ztj_ups(ji,jj,jl) ) <= 0._wp ) THEN
                        pfu_ho(ji,jj,jl)=0._wp
                        pfv_ho(ji,jj,jl)=0._wp
                     ENDIF
                     !
                     IF(  pfu_ho(ji,jj,jl) * ( pt_ups(ji,jj,jl) - pt_ups(ji-1,jj  ,jl) ) <= 0._wp .AND.  &
                        & pfv_ho(ji,jj,jl) * ( pt_ups(ji,jj,jl) - pt_ups(ji  ,jj-1,jl) ) <= 0._wp ) THEN
                        pfu_ho(ji,jj,jl)=0._wp
                        pfv_ho(ji,jj,jl)=0._wp
                     ENDIF
                     !
                  ENDIF
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'icedyn_adv_umx', pfu_ho, 'U', -1., pfv_ho, 'V', -1. )   ! lateral boundary cond.

      ENDIF

      ! Search local extrema
      ! --------------------
      ! max/min of pt & pt_ups with large negative/positive value (-/+zbig) outside ice cover
      z1_dt = 1._wp / pdt
      DO jl = 1, jpl
         
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF    ( pt(ji,jj,jl) <= 0._wp .AND. pt_ups(ji,jj,jl) <= 0._wp ) THEN
                  zbup(ji,jj) = -zbig
                  zbdo(ji,jj) =  zbig
               ELSEIF( pt(ji,jj,jl) <= 0._wp .AND. pt_ups(ji,jj,jl) > 0._wp ) THEN
                  zbup(ji,jj) = pt_ups(ji,jj,jl)
                  zbdo(ji,jj) = pt_ups(ji,jj,jl)
               ELSEIF( pt(ji,jj,jl) > 0._wp .AND. pt_ups(ji,jj,jl) <= 0._wp ) THEN
                  zbup(ji,jj) = pt(ji,jj,jl)
                  zbdo(ji,jj) = pt(ji,jj,jl)
               ELSE
                  zbup(ji,jj) = MAX( pt(ji,jj,jl) , pt_ups(ji,jj,jl) )
                  zbdo(ji,jj) = MIN( pt(ji,jj,jl) , pt_ups(ji,jj,jl) )
               ENDIF
            END DO
         END DO

         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               !
               zup  = MAX( zbup(ji,jj), zbup(ji-1,jj), zbup(ji+1,jj), zbup(ji,jj-1), zbup(ji,jj+1) )  ! search max/min in neighbourhood
               zdo  = MIN( zbdo(ji,jj), zbdo(ji-1,jj), zbdo(ji+1,jj), zbdo(ji,jj-1), zbdo(ji,jj+1) )
               !
               zpos = MAX( 0._wp, pfu_ho(ji-1,jj  ,jl) ) - MIN( 0._wp, pfu_ho(ji  ,jj  ,jl) ) &  ! positive/negative part of the flux
                  & + MAX( 0._wp, pfv_ho(ji  ,jj-1,jl) ) - MIN( 0._wp, pfv_ho(ji  ,jj  ,jl) )
               zneg = MAX( 0._wp, pfu_ho(ji  ,jj  ,jl) ) - MIN( 0._wp, pfu_ho(ji-1,jj  ,jl) ) &
                  & + MAX( 0._wp, pfv_ho(ji  ,jj  ,jl) ) - MIN( 0._wp, pfv_ho(ji  ,jj-1,jl) )
               !
               zpos = zpos - (pt(ji,jj,jl) * MIN( 0., pu(ji,jj) - pu(ji-1,jj) ) + pt(ji,jj,jl) * MIN( 0., pv(ji,jj) - pv(ji,jj-1) ) &
                  &          ) * ( 1. - pamsk )
               zneg = zneg + (pt(ji,jj,jl) * MAX( 0., pu(ji,jj) - pu(ji-1,jj) ) + pt(ji,jj,jl) * MAX( 0., pv(ji,jj) - pv(ji,jj-1) ) &
                  &          ) * ( 1. - pamsk )
               !
               !                                  ! up & down beta terms
               IF( zpos > 0._wp ) THEN ; zbetup(ji,jj,jl) = MAX( 0._wp, zup - pt_ups(ji,jj,jl) ) / zpos * e1e2t(ji,jj) * z1_dt
               ELSE                    ; zbetup(ji,jj,jl) = 0._wp ! zbig
               ENDIF
               !
               IF( zneg > 0._wp ) THEN ; zbetdo(ji,jj,jl) = MAX( 0._wp, pt_ups(ji,jj,jl) - zdo ) / zneg * e1e2t(ji,jj) * z1_dt
               ELSE                    ; zbetdo(ji,jj,jl) = 0._wp ! zbig
               ENDIF
               !
               ! if all the points are outside ice cover
               IF( zup == -zbig )   zbetup(ji,jj,jl) = 0._wp ! zbig
               IF( zdo ==  zbig )   zbetdo(ji,jj,jl) = 0._wp ! zbig            
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
               zcu = 0.5_wp + SIGN( 0.5_wp , pfu_ho(ji,jj,jl) )
               !
               zcoef = ( zcu * zau + ( 1._wp - zcu ) * zbu )
               !
               pfu_ho(ji,jj,jl) = pfu_ho(ji,jj,jl) * zcoef + pfu_ups(ji,jj,jl)
               !
            END DO
         END DO

         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zav = MIN( 1._wp , zbetdo(ji,jj,jl) , zbetup(ji,jj+1,jl) )
               zbv = MIN( 1._wp , zbetup(ji,jj,jl) , zbetdo(ji,jj+1,jl) )
               zcv = 0.5_wp + SIGN( 0.5_wp , pfv_ho(ji,jj,jl) )
               !
               zcoef = ( zcv * zav + ( 1._wp - zcv ) * zbv )
               !
               pfv_ho(ji,jj,jl) = pfv_ho(ji,jj,jl) * zcoef + pfv_ups(ji,jj,jl)
               !
            END DO
         END DO

         ! clem test
!!         DO jj = 2, jpjm1
!!            DO ji = 2, fs_jpim1   ! vector opt.
!!               zzt = ( pt(ji,jj,jl) - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) ) * pdt * r1_e1e2t(ji,jj)  &
!!                  &                           - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) ) * pdt * r1_e1e2t(ji,jj)  &
!!                  &                     + pt(ji,jj,jl) * pdt * ( pu(ji,jj) - pu(ji-1,jj) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
!!                  &                     + pt(ji,jj,jl) * pdt * ( pv(ji,jj) - pv(ji,jj-1) ) * r1_e1e2t(ji,jj) * (1.-pamsk) &
!!                  &         ) * tmask(ji,jj,1)
!!               IF( zzt < -epsi20 ) THEN
!!                  WRITE(numout,*) 'T<0 nonosc_ice',zzt
!!               ENDIF
!!            END DO
!!         END DO

      END DO
      !
   END SUBROUTINE nonosc_ice

   
   SUBROUTINE limiter_x( pdt, pu, pt, pfu_ups, pfu_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE limiter_x  ***
      !!     
      !! **  Purpose :   compute flux limiter 
      !!----------------------------------------------------------------------
      REAL(wp)                  , INTENT(in   ) ::   pdt          ! tracer time-step
      REAL(wp), DIMENSION(:,:  ), INTENT(in   ) ::   pu           ! ice i-velocity => u*e2
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pt           ! ice tracer
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pfu_ups      ! upstream flux
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pfu_ho       ! high order flux
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

               IF( kn_limiter == 3 ) THEN

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

               ELSEIF( kn_limiter == 2 ) THEN
                  IF( Rj /= 0. ) THEN
                     IF( pu(ji,jj) > 0. ) THEN   ;   Cr = Rjm / Rj
                     ELSE                        ;   Cr = Rjp / Rj
                     ENDIF
                  ELSE
                     Cr = 0.
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

   
   SUBROUTINE limiter_y( pdt, pv, pt, pfv_ups, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE limiter_y  ***
      !!     
      !! **  Purpose :   compute flux limiter 
      !!----------------------------------------------------------------------
      REAL(wp)                   , INTENT(in   ) ::   pdt          ! tracer time-step
      REAL(wp), DIMENSION (:,:  ), INTENT(in   ) ::   pv           ! ice i-velocity => u*e2
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   pt           ! ice tracer
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   pfv_ups      ! upstream flux
      REAL(wp), DIMENSION (:,:,:), INTENT(inout) ::   pfv_ho       ! high order flux
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

               IF( kn_limiter == 3 ) THEN

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

               ELSEIF( kn_limiter == 2 ) THEN

                  IF( Rj /= 0. ) THEN
                     IF( pv(ji,jj) > 0. ) THEN   ;   Cr = Rjm / Rj
                     ELSE                        ;   Cr = Rjp / Rj
                     ENDIF
                  ELSE
                     Cr = 0.
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
