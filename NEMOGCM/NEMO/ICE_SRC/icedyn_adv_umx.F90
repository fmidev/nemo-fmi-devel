MODULE icedyn_adv_umx
   !!==============================================================================
   !!                       ***  MODULE  icedyn_adv_umx  ***
   !! sea-ice : advection using the ULTIMATE-MACHO scheme
   !!==============================================================================
   !! History :  3.6  !  2014-11  (C. Rousset, G. Madec)  Original code
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv_umx   : update the tracer trend with the 3D advection trends using a TVD scheme
   !!   ultimate_x(_y)    : compute a tracer value at velocity points using ULTIMATE scheme at various orders
   !!   macho             : ???
   !!   nonosc_2d         : compute monotonic tracer fluxes by a non-oscillatory algorithm 
   !!----------------------------------------------------------------------
   USE phycst         ! physical constant
   USE dom_oce        ! ocean domain
   USE sbc_oce , ONLY : nn_fsbc   ! update frequency of surface boundary condition
   USE ice            ! sea-ice variables
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

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2017)
   !! $Id: icedyn_adv_umx.F90 4499 2014-02-18 15:14:31Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_adv_umx( k_order, kt, pu_ice, pv_ice,  &
      &                    pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pe_s, pe_i )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_adv_umx  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of 
      !!                 tracers and add it to the general trend of tracer equations
      !!                 using an "Ultimate-Macho" scheme
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74. 
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   k_order    ! order of the scheme (1-5 or 20)
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
      INTEGER  ::   initad                  ! number of sub-timestep for the advection
      REAL(wp) ::   zcfl , zusnit, zdt      !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zudy, zvdx, zcu_box, zcv_box
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_dyn_adv_umx: Ultimate-Macho advection scheme'
      !
      ALLOCATE( zudy(jpi,jpj) , zvdx(jpi,jpj) , zcu_box(jpi,jpj) , zcv_box(jpi,jpj) )
      !
      ! --- If ice drift field is too fast, use an appropriate time step for advection (CFL test for stability) --- !        
      zcfl  =            MAXVAL( ABS( pu_ice(:,:) ) * rdt_ice * r1_e1u(:,:) )
      zcfl  = MAX( zcfl, MAXVAL( ABS( pv_ice(:,:) ) * rdt_ice * r1_e2v(:,:) ) )
      IF( lk_mpp )   CALL mpp_max( zcfl )

      IF( zcfl > 0.5 ) THEN   ;   initad = 2   ;   zusnit = 0.5_wp
      ELSE                    ;   initad = 1   ;   zusnit = 1.0_wp
      ENDIF

      zdt = rdt_ice / REAL(initad)

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
      DO jt = 1, initad
         CALL adv_umx( k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pato_i(:,:) )             ! Open water area 
         DO jl = 1, jpl
            CALL adv_umx( k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pa_i(:,:,jl) )         ! Ice area
            CALL adv_umx( k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pv_i(:,:,jl) )         ! Ice  volume
            CALL adv_umx( k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, psv_i(:,:,jl) )        ! Salt content
            CALL adv_umx( k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, poa_i(:,:,jl) )        ! Age content
            DO jk = 1, nlay_i
               CALL adv_umx( k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pe_i(:,:,jk,jl) )   ! Ice  heat content
            END DO
            CALL adv_umx( k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pv_s(:,:,jl) )         ! Snow volume
            DO jk = 1, nlay_s
               CALL adv_umx( k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pe_s(:,:,jk,jl) )   ! Snow heat content
            END DO
            IF ( ln_pnd_H12 ) THEN
               CALL adv_umx( k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pa_ip(:,:,jl) )     ! Melt pond fraction
               CALL adv_umx( k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pv_ip(:,:,jl) )     ! Melt pond volume
            ENDIF
         END DO
      END DO
      !
      DEALLOCATE( zudy, zvdx, zcu_box, zcv_box )
      !
   END SUBROUTINE ice_dyn_adv_umx
   
   SUBROUTINE adv_umx( k_order, kt, pdt, puc, pvc, pubox, pvbox, ptc )
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
      INTEGER                     , INTENT(in   ) ::   k_order        ! order of the ULTIMATE scheme
      INTEGER                     , INTENT(in   ) ::   kt             ! number of iteration
      REAL(wp)                    , INTENT(in   ) ::   pdt            ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   puc  , pvc     ! 2 ice velocity components => u*e2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pubox, pvbox   ! upstream velocity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   ptc            ! tracer content field
      !
      INTEGER  ::   ji, jj           ! dummy loop indices  
      REAL(wp) ::   ztra             ! local scalar
      REAL(wp) ::   zfp_ui, zfp_vj   !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj   !   -      -
      REAL(wp), DIMENSION(jpi,jpj) ::   zfu_ups, zfu_ho, zt_u, zt_ups
      REAL(wp), DIMENSION(jpi,jpj) ::   zfv_ups, zfv_ho, zt_v, ztrd
      !!----------------------------------------------------------------------
      !
      !  upstream advection with initial mass fluxes & intermediate update
      ! --------------------------------------------------------------------
      DO jj = 1, jpjm1         ! upstream tracer flux in the i and j direction
         DO ji = 1, fs_jpim1   ! vector opt.
            zfp_ui = puc(ji,jj) + ABS( puc(ji,jj) )
            zfm_ui = puc(ji,jj) - ABS( puc(ji,jj) )
            zfp_vj = pvc(ji,jj) + ABS( pvc(ji,jj) )
            zfm_vj = pvc(ji,jj) - ABS( pvc(ji,jj) )
            zfu_ups(ji,jj) = 0.5_wp * ( zfp_ui * ptc(ji,jj) + zfm_ui * ptc(ji+1,jj  ) )
            zfv_ups(ji,jj) = 0.5_wp * ( zfp_vj * ptc(ji,jj) + zfm_vj * ptc(ji  ,jj+1) )
         END DO
      END DO
      
      DO jj = 2, jpjm1            ! total intermediate advective trends
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ztra = - (   zfu_ups(ji,jj) - zfu_ups(ji-1,jj  )   &
               &       + zfv_ups(ji,jj) - zfv_ups(ji  ,jj-1)   ) * r1_e1e2t(ji,jj)
            !
            ztrd(ji,jj) =                         ztra                         ! upstream trend [ -div(uh) or -div(uhT) ]  
            zt_ups (ji,jj) = ( ptc(ji,jj) + pdt * ztra ) * tmask(ji,jj,1)      ! guess after content field with monotonic scheme
         END DO
      END DO
      CALL lbc_lnk( zt_ups, 'T', 1. )        ! Lateral boundary conditions   (unchanged sign)
      
      ! High order (_ho) fluxes 
      ! -----------------------
      SELECT CASE( k_order )
      CASE ( 20 )                          ! centered second order
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zfu_ho(ji,jj) = 0.5 * puc(ji,jj) * ( ptc(ji,jj) + ptc(ji+1,jj) )
               zfv_ho(ji,jj) = 0.5 * pvc(ji,jj) * ( ptc(ji,jj) + ptc(ji,jj+1) )
            END DO
         END DO
         !
      CASE ( 1:5 )                      ! 1st to 5th order ULTIMATE-MACHO scheme
         CALL macho( k_order, kt, pdt, ptc, puc, pvc, pubox, pvbox, zt_u, zt_v )
         !
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zfu_ho(ji,jj) = puc(ji,jj) * zt_u(ji,jj)
               zfv_ho(ji,jj) = pvc(ji,jj) * zt_v(ji,jj)
            END DO
         END DO
         !
      END SELECT
         
      ! antidiffusive flux : high order minus low order
      ! --------------------------------------------------
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            zfu_ho(ji,jj) = zfu_ho(ji,jj) - zfu_ups(ji,jj)
            zfv_ho(ji,jj) = zfv_ho(ji,jj) - zfv_ups(ji,jj)
         END DO
      END DO
      
      ! monotonicity algorithm
      ! -------------------------
      CALL nonosc_2d( ptc, zfu_ho, zfv_ho, zt_ups, pdt )
      
      ! final trend with corrected fluxes
      ! ------------------------------------
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.  
            ztra       = ztrd(ji,jj)  - (  zfu_ho(ji,jj) - zfu_ho(ji-1,jj  )   &
               &                         + zfv_ho(ji,jj) - zfv_ho(ji  ,jj-1) ) * r1_e1e2t(ji,jj)  
            ptc(ji,jj) = ptc(ji,jj) + pdt * ztra
         END DO
      END DO
      CALL lbc_lnk( ptc, 'T',  1. )
      !
   END SUBROUTINE adv_umx


   SUBROUTINE macho( k_order, kt, pdt, ptc, puc, pvc, pubox, pvbox, pt_u, pt_v )
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
      INTEGER                     , INTENT(in   ) ::   k_order    ! order of the ULTIMATE scheme
      INTEGER                     , INTENT(in   ) ::   kt         ! number of iteration
      REAL(wp)                    , INTENT(in   ) ::   pdt        ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   ptc        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   puc, pvc   ! 2 ice velocity components
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pubox, pvbox   ! upstream velocity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pt_u, pt_v ! tracer at u- and v-points 
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) ::   zc_box    !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: zzt
      !!----------------------------------------------------------------------
      !
      IF( MOD( (kt - 1) / nn_fsbc , 2 ) == 0 ) THEN         !==  odd ice time step:  adv_x then adv_y  ==!
         !
         !                                                           !--  ultimate interpolation of pt at u-point  --!
         CALL ultimate_x( k_order, pdt, ptc, puc, pt_u )
         !
         !                                                           !--  advective form update in zzt  --!
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zzt(ji,jj) = ptc(ji,jj) - pubox(ji,jj) * pdt * ( pt_u(ji,jj) - pt_u(ji-1,jj) ) * r1_e1t(ji,jj)  &
                  &                    - ptc  (ji,jj) * pdt * ( puc (ji,jj) - puc (ji-1,jj) ) * r1_e1e2t(ji,jj)
               zzt(ji,jj) = zzt(ji,jj) * tmask(ji,jj,1)
            END DO
         END DO
         CALL lbc_lnk( zzt, 'T', 1. )
         !
         !                                                           !--  ultimate interpolation of pt at v-point  --!
         CALL ultimate_y( k_order, pdt, zzt, pvc, pt_v )
         !
      ELSE                                                  !==  even ice time step:  adv_y then adv_x  ==!
         !
         !                                                           !--  ultimate interpolation of pt at v-point  --!
         CALL ultimate_y( k_order, pdt, ptc, pvc, pt_v )
         !
         !                                                           !--  advective form update in zzt  --!
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zzt(ji,jj) = ptc(ji,jj) - pvbox(ji,jj) * pdt * ( pt_v(ji,jj) - pt_v(ji,jj-1) ) * r1_e2t(ji,jj)  &
                  &                    - ptc  (ji,jj) * pdt * ( pvc (ji,jj) - pvc (ji,jj-1) ) * r1_e1e2t(ji,jj)
               zzt(ji,jj) = zzt(ji,jj) * tmask(ji,jj,1)
            END DO
         END DO
         CALL lbc_lnk( zzt, 'T', 1. )
         !
         !                                                           !--  ultimate interpolation of pt at u-point  --!
         CALL ultimate_x( k_order, pdt, zzt, puc, pt_u )
         !      
      ENDIF      
      !
   END SUBROUTINE macho


   SUBROUTINE ultimate_x( k_order, pdt, pt, puc, pt_u )
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
      INTEGER                     , INTENT(in   ) ::   k_order   ! ocean time-step index
      REAL(wp)                    , INTENT(in   ) ::   pdt       ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   puc       ! ice i-velocity component
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pt_u      ! tracer at u-point 
      !
      INTEGER  ::   ji, jj       ! dummy loop indices
      REAL(wp) ::   zcu, zdx2, zdx4    !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: ztu1, ztu2, ztu3, ztu4
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
      SELECT CASE (k_order )
      !
      CASE( 1 )                                                   !==  1st order central TIM  ==! (Eq. 21)
         !        
         DO jj = 2, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (                               pt(ji+1,jj) + pt(ji,jj)   &
                  &                                    - SIGN( 1._wp, puc(ji,jj) ) * ( pt(ji+1,jj) - pt(ji,jj) ) )
            END DO
         END DO
         !
      CASE( 2 )                                                   !==  2nd order central TIM  ==! (Eq. 23)
         !
         DO jj = 2, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = puc(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
               pt_u(ji,jj) = 0.5_wp * umask(ji,jj,1) * (                                   pt(ji+1,jj) + pt(ji,jj)   &
                  &                                               -              zcu   * ( pt(ji+1,jj) - pt(ji,jj) ) ) 
            END DO
         END DO
         !  
      CASE( 3 )                                                   !==  3rd order central TIM  ==! (Eq. 24)
         !
         DO jj = 2, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = puc(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
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
         DO jj = 2, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = puc(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
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
         DO jj = 2, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zcu  = puc(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
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
      !
   END SUBROUTINE ultimate_x
   
 
   SUBROUTINE ultimate_y( k_order, pdt, pt, pvc, pt_v )
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
      INTEGER                     , INTENT(in   ) ::   k_order   ! ocean time-step index
      REAL(wp)                    , INTENT(in   ) ::   pdt       ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pvc       ! ice j-velocity component
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pt_v      ! tracer at v-point 
      !
      INTEGER  ::   ji, jj       ! dummy loop indices
      REAL(wp) ::   zcv, zdy2, zdy4    !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: ztv1, ztv2, ztv3, ztv4
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
      SELECT CASE (k_order )
      !
      CASE( 1 )                                                !==  1st order central TIM  ==! (Eq. 21)
         DO jj = 1, jpjm1
            DO ji = fs_2, fs_jpim1
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (                              ( pt(ji,jj+1) + pt(ji,jj) )  &
                  &                                     - SIGN( 1._wp, pvc(ji,jj) ) * ( pt(ji,jj+1) - pt(ji,jj) ) )
            END DO
         END DO
         !
      CASE( 2 )                                                !==  2nd order central TIM  ==! (Eq. 23)
         DO jj = 1, jpjm1
            DO ji = fs_2, fs_jpim1
               zcv  = pvc(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
               pt_v(ji,jj) = 0.5_wp * vmask(ji,jj,1) * (        ( pt(ji,jj+1) + pt(ji,jj) )  &
                  &                                     - zcv * ( pt(ji,jj+1) - pt(ji,jj) ) )
            END DO
         END DO
         CALL lbc_lnk( pt_v, 'V',  1. )
         !
      CASE( 3 )                                                !==  3rd order central TIM  ==! (Eq. 24)
         DO jj = 1, jpjm1
            DO ji = fs_2, fs_jpim1
               zcv  = pvc(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
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
            DO ji = fs_2, fs_jpim1
               zcv  = pvc(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
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
            DO ji = fs_2, fs_jpim1
               zcv  = pvc(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
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
      !
   END SUBROUTINE ultimate_y
   
  
   SUBROUTINE nonosc_2d( pbef, paa, pbb, paft, pdt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc  ***
      !!     
      !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
      !!       scheme and the before field by a nonoscillatory algorithm 
      !!
      !! **  Method  :   ... ???
      !!       warning : pbef and paft must be masked, but the boundaries
      !!       conditions on the fluxes are not necessary zalezak (1979)
      !!       drange (1995) multi-dimensional forward-in-time and upstream-
      !!       in-space based differencing for fluid
      !!----------------------------------------------------------------------
      REAL(wp)                     , INTENT(in   ) ::   pdt          ! tracer time-step
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   pbef, paft   ! before & after field
      REAL(wp), DIMENSION (jpi,jpj), INTENT(inout) ::   paa, pbb     ! monotonic fluxes in the 2 directions
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      INTEGER  ::   ikm1      ! local integer
      REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zsml, z1_dt   ! local scalars
      REAL(wp) ::   zau, zbu, zcu, zav, zbv, zcv, zup, zdo            !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: zbetup, zbetdo, zbup, zbdo, zmsk, zdiv
      !!----------------------------------------------------------------------
      !
      zbig = 1.e+40_wp
      zsml = 1.e-15_wp

      ! test on divergence
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.  
            zdiv(ji,jj) =  - (  paa(ji,jj) - paa(ji-1,jj  )   &
               &              + pbb(ji,jj) - pbb(ji  ,jj-1) )  
         END DO
      END DO
      CALL lbc_lnk( zdiv, 'T', 1. )        ! Lateral boundary conditions   (unchanged sign)

      ! Determine ice masks for before and after tracers 
      WHERE( pbef(:,:) == 0._wp .AND. paft(:,:) == 0._wp .AND. zdiv(:,:) == 0._wp )   ;   zmsk(:,:) = 0._wp
      ELSEWHERE                                                                       ;   zmsk(:,:) = 1._wp * tmask(:,:,1)
      END WHERE

      ! Search local extrema
      ! --------------------
      ! max/min of pbef & paft with large negative/positive value (-/+zbig) inside land
!      zbup(:,:) = MAX( pbef(:,:) * tmask(:,:,1) - zbig * ( 1.e0 - tmask(:,:,1) ),   &
!         &             paft(:,:) * tmask(:,:,1) - zbig * ( 1.e0 - tmask(:,:,1) )  )
!      zbdo(:,:) = MIN( pbef(:,:) * tmask(:,:,1) + zbig * ( 1.e0 - tmask(:,:,1) ),   &
!         &             paft(:,:) * tmask(:,:,1) + zbig * ( 1.e0 - tmask(:,:,1) )  )
      zbup(:,:) = MAX( pbef(:,:) * zmsk(:,:) - zbig * ( 1.e0 - zmsk(:,:) ),   &
         &             paft(:,:) * zmsk(:,:) - zbig * ( 1.e0 - zmsk(:,:) )  )
      zbdo(:,:) = MIN( pbef(:,:) * zmsk(:,:) + zbig * ( 1.e0 - zmsk(:,:) ),   &
         &             paft(:,:) * zmsk(:,:) + zbig * ( 1.e0 - zmsk(:,:) )  )

      z1_dt = 1._wp / pdt
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            !
            zup  = MAX(   zbup(ji,jj), zbup(ji-1,jj  ), zbup(ji+1,jj  ),   &        ! search max/min in neighbourhood
               &                       zbup(ji  ,jj-1), zbup(ji  ,jj+1)    )
            zdo  = MIN(   zbdo(ji,jj), zbdo(ji-1,jj  ), zbdo(ji+1,jj  ),   &
               &                       zbdo(ji  ,jj-1), zbdo(ji  ,jj+1)    )
               !
            zpos = MAX( 0., paa(ji-1,jj  ) ) - MIN( 0., paa(ji  ,jj  ) )   &        ! positive/negative  part of the flux
               & + MAX( 0., pbb(ji  ,jj-1) ) - MIN( 0., pbb(ji  ,jj  ) )
            zneg = MAX( 0., paa(ji  ,jj  ) ) - MIN( 0., paa(ji-1,jj  ) )   &
               & + MAX( 0., pbb(ji  ,jj  ) ) - MIN( 0., pbb(ji  ,jj-1) )
               !
            zbt = e1e2t(ji,jj) * z1_dt                                   ! up & down beta terms
            zbetup(ji,jj) = ( zup         - paft(ji,jj) ) / ( zpos + zsml ) * zbt
            zbetdo(ji,jj) = ( paft(ji,jj) - zdo         ) / ( zneg + zsml ) * zbt
         END DO
      END DO
      CALL lbc_lnk_multi( zbetup, 'T', 1., zbetdo, 'T', 1. )   ! lateral boundary cond. (unchanged sign)

      ! monotonic flux in the i & j direction (paa & pbb)
      ! -------------------------------------
      DO jj = 2, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            zau = MIN( 1._wp , zbetdo(ji,jj) , zbetup(ji+1,jj) )
            zbu = MIN( 1._wp , zbetup(ji,jj) , zbetdo(ji+1,jj) )
            zcu = 0.5  + SIGN( 0.5 , paa(ji,jj) )
            !
            paa(ji,jj) = paa(ji,jj) * ( zcu * zau + ( 1._wp - zcu) * zbu )
         END DO
      END DO
      !
      DO jj = 1, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zav = MIN( 1._wp , zbetdo(ji,jj) , zbetup(ji,jj+1) )
            zbv = MIN( 1._wp , zbetup(ji,jj) , zbetdo(ji,jj+1) )
            zcv = 0.5  + SIGN( 0.5 , pbb(ji,jj) )
            !
            pbb(ji,jj) = pbb(ji,jj) * ( zcv * zav + ( 1._wp - zcv) * zbv )
         END DO
      END DO
      !
   END SUBROUTINE nonosc_2d

#else
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_adv_umx
