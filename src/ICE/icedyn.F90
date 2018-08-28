MODULE icedyn
   !!======================================================================
   !!                     ***  MODULE  icedyn  ***
   !!   Sea-Ice dynamics : master routine for sea ice dynamics 
   !!======================================================================
   !! history :  4.0  ! 2018  (C. Rousset)  original code SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn       : dynamics of sea ice
   !!   ice_dyn_init  : initialization and namelist read
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE ice            ! sea-ice: variables
   USE icedyn_rhg     ! sea-ice: rheology
   USE icedyn_adv     ! sea-ice: advection
   USE icedyn_rdgrft  ! sea-ice: ridging/rafting
   USE icecor         ! sea-ice: corrections
   USE icevar         ! sea-ice: operations
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn        ! called by icestp.F90
   PUBLIC   ice_dyn_init   ! called by icestp.F90
   
   INTEGER ::              nice_dyn   ! choice of the type of dynamics
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_dynFULL    = 1   ! full ice dynamics               (rheology + advection + ridging/rafting + correction)
   INTEGER, PARAMETER ::   np_dynRHGADV  = 2   ! pure dynamics                   (rheology + advection) 
   INTEGER, PARAMETER ::   np_dynADV     = 3   ! only advection w prescribed vel.(rn_uvice + advection)
   !
   ! ** namelist (namdyn) **
   LOGICAL  ::   ln_dynFULL       ! full ice dynamics               (rheology + advection + ridging/rafting + correction)
   LOGICAL  ::   ln_dynRHGADV     ! no ridge/raft & no corrections  (rheology + advection)
   LOGICAL  ::   ln_dynADV        ! only advection w prescribed vel.(rn_uvice + advection)
   REAL(wp) ::   rn_uice          !    prescribed u-vel (case np_dynADV)
   REAL(wp) ::   rn_vice          !    prescribed v-vel (case np_dynADV)
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icedyn.F90 8378 2017-07-26 13:55:59Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn( kt )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE ice_dyn  ***
      !!               
      !! ** Purpose :   this routine manages sea ice dynamics
      !!
      !! ** Action : - calculation of friction in case of landfast ice
      !!             - call ice_dyn_rhg    = rheology
      !!             - call ice_dyn_adv    = advection
      !!             - call ice_dyn_rdgrft = ridging/rafting
      !!             - call ice_cor        = corrections if fields are out of bounds
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ice time step
      !!
      INTEGER ::   ji, jj, jl         ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zhmax
      !!--------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('icedyn')
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_dyn: sea-ice dynamics'
         WRITE(numout,*)'~~~~~~~'
      ENDIF

      !                      
      IF( ln_landfast ) THEN            !-- Landfast ice parameterization: define max bottom friction
         tau_icebfr(:,:) = 0._wp
         DO jl = 1, jpl
            WHERE( h_i(:,:,jl) > ht_n(:,:) * rn_gamma )   tau_icebfr(:,:) = tau_icebfr(:,:) + a_i(:,:,jl) * rn_icebfr
         END DO
         IF( iom_use('tau_icebfr') )   CALL iom_put( 'tau_icebfr', tau_icebfr )  
      ENDIF

      zhmax(:,:,:) = h_i_b(:,:,:)      !-- Record max of the surrounding 9-pts ice thick. (for CALL Hbig)
      DO jl = 1, jpl
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
!!gm use of MAXVAL here is very probably less efficient than expending the 9 values
               zhmax(ji,jj,jl) = MAX( epsi20, MAXVAL( h_i_b(ji-1:ji+1,jj-1:jj+1,jl) ) )
            END DO
         END DO
      END DO
      CALL lbc_lnk( zhmax(:,:,:), 'T', 1. )
      !
      !
      SELECT CASE( nice_dyn )           !-- Set which dynamics is running

      CASE ( np_dynFULL )          !==  all dynamical processes  ==!
         CALL ice_dyn_rhg   ( kt )                            ! -- rheology  
         CALL ice_dyn_adv   ( kt )   ;   CALL Hbig( zhmax )   ! -- advection of ice + correction on ice thickness
         CALL ice_dyn_rdgrft( kt )                            ! -- ridging/rafting 
         CALL ice_cor       ( kt , 1 )                        ! -- Corrections

      CASE ( np_dynRHGADV  )       !==  no ridge/raft & no corrections ==!
         CALL ice_dyn_rhg   ( kt )                            ! -- rheology  
         CALL ice_dyn_adv   ( kt )                            ! -- advection of ice
         CALL Hpiling                                         ! -- simple pile-up (replaces ridging/rafting)

      CASE ( np_dynADV )           !==  pure advection ==!   (prescribed velocities)
         u_ice(:,:) = rn_uice * umask(:,:,1)
         v_ice(:,:) = rn_vice * vmask(:,:,1)
         !!CALL RANDOM_NUMBER(u_ice(:,:))
         !!CALL RANDOM_NUMBER(v_ice(:,:))
         CALL ice_dyn_adv   ( kt )                            ! -- advection of ice

      END SELECT
      !
      IF( ln_timing )   CALL timing_stop('icedyn')
      !
   END SUBROUTINE ice_dyn


   SUBROUTINE Hbig( phmax )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE Hbig  ***
      !!
      !! ** Purpose : Thickness correction in case advection scheme creates
      !!              abnormally tick ice
      !!
      !! ** Method  : 1- check whether ice thickness resulting from advection is
      !!                 larger than the surrounding 9-points before advection
      !!                 and reduce it if a) divergence or b) convergence & at_i>0.8
      !!              2- bound ice thickness with hi_max (99m)
      !!
      !! ** input   : Max thickness of the surrounding 9-points
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   phmax   ! max ice thick from surrounding 9-pts
      !
      INTEGER  ::   ji, jj, jl         ! dummy loop indices
      REAL(wp) ::   zh, zdv
      !!-------------------------------------------------------------------
      !
      CALL ice_var_zapsmall                       !-- zap small areas
      !
      DO jl = 1, jpl
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF ( v_i(ji,jj,jl) > 0._wp ) THEN  !-- bound to hmax
                  !
                  zh  = v_i (ji,jj,jl) / a_i(ji,jj,jl)
                  zdv = v_i(ji,jj,jl) - v_i_b(ji,jj,jl)  
                  !
                  IF ( ( zdv >  0.0 .AND. zh > phmax(ji,jj,jl) .AND. at_i_b(ji,jj) < 0.80 ) .OR. &
                     & ( zdv <= 0.0 .AND. zh > phmax(ji,jj,jl) ) ) THEN
                     a_i (ji,jj,jl) = v_i(ji,jj,jl) / MIN( phmax(ji,jj,jl), hi_max(jpl) )   !-- bound h_i to hi_max (99 m)
                  ENDIF
                  !
               ENDIF
            END DO
         END DO
      END DO 
      !                                           !-- correct pond fraction to avoid a_ip > a_i
      WHERE( a_ip(:,:,:) > a_i(:,:,:) )   a_ip(:,:,:) = a_i(:,:,:)
      !
   END SUBROUTINE Hbig


   SUBROUTINE Hpiling
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE Hpiling  ***
      !!
      !! ** Purpose : Simple conservative piling comparable with 1-cat models
      !!
      !! ** Method  : pile-up ice when no ridging/rafting
      !!
      !! ** input   : a_i
      !!-------------------------------------------------------------------
      INTEGER ::   jl         ! dummy loop indices
      !!-------------------------------------------------------------------
      !
      CALL ice_var_zapsmall                       !-- zap small areas
      !
      at_i(:,:) = SUM( a_i(:,:,:), dim=3 )
      DO jl = 1, jpl
         WHERE( at_i(:,:) > epsi20 )
            a_i(:,:,jl) = a_i(:,:,jl) * (  1._wp + MIN( rn_amax_2d(:,:) - at_i(:,:) , 0._wp ) / at_i(:,:)  )
         END WHERE
      END DO
      !
   END SUBROUTINE Hpiling


   SUBROUTINE ice_dyn_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_init  ***
      !!
      !! ** Purpose : Physical constants and parameters linked to the ice
      !!      dynamics
      !!
      !! ** Method  :  Read the namdyn namelist and check the ice-dynamic
      !!       parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !!
      NAMELIST/namdyn/ ln_dynFULL, ln_dynRHGADV, ln_dynADV, rn_uice, rn_vice,  &
         &             rn_ishlat  ,                                            &
         &             ln_landfast, rn_gamma , rn_icebfr, rn_lfrelax
      !!-------------------------------------------------------------------
      !
      REWIND( numnam_ice_ref )         ! Namelist namdyn in reference namelist : Ice dynamics
      READ  ( numnam_ice_ref, namdyn, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn in reference namelist', lwp )
      REWIND( numnam_ice_cfg )         ! Namelist namdyn in configuration namelist : Ice dynamics
      READ  ( numnam_ice_cfg, namdyn, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn in configuration namelist', lwp )
      IF(lwm) WRITE( numoni, namdyn )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn:'
         WRITE(numout,*) '      Full ice dynamics      (rhg + adv + ridge/raft + corr)  ln_dynFULL   = ', ln_dynFULL
         WRITE(numout,*) '      No ridge/raft & No cor (rhg + adv)                      ln_dynRHGADV = ', ln_dynRHGADV
         WRITE(numout,*) '      Advection only         (rn_uvice + adv)                 ln_dynADV    = ', ln_dynADV
         WRITE(numout,*) '         with prescribed velocity given by   (u,v)_ice = (rn_uice,rn_vice) = (', rn_uice,',', rn_vice,')'
         WRITE(numout,*) '      lateral boundary condition for sea ice dynamics         rn_ishlat    = ', rn_ishlat
         WRITE(numout,*) '      Landfast: param (T or F)                                ln_landfast  = ', ln_landfast
         WRITE(numout,*) '         fraction of ocean depth that ice must reach          rn_gamma     = ', rn_gamma
         WRITE(numout,*) '         maximum bottom stress per unit area of contact       rn_icebfr    = ', rn_icebfr
         WRITE(numout,*) '         relax time scale (s-1) to reach static friction      rn_lfrelax   = ', rn_lfrelax
         WRITE(numout,*)
      ENDIF
      !                             !== set the choice of ice dynamics ==!
      ioptio = 0 
      !      !--- full dynamics                               (rheology + advection + ridging/rafting + correction)
      IF( ln_dynFULL   ) THEN   ;   ioptio = ioptio + 1   ;   nice_dyn = np_dynFULL      ;   ENDIF
      !      !--- dynamics without ridging/rafting and corr   (rheology + advection)
      IF( ln_dynRHGADV ) THEN   ;   ioptio = ioptio + 1   ;   nice_dyn = np_dynRHGADV    ;   ENDIF
      !      !--- advection only with prescribed ice velocities (from namelist)
      IF( ln_dynADV    ) THEN   ;   ioptio = ioptio + 1   ;   nice_dyn = np_dynADV       ;   ENDIF
      !
      IF( ioptio /= 1 )    CALL ctl_stop( 'ice_dyn_init: one and only one ice dynamics option has to be defined ' )
      !
      !                                      !--- Lateral boundary conditions
      IF     (      rn_ishlat == 0.                ) THEN   ;   IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  free-slip'
      ELSEIF (      rn_ishlat == 2.                ) THEN   ;   IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  no-slip'
      ELSEIF ( 0. < rn_ishlat .AND. rn_ishlat < 2. ) THEN   ;   IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  partial-slip'
      ELSEIF ( 2. < rn_ishlat                      ) THEN   ;   IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  strong-slip'
      ENDIF
      !                                      !--- NO Landfast ice : set to zero once for all
      IF( .NOT.ln_landfast )   tau_icebfr(:,:) = 0._wp
      !
      CALL ice_dyn_rdgrft_init          ! set ice ridging/rafting parameters
      CALL ice_dyn_rhg_init             ! set ice rheology parameters
      CALL ice_dyn_adv_init             ! set ice advection parameters
      !
   END SUBROUTINE ice_dyn_init

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif 

   !!======================================================================
END MODULE icedyn
