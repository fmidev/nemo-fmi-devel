MODULE icevar
   !!======================================================================
   !!                       ***  MODULE icevar ***
   !!   sea-ice:     Different sets of ice model variables 
   !!                   how to switch from one to another
   !!
   !!                 There are three sets of variables
   !!                 VGLO : global variables of the model
   !!                        - v_i (jpi,jpj,jpl)
   !!                        - v_s (jpi,jpj,jpl)
   !!                        - a_i (jpi,jpj,jpl)
   !!                        - t_s (jpi,jpj,jpl)
   !!                        - e_i (jpi,jpj,nlay_i,jpl)
   !!                        - sv_i(jpi,jpj,jpl)
   !!                        - oa_i(jpi,jpj,jpl)
   !!                 VEQV : equivalent variables sometimes used in the model
   !!                        - h_i(jpi,jpj,jpl)
   !!                        - h_s(jpi,jpj,jpl)
   !!                        - t_i(jpi,jpj,nlay_i,jpl)
   !!                        ...
   !!                 VAGG : aggregate variables, averaged/summed over all
   !!                        thickness categories
   !!                        - vt_i(jpi,jpj)
   !!                        - vt_s(jpi,jpj)
   !!                        - at_i(jpi,jpj)
   !!                        - et_s(jpi,jpj)  total snow heat content
   !!                        - et_i(jpi,jpj)  total ice thermal content 
   !!                        - sm_i(jpi,jpj)  mean ice salinity
   !!                        - tm_i(jpi,jpj)  mean ice temperature
   !!                        - tm_s(jpi,jpj)  mean snw temperature
   !!======================================================================
   !! History :   -   ! 2006-01 (M. Vancoppenolle) Original code
   !!            3.4  ! 2011-02 (G. Madec) dynamical allocation
   !!            3.5  ! 2012    (M. Vancoppenolle)  add ice_var_itd
   !!            3.6  ! 2014-01 (C. Rousset) add ice_var_zapsmall, rewrite ice_var_itd
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_var_agg       : integrate variables over layers and categories
   !!   ice_var_glo2eqv   : transform from VGLO to VEQV
   !!   ice_var_eqv2glo   : transform from VEQV to VGLO
   !!   ice_var_salprof   : salinity profile in the ice
   !!   ice_var_salprof1d : salinity profile in the ice 1D
   !!   ice_var_zapsmall  : remove very small area and volume
   !!   ice_var_itd       : convert 1-cat to jpl-cat
   !!   ice_var_itd2      : convert N-cat to jpl-cat
   !!   ice_var_bv        : brine volume
   !!   ice_var_enthalpy  : compute ice and snow enthalpies from temperature
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants (ocean directory) 
   USE sbc_oce , ONLY : sss_m
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_var_agg          
   PUBLIC   ice_var_glo2eqv      
   PUBLIC   ice_var_eqv2glo      
   PUBLIC   ice_var_salprof      
   PUBLIC   ice_var_salprof1d    
   PUBLIC   ice_var_zapsmall
   PUBLIC   ice_var_itd
   PUBLIC   ice_var_itd2
   PUBLIC   ice_var_bv           
   PUBLIC   ice_var_enthalpy           

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icevar.F90 8422 2017-08-08 13:58:05Z clem $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_var_agg( kn )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_agg  ***
      !!
      !! ** Purpose :   aggregates ice-thickness-category variables to 
      !!              all-ice variables, i.e. it turns VGLO into VAGG
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kn     ! =1 state variables only
      !                                 ! >1 state variables + others
      !
      INTEGER ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   z1_at_i, z1_vt_i, z1_vt_s
      !!-------------------------------------------------------------------
      !
      !                                      ! integrated values
      vt_i(:,:) =       SUM( v_i(:,:,:)           , dim=3 )
      vt_s(:,:) =       SUM( v_s(:,:,:)           , dim=3 )
      at_i(:,:) =       SUM( a_i(:,:,:)           , dim=3 )
      et_s(:,:)  = SUM( SUM( e_s(:,:,:,:), dim=4 ), dim=3 )
      et_i(:,:)  = SUM( SUM( e_i(:,:,:,:), dim=4 ), dim=3 )
      !
      at_ip(:,:) = SUM( a_ip(:,:,:), dim=3 ) ! melt ponds
      vt_ip(:,:) = SUM( v_ip(:,:,:), dim=3 )
      !
      ato_i(:,:) = 1._wp - at_i(:,:)         ! open water fraction  

      IF( kn > 1 ) THEN
         !
         ALLOCATE( z1_at_i(jpi,jpj) , z1_vt_i(jpi,jpj) , z1_vt_s(jpi,jpj) )
         WHERE( at_i(:,:) > epsi20 )   ;   z1_at_i(:,:) = 1._wp / at_i(:,:)
         ELSEWHERE                     ;   z1_at_i(:,:) = 0._wp
         END WHERE
         WHERE( vt_i(:,:) > epsi20 )   ;   z1_vt_i(:,:) = 1._wp / vt_i(:,:)
         ELSEWHERE                     ;   z1_vt_i(:,:) = 0._wp
         END WHERE
         WHERE( vt_s(:,:) > epsi20 )   ;   z1_vt_s(:,:) = 1._wp / vt_s(:,:)
         ELSEWHERE                     ;   z1_vt_s(:,:) = 0._wp
         END WHERE
         !
         !                          ! mean ice/snow thickness
         hm_i(:,:) = vt_i(:,:) * z1_at_i(:,:)
         hm_s(:,:) = vt_s(:,:) * z1_at_i(:,:)
         !         
         !                          ! mean temperature (K), salinity and age
         tm_su(:,:) = SUM( t_su(:,:,:) * a_i(:,:,:) , dim=3 ) * z1_at_i(:,:)
         tm_si(:,:) = SUM( t_si(:,:,:) * a_i(:,:,:) , dim=3 ) * z1_at_i(:,:)
         om_i (:,:) = SUM( oa_i(:,:,:)              , dim=3 ) * z1_at_i(:,:)
         sm_i (:,:) = SUM( sv_i(:,:,:)              , dim=3 ) * z1_vt_i(:,:)
         !
         tm_i(:,:) = 0._wp
         tm_s(:,:) = 0._wp
         DO jl = 1, jpl
            DO jk = 1, nlay_i
               tm_i(:,:) = tm_i(:,:) + r1_nlay_i * t_i (:,:,jk,jl) * v_i(:,:,jl) * z1_vt_i(:,:)
            END DO
            DO jk = 1, nlay_s
               tm_s(:,:) = tm_s(:,:) + r1_nlay_s * t_s (:,:,jk,jl) * v_s(:,:,jl) * z1_vt_s(:,:)
            END DO
         END DO
         !
         !                           ! put rt0 where there is no ice
         WHERE( at_i(:,:)<=epsi20 )
            tm_su(:,:) = rt0
            tm_si(:,:) = rt0
            tm_i (:,:) = rt0
            tm_s (:,:) = rt0
         END WHERE

         DEALLOCATE( z1_at_i , z1_vt_i , z1_vt_s )
      ENDIF
      !
   END SUBROUTINE ice_var_agg


   SUBROUTINE ice_var_glo2eqv
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_glo2eqv ***
      !!
      !! ** Purpose :   computes equivalent variables as function of  
      !!              global variables, i.e. it turns VGLO into VEQV
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   ze_i             ! local scalars
      REAL(wp) ::   ze_s, ztmelts, zbbb, zccc       !   -      -
      REAL(wp) ::   zhmax, z1_zhmax                 !   -      -
      REAL(wp) ::   zlay_i, zlay_s                  !   -      -
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z1_a_i, z1_v_i
      !!-------------------------------------------------------------------

!!gm Question 2:  It is possible to define existence of sea-ice in a common way between 
!!                ice area and ice volume ?
!!                the idea is to be able to define one for all at the begining of this routine
!!                a criteria for icy area (i.e. a_i > epsi20 and v_i > epsi20 )

      !---------------------------------------------------------------
      ! Ice thickness, snow thickness, ice salinity, ice age and ponds
      !---------------------------------------------------------------
      !                                            !--- inverse of the ice area
      WHERE( a_i(:,:,:) > epsi20 )   ;   z1_a_i(:,:,:) = 1._wp / a_i(:,:,:)
      ELSEWHERE                      ;   z1_a_i(:,:,:) = 0._wp
      END WHERE
      !
      WHERE( v_i(:,:,:) > epsi20 )   ;   z1_v_i(:,:,:) = 1._wp / v_i(:,:,:)
      ELSEWHERE                      ;   z1_v_i(:,:,:) = 0._wp
      END WHERE
      !                                           !--- ice thickness
      h_i(:,:,:) = v_i (:,:,:) * z1_a_i(:,:,:)

      zhmax    =          hi_max(jpl)
      z1_zhmax =  1._wp / hi_max(jpl)               
      WHERE( h_i(:,:,jpl) > zhmax )   ! bound h_i by hi_max (i.e. 99 m) with associated update of ice area
         h_i  (:,:,jpl) = zhmax
         a_i   (:,:,jpl) = v_i(:,:,jpl) * z1_zhmax 
         z1_a_i(:,:,jpl) = zhmax * z1_v_i(:,:,jpl)
      END WHERE
      !                                           !--- snow thickness
      h_s(:,:,:) = v_s (:,:,:) * z1_a_i(:,:,:)
      !                                           !--- ice age      
      o_i(:,:,:) = oa_i(:,:,:) * z1_a_i(:,:,:)
      !                                           !--- pond fraction and thickness      
      a_ip_frac(:,:,:) = a_ip(:,:,:) * z1_a_i(:,:,:)
      WHERE( a_ip_frac(:,:,:) > epsi20 )   ;   h_ip(:,:,:) = v_ip(:,:,:) * z1_a_i(:,:,:) / a_ip_frac(:,:,:)
      ELSEWHERE                            ;   h_ip(:,:,:) = 0._wp
      END WHERE
      !
      !                                           !---  salinity (with a minimum value imposed everywhere)     
      IF( nn_icesal == 2 ) THEN
         WHERE( v_i(:,:,:) > epsi20 )   ;   s_i(:,:,:) = MAX( rn_simin , MIN( rn_simax, sv_i(:,:,:) * z1_v_i(:,:,:) ) )
         ELSEWHERE                      ;   s_i(:,:,:) = rn_simin
         END WHERE
      ENDIF
      CALL ice_var_salprof   ! salinity profile

      !-------------------
      ! Ice temperature   [K]   (with a minimum value (rt0 - 100.))
      !-------------------
      zlay_i   = REAL( nlay_i , wp )    ! number of layers
      DO jl = 1, jpl
         DO jk = 1, nlay_i
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF ( v_i(ji,jj,jl) > epsi20 ) THEN     !--- icy area 
                     !
                     ze_i             =   e_i (ji,jj,jk,jl) * z1_v_i(ji,jj,jl) * zlay_i               ! Energy of melting e(S,T) [J.m-3]
                     ztmelts          = - sz_i(ji,jj,jk,jl) * tmut                                 ! Ice layer melt temperature [C]
                     ! Conversion q(S,T) -> T (second order equation)
                     zbbb             = ( rcp - cpic ) * ztmelts + ze_i * r1_rhoic - lfus
                     zccc             = SQRT( MAX( zbbb * zbbb - 4._wp * cpic * lfus * ztmelts , 0._wp) )
                     t_i(ji,jj,jk,jl) = MAX( -100._wp , MIN( -( zbbb + zccc ) * 0.5_wp * r1_cpic , ztmelts ) ) + rt0   ! [K] with bounds: -100 < t_i < ztmelts
                     !
                  ELSE                                !--- no ice
                     t_i(ji,jj,jk,jl) = rt0
                  ENDIF
               END DO
            END DO
         END DO
      END DO

      !--------------------
      ! Snow temperature   [K]   (with a minimum value (rt0 - 100.))
      !--------------------
      zlay_s = REAL( nlay_s , wp )
      DO jk = 1, nlay_s
         WHERE( v_s(:,:,:) > epsi20 )        !--- icy area
            t_s(:,:,jk,:) = rt0 + MAX( -100._wp ,  &
                 &                MIN( r1_cpic * ( -r1_rhosn * ( e_s(:,:,jk,:) / v_s(:,:,:) * zlay_s ) + lfus ) , 0._wp ) )
         ELSEWHERE                           !--- no ice
            t_s(:,:,jk,:) = rt0
         END WHERE
      END DO
      !
      ! integrated values 
      vt_i (:,:) = SUM( v_i, dim=3 )
      vt_s (:,:) = SUM( v_s, dim=3 )
      at_i (:,:) = SUM( a_i, dim=3 )
      !
   END SUBROUTINE ice_var_glo2eqv


   SUBROUTINE ice_var_eqv2glo
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_eqv2glo ***
      !!
      !! ** Purpose :   computes global variables as function of 
      !!              equivalent variables,  i.e. it turns VEQV into VGLO
      !!-------------------------------------------------------------------
      !
      v_i (:,:,:) = h_i (:,:,:) * a_i (:,:,:)
      v_s (:,:,:) = h_s (:,:,:) * a_i (:,:,:)
      sv_i(:,:,:) = s_i (:,:,:) * v_i (:,:,:)
      v_ip(:,:,:) = h_ip(:,:,:) * a_ip(:,:,:)
      !
   END SUBROUTINE ice_var_eqv2glo


   SUBROUTINE ice_var_salprof
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_salprof ***
      !!
      !! ** Purpose :   computes salinity profile in function of bulk salinity     
      !!
      !! ** Method  : If bulk salinity greater than zsi1, 
      !!              the profile is assumed to be constant (S_inf)
      !!              If bulk salinity lower than zsi0,
      !!              the profile is linear with 0 at the surface (S_zero)
      !!              If it is between zsi0 and zsi1, it is a
      !!              alpha-weighted linear combination of s_inf and s_zero
      !!
      !! ** References : Vancoppenolle et al., 2007
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop index
      REAL(wp) ::   zsal, z1_dS
      REAL(wp) ::   zargtemp , zs0, zs
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   z_slope_s, zalpha    ! case 2 only
      REAL(wp), PARAMETER :: zsi0 = 3.5_wp
      REAL(wp), PARAMETER :: zsi1 = 4.5_wp
      !!-------------------------------------------------------------------

!!gm Question: Remove the option 3 ?  How many years since it last use ? 

      SELECT CASE ( nn_icesal )
      !
      !               !---------------------------------------!
      CASE( 1 )       !  constant salinity in time and space  !
         !            !---------------------------------------!
         sz_i(:,:,:,:) = rn_icesal
         s_i (:,:,:)   = rn_icesal
         !
         !            !---------------------------------------------!
      CASE( 2 )       !  time varying salinity with linear profile  !
         !            !---------------------------------------------!
         !
         ALLOCATE( z_slope_s(jpi,jpj,jpl) , zalpha(jpi,jpj,jpl) )
         !
         DO jl = 1, jpl
            DO jk = 1, nlay_i
               sz_i(:,:,jk,jl)  = s_i(:,:,jl)
            END DO
         END DO
         !                                      ! Slope of the linear profile 
         WHERE( h_i(:,:,:) > epsi20 )   ;   z_slope_s(:,:,:) = 2._wp * s_i(:,:,:) / h_i(:,:,:)
         ELSEWHERE                      ;   z_slope_s(:,:,:) = 0._wp
         END WHERE
         !
         z1_dS = 1._wp / ( zsi1 - zsi0 )
         DO jl = 1, jpl
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zalpha(ji,jj,jl) = MAX(  0._wp , MIN( ( zsi1 - s_i(ji,jj,jl) ) * z1_dS , 1._wp )  )
                  !                             ! force a constant profile when SSS too low (Baltic Sea)
                  IF( 2._wp * s_i(ji,jj,jl) >= sss_m(ji,jj) )   zalpha(ji,jj,jl) = 0._wp  
               END DO
            END DO
         END DO
         !
         ! Computation of the profile
         DO jl = 1, jpl
            DO jk = 1, nlay_i
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     !                          ! linear profile with 0 surface value
                     zs0 = z_slope_s(ji,jj,jl) * ( REAL(jk,wp) - 0.5_wp ) * h_i(ji,jj,jl) * r1_nlay_i
                     zs  = zalpha(ji,jj,jl) * zs0 + ( 1._wp - zalpha(ji,jj,jl) ) * s_i(ji,jj,jl)     ! weighting the profile
                     sz_i(ji,jj,jk,jl) = MIN( rn_simax, MAX( zs, rn_simin ) )
                  END DO
               END DO
            END DO
         END DO
         !
         DEALLOCATE( z_slope_s , zalpha )
         !
         !            !-------------------------------------------!
      CASE( 3 )       ! constant salinity with a fix profile      ! (Schwarzacher (1959) multiyear salinity profile
         !            !-------------------------------------------!                                   (mean = 2.30)
         !
         s_i(:,:,:) = 2.30_wp
!!gm Remark: if we keep the case 3, then compute an store one for all time-step
!!           a array  S_prof(1:nlay_i) containing the calculation and just do:
!         DO jk = 1, nlay_i
!            sz_i(:,:,jk,:) = S_prof(jk)
!         END DO
!!gm end
         !
         DO jl = 1, jpl
            DO jk = 1, nlay_i
               zargtemp  = ( REAL(jk,wp) - 0.5_wp ) * r1_nlay_i
               sz_i(:,:,jk,jl) =  1.6_wp * (  1._wp - COS( rpi * zargtemp**(0.407_wp/(0.573_wp+zargtemp)) )  )
            END DO
         END DO
         !
      END SELECT
      !
   END SUBROUTINE ice_var_salprof


   SUBROUTINE ice_var_salprof1d
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_var_salprof1d  ***
      !!
      !! ** Purpose :   1d computation of the sea ice salinity profile
      !!                Works with 1d vectors and is used by thermodynamic modules
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jk    ! dummy loop indices
      REAL(wp) ::   zargtemp, zsal, z1_dS   ! local scalars
      REAL(wp) ::   zs, zs0              !   -      -
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   z_slope_s, zalpha   !
      REAL(wp), PARAMETER :: zsi0 = 3.5_wp
      REAL(wp), PARAMETER :: zsi1 = 4.5_wp
      !!-------------------------------------------------------------------
      !
      SELECT CASE ( nn_icesal )
      !
      !               !---------------------------------------!
      CASE( 1 )       !  constant salinity in time and space  !
         !            !---------------------------------------!
         sz_i_1d(1:npti,:) = rn_icesal
         !
         !            !---------------------------------------------!
      CASE( 2 )       !  time varying salinity with linear profile  !
         !            !---------------------------------------------!
         !
         ALLOCATE( z_slope_s(jpij), zalpha(jpij) )
         !
         !                                      ! Slope of the linear profile 
         WHERE( h_i_1d(1:npti) > epsi20 )   ;   z_slope_s(1:npti) = 2._wp * s_i_1d(1:npti) / h_i_1d(1:npti)
         ELSEWHERE                          ;   z_slope_s(1:npti) = 0._wp
         END WHERE
         
         z1_dS = 1._wp / ( zsi1 - zsi0 )
         DO ji = 1, npti
            zalpha(ji) = MAX(  0._wp , MIN(  ( zsi1 - s_i_1d(ji) ) * z1_dS , 1._wp  )  )
            !                             ! force a constant profile when SSS too low (Baltic Sea)
            IF( 2._wp * s_i_1d(ji) >= sss_1d(ji) )   zalpha(ji) = 0._wp
         END DO
         !
         ! Computation of the profile
         DO jk = 1, nlay_i
            DO ji = 1, npti
               !                          ! linear profile with 0 surface value
               zs0 = z_slope_s(ji) * ( REAL(jk,wp) - 0.5_wp ) * h_i_1d(ji) * r1_nlay_i
               zs  = zalpha(ji) * zs0 + ( 1._wp - zalpha(ji) ) * s_i_1d(ji)
               sz_i_1d(ji,jk) = MIN( rn_simax , MAX( zs , rn_simin ) )
            END DO
         END DO
         !
         DEALLOCATE( z_slope_s, zalpha )

         !            !-------------------------------------------!
      CASE( 3 )       ! constant salinity with a fix profile      ! (Schwarzacher (1959) multiyear salinity profile
         !            !-------------------------------------------!                                   (mean = 2.30)
         !
         s_i_1d(1:npti) = 2.30_wp
         !
!!gm cf remark in ice_var_salprof routine, CASE( 3 )
         DO jk = 1, nlay_i
            zargtemp  = ( REAL(jk,wp) - 0.5_wp ) * r1_nlay_i
            zsal =  1.6_wp * ( 1._wp - COS( rpi * zargtemp**( 0.407_wp / ( 0.573_wp + zargtemp ) ) ) )
            DO ji = 1, npti
               sz_i_1d(ji,jk) = zsal
            END DO
         END DO
         !
      END SELECT
      !
   END SUBROUTINE ice_var_salprof1d


   SUBROUTINE ice_var_zapsmall
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_zapsmall ***
      !!
      !! ** Purpose :   Remove too small sea ice areas and correct fluxes
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jk   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj) ::   zswitch
      !!-------------------------------------------------------------------
      !
      DO jl = 1, jpl       !==  loop over the categories  ==!
         !
         !-----------------------------------------------------------------
         ! Zap ice energy and use ocean heat to melt ice
         !-----------------------------------------------------------------
         WHERE( a_i(:,:,jl) > epsi10 )   ;   h_i(:,:,jl) = v_i(:,:,jl) / a_i(:,:,jl)
         ELSEWHERE                       ;   h_i(:,:,jl) = 0._wp
         END WHERE
         !
         WHERE( a_i(:,:,jl) < epsi10 .OR. v_i(:,:,jl) < epsi10 .OR. h_i(:,:,jl) < epsi10 )   ;   zswitch(:,:) = 0._wp
         ELSEWHERE                                                                           ;   zswitch(:,:) = 1._wp
         END WHERE
         !
         DO jk = 1, nlay_i
            DO jj = 1 , jpj
               DO ji = 1 , jpi
                  ! update exchanges with ocean
                  hfx_res(ji,jj)   = hfx_res(ji,jj) - (1._wp - zswitch(ji,jj) ) * e_i(ji,jj,jk,jl) * r1_rdtice ! W.m-2 <0
                  e_i(ji,jj,jk,jl) = e_i(ji,jj,jk,jl) * zswitch(ji,jj)
                  t_i(ji,jj,jk,jl) = t_i(ji,jj,jk,jl) * zswitch(ji,jj) + rt0 * ( 1._wp - zswitch(ji,jj) )
               END DO
            END DO
         END DO
         !
         DO jk = 1, nlay_s
            DO jj = 1 , jpj
               DO ji = 1 , jpi
                  ! update exchanges with ocean
                  hfx_res(ji,jj)   = hfx_res(ji,jj) - (1._wp - zswitch(ji,jj) ) * e_s(ji,jj,jk,jl) * r1_rdtice ! W.m-2 <0
                  e_s(ji,jj,jk,jl) = e_s(ji,jj,jk,jl) * zswitch(ji,jj)
                  t_s(ji,jj,jk,jl) = t_s(ji,jj,jk,jl) * zswitch(ji,jj) + rt0 * ( 1._wp - zswitch(ji,jj) )
               END DO
            END DO
         END DO
         !
         DO jj = 1 , jpj
            DO ji = 1 , jpi
               ! update exchanges with ocean
               sfx_res(ji,jj)  = sfx_res(ji,jj) + (1._wp - zswitch(ji,jj) ) * sv_i(ji,jj,jl)   * rhoic * r1_rdtice
               wfx_res(ji,jj)  = wfx_res(ji,jj) + (1._wp - zswitch(ji,jj) ) * v_i (ji,jj,jl)   * rhoic * r1_rdtice
               wfx_res(ji,jj)  = wfx_res(ji,jj) + (1._wp - zswitch(ji,jj) ) * v_s (ji,jj,jl)   * rhosn * r1_rdtice
               !
               !-----------------------------------------------------------------
               ! zap ice and snow volume, add water and salt to ocean
               !-----------------------------------------------------------------
               a_i  (ji,jj,jl) = a_i (ji,jj,jl) * zswitch(ji,jj)
               v_i  (ji,jj,jl) = v_i (ji,jj,jl) * zswitch(ji,jj)
               v_s  (ji,jj,jl) = v_s (ji,jj,jl) * zswitch(ji,jj)
               t_su (ji,jj,jl) = t_su(ji,jj,jl) * zswitch(ji,jj) + t_bo(ji,jj) * ( 1._wp - zswitch(ji,jj) )
               oa_i (ji,jj,jl) = oa_i(ji,jj,jl) * zswitch(ji,jj)
               sv_i (ji,jj,jl) = sv_i(ji,jj,jl) * zswitch(ji,jj)
               !
               h_i (ji,jj,jl) = h_i (ji,jj,jl) * zswitch(ji,jj)
               h_s (ji,jj,jl) = h_s (ji,jj,jl) * zswitch(ji,jj)
               !
               a_ip (ji,jj,jl) = a_ip (ji,jj,jl) * zswitch(ji,jj)
               v_ip (ji,jj,jl) = v_ip (ji,jj,jl) * zswitch(ji,jj)
               !
            END DO
         END DO
         !
      END DO 

      ! to be sure that at_i is the sum of a_i(jl)
      at_i (:,:) = SUM( a_i(:,:,:), dim=3 )
      vt_i (:,:) = SUM( v_i(:,:,:), dim=3 )

      ! open water = 1 if at_i=0
      WHERE( at_i(:,:) == 0._wp )   ato_i(:,:) = 1._wp
      !
   END SUBROUTINE ice_var_zapsmall


   SUBROUTINE ice_var_itd( zhti, zhts, zati, zh_i, zh_s, za_i )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_itd   ***
      !!
      !! ** Purpose :  converting 1-cat ice to multiple ice categories
      !!
      !!                  ice thickness distribution follows a gaussian law
      !!               around the concentration of the most likely ice thickness
      !!                           (similar as iceistate.F90)
      !!
      !! ** Method:   Iterative procedure
      !!                
      !!               1) Try to fill the jpl ice categories (bounds hi_max(0:jpl)) with a gaussian
      !!
      !!               2) Check whether the distribution conserves area and volume, positivity and
      !!                  category boundaries
      !!              
      !!               3) If not (input ice is too thin), the last category is empty and
      !!                  the number of categories is reduced (jpl-1)
      !!
      !!               4) Iterate until ok (SUM(itest(:) = 4)
      !!
      !! ** Arguments : zhti: 1-cat ice thickness
      !!                zhts: 1-cat snow depth
      !!                zati: 1-cat ice concentration
      !!
      !! ** Output    : jpl-cat 
      !!
      !!  (Example of application: BDY forcings when input are cell averaged)  
      !!-------------------------------------------------------------------
      INTEGER  :: ji, jk, jl             ! dummy loop indices
      INTEGER  :: idim, i_fill, jl0  
      REAL(wp) :: zarg, zV, zconv, zdh, zdv
      REAL(wp), DIMENSION(:),   INTENT(in)    ::   zhti, zhts, zati    ! input ice/snow variables
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   zh_i, zh_s, za_i ! output ice/snow variables
      INTEGER , DIMENSION(4)                  ::   itest
      !!-------------------------------------------------------------------
      !
      ! ----------------------------------------
      ! distribution over the jpl ice categories
      ! ----------------------------------------
      ! a gaussian distribution for ice concentration is used
      ! then we check whether the distribution fullfills
      ! volume and area conservation, positivity and ice categories bounds
      idim = SIZE( zhti , 1 )
      zh_i(1:idim,1:jpl) = 0._wp
      zh_s(1:idim,1:jpl) = 0._wp
      za_i(1:idim,1:jpl) = 0._wp
      !
      DO ji = 1, idim
         !
         IF( zhti(ji) > 0._wp ) THEN
            !
            ! find which category (jl0) the input ice thickness falls into
            jl0 = jpl
            DO jl = 1, jpl
               IF ( ( zhti(ji) >= hi_max(jl-1) ) .AND. ( zhti(ji) < hi_max(jl) ) ) THEN
                  jl0 = jl
                  CYCLE
               ENDIF
            END DO
            !
            itest(:) = 0
            i_fill   = jpl + 1                                            !------------------------------------
            DO WHILE ( ( SUM( itest(:) ) /= 4 ) .AND. ( i_fill >= 2 ) )   ! iterative loop on i_fill categories
               !                                                          !------------------------------------
               i_fill = i_fill - 1
               !
               zh_i(ji,1:jpl) = 0._wp
               za_i(ji,1:jpl) = 0._wp
               itest(:)       = 0      
               !
               IF ( i_fill == 1 ) THEN      !-- case very thin ice: fill only category 1
                  zh_i(ji,1) = zhti(ji)
                  za_i (ji,1) = zati (ji)
               ELSE                         !-- case ice is thicker: fill categories >1
                  ! thickness
                  DO jl = 1, i_fill - 1
                     zh_i(ji,jl) = hi_mean(jl)
                  END DO
                  !
                  ! concentration
                  za_i(ji,jl0) = zati(ji) / SQRT(REAL(jpl))
                  DO jl = 1, i_fill - 1
                     IF ( jl /= jl0 ) THEN
                        zarg        = ( zh_i(ji,jl) - zhti(ji) ) / ( zhti(ji) * 0.5_wp )
                        za_i(ji,jl) =   za_i (ji,jl0) * EXP(-zarg**2)
                     ENDIF
                  END DO
                  !
                  ! last category
                  za_i(ji,i_fill) = zati(ji) - SUM( za_i(ji,1:i_fill-1) )
                  zV = SUM( za_i(ji,1:i_fill-1) * zh_i(ji,1:i_fill-1) )
                  zh_i(ji,i_fill) = ( zhti(ji) * zati(ji) - zV ) / MAX( za_i(ji,i_fill), epsi10 ) 
                  !
                  ! correction if concentration of upper cat is greater than lower cat
                  !    (it should be a gaussian around jl0 but sometimes it is not)
                  IF ( jl0 /= jpl ) THEN
                     DO jl = jpl, jl0+1, -1
                        IF ( za_i(ji,jl) > za_i(ji,jl-1) ) THEN
                           zdv = zh_i(ji,jl) * za_i(ji,jl)
                           zh_i(ji,jl    ) = 0._wp
                           za_i (ji,jl    ) = 0._wp
                           za_i (ji,1:jl-1) = za_i(ji,1:jl-1) + zdv / MAX( REAL(jl-1) * zhti(ji), epsi10 )
                        END IF
                     END DO
                  ENDIF
                  !
               ENDIF
               !
               ! Compatibility tests
               zconv = ABS( zati(ji) - SUM( za_i(ji,1:jpl) ) ) 
               IF ( zconv < epsi06 )   itest(1) = 1                                        ! Test 1: area conservation
               !
               zconv = ABS( zhti(ji)*zati(ji) - SUM( za_i(ji,1:jpl)*zh_i(ji,1:jpl) ) )
               IF ( zconv < epsi06 )   itest(2) = 1                                        ! Test 2: volume conservation
               !
               IF ( zh_i(ji,i_fill) >= hi_max(i_fill-1) )   itest(3) = 1                  ! Test 3: thickness of the last category is in-bounds ?
               !
               itest(4) = 1
               DO jl = 1, i_fill
                  IF ( za_i(ji,jl) < 0._wp ) itest(4) = 0                                ! Test 4: positivity of ice concentrations
               END DO
               !                                         !----------------------------
            END DO                                       ! end iteration on categories
            !                                            !----------------------------
         ENDIF
      END DO

      ! Add Snow in each category where za_i is not 0
      DO jl = 1, jpl
         DO ji = 1, idim
            IF( za_i(ji,jl) > 0._wp ) THEN
               zh_s(ji,jl) = zh_i(ji,jl) * ( zhts(ji) / zhti(ji) )
               ! In case snow load is in excess that would lead to transformation from snow to ice
               ! Then, transfer the snow excess into the ice (different from icethd_dh)
               zdh = MAX( 0._wp, ( rhosn * zh_s(ji,jl) + ( rhoic - rau0 ) * zh_i(ji,jl) ) * r1_rau0 ) 
               ! recompute h_i, h_s avoiding out of bounds values
               zh_i(ji,jl) = MIN( hi_max(jl), zh_i(ji,jl) + zdh )
               zh_s(ji,jl) = MAX( 0._wp, zh_s(ji,jl) - zdh * rhoic * r1_rhosn )
            ENDIF
         END DO
      END DO
      !
   END SUBROUTINE ice_var_itd


   SUBROUTINE ice_var_itd2( zhti, zhts, zati, zh_i, zh_s, za_i )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_itd2   ***
      !!
      !! ** Purpose :  converting N-cat ice to jpl ice categories
      !!
      !!                  ice thickness distribution follows a gaussian law
      !!               around the concentration of the most likely ice thickness
      !!                           (similar as iceistate.F90)
      !!
      !! ** Method:   Iterative procedure
      !!                
      !!               1) Fill ice cat that correspond to input thicknesses
      !!                  Find the lowest(jlmin) and highest(jlmax) cat that are filled
      !!
      !!               2) Expand the filling to the cat jlmin-1 and jlmax+1
      !!                   by removing 25% ice area from jlmin and jlmax (resp.) 
      !!              
      !!               3) Expand the filling to the empty cat between jlmin and jlmax 
      !!                   by a) removing 25% ice area from the lower cat (ascendant loop jlmin=>jlmax)
      !!                      b) removing 25% ice area from the higher cat (descendant loop jlmax=>jlmin)
      !!
      !! ** Arguments : zhti: N-cat ice thickness
      !!                zhts: N-cat snow depth
      !!                zati: N-cat ice concentration
      !!
      !! ** Output    : jpl-cat 
      !!
      !!  (Example of application: BDY forcings when inputs have N-cat /= jpl)  
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jl, jl1, jl2             ! dummy loop indices
      INTEGER  ::   idim, icat  
      INTEGER, PARAMETER ::   ztrans = 0.25_wp
      REAL(wp), DIMENSION(:,:), INTENT(in)    ::   zhti, zhts, zati    ! input ice/snow variables
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   zh_i, zh_s, za_i    ! output ice/snow variables
      INTEGER , DIMENSION(:,:), ALLOCATABLE   ::   jlfil, jlfil2
      INTEGER , DIMENSION(:)  , ALLOCATABLE   ::   jlmax, jlmin
      !!-------------------------------------------------------------------
      !
      idim = SIZE( zhti, 1 )
      icat = SIZE( zhti, 2 )
      !
      ALLOCATE( jlfil(idim,jpl), jlfil2(idim,jpl) )       ! allocate arrays
      ALLOCATE( jlmin(idim), jlmax(idim) )

      ! --- initialize output fields to 0 --- !
      zh_i(1:idim,1:jpl) = 0._wp
      zh_s(1:idim,1:jpl) = 0._wp
      za_i(1:idim,1:jpl) = 0._wp
      !
      ! --- fill the categories --- !
      !     find where cat-input = cat-output and fill cat-output fields  
      jlmax(:) = 0
      jlmin(:) = 999
      jlfil(:,:) = 0
      DO jl1 = 1, jpl
         DO jl2 = 1, icat
            DO ji = 1, idim
               IF( hi_max(jl1-1) <= zhti(ji,jl2) .AND. hi_max(jl1) > zhti(ji,jl2) ) THEN
                  ! fill the right category
                  zh_i(ji,jl1) = zhti(ji,jl2)
                  zh_s(ji,jl1) = zhts(ji,jl2)
                  za_i(ji,jl1) = zati(ji,jl2)
                  ! record categories that are filled
                  jlmax(ji) = MAX( jlmax(ji), jl1 )
                  jlmin(ji) = MIN( jlmin(ji), jl1 )
                  jlfil(ji,jl1) = jl1
               ENDIF
            END DO
         END DO
      END DO
      !
      ! --- fill the gaps between categories --- !  
      !     transfer from categories filled at the previous step to the empty ones in between
      DO ji = 1, idim
         jl1 = jlmin(ji)
         jl2 = jlmax(ji)
         IF( jl1 > 1 ) THEN
            ! fill the lower cat (jl1-1)
            za_i(ji,jl1-1) = ztrans * za_i(ji,jl1)
            zh_i(ji,jl1-1) = hi_mean(jl1-1)
            ! remove from cat jl1
            za_i(ji,jl1  ) = ( 1._wp - ztrans ) * za_i(ji,jl1)
         ENDIF
         IF( jl2 < jpl ) THEN
            ! fill the upper cat (jl2+1)
            za_i(ji,jl2+1) = ztrans * za_i(ji,jl2)
            zh_i(ji,jl2+1) = hi_mean(jl2+1)
            ! remove from cat jl2
            za_i(ji,jl2  ) = ( 1._wp - ztrans ) * za_i(ji,jl2)
         ENDIF
      END DO
      !
      jlfil2(:,:) = jlfil(:,:) 
      ! fill categories from low to high
      DO jl = 2, jpl-1
         DO ji = 1, idim
            IF( jlfil(ji,jl-1) /= 0 .AND. jlfil(ji,jl) == 0 ) THEN
               ! fill high
               za_i(ji,jl) = ztrans * za_i(ji,jl-1)
               zh_i(ji,jl) = hi_mean(jl)
               jlfil(ji,jl) = jl
               ! remove low
               za_i(ji,jl-1) = ( 1._wp - ztrans ) * za_i(ji,jl-1)
            ENDIF
         END DO
      END DO
      !
      ! fill categories from high to low
      DO jl = jpl-1, 2, -1
         DO ji = 1, idim
            IF( jlfil2(ji,jl+1) /= 0 .AND. jlfil2(ji,jl) == 0 ) THEN
               ! fill low
               za_i(ji,jl) = za_i(ji,jl) + ztrans * za_i(ji,jl+1)
               zh_i(ji,jl) = hi_mean(jl) 
               jlfil2(ji,jl) = jl
               ! remove high
               za_i(ji,jl+1) = ( 1._wp - ztrans ) * za_i(ji,jl+1)
            ENDIF
         END DO
      END DO
      !
      DEALLOCATE( jlfil, jlfil2 )      ! deallocate arrays
      DEALLOCATE( jlmin, jlmax )
      !
   END SUBROUTINE ice_var_itd2


   SUBROUTINE ice_var_bv
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_bv ***
      !!
      !! ** Purpose :   computes mean brine volume (%) in sea ice
      !!
      !! ** Method  : e = - 0.054 * S (ppt) / T (C)
      !!
      !! References : Vancoppenolle et al., JGR, 2007
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      !!-------------------------------------------------------------------
      !
!!gm I prefere to use WHERE / ELSEWHERE  to set it to zero only where needed   <<<=== to be done
!!   instead of setting everything to zero as just below
      bv_i (:,:,:) = 0._wp
      DO jl = 1, jpl
         DO jk = 1, nlay_i
            WHERE( t_i(:,:,jk,jl) < rt0 - epsi10 )   
               bv_i(:,:,jl) = bv_i(:,:,jl) - tmut * sz_i(:,:,jk,jl) * r1_nlay_i / ( t_i(:,:,jk,jl) - rt0 )
            END WHERE
         END DO
      END DO
      WHERE( vt_i(:,:) > epsi20 )   ;   bvm_i(:,:) = SUM( bv_i(:,:,:) * v_i(:,:,:) , dim=3 ) / vt_i(:,:)
      ELSEWHERE                     ;   bvm_i(:,:) = 0._wp
      END WHERE
      !
   END SUBROUTINE ice_var_bv


   SUBROUTINE ice_var_enthalpy
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_enthalpy *** 
      !!                 
      !! ** Purpose :   Computes sea ice energy of melting q_i (J.m-3) from temperature
      !!
      !! ** Method  :   Formula (Bitz and Lipscomb, 1999)
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jk   ! dummy loop indices
      REAL(wp) ::   ztmelts  ! local scalar 
      !!-------------------------------------------------------------------
      !
      DO jk = 1, nlay_i             ! Sea ice energy of melting
         DO ji = 1, npti
            ztmelts      = - tmut  * sz_i_1d(ji,jk)
            t_i_1d(ji,jk) = MIN( t_i_1d(ji,jk), ztmelts + rt0 ) ! Force t_i_1d to be lower than melting point
                                                                !   (sometimes zdf scheme produces abnormally high temperatures)   
            e_i_1d(ji,jk) = rhoic * ( cpic * ( ztmelts - ( t_i_1d(ji,jk) - rt0 ) )           &
               &                    + lfus * ( 1._wp - ztmelts / ( t_i_1d(ji,jk) - rt0 ) )   &
               &                    - rcp  *   ztmelts )
         END DO
      END DO
      DO jk = 1, nlay_s             ! Snow energy of melting
         DO ji = 1, npti
            e_s_1d(ji,jk) = rhosn * ( cpic * ( rt0 - t_s_1d(ji,jk) ) + lfus )
         END DO
      END DO
      !
   END SUBROUTINE ice_var_enthalpy

#else
   !!----------------------------------------------------------------------
   !!   Default option         Dummy module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icevar
