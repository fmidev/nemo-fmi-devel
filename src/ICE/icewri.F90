MODULE icewri
   !!======================================================================
   !!                     ***  MODULE  icewri  ***
   !!   sea-ice : output ice variables
   !!======================================================================
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_wri       : write of the diagnostics variables in ouput file 
   !!   ice_wri_state : write for initial state or/and abandon
   !!----------------------------------------------------------------------
   USE dianam         ! build name of file (routine)
   USE phycst         ! physical constant
   USE dom_oce        ! domain: ocean
   USE sbc_oce        ! surf. boundary cond.: ocean
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE ice            ! sea-ice: variables
   USE icevar         ! sea-ice: operations
   !
   USE ioipsl         !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC ice_wri        ! called by ice_stp
   PUBLIC ice_wri_state  ! called by dia_wri_state 

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2017)
   !! $Id: icewri.F90 8409 2017-08-07 15:29:21Z clem $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_wri( kt )
      !!-------------------------------------------------------------------
      !!  This routine computes the average of some variables and write it
      !!  on the ouput files.
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time-step
      !
      INTEGER  ::   ji, jj, jk, jl  ! dummy loop indices
      REAL(wp) ::   z2da, z2db, zrho1, zrho2
      REAL(wp), DIMENSION(jpi,jpj)     ::   z2d !  2D workspace
      REAL(wp), DIMENSION(jpi,jpj)     ::   zmsk00, zmsk05, zmsk15, zmsksn ! O%, 5% and 15% concentration mask and snow mask
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zmsk00l, zmsksnl               ! cat masks
      !
      ! Global ice diagnostics (SIMIP)
      REAL(wp) ::   zdiag_area_nh, zdiag_extt_nh, zdiag_volu_nh   ! area, extent, volume
      REAL(wp) ::   zdiag_area_sh, zdiag_extt_sh, zdiag_volu_sh 
      !!-------------------------------------------------------------------

      IF( ln_timing )   CALL timing_start('icewri')

      !----------------------------------------
      ! Brine volume, switches, missing values
      !----------------------------------------

      CALL ice_var_bv      ! brine volume

      ! tresholds for outputs
      DO jj = 1, jpj
         DO ji = 1, jpi
            zmsk00(ji,jj) = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - epsi06  ) ) ! 1 if ice    , 0 if no ice
            zmsk05(ji,jj) = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - 0.05_wp ) ) ! 1 if 5% ice , 0 if less
            zmsk15(ji,jj) = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - 0.15_wp ) ) ! 1 if 15% ice, 0 if less
            zmsksn(ji,jj) = MAX( 0._wp , SIGN( 1._wp , vt_s(ji,jj) - epsi06  ) ) ! 1 if snow   , 0 if no snow
         END DO
      END DO
      DO jl = 1, jpl
         DO jj = 1, jpj
            DO ji = 1, jpi
               zmsk00l(ji,jj,jl)  = MAX( 0._wp , SIGN( 1._wp , a_i(ji,jj,jl) - epsi06 ) )
               zmsksnl(ji,jj,jl)  = MAX( 0._wp , SIGN( 1._wp , v_s(ji,jj,jl) - epsi06 ) )
            END DO
         END DO
      END DO

      !-----------------
      ! Standard outputs
      !-----------------
      zrho1 = ( rau0 - rhoic ) * r1_rau0; zrho2 = rhosn * r1_rau0
      ! masks
      IF( iom_use('icemask'  ) )   CALL iom_put( "icemask"  , zmsk00              )   ! ice mask 0%
      IF( iom_use('icemask05') )   CALL iom_put( "icemask05", zmsk05              )   ! ice mask 5%
      IF( iom_use('icemask15') )   CALL iom_put( "icemask15", zmsk15              )   ! ice mask 15%
      !
      ! general fields
      IF( iom_use('icemass'  ) )   CALL iom_put( "icemass", rhoic * vt_i * zmsk00 )   ! Ice mass per cell area 
      IF( iom_use('snwmass'  ) )   CALL iom_put( "snwmass", rhosn * vt_s * zmsksn )   ! Snow mass per cell area
      IF( iom_use('icepres'  ) )   CALL iom_put( "icepres", zmsk00                )   ! Ice presence (1 or 0) 
      IF( iom_use('iceconc'  ) )   CALL iom_put( "iceconc", at_i  * zmsk00        )   ! ice concentration
      IF( iom_use('icevolu'  ) )   CALL iom_put( "icevolu", vt_i  * zmsk00        )   ! ice volume = mean ice thickness over the cell
      IF( iom_use('icethic'  ) )   CALL iom_put( "icethic", hm_i  * zmsk00        )   ! ice thickness
      IF( iom_use('snwthic'  ) )   CALL iom_put( "snwthic", hm_s  * zmsk00        )   ! snw thickness
      IF( iom_use('icebrv'   ) )   CALL iom_put( "icebrv" , bvm_i * zmsk00 * 100. )   ! brine volume
      IF( iom_use('iceage'   ) )   CALL iom_put( "iceage" , om_i  * zmsk00 / rday )   ! ice age
      IF( iom_use('icehnew'  ) )   CALL iom_put( "icehnew", ht_i_new              )   ! new ice thickness formed in the leads
      IF( iom_use('snwvolu'  ) )   CALL iom_put( "snwvolu", vt_s  * zmsksn        )   ! snow volume
      IF( iom_use('icefrb') ) THEN
         z2d(:,:) = ( zrho1 * hm_i(:,:) - zrho2 * hm_s(:,:) )                                         
         WHERE( z2d < 0._wp )   z2d = 0._wp
                                   CALL iom_put( "icefrb" , z2d * zmsk00          )   ! Ice freeboard
      ENDIF
      !
      ! melt ponds
      IF( iom_use('iceapnd'  ) )   CALL iom_put( "iceapnd", at_ip  * zmsk00       )   ! melt pond total fraction
      IF( iom_use('icevpnd'  ) )   CALL iom_put( "icevpnd", vt_ip  * zmsk00       )   ! melt pond total volume per unit area
      !
      ! salt
      IF( iom_use('icesalt'  ) )   CALL iom_put( "icesalt", sm_i  * zmsk00        )   ! mean ice salinity
      IF( iom_use('icesalm'  ) )   CALL iom_put( "icesalm", SUM( sv_i, DIM = 3 ) * rhoic * 1.0e-3 * zmsk00 )   ! Mass of salt in sea ice per cell area

      ! heat
      IF( iom_use('icetemp'  ) )   CALL iom_put( "icetemp", ( tm_i  - rt0 ) * zmsk00 )   ! ice mean temperature
      IF( iom_use('snwtemp'  ) )   CALL iom_put( "snwtemp", ( tm_s  - rt0 ) * zmsksn )   ! snw mean temperature
      IF( iom_use('icettop'  ) )   CALL iom_put( "icettop", ( tm_su - rt0 ) * zmsk00 )   ! temperature at the ice surface
      IF( iom_use('icetbot'  ) )   CALL iom_put( "icetbot", ( t_bo  - rt0 ) * zmsk00 )   ! temperature at the ice bottom
      IF( iom_use('icetsni'  ) )   CALL iom_put( "icetsni", ( tm_si - rt0 ) * zmsk00 )   ! temperature at the snow-ice interface
      IF( iom_use('icehc'    ) )   CALL iom_put( "icehc"  ,  -et_i          * zmsk00 )   ! ice heat content
      IF( iom_use('snwhc'    ) )   CALL iom_put( "snwhc"  ,  -et_s          * zmsksn )   ! snow heat content

      ! momentum
      IF( iom_use('uice'     ) )   CALL iom_put( "uice"   , u_ice                 )   ! ice velocity u component
      IF( iom_use('vice'     ) )   CALL iom_put( "vice"   , v_ice                 )   ! ice velocity v component
      IF( iom_use('utau_ai'  ) )   CALL iom_put( "utau_ai", utau_ice * zmsk00     )   ! Wind stress term in force balance (x)
      IF( iom_use('vtau_ai'  ) )   CALL iom_put( "vtau_ai", vtau_ice * zmsk00     )   ! Wind stress term in force balance (y)

      IF( iom_use('icevel') ) THEN 
         DO jj = 2 , jpjm1
            DO ji = 2 , jpim1
               z2da  = ( u_ice(ji,jj) + u_ice(ji-1,jj) )
               z2db  = ( v_ice(ji,jj) + v_ice(ji,jj-1) )
               z2d(ji,jj) = 0.5_wp * SQRT( z2da * z2da + z2db * z2db )
           END DO
         END DO
         CALL lbc_lnk( z2d, 'T', 1. )
         IF( iom_use('icevel') )   CALL iom_put( "icevel" , z2d                   )   ! ice velocity module
      ENDIF

      ! --- category-dependent fields --- !
      IF( iom_use('icemask_cat' ) )   CALL iom_put( "icemask_cat" , zmsk00l                                                    )   ! ice mask 0%
      IF( iom_use('iceconc_cat' ) )   CALL iom_put( "iceconc_cat" , a_i * zmsk00l                                              )   ! area for categories
      IF( iom_use('icethic_cat' ) )   CALL iom_put( "icethic_cat" , h_i * zmsk00l                                              )   ! thickness for categories
      IF( iom_use('snwthic_cat' ) )   CALL iom_put( "snwthic_cat" , h_s * zmsksnl                                              )   ! snow depth for categories
      IF( iom_use('icesalt_cat' ) )   CALL iom_put( "icesalt_cat" , s_i * zmsk00l                                              )   ! salinity for categories
      IF( iom_use('iceage_cat'  ) )   CALL iom_put( "iceage_cat"  , o_i * zmsk00l / rday                                       )   ! ice age
      IF( iom_use('icetemp_cat' ) )   CALL iom_put( "icetemp_cat" , ( SUM( t_i(:,:,:,:), dim=3 ) * r1_nlay_i - rt0 ) * zmsk00l )   ! ice temperature
      IF( iom_use('snwtemp_cat' ) )   CALL iom_put( "snwtemp_cat" , ( SUM( t_s(:,:,:,:), dim=3 ) * r1_nlay_s - rt0 ) * zmsksnl )   ! snow temperature
      IF( iom_use('icettop_cat' ) )   CALL iom_put( "icettop_cat" , ( t_su - rt0 ) * zmsk00l                                   )   ! surface temperature
      IF( iom_use('icebrv_cat'  ) )   CALL iom_put( "icebrv_cat"  ,   bv_i * 100.  * zmsk00l                                   )   ! brine volume
      IF( iom_use('iceapnd_cat' ) )   CALL iom_put( "iceapnd_cat" ,   a_ip         * zmsk00l                                   )   ! melt pond frac for categories
      IF( iom_use('icehpnd_cat' ) )   CALL iom_put( "icehpnd_cat" ,   h_ip         * zmsk00l                                   )   ! melt pond frac for categories
      IF( iom_use('iceafpnd_cat') )   CALL iom_put( "iceafpnd_cat",   a_ip_frac    * zmsk00l                                   )   ! melt pond frac for categories

      !------------------
      ! Add-ons for SIMIP
      !------------------
      ! trends
      IF( iom_use('dmithd') )   CALL iom_put( "dmithd", - wfx_bog - wfx_bom - wfx_sum - wfx_sni - wfx_opw - wfx_lam - wfx_res ) ! Sea-ice mass change from thermodynamics
      IF( iom_use('dmidyn') )   CALL iom_put( "dmidyn", - wfx_dyn + rhoic * diag_trp_vi     )   ! Sea-ice mass change from dynamics(kg/m2/s)
      IF( iom_use('dmiopw') )   CALL iom_put( "dmiopw", - wfx_opw                           )   ! Sea-ice mass change through growth in open water
      IF( iom_use('dmibog') )   CALL iom_put( "dmibog", - wfx_bog                           )   ! Sea-ice mass change through basal growth
      IF( iom_use('dmisni') )   CALL iom_put( "dmisni", - wfx_sni                           )   ! Sea-ice mass change through snow-to-ice conversion
      IF( iom_use('dmisum') )   CALL iom_put( "dmisum", - wfx_sum                           )   ! Sea-ice mass change through surface melting
      IF( iom_use('dmibom') )   CALL iom_put( "dmibom", - wfx_bom                           )   ! Sea-ice mass change through bottom melting
      IF( iom_use('dmtsub') )   CALL iom_put( "dmtsub", - wfx_sub                           )   ! Sea-ice mass change through evaporation and sublimation
      IF( iom_use('dmssub') )   CALL iom_put( "dmssub", - wfx_snw_sub                       )   ! Snow mass change through sublimation
      IF( iom_use('dmisub') )   CALL iom_put( "dmisub", - wfx_ice_sub                       )   ! Sea-ice mass change through sublimation
      IF( iom_use('dmsspr') )   CALL iom_put( "dmsspr", - wfx_spr                           )   ! Snow mass change through snow fall
      IF( iom_use('dmsssi') )   CALL iom_put( "dmsssi",   wfx_sni*rhosn*r1_rhoic            )   ! Snow mass change through snow-to-ice conversion
      IF( iom_use('dmsmel') )   CALL iom_put( "dmsmel", - wfx_snw_sum                       )   ! Snow mass change through melt
      IF( iom_use('dmsdyn') )   CALL iom_put( "dmsdyn", - wfx_snw_dyn + rhosn * diag_trp_vs )   ! Snow mass change through dynamics(kg/m2/s)

      ! Global ice diagnostics
      IF( iom_use('NH_icearea') .OR. iom_use('NH_icevolu') .OR. iom_use('NH_iceextt') )   THEN   ! NH diagnostics
         !
         WHERE( ff_t > 0._wp )   ;   zmsk00(:,:) = 1.0e-12
         ELSEWHERE               ;   zmsk00(:,:) = 0.
         END WHERE 
         zdiag_area_nh = glob_sum( at_i(:,:) * zmsk00(:,:) * e1e2t(:,:) )
         zdiag_volu_nh = glob_sum( vt_i(:,:) * zmsk00(:,:) * e1e2t(:,:) )
         !
         WHERE( ff_t > 0._wp .AND. at_i > 0.15 )   ; zmsk00(:,:) = 1.0e-12
         ELSEWHERE                                 ; zmsk00(:,:) = 0.
         END WHERE 
         zdiag_extt_nh = glob_sum( zmsk00(:,:) * e1e2t(:,:) )
         !
         IF( iom_use('NH_icearea') )   CALL iom_put( "NH_icearea" ,  zdiag_area_nh )
         IF( iom_use('NH_icevolu') )   CALL iom_put( "NH_icevolu" ,  zdiag_volu_nh )
         IF( iom_use('NH_iceextt') )   CALL iom_put( "NH_iceextt" ,  zdiag_extt_nh )
         !
      ENDIF
      !
      IF( iom_use('SH_icearea') .OR. iom_use('SH_icevolu') .OR. iom_use('SH_iceextt') )   THEN   ! SH diagnostics
         !
         WHERE( ff_t < 0._wp ); zmsk00(:,:) = 1.0e-12; 
         ELSEWHERE            ; zmsk00(:,:) = 0.
         END WHERE 
         zdiag_area_sh = glob_sum( at_i(:,:) * zmsk00(:,:) * e1e2t(:,:) ) 
         zdiag_volu_sh = glob_sum( vt_i(:,:) * zmsk00(:,:) * e1e2t(:,:) )
         !
         WHERE( ff_t < 0._wp .AND. at_i > 0.15 ); zmsk00(:,:) = 1.0e-12
         ELSEWHERE                              ; zmsk00(:,:) = 0.
         END WHERE 
         zdiag_extt_sh = glob_sum( zmsk00(:,:) * e1e2t(:,:) )
         !
         IF( iom_use('SH_icearea') ) CALL iom_put( "SH_icearea", zdiag_area_sh )
         IF( iom_use('SH_icevolu') ) CALL iom_put( "SH_icevolu", zdiag_volu_sh )
         IF( iom_use('SH_iceextt') ) CALL iom_put( "SH_iceextt", zdiag_extt_sh )
         !
      ENDIF 
      !
!!CR      !     !  Create an output files (output.lim.abort.nc) if S < 0 or u > 20 m/s
!!CR      !     IF( kindic < 0 )   CALL ice_wri_state( 'output.abort' )
!!CR      !     not yet implemented
!!gm  idem for the ocean...  Ask Seb how to get read of ioipsl....
      !
      IF( ln_timing )  CALL timing_stop('icewri')
      !
   END SUBROUTINE ice_wri

 
   SUBROUTINE ice_wri_state( kt, kid, kh_i )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE ice_wri_state  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ice state and forcing fields for ice model
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! History :   4.0  !  2013-06  (C. Rousset)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt               ! ocean time-step index
      INTEGER, INTENT( in ) ::   kid , kh_i
      INTEGER               ::   nz_i, jl
      REAL(wp), DIMENSION(jpl) ::   jcat
      !!----------------------------------------------------------------------
      !
      DO jl = 1, jpl
         jcat(jl) = REAL(jl)
      END DO
      
      CALL histvert( kid, "ncatice", "Ice Categories","", jpl, jcat, nz_i, "up")

      CALL histdef( kid, "sithic", "Ice thickness"          , "m"      , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "siconc", "Ice concentration"      , "%"      , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "sitemp", "Ice temperature"        , "C"      , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "sivelu", "i-Ice speed "           , "m/s"    , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "sivelv", "j-Ice speed "           , "m/s"    , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "sistru", "i-Wind stress over ice" , "Pa"     , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "sistrv", "j-Wind stress over ice" , "Pa"     , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "sisflx", "Solar flx over ocean"   , "W/m2"   , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "sinflx", "NonSolar flx over ocean", "W/m2"   , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "snwpre", "Snow precipitation"     , "kg/m2/s", jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "sisali", "Ice salinity"           , "PSU"    , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "sivolu", "Ice volume"             , "m"      , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "sidive", "Ice divergence"         , "10-8s-1", jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "si_amp", "Melt pond fraction"     , "%"      , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "si_vmp", "Melt pond volume"       ,  "m"     , jpi,jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      !
      CALL histdef( kid, "sithicat", "Ice thickness"        , "m"      , jpi,jpj, kh_i, jpl, 1, jpl, nz_i, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "siconcat", "Ice concentration"    , "%"      , jpi,jpj, kh_i, jpl, 1, jpl, nz_i, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "sisalcat", "Ice salinity"         , ""       , jpi,jpj, kh_i, jpl, 1, jpl, nz_i, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "snthicat", "Snw thickness"        , "m"      , jpi,jpj, kh_i, jpl, 1, jpl, nz_i, 32, "inst(x)", rdt, rdt )

      CALL histend( kid, snc4set )   ! end of the file definition

      CALL histwrite( kid, "sithic", kt, hm_i          , jpi*jpj, (/1/) )    
      CALL histwrite( kid, "siconc", kt, at_i          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "sitemp", kt, tm_i - rt0    , jpi*jpj, (/1/) )
      CALL histwrite( kid, "sivelu", kt, u_ice         , jpi*jpj, (/1/) )
      CALL histwrite( kid, "sivelv", kt, v_ice         , jpi*jpj, (/1/) )
      CALL histwrite( kid, "sistru", kt, utau_ice      , jpi*jpj, (/1/) )
      CALL histwrite( kid, "sistrv", kt, vtau_ice      , jpi*jpj, (/1/) )
      CALL histwrite( kid, "sisflx", kt, qsr           , jpi*jpj, (/1/) )
      CALL histwrite( kid, "sinflx", kt, qns           , jpi*jpj, (/1/) )
      CALL histwrite( kid, "snwpre", kt, sprecip       , jpi*jpj, (/1/) )
      CALL histwrite( kid, "sisali", kt, sm_i          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "sivolu", kt, vt_i          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "sidive", kt, divu_i*1.0e8  , jpi*jpj, (/1/) )
      CALL histwrite( kid, "si_amp", kt, at_ip         , jpi*jpj, (/1/) )
      CALL histwrite( kid, "si_vmp", kt, vt_ip         , jpi*jpj, (/1/) )
      !
      CALL histwrite( kid, "sithicat", kt, h_i         , jpi*jpj*jpl, (/1/) )    
      CALL histwrite( kid, "siconcat", kt, a_i         , jpi*jpj*jpl, (/1/) )    
      CALL histwrite( kid, "sisalcat", kt, s_i         , jpi*jpj*jpl, (/1/) )    
      CALL histwrite( kid, "snthicat", kt, h_s         , jpi*jpj*jpl, (/1/) )    

      !! The file is closed in dia_wri_state (ocean routine)
      !! CALL histclo( kid )
      !
    END SUBROUTINE ice_wri_state

#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icewri
