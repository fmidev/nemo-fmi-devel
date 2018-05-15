MODULE stpctl
   !!======================================================================
   !!                       ***  MODULE  stpctl  ***
   !! Ocean run control :  gross check of the ocean time stepping
   !!======================================================================
   !! History :  OPA  ! 1991-03  (G. Madec) Original code
   !!            6.0  ! 1992-06  (M. Imbard)
   !!            8.0  ! 1997-06  (A.M. Treguier)
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            2.0  ! 2009-07  (G. Madec)  Add statistic for time-spliting
   !!            3.7  ! 2016-09  (G. Madec)  Remove solver
   !!            4.0  ! 2017-04  (G. Madec)  regroup global communications
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp_ctl      : Control the run
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE c1d             ! 1D vertical configuration
   USE diawri          ! Standard run outputs       (dia_wri_state routine)
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing
   USE wet_dry,   ONLY : ll_wd, ssh_ref    ! reference depth for negative bathy

   USE netcdf          ! NetCDF library
   IMPLICIT NONE
   PRIVATE

   PUBLIC stp_ctl           ! routine called by step.F90

   INTEGER  ::   idrun, idtime, idssh, idu, ids, istatus
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL licence (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp_ctl( kt, kindic )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl  ***
      !!                     
      !! ** Purpose :   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Print it each 50 time steps
      !!              - Stop the run IF problem encountered by setting indic=-3
      !!                Problems checked: |ssh| maximum larger than 10 m
      !!                                  |U|   maximum larger than 10 m/s 
      !!                                  negative sea surface salinity
      !!
      !! ** Actions :   "time.step" file = last ocean time-step
      !!                "run.stat"  file = run statistics
      !!                nstop indicator sheared among all local domain (lk_mpp=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER, INTENT(inout) ::   kindic   ! error indicator
      !!
      INTEGER  ::   ji, jj, jk             ! dummy loop indices
      INTEGER  ::   iih, ijh               ! local integers
      INTEGER  ::   iiu, iju, iku          !   -       -
      INTEGER  ::   iis, ijs, iks          !   -       -
      REAL(wp) ::   zzz                    ! local real 
      INTEGER , DIMENSION(3) ::   ilocu, ilocs
      INTEGER , DIMENSION(2) ::   iloch
      REAL(wp), DIMENSION(4) ::   zmax
      CHARACTER(len=20) :: clname
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'stp_ctl : time-stepping control'
         WRITE(numout,*) '~~~~~~~'
         !                                ! open time.step file
         CALL ctl_opn( numstp, 'time.step', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
         !                                ! open run.stat file
         CALL ctl_opn( numrun, 'run.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )

         IF( lwm ) THEN
            clname = 'run.stat.nc'
            IF( .NOT. Agrif_Root() )   clname = TRIM(Agrif_CFixed())//"_"//TRIM(clname)
            istatus = NF90_CREATE( TRIM(clname), NF90_CLOBBER, idrun )
            istatus = NF90_DEF_DIM( idrun, 'time', NF90_UNLIMITED, idtime )
            istatus = NF90_DEF_VAR( idrun, 'abs_ssh_max', NF90_DOUBLE, (/ idtime /), idssh )
            istatus = NF90_DEF_VAR( idrun,   'abs_u_max', NF90_DOUBLE, (/ idtime /), idu )
            istatus = NF90_DEF_VAR( idrun,       's_min', NF90_DOUBLE, (/ idtime /), ids )
            istatus = NF90_ENDDEF(idrun)
         ENDIF
         
      ENDIF
      !
      IF(lwp) THEN                        !==  current time step  ==!   ("time.step" file)
         WRITE ( numstp, '(1x, i8)' )   kt
         REWIND( numstp )
      ENDIF
      !
      !                                   !==  test of extrema  ==!
      IF( ll_wd ) THEN
         zmax(1) = MAXVAL(  ABS( sshn(:,:) + ssh_ref*tmask(:,:,1) )  )        ! ssh max 
      ELSE
         zmax(1) = MAXVAL(  ABS( sshn(:,:) )  )                               ! ssh max
      ENDIF
      zmax(2) = MAXVAL(  ABS( un(:,:,:) )  )                                  ! velocity max (zonal only)
      zmax(3) = MAXVAL( -tsn(:,:,:,jp_sal) , mask = tmask(:,:,:) == 1._wp )   ! minus salinity max
      zmax(4) = REAL( nstop , wp )                                            ! stop indicator
      !
      IF( lk_mpp ) THEN
         CALL mpp_max_multiple( zmax(:), 4 )    ! max over the global domain
         !
         nstop = INT( zmax(4) )                 ! nstop indicator sheared among all local domains
      ENDIF
      !
      IF( MOD( kt, nwrite ) == 1 .AND. lwp ) THEN
         WRITE(numout,*) ' ==>> time-step= ', kt, ' |ssh| max: ',   zmax(1), ' |U| max: ', zmax(2),   &
            &                                     ' S min: '    , - zmax(3)
      ENDIF
      !
      IF (  zmax(1) >  15._wp .OR.   &                     ! too large sea surface height ( > 10 m)
         &  zmax(2) >  10._wp .OR.   &                     ! too large velocity ( > 10 m/s)
         &  zmax(3) >=  0._wp .OR.   &                     ! negative or zero sea surface salinity
         &  ISNAN( zmax(1) + zmax(2) + zmax(3) )  ) THEN   ! NaN encounter in the tests
         IF( lk_mpp ) THEN
            CALL mpp_maxloc( ABS(sshn)        , ssmask(:,:)  , zzz, iih, ijh )
            CALL mpp_maxloc( ABS(un)          , umask (:,:,:), zzz, iiu, iju, iku )
            CALL mpp_minloc( tsn(:,:,:,jp_sal), tmask (:,:,:), zzz, iis, ijs, iks )
         ELSE
            iloch = MINLOC( ABS( sshn(:,:)   )                               )
            ilocu = MAXLOC( ABS( un  (:,:,:) )                               )
            ilocs = MINLOC( tsn(:,:,:,jp_sal) , mask = tmask(:,:,:) == 1._wp )
            iih = iloch(1) + nimpp - 1   ;   ijh = iloch(2) + njmpp - 1
            iiu = ilocu(1) + nimpp - 1   ;   iju = ilocu(2) + njmpp - 1   ;   iku = ilocu(3)
            iis = ilocs(1) + nimpp - 1   ;   ijs = ilocs(2) + njmpp - 1   ;   iks = ilocu(3)
         ENDIF
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) ' stp_ctl: |ssh| > 10 m   or   |U| > 10 m/s   or   S < 0   or   NaN encounter in the tests'
            WRITE(numout,*) ' ======= '
            WRITE(numout,9100) kt,   zmax(1), iih, ijh
            WRITE(numout,9200) kt,   zmax(2), iiu, iju, iku
            WRITE(numout,9300) kt, - zmax(3), iis, ijs, iks
            WRITE(numout,*)
            WRITE(numout,*) '          output of last computed fields in output.abort.nc file'
         ENDIF
         kindic = -3
         !
         nstop = nstop + 1                            ! increase nstop by 1 (on all local domains)
         CALL dia_wri_state( 'output.abort', kt )     ! create an output.abort file
         !
      ENDIF
9100  FORMAT (' kt=',i8,'   |ssh| max: ',1pg11.4,', at  i j  : ',2i5)
9200  FORMAT (' kt=',i8,'   |U|   max: ',1pg11.4,', at  i j k: ',3i5)
9300  FORMAT (' kt=',i8,'   S     min: ',1pg11.4,', at  i j  : ',2i5)
      !
      !                                            !==  run statistics  ==!   ("run.stat" file)
      IF(lwp) WRITE(numrun,9400) kt, zmax(1), zmax(2), - zmax(3)
      IF( lwm ) THEN
         istatus = NF90_PUT_VAR( idrun, idssh, (/ zmax(1)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,   idu, (/ zmax(2)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,   ids, (/-zmax(3)/), (/kt/), (/1/) )
         IF( MOD( kt , 100 ) == 0 ) istatus = NF90_SYNC(idrun)
         IF( kt == nitend         ) istatus = NF90_CLOSE(idrun)
      END IF
      !
9400  FORMAT(' it :', i8, '    |ssh|_max: ', D23.16, ' |U|_max: ', D23.16,' S_min: ', D23.16)
      !
   END SUBROUTINE stp_ctl

   !!======================================================================
END MODULE stpctl
