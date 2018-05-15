#!/bin/bash -f
# set -vx
# simple SETTE report generator.
#
# This version should be run in the SETTE directory. 
# The machine name will be picked up from the sette.sh script but the location of the
# validation directory needs to be set here (currently assumed to reside in the ../CONFIG directory)
#
#########################################################################################
######################### Start of function definitions #################################
##

function restfile() {
# Rebuild ice restart for SAS CONFIG, and restartability checks. Expects LONG and SHORT run directories.
# For Stand Alone Surface configuration ocean is not running, just run ice model; so no outputs ocean files.
# Compares LONG rebuild restart ice file with equivalent entry from the SHORT rebuild restart ice file.
#
# check nco module loaded, and load it if not
if [ ! $( echo $LOADEDMODULES | grep cdo ) ]; then module load cdo >& /dev/null ; fi
#
  vdir=$1
  nam=$2
  pass=$3

#
  if [ -d $vdir/$nam ]; then
    dorv=`ls -1rt $vdir/$nam/$mach/ | tail -1l `
    dorv=`echo $dorv | sed -e 's:.*/::'`
    rep1=`ls -1rt $vdir/$nam/$mach/$dorv/ | tail -2l | head -1 `
    rep2=`ls -1rt $vdir/$nam/$mach/$dorv/ | tail -1l`
    cd ${SAS_RESTART_DIR}/LONG
    #SF add here compilation of rebuild_tools to rebuild restart files, and add comparison of restart files
    cd ${TOOLS_DIR}
    ./maketools -n REBUILD_NEMO -m ${mach} > /dev/null 2>&1
    cd ${TOOLS_DIR}/REBUILD_NEMO
    #SF echo "REBUILD LONG restart SAS files, without standard output"
    ./rebuild_nemo -t 4 ../../CONFIG/ORCA2_SAS_SI3_ST/LONG/SAS_00000240_restart_ice  $NPROC > /dev/null 2>&1
    #SF echo "REBUILD SHORT restart SAS files, without standard output"
    ./rebuild_nemo -t 4 ../../CONFIG/ORCA2_SAS_SI3_ST/SHORT/SAS_00000240_restart_ice $NPROC >&-
    cd ${SAS_RESTART_DIR}/LONG
    #SF echo "COPY rebuild restart files"
    cp SAS_00000240_restart_ice.nc $vdir/$nam/$mach/$dorv/LONG/.
    cp ../SHORT/SAS_00000240_restart_ice.nc $vdir/$nam/$mach/$dorv/SHORT/.

    f1o=$vdir/$nam/$mach/$dorv/LONG/SAS_00000240_restart_ice.nc
    f2o=$vdir/$nam/$mach/$dorv/SHORT/SAS_00000240_restart_ice.nc
    if  [ ! -f $f1o ] &&  [ ! -f $f2o ] ; then
      printf "%-27s %s\n" $nam " REBUILD SAS restart ice DOES NOT exists; incomplete test";
      return;
    fi
    #
    done_oce=0
    #
  if  [  -f $f1o ] && [  -f $f2o ]; then
## Compare the two netcdf files
    cdo diffn $f1o $f2o > cdo_diff.out 2> /dev/null
## Identical if cdo_diff.out exists but has zero size
    if [ ! -s cdo_diff.out ]; then
       difi=0
    else
## Identical if first character of $dif ==0
       dif=$( grep -om1 '[0-9]* of [0-9]* records differ' cdo_diff.out )
# difi contains the first character of summary of cdo dif. if = 0, then 0 record differ between the 2 files    
       if [ -n "$dif" ]; then
           difi=`echo $dif | cut -c -1`
       fi
    fi
    \rm cdo_diff.out 

    if [ $difi == 0 ]; then	
       if [ $pass == 0 ]; then
         printf "%-27s %s %s\n" $nam  " ice restarts are IDENTICAL  passed : " $dorv
       fi
    else
       printf "%-27s %s %s\n" $nam  " ice restarts are DIFFERENT  FAILED : " $dorv 
        #
	# Offer view of differences on the second pass
	#
        if [ $pass == 1 ]; then
          echo "<return> to view restart_ice.nc differences"
          read y
          cdo -diffv $f1o $f2o
          done_oce=1
          #echo "<return> to continue"
          #read y
        fi
    fi
  else
      printf "%-27s %s\n" $nam " incomplete test";
      return;
  fi
#
fi
}

function resttest() { 
#
# Restartability checks. Expects LONG and SHORT run directories
# Compares end of LONG stat files with equivalent entries from the SHORT stat files.
#
  vdir=$1
  nam=$2
  pass=$3
#
  if [ -d $vdir/$nam ]; then
    dorv=`ls -1rt $vdir/$nam/$mach/ | tail -1l `
    dorv=`echo $dorv | sed -e 's:.*/::'`
    f1o=$vdir/$nam/$mach/$dorv/LONG/ocean.output
    f1s=$vdir/$nam/$mach/$dorv/LONG/run.stat
    f1t=$vdir/$nam/$mach/$dorv/LONG/tracer.stat
    f2o=$vdir/$nam/$mach/$dorv/SHORT/ocean.output
    f2s=$vdir/$nam/$mach/$dorv/SHORT/run.stat
    f2t=$vdir/$nam/$mach/$dorv/SHORT/tracer.stat

    if  [ ! -f $f1s ] &&  [ ! -f $f1t ] ; then 
      printf "%-27s %s\n" $nam " incomplete test";
      return; 
    fi
    if  [ ! -f $f2s ] &&  [ ! -f $f2t ] ; then 
      printf "%-27s %s\n" $nam " incomplete test";
      return; 
    fi
#
    done_oce=0

    if  [  -f $f1s ] && [  -f $f2s ]; then 
      nl=(`wc -l $f2s`)
      tail -${nl[0]} $f1s > f1.tmp$$
      cmp -s f1.tmp$$ $f2s
      if [ $? == 0 ]; then
        if [ $pass == 0 ]; then 
          printf "%-27s %s %s\n" $nam  " run.stat    restartability  passed : " $dorv
        fi
      else
        printf "%-27s %s %s\n" $nam  " run.stat    restartability  FAILED : " $dorv 
#
# Offer view of differences on the second pass
#
        if [ $pass == 1 ]; then
          echo "<return> to view run.stat differences"
          read y
          sdiff f1.tmp$$ $f2s
          echo "<return> to view ocean.output differences"
          read y
          sdiff $f1o $f2o | grep "|"
          done_oce=1
          echo "<return> to continue"
          read y
        fi
      fi
    fi
#
# Check tracer.stat files (if they exist)
#
    if  [  -f $f1t ] && [  -f $f2t ]; then
      nl=(`wc -l $f2t`)
      tail -${nl[0]} $f1t > f1.tmp$$
      cmp -s f1.tmp$$ $f2t
      if [ $? == 0 ]; then
        if [ $pass == 0 ]; then 
          printf "%-27s %s %s\n" $nam  " tracer.stat restartability  passed : " $dorv
        fi
      else
        printf "%-27s %s %s\n" $nam  " tracer.stat restartability  FAILED : " $dorv 
#
# Offer view of differences on the second pass
#
        if [ $pass == 1 ]; then
          echo "<return> to view tracer.stat differences"
          read y
          sdiff f1.tmp$$ $f2t
#
# Only offer ocean.output view if it has not been viewed previously
#
          if [ $done_oce == 0 ]; then
            echo "<return> to view ocean.output differences"
            read y
            sdiff $f1o $f2o | grep "|"
          fi
          echo "<return> to continue"
          read y
        fi
      fi
    fi
    rm f1.tmp$$
  fi
}

function reprotest(){
#
# Reproducibility checks. Expects REPRO_N_M and REPRO_I_J run directories
# Compares end of stat files from each
#
  vdir=$1
  nam=$2
  pass=$3
#
  if [ -d $vdir/$nam ]; then
    dorv=`ls -1rt $vdir/$nam/$mach/ | tail -1l `
    dorv=`echo $dorv | sed -e 's:.*/::'`
    rep1=`ls -1rt $vdir/$nam/$mach/$dorv/ | grep REPRO | tail -2l | head -1 `
    rep2=`ls -1rt $vdir/$nam/$mach/$dorv/ | grep REPRO | tail -1l`
    f1o=$vdir/$nam/$mach/$dorv/$rep1/ocean.output
    f1s=$vdir/$nam/$mach/$dorv/$rep1/run.stat
    f1t=$vdir/$nam/$mach/$dorv/$rep1/tracer.stat
    f2o=$vdir/$nam/$mach/$dorv/$rep2/ocean.output
    f2s=$vdir/$nam/$mach/$dorv/$rep2/run.stat
    f2t=$vdir/$nam/$mach/$dorv/$rep2/tracer.stat

    if  [ ! -f $f1s ] && [ ! -f $f1t ] ; then 
      printf "%-27s %s\n" $nam " incomplete test";
      return; 
    fi
    if  [ ! -f $f2s ] && [ ! -f $f2t ] ; then 
      printf "%-27s %s\n" $nam " incomplete test";
      return; 
    fi
#
    done_oce=0

    if  [ -f $f1s ] && [ -f $f2s ] ; then
      cmp -s $f1s $f2s
      if [ $? == 0 ]; then
        if [ $pass == 0 ]; then 
          printf "%-27s %s %s\n" $nam  " run.stat    reproducibility passed : " $dorv
        fi
      else
        printf "%-27s %s %s\n" $nam  " run.stat    reproducibility FAILED : " $dorv 
#
# Offer view of differences on the second pass
#
        if [ $pass == 1 ]; then
          echo "<return> to view run.stat differences"
          read y
          sdiff f1.tmp$$ $f2s
          echo "<return> to view ocean.output differences"
          read y
          sdiff $f1o $f2o | grep "|"
          done_oce=1
          echo "<return> to continue"
          read y
        fi
      fi
    fi
#
# Check tracer.stat files (if they exist)
#
    if  [ -f $f1t ] && [ -f $f2t ] ; then
      cmp -s $f1t $f2t
      if [ $? == 0 ]; then
        if [ $pass == 0 ]; then           printf "%-27s %s %s\n" $nam  " tracer.stat reproducibility passed : " $dorv
        fi
      else
        printf "%-27s %s %s\n" $nam  " tracer.stat reproducibility  FAILED : " $dorv
#
# Offer view of differences on the second pass
#
        if [ $pass == 1 ]; then
          echo "<return> to view tracer.stat differences"
          read y
          sdiff $f1t $f2t
#
# Only offer ocean.output view if it has not been viewed previously
#
          if [ $done_oce == 0 ]; then
            echo "<return> to view ocean.output differences"
            read y
            sdiff $f1o $f2o | grep "|"
          fi
          echo "<return> to continue"
          read y
        fi
      fi
    fi
  fi
}

function runtest(){
#
# Run checks.
# Check presence of E R R O R in ocean.output from each
#
  vdir=$1
  nam=$2
  pass=$3
#
  if [ -d $vdir/$nam ]; then
    dorv=`ls -1rt $vdir/$nam/$mach/ | tail -1l `
    dorv=`echo $dorv | sed -e 's:.*/::'`
    rep1=`ls -1rt $vdir/$nam/$mach/$dorv/ | grep REPRO | tail -2l | head -1 `
    f1o=$vdir/$nam/$mach/$dorv/$rep1/ocean.output
    if  [ ! -f $f1o ] ; then
      printf "%-27s %s %s\n" $nam " ocean.output is MISSING : " $dorv
      return;
    else 
      nerr=`grep 'E R R O R' $f1o | wc -l`
      if [[ $nerr > 0 ]]; then
        printf "%-27s %s %s\n" $nam " run FAILED : " $dorv
        if [ $pass == 1 ]; then
          echo "<return> to view end of ocean.output"
          read y
          tail -100 $f1o
          echo ''
          echo "full ocean.output available here: $f1o"
        fi
      fi
    fi
  else
    printf "%-27s %s %s\n" $nam  " directory is MISSING : " $dorv
  fi
}

function identictest(){
#
#  checks AGRIF does not corrupe results with no AGRIF zoom. Expects ORCA2AGUL/AGRIFNOZ and ORCA2AGUL_NAGR/AGRIFNO  run directories
# Compares solver.stat files for each
#
  vdir=$1
  dir1=$2
  dir2=$3
  pass=$4
#
  if [ -d $vdir/$dir1 ] && [ -d $vdir/$dir2 ]; then
    dorv1=`ls -1rt $vdir/$dir1/$mach/ | tail -1l `
    dorv1=`echo $dorv1 | sed -e 's:.*/::'`
    dorv2=`ls -1rt $vdir/$dir2/$mach/ | tail -1l `
    dorv2=`echo $dorv2 | sed -e 's:.*/::'`

    rep1=`ls -1rt $vdir/$dir1/$mach/$dorv1/ |  tail -1l`
#clem    rep2=`ls -1rt $vdir/$dir2/$mach/$dorv2/ |  tail -1l`
    rep2=`ls -1rt $vdir/$dir1/$mach/$dorv1/ |  tail -1l`
    f1s=$vdir/$dir1/$mach/$dorv1/$rep1/run.stat
    f2s=$vdir/$dir2/$mach/$dorv2/$rep2/run.stat

    if  [ ! -f $f1s ] && [ ! -f $f2s ] ; then 
      printf "%-27s %s\n" $dir1 $dir2 " incomplete test";
      return; 
    fi
#
    done_oce=0

    if  [ -f $f1s ] && [ -f $f2s ] ; then
      cmp -s $f1s $f2s
      if [ $? == 0 ]; then
        if [ $pass == 0 ]; then 
          printf "%-5s %s %-5s %s %s %s\n" $rep1 "AGRIF vs" $rep2 "NOAGRIF run.stat    unchanged  -    passed : " $dorv1 $dorv2
        fi
      else
        printf "%-5s %s %-5s %s %s %s\n" $rep1 "AGRIF vs" $rep2 "NOAGRIF run.stat    changed  -     FAILED : " $dorv1 $dorv2
#
# Offer view of differences on the second pass
#
        if [ $pass == 1 ]; then
          echo "<return> to view run.stat differences"
          read y
          sdiff $f1s $f2s
          done_oce=1
          echo "<return> to continue"
          read y
        fi
      fi
    fi
  else
    printf "%-27s %s\n" $dir1 $dir2 " incomplete test";
  fi
}
########################### END of function definitions #################################
##                                                                                     ##
##    Main script                                                                      ##
##                                                                                     ##
#########################################################################################
#
  mach=`grep "COMPILER=" ./sette.sh | sed -e 's/COMPILER=//'`
  NEMO_VALID=`grep "NEMO_VALIDATION_DIR=" ./param.cfg | sed -e 's/NEMO_VALIDATION_DIR=//'`
  NEMO_VALID=`eval "echo $NEMO_VALID"`
#
  if [ ! -d $NEMO_VALID ]; then
    echo "$NEMO_VALID validation directory not found"
    exit
  fi
#
# Directory to run the tests
  SETTE_DIR=$(cd $(dirname "$0"); pwd)
  MAIN_DIR=$SETTE_DIR/../..
  CONFIG_DIR0=${MAIN_DIR}/CONFIG
  TOOLS_DIR=${MAIN_DIR}/TOOLS
  COMPIL_DIR=${TOOLS_DIR}/COMPILE
  NPROC=32
  SAS_RESTART_DIR=${CONFIG_DIR0}/ORCA2_SAS_SI3_ST
#
# Show current revision tag and branch name
#
cmd="svn"
[ ! -d "$MAIN_DIR/.svn" ] && cmd="git $cmd"
echo $cmd 
echo "$MAIN_DIR/.svn"
lastchange=`$cmd info ${MAIN_DIR} | grep -i "Last Changed Rev:" | sed -e "s/ //g" | cut -d ":" -f 2`
revision=`$cmd info | grep Revision | cut -d ":" -f 2 | tr -d ' '`
branchname=`$cmd info | grep URL | rev | cut -d "/" -f 3 | rev`
echo "SETTE validation report : $branchname @ r$revision  ( last change @ r$lastchange )"
#
# The script also needs the date or revision tag. Currently this is taken from the latest sub-directory found in each directory
#  
for pass in  0 1 
do
#
 if [ $pass == 0 ]; then 
   echo "" 
   echo "!!---------------1st pass------------------!!"
 fi
 if [ $pass == 1 ]; then
    echo ""
    echo "!!---------------2nd pass------------------!!"
 fi
#

# Rebuild and restartability test for SAS
# clem: not needed anymore
# for restart_file in WORCA2_SAS_SI3_ST
# do
#   restfile $NEMO_VALID $restart_file $pass
# done
#
# Restartability test
 echo ""
 echo "   !----restart----!   "
 for restart_test in WGYRE_PISCES_ST WORCA2_SI3_PISCES_ST WORCA2_OFF_PISCES_ST WAMM12_ST WORCA2_SAS_SI3_ST WAGRIF_NORDIC_ST WSPITZ12_ST WISOMIP_ST WOVERFLOW_ST WLOCK_EXCHANGE_ST WVORTEX_ST WWAD_ST WSAS_BIPER_ST 
 do
   resttest $NEMO_VALID $restart_test $pass
 done
#
# Reproducibility tests
 echo ""
 echo "   !----repro----!   "
 for repro_test in WGYRE_PISCES_ST WORCA2_SI3_PISCES_ST WORCA2_OFF_PISCES_ST WAMM12_ST WORCA2_SAS_SI3_ST WORCA2_SI3_OBS_ST WAGRIF_NORDIC_ST WSPITZ12_ST WISOMIP_ST
 do
   reprotest $NEMO_VALID $repro_test $pass
   runtest $NEMO_VALID $repro_test $pass
 done

# AGRIF special check to ensure results are unchanged with and without key_agrif
 echo ""
 echo "   !----agrif check----!   "
 dir1=WAGRIF_NORDIC_NOAGRIF_ST
 dir2=WAGRIF_NORDIC_ST
 identictest $NEMO_VALID $dir1 $dir2 $pass 

done
#

exit
