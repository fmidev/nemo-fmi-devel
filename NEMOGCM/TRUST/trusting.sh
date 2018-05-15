#!/bin/bash


cd $( dirname $0 )

## Set defaults
##---------------------------------------------------
TRUS_DIRE=$PWD
TRUS_DBUG=0; TRUS_PROD=0; TRUS_HELP=0
## No update on SVN directories & 'FAILED' result for 'Unknown error' )
TRUS_SVNA='svn status'; TRUS_RSLT='FAILED'; TRUS_RORR=0
xios_mode='--full'; stdout_redir='>&'
rev=$( svn info | awk '(NR == 9) {print $NF}' )


## Get options (replacing initials settings)
##---------------------------------------------------
while [ $# -ne 0 ]; do
    case $1 in
	'-a'|'--archive') TRUS_TARF=$2; shift 2;; '-b'|'--branch' ) TRUS_BRAN=$2; shift 2;;
	'-d'|'--debug'  ) TRUS_DBUG=1 ; shift  ;; '-e'|'--email'  ) TRUS_MAIL=$2; shift 2;;
	'-f'|'--forcdir') TRUS_FORC=$2; shift 2;; '-j'|'--job'    ) TRUS_NPRO=$2; shift 2;;
	'-h'|'--help'   ) TRUS_HELP=1 ; shift  ;; '-m'|'--machine') TRUS_HPCC=$2; shift 2;;
	'-n'|'--newconf') TRUS_CONF=$2; shift 2;; '-r'|'--refconf') TRUS_REFE=$2; shift 2;;
	'-t'|'--time'   ) TRUS_TOUT=$2; shift 2;; '-p'|'--prod'   ) TRUS_PROD=1 ; shift  ;;
	'-u'|'--user'   ) TRUS_USER=$2; shift 2;; '-v'|'--version') TRUS_SVNV=$2; shift 2;;
	'-w'|'--workdir') TRUS_WORK=$2; shift 2;; "*"             ) TRUS_HELP=1 ; shift  ;;
    esac
done


## Initialization (HPC & user environment)
##---------------------------------------------------
if [[ ! -e config/${TRUS_USER}.cfg || ! -e config/${TRUS_HPCC}.cfg || ${TRUS_HELP} -eq 1 ]]; then
    cat ./inc/trusting_help.txt

    if [ ${TRUS_HELP} -eq 0 ]; then
	printf "\n\n\033[0;33m"
	printf "At least one configuration (arch or user) file is missing or misspelled:"
	printf "\t'%s'.cfg\t'%s'.cfg" ${TRUS_USER} ${TRUS_HPCC}
	printf "\033[0m"
    fi

    printf "\n\nContent of 'config' folder:"
    find config -name *.cfg | cut -d/ -f2 \
	| xargs -n 4 printf "%-25s\t%-25s\t%-25s\n"
    exit 1
else
    . ./inc/trusting.env && . ./inc/trusting_func.sh

    ## DEBUG option to speed up & expand verbosity of compilation
    [ ${TRUS_DBUG} -eq 1 ] && { set -vx; xios_mode=''; stdout_redir='>'; }

    ## If -v|--version option has been set, modify default SVN action on directories
    if   [ $( echo ${TRUS_SVNV} | grep  "HEAD\|up\|update"                     ) ]; then
	TRUS_SVNA='svn update -r HEAD'
    elif [ $( echo ${TRUS_SVNV} | grep -o '{[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}}' ) ]; then
	TRUS_SVNA='svn update -r     '$( echo ${TRUS_SVNV} | grep -o '{[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}}' )
    elif [ $( echo ${TRUS_SVNV} | grep -o '[0-9]*'                             ) ]; then
	TRUS_SVNA='svn update -r     '$( echo ${TRUS_SVNV} | grep -o '[0-9]*' )
    fi

fi


## Display contextual summary of trusting test
##---------------------------------------------------
echo
if [ -t 0 ]; then cat ./inc/banner.txt; else cat ./inc/banner.html; fi
echo
echo '****************************************************************************************************'
echo '*                                                                                                  *'
echo '*                           NEMO Trusting (Continuous Integration Tool)                            *'
echo "*                                             ver.$rev                                             *"
echo '*                                                                                                  *'
echo '****************************************************************************************************'
echo
printf "\t§ Testing configuration\t\t%s based on %s\n" ${TRUS_CONF} ${TRUS_REFE}
printf "\t§ SVN working copy\t\t%s/%s\n"               ${TRUS_WORK} ${TRUS_BRAN}
printf "\t§ Benchmark folder\t\t%s\n"       	       ${TRUS_STOR}
printf "\t§ (Super)Computer\t\t%s\n" 		       ${TRUS_HPCC}
printf "\t§ User installation\t\t%s\n\n" 	       ${TRUS_USER}


## Make timestamped directory with messenger files
##---------------------------------------------------
print_step 'Timestamped testing directory'
mkdir -p ${TRUS_SCRA} ${TRUS_STOR}
cd       ${TRUS_SCRA}
echo     ${TRUS_SCRA}
init_files
get_date


## Get SVN revision on XIOS & NEMO essentials directories
##---------------------------------------------------
print_step "SVN action on NEMO directories: ${TRUS_SVNA}"
get_nemo_rev


## Check softwares versions (after sourced arch environment)
##---------------------------------------------------
print_step 'Get testing environement'
get_soft_rel
cat model.log | awk '{printf "%-20s %s %s\n", $1, $2, $3}'
env | sort > env.log


## XIOS compilation from scratch
##---------------------------------------------------
print_step 'Compile XIOS'
cd ${TRUS_XIOS}
eval ./make_xios ${xios_mode} --arch ${TRUS_HPCC} --job ${TRUS_NPRO} \
    ${stdout_redir} /dev/null
[ ! -e lib/libxios.a ] && get_out 1 || echo 'Success'


## NEMO compilation from scratch
##---------------------------------------------------
print_step "Compile ${TRUS_REFE} configuration"
cd ${TRUS_NGCM}/CONFIG
[[ -d ${TRUS_CONF} && ${TRUS_DBUG} -eq 0 ]] && ./makenemo -n ${TRUS_CONF} clean_config \
    > /dev/null <<EOF
y
EOF

eval ./makenemo -n ${TRUS_CONF} -r ${TRUS_REFE} -m ${TRUS_HPCC} -j ${TRUS_NPRO} \
                ${TRUS_KEYA} ${TRUS_KEYD}                                       \
    ${stdout_redir} /dev/null
[ ! -e ${TRUS_CONF}/BLD/bin/nemo.exe ] && get_out 2 || echo 'Success'


## Get all inputs for running
##---------------------------------------------------
print_step 'Set job (copying or extracting inputs)'
cd ${TRUS_SCRA}
get_inputs
cp   ${TRUS_NGCM}/CONFIG/${TRUS_CONF}/cpp_* .
find ${TRUS_NGCM}/CONFIG/${TRUS_CONF}/EXP00 -regex '.*\(_cfg\|.in\|opa\|_ref\|.xml\)' \
                                            -exec  cp {} . \;


## Check inputs
##---------------------------------------------------
print_step 'Compare inputs'
diff_inputs


## Job submission & computation
##---------------------------------------------------
print_step 'Submit job'
cp ${TRUS_DIRE}/batch/${TRUS_JSPT} ${TRUS_SCRA} ## Copy the submitting script to testing folder
TRUS_JIDN=$( eval ${TRUS_JSUB} )
[ $? -ne 0 ] && get_out 4 || printf "Success (job ID %s)\n" ${TRUS_JIDN}
print_step 'Pending job'
job_pending
print_step 'Job finished'


## Check job state & get computation performances if succeeded
##---------------------------------------------------
print_step 'Test job state'
[[ ! -e time.step || $( grep 'E R R O R' ocean.output ) ]] && get_out 5 || echo 'Success' ## Must be reviewed
print_step 'Get job performances'
get_time
get_memy


## Check outputs
##---------------------------------------------------
TRUS_RSLT='OK' ## 'OK' by default
print_step 'Compare outputs'
diff_results
print_step 'Compare restarts'
diff_restart
[ $TRUS_RSLT == 'FAILED' ] && get_out 8


## End, at least nothing has changed ;-)
##---------------------------------------------------
get_out 0
