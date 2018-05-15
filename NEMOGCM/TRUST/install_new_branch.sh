#!/bin/bash


cd $( dirname $0 )

## Set defaults
##---------------------------------------------------
TRUS_DBUG=0; TRUS_PROD=0; TRUS_HELP=0


## Get options (for replacing initials settings)
##---------------------------------------------------
while [ $# -ne 0 ]; do
    case $1 in
	'-d'|'--debug'  ) TRUS_DBUG=1 ; shift  ;; '-j'|'--job' 	  ) TRUS_NPRO=$2; shift 2;;
	'-m'|'--machine') TRUS_HPCC=$2; shift 2;; '-h'|'--help'	  ) TRUS_HELP=1 ; shift  ;;
	'-u'|'--user'   ) TRUS_USER=$2; shift 2;; "*"             ) TRUS_HELP=1 ; shift  ;;
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

    printf "\n\nContent of 'config' folder:\n"
    find config -name *.cfg | cut -d/ -f2 \
	| xargs -n 4 printf "\t- %-25s\t- %-25s\t- %-25s\t- %-25s\n"
    exit 1
else
    . ./inc/trusting.env
    [ ${TRUS_DBUG} -eq 1 ] && set -vx
fi


## List last branches from NEMO Forge
##---------------------------------------------------
printf "\nWhat branch do you want to install in "${TRUS_WORK}" for trusting test ? "
echo 'Enter 0 to abort'
select branch in 'trunk' $( svn ls ${TRUS_SVNH}/branches/2015 | tr -d / | sort -r ); do

    if [ $REPLY -eq 0 ]; then exit 1; else export TRUS_BRAN=$branch; fi

    printf "\nInstallation of a working copy of '%s' branch in '%s'? " ${TRUS_BRAN} ${TRUS_WORK}
    printf "\nType [Y|y|yes] to confirm, if not back to branches list number\n"
    read answer
    [[ $answer == 'Y' || $answer == 'y' || $answer == 'yes' ]] && break

done

echo


## First checkout of selected branch
##---------------------------------------------------
echo 'Initial checkout of '${TRUS_BRAN}' branch'
mkdir -p ${TRUS_WORK}/${TRUS_BRAN}/NEMOGCM
cd       ${TRUS_WORK}/${TRUS_BRAN}/NEMOGCM

svn_bran=branches/2015/${TRUS_BRAN}
[ ${TRUS_BRAN} == 'trunk' ] && svn_bran=${TRUS_BRAN}

for elmt in ${TRUS_CKOT}; do
    [ $elmt == '\' ] && continue
    printf "%s " $elmt

    if [ $elmt == 'TOOLS/maketools' ]; then
	svn co -q ${TRUS_SVNH}/${svn_bran}/NEMOGCM/TOOLS --depth empty
	svn up -q $elmt
    else
	svn co -q ${TRUS_SVNH}/${svn_bran}/NEMOGCM/$elmt $elmt
    fi

done

printf "\n\n"


## Compile rebuild_nemo.exe in anticipation
##---------------------------------------------------
if [ $( find ARCH -name arch-${TRUS_HPCC}.fcm ) ]; then
    echo 'Compile NEMO rebuild tool'
    cd TOOLS && ./maketools -n REBUILD_NEMO -m ${TRUS_HPCC} -j ${TRUS_NPRO} >& /dev/null
    [ -e REBUILD_NEMO/rebuild_nemo.exe ] && printf "\033[0;32mOK\033[0m" || printf "\033[0;31mKO\033[0m"
    printf "\n\n"
else
    printf "\033[0;31mNo arch file found to compile NEMO rebuild tool\033[0m\n\n"
fi

exit 0

