#!/bin/bash


## Messenger filenames
FILE_DATE=mesg_01_date_$PATTERNAME.txt  ; FILE_RSLT=mesg_02_result_$PATTERNAME.txt
FILE_STAT=mesg_03_state_$PATTERNAME.txt ; FILE_NEMO=mesg_04_nemo_$PATTERNAME.txt
FILE_XIOS=mesg_05_xios_$PATTERNAME.txt  ; FILE_CMPF=mesg_06_compiler_$PATTERNAME.txt
FILE_LMPI=mesg_07_mpi_$PATTERNAME.txt   ; FILE_NCDF=mesg_08_netcdf_$PATTERNAME.txt
FILE_INPT=mesg_09_inputs_$PATTERNAME.txt; FILE_TIME=mesg_10_time_$PATTERNAME.txt
FILE_MEMY=mesg_11_memory_$PATTERNAME.txt; FILE_NOTE=mesg_12_comments_$PATTERNAME.txt

## Trusting timestamped logfile & archive
TRUS_FILE=trusting_${DATE}_$PATTERNAME.txt; TRUS_ARCH=trusting_${DATE}_$PATTERNAME.tgz


## Functions in order of use
print_step() {
    local char_nb=$( echo "$1" | wc -c )
    local outline=$( printf "%${char_nb}s" )

    printf "\nStep.....\n%s\n%s\n" "$1" ${outline// /-}
}

init_files() {
    echo 'Date'        	      > ${FILE_DATE}; echo 'Result'           > ${FILE_RSLT}
    echo 'State'      	      > ${FILE_STAT}; echo 'NEMOGCM rev.'     > ${FILE_NEMO}
    echo 'XIOS rev.'   	      > ${FILE_XIOS}; echo 'Fortran compiler' > ${FILE_CMPF}
    echo 'MPI libs'    	      > ${FILE_LMPI}; echo 'NetCDF libs'      > ${FILE_NCDF}
    echo 'Input files' 	      > ${FILE_INPT}; echo 'Elapsed time'     > ${FILE_TIME}
    echo 'Memory usage (P/V)' > ${FILE_MEMY}; echo 'Comments'         > ${FILE_NOTE}

    ## 'Failed' status with 'Unknown error' by default
    echo ${TRUS_RSLT}      \
	>> ${FILE_RSLT}
    echo 'Unknown error' \
	>> ${FILE_STAT}
}

get_date() {
    ## UTC time zone for timestamping
    local dat=$( date -ud "${DATE}" +"%F %R %Z" )

    echo $dat           \
	>> ${FILE_DATE}
}

get_nemo_rev() {
    local dir rev_loc
    local rev=0

    ## Loop on essential NEMO directories
    for dir in ${TRUS_CKOT} ${TRUS_XIOS}; do

	## For time being, just get revision from XIOS with no action on directory
	if [ $dir == ${TRUS_XIOS} ]; then
	    rev_loc=$( svn info $dir | awk '/Last Changed Rev/ {print $NF}' )
	    echo 'XIOS '${rev_loc} \
		>> model.log
	    echo "<a href=\"https://forge.ipsl.jussieu.fr/ioserver/changeset/${rev_loc}\" target=\"_blank\">${rev_loc}</a>" \
		>> ${FILE_XIOS}
	    continue
	fi

	echo $dir && ${TRUS_SVNA} ${TRUS_NGCM}/$dir
	rev_loc=$( svn info ${TRUS_NGCM}/$dir | awk '/Last Changed Rev/ {print $NF}' )

	## Keep last rev. nb
	[ ${rev_loc} -gt $rev ] && rev=${rev_loc}
    done

    echo 'NEMOGCM '$rev \
	>> model.log
    echo "<a href=\"https://forge.ipsl.jussieu.fr/nemo/changeset/$rev\" target=\"_blank\">$rev</a>" \
	>> ${FILE_NEMO}
}

get_soft_rel() {
    local soft_rel str

    ## Sourcing environment
    if [ -n "${TRUS_ENVI}" ]; then
	if [[  -e ${TRUS_ENVI}.env && $( declare -F | grep ' module' ) ]]; then
            ## .env file if module function is available
	    . ${TRUS_ENVI}.env
	else
            ## .path file if existing, if not the given file
	    [ -e ${TRUS_ENVI}.path ] && . ${TRUS_ENVI}.path || . ${TRUS_ENVI}
	fi
    fi

    ## Problem with `prepend-path` of modulefile that use ':' instead of ' ' as delimiter
    [ $TRUS_HPCC == 'X64_ADA' ] && WRAPPER_LDFLAGS='-L/smplocal/pub/IdrMemMPI/1.4/lib -lidrmem '${WRAPPER_LDFLAGS}

    for str in ${TRUS_CMPV} ${TRUS_MPIR} ${TRUS_CDFR} ${TRUS_CDOR}; do
	[ -z "$str" ] && continue
	soft_rel=''

	## Software release: next word after "$soft" in $PATH (case-insensitive)
	soft_rel=$( echo $PATH | sed "s#.*$str\([0-9.a-z_]*\).*#\1#i" )

	## option --version would work for main compilers (gfortran, intel, pgfortran, ...)
	[ $str == ${TRUS_CMPV} ] && soft_rel=$( $str --version | grep -m1 -oe '\<[0-9. ]*\>' )

	## Cleaning characters string to display proper soft name
	str=$( echo $str | sed 's#\\##g; s#[/-]$##' )

	echo $str ${soft_rel} \
	    >> model.log
    done

    sed -n 3p model.log \
	>> ${FILE_CMPF}
    sed -n 4p model.log \
	>> ${FILE_LMPI}
    sed -n 5p model.log \
	>> ${FILE_NCDF}
}

get_inputs() {
    ## Extract archive or copy files in case of personal inputs
    [ -z "${TRUS_TARF}" ] && get_io="cp ${TRUS_FORC}/* ." || get_io="tar -vxf ${TRUS_FORC}/${TRUS_TARF}"

    ${get_io} > /dev/null
    [ $? -ne 0 ] && get_out 3 || echo 'Success'
    [ $( find -name '*.gz' -print -quit ) ] && gunzip *.gz

    ls -lh > inputs_list.txt
}

diff_inputs() {
    local dif file
    local files_list='' mesg='Same' 

    ## Simple diff
    for file in 'inputs_list.txt' *namelist_* *.xml cpp_*; do
	dif=''

	## Continue even if input file is not in here (see after)
	if [ -e ${TRUS_STOR}/$file ]; then dif=$( diff -q $file ${TRUS_STOR}/$file ); else dif=0; fi

	## Pass over useless file omission in benckmark directory
	[[ -n "$dif" && "$dif" != '0' ]] && { mesg='Different'; echo $dif; files_list+=$file' '; }
    done

    [ $mesg == 'Same' ] && echo $mesg
    echo $mesg          \
	>> ${FILE_INPT}

    ## List different files for web comment
    [ -n "${files_list}" ] && echo 'Inputs  : '${files_list}'differ<br>' \
	>> temp_${FILE_NOTE}
}

job_pending() {
    local outline=$( printf "%100s" ) time_elapsed=0 time_increment=30

    sleep ${time_increment}

    ## Append a log file while pending
    while [[ $( eval ${TRUS_JSTA} ) && ${time_elapsed} -lt ${TRUS_TOUT} ]]; do
	printf "\n%s\n" ${outline// /#}          \
	    >> computation.log
	[ -n "${TRUS_JINF}" ] && eval ${JOB_INFO} \
	    >> computation.log
	sleep ${time_increment}
	time_elapsed=$(( ${time_elapsed} + ${time_increment} ))
    done

    sleep ${time_increment}

    ## Kill remaining job & stop the test if it's too long
    [ ${time_elapsed} -eq ${TRUS_TOUT} ] && { eval ${JOB_DELE} &> /dev/null; get_out 6; }
}

diff_results() {
    local file
    local files_list='' mesg='Same'

    ## Simple diff
    for file in 'ocean.output' *.stat; do
	## Stop if no benchmark files (ocean.output, eventual stat files)
	[ ! -e ${TRUS_STOR}/$file ] && { TRUS_RSLT='FAILED'; get_out 7; }

	diff -q $file ${TRUS_STOR}/$file

	## Continue even if it differs
	[ $? -ne 0 ] && { TRUS_RSLT='FAILED'; mesg='Different'; files_list+=$file' '; }
    done

    [ $mesg == 'Same' ] && echo $mesg

    ## List different files for web comment
    [ -n "${files_list}" ] && echo 'Results : '${files_list}'differ<br>' \
	>> temp_${FILE_NOTE}
}

diff_restart() {
    local base_name comp dif file list_comp list_tmsp nb_dom time_step tmsp
    local files_list='' dif_sum=0

    ## Stop if no benchmark files (ie time.step)
    [ ! -e ${TRUS_STOR}/time.step ] && { TRUS_RSLT='FAILED'; get_out 7; }
    time_step=$( cat ${TRUS_STOR}/time.step | tr -d [:space:] )

    ## Find all restart files to rebuild
    if [ $( find -regex ".*_restart.*[0-9]\.nc" -print -quit ) ]; then
	base_name=$( find -regex ".*_restart.*[0-9]\.nc"                       \
	             | sed "s#^\./\(.*\)_[0-9]*_restart.*#\1#"       | sort -u   )
	list_comp=$( find -regex ".*_restart.*[0-9]\.nc"                       \
	             | sed "s#^.*\(restart[a-z_]*\)_[0-9].*\.nc#\1#" | sort -u   )
	list_tmsp=$( find -regex ".*_restart.*[0-9]\.nc"                       \
	             | sed "s#^.*\([0-9]\{8\}\)_restart.*#\1#"       | sort -u   )

	## Loop on each time step
	for tmsp in ${list_tmsp}; do

	    for comp in ${list_comp}; do
		file=${base_name}_${tmsp}_${comp}
		nb_dom=$( find -name "${file}_[0-9]*.nc" | wc -l | awk '{ print $1 }' )

		if   [ ${nb_dom} -gt 1 ]; then
		    ${TRUS_NGCM}/TOOLS/REBUILD_NEMO/rebuild_nemo -t ${TRUS_NPRO} $file ${nb_dom} \
			> /dev/null

		     ## Possibility of remaining decomposed restarts (even after rebuild)
		    [ $? -eq 0 ] && rm -f ${file}_[0-9]*.nc \
                        > /dev/null

		elif [ ${nb_dom} -eq 0 ]; then
		    TRUS_RSLT='FAILED' && get_out 8
		fi

		## Compare restart files at same time step
		if [ $tmsp -eq ${time_step} ]; then

                    ## Stop if no benchmark files (restart file)
		    if [ -e ${TRUS_STOR}/$file.nc ]; then

	                ## UNIX `cmp` not suitable (timestamp in .nc file)
			dif=$( $TRUS_CDOD $file.nc ${TRUS_STOR}/$file.nc 2> /dev/null          \
			       | awk '/records/ {print $0}' | sed '2 s/^/,/' | tr -d '\n' )

			## CDO can return void stdout with no difference
			if [[ -n "$dif" && $( echo $dif | awk '{print $1}' ) -ne 0 ]]; then
			    TRUS_RSLT='FAILED'
			    files_list+=$comp' ' && let dif_sum+=$( echo $dif | awk '{print $1}' )
			    echo $file.nc': '$dif
			fi

		    else
			TRUS_RSLT='FAILED' && get_out 7
		    fi

		else
		    continue
		fi

	    done

	done

        ## List different files for web comment with sum of different records
	if [ ${dif_sum} -ne 0 ]; then
	    echo 'Restarts: '${files_list}${dif_sum}' record(s) differ<br>' \
		>> temp_${FILE_NOTE}
	else
	    echo 'Same'
	fi

    else
	TRUS_RSLT='FAILED'
    fi

}

get_time() {
    [ -z "${TRUS_JTIM}" ] && return

    ## Interest for checking unusual time computation
    local time_cpu=$( eval ${TRUS_JTIM} )

    printf "Elapsed time: "
    echo ${time_cpu} | tee -a ${FILE_TIME}
}

get_memy() {
    [[ -z "${TRUS_JPME}" && -z "${TRUS_JVME}" ]] && return

    ## Interest for checking unusual memory usage
    local memory_pmax=$( eval ${TRUS_JPME} ) memory_vmax=$( eval ${TRUS_JVME} )

    printf "Memory max usage (physical/virtual): "
    echo ${memory_pmax}' / '${memory_vmax} | tee -a ${FILE_MEMY}
}

comments() {
    local opat
    local line='' state=$1

    if [ -e ocean.output ]; then
        ## 'W A R N I N G' pattern by default
	opat="-A2 \"^ $state\""
	[ "$state" == 'E R R O R' ] && opat="-A4 \"$state\""

        ## Select first occurence for web comment
	line=$( eval grep -m1 $opat ocean.output | tr -d '\n' )
    fi

    [ -n "$line" ] && ( echo $line; printf "$line<br>" \
	>> temp_${FILE_NOTE} )
}

log_make() {
    ## Format comments for web
    [ -e temp_${FILE_NOTE} ] && cat temp_${FILE_NOTE} | tr -d '\n' | sed 's/<br>$//' \
	>> ${FILE_NOTE}

    ## Construct txt file with all messenger files
    paste -d ';' mesg_*.txt | tee ${TRUS_FILE}
}

prod_publish() {
    local cmd
    local rev=$( awk '/NEMOGCM/ {print $NF}' model.log )

    ## Production mode (-p|--prod)
    if [ ${TRUS_PROD} -eq 1 ]; then

	## Create or append trusting logfile
	if [ -f ${TRUS_STOR}/trusting_$PATTERNAME.txt ]; then cmd='tail -1'; else cmd='cat'; fi

	$cmd ${TRUS_FILE}                            \
	    >> ${TRUS_STOR}/trusting_$PATTERNAME.txt

        ## Send mail only when FAILED
	if [[ ! -z "${TRUS_MAIL}" && ${TRUS_RSLT} == 'FAILED' ]]; then

	    ## Content
	    cat <<END_MAIL      \
		> trusting.mail
Dear all,


The trusting sequence has not completed successfully on new configuration ${TRUS_CONF} based on ${TRUS_REFE}.

Here is the model summary:
`cat model.log`

First checking would be on the trusting environment files:
${TRUS_USER}.cfg & ${TRUS_HPCC}.cfg

For more details, look into the testing folder at:
${TRUS_SCRA}

An archive has been created to share the questionable configuration for further studies:
${TRUS_STOR}/${TRUS_ARCH}

END_MAIL

	    ## Send with detailed subject
	    mail -s "[NEMO Trusting][$rev][${TRUS_BRAN}][${TRUS_REFE}] ${TRUS_RSLT} ${TRUS_RORR}" ${TRUS_MAIL} \
		<  trusting.mail
	fi

    fi
}

get_out() {
    local time_step=0

    TRUS_RORR=$1

    printf "\n\nEnd of test\n"

    ## In case of compilation error
    cd ${TRUS_SCRA}

    if [ ${TRUS_RSLT} == 'FAILED' ]; then
	echo 'Failure'

        ## Error identification
	case ${TRUS_RORR} in
            ## Compilation
	    '1') TRUS_RORR='XIOS compilation failed' ;; '2') TRUS_RORR='NEMO compilation failed';;
	    ## Submission
	    '3') TRUS_RORR='Missing input files'     ;; '4') TRUS_RORR='Job submission error'   ;;
	    ## Computation
	    '5') TRUS_RORR='Crashed at time step'    ;; '6') TRUS_RORR='Exceeded time limit'    ;;
	    ## Results
	    '7') TRUS_RORR='Missing previous outputs';; '8') TRUS_RORR='New outputs differ'     ;;
	    ## Other
	    '*') TRUS_RORR='Unknown error'           ;;
	esac

    else
	echo 'Success' && TRUS_RORR='Code is reliable'
    fi

    ## Eventual comments from ocean.output
    if [ "${TRUS_RORR}" == 'Crashed at time step' ]; then
	comments 'E R R O R'
	[ -e time.step ] && time_step=$( grep -o [0-9]* time.step )
	TRUS_RORR+=' '$time_step
    else
	comments 'W A R N I N G'
	[ "${TRUS_RORR}" == 'Exceeded time limit' ] && TRUS_RORR+=' '$(( ${TRUS_TOUT}/3600 ))'h'
    fi

    ## Last messenger files
    #export TRUS_RORR
    sed -i "2 s/.*/$TRUS_RSLT/" ${FILE_RSLT}; sed -i "2 s/.*/$TRUS_RORR/" ${FILE_STAT}

    ## Save tested configuration if trusting failed in production mode (-p|--prod)
    if [[ ${TRUS_RSLT} == 'FAILED' && ${TRUS_PROD} -eq 1 ]]; then
	echo 'Creating archive '${TRUS_ARCH}' under '${TRUS_STOR}
	tar -czf ${TRUS_STOR}/${TRUS_ARCH}               *                    \
	    -C   ${TRUS_NGCM}/CONFIG/${TRUS_CONF}/MY_SRC .                    \
	    -C   ${TRUS_NGCM}/CONFIG/${TRUS_CONF}        cpp_${TRUS_CONF}.fcm
    fi

    ## Logfile construct & eventual sending of notification email
    printf "\nTrusting digest:\n----------------\n"
    log_make
    prod_publish

    exit 0
}
