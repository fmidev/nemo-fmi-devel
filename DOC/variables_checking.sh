#!/bin/sh

#set -euvx

#shared_namelist_files=$( ls ../NEMOGCM/CONFIG/SHARED/namelist_* )
shared_namelist_files=$( ls ~/Workspace/NEMO/branches/trunk/NEMOGCM/CONFIG/SHARED/namelist_* )

namelist_variables=$( awk '($1 ~ /[cdlnrsy]n_/ ) {print $1}' ${shared_namelist_files} \
								| sed 's/[=(].*//' | sort -u												)

for variable in ${namelist_variables}; do
	[[ $( grep ${variable/_/\\\\\\\\_} tex_sub/*.tex ) ]] && echo $variable ${variable//_/\\\\\\\\_} && exit 0
done

exit 0
