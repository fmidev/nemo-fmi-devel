#!/bin/sh

echo 'Namelists not included in the manual:'

for nlst in nam*; do
	[[ ! $( grep "forfile{../namelists/$nlst}" ../tex_sub/*.tex ) ]] && printf "$nlst "
done

echo

exit 0
