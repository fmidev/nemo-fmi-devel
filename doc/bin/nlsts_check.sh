#!/bin/sh

if [[ $* != '' ]]; then

    if [[ $1 = 'all' ]]; then
        models='NEMO SI3 TOP'
    else
        models=$1
    fi

else
    models='NEMO'
fi

for model in $models; do
    [[ $model =~ ^(SI3|TOP)$ ]] && continue
    echo 'Namelists not included in the '$model' manual:'

    for nlst in $( ls namelists ); do

        if [[ ! $( grep "\\nlst{$nlst}" ./latex/$model/subfiles/*.tex ) ]]; then
            printf "$nlst "
        fi

    done

    echo

done

exit 0
