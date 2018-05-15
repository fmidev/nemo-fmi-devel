#!/bin/bash

rm -f $( ls -1 tex_main/NEMO_* | egrep -v "\.(bib|cfg|ist|sty|tex)$" )
#rm -rf _minted-*
#rm -rf html*
find tex_* -regextype posix-extended -regex ".*\.(aux|log)$" -exec rm -f {} \;

exit 0
