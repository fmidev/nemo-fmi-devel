#!/bin/sh

./inc/clean.sh
./inc/build.sh

cd tex_main
pdflatex -shell-escape NEMO_manual
cd -

exit 0
