#!/bin/sh

./inc/clean.sh
./inc/build.sh

cd tex_main
pdflatex -shell-escape NEMO_manual
mv NEMO_manual.pdf ..
cd -

exit 0
