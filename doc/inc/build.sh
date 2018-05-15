#!/bin/sh

latex_file='NEMO_manual'
latex_opts='-shell-escape -interaction=nonstopmode'

cd tex_main

latex 		${latex_opts} 				${latex_file}
makeindex 	-s ${latex_file}.ist 	${latex_file}
bibtex 										${latex_file}
latex 		${latex_opts} 				${latex_file}

cd -

exit 0
