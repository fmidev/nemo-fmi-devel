#!/bin/bash

#Config ccc_msub
#MSUB -A  gen7451
#MSUB -r  NEMO_CI
#MSUB -o  NEMO_CI_%I
#MSUB -e  NEMO_CI_%I
#MSUB -oe
#MSUB -q  xlarge
#MSUB -n  128
#MSUB -N  1
#MSUB -T  1800
##MSUB -@  ntmlod@locean-ipsl.upmc.fr:begin,end
ccc_mprun ./opa
