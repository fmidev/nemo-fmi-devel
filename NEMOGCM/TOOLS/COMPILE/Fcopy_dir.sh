#!/bin/bash
#set -x
set -o posix
#set -u
#set -e
#+
#
# ============
# Fcopy_dir.sh
# ============
#
# --------------------------
# Copy a reference directory
# --------------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ Fcopy_dir.sh
#
#
# DESCRIPTION
# ===========
#
#
# When a reference configuration is set, 
# Copy NEMO sub-directories needed (OCE_SRC, TOP_SRC ...)
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./Fcopy_dir.sh ORCA2_LIM
#
#
# TODO
# ====
#
# option debug
#
#
# EVOLUTIONS
# ==========
#
# $Id$
#
#
#
#   * creation
#
#-

declare -a ZTAB
grep "$1 " ${CONFIG_DIR}/cfg.txt > ${CONFIG_DIR}/cfg.tmp
read -a ZTAB < ${CONFIG_DIR}/cfg.tmp
TAB=( ${ZTAB[@]:1} )
\rm ${CONFIG_DIR}/cfg.tmp

unset -v ZTAB
