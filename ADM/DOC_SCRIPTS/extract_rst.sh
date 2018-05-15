#! /bin/sh
set -x
#+
#
# ==========
# extract.sh
# ==========
#
# -----------------------------------------
# extract ReStructuredText from source file
# -----------------------------------------
#
# SYNOPSIS
# ========
#
# ``extract_rst.sh -i filein -l language -o fileout``
#
# DESCRIPTION
# ===========
#
# ``extract_rst.sh`` extracts ReST comments from the file given in argument
#
# -i  input file
# -o  output file (ReST)
# -l  language
#
# Comment block (start, end) identification depends on language :
#
#  *F90*
#   FORTRAN source free form
#  *fortran*
#   FORTRAN source fixed form
#  *sh*
#   shell scripts
#  *IDL*
#   IDL source
#  *xml*
#   XML and XSL
#  *dot*
#   graphviz files
#  *php*
#   PHP files
#  *matlab*
#   matlab or octave files
#
# EXAMPLES
# ========
#
# To extract ReST comments of this shell script::
#
#  $ extract_rst.sh -i extract_rst.sh -l sh -o extract_rst.sh.rst
#  iii : rst lines of extract_rst.sh are in extract_rst.sh.rst
#
# You can produce HTML file from this new file::
#
#  $ rst2html.py --input-encoding=ISO-8859-15 extract_rst.sh.rst \
#    /usr/temp/${LOGNAME}/extract_rst.sh.html
#
# You can produce PDF file from this new file::
#
#  $ rst2newlatex.py --input-encoding=ISO-8859-15 extract_rst.sh.rst \
#    /usr/temp/${LOGNAME}/extract_rst.sh.tex
#  $ pdflatex extract_rst.sh.tex
#
# Of course beware of consistency of path on links.
#
# CAUTIONS
# ========
#
# Becaue of poor implementation of Standard FORTRAN in cpp (prepocessing)
# within gfortran and g95, ReST comments might induce trouble in
# FORTRAN sources.
#
# For example following line is pointed out be gfortran with
# ``error: unterminated comment``.
# This is because ``/*`` is the beginning of a C style comment !!
# ::
#
#      !    **MEAN** = sum( *X*\ (:) )/*ntime*
#
#
#
# One can modify this ReST line with
# ::
#
#      !    **MEAN** = sum( *X*\ (:) ) / *ntime*
#
# TODO
# ====
#
# check parameters
#
# log
#
# add perl
#
# SEE ALSO
# ========
#
# ReStructuredText_
#
# .. _ReStructuredText: http://docutils.sourceforge.net/rst.html
#
# Docutils_
#
# .. _Docutils: http://docutils.sourceforge.net/
#
# EVOLUTIONS
# ==========
#
# $Id$
#
# - fplod 2009-04-20T08:13:37Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * add CAUTIONS paragraph to warn about possible FORTRAN compiling problem
#
# - fplod 2009-04-03T14:53:18Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * usage of tr instead of sed to remove ``\r``
#     due to difference between ``/sw/bin/sed`` and ``/usr/bin/sed`` (the last
#     one do not work coorectly on ``\r`` interpertation ie: remove the first occurence of
#     ``r``)
#
# - fplod 2009-02-10T10:46:23Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * add language fortran for FORTRAN source in fixed form
#
# - fplod 2009-01-05T11:41:33Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * remove \\r (CRLF) from file before awk and sed (otherwise ReST block
#     was not found in "ISO-8859 text, with CRLF line terminators" files
#
# - fplod 2008-12-22T10:37:37Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * add matlab (octave)
#
# - fplod 2008-09-17T13:40:37Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * add php language
#
# - fplod 2008-09-08T09:34:04Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * add F90 language
#
# - fplod 2008-08-08T08:28:30Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * add KML files (with other XML files)
#   * add parameters ``-i`` ``-l`` ``-o``
#
# - fplod 200807
#
#   * creation
#
#-
#
system=$(uname)
case "${system}" in
   AIX|IRIX64)
      echo " www : no specific posix checking"
   ;;
   *)
      set -o posix
   ;;
esac
unset system
#
command=$(basename ${0})
log_date=$(date -u +"%Y%m%dT%H%M%SZ")
log=/tmp/$(basename ${command} .sh).log.${log_date}
#
usage=" Usage : ${command} -i filein -l language -o fileout"
#
minargcount=6
#echo " narg ${#}"
if [ ${#} -lt ${minargcount} ]
then
   echo "eee : not enought arguments"
   echo "${usage}"
   exit 1
fi
#
# default
# n.a.
#
while [ ! -z "${1}" ]
do
   case ${1} in
      -i)
         filein=${2}
         shift
      ;;
      -o)
         fileout=${2}
         shift
      ;;
      -l)
         language=${2}
         shift
      ;;
      -h)
         echo "${usage}"
         exit 0
      ;;
      *)
         echo "eee : unknown option ${1}"
         echo "${usage}"
         exit 1
      ;;
   esac
   # next flag
   shift
done
#
set -u
#
# ++ check param
#
case "${language}" in
   fortran)
      awkblockstart="^C\+$"
      awkblockend="^C-$"
      sedblockstart="^C+$"
      sedblockend="^C-$"
      comment="^C"
   ;;
   F90)
      awkblockstart="^!\+$"
      awkblockend="^!-$"
      sedblockstart="^!+$"
      sedblockend="^!-$"
      comment="^!"
   ;;
   IDL)
      awkblockstart="^;\+$"
      awkblockend="^;-$"
      sedblockstart="^;+$"
      sedblockend="^;-$"
      comment="^;"
   ;;
   xml)
      awkblockstart="^<!--rst$"
      awkblockend="-->$"
      sedblockstart="^<!--rst$"
      sedblockend="-->$"
      comment=""
   ;;
   sh)
      # iii : awk '/^\#\+/,/^\#\-/' $file
      awkblockstart="^\#\+$"
      awkblockend="^\#\-$"
      sedblockstart="^#+"
      sedblockend="^#-"
      comment="^#"
   ;;
   dot|php)
      awkblockstart="^\/\*rst$"
      awkblockend="*\/"
      sedblockstart="^\/\*rst$"
      sedblockend="^\*\/"
      comment=""
   ;;
   matlab)
      awkblockstart="^%\+$"
      awkblockend="^%-$"
      sedblockstart="^%+$"
      sedblockend="^%-$"
      comment="^%"
   ;;
   *)
      echo "eee : ${language} not implemented"
      exit 1
   ;;
esac
#
# just in case suppress \r at the end of lines
tr -d '\r' < ${filein} > /tmp/${$}_0
#
# put rst blocks in one temporary file
#awk '/^;+/,/^;-/' a.pro | sed -e "/^;+$/d" -e "/^;-$/d" -e "s/^;//"
cmdawk="awk '/${awkblockstart}/,/${awkblockend}/' /tmp/${$}_0 > /tmp/${$}_1" #++
eval ${cmdawk}
if [ ! -s /tmp/${$}_1 ]
then
   rm /tmp/${$}_0 /tmp/${$}_1
   echo "iii : no rst comments in ${filein}"
   exit 1
fi
#
# suppress begin and end of each block
sedcmd="sed -e \"/${sedblockstart}/d\" -e \"/${sedblockend}/d\" /tmp/${$}_1 > /tmp/${$}_2"
eval ${sedcmd}
#
# suppress comment at the beginning of each line
if [ "${comment}" != "" ]
then
   sedcmd="sed -e \"s/${comment}//\" /tmp/${$}_2 > /tmp/${$}_3"
   eval ${sedcmd}
   # suppress first blank
   cp /tmp/${$}_3 /tmp/${$}_2
   sed -e "s/^ //" /tmp/${$}_2 > /tmp/${$}_3
   cp /tmp/${$}_3 ${fileout}
else
   cp /tmp/${$}_2 ${fileout}
fi
#
echo "iii : rst lines of ${filein} are in ${fileout}"
#
# clean
rm /tmp/${$}_0 /tmp/${$}_1 /tmp/${$}_2 /tmp/${$}_3 2> /dev/null
#
# exit
exit 0
