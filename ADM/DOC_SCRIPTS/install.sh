#!/bin/sh -x
#+
#
# ==========
# install.sh
# ==========
#
# ----------------------------------------------
# publication of HTML files and associated files
# ----------------------------------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ install.sh -w dirwww -p dirpublish -u urlpublish -l login
#
# DESCRIPTION
# ===========
#
# publication (rsync) of dirwww content on dirpublish given in argument
#
# If the host of publication is cerbere.locean-ipsl.upmc.fr, a specific update
# is launched.
#
# -w  input directort
# -p  output directory
# -u  output url
# -l  login used to access on output url
#
# If needed, existing directories might be removed before (to erase obsolete
# files), using ncftp tool :
#
# For example, to clean directory on LOCEAN web server ::
#
#  $ ncftp -u fplod www.locean-ipsl.upmc.fr
#  ncftp> cd fplod
#  ncftp> rm -f pageperso/.DS_Store
#  ncftp> rm -rf pageperso
#  ncftp> exit
#
# EXAMPLES
# ========
#
# EVOLUTIONS
# ==========
#
# $Id$
#
# - fplod 2008-09-16T15:24:26Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * comments in ReStructured Text
#
# - fplod 2008-06-17T09:10:19Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * add -l parameter only used in specific case at LOCEAN when user
#     parameter of persoweb must be different tthan login (ex: acmo vs fplod)
#   * replace http://www.lodyc.jussieu.fr/info_reseau/persoweb/?fastupdate=1&user=${user}" by
#     http://intranet.locean-ipsl.upmc.fr/persoweb/?fastupdate=1&user=${user}
#
# - fplod 2008-03-28T10:26:58Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * new personnal webpages policy at LOCEAN so new command and new parameter (-u)
#
# - fplod 2007-09-28T09:30:43Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#   * parametrisation and translation
#
# - smasson 2007-06-07T16:43:42Z arete.locean-ipsl.upmc.fr (Darwin)
#
#   * can give the answer with input parameters
#
# - fplod 2007-04-26T11:51:42Z aedon.locean-ipsl.upmc.fr (Darwin)
#
#-
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
usage=" Usage : ${command} -w dirwww -p dirpublish -u urlpublish -l login"
#
minargcount=4
#echo " narg ${#}"
if [ ${#} -lt ${minargcount} ]
then
   echo "eee : not enought arguments"
   echo "${usage}"
   exit 1
fi
unset minargcount
#
# default
dirpublish="none"
urlpublish="none"
login="none"
#
while [ ! -z "${1}" ]
do
   case ${1} in
      -w)
         dirwww=${2}
         shift
      ;;
      -p)
         dirpublish=${2}
         shift
      ;;
      -u)
         urlpublish=${2}
         shift
      ;;
      -l)
         login=${2}
         shift
      ;;
   esac
   # next flag
   shift
done
#
set -u
#
# ++ check directories
#
answer=${1:-" "}
case ${answer} in
   y|Y|n|N)
   ;;
   *)
      if [ "${dirpublish}" != "none" ]
      then
         echo "Do you want to install on ${dirpublish} (y|[n]) ?"
         read answer
      fi
      if [ "${urlpublish}" != "none" ]
      then
         echo "Do you want to install on ${urlpublish} (y|[n]) ?"
         read answer
      fi
   ;;
esac
#
case  ${answer} in
   y|Y)
      if [ "${dirpublish}" != "none" ]
      then
         # copy of ${dirwww} on $dirpublish
         echo "iii : update of ${dirpublish}"
         rsync -av --exclude=".DS_Store" -e ssh ${dirwww}/ ${dirpublish}
         # detect if in dirpublish following this pattern [USER@]HOST:SRC, HOST
         # is cerbere.locean-ipsl.upmc.fr. If so, a specific update is launched
         userhost=${dirpublish%%:*}
         host=${userhost##*@}
         if [ ${login} = "none" ]
         then
            user=${userhost%%@*}
         else
            user=${login}
         fi
         if [ "${host}" = "cerbere.locean-ipsl.upmc.fr" ]
         then
            wget -q "http://intranet.locean-ipsl.upmc.fr/persoweb/?fastupdate=1&user=${user}" -O /dev/null
         fi
      else
         # urlpublish=http://www.locean-ipsl.upmc.fr/~ginette/produit
         dirpublish=${urlpublish##*~}
         cd ${dirwww}
         #lftp -e "mirror -R . ${dirpublish};quit" -u ${LOGNAME} skyros.locean-ipsl.upmc.fr
         lftp -e "mirror -R . ${dirpublish};quit" -u ${LOGNAME} localhost
         # pour acmo a la main ++
         #++lftp -e 'mirror -R . acmo/nouveaux/;quit' -u fplod www.locean-ipsl.upmc.fr
         # ++ log
      fi
   ;;
   *)
      echo "no update of ${dirpublish} or ${urlpublish}"
   ;;
esac
#
# normal exit
exit 0
