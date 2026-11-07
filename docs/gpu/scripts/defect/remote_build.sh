#!/bin/bash
# usage: remote_build.sh gpu|cpu   (runs on nividia)
set -o pipefail
W=~/alf-gpu/work/fg1; export TMPDIR=$W/tmp
cd $W/src_$1 || exit 1
if [ $1 = gpu ]; then . ~/alf-gpu/env/nvhpc.sh; make clean >/dev/null 2>&1; . ./configure.sh PGI noMPI GPU=nvidia NO-INTERACTIVE > $W/logs/conf_$1.log 2>&1
elif [ $1 = cpunv ]; then . ~/alf-gpu/env/nvhpc.sh; make clean >/dev/null 2>&1; . ./configure.sh PGI noMPI NO-INTERACTIVE > $W/logs/conf_$1.log 2>&1
else . ~/alf-gpu/env/gnu.sh; make clean >/dev/null 2>&1; . ./configure.sh GNU noMPI NO-INTERACTIVE > $W/logs/conf_$1.log 2>&1; fi
make lib program > $W/logs/make_$1.log 2>&1 && cp Prog/ALF.out $W/bin/ALF-$1.out && echo BUILD_OK_$1 || { echo BUILD_FAIL_$1; tail -30 $W/logs/make_$1.log; }
