#!/bin/bash
# Ticket 207 (rebase on 212, development only), nividia: T0 of the GPU=none checks with the harness's own acc_gpunone.sh, pointed
# at a development "candidate" in dev2/fakecand (tarballs of dev2/src = rsynced working tree and dev2/base_src = git archive
# 3f8ec8fa; changed_files.txt = git diff --name-status 3f8ec8fa of the working tree, shipped by the laptop as dev2/changed_files.txt).
#   usage: t207_d2_gpunone.sh OUTNAME [--plant]   -> dev2/OUTNAME/{ppcmp.txt,ctest.log,...}
W=$HOME/alf-gpu/work/t207/dev2; export TMPDIR=$W/tmp; : "${HOME:?}"
O=$1; P=${2:-}; case "$O" in ""|*/*|*..*) echo bad name; exit 1;; esac
F=$W/fakecand
rm -rf -- "$W/fakecand"; mkdir -p $F
tar -C $W/src -cf $F/src.tar . && tar -C $W/base_src -cf $F/base_src.tar . || exit 2
mkdir $F/src $F/base_src
cp $W/changed_files.txt $F/changed_files.txt || exit 2
{ echo "sha dev2-working-tree"; echo "base_sha 3f8ec8fa"; echo "tar_sha256 $(sha256sum $F/src.tar | cut -d' ' -f1)"
  echo "base_tar_sha256 $(sha256sum $F/base_src.tar | cut -d' ' -f1)"; } > $F/.alf_candidate
sed -e 's#^C=\$HOME/alf-gpu/work/t\$TICKET/cand/\$WHO/\$S12#C='"$F"'#' -e 's#export TMPDIR=\$HOME/alf-gpu/work/t\$TICKET/tmp#export TMPDIR='"$W"'/tmp#' \
    -e 's/-j4/-j3/g' $W/src/testsuite/gpu_reference/acceptance/acc_gpunone.sh > $W/acc_gpunone_dev.sh   # -j3: a GPU run may be alongside (<= 4 cores)
grep -q "^C=$F" $W/acc_gpunone_dev.sh || { echo "path substitution failed"; exit 2; }
bash $W/acc_gpunone_dev.sh 207 dev x "$W/$O" $P
