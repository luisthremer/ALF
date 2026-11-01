#!/bin/bash
# Ticket 207 (development only), nividia: evaluation of t207_d2_t0b.sh RUN -> dev2/logs/t0b_<RUN>_eval.txt
W=$HOME/alf-gpu/work/t207/dev2; R=$1; out=$W/logs/t0b_${R}_eval.txt; : > $out
P=~/alf-gpu/venv/bin/python3; D=$W/runs/t0b$R
id() { echo "== ident $*" >> $out; $P -I $W/t207_d2_ident.py "$@" >> $out 2>&1; echo "   exit $?" >> $out; }
fd() { echo "== same seed $1 vs $2" >> $out; $P -I $W/t207_firstdiv.py $1/data.h5 $2/data.h5 | grep -E "bins:|first bin|tau-resolved|equal-time|FAIL" >> $out 2>&1; }
for b in nv_def nv_log; do
  id ${D}_${b}_6a/hub6x6_cb_fd ${D}_${b}_6b/hub6x6_cb_fd
  id ${D}_${b}_b10a/hub12x12_cb_b10_s ${D}_${b}_b10b/hub12x12_cb_b10_s
  id ${D}_${b}_6a/hub6x6_cb_fd ${D}_${b}_6lt0/hub6x6_cb_fd_lt0 --eq
  id ${D}_${b}_b10a/hub12x12_cb_b10_s ${D}_${b}_b10lt0/hub12x12_cb_b10_s_lt0 --eq
  id ${D}_base_${b}_6/hub6x6_cb_fd ${D}_${b}_6a/hub6x6_cb_fd --eq
done
fd ${D}_nv_def_6a/hub6x6_cb_fd ${D}_cpp6/hub6x6_cb_fd
fd ${D}_nv_log_6a/hub6x6_cb_fd ${D}_clog6/hub6x6_cb_fd
fd ${D}_nv_def_12/hub12x12_cb_fd ${D}_cpp12/hub12x12_cb_fd
fd ${D}_nv_log_12/hub12x12_cb_fd ${D}_clog12/hub12x12_cb_fd
for x in ${D}_*/*/info; do echo "$x: $(grep -E "Precision (Green|tau)" $x | tr -s ' ' | tr '\n' ' ')" ; done >> $out
