#!/bin/bash
cd ~/work/ALF_P2/nv_check/trees
FILES="Libraries/Modules/linalg_mod.F90 Prog/gpu_state_mod.F90 Prog/gpu_wrap_mod.F90 Prog/gpu_update_mod.F90 Prog/gpu_lowrank_mod.F90 Prog/gpu_stab_mod.F90 Prog/gpu_obs_mod.F90 Prog/gpu_tau_mod.F90"
pp() { (cd $1 && for f in $FILES; do gfortran -E -cpp -P $2 $f 2>/dev/null | grep -v "^[[:space:]]*$"; done); }
norm() { sed -E 's/(Terminate_on_error\([^,]*,"[^"]*",)[0-9]+\)/\1L)/'; }

for d in "" "-DSTABLOG" "-DMPI" "-DALF_GPU" "-DALF_GPU -DSTABLOG" "-DALF_GPU -DMPI"; do
  b=$(pp $BASE "$d" | md5sum); bn=$(pp $BASE "$d" | norm | md5sum)
  for c in $(ls | grep -v $BASE); do
    x=$(pp $c "$d" | md5sum); xn=$(pp $c "$d" | norm | md5sum)
    r=DIFF; [ "$xn" = "$bn" ] && r="same except __LINE__ of error calls"; [ "$x" = "$b" ] && r=same
    echo "[$d] ${c:0:8} $r"
  done
done
