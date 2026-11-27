#!/bin/bash
# plant.sh: build one copy per planted defect, run the K1 self-test and the e2e cases; every defect must fail something
cd ~/work/ALF_P2/p26
plant() { # name file python-replacement
  d=pl-$1; rsync -a --delete --exclude .git ~/work/ALF_P2/wt-p26/ $d/
  python3 - "$d/Prog/$2" "$3" "$4" <<"PY"
import sys; p,a,b=sys.argv[1:4]; s=open(p).read(); assert s.count(a)>=1,(p,a); open(p,"w").write(s.replace(a,b,1))
PY
  (cd $d && . ./configure.sh GNU noMPI GPU=host NO-INTERACTIVE >/dev/null 2>&1 && make -j16 lib program >../$d.make.log 2>&1) || { echo "$1 BUILD FAIL"; return; }
  rm -rf $d-run; python3 new-host/docs/ai/phase2/p25_runs/p25_run.py --exe $PWD/$d/Prog/ALF.out --start $PWD/new-host/Scripts_and_Parameters_files/Start --work $PWD/$d-run --only tV-std,tV-dly4,kondo-std,kondo-dly4,sp-std,sp-dly4 new-host/docs/ai/phase2/p26_runs/cases_p26.py | cut -f2,4,6 | tr "\n" " "; echo
  for m in tV-std sp-std; do s=$d-run/st-$m; mkdir -p $s; cp run2/$m/active/parameters run2/$m/active/seeds $s/; (cd $s && ALF_GPU_SELFTEST=2 OMP_NUM_THREADS=1 timeout 900 ../../$d/Prog/ALF.out > stdout 2>&1; echo "  $1 $m selftest2 rc=$? $(grep -h "FAILED" stdout | head -1)"); done
}
plant ntype2 gpu_wrap_mod.F90 "       call k_nd_rows(G, ld, Ndim, nf, t, n, n, 6, 1, 0)" "       continue"
plant ewrongidx gpu_wrap_mod.F90 "          w = zk(a) * w" "          w = zk(b) * w"
plant wtransp gpu_update_mod.F90 "    W = matmul(X, transpose(Tri))" "    W = transpose(matmul(X, transpose(Tri)))"
plant nosqrt gpu_wrap_mod.F90 " * sqrt(cmplx(1.d0 + aimag(f), 0.d0, kind(0.d0)))" ""
plant genwrapY gpu_lowrank_mod.F90 "    call k_nd_rows(Y_d(1,1,nf), Ndim, jc(nf), nf, t, n, n, mc, sc, 1)" "    continue"
