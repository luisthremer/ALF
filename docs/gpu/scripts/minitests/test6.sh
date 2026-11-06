cd ~/alf-gpu/work/dlsub
export TMPDIR=~/alf-gpu/work/dlsub/tmp; . ~/alf-gpu/env/nvhpc.sh
rm -rf srcB; cp -a src srcB
cd srcB/Prog
python3 - <<'PY'
t=open('gpu_lowrank_mod.F90').read()
assert '!pgi$r opt=1\n  Subroutine blk_end' in t
t=t.replace('!pgi$r opt=1\n  Subroutine blk_end','!(pin removed)\n  Subroutine blk_end')
open('gpu_lowrank_mod.F90','w').write(t)
PY
cd ..; . ./configure.sh PGI noMPI GPU=nvidia NO-INTERACTIVE >/dev/null 2>&1
BASE="$ALF_FLAGS_PROG"; cd Prog
runv() { # label, flags
  export ALF_FLAGS_PROG="$1"; rm -f gpu_lowrank_mod.o
  if ! make ALF.out > ../../make_t6.log 2>&1; then echo "$2: BUILD FAILED: $(grep -i 'NVFORTRAN-' ../../make_t6.log | head -1 | cut -c1-90)"; return; fi
  r=$(cd ~/alf-gpu/work/dlsub; VAR=srcB SCHEME=submatrix TMO=60 bash rc.sh nv t6_$3 t6.sed | head -1)
  res=$(grep -c "Unstable\|exceeding" ~/alf-gpu/work/dlsub/r_t6_$3/out.log)
  echo "$2: $( [ $res -gt 0 ] && echo FAIL-unstable || echo ok ) | $r" | cut -c1-120
}
runv "$BASE" "default (-O3 -fast), no pin" v0
runv "$(echo $BASE | sed 's/-O3/-O2/')" "-O2 -fast" v1
runv "$(echo $BASE | sed 's/-O3/-O1/; s/-fast//')" "-O1" v2
runv "$BASE -Mnovect" "-O3 -fast -Mnovect" v3
runv "$BASE -Mnounroll" "-O3 -fast -Mnounroll" v4
runv "$BASE -Mnoautoinline" "-O3 -fast -Mnoautoinline" v5
runv "$BASE -Mnolre" "-O3 -fast -Mnolre" v6
