cd ~/alf-gpu/work/dlsub
M=$1
mk() { cp t3.sed sed_$1.sed; shift; for e in "$@"; do echo "$e" >> sed_$1.sed 2>/dev/null; done; }
cp t3.sed sed_Hubbard.sed
cp t3.sed sed_Hubbard_Plain_Vanilla.sed
cp t3.sed sed_Kondo.sed; printf '%s\n' 's/^Lattice_type .*/Lattice_type = "Bilayer_square"/' 's/^Model .*/Model = "Kondo"/' >> sed_Kondo.sed
cp t3.sed sed_LRC.sed; printf '%s\n' 's/^Model .*/Model = "LRC"/' >> sed_LRC.sed
cp t3.sed sed_tV.sed; printf '%s\n' 's/^Model .*/Model = "tV"/' 's/^Nwrap .*/Nwrap = 4/' >> sed_tV.sed
cp t3.sed sed_Z2_Matter.sed; printf '%s\n' 's/^Model .*/Model = "Z2_Matter"/' >> sed_Z2_Matter.sed
cp t3.sed sed_Spin_Peierls.sed; printf '%s\n' 's/^Model .*/Model = "Spin_Peierls"/' >> sed_Spin_Peierls.sed
cp t3.sed sed_LR_Jz_Test.sed; printf '%s\n' 's/^Model .*/Model = "LR_Jz_Test"/' >> sed_LR_Jz_Test.sed
for h in Hubbard Hubbard_Plain_Vanilla Kondo LRC tV Z2_Matter Spin_Peierls LR_Jz_Test; do
  TMO=120 bash rc.sh $M t3_$h sed_$h.sed $h | head -1 | cut -c1-175
  if [ $h = LR_Jz_Test ]; then
     rm -f r_t3_$h/RUNNING
     printf '\n&VAR_LR_Jz_Test\nham_T = 1.d0\nham_U = 4.d0\nham_Jz = 1.d0\nHam_chem = 0.d0\nDtau = 0.1d0\nBeta = 2.d0\n/\n' >> r_t3_$h/parameters
     (cd r_t3_$h; . ~/alf-gpu/env/$( [ $M = nv ] && echo nvhpc || echo intel ).sh 2>/dev/null; export OMP_TARGET_OFFLOAD=MANDATORY ONEAPI_DEVICE_SELECTOR=level_zero:0; timeout 120 ~/alf-gpu/work/dlsub/src/Prog/ALF.out > out.log 2>&1; echo "   LR_Jz_Test rerun rc=$? $(grep -h 'GPU path (gpu_state)' info | cut -c1-110)")
  fi
  rc=$(grep -c . r_t3_$h/out.log)
  grep -h -i "error\|unstable\|Wrong\|exceeding" r_t3_$h/out.log | head -2 | cut -c1-150 | sed "s/^/   ! /"
done
