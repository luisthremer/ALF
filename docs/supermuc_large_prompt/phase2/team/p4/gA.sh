# Gate 3a/3c: same-seed GPU (device QR on, ALFSH_DEVQR=1) vs the same-compiler CPU build, + planted defect (tau(1) *= 1+1e-12 in a throwaway libalfdevqr.so)
source $SCRATCH/large_prompt/team/p4/jobs/env_intel_pvc.sh
L0=$SCRATCH/large_prompt/team; P=$L0/p4; cd $P/jobs
export W=$P/work_gA CPUBIN=$L0/build_cpu-txt/Prog/ALF.out GPUBIN=$P/build_txtnew/Prog/ALF.out SEEDS=$L0/src/Scripts_and_Parameters_files/Start/seeds
export OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0 ALFSH_TIMING=1
mkdir -p $W; echo "host $(hostname) $(date -Is) gpu=$GPUBIN"
run() { python3 run_cmp_sng.py "$@" 2>&1 | grep -E "wall=|PASS|FAIL|confout|acceptance|_scal "; }
echo "######## devqr ON (default lib next to the build)"
export ALFSH_DEVQR=1 ALFSH_QRNB=${NB:-32}
run L6_mzF_sub  L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
CMP_TOL=1e-10 run L6_mzT_sub  L=6 mz=T scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
run L12_mzF_sub L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72
echo "--- shim call counts of the last GPU run (geqrf_dev must be > 0):"; grep -h "ALFSH_TIMING geqrf" $W/runs/L12_mzF_sub/gpu/run.log
echo "######## 3c planted defect: tau(1) *= 1+1e-12 in the device QR library"
mkdir -p $P/plant && cd $P/plant && cp $P/src/Libraries/Modules/alf_la_sycl_devqr.* . 
sed -i 's|tau\[2 \* i\] = t.r; tau\[2 \* i + 1\] = t.i;|tau[2 * i] = t.r * (i == 0 ? 1.0 + 1e-12 : 1.0); tau[2 * i + 1] = t.i * (i == 0 ? 1.0 + 1e-12 : 1.0);|' alf_la_sycl_devqr.hpp
grep -c "1.0 + 1e-12" alf_la_sycl_devqr.hpp
icpx -fsycl -O2 -std=c++17 -fPIC -shared -I$MKLROOT/include alf_la_sycl_devqr.cpp -o libalfdevqr_plant.so -qmkl -lmkl_sycl_blas || exit 1
cd $P/jobs; export ALFSH_DEVQR_LIB=$P/plant/libalfdevqr_plant.so
run D_L6_mzF_sub  L=6 mz=F scheme=submatrix ndelay=8 nsweep=10 nbin=2 U=5.3d0 symm=T beta=6.d0 nphi=18
run D_L12_mzF_sub L=12 mz=F scheme=submatrix ndelay=32 nsweep=4 nbin=1 U=5.3d0 symm=T beta=12.d0 nphi=72
echo "done $(date -Is)"
