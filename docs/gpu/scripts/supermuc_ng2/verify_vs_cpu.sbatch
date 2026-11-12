#!/bin/bash
#SBATCH -J alf_cmp2
#SBATCH --account=pn93la
#SBATCH --time=00:40:00
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -o /hppfs/scratch/0C/di38zuw/alf/cmp2_%j.out
module load slurm_setup
source ~/alf-gpu/env251.sh
export OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1
export W=/hppfs/scratch/0C/di38zuw/alf/cmp2 CPUBIN=$HOME/alf-gpu/bin/ALF-cpuref2.out GPUBIN=$HOME/alf-gpu/bin/ALF-gpuaot.out
rm -rf $W; mkdir -p $W/tmp; ln -sfn $HOME/alf-gpu/srcm $W/src_gpu
cd ~/alf-gpu
for spec in "c1 L=6 mz=F scheme=standard" "c2 L=6 mz=T scheme=standard" "c3 L=6 mz=F scheme=submatrix ndelay=8" "c4 L=6 mz=T scheme=submatrix ndelay=8" "c5 L=8 mz=F scheme=submatrix ndelay=32 U=5.3d0 nsweep=5 nbin=3" "c6 L=12 mz=F scheme=submatrix ndelay=128 U=5.3d0 nsweep=3 nbin=2" "c7 L=12 mz=T scheme=submatrix ndelay=128 U=5.3d0 nsweep=3 nbin=2" "c8 L=12 mz=F scheme=submatrix ndelay=32 U=5.3d0 nsweep=3 nbin=2"; do
  echo "######## $spec"; python3 run_cmp.py $spec 2>&1 | grep -v joules
done
echo "######## production binary (MPI+HDF5, 1 rank) on case c1 parameters"
P=$W/prod1; mkdir -p $P; cp $W/runs/c1/gpu/parameters $W/runs/c1/gpu/seeds $P/; cd $P
mpiexec -n 1 $HOME/alf-gpu/bin/ALF-prodaot.out > run.log 2>&1; echo rc=$?; ls
