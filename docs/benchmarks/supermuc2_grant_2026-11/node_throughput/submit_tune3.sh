#!/bin/bash
# Tuning stage 3: N_delay beyond 256 and combined with ASYNCBLK=0 (Intel-opt.), beta = 4, best layouts of stage 1; L = 32 N_delay scan. Results in runs_tune3/.
set -e
A=/hppfs/work/pn73xu/di38zuw/ANTRAG_2026-11-27
export BETA=4 RUNDIR=runs_tune3 EXE_GPU=$A/build_pp_alfgpu_detfix/Prog/ALF.out BINB=$A/bin_batched_65800b82/pp
P=PINMODE=excl
sbatch --export=BETA,RUNDIR,EXE_GPU,BINB -N 6 -J nt_t3_gpub_40 -p general -t 00:50:00 multi2.sh gpub:40:7_2:6:nd256asy0:$P,NDELAY=256,ALFSH_ASYNCBLK=0 gpub:40:7_2:6:nd512:$P,NDELAY=512 gpub:40:7_2:6:nd512asy0:$P,NDELAY=512,ALFSH_ASYNCBLK=0 gpub:40:7_2:6:nd1024:$P,NDELAY=1024 gpub:40:8_1:6:nd256:$P,NDELAY=256 gpub:40:8_1:6:nd512:$P,NDELAY=512
sbatch --export=BETA,RUNDIR,EXE_GPU,BINB -N 5 -J nt_t3_gpub_56 -p general -t 01:30:00 multi2.sh gpub:56:8_1:4:nd512:$P,NDELAY=512 gpub:56:8_1:4:nd1024:$P,NDELAY=1024 gpub:56:8_1:4:nd512asy0:$P,NDELAY=512,ALFSH_ASYNCBLK=0 gpub:56:7_2:4:nd256:$P,NDELAY=256 gpub:56:7_2:4:nd512:$P,NDELAY=512
sbatch --export=BETA,RUNDIR,EXE_GPU,BINB -N 4 -J nt_t3_gpu_40 -p general -t 00:50:00 multi2.sh gpu:40:7_2:6:nd512:$P,NDELAY=512 gpu:40:7_2:6:nd1024:$P,NDELAY=1024 gpu:40:8_1:6:nd256:$P,NDELAY=256 gpu:40:8_1:6:nd512:$P,NDELAY=512
sbatch --export=BETA,RUNDIR,EXE_GPU,BINB -N 4 -J nt_t3_gpu_56 -p general -t 01:30:00 multi2.sh gpu:56:7_2:4:nd512:$P,NDELAY=512 gpu:56:7_2:4:nd1024:$P,NDELAY=1024 gpu:56:8_1:4:nd256:$P,NDELAY=256 gpu:56:8_1:4:nd512:$P,NDELAY=512
sbatch --export=BETA,RUNDIR,EXE_GPU,BINB -N 8 -J nt_t3_L32 -p general -t 00:50:00 multi2.sh gpub:32:4_2:6:nd64:$P,NDELAY=64 gpub:32:4_2:6:nd128:$P,NDELAY=128 gpub:32:4_2:6:nd256:$P,NDELAY=256 gpub:32:7_2:6:nd128:$P,NDELAY=128 gpu:32:4_3:6:nd64:$P,NDELAY=64 gpu:32:4_3:6:nd128:$P,NDELAY=128 gpu:32:4_3:6:nd256:$P,NDELAY=256 gpu:32:7_2:6:nd128:$P,NDELAY=128
