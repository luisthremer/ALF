# GNU 13.3 + OpenBLAS + Intel MPI (mpif90 wraps gfortran; ~/alf-gpu/bin/mpifort -> mpif90)
. /opt/intel/oneapi/mpi/latest/env/vars.sh
export I_MPI_F90=gfortran
export PATH=$HOME/alf-gpu/bin:$PATH
. $HOME/alf-gpu/venv/bin/activate
