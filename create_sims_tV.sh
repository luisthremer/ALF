#!/bin/bash

Theta=80
Ltau=0
Nwrap=10
NSweep=20
NBin=1
n_skip=5
N_rebin=1

fortran_double() {
    local val=$1
    if [[ "$val" == *.* ]]; then
        printf "%sd0" "$val"
    else
        printf "%s.d0" "$val"
    fi
}

for h in 1; do
    for L in 130; do
        ham_h0=${h}
        N_Phi=$(( L / 2 ))
        Beta=${L}
        ham_h0_fortran=$(fortran_double "$ham_h0")

        V_values=(0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9)
        N_params=${#V_values[@]}
        mpi_per_parameter_set=48

        dir="$ALF/Benchmark/tV/Den/tV_L${L}_Theta${Theta}_HamH${ham_h0}_PP"
        mkdir -p "$dir"

        cp $ALF/Scripts_and_Parameters_files/Start/seeds "$dir"
        cp $ALF/Scripts_and_Parameters_files/Start/out_to_in.sh "$dir"

        cat <<PAR > "$dir/parameters"
!=======================================================================================
!  Input variables for a general ALF run
!---------------------------------------------------------------------------------------
&VAR_ham_name
ham_name = "tV"
/

&VAR_lattice
L1           = 130
L2           = 1
Lattice_type = "Square"
Model        = "tV"
/

&VAR_Model_Generic
Checkerboard = .T.
Symm         = .T.
N_SUN        = 1
N_FL         = 1
Phi_X        = 0.d0
Phi_Y        = 0.d0
Bulk         = .T.
N_Phi        = 0
Dtau         = 0.1d0
Beta         = 4.d0
Projector    = .T.
Theta        = 0.1d0
/

&VAR_QMC
Nwrap                = 10
NSweep               = 20
NBin                 = 5
Ltau                 = 1
LOBS_ST              = 0
LOBS_EN              = 0
CPU_MAX              = 0.0
Propose_S0           = .F.
Global_moves         = .F.
N_Global             = 1
Global_tau_moves     = .F.
N_Global_tau         = 1
Sequential           = .T.
Nt_sequential_start  = 0
Nt_sequential_end    = -1
Langevin             = .F.
Delta_t_Langevin_HMC = 0.1
Max_Force            = 5.0
HMC                  = .F.
N_HMC_sweeps         = 1
Leapfrog_steps       = 0
Amplitude            = 1.d0
/

&VAR_errors
n_skip  = 1
N_rebin = 1
N_Cov   = 0
Extended_Zone = .F.
N_BZ_Zones = 1
/

&VAR_TEMP
N_exchange_steps      = 6
N_Tempering_frequency = 10
mpi_per_parameter_set = ${mpi_per_parameter_set}
Tempering_calc_det    = .T.
/

&VAR_tV
ham_T      = 1.d0
ham_chem   = 0.d0
ham_V      = 0.d0
ham_h0     = 1.d0
/
PAR

        for idx in "${!V_values[@]}"; do
            V=${V_values[$idx]}
            ham_V_fortran=$(fortran_double "$V")

            tdir="$dir/Temp_${idx}"
            mkdir -p "$tdir"

            cat <<PAR > "$tdir/parameters"
&VAR_ham_name
ham_name = "tV"
/

&VAR_lattice
L1           = ${L}
L2           = 1
Lattice_type = "Square"
Model        = "tV"
/

&VAR_Model_Generic
Checkerboard = .T.
Symm         = .T.
N_SUN        = 1
N_FL         = 1
Phi_X        = 0.d0
Phi_Y        = 0.d0
Bulk         = .T.
N_Phi        = 0
Dtau         = 0.1d0
Beta         = ${Beta}.d0
Projector    = .T.
Theta        = ${Theta}.d0
/

&VAR_QMC
Nwrap                = ${Nwrap}
NSweep               = ${NSweep}
NBin                 = ${NBin}
Ltau                 = ${Ltau}
LOBS_ST              = 0
LOBS_EN              = 0
CPU_MAX              = 11.9
Propose_S0           = .F.
Global_moves         = .F.
N_Global             = 1
Global_tau_moves     = .F.
N_Global_tau         = 1
Sequential           = .T.
Nt_sequential_start  = 0
Nt_sequential_end    = -1
Langevin             = .F.
Delta_t_Langevin_HMC = 0.1
Max_Force            = 5.0
HMC                  = .F.
N_HMC_sweeps         = 1
Leapfrog_steps       = 0
Amplitude            = 1.d0
/

&VAR_errors
n_skip        = ${n_skip}
N_rebin       = ${N_rebin}
N_Cov         = 0
Extended_Zone = .F.
N_BZ_Zones    = 1
/

&VAR_tV
ham_T      = 1.d0
ham_chem   = 0.d0
ham_V      = ${ham_V_fortran}
ham_h0     = ${ham_h0_fortran}
ham_T2     = 1.d0
ham_V2     = 4.d0
ham_Tperp  = 1.d0
ham_Vperp  = 1.d0
/
PAR

            cat <<PAR > "$tdir/params.json"
{
    "Pin": "Point",
    "L": ${L},
    "Theta": ${Theta},
    "HamH": ${ham_h0},
    "U": ${V}
}
PAR

        done

        N_nodes=${N_params}
        cat <<PAR > "$dir/job.sh"
#!/bin/bash
#SBATCH --job-name PP_tV_L${L}_H${ham_h0}
#SBATCH --output=out.%j.log
#SBATCH --error=err.%j.log
#SBATCH --mail-type=NONE
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn73xu
#SBATCH --partition=micro
#SBATCH --nodes=${N_nodes}
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1

module load slurm_setup
module load hdf5/1.10.7-intel21

export KMP_AFFINITY=granularity=fine,compact
export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=auto:cache
export I_MPI_PIN_ORDER=scatter
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

bash ./out_to_in.sh
mpiexec -n \$SLURM_NTASKS \$ALF/Prog/ALF.out
PAR

    done
done