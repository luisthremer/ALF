! Unit test for the start-up checks in QMC_runtime_var. Usage: 42-QMC-runtime-checks-fail <case>
! Every case sets a valid parameter set except one value and must terminate with the message that
! ctest matches (PASS_REGULAR_EXPRESSION). KNOWN_ISSUES.md: H3, M2, M4, M5, Z3.
Program Test_QMC_runtime_checks_fail

  Use QMC_runtime_var, only: set_QMC_runtime_default_var, set_NSweep, set_LOBS_ST, set_LOBS_EN, &
       &                     set_Ltau, set_Nt_sequential_start, set_Nt_sequential_end, &
       &                     set_default_values_measuring_interval, check_projector_Ltau, &
       &                     check_Nt_sequential_range, check_Langevin_HMC_variables, set_Amplitude, &
       &                     set_Langevin, set_HMC, set_Delta_t_Langevin_HMC, set_Max_Force, &
       &                     set_Leapfrog_steps, set_N_HMC_sweeps, set_NBin
#ifdef MPI
  Use mpi
#endif
  Implicit none
  Character (len=32) :: name
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
#endif

  Call set_QMC_runtime_default_var()
  Call set_NSweep(5)
  Call set_NBin(2)
  Call get_command_argument(1, name)

  ! Valid values must pass
  Call set_default_values_measuring_interval(0, 10, .false.)
  Call check_projector_Ltau(.true., 2)
  Call set_Nt_sequential_start(1);  Call set_Nt_sequential_end(4)
  Call check_Nt_sequential_range(4)
  Call set_Delta_t_Langevin_HMC(0.1d0);  Call set_Max_Force(5.d0);  Call set_Leapfrog_steps(2)
  Call set_Langevin(.true.);  Call check_Langevin_HMC_variables()
  Call set_Langevin(.false.);  Call set_HMC(.true.);  Call check_Langevin_HMC_variables()
  write(*,*) "valid values accepted"

  Select Case (trim(name))
  Case ("nsweep")        ! no measurements
     Call set_NSweep(0)
     Call set_default_values_measuring_interval(0, 10, .false.)
  Case ("lobs")          ! empty measuring interval
     Call set_LOBS_ST(5);  Call set_LOBS_EN(2)
     Call set_default_values_measuring_interval(0, 10, .false.)
  Case ("lobs_window")   ! measuring interval outside 1..Ltrot (P2.8 fuzzing)
     Call set_LOBS_ST(-5);  Call set_LOBS_EN(-1)
     Call set_default_values_measuring_interval(0, 10, .false.)
  Case ("nbin")          ! no bins (P2.8 fuzzing)
     Call set_NBin(0)
     Call set_default_values_measuring_interval(0, 10, .false.)
  Case ("projector")     ! Thtrot = 0 with time-displaced measurements
     Call set_Ltau(1)
     Call check_projector_Ltau(.true., 0)
  Case ("ntseq_start")
     Call set_Nt_sequential_start(0)
     Call check_Nt_sequential_range(4)
  Case ("ntseq_end")
     Call set_Nt_sequential_end(5)
     Call check_Nt_sequential_range(4)
  Case ("ltau")          ! Ltau other than 0 or 1
     Call set_Ltau(2)
     Call check_Langevin_HMC_variables()
  Case ("amplitude")
     Call set_Amplitude(0.d0)
     Call check_Langevin_HMC_variables()
  Case ("langevin_dt")
     Call set_HMC(.false.);  Call set_Langevin(.true.);  Call set_Delta_t_Langevin_HMC(0.d0)
     Call check_Langevin_HMC_variables()
  Case ("langevin_force")
     Call set_HMC(.false.);  Call set_Langevin(.true.);  Call set_Max_Force(0.d0)
     Call check_Langevin_HMC_variables()
  Case ("hmc_dt")
     Call set_Delta_t_Langevin_HMC(0.d0)
     Call check_Langevin_HMC_variables()
  Case ("hmc_leapfrog")
     Call set_Leapfrog_steps(0)
     Call check_Langevin_HMC_variables()
  Case ("hmc_sweeps")
     Call set_N_HMC_sweeps(0)
     Call check_Langevin_HMC_variables()
  Case Default
     write(*,*) "ERROR: unknown case ", trim(name)
  end Select

  write(*,*) "ERROR: invalid value did not terminate"
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
end Program Test_QMC_runtime_checks_fail
