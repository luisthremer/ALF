! Unit test for qmc_walker_mod: Walker_layout, Walker_alloc and Walker_dealloc
! on a small model set with Set_model_for_testing (Ndim = 4, N_FL = 2,
! Ltrot = 10, Nwrap = 3; finite temperature, no projector).
! This program should exit successfully (return code 0).
Program Test_walker
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main, only: Set_model_for_testing, Setup_flavor_map, Clear_Ham_state, N_FL_eff, nsigma
  Use QMC_runtime_var, only: set_N_walkers
  Use Operator_mod, only: Operator
  Use qmc_walker_mod, only: qmc_walker_t, Walker_layout, Walker_alloc, Walker_dealloc, &
       &                    Walker_activate, Walker_deactivate
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use iso_fortran_env, only: error_unit
  Implicit none

  Type (qmc_walker_t) :: walker, w1, w2
  Integer, allocatable :: seed1(:), seed2(:)
  Real (Kind=Kind(0.d0)) :: x, x_ref(2)
  Integer :: k
  Type (Operator), dimension(:,:), allocatable :: Op_T, Op_V
  Integer :: n_fail, nf_eff, n
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
#endif
  n_fail = 0

  allocate(Op_T(0,2), Op_V(0,2))
  Call Set_model_for_testing(4, 2, 1, 10, 0, .false., .false., Op_T, Op_V)
  Call Setup_flavor_map()

  ! Layout: 10 time slices in intervals of 3 -> NSTM = 4, Stab_nt = 0 3 6 9 10
  Call Walker_layout(walker, 10, 3)
  Call check(walker%NSTM == 4, "NSTM")
  Call check(all(walker%Stab_nt == [0, 3, 6, 9, 10]), "Stab_nt")
  Call check(lbound(walker%Stab_nt,1) == 0, "lbound Stab_nt")

  Call Walker_alloc(walker)
  Call check(all(shape(walker%GR) == [4, 4, 2]), "shape GR")
  Call check(all(shape(walker%GR_Tilde) == [4, 4, 2]), "shape GR_Tilde")
  Call check(all(shape(walker%Test) == [4, 4]), "shape Test")
  Call check(size(walker%Phase_array) == 2, "size Phase_array")
  Call check(N_FL_eff == 2, "N_FL_eff")
  Call check(size(walker%udvl) == N_FL_eff .and. size(walker%udvr) == N_FL_eff, "size udvl, udvr")
  Call check(all(shape(walker%udvst) == [4, N_FL_eff]), "shape udvst")
  do nf_eff = 1, N_FL_eff
     Call check(walker%udvl(nf_eff)%side == 'L' .or. walker%udvl(nf_eff)%side == 'l', "side udvl")
     Call check(walker%udvr(nf_eff)%side == 'R' .or. walker%udvr(nf_eff)%side == 'r', "side udvr")
     do n = 1, walker%NSTM
        Call check(allocated(walker%udvst(n,nf_eff)%U), "udvst allocated")
     enddo
  enddo

  Call Walker_dealloc(walker)
  Call check(.not. allocated(walker%GR) .and. .not. allocated(walker%udvst) .and. &
       &     .not. allocated(walker%Stab_nt) .and. .not. allocated(walker%Phase_array), "deallocated")

  ! Swap of fields and random number state between two walkers (N_walkers > 1)
  Call set_N_walkers(2)
  Call random_seed(size=k)
  allocate(seed1(k), seed2(k))
  seed1 = [(12345 + 7*n, n = 1, k)]
  seed2 = [(54321 + 3*n, n = 1, k)]
  Call random_seed(put=seed1)
  Call random_number(x_ref(1))
  Call random_number(x_ref(2))          ! reference: stream 1 without interruption

  Call nsigma%make(2, 3)
  nsigma%f = cmplx(1.d0, 0.d0, kind(0.d0)); nsigma%t = 1
  Call random_seed(put=seed1)
  Call random_number(x)                  ! walker 1 draws its first number
  Call Walker_deactivate(w1)
  Call check(.not. allocated(nsigma%f) .and. allocated(w1%f) .and. allocated(nsigma%t), "deactivate moves fields")

  allocate(nsigma%f(2, 3))
  nsigma%f = cmplx(-1.d0, 0.d0, kind(0.d0))
  Call random_seed(put=seed2)
  Call random_number(x)                  ! walker 2 uses another stream
  Call Walker_deactivate(w2)

  Call Walker_activate(w1)
  Call check(all(nsigma%f == cmplx(1.d0, 0.d0, kind(0.d0))), "fields of walker 1 back")
  Call random_number(x)
  Call check(x == x_ref(2), "random stream of walker 1 continues bitwise")
  Call Walker_deactivate(w1)
  Call Walker_activate(w2)
  Call check(all(nsigma%f == cmplx(-1.d0, 0.d0, kind(0.d0))), "fields of walker 2 back")
  Call nsigma%clear()
  Call set_N_walkers(1)

  Call Clear_Ham_state()

  If (n_fail > 0) then
     Write(error_unit,*) "walker test: ", n_fail, " failure(s)"
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  endif
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
  Write(*,*) "walker test passed"

contains

  Subroutine check(ok, label)
    Logical, intent(in) :: ok
    Character(len=*), intent(in) :: label
    If (.not. ok) then
       Write(error_unit,*) "FAIL ", label
       n_fail = n_fail + 1
    endif
  end Subroutine check

End Program Test_walker
