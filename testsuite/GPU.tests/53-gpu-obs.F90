! Test of the device observables (P2.4 step 1e), registered only when the build has -DALF_GPU.
! Usage: 53-gpu-obs eq <case> | control <1|2>; one ctest per case (see CMakeLists.txt), every run in its own directory.
!
! Random G (and random GT0, G0T, G00, GTT, which are not Green functions: the observables are algebraic in them) go through
! the host ham%Obser / ham%OBSERT (on Symm(G) when the case has Symm, as the sweep does) and through gpu_obs_measure /
! gpu_obs_tau_measure + gpu_obs_flush.  Every observable (counters, sign, scalars, Obs_Latt, Obs_Latt0, Obs_local) must agree to
! 1e-12 relative to the largest entry of that observable.  The second equal-time measurement reads the device G only (the host
! array handed over is zero), the first uploads the host G.
!
!   case   lattice            N_FL  Symm  ham_h0
!   sq1    Square 4x4          1    F     0
!   sq1s   Square 4x4          1    T     0
!   sq2    Square 4x4          2    F     0
!   sq2s   Square 4x4          2    T     0
!   sq2h   Square 4x4          2    T     0.5   (local observable)
!   bil1   Bilayer_honeycomb   1    F     0     (U = 4, U2 = 2)
!   bil2   Bilayer_honeycomb   2    T     0     (U = 4, U2 = 2)
!
!   control 1|2   negative controls on bil2: 1 = the device G + 1e-6, 2 = the weights U on every orbital (U2 ignored); the
!             disagreement must exceed the tolerance by at least 1e3 (ALF_GPU_FAULT=obs is run as an eq case that must fail).
! This program should exit successfully (return code 0) in the modes that are expected to pass.
Program Test_gpu_obs
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main, only: ham, Clear_Ham_state, Ltrot, Thtrot, N_FL, N_SUN, Ndim, Obs_scal, Obs_eq, Obs_tau, Obs_local
#ifdef MPI
  Use Hamiltonian_main, only: Set_Group_Comm
#endif
  Use QMC_runtime_var, only: set_QMC_runtime_default_var, set_ham_name, set_Nwrap, set_Ltau, set_NSweep
  Use qmc_setup_mod, only: Setup_model, Setup_fields, Setup_update_schemes, Setup_gpu, Finalize_gpu
  Use qmc_walker_mod, only: qmc_walker_t, Walker_layout, Walker_alloc, Walker_dealloc
  Use gpu_state_mod, only: gpu_active, G_d, g_dev_ok
  Use gpu_stab_mod, only: gpu_stab_on
  Use gpu_obs_mod, only: gpu_obs_setup, gpu_obs_report, gpu_obs_active, gpu_obs_tau_active, gpu_obs_measure, &
       &                 gpu_obs_tau_measure, gpu_obs_flush, gpu_obs_finalize
  Use Hop_mod, only: Hop_mod_Symm
  Use Lattices_v3, only: Lattice, Unit_cell
  Use Predefined_Lattices, only: Predefined_Latt
  Use Predefined_Hoppings, only: Hopping_Matrix_type, Set_Default_hopping_parameters_square, &
       &                         Set_Default_hopping_parameters_Bilayer_honeycomb
  Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC
  Use iso_fortran_env, only: error_unit, int64
  Implicit none

  Type (qmc_walker_t) :: walker
  Character (len=64) :: mode, cname
  Integer :: n_fail, ncontrol
  Real (Kind=Kind(0.d0)), parameter :: tol = 1.d-12
  Complex (Kind=Kind(0.d0)), Allocatable :: buf(:)
  Integer :: nbuf
  Integer :: blk(64), nblk
  Character (len=24) :: blabel(64)
  Integer(int64) :: seed
  Logical :: symm_run
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  Call Set_Group_Comm(MPI_COMM_WORLD)
#endif
  n_fail = 0
  seed = 987654321_int64
  Call get_command_argument(1, mode)
  Call get_command_argument(2, cname)
  Select case (trim(mode))
  Case ('eq')
     Call run_case(trim(cname), 0)
  Case ('control')
     Read(cname,*) ncontrol
     Call run_case('bil2', ncontrol)
  Case default
     Write(error_unit,*) 'usage: 53-gpu-obs eq <case> | control <1|2>'
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  End Select
  If (n_fail > 0) then
     Write(error_unit,*) 'gpu-obs test FAILED (', trim(mode), ' ', trim(cname), '): ', n_fail, ' failure(s)'
     Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  endif
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
  Write(*,*) 'gpu-obs test passed: ', trim(mode), ' ', trim(cname)

contains

  !> One case (see the header). control = 0: device against host, must agree. 1: device G + 1e-6. 2: weights U on every orbital.
  !> In the control modes the disagreement must be at least 1e3 times the tolerance.
  Subroutine run_case(cs, control)
    Character (len=*), intent(in) :: cs
    Integer, intent(in) :: control
    Integer :: nfl_req, L1, L2, k, nt
    Integer :: nts(3)
    Logical :: symm_c, bilayer
    Real (Kind=Kind(0.d0)) :: h0, uu, u2, w(3)
    Character (len=64) :: lat
    Complex (Kind=Kind(0.d0)), Allocatable :: Gh(:,:,:,:), Gs(:,:,:), Z0(:,:,:), Hbuf(:), Dbuf(:)
    Complex (Kind=Kind(0.d0)), Allocatable :: GT0(:,:,:), G0T(:,:,:), G00(:,:,:), GTT(:,:,:)
    Complex (Kind=Kind(0.d0)), Allocatable :: sT0(:,:,:), s0T(:,:,:), s00(:,:,:), sTT(:,:,:)
    Complex (Kind=Kind(0.d0)) :: ph(3)
    Integer :: hblk(64), nhb
    Character (len=24) :: hlabel(64)
    Type (Lattice) :: Latt
    Type (Unit_cell) :: Latt_unit
    Type (Hopping_Matrix_type), Allocatable :: hop(:)
    Integer, Allocatable :: List(:,:), Invlist(:,:), nphi(:)
    Real (Kind=Kind(0.d0)), Allocatable :: tv(:), t2v(:), tpv(:), cv(:), fx(:), fy(:)
    Integer :: Ndim_l, N_op
    Real (Kind=Kind(0.d0)) :: worst, e
    Logical :: ok

    Select case (cs)
    Case ('sq1');  lat = 'Square'; nfl_req = 1; symm_c = .false.; h0 = 0.d0
    Case ('sq1s'); lat = 'Square'; nfl_req = 1; symm_c = .true.;  h0 = 0.d0
    Case ('sq2');  lat = 'Square'; nfl_req = 2; symm_c = .false.; h0 = 0.d0
    Case ('sq2s'); lat = 'Square'; nfl_req = 2; symm_c = .true.;  h0 = 0.d0
    Case ('sq2h'); lat = 'Square'; nfl_req = 2; symm_c = .true.;  h0 = 0.5d0
    Case ('bil1'); lat = 'Bilayer_honeycomb'; nfl_req = 1; symm_c = .false.; h0 = 0.d0
    Case ('bil2'); lat = 'Bilayer_honeycomb'; nfl_req = 2; symm_c = .true.;  h0 = 0.d0
    Case default
       Write(error_unit,*) 'unknown case ', cs
       Call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
    End Select
    bilayer = (trim(lat) == 'Bilayer_honeycomb')
    symm_run = symm_c
    If (bilayer) then
       L1 = 2; L2 = 2
    else
       L1 = 4; L2 = 4
    endif
    uu = 4.d0; u2 = 2.d0
    nts = [0, 3, 7]
    Call write_parameters(trim(lat), L1, L2, nfl_req, symm_c, uu, u2, h0)
    Call set_QMC_runtime_default_var()
    Call set_ham_name('Hubbard')
    Call set_Nwrap(5)
    Call set_NSweep(1)
    Call set_Ltau(1)
    Call Setup_model(N_op)
    Call Setup_fields(N_op)
    Call ham%Alloc_obs(1)
    Call Setup_update_schemes(N_op, 'info')
    Call Walker_layout(walker, Ltrot, 5)
    Call Walker_alloc(walker)
    Call Setup_gpu(walker, 'info')
    If (.not. gpu_active .or. .not. gpu_stab_on) then
       Write(error_unit,*) 'GPU path or device stabilisation inactive'
       n_fail = n_fail + 1
       Return
    Endif

    ! the hook of Alloc_obs has set the device observables up (Setup_gpu reported them): drop them, so that ham%Obser below is the
    ! host reference; the device is set up again from the replica of the model data after it
    Call gpu_obs_finalize()
    Allocate(Gh(Ndim,Ndim,N_FL,2), Gs(Ndim,Ndim,N_FL), Z0(Ndim,Ndim,N_FL))
    Z0 = (0.d0, 0.d0)
    Do k = 1, 2
       Call fill_random(Gh(:,:,:,k))
    Enddo
    ph(1) = cmplx(1.d0, 0.d0, kind(0.d0)); ph(2) = cmplx(-0.6d0, 0.2d0, kind(0.d0)); ph(3) = cmplx(0.9d0, -0.1d0, kind(0.d0))
    w(1) = 1.d0; w(2) = 0.7d0; w(3) = 1.3d0

    ! ---- equal time, host: two measurements (the sweep hands Symm(G) to Obser)
    Call ham%Init_obs(1)
    Do k = 1, 2
       Call host_obser(Gh(:,:,:,k), ph(k), w(k))
    Enddo
    Allocate(Hbuf(1))
    Call pack_all(Hbuf, nhb, hblk, hlabel, .false.)

    ! model data for the device observables: the same lattice and hopping as Ham_Set builds (Alloc_obs hands the model's own over, see e2e)
    Call Predefined_Latt(lat, L1, L2, Ndim_l, List, Invlist, Latt, Latt_unit)
    Allocate(tv(N_FL), t2v(N_FL), tpv(N_FL), cv(N_FL), fx(N_FL), fy(N_FL), nphi(N_FL))
    tv = 1.d0; t2v = 1.d0; tpv = 1.d0; cv = 0.d0; fx = 0.d0; fy = 0.d0; nphi = 0
    If (bilayer) then
       Call Set_Default_hopping_parameters_Bilayer_honeycomb(hop, tv, t2v, tpv, cv, fx, fy, .true., nphi, N_FL, List, Invlist, Latt, Latt_unit)
    else
       Call Set_Default_hopping_parameters_square(hop, tv, cv, fx, fy, .true., nphi, N_FL, List, Invlist, Latt, Latt_unit)
    endif
    If (control == 2) then
       Call gpu_obs_setup(Latt, Latt_unit, List, Invlist, hop, N_SUN, N_FL, uu, u2, .false., h0, symm_c)   ! U2 ignored
    else
       Call gpu_obs_setup(Latt, Latt_unit, List, Invlist, hop, N_SUN, N_FL, uu, u2, bilayer, h0, symm_c)
    endif
    Call gpu_obs_report()
    If (.not. gpu_obs_active()) then
       Write(error_unit,*) 'device observables inactive'
       n_fail = n_fail + 1
       Return
    Endif

    ! ---- equal time, device: first measurement uploads the host G, the second reads the device G with a zero host array
    Call ham%Init_obs(1)
    g_dev_ok = .false.
    If (control == 1) then
       Gs = Gh(:,:,:,1) + 1.d-6*Gh(:,:,:,2)     ! perturbed host G (uploaded)
       Call gpu_obs_measure(Gs, ph(1), 0, w(1))
    else
       Call gpu_obs_measure(Gh(:,:,:,1), ph(1), 0, w(1))
    endif
    G_d = Gh(:,:,:,2)
    If (control == 1) G_d = G_d + 1.d-6*Gh(:,:,:,1)
    ! H->D: test input, the device copy of G is set by the host
    !$omp target update to(G_d)
    g_dev_ok = .true.
    Call gpu_obs_measure(Z0, ph(2), 0, w(2))
    Call gpu_obs_flush()
    Call pack_all(Dbuf, nbuf, blk, blabel, .true.)
    worst = compare(Hbuf, Dbuf, nhb, hblk, hlabel, 'equal time')

    ! ---- time displaced
    Allocate(GT0(Ndim,Ndim,N_FL), G0T(Ndim,Ndim,N_FL), G00(Ndim,Ndim,N_FL), GTT(Ndim,Ndim,N_FL))
    Allocate(sT0(Ndim,Ndim,N_FL), s0T(Ndim,Ndim,N_FL), s00(Ndim,Ndim,N_FL), sTT(Ndim,Ndim,N_FL))
    ! the matrices are used by the device kernels: present on the device, set by the host (test input)
    !$omp target enter data map(alloc: GT0, G0T, G00, GTT)
    If (.not. gpu_obs_tau_active()) then
       Write(error_unit,*) 'device tau observables inactive'
       n_fail = n_fail + 1
    Else
       Call ham%Init_obs(1)
       seed = 555555555_int64
       Do k = 1, 3
          nt = nts(k)
          Call fill_random(GT0); Call fill_random(G0T); Call fill_random(G00); Call fill_random(GTT)
          If (symm_c) then
             Call Hop_mod_Symm(sT0, GT0, 1); Call Hop_mod_Symm(s0T, G0T, 1)
             Call Hop_mod_Symm(s00, G00, 1); Call Hop_mod_Symm(sTT, GTT, 1)
             Call ham%ObserT(nt, sT0, s0T, s00, sTT, ph(k), w(k))
          else
             Call ham%ObserT(nt, GT0, G0T, G00, GTT, ph(k), w(k))
          endif
       Enddo
       Call pack_all(Hbuf, nhb, hblk, hlabel, .false.)
       Call ham%Init_obs(1)
       seed = 555555555_int64                       ! same matrices again: replay the fill
       Do k = 1, 3
          nt = nts(k)
          Call fill_random(GT0); Call fill_random(G0T); Call fill_random(G00); Call fill_random(GTT)
          If (control == 1) GT0 = GT0 + 1.d-6*G0T
          ! H->D: test input
          !$omp target update to(GT0, G0T, G00, GTT)
          Call gpu_obs_tau_measure(nt, GT0, G0T, G00, GTT, ph(k), w(k))
       Enddo
       Call gpu_obs_flush()
       Call pack_all(Dbuf, nbuf, blk, blabel, .true.)
       e = compare(Hbuf, Dbuf, nhb, hblk, hlabel, 'time displaced')
       worst = max(worst, e)
    Endif
    !$omp target exit data map(delete: GT0, G0T, G00, GTT)

    Select case (control)
    Case (0)
       ok = (worst <= tol)
       Write(*,'(A,A,A,ES10.3,A,ES9.2)') ' gpu-obs ', cs, ': largest relative disagreement ', worst, ', tolerance ', tol
       If (.not. ok) n_fail = n_fail + 1
    Case default
       ok = (worst >= 1.d3*tol)
       Write(*,'(A,I0,A,ES10.3,A,ES9.2,A,L1)') ' gpu-obs negative control ', control, ': disagreement ', worst, ', 1e3 x tolerance ', &
            & 1.d3*tol, ', exceeded ', ok
       If (.not. ok) n_fail = n_fail + 1
    End Select
    Call gpu_obs_finalize()
    Call Finalize_gpu()
    Call Walker_dealloc(walker)
    Call Clear_Ham_state()
  end Subroutine run_case

  !> host Obser on Symm(G) (as the sweep) or on G
  Subroutine host_obser(G, phase, wt)
    Complex (Kind=Kind(0.d0)), intent(in) :: G(Ndim,Ndim,N_FL)
    Complex (Kind=Kind(0.d0)), intent(in) :: phase
    Real (Kind=Kind(0.d0)), intent(in) :: wt
    Complex (Kind=Kind(0.d0)), Allocatable :: Gt(:,:,:)
    Allocate(Gt(Ndim,Ndim,N_FL))
    If (symm_run) then
       Call Hop_mod_Symm(Gt, G, 1)
       Call ham%Obser(Gt, phase, 1, wt)
    else
       Call ham%Obser(G, phase, 1, wt)
    endif
  end Subroutine host_obser

  !> random complex entries in (-0.5, 0.5)
  Subroutine fill_random(A)
    Complex (Kind=Kind(0.d0)), intent(out) :: A(:,:,:)
    Integer :: i, j, f
    Real (Kind=Kind(0.d0)) :: x, y
    Do f = 1, size(A,3)
       Do j = 1, size(A,2)
          Do i = 1, size(A,1)
             x = rnd01(seed) - 0.5d0
             y = rnd01(seed) - 0.5d0
             A(i,j,f) = cmplx(x, y, kind(0.d0))
          Enddo
       Enddo
    Enddo
  end Subroutine fill_random

  !> every counter and accumulator of the Obs_* into one vector; blocks (start index) and labels in nb, bs, bl
  Subroutine pack_all(v, nb_out, bs, bl, is_dev)
    Complex (Kind=Kind(0.d0)), Allocatable, intent(inout) :: v(:)
    Integer, intent(out) :: nb_out, bs(64)
    Character (len=24), intent(out) :: bl(64)
    Logical, intent(in) :: is_dev
    Integer :: i
    Character (len=24) :: nm
    nblk = 0
    If (allocated(buf)) deallocate(buf)
    Allocate(buf(0))
    Do i = 1, size(Obs_scal,1)
       Write(nm,'(A,I0)') 'scal ', i
       Call push_start(nm, bs, bl)
       Call push([cmplx(dble(Obs_scal(i)%N), Obs_scal(i)%Ave_sign, kind(0.d0))])
       Call push(Obs_scal(i)%Obs_vec)
    Enddo
    If (allocated(Obs_local)) then
       Call push_start('local', bs, bl)
       Call push([cmplx(dble(Obs_local(1)%N), Obs_local(1)%Ave_sign, kind(0.d0))])
       Call push(reshape(Obs_local(1)%Obs_Latt, [size(Obs_local(1)%Obs_Latt)]))
    Endif
    Do i = 1, size(Obs_eq,1)
       Write(nm,'(A,A)') 'eq ', trim(Obs_eq(i)%File_Latt)
       If (i == 3 .and. size(Obs_eq,1) == 5) nm = 'eq SpinXY'
       If (i == 4 .and. size(Obs_eq,1) == 5) nm = 'eq SpinT'
       Call push_start(nm, bs, bl)
       Call push([cmplx(dble(Obs_eq(i)%N), Obs_eq(i)%Ave_sign, kind(0.d0))])
       Call push(reshape(Obs_eq(i)%Obs_Latt, [size(Obs_eq(i)%Obs_Latt)]))
       Call push(Obs_eq(i)%Obs_Latt0)
    Enddo
    If (allocated(Obs_tau)) then
       Do i = 1, size(Obs_tau,1)
          Write(nm,'(A,I0,A,A)') 'tau ', i, ' ', trim(Obs_tau(i)%File_Latt)
          Call push_start(nm, bs, bl)
          Call push([cmplx(dble(Obs_tau(i)%N), Obs_tau(i)%Ave_sign, kind(0.d0))])
          Call push(reshape(Obs_tau(i)%Obs_Latt, [size(Obs_tau(i)%Obs_Latt)]))
          Call push(Obs_tau(i)%Obs_Latt0)
       Enddo
    Endif
    nb_out = nblk
    bs(nblk+1) = size(buf) + 1
    If (allocated(v)) deallocate(v)
    Allocate(v(size(buf)))
    v = buf
  end Subroutine pack_all

  Subroutine push_start(label, bs, bl)
    Character (len=*), intent(in) :: label
    Integer, intent(inout) :: bs(64)
    Character (len=24), intent(inout) :: bl(64)
    nblk = nblk + 1
    bs(nblk) = size(buf) + 1
    bl(nblk) = label
  end Subroutine push_start

  Subroutine push(x)
    Complex (Kind=Kind(0.d0)), intent(in) :: x(:)
    Complex (Kind=Kind(0.d0)), Allocatable :: t(:)
    Allocate(t(size(buf) + size(x)))
    t(1:size(buf)) = buf
    t(size(buf)+1:) = x
    Call move_alloc(t, buf)
  end Subroutine push

  !> largest per-block relative difference max|D-H| / max(max|H|, 1e-3); printed per block
  Real (Kind=Kind(0.d0)) Function compare(H, D, nb, bs, bl, what)
    Complex (Kind=Kind(0.d0)), intent(in) :: H(:), D(:)
    Integer, intent(in) :: nb, bs(64)
    Character (len=24), intent(in) :: bl(64)
    Character (len=*), intent(in) :: what
    Integer :: b, a, z
    Real (Kind=Kind(0.d0)) :: dm, hm, r
    compare = 0.d0
    If (size(H) /= size(D)) then
       Write(error_unit,*) 'compare: packed sizes differ ', size(H), size(D)
       compare = huge(1.d0)
       Return
    Endif
    Do b = 1, nb
       a = bs(b); z = bs(b+1) - 1
       If (z < a) cycle
       dm = maxval(abs(D(a:z) - H(a:z)))
       hm = maxval(abs(H(a:z)))
       r = dm/max(hm, 1.d-3)
       compare = max(compare, r)
       If (r > 1.d-13) Write(*,'(A,A,A,A,A,ES10.3,A,ES10.3)') ' ', what, ' ', trim(bl(b)), ': max|D-H| ', dm, ', relative ', r
    Enddo
  end Function compare

  !> uniform number in [0,1)
  Real (Kind=Kind(0.d0)) function rnd01(s)
    Integer(int64), intent(inout) :: s
    s = mod(1103515245_int64*s + 12345_int64, 2147483648_int64)
    rnd01 = real(s, kind(0.d0))/2147483648.d0
  end function rnd01

  Subroutine write_parameters(lat, L1, L2, nfl_req, symm_c, uu, u2, h0)
    Character (len=*), intent(in) :: lat
    Integer, intent(in) :: L1, L2, nfl_req
    Logical, intent(in) :: symm_c
    Real (Kind=Kind(0.d0)), intent(in) :: uu, u2, h0
    Integer :: u, i
    Open(newunit=u, file='parameters', status='replace', action='write')
    Write(u,'(A)') '&VAR_lattice'
    Write(u,'(A,I0)') 'L1 = ', L1
    Write(u,'(A,I0)') 'L2 = ', L2
    Write(u,'(A)') 'Lattice_type = "'//lat//'"'
    Write(u,'(A)') 'Model = "Hubbard"'
    Write(u,'(A)') '/'
    Write(u,'(A)') '&VAR_Model_Generic'
    Write(u,'(A)') 'Checkerboard = .T.'
    If (symm_c) then
       Write(u,'(A)') 'Symm = .T.'
    else
       Write(u,'(A)') 'Symm = .F.'
    endif
    Write(u,'(A)') 'N_SUN = 2'
    Write(u,'(A,I0)') 'N_FL = ', nfl_req
    Write(u,'(A)') 'Dtau = 0.1d0'
    Write(u,'(A)') 'Beta = 2.d0'
    Write(u,'(A)') 'Projector = .F.'
    Write(u,'(A)') '/'
    Write(u,'(A)') '&VAR_Hubbard'
    If (nfl_req == 2) then
       Write(u,'(A)') 'Mz = .T.'
    else
       Write(u,'(A)') 'Mz = .F.'
    endif
    Write(u,'(A)') 'Continuous = .F.'
    Write(u,'(A)') 'ham_T = 1.d0'
    Write(u,'(A)') 'ham_T2 = 1.d0'
    Write(u,'(A)') 'ham_Tperp = 1.d0'
    Write(u,'(A)') 'ham_chem = 0.d0'
    Write(u,'(A,F6.2)') 'ham_U = ', uu
    Write(u,'(A,F6.2)') 'ham_U2 = ', u2
    Write(u,'(A,F6.2)') 'ham_h0 = ', h0
    Write(u,'(A)') '/'
    Close(u)
    Open(newunit=u, file='seeds', status='replace', action='write')
    Do i = 1, 12
       Write(u,'(I10)') 1000003 * i + 17
    Enddo
    Close(u)
  end Subroutine write_parameters

end Program Test_gpu_obs
