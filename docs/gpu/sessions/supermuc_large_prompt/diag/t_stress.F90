! Stress test of the alf_la ONEMKL chain used by gpu_stab_mod dev_cgr / cgr_solve (NVAR=1 and NVAR=2), to catch sporadic
! garbage results (production aborts "G difference > 10", 2026-11-12). Same inputs every iteration; iteration 1 is the
! reference; every later iteration is compared on the device (max |X - Xref| / max |Xref|) and bitwise.
! env: STRESS_N (1024), STRESS_SEC (wall seconds, 600), STRESS_SYNC=1 (la_sync after every alf_la call),
!      STRESS_PLANT=k (positive control: corrupt one T entry by 1.0 in iteration k), STRESS_SEED.
! Every bad iteration prints a line 'BAD rank it what rel nbits'; the end prints 'STRESS_RESULT ...'.
program t_stress
  use alf_la_mod
  use omp_lib
  implicit none
  integer :: n, sec, plant, it, ios, rank, nbad_bit, nbad_tol, nbad_big, i, j, seed
  logical :: dosync
  complex(8), allocatable, target :: UL(:,:), UR(:,:), VL(:,:), VR(:,:), sA(:,:), sB(:,:), sC(:,:), G1(:,:), G2(:,:)
  complex(8), allocatable, target :: rA(:,:), rG1(:,:), rG2(:,:), tau(:)
  integer(4), allocatable, target :: ipiv(:), ipvt(:), ripiv(:), ripvt(:)
  complex(8), parameter :: z1 = (1.d0,0.d0), z0 = (0.d0,0.d0)
  real(8) :: t0, rel, maxrel(3)
  integer :: nb
  character(len=32) :: s
  character(len=8) :: what(3) = ['LU     ', 'G_NVAR1', 'G_NVAR2']

  n = ienv('STRESS_N', 1024); sec = ienv('STRESS_SEC', 600); plant = ienv('STRESS_PLANT', -1)
  dosync = ienv('STRESS_SYNC', 0) == 1; seed = ienv('STRESS_SEED', 12345)
  rank = ienv('PMI_RANK', ienv('SLURM_PROCID', 0))
  ! device of rank r = r mod (number of devices), as ALF does (production jobfile: 32 ranks/node, 4 per tile)
  call omp_set_default_device(mod(ienv('STRESS_DEV', rank), max(omp_get_num_devices(), 1)))
  write(*,'(A,I5,A,I3,A,I3,A,I3)') 'STRESS rank', rank, ' device', omp_get_default_device(), ' of', omp_get_num_devices(), &
       ' threads', omp_get_max_threads()
  allocate(UL(n,n), UR(n,n), VL(n,n), VR(n,n), sA(n,n), sB(n,n), sC(n,n), G1(n,n), G2(n,n), rA(n,n), rG1(n,n), rG2(n,n))
  allocate(tau(n), ipiv(n), ipvt(n), ripiv(n), ripvt(n))
  call fill(UL, seed); call fill(UR, seed+1); call fill(VL, seed+2); call fill(VR, seed+3)
  sA = z0; sB = z0; sC = z0; G1 = z0; G2 = z0; rA = z0; rG1 = z0; rG2 = z0; tau = z0; ipiv = 0; ipvt = 0
  call la_init()
  !$omp target enter data map(to: UL, UR, VL, VR, sA, sB, sC, G1, G2, rA, rG1, rG2, tau, ipiv, ipvt)
  nbad_bit = 0; nbad_tol = 0; nbad_big = 0; maxrel = 0.d0
  t0 = omp_get_wtime(); it = 0
  do
    it = it + 1
    ! ---- dev_cgr: RHS = UR^H UL (sA), T = VR VL (sB) + sA (stand-in for TPUP), LU of sA
    call la_gemm(LA_C, LA_N, n, n, n, z1, UR, n, UL, n, z0, sA, n); call sy()
    call la_gemm(LA_N, LA_N, n, n, n, z1, VR, n, VL, n, z0, sB, n); call sy()
    call k_add(sB, sA)
    if (it == plant) call k_poison(sB)
    call la_getrf(n, n, sA, n, ipiv); call sy()
    call cmp(1, it)                             ! LU factors of sA
    ! ---- cgr_solve NVAR=1 on T=sB (qrpp, ormqr L C, trsm L, gemm)
    call k_copy(sC, sB)                         ! keep T for NVAR=2 (sC = T)
    call la_qrpp(n, n, sB, n, ipvt, tau); call sy()
    call k_ct(G1, UR)                           ! RHS = UR^H
    call la_ormqr('L', LA_C, n, n, n, sB, n, tau, G1, n); call sy()
    call la_trsm('L', 'U', LA_N, 'N', n, n, z1, sB, n, G1, n); call sy()
    call la_gemm(LA_N, LA_N, n, n, n, z1, UL, n, G1, n, z0, sA, n); call sy()
    call k_copy(G1, sA)
    ! ---- NVAR=2 on T^H (qrpp, ormqr R N, trsm R C, gemm)
    call k_ct(sB, sC)
    call la_qrpp(n, n, sB, n, ipvt, tau); call sy()
    call k_copy(G2, UL)
    call la_ormqr('R', LA_N, n, n, n, sB, n, tau, G2, n); call sy()
    call la_trsm('R', 'U', LA_C, 'N', n, n, z1, sB, n, G2, n); call sy()
    call la_gemm(LA_N, LA_C, n, n, n, z1, G2, n, UR, n, z0, sA, n); call sy()
    call k_copy(G2, sA)
    call cmp(2, it); call cmp(3, it)
    if (omp_get_wtime() - t0 > sec) exit
  end do
  write(*,'(A,I5,A,I8,A,F8.3,A,3I7,A,3ES10.2)') 'STRESS_RESULT rank', rank, ' iters', it, ' s/iter', (omp_get_wtime()-t0)/it, &
       ' nbad bit/1e-8/1e-3', nbad_bit, nbad_tol, nbad_big, ' maxrel LU,G1,G2', maxrel
  call la_finalize()
contains
  integer function ienv(nm, def)
    character(len=*), intent(in) :: nm
    integer, intent(in) :: def
    character(len=32) :: v
    integer :: st
    call get_environment_variable(nm, v, status=st)
    ienv = def
    if (st == 0 .and. len_trim(v) > 0) read(v, *, iostat=st) ienv
    if (st /= 0) ienv = def
  end function
  subroutine sy()
    if (dosync) call la_sync()
  end subroutine
  subroutine fill(A, sd)
    complex(8), intent(out) :: A(:,:)
    integer, intent(in) :: sd
    integer(8) :: r
    integer :: i, j
    real(8) :: a1, a2
    r = sd
    do j = 1, size(A,2)
      do i = 1, size(A,1)
        r = mod(r * 1103515245_8 + 12345_8, 2147483648_8); a1 = dble(r) / 2147483648.d0 - 0.5d0
        r = mod(r * 1103515245_8 + 12345_8, 2147483648_8); a2 = dble(r) / 2147483648.d0 - 0.5d0
        A(i,j) = cmplx(a1, a2, 8)
      end do
      A(j,j) = A(j,j) + 4.d0          ! well conditioned
    end do
  end subroutine
  subroutine k_add(X, Y)
    complex(8) :: X(:,:), Y(:,:)
    integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
      do i = 1, n
        X(i,j) = X(i,j) + Y(i,j)
      end do
    end do
  end subroutine
  subroutine k_copy(X, Y)
    complex(8) :: X(:,:), Y(:,:)
    integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
      do i = 1, n
        X(i,j) = Y(i,j)
      end do
    end do
  end subroutine
  subroutine k_ct(X, Y)
    complex(8) :: X(:,:), Y(:,:)
    integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, n
      do i = 1, n
        X(i,j) = conjg(Y(j,i))
      end do
    end do
  end subroutine
  subroutine k_poison(X)
    complex(8) :: X(:,:)
    !$omp target
    X(n/2, n/3) = X(n/2, n/3) + (1.d0, 0.d0)
    !$omp end target
  end subroutine
  ! compare X (1: LU factors sA + ipiv, 2: G1, 3: G2) with its reference; iteration 1 stores the reference
  subroutine cmp(k, it)
    integer, intent(in) :: k, it
    real(8) :: dmax, rmax
    integer :: nbits, i, j
    if (it == 1) then
      select case (k)
      case (1); call k_copy(rA, sA)
      case (2); call k_copy(rG1, G1)
      case (3); call k_copy(rG2, G2)
      end select
      return
    end if
    dmax = 0.d0; rmax = 0.d0; nbits = 0
    select case (k)
    case (1); call dcmp(sA, rA, dmax, rmax, nbits)
    case (2); call dcmp(G1, rG1, dmax, rmax, nbits)
    case (3); call dcmp(G2, rG2, dmax, rmax, nbits)
    end select
    rel = dmax / max(rmax, tiny(1.d0))
    maxrel(k) = max(maxrel(k), rel)
    if (nbits > 0) nbad_bit = nbad_bit + 1
    if (rel > 1.d-8) nbad_tol = nbad_tol + 1
    if (rel > 1.d-3) nbad_big = nbad_big + 1
    if (nbits > 0) write(*,'(A,I5,A,I8,1X,A,A,ES10.2,A,I9,A,F10.1)') 'BAD rank', rank, ' it', it, what(k), ' rel', rel, &
         ' nbits', nbits, ' t', omp_get_wtime() - t0
  end subroutine
  subroutine dcmp(X, R, dmax, rmax, nbits)
    complex(8) :: X(:,:), R(:,:)
    real(8), intent(inout) :: dmax, rmax
    integer, intent(inout) :: nbits
    integer :: i, j
    !$omp target teams distribute parallel do collapse(2) reduction(max: dmax, rmax) reduction(+: nbits) map(tofrom: dmax, rmax, nbits)
    do j = 1, n
      do i = 1, n
        dmax = max(dmax, abs(X(i,j) - R(i,j)))
        rmax = max(rmax, abs(R(i,j)))
        if (transfer(dble(X(i,j)), 0_8) /= transfer(dble(R(i,j)), 0_8) .or. &
            transfer(aimag(X(i,j)), 0_8) /= transfer(aimag(R(i,j)), 0_8)) nbits = nbits + 1
      end do
    end do
  end subroutine
end program t_stress
