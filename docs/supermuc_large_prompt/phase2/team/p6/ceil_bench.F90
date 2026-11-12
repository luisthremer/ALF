! Phase 6 gate 1: zgemm ceiling on ONE PVC tile through ALF's alf_la (ONEMKL backend, same dispatch as the ALF wraps).
! One MPI rank = one process with its own device (rank mod ndev, as ALF). All ranks run every configuration in lockstep
! (MPI_Barrier before and after), so with 4 ranks per tile and CCS 4 the per-tile aggregate is the sum over the ranks of the tile.
! The calls are host-synchronous (omp dispatch), timing = omp_get_wtime over many repetitions + la_sync.
! Configurations (N = 576, 1024): single (back-to-back), wrap-like (gemm to scratch + copy kernel, as hopd_left),
!   LEFT  pattern C_b = H M_b : separate gemms, strided batch (sA = 0), wide gemm (n = B N)
!   RIGHT pattern C_b = M_b H : separate gemms, strided batch (sB = 0), tall gemm (layout (i,b,j), m = B N)
!   DD pattern  C_b = A_b M_b : strided batch with both operands distinct (CGR/stabilisation products)
! Bitwise check of every batched/wide/tall result against the separate-gemm result of the same inputs.
! env: CB_N (list given by CB_N1, CB_N2), CB_FLOPS (flops per config per rank, 3e12), CB_BMAX (16)
program ceil_bench
  use alf_la_mod
  use omp_lib
  use mpi
  implicit none
  integer, parameter :: nbl = 5
  integer :: bl(nbl) = [2, 4, 8, 16, 1]
  integer :: rank, nranks, ierr, nn, ib, ni, bmax, nd, ndev
  integer :: nlist(2)
  real(8) :: flops_cfg
  complex(8), parameter :: z1 = (1.d0,0.d0), z0 = (0.d0,0.d0)
  complex(8), allocatable, target :: H(:,:), Mw(:,:,:), Cw(:,:,:), Rw(:,:,:), Mt(:,:,:), Ct(:,:,:), Rt(:,:,:), Aw(:,:,:), Sc(:,:)
  complex(8), allocatable, target :: Rs1(:,:,:), Rw2(:,:,:), Rt1(:,:,:)
  character(len=32) :: s

  call MPI_Init(ierr); call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr); call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  ndev = max(omp_get_num_devices(), 1)
  call omp_set_default_device(mod(rank, ndev))
  nlist = [ienv('CB_N1', 576), ienv('CB_N2', 1024)]
  bmax = ienv('CB_BMAX', 16)
  flops_cfg = 3.d12
  call get_environment_variable('CB_FLOPS', s, status=ierr); if (ierr == 0) read(s,*) flops_cfg
  if (rank == 0) write(*,'(A,I4,A,I3,A,ES9.2)') 'CEIL start nranks', nranks, ' ndev', ndev, ' flops/cfg/rank', flops_cfg
  call la_init()
  do ni = 1, 2
    nn = nlist(ni)
    if (nn <= 0) cycle
    call run_n(nn)
  end do
  call la_finalize()
  call MPI_Finalize(ierr)
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
  subroutine fill3(A, sd)
    complex(8), intent(out) :: A(:,:,:)
    integer, intent(in) :: sd
    integer(8) :: r
    integer :: i, j, k
    real(8) :: a1, a2
    r = sd
    do k = 1, size(A,3)
      do j = 1, size(A,2)
        do i = 1, size(A,1)
          r = mod(r * 1103515245_8 + 12345_8, 2147483648_8); a1 = dble(r) / 2147483648.d0 - 0.5d0
          r = mod(r * 1103515245_8 + 12345_8, 2147483648_8); a2 = dble(r) / 2147483648.d0 - 0.5d0
          A(i,j,k) = cmplx(a1, a2, 8)
        end do
      end do
    end do
  end subroutine
  subroutine k_cpy(nel, a, b)
    integer, intent(in) :: nel
    complex(8), intent(in) :: a(*)
    complex(8), intent(inout) :: b(*)
    integer :: i
    !$omp target teams distribute parallel do
    do i = 1, nel
      b(i) = a(i)
    end do
  end subroutine
  integer(8) function nbit(X, R, nel)     ! number of elements with any differing bit (host arrays)
    complex(8), intent(in) :: X(*), R(*)
    integer, intent(in) :: nel
    integer :: i
    nbit = 0
    do i = 1, nel
      if (transfer(dble(X(i)), 0_8) /= transfer(dble(R(i)), 0_8) .or. transfer(aimag(X(i)), 0_8) /= transfer(aimag(R(i)), 0_8)) nbit = nbit + 1
    end do
  end function
  subroutine rep(cfg, N, B, reps, t, fl, nb)
    character(len=*), intent(in) :: cfg
    integer, intent(in) :: N, B, reps
    real(8), intent(in) :: t, fl
    integer(8), intent(in) :: nb
    write(*,'(A,I4,A,I5,1X,A,I3,A,I6,A,F9.5,A,F8.3,A,I9,A,F9.2)') 'CEIL rank', rank, ' N', N, trim(cfg)//' B', B, ' reps', reps, &
         ' t', t, ' TFLOPs', fl*reps/t/1.d12, ' bitdiff', nb, ' usgemm', 1.d6*t/(reps*max(B,1))
  end subroutine
  subroutine run_n(N)
    integer, intent(in) :: N
    integer :: B, bi, b1, reps, r, ib2, ki, kk
    integer :: kl(2) = [8, 32]
    real(8) :: fk
    character(len=16) :: cfg
    real(8) :: t0, t1, fl1
    integer(8) :: nb
    fl1 = 8.d0 * dble(N)**3
    allocate(Rs1(N,N,bmax), Rw2(N,N,bmax), Rt1(N*bmax,N,1))
    allocate(H(N,N), Mw(N,N,bmax), Cw(N,N,bmax), Rw(N,N,bmax), Mt(N*bmax,N,1), Ct(N*bmax,N,1), Rt(N*bmax,N,1), Aw(N,N,bmax), Sc(N,N*bmax))
    call fill3(Mw, 11 + rank); call fill3(Aw, 77 + rank); Cw = z0; Rw = z0; Ct = z0; Rt = z0; Sc = z0
    call fill1(H, 5 + rank)
    call fill3(Mt, 23 + rank)
    !$omp target enter data map(to: H, Mw, Cw, Rw, Mt, Ct, Rt, Aw, Sc, Rs1, Rw2, Rt1)
    ! ---- single, back-to-back (N^3 gemm, C = H M_1)
    reps = max(10, int(flops_cfg / fl1))
    do r = 1, 5
      call la_gemm(LA_N, LA_N, N, N, N, z1, H, N, Mw(1,1,1), N, z0, Cw(1,1,1), N)
    end do
    call la_sync()
    call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
    do r = 1, reps
      call la_gemm(LA_N, LA_N, N, N, N, z1, H, N, Mw(1,1,1), N, z0, Cw(1,1,1), N)
    end do
    call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call rep('single      ', N, 1, reps, t1 - t0, fl1, 0_8)
    ! ---- wrap-like single: gemm to scratch + copy kernel (hopd_left)
    call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
    do r = 1, reps
      call la_gemm(LA_N, LA_N, N, N, N, z1, H, N, Mw(1,1,1), N, z0, Sc(1,1), N)
      call k_cpy(N*N, Sc(1,1), Cw(1,1,1))
    end do
    call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call rep('wraplike1   ', N, 1, reps, t1 - t0, fl1, 0_8)
    ! ---- reference results (separate gemms, all bmax walkers)
    do b1 = 1, bmax
      call la_gemm(LA_N, LA_N, N, N, N, z1, H, N, Mw(1,1,b1), N, z0, Rw(1,1,b1), N)       ! left:  R_b = H M_b
    end do
    do b1 = 1, bmax
      call la_gemm(LA_N, LA_N, N, N, N, z1, Mt(1+(b1-1)*N,1,1), N*bmax, H, N, z0, Rt(1+(b1-1)*N,1,1), N*bmax)   ! right: R_b = M_b H in layout (i,b,j)
    end do
    Rs1 = z0; Rw2 = z0; Rt1 = z0
    !$omp target update to(Rs1, Rw2, Rt1)
    do b1 = 1, bmax
      call la_gemm_batched(LA_N, LA_N, N, N, N, z1, H, N, 0, Mw(1,1,b1), N, N*N, z0, Rs1(1,1,b1), N, N*N, 1)
      call la_gemm_batched(LA_N, LA_N, N, N, N, z1, Mt(1+(b1-1)*N,1,1), N*bmax, N, H, N, 0, z0, Rt1(1+(b1-1)*N,1,1), N*bmax, N, 1)
    end do
    do b1 = 1, bmax - 1, 2
      call la_gemm(LA_N, LA_N, N, N*2, N, z1, H, N, Mw(1,1,b1), N, z0, Rw2(1,1,b1), N)
    end do
    call la_sync()
    !$omp target update from(Rw, Rt, Rs1, Rw2, Rt1)
    do bi = 1, nbl
      B = bl(bi)
      if (B > bmax) cycle
      reps = max(5, int(flops_cfg / (fl1 * B)))
      ! ---- LEFT separate gemms (B back-to-back)
      call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
      do r = 1, reps
        do b1 = 1, B
          call la_gemm(LA_N, LA_N, N, N, N, z1, H, N, Mw(1,1,b1), N, z0, Cw(1,1,b1), N)
        end do
      end do
      call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call rep('L_separate  ', N, B, reps, t1 - t0, fl1*B, 0_8)
      if (B == 1) cycle
      ! ---- LEFT strided batch, shared A
      Cw = z0
      !$omp target update to(Cw)
      call la_gemm_batched(LA_N, LA_N, N, N, N, z1, H, N, 0, Mw, N, N*N, z0, Cw, N, N*N, B); call la_sync()
      !$omp target update from(Cw)
      nb = nbit(Cw, Rw, N*N*B)
      call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
      do r = 1, reps
        call la_gemm_batched(LA_N, LA_N, N, N, N, z1, H, N, 0, Mw, N, N*N, z0, Cw, N, N*N, B)
      end do
      call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call rep('L_strided   ', N, B, reps, t1 - t0, fl1*B, nb)
      call rep('L_str_vs_nb1', N, B, 1, 1.d0, 0.d0, nbit(Cw, Rs1, N*N*B))
      ! ---- LEFT wide gemm n = B N
      Cw = z0
      !$omp target update to(Cw)
      call la_gemm(LA_N, LA_N, N, N*B, N, z1, H, N, Mw(1,1,1), N, z0, Cw(1,1,1), N); call la_sync()
      !$omp target update from(Cw)
      nb = nbit(Cw, Rw, N*N*B)
      call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
      do r = 1, reps
        call la_gemm(LA_N, LA_N, N, N*B, N, z1, H, N, Mw(1,1,1), N, z0, Cw(1,1,1), N)
      end do
      call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call rep('L_wide      ', N, B, reps, t1 - t0, fl1*B, nb)
      call rep('L_wide_vs_B2', N, B, 1, 1.d0, 0.d0, nbit(Cw, Rw2, N*N*B))
      ! ---- wrap-like batched: wide gemm to scratch + one copy kernel (what the batched wrap would do)
      call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
      do r = 1, reps
        call la_gemm(LA_N, LA_N, N, N*B, N, z1, H, N, Mw(1,1,1), N, z0, Sc(1,1), N)
        call k_cpy(N*N*B, Sc(1,1), Cw(1,1,1))
      end do
      call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call rep('wraplikeWide', N, B, reps, t1 - t0, fl1*B, 0_8)
      ! ---- RIGHT separate gemms in layout (i,b,j)
      call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
      do r = 1, reps
        do b1 = 1, B
          call la_gemm(LA_N, LA_N, N, N, N, z1, Mt(1+(b1-1)*N,1,1), N*bmax, H, N, z0, Ct(1+(b1-1)*N,1,1), N*bmax)
        end do
      end do
      call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call rep('R_separate  ', N, B, reps, t1 - t0, fl1*B, 0_8)
      ! ---- RIGHT strided batch (layout (i,b,j): sA = N, lda = N*bmax), shared B
      Ct = z0
      !$omp target update to(Ct)
      call la_gemm_batched(LA_N, LA_N, N, N, N, z1, Mt, N*bmax, N, H, N, 0, z0, Ct, N*bmax, N, B); call la_sync()
      !$omp target update from(Ct)
      nb = nbit_t(Ct, Rt, N, B)
      call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
      do r = 1, reps
        call la_gemm_batched(LA_N, LA_N, N, N, N, z1, Mt, N*bmax, N, H, N, 0, z0, Ct, N*bmax, N, B)
      end do
      call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call rep('R_strided   ', N, B, reps, t1 - t0, fl1*B, nb)
      call rep('R_str_vs_nb1', N, B, 1, 1.d0, 0.d0, nbit_t(Ct, Rt1, N, B))
      ! ---- RIGHT tall gemm: rows (i,b) with m = N*B and lda = N*bmax (rows of walkers > B untouched)
      if (B == bmax) then
        Ct = z0
        !$omp target update to(Ct)
        call la_gemm(LA_N, LA_N, N*B, N, N, z1, Mt(1,1,1), N*bmax, H, N, z0, Ct(1,1,1), N*bmax); call la_sync()
        !$omp target update from(Ct)
        nb = nbit_t(Ct, Rt, N, B)
        call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
        do r = 1, reps
          call la_gemm(LA_N, LA_N, N*B, N, N, z1, Mt(1,1,1), N*bmax, H, N, z0, Ct(1,1,1), N*bmax)
        end do
        call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call rep('R_tall      ', N, B, reps, t1 - t0, fl1*B, nb)
      end if
      ! ---- DD strided batch: both operands distinct
      call la_gemm_batched(LA_N, LA_N, N, N, N, z1, Aw, N, N*N, Mw, N, N*N, z0, Cw, N, N*N, B); call la_sync()
      call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
      do r = 1, reps
        call la_gemm_batched(LA_N, LA_N, N, N, N, z1, Aw, N, N*N, Mw, N, N*N, z0, Cw, N, N*N, B)
      end do
      call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call rep('DD_strided  ', N, B, reps, t1 - t0, fl1*B, 0_8)
    end do
    ! ---- rank-k update pattern of the submatrix flush: C_b -= X_b Y_b^T (N x N, k small), all operands distinct per walker
    do ki = 1, 2
      kk = kl(ki)
      fk = 8.d0 * dble(N)**2 * kk
      do bi = 1, nbl
        B = bl(bi)
        if (B > bmax) cycle
        reps = max(20, min(20000, int(flops_cfg / (fk * B))))
        call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
        do r = 1, reps
          do b1 = 1, B
            call la_gemm(LA_N, LA_T, N, N, kk, -z1, Aw(1,1,b1), N, Mw(1,1,b1), N, z1, Cw(1,1,b1), N)
          end do
        end do
        call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
        write(cfg,'(A,I2.2)') 'K_separate_k', kk
        call rep(cfg, N, B, reps, t1 - t0, fk*B, 0_8)
        if (B == 1) cycle
        call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = omp_get_wtime()
        do r = 1, reps
          call la_gemm_batched(LA_N, LA_T, N, N, kk, -z1, Aw, N, N*N, Mw, N, N*N, z1, Cw, N, N*N, B)
        end do
        call la_sync(); t1 = omp_get_wtime(); call MPI_Barrier(MPI_COMM_WORLD, ierr)
        write(cfg,'(A,I2.2)') 'K_strided_k', kk
        call rep(cfg, N, B, reps, t1 - t0, fk*B, 0_8)
      end do
    end do
    !$omp target exit data map(delete: H, Mw, Cw, Rw, Mt, Ct, Rt, Aw, Sc, Rs1, Rw2, Rt1)
    deallocate(H, Mw, Cw, Rw, Mt, Ct, Rt, Aw, Sc, Rs1, Rw2, Rt1)
  end subroutine
  integer(8) function nbit_t(X, R, N, B)    ! layout (i,b,j) with b-stride N, leading dim N*bmax: first B walkers only
    complex(8), intent(in) :: X(*), R(*)
    integer, intent(in) :: N, B
    integer :: i, j, k
    nbit_t = 0
    do j = 1, N
      do k = 1, N*B
        i = k + (j-1)*N*bmax
        if (transfer(dble(X(i)), 0_8) /= transfer(dble(R(i)), 0_8) .or. transfer(aimag(X(i)), 0_8) /= transfer(aimag(R(i)), 0_8)) nbit_t = nbit_t + 1
      end do
    end do
  end function
  subroutine fill1(A, sd)
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
    end do
  end subroutine
end program ceil_bench
