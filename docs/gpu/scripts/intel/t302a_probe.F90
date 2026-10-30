! 302a development probe (queue decision): OpenMP interop sycl::queue, dispatch BLAS, shim LAPACK, ordering.
include "mkl_blas_omp_offload.f90"
program t302a_probe
  use iso_c_binding
  use omp_lib, only: omp_interop_kind, omp_get_interop_ptr, omp_get_interop_int, omp_ipr_fr_id, omp_ipr_targetsync, &
       omp_interop_none
  use onemkl_blas_omp_offload_lp64, only: zgemm
  implicit none
  interface
    integer(c_int) function alfsh_init(q) bind(C, name='alfsh_init')
      import :: c_ptr, c_int
      type(c_ptr), value :: q
    end function
    integer(c_int) function alfsh_describe() bind(C, name='alfsh_describe')
      import :: c_int
    end function
    integer(c_int) function alfsh_ptr_kind(p) bind(C, name='alfsh_ptr_kind')
      import :: c_ptr, c_int
      type(c_ptr), value :: p
    end function
    integer(c_int) function alfsh_wait() bind(C, name='alfsh_wait')
      import :: c_int
    end function
    integer(c_int) function alfsh_zgeqrf_ws(m, n, lda, ws) bind(C, name='alfsh_zgeqrf_ws')
      import :: c_int, c_int64_t
      integer(c_int64_t), value :: m, n, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_zgeqrf(m, n, a, lda, tau, work, lwork) bind(C, name='alfsh_zgeqrf')
      import :: c_ptr, c_int, c_int64_t
      type(c_ptr), value :: a, tau, work
      integer(c_int64_t), value :: m, n, lda, lwork
    end function
  end interface
  integer, parameter :: N = 1024
  integer(omp_interop_kind) :: iop
  type(c_ptr) :: qh
  integer :: rc, i, j, rep, mode
  integer(c_intptr_t) :: fr
  complex(8), allocatable, target :: A(:,:), B(:,:), C(:,:), Cr(:,:), tau(:), w(:), Aq(:,:), Aqr(:,:)
  complex(8), allocatable :: hw(:)
  integer(c_int64_t) :: ws
  real(8) :: err, bsum, worst, accd(24)
  complex(8) :: one = (1.d0, 0.d0), zero = (0.d0, 0.d0), wq(1)
  character(len=8) :: arg, arg2
  integer :: stage
  external :: zgeqrf
  call get_command_argument(1, arg); call get_command_argument(2, arg2)
  stage = 99
  if (len_trim(arg) > 0) read(arg,*) stage
  mode = 0
  if (trim(arg2) == 'own') mode = 1

  iop = omp_interop_none
  if (mode == 0) then
    !$omp interop init(prefer_type("sycl"), targetsync: iop)
    print *, 'iop =', iop
    fr = omp_get_interop_int(iop, omp_ipr_fr_id, rc)
    print *, 'fr_id =', fr, ' rc =', rc, ' (4 = sycl)'
    qh = omp_get_interop_ptr(iop, omp_ipr_targetsync, rc)
    print *, 'targetsync ptr null? ', .not. c_associated(qh), ' rc =', rc
  else
    qh = c_null_ptr
  end if
  print *, 'alfsh_init rc =', alfsh_init(qh)
  rc = alfsh_describe()
  if (stage < 2) stop 'stage 1 done'

  allocate(A(N,N), B(N,N), C(N,N), Cr(N,N))
  if (stage < 2) stop
  call fill(A, 1); call fill(B, 2); C = 0.d0; Cr = 0.d0
  !$omp target enter data map(to: A, B, C, Cr)
  !$omp target data use_device_addr(A)
  print *, 'ptr kind of A on device (2 = device, 3 = shared, 1 = host, 0 = unknown): ', alfsh_ptr_kind(c_loc(A(1,1)))
  !$omp end target data

  if (stage < 3) stop 'stage 2 done'
  ! (a) dispatch directly on the mapped arrays
  !$omp dispatch
  call zgemm('N', 'N', N, N, N, one, A, N, B, N, zero, C, N)
  !$omp target update from(C)
  call zgemm('N', 'N', N, N, N, one, A, N, B, N, zero, Cr, N)
  err = maxval(abs(C - Cr)) / maxval(abs(Cr))
  print *, '(a) dispatch zgemm direct: rel max err vs host =', err
  ! was it done on the device? C changed on the host only via update; put garbage on the host first
  C = 99.d0
  !$omp target update from(C)
  print *, '    after update from: max|C-Cr| / max|Cr| =', maxval(abs(C - Cr)) / maxval(abs(Cr))

  if (stage < 4) stop 'stage 3 done'
  ! (b) dispatch inside target data use_device_addr
  !$omp target data use_device_addr(A, B, C)
  !$omp dispatch
  call zgemm('N', 'N', N, N, N, one, A, N, B, N, zero, C, N)
  !$omp end target data
  C = 99.d0
  !$omp target update from(C)
  print *, '(b) dispatch in use_device_addr: rel max err =', maxval(abs(C - Cr)) / maxval(abs(Cr))

  if (stage < 5) stop 'stage 4 done'
  ! (c) ordering hazards (t_order-like), dispatch only
  bsum = 0.d0
  do i = 1, N
    bsum = bsum + dble(B(i,i))
  end do
  worst = 0.d0
  do rep = 1, 24
    call set_diag(A, dble(rep))
    if (rep <= 3) print *, 'RAW rep', rep, ' set_diag done'
    !$omp dispatch
    call zgemm('N', 'N', N, N, N, one, A, N, B, N, zero, C, N)
    if (rep <= 3) print *, 'RAW rep', rep, ' zgemm done'
    accd(rep) = sumdiag(C)
    if (rep <= 3) print *, 'RAW rep', rep, ' sumdiag done'
    worst = max(worst, abs(accd(rep) - dble(rep)*bsum) / abs(dble(rep)*bsum))
  end do
  print *, '(c) RAW worst rel err =', worst
  if (stage < 6) stop 'stage 5 done'
  worst = 0.d0
  do rep = 1, 24
    call set_diag(A, 1.d0)
    !$omp dispatch
    call zgemm('N', 'N', N, N, N, one, A, N, B, N, zero, Cr, N)
    call set_diag(A, 2.d0)
    !$omp dispatch
    call zgemm('N', 'N', N, N, N, one, A, N, B, N, zero, C, N)
    worst = max(worst, abs(sumdiag(Cr) - bsum) / abs(bsum), abs(sumdiag(C) - 2*bsum) / abs(2*bsum))
  end do
  print *, '    WAR worst rel err =', worst
  if (stage < 7) stop 'stage 6 done'
  worst = 0.d0
  call set_diag(A, 3.d0)
  do rep = 1, 24
    !$omp dispatch
    call zgemm('N', 'N', N, N, N, one, A, N, B, N, zero, C, N)
    call fill_const(C, 7.d0)
    worst = max(worst, abs(sumall(C) - 7.d0*N*N) / (7.d0*N*N))
  end do
  print *, '    WAW worst rel err =', worst
  if (stage < 8) stop 'stage 7 done'
  !$omp target exit data map(delete: A, B, C, Cr)

  if (stage < 9) stop 'stage 8 done'
  ! (d) shim geqrf on the device, on the queue chosen above
  allocate(Aq(N,N), Aqr(N,N), tau(N), hw(1))
  call fill(Aq, 3); Aqr = Aq
  call zgeqrf(N, N, Aqr, N, tau, wq, -1, rc)
  deallocate(hw); allocate(hw(int(dble(wq(1)))))
  call zgeqrf(N, N, Aqr, N, tau, hw, size(hw), rc)
  !$omp target enter data map(to: Aq) map(alloc: tau)
  print *, 'geqrf_ws rc =', alfsh_zgeqrf_ws(int(N,8), int(N,8), int(N,8), ws), ' ws =', ws, ' host lwork =', size(hw)
  allocate(w(ws))
  !$omp target enter data map(alloc: w)
  call set_diag_dummy()
  !$omp target data use_device_addr(Aq, tau, w)
  rc = alfsh_zgeqrf(int(N,8), int(N,8), c_loc(Aq(1,1)), int(N,8), c_loc(tau(1)), c_loc(w(1)), ws)
  !$omp end target data
  print *, 'shim geqrf rc =', rc
  Aq = 0.d0
  !$omp target update from(Aq)
  print *, 'shim geqrf: |A - Aref|_max / |Aref|_max =', maxval(abs(Aq - Aqr)) / maxval(abs(Aqr))
  !$omp target exit data map(delete: Aq, tau, w)
contains
  subroutine set_diag_dummy()
  end subroutine
  subroutine fill(M, s)
    complex(8), intent(out) :: M(:,:)
    integer, intent(in) :: s
    integer :: ii, jj
    integer(8) :: r
    real(8) :: a, b
    r = 12345_8 + 7919_8*s
    do jj = 1, size(M,2)
      do ii = 1, size(M,1)
        r = mod(r * 1103515245_8 + 12345_8, 2147483648_8); a = dble(r) / 2147483648.d0 - 0.5d0
        r = mod(r * 1103515245_8 + 12345_8, 2147483648_8); b = dble(r) / 2147483648.d0 - 0.5d0
        M(ii,jj) = cmplx(a, b, kind=8)
      end do
    end do
  end subroutine
  subroutine set_diag(M, v)
    complex(8), intent(inout) :: M(N,N)
    real(8), intent(in) :: v
    integer :: ii, jj
    !$omp target teams distribute parallel do collapse(2)
    do jj = 1, N
      do ii = 1, N
        if (ii == jj) then
          M(ii,jj) = cmplx(v, 0.d0, kind=8)
        else
          M(ii,jj) = (0.d0, 0.d0)
        end if
      end do
    end do
  end subroutine
  subroutine fill_const(M, v)
    complex(8), intent(inout) :: M(N,N)
    real(8), intent(in) :: v
    integer :: ii, jj
    !$omp target teams distribute parallel do collapse(2)
    do jj = 1, N
      do ii = 1, N
        M(ii,jj) = cmplx(v, 0.d0, kind=8)
      end do
    end do
  end subroutine
  real(8) function sumdiag(M)
    complex(8), intent(in) :: M(N,N)
    integer :: ii
    real(8) :: s
    s = 0.d0
    !$omp target teams distribute parallel do reduction(+: s)
    do ii = 1, N
      s = s + dble(M(ii,ii))
    end do
    sumdiag = s
  end function
  real(8) function sumall(M)
    complex(8), intent(in) :: M(N,N)
    integer :: ii, jj
    real(8) :: s
    s = 0.d0
    !$omp target teams distribute parallel do collapse(2) reduction(+: s)
    do jj = 1, N
      do ii = 1, N
        s = s + dble(M(ii,jj))
      end do
    end do
    sumall = s
  end function
end program t302a_probe
