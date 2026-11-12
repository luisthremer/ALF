include "mkl_blas_omp_offload.f90"
! Gate A' (Phase 3b): OpenMP sync baseline vs SYCL submission on the interop queue (same chain as asyncbench.F90).
program syclbench
  use omp_lib, only: omp_interop_kind, omp_interop_none, omp_get_default_device, omp_ifr_sycl, omp_get_wtime, omp_set_default_device
  use iso_c_binding
  use onemkl_blas_omp_offload_lp64, only: zgemm
  implicit none
  integer, parameter :: n = 1024, k = 32, NCH = 32, NGEMM_EVERY = 8
  integer(c_int), parameter :: omp_ipr_fr_id = -1, omp_ipr_targetsync = -9
  interface
    integer(c_intptr_t) function omp_get_interop_int(interop, property_id, ret_code) bind(C, name='omp_get_interop_int')
      import :: c_intptr_t, c_int
      integer(c_intptr_t), value :: interop
      integer(c_int), value :: property_id
      integer(c_int), intent(out) :: ret_code
    end function
    type(c_ptr) function omp_get_interop_ptr(interop, property_id, ret_code) bind(C, name='omp_get_interop_ptr')
      import :: c_intptr_t, c_int, c_ptr
      integer(c_intptr_t), value :: interop
      integer(c_int), value :: property_id
      integer(c_int), intent(out) :: ret_code
    end function
    integer(c_int) function sb_init(q) bind(C, name='sb_init')
      import :: c_ptr, c_int
      type(c_ptr), value :: q
    end function
    subroutine sb_set(v, a, b, c, nn, kk) bind(C, name='sb_set')
      import :: c_ptr, c_int
      type(c_ptr), value :: v, a, b, c
      integer(c_int), value :: nn, kk
    end subroutine
    subroutine sb_kern(c, d) bind(C, name='sb_kern')
      import :: c_double
      real(c_double), value :: c, d
    end subroutine
    subroutine sb_gemm() bind(C, name='sb_gemm')
    end subroutine
    subroutine sb_copy() bind(C, name='sb_copy')
    end subroutine
    subroutine sb_wait() bind(C, name='sb_wait')
    end subroutine
  end interface
  integer :: nrep, nhw, i, j, rc
  real(8) :: t0, t1, t2, hw, chk, cc, dd
  complex(8), allocatable, target :: Am(:,:), Bm(:,:), Cm(:,:)
  real(8), allocatable, target :: v(:)
  complex(8), parameter :: one = (1.d0,0.d0), zero = (0.d0,0.d0)
  integer(omp_interop_kind) :: iop
  type(c_ptr) :: qh
  integer(c_intptr_t) :: fr
  integer :: nops, inord
  character(len=32) :: arg
  nrep = 300; nhw = 0
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read(arg,*) nrep; end if
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read(arg,*) nhw; end if
  call omp_set_default_device(0)
  allocate(Am(n,k), Bm(k,k), Cm(n,k), v(n))
  call reset()
  !$omp target enter data map(to: Am, Bm, Cm, v)
  iop = omp_interop_none
  !$omp interop init(prefer_type("sycl"), targetsync: iop) device(omp_get_default_device())
  fr = omp_get_interop_int(iop, omp_ipr_fr_id, rc)
  qh = c_null_ptr
  if (rc == 0 .and. fr == omp_ifr_sycl) qh = omp_get_interop_ptr(iop, omp_ipr_targetsync, rc)
  if (.not. c_associated(qh)) stop 'no sycl queue'
  inord = sb_init(qh)
  write(*,*) 'sycl queue ok'; flush(6)
  !$omp target data use_device_addr(v, Am, Bm, Cm)
  call sb_set(c_loc(v), c_loc(Am), c_loc(Bm), c_loc(Cm), n, k)
  !$omp end target data
  nops = nrep*(NCH + 2*(NCH/NGEMM_EVERY))
  write(*,'(A,I0,A,I0,A,I0)') 'in-order queue=', inord, '  nrep=', nrep, '  hostwork=', nhw
  ! 1: OpenMP sync baseline (today's path)
  do i = 1, 3
    call fresh()
    t0 = omp_get_wtime()
    do j = 1, nrep; call chain_omp(); end do
    t1 = omp_get_wtime()
    call report('OMP sync (target + dispatch)', i, t0, t1, t1)
  end do
  ! 2: SYCL, wait after every op
  do i = 1, 3
    call fresh()
    t0 = omp_get_wtime()
    do j = 1, nrep; call chain_sycl(.true.); end do
    t1 = omp_get_wtime()
    call report('SYCL wait after every op', i, t0, t1, t1)
  end do
  ! 3: SYCL, submit all, one wait at the end
  do i = 1, 3
    call fresh()
    t0 = omp_get_wtime()
    do j = 1, nrep; call chain_sycl(.false.); end do
    t2 = omp_get_wtime()
    call sb_wait()
    t1 = omp_get_wtime()
    call report('SYCL async, 1 wait at end', i, t0, t1, t2)
  end do
  !$omp interop destroy(iop)
contains
  subroutine reset()
    integer :: i, j
    do j = 1, n; v(j) = 1.d0 + 1.d-3*j; end do
    do j = 1, k
      do i = 1, n; Am(i,j) = cmplx(1.d0/(1+mod(i+j,7)), 1.d-2*j, 8); end do
    end do
    Bm = zero
    do j = 1, k; Bm(j,j) = cmplx(0.5d0, 0.5d0, 8); if (j < k) Bm(j+1,j) = cmplx(0.1d0, -0.05d0, 8); end do
    Cm = zero
  end subroutine
  subroutine fresh()
    call reset()
    !$omp target update to(v, Am, Bm, Cm)
    hw = 0.d0
  end subroutine
  subroutine report(lab, it, ta, tb, tissue)
    character(*), intent(in) :: lab
    integer, intent(in) :: it
    real(8), intent(in) :: ta, tb, tissue
    !$omp target update from(v, Am)
    chk = sum(v) + sum(real(Am)) + sum(aimag(Am))
    write(*,'(A,I1,A,F8.2,A,F8.2,A,ES24.16)') lab//' #', it, ': total ', 1.d6*(tb-ta)/nops, ' us/op  issue-only ', &
         1.d6*(tissue-ta)/nops, ' us/op  chk ', chk
  end subroutine
  subroutine hostwork()
    integer :: ih
    do ih = 1, nhw
      hw = hw*0.9999999d0 + 1.d-9*sin(hw+ih)
    end do
  end subroutine
  subroutine chain_omp()
    integer :: kk
    do kk = 1, NCH
      cc = 0.999d0 - 1.d-4*kk; dd = 1.d-3*kk
      call kern(cc, dd)
      call hostwork()
      if (mod(kk, NGEMM_EVERY) == 0) then
        !$omp dispatch
        call zgemm('N','N',n,k,k,one,Am,n,Bm,k,zero,Cm,n)
        call kcopy()
      end if
    end do
  end subroutine
  subroutine chain_sycl(wait_each)
    integer :: kk
    logical, intent(in) :: wait_each
    do kk = 1, NCH
      cc = 0.999d0 - 1.d-4*kk; dd = 1.d-3*kk
      call sb_kern(cc, dd); if (wait_each) call sb_wait()
      call hostwork()
      if (mod(kk, NGEMM_EVERY) == 0) then
        call sb_gemm(); if (wait_each) call sb_wait()
        call sb_copy(); if (wait_each) call sb_wait()
      end if
    end do
  end subroutine
  subroutine kern(c, d)
    real(8), intent(in) :: c, d
    integer :: jj
    !$omp target teams distribute parallel do map(present, alloc: v)
    do jj = 1, n
      v(jj) = v(jj)*c + d*(1.d0 + 1.d-4*jj)
    end do
  end subroutine
  subroutine kcopy()
    integer :: jj, ii
    !$omp target teams distribute parallel do collapse(2) map(present, alloc: Am, Cm, v)
    do jj = 1, k
      do ii = 1, n
        Am(ii,jj) = Cm(ii,jj)*cmplx(0.9d0, 0.d0, 8) + cmplx(1.d-3*v(ii), 0.d0, 8)
      end do
    end do
  end subroutine
end program
