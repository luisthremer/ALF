! micro timings of the pieces of one tiny la_geqrf call under ONEMKL (development only)
program micro
  use alf_la_mod
  use la_test_util, only: fill_z, dev_enter_z, dev_exit_z
  use omp_lib
  use iso_c_binding
  implicit none
  interface
    integer(c_int) function alfsh_zgeqrf(m, n, a, lda, tau, work, lwork) bind(C, name='alfsh_zgeqrf')
      import :: c_int, c_int64_t, c_ptr
      integer(c_int64_t), value :: m, n, lda, lwork
      type(c_ptr), value :: a, tau, work
    end function
    integer(c_int) function alfsh_zgeqrf_ws(m, n, lda, ws) bind(C, name='alfsh_zgeqrf_ws')
      import :: c_int, c_int64_t
      integer(c_int64_t), value :: m, n, lda
      integer(c_int64_t) :: ws
    end function
    integer(c_int) function alfsh_wait() bind(C, name='alfsh_wait')
      import :: c_int
    end function
  end interface
  integer, parameter :: NR = 2000
  complex(8), allocatable, target :: small(:,:), tau(:), w(:)
  integer :: i, rc, dev, lw
  integer(c_int64_t) :: ws
  real(8) :: t0
  integer :: np
  call la_init()
  allocate(small(8,8), tau(8)); call fill_z(8, 8, small); tau = 0
  lw = la_geqrf_lwork(8, 8, small, 8)
  allocate(w(lw)); w = 0
  call dev_enter_z(small, 64); call dev_enter_z(tau, 8); call dev_enter_z(w, lw)
  dev = omp_get_default_device()
  t0 = wt()
  do i = 1, NR; call la_geqrf(8, 8, small, 8, tau); end do
  call la_sync(); call tp('la_geqrf total per call')
  t0 = wt(); np = 0
  do i = 1, NR*6; if (omp_target_is_present(c_loc(small(1,1)), dev) /= 0) np = np + 1; end do
  call tp('6 x is_present per call (per 6)', 6)
  t0 = wt()
  do i = 1, NR
    !$omp target data use_device_addr(small, tau, w)
    np = np + 1
    !$omp end target data
  end do
  call tp('target data use_device_addr')
  t0 = wt()
  do i = 1, NR
    !$omp target data use_device_addr(small, tau, w)
    rc = alfsh_zgeqrf(8_c_int64_t, 8_c_int64_t, c_loc(small(1,1)), 8_c_int64_t, c_loc(tau(1)), c_loc(w(1)), int(lw,c_int64_t))
    !$omp end target data
  end do
  call tp('data + shim zgeqrf')
  t0 = wt()
  do i = 1, NR; rc = alfsh_zgeqrf_ws(8_c_int64_t, 8_c_int64_t, 8_c_int64_t, ws); end do
  call tp('geqrf_scratchpad_size')
  t0 = wt(); np = 0
  do i = 1, NR; np = np + la_geqrf_lwork(8, 8, small, 8); end do
  call tp('la_geqrf_lwork'); print *, np
  call la_finalize()
contains
  real(8) function wt()
    integer(8) :: c, r
    call system_clock(c, r); wt = dble(c)/dble(r)
  end function
  subroutine tp(s, d)
    character(*), intent(in) :: s
    integer, intent(in), optional :: d
    real(8) :: dd
    dd = NR; if (present(d)) dd = dd*d
    write(*,'(A40,F10.1,A)') s, (wt() - t0)/dd*1.d6, ' us'
  end subroutine
end program micro
