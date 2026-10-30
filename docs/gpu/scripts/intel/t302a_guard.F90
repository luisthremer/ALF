! 302a development: check of the host-side presence guard of the ONEMKL backend. Not a ctest.
!   t302a_guard present : omp_target_is_present alone (mapped vs unmapped array); no alf_la, no dispatch, no shim call
!   t302a_guard pos     : la_gemm with all arguments mapped (must work)
!   t302a_guard neg A|B|C : la_gemm with that argument NOT mapped: must stop with "argument X is not present on the device"
!                         before anything is queued (exit code /= 0, no omptarget message, no crash)
!   t302a_guard neg_ipiv  : la_getrf with A mapped and ipiv NOT mapped (same expectation, argument ipiv)
!   t302a_guard neg_ipvt  : la_qrpp with A and tau mapped and ipvt NOT mapped (same expectation, argument ipvt)
program t302a_guard
  use iso_c_binding
  use omp_lib
  use alf_la_mod
  implicit none
  character(len=16) :: mode, which
  complex(8), allocatable, target :: A(:,:), B(:,:), C(:,:), tau(:)
  integer(c_int), allocatable, target :: ip(:)
  integer, parameter :: n = 16
  integer :: i, j
  call get_command_argument(1, mode)
  which = ''
  if (command_argument_count() > 1) call get_command_argument(2, which)
  allocate(A(n,n), B(n,n), C(n,n))
  do j = 1, n
    do i = 1, n
      A(i,j) = cmplx(i, j, 8); B(i,j) = cmplx(j, -i, 8); C(i,j) = 0
    end do
  end do
  if (trim(mode) == 'present') then
    !$omp target enter data map(to: A)
    write(*,'(A,I0)') 'present(mapped A)   = ', omp_target_is_present(c_loc(A(1,1)), omp_get_default_device())
    write(*,'(A,I0)') 'present(mapped A(n,n)) = ', omp_target_is_present(c_loc(A(n,n)), omp_get_default_device())
    write(*,'(A,I0)') 'present(unmapped B) = ', omp_target_is_present(c_loc(B(1,1)), omp_get_default_device())
    !$omp target exit data map(delete: A)
    write(*,'(A,I0)') 'present(A after delete) = ', omp_target_is_present(c_loc(A(1,1)), omp_get_default_device())
    stop 0
  end if
  call la_init()
  if (trim(mode) == 'neg_ipiv' .or. trim(mode) == 'neg_ipvt') then
    allocate(tau(n), ip(n)); tau = 0; ip = 0
    !$omp target enter data map(to: A, tau)          ! ip stays unmapped on purpose
    if (trim(mode) == 'neg_ipiv') then
      call la_getrf(n, n, A, n, ip)
    else
      call la_qrpp(n, n, A, n, ip, tau)
    end if
    call la_sync()
    write(*,'(A)') 'FAILED guard: the call returned although ipiv/ipvt is not mapped'
    call la_finalize()
    stop 0
  end if
  if (trim(mode) == 'pos' .or. trim(which) /= 'A') then
    !$omp target enter data map(to: A)
  end if
  if (trim(mode) == 'pos' .or. trim(which) /= 'B') then
    !$omp target enter data map(to: B)
  end if
  if (trim(mode) == 'pos' .or. trim(which) /= 'C') then
    !$omp target enter data map(to: C)
  end if
  call la_gemm(0, 0, n, n, n, (1.d0, 0.d0), A, n, B, n, (0.d0, 0.d0), C, n)
  call la_sync()
  if (trim(mode) == 'pos') then
    !$omp target exit data map(from: C)
    write(*,'(A,2F10.1)') 'pos: C(1,1) = ', C(1,1)
    write(*,'(A,2F10.1)') 'pos: reference sum_k A(1,k)B(k,1) = ', sum(A(1,:)*B(:,1))
  else
    write(*,'(A)') 'FAILED guard: la_gemm returned although an argument is not mapped'
  end if
  call la_finalize()
end program t302a_guard
