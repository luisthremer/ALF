! 302a development (S2): what does oneMKL zgetrf on the B580 do with an exactly singular matrix? (HOST ZGETRF tolerates
! info > 0 and la_getrf tolerates it too.) Not a ctest. All arrays mapped. Prints the outcome; must not crash.
program t302a_singular
  use iso_c_binding
  use alf_la_mod
  implicit none
  integer, parameter :: n = 8
  complex(8), allocatable, target :: A(:,:), B(:,:)
  integer(c_int), allocatable, target :: ip(:)
  integer :: i, j
  allocate(A(n,n), B(n,1), ip(n))
  do j = 1, n
    do i = 1, n
      A(i,j) = cmplx(i + j, i - j, 8)
    end do
  end do
  A(:,3) = 0                    ! an exactly zero column: singular U
  B = 1; ip = 0
  call la_init()
  !$omp target enter data map(to: A, B, ip)
  call la_getrf(n, n, A, n, ip)
  write(*,'(A)') 'singular: la_getrf returned'
  call la_sync()
  write(*,'(A)') 'singular: la_sync returned (info > 0 tolerated)'
  !$omp target exit data map(from: A, ip) 
  write(*,'(A,8I4)') 'singular: ipiv = ', ip
  write(*,'(A,2ES12.4)') 'singular: U(3,3) = ', A(3,3)
  call la_finalize()
  write(*,'(A)') 'singular: done'
end program t302a_singular
