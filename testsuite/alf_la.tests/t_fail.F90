! alf_la unit test: an illegal argument must end the run with a nonzero exit code and the alf_la message.
! HOST: relies on xerbla RETURNING after the argument error (OpenBLAS, MKL), so that la_ungqr sees info < 0; a LAPACK whose
! xerbla stops the program would still terminate nonzero but without the alf_la message.
program t_fail
  use alf_la_mod
  use la_test_util
  implicit none
  complex(8) :: A(2,5), tau(5)
  call la_init()
  call fill_z(2, 5, A); tau = (0.d0, 0.d0)
  call dev_enter_z(A, 10); call dev_enter_z(tau, 5)
  call la_ungqr(2, 5, 5, A, 2, tau)      ! m < n: illegal
  call la_sync()
  write(*,*) 'ERROR: t_fail reached the end without failing'
end program t_fail
