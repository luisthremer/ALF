! alf_exit (Libraries/Modules/runtime_error_mod.F90): the process exit status equals the argument for
! 0 and 2, and is 1 for any other value. The code is passed as the first command-line argument; the
! exit status is checked by expect_exit.cmake (ctest tests 37-alf-exit-<code>).
Program alf_exit_test
  Use runtime_error_mod, only: alf_exit
  Implicit none
  Character (len=16) :: arg
  Integer :: code, st
  call get_command_argument(1, arg)
  read(arg, *, iostat=st) code
  if (st /= 0) then
     write(*,*) 'usage: 37-alf-exit <code>'
     error stop 99
  endif
  write(*,*) 'calling alf_exit(', code, ')'
  call alf_exit(code)
  write(*,*) 'FAIL: alf_exit returned'
  error stop 98
end Program alf_exit_test
