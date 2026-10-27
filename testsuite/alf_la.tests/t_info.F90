! alf_la unit test: a solver failure must terminate the run, also without an explicit la_sync (info-slot recycling).
!   t_info recycle : a failed info slot (test hook la_debug_poison_last_info) followed by more than LA_NSLOT solver calls
!                    and NO la_sync: the run must terminate with "geqrf failed" when the slots are read back.
!   t_info sync    : the same failure must be reported by la_sync.
! Registered with PASS_REGULAR_EXPRESSION "alf_la: .*geqrf failed"; reaching the end prints a message that cannot match.
program t_info
  use alf_la_mod
  use la_test_util
  implicit none
  complex(8) :: small(8,8), tau(8)
  character(len=16) :: arg
  integer :: i
  call get_command_argument(1, arg)
  call la_init()
  call fill_z(8, 8, small); tau = (0.d0, 0.d0)
  call dev_enter_z(small, 64); call dev_enter_z(tau, 8)
  call la_geqrf(8, 8, small, 8, tau)
  call la_debug_poison_last_info()
  if (trim(arg) == 'recycle') then
    do i = 1, LA_NSLOT + 50
      call la_geqrf(8, 8, small, 8, tau)
    end do
    write(*,*) 'NOTE: LA_NSLOT+50 calls done without termination'
    call la_sync()
  else
    call la_sync()
  end if
  write(*,*) 'ERROR: t_info reached the end without failing'
end program t_info
