! Timers_mod (Prog/Timers_mod.F90), compiled into this test with -DALF_TIMERS: one bin with two sweeps,
! nested push/pop/switch and a flush, then timer_bin_end writes timers.txt (directory of "info").
! timer_bin_end itself checks the sum rule (exclusive timers = inclusive sweep time) and stops the
! program if it fails; the test then checks the block: sweeps, counters, and that every timer is >= 0.
Program timers_test
  Use Timers_mod
  Use comm_mod, only: comm_world
#ifdef MPI
  Use mpi
#endif
  Implicit none
#ifdef MPI
  Integer :: ierr
#endif
  Integer :: u, st, k, nb, nsw, nprop, nacc, nflush, ntimer
  Character (len=256) :: line
  Character (len=32)  :: tag, name
  Real (Kind=Kind(0.d0)) :: x0, xmax, xmean

#ifdef MPI
  call MPI_Init(ierr)
#endif
  open(newunit=u, file="timers.txt", status="unknown")
  close(u, status="delete")
  call timer_bin_begin()
  do k = 1, 2
     call timer_sweep_begin()
     call timer_push(T_HOP)
     call spin()
     call timer_pop(T_HOP)
     call timer_push(T_VWRAP)
     call timer_switch(T_VWRAP, T_UPDATE)
     call spin()
     call timer_count_proposal(k == 1)
     timer_flush_reason = FR_CAPACITY
     call timer_flush_begin(3)
     call spin()
     call timer_flush_end()
     call timer_switch(T_UPDATE, T_VWRAP)
     call timer_pop(T_VWRAP)
     timer_flush_reason = FR_SLICE_END
     call timer_flush_skip()
     call timer_sweep_end()
  enddo
  call timer_bin_end(1, comm_world(), "info")

  nsw = -1; nprop = -1; nacc = -1; nflush = -1; ntimer = 0
  open(newunit=u, file="timers.txt", status="old", action="read")
  do
     read(u, '(A)', iostat=st) line
     if (st /= 0) exit
     if (line(1:1) == "#") cycle
     read(line, *) tag
     select case (trim(tag))
     case ("TIMER_SWEEPS")
        read(line, *) tag, nb, nsw
     case ("TIMER")
        read(line, *) tag, nb, name, x0, xmax, xmean
        if (x0 < 0.d0) call fail("negative timer " // trim(name))
        ntimer = ntimer + 1
     case ("COUNTER")
        read(line, *) tag, nb, name
        if (trim(name) == "proposals")      read(line, *) tag, nb, name, nprop
        if (trim(name) == "accepts")        read(line, *) tag, nb, name, nacc
        if (trim(name) == "flush_capacity") read(line, *) tag, nb, name, nflush
     end select
  enddo
  close(u)
  if (nsw /= 2 .or. nprop /= 2 .or. nacc /= 1 .or. nflush /= 2 .or. ntimer /= N_TIMERS + 2) then
     write(*,*) nsw, nprop, nacc, nflush, ntimer
     call fail("unexpected timers.txt block")
  endif
  write(*,*) 'PASS: timers.txt block consistent'
#ifdef MPI
  call MPI_Finalize(ierr)
#endif

contains

  subroutine spin()
    Integer :: i
    Real (Kind=Kind(0.d0)) :: a
    a = 0.d0
    do i = 1, 100000
       a = a + sqrt(dble(i))
    enddo
    if (a < 0.d0) write(*,*) a
  end subroutine spin

  subroutine fail(msg)
    Character (len=*), intent(in) :: msg
    write(*,*) 'FAIL: ', msg
    error stop 1
  end subroutine fail

end Program timers_test
