! Trace_updates_mod (Prog/Trace_updates_mod.F90): Trace_update_write appends one line per call to
! trace_updates.txt with nt, n_op, Re(R), Im(R), Weight, ranf, toggle (0/1). The test removes the file,
! writes two lines and reads them back; values exact in binary, so the comparison is exact.
Program trace_updates_test
  Use Trace_updates_mod, only: Trace_update_write
  Implicit none
  Integer :: u, st, nt, n_op, itog
  Real (Kind=Kind(0.d0)) :: re, im, w, r
  open(newunit=u, file="trace_updates.txt", status="unknown")
  close(u, status="delete")
  call Trace_update_write(3, 7, cmplx(0.5d0, -0.25d0, kind(0.d0)), 0.75d0, 0.125d0, .true.)
  call Trace_update_write(4, 1, cmplx(2.d0, 0.d0, kind(0.d0)), 0.0625d0, 0.5d0, .false.)
  flush(97)
  open(newunit=u, file="trace_updates.txt", status="old", action="read")
  read(u, *, iostat=st) nt, n_op, re, im, w, r, itog
  if (st /= 0 .or. nt /= 3 .or. n_op /= 7 .or. re /= 0.5d0 .or. im /= -0.25d0 .or. w /= 0.75d0 &
       & .or. r /= 0.125d0 .or. itog /= 1) then
     write(*,*) 'FAIL: line 1'
     error stop 1
  endif
  read(u, *, iostat=st) nt, n_op, re, im, w, r, itog
  if (st /= 0 .or. nt /= 4 .or. n_op /= 1 .or. re /= 2.d0 .or. im /= 0.d0 .or. w /= 0.0625d0 &
       & .or. r /= 0.5d0 .or. itog /= 0) then
     write(*,*) 'FAIL: line 2'
     error stop 1
  endif
  close(u)
  write(*,*) 'PASS: trace lines written and read back'
end Program trace_updates_test
