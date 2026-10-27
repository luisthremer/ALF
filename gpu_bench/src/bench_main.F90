program bench_main
  use bench_la
  use bench_util
  use bench_sec
  implicit none
  integer :: Ns(16), nN, i, k, nc, seed(64), ns_seed
  character(len=256) :: arg, secs
  integer :: stat
  real(8) :: t0
  nN = 6; Ns(1:6) = [64, 144, 256, 400, 576, 1024]
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    nc = 1
    do i = 1, len_trim(arg)
      if (arg(i:i) == ',') nc = nc + 1
    end do
    nN = nc
    read(arg,*) Ns(1:nN)
  end if
  secs = 'cb,dense,opv,qr,upd,batch'
  call get_environment_variable('BENCH_SECTIONS', secs, status=stat)
  if (stat /= 0) secs = 'cb,dense,opv,qr,upd,batch'
  call random_seed(size=ns_seed)
  seed = 0
  do i = 1, ns_seed
    seed(i) = 12345 + 7*i
  end do
  call random_seed(put=seed(1:ns_seed))
  call la_init()
  write(*,'(a,a)') 'backend: ', trim(la_backend_name())
  write(*,'(a,a)') 'sections: ', trim(secs)
  t0 = now()
  do k = 1, nN
    write(*,'(/,a,i0)') '--- N = ', Ns(k)
    if (index(secs,'cb') > 0)    call sec_cb(Ns(k))
    if (index(secs,'dense') > 0) call sec_dense(Ns(k))
    if (index(secs,'opv') > 0)   call sec_opv(Ns(k))
    if (index(secs,'qr') > 0)    call sec_qr(Ns(k))
    if (index(secs,'upd') > 0)   call sec_upd(Ns(k))
    if (index(secs,'batch') > 0) call sec_batch(Ns(k))
  end do
  write(*,'(/,a,f10.1,a)') 'total wall time of run: ', now() - t0, ' s'
  write(*,'(a,i0,a,i0,a)') 'checks: ', nchk, ' run, ', nfail, ' FAILED'
  call print_tables(6)
  call la_finalize()
  if (nfail > 0) stop 1
end program
