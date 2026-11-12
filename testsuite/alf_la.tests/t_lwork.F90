! alf_la unit test: workspace queries la_geqrf_lwork, la_ungqr_lwork, la_ormqr_lwork, la_getrf_lwork and
! la_reserve_work (not called directly by the other alf_la tests).
! HOST: each *_lwork equals the LAPACK work query (lwork = -1) exactly; la_getrf_lwork = 0 (ZGETRF has no work array);
!       after la_reserve_work (larger and smaller than needed) la_geqrf / la_ungqr / la_ormqr stay bit-identical to
!       the direct LAPACK calls with the queried lwork.
! Other backends: only lw >= 0 (the device buffer sizes are the library's own and are not compared).
program t_lwork
  use alf_la_mod
  use la_test_util
  implicit none
  external :: zgeqrf, zungqr, zunmqr
  integer :: in, cs
  call test_begin('t_lwork')
  call la_init()
  do in = 1, size(sizes)
    do cs = 1, 2
      call run(sizes(in), cs)
    end do
  end do
  call la_finalize()
  call test_end()
contains

  subroutine run(N, cs)
    integer, intent(in) :: N, cs
    integer :: m, lda, ldc, info, lw, lq, i_side, i_tr, nr, ncol, kk
    complex(8), allocatable :: A(:,:), Aref(:,:), tau(:), tauref(:), C(:,:), Cref(:,:), w(:)
    complex(8) :: wq(1)
    character(len=1) :: sd
    character(len=64) :: lab
    integer, parameter :: trs(2) = [LA_N, LA_C]
    character(len=1), parameter :: trc(2) = ['N', 'C']

    m = N; lda = N
    if (cs == 2) then
      m = N + 5; lda = m + 3
    end if
    allocate(A(lda,N), Aref(lda,N), tau(N), tauref(N))
    call fill_z(lda, N, A)
    Aref = A

    ! ---- geqrf: query, then factorization after la_reserve_work (larger, then smaller than needed)
    call zgeqrf(m, N, Aref, lda, tauref, wq, -1, info)
    lq = int(dble(wq(1)))
    call dev_enter_z(A, size(A)); call dev_enter_z(tau, size(tau))
    lw = la_geqrf_lwork(m, N, A, lda)
    write(lab,'(A,I0,A,I0)') 'geqrf_lwork N=', N, ' case', cs
    call cmp_lw(lab, lw, lq)
    allocate(w(max(lq,1)))
    call zgeqrf(m, N, Aref, lda, tauref, w, lq, info)
    deallocate(w)
    call la_reserve_work(4*lq + 1000)
    call la_geqrf(m, N, A, lda, tau)
    call la_sync()
    call dev_from_z(A, size(A)); call dev_from_z(tau, size(tau))
    if (is_host()) then
      call check(trim(lab)//' geqrf after reserve A', maxdiff(lda, N, A, Aref), 0.d0)
      call check(trim(lab)//' geqrf after reserve tau', maxabs_diff(N, tau, tauref), 0.d0)
    end if

    ! ---- ormqr: query for every side / trans, then apply after a small la_reserve_work
    do i_side = 1, 2
      sd = 'L'; if (i_side == 2) sd = 'R'
      do i_tr = 1, 2
        if (sd == 'L') then
          nr = m; ncol = 3; kk = N
        else
          nr = 3; ncol = m; kk = N
        end if
        ldc = nr + 2
        allocate(C(ldc, ncol), Cref(ldc, ncol))
        call fill_z(ldc, ncol, C)
        Cref = C
        call zunmqr(sd, trc(i_tr), nr, ncol, kk, Aref, lda, tauref, Cref, ldc, wq, -1, info)
        lq = int(dble(wq(1)))
        call dev_enter_z(C, size(C))
        lw = la_ormqr_lwork(sd, trs(i_tr), nr, ncol, kk, A, lda, tau, C, ldc)
        write(lab,'(A,A,A,A,I0,A,I0)') 'ormqr_lwork ', sd, trc(i_tr), ' N=', N, ' case', cs
        call cmp_lw(lab, lw, lq)
        allocate(w(max(lq,1)))
        call zunmqr(sd, trc(i_tr), nr, ncol, kk, Aref, lda, tauref, Cref, ldc, w, lq, info)
        deallocate(w)
        call la_reserve_work(1)
        call la_ormqr(sd, trs(i_tr), nr, ncol, kk, A, lda, tau, C, ldc)
        call la_sync()
        call dev_exit_z(C, size(C))
        if (is_host()) call check(trim(lab)//' ormqr after reserve', maxdiff(ldc, ncol, C, Cref), 0.d0)
        deallocate(C, Cref)
      end do
    end do

    ! ---- ungqr: query, then Q
    call zungqr(m, N, N, Aref, lda, tauref, wq, -1, info)
    lq = int(dble(wq(1)))
    lw = la_ungqr_lwork(m, N, N, A, lda, tau)
    write(lab,'(A,I0,A,I0)') 'ungqr_lwork N=', N, ' case', cs
    call cmp_lw(lab, lw, lq)
    allocate(w(max(lq,1)))
    call zungqr(m, N, N, Aref, lda, tauref, w, lq, info)
    deallocate(w)
    call la_ungqr(m, N, N, A, lda, tau)
    call la_sync()
    call dev_from_z(A, size(A))
    if (is_host()) call check(trim(lab)//' ungqr Q', maxdiff(lda, N, A, Aref), 0.d0)

    ! ---- getrf: no work array with LAPACK
    lw = la_getrf_lwork(m, N, A, lda)
    write(lab,'(A,I0,A,I0)') 'getrf_lwork N=', N, ' case', cs
    call cmp_lw(lab, lw, 0)

    call dev_exit_z(A, size(A)); call dev_exit_z(tau, size(tau))
    deallocate(A, Aref, tau, tauref)
  end subroutine run

  subroutine cmp_lw(lab, lw, lq)
    character(len=*), intent(in) :: lab
    integer, intent(in) :: lw, lq
    if (is_host()) then
      call check(trim(lab)//' = LAPACK query', dble(abs(lw - lq)), 0.d0)
    else
      call check(trim(lab)//' >= 0', merge(0.d0, 1.d0, lw >= 0), 0.d0)
    end if
  end subroutine cmp_lw

end program t_lwork
