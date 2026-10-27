! alf_la unit test: la_geqrf, la_ungqr, la_ormqr (ticket 202).
! HOST: factors, Q and the ormqr results bit-identical to the direct LAPACK calls.
! CUDA: invariants, not raw factors (cuSOLVER may differ in rounding / block structure):
!   |R_jj| sequence vs LAPACK (relative), ||A - Q R|| / ||A||, ||Q^H Q - I||_F  (Q formed by host ZUNGQR from the
!   device factors), device ZUNGQR vs host ZUNGQR on the same factors, device ZUNMQR vs host ZUNMQR on the same factors.
program t_qr
  use alf_la_mod
  use la_test_util
  implicit none
  external :: zgeqrf, zungqr, zunmqr, zgemm
  integer :: in, cs
  call test_begin('t_qr')
  call la_init()
  do in = 1, size(sizes)
    do cs = 1, 3
      call run(sizes(in), cs)
    end do
  end do
  call la_finalize()
  call test_end()
contains
  character function ch(op)
    integer, intent(in) :: op
    select case (op)
    case (LA_N); ch = 'N'
    case (LA_T); ch = 'T'
    case default; ch = 'C'
    end select
  end function

  subroutine run(N, cs)
    integer, intent(in) :: N, cs
    integer :: m, n2, oi, oj, lda, info, lw, i, j, nrc, ldc, side_i, tr_i
    complex(8), allocatable :: A0(:,:), A(:,:), Aref(:,:), tau(:), tauref(:), Q(:,:), Qref(:,:), R(:,:), QR(:,:), G(:,:), &
         C(:,:), Cref(:,:), w(:)
    complex(8) :: wq(1)
    real(8) :: e, rjj, dmax
    character(len=64) :: lab
    character(len=1) :: sd
    n2 = N; m = N; oi = 1; oj = 1; lda = N
    if (cs >= 2) then
      m = N + 9; lda = m + 5
    end if
    if (cs == 3) then
      oi = 3; oj = 2; lda = m + oi - 1 + 4
    end if
    allocate(A0(lda, n2+oj), A(lda, n2+oj), Aref(lda, n2+oj), tau(n2), tauref(n2), Q(lda, n2+oj), Qref(lda, n2+oj))
    call fill_z(lda, n2+oj, A0)
    A = A0; Aref = A0
    ! ---- geqrf
    call zgeqrf(m, n2, Aref(oi,oj), lda, tauref, wq, -1, info)
    lw = int(dble(wq(1))); allocate(w(lw))
    call zgeqrf(m, n2, Aref(oi,oj), lda, tauref, w, lw, info)
    call dev_enter_z(A, size(A)); call dev_enter_z(tau, size(tau))
    call la_geqrf(m, n2, A(oi,oj), lda, tau)
    call la_sync()
    call dev_from_z(A, size(A)); call dev_from_z(tau, size(tau))
    write(lab,'(A,I0,A,I0)') 'geqrf N=', N, ' case', cs
    if (is_host()) then
      call check(trim(lab)//' A', maxdiff(lda, n2+oj, A, Aref), 0.d0)
      call check(trim(lab)//' tau', maxabs_diff(n2, tau, tauref), 0.d0)
    end if
    call check(trim(lab)//' outside', outside_diff(lda, n2+oj, A, A0, oi, oj, m, n2), 0.d0)
    dmax = 0.d0
    do j = 1, n2
      rjj = abs(Aref(oi+j-1, oj+j-1))
      dmax = fmax(dmax, abs(abs(A(oi+j-1, oj+j-1)) - rjj) / rjj)
    end do
    call check(trim(lab)//' |Rjj|', dmax, TOL_CUDA)
    ! invariants from the device factors: Q by host ZUNGQR
    Q = A
    call zungqr(m, n2, n2, Q(oi,oj), lda, tau, wq, -1, info)
    lw = int(dble(wq(1))); deallocate(w); allocate(w(lw))
    call zungqr(m, n2, n2, Q(oi,oj), lda, tau, w, lw, info)
    allocate(R(n2, n2), QR(m, n2), G(n2, n2))
    R = (0.d0, 0.d0)
    do j = 1, n2
      do i = 1, j
        R(i,j) = A(oi+i-1, oj+j-1)
      end do
    end do
    call zgemm('N', 'N', m, n2, n2, (1.d0, 0.d0), Q(oi,oj), lda, R, n2, (0.d0, 0.d0), QR, m)
    call check_inv(trim(lab)//' ||A-QR||/||A||', relfro(m, n2, QR, m, A0(oi,oj), lda), TOL_CUDA)
    call zgemm('C', 'N', n2, n2, m, (1.d0, 0.d0), Q(oi,oj), lda, Q(oi,oj), lda, (0.d0, 0.d0), G, n2)
    do i = 1, n2
      G(i,i) = G(i,i) - (1.d0, 0.d0)
    end do
    call check_inv(trim(lab)//' ||QhQ-I||', nrmf(n2, n2, G, n2), TOL_CUDA)
    deallocate(R, QR, G)
    ! ---- ungqr on the device, from the device factors (Qref: host ZUNGQR on the same factors = Q above)
    Qref = Q
    Q = A
    call dev_enter_z(Q, size(Q))
    call la_ungqr(m, n2, n2, Q(oi,oj), lda, tau)
    call la_sync()
    call dev_from_z(Q, size(Q))
    call dev_exit_z(Q, size(Q))
    write(lab,'(A,I0,A,I0)') 'ungqr N=', N, ' case', cs
    call check(trim(lab), relfro(m, n2, Q(oi,oj), lda, Qref(oi,oj), lda), TOL_CUDA)
    call check(trim(lab)//' outside', outside_diff(lda, n2+oj, Q, Qref, oi, oj, m, n2), 0.d0)
    ! ---- ormqr: side L (C is m x nrc) and R (C is nrc x m), trans N and C, device factors
    nrc = N + 2
    do side_i = 1, 2
      do tr_i = 1, 2
        sd = merge('L', 'R', side_i == 1)
        ldc = merge(m, nrc, side_i == 1) + 3
        j = merge(nrc, m, side_i == 1)
        allocate(C(ldc, j), Cref(ldc, j))
        call fill_z(ldc, j, C)
        Cref = C
        i = merge(LA_N, LA_C, tr_i == 1)
        if (side_i == 1) then
          call zunmqr(sd, ch(i), m, nrc, n2, A(oi,oj), lda, tau, Cref, ldc, wq, -1, info)
        else
          call zunmqr(sd, ch(i), nrc, m, n2, A(oi,oj), lda, tau, Cref, ldc, wq, -1, info)
        end if
        lw = int(dble(wq(1))); deallocate(w); allocate(w(lw))
        if (side_i == 1) then
          call zunmqr(sd, ch(i), m, nrc, n2, A(oi,oj), lda, tau, Cref, ldc, w, lw, info)
        else
          call zunmqr(sd, ch(i), nrc, m, n2, A(oi,oj), lda, tau, Cref, ldc, w, lw, info)
        end if
        call dev_enter_z(C, size(C))
        if (side_i == 1) then
          call la_ormqr(sd, i, m, nrc, n2, A(oi,oj), lda, tau, C, ldc)
        else
          call la_ormqr(sd, i, nrc, m, n2, A(oi,oj), lda, tau, C, ldc)
        end if
        call la_sync()
        call dev_from_z(C, size(C))
        call dev_exit_z(C, size(C))
        write(lab,'(A,I0,A,I0,A,A,A)') 'ormqr N=', N, ' case', cs, ' ', sd // ch(i), ''
        if (side_i == 1) then
          e = relfro(m, nrc, C, ldc, Cref, ldc)
          call check(trim(lab)//' outside', maxdiff_rows(ldc, size(C,2), m+1, ldc, C, Cref), 0.d0)
        else
          e = relfro(nrc, m, C, ldc, Cref, ldc)
          call check(trim(lab)//' outside', maxdiff_rows(ldc, size(C,2), nrc+1, ldc, C, Cref), 0.d0)
        end if
        call check(trim(lab), e, TOL_CUDA)
        deallocate(C, Cref)
      end do
    end do
    call dev_exit_z(A, size(A)); call dev_exit_z(tau, size(tau))
  end subroutine
end program t_qr
