! alf_la unit test: la_qrpp (column-norm pre-pivot + ZGEQRF, ticket 202; semantics of QR_PREPIVOT variant 1).
! Reference (in this program): DZNRM2 per column, stable descending insertion sort (the code of
! Prog/cgr1_mod.F90 QR_PREPIVOT variant 1), column permutation, ZGEQRF with the optimal LWORK.
! HOST: A, tau, ipvt bit-identical to the reference. CUDA: ipvt identical to the reference (the matrices have
! columns with random scales 10^(-3 p/n), p a random permutation of 0..n-1, so the column norms are well separated
! and the device norm kernel and DZNRM2 give the same order), and invariants: |R_jj| relative, ||A P - Q R||/||A||,
! ||Q^H Q - I||_F with Q from host ZUNGQR applied to the device factors.
program t_qrpp
  use alf_la_mod
  use la_test_util
  implicit none
  external :: zgeqrf, zungqr, zgemm, zgeqp3
  real(8), external :: dznrm2
  integer :: in, cs
  call test_begin('t_qrpp')
  call la_init()
  do in = 1, size(sizes)
    do cs = 1, 3
      call run(sizes(in), cs)
    end do
  end do
  call la_finalize()
  call test_end()
contains
  subroutine run(N, cs)
    integer, intent(in) :: N, cs
    integer :: m, n2, oi, oj, lda, info, lw, i, j, p, k
    complex(8), allocatable :: A0(:,:), A(:,:), Aref(:,:), tau(:), tauref(:), Q(:,:), R(:,:), QR(:,:), AP(:,:), G(:,:), w(:)
    integer(4), allocatable :: ip(:), ipref(:), perm(:)
    real(8), allocatable :: nrm(:)
    complex(8) :: wq(1)
    complex(8), allocatable :: A3(:,:), tau3(:)
    real(8), allocatable :: rw(:)
    integer(4), allocatable :: jp(:)
    real(8) :: dmax, rjj, s
    integer(4) :: kk
    character(len=64) :: lab
    n2 = N; m = N; oi = 1; oj = 1; lda = N
    if (cs >= 2) then
      m = N + 9; lda = m + 5
    end if
    if (cs == 3) then
      oi = 3; oj = 2; lda = m + oi - 1 + 4
    end if
    allocate(A0(lda, n2+oj), A(lda, n2+oj), Aref(lda, n2+oj), tau(n2), tauref(n2), Q(lda, n2+oj), ip(n2), ipref(n2), &
         perm(n2), nrm(n2))
    call fill_z(lda, n2+oj, A0)
    do j = 1, n2                                   ! random permutation of 0..n-1 for the column scales
      perm(j) = j - 1
    end do
    do j = n2, 2, -1
      p = min(max(1 + int((rnd() + 0.5d0) * dble(j)), 1), j)
      kk = perm(j); perm(j) = perm(p); perm(p) = kk
    end do
    do j = 1, n2
      s = 10.d0**(-3.d0 * dble(perm(j)) / dble(n2))
      A0(:, oj+j-1) = A0(:, oj+j-1) * s
    end do
    A = A0
    ! ---- reference
    Aref = A0
    do j = 1, n2
      nrm(j) = dznrm2(m, Aref(oi, oj+j-1), 1)
      ipref(j) = j
    end do
    do i = 2, n2
      kk = ipref(i)
      j = i - 1
      do while (j >= 1)
        if (nrm(ipref(j)) >= nrm(kk)) exit
        ipref(j+1) = ipref(j)
        j = j - 1
      end do
      ipref(j+1) = kk
    end do
    do j = 1, n2
      A(oi:oi+m-1, oj+j-1) = A0(oi:oi+m-1, oj+ipref(j)-1)     ! A temporarily holds A0*P
    end do
    Aref(oi:oi+m-1, oj:oj+n2-1) = A(oi:oi+m-1, oj:oj+n2-1)
    allocate(AP(m, n2))
    AP = A(oi:oi+m-1, oj:oj+n2-1)
    call zgeqrf(m, n2, Aref(oi,oj), lda, tauref, wq, -1, info)
    lw = int(dble(wq(1))); allocate(w(lw))
    call zgeqrf(m, n2, Aref(oi,oj), lda, tauref, w, lw, info)
    ! ---- variant 1's real path: QDRP_decompose = ZGEQP3 with all columns fixed (jpvt = 1) on the pre-permuted matrix,
    ! with ALF's LWORK query (Prog/QDRP_decompose_mod.F90). la_qrpp uses ZGEQRF; compared below (HOST).
    allocate(A3(m, n2), tau3(n2), rw(2*m), jp(n2))
    A3 = AP
    jp = 1
    call zgeqp3(m, n2, A3, m, jp, tau3, wq, -1, rw, info)
    lw = int(dble(wq(1))); deallocate(w); allocate(w(lw))
    call zgeqp3(m, n2, A3, m, jp, tau3, w, lw, rw, info)
    ! ---- la_qrpp
    A = A0
    ip = -7
    call dev_enter_z(A, size(A)); call dev_enter_z(tau, size(tau)); call dev_enter_i(ip, size(ip))
    call la_qrpp(m, n2, A(oi,oj), lda, ip, tau)
    call la_sync()
    call dev_from_z(A, size(A)); call dev_from_z(tau, size(tau)); call dev_from_i(ip, size(ip))
    call dev_exit_z(A, size(A)); call dev_exit_z(tau, size(tau)); call dev_exit_i(ip, size(ip))
    write(lab,'(A,I0,A,I0)') 'qrpp N=', N, ' case', cs
    call check_int(trim(lab)//' ipvt', n2, ip, ipref)
    call check(trim(lab)//' outside', outside_diff(lda, n2+oj, A, A0, oi, oj, m, n2), 0.d0)
    if (is_host()) then
      ! variant 1's ZGEQP3(jpvt=1) vs la_qrpp's ZGEQRF: reported, bit-identity required (see report if it fails)
      call check(trim(lab)//' A vs ZGEQP3(jpvt=1)', relfro(m, n2, A(oi,oj), lda, A3, m), 0.d0)
      call check(trim(lab)//' tau vs ZGEQP3(jpvt=1)', maxabs_diff(n2, tau, tau3), 0.d0)
      call check(trim(lab)//' A', relfro(m, n2, A(oi,oj), lda, Aref(oi,oj), lda), 0.d0)
      call check(trim(lab)//' tau', maxabs_diff(n2, tau, tauref), 0.d0)
    end if
    dmax = 0.d0
    do j = 1, n2
      rjj = abs(Aref(oi+j-1, oj+j-1))
      dmax = fmax(dmax, abs(abs(A(oi+j-1, oj+j-1)) - rjj) / rjj)
    end do
    call check(trim(lab)//' |Rjj|', dmax, TOL_CUDA)
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
    ! A0 * P with the device ipvt
    do j = 1, n2
      k = ip(j)
      AP(:, j) = A0(oi:oi+m-1, oj+k-1)
    end do
    call check_inv(trim(lab)//' ||AP-QR||/||A||', relfro(m, n2, QR, m, AP, m), TOL_CUDA)
    call zgemm('C', 'N', n2, n2, m, (1.d0, 0.d0), Q(oi,oj), lda, Q(oi,oj), lda, (0.d0, 0.d0), G, n2)
    do i = 1, n2
      G(i,i) = G(i,i) - (1.d0, 0.d0)
    end do
    call check_inv(trim(lab)//' ||QhQ-I||', nrmf(n2, n2, G, n2), TOL_CUDA)
  end subroutine
end program t_qrpp
