! Unit test of OpT_time_dependent (Prog/OpT_time_dependent.F90), the time-dependent hopping element
! A(g) = U exp(g E) U^dagger (T = U E U^dagger, acting on the sites P) with g = g_t(t), used for the
! space-time point defect (docs/gpu/TIMEDEP_HOPPING.md).
!
! Reference: exp(g T) computed independently of Op_set's eigendecomposition, by a Taylor series with
! scaling and squaring, embedded into the Ndim x Ndim identity on the sites P.
! Tolerance 1d-12 (absolute, entries of order 1): the routines are a handful of dense products of
! Ndim = 7; observed deviations are ~1d-15. A wrong sign of g, a missing symmetric average,
! a wrong time index or a wrong site map gives O(1d-2) or larger.
!
! Checks
!  1  lmult(I, t)  = A(g_t(t))   and   rmult(I, t) = A(g_t(t))         (symm = .false.)
!  2  lmultinv(lmult(X)) = X,  rmultinv(rmult(X)) = X,  lmult(lmultinv(X)) = X   (A(g) A(-g) = 1)
!  3  lmultinv(I, t) = A(-g_t(t))
!  4  symm = .true.: g = (g_t(t) + g_t(t-1))/2, t-1 = Ltrot for t = 1 (periodic)
!  5  g_t(t) = 0 (and symmetric average 0): the element is skipped, X is returned bit-identical
!  6  adjointaction(X, t1, t2) = A(g_t(t1)/2) X A(-g_t(t2)/2);  t1 = 0 means Ltrot
Program OpT_time_dependent_test
  Use Operator_mod
  Use OpT_time_dependent_mod
  Implicit none

  Integer, parameter :: Ndim = 7, N = 4, Ltrot = 5
  Real (Kind=Kind(0.d0)), parameter :: tol = 1.d-12
  Type(Operator) :: Op
  Type(OpT_time_dependent) :: el, el_s
  Complex (Kind=Kind(0.d0)) :: X(Ndim,Ndim), Y(Ndim,Ndim), Id(Ndim,Ndim), Ref(Ndim,Ndim), g
  Complex (Kind=Kind(0.d0)) :: gt(Ltrot)
  Integer :: P(N), i, j, t, tm1, nfail
  Real (Kind=Kind(0.d0)) :: d

  nfail = 0
  P = (/ 2, 5, 3, 7 /)       ! non-contiguous, unsorted site map
  ! g_t: real and complex values, one zero slice (t = 4) and one slice with g(t)+g(t-1) = 0 (t = 3)
  gt = (/ cmplx(-0.3d0, 0.d0, kind(0.d0)), cmplx(0.45d0, 0.d0, kind(0.d0)), cmplx(-0.45d0, 0.d0, kind(0.d0)), &
       &  cmplx(0.d0, 0.d0, kind(0.d0)),   cmplx(0.2d0, -0.15d0, kind(0.d0)) /)

  Call Op_make(Op, N)
  Op%O(1,1) = cmplx( 0.7d0, 0.d0, kind(0.d0))
  Op%O(2,2) = cmplx(-0.4d0, 0.d0, kind(0.d0))
  Op%O(3,3) = cmplx( 0.1d0, 0.d0, kind(0.d0))
  Op%O(4,4) = cmplx( 1.2d0, 0.d0, kind(0.d0))
  Op%O(1,2) = cmplx(-1.0d0, 0.3d0, kind(0.d0));  Op%O(2,1) = conjg(Op%O(1,2))
  Op%O(1,4) = cmplx( 0.5d0,-0.2d0, kind(0.d0));  Op%O(4,1) = conjg(Op%O(1,4))
  Op%O(2,3) = cmplx(-0.8d0, 0.0d0, kind(0.d0));  Op%O(3,2) = conjg(Op%O(2,3))
  Op%O(3,4) = cmplx( 0.0d0, 0.6d0, kind(0.d0));  Op%O(4,3) = conjg(Op%O(3,4))
  Op%P = P
  Op%g = cmplx(1.d0, 0.d0, kind(0.d0))
  Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
  Allocate(Op%g_t(Ltrot))
  Op%g_t = gt
  Call Op_set(Op)
  If (.not. Op%get_g_t_alloc()) call fail('Op_set did not register g_t')

  Call el%init(Op, .false.)
  Call el_s%init(Op, .true.)

  Id = cmplx(0.d0, 0.d0, kind(0.d0))
  do i = 1, Ndim
     Id(i,i) = cmplx(1.d0, 0.d0, kind(0.d0))
  enddo
  ! a general (non-Hermitian) test matrix
  do j = 1, Ndim
     do i = 1, Ndim
        X(i,j) = cmplx(sin(1.3d0*i + 0.7d0*j), cos(0.4d0*i*j + 0.1d0), kind(0.d0))
     enddo
  enddo

  ! ---- 1 and 3: action on the identity vs the reference exp(g T)
  do t = 1, Ltrot
     if (t == 4) cycle
     Call A_ref(gt(t), Ref)
     Y = Id; Call el%lmult(Y, t);    Call check(Y, Ref, 'lmult(I) = A(g_t(t))', t)
     Y = Id; Call el%rmult(Y, t);    Call check(Y, Ref, 'rmult(I) = A(g_t(t))', t)
     Call A_ref(-gt(t), Ref)
     Y = Id; Call el%lmultinv(Y, t); Call check(Y, Ref, 'lmultinv(I) = A(-g_t(t))', t)
     Y = Id; Call el%rmultinv(Y, t); Call check(Y, Ref, 'rmultinv(I) = A(-g_t(t))', t)
     ! product with a general matrix from the left / right
     Call A_ref(gt(t), Ref)
     Y = X; Call el%lmult(Y, t);     Call check(Y, matmul(Ref, X), 'lmult(X) = A X', t)
     Y = X; Call el%rmult(Y, t);     Call check(Y, matmul(X, Ref), 'rmult(X) = X A', t)
  enddo

  ! ---- 2: inverse consistency A(g) A(-g) = 1
  do t = 1, Ltrot
     Y = X; Call el%lmult(Y, t); Call el%lmultinv(Y, t); Call check(Y, X, 'lmultinv(lmult(X)) = X', t)
     Y = X; Call el%lmultinv(Y, t); Call el%lmult(Y, t); Call check(Y, X, 'lmult(lmultinv(X)) = X', t)
     Y = X; Call el%rmult(Y, t); Call el%rmultinv(Y, t); Call check(Y, X, 'rmultinv(rmult(X)) = X', t)
     Y = X; Call el_s%lmult(Y, t); Call el_s%lmultinv(Y, t); Call check(Y, X, 'symm lmultinv(lmult(X)) = X', t)
  enddo

  ! ---- 4: symmetric averaging, periodic in t
  do t = 1, Ltrot
     tm1 = t - 1
     if (t == 1) tm1 = Ltrot
     g = 0.5d0*(gt(t) + gt(tm1))
     Call A_ref(g, Ref)
     Y = X; Call el_s%lmult(Y, t);    Call check(Y, matmul(Ref, X), 'symm lmult(X) = A((g(t)+g(t-1))/2) X', t)
     Y = X; Call el_s%rmult(Y, t);    Call check(Y, matmul(X, Ref), 'symm rmult(X) = X A((g(t)+g(t-1))/2)', t)
     Call A_ref(-g, Ref)
     Y = X; Call el_s%lmultinv(Y, t); Call check(Y, matmul(Ref, X), 'symm lmultinv(X)', t)
     Y = X; Call el_s%rmultinv(Y, t); Call check(Y, matmul(X, Ref), 'symm rmultinv(X)', t)
  enddo

  ! ---- 5: g = 0 skips the element (bit-identical); t = 4: g_t = 0; symm, t = 3: (g(3)+g(2))/2 = 0
  Y = X; Call el%lmult(Y, 4);    If (any(Y /= X)) call fail('lmult with g_t = 0 changed X')
  Y = X; Call el%rmult(Y, 4);    If (any(Y /= X)) call fail('rmult with g_t = 0 changed X')
  Y = X; Call el%lmultinv(Y, 4); If (any(Y /= X)) call fail('lmultinv with g_t = 0 changed X')
  Y = X; Call el%rmultinv(Y, 4); If (any(Y /= X)) call fail('rmultinv with g_t = 0 changed X')
  Y = X; Call el_s%lmult(Y, 3);  If (any(Y /= X)) call fail('symm lmult with average 0 changed X')
  Y = X; Call el_s%rmult(Y, 3);  If (any(Y /= X)) call fail('symm rmult with average 0 changed X')
  ! and t = 4 is not skipped in the symmetric form ((0 + g(3))/2 /= 0): covered by check 4

  ! ---- 6: adjoint action A(g(t1)/2) X A(-g(t2)/2)
  Call adj_check(2, 2)
  Call adj_check(1, 5)
  Call adj_check(0, 3)     ! t1 = 0 means Ltrot, t2 explicit
  Call adj_check(4, 2)     ! left factor skipped (g = 0)
  ! t2 absent: t2 = t1
  Y = X; Call el%adjointaction(Y, 5)
  Call A_ref(0.5d0*gt(5), Ref)
  Call check(Y, matmul(matmul(Ref, X), A_mat(-0.5d0*gt(5))), 'adjointaction(X, t) (t2 = t1)', 5)

  If (nfail > 0) then
     write(*,*) 'FAIL: ', nfail, ' checks failed'
     error stop 1
  endif
  write(*,*) 'PASS: OpT_time_dependent (lmult/rmult/inverse/symm/skip/adjointaction)'

contains

  subroutine adj_check(t1, t2)
    Integer, intent(in) :: t1, t2
    Integer :: n1, n2
    n1 = t1; if (t1 == 0) n1 = Ltrot
    n2 = t2; if (t2 == 0) n2 = Ltrot
    Y = X; Call el%adjointaction(Y, t1, t2)
    Call check(Y, matmul(matmul(A_mat(0.5d0*gt(n1)), X), A_mat(-0.5d0*gt(n2))), 'adjointaction(X, t1, t2)', 10*t1 + t2)
  end subroutine adj_check

  function A_mat(gg) result(R)
    Complex (Kind=Kind(0.d0)), intent(in) :: gg
    Complex (Kind=Kind(0.d0)) :: R(Ndim,Ndim)
    Call A_ref(gg, R)
  end function A_mat

  !> exp(gg * O) embedded on the sites P, by Taylor series with scaling and squaring (independent of Op%U, Op%E)
  subroutine A_ref(gg, R)
    Complex (Kind=Kind(0.d0)), intent(in)  :: gg
    Complex (Kind=Kind(0.d0)), intent(out) :: R(Ndim,Ndim)
    Complex (Kind=Kind(0.d0)) :: M(N,N), E(N,N), term(N,N)
    Integer :: k, s, a, b
    Integer, parameter :: nsq = 6
    M = gg * Op%O / dble(2**nsq)
    E = cmplx(0.d0, 0.d0, kind(0.d0)); term = E
    do a = 1, N
       E(a,a) = cmplx(1.d0, 0.d0, kind(0.d0)); term(a,a) = E(a,a)
    enddo
    do k = 1, 30
       term = matmul(term, M) / dble(k)
       E = E + term
    enddo
    do s = 1, nsq
       E = matmul(E, E)
    enddo
    R = Id
    do b = 1, N
       do a = 1, N
          R(P(a), P(b)) = E(a,b)
       enddo
    enddo
  end subroutine A_ref

  subroutine check(A, B, what, t)
    Complex (Kind=Kind(0.d0)), intent(in) :: A(:,:), B(:,:)
    Character (len=*), intent(in) :: what
    Integer, intent(in) :: t
    d = maxval(abs(A - B))
    if (.not. (d <= tol)) then
       write(*,'(A,A,A,I0,A,ES10.3)') 'FAIL: ', what, '  t = ', t, '  max dev = ', d
       nfail = nfail + 1
    endif
  end subroutine check

  subroutine fail(what)
    Character (len=*), intent(in) :: what
    write(*,*) 'FAIL: ', what
    nfail = nfail + 1
  end subroutine fail

end Program OpT_time_dependent_test
