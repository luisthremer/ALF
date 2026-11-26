! Tracked determinant phase of the UDV states (udv_state_mod, det_U; CGR phase without LU, ALF_CGR_LU in cgr1_mod).
!  1. Householder_det_phase(TAU) of a ZGEQRF against the LU phase of Q (ZUNGQR), with tau = 0 columns; negative
!     control: the product without the last reflector must differ.
!  2. reset: det_U = 1; decompose (side r and l, several products of random complex matrices, N = 4, 50, 200):
!     det_U equals the LU phase of U (DET_C_PHASE); assignment copies det_U.
!  3. N = 3200 diagonal matrix with random phases: det_U equals exp(i sum theta) (the raw LU product of a unitary
!     matrix underflows near this size, KNOWN_ISSUES "Determinant phase underflow"; the tracked phase never forms it).
Program udv_det_phase

  Use Random_Wrap
  Use UDV_State_mod, only: UDV_State, Householder_det_phase
  Use MyMats, only: DET_C_PHASE, MMULT

  Implicit none

  Real (Kind=Kind(0.d0)), parameter :: TOL = 1.d-10
  Integer :: Seed_Vec(1), nsz, it, i, j, n, lwork, info, sides
  Integer, parameter :: Nlist(3) = [4, 50, 200]
  Complex (Kind=Kind(0.d0)), allocatable :: A(:,:), Q(:,:), T1(:,:), tau(:), work(:)
  Complex (Kind=Kind(0.d0)) :: z, zl, zw
  Real (Kind=Kind(0.d0)) :: th, thsum, dmax
  Type (UDV_State) :: udv, udv2
  Character :: side

  Seed_vec(1) = 4782347
  Call Ranset(Seed_vec)
  dmax = 0.d0

  ! 1. reflector formula
  n = 50
  allocate(A(n,n), Q(n,n), tau(n), work(64*n))
  call random_matrix(A, n)
  A(:,7) = cmplx(0.d0, 0.d0, kind(0.d0))          ! a zero column below the diagonal part: tau = 0 there
  A(7,7) = cmplx(2.d0, 0.d0, kind(0.d0))
  A(n,n) = cmplx(0.3d0, 0.8d0, kind(0.d0))        ! complex last element: last reflector has tau /= 0
  lwork = size(work)
  call ZGEQRF(n, n, A, n, tau, work, lwork, info)
  Q = A
  call ZUNGQR(n, n, n, Q, n, tau, work, lwork, info)
  z  = Householder_det_phase(tau, n)
  zl = DET_C_PHASE(Q, n)
  zl = zl/abs(zl)
  call check('reflectors vs LU of Q', z, zl)
  if (tau(n) == cmplx(0.d0, 0.d0, kind(0.d0))) call fail('last reflector should be non-trivial')
  zw = Householder_det_phase(tau, n-1)               ! negative control
  if (abs(zw - zl) < 1.d-6) call fail('negative control: dropping the last reflector was not seen')
  deallocate(A, Q, tau, work)

  ! 2. decompose, both sides, N = 4, 50, 200
  do nsz = 1, size(Nlist)
     n = Nlist(nsz)
     allocate(A(n,n), T1(n,n), Q(n,n))
     do sides = 1, 2
        side = merge('r', 'l', sides == 1)
        call udv%alloc(n)
        if (udv%det_ok) call fail('alloc: det_ok must be false')
        call udv%reset(side)
        if (.not. udv%det_ok .or. udv%det_U /= cmplx(1.d0, 0.d0, kind(0.d0))) call fail('reset: det_U /= 1')
        do it = 1, 4
           call random_matrix(A, n)
           call MMULT(T1, A, udv%U)
           udv%U = T1
           call udv%decompose()
           if (.not. udv%det_ok) call fail('decompose: det_ok false')
           Q = udv%U
           zl = DET_C_PHASE(Q, n)
           call check('decompose N='//trim(itoa(n))//' side '//side, udv%det_U, zl/abs(zl))
        enddo
        call udv2%alloc(n)
        udv2 = udv
        if (udv2%det_U /= udv%det_U .or. .not. udv2%det_ok) call fail('assignment does not copy det_U')
        call udv%dealloc
        call udv2%dealloc
     enddo
     deallocate(A, T1, Q)
  enddo

  ! 3. N = 3200 diagonal, analytic phase
  n = 3200
  call udv%alloc(n)
  call udv%reset('r')
  thsum = 0.d0
  udv%U = cmplx(0.d0, 0.d0, kind(0.d0))
  do i = 1, n
     th = 6.2831853071795862d0*(ranf_wrap() - 0.5d0)
     thsum = thsum + th
     udv%U(i,i) = (0.5d0 + 1.5d0*ranf_wrap())*exp(cmplx(0.d0, th, kind(0.d0)))
  enddo
  call udv%decompose()
  call check('diagonal N=3200', udv%det_U, exp(cmplx(0.d0, thsum, kind(0.d0))))
  if (abs(abs(udv%det_U) - 1.d0) > 1.d-12) call fail('diagonal N=3200: |det_U| /= 1')
  call udv%dealloc

  write(6,'(A,ES10.3)') '49-udv-det-phase: OK, max |tracked - reference| = ', dmax

contains

  subroutine random_matrix(M, n)
    Integer, intent(in) :: n
    Complex (Kind=Kind(0.d0)), intent(out) :: M(n,n)
    Integer :: i, j
    do j = 1, n
       do i = 1, n
          M(i,j) = cmplx(4.d0*(ranf_wrap() - 0.5d0), 2.d0*(ranf_wrap() - 0.5d0), kind(0.d0))
       enddo
    enddo
  end subroutine random_matrix

  subroutine check(what, a, b)
    Character (len=*), intent(in) :: what
    Complex (Kind=Kind(0.d0)), intent(in) :: a, b
    write(6,'(A,A,ES10.3)') what, ': |tracked - reference| = ', abs(a - b)
    if (.not. (abs(a - b) <= TOL)) call fail(what)
    dmax = max(dmax, abs(a - b))
  end subroutine check

  subroutine fail(what)
    Character (len=*), intent(in) :: what
    write(6,*) '49-udv-det-phase: FAILED: ', what
    error stop 1
  end subroutine fail

  function itoa(k) result(s)
    Integer, intent(in) :: k
    Character (len=12) :: s
    write(s,'(I0)') k
    s = trim(s)
  end function itoa

end Program udv_det_phase
