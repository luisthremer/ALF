! Unit test for the underflow-safe determinant phase of MyMats (KNOWN_ISSUES.md,
! "Determinant phase underflow"):
!  1. DET_C_PHASE, LU_PHASE_SAFE (with start) and INV_C(..., DET_PHASE) return bitwise the result
!     of DET_C / of the raw loops when that result is usable (no change for normal runs);
!  2. on a matrix with |det| = 1 whose running LU product dips below the floating-point range,
!     DET_C returns 0, while DET_C_PHASE returns the exact phase with modulus 1.
! The matrix of part 2 is A = P U: P a permutation, U upper triangular with diagonal moduli
! 1e-20 (first half) and 1e+20 (second half). Partial pivoting then returns U's diagonal as pivots,
! so det A = sign(P) prod U(k,k) is known exactly. (The mirrored case, an overflow of the raw
! product, is handled the same way but not tested: -ffpe-trap=overflow in the CI debug job would
! stop DET_C there.)
program test_det_phase
  use MyMats, only: DET_C, DET_C_PHASE, LU_PHASE_SAFE, INV
  implicit none
  integer, parameter :: dp = kind(0.d0), n1 = 50, n2 = 40
  complex(dp) :: A(n1,n1), B(n1,n1), Ainv(n1,n1), z1, z2, d, dp1, start
  complex(dp) :: P2(n2,n2), U2(n2,n2), ref
  integer     :: ipiv(n1), info, i, k, perm(n2), n_fail
  real(dp)    :: r, theta
  external :: zgetrf

  n_fail = 0

  ! --- 1. bitwise equal in range ------------------------------------------------------------
  call fill(A)
  B = A; z1 = DET_C(B, n1)
  B = A; z2 = DET_C_PHASE(B, n1)
  call check('DET_C_PHASE == DET_C (in range)', z1 == z2)

  B = A
  call zgetrf(n1, n1, B, n1, ipiv, info)
  start = cmplx(0.3d0, -1.7d0, dp)
  d = start                                   ! raw loop as in CGR (STAB1/2): sign, then product
  do i = 1, n1
     if (ipiv(i) /= i) d = -d
     d = d*B(i,i)
  enddo
  call check('LU_PHASE_SAFE(raw, start) == raw loop', LU_PHASE_SAFE(d, B, ipiv, n1, start) == d)

  call INV(A, Ainv, z1, DET_PHASE=dp1)
  call check('INV DET_PHASE == DET (in range)', z1 == dp1)

  ! --- 2. underflow and overflow with |det| = 1 ----------------------------------------------
  do k = 1, 1
     call make_PU(k == 2, ref)
     U2 = P2
     z1 = DET_C(U2, n2)
     call check('DET_C underflows (shows the bug)', abs(z1) < 1.d-300)
     U2 = P2
     z2 = DET_C_PHASE(U2, n2)
     call check('DET_C_PHASE has modulus 1', abs(abs(z2) - 1.d0) < 1.d-12)
     call check('DET_C_PHASE has the exact phase', abs(z2 - ref) < 1.d-12)
  enddo

  if (n_fail > 0) then
     write(*,'(A,I0,A)') '14-det-phase: ', n_fail, ' check(s) FAILED'
     error stop 1
  endif
  write(*,'(A)') '14-det-phase: all checks passed'

contains

  subroutine check(what, ok)
    character(len=*), intent(in) :: what
    logical, intent(in) :: ok
    if (ok) then
       write(*,'(A,A)') '  ok    ', what
    else
       write(*,'(A,A)') '  FAIL  ', what
       n_fail = n_fail + 1
    endif
  end subroutine check

  ! Deterministic pseudo-random entries in [-1, 1]
  real(dp) function rnd()
    integer(8), save :: s = 12345_8
    s = mod(1103515245_8*s + 12345_8, 2147483647_8)   ! 64-bit: no overflow
    rnd = 2.d0*real(s, dp)/2147483647.d0 - 1.d0
  end function rnd

  subroutine fill(M)
    complex(dp), intent(out) :: M(:,:)
    integer :: i, j
    do j = 1, size(M,2)
       do i = 1, size(M,1)
          M(i,j) = cmplx(rnd(), rnd(), dp)
       enddo
    enddo
  end subroutine fill

  ! P2 = P U with |det| = 1; ref = exact det (|ref| = 1). up_first: moduli 1e+20 then 1e-20.
  subroutine make_PU(up_first, det)
    logical, intent(in) :: up_first
    complex(dp), intent(out) :: det
    integer :: i, j, t, sgn
    U2 = cmplx(0.d0, 0.d0, dp)
    det = cmplx(1.d0, 0.d0, dp)
    do j = 1, n2
       do i = 1, j - 1
          U2(i,j) = cmplx(rnd(), rnd(), dp)
       enddo
       theta = 3.d0*rnd()
       r = 1.d-20
       if ((j > n2/2) .neqv. up_first) r = 1.d20
       U2(j,j) = r*cmplx(cos(theta), sin(theta), dp)
       det = det*cmplx(cos(theta), sin(theta), dp)
    enddo
    ! a random permutation (Fisher-Yates) and its sign
    perm = [(i, i = 1, n2)]
    sgn = 1
    do i = n2, 2, -1
       j = 1 + int(real(i, dp)*0.5d0*(rnd() + 1.d0))
       if (j > i) j = i
       if (j /= i) then
          t = perm(i); perm(i) = perm(j); perm(j) = t
          sgn = -sgn
       endif
    enddo
    do i = 1, n2
       P2(perm(i),:) = U2(i,:)
    enddo
    if (sgn == -1) det = -det
  end subroutine make_PU

end program test_det_phase
