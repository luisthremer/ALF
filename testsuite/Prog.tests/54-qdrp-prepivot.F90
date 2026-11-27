! Unit test of QDRP_prepivot (test-only pre-pivoted QR, -DQR_PREPIVOT=1/2/3) in the default build:
! the variant is a runtime argument, so all variants are tested here.
! Checks per variant and matrix: info of the downstream ZUNMQR with the returned WORK/LWORK, Q*(D*R) = A*P,
! IPVT is a permutation (identity beyond N_part); variant 0 is bitwise QDRP_decompose; variant 3 with nb = 1
! makes ZGEQP3's pivots on graded matrices.
! Matrices: seeded graded (shuffled columns), duplicated (tied) columns, Hilbert, rectangular N_part < Ndim;
! block sizes include a partial last block and nb > N_part.
Program TEST_QDRP_PREPIVOT
  Use QDRP_mod, only: QDRP_decompose, QDRP_prepivot
  Implicit none
  Integer, parameter :: dp = Kind(0.d0)
  Integer :: nfail, imat, iv, inb, Ndim, Npart, seed
  Integer, parameter :: nbs(4) = [1, 4, 7, 64]

  nfail = 0
  seed = 12345
  do imat = 1, 5
     select case (imat)
     case (1); Ndim = 48; Npart = 48    ! graded, square
     case (2); Ndim = 50; Npart = 50    ! graded with tied columns
     case (3); Ndim = 12; Npart = 12    ! Hilbert
     case (4); Ndim = 40; Npart = 25    ! rectangular graded
     case (5); Ndim = 33; Npart = 33    ! graded, odd size
     end select
     call check_variant(0, 1)
     do iv = 1, 2
        call check_variant(iv, 1)
     enddo
     do inb = 1, size(nbs)
        call check_variant(3, nbs(inb))
     enddo
  enddo
  call check_nb1_pivots()

  if (nfail > 0) then
     write(*,*) 'FAILED checks:', nfail
     stop 1
  endif
  write(*,*) 'all QDRP_prepivot checks passed'

contains

  real(dp) function rnd()
    seed = mod(seed * 16807, 2147483647)
    if (seed < 0) seed = seed + 2147483647
    rnd = dble(seed) / 2147483647.d0
  end function rnd

  subroutine make_matrix(A, kind_)
    Complex(dp), intent(out) :: A(:,:)
    Integer, intent(in) :: kind_
    Integer :: i, j, m, n, p
    Complex(dp) :: col(size(A,1))
    m = size(A,1); n = size(A,2)
    if (kind_ == 3) then
       do j = 1, n
          do i = 1, m
             A(i,j) = 1.d0 / dble(i + j)
          enddo
       enddo
       return
    endif
    do j = 1, n
       do i = 1, m
          A(i,j) = cmplx(rnd() - 0.5d0, rnd() - 0.5d0, dp) * 10.d0**(-12.d0 * dble(j-1) / dble(n))
       enddo
    enddo
    if (kind_ == 2) then   ! ties: columns 2k equal to columns 2k-1
       do j = 2, n, 2
          A(:,j) = A(:,j-1)
       enddo
    endif
    do j = n, 2, -1        ! shuffle the columns
       p = 1 + int(rnd() * dble(j))
       if (p > j) p = j
       col = A(:,j); A(:,j) = A(:,p); A(:,p) = col
    enddo
  end subroutine make_matrix

  subroutine fail(msg, v, nb, x)
    Character(len=*), intent(in) :: msg
    Integer, intent(in) :: v, nb
    Real(dp), intent(in) :: x
    write(*,'(A,A,I2,A,I3,A,I3,A,I3,A,ES10.3)') 'FAIL ', msg, imat, ' v=', v, ' nb=', nb, ' Ndim=', Ndim, ' value ', x
    nfail = nfail + 1
  end subroutine fail

  subroutine check_variant(v, nb)
    Integer, intent(in) :: v, nb
    Complex(dp), allocatable :: A(:,:), M(:,:), C(:,:), M0(:,:), D0(:), D(:)
    Complex(dp), allocatable :: TAU(:), WORK(:), TAU0(:), WORK0(:)
    Integer, allocatable :: IPVT(:), IPVT0(:)
    Integer :: LWORK, LWORK0, i, j, info
    Logical :: seen(Npart)
    Real(dp) :: err, nrmA
    Integer :: s0

    s0 = seed
    allocate(A(Ndim, Npart), M(Ndim, Npart), C(Ndim, Npart), D(Npart), TAU(Npart), IPVT(Ndim))
    call make_matrix(A, merge(imat, 1, imat <= 3))
    M = A
    IPVT = 0
    call QDRP_prepivot(Ndim, Npart, M, D, IPVT, TAU, WORK, LWORK, v, nb)
    ! permutation
    seen = .false.
    do j = 1, Npart
       if (IPVT(j) < 1 .or. IPVT(j) > Npart) then
          call fail('IPVT out of range', v, nb, dble(IPVT(j))); exit
       endif
       seen(IPVT(j)) = .true.
    enddo
    if (.not. all(seen)) call fail('IPVT not a permutation', v, nb, 0.d0)
    if (v > 0) then
       do j = Npart + 1, Ndim
          if (IPVT(j) /= j) call fail('IPVT tail not identity', v, nb, dble(IPVT(j)))
       enddo
    endif
    ! reconstruction Q * (D R) with ZUNMQR and the returned WORK/LWORK
    C = 0
    do j = 1, Npart
       do i = 1, j
          C(i,j) = D(i) * M(i,j)   ! R = D * (upper triangle of M after the D separation)
       enddo
    enddo
    call ZUNMQR('L', 'N', Ndim, Npart, Npart, M, Ndim, TAU, C, Ndim, WORK, LWORK, info)
    if (info /= 0) call fail('ZUNMQR info', v, nb, dble(info))
    nrmA = maxval(abs(A))
    err = 0
    do j = 1, Npart
       err = max(err, maxval(abs(C(:,j) - A(:,IPVT(j)))))
    enddo
    err = err / nrmA
    if (.not. (err < 1.d-13)) call fail('reconstruction', v, nb, err)
    ! variant 0 is QDRP_decompose bitwise
    if (v == 0) then
       allocate(M0(Ndim, Npart), D0(Npart), TAU0(Npart), IPVT0(Ndim))
       M0 = A; IPVT0 = 0
       call QDRP_decompose(Ndim, Npart, M0, D0, IPVT0, TAU0, WORK0, LWORK0)
       if (any(M0 /= M) .or. any(D0 /= D) .or. any(TAU0 /= TAU) .or. any(IPVT0 /= IPVT) .or. LWORK0 /= LWORK) &
            & call fail('variant 0 not bitwise QDRP_decompose', v, nb, 0.d0)
       deallocate(M0, D0, TAU0, IPVT0, WORK0)
    endif
    deallocate(A, M, C, D, TAU, IPVT, WORK)
    seed = s0   ! same matrix for every variant of this case
  end subroutine check_variant

  subroutine check_nb1_pivots()
    ! variant 3 with nb = 1 (exact residual norms) against ZGEQP3 on graded matrices
    Complex(dp), allocatable :: A(:,:), M(:,:), D(:), TAU(:), WORK(:)
    Integer, allocatable :: IPVT(:), IPVT0(:)
    Integer :: LWORK, t
    do t = 1, 3
       Ndim = 30 + 7 * t; Npart = Ndim
       allocate(A(Ndim, Npart), M(Ndim, Npart), D(Npart), TAU(Npart), IPVT(Npart), IPVT0(Npart))
       call make_matrix(A, 1)
       M = A; IPVT0 = 0
       call QDRP_decompose(Ndim, Npart, M, D, IPVT0, TAU, WORK, LWORK)
       deallocate(WORK)
       M = A; IPVT = 0
       call QDRP_prepivot(Ndim, Npart, M, D, IPVT, TAU, WORK, LWORK, 3, 1)
       if (count(IPVT /= IPVT0) > 0) then
          imat = 10 + t
          call fail('v3 nb=1 pivots differ from ZGEQP3 (count)', 3, 1, dble(count(IPVT /= IPVT0)))
       endif
       deallocate(A, M, D, TAU, IPVT, IPVT0, WORK)
    enddo
  end subroutine check_nb1_pivots

End Program TEST_QDRP_PREPIVOT
