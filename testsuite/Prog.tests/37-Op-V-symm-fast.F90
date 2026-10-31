! Test of the fast Op_V_is_symmetric against the original dense implementation (reference copy below)
! and of its cost at large Ndim.
!  Part A: random operator sets (diagonal, non-diagonal, palindromic, disjoint/commuting, repeated
!          projector indices on diagonal and non-diagonal operators, complex g, N_FL = 2), Ndim <= 24: the logical result must equal the
!          reference and the value expected for the construction.
!  Part B: Ndim = 1600 (Hubbard-like, diagonal, N_FL = 2) and Ndim = 600 (non-diagonal bonds,
!          palindromic): each call must take < 5 s (the dense reference needs hours / minutes).

Program Test_Op_V_symm_fast
  Use Operator_mod
  Implicit None

  Type (Operator), allocatable :: Op_V(:,:)
  Integer :: n_fail, n_cases, icase, kind_case, Ndim, n_op, N_FL, nf, n, i, j, Nb
  Logical :: res_new, res_ref, expected
  Real (Kind=Kind(0.d0)) :: r(4), t0, t1
  Integer (Kind=8) :: c0, c1, crate

  n_fail  = 0
  n_cases = 0
  call random_seed_fixed()

  ! ---------------- Part A ----------------
  do icase = 1, 270
     kind_case = mod(icase-1, 9) + 1
     call random_number(r)
     Ndim = 6 + int(r(1)*18.d0)
     N_FL = 1 + mod(icase, 2)
     select case (kind_case)
     case (1)   ! diagonal 1x1 on-site, random order (Hubbard)            -> true
        n_op = Ndim
        allocate (Op_V(n_op, N_FL))
        do nf = 1, N_FL
           do n = 1, n_op
              call make_op(Op_V(n,nf), 1)
              Op_V(n,nf)%P(1) = 1 + mod(n*7, Ndim)
              Op_V(n,nf)%O(1,1) = rnd_real()
              Op_V(n,nf)%g = rnd_cplx()
           enddo
        enddo
        expected = .true.
     case (2)   ! diagonal 2x2 with overlapping supports and repeated P   -> true
        n_op = 2*Ndim
        allocate (Op_V(n_op, N_FL))
        do nf = 1, N_FL
           do n = 1, n_op
              call make_op(Op_V(n,nf), 2)
              Op_V(n,nf)%P(1) = 1 + rnd_int(Ndim)
              Op_V(n,nf)%P(2) = 1 + rnd_int(Ndim)       ! may equal P(1)
              Op_V(n,nf)%O(1,1) = rnd_real()
              Op_V(n,nf)%O(2,2) = rnd_real()
              Op_V(n,nf)%g = rnd_cplx()
           enddo
        enddo
        expected = .true.
     case (3)   ! random Hermitian bonds on a chain, random order          -> false
        n_op = Ndim
        allocate (Op_V(n_op, N_FL))
        do nf = 1, N_FL
           do n = 1, n_op
              call make_bond(Op_V(n,nf), 1 + mod(n*5, Ndim), 1 + mod(n*5+1, Ndim))
           enddo
        enddo
        expected = .false.
     case (4)   ! palindromic arrangement of random bonds                   -> true
        Nb = Ndim/2
        n_op = 2*Nb + 1
        allocate (Op_V(n_op, N_FL))
        do nf = 1, N_FL
           do n = 1, Nb + 1
              call make_bond(Op_V(n,nf), 1 + mod(n, Ndim), 1 + mod(n+1, Ndim))
           enddo
           do n = 1, Nb
              call copy_op(Op_V(n,nf), Op_V(n_op+1-n,nf))
           enddo
        enddo
        expected = .true.
     case (5)   ! bonds on disjoint supports (commuting), random order     -> true
        n_op = Ndim/2
        allocate (Op_V(n_op, N_FL))
        do nf = 1, N_FL
           do n = 1, n_op
              i = n_op + 1 - n                 ! reversed order, distinct supports
              call make_bond(Op_V(n,nf), 2*i-1, 2*i)
           enddo
        enddo
        expected = .true.
     case (6)   ! mixed: flavour 1 diagonal, flavour 2 asymmetric bonds    -> false if N_FL = 2
        n_op = Ndim
        allocate (Op_V(n_op, N_FL))
        do n = 1, n_op
           call make_op(Op_V(n,1), 1)
           Op_V(n,1)%P(1) = n
           Op_V(n,1)%O(1,1) = rnd_real()
           Op_V(n,1)%g = rnd_cplx()
        enddo
        do nf = 2, N_FL
           do n = 1, n_op
              call make_bond(Op_V(n,nf), 1 + mod(n*3, Ndim), 1 + mod(n*3+2, Ndim))
           enddo
        enddo
        expected = (N_FL == 1)
     case (7)   ! one operator                                             -> true
        n_op = 1
        allocate (Op_V(n_op, N_FL))
        do nf = 1, N_FL
           call make_bond(Op_V(1,nf), 1, 2)
        enddo
        expected = .true.
     case (8)   ! dense 3x3 operators with tiny coupling: asymmetry below/above tolerance
        n_op = 3
        allocate (Op_V(n_op, N_FL))
        do nf = 1, N_FL
           do n = 1, n_op
              call make_op(Op_V(n,nf), 3)
              Op_V(n,nf)%P = (/ 1, 2, 3 /)
              do i = 1, 3
                 do j = 1, 3
                    Op_V(n,nf)%O(i,j) = cmplx(rnd_real(), rnd_real(), kind(0.d0))
                 enddo
              enddo
              Op_V(n,nf)%O = Op_V(n,nf)%O + conjg(transpose(Op_V(n,nf)%O))
              Op_V(n,nf)%g = cmplx(10.d0**(-3 - mod(icase,6)), 0.d0, kind(0.d0))   ! 1e-3 .. 1e-8
           enddo
        enddo
        ! [g O_1, g O_2] ~ g^2: above 1e-10 for g >= 1e-4, below for g <= 1e-6
        expected = (Op_V(1,1)%g%re < 1.d-5)
     case (9)   ! non-diagonal 3x3 on P = (s, s, t) (O(1,2) /= 0 although P(1) = P(2)), built so that the
                ! accumulated embedding is diagonal: O(1,3) = -O(2,3) => E(s,t) = 0. All embedded
                ! operators then commute (-> true), mixed with diagonal 1x1 operators; a wrong
                ! accumulation of repeated indices leaves E(s,t) /= 0 and breaks the symmetry.
        n_op = Ndim
        allocate (Op_V(n_op, N_FL))
        do nf = 1, N_FL
           do n = 1, n_op
              if (mod(n,3) == 0) then
                 call make_op(Op_V(n,nf), 1)
                 Op_V(n,nf)%P(1) = 1 + rnd_int(Ndim)
                 Op_V(n,nf)%O(1,1) = rnd_real()
              else
                 call make_op(Op_V(n,nf), 3)
                 i = 1 + rnd_int(Ndim)
                 j = 1 + mod(i + rnd_int(Ndim-1), Ndim)     ! j /= i
                 Op_V(n,nf)%P = (/ i, i, j /)
                 Op_V(n,nf)%O(1,1) = rnd_real()
                 Op_V(n,nf)%O(2,2) = rnd_real()
                 Op_V(n,nf)%O(3,3) = rnd_real()
                 Op_V(n,nf)%O(1,2) = cmplx(rnd_real(), rnd_real(), kind(0.d0))
                 Op_V(n,nf)%O(1,3) = cmplx(rnd_real(), rnd_real(), kind(0.d0))
                 Op_V(n,nf)%O(2,3) = - Op_V(n,nf)%O(1,3)
                 Op_V(n,nf)%O(2,1) = conjg(Op_V(n,nf)%O(1,2))
                 Op_V(n,nf)%O(3,1) = conjg(Op_V(n,nf)%O(1,3))
                 Op_V(n,nf)%O(3,2) = conjg(Op_V(n,nf)%O(2,3))
              endif
              Op_V(n,nf)%g = rnd_cplx()
           enddo
        enddo
        expected = .true.
     end select

     res_ref = ref_is_symmetric(Op_V, N_FL, Ndim)
     res_new = Op_V_is_symmetric(Op_V, N_FL, Ndim)
     n_cases = n_cases + 1
     if (res_new .neqv. res_ref) then
        n_fail = n_fail + 1
        write(*,'(A,I0,A,I0,A,L1,A,L1)') 'FAIL: case ', icase, ' kind ', kind_case, ': new = ', res_new, ', ref = ', res_ref
     endif
     if (kind_case /= 8 .or. abs(log10(Op_V(1,1)%g%re) + 5.d0) > 0.6d0) then
        if (res_ref .neqv. expected) then
           n_fail = n_fail + 1
           write(*,'(A,I0,A,I0,A,L1,A,L1)') 'FAIL: case ', icase, ' kind ', kind_case, ': ref = ', res_ref, &
                &                            ', expected = ', expected
        endif
     endif
     call free_ops(Op_V)
  enddo
  write(*,'(A,I0,A,I0)') 'Part A: cases ', n_cases, ', failures ', n_fail

  ! ---------------- Part B ----------------
  ! Hubbard-like, Ndim = 1600, N_FL = 2
  Ndim = 1600
  N_FL = 2
  allocate (Op_V(Ndim, N_FL))
  do nf = 1, N_FL
     do n = 1, Ndim
        call make_op(Op_V(n,nf), 1)
        Op_V(n,nf)%P(1) = n
        Op_V(n,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.d0))
        Op_V(n,nf)%g = cmplx(0.3d0*(3-2*nf), 0.d0, kind(0.d0))
     enddo
  enddo
  call system_clock(c0, crate)
  res_new = Op_V_is_symmetric(Op_V, N_FL, Ndim)
  call system_clock(c1)
  t0 = real(c1-c0, kind(0.d0))/real(crate, kind(0.d0))
  write(*,'(A,L1,A,F8.3,A)') 'Part B diagonal Ndim=1600: result ', res_new, ', ', t0, ' s'
  if (.not. res_new .or. t0 > 5.d0) n_fail = n_fail + 1
  call free_ops(Op_V)

  ! Non-diagonal bonds, palindromic, Ndim = 600, N_FL = 1
  Ndim = 600
  N_FL = 1
  Nb   = Ndim - 1
  n_op = 2*Nb + 1
  allocate (Op_V(n_op, N_FL))
  do n = 1, Nb + 1
     call make_bond(Op_V(n,1), 1 + mod(n-1, Ndim), 1 + mod(n, Ndim))
  enddo
  do n = 1, Nb
     call copy_op(Op_V(n,1), Op_V(n_op+1-n,1))
  enddo
  call system_clock(c0, crate)
  res_new = Op_V_is_symmetric(Op_V, N_FL, Ndim)
  call system_clock(c1)
  t1 = real(c1-c0, kind(0.d0))/real(crate, kind(0.d0))
  write(*,'(A,L1,A,F8.3,A)') 'Part B palindromic bonds Ndim=600: result ', res_new, ', ', t1, ' s'
  if (.not. res_new .or. t1 > 5.d0) n_fail = n_fail + 1
  call free_ops(Op_V)

  if (n_fail > 0) then
     write(*,'(A,I0)') 'FAILURES: ', n_fail
     stop 2
  endif
  write(*,*) 'success'

contains

  subroutine random_seed_fixed()
    Integer :: nseed
    Integer, allocatable :: seed(:)
    call random_seed(size=nseed)
    allocate (seed(nseed))
    seed = 20261031
    call random_seed(put=seed)
  end subroutine random_seed_fixed

  Real (Kind=Kind(0.d0)) function rnd_real()
    call random_number(rnd_real)
    rnd_real = 2.d0*rnd_real - 1.d0
  end function rnd_real

  Complex (Kind=Kind(0.d0)) function rnd_cplx()
    rnd_cplx = cmplx(rnd_real(), 0.3d0*rnd_real(), kind(0.d0))
  end function rnd_cplx

  Integer function rnd_int(m)
    Integer, intent(in) :: m
    Real (Kind=Kind(0.d0)) :: x
    call random_number(x)
    rnd_int = min(int(x*m), m-1)
  end function rnd_int

  subroutine make_op(Op, N)
    Type (Operator), intent(inout) :: Op
    Integer, intent(in) :: N
    call Op_make(Op, N)
    Op%type = 2
  end subroutine make_op

  subroutine make_bond(Op, i1, i2)
    Type (Operator), intent(inout) :: Op
    Integer, intent(in) :: i1, i2
    call make_op(Op, 2)
    Op%P(1) = i1
    Op%P(2) = i2
    Op%O(1,1) = cmplx(rnd_real(), 0.d0, kind(0.d0))
    Op%O(2,2) = cmplx(rnd_real(), 0.d0, kind(0.d0))
    Op%O(1,2) = cmplx(rnd_real(), rnd_real(), kind(0.d0))
    Op%O(2,1) = conjg(Op%O(1,2))
    Op%g = rnd_cplx()
  end subroutine make_bond

  subroutine copy_op(Op_in, Op_out)
    Type (Operator), intent(in)    :: Op_in
    Type (Operator), intent(inout) :: Op_out
    call make_op(Op_out, Op_in%N)
    Op_out%P = Op_in%P
    Op_out%O = Op_in%O
    Op_out%g = Op_in%g
  end subroutine copy_op

  subroutine free_ops(Op_V)
    Type (Operator), allocatable, intent(inout) :: Op_V(:,:)
    Integer :: n, nf
    do nf = 1, size(Op_V,2)
       do n = 1, size(Op_V,1)
          deallocate (Op_V(n,nf)%O, Op_V(n,nf)%P)
       enddo
    enddo
    deallocate (Op_V)
  end subroutine free_ops

  ! Reference: the original dense implementation of Op_V_is_symmetric (ALF f81b5cbc,
  ! Prog/Operator_mod.F90), copied verbatim apart from the name.
  logical function ref_is_symmetric(Op_V_arr, N_FL_in, Ndim_in)
    Implicit None
    Type(Operator), intent(in) :: Op_V_arr(:,:)
    Integer, intent(in) :: N_FL_in, Ndim_in

    Integer :: n_op, n, nf, i, j, NN
    Complex(Kind=Kind(0.d0)), allocatable :: M_fwd(:,:), M_rev(:,:), M_tmp(:,:), M_op(:,:)
    Real(Kind=Kind(0.d0)) :: diff, scale
    Real(Kind=Kind(0.d0)), parameter :: sym_tol = 1.D-10

    ref_is_symmetric = .true.
    n_op = size(Op_V_arr, 1)
    NN = Ndim_in

    if (n_op <= 1) return

    allocate(M_fwd(NN, NN), M_rev(NN, NN), M_tmp(NN, NN), M_op(NN, NN))

    do nf = 1, N_FL_in
       M_fwd = cmplx(0.d0, 0.d0, kind(0.d0))
       do i = 1, NN
          M_fwd(i,i) = cmplx(1.d0, 0.d0, kind(0.d0))
       enddo
       do n = 1, n_op
          M_op = cmplx(0.d0, 0.d0, kind(0.d0))
          do i = 1, NN
             M_op(i,i) = cmplx(1.d0, 0.d0, kind(0.d0))
          enddo
          do j = 1, Op_V_arr(n, nf)%N
             do i = 1, Op_V_arr(n, nf)%N
                M_op(Op_V_arr(n, nf)%P(i), Op_V_arr(n, nf)%P(j)) = &
                   M_op(Op_V_arr(n, nf)%P(i), Op_V_arr(n, nf)%P(j)) + &
                   Op_V_arr(n, nf)%g * Op_V_arr(n, nf)%O(i, j)
             enddo
          enddo
          M_tmp = matmul(M_fwd, M_op)
          M_fwd = M_tmp
       enddo
       M_rev = cmplx(0.d0, 0.d0, kind(0.d0))
       do i = 1, NN
          M_rev(i,i) = cmplx(1.d0, 0.d0, kind(0.d0))
       enddo
       do n = n_op, 1, -1
          M_op = cmplx(0.d0, 0.d0, kind(0.d0))
          do i = 1, NN
             M_op(i,i) = cmplx(1.d0, 0.d0, kind(0.d0))
          enddo
          do j = 1, Op_V_arr(n, nf)%N
             do i = 1, Op_V_arr(n, nf)%N
                M_op(Op_V_arr(n, nf)%P(i), Op_V_arr(n, nf)%P(j)) = &
                   M_op(Op_V_arr(n, nf)%P(i), Op_V_arr(n, nf)%P(j)) + &
                   Op_V_arr(n, nf)%g * Op_V_arr(n, nf)%O(i, j)
             enddo
          enddo
          M_tmp = matmul(M_rev, M_op)
          M_rev = M_tmp
       enddo
       diff = maxval(abs(M_fwd - M_rev))
       scale = max(1.d0, maxval(abs(M_fwd)), maxval(abs(M_rev)))
       if (diff > sym_tol * scale) then
          ref_is_symmetric = .false.
          deallocate(M_fwd, M_rev, M_tmp, M_op)
          return
       endif
    enddo

    deallocate(M_fwd, M_rev, M_tmp, M_op)
  end function ref_is_symmetric

end Program Test_Op_V_symm_fast
