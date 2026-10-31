! Ticket 207 dev5 (development only, debug aid): offline accuracy of the tau-pass CGR2_2 against a quad-precision reference.
! Reads the dumps of ALF_GPU_TAUDIAG_QREF (Prog/gpu_tau_mod.F90 qref_post: inputs and outputs of the production CGR2_2 at one
! stabilization point) and prints one line per dump with the errors, against the quad reference Gq of the same inputs, of
!   dev  the device CGR2_2 (production output), var<m> device variants (ALF_GPU_TAUDIAG_QVAR masks), nvh the host CGR2_2 of the
!   nvfortran binary (QHOST), cpu the host CGR2_2 of cgr2_2_mod compiled here (gfortran; default: -DQR_PREPIVOT variant 1,
!   LOG: -DSTABLOG ZGEQP3), n<code> rounding-order null variants of the host CGR2_2 (copy of nul_cgr22 of gpu_tau_mod), u<k> the
!   host CGR2_2 on the right state with U perturbed by ~1 ulp (sensitivity, measured against the unperturbed Gq).
! Error of a method M per block b (G00, GT0, GTT, G0T): max|M_b - Gq_b| / max|Gq_b|; eF: Frobenius norm of M - Gq over the four
! blocks / that of Gq.  All differences are formed in quad.  With the next block's device propagation B, Bi in the dump (version 2,
! not after the last point of a pass) three more numbers: the error as the next precision check receives it, max|B dTT Bi|,
! max|B dT0|, max|d0T Bi| (absolute, as Control_Precision_tau; d = M - Gq formed in quad, rounded to double, products in double);
! -1 when absent.
! The quad reference qcgr22 is CGR2_2 in real(16) / complex(16) (IEEE binary128, libquadmath): same formulae, same branch, same
! scale splitting (STABLOG), V1^-1 by LU with partial pivoting, the 2N x 2N QR by Householder reflections (LAPACK zlarfg
! convention) after the column-norm pre-pivot (default) or with full column pivoting on exact residual norms (STABLOG, as
! ZGEQP3), D3 = |R_ii| separated, the extended system solved as in solve_extended_System (rows permuted by IPVT, R^H solve, rows
! scaled by conj(1/D3), Q applied).  Check: a second quad run on the randomly row-permuted QR input (qself, every QS-th dump).
!   build (t207_d5_qbuild.sh): three parts of this file, -DQ5_PART=1 (module q5: quad reference, gfortran -O2, no fast-math),
!          2 (module q5nul: the null copy of CGR2_2, compiled with the build tree's ALF_FLAGS_PROG like cgr2_2_mod, so that code 0
!          reproduces the oracle bit for bit), 3 (the program, -O2); linked with cgr2_2_mod / udv_state_mod / QDRP_decompose_mod.
!   run:   t207_d5_qref.out QS NULCODES HOSTIN < list_of_dump_files     (NULCODES: comma list or '-'; HOSTIN 0|1)
#if Q5_PART == 1
module q5
  use iso_fortran_env, only: int64
  implicit none
  integer, parameter :: dp = kind(0.d0), qp = 16
contains

  real(dp) function rnd(s)
    integer(int64), intent(inout) :: s
    s = mod(1103515245_int64*s + 12345_int64, 2147483648_int64)
    rnd = real(s, dp)/2147483648.d0 - 0.5d0
  end function rnd

  subroutine nperm(n, p, seed)
    integer, intent(in) :: n
    integer, intent(out) :: p(n)
    integer(int64), intent(inout) :: seed
    integer :: i, j, t
    do i = 1, n
       p(i) = i
    enddo
    do i = n, 2, -1
       j = 1 + int((rnd(seed) + 0.5d0)*dble(i))
       j = min(max(j, 1), i)
       t = p(i); p(i) = p(j); p(j) = t
    enddo
  end subroutine nperm

  ! X = A^-1 in quad: LU with partial pivoting, then the solve with the identity
  subroutine qinv(n, Ain, X)
    integer, intent(in) :: n
    complex(qp), intent(in) :: Ain(n,n)
    complex(qp), intent(out) :: X(n,n)
    complex(qp), allocatable :: LU(:,:), rw(:)
    integer, allocatable :: ip(:)
    integer :: i, j, k, p
    real(qp) :: amax, a
    complex(qp) :: t
    allocate(LU(n,n), ip(n), rw(n))
    LU = Ain
    do k = 1, n
       p = k; amax = abs(LU(k,k))
       do i = k+1, n
          a = abs(LU(i,k))
          if (a > amax) then
             amax = a; p = i
          endif
       enddo
       ip(k) = p
       if (p /= k) then
          rw = LU(k,:); LU(k,:) = LU(p,:); LU(p,:) = rw
       endif
       t = 1.0_qp/LU(k,k)
       LU(k+1:n,k) = LU(k+1:n,k)*t
       !$omp parallel do private(i, t) schedule(static)
       do j = k+1, n
          t = LU(k,j)
          do i = k+1, n
             LU(i,j) = LU(i,j) - LU(i,k)*t
          enddo
       enddo
       !$omp end parallel do
    enddo
    X = (0.0_qp, 0.0_qp)
    do i = 1, n
       X(i,i) = (1.0_qp, 0.0_qp)
    enddo
    do k = 1, n
       if (ip(k) /= k) then
          rw = X(k,:); X(k,:) = X(ip(k),:); X(ip(k),:) = rw
       endif
    enddo
    !$omp parallel do private(i, k, t) schedule(static)
    do j = 1, n
       do k = 1, n                       ! L (unit lower)
          t = X(k,j)
          do i = k+1, n
             X(i,j) = X(i,j) - LU(i,k)*t
          enddo
       enddo
       do k = n, 1, -1                   ! U
          X(k,j) = X(k,j)/LU(k,k)
          t = X(k,j)
          do i = 1, k-1
             X(i,j) = X(i,j) - LU(i,k)*t
          enddo
       enddo
    enddo
    !$omp end parallel do
  end subroutine qinv

  ! CGR2_2 in quad.  islog: STABLOG inputs (S = L) else default (S = D).  br: 1 = D1(1) > D2(1) branch.  rps > 0: QR input rows
  ! randomly permuted with seed rps (self-check).  G(:,:,1..4) = G00, GT0, GTT, G0T.
  subroutine qcgr22(n, islog, br, U2, D2, L2, V2, U1, D1, L1, V1, rps, G)
    integer, intent(in) :: n, br
    logical, intent(in) :: islog
    complex(dp), intent(in) :: U2(n,n), V2(n,n), U1(n,n), V1(n,n), D2(n), D1(n)
    real(dp), intent(in) :: L2(n), L1(n)
    integer(int64), intent(in) :: rps
    complex(qp), intent(out) :: G(n,n,4)
    complex(qp), allocatable :: MYU2(:,:), V1INV(:,:), H(:,:), Y(:,:), Rh(:,:), T(:,:), tau(:), D1m(:), D2m(:)
    real(qp), allocatable :: D3(:), cn(:)
    integer, allocatable :: perm(:), rp(:)
    integer :: n2, i, j, k, p, kk
    real(qp) :: xn2, bnorm, a
    complex(qp) :: alpha, beta, w, sc, z
    integer(int64) :: s
    n2 = 2*n
    allocate(MYU2(n,n), V1INV(n,n), H(n2,n2), Y(n2,n2), Rh(n2,n2), tau(n2), D1m(n), D2m(n), D3(n2), cn(n2), perm(n2))
    MYU2 = conjg(transpose(cmplx(U2, kind=qp)))
    call qinv(n, cmplx(V1, kind=qp), V1INV)
    if (islog) then
       do j = 1, n
          if (L1(j) <= 0.d0) then
             D1m(j) = cmplx(exp(real(L1(j), qp)), 0.0_qp, qp)
          else
             D1m(j) = (1.0_qp, 0.0_qp)
             V1INV(j,:) = V1INV(j,:)*exp(-real(L1(j), qp))
          endif
          if (L2(j) <= 0.d0) then
             D2m(j) = cmplx(exp(real(L2(j), qp)), 0.0_qp, qp)
          else
             D2m(j) = (1.0_qp, 0.0_qp)
             MYU2(j,:) = MYU2(j,:)*exp(-real(L2(j), qp))
          endif
       enddo
    else
       D1m = cmplx(D1, kind=qp); D2m = cmplx(D2, kind=qp)
    endif
    ! HLPB2, then H = HLPB1 = HLPB2^H
    if (br == 1) then
       Y(1:n,1:n) = V1INV; Y(n+1:n2,n+1:n2) = MYU2
       do j = 1, n
          do i = 1, n
             Y(i  , j+n) =  D1m(i)*conjg(cmplx(U1(j,i), kind=qp))
             Y(i+n, j  ) = -D2m(i)*cmplx(V2(i,j), kind=qp)
          enddo
       enddo
    else
       Y(1:n,1:n) = MYU2; Y(n+1:n2,n+1:n2) = V1INV
       do j = 1, n
          do i = 1, n
             Y(i  , j+n) = -D2m(i)*cmplx(V2(i,j), kind=qp)
             Y(i+n, j  ) =  D1m(i)*conjg(cmplx(U1(j,i), kind=qp))
          enddo
       enddo
    endif
    H = conjg(transpose(Y))
    do j = 1, n2
       perm(j) = j
    enddo
    if (.not. islog) then
       ! pre-pivot (CGR2_2 QR_PREPIVOT variant 1 / la_qrpp): column 2-norms, stable descending insertion sort
       do j = 1, n2
          cn(j) = sqrt(sum(real(H(:,j))**2 + aimag(H(:,j))**2))
       enddo
       do i = 2, n2
          k = perm(i); j = i - 1
          do while (j >= 1)
             if (cn(perm(j)) >= cn(k)) exit
             perm(j+1) = perm(j); j = j - 1
          enddo
          perm(j+1) = k
       enddo
       Y = H
       do j = 1, n2
          H(:,j) = Y(:,perm(j))
       enddo
    endif
    if (rps > 0) then
       allocate(rp(n2)); s = rps
       call nperm(n2, rp, s)
       Y = H
       do i = 1, n2
          H(i,:) = Y(rp(i),:)
       enddo
    endif
    ! Householder QR (STABLOG: column pivoting on the exact residual norms, as ZGEQP3)
    do k = 1, n2
       if (islog) then
          p = k; bnorm = -1.0_qp
          do j = k, n2
             a = sum(real(H(k:n2,j))**2 + aimag(H(k:n2,j))**2)
             if (a > bnorm) then
                bnorm = a; p = j
             endif
          enddo
          if (p /= k) then
             Y(:,1) = H(:,k); H(:,k) = H(:,p); H(:,p) = Y(:,1)
             kk = perm(k); perm(k) = perm(p); perm(p) = kk
          endif
       endif
       alpha = H(k,k)
       xn2 = 0.0_qp
       if (k < n2) xn2 = sum(real(H(k+1:n2,k))**2 + aimag(H(k+1:n2,k))**2)
       if (xn2 == 0.0_qp .and. aimag(alpha) == 0.0_qp) then
          tau(k) = (0.0_qp, 0.0_qp)
       else
          beta = cmplx(-sign(sqrt(real(alpha)**2 + aimag(alpha)**2 + xn2), real(alpha)), 0.0_qp, qp)
          tau(k) = (beta - alpha)/beta
          sc = 1.0_qp/(alpha - beta)
          H(k+1:n2,k) = H(k+1:n2,k)*sc
          H(k,k) = beta
          z = conjg(tau(k))
          !$omp parallel do private(i, w) schedule(static)
          do j = k+1, n2
             w = H(k,j)
             do i = k+1, n2
                w = w + conjg(H(i,k))*H(i,j)
             enddo
             w = z*w
             H(k,j) = H(k,j) - w
             do i = k+1, n2
                H(i,j) = H(i,j) - H(i,k)*w
             enddo
          enddo
          !$omp end parallel do
       endif
    enddo
    ! D3 = |R_ii|, R~ = D3^-1 R; Rh = R~^H (lower)
    Rh = (0.0_qp, 0.0_qp)
    do i = 1, n2
       D3(i) = abs(H(i,i))
       do j = i, n2
          Rh(j,i) = conjg(H(i,j)/D3(i))
       enddo
    enddo
    ! RHS [[UCT,0],[0,VINV]], rows in the order of ZLAPMR(forward, perm): new row i = old row perm(i)
    allocate(T(n2,n2))
    T = (0.0_qp, 0.0_qp)
    if (br == 1) then
       T(1:n,1:n) = V1INV; T(n+1:n2,n+1:n2) = MYU2
    else
       T(1:n,1:n) = MYU2; T(n+1:n2,n+1:n2) = V1INV
    endif
    do i = 1, n2
       Y(i,:) = T(perm(i),:)
    enddo
    ! R~^H Z = Y, rows * conj(1/D3), then Q * Z (Q = H_1 ... H_n2)
    !$omp parallel do private(i, k, w) schedule(static)
    do j = 1, n2
       do k = 1, n2
          Y(k,j) = Y(k,j)/Rh(k,k)
          w = Y(k,j)
          do i = k+1, n2
             Y(i,j) = Y(i,j) - Rh(i,k)*w
          enddo
       enddo
       do i = 1, n2
          Y(i,j) = Y(i,j)/D3(i)
       enddo
       do k = n2, 1, -1
          w = Y(k,j)
          do i = k+1, n2
             w = w + conjg(H(i,k))*Y(i,j)
          enddo
          w = tau(k)*w
          Y(k,j) = Y(k,j) - w
          do i = k+1, n2
             Y(i,j) = Y(i,j) - H(i,k)*w
          enddo
       enddo
    enddo
    !$omp end parallel do
    if (rps > 0) then                    ! Q = P^T Q': rows back
       T = Y
       do i = 1, n2
          Y(rp(i),:) = T(i,:)
       enddo
    endif
    ! get_blocks
    if (br == 1) then
       G(:,:,1) = Y(1:n,1:n); G(:,:,4) = Y(1:n,n+1:n2); G(:,:,2) = Y(n+1:n2,1:n); G(:,:,3) = Y(n+1:n2,n+1:n2)
    else
       G(:,:,3) = Y(1:n,1:n); G(:,:,2) = Y(1:n,n+1:n2); G(:,:,4) = Y(n+1:n2,1:n); G(:,:,1) = Y(n+1:n2,n+1:n2)
    endif
  end subroutine qcgr22

  ! e(1:4) per block max|M - Gq| / max|Gq|, e(5) Frobenius over the four blocks
  subroutine errs(n, M, Gq, e)
    integer, intent(in) :: n
    complex(dp), intent(in) :: M(n,n,4)
    complex(qp), intent(in) :: Gq(n,n,4)
    real(dp), intent(out) :: e(5)
    real(qp) :: md, mb, d, sd, sb
    integer :: b, i, j
    sd = 0; sb = 0
    do b = 1, 4
       md = 0; mb = 0
       do j = 1, n
          do i = 1, n
             d = abs(cmplx(M(i,j,b), kind=qp) - Gq(i,j,b)); md = max(md, d); sd = sd + d*d
             d = abs(Gq(i,j,b)); mb = max(mb, d); sb = sb + d*d
          enddo
       enddo
       e(b) = real(md/max(mb, tiny(1.0_qp)), dp)
       if (.not. (abs(e(b)) <= huge(1.d0))) e(b) = huge(1.d0)
    enddo
    e(5) = real(sqrt(sd/max(sb, tiny(1.0_qp))), dp)
    if (.not. (abs(e(5)) <= huge(1.d0))) e(5) = huge(1.d0)
  end subroutine errs

  subroutine qdiff(n, Qa, Qb, e)
    integer, intent(in) :: n
    complex(qp), intent(in) :: Qa(n,n,4), Qb(n,n,4)
    real(dp), intent(out) :: e(5)
    real(qp) :: md, mb, d, sd, sb
    integer :: b, i, j
    sd = 0; sb = 0
    do b = 1, 4
       md = 0; mb = 0
       do j = 1, n
          do i = 1, n
             d = abs(Qa(i,j,b) - Qb(i,j,b)); md = max(md, d); sd = sd + d*d
             d = abs(Qb(i,j,b)); mb = max(mb, d); sb = sb + d*d
          enddo
       enddo
       e(b) = real(md/max(mb, tiny(1.0_qp)), dp)
    enddo
    e(5) = real(sqrt(sd/max(sb, tiny(1.0_qp))), dp)
  end subroutine qdiff
end module q5
#endif

#if Q5_PART == 2
module q5nul
  ! copy of nul_cgr22 of Prog/gpu_tau_mod.F90 (ticket 207 dev4): host CGR2_2 with rounding-order null variants; code c = flags +
  ! 100*nb + 100000*stream; flags 1 V1^-1 with permuted columns, 2 unblocked QR / Q, 4 row-permuted QR input, 8 compact WY (nb)
  use iso_fortran_env, only: int64
  use q5, only: nperm
  use UDV_State_mod
  implicit none
contains
  Subroutine nul_cgr22(LQ, GRT0, GR00, GRTT, GR0T, udv2, udv1, code, seed)
    Use QDRP_mod, only: QDRP_decompose
    Use cgr2_2_mod, only: solve_extended_System, get_blocks
    Use MyMats, only: INV
    Implicit none
    Integer, intent(in) :: LQ
    Complex (Kind=Kind(0.d0)), intent(inout) :: GRT0(:,:), GR00(:,:), GRTT(:,:), GR0T(:,:)
    Class (UDV_State), intent(in) :: udv2, udv1
    Integer, intent(in) :: code
    Integer(int64), intent(inout) :: seed
    Complex (Kind=Kind(0.d0)), allocatable, Dimension(:) :: D3, D1m, D2m, TAU, WORK
    Complex (Kind=Kind(0.d0)), allocatable, Dimension(:,:) :: MYU2, HLPB1, HLPB2, V1INV, TMP, TT
    Integer, allocatable :: IPVT(:), RP(:), CP(:)
    Complex (Kind=Kind(0.d0)) :: Z
    Integer :: LQ2, I, J, LWORK, fl, nb, mode
    Logical :: br1, rowp

    LQ2 = 2*LQ
    fl = mod(code, 100); nb = mod(code/100, 1000)
    rowp = iand(fl, 4) /= 0
    mode = 0
    if (iand(fl, 2) /= 0) mode = 2
    if (iand(fl, 8) /= 0) mode = 8
    ALLOCATE(MYU2(LQ, LQ), V1INV(LQ,LQ), HLPB1(LQ2, LQ2), HLPB2(LQ2, LQ2), D3(LQ2), D1m(LQ), D2m(LQ))
    Allocate(IPVT(LQ2), TAU(LQ2))
    IPVT = 0
    MYU2 = CONJG(TRANSPOSE(udv2%U))
    if (iand(fl, 1) /= 0) then
       allocate(CP(LQ), TMP(LQ,LQ), TT(LQ,LQ))
       call nperm(LQ, CP, seed)
       do I = 1, LQ
          TMP(:,I) = udv1%V(:,CP(I))
       enddo
       CALL INV(TMP, TT, Z)
       do I = 1, LQ
          V1INV(CP(I),:) = TT(I,:)
       enddo
       deallocate(CP, TMP, TT)
    else
       CALL INV(udv1%V, V1INV, Z)
    endif
#if defined(STABLOG)
    DO J=1,LQ
      If(udv1%L(J) <= 0.d0) then
        D1m(J)=cmplx(exp(udv1%L(J)),0.d0,kind(0.d0))
      else
        D1m(J)=cmplx(1.d0,0.d0,kind(0.d0))
        call zscal(LQ,cmplx(exp(-udv1%L(J)),0.d0,kind(0.d0)),V1INV(J,1),LQ)
      endif
      If(udv2%L(J) <= 0.d0) then
        D2m(J)=cmplx(exp(udv2%L(J)),0.d0,kind(0.d0))
      else
        D2m(J)=cmplx(1.d0,0.d0,kind(0.d0))
        call zscal(LQ,cmplx(exp(-udv2%L(J)),0.d0,kind(0.d0)),MYU2(J,1),LQ)
      endif
    ENDDO
    br1 = udv1%L(1) > udv2%L(1)
#else
    D1m=udv1%D
    D2m=udv2%D
    br1 = dble(udv1%D(1)) > dble(udv2%D(1))
#endif
    if (br1) then
       call zlacpy('A',LQ,LQ,V1INV(1,1) ,LQ,HLPB2(1   ,1   ),LQ2)
       call zlacpy('A',LQ,LQ,MYU2(1,1)  ,LQ,HLPB2(1+LQ,1+LQ),LQ2)
       DO J = 1,LQ
          DO I = 1,LQ
             HLPB2(I   , J+LQ ) =  D1m(I)*conjg(udv1%U(J,I))
             HLPB2(I+LQ, J    ) = -D2m(I)*udv2%V(I,J)
          ENDDO
       ENDDO
    else
       call zlacpy('A',LQ,LQ,MYU2(1,1) ,LQ,HLPB2(1   ,1   ),LQ2)
       call zlacpy('A',LQ,LQ,V1INV(1,1),LQ,HLPB2(1+LQ,1+LQ),LQ2)
       DO J = 1,LQ
          DO I = 1,LQ
             HLPB2(I   , J+LQ ) = -D2m(I)*udv2%V(I,J)
             HLPB2(I+LQ, J    ) =  D1m(I)*conjg(udv1%U(J,I))
          ENDDO
       ENDDO
    endif
    HLPB1 = conjg(transpose(HLPB2))
#if defined(QR_PREPIVOT) && !defined(STABLOG)
    block
      Real(Kind=Kind(0.d0)), allocatable :: pp_nrm(:)
      Real(Kind=Kind(0.d0)) :: DZNRM2
      Complex(Kind=Kind(0.d0)), allocatable :: pp_tmp(:,:)
      Integer, allocatable :: pp_perm(:)
      Integer :: pp_i, pp_j, pp_k
      Allocate(pp_nrm(LQ2), pp_tmp(LQ2,LQ2), pp_perm(size(IPVT)))
      do pp_j = 1, LQ2
         pp_nrm(pp_j) = DZNRM2(LQ2, HLPB1(1,pp_j), 1)
         pp_perm(pp_j) = pp_j
      enddo
      do pp_i = 2, LQ2
         pp_k = pp_perm(pp_i)
         pp_j = pp_i - 1
         do while (pp_j >= 1)
            if (pp_nrm(pp_perm(pp_j)) >= pp_nrm(pp_k)) exit
            pp_perm(pp_j+1) = pp_perm(pp_j)
            pp_j = pp_j - 1
         enddo
         pp_perm(pp_j+1) = pp_k
      enddo
      do pp_j = 1, LQ2
         pp_tmp(:,pp_j) = HLPB1(1:LQ2,pp_perm(pp_j))
      enddo
      HLPB1(1:LQ2,1:LQ2) = pp_tmp
      IPVT = 1
      if (rowp) call nul_rowperm()
      call nul_qr(.true.)
      IPVT = pp_perm
      Deallocate(pp_nrm, pp_tmp, pp_perm)
    end block
#else
    if (rowp) call nul_rowperm()
    call nul_qr(.false.)
#endif
    if (br1) then
       call nul_solve(V1INV, MYU2)
    else
       call nul_solve(MYU2, V1INV)
    endif
    if (rowp) then
       allocate(TMP(LQ2,LQ2))
       do I = 1, LQ2
          TMP(RP(I),:) = HLPB2(I,:)
       enddo
       HLPB2 = TMP
       deallocate(TMP, RP)
    endif
    if (br1) then
       call get_blocks(GR00, GR0T, GRT0, GRTT, HLPB2, LQ)
    else
       call get_blocks(GRTT, GRT0, GR0T, GR00, HLPB2, LQ)
    endif
    DEALLOCATE(MYU2, V1INV, HLPB1, HLPB2, IPVT, TAU, D3, D1m, D2m)
    if (allocated(WORK)) deallocate(WORK)
    if (allocated(TT)) deallocate(TT)
  contains
    Subroutine nul_rowperm()
      allocate(RP(LQ2), TMP(LQ2,LQ2))
      call nperm(LQ2, RP, seed)
      do I = 1, LQ2
         TMP(I,:) = HLPB1(RP(I),:)
      enddo
      HLPB1 = TMP
      deallocate(TMP)
    end subroutine nul_rowperm
    Subroutine nul_qr(fixed)
      Logical, intent(in) :: fixed
      Complex (Kind=Kind(0.d0)), allocatable :: ZW(:)
      Real (Kind=Kind(0.d0)), allocatable :: VN1(:), VN2(:)
      Real (Kind=Kind(0.d0)) :: X, DZNRM2
      Integer :: info, ii, jj
      if (mode == 0) then
         call QDRP_decompose(LQ2, LQ2, HLPB1, D3, IPVT, TAU, WORK, LWORK)
         return
      endif
      if (mode == 2) then
         allocate(ZW(LQ2))
         if (fixed) then
            call ZGEQR2(LQ2, LQ2, HLPB1, LQ2, TAU, ZW, info)
         else
            allocate(VN1(LQ2), VN2(LQ2))
            do jj = 1, LQ2
               IPVT(jj) = jj
               VN1(jj) = DZNRM2(LQ2, HLPB1(1,jj), 1)
               VN2(jj) = VN1(jj)
            enddo
            call ZLAQP2(LQ2, LQ2, 0, HLPB1, LQ2, IPVT, TAU, VN1, VN2, ZW)
            deallocate(VN1, VN2)
            info = 0
         endif
         deallocate(ZW)
      else
         allocate(TT(nb, LQ2), ZW(nb*LQ2))
         call ZGEQRT(LQ2, LQ2, nb, HLPB1, LQ2, TT, nb, ZW, info)
         deallocate(ZW)
      endif
      if (info /= 0) stop 'nul_qr info'
      do ii = 1, LQ2
         X = ABS(HLPB1(ii, ii))
         D3(ii) = X
         do jj = ii, LQ2
            HLPB1(ii, jj) = HLPB1(ii, jj) / X
         enddo
      enddo
    end subroutine nul_qr
    Subroutine nul_solve(UCT, VINV)
      Complex (Kind=Kind(0.d0)), intent(in) :: UCT(LQ,LQ), VINV(LQ,LQ)
      Complex (Kind=Kind(0.d0)), allocatable :: TMPVEC(:), ZW(:)
      Complex (Kind=Kind(0.d0)) :: zz
      Integer :: info, ii, jj
      if (mode == 0) then
         call solve_extended_System(HLPB2, UCT, VINV, HLPB1, D3, TAU, IPVT, LQ, WORK, LWORK)
         return
      endif
      zz = 0.D0
      ALLOCATE(TMPVEC(LQ2))
      TMPVEC = conjg(1.D0/D3)
      call zlaset('A', LQ2, LQ2, zz, zz, HLPB2, LQ2)
      DO ii = 1, LQ
         DO jj = 1, LQ
            HLPB2(ii   , jj    ) =  UCT(ii, jj)
            HLPB2(ii+LQ, jj+LQ ) =  VINV(ii,jj)
         ENDDO
      ENDDO
      zz = 1.D0
      CALL ZLAPMR(.true., LQ2, LQ2, HLPB2(1, 1), LQ2, IPVT(1))
      CALL ZTRSM('L', 'U', 'C', 'N', LQ2, LQ2, zz, HLPB1(1, 1), LQ2, HLPB2(1, 1), LQ2)
      DO jj = 1,LQ2
         DO ii = 1,LQ2
            HLPB2(ii, jj)  = TMPVEC(ii)*HLPB2(ii, jj)
         ENDDO
      ENDDO
      deallocate(TMPVEC)
      if (mode == 2) then
         allocate(ZW(LQ2))
         call ZUNM2R('L', 'N', LQ2, LQ2, LQ2, HLPB1, LQ2, TAU, HLPB2, LQ2, ZW, info)
      else
         allocate(ZW(nb*LQ2))
         call ZGEMQRT('L', 'N', LQ2, LQ2, LQ2, nb, HLPB1, LQ2, TT, nb, HLPB2, LQ2, ZW, info)
      endif
      deallocate(ZW)
      if (info /= 0) stop 'nul_solve info'
    end subroutine nul_solve
  end subroutine nul_cgr22
end module q5nul
#endif

#if Q5_PART == 3
program t207_d5_qref
  use iso_fortran_env, only: int64, output_unit
  use q5
  use q5nul
  use UDV_State_mod
  use cgr2_2_mod, only: CGR2_2
  implicit none
  integer :: hdr(13), u, ios, n, nv, k, qs, nfl, i, hostin, nnul, ib, nbr
  integer(int64) :: p0, fsz
  complex(dp), allocatable :: Bm(:,:), Bi(:,:)
  logical :: hasb
  integer, allocatable :: masks(:), nul(:)
  character(len=1024) :: fn, arg
  character(len=8000) :: line
  complex(dp), allocatable :: M(:,:,:), Dz(:,:)
  complex(qp), allocatable :: Gq(:,:,:), Gq2(:,:,:), Gqh(:,:,:)
  real(dp), allocatable :: Lr(:,:)
  type(UDV_State) :: st(3), sp
  real(dp) :: e(5), t0, t1, tq
  integer(int64) :: seed, rate, c0, c1
  logical :: islog, hasnvh, hasdg

  call get_command_argument(1, arg); read(arg,*) qs
  call get_command_argument(2, arg)
  nnul = 0
  allocate(nul(20))
  if (trim(arg) /= '-') then
     i = 1
     do
        k = index(arg(i:), ',')
        nnul = nnul + 1
        if (k == 0) then
           read(arg(i:),*) nul(nnul); exit
        endif
        read(arg(i:i+k-2),*) nul(nnul)
        i = i + k
     enddo
  endif
  call get_command_argument(3, arg); read(arg,*) hostin
#if defined(STABLOG)
  islog = .true.
#else
  islog = .false.
#endif
  ! the rule of nul_code_ok (gpu_tau_mod): flags <= 31 without 16 here; flag 8 (compact WY, no pivoting) default path only, with
  ! 1 <= nb and without flag 2; otherwise nb = 0 (an invalid code would leave IPVT unset: ZLAPMR does not return)
  do k = 1, nnul
     i = mod(nul(k), 100); ib = mod(nul(k)/100, 1000)
     if (nul(k) < 0 .or. nul(k) >= 10000000 .or. i > 15) stop 'invalid null code'
     if (iand(i, 8) /= 0) then
        if (islog .or. iand(i, 2) /= 0 .or. ib < 1) stop 'invalid null code (flag 8: default path only, nb >= 1, not with 2)'
     else if (ib /= 0) then
        stop 'invalid null code (nb without flag 8)'
     endif
  enddo
  nfl = 0
  do
     read(*,'(A)',iostat=ios) fn
     if (ios /= 0) exit
     open(newunit=u, file=trim(fn), access='stream', form='unformatted', status='old', action='read')
     read(u) hdr(1:12)
     if (hdr(1) /= 20791 .or. (hdr(2) /= 1 .and. hdr(2) /= 2)) stop 'bad dump header'
     hdr(13) = 0
     if (hdr(2) == 2) read(u) hdr(13)
     hasb = hdr(13) == 1
     if ((hdr(9) == 1) .neqv. islog) stop 'dump path (default/LOG) does not match this build'
     n = hdr(8); nv = hdr(10); hasnvh = hdr(11) == 1; hasdg = hdr(12) == 1
     if (allocated(masks)) deallocate(masks)
     allocate(masks(max(nv,1)))
     if (nv > 0) read(u) masks(1:nv)
     if (nfl == 0) then
        do i = 1, 3
           call st(i)%alloc(n)
        enddo
        call sp%alloc(n)
        st(1)%side = 'r'; st(2)%side = 'l'; st(3)%side = 'r'; sp%side = 'r'
        allocate(M(n,n,4), Dz(n,3), Lr(n,3), Gq(n,n,4), Gq2(n,n,4), Gqh(n,n,4), Bm(n,n), Bi(n,n))
        Dz = (0.d0, 0.d0); Lr = 0.d0
     endif
     nfl = nfl + 1
     if (hasb) then                       ! B, Bi are the last 2 n^2 complex numbers of the file
        inquire(unit=u, pos=p0, size=fsz)
        read(u, pos=fsz - 32_int64*n*n + 1) Bm, Bi
        read(u, pos=p0)
     endif
     call rstate(u, st(1), 1); call rstate(u, st(2), 2)
     ! the branch rule of CGR2_2 on these inputs
#if defined(STABLOG)
     nbr = 2; if (st(2)%L(1) > st(1)%L(1)) nbr = 1
#else
     nbr = 2; if (dble(st(2)%D(1)) > dble(st(1)%D(1))) nbr = 1
#endif
     call system_clock(c0, rate)
     call qcgr22(n, islog, nbr, st(1)%U, Dz(:,1), Lr(:,1), st(1)%V, st(2)%U, Dz(:,2), Lr(:,2), st(2)%V, 0_int64, Gq)
     call system_clock(c1); tq = dble(c1 - c0)/dble(rate)
     write(line,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F7.2)') 'QREF file ', hdr(3), ' pass ', hdr(4), ' nst ', hdr(5), ' nfe ', &
          & hdr(6), ' br ', hdr(7), ' hbr ', nbr, ' tq ', tq
     read(u) M
     call put('dev')
     do k = 1, nv
        read(u) M
        write(arg,'(A,I0)') 'var', masks(k); call put(trim(arg))
     enddo
     if (hasnvh) then
        read(u) M
        call put('nvh')
     endif
     call CGR2_2(M(:,:,2), M(:,:,1), M(:,:,3), M(:,:,4), st(1), st(2), n)
     call put('cpu')
     do k = 1, nnul
        seed = mod(1000003_int64*int(hdr(3), int64) + 7919_int64*int(nul(k), int64) + 2654435_int64, 2147483648_int64)
        call nul_cgr22(n, M(:,:,2), M(:,:,1), M(:,:,3), M(:,:,4), st(1), st(2), nul(k), seed)
        write(arg,'(A,I0)') 'n', nul(k); call put(trim(arg))
     enddo
     ! 1-ulp perturbed U of the right state (perturb_state of gpu_tau_mod), against the unperturbed Gq: input sensitivity
     do k = 1, 2
        seed = mod(20261030_int64 + 104729_int64*int(hdr(3), int64) + 7_int64*k, 2147483648_int64)
        call pstate(st(1), sp, seed)
        call CGR2_2(M(:,:,2), M(:,:,1), M(:,:,3), M(:,:,4), sp, st(2), n)
        write(arg,'(A,I0)') 'u', k; call put(trim(arg))
     enddo
     if (qs > 0) then
        if (mod(hdr(3) - 1, qs) == 0) then
           call qcgr22(n, islog, nbr, st(1)%U, Dz(:,1), Lr(:,1), st(1)%V, st(2)%U, Dz(:,2), Lr(:,2), st(2)%V, &
                & int(hdr(3), int64)*7_int64 + 99991_int64, Gq2)
           call qdiff(n, Gq2, Gq, e)
           call putq('qself', e)
           ! the other branch in quad (same matrix, rearranged): mathematically identical
           call qcgr22(n, islog, 3 - nbr, st(1)%U, Dz(:,1), Lr(:,1), st(1)%V, st(2)%U, Dz(:,2), Lr(:,2), st(2)%V, 0_int64, Gq2)
           call qdiff(n, Gq2, Gq, e)
           call putq('qbr', e)
        endif
     endif
     if (hasdg) then
        call rstate(u, st(3), 3)
        if (hostin == 1) then
           call qcgr22(n, islog, nbr, st(3)%U, Dz(:,3), Lr(:,3), st(3)%V, st(2)%U, Dz(:,2), Lr(:,2), st(2)%V, 0_int64, Gqh)
           call putqq('din', Gq, Gqh)               ! quad G of the device input vs quad G of the host chain's input (input part)
           call CGR2_2(M(:,:,2), M(:,:,1), M(:,:,3), M(:,:,4), st(3), st(2), n)
           call putqq('cpuh', cmplx(M, kind=qp), Gqh)   ! host CGR2_2 on the host input against its own quad reference
           call putqq('cpuhx', cmplx(M, kind=qp), Gq)   ! host CGR2_2 on the host input against the quad G of the device input
        endif
     endif
     close(u)
     write(output_unit,'(A)') trim(line)
     flush(output_unit)
  enddo

contains

  subroutine put(name)
    character(len=*), intent(in) :: name
    real(dp) :: ep(3)
    complex(dp), allocatable :: Dd(:,:)
    call errs(n, M, Gq, e)
    ep = -1.d0
    if (hasb) then
       allocate(Dd(n,n))
       Dd = cmplx(cmplx(M(:,:,3), kind=qp) - Gq(:,:,3), kind=dp)
       ep(1) = maxval(abs(matmul(matmul(Bm, Dd), Bi)))
       Dd = cmplx(cmplx(M(:,:,2), kind=qp) - Gq(:,:,2), kind=dp)
       ep(2) = maxval(abs(matmul(Bm, Dd)))
       Dd = cmplx(cmplx(M(:,:,4), kind=qp) - Gq(:,:,4), kind=dp)
       ep(3) = maxval(abs(matmul(Dd, Bi)))
    endif
    write(line(len_trim(line)+1:),'(A,A,8(1X,ES10.3))') ' | ', name, e, ep
  end subroutine put

  ! as put for a quad A against the quad reference R (errors of A - R; propagated numbers with B, Bi when present)
  subroutine putqq(name, A, R)
    character(len=*), intent(in) :: name
    complex(qp), intent(in) :: A(n,n,4), R(n,n,4)
    real(dp) :: ep(3)
    complex(dp), allocatable :: Dd(:,:)
    call qdiff(n, A, R, e)
    ep = -1.d0
    if (hasb) then
       allocate(Dd(n,n))
       Dd = cmplx(A(:,:,3) - R(:,:,3), kind=dp)
       ep(1) = maxval(abs(matmul(matmul(Bm, Dd), Bi)))
       Dd = cmplx(A(:,:,2) - R(:,:,2), kind=dp)
       ep(2) = maxval(abs(matmul(Bm, Dd)))
       Dd = cmplx(A(:,:,4) - R(:,:,4), kind=dp)
       ep(3) = maxval(abs(matmul(Dd, Bi)))
    endif
    write(line(len_trim(line)+1:),'(A,A,8(1X,ES10.3))') ' | ', name, e, ep
  end subroutine putqq

  subroutine putq(name, e)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: e(5)
    write(line(len_trim(line)+1:),'(A,A,5(1X,ES10.3))') ' | ', name, e
  end subroutine putq

  subroutine rstate(u, s, k)
    integer, intent(in) :: u, k
    type(UDV_State), intent(inout) :: s
#if defined(STABLOG)
    read(u) s%U, s%L, s%V
    Lr(:,k) = s%L
#else
    read(u) s%U, s%D, s%V
    Dz(:,k) = s%D
#endif
  end subroutine rstate

  subroutine pstate(a, p, seed)
    type(UDV_State), intent(in) :: a
    type(UDV_State), intent(inout) :: p
    integer(int64), intent(inout) :: seed
    integer :: i, j
    p%U = a%U; p%V = a%V
#if defined(STABLOG)
    p%L = a%L
#else
    p%D = a%D
#endif
    do j = 1, size(p%U,2)
       do i = 1, size(p%U,1)
          p%U(i,j) = p%U(i,j)*(1.d0 + 4.d-16*rnd(seed))
       enddo
    enddo
  end subroutine pstate
end program t207_d5_qref
#endif
