! p8 micro benchmark: old corr kernels (verbatim) vs candidate mappings, bitwise comparison + timing, one tile.
! build: see ../jobs/b1_build.sh.  env: OB_L (24), OB_NORB (1), OB_NREP (20)
module bm
  use omp_lib
  implicit none
  integer, parameter :: dp = kind(0.d0)
  integer :: N, Norb, NFL = 2, NSUN = 1
  integer, allocatable :: o_inv(:,:), o_jc(:,:)
  complex(dp), allocatable :: G_d(:,:,:), acc(:,:,:,:,:), acc0(:,:)
  complex(dp), allocatable :: P1(:,:), P2(:,:), P3(:,:)
  complex(dp), allocatable :: T1(:,:), T2(:,:), T3(:,:), Z0(:,:)      ! term arrays (I, m) / (m, I)
  integer, parameter :: QB = QBV, NP = 2**QB, NBK = 10-QB, NK = 2**NBK
  complex(dp), allocatable :: Q1(:,:), Q2(:,:), Q3(:,:), Q4(:,:), Q5(:,:), GT(:,:,:)
  complex(dp) :: zq
  integer :: chm = 0
  integer, parameter :: CH = 16
  integer, allocatable :: o_jct(:,:)
  complex(dp), allocatable :: ZD(:,:), ZS(:,:)
  real(dp) :: wr = 0.7d0, wi = -0.3d0, sq3 = 1.7320508075688772d0
contains
  subroutine setup(L, Norb_)
    integer, intent(in) :: L, Norb_
    integer :: I, J, x, y, xj, yj, k, no
    real(dp) :: r1, r2
    N = L*L; Norb = Norb_
    allocate(o_inv(N,Norb), o_jc(N,N), G_d(N*Norb,N*Norb,NFL), acc(N,1,Norb,Norb,5), acc0(Norb,5))
    do no = 1, Norb; do I = 1, N; o_inv(I,no) = I + (no-1)*N; enddo; enddo
    do I = 1, N
       x = mod(I-1,L); y = (I-1)/L
       do J = 1, N
          xj = mod(J-1,L); yj = (J-1)/L
          k = mod(xj-x+L,L) + 1 + L*mod(yj-y+L,L)
          o_jc(I,k) = J
       enddo
    enddo
    call random_seed()
    do k = 1, NFL
       do J = 1, N*Norb; do I = 1, N*Norb
          call random_number(r1); call random_number(r2)
          G_d(I,J,k) = cmplx(r1-0.5d0, r2-0.5d0, dp)*exp(-0.3d0*abs(real(mod(I-J+N*Norb,N*Norb))))
       enddo; enddo
    enddo
    acc = 0; acc0 = 0; zq = (0.d0, 0.d0)
    !$omp target enter data map(to: o_inv, o_jc, G_d, acc, acc0)
  end subroutine setup

  subroutine zero_acc()
    acc = 0; acc0 = 0
    !$omp target update to(acc, acc0)
  end subroutine zero_acc

  ! ---------------- OLD (verbatim from gpu_obs_mod corr) ----------------
  subroutine old_green(k)
    integer, intent(in) :: k
    integer :: m, imj, no, no2, I, J, i1, j1, nf
    complex(dp) :: s1, z, a
    real(dp) :: ns
    ns = dble(NSUN)
    !$omp target teams distribute private(imj,no,no2,s1)
    do m = 0, N*Norb*Norb - 1
       imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
       s1 = (0.d0, 0.d0)
       !$omp parallel do private(I,J,i1,j1,nf,z,a) reduction(+:s1)
       do I = 1, N
          J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
          z = (0.d0, 0.d0)
          do nf = 1, NFL
             a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
             z = z + a
          enddo
          s1 = s1 + z*ns
       enddo
       !$omp end parallel do
       acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
    enddo
  end subroutine old_green

  subroutine old_den(k)
    integer, intent(in) :: k
    integer :: m, imj, no, no2, I, J, i1, j1, nf
    complex(dp) :: s1, z, zi, zj, a, zs0
    real(dp) :: ns
    ns = dble(NSUN)
    !$omp target teams distribute private(imj,no,no2,s1)
    do m = 0, N*Norb*Norb - 1
       imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
       s1 = (0.d0, 0.d0)
       !$omp parallel do private(I,J,i1,j1,nf,z,zi,zj,a) reduction(+:s1)
       do I = 1, N
          J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
          zi = (0.d0, 0.d0); zj = (0.d0, 0.d0); z = (0.d0, 0.d0)
          do nf = 1, NFL
             zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
             zj = zj + ((1.d0, 0.d0) - G_d(j1,j1,nf))
             a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
             z = z + a*G_d(i1,j1,nf)
          enddo
          s1 = s1 + (zi*ns)*(zj*ns) + z*ns
       enddo
       !$omp end parallel do
       acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
    enddo
    !$omp target teams distribute parallel do private(I,i1,zi,nf,zs0)
    do no = 1, Norb
       zs0 = (0.d0, 0.d0)
       do I = 1, N
          i1 = o_inv(I,no); zi = (0.d0, 0.d0)
          do nf = 1, NFL
             zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
          enddo
          zs0 = zs0 + zi*ns
       enddo
       acc0(no,k) = acc0(no,k) + zs0*cmplx(wr,wi,kind(0.d0))
    enddo
  end subroutine old_den

  subroutine old_spin(k)
    integer, intent(in) :: k
    integer :: m, imj, no, no2, I, J, i1, j1
    complex(dp) :: s1, s2, s3, zxy, zz, a, b, zi, zj, zs0
    !$omp target teams distribute private(imj,no,no2,s1,s2,s3)
    do m = 0, N*Norb*Norb - 1
       imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
       s1 = (0.d0, 0.d0); s2 = (0.d0, 0.d0); s3 = (0.d0, 0.d0)
       !$omp parallel do private(I,J,i1,j1,zxy,zz,a,b,zi,zj) reduction(+:s1,s2,s3)
       do I = 1, N
          J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
          a = -G_d(j1,i1,1); if (i1 == j1) a = a + (1.d0, 0.d0)
          b = -G_d(j1,i1,2); if (i1 == j1) b = b + (1.d0, 0.d0)
          zxy = a*G_d(i1,j1,2) + b*G_d(i1,j1,1)
          zi = ((1.d0, 0.d0) - G_d(i1,i1,2)) - ((1.d0, 0.d0) - G_d(i1,i1,1))
          zj = ((1.d0, 0.d0) - G_d(j1,j1,2)) - ((1.d0, 0.d0) - G_d(j1,j1,1))
          zz = a*G_d(i1,j1,1) + b*G_d(i1,j1,2) + zi*zj
          s1 = s1 + zz
          s2 = s2 + zxy
          s3 = s3 + (2.d0*zxy + zz)/3.d0
       enddo
       !$omp end parallel do
       acc(imj,1,no,no2,k)   = acc(imj,1,no,no2,k)   + s1*cmplx(wr,wi,kind(0.d0))
       acc(imj,1,no,no2,k+1) = acc(imj,1,no,no2,k+1) + s2*cmplx(wr,wi,kind(0.d0))
       acc(imj,1,no,no2,k+2) = acc(imj,1,no,no2,k+2) + s3*cmplx(wr,wi,kind(0.d0))
    enddo
    !$omp target teams distribute parallel do private(I,i1,zi,zs0)
    do no = 1, Norb
       zs0 = (0.d0, 0.d0)
       do I = 1, N
          i1 = o_inv(I,no)
          zi = ((1.d0, 0.d0) - G_d(i1,i1,2)) - ((1.d0, 0.d0) - G_d(i1,i1,1))
          zs0 = zs0 + zi
       enddo
       acc0(no,k)   = acc0(no,k)   + zs0*cmplx(wr,wi,kind(0.d0))
       acc0(no,k+2) = acc0(no,k+2) + zs0*cmplx(wr,wi,kind(0.d0))/sq3
    enddo
  end subroutine old_spin

  ! ---------------- candidate: the old kernel with thread_limit / num_teams variations (bitwise question) ----------------
  subroutine old_green_tl(k, nth, nteam)
    integer, intent(in) :: k, nth, nteam
    integer :: m, imj, no, no2, I, J, i1, j1, nf
    complex(dp) :: s1, z, a
    real(dp) :: ns
    ns = dble(NSUN)
    !$omp target teams distribute private(imj,no,no2,s1) thread_limit(nth) num_teams(nteam)
    do m = 0, N*Norb*Norb - 1
       imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
       s1 = (0.d0, 0.d0)
       !$omp parallel do private(I,J,i1,j1,nf,z,a) reduction(+:s1)
       do I = 1, N
          J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
          z = (0.d0, 0.d0)
          do nf = 1, NFL
             a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
             z = z + a
          enddo
          s1 = s1 + z*ns
       enddo
       !$omp end parallel do
       acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
    enddo
  end subroutine old_green_tl

  ! ---------------- candidate TS: terms kernel (coalesced over m) + same-tree reduction kernel ----------------
  subroutine alloc_T()
    allocate(T1(N*Norb*Norb, N), T2(N*Norb*Norb, N), T3(N*Norb*Norb, N), Z0(N, N*Norb*Norb))
    !$omp target enter data map(alloc: T1, T2, T3, Z0)
  end subroutine alloc_T

  subroutine ts_green(k, layout)
    integer, intent(in) :: k, layout
    integer :: m, imj, no, no2, I, J, i1, j1, nf, M0
    complex(dp) :: s1, z, a
    real(dp) :: ns
    ns = dble(NSUN); M0 = N*Norb*Norb
    !$omp target teams distribute parallel do collapse(2) private(imj,no,no2,J,i1,j1,nf,z,a)
    do I = 1, N
       do m = 0, M0 - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
          z = (0.d0, 0.d0)
          do nf = 1, NFL
             a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
             z = z + a
          enddo
          T1(m+1,I) = z
       enddo
    enddo
    !$omp target teams distribute private(imj,no,no2,s1)
    do m = 0, M0 - 1
       imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
       s1 = (0.d0, 0.d0)
       !$omp parallel do private(I) reduction(+:s1)
       do I = 1, N
          s1 = s1 + T1(m+1,I)*ns
       enddo
       !$omp end parallel do
       acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
    enddo
  end subroutine ts_green

  ! ---------------- candidate SER: one thread per output, serial sum over I (order 1..N), NOT the old tree order ----------------
  subroutine ser_green(k)
    integer, intent(in) :: k
    integer :: m, imj, no, no2, I, J, i1, j1, nf
    complex(dp) :: s1, z, a
    real(dp) :: ns
    ns = dble(NSUN)
    !$omp target teams distribute parallel do private(imj,no,no2,s1,I,J,i1,j1,nf,z,a)
    do m = 0, N*Norb*Norb - 1
       imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
       s1 = (0.d0, 0.d0)
       do I = 1, N
          J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
          z = (0.d0, 0.d0)
          do nf = 1, NFL
             a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
             z = z + a
          enddo
          s1 = s1 + z*ns
       enddo
       acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
    enddo
  end subroutine ser_green

  ! ---------------- 1x1 kernels: old vs terms + serial per-output sum ----------------
  subroutine acc0_den_new(k)
    integer, intent(in) :: k
    integer :: I, no, i1, nf
    complex(dp) :: zi, zs0
    real(dp) :: ns
    ns = dble(NSUN)
    !$omp target teams distribute parallel do collapse(2) private(i1,zi,nf)
    do no = 1, Norb
       do I = 1, N
          i1 = o_inv(I,no); zi = (0.d0, 0.d0)
          do nf = 1, NFL
             zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
          enddo
          Z0(I,no) = zi
       enddo
    enddo
    !$omp target teams distribute parallel do private(I,zs0)
    do no = 1, Norb
       zs0 = (0.d0, 0.d0)
       do I = 1, N
          zs0 = zs0 + Z0(I,no)*ns
       enddo
       acc0(no,k) = acc0(no,k) + zs0*cmplx(wr,wi,kind(0.d0))
    enddo
  end subroutine acc0_den_new

  ! single kernel version: 1 team, parallel terms, then master serial sum (no second launch)
  subroutine acc0_den_new1(k)
    integer, intent(in) :: k
    integer :: I, no, i1, nf
    complex(dp) :: zi, zs0
    real(dp) :: ns
    ns = dble(NSUN)
    !$omp target teams num_teams(1) private(zs0)
    !$omp parallel do collapse(2) private(i1,zi,nf)
    do no = 1, Norb
       do I = 1, N
          i1 = o_inv(I,no); zi = (0.d0, 0.d0)
          do nf = 1, NFL
             zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
          enddo
          Z0(I,no) = zi
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(I,zs0)
    do no = 1, Norb
       zs0 = (0.d0, 0.d0)
       do I = 1, N
          zs0 = zs0 + Z0(I,no)*ns
       enddo
       acc0(no,k) = acc0(no,k) + zs0*cmplx(wr,wi,kind(0.d0))
    enddo
    !$omp end parallel do
    !$omp end target teams
  end subroutine acc0_den_new1

  ! ======================= p8 candidates: reproduce the old work-group reduction tree =======================
  ! hypothesis H: reduction(+:s) over <= 1024 work-items = halving tree on the work-item index t = I-1:
  !   A_1024(t) = x_t (0 for t >= N);  A_s(t) = A_2s(t) + A_2s(t+s);  result A_1(0)
  ! = two stages: partial over k (t = tp + 32 k, strides 512..32), then over tp (strides 16..1); each as a bit-reversal pairwise sum
  subroutine push5(st, j, v)
    !$omp declare target
    complex(dp), intent(inout) :: st(0:5)
    integer, intent(in) :: j
    complex(dp), intent(in) :: v
    complex(dp) :: w
    integer :: c, lv
    w = v; c = j; lv = 0
    do while (iand(c,1) == 1)
       w = st(lv) + w; c = ishft(c,-1); lv = lv + 1
    enddo
    st(lv) = w
  end subroutine push5
  integer function brev5(j)
    !$omp declare target
    integer, intent(in) :: j
    brev5 = ishft(iand(j,1),4) + ishft(iand(j,2),2) + iand(j,4) + ishft(iand(j,8),-2) + ishft(iand(j,16),-4)
  end function brev5

  subroutine alloc_P()
    allocate(P1(N*Norb*Norb,0:31), P2(N*Norb*Norb,0:31), P3(N*Norb*Norb,0:31))
    !$omp target enter data map(alloc: P1, P2, P3)
  end subroutine alloc_P

  subroutine tree_green(k)
    integer, intent(in) :: k
    integer :: m, imj, no, no2, I, J, i1, j1, nf, tp, M0, jj, kk
    complex(dp) :: s1, z, a, st(0:5)
    real(dp) :: ns
    ns = dble(NSUN); M0 = N*Norb*Norb
    !$omp target teams distribute parallel do collapse(2) private(imj,no,no2,I,J,i1,j1,nf,z,a,s1,st,jj,kk)
    do tp = 0, 31
       do m = 0, M0 - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          st = (0.d0, 0.d0)
          !dir$ unroll(32)
          do jj = 0, 31
             kk = brev5(jj); I = tp + 32*kk + 1
             s1 = (0.d0, 0.d0)
             if (I <= N) then
                J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
                z = (0.d0, 0.d0)
                do nf = 1, NFL
                   a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                   z = z + a
                enddo
                s1 = s1 + z*ns
             endif
             call push5(st, jj, s1)
          enddo
          P1(m+1,tp) = st(5)
       enddo
    enddo
    call stage2(k, 1)
  end subroutine tree_green

  subroutine stage2(k, which)
    integer, intent(in) :: k, which
    integer :: m, imj, no, no2, M0, j
    complex(dp) :: st(0:5), s1
    M0 = N*Norb*Norb
    !$omp target teams distribute parallel do private(imj,no,no2,st,j,s1)
    do m = 0, M0 - 1
       imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
       st = (0.d0, 0.d0)
       do j = 0, 31
          call push5(st, j, P1(m+1, brev5(j)))
       enddo
       s1 = st(5)
       acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
    enddo
  end subroutine stage2

  ! 1x1 replacement variants (b): chunked gather in one single-thread-per-output kernel; (c) terms kernel + chunk-buffered serial sum
  subroutine acc0_den_chunk(k)
    integer, intent(in) :: k
    integer :: I, I0, no, i1, nf, c
    complex(dp) :: zi, zs0, zb(0:31)
    real(dp) :: ns
    ns = dble(NSUN)
    !$omp target teams distribute parallel do private(I,I0,i1,zi,nf,zs0,zb,c)
    do no = 1, Norb
       zs0 = (0.d0, 0.d0)
       do I0 = 1, N, 32
          do c = 0, 31
             I = I0 + c
             zi = (0.d0, 0.d0)
             if (I <= N) then
                i1 = o_inv(I,no)
                do nf = 1, NFL
                   zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
                enddo
             endif
             zb(c) = zi
          enddo
          do c = 0, 31
             if (I0 + c <= N) zs0 = zs0 + zb(c)*ns
          enddo
       enddo
       acc0(no,k) = acc0(no,k) + zs0*cmplx(wr,wi,kind(0.d0))
    enddo
  end subroutine acc0_den_chunk

  subroutine old_acc0_den(k)
    integer, intent(in) :: k
    integer :: I, no, i1, nf
    complex(dp) :: zi, zs0
    real(dp) :: ns
    ns = dble(NSUN)
    !$omp target teams distribute parallel do private(I,i1,zi,nf,zs0)
    do no = 1, Norb
       zs0 = (0.d0, 0.d0)
       do I = 1, N
          i1 = o_inv(I,no); zi = (0.d0, 0.d0)
          do nf = 1, NFL
             zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
          enddo
          zs0 = zs0 + zi*ns
       enddo
       acc0(no,k) = acc0(no,k) + zs0*cmplx(wr,wi,kind(0.d0))
    enddo
  end subroutine old_acc0_den

  ! ======================= generic version (compile-time QB = log2 of the number of partial sums) =======================
  subroutine pushq(st, j, v)
    !$omp declare target
    complex(dp), intent(inout) :: st(0:10)
    integer, intent(in) :: j
    complex(dp), intent(in) :: v
    complex(dp) :: w
    integer :: c, lv
    w = v; c = j; lv = 0
    do while (iand(c,1) == 1)
       w = st(lv) + w; c = ishft(c,-1); lv = lv + 1
    enddo
    st(lv) = w
  end subroutine pushq
  integer function brevn(j, nb)
    !$omp declare target
    integer, intent(in) :: j, nb
    integer :: b
    brevn = 0
    do b = 0, nb-1
       brevn = ior(ishft(brevn,1), iand(ishft(j,-b),1))
    enddo
  end function brevn

  subroutine alloc_Q()
    allocate(o_jct(N,N)); o_jct = transpose(o_jc)
    !$omp target enter data map(to: o_jct)
    allocate(Q1(N*Norb*Norb,0:NP-1), Q2(N*Norb*Norb,0:NP-1), Q3(N*Norb*Norb,0:NP-1), Q4(N*Norb*Norb,0:NP-1), Q5(N*Norb*Norb,0:NP-1), GT(N*Norb,N*Norb,NFL))
    !$omp target enter data map(alloc: Q1, Q2, Q3, Q4, Q5, GT)
  end subroutine alloc_Q

  subroutine transpose_G()
    integer :: i, j, nf, NN
    NN = N*Norb
    !$omp target teams distribute parallel do collapse(3) private(i,j,nf)
    do nf = 1, NFL
       do j = 1, NN
          do i = 1, NN
             GT(i,j,nf) = G_d(j,i,nf)
          enddo
       enddo
    enddo
  end subroutine transpose_G

  subroutine tq_green(k)
    integer, intent(in) :: k
    integer :: m, imj, no, no2, I, J, i1, j1, nf, tp, M0, jj, kk
    complex(dp) :: s1, z, a, st(0:10)
    real(dp) :: ns
    ns = dble(NSUN); M0 = N*Norb*Norb
    !$omp target teams distribute parallel do collapse(2) private(imj,no,no2,I,J,i1,j1,nf,z,a,s1,st,jj,kk)
    do tp = 0, NP-1
       do m = 0, M0 - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          st = (0.d0, 0.d0)
          do jj = 0, NK-1
             kk = brevn(jj, NBK); I = tp + NP*kk + 1
             s1 = (0.d0, 0.d0)
             if (I <= N) then
                J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
                z = (0.d0, 0.d0)
                do nf = 1, NFL
                   a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                   z = z + a
                enddo
                s1 = s1 + z*ns
             endif
             call pushq(st, jj, s1)
          enddo
          Q1(m+1,tp) = st(NBK)
       enddo
    enddo
    call tq_stage2(k, 1)
  end subroutine tq_green

  subroutine tq_den(k, useT)
    integer, intent(in) :: k, useT
    integer :: m, imj, no, no2, I, J, i1, j1, nf, tp, M0, jj, kk
    complex(dp) :: s1, z, zi, zj, a, st(0:10)
    real(dp) :: ns
    ns = dble(NSUN); M0 = N*Norb*Norb
    !$omp target teams distribute parallel do collapse(2) private(imj,no,no2,I,J,i1,j1,nf,z,zi,zj,a,s1,st,jj,kk)
    do tp = 0, NP-1
       do m = 0, M0 - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          st = (0.d0, 0.d0)
          do jj = 0, NK-1
             kk = brevn(jj, NBK); I = tp + NP*kk + 1
             s1 = (0.d0, 0.d0)
             if (I <= N) then
                J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
                zi = (0.d0, 0.d0); zj = (0.d0, 0.d0); z = (0.d0, 0.d0)
                if (useT == 0) then
                do nf = 1, NFL
                   zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
                   zj = zj + ((1.d0, 0.d0) - G_d(j1,j1,nf))
                   a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                   z = z + a*G_d(i1,j1,nf)
                enddo
                else
                do nf = 1, NFL
                   zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
                   zj = zj + ((1.d0, 0.d0) - G_d(j1,j1,nf))
                   a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                   z = z + a*GT(j1,i1,nf)
                enddo
                endif
                s1 = s1 + (zi*ns)*(zj*ns) + z*ns
             endif
             call pushq(st, jj, s1)
          enddo
          Q1(m+1,tp) = st(NBK)
       enddo
    enddo
    call tq_stage2(k, 1)
  end subroutine tq_den

  subroutine tq_spin(k, useT)
    integer, intent(in) :: k, useT
    integer :: m, imj, no, no2, I, J, i1, j1, tp, M0, jj, kk
    complex(dp) :: s1, s2, s3, zxy, zz, a, b, zi, zj, gA, gB, st1(0:10), st2(0:10), st3(0:10)
    M0 = N*Norb*Norb
    !$omp target teams distribute parallel do collapse(2) private(imj,no,no2,I,J,i1,j1,zxy,zz,a,b,zi,zj,gA,gB,s1,s2,s3,st1,st2,st3,jj,kk)
    do tp = 0, NP-1
       do m = 0, M0 - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          st1 = (0.d0, 0.d0); st2 = (0.d0, 0.d0); st3 = (0.d0, 0.d0)
          do jj = 0, NK-1
             kk = brevn(jj, NBK); I = tp + NP*kk + 1
             s1 = (0.d0, 0.d0); s2 = (0.d0, 0.d0); s3 = (0.d0, 0.d0)
             if (I <= N) then
                J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
                a = -G_d(j1,i1,1); if (i1 == j1) a = a + (1.d0, 0.d0)
                b = -G_d(j1,i1,2); if (i1 == j1) b = b + (1.d0, 0.d0)
                if (useT == 0) then
                   gA = G_d(i1,j1,1); gB = G_d(i1,j1,2)
                else
                   gA = GT(j1,i1,1); gB = GT(j1,i1,2)
                endif
                zxy = a*gB + b*gA
                zi = ((1.d0, 0.d0) - G_d(i1,i1,2)) - ((1.d0, 0.d0) - G_d(i1,i1,1))
                zj = ((1.d0, 0.d0) - G_d(j1,j1,2)) - ((1.d0, 0.d0) - G_d(j1,j1,1))
                zz = a*gA + b*gB + zi*zj
                s1 = s1 + zz
                s2 = s2 + zxy
                s3 = s3 + (2.d0*zxy + zz)/3.d0
             endif
             call pushq(st1, jj, s1); call pushq(st2, jj, s2); call pushq(st3, jj, s3)
          enddo
          Q1(m+1,tp) = st1(NBK); Q2(m+1,tp) = st2(NBK); Q3(m+1,tp) = st3(NBK)
       enddo
    enddo
    call tq_stage2(k, 3)
  end subroutine tq_spin

  subroutine tq_stage2(k, nsum)
    integer, intent(in) :: k, nsum
    integer :: m, imj, no, no2, M0, jj
    complex(dp) :: st1(0:10), st2(0:10), st3(0:10)
    M0 = N*Norb*Norb
    !$omp target teams distribute parallel do private(imj,no,no2,st1,st2,st3,jj)
    do m = 0, M0 - 1
       imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
       st1 = (0.d0, 0.d0)
       do jj = 0, NP-1
          call pushq(st1, jj, Q1(m+1, brevn(jj,QB)))
       enddo
       acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + st1(QB)*cmplx(wr,wi,kind(0.d0))
       if (nsum == 3) then
          st2 = (0.d0, 0.d0); st3 = (0.d0, 0.d0)
          do jj = 0, NP-1
             call pushq(st2, jj, Q2(m+1, brevn(jj,QB))); call pushq(st3, jj, Q3(m+1, brevn(jj,QB)))
          enddo
          acc(imj,1,no,no2,k+1) = acc(imj,1,no,no2,k+1) + st2(QB)*cmplx(wr,wi,kind(0.d0))
          acc(imj,1,no,no2,k+2) = acc(imj,1,no,no2,k+2) + st3(QB)*cmplx(wr,wi,kind(0.d0))
       endif
    enddo
  end subroutine tq_stage2

  ! ======================= v2: leaves in a private array (static indices), halving tree in registers, opaque zero =======================
  subroutine tx_green(k)
    integer, intent(in) :: k
    integer :: m, imj, no, no2, I, Ic, J, i1, j1, nf, tp, M0, kk, lev, sst
    complex(dp) :: s1, z, a, x(0:NK-1), zq_
    real(dp) :: ns
    ns = dble(NSUN); M0 = N*Norb*Norb; zq_ = zq
    !$omp target teams distribute parallel do collapse(2) private(imj,no,no2,I,Ic,J,i1,j1,nf,z,a,s1,x,kk,lev,sst)
    do tp = 0, NP-1
       do m = 0, M0 - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          do kk = 0, NK-1
             I = tp + NP*kk + 1; Ic = min(I,N)
             J = o_jct(imj,Ic); i1 = o_inv(Ic,no); j1 = o_inv(J,no2)
             z = (0.d0, 0.d0)
             do nf = 1, NFL
                a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                z = z + a
             enddo
             s1 = zq_
             s1 = s1 + z*ns
             if (I > N) s1 = (0.d0, 0.d0)
             x(kk) = s1
          enddo
          do lev = NBK-1, 0, -1
             sst = ishft(1,lev)
             do kk = 0, sst-1
                x(kk) = x(kk) + x(kk+sst)
             enddo
          enddo
          Q1(m+1,tp) = x(0)
       enddo
    enddo
    call tq_stage2(k, 1)
  end subroutine tx_green

  subroutine tx_den(k, useT)
    integer, intent(in) :: k, useT
    integer :: m, imj, no, no2, I, Ic, J, i1, j1, nf, tp, M0, kk, lev, sst
    complex(dp) :: s1, z, zi, zj, a, x(0:NK-1), zq_, gA
    real(dp) :: ns
    ns = dble(NSUN); M0 = N*Norb*Norb; zq_ = zq
    !$omp target teams distribute parallel do collapse(2) private(imj,no,no2,I,Ic,J,i1,j1,nf,z,zi,zj,a,gA,s1,x,kk,lev,sst)
    do tp = 0, NP-1
       do m = 0, M0 - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          do kk = 0, NK-1
             I = tp + NP*kk + 1; Ic = min(I,N)
             J = o_jct(imj,Ic); i1 = o_inv(Ic,no); j1 = o_inv(J,no2)
             zi = (0.d0, 0.d0); zj = (0.d0, 0.d0); z = (0.d0, 0.d0)
             do nf = 1, NFL
                zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
                zj = zj + ((1.d0, 0.d0) - G_d(j1,j1,nf))
                a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                if (useT == 0) then
                   gA = G_d(i1,j1,nf)
                else
                   gA = GT(j1,i1,nf)
                endif
                z = z + a*gA
             enddo
             s1 = zq_
             s1 = s1 + (zi*ns)*(zj*ns) + z*ns
             if (I > N) s1 = (0.d0, 0.d0)
             x(kk) = s1
          enddo
          do lev = NBK-1, 0, -1
             sst = ishft(1,lev)
             do kk = 0, sst-1
                x(kk) = x(kk) + x(kk+sst)
             enddo
          enddo
          Q1(m+1,tp) = x(0)
       enddo
    enddo
    call tq_stage2(k, 1)
  end subroutine tx_den

  subroutine tx_spin(k, useT)
    integer, intent(in) :: k, useT
    integer :: m, imj, no, no2, I, Ic, J, i1, j1, tp, M0, kk, lev, sst
    complex(dp) :: s1, s2, s3, zxy, zz, a, b, zi, zj, gA, gB, x1(0:NK-1), x2(0:NK-1), x3(0:NK-1), zq_
    M0 = N*Norb*Norb; zq_ = zq
    !$omp target teams distribute parallel do collapse(2) private(imj,no,no2,I,Ic,J,i1,j1,zxy,zz,a,b,zi,zj,gA,gB,s1,s2,s3,x1,x2,x3,kk,lev,sst)
    do tp = 0, NP-1
       do m = 0, M0 - 1
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          do kk = 0, NK-1
             I = tp + NP*kk + 1; Ic = min(I,N)
             J = o_jct(imj,Ic); i1 = o_inv(Ic,no); j1 = o_inv(J,no2)
             a = -G_d(j1,i1,1); if (i1 == j1) a = a + (1.d0, 0.d0)
             b = -G_d(j1,i1,2); if (i1 == j1) b = b + (1.d0, 0.d0)
             if (useT == 0) then
                gA = G_d(i1,j1,1); gB = G_d(i1,j1,2)
             else
                gA = GT(j1,i1,1); gB = GT(j1,i1,2)
             endif
             zxy = a*gB + b*gA
             zi = ((1.d0, 0.d0) - G_d(i1,i1,2)) - ((1.d0, 0.d0) - G_d(i1,i1,1))
             zj = ((1.d0, 0.d0) - G_d(j1,j1,2)) - ((1.d0, 0.d0) - G_d(j1,j1,1))
             zz = a*gA + b*gB + zi*zj
             s1 = zq_; s2 = zq_; s3 = zq_
             s1 = s1 + zz
             s2 = s2 + zxy
             s3 = s3 + (2.d0*zxy + zz)/3.d0
             if (I > N) then; s1 = (0.d0, 0.d0); s2 = (0.d0, 0.d0); s3 = (0.d0, 0.d0); endif
             x1(kk) = s1; x2(kk) = s2; x3(kk) = s3
          enddo
          do lev = NBK-1, 0, -1
             sst = ishft(1,lev)
             do kk = 0, sst-1
                x1(kk) = x1(kk) + x1(kk+sst); x2(kk) = x2(kk) + x2(kk+sst); x3(kk) = x3(kk) + x3(kk+sst)
             enddo
          enddo
          Q1(m+1,tp) = x1(0); Q2(m+1,tp) = x2(0); Q3(m+1,tp) = x3(0)
       enddo
    enddo
    call tq_stage2(k, 3)
  end subroutine tx_spin

  ! ======================= fused: all kinds in one stage-1 kernel and one stage-2 kernel (incl. the acc0 serial sums) =======================
  subroutine old_sun(k)
    integer, intent(in) :: k
    integer :: m, imj, no, no2, I, J, i1, j1
    complex(dp) :: s1, a
    real(dp) :: ns
    ns = dble(NSUN)
    !$omp target teams distribute private(imj,no,no2,s1)
    do m = 0, N*Norb*Norb - 1
       imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
       s1 = (0.d0, 0.d0)
       !$omp parallel do private(I,J,i1,j1,a) reduction(+:s1)
       do I = 1, N
          J = o_jc(I,imj); i1 = o_inv(I,no); j1 = o_inv(J,no2)
          a = -G_d(j1,i1,1); if (i1 == j1) a = a + (1.d0, 0.d0)
          s1 = s1 + a*G_d(i1,j1,1)*ns
       enddo
       !$omp end parallel do
       acc(imj,1,no,no2,k) = acc(imj,1,no,no2,k) + s1*cmplx(wr,wi,kind(0.d0))
    enddo
  end subroutine old_sun

  subroutine alloc_F()
    allocate(ZD(N,Norb), ZS(N,Norb))
    !$omp target enter data map(alloc: ZD, ZS)
  end subroutine alloc_F

  !> fg: Green at k=kg (always); fd: Den at kd (+acc0); fs: SpinMz at ks (k,k+1,k+2) (+acc0); fu: SpinSUN at ku
  subroutine tf_all(kg, fd, kd, fs, ks, fu, ku)
    integer, intent(in) :: kg, kd, ks, ku
    logical, intent(in) :: fd, fs, fu
    integer :: lin, m, imj, no, no2, I, Ic, J, i1, j1, nf, tp, M0, kk, lev, sst, d, NL, NFLl, jj
    complex(dp) :: st1(0:10), st2(0:10), st3(0:10), st4(0:10), st5(0:10), zs0, zb(0:CH-1)
    integer :: I0, c, chm_l
    complex(dp) :: sg, sd, s1, s2, s3, z, zg, zi, zj, a, b, gA, gB, zxy, zz, zq_
    complex(dp) :: x1(0:NK-1), x2(0:NK-1), x3(0:NK-1), x4(0:NK-1), x5(0:NK-1)
    real(dp) :: ns
    ns = dble(NSUN); M0 = N*Norb*Norb; zq_ = zq; NL = NP*M0; NFLl = NFL; chm_l = chm
    !$omp target teams distribute parallel do private(m,imj,no,no2,I,Ic,J,i1,j1,nf,tp,kk,lev,sst,d,sg,sd,s1,s2,s3,z,zg,zi,zj,a,b,gA,gB,zxy,zz,x1,x2,x3,x4,x5)
    do lin = 0, NL + N*Norb - 1
       if (lin < NL) then
          tp = lin / M0; m = lin - tp*M0
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          do kk = 0, NK-1
             I = tp + NP*kk + 1; Ic = min(I,N)
             J = o_jct(imj,Ic); i1 = o_inv(Ic,no); j1 = o_inv(J,no2)
             ! shared pieces
             zg = (0.d0, 0.d0)
             do nf = 1, NFLl
                a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                zg = zg + a
             enddo
             sg = zq_; sg = sg + zg*ns
             if (I > N) sg = (0.d0, 0.d0)
             x1(kk) = sg
             if (fd) then
                zi = (0.d0, 0.d0); zj = (0.d0, 0.d0); z = (0.d0, 0.d0)
                do nf = 1, NFLl
                   zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
                   zj = zj + ((1.d0, 0.d0) - G_d(j1,j1,nf))
                   a = -G_d(j1,i1,nf); if (i1 == j1) a = a + (1.d0, 0.d0)
                   z = z + a*G_d(i1,j1,nf)
                enddo
                sd = zq_
                sd = sd + (zi*ns)*(zj*ns) + z*ns
                if (I > N) sd = (0.d0, 0.d0)
                x2(kk) = sd
             endif
             if (fs) then
                a = -G_d(j1,i1,1); if (i1 == j1) a = a + (1.d0, 0.d0)
                b = -G_d(j1,i1,2); if (i1 == j1) b = b + (1.d0, 0.d0)
                gA = G_d(i1,j1,1); gB = G_d(i1,j1,2)
                zxy = a*gB + b*gA
                zi = ((1.d0, 0.d0) - G_d(i1,i1,2)) - ((1.d0, 0.d0) - G_d(i1,i1,1))
                zj = ((1.d0, 0.d0) - G_d(j1,j1,2)) - ((1.d0, 0.d0) - G_d(j1,j1,1))
                zz = a*gA + b*gB + zi*zj
                s1 = zq_; s2 = zq_; s3 = zq_
                s1 = s1 + zz
                s2 = s2 + zxy
                s3 = s3 + (2.d0*zxy + zz)/3.d0
                if (I > N) then; s1 = (0.d0, 0.d0); s2 = (0.d0, 0.d0); s3 = (0.d0, 0.d0); endif
                x3(kk) = s1; x4(kk) = s2; x5(kk) = s3
             elseif (fu) then
                a = -G_d(j1,i1,1); if (i1 == j1) a = a + (1.d0, 0.d0)
                s1 = zq_
                s1 = s1 + a*G_d(i1,j1,1)*ns
                if (I > N) s1 = (0.d0, 0.d0)
                x3(kk) = s1
             endif
          enddo
          do lev = NBK-1, 0, -1
             sst = ishft(1,lev)
             do kk = 0, sst-1
                x1(kk) = x1(kk) + x1(kk+sst)
                if (fd) x2(kk) = x2(kk) + x2(kk+sst)
                if (fs .or. fu) x3(kk) = x3(kk) + x3(kk+sst)
                if (fs) then; x4(kk) = x4(kk) + x4(kk+sst); x5(kk) = x5(kk) + x5(kk+sst); endif
             enddo
          enddo
          Q1(m+1,tp) = x1(0)
          if (fd) Q2(m+1,tp) = x2(0)
          if (fs .or. fu) Q3(m+1,tp) = x3(0)
          if (fs) then; Q4(m+1,tp) = x4(0); Q5(m+1,tp) = x5(0); endif
       else
          d = lin - NL; no = d / N + 1; I = d - (no-1)*N + 1
          i1 = o_inv(I,no)
          if (fd) then
             zi = (0.d0, 0.d0)
             do nf = 1, NFLl
                zi = zi + ((1.d0, 0.d0) - G_d(i1,i1,nf))
             enddo
             ZD(I,no) = zi
          endif
          if (fs) ZS(I,no) = ((1.d0, 0.d0) - G_d(i1,i1,2)) - ((1.d0, 0.d0) - G_d(i1,i1,1))
       endif
    enddo
    ! ---- stage 2: trees over the NP partials; the extra Norb threads do the serial acc0 sums (order I = 1..N)
    !$omp target teams distribute parallel do private(m,imj,no,no2,jj,st1,st2,st3,st4,st5,zs0,I,zi,d,I0,c,zb)
    do lin = 0, M0 + 2*Norb - 1
       if (lin < M0) then
          m = lin
          imj = mod(m,N) + 1; no = mod(m/N,Norb) + 1; no2 = m/(N*Norb) + 1
          call tree_sum(Q1, m, st1)
          acc(imj,1,no,no2,kg) = acc(imj,1,no,no2,kg) + st1(QB)*cmplx(wr,wi,kind(0.d0))
          if (fd) then
             call tree_sum(Q2, m, st1)
             acc(imj,1,no,no2,kd) = acc(imj,1,no,no2,kd) + st1(QB)*cmplx(wr,wi,kind(0.d0))
          endif
          if (fu) then
             call tree_sum(Q3, m, st1)
             acc(imj,1,no,no2,ku) = acc(imj,1,no,no2,ku) + st1(QB)*cmplx(wr,wi,kind(0.d0))
          endif
          if (fs) then
             call tree_sum(Q3, m, st1)
             acc(imj,1,no,no2,ks) = acc(imj,1,no,no2,ks) + st1(QB)*cmplx(wr,wi,kind(0.d0))
             call tree_sum(Q4, m, st1)
             acc(imj,1,no,no2,ks+1) = acc(imj,1,no,no2,ks+1) + st1(QB)*cmplx(wr,wi,kind(0.d0))
             call tree_sum(Q5, m, st1)
             acc(imj,1,no,no2,ks+2) = acc(imj,1,no,no2,ks+2) + st1(QB)*cmplx(wr,wi,kind(0.d0))
          endif
       else
          d = lin - M0
          if (d < Norb) then
             no = d + 1
             if (fd) then
                zs0 = (0.d0, 0.d0)
                if (chm_l == 0) then
                   do I = 1, N
                      zs0 = zs0 + ZD(I,no)*ns
                   enddo
                else
                   do I0 = 1, N, CH
                      do c = 0, CH-1
                         zb(c) = ZD(min(I0+c,N),no)
                      enddo
                      do c = 0, CH-1
                         if (I0+c <= N) zs0 = zs0 + zb(c)*ns
                      enddo
                   enddo
                endif
                acc0(no,kd) = acc0(no,kd) + zs0*cmplx(wr,wi,kind(0.d0))
             endif
          else
             no = d - Norb + 1
             if (fs) then
                zs0 = (0.d0, 0.d0)
                if (chm_l == 0) then
                   do I = 1, N
                      zs0 = zs0 + ZS(I,no)
                   enddo
                else
                   do I0 = 1, N, CH
                      do c = 0, CH-1
                         zb(c) = ZS(min(I0+c,N),no)
                      enddo
                      do c = 0, CH-1
                         if (I0+c <= N) zs0 = zs0 + zb(c)
                      enddo
                   enddo
                endif
                acc0(no,ks)   = acc0(no,ks)   + zs0*cmplx(wr,wi,kind(0.d0))
                acc0(no,ks+2) = acc0(no,ks+2) + zs0*cmplx(wr,wi,kind(0.d0))/sq3
             endif
          endif
       endif
    enddo
  end subroutine tf_all

  subroutine tree_sum(Q, m, st)
    !$omp declare target
    complex(dp), intent(in) :: Q(:,0:)
    integer, intent(in) :: m
    complex(dp), intent(inout) :: st(0:10)
    integer :: jj
    st = (0.d0, 0.d0)
    do jj = 0, NP-1
       call pushq(st, jj, Q(m+1, brevn(jj,QB)))
    enddo
  end subroutine tree_sum
end module bm

program obsbench
  use bm
  implicit none
  integer :: L, Nb, nrep, r, i, j
  character(len=32) :: s
  real(dp) :: tm0, tm1
  complex(dp), allocatable :: ref(:,:,:,:,:), ref0(:,:), cur(:,:,:,:,:), cur0(:,:)
  integer :: nbad

  L = 24; Nb = 1; nrep = 20
  call get_environment_variable('OB_L', s); if (len_trim(s) > 0) read(s,*) L
  call get_environment_variable('OB_NORB', s); if (len_trim(s) > 0) read(s,*) Nb
  call get_environment_variable('OB_NREP', s); if (len_trim(s) > 0) read(s,*) nrep
  call setup(L, Nb); call alloc_T()
  allocate(ref(N,1,Norb,Norb,5), ref0(Norb,5), cur(N,1,Norb,Norb,5), cur0(Norb,5))
  print '(A,I0,A,I0,A,I0,A,I0)', 'OB N=', N, ' Norb=', Norb, ' NFL=', NFL, ' devices=', omp_get_num_devices()

  ! reference: old kernels twice (reproducibility), keep first
  call zero_acc(); call old_green(1); call old_spin(2); call old_den(5)
  !$omp target update from(acc, acc0)
  ref = acc; ref0 = acc0
  call zero_acc(); call old_green(1); call old_spin(2); call old_den(5)
  !$omp target update from(acc, acc0)
  print '(A,L1)', 'OB old run-to-run bitwise: ', all(acc == ref) .and. all(acc0 == ref0)



  call alloc_Q(); call alloc_F()
  call zero_acc()
  call old_green(1)
  call old_den(5)
  call old_spin(2)
  !$omp target update from(acc, acc0)
  ref = acc; ref0 = acc0
  call timeit('old_all', 17)
  do chm = 0, 1
     call zero_acc(); call tf_all(1, .true., 5, .true., 2, .false., 0)
     !$omp target update from(acc, acc0)
     print '(A,I0,A,L1,A,L1)', 'OB chunk mode ', chm, ' acc0 bitwise: ', all(acc0 == ref0), '  acc bitwise: ', all(acc == ref)
     call timeit('tf_all', 18)
  enddo
contains
  subroutine report(name, k)
    character(len=*), intent(in) :: name
    integer, intent(in) :: k
    integer :: nb2
    real(dp) :: dmax
    nb2 = count(acc(:,:,:,:,k) /= ref(:,:,:,:,k))
    dmax = maxval(abs(acc(:,:,:,:,k) - ref(:,:,:,:,k)))
    print '(A,A,A,I0,A,ES10.2,A,I0)', 'OB ', name, ' vs old: n differing = ', nb2, ' max|d| = ', dmax, ' of ', size(acc(:,:,:,:,k))
  end subroutine report
  subroutine chk_tl(nth, nteam)
    integer, intent(in) :: nth, nteam
    call zero_acc(); call old_green_tl(1, nth, nteam)
    !$omp target update from(acc, acc0)
    print '(A,I0,A,I0)', 'OB old_green thread_limit ', nth, ' num_teams ', nteam
    call report('  tl-variant', 1)
    tm0 = omp_get_wtime()
    do r = 1, nrep; call old_green_tl(1, nth, nteam); enddo
    tm1 = omp_get_wtime()
    print '(A,F10.1,A)', 'OB   time/call ', 1.d6*(tm1-tm0)/nrep, ' us'
  end subroutine chk_tl
  subroutine timeit(name, which)
    character(len=*), intent(in) :: name
    integer, intent(in) :: which
    call zero_acc()
    do r = 1, 3; call run(which); enddo
    tm0 = omp_get_wtime()
    do r = 1, nrep; call run(which); enddo
    tm1 = omp_get_wtime()
    print '(A,A,A,F10.1,A)', 'OB TIME ', name, ' ', 1.d6*(tm1-tm0)/nrep, ' us/call'
  end subroutine timeit
  subroutine run(which)
    integer, intent(in) :: which
    select case (which)
    case (1); call old_green(1)
    case (2); call old_spin(2)
    case (3); call old_den(5)
    case (4); call ts_green(1, 1)
    case (5); call ser_green(1)
    case (6); call acc0_den_new(5)
    case (7); call acc0_den_new1(5)
    case (8); call tree_green(1)
    case (9); call old_acc0_den(5)
    case (10); call acc0_den_chunk(5)
    case (17); call old_green(1); call old_den(5); call old_spin(2)
    case (18); call tf_all(1, .true., 5, .true., 2, .false., 0)
    case (11); call tx_green(1)
    case (12); call tx_den(5,0)
    case (13); call tx_den(5,1)
    case (14); call transpose_G()
    case (15); call tx_spin(2,0)
    case (16); call tx_spin(2,1)
    end select
  end subroutine run
end program obsbench
