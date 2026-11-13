! p9/devupd feasibility bench: can a device (ifx OpenMP offload) evaluation of the per-site update arithmetic of the GPU
! submatrix block path (blk_gpp, site_update ratio/decision/phase, gpu_lr_accept M update) be bitwise equal to the host?
! Same source expressions on host and device (ALF flags), random inputs, count of bit differences per operation.
! Part 2: time of a one-team device kernel doing the sequential Metropolis of one block (b sites, 2 flavours) vs host.
! env: DU_NS (samples, 200000), DU_K (pending sites in the chain tests, 24), DU_NSUN (2), DU_B (64), DU_NREP (200)
module dub
  use omp_lib
  implicit none
  integer, parameter :: dp = kind(0.d0)
contains
  integer function env_i(name, def)
    character(len=*), intent(in) :: name
    integer, intent(in) :: def
    character(len=32) :: v
    integer :: st
    env_i = def
    call get_environment_variable(name, v, status=st)
    if (st == 0 .and. len_trim(v) > 0) read(v,*) env_i
  end function env_i

  subroutine rnd(z, n, sc, sh)
    integer, intent(in) :: n
    complex(dp), intent(out) :: z(n)
    real(dp), intent(in) :: sc, sh
    real(dp) :: r1(n), r2(n)
    call random_number(r1); call random_number(r2)
    z = cmplx(sc*(r1-0.5d0) + sh, sc*(r2-0.5d0), dp)
  end subroutine rnd

  ! bit differences and max relative difference of two complex arrays
  subroutine cmpz(name, a, b, n)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    complex(dp), intent(in) :: a(n), b(n)
    integer :: i, nd
    real(dp) :: mx
    nd = 0; mx = 0.d0
    do i = 1, n
       if (transfer(real(a(i)), 0_8) /= transfer(real(b(i)), 0_8) .or. &
           transfer(aimag(a(i)), 0_8) /= transfer(aimag(b(i)), 0_8)) then
          nd = nd + 1
          mx = max(mx, abs(a(i)-b(i))/max(abs(a(i)), 1.d-300))
       endif
    enddo
    write(*,'(A,A24,A,I9,A,I9,A,ES10.3)') ' BITCMP ', name, ' differ ', nd, ' of ', n, '  max rel ', mx
  end subroutine cmpz
end module dub

program dubench
  use dub
  implicit none
  integer :: ns, k, nsun, b, nrep, i, a, m, c, it, fl, kk, rep
  complex(dp), allocatable :: z(:), d(:), v(:), ph(:), g0(:), myexp(:,:), ea(:,:), gam(:), zz(:,:)
  complex(dp), allocatable :: Mm(:,:,:), Gc(:,:), Gr(:,:), av(:,:), Mn(:,:,:)
  complex(dp), allocatable :: oh(:,:), od(:,:)
  real(dp), allocatable :: wh(:), wd(:), s0(:)
  complex(dp) :: s, gpp(2), R(2), Rt, Dm, zx
  real(dp) :: t0, th, td
  ! block kernel data
  complex(dp), allocatable :: Gb(:,:,:), Mb(:,:,:), Wb(:,:), gw(:,:), zb(:,:)
  real(dp), allocatable :: rn(:)
  integer, allocatable :: Sb(:,:), kb(:), accb(:)
  complex(dp) :: ex(2), e2(2)

  ns = env_i('DU_NS', 200000); k = env_i('DU_K', 24); nsun = env_i('DU_NSUN', 2)
  b = env_i('DU_B', 64); nrep = env_i('DU_NREP', 200)
  write(*,*) 'samples', ns, ' k', k, ' N_SUN', nsun, ' block', b, ' nrep', nrep
  allocate(z(ns), d(ns), v(ns), ph(ns), g0(ns), myexp(ns,2), ea(ns,2), gam(ns), zz(ns,2), s0(ns))
  allocate(oh(ns,8), od(ns,8), wh(ns), wd(ns))
  call random_seed()
  call rnd(z, ns, 1.0d0, 1.0d0); call rnd(d, ns, 1.0d0, -0.3d0); call rnd(v, ns, 1.0d0, 0.5d0)
  call rnd(ph, ns, 2.0d0, 0.2d0); call rnd(g0, ns, 1.0d0, 0.5d0); call rnd(gam, ns, 0.2d0, 1.0d0)
  call rnd(myexp(:,1), ns, 0.5d0, 1.0d0); call rnd(myexp(:,2), ns, 0.5d0, 1.0d0)
  call rnd(ea(:,1), ns, 0.2d0, 1.0d0); call rnd(ea(:,2), ns, 0.2d0, 1.0d0)
  call random_number(s0)
  zz = myexp - 1.d0
  allocate(Mm(k,k,ns), Gc(k,ns), Gr(k,ns), av(k,ns), Mn(k,k,ns))
  do i = 1, ns
     do a = 1, k
        call rnd(Mm(:,a,i), k, 0.4d0, 0.d0)
     enddo
     call rnd(Gc(:,i), k, 0.4d0, 0.d0); call rnd(Gr(:,i), k, 0.4d0, 0.d0); call rnd(av(:,i), k, 1.0d0, 0.d0)
  enddo

  ! ---- host ----
  do i = 1, ns
     oh(i,1) = z(i)**dble(nsun)                                   ! cpow (site_update: Ratiotot**dble(N_SUN))
     oh(i,2) = ph(i) * z(i)/sqrt(z(i)*conjg(z(i)))                ! csqrt (Phase update)
     oh(i,3) = d(i)/(1.d0 + v(i)*d(i))                            ! complex division (xv)
     oh(i,4) = exp(z(i)*v(i))                                     ! cexp (myexp, exp(Z1*alpha); precomputable on the host)
  enddo
  do i = 1, ns                                                    ! blk_gpp chain: W = -M Gcol, gpp = g0 - Grow . W
     s = (0.d0,0.d0)
     do a = 1, k
        zx = (0.d0,0.d0)
        do m = 1, k
           zx = zx - Mm(a,m,i) * Gc(m,i)
        enddo
        s = s + Gr(a,i) * zx
     enddo
     oh(i,5) = g0(i) - s
  enddo
  do i = 1, ns                                                    ! ratio + weight (site_update), gpp given
     gpp(1) = g0(i); gpp(2) = v(i)
     do fl = 1, 2
        Dm = - zz(i,fl) * gpp(fl)
        Dm = myexp(i,fl) + Dm
        R(fl) = Dm * ea(i,fl)
     enddo
     Rt = product(R)
     Rt = (Rt**dble(nsun)) * gam(i)/gam(ns+1-i)
     oh(i,6) = Rt
     wh(i) = s0(i) * abs( real(ph(i)*Rt, dp)/real(ph(i), dp) )
  enddo
  do i = 1, ns                                                    ! M update of gpu_lr_accept (only element (k,k) and bt(1) kept)
     zx = (0.d0,0.d0)
     do a = 1, k
        zx = zx + Gr(a,i) * Mm(a,1,i)
     enddo
     oh(i,7) = zx
     oh(i,8) = Mm(k,k,i) + d(i) * av(k,i) * Gc(k,i)
  enddo

  ! ---- device (same expressions) ----
  !$omp target data map(to: z, d, v, ph, g0, myexp, ea, gam, zz, s0, Mm, Gc, Gr, av) map(from: od, wd)
  !$omp target teams distribute parallel do
  do i = 1, ns
     od(i,1) = z(i)**dble(nsun)
     od(i,2) = ph(i) * z(i)/sqrt(z(i)*conjg(z(i)))
     od(i,3) = d(i)/(1.d0 + v(i)*d(i))
     od(i,4) = exp(z(i)*v(i))
  enddo
  !$omp target teams distribute parallel do private(s, zx, a, m)
  do i = 1, ns
     s = (0.d0,0.d0)
     do a = 1, k
        zx = (0.d0,0.d0)
        do m = 1, k
           zx = zx - Mm(a,m,i) * Gc(m,i)
        enddo
        s = s + Gr(a,i) * zx
     enddo
     od(i,5) = g0(i) - s
  enddo
  !$omp target teams distribute parallel do private(gpp, R, Rt, Dm, fl)
  do i = 1, ns
     gpp(1) = g0(i); gpp(2) = v(i)
     do fl = 1, 2
        Dm = - zz(i,fl) * gpp(fl)
        Dm = myexp(i,fl) + Dm
        R(fl) = Dm * ea(i,fl)
     enddo
     Rt = product(R)
     Rt = (Rt**dble(nsun)) * gam(i)/gam(ns+1-i)
     od(i,6) = Rt
     wd(i) = s0(i) * abs( real(ph(i)*Rt, dp)/real(ph(i), dp) )
  enddo
  !$omp target teams distribute parallel do private(zx, a)
  do i = 1, ns
     zx = (0.d0,0.d0)
     do a = 1, k
        zx = zx + Gr(a,i) * Mm(a,1,i)
     enddo
     od(i,7) = zx
     od(i,8) = Mm(k,k,i) + d(i) * av(k,i) * Gc(k,i)
  enddo
  !$omp end target data

  call cmpz('cpow z**dble(N_SUN)', oh(:,1), od(:,1), ns)
  call cmpz('phase csqrt', oh(:,2), od(:,2), ns)
  call cmpz('cdiv xv', oh(:,3), od(:,3), ns)
  call cmpz('cexp', oh(:,4), od(:,4), ns)
  call cmpz('gpp chain (blk_gpp)', oh(:,5), od(:,5), ns)
  call cmpz('Ratiotot', oh(:,6), od(:,6), ns)
  call cmpz('bt dot (accept)', oh(:,7), od(:,7), ns)
  call cmpz('M rank-1 elem', oh(:,8), od(:,8), ns)
  write(*,'(A,I9,A,ES10.3)') ' BITCMP weight differ ', count(transfer(wh,0_8,ns) /= transfer(wd,0_8,ns)), '  max rel ', &
       maxval(abs(wh-wd)/max(wh,1.d-300))

  ! ---- part 2: one block of b sites, sequential Metropolis, 2 flavours, on one team vs host ----
  allocate(Gb(b,b,2), Mb(b,b,2), Wb(b,2), gw(b,2), zb(b,2), rn(b), Sb(b,2), kb(2), accb(b))
  do fl = 1, 2
     do a = 1, b
        call rnd(Gb(:,a,fl), b, 0.2d0, 0.d0)
        Gb(a,a,fl) = Gb(a,a,fl) + 0.5d0
     enddo
  enddo
  call random_number(rn)
  call rnd(ex, 2, 0.2d0, 1.0d0); e2 = ex - 1.d0
  ! host
  t0 = omp_get_wtime()
  do rep = 1, nrep
     call blk_host()
  enddo
  th = (omp_get_wtime() - t0)/nrep
  write(*,'(A,I4,A,F10.2,A,I4,A,I4)') ' BLK host  b=', b, '  us/block ', th*1.d6, '  accepted ', count(accb == 1), '  k ', kb(1)
  !$omp target enter data map(to: Gb, rn, ex, e2) map(alloc: Mb, Wb, gw, Sb, kb, accb)
  call blk_dev()                                   ! warm-up (JIT / first launch)
  t0 = omp_get_wtime()
  do rep = 1, nrep
     call blk_dev()
  enddo
  td = (omp_get_wtime() - t0)/nrep
  !$omp target update from(accb, kb)
  write(*,'(A,I4,A,F10.2,A,I4,A,I4)') ' BLK dev1  b=', b, '  us/block ', td*1.d6, '  accepted ', count(accb == 1), '  k ', kb(1)

contains
  ! simplified block Metropolis (submatrix M form, 2 flavours): per site W = M (e_i - G(S,p)), gpp, ratio, decision, M update
  subroutine blk_host()
    integer :: it, fl, a, m, kk
    complex(dp) :: s, gp(2), R2, xv, bt(b)
    Mb = (0.d0,0.d0); kb = 0; accb = 0
    do it = 1, b
       do fl = 1, 2
          kk = kb(fl)
          do a = 1, kk
             s = (0.d0,0.d0)
             do m = 1, kk
                s = s - Mb(a,m,fl) * Gb(m,it,fl)
             enddo
             Wb(a,fl) = s; gw(a,fl) = Gb(it,a,fl)
          enddo
          s = (0.d0,0.d0)
          do a = 1, kk
             s = s + gw(a,fl) * Wb(a,fl)
          enddo
          gp(fl) = Gb(it,it,fl) - s
       enddo
       R2 = ((ex(1) - e2(1)*gp(1)) * (ex(2) - e2(2)*gp(2)))**2
       if (abs(real(R2)) > rn(it)) then
          accb(it) = 1
          do fl = 1, 2
             kk = kb(fl) + 1
             xv = e2(fl)/(1.d0 + (1.d0 - gp(fl))*e2(fl))
             do m = 1, kk
                s = (0.d0,0.d0)
                do a = 1, kk - 1
                   s = s + gw(a,fl) * Mb(a,m,fl)
                enddo
                bt(m) = s
             enddo
             bt(kk) = bt(kk) + 1.d0
             Wb(kk,fl) = (0.d0,0.d0)
             do m = 1, kk
                do a = 1, kk
                   Mb(a,m,fl) = Mb(a,m,fl) + xv * (1.d0 - Wb(a,fl)) * bt(m)
                enddo
             enddo
             kb(fl) = kk
          enddo
       endif
    enddo
  end subroutine blk_host

  subroutine blk_dev()
    integer :: it, fl, a, m, kk, acc
    complex(dp) :: s, gp(2), R2, xv(2), bt(b,2)
    !$omp target teams num_teams(1) thread_limit(256) private(it, fl, a, m, kk, s, gp, R2, xv, acc, bt)
    !$omp parallel num_threads(256) private(it, fl, a, m, kk, s) shared(gp, R2, xv, acc, bt)
    !$omp do collapse(3)
    do fl = 1, 2
       do m = 1, b
          do a = 1, b
             Mb(a,m,fl) = (0.d0,0.d0)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp single
    kb = 0; accb = 0
    !$omp end single
    do it = 1, b
       !$omp do collapse(2)
       do fl = 1, 2
          do a = 1, b
             if (a <= kb(fl)) then
                s = (0.d0,0.d0)
                do m = 1, kb(fl)
                   s = s - Mb(a,m,fl) * Gb(m,it,fl)
                enddo
                Wb(a,fl) = s; gw(a,fl) = Gb(it,a,fl)
             endif
          enddo
       enddo
       !$omp end do
       !$omp single
       do fl = 1, 2
          s = (0.d0,0.d0)
          do a = 1, kb(fl)
             s = s + gw(a,fl) * Wb(a,fl)
          enddo
          gp(fl) = Gb(it,it,fl) - s
       enddo
       R2 = ((ex(1) - e2(1)*gp(1)) * (ex(2) - e2(2)*gp(2)))**2
       acc = 0
       if (abs(real(R2)) > rn(it)) acc = 1
       accb(it) = acc
       do fl = 1, 2
          xv(fl) = e2(fl)/(1.d0 + (1.d0 - gp(fl))*e2(fl))
       enddo
       !$omp end single
       if (acc == 1) then
          !$omp do collapse(2)
          do fl = 1, 2
             do m = 1, b
                if (m <= kb(fl) + 1) then
                   s = (0.d0,0.d0)
                   do a = 1, kb(fl)
                      s = s + gw(a,fl) * Mb(a,m,fl)
                   enddo
                   if (m == kb(fl) + 1) s = s + 1.d0
                   bt(m,fl) = s
                endif
             enddo
          enddo
          !$omp end do
          !$omp do collapse(3)
          do fl = 1, 2
             do m = 1, b
                do a = 1, b
                   if (a <= kb(fl) + 1 .and. m <= kb(fl) + 1) then
                      if (a <= kb(fl)) then
                         Mb(a,m,fl) = Mb(a,m,fl) + xv(fl) * (1.d0 - Wb(a,fl)) * bt(m,fl)
                      else
                         Mb(a,m,fl) = Mb(a,m,fl) + xv(fl) * bt(m,fl)
                      endif
                   endif
                enddo
             enddo
          enddo
          !$omp end do
          !$omp single
          kb = kb + 1
          !$omp end single
       endif
    enddo
    !$omp end parallel
    !$omp end target teams
  end subroutine blk_dev
end program dubench
