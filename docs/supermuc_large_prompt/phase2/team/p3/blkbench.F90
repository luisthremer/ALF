include "mkl_blas_omp_offload.f90"
! Phase 3b step 2: time per submatrix block (blk_begin + blk_end device sequence of gpu_lowrank_mod, Ndim=576, nd=b=32, 2 flavours)
! (a) OpenMP, synchronous, exactly the calls of ALF (target regions + dispatch zgemm + target updates + la_sync-like wait)
! (b) SYCL on the interop queue, every op followed by a wait, (c) SYCL, one wait per block (the D2H of the block gather).
program blkbench
  use omp_lib, only: omp_interop_kind, omp_interop_none, omp_get_default_device, omp_ifr_sycl, omp_get_wtime, omp_set_default_device
  use iso_c_binding
  use onemkl_blas_omp_offload_lp64, only: zgemm
  implicit none
  integer, parameter :: ND = 32, NDIM = 576, NFL = 2, LDC = 2*ND + ND*ND
  integer(c_int), parameter :: omp_ipr_fr_id = -1, omp_ipr_targetsync = -9
  interface
    integer(c_intptr_t) function omp_get_interop_int(interop, property_id, ret_code) bind(C, name='omp_get_interop_int')
      import :: c_intptr_t, c_int
      integer(c_intptr_t), value :: interop
      integer(c_int), value :: property_id
      integer(c_int), intent(out) :: ret_code
    end function
    type(c_ptr) function omp_get_interop_ptr(interop, property_id, ret_code) bind(C, name='omp_get_interop_ptr')
      import :: c_intptr_t, c_int, c_ptr
      integer(c_intptr_t), value :: interop
      integer(c_int), value :: property_id
      integer(c_int), intent(out) :: ret_code
    end function
    integer(c_int) function bb_init(q) bind(C, name='bb_init')
      import :: c_ptr, c_int
      type(c_ptr), value :: q
    end function
    type(c_ptr) function bb_alloc_host(nbytes) bind(C, name='bb_alloc_host')
      import :: c_ptr, c_long
      integer(c_long), value :: nbytes
    end function
    subroutine bb_wait() bind(C, name='bb_wait')
    end subroutine
    subroutine bb_memcpy(dst, src, nbytes) bind(C, name='bb_memcpy')
      import :: c_ptr, c_long
      type(c_ptr), value :: dst, src
      integer(c_long), value :: nbytes
    end subroutine
    subroutine bb_gather(nfl, ld, b, nd, gbk, bsk, Gk) bind(C, name='bb_gather')
      import :: c_ptr, c_int
      integer(c_int), value :: nfl, ld, b, nd
      type(c_ptr), value :: gbk, bsk, Gk
    end subroutine
    subroutine bb_scale(nfl, ld, b, nd, ldc, Gk, bsk, cbk) bind(C, name='bb_scale')
      import :: c_ptr, c_int
      integer(c_int), value :: nfl, ld, b, nd, ldc
      type(c_ptr), value :: Gk, bsk, cbk
    end subroutine
    subroutine bb_fgather(G, ld, nrow, S, k, nf, C, Qt, slot) bind(C, name='bb_fgather')
      import :: c_ptr, c_int
      integer(c_int), value :: ld, nrow, k, nf
      type(c_ptr), value :: G, S, C, Qt, slot
    end subroutine
    subroutine bb_gemm(ta, tb, m, n, k, ar, ai, A, lda, B, ldb, br, bi, C, ldc) bind(C, name='bb_gemm')
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ta, tb, m, n, k, lda, ldb, ldc
      real(c_double), value :: ar, ai, br, bi
      type(c_ptr), value :: A, B, C
    end subroutine
  end interface
  complex(8), parameter :: cone = (1.d0,0.d0), cmone = (-1.d0,0.d0), czero = (0.d0,0.d0)
  complex(8), allocatable, target :: G_d(:,:,:), X_d(:,:,:), scr_c(:,:,:), gbd(:,:,:), cbuf(:,:), G0(:,:,:)
  integer, allocatable, target :: bsd(:,:), bS_d(:,:), slot_d(:,:), nsig_d(:,:)
  integer(omp_interop_kind) :: iop
  type(c_ptr) :: qh
  integer(c_intptr_t) :: fr
  integer :: nblk, rc, i, it, inord
  real(8) :: t0, t1, t2, chk
  real(8) :: topt(0:20) = 0.d0, tlast = 0.d0
  character(len=32) :: arg
  ! device pointers (use_device_addr)
  type(c_ptr) :: pG, pX, pS, pGb, pBs, pBS2, pCb, pSl, pNs
  nblk = 2000
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read(arg,*) nblk; end if
  call omp_set_default_device(0)
  allocate(G_d(NDIM,NDIM,NFL), G0(NDIM,NDIM,NFL), X_d(NDIM,ND,NFL), scr_c(NDIM,NDIM,2), gbd(ND,ND,2), cbuf(LDC,2))
  allocate(bsd(ND,2), bS_d(ND,2), slot_d(NDIM,2), nsig_d(NDIM,2))
  call init_data()
  !$omp target enter data map(to: G_d, X_d, scr_c, gbd, cbuf, bsd, bS_d, slot_d, nsig_d)
  iop = omp_interop_none
  !$omp interop init(prefer_type("sycl"), targetsync: iop) device(omp_get_default_device())
  fr = omp_get_interop_int(iop, omp_ipr_fr_id, rc)
  qh = c_null_ptr
  if (rc == 0 .and. fr == omp_ifr_sycl) qh = omp_get_interop_ptr(iop, omp_ipr_targetsync, rc)
  if (.not. c_associated(qh)) stop 'no sycl queue'
  inord = bb_init(qh)
  !$omp target data use_device_addr(G_d, X_d, scr_c, gbd, cbuf, bsd, bS_d, slot_d, nsig_d)
  pG = c_loc(G_d); pX = c_loc(X_d); pS = c_loc(scr_c); pGb = c_loc(gbd); pCb = c_loc(cbuf); pBs = c_loc(bsd); pBS2 = c_loc(bS_d)
  pSl = c_loc(slot_d); pNs = c_loc(nsig_d)
  !$omp end target data
  write(*,'(A,I0,A,I0)') 'in-order queue=', inord, '  blocks per run=', nblk
  do it = 1, 3
    call fresh(); t0 = omp_get_wtime()
    do i = 1, nblk; call block_omp(i); end do
    call bb_wait(); t1 = omp_get_wtime(); call rep('OMP sync (ALF calls)', it, t0, t1, t1)
  end do
  do it = 1, 3
    call fresh(); t0 = omp_get_wtime()
    do i = 1, nblk; call block_sycl(i, .true.); end do
    call bb_wait(); t1 = omp_get_wtime(); call rep('SYCL wait after every op', it, t0, t1, t1)
    write(*,'(A,10F8.1)') '   us/block per op: 1 H2D bsd, 2 gather, 3 H2D cbuf, 4 H2D bS, 5 scale, 6 fgather x2, 7 gemm(576x32x32) x2, 8 gemm(576x576x32) x2, 9 H2D nsig, 10 D2H gbd: ', 1.d6*topt(1:10)/nblk
    topt = 0.d0
  end do
  do it = 1, 3
    call fresh(); t0 = omp_get_wtime()
    do i = 1, nblk; call block_sycl(i, .false.); end do
    t2 = omp_get_wtime(); call bb_wait(); t1 = omp_get_wtime(); call rep('SYCL async, 1 wait per block', it, t0, t1, t2)
  end do
  call experiments()
  !$omp interop destroy(iop)
contains
  subroutine experiments()
    type(c_ptr) :: ph
    complex(8), pointer :: gh(:,:,:)
    integer :: it2, i2
    real(8) :: s0, s1
    ph = bb_alloc_host(int(16*ND*ND*2, c_long))
    call c_f_pointer(ph, gh, [ND, ND, 2])
    call fresh()
    call host_fill(1)
    !$omp target update to(bsd)
    do it2 = 1, 2
      s0 = omp_get_wtime()
      do i2 = 1, nblk
        call k_blk_gather(2, 1, 2, NDIM, ND, ND, 2, gbd, bsd, G_d)
        !$omp target update from(gbd)
      end do
      s1 = omp_get_wtime(); write(*,'(A,F8.1,A)') 'E1 OMP gather kernel + update from(gbd) 32 KB       : ', 1.d6*(s1-s0)/nblk, ' us'
      s0 = omp_get_wtime()
      do i2 = 1, nblk
        call k_blk_gather(2, 1, 2, NDIM, ND, ND, 2, gbd, bsd, G_d)
      end do
      s1 = omp_get_wtime(); write(*,'(A,F8.1,A)') 'E0 OMP gather kernel only                          : ', 1.d6*(s1-s0)/nblk, ' us'
      s0 = omp_get_wtime()
      do i2 = 1, nblk
        !$omp target update from(gbd(1:2,1,1))
      end do
      s1 = omp_get_wtime(); write(*,'(A,F8.1,A)') 'E0b OMP update from 32 B                           : ', 1.d6*(s1-s0)/nblk, ' us'
      s0 = omp_get_wtime()
      do i2 = 1, nblk
        !$omp target update from(gbd(1:ND,1:ND,1:1))
      end do
      s1 = omp_get_wtime(); write(*,'(A,F8.1,A)') 'E0c OMP update from 16 KB                          : ', 1.d6*(s1-s0)/nblk, ' us'
      s0 = omp_get_wtime()
      do i2 = 1, nblk
        call bb_gather(2, NDIM, ND, ND, pGb, pBs, pG)
        call bb_memcpy(c_loc(gbd), pGb, int(16*ND*ND*2, c_long)); call bb_wait()
      end do
      s1 = omp_get_wtime(); write(*,'(A,F8.1,A)') 'E2 SYCL gather + memcpy to pageable + wait          : ', 1.d6*(s1-s0)/nblk, ' us'
      s0 = omp_get_wtime()
      do i2 = 1, nblk
        call bb_gather(2, NDIM, ND, ND, pGb, pBs, pG)
        call bb_memcpy(ph, pGb, int(16*ND*ND*2, c_long)); call bb_wait()
      end do
      s1 = omp_get_wtime(); write(*,'(A,F8.1,A)') 'E3 SYCL gather + memcpy to host-USM + wait          : ', 1.d6*(s1-s0)/nblk, ' us'
      s0 = omp_get_wtime()
      do i2 = 1, nblk
        call bb_gather(2, NDIM, ND, ND, ph, pBs, pG); call bb_wait()
      end do
      s1 = omp_get_wtime(); write(*,'(A,F8.1,A)') 'E4 SYCL gather kernel writing host-USM directly + wait: ', 1.d6*(s1-s0)/nblk, ' us'
      s0 = omp_get_wtime()
      do i2 = 1, nblk
        call bb_gather(2, NDIM, ND, ND, pGb, pBs, pG); call bb_wait()
      end do
      s1 = omp_get_wtime(); write(*,'(A,F8.1,A)') 'E5 SYCL gather kernel to device + wait (no copy)     : ', 1.d6*(s1-s0)/nblk, ' us'
    end do
  end subroutine
  subroutine init_data()
    integer :: i, j, f
    do f = 1, NFL
      do j = 1, NDIM; do i = 1, NDIM
        G0(i,j,f) = cmplx(0.01d0*sin(0.37d0*i + 0.11d0*j + f), 0.01d0*cos(0.23d0*i - 0.17d0*j), 8)
      end do; end do
      do i = 1, NDIM; G0(i,i,f) = G0(i,i,f) + cmplx(0.5d0, 0.d0, 8); end do
    end do
    X_d = czero; scr_c = czero; gbd = czero; cbuf = czero; slot_d = 1; nsig_d = 1
    bsd = 1; bS_d = 1
  end subroutine
  subroutine fresh()
    G_d = G0
    !$omp target update to(G_d)
  end subroutine
  subroutine rep(lab, it, ta, tb, ti)
    character(*), intent(in) :: lab
    integer, intent(in) :: it
    real(8), intent(in) :: ta, tb, ti
    !$omp target update from(G_d, gbd)
    chk = sum(real(G_d)) + sum(aimag(G_d)) + sum(real(gbd))
    write(*,'(A,A,I1,A,F8.1,A,F8.1,A,ES24.16)') lab, ' #', it, ': ', 1.d6*(tb-ta)/nblk, ' us/block  issue-only ', 1.d6*(ti-ta)/nblk, ' us  chk ', chk
  end subroutine
  subroutine tick(op)
    integer, intent(in) :: op
    real(8) :: tn
    call bb_wait()
    tn = omp_get_wtime()
    topt(op) = topt(op) + (tn - tlast)
    tlast = tn
  end subroutine
  subroutine host_fill(blk)
    integer, intent(in) :: blk
    integer :: m, f
    do f = 1, 2
      do m = 1, ND
        bsd(m,f) = mod((m-1)*17 + blk*5 + 100*(f-1), NDIM) + 1
        bS_d(m,f) = mod((m-1)*13 + blk*7 + 50*(f-1), NDIM) + 1
      end do
      cbuf(1:ND,f) = cmplx(1.d0 + 1.d-4*mod(blk,7), 1.d-4, 8)
      cbuf(ND+1:2*ND,f) = cmplx(1.d0 - 1.d-4*mod(blk,5), -1.d-4, 8)
      do m = 1, ND*ND
        cbuf(2*ND+m,f) = cmplx(2.d-4*sin(0.1d0*m + blk), 1.d-4*cos(0.2d0*m), 8)
      end do
    end do
  end subroutine
  subroutine block_omp(blk)
    integer, intent(in) :: blk
    integer :: fl
    call host_fill(blk)
    call bb_wait()                       ! la_sync (alfsh_wait)
    !$omp target update to(bsd)
    call k_blk_gather(2, 1, 2, NDIM, ND, ND, 2, gbd, bsd, G_d)
    !$omp target update from(gbd)
    call bb_wait()                       ! la_sync at blk_end
    !$omp target update to(cbuf, bS_d)
    call k_blk_scale(2, 1, 2, NDIM, ND, ND, 2, LDC, G_d, bsd, cbuf)
    do fl = 1, 2
      call k_sub_flush_gather(G_d(1,1,fl), NDIM, NDIM, bS_d(1,fl), ND, fl, X_d(1,1,fl), scr_c(1,1,1), slot_d, 2)
      !$omp dispatch
      call zgemm('N','T',NDIM,ND,ND,cone,scr_c(1,1,1),NDIM,cbuf(2*ND+1,fl),ND,czero,scr_c(1,1,2),NDIM)
      !$omp dispatch
      call zgemm('N','T',NDIM,NDIM,ND,cmone,X_d(1,1,fl),NDIM,scr_c(1,1,2),NDIM,cone,G_d(1,1,fl),NDIM)
    end do
    !$omp target update to(nsig_d(1:ND,1))
  end subroutine
  subroutine block_sycl(blk, wait_each)
    integer, intent(in) :: blk
    logical, intent(in) :: wait_each
    integer :: fl
    call host_fill(blk)
    call bb_wait()
    tlast = omp_get_wtime()
    call bb_memcpy(pBs, c_loc(bsd), int(4*ND*2, c_long)); if (wait_each) call tick(1)
    call bb_gather(2, NDIM, ND, ND, pGb, pBs, pG); if (wait_each) call tick(2)
    call bb_memcpy(c_loc(gbd), pGb, int(16*ND*ND*2, c_long))
    if (wait_each) then; call tick(10); else; call bb_wait(); end if     ! the one unavoidable sync (host needs the block of G)
    call bb_memcpy(pCb, c_loc(cbuf), int(16*LDC*2, c_long)); if (wait_each) call tick(3)
    call bb_memcpy(pBS2, c_loc(bS_d), int(4*ND*2, c_long)); if (wait_each) call tick(4)
    call bb_scale(2, NDIM, ND, ND, LDC, pG, pBs, pCb); if (wait_each) call tick(5)
    do fl = 1, 2
      call bb_fgather(off(pG, 16_c_long*NDIM*NDIM*(fl-1)), NDIM, NDIM, off(pBS2, 4_c_long*ND*(fl-1)), ND, fl, &
           off(pX, 16_c_long*NDIM*ND*(fl-1)), pS, pSl); if (wait_each) call tick(6)
      call bb_gemm(0, 1, NDIM, ND, ND, 1.d0, 0.d0, pS, NDIM, off(pCb, 16_c_long*(2*ND + LDC*(fl-1))), ND, 0.d0, 0.d0, &
           off(pS, 16_c_long*NDIM*NDIM), NDIM); if (wait_each) call tick(7)
      call bb_gemm(0, 1, NDIM, NDIM, ND, -1.d0, 0.d0, off(pX, 16_c_long*NDIM*ND*(fl-1)), NDIM, off(pS, 16_c_long*NDIM*NDIM), NDIM, &
           1.d0, 0.d0, off(pG, 16_c_long*NDIM*NDIM*(fl-1)), NDIM); if (wait_each) call tick(8)
    end do
    call bb_memcpy(pNs, c_loc(nsig_d), int(4*ND, c_long)); if (wait_each) call tick(9)
  end subroutine
  function off(p, nbytes) result(q)
    type(c_ptr), intent(in) :: p
    integer(c_long), intent(in) :: nbytes
    type(c_ptr) :: q
    q = transfer(transfer(p, 0_c_intptr_t) + nbytes, q)
  end function
  ! ---- the three ALF kernels as in gpu_lowrank_mod (explicit-shape dummies, map(present,alloc:)) ----
  Subroutine k_blk_gather(nfl, nfa, nfb, ld, b, nd_, nfg, gbk, bsk, Gk)
    Integer, intent(in) :: nfl, nfa, nfb, ld, b, nd_, nfg
    Complex (Kind=Kind(0.d0)), intent(inout) :: gbk(nd_,nd_,2), Gk(ld,ld,nfg)
    Integer, intent(in) :: bsk(nd_,2)
    Integer :: fl, a, c, nf
    !$omp target teams distribute parallel do collapse(3) private(nf) map(present,alloc: gbk, bsk, Gk)
    do fl = 1, nfl
       do c = 1, b
          do a = 1, b
             nf = nfa
             if (fl == 2) nf = nfb
             gbk(a,c,fl) = Gk(bsk(a,fl),bsk(c,fl),nf)
          enddo
       enddo
    enddo
  end subroutine
  Subroutine k_blk_scale(nfl, nfa, nfb, ld, b, nd_, nfg, ldc_, Gk, bsk, cbk)
    Integer, intent(in) :: nfl, nfa, nfb, ld, b, nd_, nfg, ldc_
    Complex (Kind=Kind(0.d0)), intent(inout) :: Gk(ld,ld,nfg)
    Integer, intent(in) :: bsk(nd_,2)
    Complex (Kind=Kind(0.d0)), intent(in) :: cbk(ldc_,2)
    Integer :: fl, a, j, c, cj, nf, sa
    Complex (Kind=Kind(0.d0)) :: v
    !$omp target teams distribute parallel do collapse(3) private(nf, sa, c, cj, v) map(present,alloc: Gk, bsk, cbk)
    do fl = 1, nfl
       do a = 1, b
          do j = 1, ld
             nf = nfa
             if (fl == 2) nf = nfb
             sa = bsk(a,fl)
             cj = 0
             do c = 1, b
                if (bsk(c,fl) == j) cj = c
             enddo
             v = cbk(a,fl) * Gk(sa,j,nf)
             if (cj > 0) v = cbk(nd_+cj,fl) * v
             Gk(sa,j,nf) = v
             if (cj == 0) Gk(j,sa,nf) = cbk(nd_+a,fl) * Gk(j,sa,nf)
          enddo
       enddo
    enddo
  end subroutine
  Subroutine k_sub_flush_gather(G, ld, nrow, S, k, nf, C, Qt, slot, nfg)
    Integer, intent(in) :: ld, nrow, k, nf, nfg
    Integer, intent(in) :: S(k)
    Complex (Kind=Kind(0.d0)), intent(in) :: G(ld,ld)
    Complex (Kind=Kind(0.d0)), intent(inout) :: C(ld,k), Qt(ld,k)
    Integer, intent(inout) :: slot(ld,nfg)
    Integer :: i, m
    !$omp target teams distribute parallel do collapse(2) map(present,alloc: G, S, C, Qt, slot)
    do m = 1, k
       do i = 1, nrow
          C(i,m)  = G(i,S(m))
          Qt(i,m) = - G(S(m),i)
          if (i == S(m)) then
             Qt(i,m) = Qt(i,m) + 1.d0
             slot(i,nf) = 0
          endif
       enddo
    enddo
  end subroutine
end program
