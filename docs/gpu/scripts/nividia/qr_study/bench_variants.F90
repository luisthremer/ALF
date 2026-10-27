! Ticket 108: standalone check of the QR_PREPIVOT=2/3 blocks against ZGEQP3 (generated from gen_variants.py text).
module t
contains
subroutine zgeqp3_fixed(M,NC,A,D,IPVT,TAU,WORK,LWORK)
  integer :: M,NC,LWORK; integer, allocatable :: IPVT(:)
  complex(kind(0.d0)) :: A(:,:), D(:); complex(kind(0.d0)), allocatable :: TAU(:), WORK(:)
  complex(kind(0.d0)) :: z; real(kind(0.d0)) :: xx; integer :: ii,jj; real(kind(0.d0)), allocatable :: rw(:); integer :: info
  allocate(rw(2*M)); call ZGEQP3(M,NC,A(1,1),M,IPVT,TAU(1),z,-1,rw,info); LWORK=int(dble(z)); allocate(WORK(LWORK))
  call ZGEQP3(M,NC,A(1,1),M,IPVT,TAU(1),WORK(1),LWORK,rw,info)
  do ii=1,NC; xx=abs(A(ii,ii)); D(ii)=xx; do jj=ii,NC; A(ii,jj)=A(ii,jj)/xx; enddo; enddo
end subroutine
subroutine fac(variant,nb,M,NC,A,D,IPVT,TAU)
  integer :: variant,nb,M,NC; complex(kind(0.d0)) :: A(:,:), D(:); integer, allocatable :: IPVT(:)
  complex(kind(0.d0)), allocatable :: TAU(:), WORK(:); integer :: LWORK
  integer, parameter :: QR_PPNB_DUMMY=0
  if (variant==2) then

        ! Test-only (ticket 108), variant 2 (two-pass): pass 1 = static norm sort + unpivoted QR of a copy;
        ! the columns are then re-sorted by |R_jj| (residual norm in the pass-1 order), and the final unpivoted
        ! QR of the re-sorted original matrix is done by QDRP_decompose with all columns fixed (cost 2 x geqrf).
        ! IPVT has ZGEQP3 semantics (column j of A*P is column IPVT(j) of A).
        block
          Real(Kind=Kind(0.d0)), allocatable :: pp_nrm(:)
          Real(Kind=Kind(0.d0)) :: DZNRM2
          Complex(Kind=Kind(0.d0)), allocatable :: pp_tmp(:,:), pp_work(:), pp_tau(:)
          Complex(Kind=Kind(0.d0)) :: pp_z
          Integer, allocatable :: pp_perm(:), pp_perm2(:), pp_ord(:)
          Integer :: pp_i, pp_j, pp_k, pp_lw, pp_info
          Allocate(pp_nrm(NC), pp_tmp(M,NC), pp_perm(size(IPVT)), pp_perm2(size(IPVT)), pp_ord(NC), pp_tau(NC))
          do pp_j = 1, NC
             pp_nrm(pp_j) = DZNRM2(M, A(1,pp_j), 1)
             pp_perm(pp_j) = pp_j
          enddo
          do pp_j = NC+1, size(IPVT)
             pp_perm(pp_j) = pp_j
          enddo
          do pp_i = 2, NC
             pp_k = pp_perm(pp_i)
             pp_j = pp_i - 1
             do while (pp_j >= 1)
                if (pp_nrm(pp_perm(pp_j)) >= pp_nrm(pp_k)) exit
                pp_perm(pp_j+1) = pp_perm(pp_j)
                pp_j = pp_j - 1
             enddo
             pp_perm(pp_j+1) = pp_k
          enddo
          do pp_j = 1, NC
             pp_tmp(:,pp_j) = A(1:M,pp_perm(pp_j))
          enddo
          ! pass 1: unpivoted QR of the sorted copy (only |R_jj| is used)
          call ZGEQRF(M, NC, pp_tmp(1,1), M, pp_tau(1), pp_z, -1, pp_info)
          pp_lw = INT(DBLE(pp_z))
          Allocate(pp_work(pp_lw))
          call ZGEQRF(M, NC, pp_tmp(1,1), M, pp_tau(1), pp_work(1), pp_lw, pp_info)
          do pp_j = 1, NC
             pp_nrm(pp_j) = ABS(pp_tmp(pp_j,pp_j))
          enddo
          ! pass 2 order: stable descending sort of the pass-1 positions by |R_jj|
          do pp_j = 1, NC
             pp_ord(pp_j) = pp_j
          enddo
          do pp_i = 2, NC
             pp_k = pp_ord(pp_i)
             pp_j = pp_i - 1
             do while (pp_j >= 1)
                if (pp_nrm(pp_ord(pp_j)) >= pp_nrm(pp_k)) exit
                pp_ord(pp_j+1) = pp_ord(pp_j)
                pp_j = pp_j - 1
             enddo
             pp_ord(pp_j+1) = pp_k
          enddo
          pp_perm2 = pp_perm
          do pp_j = 1, NC
             pp_perm2(pp_j) = pp_perm(pp_ord(pp_j))
             pp_tmp(:,pp_j) = A(1:M,pp_perm2(pp_j))
          enddo
          A(1:M,1:NC) = pp_tmp
          IPVT = 1   ! all columns fixed: LAPACK pivoting off
          call zgeqp3_fixed(M,NC,A,D,IPVT,TAU,WORK,LWORK)
          IPVT = pp_perm2
          Deallocate(pp_nrm, pp_tmp, pp_perm, pp_perm2, pp_ord, pp_tau, pp_work)
        end block

  else

        ! Test-only (ticket 108), variant 3 (blocked re-ranking): per block of nb columns the remaining
        ! columns are re-sorted by their current residual 2-norm, the block is factored by ZGEQRF, Q^H is applied
        ! to the trailing columns by ZUNMQR, and the trailing residual norms are recomputed exactly (DZNRM2,
        ! no downdating). Total flops = one geqrf plus O(M*N^2/nb) for the norms. QRCP with block-level
        ! pivoting; nb=1 would be ZGEQP3 with exact norms. D is separated as in QDRP_decompose.
        ! IPVT has ZGEQP3 semantics (column j of A*P is column IPVT(j) of A).
        block
          Real(Kind=Kind(0.d0)), allocatable :: pp_nrm(:), pp_nrm2(:), pp_rwr(:)
          Real(Kind=Kind(0.d0)) :: DZNRM2, pp_x
          Complex(Kind=Kind(0.d0)), allocatable :: pp_tmp(:,:), pp_work(:)
          Complex(Kind=Kind(0.d0)) :: pp_z
          Integer, allocatable :: pp_perm(:), pp_perm2(:), pp_ord(:)
          Integer :: pp_i, pp_j, pp_j0, pp_k, pp_b, pp_r, pp_info, pp_lw, pp_nb
          pp_nb = nb
          Allocate(pp_nrm(NC), pp_nrm2(NC), pp_tmp(M,NC), pp_perm(size(IPVT)), pp_perm2(size(IPVT)), pp_ord(NC))
          do pp_j = 1, NC
             pp_nrm(pp_j) = DZNRM2(M, A(1,pp_j), 1)
          enddo
          do pp_j = 1, size(IPVT)
             pp_perm(pp_j) = pp_j
          enddo
          pp_b = min(pp_nb, NC)
          call ZGEQRF(M, pp_b, A(1,1), M, TAU(1), pp_z, -1, pp_info)
          pp_lw = INT(DBLE(pp_z))
          call ZUNMQR('L', 'C', M, NC, pp_b, A(1,1), M, TAU(1), A(1,1), M, pp_z, -1, pp_info)
          pp_lw = max(pp_lw, INT(DBLE(pp_z)))
          Allocate(pp_work(pp_lw))
          pp_k = 1
          do while (pp_k <= NC)
             pp_r = NC - pp_k + 1
             pp_b = min(pp_nb, pp_r)
             ! re-rank the remaining columns by current residual norm (stable, descending)
             do pp_i = 1, pp_r
                pp_ord(pp_i) = pp_k + pp_i - 1
             enddo
             do pp_i = 2, pp_r
                pp_j0 = pp_ord(pp_i)
                pp_j = pp_i - 1
                do while (pp_j >= 1)
                   if (pp_nrm(pp_ord(pp_j)) >= pp_nrm(pp_j0)) exit
                   pp_ord(pp_j+1) = pp_ord(pp_j)
                   pp_j = pp_j - 1
                enddo
                pp_ord(pp_j+1) = pp_j0
             enddo
             do pp_i = 1, pp_r
                pp_tmp(1:M,pp_i) = A(1:M,pp_ord(pp_i))
                pp_nrm2(pp_i) = pp_nrm(pp_ord(pp_i))
                pp_perm2(pp_i) = pp_perm(pp_ord(pp_i))
             enddo
             A(1:M,pp_k:NC) = pp_tmp(1:M,1:pp_r)
             pp_nrm(pp_k:NC) = pp_nrm2(1:pp_r)
             pp_perm(pp_k:NC) = pp_perm2(1:pp_r)
             ! factor the panel, apply Q^H to the trailing columns, recompute their residual norms exactly
             call ZGEQRF(M-pp_k+1, pp_b, A(pp_k,pp_k), M, TAU(pp_k), pp_work(1), pp_lw, pp_info)
             if (pp_k + pp_b <= NC) then
                call ZUNMQR('L', 'C', M-pp_k+1, NC-pp_k-pp_b+1, pp_b, A(pp_k,pp_k), M, TAU(pp_k), &
                     &      A(pp_k,pp_k+pp_b), M, pp_work(1), pp_lw, pp_info)
                do pp_j = pp_k + pp_b, NC
                   pp_nrm(pp_j) = DZNRM2(M-pp_k-pp_b+1, A(pp_k+pp_b,pp_j), 1)
                enddo
             endif
             pp_k = pp_k + pp_b
          enddo
          ! same LWORK/WORK contract as QDRP_decompose (WORK is used afterwards by ZUNMQR)
          Allocate(pp_rwr(2*M))
          call ZGEQP3(M, NC, A(1,1), M, IPVT, TAU(1), pp_z, -1, pp_rwr(1), pp_info)
          LWORK = INT(DBLE(pp_z))
          Allocate(WORK(LWORK))
          ! separate off D (as QDRP_decompose)
          do pp_i = 1, NC
             pp_x = ABS(A(pp_i, pp_i))
             D(pp_i) = pp_x
             do pp_j = pp_i, NC
                A(pp_i, pp_j) = A(pp_i, pp_j) / pp_x
             enddo
          enddo
          IPVT = pp_perm
          Deallocate(pp_nrm, pp_nrm2, pp_rwr, pp_tmp, pp_work, pp_perm, pp_perm2, pp_ord)
        end block

  endif
end subroutine
end module
program main
  use t
  implicit none
  integer, parameter :: ns(2)=(/400,1024/)
  integer :: in, M, v, nb, rep, c0, c1, cr, info, lw, j, k
  complex(kind(0.d0)), allocatable :: A0(:,:), A(:,:), D(:), TAU(:), WORK(:), TAU2(:)
  integer, allocatable :: IPVT(:)
  real(kind(0.d0)), allocatable :: x(:,:), y(:,:), rw(:)
  complex(kind(0.d0)) :: z
  real(kind(0.d0)) :: tbest, tt, tref, sc
  integer, parameter :: nbl(3)=(/16,32,64/)
  do in=1,2
    M=ns(in)
    allocate(A0(M,M),A(M,M),D(M),IPVT(M),TAU(M),x(M,M),y(M,M),TAU2(M),rw(2*M))
    call random_number(x); call random_number(y); A0=cmplx(x-0.5d0,y-0.5d0,kind(0.d0))
    do j=1,M; sc=10.d0**(-6.d0*dble(mod(j*37,M))/M); A0(:,j)=A0(:,j)*sc; enddo
    do v=0,3
      do k=1,3
        if (v==3) then; nb=nbl(k); elseif (k>1) then; cycle; else; nb=0; endif
        tbest=1d30
        do rep=1,3
          A=A0; IPVT=0
          call system_clock(c0,cr)
          if (v==0) then
            call ZGEQP3(M,M,A(1,1),M,IPVT,TAU(1),z,-1,rw,info); lw=int(dble(z)); allocate(WORK(lw))
            call system_clock(c0)
            call ZGEQP3(M,M,A(1,1),M,IPVT,TAU(1),WORK(1),lw,rw,info); deallocate(WORK)
          elseif (v==1) then
            ! static sort + unpivoted QR (ticket 104 variant)
            call ZGEQRF(M,M,A(1,1),M,TAU(1),z,-1,info); lw=int(dble(z)); allocate(WORK(lw))
            call system_clock(c0)
            call static_qr(M,A,IPVT,TAU,WORK,lw); deallocate(WORK)
          else
            call system_clock(c0)
            call fac(v,nb,M,M,A,D,IPVT,TAU)
          endif
          call system_clock(c1); tt=dble(c1-c0)/cr
          tbest=min(tbest,tt)
        enddo
        if (v==0) tref=tbest
        write(*,'(a,i5,a,i2,a,i3,a,f9.4,a,f6.2)') 'N',M,' variant',v,' nb',nb,'  best time s',tbest,'  ratio to ZGEQP3',tbest/tref
      enddo
    enddo
    deallocate(A0,A,D,IPVT,TAU,x,y,TAU2,rw)
  enddo
contains
  subroutine static_qr(M,A,IPVT,TAU,WORK,lw)
    integer :: M, lw; complex(kind(0.d0)) :: A(:,:), TAU(:), WORK(:); integer :: IPVT(:)
    real(kind(0.d0)) :: nrm(M), DZNRM2; integer :: i,j,k,info; complex(kind(0.d0)), allocatable :: tmp(:,:)
    allocate(tmp(M,M))
    do j=1,M; nrm(j)=DZNRM2(M,A(1,j),1); IPVT(j)=j; enddo
    do i=2,M
      k=IPVT(i); j=i-1
      do while (j>=1); if (nrm(IPVT(j))>=nrm(k)) exit; IPVT(j+1)=IPVT(j); j=j-1; enddo
      IPVT(j+1)=k
    enddo
    do j=1,M; tmp(:,j)=A(:,IPVT(j)); enddo
    A=tmp
    call ZGEQRF(M,M,A(1,1),M,TAU(1),WORK(1),lw,info)
  end subroutine
end program
