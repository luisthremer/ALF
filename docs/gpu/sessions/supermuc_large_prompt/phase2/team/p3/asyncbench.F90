include "mkl_blas_omp_offload.f90"
! Gate A microbenchmark: dependent chain of small device kernels + small dispatched zgemm (k=32), as in the update phase.
! Build with -DASYNC for nowait + DEPCL on every target region and dispatch, one taskwait per chain.
! Optional host work between launches: HOSTWORK=<n> iterations of a dependent flop loop (models the host block update).
#ifdef NODEP
#define DEPCL
#else
#define DEPCL depend(inout: tok)
#endif
program asyncbench
  use omp_lib
  use iso_c_binding
  use onemkl_blas_omp_offload_lp64, only: zgemm
  implicit none
  integer, parameter :: n = 1024, k = 32, NCH = 32, NGEMM_EVERY = 8
  integer :: nrep, i, ic, j, nhw, ih
  real(8) :: t0, t1, t2, hw
  complex(8), allocatable :: Am(:,:), Bm(:,:), Cm(:,:)
  real(8), allocatable :: v(:)
  complex(8), parameter :: one = (1.d0,0.d0), zero = (0.d0,0.d0)
  integer :: tok, nkern
  real(8) :: cc, dd, chk
  character(len=32) :: arg
  nrep = 300; nhw = 0
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read(arg,*) nrep; end if
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read(arg,*) nhw; end if
  call omp_set_default_device(0)
  allocate(Am(n,k), Bm(k,k), Cm(n,k), v(n))
  tok = 0
  call reset()
  !$omp target enter data map(to: Am, Bm, Cm, v)
  ! warm-up
  do ic = 1, 3
    call chain(hw)
  end do
  !$omp target update from(v)
  call reset()
  !$omp target update to(v, Am, Bm, Cm)
  t0 = omp_get_wtime()
  hw = 0.d0
  t1 = 0.d0
  do i = 1, nrep
    call chain(hw)
  end do
  t2 = omp_get_wtime()
#ifdef ASYNC
  !$omp taskwait
#endif
  t1 = omp_get_wtime()
  !$omp target update from(v, Am)
  chk = sum(v) + sum(real(Am)) + sum(aimag(Am))
  nkern = nrep*(NCH + 2*(NCH/NGEMM_EVERY))      ! kernels + (gemm + copy kernel) per group
  write(*,'(A,I0,A,I0,A,ES24.16,A,ES9.2)') 'nrep=', nrep, ' hostwork=', nhw, ' checksum=', chk, ' hw=', hw
  write(*,'(A,F8.2,A,F8.2,A)') 'RESULT total ', 1.d6*(t1-t0)/nkern, ' us/op (kernel or gemm)   issue-only ', &
       1.d6*(t2-t0)/nkern, ' us/op'
contains
  subroutine reset()
    do j = 1, n; v(j) = 1.d0 + 1.d-3*j; end do
    do j = 1, k
      do i = 1, n; Am(i,j) = cmplx(1.d0/(1+mod(i+j,7)), 1.d-2*j, 8); end do
    end do
    Bm = zero
    do j = 1, k; Bm(j,j) = cmplx(0.5d0, 0.5d0, 8); if (j < k) Bm(j+1,j) = cmplx(0.1d0, -0.05d0, 8); end do
    Cm = zero
  end subroutine
  subroutine chain(h)
    real(8), intent(inout) :: h
    integer :: kk, jj, ii
    do kk = 1, NCH
      cc = 0.999d0 - 1.d-4*kk; dd = 1.d-3*kk
      call kern(cc, dd)
      do ih = 1, nhw
        h = h*0.9999999d0 + 1.d-9*sin(h+ih)
      end do
      if (mod(kk, NGEMM_EVERY) == 0) then
#ifdef ASYNC
        !$omp dispatch nowait DEPCL
#else
        !$omp dispatch
#endif
        call zgemm('N','N',n,k,k,one,Am,n,Bm,k,zero,Cm,n)
        call kcopy()
      end if
    end do
  end subroutine
  subroutine kern(c, d)
    real(8), intent(in) :: c, d
    integer :: jj
#ifdef ASYNC
    !$omp target teams distribute parallel do nowait DEPCL map(present, alloc: v)
#else
    !$omp target teams distribute parallel do map(present, alloc: v)
#endif
    do jj = 1, n
      v(jj) = v(jj)*c + d*(1.d0 + 1.d-4*jj)
    end do
  end subroutine
  subroutine kcopy()
    integer :: jj, ii
#ifdef ASYNC
    !$omp target teams distribute parallel do collapse(2) nowait DEPCL map(present, alloc: Am, Cm, v)
#else
    !$omp target teams distribute parallel do collapse(2) map(present, alloc: Am, Cm, v)
#endif
    do jj = 1, k
      do ii = 1, n
        Am(ii,jj) = Cm(ii,jj)*cmplx(0.9d0, 0.d0, 8) + cmplx(1.d-3*v(ii), 0.d0, 8)
      end do
    end do
  end subroutine
end program
