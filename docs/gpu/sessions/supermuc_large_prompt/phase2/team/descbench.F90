! Team session: per-launch cost of an OpenMP target kernel on PVC depending on how its arrays are passed.
!  M: module allocatable arrays referenced directly in the target region (as gpu_lowrank_mod's k_blk_* kernels do)
!  D: the same kernel with the arrays passed as explicit-shape dummy arguments (no descriptors)
!  R: M, but the module arrays are 'declare target' (descriptor resident on the device)
! Same arithmetic in all three (k_blk_scale_rows shape: b = 32 rows of a 576 x 576 x 2 complex array).
module dm
  implicit none
  complex(8), allocatable :: G1(:,:,:), cb1(:,:)
  integer, allocatable :: bs1(:,:)
  complex(8), allocatable :: G2(:,:,:), cb2(:,:)
  integer, allocatable :: bs2(:,:)
  !$omp declare target(G2, cb2, bs2)
contains
  subroutine k_mod(nf, fl, ld, b)
    integer, intent(in) :: nf, fl, ld, b
    integer :: a, j
    !$omp target teams distribute parallel do collapse(2)
    do a = 1, b
      do j = 1, ld
        G1(bs1(a,fl),j,nf) = cb1(a,fl) * G1(bs1(a,fl),j,nf)
      enddo
    enddo
  end subroutine
  subroutine k_res(nf, fl, ld, b)
    integer, intent(in) :: nf, fl, ld, b
    integer :: a, j
    !$omp target teams distribute parallel do collapse(2)
    do a = 1, b
      do j = 1, ld
        G2(bs2(a,fl),j,nf) = cb2(a,fl) * G2(bs2(a,fl),j,nf)
      enddo
    enddo
  end subroutine
  subroutine k_dum(nf, fl, ld, b, nd, G, cb, bs)
    integer, intent(in) :: nf, fl, ld, b, nd
    complex(8) :: G(ld,ld,2), cb(nd,2)
    integer :: bs(nd,2)
    integer :: a, j
    !$omp target teams distribute parallel do collapse(2) map(present,alloc: G, cb, bs)
    do a = 1, b
      do j = 1, ld
        G(bs(a,fl),j,nf) = cb(a,fl) * G(bs(a,fl),j,nf)
      enddo
    enddo
  end subroutine
end module dm

program descbench
  use dm
  use omp_lib
  implicit none
  integer, parameter :: ld = 576, b = 32, nd = 64, NREP = 4000
  integer :: i, a
  real(8) :: t0, t1
  allocate(G1(ld,ld,2), cb1(nd,2), bs1(nd,2), G2(ld,ld,2), cb2(nd,2), bs2(nd,2))
  G1 = (1.d0,0.d0); cb1 = (1.d0,0.d0); G2 = G1; cb2 = cb1
  do a = 1, nd
    bs1(a,:) = a; bs2(a,:) = a
  enddo
  !$omp target enter data map(to: G1, cb1, bs1)
  !$omp target update to(G2, cb2, bs2)
  do i = 1, 50
    call k_mod(1, 1, ld, b); call k_res(1, 1, ld, b); call k_dum(1, 1, ld, b, nd, G1, cb1, bs1)
  enddo
  t0 = omp_get_wtime(); do i = 1, NREP; call k_mod(1 + mod(i,2), 1, ld, b); enddo; t1 = omp_get_wtime()
  write(*,'(A,F9.2,A)') 'RESULT M module allocatables      ', 1d6*(t1-t0)/NREP, ' us/launch'
  t0 = omp_get_wtime(); do i = 1, NREP; call k_res(1 + mod(i,2), 1, ld, b); enddo; t1 = omp_get_wtime()
  write(*,'(A,F9.2,A)') 'RESULT R declare target module    ', 1d6*(t1-t0)/NREP, ' us/launch'
  t0 = omp_get_wtime(); do i = 1, NREP; call k_dum(1 + mod(i,2), 1, ld, b, nd, G1, cb1, bs1); enddo; t1 = omp_get_wtime()
  write(*,'(A,F9.2,A)') 'RESULT D explicit-shape dummies   ', 1d6*(t1-t0)/NREP, ' us/launch'
  !$omp target update from(G1, G2)
  write(*,'(A,2ES12.4)') 'check (must be equal, 1): ', abs(G1(5,7,1)), abs(G2(5,7,1))
end program descbench
