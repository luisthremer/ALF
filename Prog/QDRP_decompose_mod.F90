!  Copyright (C) 2017, 2018 The ALF project
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.


!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This constructs a decompostion Mat = Q D R P^* using a pivoted QR decomposition
!--------------------------------------------------------------------
Module QDRP_mod
  Use linalg_mod, only: la_zgeqp3, la_zgeqrf, la_zunmqr

  !> Test-only pre-pivoted QR in CGR / CGR2_2 (default off). Compile with ALF_FLAGS_EXT="-DQR_PREPIVOT[=v]"
  !> (v = 1 static norm sort, 2 two-pass, 3 blocked re-ranking with block size -DQR_PPNB=n, default 32).
  !> The call sites use QDRP_prepivot(..., QR_PP_VARIANT, QR_PP_NB) instead of QDRP_decompose.
#if defined(QR_PREPIVOT)
#if QR_PREPIVOT < 1 || QR_PREPIVOT > 3
#error "QR_PREPIVOT must be 1, 2 or 3"
#endif
  Integer, parameter :: QR_PP_VARIANT = QR_PREPIVOT
#else
  Integer, parameter :: QR_PP_VARIANT = 0
#endif
#if defined(QR_PPNB)
  Integer, parameter :: QR_PP_NB = QR_PPNB
#else
  Integer, parameter :: QR_PP_NB = 32
#endif

Contains

!--------------------------------------------------------------------
!> @brief
!> QDRP_decompose with a pre-pivot instead of ZGEQP3's column pivoting (gpu-port tickets 104, 108).
!> @details
!> Same contract as QDRP_decompose (Mat = reflectors and R with D separated, TAU, WORK allocated here,
!> LWORK); IPVT has ZGEQP3 semantics (column j of Mat*P is column IPVT(j) of Mat), entries beyond N_part
!> are the identity. The reflectors form one ZGEQRF-format set, so ZUNMQR and the det(Q) phase work unchanged.
!>   variant 0: plain QDRP_decompose (ZGEQP3).
!>   variant 1: column 2-norms, stable descending sort, unpivoted QR (the CPU form of the device la_zqrpp).
!>   variant 2: two-pass: variant-1 order, unpivoted QR of a copy, re-sort (stable) by |R_jj|, unpivoted QR
!>              of the re-sorted original (2 x geqrf).
!>   variant 3: blocked re-ranking: per panel of nb columns the remaining columns are sorted by their exact
!>              residual 2-norm, the panel is factored (ZGEQRF), Q^H is applied to the trailing columns (ZUNMQR)
!>              and their residual norms are recomputed (DZNRM2, no downdating). nb = 1 gives ZGEQP3's pivots.
!--------------------------------------------------------------------
    SUBROUTINE QDRP_prepivot(Ndim, N_part, Mat, D, IPVT, TAU, WORK, LWORK, variant, nb)
      Implicit None
      Integer, intent(in) :: Ndim, N_part, variant, nb
      Integer, intent(inout) :: LWORK
      Integer, Dimension(:), intent(inout), Allocatable :: IPVT
      COMPLEX(Kind=Kind(0.d0)), Intent(inout) :: Mat(Ndim, N_part)   ! leading Ndim x N_part block of the caller's array
      COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(inout) :: D
      COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(inout), Allocatable :: TAU
      COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(INOUT), Allocatable :: WORK

      Real(Kind=Kind(0.d0)), external :: DZNRM2
      Real(Kind=Kind(0.d0)), allocatable :: nrm(:), RWORK(:)
      Complex(Kind=Kind(0.d0)), allocatable :: tmp(:,:), w1(:)
      Complex(Kind=Kind(0.d0)) :: wq(1)
      Integer, allocatable :: perm(:), ord(:)
      Integer :: i, j, k, b, r, info, lw, nbl
      Real(Kind=Kind(0.d0)) :: X

      if (variant < 1 .or. variant > 3) then
         call QDRP_decompose(Ndim, N_part, Mat, D, IPVT, TAU, WORK, LWORK)
         return
      endif
      allocate(nrm(N_part), tmp(Ndim, N_part), perm(size(IPVT)), ord(N_part))
      do j = 1, N_part
         nrm(j) = DZNRM2(Ndim, Mat(1,j), 1)
      enddo
      do j = 1, size(IPVT)
         perm(j) = j
      enddo

      if (variant /= 3) then
         call sort_desc(nrm, perm(1:N_part))
         if (variant == 2) then
            ! pass 1: unpivoted QR of the sorted copy; only |R_jj| is used
            do j = 1, N_part
               tmp(:,j) = Mat(1:Ndim, perm(j))
            enddo
            call la_zgeqrf(Ndim, N_part, tmp, Ndim, TAU, wq, -1, info)
            lw = int(dble(wq(1)))
            allocate(w1(lw))
            call la_zgeqrf(Ndim, N_part, tmp, Ndim, TAU, w1, lw, info)
            deallocate(w1)
            do j = 1, N_part
               nrm(j) = abs(tmp(j,j))
               ord(j) = j
            enddo
            call sort_desc(nrm, ord)
            ord = perm(ord)
            perm(1:N_part) = ord
         endif
         do j = 1, N_part
            tmp(:,j) = Mat(1:Ndim, perm(j))
         enddo
         Mat(1:Ndim, 1:N_part) = tmp
         IPVT = 1   ! all columns fixed: LAPACK pivoting off
         call QDRP_decompose(Ndim, N_part, Mat, D, IPVT, TAU, WORK, LWORK)
         IPVT = perm
         deallocate(nrm, tmp, perm, ord)
         return
      endif

      ! variant 3
      nbl = max(1, min(nb, N_part))
      call la_zgeqrf(Ndim, nbl, Mat, Ndim, TAU, wq, -1, info)
      lw = int(dble(wq(1)))
      call la_zunmqr('L', 'C', Ndim, N_part, nbl, Mat, Ndim, TAU, Mat, Ndim, wq, -1, info)
      lw = max(lw, int(dble(wq(1))))
      allocate(w1(lw))
      k = 1
      do while (k <= N_part)
         r = N_part - k + 1
         b = min(nbl, r)
         ! re-rank the remaining columns by their current residual norm (stable, descending)
         do i = 1, r
            ord(i) = k + i - 1
         enddo
         call sort_desc(nrm, ord(1:r))
         do i = 1, r
            tmp(:,i) = Mat(1:Ndim, ord(i))
         enddo
         Mat(1:Ndim, k:N_part) = tmp(:, 1:r)
         nrm(k:N_part) = nrm(ord(1:r))
         perm(k:N_part) = perm(ord(1:r))
         ! factor the panel, apply Q^H to the trailing columns, recompute their residual norms exactly
         call la_zgeqrf(Ndim-k+1, b, Mat(k,k), Ndim, TAU(k), w1, lw, info)
         if (k + b <= N_part) then
            call la_zunmqr('L', 'C', Ndim-k+1, N_part-k-b+1, b, Mat(k,k), Ndim, TAU(k), &
                 &         Mat(k,k+b), Ndim, w1, lw, info)
            do j = k + b, N_part
               nrm(j) = DZNRM2(Ndim-k-b+1, Mat(k+b,j), 1)
            enddo
         endif
         k = k + b
      enddo
      deallocate(w1)
      ! same WORK/LWORK contract as QDRP_decompose (WORK is used afterwards by ZUNMQR)
      allocate(RWORK(2*Ndim))
      call la_zgeqp3(Ndim, N_part, Mat, Ndim, IPVT, TAU, wq, -1, RWORK, info)
      LWORK = int(dble(wq(1)))
      allocate(WORK(LWORK))
      deallocate(RWORK)
      ! separate off D (as QDRP_decompose)
      do i = 1, N_part
         X = ABS(Mat(i, i))
         D(i) = X
         do j = i, N_part
            Mat(i, j) = Mat(i, j) / X
         enddo
      enddo
      IPVT = perm
      deallocate(nrm, tmp, perm, ord)

    contains

      !> Stable insertion sort of idx, descending by key(idx(:)).
      subroutine sort_desc(key, idx)
        Real(Kind=Kind(0.d0)), intent(in) :: key(:)
        Integer, intent(inout) :: idx(:)
        Integer :: ii, jj, kk
        do ii = 2, size(idx)
           kk = idx(ii)
           jj = ii - 1
           do while (jj >= 1)
              if (key(idx(jj)) >= key(kk)) exit
              idx(jj+1) = idx(jj)
              jj = jj - 1
           enddo
           idx(jj+1) = kk
        enddo
      end subroutine sort_desc
    END SUBROUTINE QDRP_prepivot

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> ! This constructs a decompostion Mat = Q D R P^* using a pivoted QR decomposition
!
!> @param Ndim[in] The size of the involved matrices
!> @param Mat[inout] The matrix that we want to decompose. The Householder reflectors and the upper
!>            triangular matrix are returned in Mat on exit
!> @param D[inout] The diagonal elements
!> @param IPVT[inout] A Pivoting vector that collects the permutations. Can be used by ?LAPMR and ?LAPMT
!> @param TAU[inout] The scalar factors for the Householder decomposition
!> @param WORK[inout] work memory. We query and allocate it in this routine. Needs to be deallocated outside.
!> @param LWORK[inout] optimal size of the work memory.
!--------------------------------------------------------------------
    SUBROUTINE QDRP_decompose(Ndim, N_part, Mat, D, IPVT, TAU, WORK, LWORK)
      Implicit None
      Integer, intent(in) :: Ndim
      Integer, intent(in) :: N_part
      Integer, intent(inout) :: LWORK
      Integer, Dimension(:), intent(inout), Allocatable :: IPVT
      COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(inout) :: Mat
      COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(inout) :: D
      COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(inout), Allocatable :: TAU
      COMPLEX(Kind=Kind(0.d0)), Dimension(:), Intent(INOUT), Allocatable :: WORK
      
      Real(Kind=Kind(0.d0)), Dimension(:), Allocatable :: RWORK
      COMPLEX(Kind=Kind(0.d0)) :: WORK_QUERY(1)
      Integer :: info, i, j
      Real(Kind=Kind(0.d0)) :: X
      
      ALLOCATE(RWORK(2*Ndim))
      ! Query optimal amount of memory
      call la_zgeqp3(Ndim, N_part, Mat, Ndim, IPVT, TAU, WORK_QUERY, -1, RWORK, INFO)
      LWORK = INT(DBLE(WORK_QUERY(1)))
      ALLOCATE(WORK(LWORK))
      ! QR decomposition of Mat with full column pivoting, Mat * P = Q * R
      call la_zgeqp3(Ndim, N_part, Mat, Ndim, IPVT, TAU, WORK, LWORK, RWORK, INFO)
      DEALLOCATE(RWORK)
      ! separate off D
      do i = 1, N_part
         ! plain diagonal entry
         X = ABS(Mat(i, i))
         !             ! a inf-norm
         !             X = TPUP(i, i+izamax(Ndim+1-i, TPUP(i, i), Ndim)-1)
         !             ! another inf-norm
         !             X = TPUP(i, i-1+izmax1(Ndim+1-i, TPUP(i, i), Ndim))
         !             ! 1-norm
         !            X = DZSUM1(N_size+1-i, TPUP(i, i), N_size)
         ! 2-norm
         !            X = DZNRM2(N_size+1-i, TPUP(i, i), N_size)
         D(i) = X
         do j = i, N_part
            Mat(i, j) = Mat(i, j) / X
         enddo
      enddo
    END SUBROUTINE QDRP_decompose
    
    SUBROUTINE Pivot_phase(Phase, IPVT, N_size)
      Implicit none
      COMPLEX(kind=kind(0.d0)), Intent(INOUT) :: Phase
      Integer, Dimension(:), Intent(IN)       :: IPVT
      Integer,               Intent(IN)       :: N_size
      
      Integer:: i, next, L, VISITED(N_size)
      
      VISITED=0
      do i = 1, N_size
         if (VISITED(i) .eq. 0) then
            next = i
            L = 0
            do while (VISITED(next) .eq. 0)
               L = L + 1
               VISITED(next) = 1
               next = IPVT(next)
            enddo
            if(MOD(L, 2) .eq. 0) then
               PHASE = -PHASE
            endif
         endif
      enddo
    END SUBROUTINE Pivot_phase
    
  End Module QDRP_mod
