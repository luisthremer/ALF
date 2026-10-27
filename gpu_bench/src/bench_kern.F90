! Device kernels: OpenMP target teams distribute parallel do on plain contiguous arrays.
! All arrays must already be present in the device data environment (mapped by the caller).
! In the HOST build the same directives run on the host (nvfortran -mp).
module bench_kern
  implicit none
contains

  ! G <- M G on rows (i,j) of every bond: M = [[c,s],[s,c]].  Bonds of one family are disjoint.
  ! Loop order: bond index fastest so neighbouring threads touch neighbouring rows i (stride-2 bonds).
  subroutine k_cb_left(G, N, nb, bi, bj, c, s)
    integer, intent(in) :: N, nb
    complex(8) :: G(N,N)
    integer :: bi(nb), bj(nb)
    real(8), intent(in) :: c, s
    integer :: l, b, i, j
    complex(8) :: gi, gj
    !$omp target teams distribute parallel do collapse(2) private(i,j,gi,gj)
    do l = 1, N
      do b = 1, nb
        i = bi(b); j = bj(b)
        gi = G(i,l); gj = G(j,l)
        G(i,l) = c*gi + s*gj
        G(j,l) = s*gi + c*gj
      end do
    end do
  end subroutine

  ! G <- G M^-1 on columns (i,j): M^-1 = [[c,-s],[-s,c]].  Row index fastest (coalesced).
  subroutine k_cb_right(G, N, nb, bi, bj, c, s)
    integer, intent(in) :: N, nb
    complex(8) :: G(N,N)
    integer :: bi(nb), bj(nb)
    real(8), intent(in) :: c, s
    integer :: k, b, i, j
    complex(8) :: gi, gj
    !$omp target teams distribute parallel do collapse(2) private(i,j,gi,gj)
    do b = 1, nb
      do k = 1, N
        i = bi(b); j = bj(b)
        gi = G(k,i); gj = G(k,j)
        G(k,i) = c*gi - s*gj
        G(k,j) = -s*gi + c*gj
      end do
    end do
  end subroutine

  ! G <- B G B^-1 with B = F1 F2 F3 F4 (family f has bonds bi(:,f), bj(:,f)); 8 launches.
  subroutine cb_apply(G, N, nb, bi, bj, c, s)
    integer, intent(in) :: N, nb
    complex(8) :: G(N,N)
    integer :: bi(nb,4), bj(nb,4)
    real(8), intent(in) :: c, s
    integer :: f
    do f = 4, 1, -1
      call k_cb_left(G, N, nb, bi(1,f), bj(1,f), c, s)
    end do
    do f = 4, 1, -1
      call k_cb_right(G, N, nb, bi(1,f), bj(1,f), c, s)
    end do
  end subroutine

  ! Op_V diagonal wrap for one site i: row i * ea, column i / ea (G(i,i) unchanged).
  subroutine k_opv_site(G, N, i, ea)
    integer, intent(in) :: N, i
    complex(8) :: G(N,N)
    real(8), intent(in) :: ea
    integer :: m
    !$omp target teams distribute parallel do
    do m = 1, N
      if (m /= i) then
        G(i,m) = G(i,m) * ea
        G(m,i) = G(m,i) / ea
      end if
    end do
  end subroutine

  ! All sites at once: G_kl <- e^{a_k} G_kl e^{-a_l}, ea(k) = e^{a_k}.
  subroutine k_opv_all(G, N, ea)
    integer, intent(in) :: N
    complex(8) :: G(N,N)
    real(8) :: ea(N)
    integer :: k, l
    !$omp target teams distribute parallel do collapse(2)
    do l = 1, N
      do k = 1, N
        G(k,l) = G(k,l) * (ea(k) / ea(l))
      end do
    end do
  end subroutine

  ! G <- G + x y^T, y has stride incy
  subroutine k_ger(G, N, x, y, incy)
    integer, intent(in) :: N, incy
    complex(8) :: G(N,N), x(N), y(incy*(N-1)+1)
    integer :: k, l
    !$omp target teams distribute parallel do collapse(2)
    do l = 1, N
      do k = 1, N
        G(k,l) = G(k,l) + x(k) * y(1+(l-1)*incy)
      end do
    end do
  end subroutine

  ! same as k_ger but asynchronous (target task, ordered by depend on G); caller does !$omp taskwait
  subroutine k_ger_nw(G, N, x, y, incy)
    integer, intent(in) :: N, incy
    complex(8) :: G(N,N), x(N), y(incy*(N-1)+1)
    integer :: k, l
    !$omp target teams distribute parallel do collapse(2) nowait depend(inout: G)
    do l = 1, N
      do k = 1, N
        G(k,l) = G(k,l) + x(k) * y(1+(l-1)*incy)
      end do
    end do
  end subroutine

  subroutine k_copy(src, dst, n)
    integer, intent(in) :: n
    complex(8) :: src(n), dst(n)
    integer :: i
    !$omp target teams distribute parallel do
    do i = 1, n
      dst(i) = src(i)
    end do
  end subroutine

  ! B = A^T (plain transpose, naive: no shared-memory tiling in OpenMP)
  subroutine k_transpose(A, B, N)
    integer, intent(in) :: N
    complex(8) :: A(N,N), B(N,N)
    integer :: i, j
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, N
      do i = 1, N
        B(j,i) = A(i,j)
      end do
    end do
  end subroutine

  ! Column 2-norms squared (sum of |a|^2); one team per column, reduction inside the team.
  subroutine k_colnorm(A, N, nrm)
    integer, intent(in) :: N
    complex(8) :: A(N,N)
    real(8) :: nrm(N)
    integer :: i, j
    real(8) :: s
    !$omp target teams distribute private(s)
    do j = 1, N
      s = 0d0
      !$omp parallel do reduction(+:s)
      do i = 1, N
        s = s + real(A(i,j))**2 + aimag(A(i,j))**2
      end do
      nrm(j) = s
    end do
  end subroutine

  ! Stable descending sort by rank counting: perm(rank(j)) = j, rank(j) = 1 + #{k : n_k > n_j or
  ! (n_k == n_j and k < j)}.  O(N^2) work, perfectly parallel, no device sort library needed.
  subroutine k_sortperm(nrm, N, perm)
    integer, intent(in) :: N
    real(8) :: nrm(N)
    integer :: perm(N)
    integer :: j, k, r
    !$omp target teams distribute parallel do private(r)
    do j = 1, N
      r = 1
      do k = 1, N
        if (nrm(k) > nrm(j) .or. (nrm(k) == nrm(j) .and. k < j)) r = r + 1
      end do
      perm(r) = j
    end do
  end subroutine

  ! B(:,jj) = A(:,perm(jj))
  subroutine k_permcols(A, N, perm, B)
    integer, intent(in) :: N
    complex(8) :: A(N,N), B(N,N)
    integer :: perm(N)
    integer :: i, jj
    !$omp target teams distribute parallel do collapse(2)
    do jj = 1, N
      do i = 1, N
        B(i,jj) = A(i,perm(jj))
      end do
    end do
  end subroutine

end module bench_kern
