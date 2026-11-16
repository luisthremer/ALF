!  Copyright (C) 2026 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
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
!
!     MODIFIED (luis_playground): not part of the original ALF distribution.

module ana_localtau_mod
!--------------------------------------------------------------------
!> @brief
!> Analysis of time-displaced *local* observables "<name>_localtau" stored in
!> data.h5 in real space (Obs_Latt(Ns, Ntau, Norb), Norb = 1).
!>
!> @details
!> Reproduces, bit for bit, the two-step workflow
!>    make_SpinZlocal_tau.py  (numpy Fourier transform -> text file <name>local_tau)
!>    ana.out <name>local_tau (read_latt + ana_tau)
!> directly from data.h5. To be bit-identical with numpy, the transform calls
!> exactly the same routines numpy uses:
!>    qr      = recip_latt @ latt.T   -> cblas_dgemm(RowMajor, N, T)  of numpy's OpenBLAS
!>    phase   = exp(-1j*qr)           -> glibc cexp
!>    spin_q  = spin @ phase.T        -> cblas_zgemm(RowMajor, N, T)  per bin
!> The OpenBLAS library is loaded at run time with dlopen (no symbol clash with MKL).
!> Its path defaults to ALF_LOCALTAU_BLAS_DEFAULT and can be overridden with the
!> environment variable ALF_LOCALTAU_BLAS. Bit-identity holds for the same OpenBLAS
!> build on the same CPU type that ran the Python script.
!>
!> Python conventions copied: sign 1.0 (taken from data.h5; warns if /= 1),
!> background 0, Channel 'PH', only the real part of obser, output prefix
!> "<name>local" (e.g. SpinZ_localtau -> SpinZlocal_*).
!--------------------------------------------------------------------
   use iso_c_binding
   use iso_fortran_env, only: output_unit, error_unit
   use Lattices_v3, only: Unit_cell, Lattice
   use ana_mod,     only: read_local_hdf5, ana_tau
   implicit none
   private
   public :: Cov_localtau

#ifndef ALF_LOCALTAU_BLAS_DEFAULT
#define ALF_LOCALTAU_BLAS_DEFAULT "/home/atuin/b133ae/b133ae29/SOFTWARE/micromamba/envs/general_python/lib/libblas.so.3"
#endif

   integer(c_int), parameter :: CblasRowMajor = 101, CblasNoTrans = 111, CblasTrans = 112
   integer(c_int), parameter :: RTLD_NOW = 2

   interface
      function dlopen(filename, mode) bind(C, name="dlopen")
         import :: c_ptr, c_char, c_int
         type(c_ptr) :: dlopen
         character(kind=c_char), dimension(*), intent(in) :: filename
         integer(c_int), value :: mode
      end function dlopen
      function dlsym(handle, symbol) bind(C, name="dlsym")
         import :: c_ptr, c_funptr, c_char
         type(c_funptr) :: dlsym
         type(c_ptr), value :: handle
         character(kind=c_char), dimension(*), intent(in) :: symbol
      end function dlsym
   end interface

   abstract interface
      subroutine cblas_dgemm_t(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) bind(C)
         import :: c_int, c_double
         integer(c_int), value :: order, transa, transb, m, n, k, lda, ldb, ldc
         real(c_double), value :: alpha, beta
         real(c_double), intent(in)    :: a(*), b(*)
         real(c_double), intent(inout) :: c(*)
      end subroutine cblas_dgemm_t
      subroutine cblas_zgemm_t(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) bind(C)
         import :: c_int, c_double_complex
         integer(c_int), value :: order, transa, transb, m, n, k, lda, ldb, ldc
         complex(c_double_complex), intent(in)    :: alpha, beta     ! const void*
         complex(c_double_complex), intent(in)    :: a(*), b(*)
         complex(c_double_complex), intent(inout) :: c(*)
      end subroutine cblas_zgemm_t
      function cexp_t(z) bind(C) result(r)
         import :: c_double_complex
         complex(c_double_complex), value :: z
         complex(c_double_complex) :: r
      end function cexp_t
   end interface

   procedure(cblas_dgemm_t), pointer, save :: p_dgemm => null()
   procedure(cblas_zgemm_t), pointer, save :: p_zgemm => null()
   procedure(cexp_t),        pointer, save :: p_cexp  => null()

contains

   subroutine load_libs()
      character(len=1024) :: path
      integer :: stat, plen
      type(c_ptr)    :: h_blas, h_m
      type(c_funptr) :: f

      if (associated(p_zgemm)) return

      call get_environment_variable("ALF_LOCALTAU_BLAS", path, plen, stat)
      if (stat /= 0 .or. plen == 0) path = ALF_LOCALTAU_BLAS_DEFAULT

      h_blas = dlopen(trim(path)//c_null_char, RTLD_NOW)
      if (.not. c_associated(h_blas)) then
         write(error_unit,'(3A)') "Cov_localtau: cannot dlopen OpenBLAS '", trim(path), "'."
         write(error_unit,'(A)')  "  Set ALF_LOCALTAU_BLAS to numpy's libopenblas (np.show_config())."
         error stop 1
      endif
      h_m = dlopen("libm.so.6"//c_null_char, RTLD_NOW)
      if (.not. c_associated(h_m)) then
         write(error_unit,'(A)') "Cov_localtau: cannot dlopen libm.so.6"
         error stop 1
      endif

      f = dlsym(h_blas, "cblas_dgemm"//c_null_char)
      if (.not. c_associated(f)) error stop "Cov_localtau: cblas_dgemm not found"
      call c_f_procpointer(f, p_dgemm)
      f = dlsym(h_blas, "cblas_zgemm"//c_null_char)
      if (.not. c_associated(f)) error stop "Cov_localtau: cblas_zgemm not found"
      call c_f_procpointer(f, p_zgemm)
      f = dlsym(h_m, "cexp"//c_null_char)
      if (.not. c_associated(f)) error stop "Cov_localtau: cexp not found"
      call c_f_procpointer(f, p_cexp)

      write(output_unit,'(2A)') " Cov_localtau: using BLAS ", trim(path)
   end subroutine load_libs

!==============================================================================
   subroutine Cov_localtau(name_obs, filename_h5)
      character(len=*), intent(in) :: name_obs
      character(len=*), intent(in) :: filename_h5

      real    (kind=kind(0.d0)), allocatable :: sgn(:)
      complex (kind=kind(0.d0)), pointer     :: bins_loc(:,:,:,:)
      complex (kind=kind(0.d0)), pointer     :: bins(:,:,:,:,:), bins0(:,:)
      type (Lattice)   :: Latt
      type (Unit_cell) :: Latt_unit
      real (kind=kind(0.d0))         :: dtau, x_p(2), t, d
      character (len=:), allocatable :: Channel_h5, Channel
      character (len=64)             :: name_out

      integer :: Ns, Ntau, Nbins, L, off, i, j, ja, jb, nb, nt, n0, nlen
      real    (kind=c_double),         allocatable :: q(:), rl(:,:), lt(:,:), qr(:,:)
      complex (kind=c_double_complex), allocatable :: ph(:,:), A(:,:), Cq(:,:)
      complex (kind=c_double_complex) :: one, zero
      real    (kind=kind(0.d0)), parameter :: pi = 3.141592653589793d0   ! == numpy.pi

      call load_libs()

      call read_local_hdf5(filename_h5, name_obs, sgn, bins_loc, Latt, Latt_unit, dtau, Channel_h5)
      Ns    = size(bins_loc, 1)
      Ntau  = size(bins_loc, 2)
      Nbins = size(bins_loc, 4)
      if (size(bins_loc, 3) /= 1) then
         write(error_unit,*) "Cov_localtau: only Norb = 1 is supported, got ", size(bins_loc, 3)
         error stop 1
      endif
      L = int(sqrt(dble(Ns)))
      if (L*L /= Ns) then
         write(error_unit,*) "Cov_localtau: Ns is not a square L*L: ", Ns
         error stop 1
      endif
      if (any(sgn /= 1.d0)) write(error_unit,*) &
           & "Cov_localtau: WARNING sign /= 1 in some bins; the Python workflow assumed sign = 1."

      ! q = 2*np.pi/L * n,  n = arange(-L//2+1, L//2+1)   (Python floor division)
      allocate(q(L))
      n0 = -((L+1)/2) + 1
      t  = (2.d0*pi) / dble(L)
      do i = 1, L
         q(i) = t * dble(n0 + i - 1)
      enddo

      ! recip_latt (row-major Ns x 2) == rl(2,Ns);  latt (row-major Ns x 2) == lt(2,Ns)
      allocate(rl(2,Ns), lt(2,Ns), qr(Ns,Ns))
      off = L/2 - 1
      do j = 0, Ns-1
         ja = j/L ; jb = mod(j,L)
         rl(1,j+1) = q(ja+1)
         rl(2,j+1) = q(jb+1)
         lt(1,j+1) = dble(j/L - off)
         lt(2,j+1) = dble(mod(j,L) - off)
      enddo

      ! Sanity: momenta must match ALF's Latt%listk ordering (as read_latt checks).
      do j = 1, Ns
         x_p = dble(Latt%listk(j,1))*Latt%b1_p + dble(Latt%listk(j,2))*Latt%b2_p
         d   = (x_p(1)-rl(1,j))**2 + (x_p(2)-rl(2,j))**2
         if (d > 0.00001d0) then
            write(error_unit,*) "Cov_localtau: momenta do not fit", j, x_p, rl(:,j)
            error stop 1
         endif
      enddo

      ! qr (row-major Nq x Nr) == qr(Nr,Nq)  =  recip_latt @ latt.T
      qr = 0.d0
      call p_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, int(Ns,c_int), int(Ns,c_int), 2_c_int, &
           &       1.d0, rl, 2_c_int, lt, 2_c_int, 0.d0, qr, int(Ns,c_int))

      ! phase = exp(-1j*qr)  (row-major Nq x Nr) == ph(Nr,Nq)
      allocate(ph(Ns,Ns))
      do j = 1, Ns
         do i = 1, Ns
            ph(i,j) = p_cexp(cmplx(0.d0, -qr(i,j), kind=c_double_complex))
         enddo
      enddo
      deallocate(qr, lt, q)

      ! spin_q[b] = spin[b] @ phase.T ;  spin[b] row-major (Ntau x Ns) == A(Ns,Ntau)
      allocate(bins(Ns,Ntau,1,1,Nbins), bins0(1,Nbins), A(Ns,Ntau), Cq(Ns,Ntau))
      one  = cmplx(1.d0, 0.d0, kind=c_double_complex)
      zero = cmplx(0.d0, 0.d0, kind=c_double_complex)
      do nb = 1, Nbins
         do nt = 1, Ntau
            do i = 1, Ns
               A(i,nt) = cmplx(dble(bins_loc(i,nt,1,nb)), 0.d0, kind=c_double_complex)
            enddo
         enddo
         call p_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, int(Ntau,c_int), int(Ns,c_int), &
              &       int(Ns,c_int), one, A, int(Ns,c_int), ph, int(Ns,c_int), zero, Cq, int(Ns,c_int))
         bins(:,:,1,1,nb) = Cq
      enddo
      bins0 = cmplx(0.d0, 0.d0, kind(0.d0))
      deallocate(A, Cq, ph, rl)
      deallocate(bins_loc)

      Channel = 'PH'
      nlen = len_trim(name_obs) - len('_localtau')
      write(name_out, '(2A)') name_obs(1:nlen), "local"
      write(output_unit,'(3A)') ' Channel is ', Channel, '   (as in make_SpinZlocal_tau.py)'
      write(output_unit,'(2A)') ' Output prefix: ', trim(name_out)

      call ana_tau(name_out, sgn, bins, bins0, Latt, Latt_unit, dtau, Channel)

      deallocate(bins, bins0, sgn)
   end subroutine Cov_localtau

end module ana_localtau_mod
