!  Copyright (C) 2026 The ALF project
!
!  This file is part of the ALF project.
!
!  The ALF project is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  The ALF project is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with ALF.  If not, see http://www.gnu.org/licenses/.
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
!  New file of the delayed/submatrix research branch (feature/delayed-submatrix, 2026): BLAS ceiling microbenchmark.

!--------------------------------------------------------------------
!> @brief
!> BLAS ceiling microbenchmark of the low-rank update study: the kernels that bound the delayed/submatrix schemes,
!> measured with the BLAS ALF is linked against. One process = one rank; testsuite/lowrank/blas_bench.sh runs one per
!> rank of a layout and writes the provenance.
!
!> Two builds: noMPI (the local launcher: ranks synchronise through ready files) and -DBLAS_BENCH_MPI (MPI_Init, so Intel MPI
!> pins the ranks as for ALF; rank from MPI_Comm_rank; MPI_Barrier). The ranks are synchronised before EVERY measurement point
!> (timer cost, each kernel x N x k), so the ranks measure the same kernel at the same time.
!
!> Usage:  blas_bench [KEY=VALUE ...] | --probe | --help
!>   N=400,1600,3600   matrix sizes (default 400,1600,3600; up to 8000 and beyond: memory 16*N^2 bytes per matrix, the shell
!>                     script drops the N that do not fit R ranks into the node memory)
!>   KS=16,64,256      inner dimensions of the ZGEMM N x N x k kernel (default 16,64,256)
!>   KERNELS=a,b,..    subset of zgemm_k, zgemm_n, zgeru, zgemv, zgeqp3 (default all five)
!>   REPS=7            timed repetitions per point (>= 5 recommended; the median is reported)
!>   WARMUP=2          untimed warm-up calls per point (the last one calibrates the batch size); ZGEQP3: always 1
!>   MINTIME=0.02      a timed repetition is a batch of calls lasting >= MINTIME seconds (not for zgeqp3: 1 call)
!>   CACHE_MB=256      the BLAS-2 kernels (zgeru, zgemv) rotate over copies of the matrix covering >= this many MB
!>                     (>= 1 copy), so that the matrix does not come from cache (blas_bench.sh sets it from the L3 size)
!>   TIMER_REPS=2000000  push/pop pairs per timer-cost sample (0 = skip the timer-cost measurement)
!>   RANK=0 THREADS=1 NRANKS=1   recorded in every line (the shell script passes them). With MPI, RANK/NRANKS must agree
!>                     with MPI_COMM_WORLD or the run stops
!>   READYDIR=DIR READY_TIMEOUT=900   noMPI, NRANKS > 1: the point barrier through files DIR/ready_<point>_<rank>; a rank
!>                     that waits longer than READY_TIMEOUT seconds stops the run with an error (never silently continues)
!>   START=FILE        wait until the wall-clock epoch (seconds, in FILE) before the first measurement
!>   --probe           one 8x8x8 ZGEMM and one 8x8 ZGEQP3, prints "probe ok", exit 0 (MKL_VERBOSE probe of the campaign)
!> Output (stdout, one line per record, whitespace separated; the header lines start with '#'):
!>   META key=value ...                     build/run facts (compiler, options, clock, MPI, timers compiled in)
!>   SYNC rank start_epoch release_epoch lag_s      with START: how late the process was released
!>   TIMERCOST rank threads what calls_per_sample samples t_median_ns t_min_ns t_max_ns
!>        what = push_pop (timer_push+timer_pop of Prog/Timers_mod.F90 compiled with -DALF_TIMERS; the cost of one
!>        start/stop pair), switch (one timer_switch), system_clock (one call, the clock of tnow) and wtime
!>        (omp_get_wtime alone; = the system_clock cost in a build without OpenMP)
!>   RESULT rank threads kernel N k reps inner t_median_s t_min_s t_max_s gflops gflops_best gbs t_release nbuf
!>        times are per call; gflops uses the median, gflops_best the minimum time; gbs = nominal minimum memory
!>        traffic / median time (not a measured bandwidth); t_release = wall-clock epoch at which this rank left the point
!>        barrier (the spread over the ranks is the release skew of the point); nbuf = rotating copies (BLAS-2 kernels)
!> Flop and byte models (complex double, 16 bytes per element):
!>   zgemm N x N x k  : 8 N^2 k flops, C := C - A B (alpha = -1, beta = 1, as the flush of Lowrank_update_mod.F90),
!>                      bytes 16 (2 N k + 2 N^2) (A, B read, C read and written)
!>   zgeru N x N      : 8 N^2 flops, A := A + alpha x y^T, bytes 16 (2 N^2 + 2 N)
!>   zgemv N x N      : 8 N^2 flops, y := A x, bytes 16 (N^2 + 2 N)
!>   zgeqp3 N x N     : (16/3) N^3 flops (4 x the real Householder count 4/3 N^3), bytes 32 N^2 (nominal); the input
!>                      is restored from a master copy before every call (warm-up included), outside the timed region
!> The buffers are filled in parallel (OpenMP, static schedule) so that the first touch follows the threads of the kernel.
!> The BLAS-2 numbers come from memory (rotating buffers); ALF's G may stay in cache at small N (a lower bound there).
!> Never uses ZDOTU/ZDOTC (Accelerate's f2c return convention).
program blas_bench
   use iso_fortran_env, only: int64, error_unit, output_unit, compiler_version, compiler_options
   use Timers_mod
#ifdef _OPENMP
   use omp_lib
#endif
#ifdef BLAS_BENCH_MPI
   use mpi
#endif
   implicit none

   Integer, parameter :: dp = kind(0.d0)
   Integer, parameter :: MAXLIST = 64
   Integer, parameter :: NKERN = 5
   Character (len=8), parameter :: kname(NKERN) = [Character(len=8) :: "zgemm_k", "zgemm_n", "zgeru", "zgemv", "zgeqp3"]

   Integer  :: ns(MAXLIST) = 0, ks(MAXLIST) = 0, nn = 0, nk = 0
   Logical  :: run_kernel(NKERN) = .true.
   Integer  :: reps = 7, warmup = 2, rank = 0, threads = 1, nranks = 1
   Integer  :: timer_reps = 2000000
   Real (Kind=dp) :: mintime = 0.02_dp, cache_mb = 256.0_dp, ready_timeout = 900.0_dp
   Character (len=1024) :: start_file = "", ready_dir = ""
   Logical  :: probe = .false.
   Integer  :: barrier_seq = 0
   Real (Kind=dp) :: t_release = 0
   Integer  :: ierr_mpi = 0

   Integer  :: iarg, ik, in, ikern

   ! state of the kernel under test (shared by bench_point and one_call)
   Complex (Kind=dp), allocatable :: a(:,:,:), b(:,:), c(:,:), x(:), y(:), work(:), tau(:)
   Real (Kind=dp), allocatable :: rwork(:)
   Integer, allocatable :: jpvt(:)
   Complex (Kind=dp), parameter :: one = (1.0_dp, 0.0_dp), zero = (0.0_dp, 0.0_dp), mone = (-1.0_dp, 0.0_dp)
   Complex (Kind=dp) :: alpha
   Integer :: cur_kern = 0, cur_n = 0, cur_k = 0, nbuf = 1, ic = 0, lwork = 1
   Character (len=1024) :: arg

   ns(1:3) = [400, 1600, 3600]; nn = 3
   ks(1:3) = [16, 64, 256];     nk = 3

   do iarg = 1, command_argument_count()
      call get_command_argument(iarg, arg)
      call parse_arg(trim(arg))
   enddo

   if (probe) then
      call run_probe()
      stop
   endif
   if (nranks < 1 .or. rank < 0 .or. rank >= nranks) call die("RANK/NRANKS out of range")

#ifdef BLAS_BENCH_MPI
   call MPI_Init(ierr_mpi)
   block
      Integer :: r_mpi, n_mpi
      call MPI_Comm_rank(MPI_COMM_WORLD, r_mpi, ierr_mpi)
      call MPI_Comm_size(MPI_COMM_WORLD, n_mpi, ierr_mpi)
      if (r_mpi /= rank .or. n_mpi /= nranks) then
         write(error_unit,'(A,4(1X,I0))') 'blas_bench: MPI rank/size /= RANK/NRANKS arguments:', r_mpi, n_mpi, rank, nranks
         call MPI_Abort(MPI_COMM_WORLD, 2, ierr_mpi)
         error stop 2
      endif
   end block
#endif

   call print_header()
   if (len_trim(start_file) > 0) call wait_for_start(trim(start_file))
   if (timer_reps > 0) then
      call point_barrier()
      call timer_cost()
   endif

   do in = 1, nn
      do ikern = 1, NKERN
         if (.not. run_kernel(ikern)) cycle
         if (ikern == 1) then
            do ik = 1, nk
               call bench_point(ikern, ns(in), ks(ik))
            enddo
         elseif (ikern == 2) then
            call bench_point(ikern, ns(in), ns(in))
         else
            call bench_point(ikern, ns(in), 0)
         endif
      enddo
   enddo
   write(output_unit,'(A)') "# done"
   flush(output_unit)
#ifdef BLAS_BENCH_MPI
   call MPI_Finalize(ierr_mpi)
#endif

contains

   !-----------------------------------------------------------------------------------------------------------
   subroutine die(msg)
      Character (len=*), intent(in) :: msg
      write(error_unit,'(A)') 'blas_bench: ' // msg
#ifdef BLAS_BENCH_MPI
      if (.not. probe) call MPI_Abort(MPI_COMM_WORLD, 2, ierr_mpi)
#endif
      error stop 2
   end subroutine die

   subroutine parse_list(s, list, cnt)
      Character (len=*), intent(in) :: s
      Integer, intent(out) :: list(MAXLIST), cnt
      Integer :: p, q, ios
      cnt = 0
      p = 1
      do while (p <= len_trim(s))
         q = index(s(p:), ',')
         if (q == 0) then
            q = len_trim(s) + 1
         else
            q = p + q - 1
         endif
         cnt = cnt + 1
         if (cnt > MAXLIST) call die("list too long: " // s)
         read(s(p:q-1), *, iostat=ios) list(cnt)
         if (ios /= 0 .or. list(cnt) < 1) call die("bad list entry in: " // s)
         p = q + 1
      enddo
      if (cnt == 0) call die("empty list")
   end subroutine parse_list

   subroutine parse_arg(a)
      Character (len=*), intent(in) :: a
      Integer :: eq, ios, p, q, i
      Character (len=len(a)) :: key, val
      Character (len=16) :: name
      Logical :: found
      if (a == "--probe") then
         probe = .true.; return
      endif
      if (a == "--help" .or. a == "-h") then
         write(output_unit,'(A)') "usage: blas_bench [N=400,1600] [KS=16,64,256] [KERNELS=zgemm_k,zgemm_n,zgeru,zgemv,zgeqp3] " // &
              & "[REPS=7] [WARMUP=2] [MINTIME=0.02] [CACHE_MB=256] [TIMER_REPS=2000000] [RANK=0] [THREADS=1] [NRANKS=1] " // &
              & "[READYDIR=DIR] [READY_TIMEOUT=900] [START=FILE] | --probe"
         stop
      endif
      eq = index(a, '=')
      if (eq < 2) call die("argument not KEY=VALUE: " // a)
      key = a(1:eq-1); val = a(eq+1:)
      do i = 1, len_trim(key)
         if (key(i:i) >= 'a' .and. key(i:i) <= 'z') key(i:i) = achar(iachar(key(i:i)) - 32)
      enddo
      ios = 0
      select case (trim(key))
      case ("N");          call parse_list(trim(val), ns, nn)
      case ("KS");         call parse_list(trim(val), ks, nk)
      case ("REPS");       read(val, *, iostat=ios) reps
      case ("WARMUP");     read(val, *, iostat=ios) warmup
      case ("MINTIME");    read(val, *, iostat=ios) mintime
      case ("CACHE_MB");   read(val, *, iostat=ios) cache_mb
      case ("TIMER_REPS"); read(val, *, iostat=ios) timer_reps
      case ("RANK");       read(val, *, iostat=ios) rank
      case ("THREADS");    read(val, *, iostat=ios) threads
      case ("NRANKS");     read(val, *, iostat=ios) nranks
      case ("READY_TIMEOUT"); read(val, *, iostat=ios) ready_timeout
      case ("READYDIR");   ready_dir = val
      case ("START");      start_file = val
      case ("KERNELS")
         run_kernel = .false.
         p = 1
         do while (p <= len_trim(val))
            q = index(val(p:), ',')
            if (q == 0) then
               q = len_trim(val) + 1
            else
               q = p + q - 1
            endif
            name = val(p:q-1)
            found = .false.
            do i = 1, NKERN
               if (trim(name) == trim(kname(i))) then
                  run_kernel(i) = .true.
                  found = .true.
               endif
            enddo
            if (.not. found) call die("unknown kernel: " // trim(name))
            p = q + 1
         enddo
      case default
         call die("unknown key: " // trim(key))
      end select
      if (ios /= 0) call die("bad value in: " // a)
      if (reps < 1 .or. warmup < 0) call die("REPS must be >= 1, WARMUP >= 0")
   end subroutine parse_arg

   !-----------------------------------------------------------------------------------------------------------
   Real (Kind=dp) function wall()
      Integer (Kind=int64) :: c, r
#ifdef _OPENMP
      wall = omp_get_wtime()
#else
      call system_clock(c, r)
      wall = dble(c)/dble(r)
#endif
   end function wall

   !> seconds since 1970-01-01 00:00:00 UTC from date_and_time (ms resolution), standard Fortran only
   Real (Kind=dp) function epoch_now()
      Integer :: v(8), y, m, era, yoe, doy, doe
      Integer (Kind=int64) :: days
      call date_and_time(values=v)
      y = v(1); m = v(2)
      if (m <= 2) y = y - 1
      if (y >= 0) then
         era = y/400
      else
         era = (y - 399)/400
      endif
      yoe = y - era*400
      if (m > 2) then
         doy = (153*(m - 3) + 2)/5 + v(3) - 1
      else
         doy = (153*(m + 9) + 2)/5 + v(3) - 1
      endif
      doe = yoe*365 + yoe/4 - yoe/100 + doy
      days = int(era, int64)*146097_int64 + doe - 719468_int64
      epoch_now = dble(days)*86400.d0 + dble(v(5))*3600.d0 + dble(v(6))*60.d0 + dble(v(7)) + dble(v(8))/1000.d0 &
           &      - dble(v(4))*60.d0
   end function epoch_now

   subroutine wait_for_start(fname)
      Character (len=*), intent(in) :: fname
      Integer :: u, ios
      Real (Kind=dp) :: t_start, t_now
      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) call die("cannot open START file " // fname)
      read(u, *, iostat=ios) t_start
      close(u)
      if (ios /= 0) call die("START file has no epoch time: " // fname)
      t_now = epoch_now()
      if (t_start - t_now > 300.d0) call die("START time is more than 300 s in the future")
      do
         t_now = epoch_now()
         if (t_now >= t_start) exit
      enddo
      write(output_unit,'(A,1X,I0,2(1X,F16.3),1X,F9.3)') "SYNC", rank, t_start, t_now, t_now - t_start
      flush(output_unit)
   end subroutine wait_for_start

   !> All ranks leave together: MPI_Barrier (MPI build) or the ready files (noMPI, NRANKS > 1; timeout = error).
   !> Called before the timer cost and before every kernel x N x k point, so that all ranks run the same point.
   subroutine point_barrier()
      Integer :: r, u, ios
      Logical :: ex, all_ready
      Real (Kind=dp) :: tw
      Character (len=1200) :: f
      barrier_seq = barrier_seq + 1
#ifdef BLAS_BENCH_MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)
#else
      if (nranks > 1) then
         if (len_trim(ready_dir) == 0) call die("NRANKS > 1 without MPI needs READYDIR (the point barrier)")
         write(f,'(A,A,I0,A,I0)') trim(ready_dir), '/ready_', barrier_seq, '_', rank
         open(newunit=u, file=trim(f), status='replace', action='write', iostat=ios)
         if (ios /= 0) call die("cannot create ready file " // trim(f))
         write(u,'(A)') 'ready'
         close(u)
         tw = wall()
         do
            all_ready = .true.
            do r = 0, nranks - 1
               write(f,'(A,A,I0,A,I0)') trim(ready_dir), '/ready_', barrier_seq, '_', r
               inquire(file=trim(f), exist=ex)
               if (.not. ex) then
                  all_ready = .false.
                  exit
               endif
            enddo
            if (all_ready) exit
            write(f,'(A,A)') trim(ready_dir), '/ABORT'
            inquire(file=trim(f), exist=ex)
            if (ex) call die("point barrier: another rank failed (ABORT file)")
            if (wall() - tw > ready_timeout) then
               call die("point barrier timeout: not all ranks arrived (a rank died or is far behind)")
            endif
            call pause_ms(5.0_dp)
         enddo
      endif
#endif
      t_release = epoch_now()
   end subroutine point_barrier

   subroutine pause_ms(ms)
      Real (Kind=dp), intent(in) :: ms
      Real (Kind=dp) :: t0
      t0 = wall()
      do while ((wall() - t0)*1.d3 < ms)
      enddo
   end subroutine pause_ms

   !-----------------------------------------------------------------------------------------------------------
   subroutine print_header()
      write(output_unit,'(A)') "# blas_bench: BLAS ceiling microbenchmark (testsuite/lowrank/blas_bench.F90)"
      write(output_unit,'(A)') "# SYNC rank start_epoch release_epoch lag_s"
      write(output_unit,'(A)') "# TIMERCOST rank threads what calls_per_sample samples t_median_ns t_min_ns t_max_ns"
      write(output_unit,'(A)') "# RESULT rank threads kernel N k reps inner t_median_s t_min_s t_max_s gflops gflops_best gbs t_release nbuf"
      write(output_unit,'(A)') "META compiler=" // trim(compiler_version())
      write(output_unit,'(A)') "META compiler_options=" // trim(compiler_options())
      write(output_unit,'(A)') "META clock=" // trim(clock_name())
#ifdef BLAS_BENCH_MPI
      write(output_unit,'(A)') "META mpi=yes (MPI_Init, MPI_Barrier before every point)"
#else
      write(output_unit,'(A)') "META mpi=no (point barrier through ready files, NRANKS > 1)"
#endif
#ifdef ALF_TIMERS
      write(output_unit,'(A)') "META timers_compiled=yes (-DALF_TIMERS)"
#else
      write(output_unit,'(A)') "META timers_compiled=no (timer calls are empty: timer cost is not meaningful)"
#endif
      write(output_unit,'(A,6(1X,I0))') "META reps_warmup_rank_threads_nranks_timerreps=", reps, warmup, rank, threads, &
           & nranks, timer_reps
      write(output_unit,'(A,F8.3,1X,A,F10.1)') "META mintime_s=", mintime, "cache_mb=", cache_mb
      flush(output_unit)
   end subroutine print_header

   Character (len=40) function clock_name()
#ifdef _OPENMP
      clock_name = "omp_get_wtime"
#else
      clock_name = "system_clock"
#endif
   end function clock_name

   !-----------------------------------------------------------------------------------------------------------
   Real (Kind=dp) function median(x, n)
      Integer, intent(in) :: n
      Real (Kind=dp), intent(in) :: x(n)
      Real (Kind=dp) :: y(n), tmp
      Integer :: i, j
      y = x
      do i = 2, n
         tmp = y(i)
         j = i - 1
         do while (j >= 1)
            if (y(j) <= tmp) exit
            y(j+1) = y(j)
            j = j - 1
         enddo
         y(j+1) = tmp
      enddo
      if (mod(n, 2) == 1) then
         median = y((n+1)/2)
      else
         median = 0.5_dp*(y(n/2) + y(n/2+1))
      endif
   end function median

   !-----------------------------------------------------------------------------------------------------------
   !> cost of the T1 timers (and of the clocks they read) in ns per call
   subroutine timer_cost()
      Integer, parameter :: NS_SAMPLES = 7
      Real (Kind=dp) :: s(NS_SAMPLES), t0, t1, sink
      Integer (Kind=int64) :: c, csum
      Integer :: i, j, nrep, what
      Character (len=12) :: wname(4) = [Character(len=12) :: "push_pop", "switch", "system_clock", "wtime"]
      nrep = timer_reps
      do what = 1, 4
         do i = 1, NS_SAMPLES + 1            ! first sample is a warm-up, dropped
            sink = 0
            csum = 0
            select case (what)
            case (1)
               call timer_sweep_begin()
               t0 = wall()
               do j = 1, nrep
                  call timer_push(T_UPDATE)
                  call timer_pop(T_UPDATE)
               enddo
               t1 = wall()
               call timer_sweep_end()
            case (2)
               call timer_sweep_begin()
               call timer_push(T_UPDATE)
               t0 = wall()
               do j = 1, nrep
                  call timer_switch(T_UPDATE, T_FLUSH)
                  call timer_switch(T_FLUSH, T_UPDATE)
               enddo
               t1 = wall()
               call timer_pop(T_UPDATE)
               call timer_sweep_end()
            case (3)
               t0 = wall()
               do j = 1, nrep
                  call system_clock(c)
                  csum = csum + iand(c, 1_int64)
               enddo
               t1 = wall()
               sink = dble(csum)
            case (4)
               ! the clock alone (omp_get_wtime; system_clock in a build without OpenMP): the timing calls around the
               ! loop are two further calls, negligible for nrep >= 1e4
               t0 = wall()
               do j = 1, nrep
#ifdef _OPENMP
                  sink = sink + omp_get_wtime()
#else
                  call system_clock(c)
                  sink = sink + dble(c)
#endif
               enddo
               t1 = wall()
            end select
            if (sink < 0) write(error_unit,*) sink
            if (i > 1) s(i-1) = (t1 - t0)*1.d9/dble(nrep)
            if (what == 2 .and. i > 1) s(i-1) = s(i-1)*0.5_dp          ! two switches per iteration
         enddo
         write(output_unit,'(A,1X,I0,1X,I0,1X,A,1X,I0,1X,I0,3(1X,F10.2))') "TIMERCOST", rank, threads, trim(wname(what)), &
              & nrep, NS_SAMPLES, median(s, NS_SAMPLES), minval(s), maxval(s)
         flush(output_unit)
      enddo
   end subroutine timer_cost

   !-----------------------------------------------------------------------------------------------------------
   !> thread-safe, deterministic fill (xorshift64 per column), parallel with a static schedule: first touch by the threads
   subroutine fill(a)
      Complex (Kind=dp), intent(out) :: a(:,:)
      Integer :: i, j
      Integer (Kind=int64) :: s
      Real (Kind=dp) :: re, im
      Real (Kind=dp), parameter :: sc = 1.0_dp/9007199254740992.0_dp
!$omp parallel do schedule(static) private(i, s, re, im)
      do j = 1, size(a,2)
         s = ieor(int(j, int64)*2654435761_int64, 88172645463325252_int64)
         do i = 1, size(a,1)
            s = ieor(s, shiftl(s, 13)); s = ieor(s, shiftr(s, 7)); s = ieor(s, shiftl(s, 17))
            re = dble(shiftr(s, 11))*sc - 0.5_dp
            s = ieor(s, shiftl(s, 13)); s = ieor(s, shiftr(s, 7)); s = ieor(s, shiftl(s, 17))
            im = dble(shiftr(s, 11))*sc - 0.5_dp
            a(i,j) = cmplx(re, im, kind=dp)
         enddo
      enddo
   end subroutine fill

   subroutine fill_vec(v)
      Complex (Kind=dp), intent(out) :: v(:)
      Complex (Kind=dp) :: m(size(v), 1)
      call fill(m)
      v = m(:,1)
   end subroutine fill_vec

   subroutine free_state()
      if (allocated(a)) deallocate(a)
      if (allocated(b)) deallocate(b)
      if (allocated(c)) deallocate(c)
      if (allocated(x)) deallocate(x)
      if (allocated(y)) deallocate(y)
      if (allocated(work)) deallocate(work)
      if (allocated(tau)) deallocate(tau)
      if (allocated(rwork)) deallocate(rwork)
      if (allocated(jpvt)) deallocate(jpvt)
   end subroutine free_state

   subroutine bench_point(ikern, n, k)
      Integer, intent(in) :: ikern, n, k
      Complex (Kind=dp), allocatable :: master(:,:)
      Complex (Kind=dp) :: wq(1)
      Real (Kind=dp), allocatable :: t(:)
      Real (Kind=dp) :: flops, bytes, t0, t1, tmed, tmin, tmax, tcal, nrm
      Integer (Kind=int64) :: matbytes
      Integer :: ib, inner, i, r, info, ierr, nwarm
      Logical :: refresh

      if (n < 1) call die("N < 1")
      cur_kern = ikern; cur_n = n; cur_k = k
      call free_state()
      allocate(t(reps))
      ic = 0
      nbuf = 1
      refresh = .false.
      select case (ikern)
      case (1, 2)
         allocate(a(n,k,1), b(k,n), c(n,n), stat=ierr)
         if (ierr /= 0) call die("out of memory for zgemm")
         call fill(a(:,:,1)); call fill(b); call fill(c)
         flops = 8.0_dp*dble(n)*dble(n)*dble(k)
         bytes = 16.0_dp*(2.0_dp*dble(n)*dble(k) + 2.0_dp*dble(n)*dble(n))
      case (3, 4)
         matbytes = 16_int64*int(n, int64)*int(n, int64)
         nbuf = max(1, int(ceiling(cache_mb*1048576.0_dp/dble(matbytes))))
         allocate(a(n,n,nbuf), x(n), y(n), stat=ierr)
         if (ierr /= 0) call die("out of memory for the rotating buffers of " // trim(kname(ikern)))
         do ib = 1, nbuf
            call fill(a(:,:,ib))
         enddo
         call fill_vec(x); call fill_vec(y)
         alpha = cmplx(1.0d-6, 0.0d0, kind=dp)
         flops = 8.0_dp*dble(n)*dble(n)
         if (ikern == 3) then
            bytes = 16.0_dp*(2.0_dp*dble(n)*dble(n) + 2.0_dp*dble(n))
         else
            bytes = 16.0_dp*(dble(n)*dble(n) + 2.0_dp*dble(n))
         endif
      case (5)
         allocate(a(n,n,1), master(n,n), tau(n), jpvt(n), rwork(2*n), stat=ierr)
         if (ierr /= 0) call die("out of memory for zgeqp3")
         call fill(master)
         a(:,:,1) = master
         jpvt = 0
         call zgeqp3(n, n, a(:,:,1), n, jpvt, tau, wq, -1, rwork, info)
         lwork = max(1, int(real(wq(1))))
         allocate(work(lwork))
         flops = 16.0_dp/3.0_dp*dble(n)**3
         bytes = 32.0_dp*dble(n)*dble(n)
         refresh = .true.
      end select

      call point_barrier()

      nwarm = warmup + 1
      if (refresh) nwarm = 1
      tcal = 0
      do i = 1, nwarm                        ! warm-up; the last call calibrates the batch size
         if (refresh) a(:,:,1) = master
         t0 = wall()
         call one_call()
         t1 = wall()
         tcal = t1 - t0
      enddo
      inner = 1
      if (.not. refresh) inner = max(1, int(ceiling(mintime/max(tcal, 1.0d-9))))

      do r = 1, reps
         if (refresh) a(:,:,1) = master
         t0 = wall()
         do i = 1, inner
            call one_call()
         enddo
         t1 = wall()
         t(r) = (t1 - t0)/dble(inner)
      enddo
      tmed = median(t, reps); tmin = minval(t); tmax = maxval(t)

      ! guard against a kernel that silently did nothing or produced garbage
      select case (ikern)
      case (1, 2);  nrm = abs(c(1,1)) + abs(c(n,n))
      case (3);     nrm = abs(a(1,1,1)) + abs(a(n,n,1))
      case (4);     nrm = abs(y(1)) + abs(y(n))
      case (5);     nrm = abs(a(1,1,1))
      end select
      if (.not. (nrm == nrm) .or. nrm > huge(nrm)) call die("non-finite result in " // trim(kname(ikern)))

      write(output_unit,'(A,2(1X,I0),1X,A,4(1X,I0),3(1X,ES13.6),2(1X,F10.3),1X,F10.3,1X,F16.3,1X,I0)') "RESULT", rank, &
           & threads, trim(kname(ikern)), n, k, reps, inner, tmed, tmin, tmax, flops/tmed*1.d-9, flops/tmin*1.d-9, &
           & bytes/tmed*1.d-9, t_release, nbuf
      flush(output_unit)
      deallocate(t)
      call free_state()
   end subroutine bench_point

   subroutine one_call()
      Integer :: info
      select case (cur_kern)
      case (1, 2)
         call zgemm('N', 'N', cur_n, cur_n, cur_k, mone, a(:,:,1), cur_n, b, cur_k, one, c, cur_n)
      case (3)
         ic = ic + 1; if (ic > nbuf) ic = 1
         call zgeru(cur_n, cur_n, alpha, x, 1, y, 1, a(:,:,ic), cur_n)
      case (4)
         ic = ic + 1; if (ic > nbuf) ic = 1
         call zgemv('N', cur_n, cur_n, one, a(:,:,ic), cur_n, x, 1, zero, y, 1)
      case (5)
         jpvt = 0
         call zgeqp3(cur_n, cur_n, a(:,:,1), cur_n, jpvt, tau, work, lwork, rwork, info)
         if (info /= 0) call die("zgeqp3 info /= 0")
      end select
   end subroutine one_call

   !-----------------------------------------------------------------------------------------------------------
   !> tiny ZGEMM and ZGEQP3 (the MKL_VERBOSE probe of campaign/job_template.sh): "probe ok", exit 0
   subroutine run_probe()
      Integer, parameter :: m = 8
      Complex (Kind=dp) :: aa(m,m), bb(m,m), cc(m,m), tt(m), ww(4*m)
      Real (Kind=dp) :: rw(2*m)
      Integer :: jp(m), info
      call fill(aa); call fill(bb)
      cc = zero
      call zgemm('N', 'N', m, m, m, one, aa, m, bb, m, zero, cc, m)
      if (.not. (sum(abs(cc)) > 0.0_dp)) then
         write(output_unit,'(A)') "probe FAILED (zgemm)"; error stop 1
      endif
      jp = 0
      call zgeqp3(m, m, aa, m, jp, tt, ww, 4*m, rw, info)
      if (info /= 0 .or. .not. (sum(abs(aa)) > 0.0_dp)) then
         write(output_unit,'(A)') "probe FAILED (zgeqp3)"; error stop 1
      endif
      write(output_unit,'(A)') "probe ok"
   end subroutine run_probe

end program blas_bench
