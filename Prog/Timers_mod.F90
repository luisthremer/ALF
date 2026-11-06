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
!  New file of the delayed/submatrix research branch (feature/delayed-submatrix, 2026): benchmark timing split.

!--------------------------------------------------------------------
!> @brief
!> Exclusive (stack) timing split of the sweep, for benchmarks (delayed/submatrix branch).
!>
!> The module is compiled in every build; its routines only do something with -DALF_TIMERS
!> (otherwise they are empty). Every call site and every `use Timers_mod` in the program is under
!> #ifdef ALF_TIMERS, so the default build does not call it at all.
!>
!> Exclusive timers: a stack of category ids. timer_push(id) charges the time since the last
!> transition to the current top and makes id the top; timer_pop(id) (id must be the top) charges
!> id and resumes its parent; timer_switch(a, b) replaces the top a by b (one clock read, used per
!> vertex between the wraps and the update). Each transition reads the clock once and charges the
!> integer count difference to exactly one id, so the charges between timer_sweep_begin and
!> timer_sweep_end telescope: sum over all ids (the sweep frame's own exclusive time is "other")
!> equals the inclusive sweep time exactly (asserted at every bin end). Outside a sweep (between
!> timer_sweep_end and the next timer_sweep_begin) push/pop/switch are no-ops.
!>
!> Categories: T_UPDATE (ham%S0, proposal ranf, Upgrade2/Upgrade_lowrank without nested flushes),
!> T_FLUSH (Lowrank_flush with something pending), T_VWRAP (Op_Wrapup/Op_Wrapdo of the vertices and the
!> low-rank hooks), T_HOP (hopping wraps incl. cross-slice hops), T_STAB (WRAPUR/WRAPUL/CGR, precision
!> control, phase), T_MEAS (Obser, Hop_mod_Symm, reconstruction, Tau_m/Tau_p), T_GLOBAL (global moves,
!> Global_tau, tempering exchange), T_LANGEVIN_HMC (whole Langevin/HMC branches), T_SWEEP (the sweep
!> frame: exclusive = "other", inclusive = the whole body of the DO NSW loop).
!>
!> Counters (sequential updates): proposals, accepts; flushes with something pending by reason
!> (timer_flush_reason, set by the caller right before calling Lowrank_flush and consumed by it):
!> capacity, slice end, ineligible vertex, read point (cross-slice key); mean pending rank at a
!> flush (max over flavours of the pending columns / sites).
!>
!> Output (timer_bin_end, collective on the group communicator, group rank 0 writes; one block per
!> bin, appended to `timers.txt`, Temp_<g>/timers.txt in tempering builds), seconds:
!>   TIMER_RATE <count_rate>  TIMER_CLOCK <system_clock|omp_get_wtime>  TIMER_RANKS <n>
!>   TIMER_TICK <measured clock resolution in s>                                      (once per run)
!>   TIMER_SWEEPS <bin> <n>
!>   TIMER <bin> <id> <rank0> <max> <mean>    id = update flush vwrap hop stab meas global
!>                                            langevin_hmc other sweep bin_overhead
!>   COUNTER <bin> <name> <rank0> <max> <mean>  (counts: rank0, max integers, mean ES; flush_rank_mean: ES)
!> bin_overhead = bin wall time (timer_bin_begin .. timer_bin_end: Init_obs, sweeps, Pr_obs, output)
!> minus the inclusive sweep times. Timers and counters are reset at every bin start.
!--------------------------------------------------------------------
module Timers_mod
   use iso_fortran_env, only: int64, error_unit
#ifdef ALF_TIMERS
!$ use omp_lib, only: omp_get_wtime
#endif
   implicit none
   private

   Integer, parameter, public :: T_UPDATE = 1, T_FLUSH = 2, T_VWRAP = 3, T_HOP = 4, T_STAB = 5, &
        &                        T_MEAS = 6, T_GLOBAL = 7, T_LANGEVIN_HMC = 8, T_SWEEP = 9
   Integer, parameter, public :: N_TIMERS = 9
   !> Flush reasons (timer_flush_reason)
   Integer, parameter, public :: FR_NONE = 0, FR_CAPACITY = 1, FR_SLICE_END = 2, FR_INELIGIBLE = 3, &
        &                        FR_READ_POINT = 4
   !> Reason of the next Lowrank_flush call; set by the caller, reset to FR_NONE by the flush
   Integer, public :: timer_flush_reason = FR_NONE

   public :: timer_push, timer_pop, timer_switch, timer_sweep_begin, timer_sweep_end, &
        &    timer_bin_begin, timer_bin_end, timer_flush_begin, timer_flush_end, timer_flush_skip, &
        &    timer_count_proposal

#ifdef ALF_TIMERS
   !> Output names: the N_TIMERS exclusive timers (T_SWEEP's exclusive time is "other"), then the inclusive
   !> sweep time and the bin overhead
   Character (len=12), parameter :: out_name(N_TIMERS+2) = [Character(len=12) :: "update", "flush", &
        &  "vwrap", "hop", "stab", "meas", "global", "langevin_hmc", "other", "sweep", "bin_overhead"]
   Integer, parameter :: MAX_DEPTH = 16
   Logical        :: initialised = .false., in_sweep = .false., use_omp = .false., header_done = .false.
   Integer (Kind=int64) :: count_rate = 0
   Integer (Kind=int64) :: tick = 0                ! measured clock resolution in counts (smallest step seen)
   Integer        :: depth = 0
   Integer        :: stk(0:MAX_DEPTH) = 0
   Integer (Kind=int64) :: t_last = 0, t_sweep0 = 0, t_bin0 = 0
   Integer (Kind=int64) :: excl(N_TIMERS) = 0      ! exclusive counts per id (T_SWEEP: "other")
   Integer (Kind=int64) :: sweep_incl = 0          ! inclusive sweep counts in this bin
   Integer (Kind=int64) :: n_sweeps = 0
   ! counters
   Integer, parameter :: C_PROP = 1, C_ACC = 2, C_FCAP = 3, C_FSLICE = 4, C_FINEL = 5, C_FREAD = 6, &
        &                C_FNONE = 7, C_FTOT = 8, N_COUNTERS = 8
   Character (len=18), parameter :: counter_name(N_COUNTERS+1) = [Character(len=18) :: "proposals", &
        &  "accepts", "flush_capacity", "flush_slice_end", "flush_ineligible", "flush_read_point", &
        &  "flush_unattributed", "flush_total", "flush_rank_mean"]
   Integer (Kind=int64) :: cnt(N_COUNTERS) = 0
   Integer (Kind=int64) :: flush_rank_sum = 0
#endif

contains

#ifdef ALF_TIMERS
   !> Clock in counts (count_rate per second)
   Integer (Kind=int64) function tnow()
      if (use_omp) then
         tnow = 0
!$       tnow = int(omp_get_wtime()*1.d9, int64)
      else
         call system_clock(tnow)
      endif
   end function tnow

   subroutine timers_init()
      Integer (Kind=int64) :: r, t0, t1
      Integer              :: i
      call system_clock(count_rate=r)
      count_rate = r
      use_omp = .false.
      if (r < 1000000_int64) then
!$       use_omp = .true.
!$       count_rate = 1000000000_int64
      endif
      ! Effective resolution (e.g. 1 us for system_clock on macOS although count_rate = 1e9): smallest
      ! nonzero step of consecutive reads. Per-segment rounding is statistically unbiased; reported only.
      tick = huge(tick)
      do i = 1, 20
         t0 = tnow()
         do
            t1 = tnow()
            if (t1 /= t0) exit
         enddo
         tick = min(tick, t1 - t0)
      enddo
      initialised = .true.
   end subroutine timers_init

   subroutine timer_fail(msg)
      Character (len=*), intent(in) :: msg
      write(error_unit,'(A)') ' Timers_mod: ' // trim(msg)
      error stop 1
   end subroutine timer_fail

   Character (len=12) function idname(id)
      Integer, intent(in) :: id
      if (id >= 1 .and. id < N_TIMERS) then
         idname = out_name(id)
      elseif (id == T_SWEEP) then
         idname = "sweep"
      else
         write(idname,'(I0)') id
      endif
   end function idname
#endif

   !> Pauses the current top and starts id (no-op outside a sweep)
   subroutine timer_push(id)
      Integer, intent(in) :: id
#ifdef ALF_TIMERS
      Integer (Kind=int64) :: t
      if (.not. in_sweep) return
      if (depth >= MAX_DEPTH) call timer_fail("timer_push: stack overflow")
      t = tnow()
      excl(stk(depth)) = excl(stk(depth)) + (t - t_last)
      depth = depth + 1
      stk(depth) = id
      t_last = t
#endif
   end subroutine timer_push

   !> Stops id (must be the top) and resumes its parent (no-op outside a sweep)
   subroutine timer_pop(id)
      Integer, intent(in) :: id
#ifdef ALF_TIMERS
      Integer (Kind=int64) :: t
      if (.not. in_sweep) return
      if (depth <= 1 .or. stk(depth) /= id) &
           & call timer_fail("timer_pop: " // trim(idname(id)) // " is not the top (top: " // &
           &                 trim(idname(stk(depth))) // ")")
      t = tnow()
      excl(id) = excl(id) + (t - t_last)
      depth = depth - 1
      t_last = t
#endif
   end subroutine timer_pop

   !> Replaces the top id_from by id_to (one clock read; no-op outside a sweep)
   subroutine timer_switch(id_from, id_to)
      Integer, intent(in) :: id_from, id_to
#ifdef ALF_TIMERS
      Integer (Kind=int64) :: t
      if (.not. in_sweep) return
      if (depth <= 1 .or. stk(depth) /= id_from) &
           & call timer_fail("timer_switch: " // trim(idname(id_from)) // " is not the top (top: " // &
           &                 trim(idname(stk(depth))) // ")")
      t = tnow()
      excl(id_from) = excl(id_from) + (t - t_last)
      stk(depth) = id_to
      t_last = t
#endif
   end subroutine timer_switch

   !> Opens the sweep frame (the stack must be empty)
   subroutine timer_sweep_begin()
#ifdef ALF_TIMERS
      if (.not. initialised) call timers_init()
      if (in_sweep .or. depth /= 0) call timer_fail("timer_sweep_begin: stack not empty")
      t_last   = tnow()
      t_sweep0 = t_last
      depth    = 1
      stk(1)   = T_SWEEP
      in_sweep = .true.
#endif
   end subroutine timer_sweep_begin

   !> Closes the sweep frame; asserts that only the sweep frame is open (stack depth 0 afterwards)
   subroutine timer_sweep_end()
#ifdef ALF_TIMERS
      Integer (Kind=int64) :: t
      if (.not. in_sweep .or. depth /= 1 .or. stk(1) /= T_SWEEP) &
           & call timer_fail("timer_sweep_end: timers still open at the end of the sweep (top: " // &
           &                 trim(idname(stk(depth))) // ")")
      t = tnow()
      excl(T_SWEEP) = excl(T_SWEEP) + (t - t_last)
      sweep_incl    = sweep_incl + (t - t_sweep0)
      depth    = 0
      in_sweep = .false.
      n_sweeps = n_sweeps + 1
#endif
   end subroutine timer_sweep_end

   !> Resets timers and counters; the bin clock starts (call before Init_obs)
   subroutine timer_bin_begin()
#ifdef ALF_TIMERS
      if (.not. initialised) call timers_init()
      if (in_sweep) call timer_fail("timer_bin_begin: inside a sweep")
      excl = 0
      sweep_incl = 0
      n_sweeps = 0
      cnt = 0
      flush_rank_sum = 0
      t_bin0 = tnow()
#endif
   end subroutine timer_bin_begin

   !> Called by Lowrank_flush when something is pending: counts the flush by timer_flush_reason
   !> (then reset), adds the pending rank and starts T_FLUSH
   subroutine timer_flush_begin(rank)
      Integer, intent(in) :: rank
#ifdef ALF_TIMERS
      if (in_sweep) then
         select case (timer_flush_reason)
         case (FR_CAPACITY);   cnt(C_FCAP)   = cnt(C_FCAP)   + 1
         case (FR_SLICE_END);  cnt(C_FSLICE) = cnt(C_FSLICE) + 1
         case (FR_INELIGIBLE); cnt(C_FINEL)  = cnt(C_FINEL)  + 1
         case (FR_READ_POINT); cnt(C_FREAD)  = cnt(C_FREAD)  + 1
         case default;         cnt(C_FNONE)  = cnt(C_FNONE)  + 1
         end select
         cnt(C_FTOT) = cnt(C_FTOT) + 1
         flush_rank_sum = flush_rank_sum + rank
      endif
      timer_flush_reason = FR_NONE
      call timer_push(T_FLUSH)
#endif
   end subroutine timer_flush_begin

   !> End of a flush started with timer_flush_begin
   subroutine timer_flush_end()
#ifdef ALF_TIMERS
      call timer_pop(T_FLUSH)
#endif
   end subroutine timer_flush_end

   !> Lowrank_flush with nothing pending: drops the caller's reason
   subroutine timer_flush_skip()
#ifdef ALF_TIMERS
      timer_flush_reason = FR_NONE
#endif
   end subroutine timer_flush_skip

   !> One sequential proposal (Upgrade2 / Upgrade_lowrank called); acc = accepted
   subroutine timer_count_proposal(acc)
      Logical, intent(in) :: acc
#ifdef ALF_TIMERS
      if (.not. in_sweep) return
      cnt(C_PROP) = cnt(C_PROP) + 1
      if (acc) cnt(C_ACC) = cnt(C_ACC) + 1
#endif
   end subroutine timer_count_proposal

   !> End of bin nbc: checks the sum rule, reduces over the group communicator comm (collective:
   !> every rank of the group must call it) and appends the block to file_name on group rank 0.
   subroutine timer_bin_end(nbc, comm, file_name)
#if defined(ALF_TIMERS) && defined(MPI)
      use mpi
#endif
      Integer,           intent(in) :: nbc, comm
      Character (len=*), intent(in) :: file_name
#ifdef ALF_TIMERS
      Integer, parameter   :: NV = N_TIMERS + 2 + N_COUNTERS + 1
      Real (Kind=Kind(0.d0)) :: v(NV), vsum(NV), vmax(NV), sec
      Integer (Kind=int64) :: t_bin, tsum
      Integer              :: irank, isize, i, unit_t
#ifdef MPI
      Integer              :: ierr
#endif

      t_bin = tnow() - t_bin0
      if (in_sweep .or. depth /= 0) call timer_fail("timer_bin_end: inside a sweep")
      ! Sum rule (exact in integer counts): categories + other = inclusive sweep time
      tsum = sum(excl)
      if (tsum /= sweep_incl) call timer_fail("timer_bin_end: sum of the exclusive timers /= sweep time")

      sec = 1.d0/dble(count_rate)
      ! v: update .. langevin_hmc, other (= exclusive T_SWEEP), sweep (inclusive), bin_overhead, counters, rank mean
      do i = 1, N_TIMERS
         v(i) = dble(excl(i))*sec
      enddo
      v(N_TIMERS+1) = dble(sweep_incl)*sec
      v(N_TIMERS+2) = dble(t_bin - sweep_incl)*sec
      do i = 1, N_COUNTERS
         v(N_TIMERS+2+i) = dble(cnt(i))
      enddo
      v(NV) = 0.d0
      if (cnt(C_FTOT) > 0) v(NV) = dble(flush_rank_sum)/dble(cnt(C_FTOT))

      irank = 0
      isize = 1
#ifdef MPI
      call MPI_Comm_rank(comm, irank, ierr)
      call MPI_Comm_size(comm, isize, ierr)
      call MPI_Reduce(v, vsum, NV, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      call MPI_Reduce(v, vmax, NV, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
#else
      vsum = v
      vmax = v
#endif
      if (irank /= 0) return

      open(newunit=unit_t, file=file_name, status="unknown", position="append", action="write")
      if (.not. header_done) then
         write(unit_t,'(A)') "# ALF timing split (-DALF_TIMERS): exclusive stack timers, seconds per bin;"
         write(unit_t,'(A)') "# TIMER <bin> <id> <rank0> <max> <mean>, COUNTER <bin> <name> <rank0> <max> <mean>;"
         write(unit_t,'(A)') "# update+flush+vwrap+hop+stab+meas+global+langevin_hmc+other = sweep (per rank, exact)"
         write(unit_t,'(A,1X,I0)') "TIMER_RATE", count_rate
         if (use_omp) then
            write(unit_t,'(A)') "TIMER_CLOCK omp_get_wtime"
         else
            write(unit_t,'(A)') "TIMER_CLOCK system_clock"
         endif
         write(unit_t,'(A,1X,I0)') "TIMER_RANKS", isize
         write(unit_t,'(A,1X,ES12.5)') "TIMER_TICK", dble(tick)*sec
         header_done = .true.
      endif
      write(unit_t,'(A,1X,I0,1X,I0)') "TIMER_SWEEPS", nbc, n_sweeps
      do i = 1, N_TIMERS + 2
         write(unit_t,'(A,1X,I0,1X,A,3(1X,ES16.9))') "TIMER", nbc, trim(out_name(i)), v(i), vmax(i), &
              &  vsum(i)/dble(isize)
      enddo
      do i = 1, N_COUNTERS
         write(unit_t,'(A,1X,I0,1X,A,2(1X,I0),1X,ES16.9)') "COUNTER", nbc, trim(counter_name(i)), &
              &  nint(v(N_TIMERS+2+i), int64), nint(vmax(N_TIMERS+2+i), int64), vsum(N_TIMERS+2+i)/dble(isize)
      enddo
      write(unit_t,'(A,1X,I0,1X,A,3(1X,ES16.9))') "COUNTER", nbc, trim(counter_name(N_COUNTERS+1)), v(NV), vmax(NV), &
           &  vsum(NV)/dble(isize)
      close(unit_t)
#endif
   end subroutine timer_bin_end

end module Timers_mod
