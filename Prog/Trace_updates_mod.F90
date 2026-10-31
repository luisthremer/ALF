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
!  New file of the delayed/submatrix research branch (feature/delayed-submatrix, 2026). See docs/algorithm-notes.md.

!> @brief
!> Debug trace of the Metropolis decisions of the sequential local update.
!>
!> The module is compiled in every build, but Trace_update_write is only called from
!> Upgrade2 under #ifdef ALF_TRACE_UPDATES. It is used by testsuite/lowrank/compare_trajectory.sh
!> to compare the accept/reject decisions of different Sequential_scheme values.
!> Each call appends one line to trace_updates.txt in the working directory:
!>   nt  n_op  Re(R)  Im(R)  Weight  ranf  toggle        (ES24.16 for the reals, toggle as 0/1)
!> where R = Ratiotot is the full fermion ratio as used for the decision.
module Trace_updates_mod
   implicit none
   private
   public :: Trace_update_write

   integer, parameter :: trace_unit = 97
   logical, save      :: trace_open = .false.

contains

   subroutine Trace_update_write(nt, n_op, Ratiotot, Weight, ranf_val, toggle)
      integer,                      intent(in) :: nt, n_op
      complex (kind=kind(0.d0)),    intent(in) :: Ratiotot
      real    (kind=kind(0.d0)),    intent(in) :: Weight, ranf_val
      logical,                      intent(in) :: toggle
      integer :: itog

      if ( .not. trace_open ) then
         open(unit=trace_unit, file="trace_updates.txt", status="unknown", position="append")
         trace_open = .true.
      endif
      itog = 0
      if ( toggle ) itog = 1
      write(trace_unit, '(2(I8,1X),4(ES24.16,1X),I2)') nt, n_op, real(Ratiotot, kind(0.d0)), &
           &  aimag(Ratiotot), Weight, ranf_val, itog
   end subroutine Trace_update_write

end module Trace_updates_mod
