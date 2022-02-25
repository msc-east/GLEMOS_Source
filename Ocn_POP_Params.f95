!*************************************************************************!
!*  Copyright (C) Meteorological Synthesizing Centre - East of EMEP, 2017 
!*
!*  Contact information:
!*  Meteorological Synthesizing Centre - East of EMEP
!*  2nd Roshchinsky proezd, 8/5
!*  115419 Moscow, Russia
!*  email: msce@msceast.org
!*  http://www.msceast.org
!*
!*  This program is free software: you can redistribute it and/or modify
!*  it under the terms of the GNU General Public License as published by
!*  the Free Software Foundation, either version 3 of the License, or
!*  (at your option) any later version.
!*
!*  This program is distributed in the hope that it will be useful,
!*  but WITHOUT ANY WARRANTY; without even the implied warranty of
!*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*  GNU General Public License for more details.
!*
!*  You should have received a copy of the GNU General Public License
!*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*************************************************************************!

!=======================================================================================================
! Ocean compartment POP parameters
! Version:      1.0
! Modified:     24.11.2014
!=======================================================================================================
#ifdef G_POP
Module Ocn_POP_Params
use GeneralParams
implicit none

  real CDWater(2,MaxSubs)        ! CDWater(f,p) - degradation coefficient for form f (1 - dissolved, 2 - particulate)
  real k_ph_rd(MaxSubs)          ! k_ph_rd(p) - ratio of concentration in dissolved form to that in particulate form

end Module Ocn_POP_Params
#endif