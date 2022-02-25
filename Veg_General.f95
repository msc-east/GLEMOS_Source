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

#ifdef M_VEG
Module Veg_General
    
use GeneralParams
use Veg_Params
use Veg_POP_Params
use Veg_POP_General
use Exch_Params

use Veg_Input
use Veg_Output

use GeneralOutput

implicit none

contains

!==============================================================
! Vegetation memory allocation/deallocation
!==============================================================
subroutine Veg_Memory(operation)

	character(*) operation

	selectcase(operation)
!---------------------------------------------------------------------------------
	case('Initial')
	  call Veg_MemAlloc
          call Veg_POP_MemAlloc
!---------------------------------------------------------------------------------
	case('Final')

	  call Veg_MemDealloc
          call Veg_POP_MemDealloc
	endselect

end subroutine Veg_Memory

!==============================================================
! Vegetation processes
!==============================================================
subroutine Veg_Process	

  real dTime, s, s1
  integer NSubs

#ifdef DEBUG_MODE
    print *, '+Entering Veg_Process... '
#endif

!   call Veg_Interaction     ! Subroutine describing interaction between pollutants in vegetation

  do NSubs = 1, NumSubs
	select case (SubsGroup(NSubs))
	  case (POP)                         
            call Veg_POPProcess(NSubs,Tstep(Veg))
	  case default
	    continue
	end select
  end do

#ifdef DEBUG_MODE
    print *, '+Exit Veg_Process'
#endif
end subroutine Veg_Process

!==============================================================
! Vegetation time step
!==============================================================
subroutine Veg_TStepCalc(dTiter)

  real dTiter		
  integer NSubs

  dTiter = dTInput
  do NSubs = 1, NumSubs
	select case (SubsGroup(NSubs))
	  case (POP)                         
		call Veg_POPTStepCalc(NSubs,dTiter)
	  case default
		continue
	end select
  end do

end subroutine Veg_TStepCalc

!==============================================================
! Vegetation exchange
!==============================================================
subroutine Veg_ExchPar           ! To be called each MP

  integer NSubs

  Veg_VExch = 0.
  do NSubs = 1, NumSubs
	select case (SubsGroup(NSubs))
	  case (POP)                         
!		call Veg_POPExchPar(NSubs)
	  case default
		continue
	end select
  end do

end subroutine Veg_ExchPar

!==============================================================
! Vegetation mass
!==============================================================
real(8) function Veg_Mass(NSubs)

  integer NSubs

  select case (SubsGroup(NSubs))
	case (POP)                          
	  Veg_Mass = Veg_POPMass(NSubs)
	case default
	  continue
  end select

end function Veg_Mass

!==============================================================
! The following subroutine should describe interactions 
! in vegetation between different pollutants includded in the model run
!==============================================================
subroutine Veg_Interaction

  continue

end subroutine Veg_Interaction


end Module Veg_General
#endif