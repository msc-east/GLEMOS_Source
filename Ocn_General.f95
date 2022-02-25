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

#ifdef M_OCN
module Ocn_General

use GeneralParams
use Atm_POP_Params
use Ocn_Params
use Ocn_POP_Params
use Ocn_POP_General
use Ocn_Input
use Ocn_Output
use Exch_Params

implicit none

 contains

!==============================================================
! Ocean processes
!==============================================================
subroutine Ocn_Process	 

  real dTime
  integer Subs

#ifdef DEBUG_MODE
    print *, '+Entering Ocn_Process... '
#endif

#if (OCEAN_RUN==1)
      call Ocn_Interaction
      call InterpolateInputOceanData
#else
      ! Switched off to run without ocean (Sep 2016 V.E.Sh.)
#endif

  dtOceanAdvDiff=Tstep(Ocn)
  kocean=kocean+1
  avg_ts = .false.
  if (mod(kocean,17).eq.0) then
      avg_ts = .true.
      kocean=0
  end if

  do Subs = 1, NumSubs
	select case (SubsGroup(Subs))
	  case (POP)                         
	    call Ocn_POPProcess(Subs, Tstep(Ocn))	 
	  case default
	    continue
	end select
  end do

#ifdef DEBUG_MODE
    print *, '+Exit Ocn_Process'
#endif
end subroutine Ocn_Process

!==============================================================
! Ocean time step
!==============================================================
subroutine Ocn_TStepCalc(dTiter)

  real dTiter
  integer NSubs

  dTiter=600.

end subroutine Ocn_TStepCalc

!==============================================================
! Ocean exchange (Currently is not called from anywhere !!!!)
!==============================================================
subroutine Ocn_ExchPar          

  integer NSubs

  Ocn_VExch = 0.
  do NSubs = 1, NumSubs
	select case (SubsGroup(NSubs))
	  case (POP)
#ifdef O_NO_EQUIL
		call Ocn_POPExchPar(NSubs)
#endif
#ifdef O_EQUIL
                ! Calculation of exchange between air and water using instanteneous equilibrium (Exch_General_POP)
                ! Code included in the Exch_General_POP.f95
#endif
	  case default
		continue
	end select
  end do

end subroutine Ocn_ExchPar

!==============================================================
! Ocean mass
!==============================================================
real(8) function Ocn_Mass(NSubs)

  integer NSubs

  select case (SubsGroup(NSubs))
	case (POP)                      
	  Ocn_Mass = Ocn_POPMass(NSubs)
	case default
	  continue
  end select

end function Ocn_Mass

!==============================================================
! The following subroutine should describe interactions in ocean between different pollutants includded in the model run
!==============================================================
subroutine Ocn_Interaction

  continue

end subroutine Ocn_Interaction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine managing memory allocation/deallocation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Ocn_Memory(operation)

	character(*) operation

	selectcase(operation)
!---------------------------------------------------------------------------------
	case('Initial')
	  call Ocn_allocation
!---------------------------------------------------------------------------------
	case('Final')
	  call Ocn_deallocation
	endselect

end subroutine Ocn_Memory


end Module Ocn_General
#endif
