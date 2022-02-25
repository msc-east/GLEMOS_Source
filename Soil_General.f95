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

#ifdef M_SOIL
Module Soil_General
    
use GeneralParams
use Soil_Params
use Soil_POP_Params
use Soil_POP_General
use Exch_Params

use Soil_Input
use Soil_Output

use GeneralOutput

implicit none

contains

!==============================================================
! Soil memory allocation/deallocation
!==============================================================
subroutine Soil_Memory(operation)

	character(*) operation

	selectcase(operation)
!---------------------------------------------------------------------------------
	case('Initial')
	  call Soil_MemAlloc
          call Soil_POP_MemAlloc
!---------------------------------------------------------------------------------
	case('Final')
	  call Soil_MemDealloc
          call Soil_POP_MemDealloc
	endselect

end subroutine Soil_Memory

!==============================================================
! Soil processes
!==============================================================
subroutine Soil_Process	

  real dTime, s, s1
  integer NSubs

#ifdef DEBUG_MODE
    print *, '+Entering Soil_Process... '
#endif

   call Soil_Emission        ! 09.04.2012 Gusev
!   call Soil_Interaction     ! Subroutine describing interaction between pollutants in soil

  do NSubs = 1, NumSubs
	select case (SubsGroup(NSubs))
	  case (POP)                         
            call Soil_POPProcess(NSubs,Tstep(Soil))
	  case default
	    continue
	end select
  end do

#ifdef DEBUG_MODE
    print *, '+Exit Soil_Process'
#endif
end subroutine Soil_Process

!==============================================================
! Soil time step
!==============================================================
subroutine Soil_TStepCalc(dTiter)

  real dTiter		
  integer NSubs

  dTiter = dTInput
  do NSubs = 1, NumSubs
	select case (SubsGroup(NSubs))
	  case (POP)                          
		call Soil_POPTStepCalc(NSubs,dTiter)
	  case default
		continue
	end select
  end do

end subroutine Soil_TStepCalc

!==============================================================
! Soil exchange
!==============================================================
subroutine Soil_ExchPar           ! To be called each MP

  integer NSubs

  Soil_VExch = 0.
  do NSubs = 1, NumSubs
	select case (SubsGroup(NSubs))
	  case (POP)                          
		call Soil_POPExchPar(NSubs)
	  case default
		continue
	end select
  end do

end subroutine Soil_ExchPar

!==============================================================
! Soil mass
!==============================================================
real(8) function Soil_Mass(NSubs)

  integer NSubs

  select case (SubsGroup(NSubs))
	case (POP)                          
	  Soil_Mass = Soil_POPMass(NSubs)
	case default
	  continue
  end select

end function Soil_Mass

!==============================================================
! The following subroutine should describe interactions 
! in soil between different pollutants includded in the model run
!==============================================================
subroutine Soil_Interaction

  continue

end subroutine Soil_Interaction


!==============================================================
! Enter pollutant emission directly to soil
! Added 09.04.2012
!==============================================================
subroutine Soil_Emission

	integer i, j, k, L, Ind, Form
	real    dT
        real*8 emis,landsum
        real*8  cmass(120,60),Mass1,Mass2,Emis1,Emis2

    dT=Tstep(Soil)
    do j=Jmin, Jmax
        do i=minI(j), maxI(j)
            do Ind=1, SoilEmisNum
! Anthropogenic emission
                Form=SoilEmisInd(Ind)
                emis = SoilEmisFlux(i,j,Form)*dt
                landsum = sum(Soil_Frac(i,j,1:Soil_NumType))
                if(landsum > 0.) then
                        do L = 1, Soil_NumType
                                do K = 1, Soil_KMax
                                        Soil_Conc(i,j,K,L,Form) = Soil_Conc(i,j,K,L,Form) + &
                                        & emis*soilprof(K)/landsum/dble(MeshArea(i,j))/dble((Poros-wl))/dble(Soil_dz(K))
                              end do
                       end do

                MassSoilEmis(Form)=MassSoilEmis(Form)+emis  ! Keep total amount for balance
                end if
           end do
        end do
    end do

end subroutine Soil_Emission

end Module Soil_General
#endif