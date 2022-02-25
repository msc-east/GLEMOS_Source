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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Deposition module
! Contains subroutines: 
!	DryDepVeloc
!	DryDepos_Part
!	WetDepos		
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

module Atm_Exchange 

  use GeneralParams
  use Exch_Params
#ifdef G_HG
  use Atm_Hg_General              
#endif
#ifdef G_Tracer
  use Atm_Tracer_General
#endif
#ifdef G_HM
  use Atm_HM_General
#endif
#ifdef G_AERO
  use Atm_AERO_General
#endif

  implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating atmospheric exchange parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_ExchPar

    integer Ind, Form, Subs, Grp

      do Ind=1, DryDepNum
        Form=DryDepInd(Ind)
        Subs=FormSubs(Atm,Form)
        Grp=SubsGroup(Subs)

        selectcase(Grp)
#ifdef G_HM
          case(HM)
            call Atm_HM_ExchPar(Ind,toDay)
            call Atm_HM_ExchPar(Ind,toMor)
#endif
#ifdef G_HG
          case(HG)
            call Atm_Hg_ExchPar(Ind,toDay)
            call Atm_Hg_ExchPar(Ind,toMor)
              
#endif
#ifdef G_POP
          case(POP)
#endif
#ifdef G_AERO
            case(AERO)
            call Atm_AERO_ExchPar(Ind)
#endif
#ifdef G_Tracer
          case(Tracer)
            call Atm_Tracer_ExchPar(Ind)
#endif
        endselect
      enddo

end subroutine Atm_ExchPar

end module Atm_Exchange