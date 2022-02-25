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

module Atm_HM_Params

#ifdef G_HM

  use GeneralParams
  use Atm_Params

  implicit none


    real, allocatable :: DustFlux(:,:,:,:,:)            ! Flux of suspended dust [kg/cell]
    real, allocatable :: HMinSoil(:,:,:)                ! HM concentration in soil [kg/kg]
    real EFsea                                          ! Seawate emission factor [kg/kg]
    integer LUresuspNum
    character(20) LUresusp(MaxSurf)                     ! Land cover types with re-suspension
    integer Npart, Part(MaxForm)                        ! Pollutant form identifier
  

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine allocating memory for the variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_HM_MemAlloc

	integer :: AllocErr=0, aErr

! Allocation of pollutant specific variables
	allocate(DustFlux(4,4,Imin:Imax,Jmin:Jmax,4),&
            &HMinSoil(Imin:Imax,Jmin:Jmax,MaxSurf), stat=aErr)

    AllocErr=AllocErr+aErr

	if(AllocErr/=0) stop 'STOP: Memory allocation error'

end subroutine Atm_HM_MemAlloc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine deallocating memory 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_HM_MemDealloc

	integer :: DeAllocErr=0, dErr

! Deallocation of pollutant specific variables
	deallocate(DustFlux,&
            &HMinSoil, stat=dErr)

    DeAllocErr=DeAllocErr+dErr

	if(DeAllocErr/=0) stop 'STOP: Memory deallocation error'

end subroutine Atm_HM_MemDealloc

#endif

end module Atm_HM_Params