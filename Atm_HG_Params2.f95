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

module Atm_Hg_Params

#ifdef G_HG

  use GeneralParams
  use Atm_Params

  implicit none

! Quantities related to mercury chemistry
  real(8), allocatable :: ReactConc(:,:,:,:)        ! Reactants concentrations in the air [molec/cm3]
  
  real(8), allocatable :: DensAirNum(:,:,:)         ! Number air density [molec/cm3]
  real(8), allocatable :: ConcCl2(:,:,:,:)          ! Cl2 concentration in the air [molec/cm3]
  real(8), allocatable :: ConcHgBr(:,:,:,:)         ! HgBr concentration in the air [kg/m3]
  real(8), allocatable :: PhotoRate(:,:,:,:)        ! Photolysis rates [1/s]
  
  real(8), allocatable :: ReactField(:,:,:,:,:,:)   ! Field of ractants mixing ratio [ppbv]
  real(8), allocatable :: dReactField(:,:,:,:,:)    ! Time derivative of reactant mixing ratio
  real(8), allocatable :: d2ReactField(:,:,:,:,:)   ! Second time derivative of reeactants mixing ratio
  
  real(8) ConcCl                                    ! [Cl-] concentration in cloud water
  real(8) CHenryHg(2)                               ! Henry constant coefficients for Hg
  real(8) CHenryO3(2)                               ! Henry constant coefficients for ozone
  real(8) CHenryOH(2)                               ! Henry constant coefficients for ozone
  real(8) CHenryHgCl2(2)                            ! Henry constant coefficients for HgCl2
  real(8) CHenryCl2(2)                              ! Henry constant coefficients for Cl2
  real(8) Clcoef(4)                                 ! Chlorine ion reaction constants
  real(8) pH                                        ! Hydrogen ion exponent
  real(8) Rsoot                                     ! Dissiolved-to-adsorbed ratio
  real(8) Kac1, Kac2, Kac3, Kba, Kcb, Kca           ! Reaction constants
  real(8) O3coef(2), Cl2coef(2), OHcoef(2)          ! Gas phase reaction constants
  integer Gas, Dis, Sulf, Chlor                     ! Pollutant form identifiers
  integer Hg0, HgBr2, HgBrOH, HgBrOOH, HgBrONO, HgBr, HOHg, HOHgONO, HOHgOOH,&
          &HgBr2_part, HgBrOH_part, HgBrOOH_part, HgBrONO_part, HgBr_part, HOHg_part,&
          &HOHgONO_part, HOHgOOH_part, HgO_part     ! Pollutant form identifiers
  integer Npart, Noxid, Part(MaxForm), Oxid(MaxForm), GasPart(2,MaxForm), Nreact                                  

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine allocating memory for the variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Atm_Hg_MemAlloc
    integer :: AllocErr=0, aErr

! Allocation of the chemical variables
    allocate(ReactConc(Imin:Imax,Jmin:Jmax,Atm_Kmax,Nreact),&
            &ReactField(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2,Nreact),&
            &dReactField(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1,Nreact),&
            &d2ReactField(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1,Nreact),&
            &DensAirNum(Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &ConcCl2(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer),&
            &ConcHgBr(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumSrc),&
            &PhotoRate(Noxid,Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &stat=aErr)
    AllocErr=AllocErr+aErr

    if(AllocErr/=0) stop 'STOP5: Memory allocation error'

end subroutine Atm_Hg_MemAlloc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine deallocating memory 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_MemDealloc

    integer :: DeAllocErr=0, dErr

! Deallocation of the chemical variables
    deallocate(ReactConc,&
            &DensAirNum,&
            &ConcCl2,&
            &ConcHgBr,&
            &ReactField,&
            &dReactField,&
            &dReactField,&
            &PhotoRate,&
            &stat=dErr)
    DeAllocErr=DeAllocErr+dErr

    if(DeAllocErr/=0) stop 'STOP: Memory deallocation error'

end subroutine Atm_Hg_MemDealloc

#endif

end module Atm_Hg_Params
