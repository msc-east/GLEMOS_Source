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
  real(8), allocatable :: ConcO3(:,:,:)             ! O3 concentration in the air [molec/cm3]
  real(8), allocatable :: ConcSO2(:,:,:)            ! SO2 concentration in the air [molec/cm3]
  real(8), allocatable :: ConcOH(:,:,:)             ! OH concentration in the air [molec/cm3]
  real(8), allocatable :: ConcBrO(:,:,:)            ! BrO concentration in the air [molec/cm3]
  real(8), allocatable :: ConcBr(:,:,:)             ! BrO concentration in the air [molec/cm3]
  real(8), allocatable :: ConcPM(:,:,:)             ! OH concentration in the air [molec/cm3]
  real(8), allocatable :: ConcNO2(:,:,:)            ! OH concentration in the air [molec/cm3]
  real(8), allocatable :: ConcHO2(:,:,:)            ! OH concentration in the air [molec/cm3]
  real(8), allocatable :: DensAirNum(:,:,:)         ! Number air density [molec/cm3]
  real(8), allocatable :: ConcCl2(:,:,:,:)          ! Cl2 concentration in the air [molec/cm3]
  real(8), allocatable :: ConcHgBr(:,:,:,:)         ! HgBr concentration in the air [kg/m3]
  real(8), allocatable :: PhotoRate(:,:,:,:)        ! Photolysis rates [1/s]
  
  real, allocatable :: O3field(:,:,:,:,:)           ! Field of O3 mixing ratio [ppbv]
  real, allocatable :: SO2field(:,:,:,:,:)          ! Field of SO2 mixing ratio [ppbv]
  real, allocatable :: OHfield(:,:,:,:,:)           ! Field of OH mixing ratio [ppbv]
  real, allocatable :: BrOfield(:,:,:,:,:)          ! Field of BrO mixing ratio [ppbv]
  real, allocatable :: Brfield(:,:,:,:,:)           ! Field of Br mixing ratio [ppbv]
  real, allocatable :: PMfield(:,:,:,:,:)           ! Field of PM mixing ratio [ppbm]
  real, allocatable :: NO2field(:,:,:,:,:)          ! Field of PM mixing ratio [ppbm]
  real, allocatable :: HO2field(:,:,:,:,:)          ! Field of PM mixing ratio [ppbm]
  real, allocatable :: dO3field(:,:,:,:)            ! Time derivative of O3 mixing ratio
  real, allocatable :: dSO2field(:,:,:,:)           ! Time derivative of SO2 mixing ratio
  real, allocatable :: dOHfield(:,:,:,:)            ! Time derivative of OH mixing ratio
  real, allocatable :: dBrfield(:,:,:,:)            ! Time derivative of Br mixing ratio
  real, allocatable :: dBrOfield(:,:,:,:)           ! Time derivative of BrO mixing ratio
  real, allocatable :: dPMfield(:,:,:,:)            ! Time derivative of PM mixing ratio
  real, allocatable :: dNO2field(:,:,:,:)           ! Time derivative of PM mixing ratio
  real, allocatable :: dHO2field(:,:,:,:)           ! Time derivative of PM mixing ratio
  real, allocatable :: d2O3field(:,:,:,:)           ! Second time derivative of O3 mixing ratio
  real, allocatable :: d2SO2field(:,:,:,:)          ! Second time derivative of SO2 mixing ratio
  real, allocatable :: d2OHfield(:,:,:,:)           ! Second time derivative of OH mixing ratio
  real, allocatable :: d2Brfield(:,:,:,:)           ! Second time derivative of Br mixing ratio
  real, allocatable :: d2BrOfield(:,:,:,:)          ! Second time derivative of BrO mixing ratio
  real, allocatable :: d2PMfield(:,:,:,:)           ! Second time derivative of PM mixing ratio
  real, allocatable :: d2NO2field(:,:,:,:)          ! Second time derivative of PM mixing ratio
  real, allocatable :: d2HO2field(:,:,:,:)          ! Second time derivative of PM mixing ratio
  
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
  integer Npart, Noxid, Part(MaxForm), Oxid(MaxForm), GasPart(2,MaxForm)                                  

 
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine allocating memory for the variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_MemAlloc

    integer :: AllocErr=0, aErr

! Allocation of the chemical variables
    allocate(ConcO3(Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &ConcSO2(Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &ConcOH(Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &ConcBr(Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &ConcBrO(Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &ConcPM(Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &ConcNO2(Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &ConcHO2(Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &DensAirNum(Imin:Imax,Jmin:Jmax,Atm_Kmax),&
            &ConcCl2(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer),&
            &ConcHgBr(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumSrc),&
            &O3field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2),&
            &SO2field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2),&
            &OHfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2),&
            &Brfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2),&
            &BrOfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2),&
            &PMfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2),&
            &NO2field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2),&
            &HO2field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2),&
            &dO3field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &dSO2field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &dOHfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &dBrfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &dBrOfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &dPMfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &dNO2field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &dHO2field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &d2O3field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &d2SO2field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &d2OHfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &d2Brfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &d2BrOfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &d2PMfield(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &d2NO2field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
            &d2HO2field(Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1),&
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
    deallocate(ConcO3,&
            &ConcSO2,&
            &ConcOH,&
            &ConcBr,&
            &ConcBrO,&
            &ConcPM,&
            &ConcNO2,&
            &ConcHO2,&
            &DensAirNum,&
            &ConcCl2,&
            &ConcHgBr,&
            &O3field,&
            &SO2field,&
            &OHfield,&
            &Brfield,&
            &BrOfield,&
            &PMfield,&
            &NO2field,&
            &HO2field,&
            &dO3field,&
            &dSO2field,&
            &dOHfield,&
            &dBrfield,&
            &dBrOfield,&
            &dPMfield,&
            &dNO2field,&
            &dHO2field,&
            &d2O3field,&
            &d2SO2field,&
            &d2OHfield,&
            &d2Brfield,&
            &d2BrOfield,&
            &d2PMfield,&
            &d2NO2field,&
            &d2HO2field,&
            &PhotoRate,&
            &stat=dErr)
    DeAllocErr=DeAllocErr+dErr

    if(DeAllocErr/=0) stop 'STOP: Memory deallocation error'

end subroutine Atm_Hg_MemDealloc

#endif

end module Atm_Hg_Params