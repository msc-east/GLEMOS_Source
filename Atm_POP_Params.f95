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

#ifdef G_POP
Module Atm_POP_Params

  use GeneralParams
  use Atm_Params
  
implicit none

    real, parameter :: RhoOrg = 824.0*1.4e9      ! Octanol density (kg/m3)
    real, parameter :: RhoEC = 1.8e12            ! Soot density (g/cm3)
    real, parameter :: AlphaSoot = 100.          ! Soot specific surface (m2/g)
    integer, parameter :: NReactMax = 10     ! Maximum number of reactents and aerosol species    05.07 2018

! **************** To be read in POP_ReadProps ****************************************************
    real HT0(MaxSubs), HT(MaxSubs)       ! Parameters of Henry's constant temperature dependence
    real DWater(MaxSubs), DAir(MaxSubs)  ! Molecular diffusion coefficients for water and air

! **************** To be read in POP_ReadProps ****************************************************
    real(8) p0(MaxSubs), pT(MaxSubs)              ! Parameters of temperature dependence for pressure over subcooled liquid
    real(8)  Koa0(MaxSubs), KoaT(MaxSubs)    ! Parameters of temperature dependence for octanol-air partition coefficient
    logical TD(MaxSubs)                                 ! Switch for temperature dependence of air degradation
    real ADegr(MaxSubs), Ea(MaxSubs)        ! Parameters of temperature dependence for air degradation
    real CDWinter(MaxSubs)         ! Average air degradation rate in winter
    real CDSpring(MaxSubs)         ! Average air degradation rate in spring/fall
    real CDSummer(MaxSubs)         ! Average air degradation rate in summer
    logical DegOH(MaxSubs)         ! Switch for photodegradation
    logical DegO3(MaxSubs)         ! Switch for photodegradation
    real PhDCoeff(MaxSubs)         ! Coefficient for photodegradation
    real Kd_OC, Kd_EC              ! Degradation coefficients for OC and EC         ! 28.03.2017
    real KO3_OC, KO3_EC              ! Degradation coefficients for OC and EC         ! 03.04.2017
    real(8) WRatioGas0(MaxSubs)       ! Ratio of concentrations in precipitation and air at 0 C
    real(8) WRatioPart
    logical WashExp(MaxSubs)       ! WashExp = .true. - usage of experimental value of washout ratio,
                                           ! WashExp = .false. - calculations via Henry's law coefficient
    real, parameter :: T0 = 283.15
    real, parameter :: HumidCoef=0.60783

    real(8), allocatable :: Atm_WDVelGas(:,:,:,:)  ! Atm_WDVelGas(i,j,k,p) - wet deposition velocity of gas form in cell (i,j,k)
    real(8), allocatable :: Atm_WDCoeff(:,:,:,:)   ! Atm_WDCoeff(i,j,k,p) - wet deposition velocity of particulate form in cell (i,j,k)

    character(8) :: PartType = 'ad'		    ! Switch for the type of partitioning in the atmosphere:
											! ad - adsorption (Yunge-Pankow approach)
											! ab - absorption (Koa approach)
											! co - conbined ad/absorption
    real(8), allocatable :: kpt(:,:,:,:,:)       ! kpt(i,j,k,t,n) - partitioning coefficient between particular and gaseous forms in
                                            ! the atmosphere in the cell i, j, k in the meteoperiod t for POP number n
    real(8), allocatable :: kDegrG(:,:,:,:,:)    ! kDegrG(i,j,k,t,n) - air degradation constant for POP number n in the cell
                                            ! i, j, k in the meteoperiod t (gaseous form)
    real(8), allocatable :: kDegrP(:,:,:,:,:)    ! kDegrP(i,j,k,t,n) - air degradation constant for POP number n in the cell
                                            ! i, j, k in the meteoperiod t (particulate form)
  real(8), allocatable :: DensAirNum(:,:,:)         ! Number air density [molec/cm3]

! Reactant and aerosol quantities related to partitioning and degradation (MOZART data)
  integer NmbAero, NmbReacts                   ! Numbers of aerosol species and reactants used for the pollutant   05.07.2018
  character*10 AName(NReactMax), RName(NReactMax)     ! Names of aerosol species and reactants, used for reading files as prefix
  character*10 AUnit(NReactMax), RUnit(NReactMax)     ! Units for aerosol species and reactants, used in initial files
  real CorrectA(NReactMax), CorrectR(NReactMax)
  real Cp0(NReactMax), CpT(NReactMax)      ! Coefficients for partitioning to given aerosol species
  real Kd0(NReactMax,0:NReactMax), KdT(NReactMax,0:NReactMax)   ! Coefficients for degradation by given reactant on
                                                                ! given aerosol species or air (with second index 0)

  real KReact(NReactMax,0:NReactMax)          ! 10-02-2020
  logical :: WaterCorrect = .false.
  real, allocatable :: Aerofield(:,:,:,:,:,:), Reactfield(:,:,:,:,:,:)   ! 05.07.2018
  real, allocatable :: dAerofield(:,:,:,:,:), dReactfield(:,:,:,:,:)   ! 05.07.2018
  real, allocatable :: d2Aerofield(:,:,:,:,:), d2Reactfield(:,:,:,:,:)   ! 05.07.2018
  real(8), allocatable :: AeroConc(:,:,:,:), ReactConc(:,:,:,:)   ! 05.07.2018

! OH concentration in the air [molec/cm3]
! OC concentration in the aerosol [mkg/m3]
! EC concentration in the aerosol [mkg/m3]
! O3 concentration in the aerosol [mkg/m3]
! Aerosol surface [m2/m3]

! Field of OH mixing ratio [ppb]
! Field of OC mixing ratio [kg/kg]
! Field of ozone mixing ratio [kg/kg]         28.03.2017
! Field of OC mixing ratio [kg/kg]
! Field of aerosol surface [m2/kg]

  ! For ROI-T scheme
  real DegrP0(10,3), DegrPT(10,3), BaseD(11,3), MaxD(11,3), XHalfD(11,3), RateD(11,3), TBound(11,3)


contains

!==============================================================
! Atmosphere POP memory allocation
!==============================================================
subroutine Atm_POP_MemAlloc

    integer AllocErr, aErr

    AllocErr = 0
    allocate(AeroConc(NmbAero,Imin:Imax,Jmin:Jmax,Atm_Kmax),&      ! 05.07.2018
            &ReactConc(NmbReacts,Imin:Imax,Jmin:Jmax,Atm_Kmax),&      ! 05.07.2018
            &DensAirNum(Imin:Imax,Jmin:Jmax,Atm_Kmax),&

            &AeroField(NmbAero,Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2), &      ! 05.07.2018
            &ReactField(NmbReacts,Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer,2), &      ! 05.07.2018
            &dAeroField(NmbAero,Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1), &      ! 05.07.2018
            &dReactField(NmbReacts,Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1), &      ! 05.07.2018
            &d2AeroField(NmbAero,Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1), &      ! 05.07.2018
            &d2ReactField(NmbReacts,Imin:Imax,Jmin:Jmax,Atm_Kmax,NumPer+1), stat = AllocErr)      ! 05.07.2018

    allocate (kpt(IMin:IMax,JMin:JMax,Atm_KMax,NumPer,gSubsNum(GroupInd(POP))), &
            &kDegrG(IMin:IMax,JMin:JMax,Atm_KMax,NumPer,gSubsNum(GroupInd(POP))), &
            &kDegrP(IMin:IMax,JMin:JMax,Atm_KMax,NumPer,gSubsNum(GroupInd(POP))), &
            &Atm_WDVelGas(IMin:IMax,JMin:JMax,Atm_KMax,gSubsNum(GroupInd(POP))),&
            &Atm_WDCoeff(IMin:IMax,JMin:JMax,Atm_KMax,gSubsNum(GroupInd(POP))), stat=aErr)

    AllocErr = AllocErr + aErr

    if (AllocErr /= 0) stop 'STOP: Memory allocation error in Atm_POP_MemAlloc'

end subroutine Atm_POP_MemAlloc

!==============================================================
! Atmosphere POP memory deallocation
!==============================================================
subroutine Atm_POP_MemDealloc      ! modified 05.07.2018


    integer DeallocErr, dErr

    DeallocErr = 0
    deallocate(AeroConc,&
            &ReactConc,&
            &Aerofield,&
            &Reactfield,&
            &DensAirNum,&

            &dAerofield,&
            &dReactfield,&

            &d2Aerofield,&
            &d2Reactfield,&
            stat = DeallocErr)

    deallocate (kpt, &
            &kDegrG, &
            &kDegrP, &
            &Atm_WDVelGas,&
            &Atm_WDCoeff, stat=dErr)

    DeallocErr = DeallocErr + dErr
    if (DeallocErr /= 0) stop 'STOP: Memory deallocation error in Atm_POP_MemDealloc'

end subroutine Atm_POP_MemDealloc

end Module Atm_POP_Params
#endif
