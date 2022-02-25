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
! Soil compartment specific parameters
! Version:      1.0
! Modified:     24.11.2014
!=======================================================================================================
#ifdef M_SOIL
Module Soil_Params
use GeneralParams
use Soil_POP_Params
implicit none

  integer, parameter :: NLSMax = 10               ! Maximum number of layers in soil (for array declarations)
  integer, parameter :: Soil_TypeMax = 10               ! 

! Common variables

  integer Soil_KMax                               ! Number of soil layers in vertical
  integer Soil_NumType                            ! Number of soil types considered
  real(8) Soil_dz(NLSMax)                            ! Soil_dz(k) - thickness of soil layer k
  real(8) DiffSoil(NLSMax)                             ! Diffusion coefficient in soil
  real(8) LC2SoilType(0:Soil_TypeMax,MaxSurf)           ! Transition matrix from LandCover to types of soil
  real(8) Soil_Disp(0:Soil_TypeMax,NumSns)
  real(8) Soil_Z0(0:Soil_TypeMax,NumSns)
  real Soil_Height(0:Soil_TypeMax)
  real(8), allocatable :: Soil_Conc(:,:,:,:,:)    ! Soil_Conc(i,j,k,L,Form); i, j - horizontal cell numbers
                                                  !                      k - vertical cell number beginning from the surface
											      !                      L - number of soil type
											      !                      Form- number of pollutant form
  real(8), allocatable :: Soil_ConcDay(:,:,:,:,:) ! Soil_ConcDay(i,j,k,L,p); counter for averaging, p - pollutant number
  real(8), allocatable :: Soil_ConcMonth(:,:,:,:,:) ! Soil_ConcMonth(i,j,k,L,p); counter for averaging
  real(8), allocatable :: Soil_ConcYear(:,:,:,:,:)  ! Soil_ConcYear(i,j,k,L,p); counter for averaging
  real(8) Soil_Degr                               ! Mass degraded in soil
  character*10 SoilType(0:Soil_TypeMax)        ! SoilType(L) - name of soil type L (for output)
  real(8), allocatable :: Soil_Frac(:,:,:)           ! Soil_Frac(i,j,L) - area fraction in cell (i,j) under soil of type L
  real(8) Soil_TStep                                 ! Maximum permissible time step
  real(8), allocatable :: Soil_VExch(:,:,:,:,:,:)    ! Soil_VExch(i,j,L,m,Dir,p) - exchange velocity in the cell (i,j)
                                                  !    Dir = 1 - from medium m to soil of type L, Dir = 2 - from soil to medium m
												  !    for pollutant p

! Water flux

  real(8), allocatable :: Ve(:,:)                    ! Velocity of soil solute
  real(8), allocatable :: QPrec(:,:)                 ! Variable for calculating the velocity of soil solute
  real(8), allocatable :: DeltaQPrec(:,:)            ! Variable for calculating the velocity of soil solute
  integer, allocatable :: Inverse(:,:)		      ! How long (in MP) inverse flux can occur
  integer, parameter :: InvMax = 20	              ! 5 days of inverse flux (im meteoperiods)
  real, parameter :: InvConst = 0.6	              ! fraction of precipitation evaporated from soil. Data from
										          ! ���������, ����������, ���������, 1991, ���. 22
										          ! ����� ���������� ������� � ������� = 0.8 �.
										          ! ��������� = 0.485 �.

! Output parameters

    integer Soil_DetalizationLevel                  ! 1 - minimum, 2 - moderate, 3 - full levels of output

    real*8    soilprof(NLSMax)                      ! Coefficients to enter emissions into soil layers
    real*8, allocatable   ::  SoilEmisFlux(:,:,:)
    real*8, allocatable   ::  MassSoilEmis(:)                     ! Mass of soil emission for balanse

    data    soilprof / 0.20690679,0.19681580,0.18265182,0.15740584,0.12284288,0.08681891,0.04655796,&
                        &0., 0., 0. /
! Coefficient defining concentration profile in soil to distribute emission
! according to this profile.
! Needs to be generalised adding calculation of this coeffs depending on the substance!!!

  character*10 :: SoilEmissionStep = 'nodef' !Yearly, Monthly, Daily   

  contains

!======================================================================================
! Allocate soil related arrays
!======================================================================================
subroutine Soil_MemAlloc

  integer Ares, g

  allocate (Soil_Conc(IMin:IMax,JMin:JMax,Soil_KMax,Soil_NumType,NumForm(Soil)), &
	& Soil_ConcDay(IMin:IMax,JMin:JMax,Soil_KMax,Soil_NumType,NumSubs), &
	& Soil_ConcMonth(IMin:IMax,JMin:JMax,Soil_KMax,Soil_NumType,NumSubs), &
	& Soil_ConcYear(IMin:IMax,JMin:JMax,Soil_KMax,Soil_NumType,NumSubs), &
	& Soil_Frac(IMin:IMax,JMin:JMax,0:Soil_NumType), &
	& Soil_VExch(IMin:IMax,JMin:JMax,Soil_NumType,NumMed,2,NumSubs), &
        & Ve(IMin:IMax,JMin:JMax), &
	& QPrec(IMin:IMax,JMin:JMax), &
        & DeltaQPrec(IMin:IMax,JMin:JMax), &
        & Inverse(IMin:IMax,JMin:JMax), &
	& SoilEmisFlux(IMin:IMax,JMin:JMax,NumForm(Soil)), &        ! 09.04.2012 Gusev
        & MassSoilEmis(NumForm(Soil)), &                            ! 09.04.2012 Gusev
        & stat = ARes)

	if (ARes /= 0) then
	  print*, 'Allocation failed'
	  stop
	end if

        Soil_Conc = 0.
        Soil_ConcDay = 0.
        Soil_ConcMonth = 0.
        Soil_ConcYear = 0.
        Soil_Frac = 0.
        Soil_VExch = 0.
        Ve = 0.
        QPrec = 0.
        DeltaQPrec = 0.
        Inverse = 0.
        SoilEmisFlux = 0.
        MassSoilEmis = 0.
        Soil_Degr = 0.

end subroutine Soil_MemAlloc


!======================================================================================
! Deallocate soil related arrays
!======================================================================================
subroutine Soil_MemDealloc

  integer g

  deallocate (Soil_Conc, Soil_ConcDay, Soil_ConcMonth, Soil_ConcYear, &
            & Soil_Frac, Soil_VExch, Ve, QPrec, DeltaQPrec, Inverse, &
            & SoilEmisFlux, MassSoilEmis)                                  ! 09.04.2012 Gusev

end subroutine Soil_MemDealloc


!======================================================================================
! Allocate soil related arrays (POPs)
!======================================================================================
subroutine Soil_POP_MemAlloc

  integer ARes, Nn

Nn = gSubsNum(GroupInd(POP))

allocate (FocEff(IMin:IMax,JMin:JMax), Fdo(IMin:IMax,JMin:JMax), &
	& PartL(IMin:IMax,JMin:JMax,Soil_KMax,NumPer,Nn), &
	& PartG(IMin:IMax,JMin:JMax,Soil_KMax,NumPer,Nn), &
	& PartS(IMin:IMax,JMin:JMax,Soil_KMax,NumPer,Nn),&
        & POPVdAero(IMin:IMax,JMin:JMax,0:Soil_NumType,Nn), &
        & CDSoil(Soil_KMax,Nn), &
        & stat = ARes)

        if (ARes /= 0) then
            print*, 'Allocation failed'
            stop
        end if

        FocEff = 0.
        PartL = 0.
        PartG = 0.
        PartS = 0.
        POPVdAero = 0.
        CDSoil = 0.

end subroutine Soil_POP_MemAlloc


!======================================================================================
! Deallocate soil related arrays (POPs)
!======================================================================================
subroutine Soil_POP_MemDealloc

  deallocate (FocEff, Fdo, PartL, PartG, PartS, POPVdAero, CDSoil)

end subroutine Soil_POP_MemDealloc


end Module Soil_Params
#endif