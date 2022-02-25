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
! Module of vegetation specific parameters
! Version:      1.0
! Modified:     24.11.2014
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifdef M_VEG
module Veg_Params
	
use GeneralParams
implicit none

	integer, parameter :: Veg_NumType = 5		! Number of vegetation types involved
												! 1 - Decid forest,2 - Conif forest, 3 - Grass 
												! 4 - Scrabs,5 - Arable	
	real,parameter :: CEFall=2.19e-9		! Constant of forest litter/Veg exchange, s-1
        real,parameter :: SLeaf = 8.e3				! Ratio of leaf surface to its volume (m2/m3)
	real,parameter :: VegDensity = 267.			! Vegetation density (kg/m3)
        real,parameter :: SLeaf2VegDens = SLeaf/VegDensity
        real,parameter :: LCovVegAccuracy = 0.001       ! Account for veg fractions in cells exceeding tenth of percent
	real,parameter :: ThroughFall = 3.3333E-01	! Part of dry particle deposition to the forest,entered to the vegetation
      
      character*2 SeasDistr					! Switch for season distribution in each cell
	
	real CDVeg				! Degradation coefficient in vegetation
	real CDFall				! Degradation coefficient in forest litter
	real Vmol				! Molar volume (cm**3/mol)

	real(8), allocatable :: Veg_Frac(:,:,:)         ! Veg_Frac(i,j,L) - area fraction in cell (i,j) under vegetation of type L
	real, allocatable :: LAI(:,:)				! LAI(i,j) is the leaf area index in the cell (m2/m2).
	real, allocatable :: LAIPriv(:,:,:)			! Private LAI
	real Alpha(Veg_NumType,12)				! Private LAI coefficients by months or seasons
	real, allocatable :: AlphaCurr(:,:,:)		! Current private LAI coefficients
	real, allocatable :: AlphaConst(:,:,:)		! Current private LAI coefficients

	real(8), allocatable ::  Veg_Conc(:,:,:,:)	! Pollutant "concentration" in vegetation (kg/m2)
	real(8), allocatable ::  Veg_ConcDay(:,:,:,:)
	real(8), allocatable ::  Veg_ConcMonth(:,:,:,:)
	real(8), allocatable ::  Veg_ConcYear(:,:,:,:)
        real(8) Veg_Degr                               ! Mass degraded in vegetation

	real(8), allocatable ::  Fall_Conc(:,:,:,:)		! Pollutant "concentration" in forest litter (kg/m2)
	real(8), allocatable ::  Fall_ConcDay(:,:,:,:)
	real(8), allocatable ::  Fall_ConcMonth(:,:,:,:)
	real(8), allocatable ::  Fall_ConcYear(:,:,:,:)

  real, allocatable :: Veg_VExch(:,:,:,:,:,:)     ! Veg_VExch(i,j,L,m,Dir,p) - exchange velocity in the cell (i,j)
                                                  !    Dir = 1 - from medium m to veg of type L, Dir = 2 - from veg to medium m
! Output parameters

!  logical :: Veg_DailyOutput = .false.            ! If DailyOutput = .true. daily avreaged fields are written
												  
contains
	
!======================================================================================
! Allocate vegetation related arrays
!======================================================================================
subroutine Veg_MemAlloc

  integer Ares

  allocate (Veg_Conc(IMin:IMax,JMin:JMax,Veg_NumType,NumForm(Veg)), &
	& Veg_ConcDay(IMin:IMax,JMin:JMax,Veg_NumType,NumForm(Veg)), &
	& Veg_ConcMonth(IMin:IMax,JMin:JMax,Veg_NumType,NumForm(Veg)), &
	& Veg_ConcYear(IMin:IMax,JMin:JMax,Veg_NumType,NumForm(Veg)), &
        & Fall_Conc(IMin:IMax,JMin:JMax,Veg_NumType,NumSubsMedia(Veg)), &
	& Fall_ConcDay(IMin:IMax,JMin:JMax,Veg_NumType,NumSubsMedia(Veg)), &
	& Fall_ConcMonth(IMin:IMax,JMin:JMax,Veg_NumType,NumSubsMedia(Veg)), &
	& Fall_ConcYear(IMin:IMax,JMin:JMax,Veg_NumType,NumSubsMedia(Veg)), &
	& Veg_Frac(IMin:IMax,JMin:JMax,Veg_NumType), &
	& LAI(IMin:IMax,JMin:JMax), &
	& LAIPriv(IMin:IMax,JMin:JMax,Veg_NumType), &
	& AlphaCurr(IMin:IMax,JMin:JMax,Veg_NumType), &
	& AlphaConst(IMin:IMax,JMin:JMax,Veg_NumType), &
	& Veg_VExch(IMin:IMax,JMin:JMax,Veg_NumType,NumMed,2,NumSubsMedia(Veg)), &
        & stat = ARes)

	if (ARes /= 0) then
	  print*, 'Allocation failed'
	  stop
	end if

        Veg_Conc = 0.
        Veg_ConcDay = 0.
        Veg_ConcMonth = 0.
        Veg_ConcYear = 0.
        Fall_Conc = 0.
        Fall_ConcDay = 0.
        Fall_ConcMonth = 0.
        Fall_ConcYear = 0.
        Veg_Frac = 0.
        LAI = 0.
        LAIPriv = 0.
        AlphaCurr = 0.
        AlphaConst = 0.
        Veg_VExch = 0.

        Veg_Degr = 0.

end subroutine Veg_MemAlloc


!======================================================================================
! Deallocate vegetation related arrays
!======================================================================================
subroutine Veg_MemDealloc

  integer :: DeAllocErr = 0

  deallocate (Veg_Conc, Veg_ConcDay, Veg_ConcMonth, Veg_ConcYear, &
            & Fall_Conc, Fall_ConcDay, Fall_ConcMonth, Fall_ConcYear, &
            & Veg_Frac, LAI, LAIPriv, AlphaCurr, AlphaConst, Veg_VExch,&
            &stat=DeAllocErr)
if(DeAllocErr/=0) stop 'STOP: Memory deallocation error'

end subroutine Veg_MemDealloc


!======================================================================================
! Allocate vegetation related arrays (POPs)
!======================================================================================
subroutine Veg_POP_MemAlloc


end subroutine Veg_POP_MemAlloc


!======================================================================================
! Deallocate vegetation related arrays (POPs)
!======================================================================================
subroutine Veg_POP_MemDealloc


end subroutine Veg_POP_MemDealloc


endmodule Veg_Params
#endif