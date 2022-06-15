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
! Module of exchangec parameters
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module Exch_Params

  use GeneralParams

  implicit none

!********************  Deposition storages  ********************!
	real(8), allocatable :: DryDepYear(:,:,:,:,:)		! Annual dry deposition field
	real(8), allocatable :: DryDepMonth(:,:,:,:,:)		! Monthly dry deposition field
	real(8), allocatable :: DryDepDay(:,:,:,:,:)		! Daily dry deposition field
        real(8), allocatable :: DryDepDayTmp(:,:,:,:,:)		! Daily dry deposition field for aggregated LU categories 16-10-2019
	real(8), allocatable :: WetDepYear(:,:,:,:,:)		! Annual wet deposition field
	real(8), allocatable :: WetDepMonth(:,:,:,:,:)		! Monthly wet deposition field
	real(8), allocatable :: WetDepDay(:,:,:,:,:)		! Daily wet deposition field
        real(8), allocatable :: WetDepDayTmp(:,:,:,:,:)		! Daily wet deposition field for aggregated LU categories 16-10-2019
	real(8), allocatable :: FreshDep(:,:,:,:)       	! Annual wet deposition field
	real(8), allocatable :: PrecipYear(:,:)				! Annual precipitation amount
	real(8), allocatable :: PrecipMonth(:,:)			! Monthly precipitation amount
	real(8), allocatable :: PrecipDay(:,:)  			! Daily precipitation amount
	real(8), allocatable :: DryRemYear(:,:,:,:,:)		! Annual dry re-emission field
	real(8), allocatable :: DryRemMonth(:,:,:,:,:)		! Monthly dry re-emission field
	real(8), allocatable :: DryRemDay(:,:,:,:,:)		! Daily dry re-emission field
        real(8), allocatable :: DryRemDayTmp(:,:,:,:,:)		! Daily dry re-emission field for aggregated LU categories 16-10-2019
#if G_POP
    real(8), allocatable :: MassDDSrc(:,:,:)            ! 18-09-2019
    real(8), allocatable :: MassDRSrc(:,:,:)            ! 18-09-2019
    real(8), allocatable :: MassWDSrc(:,:)
#endif
#if RTYPE==2
	real(8), allocatable :: TotDepMatrix(:,:)				! Deposition matrix
	real(8), allocatable :: DryDepMatrix(:,:)				! Deposition matrix
	real(8), allocatable :: WetDepMatrix(:,:)				! Deposition matrix
#endif

!****************  Deposition characteristics  *****************!
	real, allocatable :: Vd(:,:,:,:,:,:)	! Dry deposition velocity (m/sec)
	real Vdcoef(3)				            ! Dry deposition coefficients
	real Ain, Bin, Abelow, Bbelow		    ! Scavenging coefficient
	real Aeff, Asol

!****************  Dry deposition parameters  ******************!
	real Aforest, Bforest, Cforest		! 
	real Agrass, Bgrass, Cgrass, Dgrass
	real Season(Imin:Imax,Jmin:Jmax,12)
	real heightLU(MaxSurf)
	real dispLU(MaxSurf,NumSns)
	real ZoLU(MaxSurf,NumSns)
	real Bm, Gm, Bh, Gh
	real Dp, RhoP, Dfog
	real Rbrok, Rcmin

#if G_POP
    integer, parameter :: KAirMax = 20 ! 30
    integer, parameter :: KMedMax = 24 ! 40
    integer Exch_KMax, Exch_KMed
    real(8) CAtm(KAirMax+KAirMax),CMed(KMedMax)
#if RTYPE==2
    real(8) CABeg(KAirMax), CAFin(KAirMax)
    real(8) Contrib(KAirMax,MaxMatr)
#endif
    real(8)  DeltaZMed(KMedMax), DeltaZ(KAirMax+KAirMax)
    real(8) MV(KMedMax), MA
    real Exch_dT
!    real(8) WetDepVel(KAirMax), WetDepCoeff(KAirMax), WetDepMask(KAirMax)
    real(8) WetDepVel(KAirMax), WetDepCoeff(KAirMax)
    real(8) ExchVel(0:KMedMax,0:KMedMax) , DryDepVel(0:KMedMax), WetDepMask(KMedMax)
    real LC(KMedMax)
    real(8) RHP(KAirMax+KAirMax + KMedMax)                 ! RHP - right-hand parts
	real(8) DCIn(KAirMax - 1), DCOut(KAirMax - 1), DCInMed  !matrix 24.12.2014
    real DeltaT
    integer NSteps

    real(8) :: AirSoilFlux = 0.
    real(8) :: AirVegFlux  = 0.
    real(8) :: AirOcnFlux  = 0.
    real(8) :: SoilAirFlux = 0.
    real(8) :: VegAirFlux  = 0.
    real(8) :: OcnAirFlux  = 0.
    real(8) :: VegSoilFlux = 0.         ! End 19-09-2019
#endif


 contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine allocating memory for the variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Exch_MemAlloc

	integer :: AllocErr=0, aErr

! Allocation of the atmospheric variables
	allocate(DryDepYear(Imin:Imax,bJmin:bJmax,NumForm(Atm),NumSurf,NumSrc),&
			&DryDepMonth(Imin:Imax,bJmin:bJmax,NumForm(Atm),NumSurf,NumSrc),&
			&DryDepDay(Imin:Imax,bJmin:bJmax,NumForm(Atm),NumSurf,NumSrc),&
			&WetDepYear(Imin:Imax,bJmin:bJmax,NumForm(Atm),NumSurf,NumSrc),&
			&WetDepMonth(Imin:Imax,bJmin:bJmax,NumForm(Atm),NumSurf,NumSrc),&
			&WetDepDay(Imin:Imax,bJmin:bJmax,NumForm(Atm),NumSurf,NumSrc),&
			&FreshDep(Imin:Imax,bJmin:bJmax,NumSurf,NumSrc),&
            &PrecipYear(Imin:Imax,bJmin:bJmax),&
            &PrecipMonth(Imin:Imax,bJmin:bJmax),&
            &PrecipDay(Imin:Imax,bJmin:bJmax),&
            &DryRemYear(Imin:Imax,bJmin:bJmax,NumForm(Atm),NumSurf,NumSrc),&
            &DryRemMonth(Imin:Imax,bJmin:bJmax,NumForm(Atm),NumSurf,NumSrc),&
            &DryRemDay(Imin:Imax,bJmin:bJmax,NumForm(Atm),NumSurf,NumSrc),&
#if G_POP
            &DryDepDayTmp(Imin:Imax,bJmin:bJmax,NumForm(Atm),KMedMax,NumSrc),&   ! 16-10-2019
            &DryRemDayTmp(Imin:Imax,bJmin:bJmax,NumForm(Atm),KMedMax,NumSrc),&       ! 16-10-2019
            &WetDepDayTmp(Imin:Imax,bJmin:bJmax,NumForm(Atm),KMedMax,NumSrc),&   ! 16-10-2019
#endif
	    &Vd(Imin:Imax,bJmin:bJmax,NumPer,DryDepNum,NumSurf,2), stat=aErr)
	AllocErr=AllocErr+aErr

! Allocation of the matrix variables
#if RTYPE==2
	allocate(TotDepMatrix(NumSrc,NumRcp), &
        &WetDepMatrix(NumSrc,NumRcp), DryDepMatrix(NumSrc,NumRcp),&
        &stat=aErr)
	AllocErr=AllocErr+aErr
#endif

#if G_POP
    allocate(MassDDSrc(NumForm(Atm),KMedMax,NumSrc),&
            &MassDRSrc(NumForm(Atm),KMedMax,NumSrc),&
            &MassWDSrc(NumForm(Atm),NumSrc),&
            &stat=aErr)                ! 18-09-2019  Corrected 16-10-2019
	AllocErr=AllocErr+aErr

#endif

	if(AllocErr/=0) stop 'STOP: Memory allocation error'

end subroutine Exch_MemAlloc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine deallocating memory 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Exch_MemDealloc

	integer :: DeAllocErr=0, dErr

! Deallocation of the atmospheric variables
	deallocate(DryDepYear,&
	    &DryDepMonth,&
	    &DryDepDay,&
	    &WetDepYear,&
	    &WetDepMonth,&
	    &WetDepDay,&
	    &FreshDep,&
            &PrecipYear,&
            &PrecipMonth,&
            &PrecipDay,&
            &DryRemYear,&
            &DryRemMonth,&
            &DryRemDay,&
#if def
            &WetDepDayTmp,&     ! 16-10-2019
            &DryDepDayTmp,&     ! 16-10-2019
            &DryRemDayTmp,&      ! 16-10-2019
#endif
	    &Vd,&
            &stat=dErr)
	DeAllocErr=DeAllocErr+dErr

! Deallocation of the matrix variables
#if RTYPE==2
	  deallocate(TotDepMatrix,WetDepMatrix,DryDepMatrix, stat=dErr)
	  DeAllocErr=DeAllocErr+dErr
#endif

#if G_POP
	deallocate(MassDDSrc,&
            &MassDRSrc,&
            &MassWDSrc,&
            &stat=dErr)
	DeAllocErr=DeAllocErr+dErr
#endif

	if(DeAllocErr/=0) stop 'STOP1: Memory deallocation error'

end subroutine Exch_MemDealloc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine initializing atmospheric variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Exch_Initial

	DryDepYear=0.
	DryDepMonth=0.			
	DryDepDay=0.
	WetDepYear=0.
	WetDepMonth=0.			
	WetDepDay=0.
	FreshDep=0.
	PrecipYear=0.		
	PrecipMonth=0.		
	PrecipDay=0.		
	DryRemYear=0.
	DryRemMonth=0.
	DryRemDay=0.
#if G_POP
        WetDepDayTmp=0.   ! 16-10-2-19
        DryDepDayTmp=0.   ! 16-10=2019
        DryRemDayTmp=0.   ! 16-10-2019
#endif

end subroutine Exch_Initial

end module Exch_Params
