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
! Module of atmospheric parameters
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module Atm_Params

  use GeneralParams

  implicit none

!*******************  Grid characteristics  ********************!
    integer, parameter :: Atm_Kmax=A_KMAX	! Number of layers in the model
    real Ptop                   		! Upper boundary of the model domain (Pa)
    real dS(0:Atm_Kmax+1)			! Grid step in vertical direction
    real Sigma(Atm_Kmax)			! Centers of the vertical sigma layers
    real Slevel(0:Atm_Kmax)			! Boundaries the vertical sigma layers
    real, allocatable ::   Veff(:,:,:)		! Grid mesh effective volume
    character(10) VertVelocOrig                 ! Origin of vertical wind velocity

!*****************  Surface pressure variables  ****************!
    real(8), allocatable ::  PxCurr(:,:)	! P* at the current time step [Pa]
    real(8), allocatable ::  PxNext(:,:)	! P* at the next time step [Pa]
    real(8), allocatable ::  PxDay(:,:)		! Daily mean P* [Pa] (AG 09-01-2017)
    real(8), allocatable ::  PxMonth(:,:)	! Monthly mean P* [Pa] (AG 09-01-2017)
    real(8), allocatable ::  PxYear(:,:)	! Yearly mean P* [Pa] (AG 09-01-2017)
    real(8), allocatable ::  Pcalc(:,:,:)

!*******************  Current met variables  *******************!
    real, allocatable :: Ucurr(:,:,:)           ! Zonal component of wind speed at the current time step [m/s]
    real, allocatable :: Vcurr(:,:,:)           ! Meridional component of wind speed at the current time step [m/s]
    real, allocatable :: TairCurr(:,:,:)        ! Air temperature at the current time step [K]
    real, allocatable :: RTcurr(:,:,:)
    real, allocatable :: dUwind(:,:,:,:)	! Zonal component of wind velocity [m/sec]
    real, allocatable :: dVwind(:,:,:,:)	! Meridional component of wind velocity [m/sec]
    real, allocatable :: dPxdT(:,:,:)		! Time derivative of the surface pressure [Pa/sec]
    real, allocatable :: d2PxdT(:,:,:)		! Second time derivative of surface pressure
    real, allocatable :: d2Uwind(:,:,:,:)       ! Second time derivative of
    real, allocatable :: d2Vwind(:,:,:,:)       ! Second time derivative of
    real, allocatable :: dTair(:,:,:,:)         ! Time derivative of air temperature
    real, allocatable :: d2Tair(:,:,:,:)        ! Second time derivative of air temperature
    real, allocatable :: dKsigma(:,:,:,:)       ! Time derivative of eddy diffusion coefficient
    real, allocatable :: d2Ksigma(:,:,:,:)      ! Second time derivative of eddy diffusion coefficient
#if A_VTYPE==2
    real, allocatable :: Scurr(:,:,:)           ! Vertical component of wind speed at the current time step [m/s]
    real, allocatable :: dSwind(:,:,:,:)	! First time derivative of vertical wind speed
    real, allocatable :: d2Swind(:,:,:,:)       ! Second time derivative of vertical wind speed
#endif

!*********  Pollutant concentration in the atmosphere  *********!
    real(8), allocatable :: Atm_MixRatio(:,:,:,:)	! Pollutant mixing ratio [kg/kg]
    real(8), allocatable :: Atm_Conc(:,:,:,:)		! Pollutant air concentration [kg/m3]
#if RTYPE==2
    real(8), allocatable :: Atm_Contrib(:,:,:,:,:)	! Contribution of continents [kg/m3]
#endif
    real(8), allocatable :: Atm_BoundW(:,:,:,:,:)	! Western boundary mixing ratio [kg/kg]
    real(8), allocatable :: Atm_BoundE(:,:,:,:,:)	! Eastern boundary mixing ratio [kg/kg]
    real(8), allocatable :: Atm_BoundS(:,:,:,:,:)	! Southern boundary mixing ratio [kg/kg]
    real(8), allocatable :: Atm_BoundN(:,:,:,:,:)	! Nothern boundary mixing ratio [kg/kg]

!******************  Concentration averages  *******************!
    real(8), allocatable :: Atm_ConcYear(:,:,:,:,:)	! Yearly mean air concentration
    real(8), allocatable :: Atm_ConcMonth(:,:,:,:,:)	! Monthly mean air concentration
    real(8), allocatable :: Atm_ConcDay(:,:,:,:,:)	! Dayly mean air concentration
    real(8), allocatable :: Atm_MixYear(:,:,:,:,:)	! Yearly mean mixing ratio
    real(8), allocatable :: Atm_MixMonth(:,:,:,:,:)	! Monthly mean mixing ratio
    real(8), allocatable :: Atm_MixDay(:,:,:,:,:)	! Dayly mean mixing ratio
#if RTYPE==2
    real(8), allocatable :: ConcMatrix(:,:)		! Concentration matrix
#endif

!******************  Mass balance storages  ********************!
    real(8), allocatable :: MassAtmInit(:)		! Initial pollutant mass in the atmosphere
    real(8), allocatable :: MassAtmAntEmis(:)		! Anthropogenic emission mass
    real(8), allocatable :: MassAtmNatEmis(:)		! Natural emission mass
    real(8), allocatable :: MassReEmis(:)		! Re-emission mass
    real(8), allocatable :: MassAtmFin(:)		! Final pollutant mass in the atmosphere
    real(8), allocatable :: MassAtmUp(:)		! Mass incoming through the upper boundary
    real(8), allocatable :: MassAtmBnd(:)		! Mass incoming through the lateral boundary
    real(8), allocatable :: MassDryDep(:)		! Mass dry deposited
    real(8), allocatable :: MassDryRem(:)		! Mass dry deposited
    real(8), allocatable :: MassWetDep(:)		! Mass wet deposited
    real(8), allocatable :: MassChemEx(:)		! Mass of chemical exchange
    real(8), allocatable :: MassAtmDegr(:)		! Mass degradated in the atmosphare
    real(8), allocatable :: MassAtmPartit(:)		! Gas-Particle partitioning

#if RTYPE==2
    real(8), allocatable :: MatrInit(:)			! Mass of chemical exchange
    real(8), allocatable :: MatrEmis(:)			! Mass of chemical exchange
    real(8), allocatable :: MatrBndIn(:)		! Mass of chemical exchange
    real(8), allocatable :: MatrBndOut(:)		! Mass of chemical exchange
    real(8), allocatable :: MatrDep(:)			! Mass of chemical exchange
    real(8), allocatable :: MatrAtm(:)			! Mass of chemical exchange
    real(8), allocatable :: MatrBalance(:)		! Mass of chemical exchange
#endif

!*******************  Emission parameters  *********************!
    integer HemisAnt, HemisNat                      ! Number of emission layers
    real(8), allocatable :: AntEmisFlux(:,:,:,:,:)	! Emission flux from anthropogenic sources [kg/sec]
    real(8), allocatable :: NatEmisFlux(:,:,:,:)	! Emission flux from natural sources [kg/sec]
    real(8), allocatable :: ReEmisFlux(:,:,:,:)		! Remission flux [kg/sec]
    real(8), allocatable :: AntEmisMonth(:,:)		! Accumulated monthly emission from anthropogenic sources [kg/month]
    real(8), allocatable :: AntEmisYear(:,:)		! Accumulated annual emission from anthropogenic sources [kg/yr]
    real(8), allocatable :: NatEmisMonth(:,:)		! Accumulated monthly emission from natural sources [kg/month]
    real(8), allocatable :: NatEmisYear(:,:)		! Accumulated annual emission from natural sources [kg/yr]
    real(8), allocatable :: ReEmisMonth(:,:)		! Accumulated monthly emission from natural sources [kg/month]
    real(8), allocatable :: ReEmisYear(:,:)		! Accumulated annual emission from natural sources [kg/yr]

    type station
    character(300) nameSt
    character(300) indSt
    real longSt, latSt, Ast(2,2)
    integer iSt(2,2), jSt(2,2)
    real stMod(2,2,MaxForm)
    endtype station

    integer, parameter :: MaxL=200
    integer, parameter :: NstatMax=10000
    integer, parameter :: MaxB=800

    integer Nair, Nprec
    type(station) AirStat(NstatMax), PrecStat(NstatMax), ReadStat

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine allocating memory for the variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_MemAlloc

    integer :: AllocErr=0, aErr

    allocate(Veff(bImin:bImax,bJmin:bJmax,Atm_Kmax+1), stat=aErr)
    AllocErr=AllocErr+aErr

    allocate(Uwind(Imin-1:bImax,bJmin:bJmax,Atm_Kmax,NumPer,2),&
            &Vwind(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer,2),&
            &dUwind(Imin-1:bImax,bJmin:bJmax,Atm_Kmax,NumPer+1),&
            &dVwind(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer+1),&
#if A_VTYPE==2
            &Swind(bImin:bImax,bJmin:bJmax,0:Atm_Kmax,NumPer,2),&
            &dSwind(bImin:bImax,bJmin:bJmax,0:Atm_Kmax,NumPer+1),&
            &d2Swind(bImin:bImax,bJmin:bJmax,0:Atm_Kmax,NumPer+1),&
            &Scurr(bImin:bImax,bJmin:bJmax,0:Atm_Kmax),&
#endif
            &Px(bImin:bImax,bJmin:bJmax,NumPer,2),&
            &dPxdT(bImin:bImax,bJmin:bJmax,NumPer+1),&
            &d2PxdT(bImin:bImax,bJmin:bJmax,NumPer+1),&
            &Zs(bImin:bImax,bJmin:bJmax),&
            &TempAir(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer,2),&
            &HumidAir(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer,2),&
            &PrecRainConv(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer),&
            &PrecRainStrat(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer),&
            &TempSurf(bImin:bImax,bJmin:bJmax,NumPer),&
            &Roughness(bImin:bImax,bJmin:bJmax,NumPer),&
            &SoilWater(bImin:bImax,bJmin:bJmax,NumPer),&
            &SnowDepth(bImin:bImax,bJmin:bJmax,NumPer),&
            &Ksigma(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer,2),&
            &MOLength(bImin:bImax,bJmin:bJmax,NumPer,2),&
            &Ufric(bImin:bImax,bJmin:bJmax,NumPer,2),&
            &RTemp(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer,2),&
            &WaterCont(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer),&
            &LiqCont(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer),&
            &FrozCont(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer),&
            &SeaIce(bImin:bImax,bJmin:bJmax,NumPer),&
            &DensAir(bImin:bImax,bJmin:bJmax,Atm_Kmax),&
            &SolarRad(bImin:bImax,bJmin:bJmax,NumPer),&
            &BLHeight(bImin:bImax,bJmin:bJmax,NumPer),&
            &d2Uwind(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer+1),&
            &d2Vwind(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer+1),&
            &dTair(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer+1),&
            &d2Tair(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer+1),&
            &dKsigma(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer+1),&
            &d2Ksigma(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumPer+1),&
            &TairCurr(bImin:bImax,bJmin:bJmax,Atm_Kmax),&
            &RTcurr(bImin:bImax,bJmin:bJmax,Atm_Kmax),&
            &Ucurr(bImin:bImax,bJmin:bJmax,Atm_Kmax),&
            &Vcurr(bImin:bImax,bJmin:bJmax,Atm_Kmax),&
            &PxCurr(bImin:bImax,bJmin:bJmax),&
            &PxNext(bImin:bImax,bJmin:bJmax),&
            &PxDay(bImin:bImax,bJmin:bJmax),&                                       ! (AG 09-01-2017)
            &PxMonth(bImin:bImax,bJmin:bJmax),&                                     ! (AG 09-01-2017)
            &PxYear(bImin:bImax,bJmin:bJmax),&                                      ! (AG 09-01-2017)
            &Pcalc(bImin:bImax,bJmin:bJmax,Atm_Kmax), stat=aErr)
    AllocErr=AllocErr+aErr

    Uwind=0.
    Vwind=0.
    dUwind=0.
    dVwind=0.
#if A_VTYPE==2
    Swind=0.
    dSwind=0.
    d2Swind=0.
    Scurr=0.
#endif
    Px=0.
    dPxdT=0.
    Zs=0.
    TempAir=0.
    HumidAir=0.
    PrecRainConv=0.
    PrecRainStrat=0.
    TempSurf=0.
    Roughness=0.
    SoilWater=0.
    SnowDepth=0.
    Ksigma=0.
    MOLength=0.
    Ufric=0.
    RTemp=0.
    WaterCont=0.
    LiqCont=0.
    FrozCont=0.
    SeaIce=0.
    DensAir=0.
    SolarRad=0.
    BLHeight=0.
    d2PxdT=0.
    d2Uwind=0.
    d2Vwind=0.
    dTair=0.
    d2Tair=0.
    dKsigma=0.
    d2Ksigma=0.
    TairCurr=0.
    RTcurr=0.
    Ucurr=0.
    Vcurr=0.
    PxCurr=0.
    PxNext=0.
    PxDay = 0.                                                                        ! (AG 09-01-2017)
    PxMonth = 0.                                                                      ! (AG 09-01-2017)
    PxYear = 0.                                                                       ! (AG 09-01-2017)
    Pcalc=0.

! Allocation of the atmospheric variables
    allocate(AntEmisFlux(Imin:Imax,bJmin:bJmax,AntEmisNum,NumAnth,Atm_Kmax),&
            &NatEmisFlux(Imin:Imax,bJmin:bJmax,NatEmisNum,NumNat),&
            &ReEmisFlux(Imin:Imax,bJmin:bJmax,ReEmisNum,NumSrc),&
            &Atm_MixRatio(bImin:bImax,bJmin:bJmax,0:Atm_Kmax+1,NumForm(Atm)),&
            &Atm_Conc(bImin:bImax,bJmin:bJmax,0:Atm_Kmax+1,NumForm(Atm)),&
            &Atm_ConcYear(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumForm(Atm),NumSrc),&
            &Atm_MixYear(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumForm(Atm),NumSrc),&
            &Atm_ConcMonth(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumForm(Atm),NumSrc),&
            &Atm_MixMonth(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumForm(Atm),NumSrc),&
            &Atm_ConcDay(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumForm(Atm),NumSrc),&
            &Atm_MixDay(bImin:bImax,bJmin:bJmax,Atm_Kmax,NumForm(Atm),NumSrc),&
            &LandCover(Imin:Imax,bJmin:bJmax,NumSurf),&          
            &Atm_BoundW(bJmin:bJmax,Atm_Kmax,NumForm(Atm),NumPer,NumSrc),&
            &Atm_BoundE(bJmin:bJmax,Atm_Kmax,NumForm(Atm),NumPer,NumSrc),&
            &Atm_BoundS(bImin:bImax,Atm_Kmax,NumForm(Atm),NumPer,NumSrc),&
            &Atm_BoundN(bImin:bImax,Atm_Kmax,NumForm(Atm),NumPer,NumSrc),&
            &AntEmisMonth(Imin:Imax,Jmin:Jmax),&
            &AntEmisYear(Imin:Imax,Jmin:Jmax),&
            &NatEmisMonth(Imin:Imax,Jmin:Jmax),&
            &NatEmisYear(Imin:Imax,Jmin:Jmax),&
            &ReEmisMonth(Imin:Imax,Jmin:Jmax),&
            &ReEmisYear(Imin:Imax,Jmin:Jmax),&
            &MassAtmInit(NumForm(Atm)),&
            &MassAtmAntEmis(NumForm(Atm)),&
            &MassAtmNatEmis(NumForm(Atm)),&
            &MassReEmis(NumForm(Atm)),&
            &MassAtmFin(NumForm(Atm)),&
            &MassAtmUp(NumForm(Atm)),&
            &MassAtmBnd(NumForm(Atm)),&
            &MassDryDep(NumForm(Atm)),&
            &MassDryRem(NumForm(Atm)),&
            &MassWetDep(NumForm(Atm)),& 
            &MassChemEx(NumForm(Atm)),&
            &MassAtmDegr(NumForm(Atm)),&                
            &MassAtmPartit(NumForm(Atm)), stat=aErr)    
    AllocErr=AllocErr+aErr

! Allocation of the matrix variables
#if RTYPE==2
    allocate(ConcMatrix(NumSrc,NumRcp),&
            &Atm_Contrib(bImin:bImax,bJmin:bJmax,0:Atm_Kmax+1,NumForm(Atm),NumSrc),&
            &RecepPart(Imin:Imax,bJmin:bJmax,NumRcp),&
            &RecepArea(NumRcp), stat=aErr)
    AllocErr=AllocErr+aErr
    allocate(MatrInit(NumSrc),&
            &MatrEmis(NumSrc),&
            &MatrBndIn(NumSrc),&
            &MatrBndOut(NumSrc),&
            &MatrDep(NumSrc),&
            &MatrAtm(NumSrc),&
            &MatrBalance(NumSrc),&
            &stat=aErr)
    AllocErr=AllocErr+aErr
#endif

    if(AllocErr/=0) then
        print *, 'STOP2: Memory allocation error'
        stop
    endif

end subroutine Atm_MemAlloc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine deallocating memory 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_MemDealloc

    integer :: DeAllocErr=0, dErr

    deallocate(Veff, stat=dErr)
    DeAllocErr=DeAllocErr+dErr

! Deallocation of meteorological parameters
    deallocate( Uwind,&
                &Vwind,&
                &dUwind,&
                &dVwind,&
#if A_VTYPE==2
                &Swind,&
                &dSwind,&
                &d2Swind,&
                &Scurr,&
#endif
                &Px,&
                &dPxdT,&
                &Zs,&
                &TempAir,&
                &HumidAir,&
                &PrecRainConv,&
                &PrecRainStrat,&
                &TempSurf,&
                &Roughness,&
                &SoilWater,&
                &SnowDepth,&
                &Ksigma,&
                &MOLength,&
                &Ufric,&
                &RTemp,&
                &WaterCont,&
                &LiqCont,&
                &FrozCont,&
                &SeaIce,&
                &DensAir,&
                &SolarRad,&
                &BLHeight,&
                &d2PxdT,&
                &d2Uwind,&
                &d2Vwind,&
                &dTair,&
                &d2Tair,&
                &dKsigma,&
                &d2Ksigma,&
                &TairCurr,&
                &RTcurr,&
                &Ucurr,&
                &Vcurr,&
                &PxCurr,&
                &PxNext,&
                &PxDay,&
                &PxMonth,&                       
                &PxYear,&
                &Pcalc, stat=dErr)
    DeAllocErr=DeAllocErr+dErr

! Deallocation of the atmospheric variables
    deallocate(AntEmisFlux,&
                &NatEmisFlux,&
                &ReEmisFlux,&
                &Atm_MixRatio,&
                &Atm_Conc,&
                &Atm_ConcYear,&
                &Atm_MixYear,&
                &Atm_ConcMonth,&
                &Atm_MixMonth,&
                &Atm_ConcDay,&
                &Atm_MixDay,&
                &LandCover,&
                &Atm_BoundW,&
                &Atm_BoundE,&
                &Atm_BoundS,&
                &Atm_BoundN,&
                &AntEmisMonth,&
                &AntEmisYear,&
                &NatEmisMonth,&
                &NatEmisYear,&
                &ReEmisMonth,&
                &ReEmisYear,&
                &MassAtmInit,&
                &MassAtmAntEmis,&
                &MassAtmNatEmis,&
                &MassReEmis,&
                &MassAtmFin,&
                &MassAtmUp,&
                &MassAtmBnd,&
                &MassDryDep,&
                &MassDryRem,&
                &MassWetDep,&
                &MassChemEx,&
                &MassAtmDegr,&
                &MassAtmPartit, stat=dErr)
    DeAllocErr=DeAllocErr+dErr

! Deallocation of the matrix variables
#if RTYPE==2
    deallocate(ConcMatrix,&
                &Atm_Contrib,&
                &RecepPart,&
                &RecepArea, stat=dErr)
    DeAllocErr=DeAllocErr+dErr
    
    deallocate(MatrInit,&
                &MatrEmis,&
                &MatrBndIn,&
                &MatrBndOut,&
                &MatrDep,&
                &MatrAtm,&
                &MatrBalance,&
                &stat=dErr)
    DeAllocErr=DeAllocErr+dErr
#endif

    if(DeAllocErr/=0) then
        print *, 'STOP2: Memory deallocation error'
        stop
    endif

end subroutine Atm_MemDealloc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine initializing atmospheric variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Initial

    Atm_MixRatio=0.
#if RTYPE==2
    Atm_Contrib=0.
    MatrInit=0.
    MatrEmis=0.
    MatrBndIn=0
    MatrBndOut=0.
    MatrDep=0.
    MatrAtm=0.
    MatrBalance=0.    
#endif
    Atm_ConcYear=0.
    Atm_ConcMonth=0.
    Atm_ConcDay=0.
    Atm_MixYear=0.
    Atm_MixMonth=0.
    Atm_MixDay=0.
    AntEmisFlux=0.
    NatEmisFlux=0.
    ReEmisFlux=0.
    AntEmisYear=0.
    AntEmisMonth=0.
    NatEmisYear=0.
    NatEmisMonth=0.
    ReEmisYear=0.
    ReEmisMonth=0.
    Atm_BoundW=0.
    Atm_BoundE=0.
    Atm_BoundS=0.
    Atm_BoundN=0.

    MassAtmInit=0.
    MassAtmAntEmis=0.
    MassAtmNatEmis=0.
    MassReEmis=0.
    MassAtmFin=0.
    MassAtmUp=0.
    MassAtmBnd=0.
    MassDryDep=0.
    MassWetDep=0.
    MassDryRem=0.
    MassAtmDegr=0.
    MassAtmPartit=0.
    MassChemEx=0.

end subroutine Atm_Initial


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function calculating a grid mesh volume
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MeshVolume(i,j,k)

    integer i, j, k
    real MeshVolume

    MeshVolume=Veff(i,j,k)*PxCurr(i,j)/DensAir(i,j,k) 

end function MeshVolume


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function calculating the altitude of a sigma level 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Altitude(i,j,k)

    integer i, j, k, l
    real Altitude, dZ

    Altitude=Zs(i,j)
    do l=1, k-1
      dZ=RTcurr(i,j,l)/Ggrav*log((Slevel(l-1)*PxCurr(i,j)+Ptop)/(Slevel(l)*PxCurr(i,j)+Ptop))
      Altitude=Altitude+dZ
    enddo    
    Altitude=Altitude+RTcurr(i,j,k)/Ggrav*log((Slevel(k-1)*PxCurr(i,j)+Ptop)/(Sigma(k)*PxCurr(i,j)+Ptop))

end function Altitude

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function calculating height of the upper boundary of a sigma layer
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sheight(i,j,k)

    integer i, j, k, l
    real Sheight, dZ

    Sheight=Zs(i,j)
    do l=1, k
      dZ=RTcurr(i,j,l)/Ggrav*log((Slevel(l-1)*PxCurr(i,j)+Ptop)/(Slevel(l)*PxCurr(i,j)+Ptop))
      Sheight=Sheight+dZ
    enddo    

end function Sheight

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating RT variable
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RTcurrent

    integer i, j, k
    real Tair, Mv

    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        do i=minI(j), maxI(j)
          Tair=TairCurr(i,j,k)
          Mv=HumidAir(i,j,k,Period,toDay)
          RTcurr(i,j,k)=(1.+HumCoef*Mv/(1.+Mv))*Rair*Tair
        enddo
      enddo
    enddo

end subroutine RTcurrent


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating air density 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine AirDensity

    integer j, k, Form

    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        DensAir(minI(j):maxI(j),j,k)=(Sigma(k)*PxCurr(minI(j):maxI(j),j)+Ptop)&
               &/RTcurr(minI(j):maxI(j),j,k)
      enddo
    enddo

end subroutine AirDensity


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine restoring pollutant air concentration from mixing ratio 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MixRatioToConc

    integer j, k, Form

    do Form=1, NumForm(Atm)
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          Atm_Conc(minI(j):maxI(j),j,k,Form)=Atm_MixRatio(minI(j):maxI(j),j,k,Form)*&
            &DensAir(minI(j):maxI(j),j,k)
        enddo
      enddo
    enddo

end subroutine MixRatioToConc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine restoring pollutant mixing ratio from air concentration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ConcToMixRatio

    integer j, k, Form

    do Form=1, NumForm(Atm)
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          Atm_MixRatio(minI(j):maxI(j),j,k,Form)=Atm_Conc(minI(j):maxI(j),j,k,Form)/&
            &DensAir(minI(j):maxI(j),j,k)
        enddo
      enddo
    enddo

end subroutine ConcToMixRatio

end module Atm_Params
