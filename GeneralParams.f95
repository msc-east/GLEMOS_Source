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
! Module containing general parameters description
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module GeneralParams
    use netcdf
    use typesizes
    implicit none

    character(32) :: MODEL_NAME='GLEMOS' 
    character(64) :: MODEL_LONG_NAME='Global EMEP Multi-Media Modelling System' 
    character(10) :: MODEL_VERSION='2.2.2_dev'
    character(10) :: DumpVersion='2.2.2'

!*************************  Constants  *************************!
    real, parameter :: piNum=3.141592653589793          ! Pi number
    real, parameter :: pi180=1.74532925199433e-2        ! Pi/180
    real, parameter :: pi360=8.726646259971648e-3       ! Pi/360
    real(8) Zero, Eps

!********************  Physical constants  *********************!
    real, parameter :: Ggrav=9.78031                    ! Gravitational acceleration [m/sec2]
    real, parameter :: Rearth=6.371032e6                ! Earth radius [m]
    real, parameter :: Runiv=8.3145                     ! Universal Gas constant [kg*m2/sec2/mol/K]
    real, parameter :: Nav=6.02213e23                   ! Avogadro's number [molec/mol]
    real, parameter :: Rair=287.05                      ! Gas constant for dry air [m2/sec2/K]
    real, parameter :: HumCoef=0.60783                  ! HumCoef=Mdry/Mvap-1
    real, parameter :: Karman=0.4                       ! von Karman's constant
    real, parameter :: Nu=1.39e-5                       ! Air kinematic viscosity [m2/s]
    real, parameter :: RhoWat=1000.                     ! Water density [kg/m3]
    real, parameter :: Lw0=1.e-5                        ! Threshold of cloud LWC [kg/m3]
    real, parameter :: Md=0.028966                      ! Mole weight of dry air [kg/mole]
    real, parameter :: MHgCL2=0.2714                    ! Mole weight of HgCl2 [kg/mole]
    real, parameter :: Dm=4.5e-10                       ! Molecule diameter [m]
    real, parameter :: Cpd=1004.67                      ! Cp of dry air [J/kg/K]
    real, parameter :: CpCoef=0.856                     ! CpCoef=Cpvap/Cpdry-1
    real, parameter :: Prt=0.95                         ! Turbulent Prandtl number
    real, parameter :: kB=1.3807e-23                    ! Boltzmann's constant

!****************  Type of the calculation run  ****************!
    character(10) RunType, GridCode, InitCond
    integer, parameter :: RtypeCompil=RTYPE

!*******************  Grid characteristics  ********************!
    real dXstep                                         ! Grid step in zonal direction (degrees)
    real dYstep                                         ! Grid step in meridional direction (degrees)
    integer, parameter :: Imin=1                        ! Beginning of model domain in I direction
    integer, parameter :: Jmin=1                        ! Beginning of model domain in J direction
    integer, parameter :: Imax=GRIDIMAX                 ! End of model domain in I direction
    integer, parameter :: Jmax=GRIDJMAX                 ! End of model domain in J direction 
    real dLmin                                          ! Minimum grid step in zonal direction (m)
#if (REGTYPE==1)
    integer, parameter :: bImin=Imin                    ! Location of left boundary
    integer, parameter :: bImax=Imax                    ! Location of right boundary
    integer, parameter :: bJmin=Jmin                    ! Location of lower boundary (global)
    integer, parameter :: bJmax=Jmax                    ! Location of upper boundary (global)
    integer, parameter  :: iSP=1, jSP=1                 ! Indexes of South Pole
    integer, parameter  :: iNP=1, jNP=Jmax              ! Indexes of North Pole
#elif (REGTYPE==2)
    integer, parameter :: bImin=Imin-1                  ! Location of left boundary
    integer, parameter :: bImax=Imax+1                  ! Location of right boundary
    integer, parameter :: bJmin=Jmin-1                  ! Location of lower boundary (regional)
    integer, parameter :: bJmax=Jmax+1                  ! Location of upper boundary (regional)
    integer, parameter  :: iSP=bImin, jSP=bJmin         ! Indexes of South Pole
    integer, parameter  :: iNP=bImin, jNP=bJmax         ! Indexes of North Pole
#endif

    real xOrig, yOrig                                   ! Coordinates of the lower left corner of the grid
    real dXmin                                          ! Minimum grid step in zonal direction (radians)
    real dY                                             ! Grid step in meridional direction (radians)
    integer jGlob                                       ! Number of gridcells along the meridian
    real dX(bJmin:bJmax)                                ! Grid step in meridional direction (radians)
    integer iGlob(bJmin:bJmax), minI(bJmin:bJmax), maxI(bJmin:bJmax), bminI(bJmin:bJmax), bmaxI(bJmin:bJmax)
    integer iS(Imin:Imax,bJmin:bJmax,2), iN(Imin:Imax,bJmin:bJmax,2)
    real LongMesh(bImin:bImax,bJmin:bJmax)              ! Longitude of a grid mesh center (aggregated)
    real LongMeshR(bImin:bImax)                         ! Longitude of a grid mesh center (regular)
    real LatMesh(bJmin:bJmax)                           ! Latitude of a grid mesh center
    real MeshArea(bImin:bImax,bJmin:bJmax)              ! Area of a grid mesh (m2)
    real cosdY, cosY(bJmin:bJmax), ctgdY2, tgY(bJmin:bJmax) ! Trig. functions

!******************  Media characteristics  ********************!
    integer NumMed                                      ! Number of media
    integer, parameter :: MaxMed=4                      ! Maximum number of media
    integer, parameter :: Atm=1, Soil=2, Ocn=3, Veg=4   ! Media indexes
    character(6) MedmID(MaxMed)                         ! Names of media
    data MedmID /'Atm','Soil','Ocn','Veg'/              ! Defenition of media names
    integer MedmLst(MaxMed)                             ! List of substance groups in the current model run    
    real Tstep(MaxMed)                                  ! Time step in media
    integer NumStep(MaxMed)                             ! Number of time steps in media
    integer MedOrder(MaxMed)                            ! Order of media according to time step size

!****************  Pollutant characteristics  ******************!
    integer NumSubs                                     ! Number of substances
    integer, parameter :: MaxSubs=50                    ! Maximum number of substances
    character(10) SubsID(MaxSubs)                       ! Pollutant short name
    integer NumGroups                                   ! Number of substance groups
    integer, parameter :: MaxGroups=5                   ! Maximum number of substance groups
    integer, parameter :: HM=1, HG=2, POP=3, TRACER=4, AERO=5   ! Defenition of substance group indentifiers
    character(10) SubsGroupID(MaxGroups)                ! Names of substance group
    data SubsGroupID /'HM','HG','POP','TRACER','AERO'/  ! Defenition of substance group names
    integer SubsGrpLst(MaxGroups)                       ! List of substance groups in the current model run    
    integer GroupInd(MaxGroups)                         ! Index of group in the list of groups 
    integer SubsGroup(MaxSubs)                          ! Correspondens of substances to groups
    integer gSubsNum(MaxGroups)                         ! Number substances in groups
    integer gSubsInd(MaxGroups,MaxSubs)                 ! Substance indexes in whole list of substances
    integer gSubsGroupInd(MaxSubs)                      ! Substance indexes in group list of substances
    integer gSubsMediaInd(MaxMed,MaxSubs)               ! Indexes of particular substances in media 

    integer :: NumForm(MaxMed)=0                        ! Total number of substance forms in media
    integer :: NumSubsMedia(MaxMed)=0                   ! Total number of substances in media    
    integer, parameter :: MaxForm=50                    ! Maximum number of the substance form in media
    integer FormSubs(MaxMed,MaxSubs*MaxForm)            ! Correspondence of forms to substances
    integer sFormNum(MaxMed,MaxSubs)                    ! Number of particular substance forms in media
    integer sFormInd(MaxMed,MaxSubs,MaxForm)            ! Form indexes in the whole list of forms
    character(20) FormID(MaxMed,MaxSubs,MaxForm)        ! Names of the substance forms

    integer :: AtmTransInd(MaxSubs*MaxForm)=0, AtmTransNum=0   ! Atmospheric transport indexes
    integer :: AtmBndInd(MaxSubs*MaxForm)=0, AtmBndNum=0       ! Atmospheric boundary indexes
    integer :: DryDepInd(MaxSubs*MaxForm)=0, DryDepNum=0       ! Dry deposition indexes
    integer :: BCScvInd(MaxSubs*MaxForm)=0, BCScvNum=0         ! Below-cloud deposition indexes
    integer :: ICScvInd(MaxSubs*MaxForm)=0, ICScvNum=0         ! In-cloud deposition indexes
    integer :: AqFormInd(MaxSubs*MaxForm)=0, AqFormNum=0       ! Aqueous form indexes
    integer :: GasExchInd(MaxSubs*MaxForm)=0, GasExchNum=0     ! Gas exchange indexes

    integer :: AntEmisInd(MaxSubs*MaxForm)=0, AntEmisNum=0     ! Anthropogenic emission indexes
    integer :: NatEmisInd(MaxSubs*MaxForm)=0, NatEmisNum=0     ! Natural emission indexes
    integer :: ReEmisInd(MaxSubs*MaxForm)=0, ReEmisNum=0       ! Reemission indexes
    integer :: ReadAntInd(MaxSubs*MaxForm,3)=0                 ! Indexes of reading anthropogenic emission 
    integer :: ReadNatInd(MaxSubs*MaxForm,3)=0                 ! Indexes of reading natural emission 
    integer :: ReadReInd(MaxSubs*MaxForm,3)=0                  ! Indexes of reading re-emission 
    integer :: ReadBndInd(MaxSubs*MaxForm,3)=0                 ! Indexes of reading boundary
#ifdef M_SOIL
    integer :: SoilEmisInd(MaxSubs*MaxForm)=0, SoilEmisNum=0   ! Soil emission indexes 
    integer :: ReadSoilInd(MaxSubs*MaxForm,3)=0                ! Indexes of reading soil emission 
#endif

!****************  Land cover characteristics  *****************!
    integer NumSurf                                     ! Number of surface types
    integer, parameter :: NumSns=5                      ! Number of seasons
    integer, parameter :: MaxSurf=50                    ! Maximal number of surface types
    integer SurfGroup(MaxSurf)                          ! Groups of LC types
    integer, parameter :: Urban=1, Forest=2, Arable=3, Grass=4, Barren=5, Water=6    ! Land cover groups
    integer gNum(6), gInd(6,MaxSurf)                    ! Number of surfaces and surface indexes in a LC group
    real, allocatable ::  LandCover(:,:,:)              ! Distribution of land cover types 
    character(20) SurfType(MaxSurf)

!*****************  Matrix characteristics  ********************!
    integer NumSrc                                      ! Total number of sources
    integer NumAnth                                     ! Number of anthropogenic sources
    integer NumNat                                      ! Number of natural sources
    integer NumBnd                                      ! Number of boundary sources
    integer NumRcp                                      ! Number of receptors
    integer Init                                        ! Initial conditions index
    integer Reem                                        ! Re-emission index
    integer restAnt                                     ! The rest of anthropogenic sources
    integer restNat                                     ! The rest of natural sources
    integer restBnd                                     ! The rest of boundary sources
    integer, parameter :: MaxMatr=100
    integer MaxAnth, MaxNat, MaxRcp, MaxBnd
    character(8) SourcID(MaxMatr), AnthrID(MaxMatr), AnthrIDmax(MaxMatr), NaturID(MaxMatr), NaturIDmax(MaxMatr)
    character(8) RecepID(MaxMatr), RecepIDmax(MaxMatr), BoundID(MaxMatr), BoundIDmax(MaxMatr)
    integer :: NatSrcInd(MaxMatr)
    integer antSRCmode, natSRCmode, rcpSRCmode, bndSRCmode, initSRCmode, reemisSRCmode
#if RTYPE==2
    real, allocatable :: RecepPart(:,:,:)               ! Relative part of receptors in a cell
    real, allocatable :: RecepArea(:)                   !
#endif

!*******************  Time characteristics  ********************!
    integer Year, Month, Day, Period                    ! Components of the current date throughout the model run
    real DayTime                                        ! Time from the beginning of day (sec)
    integer BegDate(3), FinDate(3)                      ! Initial and final dates of the nodel run
    integer, parameter :: da=1, mn=2, yr=3
    integer, parameter :: SecInDay=86400
    integer, parameter :: NumPer=4                      ! Number of meteorological periods in day
    real :: dTinput=real(SecInDay/NumPer)               ! Time step of meteorology input, 6 hours (sec)
    real timeCalc                                       ! Model internal calculation time (sec)
    real timeSum, timeShare(2), timeRun                 ! Real calculation time (sec)
    logical climRun, climReactRun, climLCRun         ! Switch of climatic run
    integer ClimYear, climReactYear, climLCYear      ! Fixed year for climatic run
    integer, parameter :: toDay=1, toMor=2
    real timePer(8)

!********************  Reactant parameters  ********************!
        character(300) ReactSource

!*****************  Meteorological parameters  *****************!
    real, allocatable :: Uwind(:,:,:,:,:)               ! Zonal component of wind velocity [m/sec]
    real, allocatable :: Vwind(:,:,:,:,:)               ! Meridional component of wind velocity [m/sec]
#if A_VTYPE==2
    real, allocatable :: Swind(:,:,:,:,:)               ! Analog of vertical velocity [1/sec]
#endif
    real, allocatable :: Px(:,:,:,:)                    ! P*=Ps-Pt [Pa]
    real, allocatable :: Zs(:,:)                        ! Survace altitude [m]
    real, allocatable :: TempAir(:,:,:,:,:)             ! Air temperature [K]
    real, allocatable :: HumidAir(:,:,:,:,:)            ! Air humidity
    real, allocatable :: PrecRainConv(:,:,:,:)          ! Precipitation intensity [m/sec]
    real, allocatable :: PrecRainStrat(:,:,:,:)         ! Precipitation intensity [m/sec]
    real, allocatable :: CloudConv(:,:,:,:)             ! Convective cloudiness
    real, allocatable :: CloudStrat(:,:,:,:)            ! Large-scale cloudiness
    real, allocatable :: TempSurf(:,:,:)                ! Surface temperature [K]
    real, allocatable :: Roughness(:,:,:)               ! Roughness of the underlying surface [m]
    real, allocatable :: SoilWater(:,:,:)               ! Soil water
    real, allocatable :: SnowDepth(:,:,:)               ! Snow depth
    real, allocatable :: Ksigma(:,:,:,:,:)              ! Eddy diffusion coefficient [m2/sec]
    real, allocatable :: MOLength(:,:,:,:)              ! Monin-Obukhov length [m]
    real, allocatable :: Ufric(:,:,:,:)                 ! Friction velocity [m/sec]
    real, allocatable :: RTemp(:,:,:,:,:)               ! Effective temperature (RT) [m2/sec2]
    real, allocatable :: WaterCont(:,:,:,:)             ! Cloud water content [kg/kg]
    real, allocatable :: LiqCont(:,:,:,:)               ! Liquid water content [kg/m3]
    real, allocatable :: FrozCont(:,:,:,:)              ! Frozen water content [kg/m3]
    real, allocatable :: SeaIce(:,:,:)                  ! Sea ice cover [1]
    real, allocatable :: DensAir(:,:,:)                 ! Air density [kg/m3]
    real, allocatable :: SolarRad(:,:,:)                ! Accumulated solar radiation [J/m2]
     real, allocatable :: BLHeight(:,:,:)               ! Height of boundary layer [1]

!**********************  Date parameters  **********************!
    character(3) MonthName(12)
    data MonthName /'Jan','Feb','Mar','Apr','May','Jun', 'Jul','Aug','Sep','Oct','Nov','Dec'/
    integer MonthDays(12)
    data MonthDays /31,28,31,30,31,30,31,31,30,31,30,31/

!*******************  Output configuration  ********************! 
    logical OutputDumpYearly, OutputDumpMonthly, OutputDumpDaily
    logical OutputBalanceYearly, OutputBalanceMonthly, OutputBalanceDaily, OutputBalance6hourly ! Added AG 24.01.17

#ifdef M_ATM
    logical AtmOutMonitYearly, AtmOutMonitMonthly, AtmOutMonitDaily, AtmOutMonit6hourly, AtmOutMonitHourly
    logical AtmOutMatrixYearly, AtmOutMatrixMonthly, AtmOutMatrixDaily
    logical AtmOutFieldsYearly, AtmOutFieldsMonthly, AtmOutFieldsDaily, AtmOutFields6hourly, AtmOutFieldsHourly
    logical AtmOutNCFYearly, AtmOutNCFMonthly, AtmOutNCFDaily, AtmOutNCF6hourly, AtmOutNCFHourly
#endif
#ifdef M_OCN
    logical OcnOutMonitYearly, OcnOutMonitMonthly, OcnOutMonitDaily, OcnOutMonit6hourly, OcnOutMonitHourly
    logical OcnOutMatrixYearly, OcnOutMatrixMonthly, OcnOutMatrixDaily
    logical OcnOutFieldsYearly, OcnOutFieldsMonthly, OcnOutFieldsDaily, OcnOutFields6hourly, OcnOutFieldsHourly
    logical OcnOutNCFYearly, OcnOutNCFMonthly, OcnOutNCFDaily, OcnOutNCF6hourly, OcnOutNCFHourly
#endif
#ifdef M_SOIL
    logical SoilOutMonitYearly, SoilOutMonitMonthly, SoilOutMonitDaily, SoilOutMonit6hourly, SoilOutMonitHourly
    logical SoilOutMatrixYearly, SoilOutMatrixMonthly, SoilOutMatrixDaily
    logical SoilOutFieldsYearly, SoilOutFieldsMonthly, SoilOutFieldsDaily, SoilOutFields6hourly, SoilOutFieldsHourly
    logical SoilOutNCFYearly, SoilOutNCFMonthly, SoilOutNCFDaily, SoilOutNCF6hourly, SoilOutNCFHourly
#endif
#ifdef M_VEG
    logical VegOutMonitYearly, VegOutMonitMonthly, VegOutMonitDaily ,VegOutMonit6hourly, VegOutMonitHourly
    logical VegOutMatrixYearly, VegOutMatrixMonthly, VegOutMatrixDaily
    logical VegOutFieldsYearly, VegOutFieldsMonthly, VegOutFieldsDaily, VegOutFields6hourly ,VegOutFieldsHourly
    logical VegOutNCFYearly, VegOutNCFMonthly, VegOutNCFDaily, VegOutNCF6hourly, VegOutNCFHourly
#endif

!*****************  Data paths and file names  *****************!
! Input path names
    character(800) MeteoPath, EmisPath, MatrixPath, ReactPath1, ReactPath2, InBoundPath, GeoPath, PropPath, RecepPath,&
        &StatPath, DustPath, SoilPath, ConfigPath, InitPath, LandcoverPath

! Output path names
    character(800)  OutPath
    character(800)  AtmDir, OcnDir, SoilDir, VegDir, DumpDir, BalanceDir
    character(800)  MatrixDir, MonitorDir, FieldsDir, FieldsNcfDir
    character(800)  YearlyDir, MonthlyDir, DailyDir, SixHourlyDir, HourlyDir

! General file names
    character(80) RunName, GridConfName, MediaConfName, LCconfName, PropName, OrogName, LandSource, LandName, & 
            GroupName, LogName, DumpName, SeasonName, ComGeoName, RoughName, OutputConfigName, MatrixConfName

! Atmospheric filenames
    character(80) InitName, DryName, WetName, RemName, NetName, PrecipName, TotName, ConcName, MixRatName,&
        & InPrecName, ConcMonitName, PrecMonitName, FluxMonitName, PrecModelName, StatName,&
        & OzoneName, SO2Name, OHName, AntEmisName, NatEmisName,&
        & RecepArName, MatrixTotDepName, MatrixWetDepName, MatrixDryDepName, MatrixConcName, BoundName,&
        & DryMonitName, WetMonitName, InprecMonitName, NCFMonitName, NCFOutName,&
        & EmisOutAnt, EmisOutNat, EmisOutRem, SoilConcName, DustName

! Soil names
#ifdef M_SOIL
    character(80) SoilPropName, ConcSoilName
#endif
    character(80) LC2SoilTypeName, FocName, AeroSurfName

! Vegetation names
#ifdef M_VEG
    character(80) VegPropName, ConcVegName
#endif

! Ocean names
#ifdef M_OCN
    character(800) OceanPath   !
    character(80) OcnPropName, TopoOcnName, ConcOcnName
#endif
    character(80) LC2OcnName

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine controlling NetCDF operations
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine checkNC(status,i)

    integer, intent (in) :: status, i
    character(80) str_er

    if(status/=nf90_noerr) then
      print *, 'STOP: Netcdf error - ', status, i
      stop 2
    endif

end subroutine checkNC

end module GeneralParams
