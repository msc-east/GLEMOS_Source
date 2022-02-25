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
! Module containing general model procedures
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module GeneralProcs

  use GeneralParams
  use Atm_Input
#ifdef M_ATM
  use Atm_General
#endif
#ifdef M_OCN
  use Ocn_General
#endif
#ifdef M_SOIL
  use Soil_General
#endif
#ifdef M_VEG
  use Veg_General
#endif
#ifdef G_POP
  use Exch_General_POP 
#else
  use Exch_General
#endif

  implicit none

  character(800), private :: fileName, fullName
  character(4),  private :: YearNum
  character(2),  private :: DayNum, MonthNum, YearShort

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initialization of model variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Initial

#ifdef DEBUG_MODE
    print *, '>- Entering Initial...'
#endif
! Initializing mashine time variable
    call etime(timeShare,timeSum)

! Reading model run information
    call RunInfo

! Reading model run configuration
    call ReadConfig
    call LogFile('Reading configuration files',0)
      
! Definding geometrical parameters
    call GridGeometry
    call LogFile('Defining grid geometry',0)

#ifdef G_AERO
    call init_aqueous()
#endif

#ifdef DEBUG_MODE
    print *, '<- Exit Initial'
#endif

end subroutine Initial


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine managing memory allocation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Memory(operation)

    character(*) operation

#ifdef M_ATM
    call Atm_Memory(operation)
#endif
#ifdef M_OCN
    call Ocn_Memory(operation)
#endif
#ifdef M_SOIL
    call Soil_Memory(operation)
#endif
#ifdef M_VEG
    call Veg_Memory(operation)
#endif
    call Exch_Memory(operation)

    selectcase(operation)
!---------------------------------------------------------------------------------
    case('Initial')
      call LogFile('Allocating memory for variables',0)
!---------------------------------------------------------------------------------
    case('Final')
      call LogFile('Memory deallocation  ',0)
    endselect

end subroutine Memory


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading configuration files 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadConfig

! Reading model grid configuration
    call GridConfig

! Reading land cover configuration
    call LC_config

! Reading model run output configuration
    call ReadOutputConfig                       !!!!!! New Output
    
! Creating a logfile
    call LogFile('Initial')

! Reading media configurations
#ifdef M_ATM
    call Atm_ReadConfig
#endif

#ifdef M_OCN
    call Ocn_ReadConfig
#endif

#ifdef M_SOIL
    call Soil_ReadConfig
#endif

#ifdef M_VEG
    call Veg_ReadConfig
#endif
    
end subroutine ReadConfig


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine managing model input 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine InputData(Prd)

    character(*) Prd

#ifdef DEBUG_MODE
    print *, '>- Entering InputData ', Prd
#endif

!    if(Prd=='Initial') then
    selectcase(Prd)
!---------------------------------------------------------------------------------
    case('Initial')

! Reading geophysical information 
      call GeoPhysData
      call ReadLandCover
      call LogFile('Reading geophysical information',0)

! Initializing model variables
#ifdef M_ATM
      call Atm_Initial
#endif

! Reading initial conditions
      call LogFile('Reading initial conditions',0)
      selectcase(InitCond)
        case('zero')
          call Exch_Initial
        case('cond')
          call ReadInitCond
          call Exch_Initial
        case('dump')
          call ReadDump
      endselect
!---------------------------------------------------------------------------------
    case('Yearly')

      call ReadLandCover
!---------------------------------------------------------------------------------
    endselect
!    endif
  
! Media specific input
#ifdef M_ATM
    call Atm_InputData(Prd)
#endif
    call Exch_InputData(Prd)
#ifdef M_OCN
    call Ocn_InputData(Prd)
#endif
#ifdef M_SOIL
    call Soil_InputData(Prd)
#endif
#ifdef M_VEG
    call Veg_InputData(Prd)
#endif

#ifdef DEBUG_MODE
    print *, '<- Exit InputData ', Prd
#endif

end subroutine InputData


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine managing model output 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine OutputData(Prd,Medm)

    character(*) Prd
    integer, optional :: Medm

    if(present(Medm)) then
      selectcase(Medm)
#ifdef M_ATM
        case(Atm)
          call Atm_OutputData(Prd)
#endif
#ifdef M_OCN
        case(Ocn)
          call Ocn_OutputData(Prd)
#endif
#ifdef M_SOIL
        case(Soil)
          call Soil_OutputData(Prd)
#endif
#ifdef M_VEG
          case(Veg)
          call Veg_OutputData(Prd)
#endif
        case default
          print '(/,"STOP: Inappropriate medium number ''",i2,"''",/)', Medm
      endselect
    else

#ifdef M_ATM    
      call Atm_OutputData(Prd)
#endif
#ifdef M_OCN
      call Ocn_OutputData(Prd)
#endif
#ifdef M_SOIL
      call Soil_OutputData(Prd)
#endif
#ifdef M_VEG
      call Veg_OutputData(Prd)
#endif
    endif

! Calculate mass balance characteristics and store to files
    call OutputBalance(Prd)

! General output
    select case(Prd)
        case('Final')
            call LogFile('*** Finish of the calculation cycle ***')
            call etime(timeShare,timeSum)
            if(timeSum>0.) then
                timeRun=timeSum
            else
                timeRun=timeShare(1)+timeShare(2)
            endif
            call PrintScr('Final')
            call LogFile('Final')
        case default
            continue
        end select
    
end subroutine OutputData


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating the dynamical time step
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MediaTimeSteps

    integer Medm, MedMax, MedTmp
    integer Nstep

#ifdef DEBUG_MODE
    print *, '>- Entering MediaTimeSteps...'
#endif

! Calculation of timesteps
    Tstep=dTinput
    do Medm=1, NumMed
      selectcase(MedmLst(Medm))
#ifdef M_ATM
        case(Atm)
          call Atm_TstepCalc(Tstep(Medm))
          MedOrder(Medm)=Medm
#endif
#ifdef M_OCN
        case(Ocn)
          call Ocn_TstepCalc(Tstep(Medm))
          MedOrder(Medm)=Medm
#endif
#ifdef M_SOIL
        case(Soil)
          call Soil_TstepCalc(Tstep(Medm))
          MedOrder(Medm)=Medm
#endif
#ifdef M_VEG
          case(Veg)
          call Veg_TstepCalc(Tstep(Medm))
          MedOrder(Medm)=Medm
#endif
        case default
          print '(/,"STOP: Inappropriate medium number =",i2,/)', Medm
          stop
      endselect
    enddo
    
! Media ordering by time step
    do MedMax=NumMed-1, 1, -1
      do Medm=1, MedMax          
        if(Tstep(MedOrder(Medm))>Tstep(MedOrder(Medm+1))) then
          MedTmp=MedOrder(Medm)
          MedOrder(Medm)=MedOrder(Medm+1)
          MedOrder(Medm+1)=MedTmp        
        endif
      enddo
    enddo

! Adjusting time steps
    NumStep(MedOrder(NumMed))=ceiling(dTinput/Tstep(MedOrder(NumMed)))   
    Tstep(MedOrder(NumMed))=real(dTInput/NumStep(MedOrder(NumMed)))
    Nstep=NumStep(MedOrder(NumMed))
    do Medm=NumMed-1, 1, -1          
      NumStep(MedOrder(Medm))=ceiling(Tstep(MedOrder(Medm+1))/Tstep(MedOrder(Medm)))
      Nstep=Nstep*NumStep(MedOrder(Medm))
      Tstep(MedOrder(Medm))=real(Tstep(MedOrder(Medm+1))/NumStep(MedOrder(Medm)))
    enddo

! Adjusting media parameters
#ifdef M_ATM
    call Atm_Adjust
#endif
#ifdef M_OCN
!    call Ocn_Adjust
#endif
#ifdef M_SOIL
!    call Soil_Adjust
#endif    
#ifdef M_VEG
!    call Veg_Adjust
#endif

    DayTime=real(Period-1)*dTinput
   
! Printing information on the screen
    write(*,'(a14,3x,i1)') '  Meteoperiod:', Period
    write(*,'(a)') '   Medium   Steps    Timestep'
    Nstep=1
    do Medm=NumMed, 1, -1
      Nstep=Nstep*NumStep(MedOrder(Medm))
      write(*,'(3x,a4,5x,i4,4x,f6.2,1x,a3)') MedmID(MedOrder(Medm)), Nstep, Tstep(MedOrder(Medm))/60., 'min'
    enddo

#ifdef DEBUG_MODE
    print *, '<- Exit MediaTimeSteps'
#endif
end subroutine MediaTimeSteps


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine performing time integration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recursive subroutine TimeIntegration(Medm)

    integer Medm, iStep

#ifdef DEBUG_MODE
    print *, '>- Entering TimeIntegration... ', Medm
#endif

    do iStep=1, NumStep(MedOrder(Medm))
    
        if(Medm>1) call TimeIntegration(Medm-1)

        if(Medm==1) DayTime=DayTime+Tstep(MedOrder(Medm))

        call Media_Process(MedOrder(Medm))

        if(Medm==1) then

#ifdef G_POP
            call Media_Exchange_POP
#else
            call Media_Exchange
#endif

        endif
        call OutputData('Steply',MedOrder(Medm))          ! corrected 03.04.2015
    enddo

#ifdef DEBUG_MODE
    print *, '<- Exit TimeIntegration ', Medm
#endif
end subroutine TimeIntegration


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calling process procedures of different media
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Media_Process(Medm)

    integer Medm

#ifdef DEBUG_MODE
    print *, '>- Entering Media_Process... ', Medm
#endif

    selectcase(Medm)
#ifdef M_ATM
    case(Atm)
      call Atm_Process
#endif

#ifdef M_OCN
    case(Ocn)
      call Ocn_Process
#endif

#ifdef M_SOIL
    case(Soil)
      call Soil_Process
#endif

#ifdef M_VEG
    case(Veg)
      call Veg_Process
#endif
    case default
        print '(/,"STOP: Inappropriate medium number ''",i2,"''",/)', Medm
        stop
    endselect

#ifdef DEBUG_MODE
    print *, '<- Exit Media_Process ', Medm
#endif
end subroutine Media_Process


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading conditions of the calculation run
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RunInfo

    integer i, j, FileStat, lenType, lenPar, lenStr, Src, Grp, Subs, Medm, ios, ep, Lep
    real longM, latM
    character(800) strRead, strPar, strTemp
    character(90) strType, groupTmp(MaxSubs)
    character(4) dateRun(3)
    character(6) MedmTmp(MaxMed)
    character(10) :: SubsGroupDef(MaxGroups,MaxSubs)=''
    character(1) frstSymb

    character(10) epDate
    character(20) epName
    real epLat, epLong

    integer Ip, Jp, bn
    real ex, ey, M, Lproj
    logical wrongGroup

    RunName='./RunInfo.dat'
    open(1, file=trim(RunName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open RunInfo file ''",a,"''",/)', trim(RunName)
      stop
    endif

    do
      read(1,'(a)',iostat=ios) strRead
      if(ios==-1) exit

! Recognizing comments
      if(strRead(1:1)=='!'.or.strRead(1:1)==' '.or.strRead(1:1)=='	'.or.ichar(strRead(1:1))==13) cycle
      lenType=scan(strRead,'!')
      if(lenType>0) then
        strRead=strRead(:lenType-1)
      endif

! Deleting tabs and line termination symbols
      lenType=scan(strRead,'	')
      do while (lenType>0)
        strRead=strRead(:lenType-1)//strRead(lenType+1:)
        lenType=scan(strRead,'	')
      enddo      
      lenType=scan(strRead,achar(13))
      if(lenType>0) strRead=strRead(:lenType-1)

! Reading parameter type
      lenType=scan(strRead,':')
      if(lenType==0) then
        print '(/,"STOP: Wrong format of ''",a,"'' file",/)', trim(RunName)
        stop
      endif
      strType=trim(strRead(:lenType-1))
      strPar=strRead(lenType+1:)
      strPar=adjustl(strPar)

      selectcase(strType)
!------------------------------------------------------------------
        case('Pollutants number')
            read(strPar,'(i2)') NumSubs
!------------------------------------------------------------------
        case('Pollutants ID')
          read(strPar,*) (SubsID(i), i=1, NumSubs)
!------------------------------------------------------------------
        case('Number of media')
          read(strPar,'(i2)') NumMed
!------------------------------------------------------------------
        case('Media ID')
          read(strPar,*) (MedmTmp(i), i=1, NumMed)

          do i=1, NumMed
            do Medm=1, MaxMed
              if(MedmTmp(i)==MedmID(Medm)) exit
            enddo
            if(Medm<=MaxMed) then
              MedmLst(i)=Medm
            else
              print '(/,"STOP: Wrong media type ''",a,"''",/)', trim(MedmTmp(i))
              stop 
            endif
          enddo
!------------------------------------------------------------------
        case('Run type')
          RunType=strPar
          if(RunType=='field'.and.RtypeCompil==2.or.RunType=='matrix'.and.RtypeCompil==1) then
            print '(/,"STOP: Wrong Run type ''",a,"''",/)', trim(RunType)
            stop
          endif
!------------------------------------------------------------------
        case('Grid code')
          GridCode=strPar
!------------------------------------------------------------------
        case('Conditions')
          InitCond=strPar
          if(InitCond/='zero'.and.InitCond/='cond'.and.InitCond/='dump') then
            print '(/,"STOP: Wrong initial conditions ''",a,"''",/)', trim(InitCond)
            stop
          endif
!------------------------------------------------------------------
        case('Start')
          dateRun(da)=strPar(:2)
          dateRun(mn)=strPar(4:5)
          dateRun(yr)=strPar(7:)
                  
          read(dateRun(da),'(i2)') BegDate(da)
          read(dateRun(mn),'(i2)') BegDate(mn)
          read(dateRun(yr),'(i4)') BegDate(yr)
          Year = BegDate(yr)
!------------------------------------------------------------------
        case('Finish')
          dateRun(da)=strPar(:2)
          dateRun(mn)=strPar(4:5)
          dateRun(yr)=strPar(7:)
                  
          read(dateRun(da),'(i2)') FinDate(da)
          read(dateRun(mn),'(i2)') FinDate(mn)
          read(dateRun(yr),'(i4)') FinDate(yr)
!------------------------------------------------------------------
        case('Climatic meteo mode')
          selectcase(trim(strPar))
            case('yes')
              climRun=.true.
            case('no')
              climRun=.false.
            case default
              print '(/,"STOP: Incorrect type of climatic mode ''",a,"''",/)', trim(strPar)
              stop
          endselect
!------------------------------------------------------------------
        case('Climatic meteo year')
          read(strPar,'(i4)') ClimYear
!------------------------------------------------------------------
        case('Climatic reactants mode')
          selectcase(trim(strPar))
            case('yes')
              climReactRun=.true.
            case('no')
              climReactRun=.false.
            case default
              print '(/,"STOP: Incorrect type of climatic mode ''",a,"''",/)', trim(strPar)
              stop
          endselect
!------------------------------------------------------------------
        case('Climatic reactants year')
          read(strPar,'(i4)') ClimReactYear
!------------------------------------------------------------------
        case('Climatic land cover mode')
          selectcase(trim(strPar))
            case('yes')
              climLCRun=.true.
            case('no')
              climLCRun=.false.
            case default
              print '(/,"STOP: Incorrect type of LC climatic mode ''",a,"''",/)', trim(strPar)
              stop
          endselect
!------------------------------------------------------------------
        case('Climatic land cover year')
          read(strPar, '(i4)') climLCYear
!------------------------------------------------------------------
        case('Ant sources mode')
          selectcase(trim(strPar))
            case('all')
              antSRCmode=1
            case('selected')
              antSRCmode=2
            case('none')
              antSRCmode=0
            case default
              print '(/,"STOP: Incorrect type of anthropogenic sources mode ''",a,"''",/)', trim(strPar)
              stop
          endselect
!------------------------------------------------------------------
        case('Ant sources number')
          if(antSRCmode==2) read(strPar,'(i2)') NumAnth
!------------------------------------------------------------------
        case('Ant sources codes')
          if(antSRCmode==2) read(strPar,*) (AnthrID(i), i=1, NumAnth)
!------------------------------------------------------------------
        case('Nat sources mode')
          selectcase(trim(strPar))
            case('all')
              natSRCmode=1
            case('selected')
              natSRCmode=2
            case('none')
              natSRCmode=0
            case default
              print '(/,"STOP: Incorrect type of natural sources mode ''",a,"''",/)', trim(strPar)
              stop
          endselect
!------------------------------------------------------------------
        case('Nat sources number')
          if(natSRCmode==2) read(strPar,'(i2)') NumNat
!------------------------------------------------------------------
        case('Nat sources codes')
          if(natSRCmode==2) read(strPar,*) (NaturID(i), i=1, NumNat)
!------------------------------------------------------------------
        case('Bnd sources mode')
          selectcase(trim(strPar))
            case('all')
              bndSRCmode=1
            case('selected')
              bndSRCmode=2
            case('none')
              bndSRCmode=0
            case default
              print '(/,"STOP: Incorrect type of boundary sources mode ''",a,"''",/)', trim(strPar)
              stop
          endselect
!------------------------------------------------------------------
        case('Bnd sources number')
          if(bndSRCmode==2) read(strPar,'(i2)') NumBnd
!------------------------------------------------------------------
        case('Bnd sources codes')
          if(bndSRCmode==2) read(strPar,*) (BoundID(i), i=1, NumBnd)
!------------------------------------------------------------------
        case('Receptors mode')
          selectcase(trim(strPar))
            case('all')
              rcpSRCmode=1
            case('selected')
              rcpSRCmode=2
            case('none')
              rcpSRCmode=0
            case default
              print '(/,"STOP: Incorrect type of receptors mode ''",a,"''",/)', trim(strPar)
              stop
          endselect
!------------------------------------------------------------------
        case('Receptors number')
          if(rcpSRCmode==2) read(strPar,'(i2)') NumRcp
!------------------------------------------------------------------
        case('Receptors codes')
          if(rcpSRCmode==2) read(strPar,*) (RecepID(i), i=1, NumRcp)
!------------------------------------------------------------------
        case('Initial cond mode')
          selectcase(trim(strPar))
            case('single')
              initSRCmode=1
            case('multi')
              initSRCmode=2
            case default
              print '(/,"STOP: Incorrect type of initial conditions mode ''",a,"''",/)', trim(strPar)
              stop
        endselect
!------------------------------------------------------------------
        case('Re-emission mode')
          selectcase(trim(strPar))
            case('single')
              reemisSRCmode=1
            case('multi')
              reemisSRCmode=2
            case default
              print '(/,"STOP: Incorrect type of initial conditions mode ''",a,"''",/)', trim(strPar)
              stop
        endselect
!------------------------------------------------------------------
        case('Emission path')
          EmisPath=strPar
!------------------------------------------------------------------
        case('Emission dataset')
          AntEmisName=strPar
!------------------------------------------------------------------
        case('Natur emis file')
          NatEmisName=strPar
!------------------------------------------------------------------
        case('Meteo path')
          MeteoPath=strPar
!------------------------------------------------------------------
        case('Config path')
          ConfigPath=strPar
!------------------------------------------------------------------
        case('Grid config file')
          GridConfName=strPar
!------------------------------------------------------------------
        case('Media config file')
          MediaConfName=strPar
!------------------------------------------------------------------
        case('LC config file')
          LCconfName=strPar
!------------------------------------------------------------------
        case('Matrix config file')
          MatrixConfName=strPar
!------------------------------------------------------------------
        case('Output config file')
          OutputConfigName=strPar
!------------------------------------------------------------------
        case('Geodata path')
          GeoPath=strPar
!------------------------------------------------------------------
        case('Common geophys')
            ComGeoName=strPar
!------------------------------------------------------------------
        case('Orography file')
          OrogName=strPar
!------------------------------------------------------------------
        case('Land Cover file')
          LandName=strPar
!------------------------------------------------------------------
        case('LandCover source')
          LandSource=strPar
!------------------------------------------------------------------
        case('Season file')
            SeasonName=strPar
!------------------------------------------------------------------
        case('Roughness file')
            RoughName=strPar
!------------------------------------------------------------------
        case('Properties path')
          PropPath=strPar
!------------------------------------------------------------------
        case('Properties file')
          PropName=strPar
!------------------------------------------------------------------
        case('Groups file')
          GroupName=strPar
!------------------------------------------------------------------
        case('Receptors path')
          RecepPath=strPar
!------------------------------------------------------------------
        case('Receptors file')
          RecepArName=strPar
!------------------------------------------------------------------
        case('Stations path')
          StatPath=strPar
!------------------------------------------------------------------
        case('Stations file')
          StatName=strPar
!------------------------------------------------------------------
        case('Boundary path')
          InBoundPath=strPar
!------------------------------------------------------------------
        case('Boundary file')
          BoundName=strPar
!------------------------------------------------------------------
        case('Init cond path')
          InitPath=strPar
!------------------------------------------------------------------
        case('Init file')
          InitName=strPar
!------------------------------------------------------------------
        case('Dust data path')
          DustPath=strPar
!------------------------------------------------------------------
        case('Dust file')
          DustName=strPar
!------------------------------------------------------------------
        case('Soil data path')
          SoilPath=strPar
!------------------------------------------------------------------
        case('Soil conc file')
          SoilConcName=strPar
!------------------------------------------------------------------
        case('Reactant source')
          ReactSource=strPar
!------------------------------------------------------------------
        case('Reactant path 1')
          ReactPath1=strPar
!------------------------------------------------------------------
        case('Reactant path 2')
          ReactPath2=strPar
!------------------------------------------------------------------
        case('Ozone conc file')
          OzoneName=strPar
!------------------------------------------------------------------
        case('S02 conc file')
          SO2Name=strPar
!------------------------------------------------------------------
        case('OH conc file')
          OHName=strPar
!------------------------------------------------------------------
        case('LandCover path')
          LandcoverPath=strPar
!------------------------------------------------------------------
#ifdef M_OCN
        case('Ocean path')  
          OceanPath=strPar  
!------------------------------------------------------------------
        case('Ocean topo file')
          TopoOcnName=strPar
!------------------------------------------------------------------
        case('Ocean fract file')
          LC2OcnName=strPar
#endif
!------------------------------------------------------------------
#ifdef M_SOIL
        case('Soil types file')
          LC2SoilTypeName=strPar
        case('FOC file')
          FocName=strPar
#endif
!------------------------------------------------------------------
        endselect
    enddo
    close(1)

    ReactPath1=trim(ReactPath1)//trim(ReactSource)//'/'

! Defining pollutant groups
    fileName=trim(GroupName)
    fullName=trim(ConfigPath)//fileName
    open(2, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open pollutant groups file ''",a,"''",/)', trim(fullName)
      stop
    endif

    do
      read(2,'(a)',iostat=ios) strRead
      if(ios==-1) exit

! Recognizing comments
      if(strRead(1:1)=='!'.or.strRead(1:1)==' '.or.strRead(1:1)=='	') cycle
      lenType=scan(strRead,'!')
      if(lenType>0) then
        strRead=strRead(:lenType-1)
      endif

! Deleting tabs
      lenType=scan(strRead,'    ')
	  do while (lenType>0)
        strRead=strRead(:lenType-1)//strRead(lenType+1:)
        lenType=scan(strRead,'	')
      enddo      

! Reading group definition
      lenType=scan(strRead,':')
      if(lenType==0) then
        print '(/,"STOP: Wrong format of file ''",a,"''",/)', trim(fileName)
        stop 
      endif
      strType=trim(strRead(:lenType-1))

      do Grp=1, MaxGroups
        if(strType==SubsGroupID(Grp)) exit
      enddo
      if(Grp>MaxGroups) then
        print '(/,"STOP: Wrong pollutant group type ''",a,"'' in file ''",a,"''",/)', trim(strType), trim(fileName)
        stop 
      endif

      strPar=strRead(lenType+1:)
      strPar=adjustl(strPar)
      read(strPar,*,iostat=ios) (SubsGroupDef(Grp,Subs), Subs=1, MaxSubs)
    enddo

! Identifying pollutants and groups in the current model run
    NumGroups=0          
    gSubsNum=0
    do i=1, NumSubs
      do Grp=1, MaxGroups
        do Subs=1, MaxSubs      
         if(SubsID(i)==SubsGroupDef(Grp,Subs)) then
            SubsGroup(i)=Grp
            exit
          endif        
        enddo
        if(Subs<=MaxSubs) exit
      enddo
      if(Grp>MaxGroups) then
        print '(/,"STOP: Wrong pollutant type ''",a,"'' in the list",/)', trim(SubsID(i))
        stop 
      endif

      do Grp=1, NumGroups
        if(SubsGrpLst(Grp)==SubsGroup(i)) exit
      enddo
      if(Grp>NumGroups) then
        SubsGrpLst(Grp)=SubsGroup(i)
        NumGroups=Grp
        GroupInd(SubsGroup(i))=Grp                   
      endif
      gSubsNum(Grp)=gSubsNum(Grp)+1
      gSubsInd(Grp,gSubsNum(Grp))=i        
      gSubsGroupInd(i)=gSubsNum(Grp)           
    enddo

    wrongGroup=.true.
    do Grp=1, NumGroups
      selectcase(SubsGrpLst(Grp))
#ifdef G_HM
        case(HM)
          wrongGroup=.false.
#endif 
#ifdef G_HG
        case(HG)
          wrongGroup=.false.
#endif
#ifdef G_POP
        case(POP)
          wrongGroup=.false.
#endif
#ifdef G_TRACER
        case(TRACER)
          wrongGroup=.false.
#endif
#ifdef G_AERO
        case(AERO)
          wrongGroup=.false.
#endif
      endselect
      if(wrongGroup) then
        print '(/,"STOP: Wrong pollutant group of the model run",/)'
        stop
      endif
    enddo
   
end subroutine RunInfo


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading model grid configuration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GridConfig

    integer FileStat, lenType, ios, maxI, maxJ, Xres, Yres
    character(800) strRead, strPar
    character(90) strType

#ifdef DEBUG_MODE
    print *, '>- Entering GridConfig...'
#endif

    fullName=trim(ConfigPath)//trim(GridConfName)//trim(GridCode)//'.dat'
    open(1, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open grid config file ''",a,"''",/)', trim(fullName)
      stop
    endif

    do
      read(1,'(a)',iostat=ios) strRead
      if(ios==-1) exit

! Recognizing comments
      if(strRead(1:1)=='!'.or.strRead(1:1)==' '.or.strRead(1:1)=='	'.or.ichar(strRead(1:1))==13) cycle
      lenType=scan(strRead,'!')
      if(lenType>0) then
        strRead=strRead(:lenType-1)
      endif

! Deleting tabs and line termination symbols
      lenType=scan(strRead,'	')
      do while (lenType>0)
        strRead=strRead(:lenType-1)//strRead(lenType+1:)
        lenType=scan(strRead,'	')
      enddo
      lenType=scan(strRead,achar(13))
      if(lenType>0) strRead=strRead(:lenType-1)

! Reading parameter type
      lenType=scan(strRead,':')
      if(lenType==0) then
        print '(/,"STOP: Wrong format of ''",a,"'' file",/)', trim(fullName)
        stop
      endif
      strType=trim(strRead(:lenType-1))
      strPar=strRead(lenType+1:)
      strPar=adjustl(strPar)

      selectcase(strType)
!------------------------------------------------------------------
        case('Grid code')
          if(trim(strPar)/=trim(GridCode)) then
            print '(/,"STOP: Wrong grid code - ''",a,"''",/)', trim(strPar)
            stop
          endif
!------------------------------------------------------------------
        case('Projection type')
          if(trim(strPar)=='latlong'.and.PROJTYPE/=1.or.trim(strPar)=='polar'.and.PROJTYPE/=2) then
            print '(/,"STOP: Wrong projection type - ''",a,"''",/)', trim(strPar)
            stop
          endif
!------------------------------------------------------------------
        case('Region type')
          if(trim(strPar)=='global'.and.REGTYPE/=1.or.trim(strPar)=='regional'.and.REGTYPE/=2) then
            print '(/,"STOP: Wrong region type - ''",a,"''",/)', trim(strPar)
            stop
          endif
!------------------------------------------------------------------
        case('Grid size')
          read(strPar,*) maxI, maxJ             ! Number of grid cells
          if(maxI/=Imax.or.maxJ/=Jmax) then
            print '(/,"STOP: Wrong grid size - ",a,/)', trim(strPar)
            stop
          endif
!------------------------------------------------------------------
        case('Grid resolution')
          read(strPar,*) dXstep, dYstep
          dXmin=dXstep*pi180
          dY=dYstep*pi180
          jGlob=nint(180./dYstep)
!------------------------------------------------------------------
        case('Grid origin')
          read(strPar,*) xOrig, yOrig           ! Degrees
          xOrig=xOrig*pi180
          yOrig=yOrig*pi180
!------------------------------------------------------------------
        case('Min step')
          read(strPar,*) dLmin                  ! Meters
!------------------------------------------------------------------
      endselect
    enddo
    close(1)

#ifdef DEBUG_MODE
    print *, '<- Exit GridConfig'
#endif
end subroutine GridConfig


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading configuration of land cover data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LC_config

    integer FileStat, lenType, ios, Surf, Surf1, Nsurf
    character(800) strRead, strPar
    character(300) LCfillConfName
    character(90) strType, TypeSurf(MaxSurf)
    
#ifdef DEBUG_MODE
    print *, '>- Entering LC_config...'
#endif

    LCfillConfName=trim(LCconfName)//trim(LandSource)//'_config.dat'
    fullName=trim(ConfigPath)//LCfillConfName
    open(1, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open land cover config file ''",a,"''",/)', trim(fullName)
      stop
    endif

    gNum=0
    do
      read(1,'(a)',iostat=ios) strRead
      if(ios==-1) exit

! Recognizing comments
      if(strRead(1:1)=='!'.or.strRead(1:1)==' '.or.strRead(1:1)=='	'.or.ichar(strRead(1:1))==13) cycle
      lenType=scan(strRead,'!')
      if(lenType>0) then
        strRead=strRead(:lenType-1)
      endif

! Deleting tabs and line termination symbols
      lenType=scan(strRead,'	')
      do while (lenType>0)
        strRead=strRead(:lenType-1)//strRead(lenType+1:)
        lenType=scan(strRead,'	')
      enddo
      lenType=scan(strRead,achar(13))
      if(lenType>0) strRead=strRead(:lenType-1)

! Reading parameter type
      lenType=scan(strRead,':')
      if(lenType==0) then
        print '(/,"STOP: Wrong format of ''",a,"'' file",/)', trim(fullName)
        stop
      endif
      strType=trim(strRead(:lenType-1))
      strPar=strRead(lenType+1:)
      strPar=adjustl(strPar)

      selectcase(strType)
!------------------------------------------------------------------
        case('Surfaces number')
          read(strPar,'(i2)') NumSurf
!------------------------------------------------------------------
        case('LC codes')
          read(strPar,*) (SurfType(Surf), Surf=1, NumSurf)
!------------------------------------------------------------------
        case('Water')
          call ReadList(strPar,TypeSurf,Nsurf)

          do Surf1=1, Nsurf
            do Surf=1, NumSurf
              if(TypeSurf(Surf1)==SurfType(Surf)) then
                SurfGroup(Surf)=Water
                gNum(SurfGroup(Surf))=gNum(SurfGroup(Surf))+1
                gInd(SurfGroup(Surf),gNum(SurfGroup(Surf)))=Surf
                exit
              endif    
            enddo    
          enddo    
!------------------------------------------------------------------
        case('Forest')
          call ReadList(strPar,TypeSurf,Nsurf)

          do Surf1=1, Nsurf
            do Surf=1, NumSurf
              if(TypeSurf(Surf1)==SurfType(Surf)) then
                SurfGroup(Surf)=Forest
                gNum(SurfGroup(Surf))=gNum(SurfGroup(Surf))+1
                gInd(SurfGroup(Surf),gNum(SurfGroup(Surf)))=Surf
                exit
              endif
            enddo
          enddo
!------------------------------------------------------------------
        case('Grass')
          call ReadList(strPar,TypeSurf,Nsurf)

          do Surf1=1, Nsurf
            do Surf=1, NumSurf
              if(TypeSurf(Surf1)==SurfType(Surf)) then
                SurfGroup(Surf)=Grass
                gNum(SurfGroup(Surf))=gNum(SurfGroup(Surf))+1
                gInd(SurfGroup(Surf),gNum(SurfGroup(Surf)))=Surf
                exit
              endif
            enddo
          enddo
!------------------------------------------------------------------
        case('Arable')
          call ReadList(strPar,TypeSurf,Nsurf)

          do Surf1=1, Nsurf
            do Surf=1, NumSurf
              if(TypeSurf(Surf1)==SurfType(Surf)) then
                SurfGroup(Surf)=Arable
                gNum(SurfGroup(Surf))=gNum(SurfGroup(Surf))+1
                gInd(SurfGroup(Surf),gNum(SurfGroup(Surf)))=Surf
                exit
              endif
            enddo
          enddo
!------------------------------------------------------------------
        case('Urban')
          call ReadList(strPar,TypeSurf,Nsurf)

          do Surf1=1, Nsurf
            do Surf=1, NumSurf
              if(TypeSurf(Surf1)==SurfType(Surf)) then
                SurfGroup(Surf)=Urban
                gNum(SurfGroup(Surf))=gNum(SurfGroup(Surf))+1
                gInd(SurfGroup(Surf),gNum(SurfGroup(Surf)))=Surf
                exit
              endif
            enddo
          enddo
!------------------------------------------------------------------
        case('Barren')
          call ReadList(strPar,TypeSurf,Nsurf)

          do Surf1=1, Nsurf
            do Surf=1, NumSurf
              if(TypeSurf(Surf1)==SurfType(Surf)) then
                SurfGroup(Surf)=Barren
                gNum(SurfGroup(Surf))=gNum(SurfGroup(Surf))+1
                gInd(SurfGroup(Surf),gNum(SurfGroup(Surf)))=Surf
                exit
              endif
            enddo
          enddo
!------------------------------------------------------------------
      endselect
    enddo
    close(1)
    
#ifdef DEBUG_MODE
    print *, '<- Exit LC_config'
#endif
 
contains 

subroutine ReadList(str,List,Nitem)

  integer Nitem, n, pos
  character(800) str
  character(90) List(MaxSurf)

  n=0
  pos=scan(str,',')
  do while (pos>0)
    n=n+1
    List(n)=str(:pos-1)
    str=str(pos+1:)
    str=adjustl(str)
    pos=scan(str,',')
  enddo
  n=n+1
  List(n)=str
  Nitem=n
  
end subroutine ReadList    
    
end subroutine LC_config


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading model output configuration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadOutputConfig

    integer FileStat, lenType, ios
    character(800) strRead, strPar
    character(90) strType
    
#ifdef DEBUG_MODE
    print *, '>- Entering ReadOutputConfig...'
#endif

    call InitOutputConfig

    open(1, file=trim(OutputConfigName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open output config file ''",a,"''",/)', trim(OutputConfigName)
      stop
    endif

    do
      read(1,'(a)',iostat=ios) strRead
      if(ios==-1) exit

! Recognizing comments
      if(strRead(1:1)=='!'.or.strRead(1:1)==' '.or.strRead(1:1)=='	'.or.ichar(strRead(1:1))==13) cycle
      lenType=scan(strRead,'!')
      if(lenType>0) then
        strRead=strRead(:lenType-1)
      endif

! Deleting tabs and line termination symbols
      lenType=scan(strRead,'	')
      do while (lenType>0)
        strRead=strRead(:lenType-1)//strRead(lenType+1:)
    	lenType=scan(strRead,'	')
      enddo
      lenType=scan(strRead,achar(13))
      if(lenType>0) strRead=strRead(:lenType-1)

! Reading parameter type
      lenType=scan(strRead,':')
      if(lenType==0) then
        print '(/,"STOP: Wrong format of ''",a,"'' file",/)', trim(RunName)
        stop
      endif
      strType=trim(strRead(:lenType-1))
      strPar=strRead(lenType+1:)
      strPar=adjustl(strPar)

      selectcase(strType)

!--- Output directories control and log file name
        case('Output Dir')
          OutPath=trim(strPar)//'/'
        case('Log file')
          LogName=strPar
!--- Directory for output dump and dump file name
        case('Dump Dir')
          DumpDir=trim(strPar)//'/'
        case('Dump file')
          DumpName=strPar
!--- Directory for balance files
        case('Balance Dir')
          BalanceDir=trim(strPar)//'/'
!--- Dump frequency control (annual, monthly, daily)
        case('Dump frequency')
          if(strPar == 'yearly') then
            OutputDumpYearly=.true.
          else if(strPar == 'monthly') then
            OutputDumpMonthly=.true.
          else if(strPar == 'daily') then
            OutputDumpDaily=.true.
          else
            print *, 'Error. Dump output time scale is not defined. Stop.'
            stop
          endif

!--- Dump frequency control (annual, monthly, daily)
        case('Balance frequency')
          if(strPar == 'yearly') then
            OutputBalanceYearly=.true.
          else if(strPar == 'monthly') then
            OutputBalanceMonthly=.true.
          else if(strPar == 'daily') then
            OutputBalanceDaily=.true.
          else if(strPar == '6hourly') then
            OutputBalance6hourly=.true.
          else
            print *, 'Error. Balance output time scale is not defined. Stop.'
            stop
          endif
                 
!--- Media output directories and specific output files
#ifdef M_ATM
          AtmDir=trim(MedmID(Atm))//'/'
!        case('Atm Dir')
!          AtmDir=trim(strPar)//'/'
        case('Dry depos file')
          DryName=strPar
        case('Wet depos file')
          WetName=strPar
        case('Total depos file')
          TotName=strPar
        case('Net depos file')
          NetName=strPar
        case('Remobiliz file')
          RemName=strPar
        case('Mean conc file')
          ConcName=strPar
        case('Mix ratio file')
          MixRatName=strPar
        case('In precip file')
          InPrecName=strPar
        case('Precip file')
          PrecipName=strPar
        case('Matrix tot dep file')
          MatrixTotDepName=strPar
        case('Matrix wet dep file')
          MatrixWetDepName=strPar
        case('Matrix dry dep file')
          MatrixDryDepName=strPar
        case('Matrix conc file')
          MatrixConcName=strPar
        case('Out ant file')
          EmisOutAnt=strPar
        case('Out nat file')
          EmisOutNat=strPar
        case('Out re-emis file')
          EmisOutRem=strPar
        case('NCF output file')
          NCFOutName=strPar
        case('Conc monitor')
          ConcMonitName=strPar
        case('Inprec monitor')
          InprecMonitName=strPar
        case('Prec monitor')
          PrecMonitName=strPar
        case('Flux monitor')
          FluxMonitName=strPar
        case('Dry monitor')
          DryMonitName=strPar
        case('Wet monitor')
          WetMonitName=strPar
        case('NCF monitor')
          NCFMonitName=strPar
#endif
#ifdef M_OCN
        OcnDir=trim(MedmID(Ocn))//'/'
!        case('Ocean Dir')
!          OcnDir=trim(strPar)//'/'
        case('Ocn conc file')
          ConcOcnName=strPar
#endif
#ifdef M_SOIL
        SoilDir=trim(MedmID(Soil))//'/'
!        case('Soil Dir')
!          SoilDir=trim(strPar)//'/'
        case('Soil conc file')
          ConcSoilName=strPar
#endif
#ifdef M_VEG
        VegDir=trim(MedmID(Veg))//'/'
!        case('Vegetation Dir')
!          VegDir=trim(strPar)//'/'
        case('Veg conc file')
          ConcVegName=strPar
#endif

!--- Specific directories
        case('Fields Dir')
          FieldsDir=trim(strPar)//'/'
        case('Fields NCF Dir')
          FieldsNCFDir=trim(strPar)//'/'
        case('Monitor Dir')
          MonitorDir=trim(strPar)//'/'
        case('Matrix Dir')
          MatrixDir=trim(strPar)//'/'

!--- Specific temporal variations directories
        case('Monthly Dir')
          MonthlyDir=trim(strPar)//'/'
        case('Daily Dir')
          DailyDir=trim(strPar)//'/'
        case('Yearly Dir')
          YearlyDir=trim(strPar)//'/'
        case('6hourly Dir')
          SixHourlyDir=trim(strPar)//'/'
        case('Hourly Dir')
          HourlyDir=trim(strPar)//'/'

!--- Output time scale control by media

#ifdef M_ATM
!--- Monitoring output control
        case('Atm Monitor Yearly')
          if(strPar == 'on')      AtmOutMonitYearly=.true.
        case('Atm Monitor Monthly')
          if(strPar == 'on')      AtmOutMonitMonthly=.true.
        case('Atm Monitor Daily')
          if(strPar == 'on')      AtmOutMonitDaily=.true.
        case('Atm Monitor 6hourly')
          if(strPar == 'on')      AtmOutMonit6hourly=.true.
        case('Atm Monitor Hourly')
          if(strPar == 'on')      AtmOutMonitHourly=.true.
!--- Matrix output control
        case('Atm Matrix Yearly')
          if(strPar == 'on')      AtmOutMatrixYearly=.true.
        case('Atm Matrix Monthly')
          if(strPar == 'on')      AtmOutMatrixMonthly=.true.
        case('Atm Matrix Daily')
          if(strPar == 'on')      AtmOutMatrixDaily=.true.
!--- Fields output in text files control
        case('Atm Fields Yearly')
          if(strPar == 'on')      AtmOutFieldsYearly=.true.
        case('Atm Fields Monthly')
          if(strPar == 'on')      AtmOutFieldsMonthly=.true.
        case('Atm Fields Daily')
          if(strPar == 'on')      AtmOutFieldsDaily=.true.
        case('Atm Fields 6hourly')
          if(strPar == 'on')      AtmOutFields6hourly=.true.
        case('Atm Fields Hourly')
          if(strPar == 'on')      AtmOutFieldsHourly=.true.
!--- NetCDF output control
        case('Atm NCF Yearly')
          if(strPar == 'on')      AtmOutNCFYearly=.true.
        case('Atm NCF Monthly')
          if(strPar == 'on')      AtmOutNCFMonthly=.true.
        case('Atm NCF Daily')
          if(strPar == 'on')      AtmOutNCFDaily=.true.
        case('Atm NCF 6hourly')
          if(strPar == 'on')      AtmOutNCF6hourly=.true.
        case('Atm NCF Hourly')
          if(strPar == 'on')      AtmOutNCFHourly=.true.
#endif

#ifdef M_OCN
!--- Monitoring output control
        case('Ocn Monitor Yearly')
          if(strPar == 'on')      OcnOutMonitYearly=.true.
        case('Ocn Monitor Monthly')
          if(strPar == 'on')      OcnOutMonitMonthly=.true.
        case('Ocn Monitor Daily')
          if(strPar == 'on')      OcnOutMonitDaily=.true.
        case('Ocn Monitor 6hourly')
          if(strPar == 'on')      OcnOutMonit6hourly=.true.
        case('Ocn Monitor Hourly')
          if(strPar == 'on')      OcnOutMonitHourly=.true.
!--- Matrix output control
        case('Ocn Matrix Yearly')
          if(strPar == 'on')      OcnOutMatrixYearly=.true.
        case('Ocn Matrix Monthly')
          if(strPar == 'on')      OcnOutMatrixMonthly=.true.
        case('Ocn Matrix Daily')
          if(strPar == 'on')      OcnOutMatrixDaily=.true.
!--- Fields output in text files control
        case('Ocn Fields Yearly')
          if(strPar == 'on')      OcnOutFieldsYearly=.true.
        case('Ocn Fields Monthly')
          if(strPar == 'on')      OcnOutFieldsMonthly=.true.
        case('Ocn Fields Daily')
          if(strPar == 'on')      OcnOutFieldsDaily=.true.
        case('Ocn Fields 6hourly')
          if(strPar == 'on')      OcnOutFields6hourly=.true.
        case('Ocn Fields Hourly')
          if(strPar == 'on')      OcnOutFieldsHourly=.true.
!--- NetCDF output control
        case('Ocn NCF Yearly')
          if(strPar == 'on')      OcnOutNCFYearly=.true.
        case('Ocn NCF Monthly')
          if(strPar == 'on')      OcnOutNCFMonthly=.true.
        case('Ocn NCF Daily')
          if(strPar == 'on')      OcnOutNCFDaily=.true.
        case('Ocn NCF 6hourly')
          if(strPar == 'on')      OcnOutNCF6hourly=.true.
        case('Ocn NCF Hourly')
          if(strPar == 'on')      OcnOutNCFHourly=.true.
#endif

#ifdef M_SOIL
!--- Monitoring output control
        case('Soil Monitor Yearly')
          if(strPar == 'on')      SoilOutMonitYearly=.true.
        case('Soil Monitor Monthly')
          if(strPar == 'on')      SoilOutMonitMonthly=.true.
        case('Soil Monitor Daily')
          if(strPar == 'on')      SoilOutMonitDaily=.true.
        case('Soil Monitor 6hourly')
          if(strPar == 'on')      SoilOutMonit6hourly=.true.
        case('Soil Monitor Hourly')
          if(strPar == 'on')      SoilOutMonitHourly=.true.
!--- Matrix output control
        case('Soil Matrix Yearly')
          if(strPar == 'on')      SoilOutMatrixYearly=.true.
        case('Soil Matrix Monthly')
          if(strPar == 'on')      SoilOutMatrixMonthly=.true.
        case('Soil Matrix Daily')
          if(strPar == 'on')      SoilOutMatrixDaily=.true.
!--- Fields output in text files control
        case('Soil Fields Yearly')
          if(strPar == 'on')      SoilOutFieldsYearly=.true.
        case('Soil Fields Monthly')
          if(strPar == 'on')      SoilOutFieldsMonthly=.true.
        case('Soil Fields Daily')
          if(strPar == 'on')      SoilOutFieldsDaily=.true.
        case('Soil Fields 6hourly')
          if(strPar == 'on')      SoilOutFields6hourly=.true.
        case('Soil Fields Hourly')
          if(strPar == 'on')      SoilOutFieldsHourly=.true.
!--- NetCDF output control
        case('Soil NCF Yearly')
          if(strPar == 'on')      SoilOutNCFYearly=.true.
        case('Soil NCF Monthly')
         if(strPar == 'on')      SoilOutNCFMonthly=.true.
        case('Soil NCF Daily')
          if(strPar == 'on')      SoilOutNCFDaily=.true.
        case('Soil NCF 6hourly')
          if(strPar == 'on')      SoilOutNCF6hourly=.true.
        case('Soil NCF Hourly')
          if(strPar == 'on')      SoilOutNCFHourly=.true.
#endif

#ifdef M_VEG
!--- Monitoring output control
        case('Veg Monitor Yearly')
          if(strPar == 'on')      VegOutMonitYearly=.true.
        case('Veg Monitor Monthly')
          if(strPar == 'on')      VegOutMonitMonthly=.true.
        case('Veg Monitor Daily')
          if(strPar == 'on')      VegOutMonitDaily=.true.
        case('Veg Monitor 6hourly')
          if(strPar == 'on')      VegOutMonit6hourly=.true.
        case('Veg Monitor Hourly')
          if(strPar == 'on')      VegOutMonitHourly=.true.
!--- Matrix output control
        case('Veg Matrix Yearly')
          if(strPar == 'on')      VegOutMatrixYearly=.true.
        case('Veg Matrix Monthly')
          if(strPar == 'on')      VegOutMatrixMonthly=.true.
        case('Veg Matrix Daily')
          if(strPar == 'on')      VegOutMatrixDaily=.true.
!--- Fields output in text files control
        case('Veg Fields Yearly')
          if(strPar == 'on')      VegOutFieldsYearly=.true.
        case('Veg Fields Monthly')
          if(strPar == 'on')      VegOutFieldsMonthly=.true.
        case('Veg Fields Daily')
          if(strPar == 'on')      VegOutFieldsDaily=.true.
        case('Veg Fields 6hourly')
          if(strPar == 'on')      VegOutFields6hourly=.true.
        case('Veg Fields Hourly')
          if(strPar == 'on')      VegOutFieldsHourly=.true.
!--- NetCDF output control
        case('Veg NCF Yearly')
          if(strPar == 'on')      VegOutNCFYearly=.true.
        case('Veg NCF Monthly')
          if(strPar == 'on')      VegOutNCFMonthly=.true.
        case('Veg NCF Daily')
          if(strPar == 'on')      VegOutNCFDaily=.true.
        case('Veg NCF 6hourly')
          if(strPar == 'on')      VegOutNCF6hourly=.true.
        case('Veg NCF Hourly')
          if(strPar == 'on')      VegOutNCFHourly=.true.
#endif
          
      endselect
    enddo
    close(1)
    
#ifdef DEBUG_MODE
    print *, '<- Exit ReadOutputConfig'
#endif

end subroutine ReadOutputConfig

subroutine InitOutputConfig

    OutputDumpYearly  = .false.
    OutputDumpMonthly  = .false.
    OutputDumpDaily  = .false.

#ifdef M_ATM
    AtmOutMonitYearly = .false.
    AtmOutMonitMonthly = .false.
    AtmOutMonitDaily = .false.
    AtmOutMonit6hourly = .false.
    AtmOutMonitHourly = .false.
    AtmOutMatrixYearly = .false.
    AtmOutMatrixMonthly = .false.
    AtmOutMatrixDaily = .false.
    AtmOutFieldsYearly = .false.
    AtmOutFieldsMonthly = .false.
    AtmOutFieldsDaily = .false.
    AtmOutFields6hourly = .false.
    AtmOutFieldsHourly = .false.
    AtmOutNCFYearly = .false.
    AtmOutNCFMonthly = .false.
    AtmOutNCFDaily = .false.
    AtmOutNCF6hourly = .false.
    AtmOutNCFHourly = .false.
#endif
#ifdef M_OCN
    OcnOutMonitYearly = .false.
    OcnOutMonitMonthly = .false.
    OcnOutMonitDaily = .false.
    OcnOutMonit6hourly = .false.
    OcnOutMonitHourly = .false.
    OcnOutMatrixYearly = .false.
    OcnOutMatrixMonthly = .false.
    OcnOutMatrixDaily = .false.
    OcnOutFieldsYearly = .false.
    OcnOutFieldsMonthly = .false.
    OcnOutFieldsDaily = .false.
    OcnOutFields6hourly = .false.
    OcnOutFieldsHourly = .false.
    OcnOutNCFYearly = .false.
    OcnOutNCFMonthly = .false.
    OcnOutNCFDaily = .false.
    OcnOutNCF6hourly = .false.
    OcnOutNCFHourly = .false.
#endif
#ifdef M_SOIL
    SoilOutMonitYearly  = .false.
    SoilOutMonitMonthly  = .false.
    SoilOutMonitDaily  = .false.
    SoilOutMonit6hourly  = .false.
    SoilOutMonitHourly = .false.
    SoilOutMatrixYearly  = .false.
    SoilOutMatrixMonthly  = .false.
    SoilOutMatrixDaily = .false.
    SoilOutFieldsYearly  = .false.
    SoilOutFieldsMonthly  = .false.
    SoilOutFieldsDaily  = .false.
    SoilOutFields6hourly  = .false.
    SoilOutFieldsHourly = .false.
    SoilOutNCFYearly  = .false.
    SoilOutNCFMonthly  = .false.
    SoilOutNCFDaily  = .false.
    SoilOutNCF6hourly  = .false.
    SoilOutNCFHourly = .false.
#endif
#ifdef M_VEG
    VegOutMonitYearly  = .false.
    VegOutMonitMonthly  = .false.
    VegOutMonitDaily  = .false.
    VegOutMonit6hourly  = .false.
    VegOutMonitHourly = .false.
    VegOutMatrixYearly  = .false.
    VegOutMatrixMonthly  = .false.
    VegOutMatrixDaily = .false.
    VegOutFieldsYearly  = .false.
    VegOutFieldsMonthly  = .false.
    VegOutFieldsDaily  = .false.
    VegOutFields6hourly  = .false.
    VegOutFieldsHourly = .false.
    VegOutNCFYearly  = .false.
    VegOutNCFMonthly  = .false.
    VegOutNCFDaily  = .false.
    VegOutNCF6hourly  = .false.
    VegOutNCFHourly = .false.
#endif

end subroutine InitOutputConfig


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading chracteristics of the underlying surface
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GeoPhysData

    integer i, j, k, FileStat, Xscal, lenStr, ios, n, cnt, mon
    integer Rcp, RcpNum, rcpInd, Nrcp(MaxMatr)
    integer Surf, Surf1, SurfNum
    real AreaSurf(MaxSurf), sumLU, Aver(Imin:Imax), Frac
     character(20) codeGrid, IDrcp(MaxMatr), NameRcp(MaxMatr), RegRcp(MaxMatr)
     character(30) TypeSurf(MaxSurf)
    character(1) ch1, ch2
    character(800) readStr, readPar


#ifdef DEBUG_MODE
    print *, '>- Entering GeoPhysData...'
#endif
    write(YearNum,'(i4)') BegDate(yr)
    if(BegDate(yr)<2000) then
      write(YearShort,'(i2.2)') BegDate(yr)-1900
    else
      write(YearShort,'(i2.2)') BegDate(yr)-2000
    endif
    write(DayNum,'(i2.2)') BegDate(da)
!----------------------------------------------------------------------------
! Reading the surface altitude
    fullName=trim(GeoPath)//trim(OrogName)//trim(GridCode)//'.dat'
    open(5, file=trim(fullName), status='old',iostat=FileStat, action='read')
    if(FileStat>0) then
            print '(/,"STOP: Cannot open surface altitude file ''",a,"''",/)', trim(fullName)
      stop
    endif

    read(5,*) codeGrid
    lenStr=scan(codeGrid,achar(13))
    if(lenStr>0) codeGrid=codeGrid(:lenStr-1)
    if(trim(codeGrid)/=trim(GridCode)) then
      print '(/,"STOP: Incorrect grid code ''",a,"'' in file ",a,/)', trim(codeGrid), trim(fullName)
      stop
    endif

    read(5,'(/)')
    do
      read(5,'(a)',iostat=ios) readStr
      if(ios==-1) exit
      lenStr=scan(readStr,achar(13))
      if(lenStr>0) readStr=readStr(:lenStr-1)
      read(readStr,*) i, j, Zs(i,j)
    enddo
    close(5)

! Grid aggregation
    do j=Jmin, Jmax 
      if(maxI(j)==1) cycle
      Xscal=Imax/maxI(j)
      if(Xscal==1) cycle

      Aver(Imin:Imax)=Zs(Imin:Imax,j)
      call GridAggreg(j,Xscal,Aver,1)
      Zs(minI(j):maxI(j),j)=Aver(minI(j):maxI(j))
    enddo

! Reading season indicators
    fullName=trim(GeoPath)//trim(SeasonName)//trim(GridCode)//'.dat'
    open(11, file=trim(fullName), status='old', iostat=FileStat, action='read')
      if(FileStat>0) then
      print '(/,"STOP: Cannot open file with season indicators ''",a,"''",/)', trim(fullName)
      stop
    endif

    read(11,*) codeGrid
    lenStr=scan(codeGrid,achar(13))
    if(lenStr>0) codeGrid=codeGrid(:lenStr-1)
    if(trim(codeGrid)/=trim(GridCode)) then
      print '(/,"STOP: Incorrect grid code ''",a,"'' in file ",a,/)', trim(codeGrid), trim(fullName)
      stop
    endif

! Reading file header
    do n=1, 9
      read(11,*)
    enddo

    do
      read(11,'(a)',iostat=ios) readStr
      if(ios==-1) exit
      lenStr=scan(readStr,achar(13))
      if(lenStr>0) readStr=readStr(:lenStr-1)
      read(readStr,*) i, j, (Season(i,j,mon), mon=1, 12)
    enddo
    close(11)

! Grid aggregation
    do j=Jmin, Jmax
      Xscal=Imax/maxI(j)
      if(Xscal==1) cycle

      do mon=1, 12
        Aver(Imin:Imax)=Season(Imin:Imax,j,mon)
        call GridAggreg(j,Xscal,Aver,1)
        Season(minI(j):maxI(j),j,mon)=anint(Aver(minI(j):maxI(j)))
      enddo
    enddo

!***************** Matrix calculations *****************
#if RTYPE==2
! Reading receptors relative area
    fullName=trim(RecepPath)//trim(RecepArName)//trim(GridCode)//'.dat'
    open(5, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open receptors file ''",a,"''",/)', trim(fullName)
      stop
    endif

    read(5,*) codeGrid
    lenStr=scan(codeGrid,achar(13))
    if(lenStr>0) codeGrid=codeGrid(:lenStr-1)
    if(trim(codeGrid)/=trim(GridCode)) then
      print '(/,"STOP: Incorrect grid code ''",a,"'' in file ",a,/)', trim(codeGrid), trim(fullName)
      stop
    endif
    
    RecepPart=0.
    read(5,'(a)',iostat=ios) readStr
    lenStr=scan(readStr,achar(13))
    if(lenStr>0) readStr=readStr(:lenStr-1)
    read(readStr,*) RcpNum

    do Rcp=1, RcpNum
      read(5,'(a)',iostat=ios) readStr
      if(ios==-1) exit
      lenStr=scan(readStr,achar(13))
      if(lenStr>0) readStr=readStr(:lenStr-1)

      do n=1, 4
        lenStr=scan(readStr,',')
        if(lenStr>0) then
          readPar=readStr(:lenStr-1)
          readStr=readStr(lenStr+1:)
        else
          readPar=readStr
        endif
        select case(n)
        case(1)
          read(readPar,'(a)') NameRcp(Rcp)
        case(2)
          read(readPar,'(a)') IDrcp(Rcp)
        case(3)
          read(readPar,'(i)') Nrcp(Rcp)
        case(4)
          read(readPar,'(a)') RegRcp(Rcp)
        end select
      enddo

! Checking receptor with the receptors list
      rcpInd=0
      do cnt=1, MaxRcp
        if(IDrcp(Rcp)==RecepIDmax(cnt)) then
          rcpInd=cnt
          exit
        endif
      enddo
      if(rcpInd<=0) then
        print '(/,"STOP: Unknown receptor ''",a,"'' in file "a,/)', trim(IDrcp(Rcp)), trim(fullName)
        stop
      endif

! Reading receptors area fraction
      rcpInd=0
      do cnt=1, NumRcp
        if(IDrcp(Rcp)==RecepID(cnt)) then
            rcpInd=cnt
          exit
        endif
      enddo
      if(rcpInd==0) then
        do cnt=1, Nrcp(Rcp)
            read(5,'(a)',iostat=ios) readStr
            if(ios==-1) exit
            lenStr=scan(readStr,achar(13))
            if(lenStr>0) readStr=readStr(:lenStr-1)
            read(readStr,*) i, j, Frac
        enddo
      else
         do cnt=1, Nrcp(Rcp)
        read(5,'(a)',iostat=ios) readStr
        if(ios==-1) exit
            lenStr=scan(readStr,achar(13))
            if(lenStr>0) readStr=readStr(:lenStr-1)
            read(readStr,*) i, j, Frac
            RecepPart(i,j,rcpInd)=RecepPart(i,j,rcpInd)+Frac/100.  ! Percents -> fractions
        enddo
      endif
    enddo
    close(5)

! Checking availability of information for all receptors
    do Rcp=1, NumRcp
      rcpInd=0
      do cnt=1, RcpNum
        if(RecepID(Rcp)==IDrcp(cnt)) then
          rcpInd=cnt
          exit
        endif
      enddo
      if(rcpInd==0) then
        print '(/,"STOP: No information for receptor ''",a,"''",/)', trim(RecepID(Rcp))
        stop
      endif
    enddo

! Grid aggregation
    do j=Jmin, Jmax
      Xscal=Imax/maxI(j)
      if(Xscal==1) cycle

      do Rcp=1, NumRcp
        Aver(Imin:Imax)=RecepPart(Imin:Imax,j,Rcp)
        call GridAggreg(j,Xscal,Aver,1)
        RecepPart(minI(j):maxI(j),j,Rcp)=Aver(minI(j):maxI(j))
      enddo
    enddo
    
! Calculation of receptors areas
    RecepArea=0.
    do j=Jmin, Jmax
      do i=minI(j), maxI(j)
        do Rcp=1, NumRcp
          RecepArea(Rcp)=RecepArea(Rcp)+MeshArea(i,j)*RecepPart(i,j,Rcp)
        enddo
      enddo
    enddo

#endif
!***************** Matrix calculations *****************

#ifdef DEBUG_MODE
    print *, '<- Exit GeoPhysData'
#endif
end subroutine GeoPhysData


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading land cover data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadLandCover

    integer i, j, k, FileStat, Xscal, lenStr, ios, yrCur
    integer Surf, Surf1, SurfNum
    real AreaSurf(MaxSurf), sumLU, Aver(Imin:Imax)
    character(20) codeGrid
    character(30) TypeSurf(MaxSurf)
    character(1) ch1, ch2
    character(800) readStr, readPar


#ifdef DEBUG_MODE
    print *, '>- Entering ReadLandCover...'
#endif

    if(climLCRun) then
      yrCur=climLCYear
    else
      yrCur=Year
    endif
    write(YearNum,'(i4)') yrCur

!----------------------------------------------------------------------------
! Reading landcover data
    write(fileName, '(a,a,"_",a4,a4)') trim(LandName), trim(GridCode), YearNum, '.dat'
    fullName=trim(LandcoverPath)//trim(LandSource)//'/'//trim(GridCode)//'/'//trim(fileName)
    open(7, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
      stop
    endif
    read(7,*) ! skipping first row with LC product name
    read(7,*) codeGrid
    lenStr=scan(codeGrid,achar(13))
    if(lenStr>0) codeGrid=codeGrid(:lenStr-1)
    if(trim(codeGrid)/=trim(GridCode)) then
      print '(/,"STOP: Incorrect grid code ''",a,"'' in file ",a,/)', trim(codeGrid), trim(fullName)
      stop
    endif

    read(7,'(i2)') SurfNum
    if(SurfNum/=NumSurf) then
      print '(/,"STOP: Incorrect number of LC types ''",i2,"''",/)', SurfNum
      stop
    endif

    read(7,'(a)') readStr
    lenStr=scan(readStr,achar(13))
    if(lenStr>0) readStr=readStr(:lenStr-1)
    read(readStr,*) ch1, ch2, (TypeSurf(Surf), Surf=1, NumSurf)

    do Surf1=1, NumSurf
      do Surf=1, NumSurf
        if(TypeSurf(Surf1)==SurfType(Surf)) then
          exit
        endif
      enddo
      if(Surf>NumSurf) then
        print '(/,"STOP: Unknown land cover type ''",a,"''",/)', trim(TypeSurf(Surf1))
        stop
      endif
    enddo

    do
      read(7,'(a)',iostat=ios) readStr
      if(ios==-1) exit
      lenStr=scan(readStr,achar(13))
      if(lenStr>0) readStr=readStr(:lenStr-1)
      read(readStr,*) i, j, (AreaSurf(Surf), Surf=1, NumSurf)

      do Surf=1, NumSurf
        LandCover(i,j,Surf)=AreaSurf(Surf)
      enddo
    enddo
    close(7)

! Grid aggregation
    do Surf=1, NumSurf
      do j=Jmin, Jmax 
        if(maxI(j)==1) cycle
        Xscal=Imax/maxI(j)
        if(Xscal==1) cycle
        Aver(Imin:Imax)=LandCover(Imin:Imax,j,Surf)
        call GridAggreg(j,Xscal,Aver,1)
        LandCover(minI(j):maxI(j),j,Surf)=Aver(minI(j):maxI(j))
      enddo
    enddo

! Normalizing the land use data
    do j=Jmin, Jmax
      do i=minI(j), maxI(j)
        sumLU=sum(LandCover(i,j,1:NumSurf))
        if (sumLU == 0) then
          print '(/,"STOP: Incorrect LC data at gridcell (",i3,",",i3,") in file ",a,/)', i, j, trim(fullName)
          stop
        endif
        do Surf=1, NumSurf
           LandCover(i,j,Surf)=LandCover(i,j,Surf)/sumLU
        enddo
      enddo
    enddo
    
#ifdef G_POP
#ifdef M_SOIL
    call Soil_Initial ! Initialize soil
#endif

#ifdef M_OCN
    call Ocn_Initial ! Initialize ocean
#endif

#ifdef M_VEG
    call Veg_Initial ! Initialize vegetation
#endif
#endif

#ifdef DEBUG_MODE
    print *, '<- Exit ReadLandCover'
#endif
end subroutine ReadLandCover


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating time limits of the calculation run in each month and year
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CalcLimits(tmPer,perBeg,perFin)

    integer tmPer, perBeg, perFin
    
    selectcase(tmPer)
    
!------------------------------------------------------------------
      case(yr)
        perBeg=BegDate(yr)
        perFin=FinDate(yr)
!------------------------------------------------------------------
      case(mn)
        if(Year==BegDate(yr)) then 
          perBeg=BegDate(mn)
        else
          perBeg=1
        endif
        if(Year==FinDate(yr)) then 
          perFin=FinDate(mn)
        else
          perFin=12
        endif
!------------------------------------------------------------------
      case(da)
        if(Year==BegDate(yr).and.Month==BegDate(mn)) then 
          perBeg=BegDate(da)
        else
          perBeg=1
        endif
        if(Year==FinDate(yr).and.Month==FinDate(mn)) then 
          perFin=FinDate(da)
        else
          perFin=MonthDays(Month)
        endif
!------------------------------------------------------------------
    endselect 

end subroutine CalcLimits

end module GeneralProcs
