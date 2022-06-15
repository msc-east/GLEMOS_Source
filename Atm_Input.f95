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
! Module of the atmospheric information input 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module Atm_Input

  use GeneralParams
  use Atm_Output
  use Atm_HorizAdvect
  use Atm_VertAdvect
  use Atm_Exchange
#ifdef G_HG
  use Atm_Hg_General
#endif  
#ifdef G_HM
  use Atm_HM_General
#endif  
#ifdef G_POP
  use Atm_POP_General 
#endif  
#ifdef G_AERO
  use Atm_AERO_General
#endif
#ifdef G_Tracer
  use Atm_Tracer_General
#endif  
  use netcdf
  use typeSizes

  implicit none

  character(800), private :: fileName, fullName
  character(4), private :: YearNum
  character(2), private :: DayNum, MonthNum, YearShort


contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading pollutant specific configuration files
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_ReadConfig

    integer Subs, Grp

    call Atm_Config

#if RTYPE==2
    do Grp=1, NumGroups
      call Matrix_config(Grp)
    enddo
#endif

! Defining a list of sources
    call SourceList

! Reading physical and chemical properties
      do Subs=1, NumSubs
        Grp=SubsGroup(Subs)
        selectcase(Grp)
#ifdef G_HG
          case(HG)
            call Atm_Hg_Props(Subs)
#endif
#ifdef G_HM
          case(HM)
            call Atm_HM_Props(Subs)
#endif
#ifdef G_POP
          case(POP)
            call Atm_POP_Props(Subs)
#endif
#ifdef G_AERO
            case(AERO)
            call Atm_AERO_Props(Subs)
#endif
#ifdef G_Tracer
          case(Tracer)
            call Atm_Tracer_Props(Subs)
#endif
        endselect
      enddo
      call LogFile('Reading physical and chemical properties',0)
      
end subroutine Atm_ReadConfig


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading model grid configuration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Config

    integer k, FileStat, lenType, ios
    character(300) strRead, strPar
    character(30) strType
    
#ifdef DEBUG_MODE
    print *, '>- Entering Atm_Config...'
#endif

    filename='Atm'//trim(MediaConfName)//'.dat'
    fullName=trim(ConfigPath)//trim(filename)
    open(1, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open Atm config file ''",a,"''",/)', trim(fullName)
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
        case('Upper boundary')
          read(strPar,*) Ptop
!------------------------------------------------------------------
        case('Sigma bound')
          read(strPar,*) (Slevel(k), k=0, Atm_Kmax)
!------------------------------------------------------------------
      endselect
    enddo
    close(1)
    
! Definition of vertical layers
    do k=1, Atm_Kmax
      Sigma(k)=(Slevel(k)+Slevel(k-1))/2.
      dS(k)=Slevel(k-1)-Slevel(k)  
    enddo
    dS(Atm_Kmax+1)=dS(Atm_Kmax)  
    
#ifdef DEBUG_MODE
    print *, '<- Exit Atm_Config'
#endif

end subroutine Atm_Config


#if RTYPE==2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading configuration of emission sources
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Matrix_config(Grp)

    integer FileStat, lenType, ios, Src, Grp
    character(600) strRead, strPar
    character(30) strType

    filename=trim(SubsGroupID(SubsGrpLst(Grp)))//trim(MatrixConfName)//trim(GridCode)//'.dat'    
    fullName=trim(ConfigPath)//trim(filename)
    open(1, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open matrix config file ''",a,"''",/)', trim(fullName)
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
        case('Ant sources number')
          read(strPar,'(i2)') MaxAnth
!------------------------------------------------------------------
        case('Ant sources codes')
          read(strPar,*) (AnthrIDmax(Src), Src=1, MaxAnth)
!------------------------------------------------------------------
        case('Nat sources number')
          read(strPar,'(i2)') MaxNat
!------------------------------------------------------------------
        case('Nat sources codes')
          read(strPar,*) (NaturIDmax(Src), Src=1, MaxNat)
!------------------------------------------------------------------
        case('Bnd sources number')
          read(strPar,'(i2)') MaxBnd
!------------------------------------------------------------------
        case('Bnd sources codes')
          read(strPar,*) (BoundIDmax(Src), Src=1, MaxBnd)
!------------------------------------------------------------------
        case('Receptors number')
          read(strPar,'(i2)') MaxRcp
!------------------------------------------------------------------
        case('Receptors codes')
          read(strPar,*) (RecepIDmax(Src), Src=1, MaxRcp)
!------------------------------------------------------------------
      endselect
    enddo
    close(1)

end subroutine Matrix_config
#endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine defining a list of sources (matrix run)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SourceList

    integer Src, Src1, SrcInd

#if RTYPE==2
    selectcase(antSRCmode)
    case(1)
      NumAnth=maxAnth
      SourcID=AnthrIDmax
    case(2)
      SourcID=AnthrID
      NumAnth=NumAnth+1
      restAnt=NumAnth
      SourcID(restAnt)='restAnt'
    case(3)
      NumAnth=0
    endselect

    selectcase(natSRCmode)
    case(1)
      NumNat=maxNat
      NaturID=NaturIDmax
    case(2)
      NumNat=NumNat+1
      NaturID(NumNat)='restNat'
      restNat=NumNat
    case(3)
      NumNat=0
    endselect
    do Src=1, NumNat
      SourcID(NumAnth+Src)=NaturID(Src)
    enddo
    NumSrc=NumAnth+NumNat

#if (REGTYPE==2)
    selectcase(bndSRCmode)
    case(1)
      NumBnd=maxBnd
      BoundID=BoundIDmax
    case(2)
      NumBnd=NumBnd+1
      BoundID(NumBnd)='restBnd'
    case(3)
      NumBnd=0
    endselect

    SrcInd=0
    do Src=1, NumBnd
      SrcInd=0
      do Src1=1, NumSrc
        if(BoundID(Src)==SourcID(Src1)) then
          SrcInd=Src
          exit
        endif
      enddo
      if(SrcInd==0) then
        NumSrc=NumSrc+1
        SourcID(NumSrc)=BoundID(Src)
        if(trim(BoundID(Src))=='restBnd') restBnd=NumSrc 
      endif
    enddo
#endif

    selectcase(initSRCmode)
    case(1)
      NumSrc=NumSrc+1
      Init=NumSrc
      SourcID(Init)='Init'
    endselect

    selectcase(reemisSRCmode)
    case(1)
      NumSrc=NumSrc+1
      Reem=NumSrc
      SourcID(Reem)='Reemis'
    endselect

    selectcase(rcpSRCmode)
    case(1)
      NumRcp=maxRcp
      RecepID=RecepIDmax
    case(3)
      NumRcp=0
    endselect
#else
    NumAnth=1
    NumNat=1
    NumBnd=1
    NumSrc=1
    Init=1
    restAnt=1
    restNat=1
    restBnd=1
    NumRcp=1
#endif
    
end subroutine SourceList


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine supplying information for the calculation run
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_InputData(Prd)

    integer i, j, k, Grp
    integer Form       
    character(*) Prd
character(100) Note

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_InputData...', Prd
#endif

! Reading common input data
    selectcase(Prd)
!---------------------------------------------------------------------------------
    case('Initial')

      do k=1, Atm_Kmax+1
        do j=bJmin, bJmax
          do i=bminI(j), bmaxI(j)
           Veff(i,j,k)=MeshArea(i,j)*dS(k)/Ggrav
         enddo
        enddo
      enddo

! Calculation of parameters for Bott polinomials 
      call BottPolynomLat
       call BottPolynomSigma

      call Atm_InitMeteo

! Initializing the model variables
      selectcase(InitCond)
        case('zero')
          timeCalc=0.
        case('cond')
          timeCalc=0.
          call Atm_InitMass
          call MixRatioToConc
        case('dump')
          Period=1
      endselect

! Reading monitoring sites information
      call MonitorSites
      call LogFile('Reading monitoring sites information',0)
!---------------------------------------------------------------------------------
    case('Yearly')

! Checking for leap year
      if(mod(Year,4)==0) then
        MonthDays(2)=29
      else
        MonthDays(2)=28
      endif

! Reading emissions data
      call Atm_ReadEmis(yr)

! Reading boundary conditions
#if (REGTYPE==2)
      call ReadBoundary(yr)
#endif
!---------------------------------------------------------------------------------
    case('Monthly')

! Reading emissions data
      call Atm_ReadEmis(mn)

! Reading boundary conditions
#if (REGTYPE==2)
      call ReadBoundary(mn)
#endif

      call LogFile('Monthly')
!---------------------------------------------------------------------------------
    case('Daily')

! Reading emissions data
      call Atm_ReadEmis(da)

! Reading boundary conditions
#if (REGTYPE==2)
      call ReadBoundary(da)
#endif

! Reading meteorological data
      call Atm_ReadMeteo(Year,Month,Day)

! Calculation of atmospheric exchange parameters
      call Atm_ExchPar
      call LiqWatCont
!---------------------------------------------------------------------------------
    endselect

! Reading pollutant groups specific data
    do Grp=1, NumGroups
      selectcase(SubsGrpLst(Grp))
#ifdef G_HG
        case(HG)
          call Atm_Hg_Input(Prd)
#endif
#ifdef G_HM
        case(HM)
          call Atm_HM_Input(Prd)
#endif
#ifdef G_POP
        case(POP)
          call Atm_POP_Input(Prd)
#endif
#ifdef G_AERO
        case(AERO)
          call Atm_AERO_Input(Prd)
#endif
#ifdef G_TRACER
        case(TRACER)
          call Atm_Tracer_Input(Prd)
#endif
        endselect
      enddo
      
#ifdef DEBUG_MODE
    print *, '<- Exit Atm_InputData ', Prd
#endif
end subroutine Atm_InputData


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading emissions data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_ReadEmis(tm)

    integer tm, Ind, Form, Subs, Grp
    character(100) Note

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_ReadEmis... ', tm
#endif

! Reading anthropogenic emission data
      do Ind=1, AntEmisNum
        Form=ReadAntInd(Ind,tm)
        if(Form==0) cycle
        Subs=FormSubs(Atm,Form)

         call Atm_ReadAntEmis(Subs,Ind,tm,Note)

        call LogFile('Reading anthropogenic emissions data',0,Note)
      enddo
      
! Reading natural emission data (if any)
      do Ind=1, NatEmisNum
        Form=ReadNatInd(Ind,tm)
        if(Form==0) cycle
        Subs=FormSubs(Atm,Form)
        Grp=SubsGroup(Subs)
        selectcase(Grp)
#ifdef G_HG
          case(HG)
            call Atm_Hg_ReadNatEmis(Subs,Ind,tm,Note)
#endif    
#ifdef G_HM
          case(HM)
            call Atm_HM_ReadNatEmis(Subs,Ind,tm,Note)
#endif    
#ifdef G_TRACER
          case(TRACER)
            call Atm_Tracer_ReadNatEmis(Subs,Ind,tm,Note)
#endif    
        endselect
        call LogFile('Reading natural emission data',0,Note)
      enddo
      
end subroutine Atm_ReadEmis


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading anthropogenic emission
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_ReadAntEmis(Subs,Ind,tm,Note)

    integer i, j, ln, Ind, FileStat, Hs, Xscal, Form, Subs, lenStr, tm, ios, Nread
    integer Src, Src2, SrcInd
    character(20) IDsrc, IDsubs, IDform1, IDform2, TResol, prd, IDread(0:MaxMatr)
    character(300) readStr, parVal, SrcInfo
    character(30) parType, gridID, units
    real(8) EmisSrc(Atm_KMax), ConvEmisUnit
    real Aver(Imin:Imax)
    character(*) Note

!    Form=AntEmisInd(Ind)
    Form=ReadAntInd(Ind,tm)
    IDform1=FormID(Atm,Subs,Form)
    
    selectcase(tm)
      case(yr)
        write(fileName,'(i4)') Year
        fileName=trim(SubsID(Subs))//trim(IDform1)//'_'//trim(fileName)//trim(AntEmisName)//'.dat'
        ConvEmisUnit=1000./sum(MonthDays)/SecInDay      ! Convert t/y ==> kg/s
      case(mn)
        write(fileName,'(i4,i2.2)') Year, Month
        fileName=trim(SubsID(Subs))//trim(IDform1)//'_'//trim(fileName)//trim(AntEmisName)//'.dat'
        ConvEmisUnit=1000./MonthDays(Month)/SecInDay    ! Convert t/m ==> kg/s
      case(da)
        write(fileName,'(i4,i2.2,i2.2)') Year, Month, Day
        fileName=trim(SubsID(Subs))//trim(IDform1)//'_'//trim(fileName)//trim(AntEmisName)//'.dat'
        ConvEmisUnit=1000./SecInDay                     ! Convert t/day ==> kg/s
    endselect

    write(YearNum,'(i4)') Year
    fullName=trim(EmisPath)//trim(GridCode)//'/'//trim(SubsID(Subs))//trim(AntEmisName)//'/'//YearNum//'/'//fileName
    open(4, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open anthropogenic emission file ''",a,"''",/)', trim(fullName)
      stop
    endif

! Reading header
    do ln=1, 8
      read(4,'(a)') readStr
      lenStr=scan(readStr,achar(13))
      if(lenStr>0) readStr=readStr(:lenStr-1)

      lenStr=scan(readStr,':')
      if(lenStr==0) then
        print '(/,"STOP: Wrong format of file ''",a,"''",/)', trim(fullName)
        stop
      endif
      parType=trim(readStr(:lenStr-1))
      parVal=readStr(lenStr+1:)
      parVal=adjustl(parVal)

      selectcase(parType)
        case('Substance')
          IDsubs=parVal
        case('Form')
          IDform2=parVal
        case('Grid Type')
          gridID=parVal
        case('Source')
          SrcInfo=parVal
        case('TimeResol')
          TResol = ParVal
          if (tm == yr .and. trim(TResol) /= 'yearly') then
            print *, 'STOP: Invalid time resolution of emission data:', trim(TResol)
            stop
          end if
          if (tm == mn .and. trim(TResol) /= 'monthly') then
            print *, 'STOP: Invalid time resolution of emission data:', trim(TResol)
            stop
          end if
          if (tm == da .and. trim(TResol) /= 'daily') then
            print *, 'STOP: Invalid time resolution of emission data:', trim(TResol)
            stop
          end if
        case('Period')
          prd = parVal
        case('Units')
          units=parVal
        case('Layers')
          read(parVal,*) HemisAnt
        case default
          print '(/,"STOP: Unknown input parameter ''",a,"''",/)', trim(parType)
          stop
      endselect
    enddo

    write(Note,'("Emissions of ",a,a," in ",a,", ",a)') trim(IDsubs), trim(IDform2), prd, trim(units)

    read(4,*)
    read(4,*)
    read(4,*)

    AntEmisFlux(:,:,Ind,:,:)=0.
    IDread=''
    Nread=0
    do
      read(4,'(a)',iostat=ios) readStr
      if(ios==-1) exit
      lenStr=scan(readStr,achar(13))
      if(lenStr>0) readStr=readStr(:lenStr-1)
      read(readStr,*) IDsrc, i, j, (EmisSrc(Hs), Hs=1, HemisAnt)

!***************** Matrix calculations *****************
#if RTYPE==2
      do Src=1, Nread
        if(IDread(Src)==IDsrc) exit
      enddo
      if(Src>Nread) then
        Nread=Nread+1
        IDread(Nread)=IDsrc
      endif

! Checking the source with sources list
      SrcInd=0
      do Src=1, MaxAnth
        if(IDsrc==AnthrIDmax(Src)) then
          SrcInd=Src
          exit
        endif
      enddo
      if(SrcInd<=0) then
        print '(/,"STOP: Unknown anthropogenic source ''",a,"'' in file "a,/)', trim(IDsrc), trim(fullName)
        stop
      endif

! Identifying the source
      SrcInd=0
      do Src=1, NumAnth
        if(IDsrc==SourcID(Src)) then
          SrcInd=Src
          exit
        endif
      enddo
      if(SrcInd>0) then
        AntEmisFlux(i,j,Ind,SrcInd,1:HemisAnt)=AntEmisFlux(i,j,Ind,SrcInd,1:HemisAnt)+EmisSrc(1:HemisAnt)    ! t/y, t/m, t/day
      elseif(antSRCmode==2) then
        AntEmisFlux(i,j,Ind,restAnt,1:HemisAnt)=AntEmisFlux(i,j,Ind,restAnt,1:HemisAnt)+EmisSrc(1:HemisAnt)  ! t/y, t/m, t/day
      endif
#else
!*******************************************************
      AntEmisFlux(i,j,Ind,1,1:HemisAnt)=AntEmisFlux(i,j,Ind,1,1:HemisAnt)+EmisSrc(1:HemisAnt)   ! t/y, t/m, t/day
#endif
    enddo
    close(4)
    AntEmisFlux(:,:,Ind,:,:)=AntEmisFlux(:,:,Ind,:,:)*ConvEmisUnit                ! -> kg/sec

!***************** Matrix calculations *****************
#if RTYPE==2
! Checking completeness of emission data
    do Src=1, NumAnth
      if(SourcID(Src)=='restAnt') cycle
      SrcInd=0
      do Src2=1, Nread
        if(trim(IDread(Src2))==trim(SourcID(Src))) then
          SrcInd=Src2
          exit
        endif
      enddo
      if(SrcInd==0) then
        print '(/,"STOP: No  anthropogenic emission data for source ''",a,"''",/)', trim(SourcID(Src))
        stop
      endif
    enddo
#endif
!*******************************************************

! Grid aggregation
    do j=Jmin, Jmax
      if(maxI(j)==1) cycle
      Xscal=Imax/maxI(j)
      if(Xscal==1) cycle

      do Src=1, NumAnth
        do Hs=1, HemisAnt
          Aver(Imin:Imax)=AntEmisFlux(Imin:Imax,j,Ind,Src,Hs)
          call GridAggreg(j,Xscal,Aver,2)
          AntEmisFlux(minI(j):maxI(j),j,Ind,Src,Hs)=Aver(minI(j):maxI(j))
        enddo
      enddo
    enddo

end subroutine Atm_ReadAntEmis


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading meteo information in NetCDF format
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_ReadMeteo(Year,Month,Day)

    integer i, j, k, t, ncid_in, var_id_in, jV, Xscal, nDay, flag, st
    integer Year, Month, Day, Year2, Month2, Day2, yrCur
    integer start2(3), start3(4), count2(3), count3(4)
    real fld_3d_u(Imin-1:Imax,bJmin:bJmax,Atm_Kmax)
    real fld_3d_v(Imin:Imax,Jmin-1:Jmax,Atm_Kmax)
#if A_VTYPE==2
    real fld_3d_w(Imin:Imax,Jmin:Jmax,0:Atm_Kmax)
#endif
    real fld_3d(Imin:Imax,Jmin:Jmax,Atm_Kmax)
    real fld_3d_e(bImin:bImax,bJmin:bJmax,Atm_Kmax)
    real fld_2d(Imin:Imax,Jmin:Jmax)
    real fld_2d_e(bImin:bImax,bJmin:bJmax)
    real Aver(Imin:Imax), Tair, Mv, BLH, zDn, zUp
    real splRow(NumPer*2), dFdX(NumPer*2), d2FdX(NumPer*2), dFdXleft

! Check for climatic run
    if(climRun) then
      yrCur=ClimYear
    else
      yrCur=Year
    endif

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_ReadMeteo...', yrCur, Month, Day
#endif

! Open input netCDF file
    write(YearNum,'(i4)') yrCur
    write(fileName, '(i4,i2.2,i2.2,a3)') yrCur, Month, Day, '.nc'
    fullName=trim(MeteoPath)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
    st=nf90_open(trim(fullName),nf90_nowrite,ncid_in)
    if(st/=nf90_noerr) then                                                             ! Added AG 19.01.17
      print '(/,"STOP in Atm_ReadMeteo: Cannot open file ''",a,"''",/)', trim(fullName)   ! Added AG 19.01.17
      stop                                                                            ! Added AG 19.01.17
    endif                                                                               ! Added AG 19.01.17
    start2=(/1, 1, 1 /)
    start3=(/ 1, 1, 1, 1 /)
    count2=(/Imax,Jmax,1/)

    do t=1, NumPer
      start2(3)=t
      start3(4)=t
! Reading Uwind
      count3=(/Imax+1,bJmax-bJmin+1,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'U',var_id_in),2)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_3d_u,start3,count3),3)
      Uwind(Imin-1:Imax,bJmin:bJmax,1:Atm_Kmax,t,toDay)=fld_3d_u(Imin-1:Imax,bJmin:bJmax,1:Atm_Kmax)
! Reading Vwind
      count3=(/Imax,Jmax+1,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'V',var_id_in),6)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_3d_v,start3,count3),7)
      Vwind(Imin:Imax,bJmin:bJmax-1,1:Atm_Kmax,t,toDay)=fld_3d_v(Imin:Imax,bJmin:bJmax-1,1:Atm_Kmax)
#if A_VTYPE==2
! Reading Swind
      count3=(/Imax,Jmax,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'W',var_id_in),22)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_3d_w,start3,count3),23)
      Swind(Imin:Imax,Jmin:Jmax,1:Atm_Kmax,t,toDay)=fld_3d_w(Imin:Imax,Jmin:Jmax,1:Atm_Kmax)
#endif
! Reading Px
      count2=(/bImax-bImin+1,bJmax-bJmin+1,1/)
      call checkNC(nf90_inq_varid(ncid_in,'PSFC',var_id_in),22)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_2d_e,start2,count2),23)
      Px(bImin:bImax,bJmin:bJmax,t,toDay)=fld_2d_e(bImin:bImax,bJmin:bJmax)
! Reading TempAir
      count3 = (/bImax-bImin+1,bJmax-bJmin+1,Atm_Kmax,1/)
      call checkNC( nf90_inq_varid(ncid_in,'T',var_id_in),10)
      call checkNC( nf90_get_var(ncid_in,var_id_in,fld_3d_e, start3,count3),11)
      TempAir(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,t,toDay)=fld_3d_e(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)
! Reading HumidAir
      count3=(/bImax-bImin+1,bJmax-bJmin+1,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'QVAPOR',var_id_in),12)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_3d_e,start3,count3),13)
      HumidAir(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,t,toDay)=fld_3d_e(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)
! Reading PrecRainConv and PrecRainStrat
      count3=(/Imax,Jmax,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'RAINC3D',var_id_in),14)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_3d, start3, count3),15)
      PrecRainConv(Imin:Imax,Jmin:Jmax,1:Atm_Kmax,t)=fld_3d(Imin:Imax,Jmin:Jmax,1:Atm_Kmax)
      call checkNC( nf90_inq_varid(ncid_in,'RAINNC3D',var_id_in),16)
      call checkNC( nf90_get_var(ncid_in,var_id_in,fld_3d, start3,count3),17)
      PrecRainStrat(Imin:Imax,Jmin:Jmax,1:Atm_Kmax,t)=fld_3d(Imin:Imax,Jmin:Jmax,1:Atm_Kmax)
! Reading Ksigma
      count3=(/Imax,Jmax,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'EXCH_H',var_id_in),18)
      call checkNC(nf90_get_var(ncid_in, var_id_in,fld_3d,start3,count3),19)
      Ksigma(Imin:Imax,Jmin:Jmax,1:Atm_Kmax-1,t,toDay)=fld_3d(Imin:Imax,Jmin:Jmax,2:Atm_Kmax)
! Reading WaterCont
      count3=(/Imax,Jmax,Atm_Kmax,1/)
      call checkNC( nf90_inq_varid(ncid_in,'QCLOUD',var_id_in),20)
      call checkNC( nf90_get_var(ncid_in,var_id_in,fld_3d,start3,count3),21)
      WaterCont(Imin:Imax,Jmin:Jmax,1:Atm_Kmax,t)=fld_3d(Imin:Imax,Jmin:Jmax,1:Atm_Kmax)
! Reading TempSurf
      count2=(/Imax,Jmax,1/)
      call checkNC( nf90_inq_varid(ncid_in,'TSK',var_id_in),27)
      call checkNC(nf90_get_var(ncid_in, var_id_in,fld_2d,start2,count2),28)
      TempSurf(Imin:Imax,Jmin:Jmax,t)=fld_2d(Imin:Imax,Jmin:Jmax)
! Reading SnowDepth
      call checkNC(nf90_inq_varid(ncid_in,'SNOWH',var_id_in),29)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_2d,start2,count2),30)
      SnowDepth(Imin:Imax,Jmin:Jmax,t)=fld_2d(Imin:Imax,Jmin:Jmax)
! Reading MOLength
      call checkNC( nf90_inq_varid(ncid_in, 'MOL', var_id_in), 31 )
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_2d,start2,count2),32)
      MOLength(Imin:Imax,Jmin:Jmax,t,toDay) = fld_2d(Imin:Imax,Jmin:Jmax)
! Reading Ufric
      call checkNC(nf90_inq_varid(ncid_in,'UST',var_id_in),33)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_2d,start2,count2),34)
      Ufric(Imin:Imax,Jmin:Jmax,t,toDay)=fld_2d(Imin:Imax,Jmin:Jmax)
! Reading SeaIce
      call checkNC(nf90_inq_varid(ncid_in,'SEAICE',var_id_in),35)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_2d,start2,count2),36)
      SeaIce(Imin:Imax,Jmin:Jmax,t)=fld_2d(Imin:Imax,Jmin:Jmax)
! Reading SolarRad
      call checkNC(nf90_inq_varid(ncid_in,'SWDOWN',var_id_in),37)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_2d,start2,count2),38)
      SolarRad(Imin:Imax,Jmin:Jmax,t)=fld_2d(Imin:Imax,Jmin:Jmax)
! Reading BLH
      call checkNC( nf90_inq_varid(ncid_in, 'PBLH', var_id_in), 39 )
      call checkNC( nf90_get_var(ncid_in, var_id_in, fld_2d, start2, count2), 40 )
      BLHeight(Imin:Imax,Jmin:Jmax,t)=fld_2d(Imin:Imax,Jmin:Jmax)
    enddo !t
    call checkNC(nf90_close(ncid_in),110)

! Removing possible negative values
    where(PrecRainConv<Zero) PrecRainConv=0.
    where(PrecRainStrat<Zero) PrecRainStrat=0.

    if(Day<MonthDays(Month)) then
      Day2=Day+1
      Month2=Month
      Year2=yrCur
    elseif(Month<12) then
      Day2=1
      Month2=Month+1
      Year2=yrCur
    else
      Day2=1
      Month2=1
      if(yrCur<FinDate(yr)) then
        Year2=yrCur+1
      else
        Year2=yrCur
      endif
    endif

! Check for climatic run
    if(climRun) then
      Year2=ClimYear
    endif

    write(YearNum,'(i4)') Year2
    write(fileName, '(i4,i2.2,i2.2,a3)') Year2, Month2, Day2, '.nc'
    fullName=trim(MeteoPath)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
    st = nf90_open(trim(fullName),nf90_nowrite,ncid_in)
    if(st /= nf90_noerr) then                                                             ! Added AG 19.01.17
      print '(/,"STOP in Atm_ReadMeteo: Cannot open file ''",a,"''",/)', trim(fullName)   ! Added AG 19.01.17
      stop                                                                            ! Added AG 19.01.17
    endif                                                                               ! Added AG 19.01.17
    start2=(/1, 1, 1 /)
    start3=(/ 1, 1, 1, 1 /)
    count2=(/Imax,Jmax,1/)

    do t=1, NumPer
      start2(3)=t
      start3(4)=t
! Reading Uwind
      count3=(/Imax+1,bJmax-bJmin+1,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'U',var_id_in),2)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_3d_u,start3,count3),3)
      Uwind(Imin-1:Imax,bJmin:bJmax,1:Atm_Kmax,t,toMor)=fld_3d_u(Imin-1:Imax,bJmin:bJmax,1:Atm_Kmax)
! Reading Vwind
      count3=(/Imax,Jmax+1,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'V',var_id_in),6)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_3d_v,start3,count3),7)
      Vwind(Imin:Imax,bJmin:bJmax-1,1:Atm_Kmax,t,toMor)=fld_3d_v(Imin:Imax,bJmin:bJmax-1,1:Atm_Kmax)
#if A_VTYPE==2
! Reading Swind
      count3=(/Imax,Jmax,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'W',var_id_in),22)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_3d_w,start3,count3),23)
      Swind(Imin:Imax,Jmin:Jmax,1:Atm_Kmax,t,toMor)=fld_3d_w(Imin:Imax,Jmin:Jmax,1:Atm_Kmax)
#endif
! Reading Px
      count2=(/bImax-bImin+1,bJmax-bJmin+1,1/)
      call checkNC(nf90_inq_varid(ncid_in,'PSFC',var_id_in),22)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_2d_e,start2,count2),23)
      Px(bImin:bImax,bJmin:bJmax,t,toMor)=fld_2d_e(bImin:bImax,bJmin:bJmax)
! Reading TempAir
      count3 = (/bImax-bImin+1,bJmax-bJmin+1,Atm_Kmax,1/)
      call checkNC( nf90_inq_varid(ncid_in,'T',var_id_in),10)
      call checkNC( nf90_get_var(ncid_in,var_id_in,fld_3d_e, start3,count3),11)
      TempAir(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,t,toMor)=fld_3d_e(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)
! Reading HumidAir
      count3=(/bImax-bImin+1,bJmax-bJmin+1,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'QVAPOR',var_id_in),12)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_3d_e,start3,count3),13)
      HumidAir(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,t,toMor)=fld_3d_e(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)
! Reading Ksigma
      count3=(/Imax,Jmax,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'EXCH_H',var_id_in),18)
      call checkNC(nf90_get_var(ncid_in, var_id_in,fld_3d,start3,count3),19)
      Ksigma(Imin:Imax,Jmin:Jmax,1:Atm_Kmax-1,t,toMor)=fld_3d(Imin:Imax,Jmin:Jmax,2:Atm_Kmax)
! Reading MOLength
      count2=(/Imax,Jmax,1/)
      call checkNC( nf90_inq_varid(ncid_in, 'MOL', var_id_in), 31 )
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_2d,start2,count2),32)
      MOLength(Imin:Imax,Jmin:Jmax,t,toMor)=fld_2d(Imin:Imax,Jmin:Jmax)
! Reading Ufric
      call checkNC(nf90_inq_varid(ncid_in,'UST',var_id_in),33)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_2d,start2,count2),34)
      Ufric(Imin:Imax,Jmin:Jmax,t,toMor)=fld_2d(Imin:Imax,Jmin:Jmax)
    enddo !t
    call checkNC(nf90_close(ncid_in),110)

! Aggregation of meteo parameters
    do t=1, NumPer
      do k=1, Atm_Kmax

! Meridional wind velocoty
        do j=bJmin, bJmax-1
          if(maxI(j)<maxI(j+1)) then
            jV=j+1
          else
            jV=j
          endif
          Xscal=Imax/maxI(jV)
          if(Xscal==1) cycle

          do nDay=1, 2
            Aver(Imin:Imax)=Vwind(Imin:Imax,j,k,t,nDay)
            call GridAggreg(jV,Xscal,Aver,1)
            Vwind(minI(jV):maxI(jV),j,k,t,nDay)=Aver(minI(jV):maxI(jV))
          enddo
        enddo

        do j=bJmin, bJmax
          Xscal=Imax/maxI(j)
          if(Xscal==1) cycle

! Zonal velocity
          do i=minI(j), maxI(j)
            Uwind(i,j,k,t,1:2)=Uwind(i*Xscal,j,k,t,1:2)
          enddo

! Averaging 3D fields
          do nDay=1, 2
            Aver(Imin:Imax)=TempAir(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            TempAir(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
          enddo

          do nDay=1, 2
            Aver(Imin:Imax)=HumidAir(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            HumidAir(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
          enddo

          Aver(Imin:Imax)=PrecRainConv(Imin:Imax,j,k,t)
          call GridAggreg(j,Xscal,Aver,1)
          PrecRainConv(minI(j):maxI(j),j,k,t)=Aver(minI(j):maxI(j))

          Aver(Imin:Imax)=PrecRainStrat(Imin:Imax,j,k,t)
          call GridAggreg(j,Xscal,Aver,1)
          PrecRainStrat(minI(j):maxI(j),j,k,t)=Aver(minI(j):maxI(j))

          do nDay=1, 2
            Aver(Imin:Imax)=Ksigma(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            Ksigma(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
          enddo

          Aver(Imin:Imax)=WaterCont(Imin:Imax,j,k,t)
          call GridAggreg(j,Xscal,Aver,1)
          WaterCont(minI(j):maxI(j),j,k,t)=Aver(minI(j):maxI(j))

#if A_VTYPE==2
          do nDay=1, 2
            Aver(Imin:Imax)=Swind(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            Swind(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
          enddo
#endif
        enddo
      enddo
    enddo
if(yrCur==1995.and.Month==5.and.Day==14) Vwind(:,:,:,4,1)=(Vwind(:,:,:,3,1)+Vwind(:,:,:,1,2))/2

! Calculation of RT
    do nDay=1, 2
    do t=1, NumPer
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=minI(j), maxI(j)
            Tair=TempAir(i,j,k,t,nDay)
            Mv=HumidAir(i,j,k,t,nDay)
            RTemp(i,j,k,t,nDay)=(1.+HumCoef*Mv/(1.+Mv))*Rair*Tair
          enddo
        enddo
      enddo
    enddo
    enddo

! Transformatrion of BLH from z to sigma coordinates
    do t=1, NumPer
      PxCurr(:,:)=Px(:,:,t,toDay)
      do j=Jmin, Jmax
        do i=minI(j), maxI(j)
          BLH=BLHeight(i,j,t)
          do k=1, Atm_Kmax
            zDn=Sheight(i,j,k-1)-Zs(i,j)
            zUp=Sheight(i,j,k)-Zs(i,j)
            if(BLH>=zDn.and.BLH<zUp) then
              BLHeight(i,j,t)=Slevel(k-1)+(Slevel(k)-Slevel(k-1))/(zUp-zDn)*(BLH-zDn)
              exit
            endif
          enddo
        enddo
      enddo
    enddo

! Averaging 2D fields
    do t=1, NumPer
      do j=bJmin, bJmax
        Xscal=Imax/maxI(j)
        if(Xscal==1) cycle

        do nDay=1, 2
          Aver(Imin:Imax)=Px(Imin:Imax,j,t,nDay)
          call GridAggreg(j,Xscal,Aver,1)
          Px(minI(j):maxI(j),j,t,nDay)=Aver(minI(j):maxI(j))
        enddo

        Aver(Imin:Imax)=Roughness(Imin:Imax,j,t)
        call GridAggreg(j,Xscal,Aver,1)
        Roughness(minI(j):maxI(j),j,t)=Aver(minI(j):maxI(j))

        Aver(Imin:Imax)=TempSurf(Imin:Imax,j,t)
        call GridAggreg(j,Xscal,Aver,1)
        TempSurf(minI(j):maxI(j),j,t)=Aver(minI(j):maxI(j))

        Aver(Imin:Imax)=SoilWater(Imin:Imax,j,t)
        call GridAggreg(j,Xscal,Aver,1)
        SoilWater(minI(j):maxI(j),j,t)=Aver(minI(j):maxI(j))

        Aver(Imin:Imax)=SnowDepth(Imin:Imax,j,t)
        call GridAggreg(j,Xscal,Aver,1)
        SnowDepth(minI(j):maxI(j),j,t)=Aver(minI(j):maxI(j))

        do nDay=1, 2
          Aver(Imin:Imax)=MOLength(Imin:Imax,j,t,nDay)
          call GridAggreg(j,Xscal,Aver,1)
          MOLength(minI(j):maxI(j),j,t,nDay)=Aver(minI(j):maxI(j))
        enddo

        do nDay=1, 2
          Aver(Imin:Imax)=Ufric(Imin:Imax,j,t,nDay)
          call GridAggreg(j,Xscal,Aver,1)
          Ufric(minI(j):maxI(j),j,t,nDay)=Aver(minI(j):maxI(j))
        enddo

        Aver(Imin:Imax)=SeaIce(Imin:Imax,j,t)
        call GridAggreg(j,Xscal,Aver,1)
        SeaIce(minI(j):maxI(j),j,t)=Aver(minI(j):maxI(j))

        Aver(Imin:Imax)=SolarRad(Imin:Imax,j,t)
        call GridAggreg(j,Xscal,Aver,1)
        SolarRad(minI(j):maxI(j),j,t)=Aver(minI(j):maxI(j))

        Aver(Imin:Imax)=BLHeight(Imin:Imax,j,t)
        call GridAggreg(j,Xscal,Aver,1)
        BLHeight(minI(j):maxI(j),j,t)=Aver(minI(j):maxI(j))
      enddo
    enddo

! Calculation of temporal interpolation parameters
    do j=bJmin, bJmax
      do i=bminI(j), bmaxI(j)
        splRow(1:NumPer)=Px(i,j,1:NumPer,toDay)
        splRow(NumPer+1:NumPer*2)=Px(i,j,1:NumPer,toMor)
        if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
          flag=0
          dFdXleft=0.
        else
          flag=1
          dFdXleft=dPxdT(i,j,NumPer+1)
        endif
        call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
        dPxdT(i,j,1:NumPer+1)=dFdX(1:NumPer+1)
        d2PxdT(i,j,1:NumPer+1)=d2FdX(1:NumPer+1)
!----------------------------------------------------------------------------------
        do k=1, Atm_Kmax
          splRow(1:NumPer)=Uwind(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Uwind(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dUwind(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dUwind(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2Uwind(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!----------------------------------------------------------------------------------
          splRow(1:NumPer)=Vwind(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Vwind(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dVwind(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dVwind(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2Vwind(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!----------------------------------------------------------------------------------
#if A_VTYPE==2
          splRow(1:NumPer)=Swind(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Swind(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dTair(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dSwind(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2Swind(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
#endif
!----------------------------------------------------------------------------------
          splRow(1:NumPer)=TempAir(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=TempAir(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dTair(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dTair(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2Tair(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!----------------------------------------------------------------------------------
          splRow(1:NumPer)=Ksigma(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Ksigma(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dKsigma(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dKsigma(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2Ksigma(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
        enddo
      enddo
    enddo

#ifdef DEBUG_MODE
    print *, '>- Exit Atm_ReadMeteo...'
#endif

end subroutine Atm_ReadMeteo


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating liquid water content
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LiqWatCont

    integer i, j, k, t
    real RhoAir, Tair, Wcont
    real :: T0=273., Tfr=258.

    do t=1, NumPer
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=minI(j), maxI(j)

            Tair=TempAir(i,j,k,t,toDay)
            RhoAir=(Sigma(k)*Px(i,j,t,toDay)+Ptop)/Rair/Tair

            Wcont=WaterCont(i,j,k,t)*RhoAir

            LiqCont(i,j,k,t)=Wcont*max(min((Tair-Tfr)/(T0-Tfr),1.),0.)
            FrozCont(i,j,k,t)=Wcont-LiqCont(i,j,k,t)
          enddo
        enddo
      enddo
    enddo

end subroutine LiqWatCont


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating initial meteo parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_InitMeteo

    integer i, j, k, t, ncid_in, var_id_in, jV, Xscal, yrCur, st
    integer start2(3), start3(4), count2(3), count3(4)
    real fld_3d_e(bImin:bImax,bJmin:bJmax,Atm_Kmax)
    real fld_2d_e(bImin:bImax,bJmin:bJmax)
    real Aver(Imin:Imax), Tair, Mv

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_InitMeteo...'
#endif

! Check for climatic run
    if(climRun) then
      yrCur=ClimYear
    else
      yrCur=BegDate(yr)
    endif

! Open input netCDF file
    write(YearNum,'(i4)') yrCur
    write(fileName, '(i4,i2.2,i2.2,a3)') yrCur, BegDate(mn), BegDate(da), '.nc'
    fullName=trim(MeteoPath)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
    st=nf90_open(trim(fullName),nf90_nowrite,ncid_in)
    if(st/=nf90_noerr) then                                                             ! Added AG 19.01.17
      print '(/,"STOP in Atm_InitMeteo: Cannot open file ''",a,"''",/)', trim(fullName)   ! Added AG 19.01.17
      stop                                                                            ! Added AG 19.01.17
    endif                                                                               ! Added AG 19.01.17
    start2=(/1, 1, 1 /)
    start3=(/ 1, 1, 1, 1 /)

    do t=1, NumPer
      start2(3)=t
      start3(4)=t
! Px
      count2=(/bImax-bImin+1,bJmax-bJmin+1,1/)
      call checkNC(nf90_inq_varid(ncid_in,'PSFC',var_id_in),22)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_2d_e,start2,count2),23)
      Px(bImin:bImax,bJmin:bJmax,t,toDay)=fld_2d_e(bImin:bImax,bJmin:bJmax)
! TempAir
      count3 = (/bImax-bImin+1,bJmax-bJmin+1,Atm_Kmax,1/)
      call checkNC( nf90_inq_varid(ncid_in,'T',var_id_in),10)
      call checkNC( nf90_get_var(ncid_in,var_id_in,fld_3d_e, start3,count3),11)
      TempAir(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,t,toDay)=fld_3d_e(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)
! HumidAir
      count3=(/bImax-bImin+1,bJmax-bJmin+1,Atm_Kmax,1/)
      call checkNC(nf90_inq_varid(ncid_in,'QVAPOR',var_id_in),12)
      call checkNC(nf90_get_var(ncid_in,var_id_in,fld_3d_e,start3,count3),13)
      HumidAir(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,t,toDay)=fld_3d_e(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)
    end do !t
    call checkNC(nf90_close(ncid_in),110)

    do j=bJmin, bJmax
      Xscal=Imax/maxI(j)
      if(Xscal==1) cycle

      Aver(Imin:Imax)=Px(Imin:Imax,j,1,toDay)
      call GridAggreg(j,Xscal,Aver,1)
      Px(minI(j):maxI(j),j,1,toDay)=Aver(minI(j):maxI(j))

      do k=1, Atm_Kmax
        Aver(Imin:Imax)=TempAir(Imin:Imax,j,k,1,toDay)
        call GridAggreg(j,Xscal,Aver,1)
        TempAir(minI(j):maxI(j),j,k,1,toDay)=Aver(minI(j):maxI(j))

        Aver(Imin:Imax)=HumidAir(Imin:Imax,j,k,1,toDay)
        call GridAggreg(j,Xscal,Aver,1)
        HumidAir(minI(j):maxI(j),j,k,1,toDay)=Aver(minI(j):maxI(j))
      enddo
    enddo

    do j=bJmin, bJmax
      PxCurr(bminI(j):bmaxI(j),j)=Px(bminI(j):bmaxI(j),j,1,toDay)
    enddo

    Period=1
    TairCurr(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)=TempAir(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1,toDay)

    call RTcurrent

    call AirDensity

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_InitMeteo'
#endif
end subroutine Atm_InitMeteo


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading boundary conditions
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadBoundary(tm)

    integer tm, i, j, k, t, Ind, FileStat, Xscal, Subs, Form, BndNum, Src, Src1, SrcInd
    real readMix
    real Aver(Imin:Imax)
    character(8) IDbnd(MaxMatr)

    write(YearNum,'(i4)') Year

    Atm_BoundW=0.
    Atm_BoundE=0.
    Atm_BoundS=0.
    Atm_BoundN=0.
    
    do Ind=1, AtmBndNum
      Form=ReadBndInd(Ind,tm)
      if(Form==0) cycle
      Subs=FormSubs(Atm,Form)

      selectcase(tm)
      case(da)
        write(fileName,'(i4,i2.2,i2.2,a4)') Year, Month, Day, '.bin'
      case(mn)
        write(fileName,'(i4,i2.2,a4)') Year, Month, '.bin'
      case(yr)
        write(fileName,'(i4,a4)') Year, '.bin'
      endselect
      fileName=trim(SubsID(Subs))//trim(FormID(Atm,Subs,Form))//trim(BoundName)//trim(fileName)
      fullName=trim(InBoundPath)//trim(GridCode)//'/'//trim(SubsID(Subs))//'/'//YearNum//'/'//trim(fileName)
      open(5, file=trim(fullName), form='unformatted', access="stream", status='old', iostat=FileStat, action='read')
      if(FileStat>0) then
        print '(/,"STOP: Cannot open boundary file ''",a,"''",/)', trim(fullName)
        stop
      endif

      read(5) BndNum
      read(5) (IDbnd(Src), Src=1, BndNum)

!      where(IDbnd=='restAnt') IDbnd='bndAnt'
!      where(IDbnd=='restNat') IDbnd='bndNat'

      do Src=1, BndNum

!***************** Matrix calculations *****************
#if RTYPE==2
         
! Checking source with the sources list
        SrcInd=0
        do Src1=1, MaxBnd
          if(IDbnd(Src)==BoundIDmax(Src1)) then
            SrcInd=Src1
            exit
          endif
        enddo
        if(SrcInd<=0) then
          print '(/,"STOP: Unknown boundary source ''",a,"'' in file "a,/)', trim(IDbnd(Src)), trim(fullName)
          stop
        endif

! Identifying the source
        SrcInd=0
        do Src1=1, NumSrc
          if(IDbnd(Src)==SourcID(Src1)) then
            SrcInd=Src1
            exit
          endif
        enddo
        if(SrcInd==0.and.bndSRCmode==2) then
          SrcInd=restBnd
        endif
#else        
!*******************************************************
        SrcInd=1
#endif

        do t=1, NumPer
          do k=1, Atm_Kmax
            do j=bJmin, bJmax
              read(5) readMix
              Atm_BoundW(j,k,Form,t,SrcInd)=Atm_BoundW(j,k,Form,t,SrcInd)+real(readMix,8)
              read(5) readMix
              Atm_BoundE(j,k,Form,t,SrcInd)=Atm_BoundE(j,k,Form,t,SrcInd)+real(readMix,8)
            enddo
            do i=bImin, bImax
              read(5) readMix
              Atm_BoundS(i,k,Form,t,SrcInd)=Atm_BoundS(i,k,Form,t,SrcInd)+real(readMix,8)
              read(5) readMix
              Atm_BoundN(i,k,Form,t,SrcInd)=Atm_BoundN(i,k,Form,t,SrcInd)+real(readMix,8)
            enddo
          enddo     ! k
        enddo       ! t
      enddo         ! Src
      close(5)

! Grid aggregation

!***************** Matrix calculations *****************
#if RTYPE==2
      do Src=1, NumSrc
#else        
!*******************************************************
      do Src=1, 1
#endif
        do t=1, NumPer
          do k=1, Atm_Kmax
            Xscal=Imax/maxI(Jmin)
            if(Xscal>1) then
              Aver(Imin:Imax)=Atm_BoundS(Imin:Imax,k,Form,t,Src)
              call GridAggreg(Jmin,Xscal,Aver,1)
              Atm_BoundS(minI(Jmin):maxI(Jmin),k,Form,t,Src)=Aver(minI(Jmin):maxI(Jmin))
              Atm_BoundS(bmaxI(Jmin),k,Form,t,Src)=Atm_BoundS(bImax,k,Form,t,Src)
            endif
            Xscal=Imax/maxI(Jmax)
            if(Xscal>1) then
              Aver(Imin:Imax)=Atm_BoundN(Imin:Imax,k,Form,t,Src)
              call GridAggreg(Jmax,Xscal,Aver,1)
              Atm_BoundN(minI(Jmax):maxI(Jmax),k,Form,t,Src)=Aver(minI(Jmax):maxI(Jmax))
              Atm_BoundN(bmaxI(Jmax),k,Form,t,Src)=Atm_BoundN(bImax,k,Form,t,Src)
            endif
          enddo     ! k
        enddo       ! t
      enddo         ! Src
    enddo           ! Ind

end subroutine ReadBoundary


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading initial conditions from the dump file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadInitCond

    integer FormNum, SurfNum, AnthNum, NatNum, ReNum, SrcNum, RcpNum, FileStat, Src
    integer DumpYear, DumpMonth, DumpDay
    integer DumpFreq
    character(4) YearNum
    character(8) IDsourc(MaxMatr)
    character(10) DumpGrid, InputDumpVersion
    

    fullName=trim(InitPath)//trim(InitName)

    open(1, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open dump file ''",a,"''",/)', trim(fullName)
      stop
    endif

! Reading dump version
    read(1) InputDumpVersion
    if(trim(InputDumpVersion)/=trim(DumpVersion)) then
      print '(/,"WARNING: Incorrect dump version ''",a,"'' in file ",a,/)', trim(InputDumpVersion), trim(fullName)
      read(1, pos=1)
    endif
    
! Reading grid code
    read(1) DumpGrid

    if(trim(DumpGrid)/=trim(GridCode)) then
      print '(/,"STOP: Incorrect grid code ''",a,"'' in file ",a,/)', trim(DumpGrid), trim(fullName)
      stop
    endif

! Dumping date
    read(1) DumpYear, DumpMonth, DumpDay, timeCalc

! Run characteristics
    read(1) FormNum, SurfNum, AnthNum, NatNum, ReNum, SrcNum, RcpNum

    if(FormNum/=NumForm(Atm)) then
      print '(/,"STOP: Wrong number of the pollutant components&
        & in the dump file ''",a,"''",/)', trim(fileName)
      stop
    endif

! Content in media for initial conditions
        read(1) Atm_MixRatio(bImin:bImax,bJmin:bJmax,0:Atm_Kmax+1,1:NumForm(Atm))
        read(1) Atm_Conc(bImin:bImax,bJmin:bJmax,0:Atm_Kmax+1,1:NumForm(Atm))
#ifdef M_SOIL
        read(1) Soil_Conc(IMin:IMax,JMin:JMax,1:Soil_KMax,1:Soil_NumType,1:NumForm(Soil))
#endif
#ifdef M_OCN
        read(1) Ocn_conc(Imin:Imax,Jmin:Jmax,1:Ocn_Kmax,1:NumForm(Ocn),1:3)
        read(1) Water_upl_conc(IMin:IMax,JMin:JMax,1:NumForm(Ocn))
#endif
#ifdef M_VEG
        read(1) Veg_Conc(IMin:IMax,JMin:JMax,1:Veg_NumType,1:NumForm(Veg))
        read(1) Fall_Conc(IMin:IMax,JMin:JMax,1:Veg_NumType,1:NumSubsMedia(Veg))
#endif

!***************** Matrix calculations *****************
#if RTYPE==2
        selectcase(initSRCmode)
        case(1)
          Atm_Contrib=0.
          Atm_Contrib(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),Init)=1.

        case(2)
          if(SrcNum/=NumSrc) then
            print '(/,"STOP: Wrong number of sources in initial conditions ''",i3,"''",/)', SrcNum
            stop
          endif
          
          read(1) IDsourc
          
          do Src=1, NumSrc
            if(IDsourc(Src)/=SourcID(Src)) then
              print '(/,"STOP: Wrong list of sources in initial conditions",/)'
              stop
            endif
          enddo

          read(1) Atm_Contrib(Imin:Imax,bJmin:bJmax,0:Atm_Kmax+1,1:NumForm(Atm),1:SrcNum)

        endselect  
#endif
!***************** Matrix calculations *****************

        close(1)

end subroutine ReadInitCond


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating initial pollutant mass in the atmosphere
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_InitMass

    integer i, j, k
    real(8) qInit(MaxForm), mInit(MaxForm)

#ifdef DEBUG_MODE
    print *, '<- Enter Atm_InitMass'
#endif
    
    MassAtmInit = 0.

    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        do i=minI(j), maxI(j)
          qInit(1:NumForm(Atm))=Atm_MixRatio(i,j,k,1:NumForm(Atm))
          mInit(1:NumForm(Atm))=qInit(1:NumForm(Atm))*Veff(i,j,k)*PxCurr(i,j)
          MassAtmInit(1:NumForm(Atm))=MassAtmInit(1:NumForm(Atm))+mInit(1:NumForm(Atm))
        enddo
      enddo
    enddo

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_InitMass'
#endif

end subroutine Atm_InitMass


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading meteo information
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MonitorSites

    integer i, j, FileStat, st, lenStr, ios
    real longM, latM, long11, long12, long21, long22, lat11, lat12, lat21, lat22
    real Dist(2,2), Rdist, wt(2,2), wSum, Lat1, Lat2, dLon, cosR
    character(2) typeSt
    character(300) readStr

#ifdef DEBUG_MODE
    print *, '<- Enter MonitorSites'
#endif

    fileName=trim(SubsID(1))//trim(StatName)//'.dat'
    fullName=trim(StatPath)//fileName
    open(2, file=fullName, status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file with monitoring sites ''",a,"''",/)', trim(fullName)
      stop
    endif
    read(2,*)    
    read(2,*)    
    Nair=0
    Nprec=0
    do
      read(2,'(a)',iostat=ios) readStr
      if(ios==-1) exit
      lenStr=scan(readStr,achar(13))
      if(lenStr>0) readStr=readStr(:lenStr-1)
      read(readStr,*) ReadStat%nameSt, ReadStat%indSt, ReadStat%latSt, ReadStat%longSt, typeSt
      if(ReadStat%longSt<xOrig/pi180.or.ReadStat%longSt>(xOrig+real(Imax)*dXmin)/pi180.or.&
              &ReadStat%latSt<yOrig/pi180.or.ReadStat%latSt>(yOrig+real(Jmax)*dY)/pi180) cycle
      if(typeSt(:1)=='a'.or.typeSt(2:)=='a') then
        Nair=Nair+1
        AirStat(Nair)%nameSt=ReadStat%nameSt
        AirStat(Nair)%indSt=ReadStat%indSt
        AirStat(Nair)%latSt=ReadStat%latSt
        AirStat(Nair)%longSt=ReadStat%longSt
      endif
      if(typeSt(:1)=='p'.or.typeSt(2:)=='p') then
        Nprec=Nprec+1
        PrecStat(Nprec)%nameSt=ReadStat%nameSt
        PrecStat(Nprec)%indSt=ReadStat%indSt
        PrecStat(Nprec)%latSt=ReadStat%latSt
        PrecStat(Nprec)%longSt=ReadStat%longSt
      endif
    enddo
    close(2)

! Check if stations are defined for monitoring output
    if(AtmOutMonitYearly.or.AtmOutMonitMonthly.or.AtmOutMonitDaily.or.AtmOutMonit6hourly.or.AtmOutMonitHourly) then
      if(NAir == 0) then
        print '(/,"STOP: Number of air monitoring sites is zero ",/)'
        stop
      endif
      if(NPrec == 0) then
        print '(/,"STOP: Number of precip monitoring sites is zero ",/)'
        stop
      endif
    end if

    do st=1, Nair
      Dist=1.e10
      do j=Jmin, Jmax
        latM=LatMesh(j)/pi180
        do i=Imin, Imax
          longM=LongMeshR(i)/pi180
          Rdist=AngleDist(AirStat(st)%longSt,AirStat(st)%latSt,longM,latM)
          if(longM<=AirStat(st)%longSt.and.latM<=AirStat(st)%latSt.and.Rdist<Dist(1,1)) then
            AirStat(st)%iSt(1,1)=i
            AirStat(st)%jSt(1,1)=j
            Dist(1,1)=max(Rdist,Eps)
          endif
          if(longM<=AirStat(st)%longSt.and.latM>AirStat(st)%latSt.and.Rdist<Dist(1,2)) then
            AirStat(st)%iSt(1,2)=i
            AirStat(st)%jSt(1,2)=j
            Dist(1,2)=max(Rdist,Eps)
          endif
          if(longM>AirStat(st)%longSt.and.latM<=AirStat(st)%latSt.and.Rdist<Dist(2,1)) then
            AirStat(st)%iSt(2,1)=i
            AirStat(st)%jSt(2,1)=j
            Dist(2,1)=max(Rdist,Eps)
          endif
          if(longM>AirStat(st)%longSt.and.latM>AirStat(st)%latSt.and.Rdist<Dist(2,2)) then
            AirStat(st)%iSt(2,2)=i
            AirStat(st)%jSt(2,2)=j
            Dist(2,2)=max(Rdist,Eps)
          endif
        enddo
      enddo
      long11=LongMeshR(AirStat(st)%iSt(1,1))/pi180
      long12=LongMeshR(AirStat(st)%iSt(1,2))/pi180
      long21=LongMeshR(AirStat(st)%iSt(2,1))/pi180
      long22=LongMeshR(AirStat(st)%iSt(2,2))/pi180

      lat11=LatMesh(AirStat(st)%jSt(1,1))/pi180
      lat12=LatMesh(AirStat(st)%jSt(1,2))/pi180
      lat21=LatMesh(AirStat(st)%jSt(2,1))/pi180
      lat22=LatMesh(AirStat(st)%jSt(2,2))/pi180

      if(long11==long12.and.long21==long22) then
        AirStat(st)%Ast(1,1)=abs(long22-AirStat(st)%longSt)*abs(lat22-AirStat(st)%latSt)/abs(long22-long11)/abs(lat22-lat11)
        AirStat(st)%Ast(1,2)=abs(long21-AirStat(st)%longSt)*abs(lat21-AirStat(st)%latSt)/abs(long21-long12)/abs(lat21-lat12)
        AirStat(st)%Ast(2,1)=abs(long12-AirStat(st)%longSt)*abs(lat12-AirStat(st)%latSt)/abs(long12-long21)/abs(lat12-lat21)
        AirStat(st)%Ast(2,2)=abs(long11-AirStat(st)%longSt)*abs(lat11-AirStat(st)%latSt)/abs(long11-long22)/abs(lat11-lat22)
      else
        wt(1,1)=1./Dist(1,1)/Dist(1,1)
        wt(1,2)=1./Dist(1,2)/Dist(1,2)
        wt(2,1)=1./Dist(2,1)/Dist(2,1)
        wt(2,2)=1./Dist(2,2)/Dist(2,2)
        wSum=wt(1,1)+wt(1,2)+wt(2,1)+wt(2,2)
         
        AirStat(st)%Ast(1,1)=wt(1,1)/wSum
        AirStat(st)%Ast(1,2)=wt(1,2)/wSum
        AirStat(st)%Ast(2,1)=wt(2,1)/wSum
        AirStat(st)%Ast(2,2)=wt(2,2)/wSum
      endif
    enddo

    do st=1, Nprec
      Dist=1.e10
      do j=Jmin, Jmax
        latM=LatMesh(j)/pi180
        do i=Imin, Imax
          longM=LongMeshR(i)/pi180
          Rdist=AngleDist(PrecStat(st)%longSt,PrecStat(st)%latSt,longM,latM)
          if(longM<=PrecStat(st)%longSt.and.latM<=PrecStat(st)%latSt.and.Rdist<Dist(1,1)) then
            PrecStat(st)%iSt(1,1)=i
            PrecStat(st)%jSt(1,1)=j
            Dist(1,1)=max(Rdist,Eps)
          endif
          if(longM<=PrecStat(st)%longSt.and.latM>PrecStat(st)%latSt.and.Rdist<Dist(1,2)) then
            PrecStat(st)%iSt(1,2)=i
            PrecStat(st)%jSt(1,2)=j
            Dist(1,2)=max(Rdist,Eps)
          endif
          if(longM>PrecStat(st)%longSt.and.latM<=PrecStat(st)%latSt.and.Rdist<Dist(2,1)) then
            PrecStat(st)%iSt(2,1)=i
            PrecStat(st)%jSt(2,1)=j
            Dist(2,1)=max(Rdist,Eps)
          endif
          if(longM>PrecStat(st)%longSt.and.latM>PrecStat(st)%latSt.and.Rdist<Dist(2,2)) then
            PrecStat(st)%iSt(2,2)=i
            PrecStat(st)%jSt(2,2)=j
            Dist(2,2)=max(Rdist,Eps)
          endif
        enddo
      enddo
      long11=LongMeshR(PrecStat(st)%iSt(1,1))/pi180
      long12=LongMeshR(PrecStat(st)%iSt(1,2))/pi180
      long21=LongMeshR(PrecStat(st)%iSt(2,1))/pi180
      long22=LongMeshR(PrecStat(st)%iSt(2,2))/pi180

      lat11=LatMesh(PrecStat(st)%jSt(1,1))/pi180
      lat12=LatMesh(PrecStat(st)%jSt(1,2))/pi180
      lat21=LatMesh(PrecStat(st)%jSt(2,1))/pi180
      lat22=LatMesh(PrecStat(st)%jSt(2,2))/pi180

      if(long11==long12.and.long21==long22) then
        PrecStat(st)%Ast(1,1)=abs(long22-PrecStat(st)%longSt)*abs(lat22-PrecStat(st)%latSt)/abs(long22-long11)/abs(lat22-lat11)
        PrecStat(st)%Ast(1,2)=abs(long21-PrecStat(st)%longSt)*abs(lat21-PrecStat(st)%latSt)/abs(long21-long12)/abs(lat21-lat12)
        PrecStat(st)%Ast(2,1)=abs(long12-PrecStat(st)%longSt)*abs(lat12-PrecStat(st)%latSt)/abs(long12-long21)/abs(lat12-lat21)
        PrecStat(st)%Ast(2,2)=abs(long11-PrecStat(st)%longSt)*abs(lat11-PrecStat(st)%latSt)/abs(long11-long22)/abs(lat11-lat22)
      else
        wt(1,1)=1./Dist(1,1)/Dist(1,1)
        wt(1,2)=1./Dist(1,2)/Dist(1,2)
        wt(2,1)=1./Dist(2,1)/Dist(2,1)
        wt(2,2)=1./Dist(2,2)/Dist(2,2)
        wSum=wt(1,1)+wt(1,2)+wt(2,1)+wt(2,2)

        PrecStat(st)%Ast(1,1)=wt(1,1)/wSum
        PrecStat(st)%Ast(1,2)=wt(1,2)/wSum
        PrecStat(st)%Ast(2,1)=wt(2,1)/wSum
        PrecStat(st)%Ast(2,2)=wt(2,2)/wSum
      endif
    enddo

#ifdef DEBUG_MODE
    print *, '<- Exit MonitorSites'
#endif

end subroutine MonitorSites

end module Atm_Input
