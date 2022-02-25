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
! Module containing Hg specific procedures
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module Atm_Hg_General

#ifdef G_HG

  use GeneralParams
  use Atm_Hg_Chemistry

  implicit none

  character(800), private :: fileName, fullName
  character(4),  private :: YearNum
  character(2),  private :: DayNum, YearShort


contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine including mercury specific processes
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_Input(Prd)

    character(*) Prd

    selectcase(Prd)
!---------------------------------------------------------------------------------
    case('Initial')

!---------------------------------------------------------------------------------
    case('Yearly')

!---------------------------------------------------------------------------------
    case('Monthly')

!---------------------------------------------------------------------------------
    case('Daily')
       select case (ReactSource)
           case('GEOS_Chem_2x2.5')
               call Atm_ReadReactDaily_GEOSChem
           case('GEOS_Chem_4x5')
               call Atm_ReadReactDaily_GEOSChem
           case('MOZART')
               call Atm_ReadReactDaily_MOZART
           case default
               print '(/,"STOP: Invalid source of reactant data ''",a,"''",/)', trim(ReactSource)
               stop
      end select
      
    call ReactTransformDaily

!---------------------------------------------------------------------------------
    case('6hourly')

!---------------------------------------------------------------------------------
    endselect

end subroutine Atm_Hg_Input

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine including mercury specific processes
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_Process(Subs)

    integer Subs
    integer, save :: order=1

    call ReactTransformSteply

    call Atm_Hg_Partition

    if(order>0) then
      call Atm_Hg_DropChem
!      call Atm_Hg_GasChem_Br
!      call Atm_Hg_GasChem_Br2
      call Atm_Hg_GasChem_O3_OH
    else
!      call Atm_Hg_GasChem_Br
!      call Atm_Hg_GasChem_Br2
      call Atm_Hg_GasChem_O3_OH
      call Atm_Hg_DropChem
    endif
    order=-order

    call Atm_Hg_Partition

end subroutine Atm_Hg_Process


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine defining physical and chemical properties of Hg
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_Props(Subs)

    integer i, FileStat, lenType, lenType1, Subs, Form, Limit, Ind, Src, nForm0
    character(720) strRead, strPar
    character(30) strType, FormType(MaxForm)
    integer mon, ios
    real emisAmp(MaxForm)    

    fileName='Atm'//trim(PropName)//trim(SubsID(Subs))//'.dat'
    fullName=trim(PropPath)//fileName
    open(2, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
      stop
    endif

    do
      read(2,'(a)',iostat=ios) strRead
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

! Reading pollutant form
      if(strRead(1:1)=='[') then
        Limit=scan(strRead(2:),']')
        lenType=Limit
        strType=trim(strRead(2:lenType))

        Form=0
        do i=1, sFormNum(Atm,Subs)
          if(strType==FormType(i)) then
            Form=i
            exit
          endif
        enddo
        if(Form==0) then
          print '(/,"STOP: Unknown pollutant form ''",a,"''",/)', trim(strType)
          stop
        endif

! Reading parameter type
      else
        lenType=scan(strRead,':')
        if(lenType==0) then
          print '(/,"STOP: Wrong format of file ''",a,"''",/)', trim(fileName)
          stop 
        endif
        strType=trim(strRead(:lenType-1))
        strPar=strRead(lenType+1:)
        strPar=adjustl(strPar)

        selectcase(strType)
!------------------------------------------------------------------
          case('Forms number')
            read(strPar,'(i2)') sFormNum(Atm,Subs)
            if (sFormNum(Atm,Subs)>0) then
              nForm0=NumForm(Atm)
              FormSubs(Atm,nForm0+1:nForm0+sFormNum(Atm,Subs))=Subs
              NumSubsMedia(Atm)=NumSubsMedia(Atm)+1
              gSubsMediaInd(Atm,Subs)=NumSubsMedia(Atm)
            endif
!------------------------------------------------------------------
          case('Forms')
            read(strPar,*) (FormType(i), i=1, sFormNum(Atm,Subs))
!------------------------------------------------------------------
          case('Form ID')
            FormID(Atm,Subs,Form)=strPar
            sFormInd(Atm,Subs,Form)=nForm0+Form
            FormSubs(Atm,sFormInd(Atm,Subs,Form))=Subs
            NumForm(Atm)=NumForm(Atm)+1
!------------------------------------------------------------------
          case('Dry deposition')
            if(strPar=='yes') then
              DryDepNum=DryDepNum+1
              DryDepInd(DryDepNum)=sFormInd(Atm,Subs,Form)
            endif
!------------------------------------------------------------------
          case('Below-cloud scavenging')
            if(strPar=='yes') then
              BCScvNum=BCScvNum+1
              BCScvInd(BCScvNum)=sFormInd(Atm,Subs,Form)
            endif
!------------------------------------------------------------------
          case('In-cloud scavenging')
            if(strPar=='yes') then
              ICScvNum=ICScvNum+1
              ICScvInd(ICScvNum)=sFormInd(Atm,Subs,Form)
            endif
!------------------------------------------------------------------
          case('Gas exchange')
            if(strPar=='yes') then
              GasExchNum=GasExchNum+1
              GasExchInd(GasExchNum)=sFormInd(Atm,Subs,Form)
            endif
!------------------------------------------------------------------
          case('Atm transport')
            if(strPar=='yes') then
              AtmTransNum=AtmTransNum+1
              AtmTransInd(AtmTransNum)=sFormInd(Atm,Subs,Form)
            endif
!------------------------------------------------------------------
          case('Boundary conditions')
            if(strPar=='yearly') then
              AtmBndNum=AtmBndNum+1
              ReadBndInd(AtmBndNum,yr)=sFormInd(Atm,Subs,Form)
              AtmBndInd(AtmBndNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='monthly') then
              AtmBndNum=AtmBndNum+1
              ReadBndInd(AtmBndNum,mn)=sFormInd(Atm,Subs,Form)
              AtmBndInd(AtmBndNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='daily') then
              AtmBndNum=AtmBndNum+1
              ReadBndInd(AtmBndNum,da)=sFormInd(Atm,Subs,Form)
              AtmBndInd(AtmBndNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='no') then
              continue
            else
              print '(/,"STOP: Incorrect periodicity of boundary conditions ''",a,"''",/)', trim(strPar)
              stop
            endif
!------------------------------------------------------------------
          case('Anthrop emiss')
            if(strPar=='yearly') then
              AntEmisNum=AntEmisNum+1
              ReadAntInd(AntEmisNum,yr)=sFormInd(Atm,Subs,Form)
              AntEmisInd(AntEmisNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='monthly') then
               AntEmisNum=AntEmisNum+1
              ReadAntInd(AntEmisNum,mn)=sFormInd(Atm,Subs,Form)
              AntEmisInd(AntEmisNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='daily') then
              AntEmisNum=AntEmisNum+1
              ReadAntInd(AntEmisNum,da)=sFormInd(Atm,Subs,Form)
              AntEmisInd(AntEmisNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='no') then
              continue
            else
              print '(/,"STOP: Incorrect periodicity of anthropogenic emission ''",a,"''",/)', trim(strPar)
              stop
            endif
!------------------------------------------------------------------
          case('Natural emiss')
            if(strPar=='yearly') then
               NatEmisNum=NatEmisNum+1
              ReadNatInd(NatEmisNum,yr)=sFormInd(Atm,Subs,Form)
              NatEmisInd(NatEmisNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='monthly') then
              NatEmisNum=NatEmisNum+1
              ReadNatInd(NatEmisNum,mn)=sFormInd(Atm,Subs,Form)
              NatEmisInd(NatEmisNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='daily') then
              NatEmisNum=NatEmisNum+1
              ReadNatInd(NatEmisNum,da)=sFormInd(Atm,Subs,Form)
              NatEmisInd(NatEmisNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='no') then
              continue
            else
              print '(/,"STOP: Incorrect periodicity of natural emission ''",a,"''",/)', trim(strPar)
              stop
            endif
!------------------------------------------------------------------
          case('Reemission')
            if(strPar=='yearly') then
              ReEmisNum=ReEmisNum+1
              ReadReInd(ReEmisNum,yr)=sFormInd(Atm,Subs,Form)
              ReEmisInd(ReEmisNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='monthly') then
              ReEmisNum=ReEmisNum+1
              ReadReInd(ReEmisNum,mn)=sFormInd(Atm,Subs,Form)
              ReEmisInd(ReEmisNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='daily') then
              ReEmisNum=ReEmisNum+1
              ReadReInd(ReEmisNum,da)=sFormInd(Atm,Subs,Form)
              ReEmisInd(ReEmisNum)=sFormInd(Atm,Subs,Form)
            elseif(strPar=='no') then
              continue
            else
              print '(/,"STOP: Incorrect periodicity of re-emission input ''",a,"''",/)', trim(strPar)
              stop
            endif
!------------------------------------------------------------------
          case('Aqueous form')
            if(strPar=='yes') then
              AqFormNum=AqFormNum+1
              AqFormInd(AqFormNum)=sFormInd(Atm,Subs,Form)
            endif
!------------------------------------------------------------------
          case('Particle diameter, m')
            read(strPar,*) Dp
!------------------------------------------------------------------
          case('Particle density, kg/m3')
            read(strPar,*) RhoP
!------------------------------------------------------------------
          case('Vd coeffs')
            read(strPar,*) (Vdcoef(i), i=1,3)
!------------------------------------------------------------------
          case('Surface resistance')
            read(strPar,*) Rcmin
!------------------------------------------------------------------
          case('BCS coeffs')
            read(strPar,*) Abelow, Bbelow
!------------------------------------------------------------------
          case('ICS coeffs')
            read(strPar,*) Ain, Bin
!------------------------------------------------------------------
        case('Washout effectiveness coeff')
            read(strPar,*) Aeff
!------------------------------------------------------------------
        case('Aerosol solubility')
            read(strPar,*) Asol
!------------------------------------------------------------------
          case('Emis amplitude')
            read(strPar,*) emisAmp(Form)
!------------------------------------------------------------------
          case('Henry const (Hg)')
            read(strPar,*) (CHenryHg(i), i=1,2)
!------------------------------------------------------------------
          case('Henry const (O3)')
            read(strPar,*) (CHenryO3(i), i=1,2)
!------------------------------------------------------------------
          case('Henry const (HgCl2)')
            read(strPar,*) (CHenryHgCl2(i), i=1,2)
!------------------------------------------------------------------
          case('Henry const (Cl2)')
            read(strPar,*) (CHenryCl2(i), i=1,2)
!------------------------------------------------------------------
          case('Henry const (OH)')
            read(strPar,*) (CHenryOH(i), i=1,2)
!------------------------------------------------------------------
          case('Hg ion/chloride ratio')
            read(strPar,*) (Clcoef(i), i=1,4)
!------------------------------------------------------------------
          case('Chlorine ion conc')
            read(strPar,*) ConcCl
!------------------------------------------------------------------
          case('Adsorption ratio')
            read(strPar,*) Rsoot
!------------------------------------------------------------------
          case('Reaction oxid (O3)')
            read(strPar,*) (O3coef(i), i=1,2)
!------------------------------------------------------------------
          case('Reaction oxid (Cl2)')
            read(strPar,*) (Cl2coef(i), i=1,2)
!------------------------------------------------------------------
          case('Reaction oxid (OH)')
            read(strPar,*) (OHcoef(i), i=1,2)
!------------------------------------------------------------------
          case('Reaction A->C (O3)')
            read(strPar,*) Kac1
!------------------------------------------------------------------
          case('Reaction A->C (Cl2)')
            read(strPar,*) Kac2
!------------------------------------------------------------------
          case('Reaction A->C (OH)')
            read(strPar,*) Kac3
!------------------------------------------------------------------
          case('Reaction B->A')
            read(strPar,*) Kba
!------------------------------------------------------------------
          case('Reaction C->B')
            read(strPar,*) Kcb
!------------------------------------------------------------------
          case('Reaction C->A')
            read(strPar,*) Kca
!------------------------------------------------------------------
          case('pH')
            read(strPar,*) pH
!------------------------------------------------------------------
          case default
            print '(/,"STOP: Unknown input parameter ''",a,"''",/)', trim(strType)
            stop
!------------------------------------------------------------------
        endselect

      endif
    enddo

    close(2)

    if(NumForm(Atm)<nForm0+sFormNum(Atm,Subs)) then
      do i=1, sFormNum(Atm,Subs)
        if(sFormInd(Atm,Subs,i)/=nForm0+i) then
        print '(/,"STOP: No information for pollutant form ''",a,"''",/)', trim(FormType(i))
        stop
        endif
      enddo
    endif

! Identifying pollutant forms
    Npart=0
    Noxid=0
    do Ind=1, sFormNum(Atm,Subs)
      Form=sFormInd(Atm,Subs,Ind)
      selectcase(FormID(Atm,Subs,Ind))
        case('0_gas')
          Hg0=Form
          Gas=Form
        case('0_dis')
          Dis=Form
        case('_sulf')
          Sulf=Form
        case('_chlor')
          Chlor=Form
        case('Br')
          HgBr=Form
          Noxid=Noxid+1
          Oxid(Noxid)=Form
        case('Br_part')
          HgBr_part=Form
          Npart=Npart+1
          Part(Npart)=Form
        case('Br2','II_gas')
          HgBr2=Form
          Noxid=Noxid+1
          Oxid(Noxid)=Form
        case('Br2_part','_part')
          HgBr2_part=Form
          Npart=Npart+1
          Part(Npart)=Form
        case('BrOH')
          HgBrOH=Form
          Noxid=Noxid+1
          Oxid(Noxid)=Form
        case('BrOH_part')
          HgBrOH_part=Form
          Npart=Npart+1
          Part(Npart)=Form
        case('BrOOH')
          HgBrOOH=Form
          Noxid=Noxid+1
          Oxid(Noxid)=Form
        case('BrOOH_part')
          HgBrOOH_part=Form
          Npart=Npart+1
          Part(Npart)=Form
        case('BrONO')
          HgBrONO=Form
          Noxid=Noxid+1
          Oxid(Noxid)=Form
        case('BrONO_part')
          HgBrONO_part=Form
          Npart=Npart+1
          Part(Npart)=Form
        case('OH')
          HOHg=Form
          Noxid=Noxid+1
          Oxid(Noxid)=Form
        case('OH_part')
          HOHg_part=Form
          Npart=Npart+1
          Part(Npart)=Form
        case('OHOOH')
          HOHgOOH=Form
          Noxid=Noxid+1
          Oxid(Noxid)=Form
        case('OHOOH_part')
          HOHgOOH_part=Form
          Npart=Npart+1
          Part(Npart)=Form
        case('OHONO')
          HOHgONO=Form
          Noxid=Noxid+1
          Oxid(Noxid)=Form
        case('OHONO_part')
          HOHgONO_part=Form
          Npart=Npart+1
          Part(Npart)=Form
        case('O_part')
          HgO_part=Form
          Npart=Npart+1
          Part(Npart)=Form
        case default
          print '(/,"STOP: Unknown pollutant form ''",a,"''",/)', FormType(Form)
          stop
      endselect
    enddo
    
! Identifying forms for gas-particle partitioning    
    do Ind=1, Noxid
      Form=Oxid(Ind)
      if(Form==HgBr) then
        GasPart(1:2,Ind)=(/HgBr,HgBr_part/)
      elseif(Form==HgBr2) then
        GasPart(1:2,Ind)=(/HgBr2,HgBr2_part/)
      elseif(Form==HgBrOH) then
        GasPart(1:2,Ind)=(/HgBrOH,HgBrOH_part/)
      elseif(Form==HgBrOOH) then
        GasPart(1:2,Ind)=(/HgBrOOH,HgBrOOH_part/)
      elseif(Form==HgBrONO) then
        GasPart(1:2,Ind)=(/HgBrONO,HgBrONO_part/)
      elseif(Form==HOHg) then
        GasPart(1:2,Ind)=(/HOHg,HOHg_part/)
      elseif(Form==HOHgOOH) then
        GasPart(1:2,Ind)=(/HOHgOOH,HOHgOOH_part/)
      elseif(Form==HOHgONO) then
        GasPart(1:2,Ind)=(/HOHgONO,HOHgONO_part/)
      endif
    enddo

#if RTYPE==2
    if(NumAnth==0.and.AntEmisNum/=0) then
      print '(/,"STOP: No emission sources is defined for anthropogenic emissions of forms:",<AntEmisNum>(1x,a),/)',&
                            &(trim(SubsID(Subs))//trim(FormID(Atm,Subs,AntEmisInd(Ind))), Ind=1, AntEmisNum)
      stop
    elseif(NumAnth/=0.and.AntEmisNum==0) then
      print '(/,"STOP: No anthropogenic emission is defined for sources:",<NumAnth>(1x,a),/)',&
                            &(trim(SourcID(Src)), Src=1, NumAnth)
      stop
    endif    

    if(NumNat==0.and.NatEmisNum/=0) then
      print '(/,"STOP: No emission sources is defined for natural emissions of forms:",<NatEmisNum>(1x,a),/)',&
                            &(trim(SubsID(Subs))//trim(FormID(Atm,Subs,NatEmisInd(Ind))), Ind=1, NatEmisNum)
      stop
    elseif(NumNat/=0.and.NatEmisNum==0) then
      print '(/,"STOP: No natural emission is defined for sources:",<NumNat>(1x,a),/)',&
                            &(trim(NaturID(Src)), Src=1, NumNat)
      stop
    endif

#if (REGTYPE==2)
    if(NumBnd==0.and.AtmBndNum/=0) then
      print '(/,"STOP: No sources is defined for boundary concentrations of forms:",<AtmBndNum>(1x,a),/)',&
                            &(trim(SubsID(Subs))//trim(FormID(Atm,Subs,AtmBndInd(Ind))), Ind=1, AtmBndNum)
      stop
    elseif(NumBnd/=0.and.AtmBndNum==0) then
      print '(/,"STOP: No boundary concentrations is defined for sources:",<NumBnd>(1x,a),/)',&
                            &(trim(BoundID(Src)), Src=1, NumBnd)
      stop
    endif
#endif
#endif

end subroutine Atm_Hg_Props


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading emission fields
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_ReadNatEmis(Subs,Ind,tm,Note)

    integer i, j, Ind, Form, Subs, FileStat, Hs, Xscal, lenStr, SrcNat, ios, tm, ln, Nread
    integer Src, Src2, SrcInd
    real(8) Emis, ConvEmisUnit, Ecoeff(12)
    real Aver(Imin:Imax)
    character(10) IDsrc, IDsubs, IDform, TResol, prd, IDread(0:MaxMatr)
    character(30) parType, gridID, units
    character(300) readStr, parVal, SrcInfo
    data Ecoeff /1.868,1.845,1.986,2.047,1.830,1.610,1.444,1.431,1.352,1.373,1.367,1.660/
    character(*) Note

    Form=NatEmisInd(Ind)
    selectcase(tm)
      case(yr)
        write(fileName,'(i4)') Year
        fileName=trim(SubsID(Subs))//trim(FormID(Atm,Subs,Form))//'_'//trim(fileName)//trim(NatEmisName)//'.dat'
        ConvEmisUnit=1000./sum(MonthDays)/SecInDay      ! Convert t/y ==> kg/s
      case(mn)
        write(fileName,'(i4,i2.2)') Year, Month
        fileName=trim(SubsID(Subs))//trim(FormID(Atm,Subs,Form))//'_'//trim(fileName)//trim(NatEmisName)//'.dat'
!fileName=trim(SubsID(Subs))//trim(FormID(Atm,Subs,Form))//'_'//trim(NatEmisName)//MonthName(Month)//'.dat'
        ConvEmisUnit=1000./MonthDays(Month)/SecInDay    ! Convert t/m ==> kg/s
      case(da)
        write(fileName,'(i4,i2.2,i2.2)') Year, Month, Day
        fileName=trim(SubsID(Subs))//trim(FormID(Atm,Subs,Form))//'_'//trim(fileName)//trim(NatEmisName)//'.dat'
        ConvEmisUnit=1000./SecInDay                     ! Convert t/day ==> kg/s
    endselect

    write(YearNum,'(i4)') Year
    fullName=trim(EmisPath)//trim(GridCode)//'/'//trim(SubsID(Subs))//trim(AntEmisName)//'/'//YearNum//'/'//fileName
    open(4, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
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
          IDform=parVal
        case('Grid Type')
          gridID=parVal
        case('Source')
          SrcInfo=parVal
        case('TimeResol')
          TResol = ParVal
          if (tm == yr .and. trim(TResol) /= 'yearly') then
            print*, 'STOP: Invalid time resolution of emission data'
            stop
          end if
          if (tm == mn .and. trim(TResol) /= 'monthly') then
            print*, 'STOP: Invalid time resolution of emission data'
            stop
          end if
          if (tm == da .and. trim(TResol) /= 'daily') then
            print*, 'STOP: Invalid time resolution of emission data'
            stop
          end if
        case('Period')
          prd = parVal
        case('Units')
          units=parVal
        case('Layers')
          read(parVal,*) HemisNat
        case default
          print '(/,"STOP: Unknown input parameter ''",a,"''",/)', trim(parType)
          stop
      endselect
    enddo

    write(Note,'("Emissions of ",a,a," in ",a,", ",a)') trim(IDsubs), trim(IDform), prd, trim(units)

    read(4,*)
    read(4,*)
    read(4,*)

    SrcNat=0
    NatEmisFlux(:,:,Ind,:)=0.
    IDread=''
    Nread=0
    do
      read(4,'(a)',iostat=ios) readStr
      if(ios==-1) exit
      lenStr=scan(readStr,achar(13))
      if(lenStr>0) readStr=readStr(:lenStr-1)
      read(readStr,*) IDsrc, i, j, Emis

!***************** Matrix calculations *****************
#if RTYPE==2
      if(IDread(Nread)/=IDsrc) then
        Nread=Nread+1
        IDread(Nread)=IDsrc
      endif

! Identifying the source
      SrcInd=0
      do Src=1, NumNat
        if(IDsrc==NaturID(Src)) then
          SrcInd=Src
          exit
        endif
      enddo
      if(SrcInd>0) then
        NatEmisFlux(i,j,Ind,SrcInd)=NatEmisFlux(i,j,Ind,SrcInd)+Emis   ! t/mon
      elseif(natSRCmode==2) then
        NatEmisFlux(i,j,Ind,restNat)=NatEmisFlux(i,j,Ind,restNat)+Emis   ! t/mon
      endif

#else
!*******************************************************
      NatSrcInd(1)=1
      NatEmisFlux(i,j,Ind,1)=NatEmisFlux(i,j,Ind,1)+Emis               ! t/mon
#endif
    enddo
    close(4)
    NatEmisFlux(:,:,Ind,:)=NatEmisFlux(:,:,Ind,:)*ConvEmisUnit           ! t/mon -> kg/sec
    NatEmisFlux(:,:,Ind,:)=NatEmisFlux(:,:,Ind,:)

!***************** Matrix calculations *****************
#if RTYPE==2
! Checking completeness of emission data
    do Src=NumAnth+1, NumAnth+NumNat
      if(SourcID(Src)=='restNat') cycle
      SrcInd=0
      do Src2=1, Nread
        if(trim(IDread(Src2))==trim(SourcID(Src))) then
          SrcInd=Src2
          exit
        endif
      enddo
      if(SrcInd==0) then
        print '(/,"STOP: No natural emission data for source ''",a,"''",/)', trim(SourcID(Src))
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

      do Src=1, NumNat
        Aver(Imin:Imax)=NatEmisFlux(Imin:Imax,j,Ind,Src)
        call GridAggreg(j,Xscal,Aver,2)
        NatEmisFlux(minI(j):maxI(j),j,Ind,Src)=Aver(minI(j):maxI(j))
      enddo
    enddo

end subroutine Atm_Hg_ReadNatEmis


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating dry deposition velosity (Vd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Atm_Hg_ExchPar(Ind,nDay)

    integer i, j, t, Surf, Form, Ind, s, it, SurfG, nDay, k
    real Pxx, Pa, RT, Rho, Tair, Psat, Qv, RH, Mu, Lmo0, Lmo(MaxSurf), Almo, Zr, Ui, Uj, Uref, Kcond, Dh, SnowH, Tsurf
    real DpW, Kn, Cunn, Vg, DiffW, DpS, VgS, DiffS, VgF, DiffF, Ma, DiffG, Cpw, Dfor, Hfor, ZoM, ZoH, Ux0, Ux
    real PrecN(Imin:Imax,Jmin:Jmax), PrecB(Imin:Imax,Jmin:Jmax), Lcover(MaxSurf) !, LCmask(MaxSurf)
    real St, Sc, Uh, Eb, Ein, Eim, Reb, Eff, Fbrok, Pr, Tcoef, Lw, VdFog
    real Ra, Rs, Rs1, Rs2, Rb, Rc, Zr0
    logical WetGrass, WetForest

    Form=DryDepInd(Ind)

    Vd(:,:,:,Ind,:,nDay)=0.

    do t=1, NumPer
      do j=Jmin, Jmax
        do i=minI(j), maxI(j)

! Calculation of air parameters
          Lcover(1:NumSurf)=LandCover(i,j,1:NumSurf)                    ! Land cover types
          s=Season(i,j,Month)                                            ! Seasons (1-5)
          Pxx=Px(i,j,t,nDay)                                            ! Surface pressure Px=Ps-Pt
          Pa=Sigma(1)*Pxx+Ptop                                            ! Pressure at the lowest sigma level
          RT=RTemp(i,j,1,t,nDay)                                        ! RT of air at the lowest sigma level
          Rho=Pa/RT                                                        ! Air density at the lowest sigma level
          Lw=LiqCont(i,j,1,t)+FrozCont(i,j,1,t)                            ! Liquid or frozen water content
          Tair=TempAir(i,j,1,t,nDay)                                    ! Temperature at the lowest sigma level
            Tsurf=TempSurf(i,j,t)
          if(Tair>273.15) then
            Psat=6.112*exp(17.67*(Tair-273.15)/(Tair-29.65))*100.        ! Saturation vapour pressure (t>0) [Pa]    
          else
            Psat=6.11*exp(22.514-6150./Tair)*100.                        ! Saturation vapour pressure (t<0) [Pa]                
          endif
          SnowH=SnowDepth(i,j,t)
          Qv=HumidAir(i,j,1,t,nDay)                                        ! Vapour mioxing ratio
          RH=min(Qv*(1.+HumCoef)/(1.+Qv*(1.+HumCoef))*Pa/Psat,1.)        ! Relative humidity (saturation ratio)
          Mu=Nu*Rho                                                        ! Dynamic viscosity 
          Zr=RT/Ggrav*log((Pxx+Ptop)/(Sigma(1)*Pxx+Ptop))                ! Height of the lowest sigma level
          Ui=(Uwind(i,j,1,t,nDay)+Uwind(i-1,j,1,t,nDay))/2.
          if(j>Jmin) then
            Uj=(Vwind(i,j,1,t,nDay)+Vwind(iS(i,j,1),j-1,1,t,nDay))/2.                   ! Shift due to grid aggregation
          else
            Uj=Vwind(i,j,1,t,nDay)
          endif
          Uref=sqrt(Ui*Ui+Uj*Uj)                                        ! Absolute value of vind velosity
          Ux0=Ufric(i,j,t,nDay)                                            ! Mean friction velocity
          Lmo0=MOLength(i,j,t,nDay)                                        ! Mean Monin-Obuhov length
          Lmo(1:NumSurf)=Lmo0                                            ! Surface dependant MO length
          Almo=Lmo0/Ux0/Ux0/Ux0                                            ! MO length coefficient
          Kcond=0.023807+7.1128e-5*(Tair-273.15)                        ! Thermal conductivity of dry air
          Dh=Kcond/Rho/Cpd                                                ! Molecular thermal diffusion coefficient
          Ma=Md/(1.+HumCoef*Qv/(1.+Qv))                                    ! Malecular weight of total air
          DiffG=3./8./Nav/Dm/Dm/Rho*sqrt(Runiv*Tair*Ma*(Ma+MHgCL2)/MHgCL2/2./PiNum) ! Diffusion coefficient of RGM
          Cpw=Cpd*(1.+CpCoef*Qv/(1.+Qv))                                ! Specific heat of moist air
          Pr=Mu*Cpw/Kcond                                                ! Prandtl number

           DpW=WetDiameter(Dp,RH,1.)                                        ! Wet diameter of particles
          Kn=4.*Nu/sqrt(8.*kB*Tair*Nav/PiNum/Md)/DpW                    ! Knudsen number
          Cunn=1+Kn*(1.249+0.42*exp(-0.87/Kn))                            ! Cunningham correction
          Vg=DpW*DpW*RhoP*Ggrav/18./Mu*Cunn                                ! Gravitational settling velosity
          DiffW=kB*Tair/3./PiNum/DpW/Mu*Cunn                            ! Particle diffusion coefficient in moist air 

          if(sum(Lcover(gInd(Water,1:gNum(Water))))>0.)    then            ! 
            DpS=WetDiameter(Dp,0.98,1.)                                    !
            Kn=4.*Nu/sqrt(8.*kB*Tair*Nav/PiNum/Md)/DpS                    ! Particles in water surface layer
            Cunn=1+Kn*(1.249+0.42*exp(-0.87/Kn))                        !
            VgS=DpS*DpS*RhoP*Ggrav/18./Mu*Cunn                            !  
            DiffS=kB*Tair/3./PiNum/DpS/Mu*Cunn                            !
          endif                                                            !

          if(Lw>0..and.SubsID(1)=='Hg') then                                ! 
            Kn=4.*Nu/sqrt(8.*kB*Tair*Nav/PiNum/Md)/Dfog                    ! 
            Cunn=1+Kn*(1.249+0.42*exp(-0.87/Kn))                        ! Fog deposition
            VgF=Dfog*Dfog*RhoWat*Ggrav/18./Mu*Cunn                        !  
            DiffF=kB*Tair/3./PiNum/Dfog/Mu*Cunn                            !
          endif                                                            !

! Calculation of surface wetness (for grass and forest)
          PrecN(i,j)=PrecRainConv(i,j,1,t)+PrecRainStrat(i,j,1,t)        ! Current precipitation rate
          if(PrecN(i,j)>0.) then
            WetGrass=.true.                                                ! Wet grass surface
          else
            WetGrass=.false.                                            ! Dry grass surface
          endif
          if(PrecN(i,j)>0..or.PrecB(i,j)>0.) then
            WetForest=.true.                                            ! Wet forest surface
          else
            WetForest=.false.                                            ! Dry foret surface
          endif
          PrecB(i,j)=PrecN(i,j)                                            ! Previous precipitation rate

          do Surf=1, NumSurf
            if(Lcover(Surf)<=0.) cycle
            Hfor=heightLU(Surf)                                            ! Canopy height
            Dfor=dispLU(Surf,s)                                            ! Displacement height

            if(SurfGroup(Surf)==Water) then
              ZoM=0.016*Ux0*Ux0/Ggrav+Nu/9.1*Ux0                        ! Roughness
              Ux=Karman*Uref/IphiM(Lmo(Surf),Zr,ZoM)                    ! Friction velocity
              Ux = max(Ux, 0.1)    ! 09.04.2008 Ilia
              ZoH=ZoM*exp(-Karman*(13.6*Pr**(2./3.)-12.))                ! Energy roughness
            else
              ZoM=ZoLU(Surf,s)                                            ! Momentum roughness
              do it=1, 5
                Ux=max(Karman*Uref/IphiM(Lmo(Surf),Zr-Dfor,ZoM),0.1)    ! Friction velocity
                Lmo(Surf)=Almo*Ux*Ux*Ux
              enddo
              ZoH=0.135*ZoM                                                ! Energy roughness
            endif

            Ra=IphiH(Lmo(Surf),Zr-Dfor,ZoH)/Karman/Ux                    ! Aerodynamic resistance

! Calculation of dry deposition velocities
            Form=DryDepInd(Ind)
            VdFog=0.

!---------------------------------------------------------------------------------
! Particulate form
            if(any(Part(1:Npart)==Form)) then

! Low vegetation is treated as barren surface in case of snow cover height more than 100 mm (snow density - 0.2 g/m3)
              if((SnowH>0.1.and.(SurfGroup(Surf)==Grass.or.SurfGroup(Surf)==Arable))&       ! WRF
                                        &.or.(Tsurf<270..and.SurfGroup(Surf)==Water)) then
                SurfG=Barren
              else
                SurfG=SurfGroup(Surf)
              endif
                selectcase(SurfG)
                  case(Forest)
                    St=Ux*Vg/Ggrav/1.e-3
                    Sc=Nu/DiffW
                    Uh=Ux/Karman*IphiM(Lmo(Surf),Hfor-Dfor,ZoM)
                    if(WetForest.eqv..true.) then
                      Reb=1.
                    else
                      Reb=exp(-Cforest*St**0.25)
                    endif
                    Eb=Sc**(-2./3.)
                    Ein=(0.01*DpW/(DpW+2.e-5)+0.99*DpW/(DpW+2.e-3))
                    Eim=(St/(Aforest+St))**0.5
                    Eff=Bforest*(Eb+Ein+Eim)*Reb
                    Rs=Uh/Eff/Ux/Ux
                    Vd(i,j,t,Ind,Surf,nDay)=1./(Ra+Rs)+Vg

                  case(Grass,Arable)
                    St=Ux*Vg/Ggrav/1.e-3
                    Sc=Nu/DiffW
                    Uh=Ux/Karman*IphiM(Lmo(Surf),Hfor-Dfor,ZoM)
                    if(WetGrass.eqv..true.) then
                      Reb=1.
                    else
                      Reb=exp(-Cgrass*St**0.25)
                    endif
                    Eb=Sc**(-2./3.)
                    Ein=(0.01*DpW/(DpW+2.e-5)+0.99*DpW/(DpW+2.e-3))
                    Eim=(St/(Agrass+St))**0.5
                    Eff=Bgrass*(Eb+Ein+Eim)*Reb
                    if(Lmo(Surf)>=0.) then
                      Rs=Uh/Eff/Ux/Ux
                    else
                      Rs=Uh/Eff/Ux/Ux/(1.+(-Dgrass/Lmo(Surf))**(2./3.))
                    endif
                    Vd(i,j,t,Ind,Surf,nDay)=1./(Ra+Rs)+Vg

                  case(Water)
                    St=Ux*Ux*VgS/Ggrav/Nu
                    Sc=Nu/DiffS
                    Uh=Ux/Karman*IphiM(Lmo(Surf),10.,ZoM)
                    Fbrok=min(1.7e-6*Uh**3.75,1.)
                    Rs1=Karman*Uh/Ux/Ux/(10.**(-3./St)+Sc**(-0.5))
                    Rs2=Rbrok ! 10. s/m
                    Rs=1./((1.-Fbrok)/Rs1+Fbrok/Rs2)    
                    Vd(i,j,t,Ind,Surf,nDay)=(1.+Ra*Vg)*(1.+Rs*VgS)/(Ra+Rs+Ra*Rs*VgS)

                  case(Barren)
                    St=Ux*Ux*Vg/Ggrav/Nu
                    Sc=Nu/DiffW
                    Rs=Karman*Uref/Ux/Ux/(10.**(-3./St)+Sc**(-2./3.))
                    Vd(i,j,t,Ind,Surf,nDay)=1./(Ra+Rs+Ra*Rs*Vg)+Vg

                  case(Urban)
                    St=Ux*Ux*Vg/Ggrav/Nu
                    Sc=Nu/DiffW
                    Rs=Karman*Uref/Ux/Ux/(St*St/(400.+St*St)+Sc**(-2./3.))
                    Vd(i,j,t,Ind,Surf,nDay)=1./(Ra+Rs+Ra*Rs*Vg)+Vg
              endselect
            endif
!---------------------------------------------------------------------------------
! Gaseous oxidized form
            if(any(Oxid(1:Noxid)==Form)) then

              Sc=Nu/DiffG
              Rb=2./Karman/Ux*(Sc/Pr)**(2./3.)
              Rc=Rcmin !10. s/m
              Vd(i,j,t,Ind,Surf,nDay)=1./(Ra+Rb+Rc)
            endif

!---------------------------------------------------------------------------------
! Elemental gaseous form
            if(Form==Hg0) then

                Tcoef=min(max((Tair-273.)/20.,0.),1.)
                selectcase(SurfGroup(Surf))
                  case(Forest)
                    Vd(i,j,t,Ind,Surf,nDay)=Vdcoef(1)*Tcoef

                  case(Grass,Arable)
                    Vd(i,j,t,Ind,Surf,nDay)=Vdcoef(2)*Tcoef

                  case(Urban,Barren,Water)
                    Vd(i,j,t,Ind,Surf,nDay)=Vdcoef(3)
              
                endselect
              endif

!---------------------------------------------------------------------------------
! Aqueous forms
              if(Form==Dis.or.Form==Sulf.or.Form==Chlor) then 
                if(Lw>0.) then !.and.PrecN(i,j)<=0.
                  if(VdFog>0.) then
                    Vd(i,j,t,Ind,Surf,nDay)=VdFog
                  else
                    selectcase(SurfGroup(Surf))
                      case(Forest)
                        St=Ux*VgF/Ggrav/1.e-3
                        Sc=Nu/DiffF
                        Uh=Ux/Karman*IphiM(Lmo(Surf),Hfor-Dfor,ZoM)
                        Reb=1.
                        Eb=Sc**(-2./3.)
                        Ein=(0.01*Dfog/(Dfog+2.e-5)+0.99*Dfog/(Dfog+2.e-3))
                        Eim=(St/(Aforest+St))**0.5
                        Eff=Bforest*(Eb+Ein+Eim)*Reb
                        Rs=Uh/Eff/Ux/Ux
                        VdFog=1./(Ra+Rs)+VgF
                        Vd(i,j,t,Ind,Surf,nDay)=VdFog

                      case(Grass,Arable)
                        St=Ux*VgF/Ggrav/1.e-3
                        Sc=Nu/DiffF
                        Uh=Ux/Karman*IphiM(Lmo(Surf),Hfor-Dfor,ZoM)
                        Reb=1.
                        Eb=Sc**(-2./3.)
                        Ein=(0.01*Dfog/(Dfog+2.e-5)+0.99*Dfog/(Dfog+2.e-3))
                        Eim=(St/(Agrass+St))**0.5
                        Eff=Bgrass*(Eb+Ein+Eim)*Reb
                        if(Lmo(Surf)>=0.) then
                          Rs=Uh/Eff/Ux/Ux
                        else
                          Rs=Uh/Eff/Ux/Ux/(1.+(-Dgrass/Lmo(Surf))**(2./3.))
                        endif
                        VdFog=1./(Ra+Rs)+VgF                        
                        Vd(i,j,t,Ind,Surf,nDay)=VdFog

                      case(Water)
                        St=Ux*Ux*VgF/Ggrav/Nu
                        Sc=Nu/DiffF
                        Uh=Ux/Karman*IphiM(Lmo(Surf),10.,ZoM)
                        Fbrok=min(1.7e-6*Uh**3.75,1.)
                        Rs1=Karman*Uh/Ux/Ux/(10.**(-3./St)+Sc**(-0.5))
                        Rs2=Rbrok ! 10. s/m
                        Rs=1./((1.-Fbrok)/Rs1+Fbrok/Rs2)    
                        VdFog=1./(Ra+Rs+Ra*Rs*VgF)+VgF                        
                        Vd(i,j,t,Ind,Surf,nDay)=VdFog

                      case(Barren)
                        St=Ux*Ux*VgF/Ggrav/Nu
                        Sc=Nu/DiffF
                        Rs=Karman*Uref/Ux/Ux/(10.**(-3./St)+Sc**(-2./3.))
                        VdFog=1./(Ra+Rs+Ra*Rs*VgF)+VgF                        
                        Vd(i,j,t,Ind,Surf,nDay)=VdFog

                      case(Urban)
                        St=Ux*Ux*VgF/Ggrav/Nu
                        Sc=Nu/DiffF
                        Rs=Karman*Uref/Ux/Ux/(St*St/(400.+St*St)+Sc**(-2./3.))
                        VdFog=1./(Ra+Rs+Ra*Rs*VgF)+VgF                        
                        Vd(i,j,t,Ind,Surf,nDay)=VdFog
                      endselect
                  endif
                else
                  Vd(i,j,t,Ind,Surf,nDay)=0.
                endif
              endif
          enddo            ! end Surf
        enddo            ! end i
      enddo                ! end j
    enddo                ! end t


contains
!...............................................................................
real function IphiM(L,Zr,Zo)

    real L, Zr, Zo
    real(8) Kr, Ko

    if(L==0.) L=1.e-6

    if(Zr/L<1.e-6) then                     ! Neutral    
      IphiM=log(Zr/Zo)
    elseif(L>=0.) then                        ! Stable
      IphiM=log(Zr/Zo)+Bm/L*(Zr-Zo)
    else                                    ! Unstable
      Kr=(1.-dble(Gm*Zr/L))**0.25
      Ko=(1.-dble(Gm*Zo/L))**0.25
      IphiM=dlog((Kr-1.)/(Kr+1.))-dlog((Ko-1.)/(Ko+1.))+2*datan(Kr)-2.*datan(Ko)
    endif

end function IphiM
!...............................................................................
real function IphiH(L,Zr,Zo)

    real L, Zr, Zo
    real(8) Kr, Ko

    if(L==0.) L=1.e-6

    if(Zr/L<1.e-6) then                     ! Neutral    
      IphiH=Prt*log(Zr/Zo)
    elseif(L>=0.) then                        ! Stable
      IphiH=Prt*log(Zr/Zo)+Bh/L*(Zr-Zo)
    else                                    ! Unstable
      Kr=sqrt(1.-dble(Gh*Zr/L))
      Ko=sqrt(1.-dble(Gh*Zo/L))
      IphiH=Prt*(dlog((Kr-1.)/(Kr+1.))-dlog((Ko-1.)/(Ko+1.)))
    endif

end function IphiH
!...............................................................................
real function WetDiameter(Dd,S,e)

    real Dd, S, e
    real A, B, k, k1, k2

    A=1.2*exp(0.066*S/(1.058-S))
    B=exp(0.00077*S/(1.009-S))
    k1=10.2-23.7*S+14.5*S*S
    k2=-6.7+15.5*S-9.2*S*S
    k=1.-k1*(1.-e)-k2*(1.-e*e)

    WetDiameter=A*(Dd*1.e6)**B*k*1.e-6   

end function WetDiameter
!...............................................................................
end subroutine Atm_Hg_ExchPar

#endif

end module Atm_Hg_General
