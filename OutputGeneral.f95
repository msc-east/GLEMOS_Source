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
! Module of the model output
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


module GeneralOutput

  use GeneralParams
  use Geometry
  use Atm_Params
#ifdef M_SOIL
use Soil_Params
#endif
#ifdef M_OCN
  use Ocn_Params
#endif
#ifdef M_VEG
  use Veg_Params
#endif
  use Exch_Params
  use Balance

  character(800), private :: fileName, fullName
  character(4), private :: YearNum
  
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine creating and writing records to the logfile
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LogFile(recType,res,note)

    character(*), intent(in)            :: recType
    integer, intent(in), optional       :: res
    character(*), intent(in), optional  :: note

    character(12)   :: resType(0:1)=(/'..successful','unsuccessful'/)
    integer         ::  lenField=60
    integer             numP, lenRec
    integer(8)          curVal(8)
    character(600)      logRec
    integer             fs

    fullName=trim(OutPath)//trim(LogName)                                                      

    selectcase(recType)
        case('Initial')
            open(101, file=fullName, action='write', status='replace', access='sequential', iostat=fs)
            if(fs > 0) then
              print '(/,"STOP: Cannot create log file ''",a,"''",/)', trim(fullName)
              stop
            endif
            write(101,'(a)') '@@@@@@  Log file of the calculation run  '//repeat('@',lenField-41)
            call date_and_time(values=curVal)
            write(101,'(/,"Pollutant: ",a)') SubsID(1)
            write(101,'(/,"Started: ",i2,"-",a3,"-",i4," at ",i2.2,":",i2.2,":",i2.2)') curVal(3), MonthName(curVal(2)),&
                            &curVal(1), curVal(5), curVal(6), curVal(7)
            write(101,'(/,a)') repeat('*',lenField)
            close(101)

       case('Monthly')
            open(101, file=fullName, action='write', status='old', access='append')
            write(YearNum,'(i4)') Year
            logRec=MonthName(Month)//' '//YearNum
            write(101,'(a)') logRec
            close(101)

       case('Final')
            call OutputBalanceInLogFile(trim(fullName))
    
       case default
            open(101, file=fullName, action='write', status='old', access='append')
            if(present(res)) then
                lenRec=lenField-13-len(recType)
                numP=max(0,lenRec)
                logRec=recType//':'//repeat('.',numP)//resType(res)
            else
                logRec=recType
            endif
            if(present(note)) then
        write(101,'(/,a,/,a,/)') trim(logRec), note
            else
        write(101,'(/,a)') trim(logRec)
            endif
            close(101)

    endselect

end subroutine LogFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine printing information to the screen
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PrintScr(Prd)

    character(*) Prd
    character(2) DayNum
    character(4) YearNum

    select case(Prd)

        case('Daily')
            write(YearNum,'(i4)') Year
            write(DayNum,'(i2.2)') Day
            print '(a2," ",a3,", ",a4)', DayNum, MonthName(Month), YearNum

    case('Final')
            call PrintBalanceScr
    end select
    
end subroutine PrintScr


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine writing current dumping information to the dump file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WriteDump(yr, mn, dn)

    character*8 DumpDate
    integer     FileStat, yr, mn, dn

    write(DumpDate,'(i4,i2.2,i2.2)') yr, mn, dn
    fileName=trim(DumpName)//trim(SubsID(1))//'_'//DumpDate//'.bin'
    fullName=trim(adjustL(OutPath))//trim(DumpDir)//'/'//trim(fileName)
    open(1, file=fullName, form='unformatted', access="stream", action='write', iostat=FileStat)
    if(FileStat>0) then
      print '(/,"STOP: Cannot open dump file for writing: ''",a,"''",/)', trim(fullName)
      stop
    endif
    
! Writing dump version
    write(1) DumpVersion
    
! Dumping grid code 
    write(1) GridCode
    
! Dumping date and time 
    write(1) yr, mn, dn, timeCalc
    
! Run characteristics
    write(1) NumForm(Atm), NumSurf, NumAnth, NumNat, NumRe, NumSrc, NumRcp
    
! Content in media for initial conditions
    write(1) Atm_MixRatio(bImin:bImax,bJmin:bJmax,0:Atm_Kmax+1,1:NumForm(Atm))
    write(1) Atm_Conc(bImin:bImax,bJmin:bJmax,0:Atm_Kmax+1,1:NumForm(Atm))
#ifdef M_SOIL
    write(1) Soil_Conc(IMin:IMax,JMin:JMax,1:Soil_KMax,1:Soil_NumType,1:NumForm(Soil))
#endif
#ifdef M_OCN
    write(1) Ocn_conc(Imin:Imax,Jmin:Jmax,1:Ocn_Kmax,1:NumForm(Ocn),1:3)
    write(1) Water_upl_conc(IMin:IMax,JMin:JMax,1:NumForm(Ocn))
#endif
#ifdef M_VEG
    write(1) Veg_Conc(IMin:IMax,JMin:JMax,1:Veg_NumType,1:NumForm(Veg))
    write(1) Fall_Conc(IMin:IMax,JMin:JMax,1:Veg_NumType,1:NumSubsMedia(Veg))
#endif
#if RTYPE==2
    write(1) SourcID
    write(1) Atm_Contrib(Imin:Imax,bJmin:bJmax,0:Atm_Kmax+1,1:NumForm(Atm),1:NumSrc)
#endif

! Atmospheric concentrations and fluxes
    write(1) Atm_MixYear(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),1:NumSrc)
    write(1) Atm_ConcYear(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),1:NumSrc)
    write(1) DryDepYear(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
    write(1) WetDepYear(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
    write(1) PrecipYear(Imin:Imax,bJmin:bJmax)

    write(1) Atm_MixMonth(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),1:NumSrc)
    write(1) Atm_ConcMonth(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),1:NumSrc)
    write(1) DryDepMonth(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
    write(1) WetDepMonth(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
    write(1) PrecipMonth(Imin:Imax,bJmin:bJmax)

    write(1) AntEmisYear(Imin:Imax,Jmin:Jmax)
    write(1) NatEmisYear(Imin:Imax,Jmin:Jmax)
    write(1) ReEmisYear(Imin:Imax,Jmin:Jmax)
    write(1) PxYear(bImin:bImax,bJmin:bJmax)
    
    write(1) AntEmisMonth(Imin:Imax,Jmin:Jmax)
    write(1) NatEmisMonth(Imin:Imax,Jmin:Jmax)
    write(1) ReEmisMonth(Imin:Imax,Jmin:Jmax)
    write(1) PxMonth(bImin:bImax,bJmin:bJmax)
    
#ifdef G_HG
    write(1) FreshDep(Imin:Imax,bJmin:bJmax,1:NumSurf,1:NumSrc)
    write(1) ReEmisFlux(Imin:Imax,bJmin:bJmax,1:ReEmisNum,1:NumSrc)
#endif
    write(1) DryRemYear(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
    write(1) DryRemMonth(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)

! Balance masses
    write(1) MassAtmInit, MassAtmAntEmis, MassAtmNatEmis, MassReEmis, MassAtmFin, MassAtmUp,&
        &MassAtmBnd, MassDryDep, MassDryRem, MassWetDep, MassChemEx

! Atmosphere
#ifdef G_POP
    write(1) MassAtmDegr, AirOcnFlux, AirSoilFlux, AirVegFlux, MassAtmPartit
#endif

! Soil compartment
#ifdef M_SOIL
    write(1) SoilInit, MassSoilEmis
    write(1) Soil_ConcYear(IMin:IMax,JMin:JMax,1:Soil_KMax,1:Soil_NumType,1:NumSubs)
    write(1) Soil_ConcMonth(IMin:IMax,JMin:JMax,1:Soil_KMax,1:Soil_NumType,1:NumSubs)

#ifdef G_POP
    write(1) Soil_Degr,SoilAirFlux
    write(1) QPrec(IMin:IMax,JMin:JMax), DeltaQPrec(IMin:IMax,JMin:JMax), Inverse(IMin:IMax,JMin:JMax)
#endif
#endif

! Ocean compartment
#ifdef M_OCN
    write(1) Ocn_Sedim, oldtime, curtime, newtime, OcnInit
    write(1) Water_upl_concYear(IMin:IMax,JMin:JMax,1:NumSubs)
    write(1) Ocn_ConcYear(IMin:IMax,JMin:JMax,1:Ocn_KMax,1:NumSubs)
    write(1) Water_upl_concMonth(IMin:IMax,JMin:JMax,1:NumSubs)
    write(1) Ocn_ConcMonth(IMin:IMax,JMin:JMax,1:Ocn_KMax,1:NumSubs)
#ifdef G_POP
    write(1) Ocn_Degr,OcnAirFlux
#endif
#endif

! Vegetation compartment
#ifdef M_VEG
    write(1) VegInit, FallInit
    write(1) Veg_ConcYear(IMin:IMax,JMin:JMax,Veg_NumType,NumForm(Veg))
    write(1) Fall_ConcYear(IMin:IMax,JMin:JMax,Veg_NumType,NumSubsMedia(Veg))
    write(1) Veg_ConcMonth(IMin:IMax,JMin:JMax,Veg_NumType,NumForm(Veg))
    write(1) Fall_ConcMonth(IMin:IMax,JMin:JMax,Veg_NumType,NumSubsMedia(Veg))
#ifdef G_POP
    write(1) Veg_Degr,VegAirFlux, VegSoilFlux
#endif
#endif
    close(1)

end subroutine WriteDump


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading the dump file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadDump

    integer FormNum, SurfNum, AnthNum, NatNum, ReNum, SrcNum, RcpNum, FileStat
    integer DumpYear, DumpMonth, DumpDay
    character(2) DumpFreq
    character(4) YearNum
    character(8) IDsourc(MaxMatr)
    character(10) DumpGrid, InputDumpVersion


    fullName=trim(InitPath)//trim(InitName)

    open(1, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open dump file for reading: ''",a,"''",/)', trim(fullName)
      stop
    endif

! Reading dump version
    read(1) InputDumpVersion
    if(trim(InputDumpVersion)/=trim(DumpVersion)) then
      print '(/,"STOP: Incorrect dump version ''",a,"'' in file ",a,/)', trim(InputDumpVersion), trim(fullName)
      stop
    endif
    
! Reading grid code
    read(1) DumpGrid
    if(trim(DumpGrid)/=trim(GridCode)) then
      print '(/,"STOP: Incorrect grid code ''",a,"'' in file ",a,/)', trim(DumpGrid), trim(fullName)
      stop
    endif

! Dumping date 
    read(1) DumpYear, DumpMonth, DumpDay, timeCalc
    print *, 'Read dump: ',DumpYear, DumpMonth, DumpDay

    if(DumpDay<MonthDays(DumpMonth)) then
        BegDate(da)=DumpDay+1
        BegDate(mn)=DumpMonth
        BegDate(yr)=DumpYear
    elseif(DumpMonth<12) then
        BegDate(da)=1
        BegDate(mn)=DumpMonth+1
        BegDate(yr)=DumpYear
    else
        BegDate(da)=1
        BegDate(mn)=1
        BegDate(yr)=DumpYear+1
    endif

! Run characteristics
    read(1) FormNum, SurfNum, AnthNum, NatNum, ReNum, SrcNum, RcpNum

    if(FormNum/=NumForm(Atm)) then 
      print '(/,"STOP: Wrong number of the pollutant components&
        & in the dump file ''",a,"''",/)', trim(fileName)
      stop
    endif
    if(SurfNum/=NumSurf) then 
      print '(/,"STOP: Wrong number of the surface types&
        & in the dump file ''",a,"''",/)', trim(fileName)
      stop
    endif
    if(AnthNum/=NumAnth) then 
      print '(/,"STOP: Wrong number of anthropogenic sources&
        & in the dump file ''",a,"''",/)', trim(fileName)
      stop
    endif
    if(NatNum/=NumNat) then 
      print '(/,"STOP: Wrong number of anthropogenic sources&
        & in the dump file ''",a,"''",/)', trim(fileName)
      stop
    endif
    if(ReNum/=NumRe) then 
      print '(/,"STOP: Wrong number of anthropogenic sources&
        & in the dump file ''",a,"''",/)', trim(fileName)
      stop
    endif
    if(SrcNum/=NumSrc) then 
      print '(/,"STOP: Wrong number of receptors&
        & in the dump file ''",a,"''",/)', trim(fileName)
      stop
    endif
    if(RcpNum/=NumRcp) then 
      print '(/,"STOP: Wrong number of receptors&
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
#if RTYPE==2
    read(1) IDsourc
    read(1) Atm_Contrib(Imin:Imax,bJmin:bJmax,0:Atm_Kmax+1,1:NumForm(Atm),1:NumSrc)
#endif

! Atmospheric concentrations and fluxes
    read(1) Atm_MixYear(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),1:NumSrc)
    read(1) Atm_ConcYear(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),1:NumSrc)
    read(1) DryDepYear(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
    read(1) WetDepYear(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
    read(1) PrecipYear(Imin:Imax,bJmin:bJmax)

    read(1) Atm_MixMonth(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),1:NumSrc)
    read(1) Atm_ConcMonth(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),1:NumSrc)
    read(1) DryDepMonth(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
    read(1) WetDepMonth(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
    read(1) PrecipMonth(Imin:Imax,bJmin:bJmax)

!    read(1) Atm_MixDay(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),1:NumSrc)
!    read(1) Atm_ConcDay(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1:NumForm(Atm),1:NumSrc)
!    read(1) DryDepDay(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
!    read(1) WetDepDay(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
!    read(1) PrecipDay(Imin:Imax,bJmin:bJmax)
    
    read(1) AntEmisYear(Imin:Imax,Jmin:Jmax)
    read(1) NatEmisYear(Imin:Imax,Jmin:Jmax)
    read(1) ReEmisYear(Imin:Imax,Jmin:Jmax)
    read(1) PxYear(bImin:bImax,bJmin:bJmax)
    
    read(1) AntEmisMonth(Imin:Imax,Jmin:Jmax)
    read(1) NatEmisMonth(Imin:Imax,Jmin:Jmax)
    read(1) ReEmisMonth(Imin:Imax,Jmin:Jmax)
    read(1) PxMonth(bImin:bImax,bJmin:bJmax)

#ifdef G_HG
    read(1) FreshDep(Imin:Imax,bJmin:bJmax,1:NumSurf,1:NumSrc)
    read(1) ReEmisFlux(Imin:Imax,bJmin:bJmax,1:ReEmisNum,1:NumSrc)
#endif
    read(1) DryRemYear(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
    read(1) DryRemMonth(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)
!    read(1) DryRemDay(Imin:Imax,bJmin:bJmax,1:NumForm(Atm),1:NumSurf,1:NumSrc)

! Balance masses
    read(1) MassAtmInit, MassAtmAntEmis, MassAtmNatEmis, MassReEmis, MassAtmFin, MassAtmUp,&
        &MassAtmBnd, MassDryDep, MassDryRem, MassWetDep, MassChemEx

! Atmosphere
#ifdef G_POP
    read(1) MassAtmDegr, AirOcnFlux, AirSoilFlux, AirVegFlux, MassAtmPartit
#endif
! Soil compartment
#ifdef M_SOIL
    read(1) SoilInit, MassSoilEmis
    read(1) Soil_ConcYear(IMin:IMax,JMin:JMax,1:Soil_KMax,1:Soil_NumType,1:NumSubs)
    read(1) Soil_ConcMonth(IMin:IMax,JMin:JMax,1:Soil_KMax,1:Soil_NumType,1:NumSubs)
#ifdef G_POP
    read(1) Soil_Degr,SoilAirFlux
    read(1) QPrec(IMin:IMax,JMin:JMax), DeltaQPrec(IMin:IMax,JMin:JMax), Inverse(IMin:IMax,JMin:JMax)
#endif
#endif
! Ocean compartment
#ifdef M_OCN
    read(1) Ocn_Sedim, oldtime, curtime, newtime, OcnInit
    read(1) Water_upl_concYear(IMin:IMax,JMin:JMax,1:NumSubs)
    read(1) Ocn_ConcYear(IMin:IMax,JMin:JMax,1:Ocn_KMax,1:NumSubs)
    read(1) Water_upl_concMonth(IMin:IMax,JMin:JMax,1:NumSubs)
    read(1) Ocn_ConcMonth(IMin:IMax,JMin:JMax,1:Ocn_KMax,1:NumSubs)
#ifdef G_POP
    read(1) Ocn_Degr,OcnAirFlux
#endif
#endif
! Vegetation compartment
#ifdef M_VEG
    read(1) VegInit, FallInit
    read(1) Veg_ConcYear(IMin:IMax,JMin:JMax,Veg_NumType,NumForm(Veg))
    read(1) Fall_ConcYear(IMin:IMax,JMin:JMax,Veg_NumType,NumSubsMedia(Veg))
    read(1) Veg_ConcMonth(IMin:IMax,JMin:JMax,Veg_NumType,NumForm(Veg))
    read(1) Fall_ConcMonth(IMin:IMax,JMin:JMax,Veg_NumType,NumSubsMedia(Veg))
#ifdef G_POP
    read(1) Veg_Degr,VegAirFlux, VegSoilFlux
#endif
#endif
    close(1)

    if(BegDate(da)==1.and.BegDate(mn)==1.and.BegDate(yr)==DumpYear+1) then
      Atm_MixYear=0.
      Atm_ConcYear=0.
      PrecipYear=0.
      DryDepYear=0.
      WetDepYear=0.
      DryRemYear=0.
      Ocn_ConcYear=0.
      Water_upl_concYear=0.
      Soil_ConcYear=0.
      Veg_ConcYear=0.
      Fall_ConcYear=0.
      AntEmisYear=0.
      NatEmisYear=0.
      ReEmisYear=0.
      PxYear=0.
      timeCalc=0.
    endif    

end subroutine ReadDump


end module GeneralOutput

