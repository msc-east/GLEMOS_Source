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
module Atm_Output

    use GeneralParams
    use Geometry
    use Atm_Params
    use Exch_Params
#if G_POPO
	use Soil_Params                ! 18-10-2019
#endif
    use GeneralOutput
    use TextOutputProc
    use Atm_NcfOutput
    use Balance

    implicit none

    integer, parameter      :: FIELDS_NUM    = 10
    integer, parameter      :: FIELDS_BY_FORM = 7
    integer, parameter      :: MATRIX_NUM    = 4
    integer, parameter      :: MONIT_NUM     = 3
    integer, parameter      :: OUT_VARS_NUM  = FIELDS_NUM + MATRIX_NUM + MONIT_NUM
    real, parameter         :: stp = 1.e5/Rair/273.15      ! Coef to convert mix ratio to conc at standard temperature and pressure

    character(80), private  :: hdr(5,OUT_VARS_NUM)
    character(120), private :: fName(5,OUT_VARS_NUM)                            
    real, private, allocatable :: fld_2d(:,:,:,:)
    logical, private        :: fld_out(5,OUT_VARS_NUM)
    character(120), private :: filePath(OUT_VARS_NUM)
    character(80), private  :: fld_unt(OUT_VARS_NUM)
    character(80), private  :: fld_nm(OUT_VARS_NUM)

    character(800), private  :: fileName, fullName, fp_prefix
    character(40), private  :: fp_suffix, tm_pfx(5), unt_sfx(5)
    character(4), private   :: YearNum
    character(2), private   :: MonthNum, DayNum, PerNum
    character(10), private  :: DateNum

    data tm_pfx / 'Hourly','Six hourly','Daily','Monthly','Annual'/
    data unt_sfx / '/hour','/6hour','/day','/month','/year' /
    data fld_nm / &
                    &'air concentrations',&
                    &'concentrations in precipitation',&
                    &'dry deposition flux',&
                    &'wet deposition flux',&
                    &'total deposition flux',&
                    &'net deposition flux',&
                    &'re-emission flux',&
                    &'anthropogenic emissions',&
                    &'natural emissions',&
                    &'precipitation amount',&
                    &'air concentrations matrix',&
                    &'wet deposition matrix',&
                    &'dry deposition matrix',&
                    &'total deposition matrix',&
                    &'air concentrations',&
                    &'concentrations in precipitation',&
                    &'wet deposition'&
                    & /
    data fld_unt / 'ng/m3','ng/L','g/km2','g/km2','g/km2','g/km2','g/km2','g/km2','g/km2','mm',&
                    &'ng/m3','kg','kg','kg','ng/m3','ng/L','g/km2' /

    integer, parameter      :: AC_IND    = 1
    integer, parameter      :: CP_IND    = 2
    integer, parameter      :: DD_IND    = 3
    integer, parameter      :: WD_IND    = 4
    integer, parameter      :: TD_IND    = 5
    integer, parameter      :: ND_IND    = 6
    integer, parameter      :: RE_IND    = 7
    integer, parameter      :: AE_IND    = 8
    integer, parameter      :: NE_IND    = 9
    integer, parameter      :: PA_IND    = 10
    integer, parameter      :: ACM_IND   = 11
    integer, parameter      :: WDM_IND   = 12
    integer, parameter      :: DDM_IND   = 13
    integer, parameter      :: TDM_IND   = 14
    integer, parameter      :: AC_MON    = 15
    integer, parameter      :: CP_MON    = 16
    integer, parameter      :: WD_MON    = 17

    integer :: numHour=1


contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine  of the calculation run
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_OutputData(Prd)

    character(*), intent(in)    :: Prd
    integer         i,j,k, Src,n,L
    real            dT, t

    dT=Tstep(Atm)

    selectcase(Prd)

        case('Initial')

! Set up information for output of fields, matrix output, and monitor files
            call Atm_Output_Initial                                                         
            call Atm_WriteNCF_Initial                                                       

!---------------------------------------------------------------------------------
        case('Steply')
            timeCalc=timeCalc+dT
            call ConcToMixRatio

!***************** Matrix calculations *****************
#if RTYPE==2
            do n = 1, NumForm(Atm)              
                do Src=1, NumSrc
                    do j=Jmin, Jmax
                        Atm_ConcDay(minI(j):maxI(j),j,1:Atm_Kmax,n,Src) = &               
                            &Atm_ConcDay(minI(j):maxI(j),j,1:Atm_Kmax,n,Src)&             
                            & + Atm_Conc(minI(j):maxI(j),j,1:Atm_Kmax,n)*dT&              
                            &*Atm_Contrib(minI(j):maxI(j),j,1:Atm_Kmax,n,Src)             

                        Atm_MixDay(minI(j):maxI(j),j,1:Atm_Kmax,n,Src) = &                
                            &Atm_MixDay(minI(j):maxI(j),j,1:Atm_Kmax,n,Src)&              
                            & + Atm_MixRatio(minI(j):maxI(j),j,1:Atm_Kmax,n)*dT&          
                            &*Atm_Contrib(minI(j):maxI(j),j,1:Atm_Kmax,n,Src)           
                    enddo ! j
                enddo ! Src
            end do ! n
#else
            do n = 1, NumForm(Atm)              
                do j=Jmin, Jmax
                    Atm_ConcDay(minI(j):maxI(j),j,1:Atm_Kmax,n,1) = &                     
                        &Atm_ConcDay(minI(j):maxI(j),j,1:Atm_Kmax,n,1)&                   
                        & + Atm_Conc(minI(j):maxI(j),j,1:Atm_Kmax,n)*dT                   

                    Atm_MixDay(minI(j):maxI(j),j,1:Atm_Kmax,n,1) = &                      
                        &Atm_MixDay(minI(j):maxI(j),j,1:Atm_Kmax,n,1)&                    
                        & + Atm_MixRatio(minI(j):maxI(j),j,1:Atm_Kmax,n)*dT               
                enddo
            end do ! n
#endif

            t = abs(DayTime-anint(real(numHour)*(dTinput/dT)/6.)*dT)
            if(t <= 10.) then
                if(AtmOutMonitHourly) call Atm_WriteMonitor('Hourly')
                if(AtmOutFieldsHourly) call Atm_WriteFields('Hourly')
                numHour=numHour+1
            endif

!---------------------------------------------------------------------------------
	case('6hourly')

            if(AtmOutMonit6hourly)     call Atm_WriteMonitor('6hourly')
            if(AtmOutFields6hourly)    call Atm_WriteFields('6hourly')
            if(AtmOutNCF6hourly)       call Atm_WriteFieldsNCF('6hourly',numHour-1)              
        ! Accumulate daily surface pressure
            PxDay(:,:) = PxDay(:,:) + Px(:,:,Period,1)*SecInDay/NumPer                           

!---------------------------------------------------------------------------------
	case('Daily')

#if G_POP
            call DepSurf2LU       ! 18-10-2019
#endif
            Period=Period-1             ! Set Period equal to 4 (last value during the day)
            if(AtmOutMonitDaily)    call Atm_WriteMonitor('Daily')
            if(AtmOutFieldsDaily)   call Atm_WriteFields('Daily')
            if(AtmOutNCFDaily)      call Atm_WriteFieldsNCF('Daily')

#if RTYPE==2
            if(AtmOutMatrixDaily)   call Atm_WriteMatrix('Daily')
            if(AtmOutNCFDaily)      call Atm_WriteMonitoringNCF('Daily')
#endif

! Update counters atmospheric concentrations and mixing ratio
            do n = 1, NumForm(Atm)              
                do j=Jmin, Jmax
                    Atm_ConcMonth(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)=&          
                        &Atm_ConcMonth(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)+&     
                        &Atm_ConcDay(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)         

                    Atm_MixMonth(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)=&           
                        &Atm_MixMonth(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)+&      
                        &Atm_MixDay(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)          
                enddo ! j
            end do ! n
! Update counters of dry and wet depositions
            DryDepMonth=DryDepMonth+DryDepDay
            DryDepDay=0.
            DryDepDayTmp=0.        ! 16-10-2019
!#ifdef G_POP
            DryRemMonth=DryRemMonth+DryRemDay
            DryRemDay=0.
            DryRemDayTmp=0.        ! 16-10-2019
!#endif
            WetDepMonth=WetDepMonth+WetDepDay
            WetDepDay=0.
            WetDepDayTmp=0.        ! 16-1002-19
            Atm_ConcDay=0.
            Atm_MixDay=0.
            PrecipMonth=PrecipMonth+PrecipDay
            PrecipDay=0.
        ! Accumulate monthly surface pressure
            PxMonth = PxMonth + PxDay                                           
            PxDay = 0.                                                          
            numHour=1

            if(OutputDumpDaily) call WriteDump(Year, Month, Day)

!---------------------------------------------------------------------------------

	case('Monthly')
            if(AtmOutMonitMonthly)      call Atm_WriteMonitor('Monthly')
            if(AtmOutFieldsMonthly)     call Atm_WriteFields('Monthly')
            if(AtmOutNCFMonthly)        call Atm_WriteFieldsNCF('Monthly')

#if RTYPE==2
            if(AtmOutMatrixMonthly)     call Atm_WriteMatrix('Monthly')
            if(AtmOutNCFMonthly)        call Atm_WriteMonitoringNCF('Monthly')
#endif

            do n = 1, NumForm(Atm)              
                do j=Jmin, Jmax
                    Atm_ConcYear(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)=&               
                        &Atm_ConcYear(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)+&          
                        &Atm_ConcMonth(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)           

                    Atm_MixYear(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)=&                
                        &Atm_MixYear(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)+&           
                        &Atm_MixMonth(minI(j):maxI(j),j,1:Atm_Kmax,n,1:NumSrc)            
                enddo ! j
            end do ! n

            Atm_ConcMonth=0.
            Atm_MixMonth=0.
            DryDepYear=DryDepYear+DryDepMonth
            DryDepMonth=0.
!#ifdef G_POP
            DryRemYear=DryRemYear+DryRemMonth
            DryRemMonth=0.
!#endif
            WetDepYear=WetDepYear+WetDepMonth
            WetDepMonth=0.
            PrecipYear=PrecipYear+PrecipMonth
            PrecipMonth=0.
            AntEmisYear = AntEmisYear + AntEmisMonth          
            AntEmisMonth = 0.                                 
            NatEmisYear = NatEmisYear + NatEmisMonth          
            NatEmisMonth = 0.                                 
        ! Accumulate yearly surface pressure
            PxYear = PxYear + PxMonth                         
            PxMonth = 0.                                      

! Write dump only at the end of a month.
            if(Day > MonthDays(Month) .and. OutputDumpMonthly) call WriteDump(Year, Month, Day-1)

!---------------------------------------------------------------------------------
	case('Yearly')
            if(AtmOutMonitYearly)       call Atm_WriteMonitor('Yearly')
            if(AtmOutFieldsYearly)      call Atm_WriteFields('Yearly')
            if(AtmOutNCFYearly)         call Atm_WriteFieldsNCF('Yearly')

#if RTYPE==2
            if(AtmOutMatrixYearly)      call Atm_WriteMatrix('Yearly')
            if(AtmOutNCFYearly)         call Atm_WriteMonitoringNCF('Yearly')
#endif
! Write dump only at the end of a year.
            if(Month > 12 .and. OutputDumpYearly) call WriteDump(Year, Month, Day-1)

            Atm_ConcYear=0.
            Atm_MixYear = 0.
            PrecipYear=0.
            DryDepYear=0.
!#ifdef G_POP
            DryRemYear=0.
!#endif
            WetDepYear=0.
            AntEmisYear=0.
            NatEmisYear=0.
            ReEmisYear=0.
            PxYear = 0.

            timeCalc=0.
            

!---------------------------------------------------------------------------------
	case('Final')
            call Atm_Output_Final                               
            call Atm_WriteNCF_Final                             
!---------------------------------------------------------------------------------
	endselect

end subroutine Atm_OutputData
!*********************************************************************************


subroutine Atm_Output_Initial
    integer     aErr, s, itm
    character(40)   tmDir(5)

    allocate(fld_2d(Imin:Imax,Jmin:Jmax,MaxForm,OUT_VARS_NUM), stat=aErr)
    if(aErr/=0) stop 'STOP: Memory allocation error in Atm_Output_Initial'

    tmDir(1) = trim(HourlyDir)
    tmDir(2) = trim(SixHourlyDir)
    tmDir(3) = trim(DailyDir)
    tmDir(4) = trim(MonthlyDir)
    tmDir(5) = trim(YearlyDir)

    fld_out(:,AC_IND) = (/.true.,   .true., .true., .true., .true./)
    fld_out(:,CP_IND) = (/.false., .false., .true., .true., .true./)
    fld_out(:,DD_IND) = (/.false., .false., .true., .true., .true./)
    fld_out(:,WD_IND) = (/.false., .false., .true., .true., .true./)
    fld_out(:,TD_IND) = (/.false., .false., .true., .true., .true./)
    fld_out(:,AE_IND) = (/.false., .false.,.false., .true., .true./)
    fld_out(:,NE_IND) = (/.false., .false.,.false., .true., .true./)
    fld_out(:,ND_IND) = (/.false., .false., .true., .true., .true./)
    fld_out(:,RE_IND) = (/.false., .false., .true., .true., .true./)
    fld_out(:,PA_IND) = (/.false., .false., .true., .true., .true./)
    fld_out(:,ACM_IND)= (/.false., .false., .true., .true., .true./)
    fld_out(:,WDM_IND)= (/.false., .false., .true., .true., .true./)
    fld_out(:,DDM_IND)= (/.false., .false., .true., .true., .true./)
    fld_out(:,TDM_IND)= (/.false., .false., .true., .true., .true./)
    fld_out(:,AC_MON) = (/.false., .false., .true., .true., .true./)
    fld_out(:,CP_MON) = (/.false., .false., .true., .true., .true./)
    fld_out(:,WD_MON) = (/.false., .false., .true., .true., .true./)

! Fill in parts of output file names and headers
    do s = 1, NumSubsMedia(Atm)
        do itm=1, 5
        ! File names for output variable
            fp_prefix = trim(OutPath)//trim(AtmDir)//trim(FieldsDir)//trim(tmDir(itm))//trim(SubsID(s))
            fName(itm,AC_IND) = trim(fp_prefix)//trim(ConcName)
            fName(itm,CP_IND) = trim(fp_prefix)//trim(InPrecName)
            fName(itm,DD_IND) = trim(fp_prefix)//trim(DryName)
            fName(itm,WD_IND) = trim(fp_prefix)//trim(WetName)
            fName(itm,TD_IND) = trim(fp_prefix)//trim(TotName)
            fName(itm,ND_IND) = trim(fp_prefix)//trim(NetName)
            fName(itm,RE_IND) = trim(fp_prefix)//trim(RemName)
            fName(itm,AE_IND) = trim(fp_prefix)//trim(EmisOutAnt)
            fName(itm,NE_IND) = trim(fp_prefix)//trim(EmisOutNat)

            fp_prefix = trim(OutPath)//trim(AtmDir)//trim(FieldsDir)//trim(tmDir(itm))
            fName(itm,PA_IND) = trim(fp_prefix)//trim(PrecipName)

            fp_prefix = trim(OutPath)//trim(AtmDir)//trim(MatrixDir)//trim(tmDir(itm))//trim(SubsID(s))
            fName(itm,ACM_IND) = trim(fp_prefix)//trim(MatrixConcName)
            fName(itm,WDM_IND) = trim(fp_prefix)//trim(MatrixWetDepName)
            fName(itm,DDM_IND) = trim(fp_prefix)//trim(MatrixDryDepName)
            fName(itm,TDM_IND) = trim(fp_prefix)//trim(MatrixTotDepName)

            fName(itm,AC_MON) = trim(ConcMonitName)
            fName(itm,CP_MON) = trim(InprecMonitName)
            fName(itm,WD_MON) = trim(WetMonitName)

            
        ! Text headers for output variable
            hdr(itm,AC_IND)= trim(tm_pfx(itm))//' '//trim(fld_nm(AC_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(AC_IND))
            hdr(itm,CP_IND)= trim(tm_pfx(itm))//' '//trim(fld_nm(CP_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(CP_IND))
            hdr(itm,DD_IND)= trim(tm_pfx(itm))//' '//trim(fld_nm(DD_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(DD_IND))&
                                &//trim(unt_sfx(itm))
            hdr(itm,WD_IND)= trim(tm_pfx(itm))//' '//trim(fld_nm(WD_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(WD_IND))&
                                &//trim(unt_sfx(itm))
            hdr(itm,TD_IND)= trim(tm_pfx(itm))//' '//trim(fld_nm(TD_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(TD_IND))&
                                &//trim(unt_sfx(itm))
            hdr(itm,ND_IND)= trim(tm_pfx(itm))//' '//trim(fld_nm(ND_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(ND_IND))&
                                &//trim(unt_sfx(itm))
            hdr(itm,RE_IND)= trim(tm_pfx(itm))//' '//trim(fld_nm(RE_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(RE_IND))&
                                &//trim(unt_sfx(itm))
            hdr(itm,AE_IND)= trim(tm_pfx(itm))//' '//trim(fld_nm(AE_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(AE_IND))&
                                &//trim(unt_sfx(itm))
            hdr(itm,NE_IND)= trim(tm_pfx(itm))//' '//trim(fld_nm(NE_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(NE_IND))&
                                &//trim(unt_sfx(itm))
            hdr(itm,PA_IND)=trim(tm_pfx(itm))//' '//trim(fld_nm(PA_IND))//' '//trim(fld_unt(PA_IND))//trim(unt_sfx(itm))
            hdr(itm,ACM_IND)=trim(tm_pfx(itm))//' '//trim(fld_nm(ACM_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(ACM_IND))
            hdr(itm,WDM_IND)=trim(tm_pfx(itm))//' '//trim(fld_nm(WDM_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(WDM_IND))&
                                &//trim(unt_sfx(itm))
            hdr(itm,DDM_IND)=trim(tm_pfx(itm))//' '//trim(fld_nm(DDM_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(DDM_IND))&
                                &//trim(unt_sfx(itm))
            hdr(itm,TDM_IND)=trim(tm_pfx(itm))//' '//trim(fld_nm(TDM_IND))//' '//trim(SubsID(s))//' '//trim(fld_unt(TDM_IND))&
                                &//trim(unt_sfx(itm))
            hdr(itm,AC_MON)=trim(tm_pfx(itm))//' '//trim(fld_nm(AC_MON))//' '//trim(SubsID(s))//' '//trim(fld_unt(AC_MON))
            hdr(itm,CP_MON)=trim(tm_pfx(itm))//' '//trim(fld_nm(CP_MON))//' '//trim(SubsID(s))//' '//trim(fld_unt(CP_MON))
            hdr(itm,WD_MON)=trim(tm_pfx(itm))//' '//trim(fld_nm(WD_MON))//' '//trim(SubsID(s))//' '//trim(fld_unt(WD_MON))&
                                &//trim(unt_sfx(itm))
        end do ! Time level
    end do ! Substances

! Create files for output of monitoring
    if(InitCond == 'zero' .or. InitCond == 'cond') then
        if (AtmOutMonitHourly)  call Atm_WriteMonitorHdr('Hourly')
        if (AtmOutMonit6hourly)  call Atm_WriteMonitorHdr('6hourly')
        if (AtmOutMonitDaily)   call Atm_WriteMonitorHdr('Daily')
        if (AtmOutMonitMonthly) call Atm_WriteMonitorHdr('Monthly')
        if (AtmOutMonitYearly)  call Atm_WriteMonitorHdr('Yearly')
    end if

end subroutine Atm_Output_Initial

subroutine Atm_Output_Final
    integer     dErr

    deallocate(fld_2d, stat=dErr)
    if(dErr/=0) stop 'STOP: Memory deallocation error in Atm_Output_Final'

end subroutine Atm_Output_Final


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine writing 2d fields to text files
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_WriteFields(prnType)
    character(*), intent(in)    ::  prnType
    real                        ::  tp
    integer                     ::  i,j,n,m,s,f,ind,fld_no, time_lev
    character*2                 ::  HourNum

    write(YearNum,'(i4)') Year

    do s = 1, NumSubsMedia(Atm)

        select case(prnType)
            case('Hourly')
                write(HourNum,'(i2.2)') numHour
                write(DayNum,'(i2.2)') Day
                write(MonthNum,'(i2.2)') Month
                fp_suffix = '_'//YearNum//MonthNum//DayNum//'_'//HourNum//'.dat'
                time_lev = TM_HOURLY
                fld_no = 1

            ! Add suffix to file name valid for current time
                do m = 1, fld_no
                    if(fld_out(time_lev,m)) filePath(m) = trim(fName(time_lev,m))//trim(fp_suffix)  
                end do

            ! copy surface layers for all forms of current substance
                do n = 1, sFormNum(Atm, s)
                    ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    if(fld_out(time_lev,AC_IND)) then        ! Air concentrations
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
#if G_HG
                                fld_2d(i,j,n,AC_IND) = Atm_MixRatio(i,j,1,ind) * stp * 1.e12     ! ng / m3
#else
                                fld_2d(i,j,n,AC_IND) = Atm_Conc(i,j,1,ind) * 1.e12     ! ng / m3
#endif
                            end do
                        end do
                    end if
                end do

            case('6hourly')
                write(PerNum,'(i2.2)') Period
                write(DayNum,'(i2.2)') Day
                write(MonthNum,'(i2.2)') Month
                fp_suffix = '_'//YearNum//MonthNum//DayNum//'_'//PerNum//'.dat'
                time_lev = TM_6HOURLY
                fld_no = 1

            ! Add suffix to file name valid for current time
                do m = 1, fld_no
                    if(fld_out(time_lev,m)) filePath(m) = trim(fName(time_lev,m))//trim(fp_suffix)  
                end do

            ! copy surface layers for all forms of current substance
                do n = 1, sFormNum(Atm, s)
                    ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    if(fld_out(time_lev,AC_IND)) then        ! Air concentrations
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
#if G_HG
                                fld_2d(i,j,n,AC_IND) = Atm_MixRatio(i,j,1,ind) * stp * 1.e12     ! ng / m3
#else
                                fld_2d(i,j,n,AC_IND) = Atm_Conc(i,j,1,ind) * 1.e12     ! ng / m3
#endif
                            end do
                        end do
                    end if
                end do

            case('Daily')
                write(DayNum,'(i2.2)') Day
                write(MonthNum,'(i2.2)') Month
                fp_suffix = '_'//YearNum//MonthNum//DayNum//'.dat'
                time_lev = TM_DAILY

            ! Add suffix to file name valid for current time
                do m = 1, FIELDS_NUM
                    if(fld_out(time_lev,m)) filePath(m) = trim(fName(time_lev,m))//trim(fp_suffix)  
                end do

                fld_no = FIELDS_BY_FORM                    ! Number of fields with forms division

            ! copy fields
                do n = 1, sFormNum(Atm, s)
                    ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms

                    if(fld_out(time_lev,AC_IND)) then        ! Air concentrations
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
#if G_HG
                                fld_2d(i,j,n,AC_IND) = sum(Atm_MixDay(i,j,1,ind,1:NumSrc))*stp*1.e12/SecInDay     ! ng / m3
#else
                                fld_2d(i,j,n,AC_IND) = sum(Atm_ConcDay(i,j,1,ind,1:NumSrc)) * 1.e12/SecInDay     ! ng / m3
#endif
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,CP_IND)) then        ! Concentrations in precipitation
                        fld_2d(:,:,n,CP_IND) = 0.
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                if(PrecipDay(i,j) /= 0.) then
                                    fld_2d(i,j,n,CP_IND) = sum(WetDepDay(i,j,ind,1:NumSurf,1:NumSrc)) / &
                                                    & PrecipDay(i,j) / MeshArea(i,j) * 1.e9             ! ng / L
                                end if
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,DD_IND)) then        ! Dry deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,DD_IND)= sum(DryDepDay(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/day
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,WD_IND)) then        ! Wet deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,WD_IND)= sum(WetDepDay(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/day
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,TD_IND)) then        ! Total deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,TD_IND)= (sum(WetDepDay(i,j,ind,1:NumSurf,1:NumSrc))+&
                                                 &sum(DryDepDay(i,j,ind,1:NumSurf,1:NumSrc)))/MeshArea(i,j)*1.e9 ! g/km2/day
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,ND_IND)) then        ! Net deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,ND_IND)= (sum(WetDepDay(i,j,ind,1:NumSurf,1:NumSrc))+&
                                                 &sum(DryDepDay(i,j,ind,1:NumSurf,1:NumSrc))-&
                                                 &sum(DryRemDay(i,j,ind,1:NumSurf,1:NumSrc)))/MeshArea(i,j)*1.e9 ! g/km2/yr
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,RE_IND)) then        ! Re-emission
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,RE_IND)= sum(DryRemDay(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/yr
                            end do
                        end do
                    end if
                end do ! n

                if(fld_out(time_lev,PA_IND)) then        ! Precipitation amount
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,1,PA_IND)= PrecipDay(i,j) *1.e3       ! mm/day
                        end do
                    end do
                end if

            case('Monthly')
                write(MonthNum,'(i2.2)') Month
                fp_suffix = '_'//YearNum//MonthNum//'.dat'
                time_lev = TM_MONTHLY

            ! Add suffix to file name valid for current time
                do m = 1, FIELDS_NUM
                    if(fld_out(time_lev,m)) filePath(m) = trim(fName(time_lev,m))//trim(fp_suffix)  
                end do
                fld_no = FIELDS_BY_FORM                    ! Number of fields with forms division

            ! copy fields
                do n = 1, sFormNum(Atm, s)
                    ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    if(fld_out(time_lev,AC_IND)) then        ! Air concentrations
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
#if G_HG
                                fld_2d(i,j,n,AC_IND) = sum(Atm_MixMonth(i,j,1,ind,1:NumSrc)) &
                                & *stp*1.e12/(min(SecInDay*real(MonthDays(Month)),timeCalc))     ! ng / m3
#else
                                fld_2d(i,j,n,AC_IND) = sum(Atm_ConcMonth(i,j,1,ind,1:NumSrc)) &
                                & * 1.e12/(min(SecInDay*real(MonthDays(Month)),timeCalc))     ! ng / m3
#endif
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,CP_IND)) then        ! Concentrations in precipitation
                        fld_2d(:,:,n,2) = 0.
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                if(PrecipMonth(i,j) /= 0.) then
                                    fld_2d(i,j,n,CP_IND) = sum(WetDepMonth(i,j,ind,1:NumSurf,1:NumSrc)) / &
                                                    & PrecipMonth(i,j) / MeshArea(i,j) * 1.e9             ! ng / L
                                end if
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,DD_IND)) then        ! Dry deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,DD_IND)= sum(DryDepMonth(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/day
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,WD_IND)) then        ! Wet deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,WD_IND)= sum(WetDepMonth(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/day
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,TD_IND)) then        ! Total deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,TD_IND)= (sum(WetDepMonth(i,j,ind,1:NumSurf,1:NumSrc))+&
                                                 &sum(DryDepMonth(i,j,ind,1:NumSurf,1:NumSrc)))/MeshArea(i,j)*1.e9 ! g/km2/day
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,ND_IND)) then        ! Net deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,ND_IND)= (sum(WetDepMonth(i,j,ind,1:NumSurf,1:NumSrc))+&
                                                 &sum(DryDepMonth(i,j,ind,1:NumSurf,1:NumSrc))-&
                                                 &sum(DryRemMonth(i,j,ind,1:NumSurf,1:NumSrc)))/MeshArea(i,j)*1.e9 ! g/km2/yr
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,RE_IND)) then        ! Re-emission
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,RE_IND)= sum(DryRemMonth(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/yr
                            end do
                        end do
                    end if
                end do ! n

                if(fld_out(time_lev,PA_IND)) then        ! Precipitation amount
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,1,PA_IND)= PrecipMonth(i,j) *1.e3       ! mm/day
                        end do
                    end do
                end if
                if(fld_out(time_lev,AE_IND)) then        ! Anthropogenic emission
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,1,AE_IND)= AntEmisMonth(i,j)/MeshArea(i,j)*1.e9       ! g/km2/day
                        end do
                    end do
                end if
                if(fld_out(time_lev,NE_IND)) then        ! Natural emission
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,1,NE_IND)= NatEmisMonth(i,j)/MeshArea(i,j)*1.e9       ! g/km2/day
                        end do
                    end do
                end if


            case('Yearly')
                fp_suffix = '_'//YearNum//'.dat'
                time_lev = TM_YEARLY

            ! Add suffix to file name valid for current time
                do m = 1, FIELDS_NUM
                    if(fld_out(time_lev,m)) filePath(m) = trim(fName(time_lev,m))//trim(fp_suffix)  
                end do
                fld_no = FIELDS_BY_FORM                    ! Number of fields with forms division

            ! copy fields
                do n = 1, sFormNum(Atm, s)
                    ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    if(fld_out(time_lev,AC_IND)) then        ! Air concentrations
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
#if G_HG
                                fld_2d(i,j,n,AC_IND) = sum(Atm_MixYear(i,j,1,ind,1:NumSrc)) &
                                & *stp*1.e12/timeCalc     ! ng / m3
#else
                                fld_2d(i,j,n,AC_IND) = sum(Atm_ConcYear(i,j,1,ind,1:NumSrc)) &
                                & * 1.e12/timeCalc     ! ng / m3
#endif
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,CP_IND)) then        ! Concentrations in precipitation
                        fld_2d(:,:,n,CP_IND) = 0.
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                if(PrecipYear(i,j) /= 0.) then
                                    fld_2d(i,j,n,CP_IND) = sum(WetDepYear(i,j,ind,1:NumSurf,1:NumSrc)) / &
                                                    & PrecipYear(i,j) / MeshArea(i,j) * 1.e9             ! ng / L
                                end if
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,DD_IND)) then        ! Dry deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,DD_IND)= sum(DryDepYear(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/day
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,WD_IND)) then        ! Wet deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,WD_IND)= sum(WetDepYear(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/day
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,TD_IND)) then        ! Total deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,TD_IND)= (sum(WetDepYear(i,j,ind,1:NumSurf,1:NumSrc))+&
                                                 &sum(DryDepYear(i,j,ind,1:NumSurf,1:NumSrc)))/MeshArea(i,j)*1.e9 ! g/km2/day
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,ND_IND)) then        ! Net deposition
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,ND_IND)= (sum(WetDepYear(i,j,ind,1:NumSurf,1:NumSrc))+&
                                                 &sum(DryDepYear(i,j,ind,1:NumSurf,1:NumSrc))-&
                                                 &sum(DryRemYear(i,j,ind,1:NumSurf,1:NumSrc)))/MeshArea(i,j)*1.e9 ! g/km2/yr
                            end do
                        end do
                    end if
                    if(fld_out(time_lev,RE_IND)) then        ! Re-emission
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                fld_2d(i,j,n,RE_IND)= sum(DryRemYear(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/yr
                            end do
                        end do
                    end if
                end do ! n

                if(fld_out(time_lev,PA_IND)) then        ! Precipitation amount
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,1,PA_IND)= PrecipYear(i,j) *1.e3       ! mm/day
                        end do
                    end do
                end if
                if(fld_out(time_lev,AE_IND)) then        ! Anthropogenic emission
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,1,AE_IND)= AntEmisYear(i,j)/MeshArea(i,j)*1.e9       ! g/km2/day
                        end do
                    end do
                end if
                if(fld_out(time_lev,NE_IND)) then        ! Natural emission
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,1,NE_IND)= NatEmisYear(i,j)/MeshArea(i,j)*1.e9       ! g/km2/day
                        end do
                    end do
                end if

        end select
    ! write to text file
        do n = 1, fld_no
            if(fld_out(time_lev,n)) call WriteFieldTxt(trim(filePath(n)),hdr(time_lev,n),fld_2d(:,:,:,n),Atm,s)
        end do
        if(fld_out(time_lev,AE_IND)) call WriteFieldTxt(trim(filePath(AE_IND)),hdr(time_lev,AE_IND),fld_2d(:,:,:,AE_IND),Atm,s,1)
        if(fld_out(time_lev,NE_IND)) call WriteFieldTxt(trim(filePath(NE_IND)),hdr(time_lev,NE_IND),fld_2d(:,:,:,NE_IND),Atm,s,1)
        if(fld_out(time_lev,PA_IND)) call WriteFieldTxt(trim(filePath(PA_IND)),hdr(time_lev,PA_IND),fld_2d(:,:,:,PA_IND),Atm,s,1,1)
    end do ! do s - loop for substances

end subroutine Atm_WriteFields




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating continent-to-continent matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if RTYPE==2

subroutine Atm_WriteMatrix(Prd)
! parameters
    character(*), intent(in)    :: Prd
! local vars
    integer         i, j, f, sf, src, rcp, s, time_lev, ind


    write(YearNum,'(i4)') Year

    do s = 1, NumSubsMedia(Atm)
        selectcase(Prd)
        case('Daily')
            write(DayNum,'(i2.2)') Day
            write(MonthNum,'(i2.2)') Month
            fp_suffix = '('//YearNum//MonthNum//DayNum//').dat'
            time_lev = TM_DAILY

        ! Add suffix to file name valid for current time
            filePath(ACM_IND) = trim(fName(time_lev,ACM_IND))//trim(fp_suffix)  
            filePath(WDM_IND) = trim(fName(time_lev,WDM_IND))//trim(fp_suffix)
            filePath(DDM_IND) = trim(fName(time_lev,DDM_IND))//trim(fp_suffix)
            filePath(TDM_IND) = trim(fName(time_lev,TDM_IND))//trim(fp_suffix)  

            TotDepMatrix=0.
            do src=1, NumSrc
                do j=Jmin, Jmax
                    do i=minI(j), maxI(j)
                        do sf=1, NumSurf
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
                                do rcp=1, NumRcp
                                    TotDepMatrix(src,rcp)=TotDepMatrix(src,rcp)+&
                                        &(DryDepDay(i,j,ind,sf,src)+WetDepDay(i,j,ind,sf,src))*RecepPart(i,j,rcp)
                                enddo ! rcp
                            enddo ! f
                        enddo ! sf
                    enddo ! i
                enddo ! j
            enddo ! src
            WetDepMatrix=0.
            do src=1, NumSrc
                do j=Jmin, Jmax
                    do i=minI(j), maxI(j)
                        do sf=1, NumSurf
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
                                do rcp=1, NumRcp
                                    WetDepMatrix(src,rcp)=WetDepMatrix(src,rcp)+&
                                        &(WetDepDay(i,j,ind,sf,src))*RecepPart(i,j,rcp)
                                enddo ! rcp
                            enddo ! f
                        enddo ! sf
                    enddo ! i
                enddo ! j
            enddo ! src
            DryDepMatrix=0.
            do src=1, NumSrc
                do j=Jmin, Jmax
                    do i=minI(j), maxI(j)
                        do sf=1, NumSurf
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
                                do rcp=1, NumRcp
                                    DryDepMatrix(src,rcp)=DryDepMatrix(src,rcp)+&
                                        &(DryDepDay(i,j,ind,sf,src))*RecepPart(i,j,rcp)
                                enddo ! rcp
                            enddo ! f
                        enddo ! sf
                    enddo ! i
                enddo ! j
            enddo ! src

            ConcMatrix=0.
            do src=1, NumSrc
                do rcp=1, NumRcp
                    do j=Jmin, Jmax
                        do i=minI(j), maxI(j)
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
#if G_HG
                                ConcMatrix(src,rcp)=ConcMatrix(src,rcp)+&
                                    &Atm_MixDay(i,j,1,ind,src)*stp*RecepPart(i,j,rcp)*MeshArea(i,j)
#else
                                ConcMatrix(src,rcp)=ConcMatrix(src,rcp)+&
                                    &Atm_ConcDay(i,j,1,ind,src)*RecepPart(i,j,rcp)*MeshArea(i,j)
#endif
                            enddo ! f
                        enddo ! i
                    enddo ! j
                    if(RecepArea(rcp) /= 0.) then
                        ConcMatrix(src,rcp)=ConcMatrix(src,rcp)/RecepArea(rcp)
                    else
                        print *, 'Error: zero area for the receptor ',rcp
                    end if
                enddo ! rcp
            enddo ! src
            ConcMatrix = ConcMatrix / SecInDay*1.e12

        case('Monthly')
            write(MonthNum,'(i2.2)') Month
            fp_suffix = '('//YearNum//MonthNum//').dat'
            time_lev = TM_MONTHLY

        ! Add suffix to file name valid for current time
            filePath(ACM_IND) = trim(fName(time_lev,ACM_IND))//trim(fp_suffix)  
            filePath(WDM_IND) = trim(fName(time_lev,WDM_IND))//trim(fp_suffix)
            filePath(DDM_IND) = trim(fName(time_lev,DDM_IND))//trim(fp_suffix)
            filePath(TDM_IND) = trim(fName(time_lev,TDM_IND))//trim(fp_suffix)  

            TotDepMatrix=0.
            do src=1, NumSrc
                do j=Jmin, Jmax
                    do i=minI(j), maxI(j)
                        do sf=1, NumSurf
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
                                do rcp=1, NumRcp
                                    TotDepMatrix(src,rcp)=TotDepMatrix(src,rcp)+&
                                        &(DryDepMonth(i,j,ind,sf,src)+WetDepMonth(i,j,ind,sf,src))*RecepPart(i,j,rcp)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            WetDepMatrix=0.
            do src=1, NumSrc
                do j=Jmin, Jmax
                    do i=minI(j), maxI(j)
                        do sf=1, NumSurf
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
                                do rcp=1, NumRcp
                                    WetDepMatrix(src,rcp)=WetDepMatrix(src,rcp)+&
                                        &(WetDepMonth(i,j,ind,sf,src))*RecepPart(i,j,rcp)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            DryDepMatrix=0.
            do src=1, NumSrc
                do j=Jmin, Jmax
                    do i=minI(j), maxI(j)
                        do sf=1, NumSurf
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
                                do rcp=1, NumRcp
                                    DryDepMatrix(src,rcp)=DryDepMatrix(src,rcp)+&
                                        &(DryDepMonth(i,j,ind,sf,src))*RecepPart(i,j,rcp)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            ConcMatrix=0.
            do src=1, NumSrc
                do rcp=1, NumRcp
                    do j=Jmin, Jmax
                        do i=minI(j), maxI(j)
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
#if G_HG
                                ConcMatrix(src,rcp)=ConcMatrix(src,rcp)+&
                                    &Atm_MixMonth(i,j,1,ind,src)*stp*RecepPart(i,j,rcp)*MeshArea(i,j)
#else
                                ConcMatrix(src,rcp)=ConcMatrix(src,rcp)+&
                                    &Atm_ConcMonth(i,j,1,ind,src)*RecepPart(i,j,rcp)*MeshArea(i,j)
#endif
                            enddo ! f
                        enddo ! i
                    enddo ! j
                    if(RecepArea(rcp) /= 0.) then
                        ConcMatrix(src,rcp)=ConcMatrix(src,rcp)/RecepArea(rcp)
                    else
                        print *, 'Error: zero area for the receptor ',rcp
                    end if
                enddo ! rcp
            enddo ! src
            ConcMatrix = ConcMatrix / min(SecInDay*real(MonthDays(Month)),timeCalc) * 1.e12

        case('Yearly')
            fp_suffix = '('//YearNum//').dat'
            time_lev = TM_YEARLY

        ! Add suffix to file name valid for current time
            filePath(ACM_IND) = trim(fName(time_lev,ACM_IND))//trim(fp_suffix)  
            filePath(WDM_IND) = trim(fName(time_lev,WDM_IND))//trim(fp_suffix)
            filePath(DDM_IND) = trim(fName(time_lev,DDM_IND))//trim(fp_suffix)
            filePath(TDM_IND) = trim(fName(time_lev,TDM_IND))//trim(fp_suffix)  

            TotDepMatrix=0.
            do src=1, NumSrc
                do j=Jmin, Jmax
                    do i=minI(j), maxI(j)
                        do sf=1, NumSurf
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
                                do rcp=1, NumRcp
                                    TotDepMatrix(src,rcp)=TotDepMatrix(src,rcp)+&
                                        &(DryDepYear(i,j,ind,sf,src)+WetDepYear(i,j,ind,sf,src))*RecepPart(i,j,rcp)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            WetDepMatrix=0.
            do src=1, NumSrc
                do j=Jmin, Jmax
                    do i=minI(j), maxI(j)
                        do sf=1, NumSurf
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
                                do rcp=1, NumRcp
                                    WetDepMatrix(src,rcp)=WetDepMatrix(src,rcp)+&
                                        &(WetDepYear(i,j,ind,sf,src))*RecepPart(i,j,rcp)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            DryDepMatrix=0.
            do src=1, NumSrc
                do j=Jmin, Jmax
                    do i=minI(j), maxI(j)
                        do sf=1, NumSurf
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
                                do rcp=1, NumRcp
                                    DryDepMatrix(src,rcp)=DryDepMatrix(src,rcp)+&
                                        &(DryDepYear(i,j,ind,sf,src))*RecepPart(i,j,rcp)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            ConcMatrix=0.
            do src=1, NumSrc
                do rcp=1, NumRcp
                    do j=Jmin, Jmax
                        do i=minI(j), maxI(j)
                            do f= 1, sFormNum(Atm, s)
                                ind = sFormInd(Atm, s, f)
#if G_HG
                                ConcMatrix(src,rcp)=ConcMatrix(src,rcp)+&
                                    &Atm_MixYear(i,j,1,ind,src)*stp*RecepPart(i,j,rcp)*MeshArea(i,j)
#else
                                ConcMatrix(src,rcp)=ConcMatrix(src,rcp)+&
                                    &Atm_ConcYear(i,j,1,ind,src)*RecepPart(i,j,rcp)*MeshArea(i,j)
#endif
                            enddo ! f
                        enddo ! i
                    enddo ! j
                    if(RecepArea(rcp) /= 0.) then
                        ConcMatrix(src,rcp)=ConcMatrix(src,rcp)/RecepArea(rcp)
                    else
                        print *, 'Error: zero area for the receptor ',rcp
                    end if
                enddo ! rcp
            enddo ! src
            ConcMatrix = ConcMatrix / timeCalc * 1.e12

        end select
    end do

! Write matrices
    call WriteMatrix(filePath(ACM_IND), hdr(time_lev,ACM_IND), ConcMatrix)
    call WriteMatrix(filePath(WDM_IND), hdr(time_lev,WDM_IND), WetDepMatrix)
    call WriteMatrix(filePath(DDM_IND), hdr(time_lev,DDM_IND), DryDepMatrix)
    call WriteMatrix(filePath(TDM_IND), hdr(time_lev,TDM_IND), TotDepMatrix)

end subroutine Atm_WriteMatrix


#endif



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutines to write concentrations, deposition, and precipitation amount at staions locations
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_WriteMonitor(Prd)
! parameters
    character(*), intent(in)    :: Prd

    integer         i, j, n, s, ind
    character(17)   tm                  ! time string
    character(120)  fp, fn              
    real            cf

    do s = 1, NumSubsMedia(Atm)

        select case(Prd)


        case('Hourly')
                write(tm,'(1x,i4,a1,i2.2,a1,i2.2,a1,i2.2,a3)') Year,'-',Month,'-',Day,' ',numHour,':00'
                fp = trim(OutPath)//trim(AtmDir)//trim(MonitorDir)//trim(HourlyDir)      
                do n = 1, sFormNum(Atm, s)
                    ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
#if G_HG
                            fld_2d(i,j,n,AC_IND) = Atm_MixRatio(i,j,1,ind) * stp * 1.e12     ! ng / m3
#else
                            fld_2d(i,j,n,AC_IND) = Atm_Conc(i,j,1,ind) * 1.e12     ! ng / m3
#endif
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,AC_IND), 1)
                    do i = 1, NAir
                        AirStat(i)%stMod(1,1,ind)= fld_2d(AirStat(i)%iSt(1,1),AirStat(i)%jSt(1,1),n,AC_IND)
                        AirStat(i)%stMod(1,2,ind)= fld_2d(AirStat(i)%iSt(1,2),AirStat(i)%jSt(1,2),n,AC_IND)
                        AirStat(i)%stMod(2,1,ind)= fld_2d(AirStat(i)%iSt(2,1),AirStat(i)%jSt(2,1),n,AC_IND)
                        AirStat(i)%stMod(2,2,ind)= fld_2d(AirStat(i)%iSt(2,2),AirStat(i)%jSt(2,2),n,AC_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_HOURLY,AC_MON))//'hourly.dat'
                    call WriteMonitor(trim(fn),AirStat, NAir, tm, ind)                                    
                end do

        case('6hourly')
                write(tm,'(1x,i4,a1,i2.2,a1,i2.2,a1,i2.2,a3)') Year,'-',Month,'-',Day,' ',numHour,':00'
                fp = trim(OutPath)//trim(AtmDir)//trim(MonitorDir)//trim(SixHourlyDir)      
                do n = 1, sFormNum(Atm, s)
                    ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
#if G_HG
                            fld_2d(i,j,n,AC_IND) = Atm_MixRatio(i,j,1,ind) * stp * 1.e12     ! ng / m3
#else
                            fld_2d(i,j,n,AC_IND) = Atm_Conc(i,j,1,ind) * 1.e12     ! ng / m3
#endif
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,AC_IND), 1)
                    do i = 1, NAir
                        AirStat(i)%stMod(1,1,ind)= fld_2d(AirStat(i)%iSt(1,1),AirStat(i)%jSt(1,1),n,AC_IND)
                        AirStat(i)%stMod(1,2,ind)= fld_2d(AirStat(i)%iSt(1,2),AirStat(i)%jSt(1,2),n,AC_IND)
                        AirStat(i)%stMod(2,1,ind)= fld_2d(AirStat(i)%iSt(2,1),AirStat(i)%jSt(2,1),n,AC_IND)
                        AirStat(i)%stMod(2,2,ind)= fld_2d(AirStat(i)%iSt(2,2),AirStat(i)%jSt(2,2),n,AC_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_6HOURLY,AC_MON))//'6hourly.dat'
                    call WriteMonitor(trim(fn),AirStat, NAir, tm, ind)                                    
                end do

        case('Daily')
                write(tm,'(1x,i4,a1,i2.2,a1,i2.2,"      ")') Year,'-',Month,'-',Day
                fp = trim(OutPath)//trim(AtmDir)//trim(MonitorDir)//trim(DailyDir)      
                cf = 1.e12/SecInDay
                do n = 1, sFormNum(Atm, s)
           ! Air conc monitoring
                    ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
#if G_HG
                            fld_2d(i,j,n,AC_IND) = sum(Atm_MixDay(i,j,1,ind,1:NumSrc))*stp*cf     ! ng / m3
#else
                            fld_2d(i,j,n,AC_IND) = sum(Atm_ConcDay(i,j,1,ind,1:NumSrc))*cf     ! ng / m3
#endif
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,AC_IND), 1)
                    do i = 1, NAir
                        AirStat(i)%stMod(1,1,ind)= fld_2d(AirStat(i)%iSt(1,1),AirStat(i)%jSt(1,1),n,AC_IND)
                        AirStat(i)%stMod(1,2,ind)= fld_2d(AirStat(i)%iSt(1,2),AirStat(i)%jSt(1,2),n,AC_IND)
                        AirStat(i)%stMod(2,1,ind)= fld_2d(AirStat(i)%iSt(2,1),AirStat(i)%jSt(2,1),n,AC_IND)
                        AirStat(i)%stMod(2,2,ind)= fld_2d(AirStat(i)%iSt(2,2),AirStat(i)%jSt(2,2),n,AC_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_DAILY,AC_MON))//'daily.dat'
                    call WriteMonitor(trim(fn),AirStat, NAir, tm, ind)                                    
           ! Prec conc monitoring
                    fld_2d(:,:,n,CP_IND) = 0.
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            if(PrecipDay(i,j) /= 0.) then
                                fld_2d(i,j,n,CP_IND) = sum(WetDepDay(i,j,ind,1:NumSurf,1:NumSrc)) / &
                                                & PrecipDay(i,j) / MeshArea(i,j) * 1.e9             ! ng / L
                            end if
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,CP_IND), 1)
                    do i = 1, NPrec
                        PrecStat(i)%stMod(1,1,ind)= fld_2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),n,CP_IND)
                        PrecStat(i)%stMod(1,2,ind)= fld_2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),n,CP_IND)
                        PrecStat(i)%stMod(2,1,ind)= fld_2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),n,CP_IND)
                        PrecStat(i)%stMod(2,2,ind)= fld_2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),n,CP_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_DAILY,CP_MON))//'daily.dat'
                    call WriteMonitor(trim(fn),PrecStat, NPrec, tm, ind)
           ! Wet dep monitoring
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,n,WD_IND)= sum(WetDepDay(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/day
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,WD_IND), 1)
                    do i = 1, NPrec
                        PrecStat(i)%stMod(1,1,ind)= fld_2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),n,WD_IND)
                        PrecStat(i)%stMod(1,2,ind)= fld_2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),n,WD_IND)
                        PrecStat(i)%stMod(2,1,ind)= fld_2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),n,WD_IND)
                        PrecStat(i)%stMod(2,2,ind)= fld_2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),n,WD_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_DAILY,WD_MON))//'daily.dat'
                    call WriteMonitor(trim(fn),PrecStat, NPrec, tm, ind)
                end do

        case('Monthly')
                write(tm,'(1x,i4,a1,i2.2,"         ")') Year,'-',Month
                fp = trim(OutPath)//trim(AtmDir)//trim(MonitorDir)//trim(MonthlyDir)      
                cf = 1.e12/min(SecInDay*real(MonthDays(Month)),timeCalc)
                do n = 1, sFormNum(Atm, s)
                    ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
           ! Air conc monitoring
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
#if G_HG
                            fld_2d(i,j,n,AC_IND) = sum(Atm_MixMonth(i,j,1,ind,1:NumSrc))*stp*cf     ! ng / m3
#else
                            fld_2d(i,j,n,AC_IND) = sum(Atm_ConcMonth(i,j,1,ind,1:NumSrc))*cf     ! ng / m3
#endif
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,AC_IND), 1)
                    do i = 1, NAir
                        AirStat(i)%stMod(1,1,ind)= fld_2d(AirStat(i)%iSt(1,1),AirStat(i)%jSt(1,1),n,AC_IND)
                        AirStat(i)%stMod(1,2,ind)= fld_2d(AirStat(i)%iSt(1,2),AirStat(i)%jSt(1,2),n,AC_IND)
                        AirStat(i)%stMod(2,1,ind)= fld_2d(AirStat(i)%iSt(2,1),AirStat(i)%jSt(2,1),n,AC_IND)
                        AirStat(i)%stMod(2,2,ind)= fld_2d(AirStat(i)%iSt(2,2),AirStat(i)%jSt(2,2),n,AC_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_MONTHLY,AC_MON))//'monthly.dat'
                    call WriteMonitor(trim(fn),AirStat, NAir, tm, ind)                                    
           ! Prec conc monitoring
                    fld_2d(:,:,n,CP_IND) = 0.
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            if(PrecipMonth(i,j) /= 0.) then
                                fld_2d(i,j,n,CP_IND) = sum(WetDepMonth(i,j,ind,1:NumSurf,1:NumSrc)) / &
                                                & PrecipMonth(i,j) / MeshArea(i,j) * 1.e9             ! ng / L
                            end if
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,CP_IND), 1)
                    do i = 1, NPrec
                        PrecStat(i)%stMod(1,1,ind)= fld_2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),n,CP_IND)
                        PrecStat(i)%stMod(1,2,ind)= fld_2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),n,CP_IND)
                        PrecStat(i)%stMod(2,1,ind)= fld_2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),n,CP_IND)
                        PrecStat(i)%stMod(2,2,ind)= fld_2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),n,CP_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_MONTHLY,CP_MON))//'monthly.dat'
                    call WriteMonitor(trim(fn),PrecStat, NPrec, tm, ind)
           ! Wet dep monitoring
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,n,WD_IND)= sum(WetDepMonth(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/day
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,WD_IND), 1)
                    do i = 1, NPrec
                        PrecStat(i)%stMod(1,1,ind)= fld_2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),n,WD_IND)
                        PrecStat(i)%stMod(1,2,ind)= fld_2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),n,WD_IND)
                        PrecStat(i)%stMod(2,1,ind)= fld_2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),n,WD_IND)
                        PrecStat(i)%stMod(2,2,ind)= fld_2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),n,WD_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_MONTHLY,WD_MON))//'monthly.dat'
                    call WriteMonitor(trim(fn),PrecStat, NPrec, tm, ind)
                end do

        case('Yearly')
                write(tm,'(1x,i4,"            ")') Year
                fp = trim(OutPath)//trim(AtmDir)//trim(MonitorDir)//trim(YearlyDir)      
                cf = 1.e12/timeCalc
                do n = 1, sFormNum(Atm, s)
                    ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
           ! Air conc monitoring
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
#if G_HG
                            fld_2d(i,j,n,AC_IND) = sum(Atm_MixYear(i,j,1,ind,1:NumSrc))*stp*cf     ! ng / m3
#else
                            fld_2d(i,j,n,AC_IND) = sum(Atm_ConcYear(i,j,1,ind,1:NumSrc))*cf     ! ng / m3
#endif
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,AC_IND), 1)
                    do i = 1, NAir
                        AirStat(i)%stMod(1,1,ind)= fld_2d(AirStat(i)%iSt(1,1),AirStat(i)%jSt(1,1),n,AC_IND)
                        AirStat(i)%stMod(1,2,ind)= fld_2d(AirStat(i)%iSt(1,2),AirStat(i)%jSt(1,2),n,AC_IND)
                        AirStat(i)%stMod(2,1,ind)= fld_2d(AirStat(i)%iSt(2,1),AirStat(i)%jSt(2,1),n,AC_IND)
                        AirStat(i)%stMod(2,2,ind)= fld_2d(AirStat(i)%iSt(2,2),AirStat(i)%jSt(2,2),n,AC_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_YEARLY,AC_MON))//'yearly.dat'
                    call WriteMonitor(trim(fn),AirStat, NAir, tm, ind)                                    
           ! Prec conc monitoring
                    fld_2d(:,:,n,CP_IND) = 0.
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            if(PrecipYear(i,j) /= 0.) then
                                fld_2d(i,j,n,CP_IND) = sum(WetDepYear(i,j,ind,1:NumSurf,1:NumSrc)) / &
                                                & PrecipYear(i,j) / MeshArea(i,j) * 1.e9             ! ng / L
                            end if
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,CP_IND), 1)
                    do i = 1, NPrec
                        PrecStat(i)%stMod(1,1,ind)= fld_2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),n,CP_IND)
                        PrecStat(i)%stMod(1,2,ind)= fld_2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),n,CP_IND)
                        PrecStat(i)%stMod(2,1,ind)= fld_2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),n,CP_IND)
                        PrecStat(i)%stMod(2,2,ind)= fld_2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),n,CP_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_YEARLY,CP_MON))//'yearly.dat'
                    call WriteMonitor(trim(fn),PrecStat, NPrec, tm, ind)
           ! Wet dep monitoring
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,n,WD_IND)= sum(WetDepYear(i,j,ind,1:NumSurf,1:NumSrc))/MeshArea(i,j)*1.e9 ! g/km2/day
                        end do
                    end do
                    call GridDisAggreg(fld_2d(:,:,n,WD_IND), 1)
                    do i = 1, NPrec
                        PrecStat(i)%stMod(1,1,ind)= fld_2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),n,WD_IND)
                        PrecStat(i)%stMod(1,2,ind)= fld_2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),n,WD_IND)
                        PrecStat(i)%stMod(2,1,ind)= fld_2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),n,WD_IND)
                        PrecStat(i)%stMod(2,2,ind)= fld_2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),n,WD_IND)
                    end do
                    fn=trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_YEARLY,WD_MON))//'yearly.dat'
                    call WriteMonitor(trim(fn),PrecStat, NPrec, tm, ind)
                end do

        end select
    end do

end subroutine Atm_WriteMonitor



subroutine Atm_WriteMonitorHdr(Prd)
! parameters
    character(*), intent(in)    :: Prd

    integer         s, n, ind
    character(120)  fp,fn               

    do s = 1, NumSubsMedia(Atm)

        select case(Prd)

        case('Hourly')
            fp = trim(OutPath)//trim(AtmDir)//trim(MonitorDir)//trim(HourlyDir)                    
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_HOURLY,AC_MON))//'hourly.dat'
                call WriteMonitorHeader(trim(fn),hdr(TM_HOURLY, AC_MON), AirStat, NAir)                 
            end do

        case('6hourly')
            fp = trim(OutPath)//trim(AtmDir)//trim(MonitorDir)//trim(SixHourlyDir)                    
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_6HOURLY,AC_MON))//'6hourly.dat'
                call WriteMonitorHeader(trim(fn),hdr(TM_6HOURLY, AC_MON), AirStat, NAir)                 
            end do

        case('Daily')
            fp = trim(OutPath)//trim(AtmDir)//trim(MonitorDir)//trim(DailyDir)                    
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_DAILY,AC_MON))//'daily.dat'
                call WriteMonitorHeader(trim(fn), hdr(TM_DAILY, AC_MON), AirStat, NAir)                   
            end do
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_DAILY,CP_MON))//'daily.dat'
                call WriteMonitorHeader(trim(fn), hdr(TM_DAILY, CP_MON), PrecStat, NPrec)               
            end do
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_DAILY,WD_MON))//'daily.dat'
                call WriteMonitorHeader(trim(fn), hdr(TM_DAILY, WD_MON), PrecStat, NPrec)               
            end do

        case('Monthly')
            fp = trim(OutPath)//trim(AtmDir)//trim(MonitorDir)//trim(MonthlyDir)                    
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_MONTHLY,AC_MON))//'monthly.dat'
                call WriteMonitorHeader(trim(fn), hdr(TM_MONTHLY, AC_MON), AirStat, NAir)               
            end do
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_MONTHLY,CP_MON))//'monthly.dat'
                call WriteMonitorHeader(trim(fn), hdr(TM_MONTHLY, CP_MON), PrecStat, NPrec)             
            end do
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_MONTHLY,WD_MON))//'monthly.dat'
                call WriteMonitorHeader(trim(fn), hdr(TM_MONTHLY, WD_MON), PrecStat, NPrec)               
            end do

        case('Yearly')
            fp = trim(OutPath)//trim(AtmDir)//trim(MonitorDir)//trim(YearlyDir)                    
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_YEARLY,AC_MON))//'yearly.dat'
                call WriteMonitorHeader(trim(fn),hdr(TM_YEARLY, AC_MON), AirStat, NAir)                 
            end do
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_YEARLY,CP_MON))//'yearly.dat'
                call WriteMonitorHeader(trim(fn),hdr(TM_YEARLY, CP_MON), PrecStat, NPrec)               
            end do
            do n = 1, sFormNum(Atm, s)
                ind = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                fn = trim(fp)//trim(SubsID(s))//trim(FormID(Atm,s,ind))//trim(fName(TM_YEARLY,WD_MON))//'yearly.dat'
                call WriteMonitorHeader(trim(fn), hdr(TM_YEARLY, WD_MON), PrecStat, NPrec)               
            end do

        end select
    end do

end subroutine Atm_WriteMonitorHdr


#if G_POP

subroutine DepSurf2LU       ! 17-10-2019

  real*8 Coeff(NumSurf), Ss(NumSrc), Norm
  integer Surf,L,i,j,f,Src

  DryDepDay = 0.
  WetDepDay = 0.
  DryRemDay = 0.
  do f=1, NumForm(Atm)
         do Src = 1, NumSrc
    do j = JMin, JMax
      do i = MinI(j), MaxI(j)
	    do Surf = 0, Soil_NumType
		  do L = 1, NumSurf
		    Coeff(L) = LC2SoilType(Surf,L) * LandCover(i,j,L)
		  end do
		  Norm = sum(Coeff)
		  if (Norm > 0.) then
		    Coeff = Coeff / Norm
! Dry deposition
		    if (Surf == 0) Ss(Src) = DryDepDayTmp(i,j,f,Soil_NumType+3*Veg_NumType+1,Src)
		    
		    if (Surf > 0 .and. Surf <= Veg_NumType) Ss(Src) = DryDepDayTmp(i,j,f,Surf,Src) + &
		        & DryDepDayTmp(i,j,f,Soil_NumType+Surf,Src) + &
		        & DryDepDayTmp(i,j,f,Soil_NumType+Veg_NumType+Surf,Src)+&
		        & DryDepDayTmp(i,j,f,Soil_NumType+2*Veg_NumType+Surf,Src)
		    
	            if (Surf > Veg_NumType) Ss(Src) = DryDepDayTmp(i,j,f,Surf,Src) 
	            
		    do L = 1, NumSurf
		      DryDepDay(i,j,f,L,Src) = DryDepDay(i,j,f,L,Src) + Ss(Src) * Coeff(L)
		    end do
! Wet deposition
		    if (Surf == 0) Ss(Src) = WetDepDayTmp(i,j,f,Soil_NumType+3*Veg_NumType+1,Src) 
		    
		    if (Surf > 0 .and. Surf <= Veg_NumType) Ss(Src) = WetDepDayTmp(i,j,f,Surf,Src)+&
		        & WetDepDayTmp(i,j,f,Soil_NumType+Surf,Src) + &
		        & WetDepDayTmp(i,j,f,Soil_NumType+Veg_NumType+Surf,Src)+&
		        & WetDepDayTmp(i,j,f,Soil_NumType+2*Veg_NumType+Surf,Src)
		          
	            if (Surf > Veg_NumType) Ss(Src) = WetDepDayTmp(i,j,f,Surf,Src) 
	            
		    do L = 1, NumSurf
		      WetDepDay(i,j,f,L,Src) = WetDepDay(i,j,f,L,Src) + Ss(Src) * Coeff(L)
		    end do
! Re-emission
		    if (Surf == 0) Ss(Src) = DryRemDayTmp(i,j,f,Soil_NumType+3*Veg_NumType+1,Src) 
		    
		    if (Surf > 0 .and. Surf <= Veg_NumType) Ss(Src) = DryRemDayTmp(i,j,f,Surf,Src)+&
		          & DryRemDayTmp(i,j,f,Soil_NumType+Surf,Src) + &
		          & DryRemDayTmp(i,j,f,Soil_NumType+Veg_NumType+Surf,Src)+&
		          & DryRemDayTmp(i,j,f,Soil_NumType+2*Veg_NumType+Surf,Src)
		          
	            if (Surf > Veg_NumType) Ss(Src) = DryRemDayTmp(i,j,f,Surf,Src) 
	            
		    do L = 1, NumSurf
		      DryRemDay(i,j,f,L,Src) = DryRemDay(i,j,f,L,Src) + Ss(Src) * Coeff(L)
		    end do
		  end if
	    end do
      end do ! i
    end do ! j
   end do ! Src
  end do ! f
  
end subroutine DepSurf2LU

#endif

end module Atm_Output

