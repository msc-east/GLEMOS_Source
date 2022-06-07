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
! Module of the ocean  information input
! Version:  1.0
! Modified: 03.12.2014
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifdef M_OCN
module Ocn_Input

    use GeneralOutput
    use Ocn_Params
    use GeneralParams
    use Balance


  implicit none

  character(160), private :: UFileName,VFileName,HFileName
  character(160), private :: UFile,VFile,HFile
  character(160), private :: DirName
  character(80), private :: TimeString0,TimeString1
  character(800), private :: fileName, fullName

 contains

!========================================================================
! Initialize oceanic arrays
!========================================================================
subroutine Ocn_Initial

  integer i, j, k, FileStat, Surf

#ifdef DEBUG_MODE
    print *, '+Entering Ocn_Initial... ', Period
#endif

! Calculating water fractions
        fileName=trim(LC2OcnName)
        fullName=trim(GeoPath)//trim(LC2OcnName)
        open(10, file=trim(fullName), status = 'old', iostat=FileStat, action='read')
        if(FileStat>0) then
            print '(/,"STOP: Cannot open file LC2Ocn ''",a,"''",/)', trim(fullName)
	    stop
        endif
        read(10,*) LC2Ocn(1:NumSurf)
        close(10)

        Ocn_Frac = 0.
        do j = JMin, Jmax
            do i = MinI(j), MaxI(j)
                do Surf = 1, NumSurf
                    Ocn_Frac(i,j) = Ocn_Frac(i,j) + LC2Ocn(Surf) * LandCover(i,j,Surf)
                end do
            end do
        end do

! Reading ocean config file
        filename='Ocn'//trim(MediaConfName)//'.dat'
        fullName=trim(ConfigPath)//trim(filename)
        open(10, file=trim(fullName), status = 'old', iostat=FileStat, action='read')
        if(FileStat>0) then
            print '(/,"STOP in Ocn_Config: Cannot open file ''",a,"''",/)', trim(fullName)
            stop
        endif

        read(10,*) dzwm0
        read(10,*) dhf
        read(10,*) vdc_const        ! cm/s2
        read(10,*) ah               ! cm/s2
        read(10,*) SedimCoeff
        read(10,*) RPart
! Read vertical grid structure
        read(10,*)
        do k=1,Ocn_Kmax
            read(10,*) Ocn_dz(k)
        end do
        close(10)

! Set upper water layer extent
        Water_upl_dz=Ocn_dz(1)/100.

! Calculating auxiliary parameters

!	  RSeaAdd = Water_upl_dz / (kdv + kdvUML)
!	  WSed = SedimCoeff * RPart * RPart
        kdvUML=100.
        WSed=9.0e-04                                        !2011: cm/s
	RSeaAdd = Water_upl_dz / (vdc_const/100. *kdvUML)   !2011

      call initialize_ocean_advdiff

#ifdef DEBUG_MODE
    print *, '+Exit Ocn_Initial ', Period
#endif
end subroutine Ocn_Initial

!========================================================================
! Reading ocean input data
!========================================================================
subroutine Ocn_InputData(Period)

  character(*) Period
  integer g

#ifdef DEBUG_MODE
    print *, '+Entering Ocn_InputData... ', Period
#endif

  select case (Period)

    case ('Initial')

! Initializating arrays
        Ocn_VExch = 0.
        Water_upl_concDay = 0.
        Water_upl_concMonth = 0.
        Ocn_concDay = 0.
        Ocn_concMonth = 0.

        if(InitCond == 'zero') then
            Ocn_conc = 0.
            Water_upl_conc = 0.
            Water_upl_concYear = 0.
            Ocn_concYear = 0.
            Ocn_Degr = 0.
            Ocn_Sedim = 0.
        elseif(InitCond == 'cond') then
            Water_upl_concYear = 0.
            Ocn_concYear = 0.
            Ocn_Degr = 0.
            Ocn_Sedim = 0.
        end if
        
	  call Ocn_Initial

	  do g = 1, NumGroups
		  select case (SubsGrpLst(g))
	        case (POP)                          
		      call Ocn_POPInput('Initial')
	        case default
		      continue
	      end select
	  end do

	case ('Steply')

      do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                        
		      call Ocn_POPInput('Steply')
	        case default
		      continue
	      end select
      end do

	case ('Daily')

	  do g = 1, NumGroups
 	      select case (SubsGrpLst(g))
	        case (POP)                          
		      call Ocn_POPInput('Daily')
	        case default
		      continue
	      end select
      end do

	case ('Monthly')

      do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                          
		      call Ocn_POPInput('Monthly')
	        case default
		      continue
	      end select
      end do

	case ('Yearly')

      do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                          
		      call Ocn_POPInput('Yearly')
	        case default
		      continue
	      end select
      end do

	case ('6hourly')

#if (OCEAN_RUN==1)
    call ReadOcnInputData
#else
    ! Swithced off for run without ocean Sep 2016 V. Sh
#endif

  end select

#ifdef DEBUG_MODE
    print *, '+Exit Ocn_InputData ', Period
#endif
end subroutine Ocn_InputData

!========================================================================
! Reading ocean input U- V- and H-fields
!========================================================================
  subroutine ReadOcnInputData

     integer i,j,k,ii,jj,kk

     character cchar

#ifdef DEBUG_MODE
    print *, '+Entering ReadOcnInputData...'
#endif

!    chYear='2009'
    iMonth=Month
    iDay=Day
    iHour=Period

!                  write(TimeString0,'(i4,i2.2,i2.2,a2)')Year,iMonth,iDay,chTime(iHour)
    write(TimeString0,'(a4,i2.2,i2.2,a2)')'2009',iMonth,iDay,chTime(iHour)

    if (iHour.ne.4) then
        iHour1=iHour+1
        iDay1=iDay
        iMonth1=iMonth
    else
        iHour1=1
        if (iDay.ne.MonthDays(iMonth)) then
            iDay1=iDay+1
            iMonth1=iMonth
        else
            iDay1=1
            if (iMonth.ne.12) then
                iMonth1=iMonth+1
            else
                iMonth1=1
            end if
        end if
    end if

!                 write(TimeString1,'(i4,i2.2,i2.2,a2)')Year,iMonth1,iDay1,chTime(iHour1)
    write(TimeString1,'(a4,i2.2,i2.2,a2)')'2009',iMonth1,iDay1,chTime(iHour1)

    UFileName='UVEL'//trim(TimeString0)//'.dat'
    VFileName='VVEL'//trim(TimeString0)//'.dat'
    HFileName='HEIG'//trim(TimeString0)//'.dat'

    UFile=trim(OceanPath)//trim(UFileName)
    VFile=trim(OceanPath)//trim(VFileName)
    HFile=trim(OceanPath)//trim(HFileName)

    open (unit=1233,file=UFile)
    open (unit=1244,file=VFile)
    open (unit=1255,file=HFile)

    do i=1,52
        read(1233,*)cchar
        read(1244,*)cchar
        read(1255,*)cchar
    end do

    read(1233,*)UVEL0
    read(1244,*)VVEL0
    read(1255,*)HEIG0

    close(1233)
    close(1244)
    close(1255)

    UFileName='UVEL'//trim(TimeString1)//'.dat'
    VFileName='VVEL'//trim(TimeString1)//'.dat'
    HFileName='HEIG'//trim(TimeString1)//'.dat'

    UFile=trim(OceanPath)//trim(UFileName)
    VFile=trim(OceanPath)//trim(VFileName)
    HFile=trim(OceanPath)//trim(HFileName)

    open (unit=1233,file=UFile)
    open (unit=1244,file=VFile)
    open (unit=1255,file=HFile)

    do i=1,52
        read(1233,*)cchar
        read(1244,*)cchar
        read(1255,*)cchar
    end do

    read(1233,*)UVEL1
    read(1244,*)VVEL1
    read(1255,*)HEIG1

    close(1233)
    close(1244)
    close(1255)

#ifdef DEBUG_MODE
    print *, '+Exit ReadOcnInputData'
#endif
end subroutine ReadOcnInputData

!========================================================================
! Read ocean config from file
!========================================================================
subroutine Ocn_ReadConfig

    integer Subs, Grp

! Reading physical and chemical properties
      do Subs=1, NumSubs
        Grp=SubsGroup(Subs)
        selectcase(Grp)
          case(HM)
 !           call Ocn_HM_Props(Subs)
          case(Hg)
 !           call Ocn_Hg_Props(Subs)
          case(POP)
            call Ocn_POP_Props(Subs)
          case(Tracer)
!            call Ocn_Tracer_Props(Subs)
        endselect
      enddo
      call LogFile('Reading physical and chemical properties (ocean)',0)

 end subroutine Ocn_ReadConfig

!========================================================================
! Reading ocean input data (POPs)
!========================================================================
subroutine Ocn_POPInput(Period)

  character(*) Period

  select case (Period)

    case('Initial')

        if(InitCond == 'cond') then
            OcnInit = Ocn_POPMass(1)          ! Initial mass in ocean in case of 'cond' parameter
        elseif(InitCond == 'zero') then
            OcnInit = 0.
        end if

      continue

    case ('Daily')

      continue

  end select

end subroutine Ocn_POPInput

!========================================================================
! Linear interpolation of ocean input fields
!========================================================================
subroutine InterpolateInputOceanData

     integer i,j,k

     PeriodTime = DayTime-dTInput*(Period-1)

      if (kocean.ne.0) then

     UVEL(:,:,:,oldtime)=UVEL(:,:,:,curtime)
     UVEL(:,:,:,curtime)=UVEL(:,:,:,newtime)
     VVEL(:,:,:,oldtime)=VVEL(:,:,:,curtime)
     VVEL(:,:,:,curtime)=VVEL(:,:,:,newtime)
     DH(:,:,oldtime)=DH(:,:,curtime)
     DH(:,:,curtime)=DH(:,:,newtime)

     UVEL(:,:,:,newtime)=(UVEL0(:,:,:)*(dTInput-PeriodTime)+UVEL1(:,:,:)*PeriodTime)/dTInput
     VVEL(:,:,:,newtime)=(VVEL0(:,:,:)*(dTInput-PeriodTime)+VVEL1(:,:,:)*PeriodTime)/dTInput
     DH(:,:,newtime)=(HEIG0(:,:)*(dTInput-PeriodTime)+HEIG1(:,:)*PeriodTime)/dTInput

     else

     UVEL(:,:,:,newtime)=(UVEL0(:,:,:)*(dTInput-PeriodTime)+UVEL1(:,:,:)*PeriodTime)/dTInput
     VVEL(:,:,:,newtime)=(VVEL0(:,:,:)*(dTInput-PeriodTime)+VVEL1(:,:,:)*PeriodTime)/dTInput
     DH(:,:,newtime)=(HEIG0(:,:)*(dTInput-PeriodTime)+HEIG1(:,:)*PeriodTime)/dTInput

     UVEL(:,:,:,curtime)=UVEL(:,:,:,newtime)
     UVEL(:,:,:,oldtime)=UVEL(:,:,:,curtime)
     VVEL(:,:,:,curtime)=VVEL(:,:,:,newtime)
     VVEL(:,:,:,oldtime)=VVEL(:,:,:,curtime)
     DH(:,:,curtime)=DH(:,:,newtime)
     DH(:,:,oldtime)=DH(:,:,curtime)

     end if

end subroutine InterpolateInputOceanData


  end module Ocn_Input
#endif
