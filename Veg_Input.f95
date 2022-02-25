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
! Module of the vegetation information input
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifdef M_VEG
module Veg_Input

    use GeneralOutput
    use Veg_Params
    use GeneralParams
    use Geometry
    use Veg_POP_General
    use Balance

    implicit none

    character(800), private :: fileName, fullName

    contains

!========================================================================
! Reading vegetation input data
!========================================================================
 subroutine Veg_InputData(Period)

  character(*) Period
  integer g

#ifdef DEBUG_MODE
    print *, '+Entering Veg_InputData... ', Period
#endif

  select case (Period)

    case ('Initial')

	  call Veg_Initial

	  do g = 1, NumGroups
		  select case (SubsGrpLst(g))
	        case (POP)                         
		      call Veg_POPInput('Initial')
	        case default
		      continue
	      end select
	  end do

	case ('Steply')

      do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                          
		      call Veg_POPInput('Steply')
	        case default
		      continue
	      end select
      end do

	case ('Daily')

!      call Veg_Daily

	  do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                          
		      call Veg_POPInput('Daily')
	        case default
		      continue
	      end select
      end do

	case ('Monthly')

      do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                          
		      call Veg_POPInput('Monthly')
	        case default
		      continue
	      end select
      end do

	case ('Yearly')

      do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                          
		      call Veg_POPInput('Yearly')
	        case default
		      continue
	      end select
      end do

  end select

#ifdef DEBUG_MODE
    print *, '+Exit Veg_InputData ', Period
#endif
end subroutine Veg_InputData

!========================================================================
! Read vegetation configuration from file
!========================================================================
subroutine Veg_ReadConfig

  integer FileStat, ioRes, NP
  integer lenType, k
  character(80) strRead, strPar
  character(20) strType
  character(160) fullName
  integer Subs,Grp

  filename='Veg'//trim(MediaConfName)//'.dat'
  fullName=trim(ConfigPath)//trim(filename)
  open(2, file = trim(fullName), status='old', iostat=FileStat, action='read')
	if(FileStat>0) then
            print '(/,"STOP in Veg_Config: Cannot open file ''",a,"''",/)', trim(fullName)
	  stop
	endif

    do while(.true.)
	  read(2,'(a)', iostat = ioRes) strRead
	  if (ioRes /= 0) exit

! Recognizing comments
	  if(strRead(1:1)=='!'.or.strRead(1:1)==' '.or.strRead(1:1)=='	') cycle

! Deleting tabs
	  lenType=scan(strRead,'	')
	  do while (lenType>0)
	    strRead=strRead(:lenType-1)//strRead(lenType+1:)
		lenType=scan(strRead,'	')
	  enddo
	  lenType=scan(strRead,achar(13))
      if(lenType>0) strRead=strRead(:lenType-1)


            ! Reading parameter type
	  lenType=scan(strRead,':')
	  if(lenType/=0) then
	  strType=trim(strRead(:lenType-1))
	  strPar=strRead(lenType+1:)
	  strPar=adjustl(strPar)
	  select case (strType)
!		case('VegTypes')
!		  read(strPar,*) Veg_NumType
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! NOW VEG_NUMTYPE = 5 -fixed in VEG_PARAMS module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		case('DailyOutput')
!	        Veg_DailyOutput = (trim(adjustL(strPar)) == 'on')
!		case('OutputDetail')
!		  read(strPar,*) Veg_DetalizationLevel
		case default
		  print '(/,"STOP: Unknown parameter ''",a,"''",/)', trim(strType)
		  stop
	  end select
      end if
    end do
  close(2)

      do Subs=1, NumSubs
        Grp=SubsGroup(Subs)
        selectcase(Grp)
          case(HM)
 !           call Veg_HM_Props(Subs)
          case(Hg)
 !           call Veg_Hg_Props(Subs)
          case(POP)
           call Veg_POP_Props(Subs)
          case(Tracer)
!            call Veg_Tracer_Props(Subs)
        endselect
      enddo
     call LogFile('Reading physical and chemical properties (Veg)',0)

end subroutine Veg_ReadConfig

!========================================================================
! Read vegetation input data
!========================================================================
subroutine Veg_POPInput(Period)

	real ArrTmp(Veg_NumType)
	character(*) Period
	character(20) SurfTypeTmp(Veg_NumType)
        character(80) FileName
	integer i, j, k, NN(Veg_NumType)

  select case (Period)

    case('Initial')

        if(InitCond == 'cond') then
            VegInit = Veg_POPMass(1)          ! Initial mass in Veg in case of 'cond' parameter
            FallInit = Fall_POPMass(1)          ! Initial mass in Fall in case of 'cond' parameter
        elseif(InitCond == 'zero') then
            VegInit = 0.
            FallInit = 0.
        end if

!	Reading private LAI regression coefficients (Alpha)      
    filename='AlphaMonth.txt'
    fullName=trim(GeoPath)//trim(filename)
	open (8,file=trim(fullName),STATUS='OLD',action='READ')
	read(8,*) SurfTypeTmp(1:Veg_NumType)
	NN = 0
	do i = 1, Veg_NumType
	  do j = 1, Veg_NumType
	    if (trim(SurfTypeTmp(i)) == trim(SoilType(j))) NN(i) = j
	  end do
	end do
	do i = 1, Veg_NumType
	  if (NN(i) == 0) then
	    print*, 'Error in the file ', trim(fullName)
		stop
	  end if
	end do
!	do while(.not. EOF(8))
116	  read(8,*,end=115) k, ArrTmp(1:Veg_NumType)
	  do i = 1, Veg_NumType
	    Alpha(NN(i),k) = ArrTmp(i)
	  end do
          goto 116
!	enddo
115	close(8)

!-
	case('Monthly')
! Assigning initial values of LAIPriv

	    do i = 1, IMax
	      do j = 1, JMax
		    LAIPriv(i,j,:) = Alpha(:,Month)
		  end do
	    end do


      if (Month == 12) then
	    do i = 1,IMax
	      do j = 1, JMax
	        AlphaConst(i,j,:) = (Alpha(:,1) - Alpha(:,Month)) / (MonthDays(Month) * SecInDay)
	      end do
	    end do
	  else
	    do i = 1,IMax
	      do j = 1, JMax
	        AlphaConst(i,j,:) = (Alpha(:,Month + 1) - Alpha(:,Month)) / (MonthDays(Month) * &
		                                                                         & SecInDay)
	      end do
	    end do
	  end if

    case ('Daily')

  end select

end subroutine Veg_POPInput
   
!========================================================================
! Initialize vegetation data
!========================================================================
subroutine Veg_Initial

  integer i, j, L, Surf, s,FileStat, ioRes
  character(20) IDsurf1(50)
  real nps

	  do L = 1, Veg_NumType
	    do j = JMin, Jmax
	      do i = MinI(j), MaxI(j)
			  Veg_Frac(i,j,L) = Soil_Frac(i,j,L)
		  end do
	    end do
	  end do

end subroutine Veg_Initial

end module Veg_Input
#endif