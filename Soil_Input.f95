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
! Module of the soil  information input
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifdef M_SOIL
module Soil_Input

  use GeneralOutput
  use Soil_Params
  use GeneralParams
  use Geometry
  use Soil_POP_General
  use Balance

  implicit none

  character(800), private :: fileName, fullName

 contains

!========================================================================
! Reading soil input data
!========================================================================
 subroutine Soil_InputData(Period)

  character(*) Period
  integer g

#ifdef DEBUG_MODE
    print *, '+Entering Soil_InputData... ', Period
#endif

  select case (Period)

    case ('Initial')

	  call Soil_Initial

	  do g = 1, NumGroups
		  select case (SubsGrpLst(g))
	        case (POP)                          
		      call Soil_POPInput('Initial')
	        case default
		      continue
	      end select
	  end do

	case ('Steply')

      do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                          
		      call Soil_POPInput('Steply')
	        case default
		      continue
	      end select
      end do

	case ('Daily')

! Reading emissions data
      if (SoilEmissionStep=='daily') then
        call Soil_ReadEmis(da)
      end if

      call Soil_Daily

	  do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                         
		      call Soil_POPInput('Daily')
	        case default
		      continue
	      end select
      end do

	case ('Monthly')
! Reading emissions data
      if (SoilEmissionStep=='monthly') then
      call Soil_ReadEmis(mn)
      end if

      do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                          
		      call Soil_POPInput('Monthly')
	        case default
		      continue
	      end select
      end do

	case ('Yearly')

! Reading emissions data
      if (SoilEmissionStep=='yearly') then
      call Soil_ReadEmis(yr)
      end if

      do g = 1, NumGroups
	      select case (SubsGrpLst(g))
	        case (POP)                          
		      call Soil_POPInput('Yearly')
	        case default
		      continue
	      end select
      end do

  end select

#ifdef DEBUG_MODE
    print *, '+Exit Soil_InputData ', Period
#endif
end subroutine Soil_InputData

!========================================================================
! Read soil config from file
!========================================================================
subroutine Soil_ReadConfig

  integer FileStat, ioRes, NP
  integer lenType, k
  character(80) strRead, strPar
  character(20) strType
  character(160) fullName
  integer Subs,Grp

  filename='Soil'//trim(MediaConfName)//'.dat'
  fullName=trim(ConfigPath)//trim(filename)
  open(2, file = trim(fullName), status='old', iostat=FileStat, action='read')
	if(FileStat>0) then
            print '(/,"STOP in Soil_Config: Cannot open file ''",a,"''",/)', trim(fullName)
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
		case('Layers')
		  read(strPar,*) Soil_KMax, (Soil_dz(k), k = 1, Soil_KMax)
		case('SoilTypes')
		  read(strPar,*) Soil_NumType
!		case('DailyOutput')
!	      Soil_DailyOutput = (trim(adjustL(strPar)) == 'on')
		case('OutputDetail')
		  read(strPar,*) Soil_DetalizationLevel
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
 !           call Ocn_HM_Props(Subs)
          case(Hg)
 !           call Ocn_Hg_Props(Subs)
          case(POP)
           call Soil_POP_Props(Subs)
          case(Tracer)
!            call Ocn_Tracer_Props(Subs)
        endselect
      enddo
     call LogFile('Reading physical and chemical properties (soil)',0)

end subroutine Soil_ReadConfig

!========================================================================
! Reading soil input data (POPs)
!========================================================================
subroutine Soil_POPInput(Period)

  character(*) Period
  character*128 Strr                                  ! Added Sep 2016 V. Sh
  real HenryEff, Foc, Kow, Koc, Kdoc, TSurf, FocE, KocE, Depth
  integer NPOP, LL                                    ! Modified Sep 2016 V. Sh
  integer DayTime, i, j, k, FileStat, ioRes
	real Aver(Imin:Imax)
    integer XScal
    real*8 PartNorm

#ifdef DEBUG_MODE
    print *, '+Entering Soil_POPInput ', Period
#endif

! Calculating degradation parameters
  select case (Period)

    case('Initial')

        if(InitCond == 'cond') then
            SoilInit = Soil_POPMass(1)          ! Initial mass in soil in case of 'cond' parameter
        elseif(InitCond == 'zero') then
            SoilInit = 0.
        end if

        Depth = 0.
	  do k = 1, Soil_KMax
        Depth = Depth + Soil_dz(k)
	    if (Depth <= AverDepth) then
          do NPOP = 1, gSubsNum(GroupInd(POP))
            CDSoil(k,NPOP) = CDSoilTop(NPOP) + CDSoilBio(NPOP)
	      end do
	    else
          do NPOP = 1, gSubsNum(GroupInd(POP))
	        CDSoil(k,NPOP) = CDSoilBottom(NPOP) + CDSoilBio(NPOP)
	      end do
	    end if
      end do
      do k = 1, Soil_KMax
        do NPOP = 1, gSubsNum(GroupInd(POP))
          if (CDSoil(k,NPOP) /= 0.) Soil_TStep = min(0.5 / CDSoil(k,NPOP), Soil_TStep)
	    end do
      end do
! Correcting diffusion coefficients
	  do NPOP = 1, gSubsNum(GroupInd(POP))
	    DWater(NPOP) = DWater(NPOP) * wl ** 0.3333 / Poros / Poros
	    Dair(NPOP) = Dair(NPOP) * (Poros-wl) ** 0.3333 / Poros / Poros
	  end do
! Reading Foc data
      fileName=trim(FocName)//trim(GridCode)//'.dat'
      fullName=trim(GeoPath)//trim(fileName)
	  open(10, file=trim(fullName), status='old', iostat=FileStat, action='read')
	    if(FileStat>0) then
	      print '(/,"STOP: Cannot open FOC file ''",a,"''",/)', trim(fullName)
	      stop
	    endif
112	    read(10,'(a128)',end=111) Strr                  ! Modified Sep 2016 V. Sh
	    LL=scan(Strr,achar(13))
            if(LL>0) strr=strr(:LL-1)
            read(Strr,*) i, j, Foc
            if (Foc.lt.0.0001) Foc=0.0001
            FocEff(i, j) = Foc * Facc
	    Fdo(i,j) = FdoRel * rhosol * (1. - Poros) * Foc / 1000.
            goto 112
111	  close(10)

      do j=bJmin, bJmax
	    if(j==jSP.or.j==jNP) cycle
	    Xscal=Imax/maxI(j)
	    if(Xscal==1) cycle
	    Aver(Imin:Imax)=FocEff(Imin:Imax, j)
	    call GridAggreg(j,Xscal,Aver,1)
	    FocEff(minI(j):maxI(j), j)=Aver(minI(j):maxI(j))

 	    Aver(Imin:Imax)=Fdo(Imin:Imax, j)
	    call GridAggreg(j,Xscal,Aver,1)
	    Fdo(minI(j):maxI(j), j)=Aver(minI(j):maxI(j))
     enddo

	case ('Daily')

            do DayTime = 1, NumPer
!#ifdef DEBUG_MODE
!    print *, '++DayTime start ', DayTime
!#endif
                do j = JMin, JMax
                    do i = MinI(j), MaxI(j)
                        Tsurf = TempSurf(i, j, DayTime)
                        FocE = FocEff(i, j)
#ifdef DEBUG_MODE
    if(FocE == 0.) then
        print *, '++++ i,j,Tsurf, FocE = 0! ', i,j,Tsurf,FocE
    end if
#endif
                        do NPOP = 1, gSubsNum(GroupInd(POP))
                            Kow = Kow0(NPOP) * exp(KowT(NPOP) * (1 / Tsurf - 1 / T0))
                            Koc = Kow * 0.41 / 1000.

                            if (trim(SubsID(NPOP))=='BaP'.or.trim(SubsID(NPOP))=='BbF' &
                              & .or.trim(SubsID(NPOP))=='BkF' .or.trim(SubsID(NPOP))=='IP') then
                                Kdoc = Kow ** 0.98 * 10 ** (-0.39) ! PAHs
                            else
                                Kdoc = Kow ** 0.93 * 10 ** (-0.54) ! PCB and others
                            endif

                            KocE = Koc / (1 + Kdoc * Fdo(i,j))
!#ifdef DEBUG_MODE
!    print *, '++ Kow, Koc, Kdoc, KocE ', Kow, Koc, Kdoc, KocE
!#endif

                            do k = 1, Soil_KMax
                                HenryEff = HT0(NPOP) / (Runiv * Tsurf * (1 + Fdo(i,j) * Kdoc)) * exp(- HT(NPOP) &
                                           & * (1 / Tsurf - 1 / T0))
                                PartL(i, j, k, DayTime, NPOP) = rhosol * FocE * KocE + wl + (Poros - wl) &
                                                               & * HenryEff
                                PartG(i, j, k, DayTime, NPOP) = PartL(i, j, k, DayTime, NPOP) / HenryEff
                                PartS(i, j, k, DayTime, NPOP) = PartL(i, j, k, DayTime, NPOP)  &
                                                             & / (rhosol * FocE * KocE)

                                PartNorm = (Poros-wl)/PartG(i, j, k, DayTime, NPOP) &
                                &                + wl/PartL(i, j, k, DayTime, NPOP) &
                                &                + 1./PartS(i, j, k, DayTime, NPOP)

                                PartL(i, j, k, DayTime, NPOP) = PartL(i, j, k, DayTime, NPOP) * PartNorm
                                PartG(i, j, k, DayTime, NPOP) = PartG(i, j, k, DayTime, NPOP) * PartNorm
                                PartS(i, j, k, DayTime, NPOP) = PartS(i, j, k, DayTime, NPOP) * PartNorm
                            end do
                        end do ! NPOP
                    end do ! i
                end do ! j
!#ifdef DEBUG_MODE
!    print *, '++DayTime end ', DayTime
!#endif
            end do ! DayTime

  end select

#ifdef DEBUG_MODE
    print *, '+Exit Soil_POPInput ', Period
#endif
end subroutine Soil_POPInput
   
!========================================================================
! Initialize soil arrays
!========================================================================
subroutine Soil_Initial

  integer i, j, L, Surf, s,FileStat, ioRes
  character(20) IDsurf1(50)
  real nps


! Defining soil type fractions
        fileName=trim(LC2SoilTypeName)
        fullName=trim(GeoPath)//trim(LC2SoilTypeName)
	    open(10, file=trim(fullName), status = 'old', iostat=FileStat, action='read')
		if(FileStat>0) then
	      print '(/,"STOP: Cannot open file of soil type fractions ''",a,"''",/)', trim(fullName)
	      stop
	    endif
	    read(10, *)	(IDsurf1(i), i = 1, NumSurf)
	    do i = 1, NumSurf
	      if (trim(SurfType(i)) /= trim(IDSurf1(i))) then
	        print*, 'Files ', trim(LandName), ' and ', trim(LC2SoilTypeName), ' are incompatible'
		    stop
	      end if
	    end do
	    L=0
	    do while(.true.)
	      read(10,*, iostat = ioRes) SoilType(L), LC2SoilType(L, 1:NumSurf)
		  if (ioRes /= 0) exit
	      L = L + 1
	    enddo
    close(10)
	  if ((L-1) /= Soil_NumType) then
	    print*, 'Wrong number of soil types in '//trim(adjustL(LC2SoilTypeName))
		stop
	  end if

	  Soil_Frac = 0.
	  do L = 0, Soil_NumType
	    do j = JMin, Jmax
	      do i = MinI(j), MaxI(j)
		    do Surf = 1, NumSurf
			  Soil_Frac(i,j,L) = Soil_Frac(i,j,L) + LC2SoilType(L,Surf) * LandCover(i,j,Surf)
		    end do
		  end do
	    end do
	  end do

      Soil_Disp = 0.
      Soil_Z0 = 0.
      do s=1,NumSns
          do L = 0, Soil_NumType
              nps=0.
		     do Surf = 1, NumSurf
			   Soil_Disp(L,s) = Soil_Disp(L,s) + LC2SoilType(L,Surf) * dispLU(Surf,s)
			   Soil_Z0(L,s) = Soil_Z0(L,s) + LC2SoilType(L,Surf) * ZoLU(Surf,s)
               nps = nps+LC2SoilType(L,Surf)
             end do
             Soil_Disp(L,s)=Soil_Disp(L,s)/nps
             Soil_Z0(L,s)=Soil_Z0(L,s)/nps
          end do
      end do

      Soil_Height = 0.
          do L = 0, Soil_NumType
              nps=0.
		     do Surf = 1, NumSurf
			   Soil_Height(L) = Soil_Height(L) + LC2SoilType(L,Surf) * heightLU(Surf)
               nps = nps+LC2SoilType(L,Surf)
             end do
             Soil_Height(L)=Soil_Height(L)/nps
          end do


end subroutine Soil_Initial

!========================================================================
! Updating daily soil data 
!========================================================================
subroutine Soil_Daily

  integer DayTime, i, j

#ifdef DEBUG_MODE
    print *, '+Entering Soil_Daily'
#endif

  do DayTime = 1, NumPer
    do j = JMin, JMax
	  do i = MinI(j), MaxI(j)

		if ((PrecRainConv(i,j,1,DayTime)+PrecRainStrat(i,j,1,DayTime)) .GT. 0.) then
		  Inverse(i, j) =  InvMax
		  Ve(i, j) = (PrecRainConv(i,j,1,DayTime)+PrecRainStrat(i,j,1,DayTime))
		  QPrec(i, j) = QPrec(i, j) + (PrecRainConv(i,j,1,DayTime)+PrecRainStrat(i,j,1,DayTime))
		else
		  if (Inverse(i, j) == InvMax) DeltaQPrec(i, j) = QPrec(i, j) / InvMax
		  if (Inverse(i, j) /= 0) then
			Inverse(i, j) = Inverse(i, j) - 1
			QPrec(i, j) = QPrec(i, j) - DeltaQPrec(i, j)
			Ve(i, j) = - InvConst * DeltaQPrec(i, j)
		  else
			QPrec(i, j) = 0.
			Ve(i, j) = 0.
		  endif
		endif

      end do
    end do
  end do

#ifdef DEBUG_MODE
    print *, '+Exit Soil_Daily'
#endif
end subroutine Soil_Daily

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading soil emissions data
! Added 09.04.2012
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Soil_ReadEmis(tm)

    integer tm, Ind, Form, Subs, Grp
	character(100) Note

#ifdef DEBUG_MODE
    print *, '+Entering Soil_ReadEmis... ', tm
#endif

! Reading soil anthropogenic emission data
      do Ind=1, SoilEmisNum
        Form=ReadSoilInd(Ind,tm)
        if(Form==0) cycle
        Subs=FormSubs(Soil,Form)
        Grp=SubsGroup(Subs)
        selectcase(Grp)
          case(POP)
 	        call Soil_POP_ReadEmis(Subs,Ind,tm,Note)
        endselect

        call LogFile('Reading soil emissions data',0,Note)
      enddo

#ifdef DEBUG_MODE
    print *, '+Exit Soil_ReadEmis ', tm
#endif
end subroutine Soil_ReadEmis


end module Soil_Input
#endif