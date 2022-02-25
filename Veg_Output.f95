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
! Module of vegetation information output
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifdef M_VEG
module Veg_Output

    use GeneralParams
    use Veg_Params
    use TextOutputProc
    implicit none

    character(800), private  :: fileName, fullName
    character(4), private   :: YearNum
    character(2), private   :: MonthNum, DayNum

    contains

!========================================================================
! Vegetation output (main procedure)
!========================================================================
subroutine Veg_OutputData(Period)

  character(*) Period
  integer NSubs

#ifdef DEBUG_MODE_OUTPUT
     print *, 'Enter Veg_OutputData...'
#endif

  select case (Period)

    case('Initial')

    case ('Hourly')

    case ('6hourly')

    case ('Daily')
        Veg_ConcMonth = Veg_ConcMonth + Veg_ConcDay
        Fall_ConcMonth = Fall_ConcMonth + Fall_ConcDay
        if (VegOutFieldsDaily) call Veg_WriteFields('Daily')
        Veg_ConcDay = 0.
        Fall_ConcDay = 0.

    case ('Monthly')
        Veg_ConcYear = Veg_ConcYear + Veg_ConcMonth
        Fall_ConcYear = Fall_ConcYear + Fall_ConcMonth
        if (VegOutFieldsMonthly) call Veg_WriteFields('Monthly')
        Veg_ConcMonth = 0.
        Fall_ConcMonth = 0.

    case ('Yearly')
        if(VegOutFieldsYearly) call Veg_WriteFields('Yearly')
        Veg_ConcYear = 0.
        Fall_ConcYear = 0.

  end select

#ifdef DEBUG_MODE_OUTPUT
     print *, 'Exit Veg_OutputData.'
#endif

end subroutine Veg_OutputData


subroutine Veg_WriteFields(prd)

    character(*), intent(in)    :: prd

    integer         i, j, s, n, ind
    character*80    header
    real            fld_2d(Imin:Imax,Jmin:Jmax,MaxForm), tp
    real            Lcover(1:Veg_NumType), LCovVeg

    write(YearNum,'(I4)') Year
    write(MonthNum, '(i2.2)') Month
    write(DayNum, '(i2.2)') Day

    do s = 1, NumSubsMedia(Veg)

        select case (prd)
            case('Hourly')

            case('6hourly')

            case('Daily')
                header="Daily mean concentrations of "//trim(SubsID(s))//" in vegetation (ng/g)"
                fullName=trim(OutPath)//trim(VegDir)//trim(FieldsDir)//trim(DailyDir)//&
                        &trim(SubsID(s))//trim(ConcVegName)//'('//YearNum//MonthNum//DayNum//').dat'
            ! copy veg concentrations to temp array by forms and convert units
                do n = 1, sFormNum(Veg, s)
                    ind = sFormInd(Veg, s, n)           ! Current form index in the arrays of forms
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,n) = 0.
                            LCover(1:Veg_NumType)=Veg_Frac(i,j,1:Veg_NumType)
                            LCovVeg = sum(LCover(1:Veg_NumType))
                            if(LCovVeg >= LCovVegAccuracy) then
                                fld_2d(i,j,n) = sum(Veg_ConcDay(i,j,1:Veg_NumType,ind)*&
                                                &LCover(1:Veg_NumType))/LCovVeg*&
                                                &SLeaf2VegDens*1.e9/SecInDay   ! ng / g
                            end if
                        end do
                    end do
                end do

            case('Monthly')
                header="Monthly mean concentrations of "//trim(SubsID(s))//" in vegetation (ng/g)"
                fullName=trim(OutPath)//trim(VegDir)//trim(FieldsDir)//trim(MonthlyDir)//&
                        &trim(SubsID(s))//trim(ConcVegName)//'('//YearNum//MonthNum//').dat'
                tp=min(SecInDay*real(MonthDays(Month)),timeCalc)
            ! copy veg concentrations to temp array by forms and convert units
                do n = 1, sFormNum(Veg, s)
                    ind = sFormInd(Veg, s, n)           ! Current form index in the arrays of forms
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,n) = 0.
                            LCover(1:Veg_NumType)=Veg_Frac(i,j,1:Veg_NumType)
                            LCovVeg = sum(LCover(1:Veg_NumType))
                            if(LCovVeg >= LCovVegAccuracy) then
                                fld_2d(i,j,n) = sum(Veg_ConcMonth(i,j,1:Veg_NumType,ind)*&
                                                &LCover(1:Veg_NumType))/LCovVeg*&
                                                &SLeaf2VegDens*1.e9/tp   ! ng / g
                            end if
                        end do
                    end do
                end do

            case('Yearly')
                header="Annual mean concentrations of "//trim(SubsID(s))//" in vegetation (ng/g)"
                fullName=trim(OutPath)//trim(VegDir)//trim(FieldsDir)//trim(YearlyDir)//&
                        &trim(SubsID(s))//trim(ConcVegName)//'('//YearNum//').dat'
                tp=SecInDay*365.
            ! copy veg concentrations to temp array by forms and convert units
                do n = 1, sFormNum(Veg, s)
                    ind = sFormInd(Veg, s, n)           ! Current form index in the arrays of forms
                    do j = Jmin, Jmax
                        do i = minI(j), maxI(j)
                            fld_2d(i,j,n) = 0.
                            LCover(1:Veg_NumType)=Veg_Frac(i,j,1:Veg_NumType)
                            LCovVeg = sum(LCover(1:Veg_NumType))
                            if(LCovVeg >= LCovVegAccuracy) then
                                fld_2d(i,j,n) = sum(Veg_ConcYear(i,j,1:Veg_NumType,ind)*&
                                                &LCover(1:Veg_NumType))/LCovVeg*&
                                                &SLeaf2VegDens*1.e9/tp   ! ng / g
                            end if
                        end do
                    end do
                end do

        end select

    ! write to text file
        call WriteFieldTxt(trim(fullName),header,fld_2d,Veg,s)

    end do

end subroutine Veg_WriteFields



  end module Veg_Output
#endif