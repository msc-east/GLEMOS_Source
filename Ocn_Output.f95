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
! Module of the ocean  information output
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifdef M_OCN
module Ocn_Output

  use GeneralParams
  use Ocn_Params
  use TextOutputProc

  implicit none

    character(800), private  :: fileName, fullName
    character(4), private   :: YearNum
    character(2), private   :: MonthNum, DayNum


  contains

!========================================================================
! Ocean output (main procedure)
!========================================================================
subroutine Ocn_OutputData(Period)

  character(*) Period

#ifdef DEBUG_MODE_OUTPUT
     print *, 'Enter Ocn_OutputData...'
#endif

  select case (Period)

    case('Initial')

    case ('Hourly')

    case ('6hourly')

    case ('Daily')

        Water_upl_concMonth = Water_upl_concMonth + Water_upl_concDay
        if (OcnOutFieldsDaily) call Ocn_WriteFields(Period)
        Water_upl_concDay = 0.

    case ('Monthly')

        Water_upl_concYear = Water_upl_concYear + Water_upl_concMonth
        if (OcnOutFieldsMonthly) call Ocn_WriteFields(Period)
        Water_upl_concMonth = 0.

    case ('Yearly')

        if (OcnOutFieldsYearly) call Ocn_WriteFields(Period)
        Water_upl_concYear = 0.

  end select

#ifdef DEBUG_MODE_OUTPUT
     print *, 'Exit Ocn_OutputData.'
#endif

end subroutine Ocn_OutputData


subroutine Ocn_WriteFields(prd)

    character(*), intent(in)    :: prd

    integer         i, j, s, n, ind
    character*80    header
    real            fld_2d(Imin:Imax,Jmin:Jmax,MaxForm), tp

    write(YearNum,'(I4)') Year
    write(MonthNum, '(i2.2)') Month
    write(DayNum, '(i2.2)') Day

    do s = 1, NumSubsMedia(Ocn)
        fld_2d = 0.
        select case (prd)
            case('Hourly')

            case('6hourly')

            case('Daily')
                header="Daily mean concentrations of "//trim(SubsID(s))//" in surface ocean layer (ng/L)"
                fullName=trim(OutPath)//trim(OcnDir)//trim(FieldsDir)//trim(DailyDir)//&
                        &trim(SubsID(s))//trim(ConcOcnName)//'('//YearNum//MonthNum//DayNum//').dat'
            ! copy fields
                do j = Jmin, Jmax
                    do i = minI(j), maxI(j)
                        if(Ocn_Frac(i,j) > 0.) then
                            fld_2d(i,j,s) = Water_upl_concDay(i,j,s)* 1.e9/SecInDay     ! kg/m3 => ng/L/sec
                        end if
                    end do
                end do

            case('Monthly')
                header="Monthly mean concentrations of "//trim(SubsID(s))//" in surface ocean layer (ng/L)"
                fullName=trim(OutPath)//trim(OcnDir)//trim(FieldsDir)//trim(MonthlyDir)//&
                        &trim(SubsID(s))//trim(ConcOcnName)//'('//YearNum//MonthNum//').dat'
                tp=min(SecInDay*real(MonthDays(Month)),timeCalc)
            ! copy fields
                do j = Jmin, Jmax
                    do i = minI(j), maxI(j)
                        if(Ocn_Frac(i,j) > 0.) then
                            fld_2d(i,j,s) = Water_upl_concMonth(i,j,s)* 1.e9 / tp     ! kg/m3 => ng/L/sec
                        end if
                    end do
                end do

            case('Yearly')
                header="Annual mean concentrations of "//trim(SubsID(s))//" in surface ocean layer (ng/L)"
                fullName=trim(OutPath)//trim(OcnDir)//trim(FieldsDir)//trim(YearlyDir)//&
                        &trim(SubsID(s))//trim(ConcOcnName)//'('//YearNum//').dat'
                tp=SecInDay*365.
            ! copy fields
                do j = Jmin, Jmax
                    do i = minI(j), maxI(j)
                        if(Ocn_Frac(i,j) > 0.) then
                            fld_2d(i,j,s) = Water_upl_concYear(i,j,s)*1.e9/tp   ! kg/m3 => ng/L/sec
                        end if
                    end do
                end do

        end select

    ! write to text file
        call WriteFieldTxt(trim(fullName),header,fld_2d,Ocn,s,1)

    end do

end subroutine Ocn_WriteFields

  end module Ocn_Output
#endif