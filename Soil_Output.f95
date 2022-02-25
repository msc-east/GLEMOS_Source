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
! Module of the soil  information output
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifdef M_SOIL
module Soil_Output

    use GeneralParams
    use Soil_Params
    use TextOutputProc
    implicit none

    character(800), private  :: fileName, fullName
    character(4), private   :: YearNum
    character(2), private   :: MonthNum, DayNum

  contains

!========================================================================
! Soil output (main procedure)
!========================================================================
subroutine Soil_OutputData(Period)

  character(*) Period

#ifdef DEBUG_MODE_OUTPUT
     print *, 'Enter Soil_OutputData...'
#endif

  select case (Period)

    case('Initial')

    case ('Hourly')

    case ('6hourly')

    case ('Daily')
        Soil_ConcMonth = Soil_ConcMonth + Soil_ConcDay
        if (SoilOutFieldsDaily) call Soil_WriteFields('Daily')
        Soil_ConcDay = 0.

    case ('Monthly')
        Soil_ConcYear = Soil_ConcYear + Soil_ConcMonth
        if (SoilOutFieldsMonthly) call Soil_WriteFields('Monthly')
        Soil_ConcMonth = 0.

    case ('Yearly')
        if(SoilOutFieldsYearly) call Soil_WriteFields('Yearly')
        Soil_ConcYear = 0.

  end select

#ifdef DEBUG_MODE_OUTPUT
     print *, 'Exit Soil_OutputData.'
#endif

end subroutine Soil_OutputData


subroutine Soil_WriteFields(prd)

    character(*), intent(in)    :: prd

    integer         i, j, s, n, ind
    character*80    header
    real            fld_2d(Imin:Imax,Jmin:Jmax,MaxForm), tp, aver

    write(YearNum,'(I4)') Year
    write(MonthNum, '(i2.2)') Month
    write(DayNum, '(i2.2)') Day

    do s = 1, NumSubsMedia(Soil)

        select case (prd)
            case('Hourly')

            case('6hourly')

            case('Daily')
                header="Daily mean concentrations of "//trim(SubsID(s))//" in upper soil layer (ng/g)"
                fullName=trim(OutPath)//trim(SoilDir)//trim(FieldsDir)//trim(DailyDir)//&
                        &trim(SubsID(s))//trim(ConcSoilName)//'('//YearNum//MonthNum//DayNum//').dat'
            ! averaging soil concentrations for specified depth
                do j = Jmin, Jmax
                    do i = minI(j), maxI(j)
                        call UpperSoilAverageConc(i,j,aver, Soil_ConcDay(i,j,1:Soil_KMax,1:Soil_NumType,s))
                        fld_2d(i,j,s) = aver * 1.e9 / rhoSol / SecInDay  ! ng/g
                    end do
                end do

            case('Monthly')
                header="Monthly mean concentrations of "//trim(SubsID(s))//" in upper soil layer (ng/g)"
                fullName=trim(OutPath)//trim(SoilDir)//trim(FieldsDir)//trim(MonthlyDir)//&
                        &trim(SubsID(s))//trim(ConcSoilName)//'('//YearNum//MonthNum//').dat'
                tp=min(SecInDay*real(MonthDays(Month)),timeCalc)
            ! averaging soil concentrations for specified depth
                do j = Jmin, Jmax
                    do i = minI(j), maxI(j)
                        call UpperSoilAverageConc(i,j,aver, Soil_ConcMonth(i,j,1:Soil_KMax,1:Soil_NumType,s))
                        fld_2d(i,j,s) = aver * 1.e9 / rhoSol / tp  ! ng/g
                    end do
                end do

            case('Yearly')
                header="Annual mean concentrations of "//trim(SubsID(s))//" in upper soil layer (ng/g)"
                fullName=trim(OutPath)//trim(SoilDir)//trim(FieldsDir)//trim(YearlyDir)//&
                        &trim(SubsID(s))//trim(ConcSoilName)//'('//YearNum//').dat'
                tp=SecInDay*365.
            ! averaging soil concentrations for specified depth
                do j = Jmin, Jmax
                    do i = minI(j), maxI(j)
                        call UpperSoilAverageConc(i,j,aver, Soil_ConcYear(i,j,1:Soil_KMax,1:Soil_NumType,s))
                        fld_2d(i,j,s) = aver * 1.e9 / rhoSol / tp  ! ng/g
                    end do
                end do

        end select

    ! write to text file
        call WriteFieldTxt(trim(fullName),header,fld_2d,Soil,s,1)

    end do

end subroutine Soil_WriteFields


subroutine UpperSoilAverageConc(i,j,CellAver, arr)

    real, intent(inout)     ::  CellAver
    real(8), intent(in)     ::  arr(1:Soil_KMax,1:Soil_NumType)
    integer, intent(in)     ::  i, j

    real(8)         AverConc(1:Soil_NumType)
    real            Depth, SoilF
    integer         k, l


    AverConc = 0.
    CellAver = 0.

    SoilF = 0.
    do l = 1, Soil_NumType
        SoilF = SoilF + Soil_Frac(i,j,l)
    end do

    if (SoilF >= 0.01) then
        do l = 1, Soil_NumType
            Depth = 0.
            do k = 1, Soil_KMax
                if (Depth + Soil_dz(k) < AverDepth) then
                    AverConc(l) = AverConc(l) + arr(k,l) * Soil_dz(k)
                    Depth = Depth + Soil_dz(k)
                else
                    AverConc(l) = AverConc(l) + arr(k,l) * (AverDepth - Depth)
                    exit
                end if
            end do
            AverConc(l) = AverConc(l) / AverDepth
            CellAver = CellAver + AverConc(l) * Soil_Frac(i,j,l)
        end do
        CellAver = CellAver / SoilF
    end if

end subroutine UpperSoilAverageConc


  end module Soil_Output
#endif