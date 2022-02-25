!*************************************************************************!
!*  Copyright (C) Meteorological Synthesizing Centre - East of EMEP, 2021 
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
! Global EMEP Multi-media Modelling System (GLEMOS)
! (c) EMEP/MSC-E, 2021
! Version 2.2.2
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

program MainProgram

    use GeneralParams
    use GeneralProcs

    implicit none

    integer YrBeg, YrFin, MonBeg, MonFin, DayBeg, DayFin

    call Initial
    call Memory('Initial')
    call InputData('Initial')
    call OutputData('Initial')

    ! Year calculation loop
    call LogFile('*** Start of the calculation cycle ***')
    call CalcLimits(yr,YrBeg,YrFin)
    
    do Year=YrBeg, YrFin

      ! Reading yearly input data
      call InputData('Yearly')
      ! Month calculation loop
      call CalcLimits(mn,MonBeg,MonFin)

      do Month=MonBeg, MonFin

        ! Reading monthly input data
        call InputData('Monthly')
        ! Day calculation loop
        call CalcLimits(da,DayBeg,DayFin)

        do Day=DayBeg, DayFin

          call PrintScr('Daily')
          ! Reading daily input data
          call InputData('Daily')

            ! Meteo period calculation loop 

            do Period=1, NumPer
              
              ! Reading input data per meteo period
              call InputData('6hourly')
              ! Calculation of time steps in different media
              call MediaTimeSteps
              ! Time integration procedure
              call TimeIntegration(NumMed)
              call OutputData('6hourly')

            enddo

            call OutputData('Daily')
          
          enddo

          call OutputData('Monthly')
      
        enddo

      call OutputData('Yearly')

    enddo

    call OutputData('Final')
    call Memory('Final')

end program MainProgram
