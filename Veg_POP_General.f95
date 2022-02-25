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
! Module containing POP specific procedures
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifdef G_POP
module Veg_POP_General

  use GeneralParams
  use Geometry
  use Atm_POP_Params
  use Veg_Params
  use Veg_POP_Params
  use Exch_Params

  implicit none

  character(800), private :: fileName, fullName

  contains


!========================================================================
! Vegetation related POP processes
!========================================================================
  subroutine Veg_POPProcess(NSubs, dTime)

  integer NSubs, i, j, k, f, L
  real(8) CTot
  real dTime
  integer Nair,iair
  real dTair

#ifdef DEBUG_MODE
    print *, '+Entering Veg_POPProcess... '
#endif

    call Veg_POPdegradation(NSubs, dTime)

! Counters for averaging
    do L = 1, Veg_NumType
        do j = JMin, JMax
          do i = MinI(j), MaxI(j)
!		CTot=0.
		    do f = 1, NumForm(Veg)
		    Veg_ConcDay(i,j,L,sFormInd(Veg,NSubs,f)) = Veg_ConcDay(i,j,L,sFormInd(Veg,NSubs,f)) &
                    & + Veg_Conc(i,j,L,sFormInd(Veg,NSubs,f)) * dTime
			end do
!		    Veg_ConcDay(i,j,L,sFormInd(Veg,NSubs,f)) = Veg_ConcDay(i,j,L,sFormInd(Veg,NSubs,f)) + CTot * dTime
		    Fall_ConcDay(i,j,L,gSubsMediaInd(Veg,NSubs))= Fall_ConcDay(i,j,L,gSubsMediaInd(Veg,NSubs)) &
                    & + Fall_Conc(i,j,L,gSubsMediaInd(Veg,NSubs))  * dTime
	    end do
      end do
    end do

#ifdef DEBUG_MODE
    print *, '+Exit Veg_POPProcess'
#endif
end subroutine Veg_POPProcess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine defining physical and chemical properties of a POP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Veg_POP_Props(Subs)

	integer i, FileStat, lenType, Subs, Form, Limit, Ind
	character(80) strRead, strPar
	character(30) strType, FormType(MaxForm)
	integer mon, ios

	fileName='Veg'//trim(PropName)//trim(SubsID(Subs))//'.dat'
	fullName=trim(PropPath)//fileName
	open(2, file=trim(fullName), status='old', iostat=FileStat, action='read')
	if(FileStat>0) then
	  print '(/,"STOP: Cannot open file of veg POP properties ''",a,"''",/)', trim(fullName)
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
		do i=1, sFormNum(Veg,Subs)
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
		  case('Forms number')
			read(strPar,'(i2)') sFormNum(Veg,Subs)
                        FormSubs(Veg,NumForm(Veg)+1:NumForm(Veg)+sFormNum(Veg,Subs))=Subs
			if (sFormNum(Veg,Subs).gt.0) then	
			NumSubsMedia(Veg)=NumSubsMedia(Veg)+1
			gSubsMediaInd(Veg,Subs)=NumSubsMedia(Veg)
			end if
!------------------------------------------------------------------
		  case('Forms')
			read(strPar,*) (FormType(i), i=1, sFormNum(Veg,Subs))
!------------------------------------------------------------------
		  case('Form ID')
            FormID(Veg,Subs,Form)=strPar
            sFormInd(Veg,Subs,Form)=NumForm(Veg)+1
            FormSubs(Veg,sFormInd(Veg,Subs,Form))=Subs
            NumForm(Veg)=NumForm(Veg)+1
!------------------------------------------------------------------
		  case('VegDegr')
             read(strPar,*) CDVeg
!------------------------------------------------------------------
		  case('FallDegr')
             read(strPar,*) CDFall
!------------------------------------------------------------------
		  case('MolarVol')
             read(strPar,*) Vmol
!------------------------------------------------------------------
	endselect

	  endif
	enddo

	close(2)

end subroutine Veg_POP_Props

!========================================================================
! Calculation of time step for Veg
!========================================================================
subroutine Veg_POPTStepCalc(NSubs,dTiter)

  integer NSubs
  real dTiter

dTiter=Tstep(Soil)

end subroutine Veg_POPTStepCalc

!========================================================================
! Degradation of a POP in vegetation
!========================================================================
subroutine Veg_POPdegradation(NSubs, dTime)

    integer NSubs
    integer i, j, f, L
    real DTime
    real*8 SPM1,SPM2

    SPM1 = Veg_POPMass(NSubs)
    SPM1 = SPM1 + Fall_POPMass(NSubs)

    do L = 1, Veg_NumType
        do j = JMin, JMax
            do i = MinI(j), MaxI(j)
                do f = 1, NumForm(Veg)
                    Veg_Conc(i,j,L,sFormInd(Veg,NSubs,f))=Veg_Conc(i,j,L,sFormInd(Veg,NSubs,f))*dble(1.-CDVeg*dTime)
                end do
                Fall_Conc(i,j,L,gSubsMediaInd(Veg,NSubs))=Fall_Conc(i,j,L,gSubsMediaInd(Veg,NSubs))*dble(1.-CDFall*dTime)
            end do
        end do
    end do

    SPM2 = Veg_POPMass(NSubs)
    SPM2 = SPM2 + Fall_POPMass(NSubs)
    Veg_Degr = Veg_Degr + SPM1 - SPM2

end subroutine Veg_POPdegradation

!========================================================================
! Mass in vegetation
!========================================================================
real(8) function Veg_POPMass(NSubs)

  integer i, j, k, f, L, NSubs

    Veg_POPMass = 0.
    do L = 1, Veg_NumType
        do j =  JMin, JMax
          do i = MinI(j), MaxI(j)
		do f = 1, NumForm(Veg)
 	        Veg_POPMass = Veg_POPMass + Veg_Conc(i,j,L,sFormInd(Veg,NSubs,f))* &
				& dble(MeshArea(i,j)) * Veg_Frac(i,j,L)
		end do
          end do
        end do
    end do
end function Veg_POPMass

real(8) function Veg_POPMassByForm(NSubs,f)
    integer, intent(in) :: f
    integer i, j, k, L, NSubs

    Veg_POPMassByForm = 0.
    do L = 1, Veg_NumType
        do j =  JMin, JMax
          do i = MinI(j), MaxI(j)
 	        Veg_POPMassByForm = Veg_POPMassByForm + Veg_Conc(i,j,L,sFormInd(Veg,NSubs,f))* &
				& dble(MeshArea(i,j)) * Veg_Frac(i,j,L)
          end do
        end do
    end do
end function Veg_POPMassByForm

!========================================================================
! Mass in fall
!========================================================================
real(8) function Fall_POPMass(NSubs)

  integer i, j, k, f, L, NSubs

    Fall_POPMass = 0.
    do L = 1, Veg_NumType
        do j =  JMin, JMax
          do i = MinI(j), MaxI(j)
 	        Fall_POPMass = Fall_POPMass + Fall_Conc(i,j,L,gSubsMediaInd(Veg,NSubs))* &
				& dble(MeshArea(i,j)) * Veg_Frac(i,j,L)
          end do
        end do
    end do

end function Fall_POPMass


end module Veg_POP_General
#endif