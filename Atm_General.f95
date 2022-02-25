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
! General module of the atmosphere medium
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module Atm_General

  use GeneralParams
  use Atm_HorizAdvect
  use Atm_VertAdvect
  use Atm_Diffusion
#ifdef G_HG
  use Atm_Hg_General
#endif
#ifdef G_HM
  use Atm_HM_Params
  use Atm_HM_General
#endif
#ifdef G_POP
  use Atm_POP_General
#endif
#ifdef G_AERO
  use Atm_AERO_General
#endif
#ifdef G_Tracer
  use Atm_Tracer_General
#endif

  implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine managing memory allocation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Memory(operation)

	character(*) operation

	selectcase(operation)
!---------------------------------------------------------------------------------
	case('Initial')
	  call Atm_MemAlloc
#ifdef G_HG
      call Atm_Hg_MemAlloc
#endif
#ifdef G_HM
      call Atm_HM_MemAlloc
#endif
#ifdef G_POP
      call Atm_POP_MemAlloc
#endif
#ifdef G_AERO
      call Atm_AERO_MemAlloc
#endif
#ifdef G_Tracer
      call Atm_Tracer_MemAlloc
#endif
!---------------------------------------------------------------------------------
	case('Final')
	  call Atm_MemDealloc
#ifdef G_HG
      call Atm_Hg_MemDealloc
#endif
#ifdef G_HM
      call Atm_HM_MemDealloc
#endif
#ifdef G_POP
      call Atm_POP_MemDealloc
#endif
#ifdef G_AERO
      call Atm_AERO_MemDealloc
#endif
#ifdef G_Tracer
      call Atm_Tracer_MemDealloc
#endif
	endselect

end subroutine Atm_Memory


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main subroutine of the atmosphere medium
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Process

    integer i, j, k, Subs
    real splRow(NumPer*2), d2FdX(NumPer*2), splVal

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_Process... '
#endif

! Temporal interpolation of air temperature
	do k=1, Atm_Kmax
      do j=bJmin, bJmax
	    do i=bminI(j), bmaxI(j)
          splRow(1:NumPer)=TempAir(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=TempAir(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2Tair(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          TairCurr(i,j,k)=splVal
        enddo
      enddo
	enddo
    call RTcurrent

! Transformation of volume cancentration to mass mixing ration
    call ConcToMixRatio

! Atmospheric emission calculation
	call Atm_Emission

! Boundary conditions
    call Atm_Boundary

! Atmospheric transport calculation
	call Atm_Transport

! Calculation of air density
    call AirDensity

! Transformation of mass mixing ration to volume cancentration
	call MixRatioToConc

! Pollutant specific processes
    do Subs=1, NumSubs
      selectcase(SubsGroup(Subs))
#ifdef G_HG
        case(HG)
          call Atm_Hg_Process(Subs)
#endif
#ifdef G_HM
        case(HM)
          call Atm_HM_Process(Subs)
#endif
#ifdef G_POP
        case(POP)
          call Atm_POP_Process(Subs)
#endif
#ifdef G_AERO
          case(AERO)
          call Atm_AERO_Process(Subs)
#endif
#ifdef G_Tracer
          case(TRACER)
          call Atm_Tracer_Process(Subs)
#endif
        case default
	      print '(/,"STOP: Unknown pollutant group ''",a,"''",/)', trim(SubsGroupID(SubsGroup(Subs)))
	      stop
	  endselect
    enddo    

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_Process'
#endif
    
end subroutine Atm_Process


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating atmospheric transport
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Transport

	integer i, j, k, Ind, Form, BegForm, FinForm, Iin, Nin, Nrec
	real(8) Pbeg(bImin:bImax,bJmin:bJmax,Atm_Kmax)			
	real(8) Pfin(bImin:bImax,bJmin:bJmax,Atm_Kmax), CourMaxSigma			
    real splRow(NumPer*2), d2FdX(NumPer*2), splVal

! Temporal interpolation of Px
    do j=bJmin, bJmax
	  do i=bminI(j), bmaxI(j)
        splRow(1:NumPer)=Px(i,j,1:NumPer,toDay)
        splRow(NumPer+1:NumPer*2)=Px(i,j,1:NumPer,toMor)
        d2FdX(1:NumPer+1)=d2PxdT(i,j,1:NumPer+1)
        splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
        PxNext(i,j)=splVal
      enddo
	enddo

! Calculation of diffusion cefficents
	call DiffCoefs

! Calculation of vertical velocity in sigma-coordinates
	Nin=1
	call Sigma_veloc(Nin)

! Advective horizontal transport
	do Ind=1, AtmTransNum
     Form=AtmTransInd(Ind)
	  call Horiz_Advect(Form)
	enddo

! Iterative vertical advection
#if A_VTYPE==1
! Dynamic calculation of vertical velocity
    Pbeg=Pcalc
	Nrec=0
	do while(Nin>1)
	  Nrec=Nrec+1
	  if(Nrec>5) then
		print '(/,"STOP: Unable to achieve an appropriate Courant number (Cmin=",f7.3,")",/)', CourMaxSigma
		stop
	  endif
	  do Iin=1, Nin-1
		call PressVert(Pbeg,Pfin)
	    do Ind=1, AtmTransNum
         Form=AtmTransInd(Ind)
		  call Vert_Advect(Form,Pbeg,Pfin)
		enddo
		Pbeg=Pfin
	  enddo
	  call CourNumSigma(Pbeg,PxNext,Nin)
	enddo
	do k=1, Atm_Kmax
	  Pfin(:,:,k)=PxNext(:,:)
	enddo
#elif A_VTYPE==2
! Static vertical velocity from met data
    Pbeg=Pcalc
    do Iin=1, Nin-1
      call PressVert(Pbeg,Pfin)
      do Ind=1, AtmTransNum
        Form=AtmTransInd(Ind)
        call Vert_Advect(Form,Pbeg,Pfin)
      enddo
    enddo
    call PressVert(Pbeg,Pfin)
#else
	print '(/,"STOP: Unknown method of vertical velocity calculation: ",i2,/)', A_VTYPE
	stop
#endif

! Final vertical advection and diffusion
    do Ind=1, AtmTransNum
      Form=AtmTransInd(Ind)
	  call Vert_Advect(Form,Pbeg,Pfin)
	  call Vert_Diffus(Form)
	enddo

	do j=bJmin, bJmax
	  PxCurr(bminI(j):bmaxI(j),j)=PxNext(bminI(j):bmaxI(j),j)
    enddo
! Temporal interpolation of wind velocity
	do k=1, Atm_Kmax
      do j=bJmin, bJmax
	    do i=bminI(j), bmaxI(j)
          splRow(1:NumPer)=Uwind(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Uwind(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2Uwind(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          Ucurr(i,j,k)=splVal

          splRow(1:NumPer)=Vwind(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Vwind(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2Vwind(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          Vcurr(i,j,k)=splVal

#if A_VTYPE==2
          splRow(1:NumPer)=Swind(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Swind(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2Swind(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          Scurr(i,j,k)=splVal
#endif
        enddo
      enddo
	enddo

! Calculation of next Courant numbers
	call CourNumsLong(Tstep(Atm),CourMax)
	call CourNumsLat(Tstep(Atm),CourMax)
    where(CourW>1) CourW=1
    where(CourE>1) CourE=1
    where(CourN>1) CourN=1
    where(CourS>1) CourS=1

end subroutine Atm_Transport


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating the dynamical time step in the atmosphere
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_TstepCalc(dTiter)

	integer j !, Niter
	real dTiter, CourMaxCurr, CourMaxNext

! Calculating local Courant numbers
! CourMax - maximal local Courant number
    if(Period<NumPer) then
      Ucurr(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)=Uwind(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,Period+1,toDay)
	  Vcurr(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)=Vwind(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,Period+1,toDay)
    else
      Ucurr(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)=Uwind(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1,toMor)
	  Vcurr(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)=Vwind(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,1,toMor)
    endif
    CourMaxNext=0.
	call CourNumsLong(dTinput,CourMaxNext)
	call CourNumsLat(dTinput,CourMaxNext)

	Ucurr(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)=Uwind(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,Period,toDay)
	Vcurr(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)=Vwind(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,Period,toDay)
    CourMaxCurr=0.
	call CourNumsLong(dTinput,CourMaxCurr)
	call CourNumsLat(dTinput,CourMaxCurr)

    CourMax=max(CourMaxCurr,CourMaxNext)

! Checking CFL condition and correcting dT (if necessary)
	if(CourMax<=1.) then
	  dTiter=dTinput
	else
	  dTiter=dTinput/CourMax
	endif

	do j=bJmin, bJmax
	  PxCurr(bminI(j):bmaxI(j),j)=Px(bminI(j):bmaxI(j),j,Period,toDay)
	enddo

! Definition of current air temperature
	TairCurr(bImin:bImax,bJmin:bJmax,1:Atm_Kmax)=TempAir(bImin:bImax,bJmin:bJmax,1:Atm_Kmax,Period,toDay)
    call RTcurrent

end subroutine Atm_TstepCalc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine adjusting courant numbers to calculated time step
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Adjust

	  CourW=CourW/dTinput*Tstep(Atm)
	  CourE=CourE/dTinput*Tstep(Atm)
	  CourN=CourN/dTinput*Tstep(Atm)
	  CourS=CourS/dTinput*Tstep(Atm)
	  CourMax=CourMax/dTinput*Tstep(Atm)

end subroutine Atm_Adjust


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating the model emission input
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Emission

	integer i, j, Form, Ind, Src, ht
    real(8) dT
    real(8) emiSum(Atm_Kmax), emis(Atm_Kmax,MaxMatr), trans(Atm_Kmax), q(Atm_Kmax)

    dT=Tstep(Atm)
	do j=Jmin, Jmax
	  do i=minI(j), maxI(j)

        do ht=1, HemisAnt
          trans(ht)=Veff(i,j,ht)*PxCurr(i,j)
        enddo

! Anthropogenic emission
		do Ind=1, AntEmisNum	
		  Form=AntEmisInd(Ind)
		  emiSum=0.
		  do Src=1, NumAnth
		    emis(1:HemisAnt,Src)=AntEmisFlux(i,j,Ind,Src,1:HemisAnt)*dT
		    emiSum(1:HemisAnt)=emiSum(1:HemisAnt)+emis(1:HemisAnt,Src)
		  enddo

            q(1:HemisAnt)=Atm_MixRatio(i,j,1:HemisAnt,Form)
            Atm_MixRatio(i,j,1:HemisAnt,Form)=q(1:HemisAnt)+emiSum(1:HemisAnt)/trans(1:HemisAnt)
            MassAtmAntEmis(Form)=MassAtmAntEmis(Form)+sum(emiSum(1:HemisAnt))
            AntEmisMonth(i,j)=AntEmisMonth(i,j)+sum(emiSum(1:HemisAnt))

!***************** Matrix calculations *****************
#if RTYPE==2
		  do ht=1, HemisAnt
            if(emiSum(ht)>0.) then
		      do Src=1, NumAnth
			    Atm_Contrib(i,j,ht,Form,Src)=(Atm_Contrib(i,j,ht,Form,Src)*q(ht)+&
					&emis(ht,Src)/trans(ht))/(q(ht)+emiSum(ht)/trans(ht))
			  enddo
			  do Src=NumAnth+1, NumSrc
				Atm_Contrib(i,j,ht,Form,Src)=Atm_Contrib(i,j,ht,Form,Src)*q(ht)/(q(ht)+emiSum(ht)/trans(ht))
			  enddo
		    endif
          enddo
#endif
!*******************************************************
		enddo

! Natural emission
		do Ind=1, NatEmisNum
		  Form=NatEmisInd(Ind)
		  emiSum(1)=0.
		  do Src=1, NumNat
		    emis(1,Src)=NatEmisFlux(i,j,Ind,Src)*dT
		    emiSum(1)=emiSum(1)+emis(1,Src)
		  enddo

                q(1)=Atm_MixRatio(i,j,1,Form)
                Atm_MixRatio(i,j,1,Form)=q(1)+emiSum(1)/trans(1)
                MassAtmNatEmis(Form)=MassAtmNatEmis(Form)+emiSum(1)
                NatEmisMonth(i,j)=NatEmisMonth(i,j)+emiSum(1)

!***************** Matrix calculations *****************
#if RTYPE==2
		  if(emiSum(1)>0.) then
            do Src=1, NumAnth
			  Atm_Contrib(i,j,1,Form,Src)=Atm_Contrib(i,j,1,Form,Src)*q(1)/(q(1)+emiSum(1)/trans(1))
			enddo
			do Src=1, NumNat
				Atm_Contrib(i,j,1,Form,NumAnth+Src)=(Atm_Contrib(i,j,1,Form,NumAnth+Src)*q(1)+&
					&emis(1,Src)/trans(1))/(q(1)+emiSum(1)/trans(1))
			enddo
			do Src=NumAnth+NumNat+1, NumSrc
			  Atm_Contrib(i,j,1,Form,Src)=Atm_Contrib(i,j,1,Form,Src)*q(1)/(q(1)+emiSum(1)/trans(1))
			enddo
		  endif
#endif
!*******************************************************
		enddo

! Reemission
		do Ind=1, ReEmisNum
		  Form=ReEmisInd(Ind)

		  emiSum(1)=0.
		  do Src=1, NumSrc
		    emis(1,Src)=ReEmisFlux(i,j,Ind,Src)      !*dT
		    emiSum(1)=emiSum(1)+emis(1,Src)
		  enddo

          q(1)=Atm_MixRatio(i,j,1,Form)
		  Atm_MixRatio(i,j,1,Form)=q(1)+emiSum(1)/trans(1)
		  MassReEmis(Form)=MassReEmis(Form)+emiSum(1)
          ReEmisMonth(i,j)=ReEmisMonth(i,j)+emiSum(1)

!***************** Matrix calculations *****************
#if RTYPE==2
		  if(emiSum(1)>0.) then
            if(reemisSRCmode==1) then
              emis(1,Src)=0.
              emis(1,Reem)=emiSum(1)
            endif

			do Src=1, NumSrc
			  Atm_Contrib(i,j,1,Form,Src)=(Atm_Contrib(i,j,1,Form,Src)*q(1)+&
					&emis(1,Src)/trans(1))/(q(1)+emiSum(1)/trans(1))
			enddo
		  endif
#endif
!*******************************************************
		enddo
	  enddo
	enddo

end subroutine Atm_Emission


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine setting boundary conditions
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Boundary

	integer i, j, k, Src, Ind, Form
	real(8) qW(MaxMatr), qE(MaxMatr), qS(MaxMatr), qN(MaxMatr), sumW, sumE, sumS, sumN, contB

    do Ind=1, AtmBndNum
      Form=AtmBndInd(Ind)

! Lateral boundary
#if (REGTYPE==2)
      do k=1, Atm_Kmax
        do j=bJmin, bJmax
          qW(1:NumSrc)=Atm_BoundW(j,k,Form,Period,1:NumSrc)
          qE(1:NumSrc)=Atm_BoundE(j,k,Form,Period,1:NumSrc)
          Atm_MixRatio(bminI(j),j,k,Form)=sum(qW(1:NumSrc))
          Atm_MixRatio(bmaxI(j),j,k,Form)=sum(qE(1:NumSrc))

!***************** Matrix calculations *****************
#if RTYPE==2
          sumW=max(0.,sum(qW(1:NumSrc)))
          sumE=max(0.,sum(qE(1:NumSrc)))
          do Src=1, NumSrc
            contB=0.
            if(sumW>0.) contB=qW(Src)/sumW
            Atm_Contrib(bminI(j),j,k,Form,Src)=contB

            contB=0.
            if(sumE>0.) contB=qE(Src)/sumE
            Atm_Contrib(bmaxI(j),j,k,Form,Src)=contB
          enddo
#endif
!*******************************************************
        enddo
        do i=bminI(bJmin), bmaxI(bJmin)
          qS(1:NumSrc)=Atm_BoundS(i,k,Form,Period,1:NumSrc)
          Atm_MixRatio(i,bJmin,k,Form)=sum(qS(1:NumSrc))

!***************** Matrix calculations *****************
#if RTYPE==2
          sumS=max(0.,sum(qS(1:NumSrc)))
          do Src=1, NumSrc
            contB=0.
            if(sumS>0.) contB=qS(Src)/sumS
            Atm_Contrib(i,bJmin,k,Form,Src)=contB
          enddo
#endif
!*******************************************************
        enddo
        do i=bminI(bJmax), bmaxI(bJmax)
          qN(1:NumSrc)=Atm_BoundN(i,k,Form,Period,1:NumSrc)
          Atm_MixRatio(i,bJmax,k,Form)=sum(qN(1:NumSrc))

!***************** Matrix calculations *****************
#if RTYPE==2
          sumN=max(0.,sum(qN(1:NumSrc)))
          do Src=1, NumSrc
            contB=0.
            if(sumN>0.) contB=qN(Src)/sumN
            Atm_Contrib(i,bJmax,k,Form,Src)=contB
          enddo
#endif
!*******************************************************
        enddo
      enddo
#endif

! Upper boundary
      do j=Jmin, Jmax
	    do i=minI(j), maxI(j)
          Atm_MixRatio(i,j,Atm_Kmax+1,Form)=0.

!***************** Matrix calculations *****************
#if RTYPE==2
		  Atm_Contrib(i,j,Atm_Kmax+1,Form,:)=0.
!		  if(Atm_MixRatio(i,j,Atm_Kmax+1,Form)>0.) Atm_Contrib(i,j,Atm_Kmax+1,Form,Bnd)=1.
#endif
!*******************************************************
	    enddo
	  enddo
    enddo

end subroutine Atm_Boundary

end module Atm_General