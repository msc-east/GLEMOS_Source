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
! Diffusion module (version 4.0)
! Contains subroutines: 
!		Vert_Diffus(Ind)
!		DiffCoefs 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

module Atm_Diffusion 

  use GeneralParams
  use Atm_Params
  use Geometry

  implicit none

  real, private :: Dsigma(0:Imax,0:Jmax,Atm_Kmax)

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating vertical diffusion in sigma coordinates (implicit)
! Version 1.1 - final dry deposition flux as lower boundary condition
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Vert_Diffus(Form)

	integer Form, i, j, k, Src
	real(8) Asweep(Atm_Kmax+1), Bsweep(Atm_Kmax+1), Dsweep(Atm_Kmax), Spec(Atm_Kmax+1), Qspec
	real(8) Bcont(Atm_Kmax+1,MaxMatr)
	real(8) MixBnd, MixUp, MassUp, Dsig(0:Atm_Kmax)
	real(8) aD(Atm_Kmax+1,MaxMatr), qCont(MaxMatr), qSum
    real dT

    dT=Tstep(Atm)

	do j=Jmin, Jmax
	  do i=minI(j), maxI(j)

		Spec(1:Atm_Kmax+1)=Atm_MixRatio(i,j,1:Atm_Kmax+1,Form)*PxNext(i,j)
		Dsig(1:Atm_Kmax)=Dsigma(i,j,1:Atm_Kmax)

! Forward sweep
		Dsweep(1)=Dsig(1)+dS(1)/dT
		Asweep(2)=Dsig(1)/Dsweep(1)
		Bsweep(2)=dS(1)/dT*Spec(1)/Dsweep(1)
		do k=2, Atm_Kmax
		  Dsweep(k)=(1.-Asweep(k))*Dsig(k-1)+Dsig(k)+dS(k)/dT
		  Asweep(k+1)=Dsig(k)/Dsweep(k)
		  Bsweep(k+1)=(Dsig(k-1)*Bsweep(k)+dS(k)/dT*Spec(k))/Dsweep(k)
		enddo

!***************** Matrix calculations *****************
#if RTYPE==2
        aD(1:Atm_Kmax,1:NumSrc)=Atm_Contrib(i,j,1:Atm_Kmax,Form,1:NumSrc)
		aD(Atm_Kmax+1,1:NumSrc)=Atm_Contrib(i,j,Atm_Kmax+1,Form,1:NumSrc)
		Bcont(2,1:NumSrc)=dS(1)/dT*Spec(1)*aD(1,1:NumSrc)/Dsweep(1)
		do k=2, Atm_Kmax
		  Bcont(k+1,1:NumSrc)=(Dsig(k-1)*Bcont(k,1:NumSrc)+dS(k)/dT*Spec(k)*aD(k,1:NumSrc))/Dsweep(k)
		enddo
#endif
!*******************************************************

! Backward sweep
		Qspec=Atm_MixRatio(i,j,Atm_Kmax+1,Form)*PxNext(i,j)

!***************** Matrix calculations *****************
#if RTYPE==2
		qCont(1:NumSrc)=Qspec*aD(Atm_Kmax+1,1:NumSrc)
#endif
!*******************************************************

        do k=Atm_Kmax, 1, -1
		  Qspec=Asweep(k+1)*Qspec+Bsweep(k+1)
		  Atm_MixRatio(i,j,k,Form)=Qspec/PxNext(i,j)

!***************** Matrix calculations *****************
#if RTYPE==2
		  qSum=0.
		  do Src=1, NumSrc
		    qCont(Src)=Asweep(k+1)*qCont(Src)+Bcont(k+1,Src)
		    qSum=qSum+qCont(Src)
		  enddo
		  if(qSum>0.) then
			Atm_Contrib(i,j,k,Form,1:NumSrc)=qCont(1:NumSrc)/qSum
		  else
		    Atm_Contrib(i,j,k,Form,1:NumSrc)=0.
		  endif
#endif
!*******************************************************
		enddo

		MixBnd=Atm_MixRatio(i,j,Atm_Kmax+1,Form)
		MixUp=Atm_MixRatio(i,j,Atm_Kmax,Form)
        MassUp=dT/dS(Atm_Kmax)*Dsigma(i,j,Atm_Kmax)*(MixUp-MixBnd)*&
				&Veff(i,j,Atm_Kmax)*PxNext(i,j)
		MassAtmUp(Form)=MassAtmUp(Form)-MassUp
	  enddo
	enddo

end subroutine Vert_Diffus


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating coefficients of the diffusion sweep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DiffCoefs

	integer i, j, k
	real Dcoef1(Atm_Kmax), Dcoef2(Atm_Kmax)
    real splRow(NumPer*2), d2FdX(NumPer*2), splVal, KzCurr

	do j=Jmin, Jmax
	  do i=minI(j), maxI(j)
		do k=1, Atm_Kmax
        Dcoef1(k)=(Ggrav/RTcurr(i,j,k)*(Sigma(k)+Ptop/PxCurr(i,j)))**2
        enddo

		do k=Atm_Kmax-1, 1, -1
          splRow(1:NumPer)=Ksigma(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Ksigma(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2Ksigma(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          KzCurr=max(splVal,0.)

          Dcoef2(k)=2.*(Dcoef1(k)*dS(k+1)+Dcoef1(k+1)*dS(k))/(dS(k)+dS(k+1))/(dS(k)+dS(k+1))
		  Dsigma(i,j,k)=KzCurr*Dcoef2(k)
		enddo	
		Dsigma(i,j,Atm_Kmax)=0.
	  enddo	
	enddo

end subroutine DiffCoefs

end module Atm_Diffusion