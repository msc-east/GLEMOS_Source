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
! Module of the vertical advection in sigma coordinates (version 4.0)
! Contains subroutines: 
!		Vert_Advect(Form)
!		CoefsSigma(Spec,Coef)
!		FluxSigma(i,j,Coef,Spec,FluxDn,FluxUp) 
!		CourNumSigma(dTime,CourM)
!		BottPolynomSigma
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

module Atm_VertAdvect 

  use GeneralParams
  use Atm_HorizAdvect

  implicit none

  real CourDn(Imin:Imax,Jmin:Jmax,Atm_Kmax+1), CourUp(Imin:Imax,Jmin:Jmax,Atm_Kmax) 
  real, private :: PolySig(3,0:Atm_Kmax+1,0:2)

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating vertical advection in sigma coordinates using 2nd order
! Bott scheme
! Version 2.0 - 2D Courant numbers, dynamic time step, accelerated
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Vert_Advect(Form,Pbeg,Pfin)

	integer Form, i, j, k, Src
	real(8) CoefSigma(Atm_Kmax+1,0:2)
	real(8) Jdown(Atm_Kmax+1), Jup(Atm_Kmax) !, Jbound(0:Imax,0:Jmax)
	real(8) SpecSigma(Atm_Kmax+1), Qspec
	real(8) aS(Atm_Kmax+1,MaxMatr), qCont
	real(8) Pbeg(bImin:bImax,bJmin:bJmax,Atm_Kmax)			
	real(8) Pfin(bImin:bImax,bJmin:bJmax,Atm_Kmax)			
    real(8) qCont1(MaxMatr)

#ifdef DEBUG_MODE
    print *, '<- Enter Vert_Advect'
#endif
!.................................................................................
!   Running in sigma-direction
!.................................................................................
	do j=Jmin, Jmax
	  do i=minI(j), maxI(j)

		SpecSigma(1:Atm_Kmax)=Atm_MixRatio(i,j,1:Atm_Kmax,Form)*Pbeg(i,j,1:Atm_Kmax)
		SpecSigma(Atm_Kmax+1)=Atm_MixRatio(i,j,Atm_Kmax+1,Form)*PxCurr(i,j)


!***************** Matrix calculations *****************
#if RTYPE==2
		aS(1:Atm_Kmax,1:NumSrc)=Atm_Contrib(i,j,1:Atm_Kmax,Form,1:NumSrc)
		aS(Atm_Kmax+1,1:NumSrc)=Atm_Contrib(i,j,Atm_Kmax+1,Form,1:NumSrc)
#endif
!*******************************************************

!	  Calculating Bott's polinomial coefficients
		call CoefsSigma(SpecSigma,CoefSigma)

!	  Calculation of fluxes through gridbox boundaries
		call FluxSigma(CourDn(i,j,1:Atm_Kmax+1),CourUp(i,j,1:Atm_Kmax),CoefSigma,SpecSigma,Jdown,Jup)

		Qspec=SpecSigma(1)-Jup(1)-Jdown(1)+dS(2)/dS(1)*Jdown(2) 
		Atm_MixRatio(i,j,1,Form)=max(Qspec/Pfin(i,j,1),0.) 

!***************** Matrix calculations *****************
#if RTYPE==2
		if(Qspec>0.) then
		  do Src=1, NumSrc
			qCont=(SpecSigma(1)-Jup(1)-Jdown(1))*aS(1,Src)&
					&+dS(2)/dS(1)*Jdown(2)*aS(2,Src)
!			Atm_Contrib(i,j,1,Form,Src)=qCont/Qspec
            qCont1(Src)=max(qCont/Qspec,0.)
		  enddo

! Normalizing the Atm_Contrib() array
        if(sum(qCont1(1:NumSrc))>0.) then
          Atm_Contrib(i,j,1,Form,1:NumSrc)=qCont1(1:NumSrc)/sum(qCont1(1:NumSrc))
        else
          Atm_Contrib(i,j,1,Form,1:NumSrc)=0.
        endif
		else
		  Atm_Contrib(i,j,1,Form,:)=0.
		endif
#endif
!*******************************************************

		do k=2, Atm_Kmax
		  Qspec=SpecSigma(k)-Jup(k)-Jdown(k)+dS(k-1)/dS(k)*Jup(k-1)&
				&+dS(k+1)/dS(k)*Jdown(k+1)
		  Atm_MixRatio(i,j,k,Form)=max(Qspec/Pfin(i,j,k),0.)

!***************** Matrix calculations *****************
#if RTYPE==2
		  if(Qspec>0.) then
			do Src=1, NumSrc
			  qCont=(SpecSigma(k)-Jup(k)-Jdown(k))*aS(k,Src)&
					&+dS(k-1)/dS(k)*Jup(k-1)*aS(k-1,Src)&
					&+dS(k+1)/dS(k)*Jdown(k+1)*aS(k+1,Src)
!			  Atm_Contrib(i,j,k,Form,Src)=qCont/Qspec
              qCont1(Src)=max(qCont/Qspec,0.)
			enddo

! Normalizing the Atm_Contrib() array
            if(sum(qCont1(1:NumSrc))>0.) then
              Atm_Contrib(i,j,k,Form,1:NumSrc)=qCont1(1:NumSrc)/sum(qCont1(1:NumSrc))
            else
              Atm_Contrib(i,j,k,Form,1:NumSrc)=0.
            endif
		  else
			Atm_Contrib(i,j,k,Form,:)=0.
		  endif
#endif
!*******************************************************
		enddo
	
		MassAtmUp(Form)=MassAtmUp(Form)+Jdown(Atm_Kmax+1)*Veff(i,j,Atm_Kmax+1) 
		MassAtmUp(Form)=MassAtmUp(Form)-Jup(Atm_Kmax)*Veff(i,j,Atm_Kmax)
      enddo
	enddo

#ifdef DEBUG_MODE
    print *, '<- Exit Vert_Advect'
#endif
	
end subroutine Vert_Advect


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating vertical velocity in sigma coordinates using 2nd order 
! Bott scheme 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Sigma_veloc(Nvert)		

	integer i, j, k, Nvert, s1, s2, n1, n2
	real CoefLong(bImin:bImax,0:2), CoefLat(Imin:Imax,bJmin:bJmax,0:2)
	real(8) Jwest(Imin:Imax+1), Jeast(Imin-1:Imax), Jnorth(Imin:Imax,bJmin:Jmax), Jsouth(Imin:Imax,Jmin:bJmax)
	real(8) SpecLong(bImin:bImax), SpecLat(Imin:Imax,bJmin:bJmax), SpecLatNew(Imin:Imax,bJmin:bJmax)  
	real(8) Phoriz(bImin:bImax,bJmin:bJmax,2), Qtmp, PoleSpecSum
	real(8) CourMaxSigma
	real(8) :: corr=1.e-12_8

! Calculation of horizoltal transport of the air
	do k=1, Atm_Kmax
      do j=bJmin, bJmax
		Phoriz(bminI(j):bmaxI(j),j,1)=PxCurr(bminI(j):bmaxI(j),j)*corr 
	  enddo

!.................................................................................
!   Running in zonal direction
!.................................................................................
	  do j=bJmin, bJmax
	    if(maxI(j)==1) cycle

		SpecLong(bminI(j):bmaxI(j))=Phoriz(bminI(j):bmaxI(j),j,1)

!     Calculating Bott's polinomial coefficients
		call CoefsLong(j,SpecLong,CoefLong)

!     Calculation of fluxes through gridbox boundaries
	    call FluxLong(j,k,CoefLong,SpecLong,Jwest,Jeast)

	    do i=minI(j), maxI(j)
		  Phoriz(i,j,2)=SpecLong(i)-Jwest(i)-Jeast(i)+Jeast(i-1)+Jwest(i+1)
	    enddo
	  enddo

!.................................................................................
!   Running in meridional direction
!.................................................................................

	  do j=bJmin+1, bJmax-1
	    SpecLat(minI(j):maxI(j),j)=Phoriz(minI(j):maxI(j),j,1) 
	    SpecLatNew(minI(j):maxI(j),j)=Phoriz(minI(j):maxI(j),j,2) 
	  enddo

	  if(maxI(bJmin)==1) then
		SpecLat(minI(jSP+1):maxI(jSP+1),jSP)=Phoriz(iSP,jSP,1)
		SpecLatNew(minI(jSP+1):maxI(jSP+1),jSP)=Phoriz(iSP,jSP,1)
	  else
		SpecLat(minI(bJmin):maxI(bJmin),bJmin)=Phoriz(minI(bJmin):maxI(bJmin),bJmin,1)
		SpecLatNew(minI(bJmin):maxI(bJmin),bJmin)=Phoriz(minI(bJmin):maxI(bJmin),bJmin,1)
	  endif

	  if(maxI(bJmax)==1) then
		SpecLat(minI(jNP-1):maxI(jNP-1),jNP)=Phoriz(iNP,jNP,1)
		SpecLatNew(minI(jNP-1):maxI(jNP-1),jNP)=Phoriz(iNP,jNP,1)
	  else
		SpecLat(minI(bJmax):maxI(bJmax),bJmax)=Phoriz(minI(bJmax):maxI(bJmax),bJmax,1)
		SpecLatNew(minI(bJmax):maxI(bJmax),bJmax)=Phoriz(minI(bJmax):maxI(bJmax),bJmax,1)
	  endif

!	Calculating Bott's polinomial coefficients
	  call CoefsLat(SpecLat,CoefLat)

!	Calculation of fluxes through gridbox boundaries
	  call FluxLat(k,CoefLat,SpecLatNew,Jnorth,Jsouth)

	  do j=bJmin+1, bJmax-1
		do i=minI(j), maxI(j)
		  s1=iS(i,j,1)
		  s2=iS(i,j,2)
		  n1=iN(i,j,1)
		  n2=iN(i,j,2)

		  if(n1==n2) then
			Qtmp=SpecLatNew(i,j)-dY/2.*(ctgdY2-tgY(j))*(Jnorth(i,j)-Jsouth(i,j+1))
		  else
			Qtmp=SpecLatNew(i,j)-dY/2.*(ctgdY2-tgY(j))*(Jnorth(n1,j)-Jsouth(n1,j+1)&
					&+Jnorth(n2,j)-Jsouth(n2,j+1))/2.
		  endif
		  if(s1==s2) then
			Qtmp=Qtmp+dY/2.*(ctgdY2+tgY(j))*(Jnorth(i,j-1)-Jsouth(i,j))
		  else
			Qtmp=Qtmp+dY/2.*(ctgdY2+tgY(j))*(Jnorth(s1,j-1)-Jsouth(s1,j)&
					&+Jnorth(s2,j-1)-Jsouth(s2,j))/2.
		  endif
		  Pcalc(i,j,k)=Qtmp/corr
	    enddo
	  enddo

!------------------ Global scale only -------------------
      if(maxI(bJmin)==1) then
        PoleSpecSum=0.
        do i=minI(jSP+1), maxI(jSP+1)
		  Qtmp=SpecLatNew(i,jSP)-dY/2.*(ctgdY2-tgY(jSP))*(Jnorth(i,jSP)-Jsouth(i,jSP+1))
          PoleSpecSum=PoleSpecSum+Qtmp
	    enddo
	    Pcalc(iSP,jSP,k)=PoleSpecSum/maxI(jSP+1)/corr
      endif
      if(maxI(bJmax)==1) then
        PoleSpecSum=0.
        do i=minI(jNP-1), maxI(jNP-1)
		  Qtmp=SpecLatNew(i,jNP)+dY/2.*(ctgdY2+tgY(j))*(Jnorth(i,jNP-1)-Jsouth(i,jNP))
		  PoleSpecSum=PoleSpecSum+Qtmp
        enddo
        Pcalc(iNP,jNP,k)=PoleSpecSum/maxI(jNP-1)/corr
      endif
!------------------ Global scale only -------------------
	enddo

! Calculation of vertical Courant numbers
#if A_VTYPE==1
	call CourNumSigma(Pcalc,PxNext,Nvert)
#elif A_VTYPE==2
    call CourNumSigma(Tstep(Atm),CourMaxSigma)

	if(CourMaxSigma>1.) then
	  Nvert=ceiling(CourMaxSigma)
	  CourDn=CourDn/Nvert
	  CourUp=CourUp/Nvert
	endif
#else
	print '(/,"STOP: Unknown method of vertical velocity calculation: ",i2,/)', A_VTYPE
	stop
#endif

end subroutine Sigma_veloc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating Bott polinomial coefficients in vertical direction
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CoefsSigma(Spec,Coef)

    integer k
	real(8) Spec(Atm_Kmax+1), SpDn, Sp0, SpUp
	real(8) Coef(Atm_Kmax+1,0:2)

	do k=1, Atm_Kmax+1
	  if(k==1) then				! Surface boundary conditions
	    SpDn=Spec(1)  
	  else
	    SpDn=Spec(k-1)
	  endif
	  Sp0=Spec(k)
	  if(k>Atm_Kmax) then
	    SpUp=Spec(Atm_Kmax+1)
	  else
	    SpUp=Spec(k+1)
	  endif

	  Coef(k,0)=PolySig(1,k,0)*SpDn+PolySig(2,k,0)*Sp0+PolySig(3,k,0)*SpUp
	  Coef(k,1)=PolySig(1,k,1)*SpDn+PolySig(2,k,1)*Sp0+PolySig(3,k,1)*SpUp
	  Coef(k,2)=PolySig(1,k,2)*SpDn+PolySig(2,k,2)*Sp0+PolySig(3,k,2)*SpUp
	enddo  
end subroutine CoefsSigma


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine FluxSigma(cDn,cUp,Coef,Spec,FluxDn,FluxUp)	

	integer k, l
	real cDn(Atm_Kmax+1), cUp(Atm_Kmax)
	real(8) FluxDn(Atm_Kmax+1), FluxUp(Atm_Kmax), fUp, fDn, Weight, Coef(Atm_Kmax+1,0:2)
	real(8) Spec(Atm_Kmax+1)
	real Ktwo, Ksign, Kc 

	do k=1, Atm_Kmax
	  fDn=0.
	  if(cDn(k)>Zero) then
	    Ktwo=2.
		Kc=1.-2.*cDn(k)
	    do l=0, 2
		  fDn=fDn+Coef(k,l)/(real(l)+1.)/Ktwo*(1.-Kc)
	      Ktwo=Ktwo*2.
		  Kc=Kc*(1.-2.*cDn(k))
	    enddo
	    fDn=max(0.,fDn)				! Flux limitation 
	  else
	    fDn=0.
	  endif
	  
	  fUp=0.
	  if(cUp(k)>Zero) then
	    Ktwo=2.
	    Ksign=1.
	    Kc=1.-2.*cUp(k)
		do l=0, 2
		  fUp=fUp+Coef(k,l)/(real(l)+1.)/Ktwo*Ksign*(1.-Kc)
	      Ktwo=Ktwo*2.
		  Ksign=-Ksign
		  Kc=Kc*(1.-2.*cUp(k))
	    enddo
	    fUp=max(0.,fUp)				! Flux limitation
	  else
	    fUp=0.
	  endif

	  Weight=min(1.,Spec(k)/(fDn+fUp+Zero))
	  FluxDn(k)=fDn*Weight
	  FluxUp(k)=fUp*Weight
	enddo

	fDn=0.
	if(cDn(Atm_Kmax+1)>Zero) then
 	  Ktwo=2.
	  Kc=1.-2.*cDn(Atm_Kmax+1)
	  do l=0, 2
	    fDn=fDn+Coef(Atm_Kmax+1,l)/(real(l)+1.)/Ktwo*(1.-Kc)
	    Ktwo=Ktwo*2
	    Kc=Kc*(1.-2.*cDn(Atm_Kmax+1))
	  enddo
	  fDn=max(0.,fDn)		! Flux limitation 
	else
	  fDn=0.
	endif
	FluxDn(Atm_Kmax+1)=fDn*min(1.,Spec(Atm_Kmax+1)/(fDn+Zero))

  end subroutine FluxSigma	


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating parameters of the Bott polinomials
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine BottPolynomSigma

	integer k
	real Kup, Kdn, Lup, Ldn, Delta, Stretch

	do k=1, Atm_Kmax+1
	  if(k==1) then
	    Stretch=dS(2)/dS(1)
		Kup=(1.+Stretch)/2.
	    Lup=3.+6.*Stretch+4.*Stretch*Stretch

		Kdn=1.
		Ldn=13.

	  elseif(k<Atm_Kmax+1) then
	    Stretch=dS(k+1)/dS(k)
		Kup=(1.+Stretch)/2.
	    Lup=3.+6.*Stretch+4.*Stretch*Stretch

		Stretch=dS(k-1)/dS(k)
		Kdn=(1.+Stretch)/2.
		Ldn=3.+6.*Stretch+4.*Stretch*Stretch

	  else
		Kup=1.
		Lup=13.

		Kdn=1.
		Ldn=13.
	  endif

	  Delta=Kdn*(Lup-1.)+Kup*(Ldn-1.)

	  PolySig(1,k,0)=-Kup/Delta
	  PolySig(2,k,0)=(Kdn*Lup+Kup*Ldn)/Delta
	  PolySig(3,k,0)=-Kdn/Delta
	  PolySig(1,k,1)=(Lup-1.)/Delta
	  PolySig(2,k,1)=(Ldn-Lup)/Delta
	  PolySig(3,k,1)=-(Ldn-1.)/Delta
	  PolySig(1,k,2)=12.*Kup/Delta
	  PolySig(2,k,2)=-12.*(Kup+Kdn)/Delta
	  PolySig(3,k,2)=12.*Kdn/Delta
	enddo

end subroutine BottPolynomSigma


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating local Courant numbers in vertical direction
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if A_VTYPE==1
subroutine CourNumSigma(Pbeg,Pfin,Nv)		

        integer i, j, k, Nv
        real(8) Pbeg(bImin:bImax,bJmin:bJmax,Atm_Kmax), Pfin(bImin:bImax,bJmin:bJmax)
        real(8) Pvert(Atm_Kmax+1), Pnew(Atm_Kmax)
        real(8) CoefSigma(Atm_Kmax+1,0:2), Jdown(Atm_Kmax+1), Jup(Atm_Kmax)
        real(8) cDn, cUp, CourM

        integer Info
        real(8) A, B, C, D
        logical run

        Nv=1
        run=.true.
        Info=1
      do while(run)
! Calculation of vertical Courant numbers
        CourM=0.
        do j=Jmin, Jmax
          do i=minI(j), maxI(j)

! Calculation of vertical fluxes
            Pvert(1:Atm_Kmax)=Pbeg(i,j,1:Atm_Kmax)
            Pvert(Atm_Kmax+1)=PxCurr(i,j)
            Pnew(1:Atm_Kmax)=Pvert(1:Atm_Kmax)+(Pfin(i,j)-Pvert(1:Atm_Kmax))/real(Nv,8)
            Jdown(1)=0.
            if(Pnew(1)-Pvert(1)<0.) then
              Jup(1)=Pvert(1)-Pnew(1)
              Jdown(2)=0.
            else
              Jup(1)=0.
              Jdown(2)=dS(1)/dS(2)*(Pnew(1)-Pvert(1))
            endif
            do k=2, Atm_Kmax
              if(Pnew(k)-Pvert(k)+Jdown(k)-dS(k-1)/dS(k)*Jup(k-1)<0.) then
                Jup(k)=Pvert(k)-Pnew(k)-Jdown(k)+dS(k-1)/dS(k)*Jup(k-1)
                Jdown(k+1)=0.
              else
                Jup(k)=0.
                Jdown(k+1)=dS(k)/dS(k+1)*(Pnew(k)-Pvert(k)+Jdown(k)-dS(k-1)/dS(k)*Jup(k-1))
              endif
            enddo

! Calculating Bott's polinomial coefficients
            call CoefsSigma(Pvert,CoefSigma)

! Calculation of Courant numbers
            if(Jup(1)>0.) then
              A=CoefSigma(1,2)/3.
              B=(CoefSigma(1,1)-CoefSigma(1,2))/2.
              C=CoefSigma(1,0)-CoefSigma(1,1)/2.+CoefSigma(1,2)/4.
              D=-Jup(1)
              cUp=PolinomRoot(A,B,C,D,Eps,Info)
              if(Info<0) then
                if(Nv>1) then
                  print '(/,"STOP: Accuracy of Sw calculation cannot be achieved (Info=",i3,")",/)', Info
                  print *, 'N=', Nv, 'i=', i, 'j=', j, 'k=', 1
                  stop 1
                endif
                exit
              endif
              if(cUp<0.) then
                print '(/,"STOP: Negative Courant number (Cup=",i3,")",/)', cUp
                stop
              endif
                CourUp(i,j,1)=cUp
              else
                CourUp(i,j,1)=0.
              endif		
              CourM=max(CourM,CourUp(i,j,1))
              do k=2, Atm_Kmax
                if(Jdown(k)>0.) then
                  A=CoefSigma(k,2)/3.
                  B=-(CoefSigma(k,1)+CoefSigma(k,2))/2.
                  C=CoefSigma(k,0)+CoefSigma(k,1)/2.+CoefSigma(k,2)/4.
                  D=-Jdown(k)
                  cDn=PolinomRoot(A,B,C,D,Eps,Info)
                  if(Info<0) then
                    if(Nv>1) then
                      print '(/,"STOP: Accuracy of Sw calculation cannot be achieved (Info=",i3,")",/)', Info
                      print *, 'N=', Nv, 'i=', i, '    j=', j, '    k=', k
                      stop 2
                    endif
                    exit
                  endif
                  if(cDn<0.) then
                    print '(/,"STOP: Negative Courant number (Cdn=",i3,")",/)', cDn
                    stop
                  endif
                  CourDn(i,j,k)=cDn
                else
                  CourDn(i,j,k)=0.
                endif		
                if(Jup(k)>0.) then
                  A=CoefSigma(k,2)/3.
                  B=(CoefSigma(k,1)-CoefSigma(k,2))/2.
                  C=CoefSigma(k,0)-CoefSigma(k,1)/2.+CoefSigma(k,2)/4.
                  D=-Jup(k)
                  cUp=PolinomRoot(A,B,C,D,Eps,Info)
                  if(Info<0) then
                    if(Nv>1) then
                      print '(/,"STOP: Accuracy of Sw calculation cannot be achieved (Info=",i3,")",/)', Info
                      print *, 'N=', Nv, 'i=', i, '    j=', j, '    k=', k
                      stop 3
                    endif
                    exit
                  endif
                  if(cUp<0.) then
                    print '(/,"STOP: Negative Courant number (Cup=",i3,")",/)', cUp
                    stop
                  endif
                  CourUp(i,j,k)=cUp
                else
                  CourUp(i,j,k)=0.
                endif		
                CourM=max(CourM,CourDn(i,j,k)+CourUp(i,j,k))
              enddo
              if(Info<0) exit
 
              if(Jdown(Atm_Kmax+1)>0.) then
                A=CoefSigma(Atm_Kmax+1,2)/3.
                B=-(CoefSigma(Atm_Kmax+1,1)+CoefSigma(Atm_Kmax+1,2))/2.
                C=CoefSigma(Atm_Kmax+1,0)+CoefSigma(Atm_Kmax+1,1)/2.+CoefSigma(Atm_Kmax+1,2)/4.
                D=-Jdown(Atm_Kmax+1)
                cDn=PolinomRoot(A,B,C,D,Eps,Info)
                if(Info<0) then
                  if(Nv>1) then
                    print '(/,"STOP: Accuracy of Sw calculation cannot be achieved (Info=",i3,")",/)', Info
                    print *, 'N=', Nv, 'i=', i, '    j=', j, '    k=', Atm_Kmax+1
                    stop 4
                  endif
                  exit
                endif
                if(cDn<0.) then
                  print '(/,"STOP: Negative Courant number (Cdn=",i3,")",/)', cDn
                  stop
                endif
                CourDn(i,j,Atm_Kmax+1)=cDn
              else
                CourDn(i,j,Atm_Kmax+1)=0.
              endif		
              CourM=max(CourM,CourDn(i,j,Atm_Kmax+1))
	  enddo
          if(Info<0) exit
	enddo

        if(Info<0) then
          Nv=Nv*ceiling(CourM)
        else
          run=.false.
        endif   
      enddo

	if(CourM>1.) then
	  Nv=Nv*ceiling(CourM)
	  CourDn=CourDn/Nv
	  CourUp=CourUp/Nv
	endif


end subroutine CourNumSigma

#elif A_VTYPE==2
subroutine CourNumSigma(dTime,CourM)

	integer i, j, k
	real dTime, Sw
    real(8) CourM

    do j=Jmin, Jmax
	  do i=minI(j), maxI(j)
        do k=1, Atm_Kmax			! Loop starts from 1, k-1=0
          Sw=Scurr(i,j,k-1)
		  if(Sw>0.) then
			CourDn(i,j,k)=Sw/dS(k)*dTime
		  else
			CourDn(i,j,k)=0.
		  endif

		  Sw=Scurr(i,j,k)
		  if(Sw<0.) then
			CourUp(i,j,k)=-Sw/dS(k)*dTime
		  else
			CourUp(i,j,k)=0.
		  endif
  		  CourM=max(CourM,CourUp(i,j,k)+CourDn(i,j,k))
		enddo

        Sw=Scurr(i,j,Atm_Kmax)
		if(Sw>0.) then
		  CourDn(i,j,Atm_Kmax+1)=Sw/dS(Atm_Kmax+1)*dTime
		else
		  CourDn(i,j,Atm_Kmax+1)=0.
		endif
		CourM=max(CourM,CourDn(i,j,Atm_Kmax+1))
	  enddo
	enddo

end subroutine CourNumSigma
#endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real(8) function PolinomRoot(a,b,c,d,eps,iFlag)

        integer it, iFlag
        real(8) a, b, c, d, x, x1, x2, F, F1, F2, xP, FP, dF, d2F, eps, Delta, x21, x22, F21, F22

        if(a==0..and.b==0.) then
          PolinomRoot=-d/c
          iFlag=1
          return
        endif

! Defuning the left-hand limit
        x1=0.
        F1=d

! Defuning the right-hand limit
        x2=0.
        F2=d
        Delta=b*b-3.*a*c
        if(a==0.) then
          x21=-c/2./b
          F21=a*x21*x21*x21+b*x21*x21+c*x21+d
          if(x21>0..and.x21<100..and.F1*F21<=0.) then
            x2=x21
            F2=F21
          endif
        elseif(Delta==0.) then
          x21=-b/3./a
          F21=a*x21*x21*x21+b*x21*x21+c*x21+d
          if(x21>0..and.x21<100..and.F1*F21<=0.) then
            x2=x21
            F2=F21
          endif
        elseif(Delta>0.) then
          x21=(-b-sqrt(Delta))/3./a
          F21=a*x21*x21*x21+b*x21*x21+c*x21+d
          x22=(-b+sqrt(Delta))/3./a
          F22=a*x22*x22*x22+b*x22*x22+c*x22+d

          if(x21>0..and.x21<100..and.F1*F21<=0.) then
            x2=x21
          endif
          if(x22>0..and.x22<100..and.F1*F22<=0..and.x22<x21) then
            x2=x22
          endif
          F2=a*x2*x2*x2+b*x2*x2+c*x2+d
        endif

        if(x2==0.) then
          x2=1.
          F2=a+b+c+d
          do while(F1*F2>0.)
            x2=x2+1.
            if(x2>100.) then
              iFlag=-1
              return
            endif
            F2=a*x2*x2*x2+b*x2*x2+c*x2+d
          enddo
        endif

! Searching for a flex point
        if(a/=0.) then
          xP=-b/3./a
          if(xP>x1.and.xP<x2) then
            FP=a*xP*xP*xP+b*xP*xP+c*xP+d
            if(FP*F1<0.) then
              x2=xP
              F2=FP
            else
              x1=xP
              F1=FP
            endif
          endif
        endif    

! Newton's method
        d2F=3.*a*(x1+x2)+2.*b
        if(d2F>0.) then
          F=F2
          x=x2
        else
          F=F1
          x=x1
        endif
        it=0
        do while(dabs(F/d)>eps)
          it=it+1
          if(it>10) then
            iFlag=-1
            return
          endif
          dF=3.*a*x*x+2.*b*x+c
          x=x-F/dF
          F=a*x*x*x+b*x*x+c*x+d
        enddo

        PolinomRoot=x
        iFlag=1

end function PolinomRoot


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PressVert(Pbeg,Pfin)

	integer i, j, k 
	real(8) Pbeg(bImin:bImax,bJmin:bJmax,Atm_Kmax)			
	real(8) Pfin(bImin:bImax,bJmin:bJmax,Atm_Kmax)			
	real(8) CoefSigma(Atm_Kmax+1,0:2)
	real(8) SpecSigma(Atm_Kmax+1), Jdown(Atm_Kmax+1), Jup(Atm_Kmax)

	do j=Jmin, Jmax
	  do i=minI(j), maxI(j)

		SpecSigma(1:Atm_Kmax)=Pbeg(i,j,1:Atm_Kmax)
		SpecSigma(Atm_Kmax+1)=PxCurr(i,j)

!  Calculating Bott's polinomial coefficients
		call CoefsSigma(SpecSigma,CoefSigma)

!  Calculation of fluxes through gridbox boundaries
		call FluxSigma(CourDn(i,j,:),CourUp(i,j,:),CoefSigma,SpecSigma,Jdown,Jup)

		Pfin(i,j,1)=SpecSigma(1)-Jup(1)-Jdown(1)+dS(2)/dS(1)*Jdown(2) 
		do k=2, Atm_Kmax
		  Pfin(i,j,k)=SpecSigma(k)-Jup(k)-Jdown(k)+dS(k-1)/dS(k)*Jup(k-1)&
			&+dS(k+1)/dS(k)*Jdown(k+1)
		enddo
	  enddo
	enddo

end subroutine PressVert


end module Atm_VertAdvect