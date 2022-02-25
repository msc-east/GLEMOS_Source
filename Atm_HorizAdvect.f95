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
! Module of the horizontal advection
! Contains subroutines: 
!		Horiz_Advect(Form)
!		CoefsLong(Spec,Coef)
!		CoefsLat(Spec,Coef)
!		FluxLong(j,k,Coef,Spec,FluxW,FluxE)
!		FluxLat(i,k,Coef,Spec,FluxN,FluxS)
!		BottPolynomLat
!		CourNumsLong(dTime,CourM)
!		CourNumsLat(dTime,CourM)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

module Atm_HorizAdvect 


  use GeneralParams
  use Atm_Params

  implicit none

  real :: CourW(Imin:bImax,bJmin:bJmax,Atm_Kmax), CourE(bImin:Imax,bJmin:bJmax,Atm_Kmax), CourN(Imin:Imax,bJmin:Jmax,Atm_Kmax),&
		&CourS(Imin:Imax,Jmin:bJmax,Atm_Kmax), CourMax
  real, private :: PolyLat(3,bJmin:bJmax,0:2)

contains
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating horizontal advection in spherical coordinates using 3nd 
! order Bott scheme
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Horiz_Advect(Form)

    integer :: Form, i, j, k, Src, s1, s2, n1, n2
    real CoefLong(bImin:bImax,0:2), CoefLat(Imin:Imax,bJmin:bJmax,0:2)
    real(8) Jwest(Imin:Imax+1), Jeast(Imin-1:Imax), Jnorth(Imin:Imax,bJmin:Jmax), Jsouth(Imin:Imax,Jmin:bJmax)
    real(8) SpecLong(bImin:bImax), SpecLat(Imin:Imax,bJmin:bJmax), SpecLatNew(Imin:Imax,bJmin:bJmax)
    real(8) Qspec(bImin:bImax,bJmin:bJmax,2), Qtemp, PoleSpecSum, PoleContSum
    real(8) aLn(Imin-1:Imax+1,MaxMatr), aLt(Imin:Imax,bJmin:bJmax,MaxMatr), Contr(Imin:Imax,bJmin:bJmax,MaxMatr), qCont
real(8) qCont1(MaxMatr)

    do k=1, Atm_Kmax
      do j=bJmin, bJmax
        Qspec(bminI(j):bmaxI(j),j,1)=Atm_MixRatio(bminI(j):bmaxI(j),j,k,Form)*PxCurr(bminI(j):bmaxI(j),j) 
      enddo

!.................................................................................
!   Running in zonal direction
!.................................................................................
	  do j=bJmin, bJmax
	    if(maxI(j)==1) cycle

		SpecLong(bminI(j):bmaxI(j))=Qspec(bminI(j):bmaxI(j),j,1)

!***************** Matrix calculations *****************
#if RTYPE==2
		aLn(minI(j):maxI(j),1:NumSrc)=Atm_Contrib(minI(j):maxI(j),j,k,Form,1:NumSrc)
		if(bminI(j)==minI(j).and.maxI(j)==iGlob(j)) then
		  aLn(minI(j)-1,1:NumSrc)=Atm_Contrib(maxI(j),j,k,Form,1:NumSrc)
		else
		  aLn(minI(j)-1,1:NumSrc)=Atm_Contrib(bminI(j),j,k,Form,1:NumSrc)
		endif
		if(bmaxI(j)==maxI(j).and.maxI(j)==iGlob(j)) then
		  aLn(maxI(j)+1,1:NumSrc)=Atm_Contrib(1,j,k,Form,1:NumSrc)
		else
		  aLn(maxI(j)+1,1:NumSrc)=Atm_Contrib(bmaxI(j),j,k,Form,1:NumSrc)
		endif
#endif
!*******************************************************

!     Calculating Bott's polinomial coefficients
		call CoefsLong(j,SpecLong,CoefLong)

!     Calculation of fluxes through gridbox boundaries
	    call FluxLong(j,k,CoefLong,SpecLong,Jwest,Jeast)

	    do i=minI(j), maxI(j)
		  Qspec(i,j,2)=SpecLong(i)-Jwest(i)-Jeast(i)+Jeast(i-1)+Jwest(i+1)

!***************** Matrix calculations *****************
#if RTYPE==2
		  if(Qspec(i,j,2)>0.) then
			do Src=1, NumSrc
			  qCont=(SpecLong(i)-Jwest(i)-Jeast(i))*aLn(i,Src)+Jeast(i-1)*aLn(i-1,Src)+Jwest(i+1)*aLn(i+1,Src)
			  Contr(i,j,Src)=qCont/Qspec(i,j,2)
            enddo
		  else
			Contr(i,j,:)=0.
		  endif
#endif
!*******************************************************
	    enddo
		if(j>bJmin.and.j<bJmax.and.maxI(j)<iGlob(j)) then
		  MassAtmBnd(Form)=MassAtmBnd(Form)+(Jeast(minI(j)-1)-Jwest(minI(j)))*Veff(minI(j),j,k)
		  MassAtmBnd(Form)=MassAtmBnd(Form)+(Jwest(maxI(j)+1)-Jeast(maxI(j)))*Veff(maxI(j),j,k)
        endif
      enddo

!.................................................................................
!   Running in meridional direction
!.................................................................................

      do j=bJmin+1, bJmax-1
        SpecLat(minI(j):maxI(j),j)=Qspec(minI(j):maxI(j),j,1) 
        SpecLatNew(minI(j):maxI(j),j)=Qspec(minI(j):maxI(j),j,2) 
      enddo

      if(maxI(bJmin)==1) then
        SpecLat(minI(jSP+1):maxI(jSP+1),jSP)=Qspec(iSP,jSP,1)
        SpecLatNew(minI(jSP+1):maxI(jSP+1),jSP)=Qspec(iSP,jSP,1)
      else
        SpecLat(minI(bJmin):maxI(bJmin),bJmin)=Qspec(minI(bJmin):maxI(bJmin),bJmin,1)
        SpecLatNew(minI(bJmin):maxI(bJmin),bJmin)=Qspec(minI(bJmin):maxI(bJmin),bJmin,1)
      endif

      if(maxI(bJmax)==1) then
        SpecLat(minI(jNP-1):maxI(jNP-1),jNP)=Qspec(iNP,jNP,1)
        SpecLatNew(minI(jNP-1):maxI(jNP-1),jNP)=Qspec(iNP,jNP,1)
      else
        SpecLat(minI(bJmax):maxI(bJmax),bJmax)=Qspec(minI(bJmax):maxI(bJmax),bJmax,1)
        SpecLatNew(minI(bJmax):maxI(bJmax),bJmax)=Qspec(minI(bJmax):maxI(bJmax),bJmax,1)
      endif

!***************** Matrix calculations *****************
#if RTYPE==2
      do j=bJmin+1, bJmax-1  
        aLt(minI(j):maxI(j),j,1:NumSrc)=Contr(minI(j):maxI(j),j,1:NumSrc)
      enddo

! South pole
      if(maxI(bJmin)==1) then
        do i=minI(jSP+1), maxI(jSP+1)
          aLt(i,jSP,1:NumSrc)=Atm_Contrib(iSP,jSP,k,Form,1:NumSrc)
        enddo
      else
        aLt(minI(bJmin):maxI(bJmin),bJmin,1:NumSrc)=Atm_Contrib(minI(bJmin):maxI(bJmin),bJmin,k,Form,1:NumSrc)
      endif

! North pole
      if(maxI(bJmax)==1) then
        do i=minI(jNP-1), maxI(jNP-1)
          aLt(i,jNP,1:NumSrc)=Atm_Contrib(iNP,jNP,k,Form,1:NumSrc)
        enddo
      else
        aLt(minI(bJmax):maxI(bJmax),bJmax,1:NumSrc)=Atm_Contrib(minI(bJmax):maxI(bJmax),bJmax,k,Form,1:NumSrc)
      endif
#endif
!*******************************************************

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

! Northern facet
          if(n1==n2) then
            Qtemp=SpecLatNew(i,j)-dY/2.*(ctgdY2-tgY(j))*(Jnorth(i,j)-Jsouth(i,j+1))
          else
            Qtemp=SpecLatNew(i,j)-dY/2.*(ctgdY2-tgY(j))*(Jnorth(n1,j)-Jsouth(n1,j+1)&
                    &+Jnorth(n2,j)-Jsouth(n2,j+1))/2.
          endif
! Southern facet
          if(s1==s2) then
            Qtemp=Qtemp+dY/2.*(ctgdY2+tgY(j))*(Jnorth(i,j-1)-Jsouth(i,j))
          else
            Qtemp=Qtemp+dY/2.*(ctgdY2+tgY(j))*(Jnorth(s1,j-1)-Jsouth(s1,j)&
                    &+Jnorth(s2,j-1)-Jsouth(s2,j))/2.
          endif
          Atm_MixRatio(i,j,k,Form)=max(Qtemp/Pcalc(i,j,k),0.)

!***************** Matrix calculations *****************
#if RTYPE==2
          if(Qtemp>0.) then
            do Src=1, NumSrc
              if(n1==n2) then
!                qCont=SpecLatNew(i,j)*aLt(i,j,Src)-dY/2.*(ctgdY2-tgY(j))*(Jnorth(i,j)*aLt(i,j,Src)-Jsouth(i,j+1)*aLt(i,j+1,Src))  ! Check indecies of aLt !!!!!!!!!!!!!!!!!!!!!!!!!!
                qCont=SpecLatNew(i,j)*aLt(i,j,Src)-dY/2.*(ctgdY2-tgY(j))*(Jnorth(i,j)*aLt(i,j,Src)-Jsouth(i,j+1)*aLt(n1,j+1,Src))  ! Check indecies of aLt !!!!!!!!!!!!!!!!!!!!!!!!!!
              else
!                qCont=SpecLatNew(i,j)*aLt(i,j,Src)-dY/2.*(ctgdY2-tgY(j))*(Jnorth(n1,j)*aLt(n1,j,Src)&
!                    &-Jsouth(n1,j+1)*aLt(n1,j+1,Src)+Jnorth(n2,j)*aLt(n2,j,Src)-Jsouth(n2,j+1)*aLt(n2,j+1,Src))/2.
                qCont=SpecLatNew(i,j)*aLt(i,j,Src)-dY/2.*(ctgdY2-tgY(j))*(Jnorth(n1,j)*aLt(i,j,Src)&
                    &-Jsouth(n1,j+1)*aLt(n1,j+1,Src)+Jnorth(n2,j)*aLt(i,j,Src)-Jsouth(n2,j+1)*aLt(n2,j+1,Src))/2.
              endif
              if(s1==s2) then
!                qCont=qCont+dY/2.*(ctgdY2+tgY(j))*(Jnorth(i,j-1)*aLt(i,j-1,Src)-Jsouth(i,j)*aLt(i,j,Src))  ! Check indecies of aLt !!!!!!!!!!!!!!!!!!!!!!!!!!
                qCont=qCont+dY/2.*(ctgdY2+tgY(j))*(Jnorth(i,j-1)*aLt(s1,j-1,Src)-Jsouth(i,j)*aLt(i,j,Src))  ! Check indecies of aLt !!!!!!!!!!!!!!!!!!!!!!!!!!
              else
!                qCont=qCont+dY/2.*(ctgdY2+tgY(j))*(Jnorth(s1,j-1)*aLt(s1,j-1,Src)-Jsouth(s1,j)*aLt(s1,j,Src)&
!                    &+Jnorth(s2,j-1)*aLt(s2,j-1,Src)-Jsouth(s2,j)*aLt(s2,j,Src))/2.
                qCont=qCont+dY/2.*(ctgdY2+tgY(j))*(Jnorth(s1,j-1)*aLt(s1,j-1,Src)-Jsouth(s1,j)*aLt(i,j,Src)&
                    &+Jnorth(s2,j-1)*aLt(s2,j-1,Src)-Jsouth(s2,j)*aLt(i,j,Src))/2.
              endif
              qCont1(Src)=max(qCont/Qtemp,0.)
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
      enddo

      if(bJmin<Jmin) then
        do i=minI(Jmin), maxI(Jmin)
          MassAtmBnd(Form)=MassAtmBnd(Form)+dY/2.*(ctgdY2+tgY(Jmin))*(Jnorth(i,bJmin)-Jsouth(i,Jmin))*Veff(i,Jmin,k)
        enddo
      endif
      if(bJmax>Jmax) then
        do i=minI(Jmax), maxI(Jmax)
          MassAtmBnd(Form)=MassAtmBnd(Form)+dY/2.*(ctgdY2-tgY(Jmax))*(Jsouth(i,bJmax)-Jnorth(i,Jmax))*Veff(i,Jmax,k)
        enddo
      endif

!------------------ Global scale only -------------------
#if REGTYPE==1
      if(maxI(bJmin)==1) then
        PoleSpecSum=0.
	    do i=minI(jSP+1), maxI(jSP+1)
		  Qtemp=SpecLatNew(i,jSP)-dY/2.*(ctgdY2-tgY(jSP))*(Jnorth(i,jSP)-Jsouth(i,jSP+1))
		  PoleSpecSum=PoleSpecSum+Qtemp
	    enddo
	    Atm_MixRatio(iSP,jSP,k,Form)=max(PoleSpecSum/maxI(jSP+1)/Pcalc(iSP,jSP,k),0.)

!***************** Matrix calculations *****************
#if RTYPE==2
		if(PoleSpecSum>0.) then		  
		  do Src=1, NumSrc
			PoleContSum=0.
			do i=minI(jSP+1), maxI(jSP+1)
			  qCont=SpecLatNew(i,jSP)*aLt(i,jSP,Src)-dY/2.*(ctgdY2-tgY(jSP))*(Jnorth(i,jSP)&
                                    &*aLt(i,jSP,Src)-Jsouth(i,jSP+1)*aLt(i,jSP+1,Src))
			  PoleContSum=PoleContSum+qCont
			enddo
			Atm_Contrib(iSP,jSP,k,Form,Src)=PoleContSum/PoleSpecSum
		  enddo
		else
		  Atm_Contrib(iSP,jSP,k,Form,:)=0.
		endif
#endif
!*******************************************************
      endif
      if(maxI(bJmax)==1) then
	    PoleSpecSum=0.
	    do i=minI(jNP-1), maxI(jNP-1)
          Qtemp=SpecLatNew(i,jNP)+dY/2.*(ctgdY2+tgY(jNP))*(Jnorth(i,jNP-1)-Jsouth(i,jNP))
		  PoleSpecSum=PoleSpecSum+Qtemp
        enddo
	    Atm_MixRatio(iNP,jNP,k,Form)=max(PoleSpecSum/maxI(jNP-1)/Pcalc(iNP,jNP,k),0.)

!***************** Matrix calculations *****************
#if RTYPE==2
		if(PoleSpecSum>0.) then		  
		  do Src=1, NumSrc
			PoleContSum=0.
            do i=minI(jNP-1), maxI(jNP-1)
              qCont=SpecLatNew(i,jNP)*aLt(i,jNP,Src)+dY/2.*(ctgdY2+tgY(jNP))*(Jnorth(i,jNP-1)*aLt(i,jNP-1,Src)&
                                    &-Jsouth(i,jNP)*aLt(i,jNP,Src))
			  PoleContSum=PoleContSum+qCont
			enddo
			Atm_Contrib(iNP,jNP,k,Form,Src)=PoleContSum/PoleSpecSum
          enddo
		else
		  Atm_Contrib(iNP,jNP,k,Form,:)=0.
		endif
#endif
!*******************************************************
      endif
#endif
!------------------ Global scale only -------------------
	enddo

end subroutine Horiz_Advect


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CoefsLong(j,Spec,Coef)

    integer i, j
	real(8) Spec(bImin:bImax), SpL, Sp0, SpR
	real Coef(bImin:bImax,0:2)

	do i=bminI(j), bmaxI(j)
	  if(i==1.and.maxI(j)==iGlob(j)) then		! Periodic boundary conditions
	    SpL=Spec(maxI(j))
	  elseif(i<minI(j)) then				! Left regional boundary
		SpL=2.*Spec(i)-Spec(i+1)
	  else								! Normal case
	    SpL=Spec(i-1)
	  endif
	  Sp0=Spec(i)
	  if(i==iGlob(j).and.minI(j)==1) then		! Periodic boundary conditions
	    SpR=Spec(1)
	  elseif(i>maxI(j)) then				! Right regional boundary
		SpR=2.*Spec(i)-Spec(i-1)
	  else								! Normal case
	    SpR=Spec(i+1)
	  endif

	  Coef(i,0)=-(SpR-26.*Sp0+SpL)/24.
	  Coef(i,1)=(SpR-SpL)/2.
	  Coef(i,2)=(SpR-2.*Sp0+SpL)/2.
	enddo  
end subroutine CoefsLong


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CoefsLat(Spec,Coef)

    integer i, j, s1, s2, n1, n2, minIadv, maxIadv
    real(8) Spec(Imin:Imax,bJmin:bJmax), Sp0, SpS1, SpS2, SpN1, SpN2
    real Coef(Imin:Imax,bJmin:bJmax,0:2)

    do j=bJmin, bJmax
      minIadv=minI(j)
      maxIadv=maxI(j)
      if(maxI(j)==1.and.j==jNP) maxIadv=maxI(jNP-1)
      if(maxI(j)==1.and.j==jSP) maxIadv=maxI(jSP+1)

      do i=minIadv, maxIadv
        s1=iS(i,j,1)
        s2=iS(i,j,2)
        n1=iN(i,j,1)
        n2=iN(i,j,2)

        Sp0=Spec(i,j)
        if(j>bJmin) then
          SpS1=Spec(s1,j-1)
          SpS2=Spec(s2,j-1)
        else                                        ! Regional boundary
!          SpS1=2.*Sp0-(Spec(n1,j+1)+Spec(n2,j+1))/2.
!          SpS2=2.*Sp0-(Spec(n1,j+1)+Spec(n2,j+1))/2.
          SpS1=Sp0
          SpS2=Sp0
        endif
        if(j<bJmax) then
          SpN1=Spec(n1,j+1)
          SpN2=Spec(n2,j+1)
        else                                        ! Regional boundary
!          SpN1=2.*Sp0-(Spec(s1,j-1)+Spec(s2,j-1))/2.
!          SpN2=2.*Sp0-(Spec(s1,j-1)+Spec(s2,j-1))/2.
          SpN1=Sp0
          SpN2=Sp0
        endif

        if(s1/=s2.and.j>bJmin) then
          Coef(s1,j,0:2)=PolyLat(1,j,0:2)*SpS1+PolyLat(2,j,0:2)*Sp0+PolyLat(3,j,0:2)*SpN1
          Coef(s2,j,0:2)=PolyLat(1,j,0:2)*SpS2+PolyLat(2,j,0:2)*Sp0+PolyLat(3,j,0:2)*SpN1
        elseif(n1/=n2.and.j<bJmax) then
          Coef(n1,j,0:2)=PolyLat(1,j,0:2)*SpS1+PolyLat(2,j,0:2)*Sp0+PolyLat(3,j,0:2)*SpN1
          Coef(n2,j,0:2)=PolyLat(1,j,0:2)*SpS1+PolyLat(2,j,0:2)*Sp0+PolyLat(3,j,0:2)*SpN2
        else
          Coef(i,j,0:2)=PolyLat(1,j,0:2)*SpS1+PolyLat(2,j,0:2)*Sp0+PolyLat(3,j,0:2)*SpN1
        endif
      enddo
    enddo

end subroutine CoefsLat


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine FluxLong(j,k,Coef,Spec,FluxW,FluxE)	

	integer i, j, k, l 
	real Coef(bImin:bImax,0:2), cW, cE
	real(8) Spec(bImin:bImax) 
	real(8) FluxW(Imin:Imax+1), FluxE(Imin-1:Imax), fW, fE, Weight
	real Ktwo, Ksign, Kc

	do i=minI(j), maxI(j)
	  fW=0.
	  cW=CourW(i,j,k)
	  if(cW>Zero) then
	    Ktwo=2.
	    Ksign=1.
		Kc=1.-2.*cW
	    do l=0, 2
		  fW=fW+Coef(i,l)/(real(l)+1.)/Ktwo*Ksign*(1.-Kc)
	      Ktwo=Ktwo*2.
		  Ksign=-Ksign
		  Kc=Kc*(1.-2.*cW)
	    enddo
	    fW=max(0.,fW)			! Flux limitation 
	  else
	    fW=0.
	  endif
	  
	  fE=0.
	  cE=CourE(i,j,k)
	  if(cE>Zero) then
	    Ktwo=2.
	    Ksign=1.
	    Kc=1.-2.*cE
		do l=0, 2
		  fE=fE+Coef(i,l)/(real(l)+1.)/Ktwo*(1.-Kc)
	      Ktwo=Ktwo*2.
		  Ksign=-Ksign
		  Kc=Kc*(1.-2.*cE)
	    enddo
	    fE=max(0.,fE)			! Flux limitation
	  else
	    fE=0.
	  endif

	  Weight=min(1.,Spec(i)/(fW+fE+Zero))
	  FluxW(i)=fW*Weight
	  FluxE(i)=fE*Weight
	enddo

! Left boundary
    i=minI(j)-1
	if(bminI(j)==minI(j).and.maxI(j)==iGlob(j)) then
	  FluxE(i)=FluxE(maxI(j))
	else
	  fE=0.
	  cE=CourE(i,j,k)
	  if(cE>Zero) then
	    Ktwo=2.
	    Ksign=1.
	    Kc=1.-2.*cE
		do l=0, 2
		  fE=fE+Coef(i,l)/(real(l)+1.)/Ktwo*(1.-Kc)
	      Ktwo=Ktwo*2.
		  Ksign=-Ksign
		  Kc=Kc*(1.-2.*cE)
	    enddo
	    fE=max(0.,fE)			! Flux limitation
	  else
	    fE=0.
	  endif
	  FluxE(i)=min(fE,Spec(i))
	endif

! Right boundary
    i=maxI(j)+1
	if(bmaxI(j)==maxI(j).and.maxI(j)==iGlob(j)) then
	  FluxW(i)=FluxW(1)
	else
	  fW=0.
	  cW=CourW(i,j,k)
	  if(cW>Zero) then
	    Ktwo=2.
	    Ksign=1.
		Kc=1.-2.*cW
	    do l=0, 2
		  fW=fW+Coef(i,l)/(real(l)+1.)/Ktwo*Ksign*(1.-Kc)
	      Ktwo=Ktwo*2.
		  Ksign=-Ksign
		  Kc=Kc*(1.-2.*cW)
	    enddo
	    fW=max(0.,fW)			! Flux limitation 
	  else
	    fW=0.
	  endif
	  FluxW(i)=min(fW,Spec(i))
	endif

end subroutine FluxLong	


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine FluxLat(k,Coef,Spec,FluxN,FluxS)	

	integer i, j, k, l, s1, s2, n1, n2, minIadv, maxIadv
	real Coef(Imin:Imax,bJmin:bJmax,0:2), cN, cS
	real(8) Spec(Imin:Imax,bJmin:bJmax)
	real(8) FluxN(Imin:Imax,bJmin:Jmax), FluxS(Imin:Imax,Jmin:bJmax), fN1, fN2, fS1, fS2, Weight
	real Ktwo, Ksign, Kc

	do j=bJmin+1, bJmax-1
	  do i=minI(j), maxI(j)
		s1=iS(i,j,1)
		s2=iS(i,j,2)
		n1=iN(i,j,1)
		n2=iN(i,j,2)

		if(s1/=s2) then					! Grid seggregation southward 
		  fS1=0.
		  cS=CourS(s1,j,k)
		  if(cS>Zero) then
			Ktwo=2.
			Ksign=1.
			Kc=1.-2.*cS
			do l=0, 2
			  fS1=fS1+Coef(s1,j,l)/(real(l)+1.)/Ktwo*Ksign*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Ksign=-Ksign
			  Kc=Kc*(1.-2.*cS)
			enddo
			fS1=max(0.,fS1)	! Flux limitation 
		  endif

		  fS2=0.
		  cS=CourS(s2,j,k)
		  if(cS>Zero) then
			Ktwo=2.
			Ksign=1.
			Kc=1.-2.*cS
			do l=0, 2
			  fS2=fS2+Coef(s2,j,l)/(real(l)+1.)/Ktwo*Ksign*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Ksign=-Ksign
			  Kc=Kc*(1.-2.*cS)
			enddo
			fS2=max(0.,fS2)	! Flux limitation 
		  endif

		  fN1=0.
		  cN=CourN(i,j,k)
		  if(cN>Zero) then
			Ktwo=2.
			Kc=1.-2.*cN
			do l=0, 2
			  fN1=fN1+(Coef(s1,j,l)+Coef(s2,j,l))/2./(real(l)+1.)/Ktwo*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Kc=Kc*(1.-2.*cN)
			enddo
			fN1=max(0.,fN1)	! Flux limitation
		  endif
		  fN2=fN1

		elseif(n1/=n2) then				! Grid seggregation northward
		  fN1=0.
		  cN=CourN(n1,j,k)
		  if(cN>Zero) then
			Ktwo=2.
			Kc=1.-2.*cN
			do l=0, 2
			  fN1=fN1+Coef(n1,j,l)/(real(l)+1.)/Ktwo*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Kc=Kc*(1.-2.*cN)
			enddo
			fN1=max(0.,fN1)	! Flux limitation
		  endif

		  fN2=0.
		  cN=CourN(n2,j,k)
		  if(cN>Zero) then
			Ktwo=2.
			Kc=1.-2.*cN
			do l=0, 2
			  fN2=fN2+Coef(n2,j,l)/(real(l)+1.)/Ktwo*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Kc=Kc*(1.-2.*cN)
			enddo
			fN2=max(0.,fN2)	! Flux limitation
		  endif

		  fS1=0.
		  cS=CourS(i,j,k)
		  if(cS>Zero) then
			Ktwo=2.
			Ksign=1.
			Kc=1.-2.*cS
			do l=0, 2
			  fS1=fS1+(Coef(n1,j,l)+Coef(n2,j,l))/2./(real(l)+1.)/Ktwo*Ksign*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Ksign=-Ksign
			  Kc=Kc*(1.-2.*cS)
			enddo
			fS1=max(0.,fS1)	! Flux limitation 
		  endif
		  fS2=fS1

		else							! Normal case
		  fS1=0.
		  cS=CourS(i,j,k)
		  if(cS>Zero) then
			Ktwo=2.
			Ksign=1.
			Kc=1.-2.*cS
			do l=0, 2
			  fS1=fS1+Coef(i,j,l)/(real(l)+1.)/Ktwo*Ksign*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Ksign=-Ksign
			  Kc=Kc*(1.-2.*cS)
			enddo
			fS1=max(0.,fS1)	! Flux limitation 
		  endif
		  fS2=fS1

		  fN1=0.
		  cN=CourN(i,j,k)
		  if(cN>Zero) then
			Ktwo=2.
			Kc=1.-2.*cN
			do l=0, 2
			  fN1=fN1+Coef(i,j,l)/(real(l)+1.)/Ktwo*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Kc=Kc*(1.-2.*cN)
			enddo
			fN1=max(0.,fN1)	! Flux limitation
		  endif
		  fN2=fN1
		endif

		Weight=min(1.,Spec(i,j)/(dY/2.*((ctgdY2-tgY(j))*(fN1+fN2)/2.+(ctgdY2+tgY(j))*(fS1+fS2)/2.)+Zero))
		if(s1==s2) then
		  FluxS(i,j)=fS1*Weight
		else
		  FluxS(s1,j)=fS1*Weight
		  FluxS(s2,j)=fS2*Weight
		endif
		if(n1==n2) then
		  FluxN(i,j)=fN1*Weight
		else
		  FluxN(n1,j)=fN1*Weight
		  FluxN(n2,j)=fN2*Weight
		endif
	  enddo
	enddo

! Southern regional boundary
	j=bJmin
	minIadv=minI(j)
	maxIadv=maxI(j)
	if(maxI(j)==1.and.j==jSP) maxIadv=maxI(jSP+1)

	do i=minIadv, maxIadv
		n1=iN(i,j,1)
		n2=iN(i,j,2)

		if(n1/=n2) then					! Grid seggregation northward 
		  fN1=0.
		  cN=CourN(n1,j,k)
		  if(cN>Zero) then
			Ktwo=2.
			Kc=1.-2.*cN
			do l=0, 2
			  fN1=fN1+Coef(n1,j,l)/(real(l)+1.)/Ktwo*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Kc=Kc*(1.-2.*cN)
			enddo
			fN1=max(0.,fN1)	! Flux limitation
		  endif

		  fN2=0.
		  cN=CourN(n2,j,k)
		  if(cN>Zero) then
			Ktwo=2.
			Kc=1.-2.*cN
			do l=0, 2
			  fN2=fN2+Coef(n2,j,l)/(real(l)+1.)/Ktwo*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Kc=Kc*(1.-2.*cN)
			enddo
			fN2=max(0.,fN2)	! Flux limitation
		  endif

		else							! Normal case
		  fN1=0.
		  cN=CourN(i,j,k)
		  if(cN>Zero) then
			Ktwo=2.
			Kc=1.-2.*cN
			do l=0, 2
			  fN1=fN1+Coef(i,j,l)/(real(l)+1.)/Ktwo*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Kc=Kc*(1.-2.*cN)
			enddo
			fN1=max(0.,fN1)	! Flux limitation
		  endif
		  fN2=fN1
		endif

		Weight=min(1.,Spec(i,j)/(dY/2.*(ctgdY2-tgY(j))*(fN1+fN2)/2.+Zero))
		if(n1==n2) then
		  FluxN(i,j)=fN1*Weight
		else
		  FluxN(n1,j)=fN1*Weight
		  FluxN(n2,j)=fN2*Weight
		endif
	enddo

! Northern regional boundary
	j=bJmax
	minIadv=minI(j)
	maxIadv=maxI(j)
	if(maxI(j)==1.and.j==jNP) maxIadv=maxI(jNP-1)

	do i=minIadv, maxIadv
		s1=iS(i,j,1)
		s2=iS(i,j,2)

		if(s1/=s2) then					! Grid seggregation southward 
		  fS1=0.
		  cS=CourS(s1,j,k)
		  if(cS>Zero) then
			Ktwo=2.
			Ksign=1.
			Kc=1.-2.*cS
			do l=0, 2
		      fS1=fS1+Coef(s1,j,l)/(real(l)+1.)/Ktwo*Ksign*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Ksign=-Ksign
			  Kc=Kc*(1.-2.*cS)
			enddo
			fS1=max(0.,fS1)	! Flux limitation 
		  endif

		  fS2=0.
		  cS=CourS(s2,j,k)
		  if(cS>Zero) then
			Ktwo=2.
			Ksign=1.
			Kc=1.-2.*cS
			do l=0, 2
			  fS2=fS2+Coef(s2,j,l)/(real(l)+1.)/Ktwo*Ksign*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Ksign=-Ksign
			  Kc=Kc*(1.-2.*cS)
			enddo
			fS2=max(0.,fS2)	! Flux limitation 
		  endif

		else							! Normal case
		  fS1=0.
		  cS=CourS(i,j,k)
		  if(cS>Zero) then
			Ktwo=2.
			Ksign=1.
			Kc=1.-2.*cS
			do l=0, 2
		      fS1=fS1+Coef(i,j,l)/(real(l)+1.)/Ktwo*Ksign*(1.-Kc)
			  Ktwo=Ktwo*2.
			  Ksign=-Ksign
			  Kc=Kc*(1.-2.*cS)
			enddo
			fS1=max(0.,fS1)	! Flux limitation 
		  endif
		  fS2=fS1
		endif

		Weight=min(1.,Spec(i,j)/(dY/2.*(ctgdY2+tgY(j))*(fS1+fS2)/2.+Zero))
		if(s1==s2) then	
		  FluxS(i,j)=fS1*Weight
		else
		  FluxS(s1,j)=fS1*Weight
		  FluxS(s2,j)=fS2*Weight
		endif
	enddo

end subroutine FluxLat	


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating parameters of the Bott polinomials
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine BottPolynomLat

	integer j
	real betta, K, L, M

	do j=bJmin, bJmax
	  betta=dY*tgY(j)
	  M=6.-7.*betta*betta
	  K=betta/M
	  L=3.*(1.-betta*betta)/M

	  PolyLat(1,j,0)=-(1.-betta)/24.
	  PolyLat(2,j,0)=26./24.
	  PolyLat(3,j,0)=-(1.+betta)/24.
	  PolyLat(1,j,1)=-(1.-K)/2.
	  PolyLat(2,j,1)=-K
	  PolyLat(3,j,1)=(1.+K)/2.
	  PolyLat(1,j,2)=L
	  PolyLat(2,j,2)=-2.*L
	  PolyLat(3,j,2)=L
	enddo

end subroutine BottPolynomLat


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating Courant numbers in zonal direction
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CourNumsLong(dTime,CourM)

	integer i, j, k
	real dTime, CourM, Uw

    do k=1, Atm_Kmax
	  do j=bJmin, bJmax
	    if(maxI(j)==1) cycle

		do i=bminI(j), bmaxI(j)

		  if(i==1.and.maxI(j)==iGlob(j)) then					! Periodic boundary conditions
			Uw=Ucurr(maxI(j),j,k)
			if(Uw<0.) then
			  CourW(i,j,k)=-Uw*dTime/(Rearth*dX(j)*cosY(j))
			else
			  CourW(i,j,k)=0.
			endif
			Uw=Ucurr(i,j,k)
			if(Uw>0.) then  
			  CourE(i,j,k)=Uw*dTime/(Rearth*dX(j)*cosY(j))
			else
			  CourE(i,j,k)=0.
			endif
			CourM=max(CourM,CourE(i,j,k)+CourW(i,j,k))
		  elseif(i<minI(j)) then								! Left regional boundary
			Uw=Ucurr(i,j,k)
			if(Uw>0.) then  
			  CourE(i,j,k)=Uw*dTime/(Rearth*dX(j)*cosY(j))
			else
			  CourE(i,j,k)=0.
			endif
			CourM=max(CourM,CourE(i,j,k))
		  elseif(i>maxI(j)) then								! Right regional boundary
			Uw=Ucurr(i-1,j,k)
			if(Uw<0.) then
			  CourW(i,j,k)=-Uw*dTime/(Rearth*dX(j)*cosY(j))
			else
			  CourW(i,j,k)=0.
			endif
			CourM=max(CourM,CourW(i,j,k))
		  else													! Normal case
			Uw=Ucurr(i-1,j,k)
			if(Uw<0.) then
			  CourW(i,j,k)=-Uw*dTime/(Rearth*dX(j)*cosY(j))
			else
			  CourW(i,j,k)=0.
			endif
			Uw=Ucurr(i,j,k)
			if(Uw>0.) then  
			  CourE(i,j,k)=Uw*dTime/(Rearth*dX(j)*cosY(j))
			else
			  CourE(i,j,k)=0.
			endif
			CourM=max(CourM,CourE(i,j,k)+CourW(i,j,k))
		  endif
		enddo
	  enddo
	enddo

end subroutine CourNumsLong


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating Courant numbers in meridional direction
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CourNumsLat(dTime,CourM)

	integer i, j, k, s1, s2, n1, n2, minIadv, maxIadv
	real dTime, CourM, Vw

	do k=1, Atm_Kmax
	  do j=bJmin+1, bJmax-1
		do i=minI(j), maxI(j)
		  s1=iS(i,j,1)
		  s2=iS(i,j,2)
		  n1=iN(i,j,1)
		  n2=iN(i,j,2)
		  
		  if(s1==s2) then					! Normal case
			Vw=Vcurr(i,j-1,k)
			if(Vw<0.) then
			  CourS(i,j,k)=-Vw*dTime/(Rearth*dY)
			else
			  CourS(i,j,k)=0.
			endif
		  else								! Grid seggregation southward
		    Vw=Vcurr(s1,j-1,k)
			if(Vw<0.) then
			  CourS(s1,j,k)=-Vw*dTime/(Rearth*dY)
			else
			  CourS(s1,j,k)=0.
		    endif

			Vw=Vcurr(s2,j-1,k)
			if(Vw<0.) then
			  CourS(s2,j,k)=-Vw*dTime/(Rearth*dY)
			else
			  CourS(s2,j,k)=0.
			endif
		  endif

		  if(n2==n1) then					! Normal case
		    Vw=Vcurr(i,j,k)
			if(Vw>0.) then  
			  CourN(i,j,k)=Vw*dTime/(Rearth*dY)
			else
			  CourN(i,j,k)=0.
		    endif
		  else								! Grid seggregation  northward
		    Vw=Vcurr(n1,j,k)
			if(Vw>0.) then  
			  CourN(n1,j,k)=Vw*dTime/(Rearth*dY)
			else
			  CourN(n1,j,k)=0.
		    endif

			Vw=Vcurr(n2,j,k)
			if(Vw>0.) then  
			  CourN(n2,j,k)=Vw*dTime/(Rearth*dY)
			else
			  CourN(n2,j,k)=0.
			endif
		  endif
		  CourM=max(CourM,(CourN(n1,j,k)+CourN(n2,j,k)+CourS(s1,j,k)+CourS(s2,j,k))/2.)
		enddo
	  enddo

! South regional boundary
	  j=bJmin
	  minIadv=minI(j)
	  maxIadv=maxI(j)
	  if(maxI(j)==1.and.j==jSP) maxIadv=maxI(jSP+1)

	  do i=minIadv, maxIadv
		n1=iN(i,j,1)
		n2=iN(i,j,2)

		if(n2==n1) then						! Normal case
			Vw=Vcurr(i,j,k)
			if(Vw>0.) then  
			  CourN(i,j,k)=Vw*dTime/(Rearth*dY)
			else
			  CourN(i,j,k)=0.
			endif
		else								! Grid seggregation  northward
			Vw=Vcurr(n1,j,k)
			if(Vw>0.) then  
			  CourN(n1,j,k)=Vw*dTime/(Rearth*dY)
			else
			  CourN(n1,j,k)=0.
			endif

			Vw=Vcurr(n2,j,k)
			if(Vw>0.) then  
			  CourN(n2,j,k)=Vw*dTime/(Rearth*dY)
			else
			  CourN(n2,j,k)=0.
			endif
		endif		
		CourM=max(CourM,(CourN(n1,j,k)+CourN(n2,j,k))/2.)
	  enddo

! North regional boundary
	  j=bJmax											
	  minIadv=minI(j)
	  maxIadv=maxI(j)
	  if(maxI(j)==1.and.j==jNP) maxIadv=maxI(jNP-1)

	  do i=minIadv, maxIadv
		s1=iS(i,j,1)
		s2=iS(i,j,2)

		if(s2==s1) then						! Normal case
			Vw=Vcurr(i,j-1,k)
			if(Vw<0.) then
			  CourS(i,j,k)=-Vw*dTime/(Rearth*dY)
			else
			  CourS(i,j,k)=0.
			endif
		else								! Grid seggregation southward
			Vw=Vcurr(s1,j-1,k)
			if(Vw<0.) then
			  CourS(s1,j,k)=-Vw*dTime/(Rearth*dY)
			else
			  CourS(s1,j,k)=0.
			endif

			Vw=Vcurr(s2,j-1,k)
			if(Vw<0.) then
			  CourS(s2,j,k)=-Vw*dTime/(Rearth*dY)
			else
			  CourS(s2,j,k)=0.
			endif
		endif
		CourM=max(CourM,(CourS(s1,j,k)+CourS(s2,j,k))/2.)
	  enddo
	enddo

end subroutine CourNumsLat


end module Atm_HorizAdvect