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
! Module containing Hg specific procedures
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module Atm_Hg_Chemistry

#ifdef G_HG

  use GeneralParams
  use Geometry
  use Atm_Hg_Params
  use Exch_Params

  implicit none

  character(800), private :: fileName, fullName
  character(4),  private :: YearNum
  character(2),  private :: DayNum, YearShort
!  integer, parameter :: Nform=9
!  real(8) JacobChem(Nform,Nform)
  real(8) JacobChem(50,50)
  real(8) Krate(50)

contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating mercury chemistry in gaseous phase (old O3/OH scheme)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_GasChem_O3_OH    

    integer i, j, k, Src
    real(8) Cgas, Cpart, Coxid, c_O3, c_Cl2, c_OH, c_BrO, c_Br, c_HgBr
    real(8) dCgas, dCpart, dCoxid, dChgbr, Agas, Bgas, Ahgbr, Bhgbr, Aoxid, Boxid, Apart, Bpart
    real(8) T, meshV, dCbr, RhoAir, cExch, aG(MaxMatr), sumSrc
    real(8) Ko3, Kcl2, Koh, Kbr1, Kbr2, Kbr2a, Kbr2b, Kbr3, Kbro, Kox, Kgom, Kpbm, Acoef, Bcoef, Dter, L1, L2, C1, C2
    real cosSol, solRad, cosMean, Lw, Fw, LwMax, dT, tDay

    dT=Tstep(Atm)
    tDay=DayTime
    LwMax=0.
    do k=Atm_Kmax, 1, -1
      do j=Jmin, Jmax
        do i=minI(j), maxI(j)
          meshV=MeshVolume(i,j,k)
          Lw=LiqCont(i,j,k,Period)
          Fw=FrozCont(i,j,k,Period)
          LwMax=max(Lw+Fw,LwMax)
          RhoAir=DensAirNum(i,j,k)
        
          Cgas=Atm_Conc(i,j,k,Gas)
          Cpart=Atm_Conc(i,j,k,Part(1))
          Coxid=Atm_Conc(i,j,k,Oxid(1))
          c_O3=ConcO3(i,j,k)
          c_OH=ConcOH(i,j,k)
          c_BrO=0. !ConcBrO(i,j,k)                       ! no BrO chemistry
          c_Br=ConcBr(i,j,k)                             ! no Br chemistry
          c_HgBr=sum(ConcHgBr(i,j,k,1:NumSrc))           ! no Br chemistry

          solRad=SolarRad(i,j,Period)
          cosSol=SolarAngle(LongMesh(i,j),LatMesh(j),tDay,Day,Month,Year,cosMean)
          if(solRad>0..and.cosSol>0.) then
            c_Cl2=ConcCl2(i,j,k,Period)/10.    ! 10 ppt
          else
            c_Cl2=ConcCl2(i,j,k,Period)        ! 100 ppt
          endif
          if(LwMax>Lw0) c_OH=c_OH/10.            ! Below clouds

          T=TairCurr(i,j,k)
          Ko3=O3coef(1)*dexp(-O3coef(2)/Runiv/T)*c_O3
          Kcl2=Cl2coef(1)*c_Cl2
          Koh=OHcoef(1)*c_OH

          Kbr1=1.5e-32*(T/298.)**(-1.86)*c_Br*RhoAir          ! Donohoue et al. (2006)
          Kbr2a=4.e9*dexp(-7292/T)                            ! Goodsite et al. (2012)
!          Kbr2a=1.6e-9*dexp(-7801./T)*RhoAir                  ! Dibble et al. (2012)
          Kbr2b=3.9e-11*c_Br                                  ! Balabanov et al. (2005)
          Kbr2=Kbr2a+Kbr2b
          Kbr3=2.5e-10*(T/298.)**(-0.57)*(c_Br+c_OH)          ! Goodsite et al. (2004)
          Kbro=1.e-14*c_BrO                                   ! Sumner & Spicer (2005)

          Kgom=Kbro+Kcl2
          Kpbm=Ko3+Koh
          Kox=Kbr1+Kgom+Kpbm
          Acoef=Kox+Kbr2+Kbr3
          Bcoef=Kox*(Kbr2+Kbr3)-Kbr1*Kbr2
          Dter=Acoef*Acoef-4.*Bcoef
          if(Dter>Zero) then
            L1=-0.5*(Acoef+dsqrt(Dter))
            L2=-0.5*(Acoef-dsqrt(Dter))
            C1=(Kbr2*c_HgBr-Cgas*(L2+Kox))/(L1-L2)
            C2=-(Kbr2*c_HgBr-Cgas*(L1+Kox))/(L1-L2)
          else
            L1=0.
            L2=0.
            C1=0.
            C2=0.
          endif

          if(dabs(L1)>Zero) then
            Agas=C1*dexp(L1*dT)
            Ahgbr=C1/Kbr2*(L1+Kox)*dexp(L1*dT)
            Aoxid=C1/L1*(Kgom+Kbr3/Kbr2*(L1+Kox))*(dexp(L1*dT)-1.)
            Apart=C1/L1*(dexp(L1*dT)-1.)
          else
            Agas=C1
            Ahgbr=C1/Kbr2*Kox
            Aoxid=0.
            Apart=0.
          endif
          if(dabs(L2)>Zero) then
            Bgas=C2*dexp(L2*dT)
            Bhgbr=C2/Kbr2*(L2+Kox)*dexp(L2*dT)
            Boxid=C2/L2*(Kgom+Kbr3/Kbr2*(L2+Kox))*(dexp(L2*dT)-1.)
            Bpart=C2/L2*(dexp(L2*dT)-1.)
          else
            Bgas=C2
            Bhgbr=C2/Kbr2*Kox
            Boxid=0.
            Bpart=0.
          endif    
          dCgas=Agas+Bgas-Cgas
          dChgbr=Ahgbr+Bhgbr-c_HgBr
          dCoxid=Aoxid+Boxid
          dCpart=Kpbm*(Apart+Bpart)

          Atm_Conc(i,j,k,Gas)=max(Cgas+dCgas,0.)
          Atm_Conc(i,j,k,Part(1))=max(Cpart+dCpart,0.)
          Atm_Conc(i,j,k,Oxid(1))=max(Coxid+dCoxid,0.)

          MassChemEx(Gas)=MassChemEx(Gas)+dCgas*meshV
          MassChemEx(Part(1))=MassChemEx(Part(1))+dCpart*meshV
          MassChemEx(Oxid(1))=MassChemEx(Oxid(1))+dCoxid*meshV

!***************** Matrix calculations *****************
#if RTYPE==2
          cExch=0.
          aG=0.
          if(dCgas<0.) then
            cExch=cExch-dCgas
            do Src=1, NumSrc
              aG(Src)=aG(Src)-dCgas*Atm_Contrib(i,j,k,Gas,Src)
            enddo
          endif
          if(dChgbr<0.) then
            cExch=cExch-dChgbr
            sumSrc=sum(ConcHgBr(i,j,k,1:NumSrc))
            do Src=1, NumSrc
              if(sumSrc>0.) aG(Src)=aG(Src)-dChgbr*ConcHgBr(i,j,k,Src)/sumSrc
            enddo
          endif
          where(aG<0.) aG=0.
          if(cExch>0.) then
            aG=aG/cExch
          else
             aG=0.
          endif
          if(sum(aG)>0.) aG=aG/sum(aG)

          if(dCgas>0.) then 
            do Src=1, NumSrc
              Atm_Contrib(i,j,k,Gas,Src)=(Cgas*Atm_Contrib(i,j,k,Gas,Src)+dCgas*aG(Src))/(Cgas+dCgas)
            enddo
          endif
          if(dCpart>0.) then 
            do Src=1, NumSrc
              Atm_Contrib(i,j,k,Part(1),Src)=(Cpart*Atm_Contrib(i,j,k,Part(1),Src)+dCpart*aG(Src))/(Cpart+dCpart)
            enddo
          endif
          if(dCoxid>0.) then 
            do Src=1, NumSrc
              Atm_Contrib(i,j,k,Oxid(1),Src)=(Coxid*Atm_Contrib(i,j,k,Oxid(1),Src)+dCoxid*aG(Src))/(Coxid+dCoxid)
            enddo
          endif
          if(dChgbr>0.) then 
            do Src=1, NumSrc
              ConcHgBr(i,j,k,Src)=ConcHgBr(i,j,k,Src)+dChgbr*aG(Src)
            enddo
          else  
            sumSrc=sum(ConcHgBr(i,j,k,1:NumSrc))
            do Src=1, NumSrc
              if(sumSrc>0.) ConcHgBr(i,j,k,Src)=max(ConcHgBr(i,j,k,Src)*(1.+dChgbr/sumSrc),0.)
            enddo
          endif
!*******************************************************
#else
          ConcHgBr(i,j,k,1)=max(c_HgBr+dChgbr,0.)
#endif

        enddo
      enddo
    enddo

end subroutine Atm_Hg_GasChem_O3_OH   


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating mercury chemistry in aqueous phase
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_DropChem    

    integer i, j, k, Ind, Form, Src, Naq, iOxid, fOxid, iPart, fPart
    integer FormAq(MaxForm)
    real(8) Ao, Bo, Co, Keff
    real(8) T, HenryHg, HenryO3, HenryHgCl2, HenryCl2, HenryOH !, HenryHO2
    real(8) c_O3, c_SO2, c_Cl, c_Cl2, c_OH, Lw, Fw !, c_HO2
    real(8) Rac, Rba, Rcb, Rca, alpha, betta, gamma, delta, epsil, r1, r2, r3
    real(8) L2, L3, p, q, X, Y, C1, C2, C3, Xa, Xb, Xc
    real(8) Conc(MaxForm), dC(MaxForm)
    real(8) meshV, cExch, aAq(MaxMatr)
    real cosSol, cosMean, dT, tDay

    Naq=4+Noxid+Npart
    iOxid=5
    fOxid=4+Noxid
    iPart=5+Noxid
    fPart=4+Noxid+Npart
    
    FormAq(1:4)=(/Hg0,Dis,Sulf,Chlor/)
    FormAq(iOxid:fOxid)=Oxid(1:Noxid)
    FormAq(iPart:fPart)=Part(1:Npart)

    dT=Tstep(Atm)
    tDay=DayTime
    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        do i=minI(j), maxI(j)
          meshV=MeshVolume(i,j,k)

          Conc(FormAq(1:Naq))=Atm_Conc(i,j,k,FormAq(1:Naq))

! Checking for cloud evaporation
          Lw=LiqCont(i,j,k,Period)
          Fw=FrozCont(i,j,k,Period)
          if(Lw+Fw<=Lw0) then
            Atm_Conc(i,j,k,Gas)=Conc(Gas)+Conc(Dis)
            Atm_Conc(i,j,k,Part(1:Npart))=Conc(Part(1:Npart))+(Conc(Sulf)+Conc(Chlor))/real(Npart,8)
            Atm_Conc(i,j,k,Dis)=0.
            Atm_Conc(i,j,k,Sulf)=0.
            Atm_Conc(i,j,k,Chlor)=0.

            dC(Gas)=Conc(Dis)
            dC(Dis)=-Conc(Dis)
            dC(Sulf)=-Conc(Sulf)
            dC(Chlor)=-Conc(Chlor)
            dC(Part(1:Npart))=(Conc(Sulf)+Conc(Chlor))/real(Npart,8)

!***************** Matrix calculations *****************
#if RTYPE==2
            if(dC(Gas)>0.) then
                Atm_Contrib(i,j,k,Gas,1:NumSrc)=(Conc(Gas)*Atm_Contrib(i,j,k,Gas,1:NumSrc)+&
                    &dC(Gas)*Atm_Contrib(i,j,k,Dis,1:NumSrc))/(Conc(Gas)+dC(Gas))
            endif
            if(sum(dC(Part(1:Npart)))>0.) then
              aAq(1:NumSrc)=(dC(Sulf)*Atm_Contrib(i,j,k,Sulf,1:NumSrc)+&
                    &dC(Chlor)*Atm_Contrib(i,j,k,Chlor,1:NumSrc))/(dC(Sulf)+dC(Chlor))
              do Src=1, NumSrc
                  Atm_Contrib(i,j,k,Part(1:Npart),Src)=(Conc(Part(1:Npart))*Atm_Contrib(i,j,k,Part(1:Npart),Src)+&
                    &dC(Part(1:Npart))*aAq(Src))/(Conc(Part(1:Npart))+dC(Part(1:Npart)))
              enddo
            endif
#endif
!*******************************************************

            do Ind=1, 4
              Form=FormAq(Ind)  
              MassChemEx(Form)=MassChemEx(Form)+dC(Form)*meshV
            enddo
            MassChemEx(Part(1:Npart))=MassChemEx(Part(1:Npart))+dC(Part(1:Npart))*meshV
            cycle
          endif

      if(Lw<=Lw0) cycle
          Keff=Lw/(Lw+Aeff)

! Initial conditions
          Ao=Conc(Gas)+Conc(Dis)*Lw/(Lw+Fw)
          Bo=Conc(Sulf)*Lw/(Lw+Fw)
          Co=Conc(Chlor)*Lw/(Lw+Fw)+sum(Conc(Oxid(1:Noxid)))+sum(Conc(Part(1:Npart)))*Asol*Keff

! Chemichal parameters
          cosSol=SolarAngle(LongMesh(i,j),LatMesh(j),tDay,Day,Month,Year,cosMean)
          if(cosSol>0.) then
            c_Cl2=ConcCl2(i,j,k,Period)/10.    ! 10 ppt
          else
            c_Cl2=ConcCl2(i,j,k,Period)        ! 100 ppt
          endif
          T=TairCurr(i,j,k)
          c_O3=ConcO3(i,j,k)
          c_OH=ConcOH(i,j,k)/10.           ! in clouds
          c_SO2=ConcSO2(i,j,k)
          c_Cl=ConcCl
          HenryHg=CHenryHg(1)*T*dexp(CHenryHg(2)*(1./T-1./298.))
          HenryO3=CHenryO3(1)*T*dexp(CHenryO3(2)*(1./T-1./298.))
          HenryOH=CHenryOH(1)*T*dexp(CHenryOH(2)*(1./T-1./298.))
          HenryHgCl2=CHenryHgCl2(1)*T*dexp(CHenryHgCl2(2)*(1./T-1./298.))
          HenryCl2=CHenryCl2(1)

          r1=c_Cl/Clcoef(1)+c_Cl*c_Cl/Clcoef(2)+c_Cl*c_Cl*c_Cl/Clcoef(3)&
                &+c_Cl*c_Cl*c_Cl*c_Cl/Clcoef(4)
          r2=Rsoot
          r3=HenryHgCl2*Lw/(RhoWat+HenryHgCl2*Lw)

          alpha=HenryHg*Lw/(RhoWat+HenryHg*Lw)         ! Fraction of [Hg0]dis in [A]
          betta=r2/(1.+r2)                             ! Fraction of [Hg(So3)2-]aq in [B] 
          gamma=r2*r3/(r1*r3+r2*r3+r1*r2)              ! Fraction of [Hg2+]aq in [C]
          delta=gamma*(1.+r1)                          ! Fraction of [Hg2+]aq and [HgnClm]dis in [C]
          epsil=r1*r2*(1.-r3)/(r1*r3+r2*r3+r1*r2)      ! Fraction of HgCL2(gas) in [C]

          Rac=alpha*(Kac1*HenryO3*c_O3+Kac2*HenryCl2*c_Cl2+Kac3*HenryOH*c_OH)
          Rba=betta*Kba
!          Rba=betta*T*dexp((31.971*T-12595)/T)                             ! Van Loon et al (2000)
          Rcb=gamma*Kcb*c_SO2*c_SO2*10.**(4.*pH)
          Rca=0.
!          if(dabs(Rca-Rba)/(Rca+Rba)<Eps) Rca=Rca*0.99      ! Singular case

          if(Rca+Rcb>dsqrt(Zero)) then
            p=Rac+Rba+Rcb+Rca
            q=Rac*Rba+Rcb*Rac+Rca*Rba+Rcb*Rba
            if(p*p-4*q<=0.) cycle

            L2=-0.5*(p+dsqrt(p*p-4*q))
            L3=-0.5*(p-dsqrt(p*p-4*q))
            X=(Rba-Rca)*Rcb*Co+Rca*Rac*Ao-Rba*Rba*Bo-Rca*Rca*Co
            Y=Rac*Ao-Rba*Bo-Rca*Co

            if(L2==0..or.L3==0.) cycle            ! Singular case

            C1=Ao+(X+(L2+L3+Rac)*Y)/L2/L3
            C2=(X+(L3+Rac)*Y)/L2/(L2-L3)
            C3=(X+(L2+Rac)*Y)/L3/(L3-L2)

            Xa=C1+C2*dexp(L2*dT)+C3*dexp(L3*dT)
            Xb=(Rac*Rcb*C1+(Rac*Rcb+(Rcb-Rca*(Rca+Rac)/(Rba-Rca))*L2-&
                &Rca*L2*L2/(Rba-Rca))*C2*dexp(L2*dT)+&
                &(Rac*Rcb+(Rcb-Rca*(Rca+Rac)/(Rba-Rca))*L3-&
                &Rca*L3*L3/(Rba-Rca))*C3*dexp(L3*dT))/Rba/(Rca+Rcb)
            Xc=(Rac*(Rba-Rca)*C1+(Rac*(Rba-Rca)+(Rac+Rba)*L2+L2*L2)*C2*dexp(L2*dT)+&
                &(Rac*(Rba-Rca)+(Rac+Rba)*L3+L3*L3)*C3*dexp(L3*dT))/(Rca+Rcb)/(Rba-Rca)
          else
            if(Rac/=Rba) then
              C1=Ao+Bo+Co
              C2=Ao-Rba*Bo/(Rac-Rba)
              C3=Rba*Bo/(Rac-Rba)
            else
              C1=Ao+Co
              C2=Ao
              C3=0.
              cycle
            endif

            Xa=C2*dexp(-Rac*dT)+C3*dexp(-Rba*dT)
            Xb=(Rac-Rba)/Rba*C3*dexp(-Rba*dT)
            Xc=C1-C2*dexp(-Rac*dT)-Rac/Rba*C3*dexp(-Rba*dT)
          endif

          Atm_Conc(i,j,k,Dis)=alpha*Xa+Conc(Dis)*Fw/(Lw+Fw)
          Atm_Conc(i,j,k,Gas)=(1.-alpha)*Xa
          Atm_Conc(i,j,k,Sulf)=Xb+Conc(Sulf)*Fw/(Lw+Fw)
          Atm_Conc(i,j,k,Chlor)=(1.-epsil)*Xc+Conc(Chlor)*Fw/(Lw+Fw)
          Atm_Conc(i,j,k,Oxid(1:Noxid))=epsil*Xc/real(Noxid,8)
          Atm_Conc(i,j,k,Part(1:Npart))=Conc(Part(1:Npart))*(1.-Asol*Keff)

          dC(Dis)=alpha*Xa-Conc(Dis)*Lw/(Lw+Fw)    
          dC(Gas)=(1.-alpha)*Xa-Conc(Gas)
          dC(Sulf)=Xb-Conc(Sulf)*Lw/(Lw+Fw)
          dC(Chlor)=(1.-epsil)*Xc-Conc(Chlor)*Lw/(Lw+Fw)
          dC(Oxid(1:Noxid))=epsil*Xc/real(Noxid,8)-Conc(Oxid(1:Noxid))
          dC(Part(1:Npart))=-Conc(Part(1:Npart))*Asol*Keff

          MassChemEx(FormAq(1:Naq))=MassChemEx(FormAq(1:Naq))+dC(FormAq(1:Naq))*meshV
          
! !***************** Matrix calculations *****************
#if RTYPE==2
          cExch=0.
          aAq=0.
          do Ind=1, Naq
            Form=FormAq(Ind)  
            if(dC(Form)<0.) then
              cExch=cExch-dC(Form)
              do Src=1, NumSrc
                aAq(Src)=aAq(Src)-dC(Form)*Atm_Contrib(i,j,k,Form,Src)
              enddo
            endif
          enddo
          if(cExch>0.) then
            aAq=aAq/cExch
          else
             aAq=0.
          endif

          do Ind=1, Naq
            Form=FormAq(Ind) 
            if(dC(Form)>0.) then 
              do Src=1, NumSrc
                Atm_Contrib(i,j,k,Form,Src)=(Conc(Form)*Atm_Contrib(i,j,k,Form,Src)+&
                    &dC(Form)*aAq(Src))/(Conc(Form)+dC(Form))
              enddo
            endif
          enddo
#endif
!*******************************************************
        enddo
      enddo
    enddo

end subroutine Atm_Hg_DropChem


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating mercury chemistry in gaseous phase (new Br cheme)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_GasChem_Br    

    integer i, j, k, Src, Form, iFlag, iForm(50), Nform
    real(8) Cform0(50), Cform(50), dCform(50)
    real(8) c_O3, c_Cl2, c_OH, c_BrO, c_Br, c_HgBr, c_NO2, c_HO2
    real(8) Agas, Bgas, Ahgbr, Bhgbr, Aoxid, Boxid, Apart, Bpart
    real(8) T, meshV, dCbr, RhoAirNum, cExch, aG(MaxMatr), sumSrc
    real(8) Rphoto(50)
    real(8) Lw, Fw, LwMax, dT, tDay

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_Hg_GasChem...'
#endif

    Nform=Noxid+1
    iForm(1:9)=(/Hg0,HgBr,HgBr2,HgBrOH,HgBrOOH,HgBrONO,HOHg,HOHgOOH,HOHgONO/)    
    
    dT=Tstep(Atm)
    tDay=DayTime
    LwMax=0.
    do k=Atm_Kmax, 1, -1
      do j=Jmin, Jmax
        do i=minI(j), maxI(j)
          meshV=MeshVolume(i,j,k)
          Lw=LiqCont(i,j,k,Period)
          Fw=FrozCont(i,j,k,Period)
          LwMax=max(Lw+Fw,LwMax)
          RhoAirNum=DensAirNum(i,j,k)
        
          Cform0(1:Nform)=Atm_Conc(i,j,k,iForm(1:Nform))
 
          c_Br=ConcBr(i,j,k)
          c_OH=ConcOH(i,j,k)
          c_NO2=ConcNO2(i,j,k)
          c_HO2=ConcHO2(i,j,k)
!          c_O3=ConcO3(i,j,k)
          T=TairCurr(i,j,k)
          
          Rphoto(1:Noxid)=PhotoRate(1:Noxid,i,j,k)
          
          Krate(1)=1.5e-32*(T/298.)**(-1.86)*c_Br*RhoAirNum         ! Hg0+Br+M->HgBr+M    Donohoue et al. (2006)
          Krate(2)=1.6e-9*dexp(-7801./T)*RhoAirNum                  ! HgBr+M->Hg0+Br+M    Dibble et al. (2012)
          Krate(3)=3.9e-11*c_Br                                     ! HgBr+Br->Hg0+Br2    Balabanov et al. (2005)
          Krate(4)=Rphoto(HgBr)                                     ! HgBr->Hg0+prod      Saiz-Lopez et al. (2019)
          Krate(5)=2.5e-10*(T/298.)**(-0.57)*c_Br                   ! HgBr+Br->HgBr2      Goodsite et al. (2004)
          Krate(6)=2.5e-10*(T/298.)**(-0.57)*c_OH                   ! HgBr+OH->HgBrOH     Goodsite et al. (2004)
          Krate(7)=rate_HO2(T,RhoAirNum)*c_HO2                      ! HgBr+HO2->HgBrOOH   Jiao and Dibble (2017)  
          Krate(8)=rate_NO2(T,RhoAirNum)*c_NO2                      ! HgBr+NO2->HgBrONO   Jiao and Dibble (2017)
          Krate(9)=Rphoto(HgBr2)                                    ! HgBr2->prod         Saiz-Lopez et al. (2018)
          Krate(10)=Rphoto(HgBrOH)                                  ! HgBrOH->prod        Saiz-Lopez et al. (2018)
          Krate(11)=Rphoto(HgBrOOH)                                 ! HgBrOOH>prod        Saiz-Lopez et al. (2018)
          Krate(12)=Rphoto(HgBrONO)                                 ! HgBrONO->prod       Saiz-Lopez et al. (2018)
          Krate(13)=3.34e-33*dexp(43./T)*c_OH*RhoAirNum             ! Hg0+OH+M->HOHg+M    Dibble et al. (2019)
          Krate(14)=1.22e-9*dexp(-5720./T)*RhoAirNum                ! HOHg+M->Hg0+OH+M    Dibble et al. (2019)
          Krate(15)=Rphoto(HOHg)                                    ! HOHg->Hg0+prod      Saiz-Lopez et al. (2019)
          Krate(16)=rate_HO2(T,RhoAirNum)*c_HO2                     ! HOHg+HO2->HOHgOOH   Jiao and Dibble (2017) 
          Krate(17)=rate_NO2(T,RhoAirNum)*c_NO2                     ! HOHg+NO2->HOHgONO   Jiao and Dibble (2017)
          Krate(18)=Rphoto(HgBrOOH)                                 ! HOHgOOH->HOHg+prod  Saiz-Lopez et al. (2018)
          Krate(19)=Rphoto(HgBrONO)                                 ! HOHgONO->HOHg+prod  Saiz-Lopez et al. (2018)
          
          Cform(1:Nform)=Cform0(1:Nform)

          call RoDas_ODE(Nform,0.D0,Cform,dT,Chem_Deriv,Jacob_Br_Chem,iFlag)
          if(iFlag<0) then
            print *, 'STOP: Computation of chemical ODE was not successful, IFlag=', iFlag
            stop
          endif  
!         iFlag= 1  Computation successful,
!         iFlag=-2  Larger Nnax is needed,
!         iFlag=-3  Step size becomes too small,
!         iFlag=-4  Matrix is repeatedly singular.
          
          where(Cform<0.) Cform=0.D0
          dCform(1:Nform)=Cform(1:Nform)-Cform0(1:Nform)
          if(dabs(sum(dCform(1:Nform)))>Zero) then
            Cform(1)=Cform(1)-sum(dCform(1:Nform))
            dCform(1)=Cform(1)-Cform0(1)
          endif
          
          Atm_Conc(i,j,k,iForm(1:Nform))=max(Cform(1:Nform),0.D0)
          MassChemEx(iForm(1:Nform))=MassChemEx(iForm(1:Nform))+dCform(1:Nform)*meshV
         
!***************** Matrix calculations *****************
#if RTYPE==2
          cExch=0.
          aG=0.
          do Form=1, Nform
            if(dCform(Form)<0.) then
              cExch=cExch-dCform(Form)
              do Src=1, NumSrc
                aG(Src)=aG(Src)-dCform(Form)*Atm_Contrib(i,j,k,iForm(Form),Src)
              enddo
            endif
          enddo
          where(aG<0.) aG=0.
          if(cExch>0.) then
            aG=aG/cExch
          else
            aG=0.
          endif
          if(sum(aG)>0.) aG=aG/sum(aG)
          
          do Form=1, Nform
            if(dCform(Form)>0.) then 
              do Src=1, NumSrc
                Atm_Contrib(i,j,k,iForm(Form),Src)=(Cform(Form)*Atm_Contrib(i,j,k,iForm(Form),Src)&
                          &+dCform(Form)*aG(Src))/(Cform(Form)+dCform(Form))
              enddo
            endif
          enddo
#endif
!*******************************************************
        enddo
      enddo
    enddo

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_Hg_GasChem'
#endif

  contains

  function rate_NO2(T,rhoNum)

    integer n
    real(8) T, rhoNum, rate_NO2    
    real(8) Temp(5), k0(5), k00(5), k1, k11, p
    data Temp /220., 260., 280., 298., 320./
    data k0 /27.4e-29, 13.5e-29, 9.52e-29, 7.1e-29, 5.09e-29/
    data k00 /22.e-11, 14.2e-11, 12.8e-11, 11.8e-11, 10.9e-11/
  
    do n=1, 5
      if(Temp(n)>T) exit
    enddo

    if(n==1) then
      k1=k0(1)
      k11=k00(1)
    elseif(n==6) then
      k1=k0(5)
      k11=k00(5)
    else
      k1=k0(n-1)+(k0(n)-k0(n-1))/(Temp(n)-Temp(n-1))*(Temp(n)-T)
      k11=k00(n-1)+(k00(n)-k00(n-1))/(Temp(n)-Temp(n-1))*(Temp(n)-T)
    endif

    p=1./(1.+dlog10(k1*rhoNum/k11)*dlog10(k1*rhoNum/k11))
    rate_NO2=k1*rhoNum/(1.+k1*rhoNum/k11)*0.6**p

  end function rate_NO2

  function rate_HO2(T,rhoNum)

    integer n
    real(8) T, rhoNum, rate_HO2    
    real(8) Temp(5), k0(5), k00(5), k1, k11, p
    data Temp /220., 260., 280., 298., 320./
    data k0 /8.4e-29, 4.28e-29, 3.01e-29, 2.27e-29, 1.64e-29/
    data k00 /14.8e-11, 9.1e-11, 7.55e-11, 6.99e-11, 6.11e-11/
  
    do n=1, 5
      if(Temp(n)>T) exit
    enddo

    if(n==1) then
      k1=k0(1)
      k11=k00(1)
    elseif(n==6) then
      k1=k0(5)
      k11=k00(5)
    else
      k1=k0(n-1)+(k0(n)-k0(n-1))/(Temp(n)-Temp(n-1))*(Temp(n)-T)
      k11=k00(n-1)+(k00(n)-k00(n-1))/(Temp(n)-Temp(n-1))*(Temp(n)-T)
    endif

    p=1./(1.+dlog10(k1*rhoNum/k11)*dlog10(k1*rhoNum/k11))
    rate_HO2=k1*rhoNum/(1.+k1*rhoNum/k11)*0.6**p

  end function rate_HO2

end subroutine Atm_Hg_GasChem_Br    


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating mercury chemistry in gaseous phase (old Br cheme)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_GasChem_Br2    

    integer i, j, k, Src, Form, iFlag, iForm(50), Nform
    real(8) Cform0(50), Cform(50), dCform(50)
    real(8) c_O3, c_Cl2, c_OH, c_BrO, c_Br, c_HgBr, c_NO2, c_HO2
    real(8) Agas, Bgas, Ahgbr, Bhgbr, Aoxid, Boxid, Apart, Bpart
    real(8) T, meshV, dCbr, RhoAirNum, cExch, aG(MaxMatr), sumSrc
    real(8) Rphoto(50)
    real(8) Lw, Fw, LwMax, dT, tDay

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_Hg_GasChem...'
#endif

    Nform=Noxid+1
    iForm(1:9)=(/Hg0,HgBr,HgBr2,HgBrOH,HgBrOOH,HgBrONO,HOHg,HOHgOOH,HOHgONO/)    
    
    dT=Tstep(Atm)
    tDay=DayTime
    LwMax=0.
    do k=Atm_Kmax, 1, -1
      do j=Jmin, Jmax
        do i=minI(j), maxI(j)
          meshV=MeshVolume(i,j,k)
          Lw=LiqCont(i,j,k,Period)
          Fw=FrozCont(i,j,k,Period)
          LwMax=max(Lw+Fw,LwMax)
          RhoAirNum=DensAirNum(i,j,k)
        
          Cform0(1:Nform)=Atm_Conc(i,j,k,iForm(1:Nform))
 
          c_Br=ConcBr(i,j,k)
          c_OH=ConcOH(i,j,k)
          T=TairCurr(i,j,k)
          
          Rphoto(1:Noxid)=PhotoRate(1:Noxid,i,j,k)
          
          Krate(1)=1.5e-32*(T/298.)**(-1.86)*c_Br*RhoAirNum         ! Hg0+Br+M->HgBr+M    Donohoue et al. (2006)
          Krate(2)=1.6e-9*dexp(-7801./T)*RhoAirNum                  ! HgBr+M->Hg0+Br+M    Dibble et al. (2012)
          Krate(3)=3.9e-11*c_Br                                     ! HgBr+Br->Hg0+Br2    Balabanov et al. (2005)
          Krate(4)=2.5e-10*(T/298.)**(-0.57)*(c_Br+c_OH)            ! HgBr+Br->HgBr2      Goodsite et al. (2004)
                                                                    ! HgBr+OH->HgBrOH     Goodsite et al. (2004)
          
          Cform(1:Nform)=Cform0(1:Nform)

          call RoDas_ODE(Nform,0.D0,Cform,dT,Chem_Deriv,Jacob_Br_Chem2,iFlag)
          if(iFlag<0) then
            print *, 'STOP: Computation of chemical ODE was not successful, IFlag=', iFlag
            stop
          endif  
!         iFlag= 1  Ccomputation successful,
!         iFlag=-2  Larger Nnax is needed,
!         iFlag=-3  Step size becomes too small,
!         iFlag=-4  Matrix is repeatedly singular.
          
          where(Cform<0.) Cform=0.D0
          dCform(1:Nform)=Cform(1:Nform)-Cform0(1:Nform)
          if(dabs(sum(dCform(1:Nform)))>Zero) then
            Cform(1)=Cform(1)-sum(dCform(1:Nform))
            dCform(1)=Cform(1)-Cform0(1)
          endif
          
          Atm_Conc(i,j,k,iForm(1:Nform))=max(Cform(1:Nform),0.D0)
          MassChemEx(iForm(1:Nform))=MassChemEx(iForm(1:Nform))+dCform(1:Nform)*meshV
         
!***************** Matrix calculations *****************
#if RTYPE==2
          cExch=0.
          aG=0.
          do Form=1, Nform
            if(dCform(Form)<0.) then
              cExch=cExch-dCform(Form)
              do Src=1, NumSrc
                aG(Src)=aG(Src)-dCform(Form)*Atm_Contrib(i,j,k,iForm(Form),Src)
              enddo
            endif
          enddo
          where(aG<0.) aG=0.
          if(cExch>0.) then
            aG=aG/cExch
          else
            aG=0.
          endif
          if(sum(aG)>0.) aG=aG/sum(aG)
          
          do Form=1, Nform
            if(dCform(Form)>0.) then 
              do Src=1, NumSrc
                Atm_Contrib(i,j,k,iForm(Form),Src)=(Cform(Form)*Atm_Contrib(i,j,k,iForm(Form),Src)&
                          &+dCform(Form)*aG(Src))/(Cform(Form)+dCform(Form))
              enddo
            endif
          enddo
#endif
!*******************************************************
        enddo
      enddo
    enddo

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_Hg_GasChem'
#endif

  contains

  function rate_NO2(T,rhoNum)

    integer n
    real(8) T, rhoNum, rate_NO2    
    real(8) Temp(5), k0(5), k00(5), k1, k11, p
    data Temp /220., 260., 280., 298., 320./
    data k0 /27.4e-29, 13.5e-29, 9.52e-29, 7.1e-29, 5.09e-29/
    data k00 /22.e-11, 14.2e-11, 12.8e-11, 11.8e-11, 10.9e-11/
  
    do n=1, 5
      if(Temp(n)>T) exit
    enddo

    if(n==1) then
      k1=k0(1)
      k11=k00(1)
    elseif(n==6) then
      k1=k0(5)
      k11=k00(5)
    else
      k1=k0(n-1)+(k0(n)-k0(n-1))/(Temp(n)-Temp(n-1))*(Temp(n)-T)
      k11=k00(n-1)+(k00(n)-k00(n-1))/(Temp(n)-Temp(n-1))*(Temp(n)-T)
    endif

    p=1./(1.+dlog10(k1*rhoNum/k11)*dlog10(k1*rhoNum/k11))
    rate_NO2=k1*rhoNum/(1.+k1*rhoNum/k11)*0.6**p

  end function rate_NO2

  function rate_HO2(T,rhoNum)

    integer n
    real(8) T, rhoNum, rate_HO2    
    real(8) Temp(5), k0(5), k00(5), k1, k11, p
    data Temp /220., 260., 280., 298., 320./
    data k0 /8.4e-29, 4.28e-29, 3.01e-29, 2.27e-29, 1.64e-29/
    data k00 /14.8e-11, 9.1e-11, 7.55e-11, 6.99e-11, 6.11e-11/
  
    do n=1, 5
      if(Temp(n)>T) exit
    enddo

    if(n==1) then
      k1=k0(1)
      k11=k00(1)
    elseif(n==6) then
      k1=k0(5)
      k11=k00(5)
    else
      k1=k0(n-1)+(k0(n)-k0(n-1))/(Temp(n)-Temp(n-1))*(Temp(n)-T)
      k11=k00(n-1)+(k00(n)-k00(n-1))/(Temp(n)-Temp(n-1))*(Temp(n)-T)
    endif

    p=1./(1.+dlog10(k1*rhoNum/k11)*dlog10(k1*rhoNum/k11))
    rate_HO2=k1*rhoNum/(1.+k1*rhoNum/k11)*0.6**p

  end function rate_HO2

end subroutine Atm_Hg_GasChem_Br2


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating the right-hand side of the chemical ODE system
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Chem_Deriv(N,Y,F)

    implicit none

    integer N, i, j
    real(8) Y(N), F(N)

    F=0.
    do i=1, N
      do j=1, N
        F(i)=F(i)+JacobChem(i,j)*Y(j)
      enddo
    enddo

end subroutine Chem_Deriv


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating the Jacobian of the chemical ODE system
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Jacob_Br_Chem(N,Y,dFy,LdFy)

    implicit none

    integer N, LdFy
    real(8) Y(N), dFy(LdFy,N)

    JacobChem=0.
    JacobChem(1,1)=-(Krate(1)+Krate(13))
    JacobChem(2,1)=Krate(1)
    JacobChem(7,1)=Krate(13)
    JacobChem(1,2)=Krate(2)+Krate(3)+Krate(4)
    JacobChem(2,2)=-(Krate(2)+Krate(3)+Krate(4))-(Krate(5)+Krate(6)+Krate(7)+Krate(8))
    JacobChem(3,2)=Krate(5)
    JacobChem(4,2)=Krate(6)
    JacobChem(5,2)=Krate(7)
    JacobChem(6,2)=Krate(8)
    JacobChem(1,3)=0.58*Krate(9)
    JacobChem(2,3)=0.42*Krate(9)
    JacobChem(3,3)=-Krate(9)
    JacobChem(1,4)=0.5*Krate(10)
    JacobChem(2,4)=0.15*Krate(10)
    JacobChem(4,4)=-Krate(10)
    JacobChem(7,4)=0.35*Krate(10)
    JacobChem(1,5)=0.66*Krate(11)
    JacobChem(2,5)=0.03*Krate(11)
    JacobChem(4,5)=0.31*Krate(11)
    JacobChem(5,5)=-Krate(11)
    JacobChem(2,6)=0.1*Krate(12)
    JacobChem(4,6)=0.9*Krate(12)
    JacobChem(6,6)=-Krate(12)
    JacobChem(1,7)=Krate(14)+Krate(15)
    JacobChem(7,7)=-(Krate(14)+Krate(15))-(Krate(16)+Krate(17))
    JacobChem(8,7)=Krate(16)
    JacobChem(9,7)=Krate(17)
    JacobChem(7,8)=Krate(18)
    JacobChem(8,8)=-Krate(18)
    JacobChem(7,9)=Krate(19)
    JacobChem(9,9)=-Krate(19)

    dFy(1:N,1:N)=JacobChem(1:N,1:N)

end subroutine Jacob_Br_Chem


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating the Jacobian of the chemical ODE system
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Jacob_Br_Chem2(N,Y,dFy,LdFy)

    implicit none

    integer N, LdFy
    real(8) Y(N), dFy(LdFy,N)

    JacobChem=0.
    JacobChem(1,1)=-Krate(1)
    JacobChem(2,1)=Krate(1)
    JacobChem(1,2)=Krate(2)+Krate(3)
    JacobChem(2,2)=-(Krate(2)+Krate(3))-(Krate(5)+Krate(6))
    JacobChem(3,2)=Krate(5)+Krate(6)

    dFy(1:N,1:N)=JacobChem(1:N,1:N)

end subroutine Jacob_Br_Chem2


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading monthly fields of chemical reactants from GEOSChem
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_ReadReactDaily_GEOSChem

    integer status, ncid_in, var_id_in, start3(4), count3(4), t
    real qreact_buf( Imax, Jmax, Atm_Kmax, 1:4)
    data start3 /1,1,1,1/
    data count3 /Imax, Jmax, Atm_Kmax, 4/

    integer(2) ii, jj
    integer FileStat, i, j, k, t, Xscal, daF(2), mnF(2), yrF(2), nDay, flag, yrCur
    real qO3(NumPer), qSO2(NumPer), qOH(NumPer), qBr(NumPer), qBrO(NumPer), qPM(NumPer)
    real splRow(NumPer*2), Aver(Imin:Imax), dFdX(NumPer*2), d2FdX(NumPer*2), dFdXleft

! Check for climatic run
    if(climReactRun) then
      yrCur=ClimReactYear
    else
      yrCur=Year
    endif
    
    write(YearNum,'(i4)') yrCur
    write(DayNum,'(i2.2)') Day

    daF(1)=Day
    mnF(1)=Month
    yrF(1:2)=yrCur
    if(Day<MonthDays(Month)) then
      daF(2)=Day+1
      mnF(2)=Month
    elseif(Month<12) then
      daF(2)=1
      mnF(2)=Month+1
    else
      daF(2)=1
      mnF(2)=1
    endif

    do nDay=1, 2

! Reading O3 distribution
      write(fileName,'(a,i4,i2.2,i2.2,a4)') trim(OzoneName), yrF(nDay), mnF(nDay), daF(nDay), '.nc4'
      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      status=nf90_open(trim(fullName), nf90_nowrite, ncid_in)
      if(status/=nf90_noerr) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      
      call checkNC(nf90_inq_varid(ncid_in,'O3',var_id_in),8)
      call checkNC(nf90_get_var(ncid_in,var_id_in,qReact_buf,start3,count3), 9)
      
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            do t = 1, NumPer
              O3field(i,j,k,t,nDay)=qReact_buf( i, j, k, t ) 
            enddo
          enddo
        enddo
      enddo
      where(O3field<=0.) O3field=Zero
      call checkNC(nf90_close(ncid_in), 10)

! Reading OH distribution
      write(fileName,'(a,i4,i2.2,i2.2,a4)') trim(OHName), yrF(nDay), mnF(nDay), daF(nDay), '.nc4'
      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      status=nf90_open(trim(fullName), nf90_nowrite, ncid_in)
      if(status/=nf90_noerr) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)              
        stop
      endif
      
      call checkNC(nf90_inq_varid(ncid_in,'OH',var_id_in),8)
      call checkNC(nf90_get_var(ncid_in,var_id_in,qReact_buf,start3,count3), 9)
      
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            do t = 1, NumPer
              OHfield(i,j,k,t,nDay)=qReact_buf( i, j, k, t )
            enddo
          enddo
        enddo
      enddo
      where(OHfield<=0.) OHfield=Zero
      call checkNC(nf90_close(ncid_in), 11)

! Reading SO2 distribution
      write(fileName,'(a,i4,i2.2,i2.2,a4)') trim(SO2Name), yrF(nDay), mnF(nDay), daF(nDay), '.nc4'
      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      status=nf90_open(trim(fullName), nf90_nowrite, ncid_in)
      if(status/=nf90_noerr) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      
      call checkNC(nf90_inq_varid(ncid_in,'SO2',var_id_in),8)
      call checkNC(nf90_get_var(ncid_in,var_id_in,qReact_buf,start3,count3), 9)
      
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            do t = 1, NumPer
              SO2field(i,j,k,t,nDay)=qReact_buf( i, j, k, t )
            enddo
          enddo
        enddo
      enddo
      where(SO2field<=0.) SO2field=Zero
      call checkNC(nf90_close(ncid_in), 12)

! Reading PM2.5 distribution
      write(fileName,'(a,i4,i2.2,i2.2,a4)') 'PM2.5_', yrF(nDay), mnF(nDay), daF(nDay), '.nc4'
      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      status=nf90_open(trim(fullName), nf90_nowrite, ncid_in)
      if(status/=nf90_noerr) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      
      call checkNC(nf90_inq_varid(ncid_in,'PM2.5',var_id_in),8)
      call checkNC(nf90_get_var(ncid_in,var_id_in,qReact_buf,start3,count3), 9)
      
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            do t = 1, NumPer
              PMfield(i,j,k,t,nDay)=qReact_buf( i, j, k, t )
            enddo
          enddo
        enddo
      enddo
      where(PMfield<=0.) PMfield=Zero
      call checkNC(nf90_close(ncid_in), 13)
    
! Reading NO2 distribution
      write(fileName,'(a,i4,i2.2,i2.2,a4)') 'NO2_', yrF(nDay), mnF(nDay), daF(nDay), '.nc4'
      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      status=nf90_open(trim(fullName), nf90_nowrite, ncid_in)
      if(status/=nf90_noerr) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      
      call checkNC(nf90_inq_varid(ncid_in,'NO2',var_id_in),8)
      call checkNC(nf90_get_var(ncid_in,var_id_in,qReact_buf,start3,count3), 9)
      
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            do t = 1, NumPer
              NO2field(i,j,k,t,nDay)=qReact_buf( i, j, k, t )
            enddo
          enddo
        enddo
      enddo
      where(SO2field<=0.) NO2field=Zero
      call checkNC(nf90_close(ncid_in), 12)

! Reading HO2 distribution
      write(fileName,'(a,i4,i2.2,i2.2,a4)') 'HO2_', yrF(nDay), mnF(nDay), daF(nDay), '.nc4'
      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      status=nf90_open(trim(fullName), nf90_nowrite, ncid_in)
      if(status/=nf90_noerr) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      
      call checkNC(nf90_inq_varid(ncid_in,'HO2',var_id_in),8)
      call checkNC(nf90_get_var(ncid_in,var_id_in,qReact_buf,start3,count3), 9)
      
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            do t = 1, NumPer
              HO2field(i,j,k,t,nDay)=qReact_buf( i, j, k, t )
            enddo
          enddo
        enddo
      enddo
      where(SO2field<=0.) HO2field=Zero
      call checkNC(nf90_close(ncid_in), 12)

! Reading Br distribution 
!      write(fileName,'(a,i4,i2.2,i2.2,a4)') 'Br_', yrF(nDay), mnF(nDay), daF(nDay), '.nc4'
!      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
!      status=nf90_open(trim(fullName), nf90_nowrite, ncid_in)
!      if(status/=nf90_noerr) then
!        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
!        stop
!      endif
!       
!      call checkNC(nf90_inq_varid(ncid_in,'Br',var_id_in),8)
!      call checkNC(nf90_get_var(ncid_in,var_id_in,qReact_buf,start3,count3), 9)
!      
!      do k=1, Atm_Kmax
!        do j=Jmin, Jmax
!          do i=Imin, Imax
!            do t = 1, NumPer
!              Brfield(i,j,k,t,nDay)=qReact_buf( i, j, k, t )
!            enddo
!          enddo
!        enddo
!      enddo
!      where(Brfield<=0.) Brfield=Zero
!      call checkNC(nf90_close(ncid_in), 13)
      
! ----------------------------- Experimental -------------------------------------------
      write(fileName,'(a,i4,i2.2,i2.2,a4)') 'Br_', 2010, mnF(nDay), daF(nDay), '.bin'
      fullName=trim(ReactPath2)//trim(GridCode)//'/'//'2010/'//trim(fileName)
      open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
      if(FileStat>0) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            read(10) (qBr(t), t=1, NumPer)          ! ppbv
            if(k>6) cycle                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Brfield(i,j,k,1:NumPer,nDay)=qBr(1:NumPer)
          enddo
        enddo
      enddo
      where(Brfield<=0.) Brfield=Zero
      close(10)

      do j=Jmin, Jmax                                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=Imin, Imax                                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(sum(Brfield(i,j,1,1:NumPer,nDay))/NumPer<4.e-4) then  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Brfield(i,j,1:6,1:NumPer,nDay)=0.                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          endif                                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo                                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo                                               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ----------------------------- Experimental -------------------------------------------      

    enddo
   
! Grid aggregation
    do j=Jmin, Jmax
      if(maxI(j)==1) cycle
      Xscal=Imax/maxI(j)
      if(Xscal==1) cycle

      do nDay=1, 2
        do t=1, NumPer
          do k=1, Atm_Kmax
            Aver(Imin:Imax)=O3field(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            O3field(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))

            Aver(Imin:Imax)=SO2field(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            SO2field(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))

            Aver(Imin:Imax)=OHfield(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            OHfield(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
            
            Aver(Imin:Imax)=Brfield(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            Brfield(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))

            Aver(Imin:Imax)=NO2field(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            NO2field(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))

            Aver(Imin:Imax)=HO2field(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            HO2field(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))

            Aver(Imin:Imax)=PMfield(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            PMfield(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
          enddo
        enddo
      enddo
    enddo

    do j=Jmin, Jmax
      do i=minI(j), maxI(j)
        do k=1, Atm_Kmax
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=O3field(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=O3field(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dO3field(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dO3field(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2O3field(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=SO2field(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=SO2field(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dSO2field(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dSO2field(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2SO2field(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=OHfield(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=OHfield(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dOHfield(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dOHfield(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2OHfield(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=Brfield(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Brfield(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dBrfield(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dBrfield(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2Brfield(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=NO2field(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=NO2field(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dNO2field(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dNO2field(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2NO2field(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=HO2field(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=HO2field(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dHO2field(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dHO2field(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2HO2field(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=PMfield(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=PMfield(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dPMfield(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dPMfield(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2PMfield(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
        enddo
      enddo
    enddo

end subroutine Atm_ReadReactDaily_GEOSChem


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading distribution monthly fields of chemical reactants
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_ReadReactDaily_MOZART

    integer(2) ii, jj
    integer FileStat, i, j, k, t, Xscal, daF(2), mnF(2), yrF(2), nDay, flag, yrCur
    real qO3(NumPer), qSO2(NumPer), qOH(NumPer), qBr(NumPer), qBrO(NumPer), qPM(NumPer)
    real splRow(NumPer*2), Aver(Imin:Imax), dFdX(NumPer*2), d2FdX(NumPer*2), dFdXleft

! Check for climatic run
    if(climReactRun) then
      yrCur=ClimReactYear
    else
      yrCur=Year
    endif

    write(YearNum,'(i4)') yrCur
    write(DayNum,'(i2.2)') Day

    daF(1)=Day
    mnF(1)=Month
    yrF(1:2)=yrCur
    if(Day<MonthDays(Month)) then
      daF(2)=Day+1
      mnF(2)=Month
    elseif(Month<12) then
      daF(2)=1
      mnF(2)=Month+1
    else
      daF(2)=1
      mnF(2)=1
    endif

    do nDay=1, 2

! Reading O3 distribution
      write(fileName,'(a,i4,i2.2,i2.2,a4)') trim(OzoneName), yrF(nDay), mnF(nDay), daF(nDay), '.bin'
      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
      if(FileStat>0) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            read(10) (qO3(t), t=1, NumPer)          ! ppbv
            O3field(i,j,k,1:NumPer,nDay)=qO3(1:NumPer)
          enddo
        enddo
      enddo
      where(O3field<=0.) O3field=Zero
      close(10)

! Reading OH distribution
      write(fileName,'(a,i4,i2.2,i2.2,a4)') trim(OHName), yrF(nDay), mnF(nDay), daF(nDay), '.bin'
      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
      if(FileStat>0) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            read(10) (qOH(t), t=1, NumPer)          ! ppbv
            OHfield(i,j,k,1:NumPer,nDay)=qOH(1:NumPer)
          enddo
        enddo
      enddo
      where(OHfield<=0.) OHfield=Zero
      close(10)

! Reading SO2 distribution
      write(fileName,'(a,i4,i2.2,i2.2,a4)') trim(SO2Name), yrF(nDay), mnF(nDay), daF(nDay), '.bin'
      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
      if(FileStat>0) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            read(10) (qSO2(t), t=1, NumPer)          ! ppbv
            SO2field(i,j,k,1:NumPer,nDay)=qSO2(1:NumPer)
          enddo
        enddo
      enddo
      where(SO2field<=1.e-5) SO2field=1.e-5
      close(10)

! Reading PM2.5 distribution
      write(fileName,'(a,i4,i2.2,i2.2,a4)') 'PM2.5_', yrF(nDay), mnF(nDay), daF(nDay), '.bin'
      fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
      if(FileStat>0) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            read(10) (qPM(t), t=1, NumPer)          ! ppbm
            PMfield(i,j,k,1:NumPer,nDay)=qPM(1:NumPer)
          enddo
        enddo
      enddo
      where(PMfield<=0.) PMfield=Zero
      close(10)

! Reading Br distribution
! ----------------------------- Experimental -------------------------------------------
!      write(fileName,'(a,i4,i2.2,i2.2,a4)') 'Br_', yrF(nDay), mnF(nDay), daF(nDay), '.bin'
!      fullName=trim(ReactPath2)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
      write(fileName,'(a,i4,i2.2,i2.2,a4)') 'Br_', 2010, mnF(nDay), daF(nDay), '.bin'
      fullName=trim(ReactPath2)//trim(GridCode)//'/'//'2010/'//trim(fileName)
      open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
      if(FileStat>0) then
        print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
        stop
      endif
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
          do i=Imin, Imax
            read(10) (qBr(t), t=1, NumPer)          ! ppbv
            if(k>6) cycle                                       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Brfield(i,j,k,1:NumPer,nDay)=qBr(1:NumPer)
          enddo
        enddo
      enddo
      where(Brfield<=0.) Brfield=Zero
      close(10)

      do j=Jmin, Jmax                                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=Imin, Imax                                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(sum(Brfield(i,j,1,1:NumPer,nDay))/NumPer<4.e-4) then  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Brfield(i,j,1:6,1:NumPer,nDay)=0.                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          endif                                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo                                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo                                               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ----------------------------- Experimental -------------------------------------------      

    enddo

! Grid aggregation
    do j=Jmin, Jmax
      if(maxI(j)==1) cycle
      Xscal=Imax/maxI(j)
      if(Xscal==1) cycle

      do nDay=1, 2
        do t=1, NumPer
          do k=1, Atm_Kmax
            Aver(Imin:Imax)=O3field(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            O3field(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))

            Aver(Imin:Imax)=SO2field(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            SO2field(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))

            Aver(Imin:Imax)=OHfield(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            OHfield(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))

            Aver(Imin:Imax)=Brfield(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            Brfield(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))

            Aver(Imin:Imax)=PMfield(Imin:Imax,j,k,t,nDay)
            call GridAggreg(j,Xscal,Aver,1)
            PMfield(minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
          enddo
        enddo
      enddo
    enddo

    do j=Jmin, Jmax
      do i=minI(j), maxI(j)
        do k=1, Atm_Kmax
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=O3field(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=O3field(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dO3field(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dO3field(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2O3field(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=SO2field(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=SO2field(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dSO2field(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dSO2field(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2SO2field(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=OHfield(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=OHfield(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dOHfield(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dOHfield(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2OHfield(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=Brfield(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Brfield(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dBrfield(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dBrfield(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2Brfield(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
          splRow(1:NumPer)=PMfield(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=PMfield(i,j,k,1:NumPer,toMor)
          if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
            flag=0
            dFdXleft=0.
          else
            flag=1
            dFdXleft=dPMfield(i,j,k,NumPer+1)
          endif
          call SplineParams(NumPer*2,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
          dPMfield(i,j,k,1:NumPer+1)=dFdX(1:NumPer+1)
          d2PMfield(i,j,k,1:NumPer+1)=d2FdX(1:NumPer+1)
!-----------------------------------------------------------------------------------
        enddo
      enddo
    enddo

end subroutine Atm_ReadReactDaily_MOZART


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading distribution monthly fields of chemical reactants
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadReactMonthly

    integer(2) ii, jj
    integer FileStat, i, j, k, Xscal
    real Ph_rate
    real Aver(Imin:Imax)


! Reading photolysis rate for HgBr2
    write(fileName,'(a,i4,i2.2,a4)') 'Khgbr2_', 2013, Month, '.bin'
    fullName='/home/alex/InputData/Reactants/Photolysis/'//trim(GridCode)//'/'//'2013'//'/'//trim(fileName)
    open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
      stop
    endif
    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        do i=Imin, Imax
          read(10) Ph_rate                        ! 1/s
          PhotoRate(HgBr2,i,j,k)=max(real(Ph_rate,8),real(Zero,8))
        enddo
      enddo
    enddo
    close(10)

! Reading photolysis rate for HgBrOH
    write(fileName,'(a,i4,i2.2,a4)') 'Khgbroh_', 2013, Month, '.bin'
    fullName='/home/alex/InputData/Reactants/Photolysis/'//trim(GridCode)//'/'//'2013'//'/'//trim(fileName)
    open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
      stop
    endif
    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        do i=Imin, Imax
          read(10) Ph_rate                        ! 1/s
          PhotoRate(HgBrOH,i,j,k)=max(real(Ph_rate,8),real(Zero,8))
        enddo
      enddo
    enddo
    close(10)

! Reading photolysis rate for HgBrOOH
    write(fileName,'(a,i4,i2.2,a4)') 'Khgbrooh_', 2013, Month, '.bin'
    fullName='/home/alex/InputData/Reactants/Photolysis/'//trim(GridCode)//'/'//'2013'//'/'//trim(fileName)
    open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
      stop
    endif
    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        do i=Imin, Imax
          read(10) Ph_rate                        ! 1/s
          PhotoRate(HgBrOOH,i,j,k)=max(real(Ph_rate,8),real(Zero,8))
        enddo
      enddo
    enddo
   close(10)

! Reading photolysis rate for HgBrONO
!    write(fileName,'(a,i4,i2.2,a4)') 'Khgbrno2_', 2013, Month, '.bin'
    write(fileName,'(a,i4,i2.2,a4)') 'Khgbrono_syn_', 2013, Month, '.bin'
    fullName='/home/alex/InputData/Reactants/Photolysis/'//trim(GridCode)//'/'//'2013'//'/'//trim(fileName)
    open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
      stop
    endif
    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        do i=Imin, Imax
          read(10) Ph_rate                        ! 1/s
          PhotoRate(HgBrONO,i,j,k)=max(real(Ph_rate,8),real(Zero,8))
        enddo
      enddo
    enddo
    close(10)

! Reading photolysis rate for HgBr
    write(fileName,'(a,i4,i2.2,a4)') 'Khgbr_', 2013, Month, '.bin'
    fullName='/home/alex/InputData/Reactants/Photolysis/'//trim(GridCode)//'/'//'2013'//'/'//trim(fileName)
    open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
      stop
    endif
    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        do i=Imin, Imax
          read(10) Ph_rate                        ! 1/s
          PhotoRate(HgBr,i,j,k)=max(real(Ph_rate,8),real(Zero,8))
        enddo
      enddo
    enddo
    close(10)
    
! Reading photolysis rate for HOHg
    write(fileName,'(a,i4,i2.2,a4)') 'Khgoh_', 2013, Month, '.bin'
    fullName='/home/alex/InputData/Reactants/Photolysis/'//trim(GridCode)//'/'//'2013'//'/'//trim(fileName)
    open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
      stop
    endif
    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        do i=Imin, Imax
          read(10) Ph_rate                        ! 1/s
          PhotoRate(HOHg,i,j,k)=max(real(Ph_rate,8),real(Zero,8))
        enddo
      enddo
    enddo
    close(10)
    
! Reading photolysis rate for HgO
!    write(fileName,'(a,i4,i2.2,a4)') 'Khgo_', 2013, Month, '.bin'
!    fullName='/home/alex/InputData/Reactants/Photolysis/'//trim(GridCode)//'/'//'2013'//'/'//trim(fileName)
!    open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read') !
!    if(FileStat>0) then
!      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
!      stop
!    endif
!    do k=1, Atm_Kmax
!      do j=Jmin, Jmax
!        do i=Imin, Imax
!          read(10) Ph_rate                        ! 1/s
!          PhotoRate(HgO,i,j,k)=max(real(Ph_rate,8),real(Zero,8))
!        enddo
!      enddo
!    enddo
!    close(10)
    
! Grid aggregation
    do j=Jmin, Jmax
      if(maxI(j)==1) cycle
      Xscal=Imax/maxI(j)
      if(Xscal==1) cycle

      do k=1, Atm_Kmax
        Aver(Imin:Imax)=PhotoRate(HgBr2,Imin:Imax,j,k)
        call GridAggreg(j,Xscal,Aver,1)
        PhotoRate(HgBr2,minI(j):maxI(j),j,k)=Aver(minI(j):maxI(j))

        Aver(Imin:Imax)=PhotoRate(HgBrOH,Imin:Imax,j,k)
        call GridAggreg(j,Xscal,Aver,1)
        PhotoRate(HgBrOH,minI(j):maxI(j),j,k)=Aver(minI(j):maxI(j))

        Aver(Imin:Imax)=PhotoRate(HgBrOOH,Imin:Imax,j,k)
        call GridAggreg(j,Xscal,Aver,1)
        PhotoRate(HgBrOOH,minI(j):maxI(j),j,k)=Aver(minI(j):maxI(j))

        Aver(Imin:Imax)=PhotoRate(HgBrONO,Imin:Imax,j,k)
        call GridAggreg(j,Xscal,Aver,1)
        PhotoRate(HgBrONO,minI(j):maxI(j),j,k)=Aver(minI(j):maxI(j))

        Aver(Imin:Imax)=PhotoRate(HgBr,Imin:Imax,j,k)
        call GridAggreg(j,Xscal,Aver,1)
        PhotoRate(HgBr,minI(j):maxI(j),j,k)=Aver(minI(j):maxI(j))
        
        Aver(Imin:Imax)=PhotoRate(HOHg,Imin:Imax,j,k)
        call GridAggreg(j,Xscal,Aver,1)
        PhotoRate(HOHg,minI(j):maxI(j),j,k)=Aver(minI(j):maxI(j))
        
!        Aver(Imin:Imax)=PhotoRate(HgO,Imin:Imax,j,k)
!        call GridAggreg(j,Xscal,Aver,1)
!        PhotoRate(HgO,minI(j):maxI(j),j,k)=Aver(minI(j):maxI(j))
      enddo
    enddo

end subroutine ReadReactMonthly


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine transforming reactants form [ppb] to [molec/cm3]
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReactTransformSteply

    integer i, j, k, t
    real(8) RhoAirNum, RhoAir
    real Tair, Psur
    real splRow(NumPer*2), d2FdX(NumPer*2), splVal, ReactMixR

    do j=Jmin, Jmax
      do i=minI(j), maxI(j)
        Psur=PxCurr(i,j)
        do k=1, Atm_Kmax
          RhoAir=DensAir(i,j,k)                                            ! kg/m3
          Tair=TairCurr(i,j,k)
          RhoAirNum=Nav/Runiv*(Sigma(k)*Psur+Ptop)/Tair*1.e-6              ! molec/cm3
          DensAirNum(i,j,k)=RhoAirNum

          splRow(1:NumPer)=O3field(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=O3field(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2O3field(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          ReactMixR=max(splVal,min(splRow(Period),splRow(Period+1)))
          ConcO3(i,j,k)=real(ReactMixR*1.e-9*RhoAirNum,8)                  ! ppbv -> molec/cm3

          splRow(1:NumPer)=SO2field(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=SO2field(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2SO2field(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          ReactMixR=max(splVal,0.00001)
          ConcSO2(i,j,k)=real(ReactMixR*1.e-9*RhoAirNum,8)                 ! ppbv -> molec/cm3

          splRow(1:NumPer)=OHfield(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=OHfield(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2OHfield(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          ReactMixR=max(splVal,min(splRow(Period),splRow(Period+1)))
          ConcOH(i,j,k)=real(ReactMixR*1.e-9*RhoAirNum,8)                  ! ppbv -> molec/cm3

          splRow(1:NumPer)=Brfield(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=Brfield(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2Brfield(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          ReactMixR=max(splVal,min(splRow(Period),splRow(Period+1)))
          ConcBr(i,j,k)=real(ReactMixR*1.e-9*RhoAirNum,8)                  ! ppbv -> molec/cm3

          splRow(1:NumPer)=NO2field(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=NO2field(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2NO2field(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          ReactMixR=max(splVal,min(splRow(Period),splRow(Period+1)))
          ConcNO2(i,j,k)=real(ReactMixR*1.e-9*RhoAirNum,8)                 ! ppbv -> molec/cm3

          splRow(1:NumPer)=HO2field(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=HO2field(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2HO2field(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          ReactMixR=max(splVal,min(splRow(Period),splRow(Period+1)))
          ConcHO2(i,j,k)=real(ReactMixR*1.e-9*RhoAirNum,8)                 ! ppbv -> molec/cm3

          splRow(1:NumPer)=PMfield(i,j,k,1:NumPer,toDay)
          splRow(NumPer+1:NumPer*2)=PMfield(i,j,k,1:NumPer,toMor)
          d2FdX(1:NumPer+1)=d2PMfield(i,j,k,1:NumPer+1)
          splVal=SplineInterpol(NumPer*2,timePer,splRow,d2FdX,Period,DayTime)
          ReactMixR=max(splVal,min(splRow(Period),splRow(Period+1)))
          ConcPM(i,j,k)=real(ReactMixR*1.e-9*RhoAir,8)                      ! ppbm -> kg/m3
        enddo
      enddo
    enddo

end subroutine ReactTransformSteply


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine transforming reactants form [ppb] to [molec/cm3]
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReactTransformDaily

    integer i, j, k, t
    real(8) RhoAir, Tsurf, Tair1, Psur, SeaRelArea

    do t=1, NumPer
      do j=Jmin, Jmax
        do i=minI(j), maxI(j)
          Psur=Px(i,j,t,toDay)
          Tair1=TempAir(i,j,1,t,toDay)
          RhoAir=Nav/Runiv*(Sigma(1)*Psur+Ptop)/Tair1*1.e-15            ! molec/cm3
          SeaRelArea=sum(LandCover(i,j,gInd(Water,1:gNum(Water))))

          Tsurf=TempSurf(i,j,t)
          if(Tsurf>273.) then
            ConcCl2(i,j,1,t)=0.1*SeaRelArea*RhoAir                      ! ppb -> molec/cm3
          else
            ConcCl2(i,j,1,t)=0.
          endif
          ConcCl2(i,j,2:Atm_Kmax,t)=0.
        enddo
      enddo
    enddo

end subroutine ReactTransformDaily


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating Hg gas-particle partitioning
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_Hg_Partition

    integer i, j, k, Ind, fOxid, fPart
    real(8) Coxid, Cpart, dC, Apart, Kpart, T, meshV

    do k=1, Atm_Kmax
      do j=Jmin, Jmax
        do i=minI(j), maxI(j)
          meshV=MeshVolume(i,j,k)
          T=TairCurr(i,j,k)
          Kpart=10._8**(-10._8+2500._8/T)
          Apart=Kpart*ConcPM(i,j,k)*1.e9
          
          do Ind=1, Noxid
            fOxid=GasPart(1,Ind)          
            fPart=GasPart(2,Ind)
          
            Coxid=Atm_Conc(i,j,k,fOxid)
            Cpart=Atm_Conc(i,j,k,fPart)
          
            Atm_Conc(i,j,k,fOxid)=(Coxid+Cpart)/(1._8+Apart)
            Atm_Conc(i,j,k,fPart)=(Coxid+Cpart)*Apart/(1._8+Apart)
            
            dC=(Cpart-Coxid*Apart)/(1._8+Apart)
            MassChemEx(fOxid)=MassChemEx(fOxid)+dC*meshV
            MassChemEx(fPart)=MassChemEx(fPart)-dC*meshV
          enddo
        enddo
      enddo
    enddo

end subroutine Atm_Hg_Partition


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
!     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y).
!     THIS IS AN EMBEDDED ROSENBROCK METHOD OF ORDER (3)4
!     (WITH STEP SIZE CONTROL).
!     !.F. SECTIONS IV.7  AND VI.3
!
!     AUTHORS: E. HAIRER AND G. WANNER
!              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
!              CH-1211 GENEVE 24, SWITZERLAND
!              E-MAIL:  Ernst.Hairer@math.unige.ch
!                       Gerhard.Wanner@math.unige.ch
!
!     THIS CODE IS PART OF THE BOOK:
!         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
!         SPRINGER-VERLAG 1991, SECOND EDITION 1996.
!
!     VERSION OF OCTOBER 28, 1996 (simplified on December 17, 2020)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RoDas_ODE(N,X,Y,Xend,FcN,Jac,Idid)

    implicit none
    
    integer N, Itol, IJac, MLJac, MUJac, LWork, Idid
    real(8) X, Xend, H
    logical Jband, Arret, Pred
    integer Nfcn, Naccpt, Nrejct, Nstep, Njac, Ndec, Nsol, Nmax, Meth, m1, m2, nm1
    integer i, LDJac, LDE, LDMas, Ijob, LDMas2, Istore, IEIP
    integer IEYNEW, IEDY1, IEDY, IEAk1, IEAk2, IEAk3, IEAk4, IEAk5, IEAk6, IEFX, IECON, IEJAC, IEMAS, Iee
    real(8) Y(N), Atol(N), Rtol(N), Uround, Hmax, Fac1, Fac2, Safe
    integer MBJac, Mdiag, MLE
    real(8) Work(2*N*N+14*N+20)
!    real(8), allocatable :: Work(:)

    external FcN, Jac
 
    Uround=1.D-16            ! Smallest number satisfying 1.D0+Uround>1.D0  
    Itol=0                   ! Rtol and Atol are scalars
    Rtol=1.0D-6 !1.0D-4      ! Relative error tolerance (Rtol(1)>10.*Uround)
    Atol=1.0D-6*Rtol         ! Absolute error tolerance (Atol(1)>0.)
    h=1.0D-6  !1.0D-6         ! Initial step size
    Ijac=1                   ! Compute the Jacobian analytically
    MLJac=N                  ! Jacobian is a full matrix
    MUJac=0                  ! Need not be defined if MLJac=N
    Nmax=100000              ! Maximal number of steps
    Meth=1                   ! Coefficients of the method
    Pred=.true.              ! Predictive step size control (GUSTAFSSON)
    Fac1=5.D0                ! PARAMETERS FOR STEP SIZE SELECTION
    Fac2=1.D0/6.0D0          !
    Safe=0.9D0               ! SAFETY FACTOR IN STEP SIZE PREDICTION
    Hmax=Xend-X              ! MAXIMAL STEP SIZE
    
! Zeroing working arrays
!    LWork=2*N*N+14*N+20
!    allocate(Work(LWork))
    Work(1:20)=0.D0

    m1=0.                    ! PARAMETER FOR SECOND ORDER EQUATIONS
    m2=0.                    !
    nm1=N-m1
    if(m1==0) m2=N
    if(m2==0) m2=m1
    
! COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS
    Jband=(MLJac<nm1)
    if(Jband) then
      LDJac=MLJac+MUJac+1
      LDE=MLJac+LDJac
    else
      MLJac=nm1
      MUJac=nm1
      LDJac=nm1
      LDE=nm1
    endif
    LDMas=0
    if(Jband)  then
      Ijob=2
    else
      Ijob=1
    endif
    LDMas2=max(1,LDMas)
      
! Setting the parameters
    Nfcn=0
    Naccpt=0
    Nrejct=0
    Nstep=0
    Njac=0
    Ndec=0
    Nsol=0
    Arret=.false.
      
! PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN Work
    IEYNEW=21
    IEDY1=IEYNEW+N
    IEDY=IEDY1+N
    IEAk1=IEDY+N
    IEAk2=IEAk1+N
    IEAk3=IEAk2+N
    IEAk4=IEAk3+N
    IEAk5=IEAk4+N
    IEAk6=IEAk5+N
    IEFX =IEAk6+N
    IECON=IEFX+N
    IEJAC=IECON+4*N
    IEMAS=IEJAC+N*LDJac
    Iee  =IEMAS+nm1*LDMas
      
! -------- call TO CORE INTEGRATOR ------------
    call RosCore(N,FcN,X,Y,Xend,Hmax,H,Rtol,Atol,Itol,Jac,IJac,MLJac,MUJac,Idid,&
                & Nmax,Uround,Meth,Ijob,Fac1,Fac2,Safe,Jband,Pred,LDJac,&
                & LDE,Work(IEYNEW),Work(IEDY1),Work(IEDY),Work(IEAk1),&
                & Work(IEAk2),Work(IEAk3),Work(IEAk4),Work(IEAk5),Work(IEAk6),&
                & Work(IEFX),Work(IEJAC),Work(Iee),Work(IECON),&
                & m1,m2,nm1,Nfcn,Njac,Nstep,Naccpt,Nrejct,Ndec,Nsol)
     
!    deallocate(Work)

 contains      
      
!......................................................................................................      
! ----------------------------------------------------------
!     CORE INTEGRATOR FOR RoDas_ODE
!     PARAMETERS SAME AS IN RoDas_ODE WITH WORKSPACE ADDED
! ---------------------------------------------------------- 
  subroutine RosCore(N,FcN,X,Y,Xend,Hmax,H,Rtol,Atol,Itol,Jac,IJac,MLJac,MUJac,Idid,&
                 & Nmax,Uround,Meth,Ijob,Fac1,Fac2,Safe,Banded,Pred,LDJac,&
                 & LDE,Ynew,DY1,DY,Ak1,Ak2,Ak3,Ak4,Ak5,Ak6,Fx,FJac,E,CONT,&
                 & m1,m2,nm1,Nfcn,Njac,Nstep,Naccpt,Nrejct,Ndec,Nsol)

    implicit none

    integer N, Itol, IJac, MLJac, MUJac, Idid, Nmax, Meth, Ijob
    integer LDJac, LDE,  m1, m2, nm1, Nfcn, Njac, Nstep, Naccpt, Nrejct, Ndec, Nsol, Ip(nm1)
    real(8) X, Xend, Hmax, H, Uround, Fac1, Fac2, Safe 
    real(8) Ak1(N), Ak2(N), Ak3(N), Ak4(N), Ak5(N), Ak6(N), Fx(N), FJac(LDJac,N), E(LDE,nm1)
    real(8) Y(N), Rtol(*), Atol(*), Ynew(N), DY1(N), DY(N), CONT(4*N)
    logical Reject, Banded, Last, Pred
    integer NN, NN2, NN3, LRC, Nsing, Irtrn, MLE, MUE, MUJacP, MUJacJ, MD, mm, k, j, j1, i
    integer Lbeg, Lend, L, Ier, Nsk
    real(8) A21,A31,A32,A41,A42,A43,A51,A52,A53,A54,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61
    real(8) C62,C63,C64,C65,Gamma,C2,C3,C4,D1,D2,D3,D4,D21,D22,D23,D24,D25,D31,D32,D33,D34,D35
    real(8) PosNeg, HmaxN, Hopt, Hout, Xold, Ysafe, Delt, Fac, Err, SK, Hnew, FacGus, Hacc, ErrAcc
    real(8) HC21, HC31, HC32, HC41, HC42, HC43, HC51, HC52, HC53, HC54, HC61, HC62, HC63, HC64, HC65
      
    NN=N 
    NN2=2*N
    NN3=3*N
    LRC=4*N

! ------ SET THE PARAMETERS OF THE METHOD -----
    call RoCoef(Meth,A21,A31,A32,A41,A42,A43,A51,A52,A53,A54,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,&
           &C62,C63,C64,C65,Gamma,C2,C3,C4,D1,D2,D3,D4,D21,D22,D23,D24,D25,D31,D32,D33,D34,D35)

! --- INITIAL PREPARATIONS
    if(m1>0) Ijob=Ijob+10
    PosNeg=sign(1.D0,Xend-X)
    HmaxN=min(abs(Hmax),abs(Xend-X))
    if(abs(H)<=10.D0*Uround) H=1.0D-6
    H=min(abs(H),HmaxN)
    H=sign(H,PosNeg)
    Reject=.false.
    Last=.false.
    Nsing=0
    Irtrn=1

! -------- PREPARE BAND-WIDTHS --------
    if(Banded)  then
      MLE=MLJac
      MUE=MUJac
      MBJac=MLJac+MUJac+1
      Mdiag=MLE+MUE+1
    endif

! --- BASIC INTEGRATION STEP
    do
      if(Nstep>Nmax) then
        print '(" Exit of RoDas_ODE at X=",E18.4)', X
        print *,  ' More than Nmax =',Nmax,'steps are needed'
        Idid=-2
        return
      endif
      if(0.1D0*abs(H)<=abs(X)*Uround) then
        print '(" Exit of RoDas_ODE at X=",E18.4)', X
        print *,  ' Step size too small, H=', H
        Idid=-3
        return
      endif
      if(Last) then
        H=Hopt
        Idid=1
        return
      endif
      Hopt=H
      if((X+H*1.0001D0-Xend)*PosNeg>=0.D0) then
        H=Xend-X
        Last=.true.
      endif

!*********************************************************************    
!  COMPUTATION OF THE JACOBIAN
      call FcN(N,Y,DY1)
      Nfcn=Nfcn+1
      Njac=Njac+1
      if(IJac==0) then
! --- COMPUTE JACOBIAN MATRIX NUMERICALLY
        if(Banded)  then
! --- JACOBIAN IS Banded
          MUJacP=MUJac+1
          MD=min(MBJac, N)
          do mm=1, m1/m2+1
            do k=1, MD
              j=k+(mm-1)*m2
              do
                Ak2(j)=Y(j)
                Ak3(j)=dsqrt(Uround*max(1.D-5,abs(Y(j))))
                Y(j)=Y(j)+Ak3(j)
                j=j+MD
                if(j>mm*m2) exit
              enddo
              call FcN(N,Y,Ak1)
              j=k+(mm-1)*m2
              j1=k
              Lbeg=max(1,j1-MUJac)+m1
              do
                Lend=min(m2,j1+MLJac)+m1
                Y(j)=Ak2(j)
                MUJacJ=MUJacP-j1-m1
                do L=Lbeg, Lend
                  FJac(L+MUJacJ,j)=(Ak1(L)-DY1(L))/Ak3(j)
                enddo
                j=j+MD
                j1=j1+MD
                Lbeg=Lend+1
                if(j>mm*m2) exit
              enddo
            enddo
          enddo
        else
! --- JACOBIAN IS FULL
          do i=1, N
            Ysafe=Y(i)
            Delt=dsqrt(Uround*max(1.D-5,abs(Ysafe)))
            Y(i)=Ysafe+Delt
            call FcN(N,Y,Ak1)
            do j=m1+1, N
              FJac(j-m1,i)=(Ak1(j)-DY1(j))/Delt
            enddo
            Y(i)=Ysafe
          enddo
        endif
      else
! --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
        call Jac(N,Y,FJac,LDJac)
      endif
!*********************************************************************    
    
!  COMPUTE THE STAGES
    do
      Fac=1.D0/(H*Gamma)
      call DeComr(N,FJac,LDJac,m1,m2,nm1,Fac,E,LDE,Ip,Ier,Ijob)
      if(Ier/=0) then
        Nsing=Nsing+1
        if(Nsing>=5) then
          print '(" Exit of RoDas_ODE at X=",E18.4)', X
          print *,  ' MATRIX IS REPEATEDLY SINGULAR, Ier=', Ier
          Idid=-4
          return
        endif
        H=H*0.5D0
        Reject=.true.
        Last=.false.
        cycle
      endif
    
      Ndec=Ndec+1
! --- PREPARE FOR THE COMPUTATION OF THE 6 STAGES
      HC21=C21/H
      HC31=C31/H
      HC32=C32/H
      HC41=C41/H
      HC42=C42/H
      HC43=C43/H
      HC51=C51/H
      HC52=C52/H
      HC53=C53/H
      HC54=C54/H
      HC61=C61/H
      HC62=C62/H
      HC63=C63/H
      HC64=C64/H
      HC65=C65/H
      
! --- THE STAGES
      call SlvRod(N,FJac,LDJac,MLJac,MUJac,m1,m2,nm1,Fac,E,LDE,Ip,DY1,Ak1,Fx,Ynew,Ijob,.false.)
      do i=1, N
        Ynew(i)=Y(i)+A21*Ak1(i)
      enddo
      call FcN(N,Ynew,DY)
      do i=1, N
        Ynew(i)=HC21*Ak1(i)
      enddo
      call SlvRod(N,FJac,LDJac,MLJac,MUJac,m1,m2,nm1,Fac,E,LDE,Ip,DY,Ak2,Fx,Ynew,Ijob,.true.)
      do i=1, N
        Ynew(i)=Y(i)+A31*Ak1(i)+A32*Ak2(i)
      enddo
      call FcN(N,Ynew,DY)
      do i=1, N
        Ynew(i)=HC31*Ak1(i)+HC32*Ak2(i)
      enddo
      call SlvRod(N,FJac,LDJac,MLJac,MUJac,m1,m2,nm1,Fac,E,LDE,Ip,DY,Ak3,Fx,Ynew,Ijob,.true.)
      do i=1, N
        Ynew(i)=Y(i)+A41*Ak1(i)+A42*Ak2(i)+A43*Ak3(i)
      enddo
      call FcN(N,Ynew,DY)
      do i=1, N
        Ynew(i)=HC41*Ak1(i)+HC42*Ak2(i)+HC43*Ak3(i)
      enddo
      call SlvRod(N,FJac,LDJac,MLJac,MUJac,m1,m2,nm1,Fac,E,LDE,Ip,DY,Ak4,Fx,Ynew,Ijob,.true.)
      do i=1, N
        Ynew(i)=Y(i)+A51*Ak1(i)+A52*Ak2(i)+A53*Ak3(i)+A54*Ak4(i)
      enddo
      call FcN(N,Ynew,DY)
      do i=1, N
        Ak6(i)=HC52*Ak2(i)+HC54*Ak4(i)+HC51*Ak1(i)+HC53*Ak3(i) 
      enddo
      call SlvRod(N,FJac,LDJac,MLJac,MUJac,m1,m2,nm1,Fac,E,LDE,Ip,DY,Ak5,Fx,Ak6,Ijob,.true.)
      
! ------------ EMBEDDED SOLUTION ---------------
      do i=1, Nsk
        Ynew(i)=Ynew(i)+Ak5(i)
      enddo
      call FcN(N,Ynew,DY)
      do i=1, N
        CONT(i)=HC61*Ak1(i)+HC62*Ak2(i)+HC65*Ak5(i)+HC64*Ak4(i)+HC63*Ak3(i) 
      enddo
      call SlvRod(N,FJac,LDJac,MLJac,MUJac,m1,m2,nm1,Fac,E,LDE,Ip,DY,Ak6,Fx,CONT,Ijob,.true.)
     
! ------------ NEW SOLUTION ---------------
      do  i=1, N
        Ynew(i)=Ynew(i)+Ak6(i)
      enddo
      Nsol=Nsol+6
      Nfcn=Nfcn+5 
      
!  ERROR ESTIMATION  
      Nstep=Nstep+1
      Err=0.D0
      do i=1, N
         if(Itol==0)  then
            SK=Atol(1)+Rtol(1)*max(abs(Y(i)),abs(Ynew(i)))
         else
            SK=Atol(i)+Rtol(i)*max(abs(Y(i)),abs(Ynew(i)))
         endif
         Err=Err+(Ak6(i)/SK)**2
      enddo
    
      Err=sqrt(Err/N)
! --- COMPUTATION OF Hnew
! --- WE REQUIRE .2<=Hnew/H<=6.
      Fac=max(Fac2,min(Fac1,(Err)**0.25D0/Safe))
      Hnew=H/Fac

!  IS THE ERROR SMALL ENOUGH ?
      if(Err<=1.D0)  then                      ! --- STEP IS ACCEPTED
        Naccpt=Naccpt+1
        if(Pred)  then
!--- PREDICTIVE CONTROLLER OF GUSTAFSSON
          if(Naccpt>1)  then
            FacGus=(Hacc/H)*(Err**2/ErrAcc)**0.25D0/Safe
            FacGus=max(Fac2,min(Fac1,FacGus))
            Fac=max(Fac,FacGus)
            Hnew=H/Fac
          endif
          Hacc=H
          ErrAcc=max(1.0D-2,Err)
        endif
        do i=1, N 
          Y(i)=Ynew(i)
        enddo
        Xold=X 
        X=X+H
        if(abs(Hnew)>HmaxN) Hnew=PosNeg*HmaxN
        if(Reject) Hnew=PosNeg*min(abs(Hnew),abs(H))
        Reject=.false.
        H=Hnew
        exit
      else                                     ! --- STEP IS REJECTED 
        Reject=.true.
        Last=.false.
        H=Hnew
        if(Naccpt>=1) Nrejct=Nrejct+1
        cycle
      endif

    enddo

    enddo
     
  end subroutine RosCore

!......................................................................................................      
  subroutine RoCoef(Meth,A21,A31,A32,A41,A42,A43,A51,A52,A53,A54,C21,C31,C32,C41,C42,C43,C51,C52,C53,&
            &C54,C61,C62,C63,C64,C65,Gamma,C2,C3,C4,D1,D2,D3,D4,D21,D22,D23,D24,D25,D31,D32,D33,D34,D35)
     
    implicit none

    integer Meth
    real(8) A21, A31, A32, A41, A42, A43, A51, A52, A53, A54, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, C61
    real(8) C62, C63, C64, C65, Gamma, C2, C3, C4, D1, D2, D3, D4, D21, D22, D23, D24, D25, D31, D32, D33, D34, D35
    real(8) BET2P, BET3P, BET4P

      select case(Meth)
        case(1)
          C2=0.386D0
          C3=0.21D0 
          C4=0.63D0
          BET2P=0.0317D0
          BET3P=0.0635D0
          BET4P=0.3438D0 
          D1= 0.2500000000000000D+00
          D2=-0.1043000000000000D+00
          D3= 0.1035000000000000D+00
          D4=-0.3620000000000023D-01
          A21= 0.1544000000000000D+01
          A31= 0.9466785280815826D+00
          A32= 0.2557011698983284D+00
          A41= 0.3314825187068521D+01
          A42= 0.2896124015972201D+01
          A43= 0.9986419139977817D+00
          A51= 0.1221224509226641D+01
          A52= 0.6019134481288629D+01
          A53= 0.1253708332932087D+02
          A54=-0.6878860361058950D+00
          C21=-0.5668800000000000D+01
          C31=-0.2430093356833875D+01
          C32=-0.2063599157091915D+00
          C41=-0.1073529058151375D+00
          C42=-0.9594562251023355D+01
          C43=-0.2047028614809616D+02
          C51= 0.7496443313967647D+01
          C52=-0.1024680431464352D+02
          C53=-0.3399990352819905D+02
          C54= 0.1170890893206160D+02
          C61= 0.8083246795921522D+01
          C62=-0.7981132988064893D+01
          C63=-0.3152159432874371D+02
          C64= 0.1631930543123136D+02
          C65=-0.6058818238834054D+01
          Gamma= 0.2500000000000000D+00  
          D21= 0.1012623508344586D+02
          D22=-0.7487995877610167D+01
          D23=-0.3480091861555747D+02
          D24=-0.7992771707568823D+01
          D25= 0.1025137723295662D+01
          D31=-0.6762803392801253D+00
          D32= 0.6087714651680015D+01
          D33= 0.1643084320892478D+02
          D34= 0.2476722511418386D+02
          D35=-0.6594389125716872D+01

        case(2)   
          C2=0.3507221D0
          C3=0.2557041D0 
          C4=0.6817790D0
          BET2P=0.0317D0
          BET3P=0.0047369D0
          BET4P=0.3438D0 
          D1= 0.2500000000000000D+00
          D2=-0.6902209999999998D-01
          D3=-0.9671999999999459D-03
          D4=-0.8797900000000025D-01
          A21= 0.1402888400000000D+01
          A31= 0.6581212688557198D+00
          A32=-0.1320936088384301D+01
          A41= 0.7131197445744498D+01
          A42= 0.1602964143958207D+02
          A43=-0.5561572550509766D+01
          A51= 0.2273885722420363D+02
          A52= 0.6738147284535289D+02
          A53=-0.3121877493038560D+02
          A54= 0.7285641833203814D+00
          C21=-0.5104353600000000D+01
          C31=-0.2899967805418783D+01
          C32= 0.4040399359702244D+01
          C41=-0.3264449927841361D+02
          C42=-0.9935311008728094D+02
          C43= 0.4999119122405989D+02
          C51=-0.7646023087151691D+02
          C52=-0.2785942120829058D+03
          C53= 0.1539294840910643D+03
          C54= 0.1097101866258358D+02
          C61=-0.7629701586804983D+02
          C62=-0.2942795630511232D+03
          C63= 0.1620029695867566D+03
          C64= 0.2365166903095270D+02
          C65=-0.7652977706771382D+01
          Gamma= 0.2500000000000000D+00  
          D21=-0.3871940424117216D+02
          D22=-0.1358025833007622D+03
          D23= 0.6451068857505875D+02
          D24=-0.4192663174613162D+01
          D25=-0.2531932050335060D+01
          D31=-0.1499268484949843D+02
          D32=-0.7630242396627033D+02
          D33= 0.5865928432851416D+02
          D34= 0.1661359034616402D+02
          D35=-0.6758691794084156D+00

! Coefficients for RoDas_ODE with order 4 for linear parabolic problems
! Gerd Steinebach (1993)
        case(3)   
          Gamma=0.25D0
          C2=3.d0*Gamma
          C3=0.21D0 
          C4=0.63D0
          BET2P=0.D0
          BET3P=c3*c3*(c3/6.d0-Gamma/2.d0)/(Gamma*Gamma)
          BET4P=0.3438D0 
          D1= 0.2500000000000000D+00
          D2=-0.5000000000000000D+00
          D3=-0.2350400000000000D-01
          D4=-0.3620000000000000D-01
          A21= 0.3000000000000000D+01
          A31= 0.1831036793486759D+01
          A32= 0.4955183967433795D+00
          A41= 0.2304376582692669D+01
          A42=-0.5249275245743001D-01
          A43=-0.1176798761832782D+01
          A51=-0.7170454962423024D+01
          A52=-0.4741636671481785D+01
          A53=-0.1631002631330971D+02
          A54=-0.1062004044111401D+01
          C21=-0.1200000000000000D+02
          C31=-0.8791795173947035D+01
          C32=-0.2207865586973518D+01
          C41= 0.1081793056857153D+02
          C42= 0.6780270611428266D+01
          C43= 0.1953485944642410D+02
          C51= 0.3419095006749676D+02
          C52= 0.1549671153725963D+02
          C53= 0.5474760875964130D+02
          C54= 0.1416005392148534D+02
          C61= 0.3462605830930532D+02
          C62= 0.1530084976114473D+02
          C63= 0.5699955578662667D+02
          C64= 0.1840807009793095D+02
          C65=-0.5714285714285717D+01
          D21= 0.2509876703708589D+02
          D22= 0.1162013104361867D+02
          D23= 0.2849148307714626D+02
          D24=-0.5664021568594133D+01
          D25= 0.0000000000000000D+00
          D31= 0.1638054557396973D+01
          D32=-0.7373619806678748D+00
          D33= 0.8477918219238990D+01
          D34= 0.1599253148779520D+02
          D35=-0.1882352941176471D+01
      endselect
  end subroutine RoCoef

!......................................................................................................      
  subroutine DeComr(N,FJac,LDJac,m1,m2,nm1,Fac1,E1,LDE1,IP1,Ier,Ijob)
     
    implicit none
      
    integer N, i, j, k, IP1(nm1), Ier, LDJac, m1, m2, mm, nm1, LDE1, IP1, Ijob, Jm1
    real(8) FJac(LDJac,N), E1(LDE1,nm1), Fac1, Sum
     
    selectcase(Ijob)
      case(1)           ! ---  B=IDENTITY, JACOBIAN A FULL MATRIX
        do j=1, N
          do i=1, N
            E1(i,j)=-FJac(i,j)
          enddo
          E1(j,j)=E1(j,j)+Fac1
        enddo
        call Dec(N,LDE1,E1,IP1,Ier)
        return
          
      case(11)           ! ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
        do j=1,nm1
          Jm1=j+m1
          do i=1,nm1
            E1(i,j)=-FJac(i,Jm1)
          enddo
          E1(j,j)=E1(j,j)+Fac1
        enddo
        mm=m1/m2
        do j=1,m2
          do i=1,nm1
            Sum=0.D0
            do k=0,mm-1
              Sum=(Sum+FJac(i,j+k*m2))/Fac1
            enddo
            E1(i,j)=E1(i,j)-Sum
          enddo
        enddo
        call Dec(nm1,LDE1,E1,IP1,Ier)
        return

      case(2)           ! ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
        do j=1, N
          do i=1, MBJac
            E1(i+MLE,j)=-FJac(i,j)
          enddo
          E1(Mdiag,j)=E1(Mdiag,j)+Fac1
        enddo
!        call DECB (N,LDE1,E1,MLE,MUE,IP1,Ier)
        return

      case(12)           ! ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
        do j=1, nm1
          Jm1=j+m1
          do i=1, MBJac
            E1(i+MLE,j)=-FJac(i,Jm1)
          enddo
          E1(Mdiag,j)=E1(Mdiag,j)+Fac1
        enddo
        mm=m1/m2
        do j=1,m2
          do i=1, MBJac
            Sum=0.D0
            do k=0, mm-1
              Sum=(Sum+FJac(i,j+k*m2))/Fac1
            enddo
            E1(i+MLE,j)=E1(i+MLE,j)-Sum
          enddo
        enddo
!        call DECB (nm1,LDE1,E1,MLE,MUE,IP1,Ier)
        return
        
    endselect

  end subroutine DeComr
  
!......................................................................................................      
  subroutine SlvRod(N,FJac,LDJac,MLJac,MUJac,m1,m2,nm1,Fac1,E,LDE,Ip,DY,Ak,Fx,YNew,Ijob,Stage1)
     
    implicit none
      
    integer N, j, k, Ip(nm1), LDJac, MUJac, MLJac, m1, m2, nm1, mm, LDE, Ip, Ijob, Jkm, Im1
    real(8) FJac(LDJac,N), E(LDE,nm1), DY(N), Ak(N), Fx(N), YNew(N), Fac1, Sum
    logical Stage1

    do i=1,N
      Ak(i)=DY(i)
    enddo

    selectcase(Ijob)
      case(1)           ! ---  B=IDENTITY, JACOBIAN A FULL MATRIX
        if(Stage1) then
          do i=1,N
            Ak(i)=Ak(i)+YNew(i)
          enddo
        endif
        call Sol (N,LDE,E,Ak,Ip)
        return

      case(11)           ! ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
        if(Stage1) then
          do i=1,N
            Ak(i)=Ak(i)+YNew(i)
          enddo
        endif
        mm=m1/m2
        do j=1,m2
          Sum=0.D0
          do k=mm-1,0,-1
            Jkm=j+k*m2
            Sum=(Ak(Jkm)+Sum)/Fac1
            do i=1,nm1
              Im1=i+m1
              Ak(Im1)=Ak(Im1)+FJac(i,Jkm)*Sum
            enddo
          enddo
        enddo
        call Sol (nm1,LDE,E,Ak(m1+1),Ip)
        do i=m1,1,-1
          Ak(i)=(Ak(i)+Ak(m2+i))/Fac1
        enddo
        return

      case(2)           ! ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
        if(Stage1) then
          do i=1,N
            Ak(i)=Ak(i)+YNew(i)
          enddo
        endif
!        call SOLB (N,LDE,E,MLE,MUE,Ak,Ip)
        return

      case(12)           ! ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
        if(Stage1) then
          do i=1,N
            Ak(i)=Ak(i)+YNew(i)
          enddo
        endif
        mm=m1/m2
        do j=1,m2
          Sum=0.D0
          do k=mm-1,0,-1
            Jkm=j+k*m2
            Sum=(Ak(Jkm)+Sum)/Fac1
            do i=max(1,j-MUJac),min(nm1,j+MLJac)
              Im1=i+m1
              Ak(Im1)=Ak(Im1)+FJac(i+MUJac+1-j,Jkm)*Sum
            enddo
          enddo
        enddo
!        call SOLB (nm1,LDE,E,MLE,MUE,Ak(m1+1),Ip)
        do i=m1,1,-1
          Ak(i)=(Ak(i)+Ak(m2+i))/Fac1
        enddo
        return
        
    endselect

  end subroutine SlvRod
      
!-----------------------------------------------------------------------
!  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
!  INPUT..
!     N = ORDER OF MATRIX.
!     NDIM = DECLARED DIMENSION OF ARRAY  A .
!     A = MATRIX TO BE TRIANGULARIZED.
!  OUTPUT..
!     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
!     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
!     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
!     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
!     Ier = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
!           SINGULAR at STAGE K.
!  USE  Sol  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
!  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
!  IF IP(N)=O, A IS SINGULAR, Sol WILL DIVIDE BY ZERO.
!
!  REFERENCE..
!     !. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!     !.A.!.M. 15 (1972), P. 274.
!-----------------------------------------------------------------------
  subroutine Dec(N,Ndim,A,Ip,Ier)
          
    implicit none
      
    integer N, Ndim, Ip(N), Ier, nm1, k, Kp1, m, i, j
    real(8) A(Ndim,N), T

    Ier=0
    Ip(N)=1
    
    if(N>1) then
      nm1=N-1
      do k=1, nm1
        Kp1=k+1
        m=k
        do i=Kp1, N
          if(dabs(A(i,k))>dabs(A(m,k))) m=i  
        enddo
        Ip(k)=m
        T=A(m,k)
        if(m/=k) then
          Ip(N)=-Ip(N)
          A(m,k)=A(k,k)
          A(k,k)=T
        endif
!        if(dabs(T)<1.d-20) then
        if(T==0.d0) then
          Ier=k
          Ip(N)=0
          return
        endif
        T=1.D0/T
        do i=Kp1, N
          A(i,k)=-A(i,k)*T
        enddo
        do j=Kp1, N
          T=A(m,j)
          A(m,j)=A(k,j)
          A(k,j)=T
!          if(dabs(T)>1.d-20) then
          if(T/=0.d0) then
            do i=Kp1, N
              A(i,j)=A(i,j)+A(i,k)*T
            enddo
          endif
        enddo
      enddo
    endif
!    if(dabs(A(N,N))<1.D-12) then
    if(A(N,N)==0.d0) then
      Ier=N
      Ip(N)=0
    endif
    return
      
  end subroutine Dec

!-----------------------------------------------------------------------
!  SOLUTION OF LINEAR SYSTEM, A*X = B .
!  INPUT..
!    N = ORDER OF MATRIX.
!    NDIM = DECLARED DIMENSION OF ARRAY  A .
!    A = TRIANGULARIZED MATRIX OBTAINED FROM Dec.
!    B = RIGHT HAND SIDE VECTOR.
!    IP = PIVOT VECTOR OBTAINED FROM Dec.
!  DO NOT USE IF Dec HAS SET Ier .NE. 0.
!  OUTPUT..
!    B = SOLUTION VECTOR, X .
!-----------------------------------------------------------------------
  subroutine Sol(N,Ndim,A,B,Ip)

    implicit none
      
    integer N, Ndim, Ip(N), nm1, k, Kp1, m, i, Kb, Km1, Kb
    real(8) A(Ndim,N), B(N), T

    if(N>1) then
      nm1=N-1
      do k=1, nm1
        Kp1=k+1
        m=Ip(k)
        T=B(m)
        B(m)=B(k)
        B(k)=T
        do i=Kp1, N
          B(i)=B(i)+A(i,k)*T
        enddo
      enddo
      do Kb=1, nm1
        Km1=N-Kb
        k=Km1+1
!        if(abs(A(k,k))>1.d-20) then
        if(A(k,k)/=0.d0) then
          B(k)=B(k)/A(k,k)
        else  
          print *, 'STOP: Zero diagonal element'
          stop  
        endif    
        T=-B(k)
        do i=1, Km1
          B(i)=B(i)+A(i,k)*T
        enddo
      enddo
    endif
    B(1)=B(1)/A(1,1)
    return
      
  end subroutine Sol
!-----------------------------------------------------------------------

end subroutine RoDas_ODE

#endif

end module Atm_Hg_Chemistry

