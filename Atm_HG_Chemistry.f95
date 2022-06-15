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
          Rcb=gamma*Kcb*c_SO2*c_SO2*10.**(4.*pH)
          Rca=0.

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
!
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
            if(k>6) cycle
            Brfield(i,j,k,1:NumPer,nDay)=qBr(1:NumPer)
          enddo
        enddo
      enddo
      where(Brfield<=0.) Brfield=Zero
      close(10)

      do j=Jmin, Jmax                                    
        do i=Imin, Imax                                  
          if(sum(Brfield(i,j,1,1:NumPer,nDay))/NumPer<4.e-4) then  
            Brfield(i,j,1:6,1:NumPer,nDay)=0.                  
          endif                                        
        enddo                                             
      enddo                                               

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
            if(k>6) cycle                                       
            Brfield(i,j,k,1:NumPer,nDay)=qBr(1:NumPer)
          enddo
        enddo
      enddo
      where(Brfield<=0.) Brfield=Zero
      close(10)

      do j=Jmin, Jmax                                     
        do i=Imin, Imax                                   
          if(sum(Brfield(i,j,1,1:NumPer,nDay))/NumPer<4.e-4) then  
            Brfield(i,j,1:6,1:NumPer,nDay)=0.                
          endif                                           
        enddo                                             
      enddo                                              
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


#endif

end module Atm_Hg_Chemistry

