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
module Soil_POP_General

  use GeneralParams
  use Geometry
  use Atm_POP_Params
  use Atm_POP_General
#ifdef M_SOIL
  use Soil_Params
  use Soil_POP_Params
#endif
  use Exch_Params
#ifdef M_VEG
  use Veg_Params
#endif

  implicit none

  character(1024), private :: fileName, fullName

  contains


!========================================================================
! Soil related POP processes
!========================================================================
  subroutine Soil_POPProcess(NSubs, dTime)

  integer NSubs, i, j, k, f, L
  real(8) CSoil(NLSMax), CSoil1(NLSMax), CTot
  real dTime
  integer Nair,iair
  real dTair

#ifdef DEBUG_MODE
    print *, '+Entering Soil_Process... '
#endif

    call Soil_POPpartitioning(NSubs)
    do L = 1, Soil_NumType
      do j = JMin, JMax
        do i = MinI(j), MaxI(j)
            CSoil(1:Soil_KMax) = Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,3))
            CSoil1(1:Soil_KMax) = Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,4))
            call Soil_DeepAccExch(CSoil,CSoil1,dTime)
            Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,3)) = CSoil(1:Soil_KMax)
            Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,4)) = CSoil1(1:Soil_KMax)

            CSoil(1:Soil_KMax)= Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,1)) * dble((Poros - wl)) + &
                & Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,2)) * dble(wl) + &
                & Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,3))
            DiffSoil(1:Soil_KMax) = Dair(gSubsGroupInd(NSubs)) / PartG(i, j, 1:Soil_KMax, Period, gSubsGroupInd(NSubs)) + &
                & Dwater(gSubsGroupInd(NSubs)) / PartL(i, j, 1:Soil_KMax, Period, gSubsGroupInd(NSubs)) + Difs
            call Soil_AdvDif(CSoil, Ve(i,j), DiffSoil, dTime)
            Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,1)) = &
                &CSoil(1:Soil_KMax) / dble(PartG(i, j, 1:Soil_KMax, Period, gSubsGroupInd(NSubs)))
            Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,2)) = &
                &CSoil(1:Soil_KMax) / dble(PartL(i, j, 1:Soil_KMax, Period, gSubsGroupInd(NSubs)))
            Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,3)) = &
                &CSoil(1:Soil_KMax) - Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,1)) * dble((Poros - wl)) - &
                & Soil_Conc(i,j,1:Soil_KMax,L,sFormInd(Soil,NSubs,2)) * dble(wl)
        end do
      end do
    end do
    call Soil_POPdegradation(NSubs, dTime)
    call Soil_POPpartitioning(NSubs)

! Counters for averaging
    do L = 1, Soil_NumType
	  do k = 1, Soil_KMax
        do j = JMin, JMax
          do i = MinI(j), MaxI(j)
             CTot = Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,1))*(Poros-wl)&
                    &+Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,2))*wl &
                    &+Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,3)) &
                    &+Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,4))
		    Soil_ConcDay(i,j,k,L,NSubs) = Soil_ConcDay(i,j,k,L,NSubs) + CTot*dTime
	    end do
      end do
    end do
    end do


#ifdef DEBUG_MODE
    print *, '+Exit Soil_Process'
#endif
end subroutine Soil_POPProcess

!========================================================================
! Advection/diffusion in soil
!========================================================================
subroutine Soil_AdvDif(CSoil, Vel, DCoef, dTime)

  integer k
  real(8) CSoil(NLSMax), Flux(NLSMax)
  real(8) Vel, DCoef(NLSMax)
  real    dTime
  real(8) CS,CS1

  Flux=0.

  do k = 2, Soil_KMax
    Flux(k) = dble(dTime * 2.) * (CSoil(k - 1) - CSoil(k)) * dble(DCoef (k)/ (Soil_dz(k - 1) + Soil_dz(k)))
	if (Vel > 0.) Flux(k) = Flux(k) + dble(Vel) * CSoil(k - 1)
	if (Vel < 0.) Flux(k) = Flux(k) + dble(Vel) * CSoil(k)
  end do
  do k = 1, Soil_KMax - 1
    CSoil(k) = CSoil(k) + (Flux(k) - Flux(k + 1)) / dble(Soil_dz(k))
  end do
  CSoil(Soil_KMax) = CSoil(Soil_KMax) + Flux(Soil_KMax) / dble(Soil_dz(Soil_KMax))

end subroutine Soil_AdvDif

subroutine Soil_DeepAccExch(CSoil,CSoil1,dTime)

  integer k
  real(8) CSoil(NLSMax), CSoil1(NLSMax), Flux
  real dTime

  do k = 1, Soil_KMax
    Flux = dTime * RateConst * (Facc * CSoil1(k) - (1 - Facc) * CSoil(k))
	CSoil1(k) = CSoil1(k) - Flux
	CSoil(k) = CSoil(k) + Flux
  end do

end subroutine Soil_DeepAccExch


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine defining physical and chemical properties of a POP in soil
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Soil_POP_Props(Subs)

	integer i, FileStat, lenType, Subs, Form, Limit, Ind
	character(80) strRead, strPar
	character(30) strType, FormType(MaxForm)
	integer mon, ios

    fileName='Soil'//trim(PropName)//trim(SubsID(Subs))//'.dat'
	fullName=trim(PropPath)//fileName
	open(2, file=trim(fullName), status='old', iostat=FileStat, action='read')
	if(FileStat>0) then
	  print '(/,"STOP: Cannot open file of soil POP properties ''",a,"''",/)', trim(fullName)
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
		do i=1, sFormNum(Soil,Subs)
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
!------------------------------------------------------------------ Added 09.04.2012 A.G.
                    case('Anthrop emiss')
		    if(strPar=='yearly') then
                        SoilEmissionStep='yearly'
			SoilEmisNum=SoilEmisNum+1
                        ReadSoilInd(SoilEmisNum,yr)=sFormInd(Soil,Subs,Form)
			SoilEmisInd(SoilEmisNum)=sFormInd(Soil,Subs,Form)
		    elseif(strPar=='monthly') then
                        SoilEmissionStep='monthly'
			SoilEmisNum=SoilEmisNum+1
                        ReadSoilInd(SoilEmisNum,mn)=sFormInd(Soil,Subs,Form)
			SoilEmisInd(SoilEmisNum)=sFormInd(Soil,Subs,Form)
		    elseif(strPar=='daily') then
                        SoilEmissionStep='daily'
			SoilEmisNum=SoilEmisNum+1
                        ReadSoilInd(SoilEmisNum,da)=sFormInd(Soil,Subs,Form)
			SoilEmisInd(SoilEmisNum)=sFormInd(Soil,Subs,Form)
                    else
                    endif
!------------------------------------------------------------------
		  case('Forms number')
			read(strPar,'(i2)') sFormNum(Soil,Subs)
                        FormSubs(Soil,NumForm(Soil)+1:NumForm(Soil)+sFormNum(Soil,Subs))=Subs
			if (sFormNum(Soil,Subs).gt.0) then	
			NumSubsMedia(Soil)=NumSubsMedia(Soil)+1
			gSubsMediaInd(Soil,Subs)=NumSubsMedia(Soil)
			end if
#ifdef DEBUG_MODE
			print *, 'Forms: ',sFormNum(Soil,Subs)
#endif
!------------------------------------------------------------------
		  case('Forms')
			read(strPar,*) (FormType(i), i=1, sFormNum(Soil,Subs))
#ifdef DEBUG_MODE
			print *, 'Soil forms: ',FormType
#endif
!------------------------------------------------------------------
		  case('Form ID')
		    FormID(Soil,Subs,Form)=strPar
                    lenType=scan(strPar,')')
                        sFormInd(Soil,Subs,Form)=NumForm(Soil)+1
                        FormSubs(Soil,sFormInd(Soil,Subs,Form))=Subs
                        NumForm(Soil)=NumForm(Soil)+1
#ifdef DEBUG_MODE
                        print *, 'Form: ',Form
                        print *, 'sFormInd(Soil,Subs,Form): ', sFormInd(Soil,Subs,Form)
                        print *, 'FormSubs(Soil,sFormInd(Soil,Subs,Form)): ', FormSubs(Soil,sFormInd(Soil,Subs,Form))
                        print *, 'NumForm(Soil): ',NumForm(Soil)
                        print *, ' '
#endif
!------------------------------------------------------------------
		  case('SoilDegr')
             read(strPar,*) CDSoilTop(gSubsGroupInd(Subs)),CDSoilBottom(gSubsGroupInd(Subs)),&
                                   &CDSoilBio(gSubsGroupInd(Subs))
!------------------------------------------------------------------
		  case('Kow')
              read(strPar,*)Kow0(gSubsGroupInd(Subs)),KowT(gSubsGroupInd(Subs))
!------------------------------------------------------------------
	endselect

	  endif
	enddo

	close(2)

end subroutine Soil_POP_Props

!========================================================================
! Calculation of time step for soil
!========================================================================
subroutine Soil_POPTStepCalc(NSubs,dTiter)

  integer i, j, k, NSubs
  real TmpVal,DCoeff
  real dTiter

    do j = JMin, JMax
        do i = MinI(j), MaxI(j)
            do k = 1, Soil_KMax
                if (k == 1) then
                    DCoeff = Dair(gSubsGroupInd(NSubs)) / PartG(i, j, k, Period, gSubsGroupInd(NSubs)) + &
		           & Dwater(gSubsGroupInd(NSubs)) / PartL(i, j, k, Period, gSubsGroupInd(NSubs)) + Difs
                    TmpVal = abs(Ve(i,j)) + 2.* DCoeff / (Soil_dz(k) + Soil_dz(k + 1))
                elseif(k < Soil_KMax) then
		    TmpVal = abs(Ve(i,j)) + 2.* DCoeff * (1. / (Soil_dz(k - 1) + Soil_dz(k)) + 1. &
			                                           & / (Soil_dz(k) + Soil_dz(k + 1)))
                else
		    TmpVal = abs(Ve(i,j)) + 2.* DCoeff / (Soil_dz(k - 1) + Soil_dz(k))
                end if
                dTiter = min(0.5 * Soil_dz(k)/TmpVal, dTiter)
            end do
        end do
    end do

end subroutine Soil_POPTStepCalc

!========================================================================
! Partitioning between forms in soil
!========================================================================
subroutine Soil_POPpartitioning(NSubs)

  integer NSubs
  integer i, j, k, L
  real(8) CTot

    do L = 1, Soil_NumType
        do k = 1, Soil_KMax
            do j = JMin, JMax
                do i = MinI(j), MaxI(j)
                    CTot = Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,1)) * dble((Poros - wl)) + &
		            & Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,2)) * dble(wl) + &
                            & Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,3))
                    Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,1)) = CTot / dble(PartG(i, j, k, Period, gSubsGroupInd(NSubs)))
		    Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,2)) = CTot / dble(PartL(i, j, k, Period, gSubsGroupInd(NSubs)))
		    Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,3)) = CTot - &
                            & Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,1)) * dble((Poros - wl)) - &
		            & Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,2)) * dble(wl)
	    end do
	  end do
    end do
  end do

end subroutine Soil_POPpartitioning

!========================================================================
! Degradation of a POP in soil
!========================================================================
subroutine Soil_POPdegradation(NSubs, dTime)

  integer NSubs
  integer i, j, k, f, L
  real DTime
  real*8 SPM1,SPM2

    SPM1 = Soil_POPMass(NSubs)

    do L = 1, Soil_NumType
        do k = 1, Soil_KMax
            do j = JMin, JMax
                do i = MinI(j), MaxI(j)
                    do f = 1, 3
                        Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,f)) = &
			      & Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,f)) * dble((1. - CDSoil(k,gSubsGroupInd(NSubs)) * dTime))
                    end do
		    Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,4)) = &
		               & Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,4)) * dble((1. - CDSoil1 * dTime))
                end do
            end do
        end do
    end do

    SPM2 = Soil_POPMass(NSubs)
    Soil_Degr = Soil_Degr + SPM1 - SPM2

end subroutine Soil_POPdegradation

!========================================================================
! Mass in soil
!========================================================================
real(8) function Soil_POPMass(NSubs)

  integer i, j, k, L, NSubs

    Soil_POPMass = 0.
    do L = 1, Soil_NumType
      do k = 1, Soil_KMax
        do j =  JMin, JMax
          do i = MinI(j), MaxI(j)
 	        Soil_POPMass = Soil_POPMass + (Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,1))* dble(Poros - wl) &
                            & + Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,2)) * dble(wl) &
                            & + Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,3)) &
                            & + Soil_Conc(i,j,k,L,sFormInd(Soil,NSubs,4))) &
                            & * dble(MeshArea(i,j)) * dble(Soil_dz(k)) * Soil_Frac(i,j,L)
          end do
        end do
      end do
    end do

end function Soil_POPMass

!========================================================================
! Exchange with soil
!========================================================================
subroutine Soil_POPExchPar(NSubs)

  integer i, j, NSubs, L,n,s
  real Deff, RS, SVE, pg, pl
  real PartConst        ! Partitioning constant between soil (total conc) and air (Cair/Csoil)
  real RSoil            ! Resistance of the first soil layer
  real Ra,Rb,RT,Pxx,Zr,Ux,ZoM,Ux0,Lmo0,Ui,Uj,Uref,Almo,disp,U10
#ifdef M_VEG
  real Kav01,Kav02,BCF,Koa,ResA,RTot,VdGVeg,Tair,Kow,Kcw,Pcw,Henry,Kaw,GCut,RCut
  real Ueff, LeafLength, DeltaBL, ResBL
  real KoaExp
#endif


#ifdef DEBUG_MODE
      print *, '-> Entering Soil_POPExchPar...'
#endif
    do j = JMin, JMax
        do i = MinI(j), MaxI(j)
            s = Season(i,j,Month)								! Seasons (1-5)
            pg = PartG(i, j, 1, Period, gSubsGroupInd(NSubs))
            pl = PartL(i, j, 1, Period, gSubsGroupInd(NSubs))
	    Deff = Dair(gSubsGroupInd(NSubs))/pg+Dwater(gSubsGroupInd(NSubs))/pl+Difs
	    if (Ve(i,j) > 0.) then
                RS = 1. / (2. * Deff/Soil_dz(1) + Ve(i, j))
                PartConst = RS * (2. * Deff/Soil_dz(1))
	    else
                RS = 1. / (2. * Deff/Soil_dz(1))
                PartConst = RS * (2. * Deff/Soil_dz(1) - Ve(i, j))
	    end if
	    RSoil = RS / pg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            PartConst = PartConst / PartG(i, j, 1, Period,gSubsGroupInd(NSubs))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            Ux0=Ufric(i,j,Period,toDay)												! Mean friction velocity
            Ui=(Uwind(i,j,1,Period,toDay)+Uwind(i-1,j,1,Period,toDay))/2.
            if(j == 1) then
                Uj=Vwind(i,j,1,Period,toDay)
            else
                Uj=(Vwind(i,j,1,Period,toDay)+Vwind(i,j-1,1,Period,toDay))/2.
            end if
            Uref=sqrt(Ui*Ui+Uj*Uj)						! Absolute value of wind velosity
            Pxx=Px(i,j,Period,toDay)						! Surface pressure Px=Ps-Pt
            RT=RTemp(i,j,1,Period,toDay)					! RT of air at the lowest sigma level
            Zr=RT/Ggrav*log((Pxx+Ptop)/(Sigma(1)*Pxx+Ptop))			! Height of the lowest sigma level

   	    do L = 1, Soil_NumType
                Lmo0=MOLength(i,j,Period,toDay)
                ZoM=Soil_Z0(L,s)						! Roughness
!                if(trim(adjustL(SoilType(L)))  == 'Water') then
!                    ZoM=0.016*Ux0*Ux0/Ggrav+Nu/9.1*Ux0
!                    Ux=Karman*Uref/IphiM(Lmo0,Zr,ZoM)
!                    Ux = max(Ux, 0.1)	! 09.04.2008 Ilia
!                else
                    Ux=Karman*Uref/IphiM(Lmo0,Zr,ZoM)				! Friction velocity
                    Ux = max(Ux, 0.1)	! 09.04.2008 Ilia
                    Almo = Lmo0 / Ux / Ux / Ux
                    disp=Soil_Disp(L,s)
                    do n=1, 5
                        Ux=max(Karman*Uref/IphiM(Lmo0,Zr-disp,ZoM),0.1)
                        Lmo0=Almo*Ux*Ux*Ux
                    enddo
!                end if
            
                Ra = 0.74*(alog(Zr/ZoM)-PsiH(Zr/Lmo0)+PsiH(ZoM/Lmo0))/Karman/Ux
                Rb = 2./Karman/Ux*(1.5e-5/Dair(gSubsGroupInd(NSubs))/0.71)**(2./3.)

!            if(i == 120 .and. j == 50) then
!            write(107, '(9e14.6)') Ra, Rb, Rsoil, Deff, pg, pl, Difs, PartConst, Ve(i,j)
!            end if

#ifdef M_VEG

!
! -------------- Gas exchange with vegetation (old scheme) --------------------------
!
                if (L >= 1 .and. L<= Veg_NumType) then
                    ResA= 0.74*2.5/Ux*(alog(Zr/1.)-PsiH(Zr/Lmo0)+PsiH(Soil_Height(L)/Lmo0))
                    if(ResA.lt.0.) then
                        ResA = 0.
                    end if
                    KoaExp = exp(KoaT(NSubs) * ((1. / TempAir(i,j,1,Period,toDay)) - (1. / T0)))
                    Koa = Koa0(NSubs)*KoaExp

                    if (trim(adjustL(SoilType(L))) == 'Decid') then
                        Kav01=30.
                        BCF=22.91*Koa**0.445
                        RTot = ResA + SLeaf/(Kav01*KoaExp)
                    end if
                    if (trim(adjustL(SoilType(L))) == 'Conif') then
                        Kav02=4.6
                        BCF=38.02*Koa**0.69
                        RTot = ResA + SLeaf/(Kav02*KoaExp)
                    end if
                    if((trim(adjustL(SoilType(L))) == 'Grass') .or. &
                        &(trim(adjustL(SoilType(L))) == 'Scrabs') .or. &
                        &(trim(adjustL(SoilType(L))) == 'Arable')) then
                    BCF=14.12*Koa**0.76
                    RTot = ResA+(Vmol**5.50)*(10.**(-1.43))/Koa
                    end if

                    VdGVeg = 1./RTot
                    Veg_VExch(i, j, L, Atm, 1, NSubs) = VdGVeg * LAIPriv(i,j,L)
                    Veg_VExch(i, j, L, Atm, 2, NSubs) = VdGVeg * SLeaf / BCF
                end if
#endif


#ifdef M_VEG_NEW
!
! -------------- Gas exchange with vegetation (new scheme) --------------------------
!
                if (L >= 1 .and. L<= Veg_NumType) then
                    ResA= 0.74*2.5/Ux*(alog(Zr/1.)-PsiH(Zr/Lmo0)+PsiH(Soil_Height(L)/Lmo0))
                    if(ResA.lt.0.) then
                        ResA = 0.
                    end if
		    TAir = TempAir(i,j,1,Period,toDay)
                    Koa = Koa0(NSubs)*exp(KoaT(NSubs)*((1./TAir)-(1./ T0)))
		    Kow = Kow0(NSubs) * exp(KowT(NSubs) * (1 / TAir - 1 / T0))
		    Kcw = 1.14025*Kow**0.97    ! log Kcw = 0.057 + 0.97*Kow
		    Pcw = 5.49541E-12*Kcw**0.734   ! log Pcw = 0.734 * log Kcw - 11.26
		    Henry = HT0(NSubs) / (RUniv * TAir) * exp(- HT(NSubs) * (1 / TAir - 1 / T0))
		    Kaw = Henry/RUniv/TAir
		    Gcut = Pcw / Kaw
		    Rcut = 1. / GCut

                    if (trim(adjustL(SoilType(L))) == 'Decid') then
                        BCF=14.*Koa**0.76
		        Ueff = Uref * 0.6                        ! Wind speed inside forest
		        LeafLength = 0.06                        ! Mean leaf length, m
		        DeltaBL = 4.* 0.001 * sqrt(LeafLength / Ueff)  ! Thickness of leaf boundary layer
		        ResBL = DeltaBL / DAir(NSubs)      ! Resistance of boundary layer near leafs
                        RTot = ResA + ResBL + Rcut
                    end if
                    if (trim(adjustL(SoilType(L))) == 'Conif') then
                        BCF=38.*Koa**0.69
		        Ueff = Uref * 0.6                        ! Wind speed inside forest
		        LeafLength = 0.001                        ! Mean leaf length, m
		        DeltaBL = 4.* 0.001 * sqrt(LeafLength / Ueff)  ! Thickness of leaf boundary layer
		        ResBL = DeltaBL / DAir(NSubs)      ! Resistance of boundary layer near leafs
                        RTot = ResA + ResBL + Rcut
                    end if
                    if((trim(adjustL(SoilType(L))) == 'Grass') .or. &
                        &(trim(adjustL(SoilType(L))) == 'Scrabs') .or. &
                        &(trim(adjustL(SoilType(L))) == 'Arable')) then
                        BCF=0.000275*Koa**1.15
		        LeafLength = 0.006                        ! Mean leaf length, m
		        DeltaBL = 4.* 0.001 * sqrt(LeafLength / Uref)  ! Thickness of leaf boundary layer
		        ResBL = DeltaBL / DAir(NSubs)      ! Resistance of boundary layer near leafs
                        RTot = ResA + ResBL + Rcut
                    end if

                    VdGVeg = 1./RTot
                    Veg_VExch(i, j, L, Atm, 1, NSubs) = VdGVeg * LAIPriv(i,j,L)
                    Veg_VExch(i, j, L, Atm, 2, NSubs) = VdGVeg * SLeaf / BCF
                end if
#endif
                SVE = 1. / (Ra + Rb + RSoil)
                Soil_VExch(i, j, L, Atm, 1, NSubs) = SVE
                Soil_VExch(i, j, L, Atm, 2, NSubs) = SVE * PartConst
            end do ! L
        end do ! i
    end do ! j

#ifdef DEBUG_MODE
      print *, '-> Exit Soil_POPExchPar'
#endif
end subroutine Soil_POPExchPar


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading soil emission fields
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Soil_POP_ReadEmis(Subs,Ind,tm,Note)

	integer i, j, ln, Ind, FileStat, Hs, Xscal, Form, Subs, lenStr
	integer Src, SrcInd
        integer tm
	real EmisSrc(2), ConvEmisUnit, Aver(Imin:Imax)
	character(30) parType, gridID, units
        character(4)  YearNum
        character(300) readStr, parVal, SrcInfo
        character(10) IDsrc, IDsubs, IDform, TResol,Period
        character(*) Note
        integer ios

        Form=SoilEmisInd(Ind)

    selectcase(tm)                   ! Modified Sep 2016 V. Sh
      case(yr)
        write(fileName,'(i4)') Year
	fileName=trim(SubsID(Subs))//trim(FormID(Soil,Subs,Form))//'_'//trim(fileName)//trim(AntEmisName)//'.dat'
        ConvEmisUnit=1000./sum(MonthDays)/SecInDay	! Convert t/y ==> kg/s
      case(mn)
        write(fileName,'(i4,i2.2)') Year, Month
	fileName=trim(SubsID(Subs))//trim(FormID(Soil,Subs,Form))//'_'//trim(fileName)//trim(AntEmisName)//'.dat'
        ConvEmisUnit=1000./MonthDays(Month)/SecInDay	! Convert t/m ==> kg/s
      case(da)
        write(fileName,'(i4,i2.2,i2.2)') Year, Month, Day
	fileName=trim(SubsID(Subs))//trim(FormID(Soil,Subs,Form))//'_'//trim(fileName)//trim(AntEmisName)//'.dat'
        ConvEmisUnit=1000./SecInDay	! Convert t/day ==> kg/s
    endselect

    write(YearNum,'(i4)') Year
    fullName=trim(EmisPath)//trim(GridCode)//'/'//trim(SubsID(Subs))//trim(AntEmisName)//'/'//YearNum//'/'//trim(fileName)
	open(4, file=fullName, status='old', iostat=FileStat, action='read')
	if(FileStat>0) then
	  print '(/,"STOP: Cannot open file with soil emission ''",a,"''",/)', trim(fullName)
	  stop
	endif

! Reading of header
    do ln=1, 8
      read(4,'(a)') readStr
      lenStr=scan(readStr,achar(13))
      if(lenStr>0) readStr=readStr(:lenStr-1)

      lenStr=scan(readStr,':')
      if(lenStr==0) then
        print '(/,"STOP: Wrong format of file ''",a,"''",/)', trim(fullName)
        stop
      endif
      parType=trim(readStr(:lenStr-1))
      parVal=readStr(lenStr+1:)
      parVal=adjustl(parVal)

      selectcase(parType)
	case('Substance')
          IDsubs=parVal
	case('Form')
          IDform=parVal
	case('Grid Type')
          gridID=parVal
	case('Source')
          SrcInfo=parVal
        case('TimeResol')
          TResol = ParVal
          if (tm == yr .and. trim(TResol) /= 'yearly') then
            print*, 'Invalid period of emission data'
            stop
          end if
          if (tm == mn .and. trim(TResol) /= 'monthly') then
            print*, 'Invalid period of emission data'
            stop
          end if
          if (tm == da .and. trim(TResol) /= 'daily') then
            print*, 'Invalid period of emission data'
            stop
          end if
        case('Period')
          Period = parVal
        case('Units')
          units=parVal
	case('Layers')
          read(parVal,*) HemisAnt
	case default
	  print '(/,"STOP: Unknown input parameter ''",a,"''",/)', trim(parType)
	  stop
      endselect
    enddo

        write(Note,'("Emissions of ",a,a," in ",a,", ",a)') trim(IDsubs), trim(IDform), Period, trim(units)

	read(4,*)
	read(4,*)
	read(4,*)

	SoilEmisFlux(:,:,Ind)=0.
	do
	  read(4,'(a)',iostat=ios) readStr
	  if(ios==-1) exit
            lenStr=scan(readStr,achar(13))
            if(lenStr>0) readStr=readStr(:lenStr-1)
            read(readStr,*) IDsrc, i, j, EmisSrc(1)
            SoilEmisFlux(i,j,Ind) = SoilEmisFlux(i,j,Ind) + EmisSrc(1)   ! t/y, t/m, t/d
	enddo
	close(4)
	SoilEmisFlux(:,:,Ind)=SoilEmisFlux(:,:,Ind)*ConvEmisUnit		! t/y -> kg/sec

! Grid aggregation
	do j=Jmin, Jmax
            if(j==jSP.or.j==jNP) cycle
            Xscal=Imax/maxI(j)
            if(Xscal==1) cycle

            Aver(Imin:Imax)=SoilEmisFlux(Imin:Imax,j,Ind)
            call GridAggreg(j,Xscal,Aver,2)
            SoilEmisFlux(minI(j):maxI(j),j,Ind)=Aver(minI(j):maxI(j))
	enddo

end subroutine Soil_POP_ReadEmis



end module Soil_POP_General
#endif