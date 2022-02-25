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
! Media exchange module
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifdef G_POP
Module Exch_General_POP

use GeneralParams
use Atm_Params
use Exch_Params
use Atm_Output
use GeneralOutput
use Atm_General
use Soil_General
use Ocn_General
use Veg_General

implicit none

  character(800), private :: fileName, fullName
  integer, private      :: CMedMax, iOcn, CMedSoilMax, CMedVegMax, Exch_KMax_2
 
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine managing memory allocation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Exch_Memory(operation)

	character(*) operation

	selectcase(operation)
!---------------------------------------------------------------------------------
	case('Initial')
            call Exch_MemAlloc
            CMedMax = Soil_NumType + 3*Veg_NumType + 1
            CMedSoilMax = Soil_NumType+3*Veg_NumType
            CMedVegMax = Soil_NumType+2*Veg_NumType
            Exch_KMax = Atm_Kmax
            Exch_KMax_2 = Exch_KMax + Exch_KMax
            iOcn = CMedMax
!---------------------------------------------------------------------------------
	case('Final')
	  call Exch_MemDealloc
	endselect

end subroutine Exch_Memory

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine supplying information for the calculation run
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Exch_InputData(Prd)

	character(*) Prd
	character(100) Note
        integer NSubs, i_st, j_st, k

#ifdef DEBUG_MODE
    print *, '-> Entering Exch_InputData'
#endif
	selectcase(Prd)
!---------------------------------------------------------------------------------
	case('Initial')

! Reading dry deposition params (coeffs. for ABL, z0 for seasons and l-uses etc.):
            call DryDepParams
            call LogFile('Reading dry deposition parameters',0)
            CMed = 0.
            Exch_KMed = CMedMax
            DeltaZMed(1:Soil_NumType) = Soil_dz(1)
            DeltaZMed(Soil_NumType+1:CMedSoilMax) = 1.
!---------------------------------------------------------------------------------
	case('6hourly')

  do NSubs = 1, NumSubs
	select case (SubsGroup(NSubs))
	case (POP)                       
            call Atm_CalcPOPWDVelGas(NSubs)
            call Atm_CalcPOPWDCoeff(NSubs)
            call Soil_POPExchPar(NSubs)

#ifdef O_EQUIL
                ! Calculation of exchange between air and water using instanteneous equilibrium
                ! Code included in the Exch_General_POP.f95
#else
            call Ocn_POPExchPar(NSubs)
#endif
            call Atm_CalcPOPDDVelPart(NSubs)
	case default
            continue
	end select
  end do

!---------------------------------------------------------------------------------
	endselect
	
#ifdef DEBUG_MODE
    print *, '-> Exit Exch_InputData'
#endif

end subroutine Exch_InputData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading adjusting parameters for ABL, some common geophys. data etc.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DryDepParams

integer i, FileStat, lenType, Surf, s, SurfInd(MaxSurf), lenStr
integer ios
real zLU(MaxSurf), dLU(MaxSurf), hLU(MaxSurf)
character(300) strRead, strPar
character(40) strType
character(20) IDsurf(MaxSurf), blank

#ifdef DEBUG_MODE
    print *, '-> Entering DryDepParams...'
#endif
	fullName=trim(GeoPath)//trim(ComGeoName)
	open(2, file=trim(fullName), status='old', iostat=FileStat, action='read')
	if(FileStat>0) then
	  print '(/,"STOP: Cannot open file with dry deps parameters ''",a,"''",/)', trim(fullName)
	  stop
	endif

	do
	  read(2,'(a)',iostat=ios) strRead
	  if(ios==-1) exit

! Recognizing comments
	  if(strRead(1:1)=='!'.or.strRead(1:1)==' '.or.strRead(1:1)=='	'.or.ichar(strRead(1:1))==13) cycle

! Deleting tabs and line termination symbols
	  lenType=scan(strRead,'	')
	  do while (lenType>0)
	    strRead=strRead(:lenType-1)//strRead(lenType+1:)
		lenType=scan(strRead,'	')
	  enddo
	  lenType=scan(strRead,achar(13))
          if(lenType>0) strRead=strRead(:lenType-1)

! Reading parameter type
		lenType=scan(strRead,':')
		if(lenType==0) then
	      print '(/,"STOP: Wrong format of file ''",a,"''",/)', trim(fileName)
	      stop
		endif
		strType=trim(strRead(:lenType-1))
	        strPar=strRead(lenType+1:)
	        strPar=adjustl(strPar)

		selectcase(strType)
!------------------------------------------------------------------
			case ('Aforest')
				read(StrPar, *) Aforest
!------------------------------------------------------------------
			case ('Bforest')
				read(StrPar, *) Bforest
!------------------------------------------------------------------
			case ('Cforest')
				read(StrPar, *) Cforest
!------------------------------------------------------------------
			case ('Agrass')
				read(StrPar, *) Agrass
!------------------------------------------------------------------
			case ('Bgrass')
				read(StrPar, *) Bgrass
!------------------------------------------------------------------
			case ('Cgrass')
				read(StrPar, *) Cgrass
!------------------------------------------------------------------
			case ('Dgrass')
				read(StrPar, *) Dgrass
!------------------------------------------------------------------
			case ('Beta for momentum')
				read (StrPar, *) Bm
!------------------------------------------------------------------
			case ('Beta for heat')
				read (StrPar, *) Bh
!------------------------------------------------------------------
			case ('Gamma for momentum')
				read (StrPar, *) Gm
!------------------------------------------------------------------
			case ('Gamma for heat')
				read (StrPar, *) Gh
!------------------------------------------------------------------
			case ('R of broken water surface, s/m')
				read (StrPar, *) Rbrok
!------------------------------------------------------------------
			case ('Fog droplet diameter, m')
				read (StrPar, *) Dfog
!------------------------------------------------------------------
		end select
	enddo
	close(2)

! Reading roughness, displacement and top canopy heights
	fullName=trim(GeoPath)//trim(RoughName)
	open(3, file=trim(fullName), status='old',iostat=FileStat, action='read')
	if(FileStat>0) then
	  print '(/,"STOP: Cannot open file with roughness ''",a,"''",/)', trim(fullName)
	  stop
	endif

	read(3,'(a)') strRead
        lenStr=scan(strRead,achar(13))
        if(lenStr>0) strRead=strRead(:lenStr-1)
        read(strRead,*) blank, (IDsurf(Surf), Surf=1, NumSurf)

        do Surf=1, NumSurf
	  SurfInd(Surf)=0
	  do i=1, NumSurf
		if(SurfType(Surf)==IDsurf(i)) then
		  SurfInd(Surf)=i
		  exit
		endif
	  enddo
	  if(SurfInd(Surf)==0) then
		print '(/,"STOP: No roughness information for surface type ''",a,"''",/)', trim(SurfType(Surf))
		stop
	  endif
	enddo

	read(3,*)
	do i=1, NumSns
	  read(3,'(a)') strRead
           lenStr=scan(strRead,achar(13))
           if(lenStr>0) strRead=strRead(:lenStr-1)
           read(strRead,*) s, (zLU(Surf), Surf=1, NumSurf)
           ZoLU(1:NumSurf,s)=zLU(SurfInd(1:NumSurf))
	enddo

	read(3,*)
	do i=1, NumSns
	  read(3,'(a)') strRead
          lenStr=scan(strRead,achar(13))
          if(lenStr>0) strRead=strRead(:lenStr-1)
          read(strRead,*) s, (dLU(Surf), Surf=1, NumSurf)
          dispLU(1:NumSurf,s)=dLU(SurfInd(1:NumSurf))
	enddo

	read(3,*)
	read(3,'(a)') strRead
        lenStr=scan(strRead,achar(13))
        if(lenStr>0) strRead=strRead(:lenStr-1)
        read(strRead,*) (hLU(Surf), Surf=1, NumSurf)
        heightLU(1:NumSurf)=hLU(SurfInd(1:NumSurf))
	close(3)

#ifdef DEBUG_MODE
    print *, '-> E[it DryDepParams'
#endif
end subroutine DryDepParams

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Exchange between media - main procedure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Media_Exchange_POP

  integer NSubs

#ifdef DEBUG_MODE
    print *, '+Entering Media_Exchange_POP...'
#endif

!  Exch_KMax = Atm_Kmax
  Exch_dT = minval(Tstep)

  do NSubs = 1, NumSubs
	select case (SubsGroup(NSubs))
	  case (POP)                          
	    call Exch_POP(NSubs)
	  case default
	    continue
	end select
  end do
  
#ifdef DEBUG_MODE
    print *, '+Exit Media_Exchange_POP'
#endif
end subroutine Media_Exchange_POP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Integration time step
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TimeStep

  real(8) TmpCoeff, TmpCoeff2
  integer i, j

#ifdef DEBUG_MODE
    print *, '-> TimeStep'
#endif

    DeltaT = Exch_dT
    do i = 2, Exch_KMax         ! Loop through Atm levels: 2 - 20
        TmpCoeff2 = WetDepVel(i) / DeltaZ(i) + WetDepCoeff(i)
        if(TmpCoeff2.gt.0) DeltaT = min(DeltaT, 0.5 / TmpCoeff2)
    end do

    TmpCoeff = 0.
    do i = 1, Exch_KMed
        TmpCoeff = TmpCoeff + (ExchVel(0,i)+DryDepVel(i)) * LC(i)
    end do

    TmpCoeff2 = (WetDepVel(1) + TmpCoeff) / DeltaZ(1) + WetDepCoeff(1)
    if (TmpCoeff2.gt.0) then 
       DeltaT = min(DeltaT, 0.5 / TmpCoeff2)
    end if

    do i = 1, Exch_KMed
        TmpCoeff = 0.
	do j = 0, Exch_KMed
            TmpCoeff = TmpCoeff + ExchVel(i,j)
	end do
        if(TmpCoeff.gt.0) then	
            DeltaT = min(DeltaT, 0.5 * DeltaZMed(i) / TmpCoeff)
        end if
    end do

    NSteps = ceiling(Exch_dT / DeltaT)
    DeltaT = Exch_dT / real(NSteps)
    
#ifdef DEBUG_MODE
    print *, '<- Exit TimeStep'
#endif
end subroutine TimeStep

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Right-hand parts of exchange differential equations
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CalcRHP(ic,jc)
  integer   ic, jc

  integer i, j, k, Src                ! 03.02.2015
  real(8) DMass, RMass                ! 04.02.2015 19-09-2019
  real(8)   AirMedFlux, MedAirFlux    ! 19-09-2019

#ifdef DEBUG_MODE
    print *, '>- Entering CalcRHP'
#endif
! Wet deposition --------------------------------------------------------------------------------
  RHP(Exch_KMax) = - WetDepVel(Exch_KMax) / DeltaZ(Exch_KMax)  * CAtm(Exch_KMax)
  RHP(Exch_KMax_2) = -  WetDepCoeff(Exch_KMax)* CAtm(Exch_KMax_2)
! Wet deposition for upper atm layer (Matrix calculations 03.02.2015)
  do Src = 1, NumSrc               ! 03.02.2015
#if RTYPE==2
    DMass = WetDepCoeff(Exch_KMax)* CAtm(Exch_KMax_2)*deltaT*MV(Exch_KMax)* Contrib(Exch_KMax,Src)   ! 04.02.2015
#endif
#if RTYPE==1
    DMass = WetDepCoeff(Exch_KMax)* CAtm(Exch_KMax_2)*deltaT*MV(Exch_KMax)   ! 04.02.2015
#endif
    MassWDSrc(2,Src)=MassWDSrc(2,Src)+ DMass   ! 04.02.2015
    MassWetDep(2)=MassWetDep(2) + DMass   ! 04.02.2015
  end do ! Src

  do i = 1, Exch_KMax - 1
    RHP(i) = - WetDepVel(i) / DeltaZ(i) * CAtm(i)+ WetDepVel(i + 1) / DeltaZ(i) * CAtm(i + 1)
    RHP(i+Exch_KMax) = -  WetDepCoeff(i) * CAtm(i+Exch_KMax)
#if RTYPE==2
! Wet deposition for internal layers (Matrix calculations 24.12.2014)
    DCOut(i) = WetDepVel(i) / DeltaZ(i) * CAtm(i) + WetDepCoeff(i) * CAtm(i+Exch_KMax)
    DCIn(i) = WetDepVel(i + 1) / DeltaZ(i) * CAtm(i + 1)
#endif
    do Src = 1, NumSrc               ! 03.02.2015
#if RTYPE==2
     DMass = WetDepCoeff(i)* CAtm(i+Exch_KMax) * deltaT * MV(i) * Contrib(i,Src)   ! 04.02.2015
#endif
#if RTYPE==1
     DMass = WetDepCoeff(i) * CAtm(i+Exch_KMax) * deltaT * MV(i)   ! 04.02.2015
#endif
      MassWDSrc(2,Src)=MassWDSrc(2,Src)+DMass   ! 04.02.2015
      MassWetDep(2)=MassWetDep(2) + DMass   ! 04.02.2015
    end do ! Src
  end do ! i

! Wet deposition for lower atm layer (Matrix calculations 03.02.2015)
  do Src = 1, NumSrc               ! 03.02.2015
#if RTYPE==2
    DMass = WetDepVel(1) * CAtm(1)*deltaT * Contrib(1,Src) * MA
#endif
#if RTYPE==1
    DMass = WetDepVel(1) * CAtm(1) * deltaT * MA ! 04.02.2015
#endif
    MassWDSrc(1,Src)=MassWDSrc(1,Src)+ Dmass ! 04.02.2015
    MassWetDep(1)=MassWetDep(1) + Dmass ! 04.02.2015
  end do

#if RTYPE==2
  DCInMed = 0.  ! Matrix calculations 24.12.2014
#endif

! Dry deposition --------------------------------------------------------------------------------
  do i = 1, Exch_KMed
    do Src = 1, NumSrc
#if RTYPE==2
!      DMass = -(ExchVel(i,0)*CMed(i)-ExchVel(0,i)*CAtm(1))*LC(i)*deltaT*MA*Contrib(1,Src)
       DMass = ExchVel(0,i)*CAtm(1)*LC(i)*deltaT*MA*Contrib(1,Src)
       if(Src == 1) RMass = ExchVel(i,0)*CMed(i)*LC(i)*deltaT*MA
#endif
#if RTYPE==1
      DMass = ExchVel(0,i) * CAtm(1) * LC(i) *deltaT * MA
      RMass = ExchVel(i,0) * CMed(i) * LC(i) *deltaT * MA
#endif
      MassDDSrc(1,i,Src)=MassDDSrc(1,i,Src) + DMass            ! 18-09-2019
      MassDryDep(1) = MassDryDep(1) + Dmass
      if(Src == 1) then
            MassDRSrc(1,i,Src)=MassDRSrc(1,i,Src) + RMass      ! 18-09-2019
            MassDryRem(1) = MassDryRem(1) + Rmass
      end if
#if RTYPE==2
      DMass = DryDepVel(i) * CAtm(1+Exch_KMax) * LC(i) *deltaT*MA * Contrib(1,Src)
#endif
#if RTYPE==1
      DMass = DryDepVel(i) * CAtm(1+Exch_KMax) * LC(i) *deltaT*MA
#endif
      MassDDSrc(2,i,Src)=MassDDSrc(2,i,Src) + DMass               ! 18-09-2019
      MassDryDep(2) = MassDryDep(2) + DMass
    end do
!---------------------------------------------------------------------------------------------------
    RHP(1) = RHP(1) + (ExchVel(i,0) * CMed(i) - ExchVel(0,i) * CAtm(1)) * LC(i) / DeltaZ(1)
    RHP(1+Exch_KMax) = RHP(1+Exch_KMax) - DryDepVel(i) * CAtm(1+Exch_KMax)* LC(i) / DeltaZ(1)
#if RTYPE==2
    DCOut(1) = DCOut(1) + (ExchVel(0,i) * CAtm(1) + DryDepVel(i) * CAtm(1+Exch_KMax))* LC(i) / DeltaZ(1)
    DCInMed = DCInMed + ExchVel(i,0) * CMed(i) * LC(i) / DeltaZ(1)
#endif
  end do
  do i = 1, Exch_KMed
!  where(CMed < 0.) CMed = 0.           ! 03.02.2015
    RHP(Exch_KMax_2 + i) = (- ExchVel(i,0) * CMed(i) + &
        &(WetDepVel(1) * WetDepMask(i) + ExchVel(0,i)) * CAtm(1) + &
        &DryDepVel(i) * CAtm(1+Exch_KMax)) / DeltaZMed(i)

    do k = 1, Exch_KMax
      RHP(Exch_KMax_2 + i) = RHP(Exch_KMax_2+  i)&
      &+ WetDepCoeff(k) * WetDepMask(i) * DeltaZ(k) * CAtm(k+Exch_KMax)&
      &/ DeltaZMed(i)
    end do
   end do        ! 19-09-2019

!---------------------------------------------------------------------------------------------------
! Counters for fluxes from/to media
!---------------------------------------------------------------------------------------------------
! Soil
  AirMedFlux = 0.
  MedAirFlux = 0.
  do i = 1, Soil_NumType
    AirMedFlux = AirMedFlux + ((WetDepVel(1)*WetDepMask(i)&
                              &+ExchVel(0,i))*CAtm(1)&
                              &+DryDepVel(i)*CAtm(1+Exch_KMax))* LC(i)
    do k = 1, Exch_KMax
      AirMedFlux = AirMedFlux + WetDepCoeff(k)*WetDepMask(i)*DeltaZ(k)*CAtm(k+Exch_KMax)*LC(i)
    end do
    MedAirFlux = MedAirFlux - ExchVel(i,0) * CMed(i) * LC(i)
  end do
  AirSoilFlux = AirSoilFlux + AirMedFlux*deltaT*MA
  SoilAirFlux = SoilAirFlux + MedAirFlux*deltaT*MA
!---------------------------------------------------------------------------------------------------
! Vegetation
  AirMedFlux = 0.
  MedAirFlux = 0.
  do i = Soil_NumType+1, CMedSoilMax
    AirMedFlux = AirMedFlux + ((WetDepVel(1)*WetDepMask(i)&
                              &+ExchVel(0,i))*CAtm(1)&
                              &+DryDepVel(i)*CAtm(1+Exch_KMax))* LC(i)
    do k = 1, Exch_KMax
      AirMedFlux = AirMedFlux+WetDepCoeff(k)*WetDepMask(i)*DeltaZ(k)*CAtm(k+Exch_KMax)*LC(i)
    end do
    MedAirFlux = MedAirFlux - ExchVel(i,0) * CMed(i) * LC(i)
  end do
  AirVegFlux = AirVegFlux + AirMedFlux*deltaT*MA
  VegAirFlux = VegAirFlux + MedAirFlux*deltaT*MA
!---------------------------------------------------------------------------------------------------
! Ocean
  i = Exch_KMed
  AirMedFlux = 0.
  MedAirFlux = 0.
  AirMedFlux = AirMedFlux + ((WetDepVel(1)*WetDepMask(i)+ExchVel(0,i))*CAtm(1)+DryDepVel(i)*CAtm(1+Exch_KMax))* LC(i)
  do k = 1, Exch_KMax
    AirMedFlux = AirMedFlux + WetDepCoeff(k)*WetDepMask(i)*DeltaZ(k)*CAtm(k+Exch_KMax)*LC(i)
  end do
  MedAirFlux = MedAirFlux - ExchVel(i,0) * CMed(i) * LC(i) 
  AirOcnFlux = AirOcnFlux + AirMedFlux*deltaT*MA
  OcnAirFlux = OcnAirFlux + MedAirFlux*deltaT*MA  
! End counters for fluxes from/to media 19-09-2019
!---------------------------------------------------------------------------------------------------
  do i = 1, Exch_KMed
   do j = 1, Exch_KMed
    RHP(Exch_KMax_2 + i) = RHP(Exch_KMax_2 + i) +&
        &(ExchVel(j,i) * CMed(j)- ExchVel(i,j) * CMed(i))/ DeltaZMed(i)

    if ((i.le.Soil_NumType).and.(j.gt.Soil_NumType.and.j.le.(CMedSoilMax))) &
        & VegSoilFlux = VegSoilFlux + (ExchVel(j,i) * CMed(j)- ExchVel(i,j) * CMed(i)) * MA * LC(i) * deltaT
   end do
  end do

#ifdef DEBUG_MODE
    print *, '<- Exit CalcRHP'
#endif
end subroutine CalcRHP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Gauss method for exchange differential equations
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ExchCalc(ic,jc)
  integer   ic, jc
  integer   i, Step
  integer   Src, Cont

#ifdef DEBUG_MODE
    print *, '-> Entering ExchCalc ', ic,jc
#endif


  call TimeStep

  do Step = 1, NSteps
#if RTYPE==2
    do i = 1, Exch_KMax - 1
      CABeg(i) = CAtm(i) + CAtm(Exch_KMax + i)
    end do
#endif
    call CalcRHP(ic,jc)
    do i = 1, Exch_KMax_2
      CAtm(i) = CAtm(i) + RHP(i) * DeltaT
    end do
    do i = 1, Exch_KMed
      CMed(i) = CMed(i) + RHP(Exch_KMax_2+ i) * DeltaT
    end do
#if RTYPE==2
    do i = 1, Exch_KMax - 1
      CAFin(i) = CAtm(i) + CAtm(Exch_KMax + i)
      if(CAFin(i) /= 0.) then
         ! update contribution of anthropogenic sources
          do Src = 1, NumSrc
            if( Src == Reem) then
                  if (i /= 1) Contrib(i,Reem) = ((CABeg(i) - DCOut(i)*DeltaT) * Contrib(i,Reem)&
                    &+ DCIn(i)*DeltaT * Contrib(i + 1,Reem)) / CAFin(i)
                  if (i == 1) Contrib(i,Reem) = ((CABeg(i) - DCOut(i)*DeltaT) * Contrib(i,Reem)&
                    &+ DCIn(i)*DeltaT * Contrib(i + 1,Reem) + DCInMed*DeltaT) / CAFin(i)
            else
            Contrib(i,Src) = ((CABeg(i) - DCOut(i)*DeltaT) * Contrib(i,Src) + DCIn(i)*DeltaT * Contrib(i + 1,Src)) / CAFin(i)
            end if
          end do
      end if
    end do
#endif
  end do

#ifdef DEBUG_MODE
    print *, '<- Exit ExchCalc'
#endif
end subroutine ExchCalc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Exchange between media - percistent organic pollutants
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Exch_POP(NSubs)

    integer Xscal,l,iB,iE
    integer NSubs
    integer i, j, k, n, f, m, surf, src           ! 18-09-2019
    integer ind1,ind2
    integer k1, k2, ind1, ind2

    real(8) Conc_old, Conc_diff, TASUM,wupl, t1
    real(8) Summ1, Summ2, TSurf, Henry
    real(8) atmp(Imin:Imax, Jmin:Jmax), P1(Imin:Imax)


#ifdef DEBUG_MODE
    print *, '+Entering Exch_POP ', NSubs
#endif

    ExchVel=0.
    WetDepVel=0.
    WetDepCoeff=0.
    DryDepVel=0.

    DeltaZMed(iOcn) = Water_upl_dz
    n = sFormInd(Ocn,NSubs,1)
    Water_upl_conc = 0.
    ind1 = sFormInd(Ocn,NSubs,1)
    ind2 = sFormInd(Ocn,NSubs,2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OCEAN 3-D -> 1-D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if ((REGTYPE == 1) && (GRIDIMAX == 360))
! Copy Ocn_conc into temp array with 0.5 grid cell shift
    do j = JMin, JMax
        do i = IMin, IMax
            k1 = i
            if(i < IMax) then
                k2 = i + 1
            else
                k2 = 1
            end if
! Here Water_upl_conc is used for collecting pollutant amount per 1 m of depth
            Water_upl_conc(i,j,ind1) = (Ocn_conc(k1,j,1,n,curtime)*TAREA(k1,j) + Ocn_conc(k2,j,1,n,curtime)*TAREA(k2,j)) / 20000.
        end do
    end do
#else
    do j = JMin, JMax
        do i = IMin, IMax
! Here Water_upl_conc is used for collecting pollutant amount per 1m of depth
            Water_upl_conc(i,j,ind1) = Ocn_conc(i,j,1,n,curtime)*TAREA(i,j) / 10000.
        end do
    end do
#endif
! Grid aggregation
    do j = Jmin, Jmax
        Xscal = Imax/maxI(j)
        if(Xscal > 1) then
            do i = minI(j), maxI(j)
                t1 = 0.
                do l = 1, Xscal
                    t1 = t1 + Water_upl_conc(i*Xscal-l+1,j,ind1)
                end do
                Water_upl_conc(i,j,ind1) = t1
            end do ! i
        end if

          do i = minI(j), maxI(j)
              if(Ocn_Frac(i,j) > 0.) then
! From amounts to concentrations calculated similar to air
                Water_upl_conc(i,j,ind1) = Water_upl_conc(i,j,ind1)/Ocn_Frac(i,j)/MeshArea(i,j)
              end if
          end do
    enddo ! j

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OCEAN 3-D -> 1-D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call Ocn_POPpartitioning(NSubs)

!%%%%%%%%%%%% PREPARATION OF ARRAYS (CATM, CMED, EXCHVEL, ...) %%%%%%%%%%%%%%%%%%
!%%%
!%%%   CMED indices:
!%%%     Atm:           0
!%%%     Soil:          1                              -  Soil_NumType
!%%%     Vegetation:    1+Soil_NumType                 -  Soil_NumType+2*Veg_NumType
!%%%     Forest Litter: 1+Soil_NumType+2*Veg_NumType   -  Soil_NumType+3*Veg_NumType
!%%%     Ocean:         1+Soil_NumType+3*Veg_NumType   -  Soil_NumType+3*Veg_NumType+1
    
    do j = JMin, JMax
    do i = MinI(j), MaxI(j)

      MassDDSrc = 0.
      MassDRSrc = 0.
      MassWDSrc = 0.

       MA = MeshArea(i,j)
       do k=1,Exch_KMax
        MV(k)=MeshVolume(i,j,k)
        DeltaZ(k) = MV(k)/MA
       end do

	  LC(1:Soil_NumType) = Soil_Frac(i,j,1:Soil_NumType)
          WetDepMask(1:Soil_NumType) = 1.
	  LC(Soil_NumType + 3*Veg_NumType +1) = Ocn_Frac(i,j)
          WetDepMask(Soil_NumType + 3*Veg_NumType +1) = 1.
          WetDepVel(1:Exch_KMax) = Atm_WDVelGas(i,j,1:Exch_KMax,NSubs)
          WetDepCoeff(1:Exch_KMax) = Atm_WDCoeff(i,j,1:Exch_KMax,NSubs)

          do k = 1, Soil_NumType-2
	    ExchVel(0,k) = Soil_VExch(i,j,k,Atm,1,NSubs)            ! From Atm to Media - 1
	    ExchVel(k,0) = Soil_VExch(i,j,k,Atm,2,NSubs)            ! From Media to Atm - 2
	  end do
          do k = 1, Veg_NumType
	    ExchVel(0,Soil_NumType+k) = Veg_VExch(i,j,k,Atm,1,NSubs)
	    ExchVel(Soil_NumType+k,0) = Veg_VExch(i,j,k,Atm,2,NSubs)
	  end do
!!!-Vegetation-!!!
          do k = 1, Veg_NumType
            LAIPriv(i,j,k) = LAIPriv(i,j,k) + AlphaConst(i,j,k) * Exch_dT
            if((AlphaConst(i,j,k) < 0.) .and. (LAIPriv(i,j,k) > 0.)) then
            ExchVel(Soil_NumType+k,CMedVegMax+k)=-AlphaConst(i,j,k)/LAIPriv(i,j,k)
            ExchVel(Soil_NumType+Veg_NumType+k,CMedVegMax+k)=-AlphaConst(i,j,k)/LAIPriv(i,j,k)
            ExchVel(CMedVegMax+k,k)=CEFall
#ifdef WO_DRYDEPGASVEG
            ExchVel(Soil_NumType+k,CMedVegMax+k) = 0.
            ExchVel(Soil_NumType+Veg_NumType+k,CMedVegMax+k) = 0.
#endif
            end if
	  end do
!!!-!!!
          do k = 1, Soil_NumType
            DryDepVel (k)= POPVdAero(i,j,k,NSubs)*(1.-Throughfall)
          end do
          do k = 1, Veg_NumType
            DryDepVel (Soil_NumType + Veg_NumType + k)= POPVdAero(i,j,k,NSubs) * Throughfall
          end do
          DryDepVel(iOcn) = POPVdAero(i,j,0,NSubs)
	  ExchVel(0,iOcn) = Ocn_VExch(i,j,Atm,1,NSubs)  ! *******
	  ExchVel(iOcn,0) = Ocn_VExch(i,j,Atm,2,NSubs)    ! *******
!!!-Atmosphere-!!!
          CAtm(1:Exch_KMax) = Atm_Conc(i,j,1:Exch_KMax,sFormInd(Atm,NSubs,1))
          CAtm(Exch_KMax+1:Exch_KMax_2) = Atm_Conc(i,j,1:Exch_KMax,sFormInd(Atm,NSubs,2))
#if RTYPE==2
          do k= 1, Exch_KMax
            Contrib(k,1:NumSrc) = Atm_Contrib(i,j,k,sFormInd(Atm,NSubs,1),1:NumSrc)
          end do
#endif
!!!-Soil-!!!
	  CMed(1:Soil_NumType) = Soil_Conc(i,j,1,1:Soil_NumType,sFormInd(Soil,NSubs,1)) * (Poros - wl) + &
          & Soil_Conc(i,j,1,1:Soil_NumType,sFormInd(Soil,NSubs,2)) * wl + Soil_Conc(i,j,1,1:Soil_NumType,sFormInd(Soil,NSubs,3))
!!!-Ocean-!!!
          if(Ocn_Frac(i,j) > 0.) then
            CMed(iOcn) = Water_upl_conc(i,j,sFormInd(Ocn,NSubs,1))
          else
              CMed(iOcn) = 0.
          end if
!!!-Vegetation-!!!
	  do f = 1, NumForm(Veg)
              ind1=Soil_NumType + (f-1)*Veg_NumType + 1
              ind2=Soil_NumType + f*Veg_NumType
              CMed(ind1:ind2) = Veg_Conc(i,j,1:Veg_NumType,sFormInd(Veg,NSubs,f))
              LC(ind1:ind2) = Veg_Frac(i,j,1:Veg_NumType)
              WetDepMask(ind1:ind2) = 0.
	  end do
!!!-ForestLitter-!!!
	  ind1=Soil_NumType + NumForm(Veg)*Veg_NumType + 1
	  ind2=Soil_NumType + (NumForm(Veg)+1)*Veg_NumType
	  CMed(ind1:ind2) = Fall_Conc(i,j,1:Veg_NumType,gSubsMediaInd(Veg,NSubs))
	  LC(ind1:ind2) = Veg_Frac(i,j,1:Veg_NumType)
          WetDepMask(ind1:ind2) = 0.
!!!-!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

#ifdef WO_DRYDEPGAS
          ExchVel = 0.
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN CYCLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  call ExchCalc(i,j)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%% UPDATING OF CONC AND DEP ARRAYS %%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef O_EQUIL
          Tsurf = TempSurf(i, j, Period)
          Henry = HT0(gSubsGroupInd(NSubs)) / (RUniv * TSurf) * exp(- HT(gSubsGroupInd(NSubs)) * (1 / TSurf - 1 / T0))
          Summ1 = CAtm(1) * DeltaZ(1) + CMed(iOcn) * water_upl_dz * Ocn_Frac(i,j)
          Summ2 = Henry * DeltaZ(1) + water_upl_dz *Ocn_Frac(i,j)
          CMed(iOcn) = Summ1 / Summ2
          CAtm(1) = Henry * CMed(iOcn)
#endif
!!!-Atmosphere-!!!
          do k=1,Exch_KMax
              Atm_Conc(i,j,k,sFormInd(Atm,NSubs,2)) = CAtm(Exch_KMax+k)     ! Update part phase
              Atm_Conc(i,j,k,sFormInd(Atm,NSubs,1)) = CAtm(k)               ! Update gas phase
          end do
#if RTYPE==2
          Atm_Contrib(i,j,1:Exch_KMax-1,sFormInd(Atm,NSubs,1),1:NumSrc)= Contrib(1:Exch_KMax-1,1:NumSrc)
          Atm_Contrib(i,j,1:Exch_KMax-1,sFormInd(Atm,NSubs,2),1:NumSrc)= Contrib(1:Exch_KMax-1,1:NumSrc)
#endif

!!!-Soil-!!!
	  Soil_Conc(i,j,1,1:Soil_NumType,sFormInd(Soil,NSubs,1)) = CMed(1:Soil_NumType) / PartG(i, j, 1, Period, gSubsGroupInd(NSubs))
	  Soil_Conc(i,j,1,1:Soil_NumType,sFormInd(Soil,NSubs,2)) = CMed(1:Soil_NumType) / PartL(i, j, 1, Period, gSubsGroupInd(NSubs))
	  Soil_Conc(i,j,1,1:Soil_NumType,sFormInd(Soil,NSubs,3)) = &
          & CMed(1:Soil_NumType) - Soil_Conc(i,j,1,1:Soil_NumType,sFormInd(Soil,NSubs,1)) * (Poros - wl) - &
	  & Soil_Conc(i,j,1,1:Soil_NumType,sFormInd(Soil,NSubs,2)) * wl
!!!-Ocean-!!!
          if(Ocn_Frac(i,j) > 0.) then
            Water_upl_conc(i,j,sFormInd(Ocn,NSubs,1)) = CMed(iOcn)
          end if
!!!-Vegetation-!!!
	  do f = 1, NumForm(Veg)
              ind1=Soil_NumType + (f-1)*Veg_NumType + 1
              ind2=Soil_NumType + f*Veg_NumType
              Veg_Conc(i,j,1:Veg_NumType,sFormInd(Veg,NSubs,f)) = CMed(ind1:ind2)
	  end do

!!!-ForestLitter-!!!
	  ind1=Soil_NumType + NumForm(Veg)*Veg_NumType + 1
	  ind2=Soil_NumType + (NumForm(Veg)+1)*Veg_NumType
	  Fall_Conc(i,j,1:Veg_NumType,gSubsMediaInd(Veg,NSubs)) = CMed(ind1:ind2)
!!!-DepositionFluxes_!!!
          do f=1, NumForm(Atm)
            do Surf = 1, Exch_KMed           ! 18-09-2019   Corrected 16-10-2019
		DryRemDayTmp(i,j,f,surf,1)=DryRemDayTmp(i,j,f,surf,1)+MassDRSrc(f,surf,1)
		do Src=1, NumSrc
                 DryDepDayTmp(i,j,f,surf,Src)=DryDepDayTmp(i,j,f,surf,Src)+MassDDSrc(f,surf,Src)
                end do
	    end do
	    do Surf = 1, Soil_NumType
	      do Src=1, NumSrc
	        WetDepDayTmp(i,j,f,surf,Src)=WetDepDayTmp(i,j,f,surf,Src)+MassWDSrc(f,Src) * LC(surf)
              end do
	    end do
	      do Src=1, NumSrc
	        WetDepDayTmp(i,j,f,Exch_KMed,Src)=WetDepDayTmp(i,j,f,Exch_KMed,Src)+&
	        &MassWDSrc(f,Src)*LC(Exch_KMed)
              end do
          end do
          PrecipDay(i,j)=PrecipDay(i,j)+(PrecRainConv(i,j,1,Period)+PrecRainStrat(i,j,1,Period))*deltaT
!!!-!!!
	end do
  end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OCEAN 1-D -> 3-D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind1 = sFormInd(Ocn,NSubs,1)
    ind2 = sFormInd(Ocn,NSubs,2)

! Convert from concentrations to amount of the pollutant (kg/m)
    do j = JMin, JMax
        do i = minI(j), maxI(j)
            if(Ocn_Frac(i,j) > 0.) then
                atmp(i,j) = (Water_upl_conc(i,j,ind1)+Water_upl_conc(i,j,ind2))*Ocn_Frac(i,j)*MeshArea(i,j)
            else
                atmp(i,j) = Water_upl_conc(i,j,ind1)+Water_upl_conc(i,j,ind2)
            end if
        end do
    end do

! Disaggregate
    do j = Jmin, Jmax
        Xscal = Imax/maxI(j)
        do i = minI(j), maxI(j)
	    iB=(i-1)*Xscal+1
	    iE=iB+Xscal-1
	    do l=iB, iE
                P1(l)=atmp(i,j)/real(Xscal)
            enddo
        enddo
        atmp(Imin:Imax,j)=P1(Imin:Imax)
    end do ! j

! Update upper layer concentrations in the Ocn_conc array and convert back to kg/m3
#if ((REGTYPE == 1) && (GRIDIMAX == 360))
    do j = JMin, JMax
        t1 = (atmp(IMin, j) + atmp(IMax, j))/2. / TArea(IMin,j) * 10000.! from amounts to concentrations (kg/m3)
        Ocn_conc(IMin,j,1,n,curtime)= t1
        Ocn_conc(IMin,j,1,n,oldtime)= t1
        do i = IMin+1, IMax
            t1 = (atmp(i-1,j) + atmp(i,j))/2. / TArea(i,j) * 10000.  ! from amounts to concentrations (kg/m3)
            Ocn_conc(i,j,1,n,curtime)= t1
            Ocn_conc(i,j,1,n,oldtime)= t1
        end do ! i
    end do ! j
#else
    do j = JMin, JMax
        do i = IMin, IMax
            t1 = atmp(i,j) / TArea(i,j) * 10000.  ! from amounts to concentrations
            Ocn_conc(i,j,1,n,curtime)= t1
            Ocn_conc(i,j,1,n,oldtime)= t1
        end do ! i
    end do ! j
#endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifdef DEBUG_MODE
    print *, '+Exit Exch_POP ', NSubs
#endif
end subroutine Exch_POP

end Module Exch_General_POP
#endif