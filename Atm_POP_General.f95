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
module Atm_POP_General

  use GeneralParams
  use Atm_Params
  use Geometry
  use Atm_POP_Params
  use Exch_Params
  use Soil_Params
  use typeSizes
  use netcdf
  implicit none


  character(800), private :: fileName, fullName
  character(4),  private :: YearNum
  character(2),  private :: DayNum, YearShort
  logical, private :: ROI_T = .false.
  integer O3Nmb

  real, private :: ArrTmp(IMin:IMax,JMin:JMax,Atm_KMax,NumPer)
  real cosSol, cosMean,  Fw

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine POP specific processes
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_POP_Process(Subs)

    integer Subs

#if DEBUG_MODE
    print *, '>- Entering Atm_POP_Process '
#endif

    call Atm_POP_ReactTransformSteply
    call Atm_POP_PartitDegrCoefs
    call Atm_POP_Partitioning(Subs)
    call Atm_POP_Degradation(Subs)
    
#if DEBUG_MODE
    print *, '>- Exit Atm_POP_Process '
#endif

end subroutine Atm_POP_Process


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine defining physical and chemical properties of POP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_POP_Props(Subs)

	integer i, FileStat, lenType, lenType1, Subs, Form, Limit, Ind, Src
	character(80) strRead, strPar
	character(30) strType, FormType(MaxForm)
	integer mon, ios, RN, AN

    fileName='Atm'//trim(PropName)//trim(SubsID(Subs))//'.dat'
	fullName=trim(PropPath)//fileName
	open(2, file=trim(fullName), status='old', iostat=FileStat, action='read')
	if(FileStat>0) then
	  print '(/,"STOP: Cannot open file with POP properties ''",a,"''",/)', trim(fileName)
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
		do i=1, sFormNum(Atm,Subs)
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
!------------------------------------------------------------------
		  case('Aerosols')       ! 05.07.2018
		    read(strPar,*) NmbAero
		    do AN = 1, NmbAero
		      read(2,*) AName(AN), AUnit(AN), Cp0(AN), CpT(AN), CorrectA(AN)
		    end do
		  case('ROI_T')       ! 05.07.2018
		    if (trim(strPar) == 'yes') then
		      ROI_T = .true.
		      call ReadROI_TParams
		    end if
		  case('Water_Correction')       ! 17.12.2019
		    if (trim(strPar) == 'yes') then
		      WaterCorrect = .true.
		    end if
		  case('Reactants')       ! 05.07.2018
		    read(strPar,*) NmbReacts
                    O3Nmb = 0
		    do RN = 1, NmbReacts
		      read(2,*) RName(RN), RUnit(RN), Kd0(RN,0), KdT(RN,0), & 
				&((Kd0(RN,AN),KdT(RN,AN), KReact(RN,AN)), AN = 1,NmbAero)
                      if (trim(RName(RN)) == 'O3_') then
                        O3Nmb = RN
                      end if
		    end do
		  case('Forms number')
			read(strPar,'(i2)') sFormNum(Atm,Subs)
            FormSubs(Atm,NumForm(Atm)+1:NumForm(Atm)+sFormNum(Atm,Subs))=Subs
			if (sFormNum(Atm,Subs).gt.0) then
			NumSubsMedia(Atm)=NumSubsMedia(Atm)+1
			gSubsMediaInd(Atm,Subs)=NumSubsMedia(Atm)
			end if
#ifdef DEBUG_MODE
			print *, 'sFormNum(Atm,Subs): ',sFormNum(Atm,Subs)
#endif
!------------------------------------------------------------------
		  case('Forms')
			read(strPar,*) (FormType(i), i=1, sFormNum(Atm,Subs))
#ifdef DEBUG_MODE
			print *, 'Atm form types: ',FormType
#endif
!------------------------------------------------------------------
		  case('Form ID')
                        FormID(Atm,Subs,Form)=strPar
                        sFormInd(Atm,Subs,Form)=NumForm(Atm)+1
                        FormSubs(Atm,sFormInd(Atm,Subs,Form))=Subs
                        NumForm(Atm)=NumForm(Atm)+1

#ifdef DEBUG_MODE
                        print *, 'Form: ', Form
                        print *, 'FormID(Atm,Subs,Form): ',FormID(Atm,Subs,Form)
                        print *, 'sFormInd(Atm,Subs,Form): ', sFormInd(Atm,Subs,Form)
                        print *, 'FormSubs(Atm,sFormInd(Atm,Subs,Form)): ', FormSubs(Atm,sFormInd(Atm,Subs,Form))
                        print *, 'NumForm(Atm): ', NumForm(Atm)
                        print *, ' '
#endif
!------------------------------------------------------------------
		  case('Dry deposition')
                        if(strPar=='yes') then
			  DryDepNum=DryDepNum+1
			  DryDepInd(DryDepNum)=sFormInd(Atm,Subs,Form)
			endif
!------------------------------------------------------------------
		  case('Below-cloud scavenging')
                        if(strPar=='yes') then
			  BCScvNum=BCScvNum+1
			  BCScvInd(BCScvNum)=sFormInd(Atm,Subs,Form)
			endif
!------------------------------------------------------------------
		  case('In-cloud scavenging')
                        if(strPar=='yes') then
			  ICScvNum=ICScvNum+1
			  ICScvInd(ICScvNum)=sFormInd(Atm,Subs,Form)
			endif
!------------------------------------------------------------------
		  case('Gas exchange')
                        if(strPar=='yes') then
			  GasExchNum=GasExchNum+1
			  GasExchInd(GasExchNum)=sFormInd(Atm,Subs,Form)
			endif
!------------------------------------------------------------------
		  case('Atm transport')
                        if(strPar=='yes') then
			  AtmTransNum=AtmTransNum+1
			  AtmTransInd(AtmTransNum)=sFormInd(Atm,Subs,Form)
			endif
!------------------------------------------------------------------
                  case('Boundary conditions')               ! Modified Aug 2016 V. Sh
                        if(strPar=='yearly') then
			  AtmBndNum=AtmBndNum+1
			  AtmBndInd(AtmTransNum)=sFormInd(Atm,Subs,Form)
			  ReadBndInd(AtmTransNum,yr)=sFormInd(Atm,Subs,Form)
                        elseif(strPar=='monthly') then
			  AtmBndNum=AtmBndNum+1
			  AtmBndInd(AtmTransNum)=sFormInd(Atm,Subs,Form)
			  ReadBndInd(AtmTransNum,mn)=sFormInd(Atm,Subs,Form)
                    elseif(strPar=='daily') then
			  AtmBndNum=AtmBndNum+1
			  AtmBndInd(AtmTransNum)=sFormInd(Atm,Subs,Form)
			  ReadBndInd(AtmTransNum,da)=sFormInd(Atm,Subs,Form)
                    else
                    endif
!------------------------------------------------------------------
		  case('Anthrop emiss')
		    if(strPar=='yearly') then
			AntEmisNum=AntEmisNum+1
                        ReadAntInd(AntEmisNum,yr)=sFormInd(Atm,Subs,Form)
                        AntEmisInd(AntEmisNum)=sFormInd(Atm,Subs,Form)
                    elseif(strPar=='monthly') then
 			AntEmisNum=AntEmisNum+1
                        ReadAntInd(AntEmisNum,mn)=sFormInd(Atm,Subs,Form)
			AntEmisInd(AntEmisNum)=sFormInd(Atm,Subs,Form)
                    elseif(strPar=='daily') then
			AntEmisNum=AntEmisNum+1
                        ReadAntInd(AntEmisNum,da)=sFormInd(Atm,Subs,Form)
			AntEmisInd(AntEmisNum)=sFormInd(Atm,Subs,Form)
                    else
		    endif
!------------------------------------------------------------------
		  case('Natural emiss')
		    if(strPar=='yearly') then
 			NatEmisNum=NatEmisNum+1
                        ReadNatInd(NatEmisNum,yr)=sFormInd(Atm,Subs,Form)
			NatEmisInd(NatEmisNum)=sFormInd(Atm,Subs,Form)
		    elseif(strPar=='monthly') then
			NatEmisNum=NatEmisNum+1
                        ReadNatInd(NatEmisNum,mn)=sFormInd(Atm,Subs,Form)
			NatEmisInd(NatEmisNum)=sFormInd(Atm,Subs,Form)
		    elseif(strPar=='daily') then
			NatEmisNum=NatEmisNum+1
                        ReadNatInd(NatEmisNum,da)=sFormInd(Atm,Subs,Form)
			NatEmisInd(NatEmisNum)=sFormInd(Atm,Subs,Form)
                    else
		    endif
!------------------------------------------------------------------
		  case('Reemission')
		    if(strPar=='yearly') then
			ReEmisNum=ReEmisNum+1
                        ReadReInd(ReEmisNum,yr)=sFormInd(Atm,Subs,Form)
			ReEmisInd(ReEmisNum)=sFormInd(Atm,Subs,Form)
		    elseif(strPar=='monthly') then
			ReEmisNum=ReEmisNum+1
                        ReadReInd(ReEmisNum,mn)=sFormInd(Atm,Subs,Form)
			ReEmisInd(ReEmisNum)=sFormInd(Atm,Subs,Form)
		    elseif(strPar=='daily') then
			ReEmisNum=ReEmisNum+1
                        ReadReInd(ReEmisNum,da)=sFormInd(Atm,Subs,Form)
			ReEmisInd(ReEmisNum)=sFormInd(Atm,Subs,Form)
                    else
                    endif
!------------------------------------------------------------------
                case('Henry')
		    read(strPar,*) HT0(gSubsGroupInd(Subs)),HT(gSubsGroupInd(Subs))
!------------------------------------------------------------------
                case('Pol')
		    read(strPar,*) p0(gSubsGroupInd(Subs)),pT(gSubsGroupInd(Subs))
!------------------------------------------------------------------
                case('Koa')
		    read(strPar,*) Koa0(gSubsGroupInd(Subs)),KoaT(gSubsGroupInd(Subs))
!------------------------------------------------------------------
		case('DMolA')
		    read(strPar,*) DAir(gSubsGroupInd(Subs))
!------------------------------------------------------------------
		case('DMolW')
		    read(strPar,*) DWater(gSubsGroupInd(Subs))
!------------------------------------------------------------------
		case('Particle diameter, m')
			read(strPar,*) Dp
!------------------------------------------------------------------
		case('Particle density, kg/m3')
			read(strPar,*) RhoP
!------------------------------------------------------------------
		case('BCS coeffs')
			read(strPar,*) Abelow, Bbelow
!------------------------------------------------------------------
		case('ICS coeffs')
			read(strPar,*) Ain, Bin
!------------------------------------------------------------------
		case('Washout effectiveness coeff')
			read(strPar,*) Aeff
!------------------------------------------------------------------
		case('Aerosol solubility')
			read(strPar,*) Asol
!------------------------------------------------------------------
                case('Particle washout ratio')
			read(strPar,*) WRatioPart
!------------------------------------------------------------------
		case default
		    print '(/,"STOP: Unknown input parameter ''",a,"''",/)', trim(strType)
!		    stop
		endselect

	  endif
	enddo

	close(2)

#if RTYPE==2
    if(NumAnth==0.and.AntEmisNum/=0) then
	  print '(/,"STOP: No emission sources is defined for anthropogenic emissions of forms:",<AntEmisNum>(1x,a),/)',&
                            &(trim(SubsID(Subs))//trim(FormID(Atm,Subs,AntEmisInd(Ind))), Ind=1, AntEmisNum)
	  stop
    elseif(NumAnth/=0.and.AntEmisNum==0) then
	  print '(/,"STOP: No anthropogenic emission is defined for sources:",<NumAnth>(1x,a),/)',&
                            &(trim(SourcID(Src)), Src=1, NumAnth)
	  stop
    endif

    if(NumNat==0.and.NatEmisNum/=0) then
	  print '(/,"STOP: No emission sources is defined for natural emissions of forms:",<NatEmisNum>(1x,a),/)',&
                            &(trim(SubsID(Subs))//trim(FormID(Atm,Subs,NatEmisInd(Ind))), Ind=1, NatEmisNum)
	  stop
    elseif(NumNat/=0.and.NatEmisNum==0) then
	  print '(/,"STOP: No natural emission is defined for sources:",<NumNat>(1x,a),/)',&
                            &(trim(NaturID(Src)), Src=1, NumNat)
	  stop
    endif

#if (REGTYPE==2)
    if(NumBnd==0.and.AtmBndNum/=0) then
	  print '(/,"STOP: No sources is defined for boundary concentrations of forms:",<AtmBndNum>(1x,a),/)',&
                            &(trim(SubsID(Subs))//trim(FormID(Atm,Subs,AtmBndInd(Ind))), Ind=1, AtmBndNum)
	  stop
    elseif(NumBnd/=0.and.AtmBndNum==0) then
	  print '(/,"STOP: No boundary concentrations is defined for sources:",<NumBnd>(1x,a),/)',&
                            &(trim(BoundID(Src)), Src=1, NumBnd)
	  stop
    endif
#endif
#endif

end subroutine Atm_POP_Props


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading/defining POP input fields
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_POP_Input(Prd)

    character(*) Prd
    real Pol, Koa, Kad, Kab, TAir(Atm_KMax), OHTrans
    real(8) RhoAir, Psur
    integer NPOP, i, j, k,t, NTHREADS
    real Lw,Sw,Lwmax,c_OH

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_POP_Input ', Prd
#endif

    select case (Prd)

    case ('Initial')           ! After POP properties are read

        kDegrG = 0.
        kDegrP = 0.
        do NPOP = 1, gSubsNum(GroupInd(POP))
            WRatioGas0(NPOP) = 2./((HT0(NPOP)/(Runiv * 273.15)) * exp(-HT(NPOP) * ((1./273.15)-(1./T0))))
        end do

    case ('Daily')
       select case (ReactSource)
           case('GEOS_Chem_2x2.5')
               call Atm_POP_ReadReactDaily_GEOSChem
           case('GEOS_Chem_4x5')
               call Atm_POP_ReadReactDaily_GEOSChem
           case('MOZART')
               call Atm_POP_ReadReactDaily_MOZART
           case default
               print '(/,"STOP: Invalid source of reactant data ''",a,"''",/)', trim(ReactSource)
               stop
       end select
    case ('Monthly')
        do NPOP = 1, gSubsNum(GroupInd(POP))
            if(.not. (TD(NPOP))) then
                select case (Month)
                case (1,2,12)
                    kDegrG(:,:,:,:,NPOP) = CDWinter(NPOP)
                case (3,4,5,9,10,11)
                    kDegrG(:,:,:,:,NPOP) = CDSpring(NPOP)
                case (6,7,8)
                    kDegrG(:,:,:,:,NPOP) = CDSummer(NPOP)
                end select
            end if
        end do
    case ('Yearly')

    end select

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_POP_Input ', Prd
#endif
end subroutine Atm_POP_Input

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine describing POP partitioning in the atmosphere
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_POP_partitioning(Subs)

    integer, intent(in)   ::  Subs
    real*8 ConcTot, MV, t1, t2, t3
    real*8 ConcGas, ConcPart              ! Matrix calculations 24.12.2014
    integer i, j, k,Ind,GasInd,PartInd

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_POP_partitioning'
#endif

!    if (SubsGroup(Subs) == POP) then
    do k = 1, Atm_KMax
        do j = JMin, Jmax
            do i = MinI(j), MaxI(j)
                MV = MeshVolume(i,j,k)
                ConcGas = Atm_Conc(i,j,k,sFormInd(Atm,Subs,1))   ! Matrix calculations 24.12.2014
                ConcPart = Atm_Conc(i,j,k,sFormInd(Atm,Subs,2))  ! Matrix calculations 24.12.2014
                ConcTot = ConcGas + ConcPart                     ! Matrix calculations 24.12.2014
#ifdef DEBUG_MODE
                if(i == 20 .and. j == 30 .and. k == 1) then
                    print *, 'CG, CP: ', ConcGas, ConcPart, kpt(i,j,k,Period,gSubsGroupInd(Subs))
                end if
#endif

!           Mass of POP partitioned into part phase
                MassAtmPartit(sFormInd(Atm,Subs,2)) = MassAtmPartit(sFormInd(Atm,Subs,2))&
                &+( ConcTot * kpt(i,j,k,Period,gSubsGroupInd(Subs))-Atm_Conc(i,j,k,sFormInd(Atm,Subs,2)))*MV
!           Air concentration in part phase at current step
                Atm_Conc(i,j,k,sFormInd(Atm,Subs,2)) = ConcTot * kpt(i,j,k,Period,gSubsGroupInd(Subs))
!           Mass of POP partitioned into gas phase
                MassAtmPartit(sFormInd(Atm,Subs,1)) = MassAtmPartit(sFormInd(Atm,Subs,1))&
                &+ (ConcTot - Atm_Conc(i,j,k,sFormInd(Atm,Subs,2))-Atm_Conc(i,j,k,sFormInd(Atm,Subs,1)))*MV
!           Air concentration in gas phase at current step
                Atm_Conc(i,j,k,sFormInd(Atm,Subs,1)) = ConcTot - Atm_Conc(i,j,k,sFormInd(Atm,Subs,2))
#if RTYPE==2
!           Matrix calculations
                if(ConcTot /= 0.) then

                Atm_Contrib(i,j,k,sFormInd(Atm,Subs,1),1:NumSrc) = (Atm_Contrib(i,j,k,sFormInd(Atm,Subs,1),1:NumSrc)*&
                    & ConcGas + Atm_Contrib(i,j,k,sFormInd(Atm,Subs,2),1:NumSrc) * ConcPart) / ConcTot
                Atm_Contrib(i,j,k,sFormInd(Atm,Subs,2),1:NumSrc) = Atm_Contrib(i,j,k,sFormInd(Atm,Subs,1),1:NumSrc)

                end if

#endif
 		enddo	!k
	  enddo		!j
	enddo		!i

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_POP_partitioning'
#endif
end subroutine Atm_POP_partitioning

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine describing POP degradation in the atmosphere
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_POP_degradation(Subs)

    integer Subs
    real dTime
    real*8 DMass, AConc11,AConc21,AConc12,AConc22
    integer i, j, k, Form
    real*8 APM1(NumForm(Atm)),APM2(NumForm(Atm))

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_POP_degradation'
#endif

! Mass of POP in air before degradation
    APM1=0.
    do Form=1, NumForm(Atm)
        do k=1, Atm_Kmax
            do j=Jmin, Jmax
                do i=minI(j), maxI(j)
                    APM1(Form)=APM1(Form)+Atm_Conc(i,j,k,Form)*MeshVolume(i,j,k)
                enddo
            enddo
        enddo
    enddo

    do k = 1, Atm_KMax
        do j = JMin, Jmax
            do i =MinI(j), MaxI(j)
                AConc11 = Atm_Conc(i,j,k,sFormInd(Atm,Subs,1))                            ! Gas phase
                AConc12 = Atm_Conc(i,j,k,sFormInd(Atm,Subs,2))                            ! Part phase
                AConc21 = AConc11 * exp(-kDegrG(i,j,k,Period,gSubsGroupInd(Subs))*Tstep(Atm))
                AConc22 = AConc12 * exp(-kDegrP(i,j,k,Period,gSubsGroupInd(Subs))*Tstep(Atm))
                Atm_Conc(i,j,k,sFormInd(Atm,Subs,1)) = AConc21
                Atm_Conc(i,j,k,sFormInd(Atm,Subs,2)) = AConc22
#ifdef DEGR_OUT
                AtmDegrMonthGas(i,j,k) = AtmDegrMonthGas(i,j,k) + (AConc11 - AConc21)*MeshVolume(i,j,k)
                AtmDegrMonthPart(i,j,k) = AtmDegrMonthPart(i,j,k) + (AConc12 - AConc22)*MeshVolume(i,j,k)
#endif
            enddo	!k
        enddo		!j
    enddo           !i


! Mass of POP in air after degradation
    APM2=0.
    do Form=1, NumForm(Atm)
        do k=1, Atm_Kmax
            do j=Jmin, Jmax
                do i=minI(j), maxI(j)
                    APM2(Form)=APM2(Form)+Atm_Conc(i,j,k,Form)*MeshVolume(i,j,k)
                enddo
            enddo
        enddo
    enddo

! Mass counters of POP degraded in gas and part phases
    MassAtmDegr(sFormInd(Atm,Subs,1)) = MassAtmDegr(sFormInd(Atm,Subs,1)) + APM1(1) - APM2(1)
    MassAtmDegr(sFormInd(Atm,Subs,2)) = MassAtmDegr(sFormInd(Atm,Subs,2)) + APM1(2) - APM2(2)

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_POP_degradation'
#endif
end subroutine Atm_POP_degradation


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating wet deposition velocities (gas-phase)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_CalcPOPWDVelGas(NSubs)

  integer NSubs
  integer i, j, k, ind
  real(8) Vel1, Vel2, TAir
  real(8) rc, rnc, prcp

#ifdef DEBUG_MODE
    print *, '<- Entering Atm_CalcPOPWDVelGas'
#endif

    do k = 1, Atm_KMax
        do j = JMin, Jmax
            do i =MinI(j), MaxI(j)
                Tair = TempAir(i, j, k, Period,toDay)
                rc = PrecRainConv(i,j,k,Period)
                rnc = PrecRainStrat(i,j,k,Period)
                ind = gSubsGroupInd(NSubs)
                prcp = rc + rnc
                if(.not. (WashExp(ind))) then
                    Vel1 = 1./((HT0(ind)/(Runiv * Tair)) * exp(-HT(ind) * (1./Tair-1./T0))) * prcp
                    Atm_WDVelGas(i,j,k,NSubs) = Vel1
                    if ((Tair < 273.15) .and. (Tair > 263.15)) then
                        Vel2 = WRatioGas0(ind) * exp(0.72 * pT(ind) * (1./Tair - 1/273.15)) * prcp
                        Atm_WDVelGas(i,j,k,NSubs) = ((Tair - 263.15) * Vel1 + (273.15 - Tair) * Vel2) / 10.
                    elseif (Tair <= 263.15) then
                        Atm_WDVelGas(i,j,k,NSubs) = WRatioGas0(ind) * exp(0.72 * pT(ind) * (1./Tair - 1/273.15)) * prcp
                    end if
                else
                    if ((Tair < 273.15) .and. (Tair > 263.15)) then
                        Atm_WDVelGas(i,j,k,NSubs) = (2. - (Tair - 263.15) / 10.) * WRatioGas0(ind) * prcp
                    elseif (Tair <= 263.15) then
                        Atm_WDVelGas(i,j,k,NSubs) = 2. * WRatioGas0(ind) * prcp
                    else
                        Atm_WDVelGas(i,j,k,NSubs) = WRatioGas0(ind) * prcp
                    end if
                end if
            enddo	!i
        enddo		!j
    enddo		!k

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_CalcPOPWDVelGas'
#endif
end subroutine Atm_CalcPOPWDVelGas


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating wet deposition coefficients
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_CalcPOPWDCoeff(NSubs)

    integer NSubs
    integer i, j, k
    real Iprec(Atm_Kmax), Ip
    real Lw, Fw, LFw

    real Tair,WRatioPrt,Koef,buf,Teff,dz(Atm_Kmax)
    real HumidCoef
    real(8)Conc(Atm_Kmax)
    real Keff,beta

#ifdef DEBUG_MODE
    print *, '-> Entering Atm_CalcPOPWDCoeff... ', NSubs
#endif

    Atm_WDCoeff = 0.

    do k = 1, Atm_Kmax
     do j = JMin, Jmax
      do i = MinI(j), MaxI(j)
          Iprec(k)=PrecRainConv(i,j,k,Period)+PrecRainStrat(i,j,k,Period)
	  Ip=Iprec(k)
	  Lw=LiqCont(i,j,k,Period)
	  Fw=FrozCont(i,j,k,Period)

          LFw = Lw + Fw
          if(LFw > Lw0) then					! In-cloud scavenging
	    if(Ip > 0.) then
	      Keff = LFw / (LFw + Aeff)
              Atm_WDCoeff (i,j,k,NSubs) = Ain * (Ip*3.6e6)**Bin * Keff * (1.-Asol)/(1.-Asol * Keff)
	    endif
	  else							! Below-cloud scavenging
	    if(Ip>0.) then
              Atm_WDCoeff (i,j,k,NSubs) = Abelow*(Ip*3.6e6)**Bbelow
	    endif
	  endif
        end do ! i
      end do   ! j
    end do     ! k

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_CalcPOPWDCoef'
#endif
end subroutine Atm_CalcPOPWDCoeff


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating dry particle phase deposition velocity
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_CalcPOPDDVelPart(NSubs)

    integer i, j, j1, NSubs
    integer L, t, Surf, Ind, s, it
    character*20 SurfG
    real Pxx, Pa, RT, Rho, Tair, Psat, Qv, RH, Mu, Lmo0, Lmo(0:Soil_NumType), Almo, Zr, Ui, Uj, Uref, Kcond, Dh, SnowH, Tsurf
    real DpW, Kn, Cunn, Vg, DiffW, DpS, VgS, DiffS, VgF, DiffF, Ma, DiffG, Cpw, Dfor, Hfor, ZoM, ZoH, Ux0, Ux
    real PrecN(Imin:Imax,Jmin:Jmax), PrecB(Imin:Imax,Jmin:Jmax), Lcover(0:Soil_NumType)
    real St, Sc, Uh, Eb, Ein, Eim, Reb, Eff, Fbrok, Pr, Tcoef, Lw, VdFog
    real Ra, Rs, Rs1, Rs2
    logical WetGrass, WetForest

#ifdef DEBUG_MODE
      print *, '-> Entering Atm_CalcPOPDDVelPart...'
#endif
    POPVdAero=0.
    do j=Jmin, Jmax
        do i=minI(j), maxI(j)

! Calculation of air parameters
            Lcover(0:Soil_NumType)=Soil_Frac(i,j,0:Soil_NumType)		! Land cover types
            s=Season(i,j,Month)							! Seasons (1-5)
            Pxx=Px(i,j,Period,toDay)						! Surface pressure Px=Ps-Pt
            Pa=Sigma(1)*Pxx+Ptop						! Pressure at the lowest sigma level
            RT=RTemp(i,j,1,Period,toDay)					! RT of air at the lowest sigma level
            Rho=Pa/RT								! Air density at the lowest sigma level
            Lw=LiqCont(i,j,1,Period)+FrozCont(i,j,1,Period)			! Liquid or frozen water content
            Tair=TempAir(i,j,1,Period,toDay)					! Temperature at the lowest sigma level
            Tsurf=TempSurf(i,j,Period)

            if(Tair>273.15) then
                Psat=6.112*exp(17.67*(Tair-273.15)/(Tair-29.65))*100.		! Saturation vapour pressure (t>0) [Pa]
            else
                Psat=6.11*exp(22.514-6150./Tair)*100.				! Saturation vapour pressure (t<0) [Pa]
            endif
            SnowH=SnowDepth(i,j,Period)
            Qv=HumidAir(i,j,1,Period,toDay)					! Vapour mixing ratio
            RH=min(Qv*(1.+HumCoef)/(1.+Qv*(1.+HumCoef))*Pa/Psat,1.)		! Relative humidity (saturation ratio)
            Mu=Nu*Rho								! Dynamic viscosity
            Zr=RT/Ggrav*log((Pxx+Ptop)/(Sigma(1)*Pxx+Ptop))			! Height of the lowest sigma level
            Ui=(Uwind(i,j,1,Period,toDay)+Uwind(i-1,j,1,Period,toDay))/2.
            if(j == Jmin) then
                Uj=Vwind(i,j,1,Period,toDay)
            else
                j1=j-1
                Uj=(Vwind(i,j,1,Period,toDay)+&
                        &Vwind(iS(i,j,1),j1,1,Period,toDay))/2.                ! Corrected 31.10.2018
            end if
            Uref=sqrt(Ui*Ui+Uj*Uj)						! Absolute value of vind velosity
            Ux0=Ufric(i,j,Period,toDay)						! Mean friction velocity
            Lmo0=MOLength(i,j,Period,toDay)					! Mean Monin-Obuhov length
            Lmo(0:Soil_NumType)=Lmo0						! Surface dependant MO length
            Almo=Lmo0/Ux0/Ux0/Ux0						! MO length coefficient
            Kcond=0.023807+7.1128e-5*(Tair-273.15)				! Thermal conductivity of dry air
            Dh=Kcond/Rho/Cpd							! Molecular thermal diffusion coefficient
            Ma=Md/(1.+HumCoef*Qv/(1.+Qv))					! Malecular weight of total air
            DiffG=3./8./Nav/Dm/Dm/Rho*&
                    &sqrt(Runiv*Tair*Ma*(Ma+MHgCL2)/MHgCL2/2./PiNum)            ! Diffusion coefficient of RGM
            Cpw=Cpd*(1.+CpCoef*Qv/(1.+Qv))					! Specific heat of moist air
            Pr=Mu*Cpw/Kcond							! Prandtl number

            DpW=WetDiameter(Dp,RH,1.)						! Wet diameter of particles
            Kn=4.*Nu/sqrt(8.*kB*Tair*Nav/PiNum/Md)/DpW				! Knudsen number
            Cunn=1+Kn*(1.249+0.42*exp(-0.87/Kn))				! Cunningham correction
            Vg=DpW*DpW*RhoP*Ggrav/18./Mu*Cunn					! Gravitational settling velosity
            DiffW=kB*Tair/3./PiNum/DpW/Mu*Cunn					! Particle diffusion coefficient in moist air

!		  if(Soil_Frac(i,j,0)>0)	then			! Commented 19.01.2018, caused DiffS=0 and division by zero
! In the original version for HMs the following lines related to water surfaces
            DpS=WetDiameter(Dp,0.98,1.)
            Kn=4.*Nu/sqrt(8.*kB*Tair*Nav/PiNum/Md)/DpS                          ! Particles in water surface layer
            Cunn=1+Kn*(1.249+0.42*exp(-0.87/Kn))
            VgS=DpS*DpS*RhoP*Ggrav/18./Mu*Cunn
            DiffS=kB*Tair/3./PiNum/DpS/Mu*Cunn
!		  endif

! Calculation of surface wetness (for grass and forest)
            PrecN(i,j)=PrecRainConv(i,j,1,Period)+PrecRainStrat(i,j,1,Period)	! Current precipitation rate
            if(PrecN(i,j)>0.) then
                WetGrass=.true.                                                 ! Wet grass surface
            else
                WetGrass=.false.						! Dry grass surface
            endif
            if(PrecN(i,j)>0..or.PrecB(i,j)>0.) then
                WetForest=.true.						! Wet forest surface
            else
                WetForest=.false.						! Dry forest surface
            endif
            PrecB(i,j)=PrecN(i,j)						! Previous precipitation rate

            do L=0,Soil_NumType
                Hfor=Soil_Height(L)						! Canopy height
                Dfor=Soil_Disp(L,s)						! Displacement height

                if(L==0) then                                                   !Water
                    ZoM=0.016*Ux0*Ux0/Ggrav+Nu/9.1*Ux0				! Roughness
                    Ux=Karman*Uref/IphiM(Lmo(L),Zr,ZoM)				! Friction velocity
                    Ux = max(Ux, 0.1)	! 09.04.2008 Ilia
                    ZoH=ZoM*exp(-Karman*(13.6*Pr**(2./3.)-12.))			! Energy roughness
                else
                    ZoM=Soil_Z0(L,s)						! Momentum roughness
                    do it=1, 5
                        Ux=max(Karman*Uref/IphiM(Lmo(L),Zr-Dfor,ZoM),0.1)	! Friction velocity
                        Lmo(L)=Almo*Ux*Ux*Ux
                    enddo
                    ZoH=0.135*ZoM						! Energy roughness
                endif

                Ra=IphiH(Lmo(L),Zr-Dfor,ZoH)/Karman/Ux				! Aerodynamic resistance

                if((SnowH>0.1.and.(trim(adjustL(SoilType(L)))  == 'Grass'.or.trim(adjustL(SoilType(L)))  == 'Arable')).or.&
                  &(Tsurf<270..and.trim(adjustL(SoilType(L)))  == 'Water')) then  !  Changed
                    SurfG='Bare'
		else
                    SurfG=trim(adjustL(SoilType(L)))
		endif
		selectcase(trim(adjustL(SurfG)))
		  case('Decid','Conif')
                    St=Ux*Vg/Ggrav/1.e-3
                    Sc=Nu/DiffW
                    Uh=Ux/Karman*IphiM(Lmo(L),Hfor-Dfor,ZoM)
                    if(WetForest.eqv..true.) then
                        Reb=1.
                    else
                        Reb=exp(-Cforest*St**0.25)
                    endif
                    Eb=Sc**(-2./3.)
                    Ein=(0.01*DpW/(DpW+2.e-5)+0.99*DpW/(DpW+2.e-3))
                    Eim=(St/(Aforest+St))**0.5
                    Eff=Bforest*(Eb+Ein+Eim)*Reb
                    Rs=Uh/Eff/Ux/Ux
                    POPVdAero(i,j,L,NSubs)=1./(Ra+Rs)+Vg

                  case('Grass','Arable','Scrabs')
                    St=Ux*Vg/Ggrav/1.e-3
                    Sc=Nu/DiffW
                    Uh=Ux/Karman*IphiM(Lmo(L),Hfor-Dfor,ZoM)
                    if(WetGrass.eqv..true.) then
                        Reb=1.
                    else
                        Reb=exp(-Cgrass*St**0.25)
                    endif
                    Eb=Sc**(-2./3.)
                    Ein=(0.01*DpW/(DpW+2.e-5)+0.99*DpW/(DpW+2.e-3))
                    Eim=(St/(Agrass+St))**0.5
                    Eff=Bgrass*(Eb+Ein+Eim)*Reb
                    if(Lmo(Surf)>=0.) then
                        Rs=Uh/Eff/Ux/Ux
                    else
                        Rs=Uh/Eff/Ux/Ux/(1.+(-Dgrass/Lmo(L))**(2./3.))
                    endif
                    POPVdAero(i,j,L,NSubs)=1./(Ra+Rs)+Vg

		  case('Water')
                    St=Ux*Ux*VgS/Ggrav/Nu
                    Sc=Nu/DiffS
                    Uh=Ux/Karman*IphiM(Lmo(L),10.,ZoM)
                    if (Uh.gt.100.) then
                        Fbrok=1.
                    else
                        Fbrok=min(1.7e-6*Uh**3.75,1.)
                    end if
                    Rs1=Karman*Uh/Ux/Ux/(10.**(-3./St)+Sc**(-0.5))
                    Rs2=Rbrok ! 10. s/m
                    Rs=1./((1.-Fbrok)/Rs1+Fbrok/Rs2)
                    POPVdAero(i,j,L,NSubs)=(1.+Ra*Vg)*(1.+Rs*VgS)/(Ra+Rs+Ra*Rs*VgS)

        	  case('Bare','Glacier')
                    St=Ux*Ux*Vg/Ggrav/Nu
                    Sc=Nu/DiffW
                    Rs=Karman*Uref/Ux/Ux/(10.**(-3./St)+Sc**(-2./3.))
                    POPVdAero(i,j,L,NSubs)=1./(Ra+Rs+Ra*Rs*Vg)+Vg

		  case('Urban')
                    St=Ux*Ux*Vg/Ggrav/Nu
                    Sc=Nu/DiffW
                    Rs=Karman*Uref/Ux/Ux/(St*St/(400.+St*St)+Sc**(-2./3.))
                    POPVdAero(i,j,L,NSubs)=1./(Ra+Rs+Ra*Rs*Vg)+Vg
		endselect
	  enddo			! end L
	enddo			! end i
  enddo				! end j


#ifdef DEBUG_MODE
      print *, '-> Exit Atm_CalcPOPDDVelPart'
#endif
end subroutine Atm_CalcPOPDDVelPart


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading reactants and aerosol data fields from GEOS-Chem output data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_POP_ReadReactDaily_GEOSChem

    integer status, ncid_in, var_id_in, start3(4), count3(4), t
    real qreact_buf( Imax, Jmax, Atm_Kmax, 1:4)
    data start3 /1,1,1,1/
    data count3 /Imax, Jmax, Atm_Kmax, 4/

    integer     ii, jj, AN, RN
    integer     FileStat, i, j, k, t, Xscal, daF(2), mnF(2), yrF(2), nDay, flag, yrCur
	real        qReact(4), Aver(Imin:Imax)
!   real        qOC(4), qO3(4), qAS(4), qOH(4)
    real        splRow(8), dFdX(8), d2FdX(8), dFdXleft

! Check for climatic run
    if(climReactRun) then
      yrCur=ClimReactYear
    else
      yrCur=Year
    endif

    write(YearNum,'(i4)') yrCur
    write(DayNum,'(i2.2)') Day

#ifdef DEBUG_MODE_ATM
    print *, '>- Entering Atm_POP_ReadReaictDaily_GEOSChem'
#endif
!------------------------------------------------------------------------------------------
! Read reactants and aerosol data

    daF(1)=Day
    mnF(1)=Month
    yrF(1:2)=yrCur
    if(Day<MonthDays(Month)) then
          daF(2)=Day+1
          mnF(2)=Month
    elseif(Month<12) then
          daF(2)=1
          mnF(2)=Month+1
!    else
!          daF(2)=1
!          mnF(2)=1
    else
        daF(2)=1
        mnF(2)=1
	if(yrCur<FinDate(yr)) then
	    yrF(2) = yrCur+1
	else
	    yrF(2) = yrCur
	end if
    endif

    do nDay=1, 2
! Reading aerosol 3D data
	write(YearNum,'(i4)') yrF(nDay)
        do AN = 1, NmbAero                     ! NmbAero should be defined from Properties file
            write(fileName,'(a,i4,i2.2,i2.2,a4)') trim(AName(AN)), yrF(nDay), mnF(nDay), daF(nDay), '.nc4'
                  ! AName should be read from Properties file: AName = 'OC_' for OC, AName = 'ADens_' for aerosol surface, etc
            fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read')
!            if(FileStat>0) then
!                    print '(/,"STOP: Cannot open data file ''",a,"''",/)', trim(fullName)
!                    stop
!            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            status=nf90_open(trim(fullName), nf90_nowrite, ncid_in)
            if(status/=nf90_noerr) then
                print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
                stop
            endif
            call checkNC(nf90_inq_varid(ncid_in, AName(AN)(1:len_trim(AName(AN))-1), var_id_in), 8)
            call checkNC(nf90_get_var(ncid_in, var_id_in, qreact_buf, start3, count3), 9)

            do k=1, Atm_Kmax
                do j=Jmin, Jmax
                    do i=Imin, Imax
                        do t = 1, 4
                            Aerofield(AN,i,j,k,t,nDay)=qReact_buf( i, j, k, t )   ! AeroField should be defined instead of OCField, ...
                        enddo
                    enddo
                enddo
            enddo
            where(Aerofield <= 0.) Aerofield = Zero
            call checkNC(nf90_close(ncid_in), 10)
        end do
    enddo  ! nDay

!------------------------------------------------------------------------------------------
! Grid aggregation of the data read
    do AN = 1, NmbAero                 ! For all aerosol species
        do j = Jmin, Jmax
            if(maxI(j) == 1) cycle
            Xscal = Imax/maxI(j)
            if(Xscal == 1) cycle
            do nDay = 1, 2
                do t = 1, 4
                    do k = 1, Atm_Kmax
                        Aver(Imin:Imax)=Aerofield(AN,Imin:Imax,j,k,t,nDay)
                        call GridAggreg(j,Xscal,Aver,1)
                        Aerofield(AN,minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
                    enddo ! k
                enddo ! t
            enddo ! nDay
        enddo ! j

        do j=Jmin, Jmax
            do i=minI(j), maxI(j)
                do k=1, Atm_Kmax
                    splRow(1:4)=Aerofield(AN,i,j,k,1:4,toDay)
                    splRow(5:8)=Aerofield(AN,i,j,k,1:4,toMor)
                    if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
                        flag=0
                        dFdXleft=0.
                    else
                        flag=1
                        dFdXleft=dAerofield(AN,i,j,k,5)
                    endif
                    call SplineParams(8,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
                    dAerofield(AN,i,j,k,1:5)  = dFdX(1:5)     ! dAeroField should be defined instead of dOCField, ...
                    d2Aerofield(AN,i,j,k,1:5) = d2FdX(1:5)    ! d2AeroField should be defined instead of d2OCField, ...
                enddo ! k
            enddo ! i
        enddo ! j
	end do   ! AN
! Reading reactants 3D data
    do nDay=1, 2
	write(YearNum,'(i4)') yrF(nDay)
        do RN = 1, NmbReacts                     ! NmbReacts should be defined from Properties file
																   ! If ROI_T NmbReacts should be put to 1 and RName(1) - to O3_ (corresponding to ozone file)
																   ! Correct(1) can be set to correct ozone concentations
            write(fileName,'(a,i4,i2.2,i2.2,a4)') trim(RName(RN)), yrF(nDay), mnF(nDay), daF(nDay), '.nc4'
                  ! RName should be read from Properties file: RName = 'OH_' for OH,  etc
            fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!            open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read')
!            if(FileStat>0) then
!                    print '(/,"STOP: Cannot open data file ''",a,"''",/)', trim(fullName)
!                    stop
!            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            status=nf90_open(trim(fullName), nf90_nowrite, ncid_in)
            if(status/=nf90_noerr) then
                print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
                stop
            endif
            call checkNC(nf90_inq_varid(ncid_in, RName(RN)(1:len_trim(RName(RN))-1), var_id_in), 8)
            call checkNC(nf90_get_var(ncid_in, var_id_in, qreact_buf, start3, count3), 9)

            do k=1, Atm_Kmax
                do j=Jmin, Jmax
                    do i=Imin, Imax
                        do t = 1, 4
                            Reactfield(RN,i,j,k,t,nDay)=qReact_buf( i, j, k, t )   ! AeroField should be defined instead of OCField, ...
                        enddo
                    enddo
                enddo
            enddo
            where(Aerofield <= 0.) Aerofield = Zero
            call checkNC(nf90_close(ncid_in), 10)

            !####
!            do k=1, Atm_Kmax
!                    do j=Jmin, Jmax
!                    do i=Imin, Imax
!                        read(10) (qReact(t), t=1, 4)          ! Different units for different reactants; correction - specific coefficients (see below)
!                        Reactfield(RN,i,j,k,1:4,nDay)=qReact(1:4)   ! ReactField should be defined instead of OHField, ...
!                    enddo
!                enddo
!            enddo
!            where(Reactfield <= 0.) Reactfield = Zero
!            close(10)
        end do
    enddo  ! nDay

!------------------------------------------------------------------------------------------
! Grid aggregation of the data read
    do RN = 1, NmbReacts                 ! For all reactants
        do j = Jmin, Jmax
            if(maxI(j) == 1) cycle
            Xscal = Imax/maxI(j)
            if(Xscal == 1) cycle
            do nDay = 1, 2
                do t = 1, 4
                    do k = 1, Atm_Kmax
                        Aver(Imin:Imax)=Reactfield(RN,Imin:Imax,j,k,t,nDay)
                        call GridAggreg(j,Xscal,Aver,1)
                        Reactfield(RN,minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
                    enddo ! k
                enddo ! t
            enddo ! nDay
        enddo ! j

        do j=Jmin, Jmax
            do i=minI(j), maxI(j)
                do k=1, Atm_Kmax
                    splRow(1:4)=Reactfield(RN,i,j,k,1:4,toDay)
                    splRow(5:8)=Reactfield(RN,i,j,k,1:4,toMor)
                    if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
                        flag=0
                        dFdXleft=0.
                    else
                        flag=1
                        dFdXleft=dReactfield(RN,i,j,k,5)
                    endif
                    call SplineParams(8,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
                    dReactfield(RN,i,j,k,1:5)  = dFdX(1:5)     ! dReactField should be defined instead of dOCField, ...
                    d2Reactfield(RN,i,j,k,1:5) = d2FdX(1:5)    ! d2ReactField should be defined instead of d2OCField, ...
                enddo ! k
            enddo ! i
        enddo ! j
    end do       ! RN

#ifdef DEBUG_MODE_ATM
    print *, '<- Exit Atm_POP_ReadReactDaily_GEOSChem'
#endif
end subroutine Atm_POP_ReadReactDaily_GEOSChem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading reactants and aerosol data fields from MOZART output data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Atm_POP_ReadReactDaily_MOZART

    integer     ii, jj, AN, RN
    integer     FileStat, i, j, k, t, Xscal, daF(2), mnF(2), yrF(2), nDay, flag, yrCur
	real        qReact(4), Aver(Imin:Imax)
!    real        qOC(4), qO3(4), qAS(4), qOH(4)
    real        splRow(8), dFdX(8), d2FdX(8), dFdXleft

! Check for climatic run
    if(climReactRun) then
      yrCur=ClimReactYear
    else
      yrCur=Year
    endif

    write(YearNum,'(i4)') yrCur
    write(DayNum,'(i2.2)') Day

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_POP_ReadReactDaily_MOZART'
#endif
!------------------------------------------------------------------------------------------
! Read reactants and aerosol data

    daF(1)=Day
    mnF(1)=Month
    yrF(1:2)=yrCur
    if(Day<MonthDays(Month)) then
          daF(2)=Day+1
          mnF(2)=Month
    elseif(Month<12) then
          daF(2)=1
          mnF(2)=Month+1
!    else
!          daF(2)=1
!          mnF(2)=1
    else
        daF(2)=1
        mnF(2)=1
	if(yrCur<FinDate(yr)) then
	    yrF(2) = yrCur+1
	else
	    yrF(2) = yrCur
	end if
    endif

    do nDay=1, 2
! Reading aerosol 3D data
	write(YearNum,'(i4)') yrF(nDay)
        do AN = 1, NmbAero                     ! NmbAero should be defined from Properties file
            write(fileName,'(a,i4,i2.2,i2.2,a4)') trim(AName(AN)), yrF(nDay), mnF(nDay), daF(nDay), '.bin'
                  ! AName should be read from Properties file: AName = 'OC_' for OC, AName = 'ADens_' for aerosol surface, etc
            fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
            open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read')
            if(FileStat>0) then
                    print '(/,"STOP: Cannot open data file ''",a,"''",/)', trim(fullName)
                    stop
            endif
            do k=1, Atm_Kmax
                do j=Jmin, Jmax
                    do i=Imin, Imax
                        read(10) (qReact(t), t=1, 4)          ! Different units for different aerosol species; correction - specific coefficients (see below)
                        Aerofield(AN,i,j,k,1:4,nDay)=qReact(1:4)   ! AeroField should be defined instead of OCField, ...
                    enddo
                enddo
            enddo
            where(Aerofield <= 0.) Aerofield = Zero
            close(10)
        end do
    enddo  ! nDay

!------------------------------------------------------------------------------------------
! Grid aggregation of the data read
    do AN = 1, NmbAero                 ! For all aerosol species
        do j = Jmin, Jmax
            if(maxI(j) == 1) cycle
            Xscal = Imax/maxI(j)
            if(Xscal == 1) cycle
            do nDay = 1, 2
                do t = 1, 4
                    do k = 1, Atm_Kmax
                        Aver(Imin:Imax)=Aerofield(AN,Imin:Imax,j,k,t,nDay)
                        call GridAggreg(j,Xscal,Aver,1)
                        Aerofield(AN,minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
                    enddo ! k
                enddo ! t
            enddo ! nDay
        enddo ! j

        do j=Jmin, Jmax
            do i=minI(j), maxI(j)
                do k=1, Atm_Kmax
                    splRow(1:4)=Aerofield(AN,i,j,k,1:4,toDay)
                    splRow(5:8)=Aerofield(AN,i,j,k,1:4,toMor)
                    if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
                        flag=0
                        dFdXleft=0.
                    else
                        flag=1
                        dFdXleft=dAerofield(AN,i,j,k,5)
                    endif
                    call SplineParams(8,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
                    dAerofield(AN,i,j,k,1:5)  = dFdX(1:5)     ! dAeroField should be defined instead of dOCField, ...
                    d2Aerofield(AN,i,j,k,1:5) = d2FdX(1:5)    ! d2AeroField should be defined instead of d2OCField, ...
                enddo ! k
            enddo ! i
        enddo ! j
	end do   ! AN
! Reading reactants 3D data
    do nDay=1, 2
	write(YearNum,'(i4)') yrF(nDay)
        do RN = 1, NmbReacts                     ! NmbReacts should be defined from Properties file
																   ! If ROI_T NmbReacts should be put to 1 and RName(1) - to O3_ (corresponding to ozone file)
																   ! Correct(1) can be set to correct ozone concentations
            write(fileName,'(a,i4,i2.2,i2.2,a4)') trim(RName(RN)), yrF(nDay), mnF(nDay), daF(nDay), '.bin'
                  ! RName should be read from Properties file: RName = 'OH_' for OH,  etc
            fullName=trim(ReactPath1)//trim(GridCode)//'/'//YearNum//'/'//trim(fileName)
            open(10, file=fullName, form='unformatted', access="stream", status='old', iostat=FileStat, action='read')
            if(FileStat>0) then
                    print '(/,"STOP: Cannot open data file ''",a,"''",/)', trim(fullName)
                    stop
            endif
                do k=1, Atm_Kmax
        do j=Jmin, Jmax
                    do i=Imin, Imax
                        read(10) (qReact(t), t=1, 4)          ! Different units for different reactants; correction - specific coefficients (see below)
                        Reactfield(RN,i,j,k,1:4,nDay)=qReact(1:4)   ! ReactField should be defined instead of OHField, ...
                    enddo
                enddo
            enddo
            where(Reactfield <= 0.) Reactfield = Zero
            close(10)
        end do
    enddo  ! nDay

!------------------------------------------------------------------------------------------
! Grid aggregation of the data read
    do RN = 1, NmbReacts                 ! For all reactants
        do j = Jmin, Jmax
            if(maxI(j) == 1) cycle
            Xscal = Imax/maxI(j)
            if(Xscal == 1) cycle
            do nDay = 1, 2
                do t = 1, 4
                    do k = 1, Atm_Kmax
                        Aver(Imin:Imax)=Reactfield(RN,Imin:Imax,j,k,t,nDay)
                        call GridAggreg(j,Xscal,Aver,1)
                        Reactfield(RN,minI(j):maxI(j),j,k,t,nDay)=Aver(minI(j):maxI(j))
                    enddo ! k
                enddo ! t
            enddo ! nDay
        enddo ! j

        do j=Jmin, Jmax
            do i=minI(j), maxI(j)
                do k=1, Atm_Kmax
                    splRow(1:4)=Reactfield(RN,i,j,k,1:4,toDay)
                    splRow(5:8)=Reactfield(RN,i,j,k,1:4,toMor)
                    if(yrCur==BegDate(yr).and.Month==BegDate(mn).and.Day==BegDate(da)) then
                        flag=0
                        dFdXleft=0.
                    else
                        flag=1
                        dFdXleft=dReactfield(RN,i,j,k,5)
                    endif
                    call SplineParams(8,timePer,splRow,flag,dFdXleft,0,0.,dFdX,d2FdX)
                    dReactfield(RN,i,j,k,1:5)  = dFdX(1:5)     ! dReactField should be defined instead of dOCField, ...
                    d2Reactfield(RN,i,j,k,1:5) = d2FdX(1:5)    ! d2ReactField should be defined instead of d2OCField, ...
                enddo ! k
            enddo ! i
        enddo ! j
    end do       ! RN

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_POP_ReadReactDaily_MOZART'
#endif
end subroutine Atm_POP_ReadReactDaily_MOZART


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine transforming reactants form [ppb] to [molec/cm3] or to [ug/m3]
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_POP_ReactTransformSteply

    integer i, j, k, t, AN, RN, iu
    real(8) RhoAir, Coeff
    real Tair, Psur
    real splRow(8), d2FdX(8), splVal, ReactMixR, dair

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_POP_ReactTransformSteply'
#endif

    do k=1, Atm_Kmax
        do j=Jmin, Jmax
            do i=minI(j), maxI(j)
                Psur=PxCurr(i,j)
                do AN = 1, NmbAero                 ! For all aerosol species
                    if (trim(AUnit(AN)) == 'ppbv') then
                    Tair=TairCurr(i,j,k)
                    RhoAir=Nav/Runiv*(Sigma(k)*Psur+Ptop)/Tair*1.e-6              ! molec/cm3
                      DensAirNum(i,j,k)=RhoAir                     ! ? For what?
                      Coeff = 1.e-9*RhoAir
                    elseif(trim(AUnit(AN)) == 'kg/kg' .or. trim(AUnit(AN)) == 'm2/kg') then
                      Coeff = DensAir(i,j,k)
                    else
                      print*, 'Aerosols: unknown units for ', trim(AName(AN)), ': ',trim(AUnit(AN))
                      stop
                    end if
                    Coeff = Coeff * CorrectA(AN)             ! Define Correct as global variable
                    splRow(1:4)=Aerofield(AN,i,j,k,1:4,toDay)
                    splRow(5:8)=Aerofield(AN,i,j,k,1:4,toMor)
                    d2FdX(1:5)=d2Aerofield(AN,i,j,k,1:5)
                    splVal=SplineInterpol(8,timePer,splRow,d2FdX,Period,DayTime)
                    ReactMixR=max(splVal,min(splRow(Period),splRow(Period+1)))
                    AeroConc(AN,i,j,k)=real(ReactMixR*Coeff,8)                  ! Unit correction. AeroConc to be defined instead of ConcOC,...
                enddo   ! AN
            enddo    ! j
        enddo    ! i
    end do   ! k

                do k=1, Atm_Kmax
        do j=Jmin, Jmax
            do i=minI(j), maxI(j)
                Psur=PxCurr(i,j)
                do RN = 1, NmbReacts                 ! For all reactants
                    if (ROI_T .and. O3Nmb == RN) then
                      Coeff = 1.             ! if units in the input file are ppvb
                    else
                      if (trim(RUnit(RN)) == 'ppbv') then
                    Tair=TairCurr(i,j,k)
                    RhoAir=Nav/Runiv*(Sigma(k)*Psur+Ptop)/Tair*1.e-6              ! molec/cm3
                        Coeff = 1.e-9*RhoAir
                      elseif(trim(RUnit(RN)) == 'kg/kg' .or. trim(RUnit(RN)) == 'm2/kg') then
                        Coeff = DensAir(i,j,k)
                      else
                        print*, 'Reactants: unknown units for ', trim(RName(RN)), ': ', trim(RUnit(RN))
                        stop
    end if
    end if
                    splRow(1:4)=Reactfield(RN,i,j,k,1:4,toDay)
                    splRow(5:8)=Reactfield(RN,i,j,k,1:4,toMor)
                    d2FdX(1:5)=d2Reactfield(RN,i,j,k,1:5)
                    splVal=SplineInterpol(8,timePer,splRow,d2FdX,Period,DayTime)
                    ReactMixR=max(splVal,min(splRow(Period),splRow(Period+1)))
                    ReactConc(RN,i,j,k)=real(ReactMixR*Coeff,8)                  ! Unit correction. ReactConc to be defined instead of ConcOH,...
                enddo   ! RN
            enddo    ! j
        enddo    ! i
    end do   ! k

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_POP_ReactTransformSteply'
#endif
end subroutine Atm_POP_ReactTransformSteply

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine transforming reactants form [ppb] to [molec/cm3] or to [ug/m3]
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Atm_POP_PartitDegrCoefs

    integer i, j, k, t, NPOP, AN, RN
!    real    Tair, Pol, Koa, Kad, Kab, Kec, Kair_ec, kd_OC_curr, kd_EC_curr
    real(8) RhoAir, Coeff, PSur, WaterConc
    real(8) Tair, Kp(NReactMax), KdG(NReactMax), KdP(NReactMax), Summ, Summ1, Qv      ! NReactMax should be defined as a global variable
    real    TCels, RH, O3Curr, Pxx, Pa, PSat, DegrPCoeff,   CorrCoeff

#ifdef DEBUG_MODE
    print *, '>- Entering Atm_POP_PartitDegrCoefs'
#endif

    do NPOP = 1, gSubsNum(GroupInd(POP))
      do k=1, Atm_Kmax
        do j=Jmin, Jmax
            do i=minI(j), maxI(j)
                if (ROI_T) Psur = PxCurr(i,j)
                    Tair=TairCurr(i,j,k)
! Coefficients for GPPartitioning
                Summ = 0.
                do AN = 1, NmbAero
                  Kp(AN) = Cp0(AN) * exp(CpT(AN)*(1./Tair - 1./T0)) * AeroConc(AN,i,j,k)
                  Summ = Summ + Kp(AN)
                end do
                kpt(i,j,k,Period,NPOP) = Summ / (1. + Summ)
! Coefficients for gas phase degradation
		Summ1 = 0.
		do RN = 1, NmbReacts
		  if (ROI_T .and. O3Nmb == RN) then
		    RhoAir=Nav/Runiv*(Sigma(k)*Psur+Ptop)/Tair*1.e-6              ! molec/cm3
            	    Coeff = 1.e-9*RhoAir
		  else
		    Coeff = 1.
		  end if
		  Summ1 = Summ1 + Kd0(RN,0) * exp(KdT(RN,0)*(1./Tair - 1./T0)) * ReactConc(RN,i,j,k) * Coeff    ! Kd0 and KdT defined for reactant RN and gas phase (0)
		end do
		kDegrG(i,j,k,Period,NPOP) = Summ1
! Coefficients for particle phase degradation

		if (ROI_T) then       ! Logical variable ROI_T should be read from Properties file or RunInfo?
		  Pxx=Px(i,j,Period,toDay)					! Surface pressure Px=Ps-Pt
		  Pa=Sigma(1)*Pxx+Ptop
		  if(Tair>273.15) then
		    Psat=6.112*exp(17.67*(Tair-273.15)/(Tair-29.65))*100.	 ! Saturation vapour pressure (t>0) [Pa]
		  else
		    Psat=6.11*exp(22.514-6150./Tair)*100.			 ! Saturation vapour pressure (t<0) [Pa]
		  endif
	          Qv=HumidAir(i,j,1,Period,toDay)				! Vapour mioxing ratio
		  RH=min(Qv*(1.+HumCoef)/(1.+Qv*(1.+HumCoef))*Pa/Psat,1.)	 ! Relative humidity (saturation ratio)
		  TCels = TAir - 273.15
		  O3Curr = ReactConc(O3Nmb,i,j,k)
		  call CalcDegrP_ROI_T(RH,TCels,O3Curr, DegrPCoeff)
		  kDegrP(i,j,k,Period,NPOP) = DegrPCoeff
		else  					! This part was modified 10-02-2020
		  Summ1 = 0.
		   if (WaterCorrect) then
		     WaterConc = HumidAir(i,j,k,Period,toDay) / 18.0153 * Nav / 1000.       ! molecules/cm3, 18.0153 g/mol - Water molar mass
		   end if
		   do AN = 1, NmbAero
		    do RN = 1, NmbReacts
		     CorrCoeff = 1. + KReact(RN,AN) * ReactConc(RN,i,j,k)
		     if (WaterCorrect) CorrCoeff = CorrCoeff + 2.1e-17 * WaterConc    ! 2.1e-17 cm3 = KH2O
		     Summ1 = Summ1 + Kd0(RN,AN) * exp(KdT(RN,AN)*(1./Tair - 1/T0)) * &
			         & ReactConc(RN,i,j,k) * Kp(AN) / CorrCoeff
		     end do
		   end do
		   if(Summ /= 0.) then
		    kDegrP(i,j,k,Period,NPOP) = Summ1 / Summ
		    else
		        print *, 'In Atm_POP_PartitDegrCoefs, Summ = 0 (i,j,k): ',i,j,k
		        stop
		    end if
	        end if ! if ROI_T
              enddo ! i
            enddo ! j
        enddo ! k
    end do ! NPOP

#ifdef DEBUG_MODE
    print *, '<- Exit Atm_POP_PartitDegrCoefs'
#endif
end subroutine Atm_POP_PartitDegrCoefs

subroutine CalcDegrP_ROI_T(RH,TCels,O3Curr,DegrPCoeff)

    real    TCels, RH, O3Curr, DegrPCoeff, DegrP1, DegrP2
    integer N_RH, N_O3, N_Temp, TypeAppr         ! Chanhed 14-08-2019

#ifdef DEBUG_MODE
    print *, '<- Entering CalcDegrP_ROI_T'
#endif
    TypeAppr = 1         ! Added 14-08-2019, use correct ROI_T_data!
	if (RH > 0.6) then
      N_RH = 1              ! Data for RH = 70%
    elseif (RH > 0.25) then
      N_RH = 2              ! Data for RH = 50%
                    else
      N_RH = 3              ! Data for RH = 0% (Dry)
                    end if
!	N_RH = 1         ! Added 14-08-2019
    if (TCels < -15.) then      ! "negative temperatures"
      if (O3Curr >= 100.) then     ! extrapolation wrt O3
!        DegrP1 = DegrP0(9,N_RH) * exp(DegrPT(9,N_RH) * TCels)           ! For ConcO3 = 90 ppb and TCels    ! Chanhed 14-08-2019
!        DegrP2 = DegrP0(10,N_RH) * exp(DegrPT(10,N_RH) * TCels)           ! For ConcO3 = 100 ppb and TCels
        DegrP1 = Approx(TypeAppr,DegrP0(9,N_RH),DegrPT(9,N_RH),TCels)           ! For ConcO3 = 90 ppb and TCels    ! Added 14-08-2019
        DegrP2 = Approx(TypeAppr,DegrP0(10,N_RH),DegrPT(10,N_RH),TCels)           ! For ConcO3 = 100 ppb and TCels
        DegrPCoeff = DegrP2 + (DegrP2 - DegrP1) / 10. * (O3Curr - 100.)    ! Read and define DegrP0 and DegrPT
      elseif(O3Curr > 10. .and. O3Curr < 100.) then   ! interpolation wrt O3
        N_O3 = floor(O3Curr / 10.)
!        DegrP1 = DegrP0(N_O3,N_RH) * exp(DegrPT(N_O3,N_RH) * TCels)           ! For TCels and lower ConcO3    ! Chanhed 14-08-2019
!        DegrP2 = DegrP0(N_O3 + 1,N_RH) * exp(DegrPT(N_O3 + 1,N_RH) * TCels)           ! For TCels and higher ConcO3
        DegrP1 = Approx(TypeAppr,DegrP0(N_O3,N_RH),DegrPT(N_O3,N_RH),TCels)           ! For ConcO3 = 90 ppb and TCels    ! Added 14-08-2019
        DegrP2 = Approx(TypeAppr,DegrP0(N_O3 + 1,N_RH),DegrPT(N_O3 + 1,N_RH),TCels)           ! For ConcO3 = 100 ppb and TCels
        DegrPCoeff = DegrP1 + (DegrP2 - DegrP1) / 10. *(O3Curr - real(N_O3) * 10.)
      else     ! O3Conc < 10, using data for O3Conc = 10.
        DegrPCoeff = DegrP0(1,N_RH) * exp(DegrPT(1,N_RH) * TCels)    ! Chanhed 14-08-2019
		DegrPCoeff = Approx(TypeAppr,DegrP0(1,N_RH),DegrPT(1,N_RH),TCels)
                    end if
! Read and define: BaseD, MaxD, RateD, XHalfD, TBound
    elseif (TCels >= -15. .and. TCels < 40.) then     ! interpolation wrt T
      N_Temp = floor((TCels + 15.) / 5.) + 1      ! if <= 7
      if (N_Temp > 7) N_Temp = N_Temp - 1        ! 20 omitted
      DegrP1 = BaseD(N_Temp,N_RH) + (MaxD(N_Temp,N_RH) - BaseD(N_Temp,N_RH)) / (1. + (XHalfD(N_Temp,N_RH) / O3Curr) &
                                                                                        & ** RateD(N_Temp,N_RH))     ! For lower TCels
      DegrP2 = BaseD(N_Temp + 1,N_RH) + (MaxD(N_Temp + 1,N_RH) - BaseD(N_Temp + 1,N_RH)) / (1. + (XHalfD(N_Temp &
                                                       & + 1,N_RH) / O3Curr) ** RateD(N_Temp,N_RH))     ! For higher TCels
      DegrPCoeff = DegrP1 + (DegrP2 - DegrP1) / (TBound(N_Temp + 1,N_RH) - TBound(N_Temp,N_RH)) &
                                                       & * (TCels - TBound(N_Temp,N_RH))
    else           ! TCels >= 40, using data for TCels = 40.
      DegrPCoeff = BaseD(11,N_RH) + (MaxD(11,N_RH) - BaseD(11,N_RH)) / (1. + (XHalfD(11,N_RH) / O3Curr) &
                                                                                        & ** RateD(11,N_RH))     ! For lower TCels
                    end if

#ifdef DEBUG_MODE
    print *, '<- Exit CalcDegrP_ROI_T'
#endif
end subroutine CalcDegrP_ROI_T

subroutine ReadROI_TParams

  integer NN, ii, DummyN

#ifdef DEBUG_MODE
    print *, '<- Entering ReadROI_TParams'
#endif
  open(3, file = 'ROI_T_Data.txt')
    do NN = 1, 3
      read(3,*)
      read(3,*)
      do ii = 1, 11
        read(3,*) TBound(ii,NN), BaseD(ii,NN), MaxD(ii,NN), RateD(ii,NN), XHalfD(ii,NN)
      end do
    end do
    do NN = 1, 3
      read(3,*)
      read(3,*)
      do ii = 1, 10
        read(3,*) DummyN, DegrP0(ii,NN), DegrPT(ii,NN)
      end do
    end do
  close(3)
#ifdef DEBUG_MODE
    print *, '<- Exit ReadROI_TParams'
#endif

end subroutine ReadROI_TParams

real function IphiM(L,Zr,Zo)

	real L, Zr, Zo
	real(8) Kr, Ko

    if(L==0.) L=1.e-6

    if(Zr/L<1.e-6) then                     ! Neutral
	  IphiM=log(Zr/Zo)
    elseif(L>=0.) then						! Stable
	  IphiM=log(Zr/Zo)+Bm/L*(Zr-Zo)
	else									! Unstable
	  Kr=(1.-dble(Gm*Zr/L))**0.25
	  Ko=(1.-dble(Gm*Zo/L))**0.25
	  IphiM=dlog((Kr-1.)/(Kr+1.))-dlog((Ko-1.)/(Ko+1.))+2*datan(Kr)-2.*datan(Ko)
	endif

end function IphiM

real function PsiH(x)
	real x
	if (x>0.) then		! stable
	  PsiH=-6.*x
	else				! non stable
	  PsiH=2.*alog(0.5*(1+Sqrt(1-9.*x)))
	endif
end function PsiH

real function IphiH(L,Zr,Zo)

	real L, Zr, Zo
	real(8) Kr, Ko

    if(L==0.) L=1.e-6

    if(Zr/L<1.e-6) then                     ! Neutral
	  IphiH=Prt*log(Zr/Zo)
	elseif(L>=0.) then						! Stable
	  IphiH=Prt*log(Zr/Zo)+Bh/L*(Zr-Zo)
	else									! Unstable
	  Kr=sqrt(1.-dble(Gh*Zr/L))
	  Ko=sqrt(1.-dble(Gh*Zo/L))
	  IphiH=Prt*(dlog((Kr-1.)/(Kr+1.))-dlog((Ko-1.)/(Ko+1.)))
	endif

end function IphiH
!...............................................................................
real function WetDiameter(Dd,S,e)

	real Dd, S, e
	real A, B, k, k1, k2

	A=1.2*exp(0.066*S/(1.058-S))
	B=exp(0.00077*S/(1.009-S))
	k1=10.2-23.7*S+14.5*S*S
	k2=-6.7+15.5*S-9.2*S*S
	k=1.-k1*(1.-e)-k2*(1.-e*e)

	WetDiameter=A*(Dd*1.e6)**B*k*1.e-6

end function WetDiameter
!...............................................................................

real function Approx(Tpe,Dg0,DgT,Tpr)        ! Added 14-08-2019

  integer Tpe
  real Dg0, DgT, Tpr

  if (Tpe == 1) then
    Approx = Dg0 * exp(DgT * Tpr)
  elseif(Tpe == 2) then
    Approx = Dg0 * exp(DgT * (1./(Tpr + 273.16)) - 1./273.15)
  else
    print*, 'Stop: invalid typr in function Approx'
	stop
  end if

end function Approx

end module Atm_POP_General
#endif
