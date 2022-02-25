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
! Module of the model balance output
! Updated: 19.10.2021
!  Output of increments of balance counters to Balance.dat file has been inserted 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

module Balance

    use GeneralParams
    use Geometry
    use Atm_Params
    use Exch_Params
#ifdef M_SOIL
    use Soil_Params
    use Soil_POP_General
#endif
#ifdef M_OCN
    use Ocn_Params
    use Ocn_POP_General
#endif
#ifdef M_VEG
    use Veg_Params
    use Veg_POP_General
#endif

    implicit none

    character(800), private  :: fileName, fullName

! Total mass in model domain
    real(8)       TotInit, TotAntEmis, TotNatEmis, TotUp
    real(8)       TotDegr, TotFin, TotInput, TotSum, TotBnd, Sedim

! Atmospheric compartment mass variables
    real(8)       AtmInit, AtmAntEmis, AtmNatEmis, ReEmis, DryDep, WetDep, AtmUp, AtmChemEx
    real(8)       AtmDegr, DryRem, AtmFin, AtmInput, AtmBnd, AtmPartit
    real(8)       AtmInput(MaxForm), AtmSum(MaxForm)

! Ocean compartment mass variables
    real(8)       OcnInit, OcnAntEmis, OcnNatEmis, OcnSedim
    real(8)       OcnDegr, OcnFin, OcnInput, OcnSum, OcnBnd, OcnPartit

! Soil compartment mass variables
    real(8)       SoilInit, SoilAntEmis, SoilNatEmis
    real(8)       SoilDegr, SoilFin, SoilInput, SoilSum, SoilPartit

! Veg compartment mass variables
    real(8)       FallInit, FallDegr, FallFin, FallSum
    real(8)       VegInit, VegDegr, VegFin, SoilInput, SoilSum, SoilPartit

! Balance output file variables

    integer, parameter    :: AtmBalVarNum = 12
    integer, parameter    :: AtmMassNo    = 1
    integer, parameter    :: AtmAntEmisNo = 2
    integer, parameter    :: AtmNatEmisNo = 3
    integer, parameter    :: AtmReEmisNo  = 4
    integer, parameter    :: AtmDryDepNo  = 5
    integer, parameter    :: AtmDryRemNo  = 6
    integer, parameter    :: AtmWetDepNo  = 7
    integer, parameter    :: AtmBndFluxNo = 8
    integer, parameter    :: AtmUpFluxNo  = 9
    integer, parameter    :: AtmChemExNo  = 10
    integer, parameter    :: AtmDegradNo  = 11
    integer, parameter    :: AtmInitNo  = 12
    real(8), private      :: AtmBalArr(AtmBalVarNum)
    real(8), private      :: AtmBalArrPrev(AtmBalVarNum)
    character(12),private :: AtmBalVarHdr(AtmBalVarNum)
    data AtmBalVarHdr / '     AtmMass','    AtmAntEm','    AtmNatEm','   AtmReemis','   AtmDryDep','   AtmDryRem','   AtmWetDep', &
                        '    AtmBndFl','     AtmUpFl','   AtmChemEx','   AtmDegrad','     AtmInit' /

#ifdef M_SOIL
    integer, parameter    :: SoilBalVarNum = 7
    integer, parameter    :: SoilMassNo    = 1
    integer, parameter    :: SoilAntEmisNo = 2
    integer, parameter    :: VegSoilFluxNo = 3
    integer, parameter    :: AtmSoilFluxNo = 4
    integer, parameter    :: SoilAtmFluxNo = 5
    integer, parameter    :: SoilDegradNo  = 6
    integer, parameter    :: SoilInitNo    = 7
    real(8), private      :: SoilBalArr(SoilBalVarNum)
    real(8), private      :: SoilBalArrPrev(SoilBalVarNum)
    character(12),private :: SoilBalVarHdr(SoilBalVarNum)
    data SoilBalVarHdr / '    SoilMass',' SoilAntEmis',' VegSoilFlux',' AtmSoilFlux',' SoilAtmFlux','  SoilDegrad','    SoilInit' /
#endif

#ifdef M_OCN
    integer, parameter    :: OcnBalVarNum = 7
    integer, parameter    :: OcnMassNo    = 1
    integer, parameter    :: OcnAntEmisNo = 2
    integer, parameter    :: AtmOcnFluxNo = 3
    integer, parameter    :: OcnAtmFluxNo = 4
    integer, parameter    :: OcnSedFluxNo = 5
    integer, parameter    :: OcnDegradNo  = 6
    integer, parameter    :: OcnInitNo    = 7
    real(8), private      :: OcnBalArr(OcnBalVarNum)
    real(8), private      :: OcnBalArrPrev(OcnBalVarNum)
    character(12),private :: OcnBalVarHdr(OcnBalVarNum)
    data OcnBalVarHdr / '     OcnMass','  OcnAntEmis','  AtmOcnFlux','  OcnAtmFlux','  OcnSedFlux','   OcnDegrad','     OcnInit' /
#endif

#ifdef M_VEG
    integer, parameter    :: VegBalVarNum = 9
    integer, parameter    :: VegMassNo    = 1
    integer, parameter    :: FallMassNo   = 2
    integer, parameter    :: AtmVegFluxNo = 3
    integer, parameter    :: VegAtmFluxNo = 4
    integer, parameter    :: VegFallFluxNo= 5
    integer, parameter    :: VegDegradNo  = 6
    integer, parameter    :: FallDegradNo = 7
    integer, parameter    :: VegInitNo    = 8
    integer, parameter    :: FallInitNo   = 9
    real(8), private      :: VegBalArr(VegBalVarNum)
    real(8), private      :: VegBalArrPrev(VegBalVarNum)
    character(12),private :: VegBalVarHdr(VegBalVarNum)
    data VegBalVarHdr / '     VegMass','    FallMass',' AtmVegFlux',' VegAtmFlux',' VegFallFlux','   VegDegrad',&
                       &'  FallDegrad','     VegInit','    FallInit'/
#endif

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutines calculating pollutant mass balance
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine OutputBalance(Prd)

    character(*), intent(in)    :: Prd          ! Temporal Frequency
    integer	Form

    select case(Prd)
        case('Initial')
	    if(InitCond == 'zero' .or. InitCond == 'cond') then             ! Added by OT 31.05.2017
              call WriteBalanceFileInit
	    endif                                                           ! Added by OT 31.05.2017

	    if(InitCond == 'cond') then             ! Added by AG 13.12.2021
		do Form = 1, NumForm(Atm)
		    AtmInit = AtmInit   + MassAtmInit(Form)
		end do
		AtmBalArrPrev(AtmMassNo)   = AtmInit
#ifdef M_SOIL
		SoilBalArrPrev(SoilMassNo) = SoilInit
#endif
#ifdef M_OCN
		OcnBalArrPrev(OcnMassNo)   = OcnInit
#endif
#ifdef M_VEG
		VegBalArrPrev(VegMassNo)   = VegInit
		VegBalArrPrev(FallMassNo)  = FallInit
#endif
	    endif                                   	! Added by AG 13.12.2021
	    
	    if(InitCond == 'dump') then			! Added by AG 13.12.2021
	    AtmBalArrPrev = 0.
	    do Form = 1, NumForm(Atm)
	    AtmBalArrPrev(1) = AtmBalArrPrev(1)+MassAtmFin(Form)
	    AtmBalArrPrev(2) = AtmBalArrPrev(2)+MassAtmAntEmis(Form)
	    AtmBalArrPrev(3) = AtmBalArrPrev(3)+MassAtmNatEmis(Form)
	    AtmBalArrPrev(4) = AtmBalArrPrev(4)+MassReEmis(Form)
	    AtmBalArrPrev(5) = AtmBalArrPrev(5)-MassDryDep(Form)
	    AtmBalArrPrev(6) = AtmBalArrPrev(6)+MassDryRem(Form)
	    AtmBalArrPrev(7) = AtmBalArrPrev(7)-MassWetDep(Form)
	    AtmBalArrPrev(8) = AtmBalArrPrev(8)+MassAtmBnd(Form)
	    AtmBalArrPrev(9) = AtmBalArrPrev(9)+MassAtmUp(Form)
	    AtmBalArrPrev(10)= AtmBalArrPrev(10)+MassChemEx(Form)
	    AtmBalArrPrev(11)= AtmBalArrPrev(11)-MassAtmDegr(Form)
	    AtmBalArrPrev(12)= AtmInit
	    end do
#ifdef M_SOIL
	    SoilBalArrPrev(SoilMassNo) = Soil_POPMass(1)
	    SoilBalArrPrev(2) = MasssoilEmis(1)
	    SoilBalArrPrev(3) = VegSoilFlux
	    SoilBalArrPrev(4) = AirSoilFlux
	    SoilBalArrPrev(5) = SoilAirFlux
	    SoilBalArrPrev(6) = -Soil_Degr
	    SoilBalArrPrev(7) = SoilInit
#endif
#ifdef M_OCN
	    OcnBalArrPrev(OcnMassNo)   = Ocn_POPMass(1)
	    OcnBalArrPrev(2) = 0.
	    OcnBalArrPrev(3) = AirOcnFlux
	    OcnBalArrPrev(4) = OcnAirFlux
	    OcnBalArrPrev(5) = -Ocn_Sedim
	    OcnBalArrPrev(6) = -Ocn_Degr
	    OcnBalArrPrev(7) = OcnInit
#endif
#ifdef M_VEG
	    VegBalArrPrev(VegMassNo)   = Veg_POPMass(1)
	    VegBalArrPrev(FallMassNo)  = Fall_POPMass(1)
	    VegBalArrPrev(3) = AirVegFlux
	    VegBalArrPrev(4) = VegAirFlux
	    VegBalArrPrev(5) = -VegSoilFlux
	    VegBalArrPrev(6) = -Veg_Degr
	    VegBalArrPrev(7) = 0.
	    VegBalArrPrev(8) = VegInit
	    VegBalArrPrev(9) = FallInit
#endif
	    endif                                      ! Added by AG 13.12.2021

        case('Hourly')
            ! Not yet implemented
        case('6hourly')
            if(outputBalance6hourly) call CalcBalance
        case('Daily')
            if(outputBalanceDaily)   call CalcBalance
        case('Monthly')
            if(outputBalanceMonthly) call CalcBalance
        case('Yearly')
            if(outputBalanceYearly)  call CalcBalance
    end select
    
end subroutine OutputBalance


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Calculate current mass balance
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine CalcBalance

	integer     i, j, k, Form

#ifdef DEBUG_MODE_BALANCE
     print *, 'Enter CalcBalance...'
#endif

! Calculate current mass in the atmosphere
	MassAtmFin=0.
	do Form=1, NumForm(Atm)
	  do k=1, Atm_Kmax
		do j=Jmin, Jmax
		  do i=minI(j), maxI(j)
			MassAtmFin(Form)=MassAtmFin(Form)+Atm_Conc(i,j,k,Form)*MeshVolume(i,j,k)
		  enddo
		enddo
	  enddo
	enddo

	TotInit=0.
	AtmAntEmis=0.
	AtmNatEmis=0.
	ReEmis=0.
	DryDep=0.
	DryRem=0.
	WetDep=0.
	TotUp=0.
	AtmBnd=0.
        TotBnd=0.
	AtmChemEx=0.
	TotInput=0.
	TotSum=0.
#ifdef G_POP
	AtmDegr=0.
        TotDegr=0.
#endif
        AtmSum=0.
	AtmFin=0.
        AtmInput = 0.
        AtmInit = 0.

	do Form = 1, NumForm(Atm)
	  AtmInit = AtmInit   + MassAtmInit(Form)
	  AtmAntEmis = AtmAntEmis + MassAtmAntEmis(Form)
	  AtmNatEmis = AtmNatEmis + MassAtmNatEmis(Form)
	  ReEmis  = ReEmis  + MassReEmis(Form)
	  DryDep  = DryDep  - MassDryDep(Form)
	  DryRem  = DryRem  + MassDryRem(Form)
	  WetDep  = WetDep  - MassWetDep(Form)
#ifdef G_POP
	  AtmDegr = AtmDegr - MassAtmDegr(Form)
#endif
          AtmBnd  = AtmBnd + MassAtmBnd(Form)
	  TotUp   = TotUp  + MassAtmUp(Form)
	  AtmFin  = AtmFin + MassAtmFin(Form)
#ifdef G_POP
	  AtmInput(Form)=AtmInput(Form) + MassAtmInit(Form) + MassAtmAntEmis(Form) + MassAtmNatEmis(Form)
	  AtmSum(Form) = MassAtmInit(Form) + MassAtmAntEmis(Form) + MassAtmNatEmis(Form) + &
                        &MassAtmUp(Form) - MassDryDep(Form) + MassDryRem(Form) - MassWetDep(Form) +&
			&MassAtmPartit(Form) - MassAtmDegr(Form) + MassAtmBnd(Form)
#else
	  AtmInput(Form)=AtmInput(Form) + MassAtmInit(Form) + MassAtmAntEmis(Form) + MassAtmNatEmis(Form) + MassReEmis(Form)
	  AtmSum(Form) = MassAtmInit(Form) + MassAtmAntEmis(Form) + MassAtmNatEmis(Form) + MassReEmis(Form) + &
			&MassAtmUp(Form) - MassDryDep(Form) - MassWetDep(Form) + MassChemEx(Form) + MassAtmBnd(Form)
	  AtmChemEx=AtmChemEx+MassChemEx(Form)
#endif
	enddo ! Form

        TotInit = AtmInit
        TotFin = AtmFin
        TotAntEmis = AtmAntEmis
        TotNatEmis = AtmNatEmis
        TotDegr = AtmDegr
        TotBnd = AtmBnd

#ifdef M_SOIL
        TotInit = TotInit + SoilInit
        TotFin  = TotFin + Soil_POPMass(1)
        TotDegr = TotDegr - Soil_Degr
        TotAntEmis = TotAntEmis + sum(MassSoilEmis(:))
#endif
#ifdef M_OCN
        TotInit = TotInit + OcnInit
        TotFin  = TotFin + Ocn_POPMass(1)
        TotDegr = TotDegr - Ocn_Degr
#endif
#ifdef M_VEG
        TotInit = TotInit + VegInit + FallInit
        TotFin  = TotFin + Veg_POPMass(1) + Fall_POPMass(1)
        TotDegr = TotDegr - Veg_Degr
#endif

#ifdef G_POP
	TotInput = TotInit + TotAntEmis + TotNatEmis
	TotSum = TotInit + TotAntEmis + TotNatEmis + TotBnd + TotUp - Ocn_Sedim + TotDegr
#else
	TotInput = TotInit + TotAntEmis + TotNatEmis + ReEmis
	TotSum = TotInit + TotAntEmis + TotNatEmis + ReEmis + TotBnd + TotUp + DryDep + AtmChemEx + WetDep
#endif

        call WriteBalanceStringToFile
        
	AtmBalArrPrev  = AtmBalArr
#ifdef M_SOIL
	SoilBalArrPrev = SoilBalArr
#endif
#ifdef M_OCN
	OcnBalArrPrev  = OcnBalArr
#endif
#ifdef M_VEG
	VegBalArrPrev  = VegBalArr
#endif

#ifdef DEBUG_MODE_BALANCE
     print *, 'Exit CalcBalance.'
#endif

end subroutine CalcBalance


subroutine WriteBalanceStringToFile

    character(120)      fmt
    character(2)        rep
    integer             num, i, j, k, l, st


#ifdef DEBUG_MODE_BALANCE
     print *, 'Enter WriteBalanceStringToFile...'
#endif

    AtmBalArr(1) = AtmFin	! Current mass in the atmosphere
    AtmBalArr(2) = AtmAntEmis
    AtmBalArr(3) = AtmNatEmis
    AtmBalArr(4) = ReEmis
    AtmBalArr(5) = DryDep
    AtmBalArr(6) = DryRem
    AtmBalArr(7) = WetDep
    AtmBalArr(8) = AtmBnd
    AtmBalArr(9) = TotUp
    AtmBalArr(10) = AtmChemEx
    AtmBalArr(11)= AtmDegr
    AtmBalArr(12)= AtmInit
#ifdef M_SOIL
    SoilBalArr(1) = Soil_POPMass(1)
    SoilBalArr(2) = MasssoilEmis(1)
    SoilBalArr(3) = VegSoilFlux
    SoilBalArr(4) = AirSoilFlux
    SoilBalArr(5) = SoilAirFlux
    SoilBalArr(6) = -Soil_Degr
    SoilBalArr(7) = SoilInit
#endif
#ifdef M_OCN
    OcnBalArr(1) = Ocn_POPMass(1)
    OcnBalArr(2) = 0.
    OcnBalArr(3) = AirOcnFlux
    OcnBalArr(4) = OcnAirFlux
    OcnBalArr(5) = -Ocn_Sedim
    OcnBalArr(6) = -Ocn_Degr
    OcnBalArr(7) = OcnInit
#endif
#ifdef M_VEG
    VegBalArr(1) = Veg_POPMass(1)
    VegBalArr(2) = Fall_POPMass(1)
    VegBalArr(3) = AirVegFlux
    VegBalArr(4) = VegAirFlux
    VegBalArr(5) = -VegSoilFlux
    VegBalArr(6) = -Veg_Degr
    VegBalArr(7) = 0.
    VegBalArr(8) = VegInit
    VegBalArr(9) = FallInit
#endif

    num = AtmBalVarNum
#ifdef M_SOIL
    num = num + SoilBalVarNum
#endif
#ifdef M_OCN
    num = num + OcnBalVarNum
#endif
#ifdef M_VEG
    num = num + VEGBalVarNum
#endif

    write(rep,'(i2)') num
    fmt = '(i4,i4.2,i4.2,'//trim(rep)//'e12.4,4e12.4)'
    fullName=trim(OutPath)//trim(BalanceDir)//'Balance.dat'
    open(103, file=fullName, action='write', status='old', access='append', iostat = st)
    if(st /= 0) then
        print '(/,"STOP in WriteBalanceStringToFile: Cannot create file ''",a,"''",/)', trim(fullName)
        stop
    endif

    write(103, fmt) Year,Month,Day, ((AtmBalArr(i)-AtmBalArrPrev(i),i=1,AtmBalVarNum) &
#ifdef M_SOIL
                    &,(SoilBalArr(j)-SoilBalArrPrev(j),j=1,SoilBalVarNum) &
#endif
#ifdef M_OCN
                    &,(OcnBalArr(k) - OcnBalArrPrev(k),k=1,OcnBalVarNum) &
#endif
#ifdef M_VEG
                    &,(VegBalArr(l) - VegBalArrPrev(l),l=1,VegBalVarNum) &
#endif
    &,TotSum, TotFin, TotSum/(TotInput+Zero)*100.,TotFin/(TotInput+Zero)*100.)
    close(103)

#ifdef DEBUG_MODE_BALANCE
     print *, 'Exit WriteBalanceStringToFile...'
#endif

end subroutine WriteBalanceStringToFile



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Start writing balance files from scratch
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WriteBalanceFileInit

    character(600)      tstr
    character(200)      tstr1
    character(120)      fmt
    character(2)        rep
    integer             i, st

    fullName=trim(OutPath)//trim(BalanceDir)//'Balance.dat'
    open(103, file=fullName, action='write', status='replace', access='sequential', iostat=st)
    if(st /= 0) then
        print '(/,"STOP in WriteBalanceFileInit: Cannot create file ''",a,"''",/)', trim(fullName)
        stop
    endif
    write(103,*)'Mass balance in media, kg'

    write(rep,'(i2)') AtmBalVarNum
    fmt = '(a12,'//trim(rep)//'a12)'
    write(tstr1, fmt) 'Year Mon Day',(AtmBalVarHdr(i),i=1,AtmBalVarNum)
    tstr = tstr1
#ifdef M_SOIL
    write(rep,'(i2)') SoilBalVarNum
    fmt = '('//trim(rep)//'a12)'
    write(tstr1, fmt) (SoilBalVarHdr(i),i=1,SoilBalVarNum)
    tstr = trim(tstr)//trim(tstr1)
#endif
#ifdef M_OCN
    write(rep,'(i2)') OcnBalVarNum
    fmt = '('//trim(rep)//'a12)'
    write(tstr1, fmt) (OcnBalVarHdr(i),i=1,OcnBalVarNum)
    tstr = trim(tstr)//trim(tstr1)
#endif
#ifdef M_VEG
    write(rep,'(i2)') VegBalVarNum
    fmt = '('//trim(rep)//'a12)'
    write(tstr1, fmt) (VegBalVarHdr(i),i=1,VegBalVarNum)
    tstr = trim(tstr)//trim(tstr1)
#endif
    write(103, '(a,4a12)') trim(tstr),'     MassBal','     MassFin','   MassBal,%','   MassFin,%'
    close(103)

end subroutine WriteBalanceFileInit


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine printing information to the screen
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PrintBalanceScr

#ifdef DEBUG_MODE_BALANCE
     print *, 'Enter PrintBalanceScr...'
#endif

    print 201, TotInit,    abs(TotInit/(TotInput+Zero))*100.,&
            &  TotAntEmis, abs(TotAntEmis/(TotInput+Zero))*100.,&
            &  TotNatEmis, abs(TotNatEmis/(TotInput+Zero))*100.,&
            &  ReEmis,     abs(ReEmis/(TotInput+Zero))*100.,&
            &  DryDep,     abs(DryDep/(TotInput+Zero))*100.,&
            &  WetDep,     abs(WetDep/(TotInput+Zero))*100.,&
            &  TotDegr,    abs(TotDegr/(TotInput+Zero))*100.,&
            &  TotBnd,     abs(TotBnd/(TotInput+Zero))*100.,&
            &  TotUp,      abs(TotUp/(TotInput+Zero))*100.,&
            &  TotSum,     abs(TotSum/(TotInput+Zero))*100.,&
            &  TotFin,     abs(TotFin/(TotInput+Zero))*100.
    print 203, timeRun/60.

201	format(/,'------------------  MASS BALANCE  ------------------',/,/,&
            &'Mass input:',/,&
            &e20.9,' kg (',f6.1,'%) - initial mass',/,&
            &e20.9,' kg (',f6.1,'%) - anthropogenic mass emitted',/,&
            &e20.9,' kg (',f6.1,'%) - natural mass emitted',/,&
            &e20.9,' kg (',f6.1,'%) - re-emitted mass',/,/,&
            &'Mass fluxes (% of input mass):',/,&
            &e20.9,' kg (',f6.1,'%) - mass dry deposited',/,&
            &e20.9,' kg (',f6.1,'%) - mass wet deposited',/,&
            &e20.9,' kg (',f6.1,'%) - mass degraded',/,&
            &e20.9,' kg (',f6.1,'%) - mass through lateral boundary',/,&
            &e20.9,' kg (',f6.1,'%) - mass through upper boundary',/,/,&
            &'Mass balance (% of input mass):',/,&
            &e20.9,' kg (',f6.1,'%) - mass balance',/,&
            &e20.9,' kg (',f6.1,'%) - final mass',/)
203	format(f12.5,' (min) - calculation time',/)

#ifdef DEBUG_MODE_BALANCE
     print *, 'Exit PrintBalanceScr.'
#endif

end subroutine PrintBalanceScr


subroutine OutputBalanceInLogFile(fp)

    character(*), intent(in)     :: fp              ! Path to log file

    integer             i, f, fs
    integer         ::  lenField=60
    character(600)      frmt1, frmt11, frmt2, frmt22,frmt3
    character(2)        rep
    real(8)             s, om, sm, em, vm, fm, im, percent
    integer(8)          curVal(8)

#ifdef DEBUG_MODE_BALANCE
     print *, 'Enter OutputBalanceInLogFile...'
#endif

    open(101, file=trim(fp), action='write', status='old', access='append', iostat=fs)
    if(fs > 0) then
      print '(/,"STOP: Cannot open log file ''",a,"''",/)', trim(fp)
      stop
    endif

      frmt1='(/,"---------------------  MASS BALANCE  ----------------------",//,&
            &"Mass input:",/,&
            &e20.9," kg (",f5.1,"%) - initial mass",/,&
            &e20.9," kg (",f5.1,"%) - anthropogenic mass emitted",/,&
            &e20.9," kg (",f5.1,"%) - natural mass emitted",/,&
            &e20.9," kg (",f5.1,"%) - re-emitted mass",/,/)'

      frmt11='("Mass fluxes (% of input mass):",/,&
            &e20.9," kg (",f5.1,"%) - mass dry deposited",/,&
            &e20.9," kg (",f5.1,"%) - mass wet deposited",/,&
            &e20.9," kg (",f5.1,"%) - mass degradated",/,&
            &e20.9," kg (",f5.1,"%) - mass through lateral boundary",/,&
            &e20.9," kg (",f5.1,"%) - mass through upper boundary",/,/,&
            &"Mass balance (% of input mass):",/,&
            &e20.9," kg (",f5.1,"%) - mass balance",/,&
            &e20.9," kg (",f5.1,"%) - final mass",/)'

        write(101,frmt1) TotInit, abs(TotInit/(TotInput+Zero))*100.,&
                & TotAntEmis, abs(TotAntEmis/(TotInput+Zero))*100.,&
                & TotNatEmis, abs(TotNatEmis/(TotInput+Zero))*100.,&
                & ReEmis, abs(ReEmis/(TotInput+Zero))*100.

        write(101,frmt11) DryDep, abs(DryDep/(TotInput+Zero))*100.,&
                & WetDep,  abs(WetDep/(TotInput+Zero))*100.,&
                & TotDegr, abs(TotDegr/(TotInput+Zero))*100.,&
                & TotBnd, abs(TotBnd/(TotInput+Zero))*100.,&
                & TotUp, abs(TotUp/(TotInput+Zero))*100.,&
                & TotSum, abs(TotSum/(TotInput+Zero))*100.,&
                & TotFin, abs(TotFin/(TotInput+Zero))*100.

        write(rep,'(i2)') NumForm(Atm)
        frmt2='(/,"-------------  ATMOSPHERE: DETAILED BALANCE  ----------------",/,/,&
			&24x,'//trim(adjustl(rep))//'(a10,3x),/,&
			&"Mass input, kg",/,&
			&4x,"Initial mass:      ",'//trim(adjustl(rep))//'(e11.4,2x),/,&
			&4x,"Anthrop. emission: ",'//trim(adjustl(rep))//'(e11.4,2x),/,&
			&4x,"Natural emission:  ",'//trim(adjustl(rep))//'(e11.4,2x),/,&
			&4x,"Re-emission:       ",'//trim(adjustl(rep))//'(e11.4,2x),/,/)'
        frmt22='("Mass fluxes, kg",/,&
			&4x,"Chem. exchange:    ",'//trim(adjustl(rep))//'(e11.4,2x),/,&
			&4x,"Partitioning:      ",'//trim(adjustl(rep))//'(e11.4,2x),/,&
			&4x,"Dry deposition:    ",'//trim(adjustl(rep))//'(e11.4,2x),/,&
			&4x,"Wet deposition:    ",'//trim(adjustl(rep))//'(e11.4,2x),/,&
			&4x,"Remobilization:    ",'//trim(adjustl(rep))//'(e11.4,2x),/,&
			&4x,"Degradation:       ",'//trim(adjustl(rep))//'(e11.4,2x),/,&
			&4x,"Lateral boundary:  ",'//trim(adjustl(rep))//'(e11.4,2x),/,/,&
			&4x,"Upper boundary:    ",'//trim(adjustl(rep))//'(e11.4,2x),/,/,&
			&"Mass balance, kg",/,&
			&4x,"Mass balance:      ",'//trim(adjustl(rep))//'(e11.4,2x),/,&
			&4x,"Final mass:        ",'//trim(adjustl(rep))//'(e11.4,2x),/)'
	write(101,frmt2) (trim(SubsID(1))//trim(FormID(1,1,f)), f=1, NumForm(Atm)),&
			& (MassAtmInit(f), f=1, NumForm(Atm)),&
			& (MassAtmAntEmis(f), f=1, NumForm(Atm)),&
			& (MassAtmNatEmis(f), f=1, NumForm(Atm)),&
			& (MassReEmis(f), f=1, NumForm(Atm))

	write(101,frmt22) (MassChemEx(f), f=1, NumForm(Atm)),&
			& (MassAtmPartit(f), f=1, NumForm(Atm)),&
			& (-MassDryDep(f), f=1, NumForm(Atm)),&
			& (-MassWetDep(f), f=1, NumForm(Atm)),&
			& (MassDryRem(f), f=1, NumForm(Atm)),&
			& (-MassAtmDegr(f), f=1, NumForm(Atm)),&
			& (MassAtmBnd(f), f=1, NumForm(Atm)),&
			& (MassAtmUp(f), f=1, NumForm(Atm)),&
			& (AtmSum(f), f=1, NumForm(Atm)),&
			& (MassAtmFin(f), f=1, NumForm(Atm))

#ifdef G_POP

#ifdef M_OCN
      write(101,'(a55)')'-----------------  OCEAN: BALANCE  --------------------'
      write(101,*)
    write(101,'(a19,e13.6)')'Initial mass:      ', OcnInit
    write(101,'(a19,e13.6)')'Flux from air:     ', AirOcnFlux
    write(101,'(a19,e13.6)')'Flux to air:       ', OcnAirFlux
    write(101,'(a19,e13.6)')'Degradation:       ', -Ocn_Degr
    write(101,'(a19,e13.6)')'Sedimentation:     ', -Ocn_Sedim
    write(101,'(a19,e13.6)')'Mass balance:      ', OcnInit+AirOcnFlux+OcnAirFlux-Ocn_Sedim-Ocn_Degr
    write(101,'(a19,e13.6)')'Final mass:        ', Ocn_POPMass(1)
#endif
#ifdef M_SOIL
      write(101,*)
      write(101,'(a55)')'-----------------  SOIL: BALANCE  ---------------------'
      write(101,*)
    write(101,'(a19,e13.6)')'Initial mass:      ', SoilInit
    write(101,'(a19,e13.6)')'Emission:          ', sum(MassSoilEmis)
    write(101,'(a19,e13.6)')'Flux from air:     ', AirSoilFlux
    write(101,'(a19,e13.6)')'Flux to air:       ', SoilAirFlux
    write(101,'(a19,e13.6)')'Flux from veg:     ', VegSoilFlux
    write(101,'(a19,e13.6)')'Degradation:       ', -Soil_Degr
    write(101,'(a19,e13.6)')'Mass balance:      ', SoilInit+sum(MassSoilEmis)+AirSoilFlux+SoilAirFlux+VegSoilFlux-Soil_Degr
    write(101,'(a19,e13.6)')'Final mass:        ', Soil_POPMass(1)
#endif
#ifdef M_VEG
      write(101,*)
      write(101,'(a55)')'--------------  VEGETATION: BALANCE  ------------------'
      write(101,*)
    write(101,'(a19,e13.6)')'Initial mass:      ', VegInit+FallInit
    write(101,'(a19,e13.6)')'Flux from air:     ', AirVegFlux
    write(101,'(a19,e13.6)')'Flux to air:       ', VegAirFlux
    write(101,'(a19,e13.6)')'Flux to soil:      ', -VegSoilFlux
    write(101,'(a19,e13.6)')'Degradation:       ', -Veg_Degr
    write(101,'(a19,e13.6)')'Mass balance:      ', VegInit+FallInit+AirVegFlux+VegAirFlux-VegSoilFlux-Veg_Degr
    write(101,'(a19,e13.6)')'Final mass:        ', Veg_POPMass(1)+Fall_POPMass(1)
#endif


#endif

    frmt3='(f12.5," (min) - calculation time")'
    write(101,frmt3) timeRun/60.
    write(101,'(/,a)') repeat('*',lenField)
    call date_and_time(values=curVal)
    write(101,'(/,"Finished: ",i2,"-",a3,"-",i4," at ",i2.2,":",i2.2,":",i2.2)') &
    & curVal(3), MonthName(curVal(2)),&
    & curVal(1), curVal(5), curVal(6), curVal(7)

    close(101)

#ifdef DEBUG_MODE_BALANCE
     print *, 'Exit OutputBalanceInLogFile...'
#endif

end subroutine OutputBalanceInLogFile



end module Balance
