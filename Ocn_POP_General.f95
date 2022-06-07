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
module Ocn_POP_General

  use GeneralParams
  use Ocn_Params
  use Ocn_POP_Params
  use Atm_POP_Params
  use Exch_Params
!  use Exch_Params_POP

  use OcnAdvDiff     !2011
  
  implicit none

  character(800), private :: fileName, fullName

  contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine defining physical and chemical properties of POP in the ocean
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Ocn_POP_Props(Subs)

    integer i, FileStat, lenType, Subs, Form, Limit, Ind
    character(80) strRead, strPar
    character(30) strType, FormType(MaxForm)
    integer mon, ios
        
    fileName='Ocn'//trim(PropName)//trim(SubsID(Subs))//'.dat'
    fullName=trim(PropPath)//fileName
    open(2, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file of ocean POP properties ''",a,"''",/)', trim(fullName)
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

! Deleting tabs and line termination symbols; codes: 9=\t, 32=space, 13=\n
        lenType=scan(strRead, achar(9))
        do while (lenType>0)
          strRead=strRead(:lenType-1)//strRead(lenType+1:)
          lenType=scan(strRead, achar(9)) 
        enddo
        lenType=scan(strRead,achar(13))
        if(lenType>0) strRead=strRead(:lenType-1)

! Reading pollutant form
        if(strRead(1:1)=='[') then
          Limit=scan(strRead(2:),']')
          lenType=Limit
          strType=trim(strRead(2:lenType))

          Form=0
          do i=1, sFormNum(Ocn,Subs)
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
          case('Forms number')
            read(strPar,'(i2)') sFormNum(Ocn,Subs)
            FormSubs(Ocn,NumForm(Ocn)+1:NumForm(Ocn)+sFormNum(Ocn,Subs))=Subs
            if (sFormNum(Ocn,Subs).gt.0) then
              NumSubsMedia(Ocn)=NumSubsMedia(Ocn)+1
              gSubsMediaInd(Ocn,Subs)=NumSubsMedia(Ocn)
            endif
!------------------------------------------------------------------
          case('Forms')
            read(strPar,*) (FormType(i), i=1, sFormNum(Ocn,Subs))
!------------------------------------------------------------------
          case('Form ID')
            FormID(Ocn,Subs,Form)=strPar
            sFormInd(Ocn,Subs,Form)=NumForm(Ocn)+1
            FormSubs(Ocn,sFormInd(Ocn,Subs,Form))=Subs
            NumForm(Ocn)=NumForm(Ocn)+1
!------------------------------------------------------------------
          case('SeaDegr')
            read(strPar,*) CDWater(Form,gSubsGroupInd(Subs))
!------------------------------------------------------------------
          case('SeaRedistr')
            read(strPar,*)k_ph_rd(gSubsGroupInd(Subs))
!------------------------------------------------------------------
        endselect
      endif
    enddo 

    close(2)

end subroutine Ocn_POP_Props

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main POP processes in the ocean
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine Ocn_POPProcess(NSubs, dTime)

    integer NSubs,i,j,k,n,kvert,ind1,ind2,k1,k2,l
    real(8) CTot, neg_conc, t1
    real    dTime
    integer Xscal

#ifdef DEBUG_MODE
    print *, '+Entering Ocn_Process... '
#endif

#if (OCEAN_RUN==1)
    call Ocn_adv_diff_step(NSubs)
#else
    ! Swithed off for run without ocean Sep 2016 V. Sh
#endif
    
!    where(Ocn_Conc < 0.) Ocn_Conc = 0.             ! Cut off negative concentrations
    call Ocn_POPdegradation(NSubs, dTime)
    call Ocn_POPSedimentation(NSubs, dTime)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OCEAN 3-D -> 1-D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AG 16022018
    n = sFormInd(Ocn,NSubs,1)
    Water_upl_conc = 0.
    ind1 = sFormInd(Ocn,NSubs,1)
    ind2 = sFormInd(Ocn,NSubs,2)

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
! From amounts to concentrations (kg/m3) calculated similar to air
                Water_upl_conc(i,j,ind1) = Water_upl_conc(i,j,ind1)/Ocn_Frac(i,j)/MeshArea(i,j)
              end if
          end do
    enddo ! j

    do j = JMin, JMax
        do i = MinI(j), MaxI(j)
            CTot = Water_upl_conc(i,j,sFormInd(Ocn,NSubs,1)) + Water_upl_conc(i,j,sFormInd(Ocn,NSubs,2))
            Water_upl_concDay(i,j,NSubs) = Water_upl_concDay(i,j,NSubs) + CTot * dTime
        end do
    end do

#ifdef DEBUG_MODE
    print *, '+Exit Ocn_Process'
#endif
end subroutine Ocn_POPProcess

!===========================================================
! Time step for POPs
!===========================================================
subroutine Ocn_POPTStepCalc(NSubs)

  integer NSubs

  continue

end subroutine Ocn_POPTStepCalc

!===========================================================
! Partitioning process for POPs
!===========================================================
subroutine Ocn_POPpartitioning(NSubs)

  integer NSubs
  integer i, j
  real(8) CTot, c1, c2

#ifdef DEBUG_MODE
    print *, '-> Entering Ocn_POPpartitioning... '
#endif

    do j = JMin, JMax
        do i = MinI(j), MaxI(j)
            c1 = Water_upl_conc(i,j,sFormInd(Ocn,NSubs,1))
            c2 = Water_upl_conc(i,j,sFormInd(Ocn,NSubs,2))
            CTot = c1 + c2
            Water_upl_conc(i,j,sFormInd(Ocn,NSubs,1))=CTot /(1.+k_ph_rd(gSubsGroupInd(NSubs)))
            Water_upl_conc(i,j,sFormInd(Ocn,NSubs,2))=k_ph_rd(gSubsGroupInd(NSubs))*CTot/(1.+k_ph_rd(gSubsGroupInd(NSubs)))
        end do
    end do

#ifdef DEBUG_MODE
    print *, '<- Exit Ocn_POPpartitioning'
#endif
end subroutine Ocn_POPpartitioning

!===========================================================
! Degradation process for POPs
!===========================================================
subroutine Ocn_POPdegradation(NSubs, dTime)

  integer NSubs
  integer i, j, k,n
  real DTime, cdw
  real*8 OPM1,OPM2

  OPM1 = Ocn_POPMass(NSubs)

      n=sFormInd(Ocn,NSubs,1)
      cdw = (CDWater(1,1)+CDWater(2,1)*k_ph_rd(1))/(1+k_ph_rd(1))   ! Temporarely defined. Needs to be replaced lately!!!

      do i=1,Imax
        do j=1,Jmax
            do k=1,Ocn_Kmax
                Ocn_conc(i,j,k,n,oldtime)=Ocn_conc(i,j,k,n,oldtime)-Ocn_conc(i,j,k,n,curtime)*dble(cdw*dTime)
                Ocn_conc(i,j,k,n,curtime)=Ocn_conc(i,j,k,n,curtime)*dble((1.-cdw*dTime))
           end do
        end do
     end do

  OPM2 = Ocn_POPMass(NSubs)
  Ocn_Degr = Ocn_Degr + OPM1 - OPM2

end subroutine Ocn_POPdegradation

!===========================================================
! Sedimentation process for POPs
!===========================================================
subroutine Ocn_POPSedimentation(NSubs, dTime)

  integer NSubs
  integer i, j, k,n
  real dTime, wsd
  real*8 OPM1,OPM2

  OPM1 = Ocn_POPMass(NSubs)

     n=sFormInd(Ocn,NSubs,1)
     wsd = WSed*(k_ph_rd(1)/(1+k_ph_rd(1)))         ! Temporarely defined. Needs to be replaced lately!!!
     do i = 1, Imax
        do j=1,Jmax
            Ocn_conc(i,j,1,n,oldtime) = Ocn_conc(i,j,1,n,oldtime)-Ocn_conc(i,j,1,n,curtime)*dble(wsd*dTime/Ocn_dz(1))
            Ocn_conc(i,j,1,n,curtime)=Ocn_conc(i,j,1,n,curtime)*dble((1.-wsd*dTime/Ocn_dz(1)))
        end do
     end do

     do i=1,Imax
        do j=1,Jmax
            do k=2,KMT(i,j)
                Ocn_conc(i,j,k,n,oldtime) = Ocn_conc(i,j,k,n,oldtime)-Ocn_conc(i,j,k,n,curtime)*dble(wsd*dTime/Ocn_dz(k))+&
                    Ocn_conc(i,j,k-1,n,curtime)*dble(wsd*dTime/Ocn_dz(k-1))       ! AG16022018 - was used Ocn_dz(k)
                Ocn_conc(i,j,k,n,curtime)=Ocn_conc(i,j,k,n,curtime)*dble((1.-wsd*dTime/Ocn_dz(k)))+&
                    Ocn_conc(i,j,k-1,n,curtime)*dble(wsd*dTime/Ocn_dz(k-1))       ! AG16022018 - was used Ocn_dz(k)
            end do
        end do
     end do

  OPM2 = Ocn_POPMass(NSubs)
  Ocn_Sedim = Ocn_Sedim + OPM1 - OPM2

end subroutine Ocn_POPSedimentation


!===========================================================
! POP mass calculation
!===========================================================
real(8) function Ocn_POPMass(NSubs)

  integer i, j,k,n, NSubs

    n=sFormInd(Ocn,NSubs,1)
    Ocn_POPMass = 0.
    do k = 1, Ocn_Kmax
        do j =  JMin, JMax
            do i = IMin, IMax
                Ocn_POPMass = Ocn_POPMass + (Ocn_conc(i,j,k,n,curtime)+Ocn_conc(i,j,k,n,oldtime)) &
                & *dble(TAREA(i,j)*Ocn_dz(k) /1000000./2.)
            end do
        end do
    end do

end function Ocn_POPMass


!===========================================================
! Exchange with ocean
!===========================================================
subroutine Ocn_POPExchPar(NSubs)

  integer i, j, j1, NSubs
  real OVE
  real Henry, TSurf, AbsVel, a1, a2, dzwm
  real RSea            ! Resistance of the first water layer
  real Rab, RT, Pxx, Zr, Ux, ZoM, Ux0, Lmo0, Ui, Uj, Uref

#ifdef DEBUG_MODE
      print *, '-> Entering Ocn_POPExchPar...'
#endif

    do j = Jmin, Jmax
        do i = minI(j), maxI(j)
            Ux0=Ufric(i,j,Period,toDay)                                         ! Mean friction velocity
            Lmo0=MOLength(i,j,Period,toDay)
            Pxx=Px(i,j,Period,toDay)                                            ! Surface pressure Px=Ps-Pt
            RT=RTemp(i,j,1,Period,toDay)                                        ! RT of air at the lowest sigma level
            ZoM=0.016*Ux0*Ux0/Ggrav+Nu/9.1*Ux0                                  ! Roughness
            Zr=RT/Ggrav*log((Pxx+Ptop)/(Sigma(1)*Pxx+Ptop))                     ! Height of the lowest sigma level
            Ui=(Uwind(i,j,1,Period,toDay)+Uwind(i-1,j,1,Period,toDay))/2.
            if(j>Jmin) then
              j1=j-1
              Uj=(Vwind(i,j,1,Period,toDay)+Vwind(iS(i,j,1),j1,1,Period,toDay))/2. ! Shift due to grid aggregation
            else
              Uj=Vwind(i,j,1,Period,toDay)
            endif
            Uref=sqrt(Ui*Ui+Uj*Uj)                                              ! Absolute value of vind velosity
            Ux=Karman*Uref/IphiM(Lmo0,Zr,ZoM)                                   ! Friction velocity
            Ux = max(Ux, 0.1)                                                   ! 09.04.2008 
            Rab = 0.74*2.5/Ux*(alog(Zr/ZoM)-PsiH(Zr/Lmo0)+PsiH(ZoM/Lmo0)+&
                &2*(1.5e-5/Dair(gSubsGroupInd(NSubs))/0.71)**(2./3.))
            Tsurf = TempSurf(i, j, Period)
            Henry = HT0(gSubsGroupInd(NSubs)) / (RUniv * TSurf) * exp(- HT(gSubsGroupInd(NSubs)) * (1 / TSurf - 1 / T0))
            AbsVel = sqrt(Uwind(i, j, 1, Period,toDay)**2 + Vwind(i, j, 1, Period,toDay)**2)
            a1 = 1.75 - 0.75 * exp( - 0.18 * AbsVel)
            a2 = min(1. - exp( - 0.01 * AbsVel), 0.8)
            dzwm = dzwm0 * exp(- 0.15 * AbsVel)
            RSea =  Henry / a1 / ((1 - a2) * DWater(gSubsGroupInd(NSubs)) / dzwm + a2 * Henry * dhf)
            OVE = 1. / (Rab + RSea + RSeaAdd)
            Ocn_VExch(i,j,Atm,1,NSubs) = OVE
            Ocn_VExch(i,j,Atm,2,NSubs) = OVE * Henry
        end do
    end do

#ifdef DEBUG_MODE
      print *, '-> Exit Ocn_POPExchPar'
#endif
end subroutine Ocn_POPExchPar

!...............................................................................
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

!...............................................................................
real function PsiH(x)
	real x
	if (x>0.) then		! stable
	  PsiH=-6.*x
	else				! non stable
	  PsiH=2.*alog(0.5*(1+Sqrt(1-9.*x)))
	endif
end function PsiH


end module Ocn_POP_General
#endif
