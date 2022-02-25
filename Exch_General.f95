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
! Contains subroutines: 
! Media_Exchange
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#if G_HM .or. G_HG

module Exch_General 

  use GeneralParams
  use Atm_Output
  use Exch_Params
#ifdef G_HG
  use Atm_Hg_Params
#endif
#ifdef G_HM
  use Atm_HM_Params
#endif

  implicit none

  character(800), private :: fileName, fullName


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
!---------------------------------------------------------------------------------
    case('Final')
      call Exch_MemDealloc
    endselect

end subroutine Exch_Memory


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calling process procedures of different media
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Media_Exchange

    integer, save :: order=1

    ReEmisFlux=0.

    if(order>0) then
      call DryDepos
      call WetDepos
      order=-order
    else
      call WetDepos
      call DryDepos
      order=-order
    endif

end subroutine Media_Exchange


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine supplying information for the calculation run
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Exch_InputData(Prd)

    character(*) Prd
    character(100) Note

    selectcase(Prd)
!---------------------------------------------------------------------------------
    case('Initial')

! Reading dry deposition params (coeffs. for ABL, z0 for seasons and l-uses etc.): 
      call DryDepParams
      call LogFile('Reading dry deposition parameters',0)

!---------------------------------------------------------------------------------
    endselect

end subroutine Exch_InputData


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine reading adjusting parameters for ABL, some common geophys. data etc.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DryDepParams

    integer i, FileStat, lenType, Surf, s, SurfInd(MaxSurf), lenStr, ios
    real zLU(MaxSurf), dLU(MaxSurf), hLU(MaxSurf)
    character(300) strRead, strPar
    character(40) strType
    character(20) IDsurf(MaxSurf), blank
    
    fullName=trim(GeoPath)//trim(ComGeoName)
    open(2, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fileName)
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
        case('Beta for momentum')
          read (StrPar, *) Bm
!------------------------------------------------------------------
        case('Beta for heat')
          read (StrPar, *) Bh
!------------------------------------------------------------------
        case('Gamma for momentum')
          read (StrPar, *) Gm
!------------------------------------------------------------------
        case('Gamma for heat')
          read (StrPar, *) Gh
!------------------------------------------------------------------
        case('R of broken water surface, s/m')
          read (StrPar, *) Rbrok
!------------------------------------------------------------------
        case('Fog droplet diameter, m')
          read (StrPar, *) Dfog
!------------------------------------------------------------------
      end select
    enddo
    close(2)

! Reading roughness, displacement and top canopy heights
    fullName=trim(GeoPath)//trim(RoughName)
    open(3, file=trim(fullName), status='old', iostat=FileStat, action='read')
    if(FileStat>0) then
      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(fullName)
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

end subroutine DryDepParams


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating dry deposition 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DryDepos

    integer i, j, Ind, Form, Surf, Src, Subs, Gas
    real meshV, expConc, expDep, vDep
    real(8) Conc, MassDep, aD(MaxMatr)
    real Teff, cosSol, cosMean, Lcover(MaxSurf), dT, tDay, solR, solRmonth
    real(8) Tage, Tred, expAge, expRed, MassFresh, MassAge, MassRed
    real splRow(NumPer*2)
 
! Definition of elemental gaseous form
!     do Form=1, NumForm(Atm)
!       Subs=FormSubs(Atm,Form)
!       if(FormID(Atm,Subs,Form)=='0_gas') then
!         Gas=Form
!         exit
!       endif
!     enddo
  
    dT=Tstep(Atm)
    tDay=DayTime
    do j=Jmin, Jmax
      do i=minI(j), maxI(j)

        Teff=RTcurr(i,j,1)/(Sigma(1)+Ptop/PxCurr(i,j))
        meshV=MeshVolume(i,j,1)
        Lcover(1:NumSurf)=LandCover(i,j,1:NumSurf)
        cosSol=max(SolarAngle(LongMesh(i,j),LatMesh(j),tDay,Day,Month,Year,cosMean),0.)
        do Ind=1, DryDepNum
          Form=DryDepInd(Ind)
          Subs=FormSubs(Atm,Form)

!***************** Matrix calculations *****************
#if RTYPE==2
          aD(1:NumSrc)=Atm_Contrib(i,j,1,Form,1:NumSrc)
#else
!*******************************************************
          aD(1)=1.
#endif
          Conc=Atm_Conc(i,j,1,Form)
          expConc=0.
          do Surf=1, NumSurf
            if(Lcover(Surf)<=0.) cycle

            splRow(1:NumPer)=Vd(i,j,1:NumPer,Ind,Surf,toDay)
            splRow(NumPer+1:NumPer*2)=Vd(i,j,1:NumPer,Ind,Surf,toMor)
            vDep=LinearInterpol(NumPer*2,timePer,splRow,Period,DayTime)

            if(FormID(Atm,Subs,Form)=='0_gas') vDep=vDep*cosSol                !????????????????
            expDep=exp(-vDep*dT*Ggrav/Teff/dS(1))
            
! Accumulation of the dry deposited mass
            MassDep=Conc*(1.-expDep)*Lcover(Surf)*meshV
            expConc=expConc+Lcover(Surf)*expDep

            do Src=1, NumSrc
!%%%%%%%%%%%%%%%%% Re-emission %%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef G_HG
              MassFresh=FreshDep(i,j,Surf,Src)
              if(SurfType(Surf)=='Snow_Ice'.or.SurfType(Surf)/='Water'.and.SnowDepth(i,j,Period)>0.01.or.&    ! 1 cm WRF
                                    &SurfType(Surf)=='Water'.and.SeaIce(i,j,Period)>0.1) then
                Tage=10.*real(SecInDay)       ! 10 days
                Tred=1.*real(SecInDay)        ! 1 day
                expAge=exp(-dT/Tage)
                expRed=exp(-dT/Tred)
                MassFresh=MassFresh+MassDep*aD(Src)
                MassAge=MassFresh*(1.-expAge)
                MassRed=MassFresh*expAge*(1.-expRed)
                MassFresh=MassFresh*expAge*expRed

                FreshDep(i,j,Surf,Src)=MassFresh
                ReEmisFlux(i,j,1,Src)=ReEmisFlux(i,j,1,Src)+MassRed           ! kg
DryRemMonth(i,j,Hg0,Surf,Src)=DryRemMonth(i,j,Hg0,Surf,Src)+MassRed
              elseif(MassFresh>0.) then
                FreshDep(i,j,Surf,Src)=0.
              endif
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!            do Src=1, NumSrc
!              DryDepMonth(i,j,Form,Surf,Src)=DryDepMonth(i,j,Form,Surf,Src)+MassDep*aD(Src)
              DryDepDay(i,j,Form,Surf,Src)=DryDepDay(i,j,Form,Surf,Src)+MassDep*aD(Src)
            enddo
            MassDryDep(Form)=MassDryDep(Form)+MassDep
          enddo

! Correction of the pollutant air concentration
          Atm_Conc(i,j,1,Form)=Conc*expConc
        enddo
      enddo
    enddo

end subroutine DryDepos


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating in-cloud and below-cloud wet deposition 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WetDepos

    integer i, j, k, Ind, Src, Surf, Form, Subs, Gas
    real Iprec(Atm_Kmax), Ip, dIprec, dT
    real Lw, Fw, meshV, mRatio, Lcover(MaxSurf)
    real Sfrac, Keff, expDep, WashOut(MaxForm), Wout, solR, solRmonth
    real(8) Conc, dConc, dCf(MaxForm), dCs(MaxMatr), dMass, aW(MaxMatr), MassDep(MaxForm,MaxMatr), PrecM
    real(8) Tage, Tred, expAge, expRed, MassFresh, MassAge, MassRed


! Definition of particulate and elemental gaseous forms
!      do Form=1, NumForm(Atm)
!        if(FormID(Atm,Subs,Form)=='0_gas') then
!          Gas=Form
!          exit
!       endif
!      enddo
   
    dT=Tstep(Atm)
    do j=Jmin, Jmax
      do i=minI(j), maxI(j)
        Iprec(1:Atm_Kmax)=PrecRainConv(i,j,1:Atm_Kmax,Period)+PrecRainStrat(i,j,1:Atm_Kmax,Period)
        Lcover(1:NumSurf)=LandCover(i,j,1:NumSurf)
        MassDep=0.
        do k=Atm_Kmax, 1, -1
          Ip=Iprec(k)
          if(k==Atm_Kmax) then
            dIprec=Iprec(k)
          else
            dIprec=Iprec(k)-Iprec(k+1)
          endif
          if(Ip==0..and.dIprec==0.) cycle

          Lw=LiqCont(i,j,k,Period)
          Fw=FrozCont(i,j,k,Period)
          meshV=MeshVolume(i,j,k)

          if(Lw+Fw>Lw0) then                    ! In-cloud scavenging
            if(Ip>0.) then
              Sfrac=1.                          ! Stratiform precipitation
              expDep=exp(-Ain*(Ip*3.6e6/Sfrac)**Bin*dT)
!              Keff=Lw/Sfrac/(Lw/Sfrac+Aeff)
              Keff=(Lw+Fw)/Sfrac/((Lw+Fw)/Sfrac+Aeff)
              WashOut(1:NumForm(Atm))=Sfrac*(1.-expDep)
              WashOut(Part(1:Npart))=Sfrac*(1.-expDep)*Keff*(1.-Asol)/(1.-Asol*Keff)
            else
              WashOut(1:NumForm(Atm))=0.
            endif

            do Ind=1, ICScvNum
              Form=ICScvInd(Ind)
              Wout=WashOut(Form)

! Correction of the pollutant concentration
              Conc=Atm_Conc(i,j,k,Form)
              Atm_Conc(i,j,k,Form)=Conc*(1.-Wout)

              dMass=Conc*Wout*meshV
!***************** Matrix calculations *****************
#if RTYPE==2
              do Src=1, NumSrc
                MassDep(Form,Src)=MassDep(Form,Src)+dMass*Atm_Contrib(i,j,k,Form,Src)
              enddo
#else
!*******************************************************
              MassDep(Form,1)=MassDep(Form,1)+dMass
#endif
            enddo
          else                                    ! Below-cloud scavenging
            if(Ip>0.) then                        
              Sfrac=1.                          ! Stratiform precipitation
              expDep=exp(-Abelow*(Ip*3.6e6/Sfrac)**Bbelow*dT)
              WashOut(1:NumForm(Atm))=Sfrac*(1.-expDep)
            else
              WashOut(1:NumForm(Atm))=0.
            endif

            do Ind=1, BCScvNum
              Form=BCScvInd(Ind)
              Wout=WashOut(Form)

! Correction of the pollutant concentration
              Conc=Atm_Conc(i,j,k,Form)
              Atm_Conc(i,j,k,Form)=Conc*(1.-Wout)

              dMass=Conc*Wout*meshV

!***************** Matrix calculations *****************
#if RTYPE==2
              do Src=1, NumSrc
                MassDep(Form,Src)=MassDep(Form,Src)+dMass*Atm_Contrib(i,j,k,Form,Src)
              enddo
#else
!*******************************************************
              MassDep(Form,1)=MassDep(Form,1)+dMass
#endif
            enddo
          endif

! Evaporation of precipitation
          if(dIprec<0..and.maxval(MassDep(1:NumForm(Atm),1:NumSrc))>Zero) then    
            mRatio=abs(dIprec)/Iprec(k+1)

            dCf=0.
            dConc=0.
            do Src=1, NumSrc
              dCf(1:NumForm(Atm))=dCf(1:NumForm(Atm))+MassDep(1:NumForm(Atm),Src)*mRatio/meshV
              dCs(Src)=sum(MassDep(1:NumForm(Atm),Src))*mRatio/meshV
              dConc=dConc+dCs(Src)
            enddo
            
            do Form=1, Npart
              Conc=Atm_Conc(i,j,k,Part(Form))
              Atm_Conc(i,j,k,Part(Form))=Conc+dConc/Npart
              MassChemEx(Part(Form))=MassChemEx(Part(Form))+dConc/Npart*meshV

!***************** Matrix calculations *****************
#if RTYPE==2
              aW=0.
              do Src=1, NumSrc
                aW(Src)=Atm_Contrib(i,j,k,Part(Form),Src)*Conc+dCs(Src)
              enddo
              aW=aW/(Conc+dConc)
              Atm_Contrib(i,j,k,Part(Form),1:NumSrc)=aW(1:NumSrc)
#endif
!*******************************************************
            enddo
            MassChemEx(1:NumForm(Atm))=MassChemEx(1:NumForm(Atm))-dCf(1:NumForm(Atm))*meshV
            MassDep(1:NumForm(Atm),1:NumSrc)=MassDep(1:NumForm(Atm),1:NumSrc)*(1.-mRatio)
          endif
        enddo

! Accumulation of the wet deposited mass
        do Form=1, NumForm(Atm)
          do Src=1, NumSrc
            if(MassDep(Form,Src)<=0.) cycle

            do Surf=1, NumSurf
              if(Lcover(Surf)<=0.) cycle

!%%%%%%%%%%%%%%%%% Re-emission %%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef G_HG
              MassFresh=FreshDep(i,j,Surf,Src)
              if(SurfType(Surf)=='Snow_Ice'.or.SurfType(Surf)/='Water'.and.SnowDepth(i,j,Period)>0.01.or.&    ! 1 cm WRF
                                    &SurfType(Surf)=='Water'.and.SeaIce(i,j,Period)>0.1) then
                Tage=10.*real(SecInDay)       ! 10 days
                Tred=1.*real(SecInDay)        ! 1 day
                expAge=exp(-dT/Tage)
                expRed=exp(-dT/Tred)
                MassFresh=MassFresh+MassDep(Form,Src)*Lcover(Surf)
                MassAge=MassFresh*(1.-expAge)
                MassRed=MassFresh*expAge*(1.-expRed)
                MassFresh=MassFresh*expAge*expRed

                FreshDep(i,j,Surf,Src)=MassFresh
                ReEmisFlux(i,j,1,Src)=ReEmisFlux(i,j,1,Src)+MassRed
DryRemMonth(i,j,Hg0,Surf,Src)=DryRemMonth(i,j,Hg0,Surf,Src)+MassRed
              elseif(MassFresh>0.) then
                FreshDep(i,j,Surf,Src)=0.
              endif
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!              WetDepMonth(i,j,Form,Surf,Src)=WetDepMonth(i,j,Form,Surf,Src)&
!                    &+MassDep(Form,Src)*Lcover(Surf)
              WetDepDay(i,j,Form,Surf,Src)=WetDepDay(i,j,Form,Surf,Src)+MassDep(Form,Src)*Lcover(Surf)
            enddo
            MassWetDep(Form)=MassWetDep(Form)+MassDep(Form,Src)
          enddo
        enddo

! Accumulation of precipitation
!        PrecipMonth(i,j)=PrecipMonth(i,j)+Iprec(1)*dT
        PrecipDay(i,j)=PrecipDay(i,j)+Iprec(1)*dT
      enddo
    enddo

end subroutine WetDepos

end module Exch_General

#endif
