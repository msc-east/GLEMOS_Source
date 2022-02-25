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
! Module of the model model domain geometry characteristics
! Contains subroutines: 
! GridGeometry
! MeshVolume(i,j,k)
! Altitude(i,j,k)
! Sheight(i,j,k)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module Geometry

  use GeneralParams

  implicit none

contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine calculating the main grid characteristics
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GridGeometry

    integer i, j, k, Xscal
    real precTest

! Definition of horizontal grid
    cosdY=cos(dY)
    ctgdY2=1./tan(dY/2.)
    minI=Imin
    bminI=bImin
    do j=bJmin, bJmax
      if(j==bJmin.and.bJmin==Jmin.or.j==bJmax.and.bJmax==Jmax) cycle

      LatMesh(j)=dY*(real(j)-0.5)+yOrig
      tgY(j)=tan(LatMesh(j))   
      cosY(j)=cos(LatMesh(j))
      Xscal=1
      do while(Xscal*dXmin*Rearth*cosY(j)<dLmin)
!        if(aint(360./dXstep/real(Xscal*2))*dXstep*real(Xscal*2)<360.) exit
        if(aint(real(Imax)/real(Xscal*2))*real(Xscal*2)<real(Imax)) exit
        Xscal=Xscal*2
      enddo
        dX(j)=dXmin*real(Xscal)
        iGlob(j)=nint(2.*piNum/dX(j))
        maxI(j)=Imax/Xscal
      if(bImax>Imax) then
        bmaxI(j)=maxI(j)+1
      else
        bmaxI(j)=iGlob(j)
      endif

      do i=bminI(j), bmaxI(j)
        LongMesh(i,j)=dX(j)*(real(i)-0.5)+xOrig
      enddo
      do i=bImin, bImax
        LongMeshR(i)=dXmin*(real(i)-0.5)+xOrig
      enddo
    enddo

    do j=bJmin, bJmax
      do i=bminI(j), bmaxI(j)
        MeshArea(i,j)=2.*Rearth*Rearth*dX(j)*sin(dY/2.)*cosY(j)
      enddo
    enddo

    maxI(jSP)=maxI(jSP+1)
    maxI(jNP)=maxI(jNP-1)
    do j=bJmin+1, bJmax-1
      do i=minI(j), maxI(j)
        if(maxI(j)<maxI(j+1)) then
          iN(i,j,1)=2*i-1
          iN(i,j,2)=2*i
        elseif(maxI(j)>maxI(j+1)) then
          iN(i,j,1)=ceiling(real(i)/2.)
          iN(i,j,2)=ceiling(real(i)/2.)
        else
          iN(i,j,1)=i
          iN(i,j,2)=i
        endif
        if(maxI(j)<maxI(j-1)) then
          iS(i,j,1)=2*i-1
          iS(i,j,2)=2*i
        elseif(maxI(j)>maxI(j-1)) then
          iS(i,j,1)=ceiling(real(i)/2.)
          iS(i,j,2)=ceiling(real(i)/2.)
        else
          iS(i,j,1)=i
          iS(i,j,2)=i
        endif
      enddo
    enddo

    j=bJmin
    do i=minI(j), maxI(j)
      if(maxI(j)<maxI(j+1)) then
        iN(i,j,1)=2*i-1
        iN(i,j,2)=2*i
      else
        iN(i,j,1)=i
        iN(i,j,2)=i
      endif
    enddo
    j=bJmax
    do i=minI(j), maxI(j)
      if(maxI(j)<maxI(j-1)) then
        iS(i,j,1)=2*i-1
        iS(i,j,2)=2*i
      else
        iS(i,j,1)=i
        iS(i,j,2)=i
      endif
    enddo

!------------------ The Poles -------------------
if(bJmin==Jmin) then
    dX(jSP)=2.*piNum
    iGlob(jSP)=1
    maxI(jSP)=1
    bmaxI(jSP)=1
    LongMesh(iSP,jSP)=0.
    LatMesh(jSP)=-piNum/2.
    MeshArea(iSP,jSP)=2.*piNum*Rearth*Rearth*(1.-cosdY)
endif
if(bJmax==Jmax) then
    dX(jNP)=2.*piNum
    iGlob(jNP)=1
    maxI(jNP)=1
    bmaxI(jNP)=1
    LongMesh(iNP,jNP)=0.
    LatMesh(jNP)=piNum/2.
    MeshArea(iNP,jNP)=2.*piNum*Rearth*Rearth*(1.-cosdY)
endif
!------------------ The Poles -------------------

    do i=1, NumPer*2
      timePer(i)=real(i-1)*dTinput
    enddo

! Calculation of the machine smallest number
    Zero=tiny(precTest) 
    Eps=epsilon(precTest)

end subroutine GridGeometry


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function calculating 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GridAggreg(j,Xs,P,flag)

! flag=1 - averaging
! flag=2 - summing

    integer i, j, l, Xs, flag
    real P(Imin:Imax), Aver

    do i=minI(j), maxI(j)
      Aver=0.
      do l=1, Xs
        Aver=Aver+P(i*Xs-l+1)
      enddo
!      if(flag==1) P(i)=Aver/real(Xs)
      selectcase(flag)
      case(1)
        P(i)=Aver/real(Xs)
      case(2)
        P(i)=Aver
      case default
        print '(/,"STOP: Unknown flag of grid aggregation",/)', flag
        stop
      endselect
    enddo
 
end subroutine GridAggreg


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function calculating 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GridDisAggreg(P,flag)

! flag=1 - averaging
! flag=2 - total concerving

    integer i, j, l, iB, iE, Xs, flag
    real P(Imin:Imax,Jmin:Jmax), P1(Imin:Imax)

    do j=Jmin, Jmax
      Xs=Imax/maxI(j)
      if(Xs==1) cycle
      do i=minI(j), maxI(j)
        iB=(i-1)*Xs+1
        iE=iB+Xs-1
        do l=iB, iE
          selectcase(flag)
          case(1)
            P1(l)=P(i,j)
          case(2)
            P1(l)=P(i,j)/real(Xs)
          case default
            print '(/,"STOP: Unknown flag of grid disaggregation",/)', flag
            stop
          endselect
        enddo
      enddo
      P(Imin:Imax,j)=P1(Imin:Imax)
    enddo    
    
end subroutine GridDisAggreg


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function calculating cosine of the solar zenith angle
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real function SolarAngle(Long,Lat,Tm,Da,Mn,Yr,meanCos)

    integer Da, Mn, Yr
    real Long, Lat, Tm, Njd, A, B
    real Lm, Gm, EclLong, EclObl, sinDecl, cosDecl, Ha, meanCos

    Njd=DaysFrom2000(Da,Mn,Yr)
    Lm=(280.460+0.9856474*Njd)*pi180
    Gm=(357.528+0.9856003*Njd)*pi180
    EclLong=Lm+(1.915*sin(Gm)+0.02*sin(2*Gm))*pi180
    EclObl=(23.439-4.e-7*Njd)*pi180
    sinDecl=sin(EclObl)*sin(EclLong)
    cosDecl=sqrt(1.-sinDecl*sinDecl)
    Ha=2.*piNum*Tm/real(SecInDay)+Long-piNum
    A=sin(Lat)*sinDecl
    B=cos(Lat)*cosDecl

    SolarAngle=A+B*cos(Ha)
    if(A<=-B) then
      meanCos=0.
    elseif(A>=B) then
      meanCos=A
    else
      meanCos=(A*acos(-A/B)+sqrt(B*B-A*A))/PiNum         ! Average over a whole day
!      meanCos=A+sqrt(B*B-A*A)/acos(-A/B)                ! Average over daytime
    endif

contains

real function DaysFrom2000(Da,Mn,Yr)

    integer Da, Mn, Yr, i
    real Dj

    Dj=0.
    do i=1, Mn-1
      Dj=Dj+real(MonthDays(i))
    enddo
    Dj=Dj+real(Da)

    if(Yr>=1985.and.Yr<1989) then
      DaysFrom2000=-4384.5+(Yr-1988)*365+Dj
    elseif(Yr<1993) then
      DaysFrom2000=-4018.5+(Yr-1989)*365+Dj
    elseif(Yr<1997) then
      DaysFrom2000=-2557.5+(Yr-1993)*365+Dj
    elseif(Yr<2001) then
      DaysFrom2000=-1096.5+(Yr-1997)*365+Dj
    elseif(Yr<2005) then
      DaysFrom2000=364.5+(Yr-2001)*365+Dj
    elseif(Yr<2009) then
      DaysFrom2000=1825.5+(Yr-2005)*365+Dj
    elseif(Yr<2013) then
      DaysFrom2000=3286.5+(Yr-2009)*365+Dj
    elseif(Yr<2017) then
      DaysFrom2000=4747.5+(Yr-2013)*365+Dj
    elseif(Yr<2021) then
      DaysFrom2000=6208.5+(Yr-2017)*365+Dj
    elseif(Yr<2025) then
      DaysFrom2000=7669.5+(Yr-2021)*365+Dj
    elseif(Yr<2029) then
      DaysFrom2000=9130.5+(Yr-2025)*365+Dj
    elseif(Yr<2033) then
      DaysFrom2000=10591.5+(Yr-2029)*365+Dj
    elseif(Yr<2037) then
      DaysFrom2000=12052.5+(Yr-2033)*365+Dj
    elseif(Yr<2041) then
      DaysFrom2000=13513.5+(Yr-2037)*365+Dj
    elseif(Yr<2045) then
      DaysFrom2000=14974.5+(Yr-2041)*365+Dj
    elseif(Yr<2049) then
      DaysFrom2000=16435.5+(Yr-2045)*365+Dj
    elseif(Yr<2053) then
      DaysFrom2000=17896.5+(Yr-2049)*365+Dj
    else
      print '(/,"STOP: Year ",i4," is out of the range",/)', Yr
      stop
    endif

end function DaysFrom2000

end function SolarAngle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function calculating parameters of cubic spline
! Returns tabulated the first (Y1val) and the second (Y2val) derivatives
! flagL, flagR - flags indicating types of left and right boundary conditions:
! 0 - natural conditions (y''=0); 1 - specified first derivative (Y1L or Y1R)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SplineParams(N,Xval,Yval,flagL,Y1L,flagR,Y1R,Y1val,Y2val)

    integer N, flagL, flagR
    real Xval(N), Yval(N), Y2val(N), Y1val(N), Y1L, Y1R
    integer i, i
    real p, ratio, q(N), r(N), qR, rR, dX(N-1), dY(N-1)

! Defining left boundary condition
    if(flagL>0) then        ! Natural conditions
      q(1)=-0.5
      r(1)=(3./(Xval(2)-Xval(1)))*((Yval(2)-Yval(1))/(Xval(2)-Xval(1))-Y1L)
    else                    ! Specified first derivative
      q(1)=0.
      r(1)=0.
    endif

! Decomposition loop
    dX(1)=Xval(2)-Xval(1)
    dY(1)=Yval(2)-Yval(1)
    do i=2, N-1
      dX(i)=Xval(i+1)-Xval(i)
      dY(i)=Yval(i+1)-Yval(i)
      ratio=dX(i-1)/(dX(i)+dX(i-1))
      p=ratio*q(i-1)+2.
      q(i)=(ratio-1.)/p
      r(i)=(6.*(dY(i)/dX(i)-dY(i-1)/dX(i-1))/(dX(i)+dX(i-1))-ratio*r(i-1))/p
    enddo

! Defining left boundary condition
    if(flagR>0) then        ! Natural conditions
      qR=0.5
      rR=(3./(Xval(N)-Xval(N-1)))*(Y1R-(Yval(N)-Yval(N-1))/(Xval(N)-Xval(N-1)))
    else                    ! Specified first derivative
      qR=0.
      rR=0.
    endif
    Y2val(N)=(rR-qR*r(N-1))/(qR*q(N-1)+1.)

! Backsubstitution loop
    do i=N-1,1,-1
      Y2val(i)=q(i)*Y2val(i+1)+r(i)
      Y1val(i)=dY(i)/dX(i)-dX(i)*(2.*Y2val(i)+Y2val(i+1))/6.
    enddo 
    Y1val(N)=dY(N-1)/dX(N-1)+dX(N-1)*(Y2val(N-1)+2.*Y2val(N))/6.

end subroutine SplineParams

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function performing spline interpolation of a tabular function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SplineInterpol(N,Xval,Yval,Y2val,point,X)

    integer N, point
    real X, Xval(N), Yval(N), Y2val(N)
    integer i, iL, iR
    real A, B, dXval
    real SplineInterpol

! Deterninating spline interval
!    iL=1
!    iR=N
!    do while(iR-iL>1)
!      i=int((iR+iL)/2)
!      if(Xval(i)>X)then
!        iR=i
!      else
!        iL=i
!      endif
!    enddo

    iL=point
    iR=point+1

! Calculating spline value
    dXval=Xval(iR)-Xval(iL)
    A=(Xval(iR)-X)/dXval 
    B=(X-Xval(iL))/dXval
    SplineInterpol=A*Yval(iL)+B*Yval(iR)+((A*A*A-A)*Y2val(iL)+(B*B*B-B)*Y2val(iR))*dXval*dXval/6.

end function SplineInterpol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function performing linear interpolation of a tabular function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LinearInterpol(N,Xval,Yval,point,X)

    integer N, point
    real X, Xval(N), Yval(N)
    integer i, iL, iR
    real A, B, dXval
    real LinearInterpol

    iL=point
    iR=point+1

! Calculating spline value
    dXval=Xval(iR)-Xval(iL)
    A=(Xval(iR)-X)/dXval
    B=(X-Xval(iL))/dXval
    LinearInterpol=A*Yval(iL)+B*Yval(iR)

end function LinearInterpol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function calculating angular distance between 2 points with coordinates
! (lon1,lat1) and (lon2,lat2)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AngleDist(lon1,lat1,lon2,lat2)

    real lon1, lat1, lon2, lat2
    real dLon, cosR, AngleDist

    dLon=(lon1-lon2)*pi180
    cosR=sin(Lat1*pi180)*sin(Lat2*pi180)+cos(Lat1*pi180)*cos(Lat2*pi180)*cos(dLon)
    if(cosR>-1.and.cosR<1) then
      AngleDist=acos(cosR)/pi180
    elseif(cosR>=1) then
      AngleDist=0.
    elseif(cosR<=-1) then
      AngleDist=180.
    endif

end function AngleDist

end module Geometry
