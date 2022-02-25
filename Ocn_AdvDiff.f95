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
! Module of advection and diffusion in the ocean.

! The algorythm and most of the subroutines are derived from
! POP (Parallel Ocean Program) model
! www.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#ifdef M_OCN
module OcnAdvDiff

  use Ocn_Params
  use GeneralParams

  implicit none

  character(800), private :: fileName, fullName

contains

!***********************************************************************

subroutine initialize_ocean_advdiff
#ifdef DEBUG_MODE
    print *, '-> Entering initialize_ocean_advdiff... '
#endif

   call horiz_grid_internal
   call init_grid
   call init_del2t
   call init_vertical_mix
   call init_prognostic

   kocean=0

#ifdef DEBUG_MODE
    print *, '<- Exit initialize_ocean_advdiff'
#endif

 end subroutine initialize_ocean_advdiff

!***********************************************************************

subroutine Ocn_adv_diff_step(NSubs)

 integer NSubs

 integer ::     k, n, nt

     n=sFormInd(Ocn,NSubs,1)
     nt=NumForm(Ocn)

     c2dtt = c2*dtOceanadvDiff

     call baroclinic_driver(n,DH)

     call baroclinic_correct_adjust(n,nt)

   if (avg_ts ) then     ! averaging step

            do k=2,Ocn_Kmax
               Ocn_conc(:,:,k,n,oldtime) =                &
                          p5*(Ocn_conc(:,:,k,n,oldtime) + &
                              Ocn_conc(:,:,k,n,curtime))
               Ocn_conc(:,:,k,n,curtime) =                &
                          p5*(Ocn_conc(:,:,k,n,curtime) + &
                              Ocn_conc(:,:,k,n,newtime))
            end do

               Ocn_conc(:,:,1,n,oldtime) =                   &
                   p5*((Ocn_dz(1) + DH(:,:,oldtime))*  &
                       Ocn_conc(:,:,1,n,oldtime) +           &
                       (Ocn_dz(1) + DH(:,:,curtime))*  &
                       Ocn_conc(:,:,1,n,curtime) )
               Ocn_conc(:,:,1,n,curtime) =                   &
                   p5*((Ocn_dz(1) + DH(:,:,curtime))*  &
                       Ocn_conc(:,:,1,n,curtime) +           &
                       (Ocn_dz(1) + DH(:,:,newtime))*  &
                       Ocn_conc(:,:,1,n,newtime) )

               Ocn_conc(:,:,1,n,oldtime) =                   &
               Ocn_conc(:,:,1,n,oldtime)/(Ocn_dz(1) +            &
                                      DH(:,:,oldtime))
               Ocn_conc(:,:,1,n,curtime) =                   &
               Ocn_conc(:,:,1,n,curtime)/(Ocn_dz(1) +            &
                                      DH(:,:,curtime))

   else  ! non-averaging step

      tmptime = oldtime
      oldtime = curtime
      curtime = newtime
      newtime = tmptime

   endif

    end subroutine Ocn_adv_diff_step

!***********************************************************************

   subroutine baroclinic_driver(n,DH)

   real, dimension(Imax,Jmax,3) :: &
      DH     ! change in surface height at T point

      real VELMAX

   integer  ::  &
      k,                &  ! indices for vertical level
      kp1,km1              ! level index for k+1, k-1 levels

   integer n,nt
   integer i,j,kk

   real, dimension(Imax,Jmax) :: &
      WTK                  ! vertical velocity at top of T box

      nt=NumForm(Ocn)

      do k = 1,Ocn_Kmax

         VDC=vdc_const

!-----------------------------------------------------------------------
      call Ocn_conc_update(k, n, nt, WTK,                             &
                               Ocn_conc (:,:,:,:,newtime), &
                               Ocn_conc (:,:,:,:,oldtime), &
                               Ocn_conc (:,:,:,:,oldtime), &
                               Ocn_conc (:,:,:,:,curtime), &
                               UVEL   (:,:,:,curtime), &
                               VVEL   (:,:,:,curtime), &
                               UVEL   (:,:,:,oldtime), &
                               VVEL   (:,:,:,oldtime), &
                               STF    (:,:,:), &
                               DH     (:,:,:))
!-----------------------------------------------------------------------

      enddo  ! k loop

 end subroutine baroclinic_driver
 
!***********************************************************************

 subroutine baroclinic_correct_adjust(n,nt)

   integer  ::  &
      k,                  &! vertical level index
      n , nt                     ! Ocn_conc index

                 where (KMT(:,:) > 0)
                     Ocn_conc(:,:,1,n,newtime) =               &
                                   Ocn_conc(:,:,1,n,newtime) - &
                                   Ocn_conc(:,:,1,n,oldtime)   &
                                   *(DH(:,:,newtime) -    &
                                     DH(:,:,oldtime))/    &
                                     Ocn_dz(1)
                  endwhere

               call impvmixt(Ocn_conc(:,:,:,:,newtime), &
                             Ocn_conc(:,:,:,:,oldtime), &
                             DH (:,:,    newtime), &
                             n, nt)

 end subroutine baroclinic_correct_adjust

!***********************************************************************

 subroutine Ocn_conc_update(k, n, nt,WTK, TNEW, TOLD, TMIX, TCUR,    &
                             UCUR, VCUR, UMIX, VMIX,                 &
                             STF_IN, DH_IN)

   integer, intent(in) :: k     ! depth level index

   integer :: &
      n             ! 
   integer :: &
      nt            ! 

     real*8, dimension(Imax,Jmax,Ocn_Kmax,nt), intent(in) :: &
      TCUR,                 &! TRACERs at current time level
      TOLD,                 &! TRACERs at old     time level
      TMIX                   ! TRACERs at mix     time level

   real, dimension(Imax,Jmax,Ocn_Kmax), intent(in) :: &
      UCUR, VCUR,           &! U,V  at current time
      UMIX, VMIX             ! U,V at mix time level

   real*8, dimension(Imax,Jmax,nt), intent(in) :: &
      STF_IN               ! surface TRACER fluxes

   real, dimension(Imax,Jmax,3), intent(in) :: &
      DH_IN                ! sfc height change at TRACER points

   real, dimension(Imax,Jmax), intent(inout) :: &
      WTK          ! on  input, vertical velocity at top    of T box
                   ! on output, vertical velocity at bottom of T box

   real*8, dimension(Imax,Jmax,Ocn_Kmax,nt), intent(inout) :: &
      TNEW                   ! TRACERs at new time level

      real*8, dimension(Imax,Jmax,nt) :: &
      FT,                &! sum of terms in dT/dt for the nth TRACER
      WORKN               ! work array used for various dT/dt terms

      integer kk,j,i

   FT    = c0

   call hdifft(k, n, nt, WORKN, TMIX)

 FT = FT + WORKN

   if (k == 1) WTK(:,:) = (DH_IN(:,:,curtime)-DH_IN(:,:,oldtime))/dtOceanadvDiff

   call advt(k,n,nt,WORKN,WTK,TCUR,UCUR,VCUR)

 FT = FT - WORKN   ! advt returns WORKN = +L(T)

   call vdifft(k, n, nt, WORKN, TOLD, STF_IN)

 FT = FT + WORKN


      do i=1,Imax
         do j=1,jmax
            TNEW(i,j,k,n)=c2dtt*FT(i,j,n)
            if (k > KMT(i,j)) TNEW(i,j,k,n)=TMIX(i,j,k,n)
         end do
      end do

 end subroutine Ocn_conc_update

!***********************************************************************

 subroutine vdifft(k, n, nt, VDTK, TOLD, STF)

   integer, intent(in) :: k   ! vertical level

 integer :: &
      n,                  &! dummy loop counter for TRACER number
      kp1,                &! k+1
      nt                            

    real*8, dimension(Imax,Jmax,Ocn_Kmax,nt), intent(in) :: &
      TOLD                ! TRACERs at old time level

   real*8, dimension(Imax,Jmax,nt), intent(in) :: &
      STF                 ! surface forcing for all TRACERs

  real*8, dimension(Imax,Jmax,nt), intent(out) :: &
      VDTK                ! returns VDIFF(TRACER(:,:,k,n))

    real, dimension(Imax,Jmax) :: &
      VTFB  ! vertical TRACER flux across bottom of T-box

   if (k  <  Ocn_Kmax) then
      kp1 = k + 1
   else
      kp1 = Ocn_Kmax
   endif

      if (k == 1) then
         VTF(:,:,n) = merge(STF(:,:,n), c08, KMT(:,:) >= 1)
      endif


          VTFB = merge(VDC(:,:,k)*                      &
                      (TOLD(:,:,k  ,n) - TOLD(:,:,kp1,n))*Ocn_dzwr(k) &
                      ,c08, KMT(:,:) > k)


         VDTK(:,:,n) = merge((VTF(:,:,n) - VTFB)*Ocn_dzr(k), &
                             c08, k <= KMT(:,:))

      VTF(:,:,n) = VTFB

end subroutine vdifft

!***********************************************************************

  subroutine init_del2t

   integer  ::  &
      i,j               ! dummy loop indices

   real, dimension (:,:), allocatable :: &
      WORK1         ! temporary work space

   allocate(WORK1 (Imax,Jmax))

     WORK1 = HTN(:,:)/HUW(:,:)

      DTN(:,:) = WORK1*TAREA_R(:,:)
      DTS(:,:) = eoshift(WORK1,dim=2,shift=-1)* &
                        TAREA_R(:,:)

      WORK1 = HTE(:,:)/HUS(:,:)

      DTE(:,:) = WORK1*TAREA_R(:,:)
      DTW(:,:) = eoshift(WORK1,dim=1,shift=-1)* &
                        TAREA_R(:,:)

   deallocate(WORK1)

  end subroutine init_del2t

!***********************************************************************

  subroutine init_vertical_mix

   integer  ::  &
      k                  ! vertical level index

   aidif = c1

      do k=1,Ocn_Kmax
         afac_t(k) = aidif*Ocn_dzwr(k)
      enddo

  end subroutine init_vertical_mix

!***********************************************************************

 subroutine hdifft(k,n, nt, HDTK,TMIX)

   integer, intent(in) :: k  ! depth level index

   real*8, dimension(Imax,Jmax,Ocn_Kmax,nt), intent(in) :: &
      TMIX                ! TRACERs at mix time level

   real*8, dimension(Imax,Jmax,nt), intent(out) ::  &
      HDTK                ! HDIFF(T) for TRACER n at level k

   integer  ::  &
      i,j,ip1,im1,jp1,jm1,n ,&! dummy TRACER index
      & nt

   real, dimension(Imax,Jmax) :: &
      CC,CN,CS,CE,CW! coeff of 5pt stencil for Del**2

       CN = merge(DTN(:,:), c0, (k <= KMTN(:,:)) .and. &
                                   (k <= KMT (:,:)))
      CS = merge(DTS(:,:), c0, (k <= KMTS(:,:)) .and. &
                                   (k <= KMT (:,:)))
      CE = merge(DTE(:,:), c0, (k <= KMTE(:,:)) .and. &
                                   (k <= KMT (:,:)))
      CW = merge(DTW(:,:), c0, (k <= KMTW(:,:)) .and. &
                                   (k <= KMT (:,:)))

   CC = -(CN + CS + CE + CW)  ! central coefficient

   HDTK = c08

   do j=1,Jmax
   do i=1,Imax

       ip1=i+1
       if (ip1.gt.Imax) ip1=1
       im1=i-1
       if (im1.lt.1) im1=Imax
       jp1=j+1
       if (jp1.gt.Jmax) jp1=Jmax
       jm1=j-1
       if (jm1.lt.1) jm1=1


      HDTK(i,j,n) = dble(ah*(CC(i,j)*TMIX(i  ,j  ,k,n) + &
                        CN(i,j)*TMIX(i  ,jp1,k,n) + &
                        CS(i,j)*TMIX(i  ,jm1,k,n) + &
                        CE(i,j)*TMIX(ip1,j  ,k,n) + &
                        CW(i,j)*TMIX(im1,j  ,k,n)))
   enddo
   enddo

 end subroutine hdifft

!***********************************************************************

  subroutine impvmixt(TNEW, TOLD, DH, n,nt)

   real*8, dimension(Imax,Jmax,Ocn_Kmax,nt), intent(inout) :: &
      TNEW         ! on input, contains right hand side
                   ! on output, contains updated TRACERs at new time

   real*8, dimension(Imax,Jmax,Ocn_Kmax,nt), intent(in) :: &
      TOLD         ! old TRACER to update with del(TRACER)

   real, dimension(Imax,Jmax), intent(in) :: &
      DH         ! surface pressure for use in determining
                   ! variable thickness surface layer

   integer :: &
      i,j,k,n,           &! dummy loop indices
      mt2 ,             & ! index for selecting TRACER coefficient
      nt

   real  ::          &
      a,b,c,d             ! various temporaries

   real, dimension(Ocn_Kmax) :: &
      e,f,               &! various work arrays
      hfac_t

   real, dimension(Imax,Jmax) :: &
      H1            ! factor containing full thickness of sfc layer

   do k=1,Ocn_Kmax
      hfac_t(k) = Ocn_dz(k)/c2dtt
   end do

       H1 = hfac_t(1) + DH/c2dtt

    do j=1,Jmax
      do i=1,Imax

        a = afac_t(1)*VDC(i,j,1)
         d = H1(i,j) + a
         e(1) = a/d
         b = H1(i,j)*e(1)
         f(1) = hfac_t(1)*TNEW(i,j,1,n)/d
         f(KMT(i,j)+1:Ocn_Kmax) = c0

         do k=2,KMT(i,j)

            c = a

            a = afac_t(k)*VDC(i,j,k)

            if (k == KMT(i,j)) then
               d = hfac_t(k)+b
            else
               d = hfac_t(k)+a+b
            endif

            e(k) = a/d
            b = (hfac_t(k) + b)*e(k)

            f(k) = (hfac_t(k)*TNEW(i,j,k,n) + c*f(k-1))/d

         end do

         !*** back substitution

         do k=KMT(i,j)-1,1,-1
            f(k) = f(k) + e(k)*f(k+1)
         end do

         do k = 1,Ocn_Kmax
            TNEW(i,j,k,n) = TOLD(i,j,k,n) + f(k)
         end do

      end do ! end of i-j loops
      end do

 end subroutine impvmixt

!***********************************************************************

 subroutine advt(k,n,nt,LTK,WTK,TRCR,UUU,VVV)

   integer, intent(in) :: &
      k  ! depth level index

      integer nt,kk

   real*8, dimension(Imax,Jmax,Ocn_Kmax,nt), intent(in) :: &
      TRCR               ! TRACERs at current time
   
   real*8, dimension(Imax,Jmax,Ocn_Kmax,nt) :: &
      TRCR1              ! TRACERs at current time

   real, dimension(Imax,Jmax,Ocn_Kmax), intent(in) :: &
      UUU, VVV           ! U,V at current time

   real, dimension(Imax,Jmax), intent(inout) :: &
      WTK            ! on  input flux velocity at top    of T box
                     ! on output flux velocity at bottom of T box
   real*8, dimension(Imax,Jmax,nt), intent(out) :: &
      LTK            ! returned as L(T) for nth TRACER at level k

   integer :: &
     i,j,im1,jm1,n          ! dummy loop indices

   real, dimension(Imax,Jmax) :: &
     UTE,UTW,VTN,VTS,  &    ! TRACER flux velocities across E,W,N,S faces
     WTKB                   ! vertical velocity at bottom of T box

   UTE = c0
   UTW = c0
   VTN = c0
   VTS = c0
   WTKB = c0 

   j=1
      do i=1,Imax


       im1=i-1
       if (im1.lt.1) im1=Imax
       jm1=j-1

         UTE(i,j) = p5*(UUU(i  ,j  ,k)*DYU(i  ,j))
         UTW(i,j) = p5*(UUU(im1,j  ,k)*DYU(im1,j))

         VTN(i,j) = p5*(VVV(i  ,j  ,k)*DXU(i  ,j) + &
                        VVV(im1,j  ,k)*DXU(im1,j))
         VTS(i,j) = 0.

      end do

      do j=2,Jmax
      do i=1,Imax

       im1=i-1
       if (im1.lt.1) im1=Imax
       jm1=j-1

         UTE(i,j) = p5*(UUU(i  ,j  ,k)*DYU(i  ,j) + &
                        UUU(i  ,jm1,k)*DYU(i  ,jm1))
         UTW(i,j) = p5*(UUU(im1,j  ,k)*DYU(im1,j) + &
                        UUU(im1,jm1,k)*DYU(im1,jm1))

         VTN(i,j) = p5*(VVV(i  ,j  ,k)*DXU(i  ,j) + &
                        VVV(im1,j  ,k)*DXU(im1,j))
         VTS(i,j) = p5*(VVV(i  ,jm1,k)*DXU(i  ,jm1) + &
                        VVV(im1,jm1,k)*DXU(im1,jm1))

      end do
      end do

      FC = p5*(VTN - VTS + UTE - UTW)*TAREA_R(:,:)
      WTKB = merge(WTK+c2Ocn_dz(k)*FC, c0, k < KMT(:,:))


      do kk=1,Ocn_Kmax
      TRCR1(:,:,kk,nt) = merge(TRCR(:,:,kk,nt), c08, kk <= KMT(:,:))           
      end do

    call advt_centered(k,n,nt,LTK,TRCR1,WTK,WTKB,UTE,VTN)

   WTK = WTKB

end subroutine advt

!***********************************************************************

 subroutine advt_centered(k,n,nt,LTK,TRCR,WTK,WTKB,UTE,VTN)

   integer, intent(in) :: k  ! depth level index

   integer nt

   real*8, dimension(Imax,Jmax,Ocn_Kmax,nt), intent(in) :: &
      TRCR                ! TRACERs at current time

   real, dimension(Imax,Jmax), intent(in) :: &
      UTE,VTN,         &! TRACER flux velocities across E,N faces
      WTK,             &! vert velocity at top of level k T box
      WTKB              ! vert velocity at bottom of level k T box

   real*8, dimension(Imax,Jmax,nt), intent(out) :: &
      LTK                 ! returned as L(T) for nth TRACER at level k

   integer :: &
      i,j,im1,jm1,ip1,jp1,n             ! TRACER loop index

      LTK = c0

      j=1
      do i=1,Imax

       ip1=i+1
       if (ip1.gt.Imax) ip1=1
       im1=i-1
       if (im1.lt.1) im1=Imax
       jp1=j+1

         LTK(i,j,n) = dble(p5*((VTN(i,j)+UTE(i,j)-UTE(im1,j))  &
                                      *TRCR(i  ,j  ,k,n) +           &
                          VTN(i  ,j  )*TRCR(i  ,jp1,k,n) -           &
                          UTE(i  ,j  )*TRCR(ip1,j  ,k,n) -           &
                          UTE(im1,j  )*TRCR(im1,j  ,k,n))*           &
                          TAREA_R(i,j))
      end do

      do j=2,Jmax-1
      do i=1,Imax

       ip1=i+1
       if (ip1.gt.Imax) ip1=1
       im1=i-1
       if (im1.lt.1) im1=Imax
       jp1=j+1
       jm1=j-1

         LTK(i,j,n) = dble(p5*((VTN(i,j)-VTN(i,jm1)+UTE(i,j)-UTE(im1,j))  &
                                      *TRCR(i  ,j  ,k,n) +           &
                          VTN(i  ,j  )*TRCR(i  ,jp1,k,n) -           &
                          VTN(i  ,jm1)*TRCR(i  ,jm1,k,n) +           &
                          UTE(i  ,j  )*TRCR(ip1,j  ,k,n) -           &
                          UTE(im1,j  )*TRCR(im1,j  ,k,n))*           &
                          TAREA_R(i,j))
      end do
      end do

      j=Jmax
       do i=1,Imax

       ip1=i+1
       if (ip1.gt.Imax) ip1=1
       im1=i-1
       if (im1.lt.1) im1=Imax
       jm1=j-1

         LTK(i,j,n) = dble(p5*((VTN(i,j)-VTN(i,jm1)+UTE(i,j)-UTE(im1,j))  &
                                      *TRCR(i  ,j  ,k,n) +           &
                          VTN(i  ,jm1)*TRCR(i  ,jm1,k,n) +           &
                          UTE(i  ,j  )*TRCR(ip1,j  ,k,n) -           &
                          UTE(im1,j  )*TRCR(im1,j  ,k,n))*           &
                          TAREA_R(i,j))
      end do

      if (k == 1) then
            LTK(:,:,n) = LTK(:,:,n) + &
                         dble(Ocn_dzr(k)*WTK)*TRCR(:,:,k,n)
      else
            LTK(:,:,n) = LTK(:,:,n) + dble(Ocn_dz2r(k)*WTK)*  &
                         (TRCR(:,:,k-1,n) + TRCR(:,:,k  ,n))
      endif

       if (k < Ocn_Kmax) then
            LTK(:,:,n) = LTK(:,:,n) - dble(Ocn_dz2r(k)*WTKB)* &
                         (TRCR(:,:,k,n) + TRCR(:,:,k+1,n))
      endif

 end subroutine advt_centered

!***********************************************************************

 subroutine horiz_grid_internal


   integer ::  i,j,ig,jg,jm1,n    ! dummy counters

   real ::&
      dlat, dlon,       &! lat/lon spacing for idealized grid
      lathalf            ! lat at T points

#if (REGTYPE==1)

      dlon = dXstep
      dlat = dYstep

      allocate (ULAT_G(Imax, 0:Jmax), &
                ULAT(Imax, 0:Jmax))

      do j = 0,Jmax
         ULAT_G(:,j)  = yOrig + j*dlat*pi180
      enddo

      do j=1,Jmax
            jg = j
            jm1 = jg - 1

            do i=1,Imax

                HTN(i,j) = dlon*Rearth*100.*pi180  ! convert to cm
                HTE(i,j) = dlat*Rearth*100.*pi180  ! convert to cm
                HUS(i,j) = dlon*Rearth*100.*pi180  ! convert to cm
                HUW(i,j) = dlat*Rearth*100.*pi180  ! convert to cm
                DYT(i,j) = dlat*Rearth*100.*pi180  ! convert to cm
                DYU(i,j) = dlat*Rearth*100.*pi180  ! convert to cm

               ig = i
               if (ig > 0 .and. jg > 0) then
                  ULAT(i,j) = ULAT_G(ig,jg)
                  HTN (i,j) = HTN(i,j)*cos(ULAT(i,j))
                  DXU (i,j) = HTN(i,j)
                  lathalf = (yOrig + (jg-p5)*dlat)*pi180
                  HUS (i,j) = HUS(i,j)*cos(lathalf)
                  DXT (i,j) = dlon*Rearth*100.*pi180 *  &
                              p5*(cos(ULAT_G(ig,jg )) + &
                                  cos(ULAT_G(ig,jm1)))
               else
                  ULAT(i,j) = c0
                  HTN (i,j) = c1 ! to prevent divide by zero
                  HUS (i,j) = c1 ! to prevent divide by zero
                  DXU (i,j) = c1 ! fixed up later
                  DXT (i,j) = c1 ! fixed up later
               endif

            end do
         enddo
      deallocate(ULAT_G,ULAT)

#else
      dlon = dXstep                   ! Added Sep 2016 V. Sh.
      dlat = dYstep

      allocate (ULAT_G(Imax, JMin - 1:Jmax), &
                ULAT(Imax, JMin - 1:Jmax))

!      do j = JMin,Jmax                             ! Creates NaN for ULAT(1,0) and DXT(1,1). Fixed
      do j = JMin-1,Jmax
         ULAT_G(:,j)  = yOrig + j*dlat*pi180
      enddo

      do j=JMin,Jmax
            jg = j
            jm1 = jg - 1

            do i=IMin,Imax

                HTN(i,j) = dlon*Rearth*100.*pi180  ! convert to cm
                HTE(i,j) = dlat*Rearth*100.*pi180  ! convert to cm
                HUS(i,j) = dlon*Rearth*100.*pi180  ! convert to cm
                HUW(i,j) = dlat*Rearth*100.*pi180  ! convert to cm
                DYT(i,j) = dlat*Rearth*100.*pi180  ! convert to cm
                DYU(i,j) = dlat*Rearth*100.*pi180  ! convert to cm

               ig = i
               if (ig > 0 .and. jg > 0) then
                  ULAT(i,j) = ULAT_G(ig,jg)
                  HTN (i,j) = HTN(i,j)*cos(ULAT(i,j))
                  DXU (i,j) = HTN(i,j)
                  lathalf = (yOrig + (jg-p5)*dlat)*pi180
                  HUS (i,j) = HUS(i,j)*cos(lathalf)
                  DXT (i,j) = dlon*Rearth*100.*pi180 *  &
                              p5*(cos(ULAT_G(ig,jg )) + &
                                  cos(ULAT_G(ig,jm1)))
               else
                  ULAT(i,j) = c0
                  HTN (i,j) = c1 ! to prevent divide by zero
                  HUS (i,j) = c1 ! to prevent divide by zero
                  DXU (i,j) = c1 ! fixed up later
                  DXT (i,j) = c1 ! fixed up later
               endif

            end do
         enddo
      deallocate(ULAT_G,ULAT)
!      print*,'OCEAN: this grid type is not supported yet'
#endif

end subroutine horiz_grid_internal

!***********************************************************************

 subroutine init_grid

   integer :: &
      i,j,k,ip1,jp1    ! dummy loop index variables


      TAREA(:,:) = DXT(:,:)*DYT(:,:)
      TAREA_R(:,:) = c1/TAREA(:,:)

!   call read_vert_grid                     ! Commented 24.11.2014. See GLEMOS_modifications.txt

   Ocn_dzw(0)  = p5*Ocn_dz(1)
   Ocn_dzw(Ocn_Kmax) = p5*Ocn_dz(Ocn_Kmax)
   Ocn_dzwr(0) = c1/Ocn_dzw(0)

   do k = 1,Ocn_Kmax-1
      Ocn_dzw(k) = p5*(Ocn_dz(k) + Ocn_dz(k+1))
   enddo

   do k = 1,Ocn_Kmax
      c2Ocn_dz(k) = c2*Ocn_dz(k)
      Ocn_dzr(k)  = c1/Ocn_dz(k)
      Ocn_dz2r(k) = c1/c2Ocn_dz(k)
      Ocn_dzwr(k) = c1/Ocn_dzw(k)
   enddo

#if (OCEAN_RUN==1)
      call read_topography
#else
      ! Reading topography is switched off to run without ocean (Sep 2016 V.E.Sh.)
#endif

   do j=1,Jmax
   do i=1,Imax
       ip1=i+1
       if (ip1.gt.Imax) ip1=1
       jp1=j+1
       if (jp1.gt.Jmax) jp1=Jmax

      KMU(i,j) = min(KMT(i,j ),KMT(ip1,j ), &
                       KMT(i,jp1),KMT(ip1,jp1))
   end do
   end do

   call landmasks

 end subroutine init_grid

!***********************************************************************

  subroutine init_prognostic

      integer i,j,k,ii,jj,kk

	  if(InitCond/='dump'.and.InitCond/='dum3'.and.InitCond/='dum5') then  

            oldtime = 1
            curtime = 2
            newtime = 3

          end if

      UVEL=0.
      VVEL=0.
      DH=0.
      STF=0.

 end subroutine init_prognostic

 !***********************************************************************

 subroutine read_topography

   integer i,j,ii,jj
   integer FileStat

      fileName=trim(TopoOcnName)//trim(GridCode)//'.dat'
      fullName=trim(GeoPath)//trim(fileName)
      open(10, file=trim(fullName), status = 'old', iostat=FileStat, action='read')
		if(FileStat>0) then
	      print '(/,"STOP: Cannot open Ocn topography file ''",a,"''",/)', trim(fullName)
	      stop
	    endif

      do i=1,Imax
      do j=1,Jmax

      read(10,*)ii,jj,KMT(ii,jj)

      end do
      end do

      close(10)

 end subroutine read_topography
 
!***********************************************************************

  subroutine landmasks

   KMTN = eoshift(KMT,dim=2,shift=+1)
   KMTS = eoshift(KMT,dim=2,shift=-1)
   KMTE = eoshift(KMT,dim=1,shift=+1)
   KMTW = eoshift(KMT,dim=1,shift=-1)

 end subroutine landmasks

!***********************************************************************

  subroutine read_vert_grid

   integer  :: &
      k                 ! vertical level index
   integer FileStat
  character sym

! 	  open(10, file=trim(ConfigPath)//'Ocn/'//trim(VertOcnName), status = 'old', iostat=FileStat, action='read')
!		if(FileStat>0) then
!	      print '(/,"STOP: Cannot open file ''",a,"''",/)', trim(VertOcnName)
!	      stop
!	    endif

 	  open(10, file=trim(ConfigPath)//'Ocn_Config.dat', status = 'old', iostat=FileStat, action='read')
		if(FileStat>0) then
	      print '(/,"STOP: Cannot open Vert grid file ''",a,"''",/)', 'Ocn_Config.dat'
	      stop
	    endif

      do k=1,7
      read(10,*) sym
      end do
      do k=1,Ocn_Kmax
      read(10,*) Ocn_dz(k)
      end do

      close(10)

      Water_upl_dz=Ocn_dz(1)/100.

 end subroutine read_vert_grid

!***********************************************************************

 end module OcnAdvDiff
#endif
