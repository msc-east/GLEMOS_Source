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

!=======================================================================================================
! Ocean compartment parameters
! Version:      1.0
! Modified:     24.11.2014
!=======================================================================================================
#ifdef M_OCN
module Ocn_Params

use GeneralParams
use Ocn_POP_params

implicit none

  integer, parameter :: Ocn_Kmax=15

  real(8), allocatable :: Water_upl_conc(:,:,:)         ! Water_upl_conc(i,j,Form); i, j - horizontal cell number, Form - form number                      Form- number of pollutant form
  real(8), allocatable :: Water_upl_concDay(:,:,:)      ! Water_upl_concDay(i,j,p); counter for averaging, p - pollutant number
  real(8), allocatable :: Water_upl_concMonth(:,:,:)    ! Water_upl_concMonth(i,j,p); counter for averaging
  real(8), allocatable :: Water_upl_concYear(:,:,:)     ! Water_upl_concYear(i,j,p); counter for averaging
  real(8) Ocn_Degr                                ! Mass degraded in the ocean
  real(8) Ocn_Sedim                               ! Mass sedimented from the ocean
  real :: Ocn_TStep = 21600.                      ! Maximum permissible time step
  real, allocatable :: Ocn_VExch(:,:,:,:,:)       ! Ocn_VExch(i,j,m,Dir,p) - exchange velocity in the cell (i,j)

  real Water_upl_dz                               ! Thickness of the ocean layer
  real, allocatable :: Ocn_Frac(:,:)              ! Ocn_Frac(i,j) - water fraction in the cell (i,j)
  real, allocatable :: LC2Ocn(:)                  ! LC2Ocn(n) - fraction of land type n included to model water
  real dzwm0                                      ! Thickness of molecular diffusion layer in water
  real dhf                                        ! speed of foam degradation
  real RSeaAdd                                    ! additional resistance in UML
  real kdvUML                                     ! Additional vertical diffusion coefficient in Upper Mixed Layer
  real vdc_const                                  ! Vertical diffusion coefficient (cm/s2)
  real ah                                         ! Horizontal diffusion coefficient (cm/s2)
  real SedimCoeff                                 ! Coefficient for computation of sedimentation speed,1/m/s
  real RPart                                      ! Average radius of sea particles
  real WSed                                       ! Sedimentation velocity
!  logical :: Ocn_DailyOutput = .false.            ! If DailyOutput = .true. daily avreaged fields are written

  character*2 chTime(4)
  data chTime/'00','06','12','18'/

  real,  parameter :: &
    c0     =    0.0  ,&
    c1     =    1.0  ,&
    c2     =    2.0  ,&
    c4     =    4.0  ,&
    p5     =    0.5

  real*8,  parameter :: &
    c08    =    0.0

  integer kocean
  integer kkk
  integer iMonth,iDay,iHour
  integer iMonth1,iDay1,iHour1

  real TRintegr,TRemis,TRDegr,TRSedim
  real PeriodTime

  character*4 chYear

  real(8), allocatable :: Ocn_conc (:,:,:,:,:)      ! 3d TRACER fields at 3 time levels
  real(8), allocatable :: Ocn_ConcDay(:,:,:,:)      ! counter for averaging
  real(8), allocatable :: Ocn_ConcMonth(:,:,:,:)    ! counter for averaging
  real(8), allocatable :: Ocn_ConcYear(:,:,:,:)     ! counter for averaging

  real, dimension (:,:,:,:),allocatable :: &
      UVEL ,    &                                   ! 3d horizontalu-velocity at 3 time lvls
      VVEL                                          ! 3d horizontal v-velocity at 3 time lvls

  real, dimension(:,:,:),allocatable :: &
      UVEL0,UVEL1,     &
      VVEL0,VVEL1

  real, dimension(:,:),allocatable :: &
      HEIG0,HEIG1

  real, dimension(:,:),allocatable:: &
      PSURF    ! surface pressure

  real, dimension(:,:,:), allocatable :: &
      VDC                 ! diffusivity

  real, dimension(:,:,:),allocatable :: &
      DH

  real*4, dimension(:,:),allocatable ::&
      HUS, HUW            ,&                         ! cell widths on {S,W} sides of U cell
      TAREA               ,&                         ! area of T cell
      TAREA_R                                        ! reciprocal of area of T cell

  real*8, dimension(:,:),allocatable :: DXU, DYU       ! {x,y} spacing centered at U points
  real*8, dimension(:,:),allocatable :: DXT, DYT       ! {x,y} spacing centered at T points


  real*4, dimension (:,:), allocatable :: &
      DTN,DTS,DTE,DTW

  real*8, dimension(:,:,:),allocatable :: &
      VTF                                            ! vertical TRACER flux at top of T-box

  integer ::             &                           ! time indices for prognostic arrays:
      curtime,           &                           ! current time level  (n)
      newtime,           &                           ! next time level     (n+1)
      oldtime,           &                           ! previous time level (n-1)
      tmptime

  real :: &
     dtOceanAdvDiff,&                                ! time step at each level
     c2dtt                                           ! 

  real , dimension(:,:), allocatable :: &
      ULAT_G, ULAT                                   ! latitude of U points

  real, dimension(:,:), allocatable :: &
      HTN, HTE                                       ! cell widths on {N,E} sides of T cell

  logical    :: &
      avg_ts                                         !   an averaging timestep

  real, dimension(Ocn_Kmax) :: &
      Ocn_dz                ,&                       ! thickness of layer k
      c2Ocn_dz              ,&                       ! 2*Ocn_dz
      Ocn_dzr, Ocn_dz2r                              ! reciprocals of Ocn_dz, c2Ocn_dz

  real, dimension(0:Ocn_Kmax), public :: &
      Ocn_dzw, Ocn_dzwr                              ! midpoint of k to midpoint of k+1
                                                     ! and its reciprocal

  integer , dimension(:,:), allocatable ::&
      KMT,KMU                                        ! k index of deepest grid cell on T grid

  real*8, dimension(:,:,:),allocatable :: &
      STF                                            ! surface TRACER fluxes

  real ::          &
      aidif                                          ! time-centering parameter for implicit vmix

  real, dimension(Ocn_Kmax) :: &
      afac_t

  integer , dimension(:,:), allocatable:: &
      KMTN,KMTS,KMTE,KMTW                            ! KMT field at neighbor points

  real, dimension(:,:), allocatable :: &
     FC                                              ! local temp space


  contains

!======================================================================================
! Allocate ocean related arrays
!======================================================================================
  subroutine Ocn_allocation

  integer Ares, g

  allocate ( LC2Ocn(NumSurf), &
    & Water_upl_conc(IMin:IMax,JMin:JMax,NumForm(Ocn)), &
    & Water_upl_concDay(IMin:IMax,JMin:JMax,NumSubs), &
    & Water_upl_concMonth(IMin:IMax,JMin:JMax,NumSubs), &
    & Water_upl_concYear(IMin:IMax,JMin:JMax,NumSubs), &
    & Ocn_Frac(IMin:IMax,JMin:JMax), &
    & Ocn_VExch(IMin:IMax,JMin:JMax,NumMed,2,NumSubs),&
    & Ocn_ConcDay(IMin:IMax,JMin:JMax,Ocn_KMax,NumSubs), &
    & Ocn_ConcMonth(IMin:IMax,JMin:JMax,Ocn_KMax,NumSubs), &
    & Ocn_ConcYear(IMin:IMax,JMin:JMax,Ocn_KMax,NumSubs), &
    & DTN(Imax,Jmax),&
    & DTS(Imax,Jmax),&
    & DTE(Imax,Jmax),&
    & DTW(Imax,Jmax),&
    & VDC(Imax,Jmax,Ocn_Kmax),&
    & Ocn_conc(Imax,Jmax,Ocn_Kmax,NumForm(Ocn),3),&
    & UVEL(Imax,Jmax,Ocn_Kmax,3),VVEL(Imax,Jmax,Ocn_Kmax,3),&
    & UVEL0(Imax,Jmax,Ocn_Kmax),VVEL0(Imax,Jmax,Ocn_Kmax),&
    & UVEL1(Imax,Jmax,Ocn_Kmax),VVEL1(Imax,Jmax,Ocn_Kmax),&
    & HEIG0(Imax,Jmax),HEIG1(Imax,Jmax),&
    & PSURF(Imax,Jmax),DH(Imax,Jmax,3),&
    & HUS(Imax,Jmax),HUW(Imax,Jmax),TAREA(Imax,Jmax),TAREA_R(Imax,Jmax),&
    & DXU(Imax,Jmax),DYU(Imax,Jmax),DXT(Imax,Jmax),DYT(Imax,Jmax),&
    & HTN(Imax,Jmax),HTE(Imax,Jmax),KMT(Imax,Jmax),KMU(Imax,Jmax),&
    & KMTN(Imax,Jmax),KMTS(Imax,Jmax),KMTE(Imax,Jmax),KMTW(Imax,Jmax),&
    & VTF(Imax,Jmax,NumForm(Ocn)),STF(Imax,Jmax,NumForm(Ocn)),FC(Imax,Jmax)&
    , stat = ARes)
	if (ARes /= 0) then
	  print*, 'Allocation failed'
	  stop
	end if

#ifdef G_POP
                call Ocn_POPallocation
#endif

end subroutine Ocn_allocation

!======================================================================================
! Deallocate ocean related arrays
!======================================================================================
subroutine Ocn_deallocation

  integer g

  deallocate (LC2Ocn, Water_upl_conc, Water_upl_concDay, Water_upl_concMonth, Water_upl_concYear, Ocn_Frac, Ocn_VExch,&
  &DTN,DTS,DTE,DTW,VDC,Ocn_conc,Ocn_ConcDay,Ocn_ConcMonth,Ocn_ConcYear,UVEL,VVEL,UVEL0,VVEL0,UVEL1,VVEL1,HEIG0,HEIG1,PSURF,DH,&
  &HUS,HUW,TAREA,TAREA_R,DXU,DYU,DXT,DYT,VTF,HTN,HTE,KMT,KMU,STF,FC,KMTN,KMTE,KMTS,KMTW)

#ifdef G_POP
                call Ocn_POPdeallocation
#endif

end subroutine Ocn_deallocation

!======================================================================================
! Allocate soil related arrays (POPs)
!======================================================================================
subroutine Ocn_POPallocation

end subroutine Ocn_POPallocation

!======================================================================================
! Deallocate soil related arrays (POPs)
!======================================================================================
subroutine Ocn_POPdeallocation

end subroutine Ocn_POPdeallocation

end Module Ocn_Params
#endif