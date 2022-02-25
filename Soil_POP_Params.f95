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
! Soil compartment POP specific parameters
! Version:      1.0
! Modified:     24.11.2014
!=======================================================================================================
#ifdef G_POP
module Soil_POP_Params
use GeneralParams
implicit none

  real, parameter :: Facc= 0.3		              ! Accessible organic carbon fraction
  real, parameter :: FdoRel = 0.005               ! Relative fraction of organic carbon dissolved in soil solute
  real, parameter :: rhosol=1350.		          ! Soil density, kg/m3
  real, parameter :: Poros =0.5		              ! Soil porosity
  real, parameter :: wl=0.3			              ! Volumetric water fraction in soil
  real, allocatable :: FocEff(:,:)                ! Effective fraction of organic carbon in soil
  real, allocatable :: Fdo(:,:)                   ! Fraction of dissolved organics in soil solute
  real, parameter :: RateConst=2.2E-8             ! Exchange rate between solid phases
  real, parameter  :: AverDepth = 0.05		      ! Depth of soil concentration averaging, m
  real, allocatable :: CDSoil(:,:)                ! Degradation coefficients in soil layers
  real, parameter :: CDSoil1 = 8.8E-10	          ! Degradation rate for deep fraction
  real, allocatable :: PartL(:,:,:,:,:)           ! Inverse for liquid fraction in soil, last index is a pollutant number
  real, allocatable :: PartG(:,:,:,:,:)           ! Inverse for gas fraction in soil
  real, allocatable :: PartS(:,:,:,:,:)           ! Inverse for solid fraction in soil (accessible)
  real, allocatable :: POPVdAero(:,:,:,:)           ! Dry deposition velocity of aerosol phase
  real, parameter :: Difs = 6.e-12	              ! Effective diffusion coefficient due to bioturbation
  real CDSoilTop(MaxSubs)                         ! Degradation coefficients in soil
  real CDSoilBottom(MaxSubs)                     ! Degradation coefficients in soil
  real CDSoilBio(MaxSubs)                           ! Degradation coefficients in soil
  real Kow0(MaxSubs), KowT(MaxSubs)    ! Parameters of Kow temperature dependence
!  character*80 FocName

  contains

end module Soil_POP_Params
#endif