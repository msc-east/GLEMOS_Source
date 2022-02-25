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

module Atm_NcfOutput

    use GeneralParams
    use Atm_Params
    use Exch_Params
    use netcdf
    use Geometry

    implicit none

!-----------------------------------------------------------------
! netCDF attributes
    character (*), parameter :: UNITS      = "units"
    character (*), parameter :: CALENDAR   = "calendar"
    character (*), parameter :: AXIS       = "axis"
    character (*), parameter :: LNG_NAME   = "long_name"
    character (*), parameter :: STD_NAME   = "standard_name"
    character (*), parameter :: TERMS      = "formula_terms"
!-----------------------------------------------------------------
! Coordinate parameters
    character (*), parameter :: TIME_NAME  = "time"
    character (*), parameter :: TIME_CLD   = "julian"
    character (*), parameter :: TIME_AXS   = "T"
    character (*), parameter :: TIME_STD   = "time"
    character (*), parameter :: TIME_LNG   = "time"

    character (*), parameter :: LAT_NAME   = "lat"
    character (*), parameter :: LAT_UNITS  = "degrees_north"
    character (*), parameter :: LAT_AXS    = "Y"
    character (*), parameter :: LAT_STD    = "latitude"
    character (*), parameter :: LAT_LNG    = "latitude"

    character (*), parameter :: LON_NAME   = "lon"
    character (*), parameter :: LON_UNITS  = "degrees_east"
    character (*), parameter :: LON_AXS    = "X"
    character (*), parameter :: LON_STD    = "longitude"
    character (*), parameter :: LON_LNG    = "longitude"

    character (*), parameter :: LEV_NAME   = "lev"
    character (*), parameter :: LEV_UNITS  = "1"
    character (*), parameter :: LEV_AXS    = "Z"
    character (*), parameter :: LEV_STD    = "atmosphere_standard_sigma_coordinate"
    character (*), parameter :: LEV_LNG    = "Standard sigma pressure coordinate"
    character (*), parameter :: LEV_TERM   = "sigma: lev ptop: ptop ps: ps"

    character (*), parameter :: PTOP_NAME   = "ptop"
    character (*), parameter :: PTOP_UNITS  = "Pa"
    character (*), parameter :: PTOP_STD    = "atmosphere_standard_sigma_coordinate"
    character (*), parameter :: PTOP_LNG    = "Pressure at the top of the model"

!-----------------------------------------------------------------
! Concentrations and deposition parameters
    character (*), parameter :: SRF_NAME   = "srf"
    character (*), parameter :: SRF_UNITS  = "surface"
    character (*), parameter :: SRF_STD    = "srf"
    character (*), parameter :: SRF_LNG    = "underlying surface"

    character (*), parameter :: SRC_NAME   = "src"
    character (*), parameter :: SRC_UNITS  = "source"
    character (*), parameter :: SRC_STD    = "src"
    character (*), parameter :: SRC_LNG    = "emission source"

    character (*), parameter :: SRC_ARR   = "src_nm"
    character (*), parameter :: SRF_ARR   = "srf_nm"

!-----------------------------------------------------------------
! Concentrations and deposition parameters for stations
    character (*), parameter :: STA_NAME   = "air_stat"
    character (*), parameter :: STA_UNITS  = "station"
    character (*), parameter :: STA_STD    = "air_stat"
    character (*), parameter :: STA_LNG    = "air conc monitoring stations"

    character (*), parameter :: STP_NAME   = "prec_stat"
    character (*), parameter :: STP_UNITS  = "station"
    character (*), parameter :: STP_STD    = "prec_stat"
    character (*), parameter :: STP_LNG    = "prec conc monitoring stations"

    character (*), parameter :: STPR_NAME   = "precip_stat"
    character (*), parameter :: STPR_UNITS  = "station"
    character (*), parameter :: STPR_STD    = "precip_stat"
    character (*), parameter :: STPR_LNG    = "precip monitoring stations"

    character (*), parameter :: STA_ARR   = "air_stat_nm"
    character (*), parameter :: STP_ARR   = "prec_stat_nm"
    
!-----------------------------------------------------------------
! Meteo parameters
    character (*), parameter :: PSUR_NAME   = "ps"
    character (*), parameter :: PSUR_UNITS  = "Pa"
    character (*), parameter :: PSUR_STD    = "surface_air_pressure"
    character (*), parameter :: PSUR_LNG    = "surface air pressure"

    character (*), parameter :: PREC_NAME   = "precip"
    character (*), parameter :: PREC_UNITS  = "m"
    character (*), parameter :: PREC_STD    = "precipitation"
    character (*), parameter :: PREC_LNG    = "precipitation amount"

    real, parameter :: P0=1.e3           ! Pa

!-----------------------------------------------------------------
! Output time periods
    integer, parameter  :: TM_HOURLY = 1
    integer, parameter  :: TM_6HOURLY = 2
    integer, parameter  :: TM_DAILY = 3
    integer, parameter  :: TM_MONTHLY = 4
    integer, parameter  :: TM_YEARLY = 5
    
    integer, private	:: times_6h(4)
    data times_6h/0,6,12,18/

!-----------------------------------------------------------------
! netCDF variables
    integer, parameter  :: NCF_VARS_NO  = 7
    integer, parameter  :: AC_3D    = 1
    integer, parameter  :: MR_3D    = 2

    integer, parameter  :: DD_2D    = 3
    integer, parameter  :: WD_2D    = 4
    integer, parameter  :: RE_2D    = 5

    integer, parameter  :: AC_ST    = 6
    integer, parameter  :: WD_ST    = 7

    logical, private         :: fld_ncf_out(5,NCF_VARS_NO)
    character(80), private   :: v_hdr(5,NCF_VARS_NO)
    character(200), private  :: ncfPath(5,NCF_VARS_NO)

    real levs(Atm_KMax)
    character(80), private  :: fld_unt(NCF_VARS_NO)
    character(80), private  :: fld_nm(NCF_VARS_NO)
    character(20), private  :: v_nm(NCF_VARS_NO)
    integer, private        :: v_id(NCF_VARS_NO)
    integer, private        :: v_bc, psid_bc
    data fld_nm / &
                    &'air concentrations',&
                    &'mixing ratio',&
                    &'dry deposition flux',&
                    &'wet deposition flux',&
                    &'re-emission flux',&
                    &'air concentrations at monitoring stations',&
                    &'wet deposition flux at monitoring stations'&
                    & /
    data v_nm / &
                    &'air_conc',&
                    &'mix_ratio',&
                    &'dry_dep_flux',&
                    &'wet_dep_flux',&
                    &'re-emis_flux',&
                    &'air_conc_station',&
                    &'wet_dep_flux_station'&
                    & /
    data fld_unt / 'kg m-3','kg kg-1','kg','kg','kg','kg m-3','kg' /

    real, private, allocatable  :: var3d(:,:,:,:,:,:)
    real, private, allocatable  :: var2ds(:,:,:,:,:,:)
    real, private, allocatable  :: var2dlc(:,:,:)
    integer, private, allocatable   :: src_ind(:)
    integer, private, allocatable   :: srf_ind(:)
    real, private               :: var2d(Imin:Imax, Jmin:Jmax, 2)

!-----------------------------------------------------------------
! stations variables
    character(30), allocatable   :: sta_nm(:), stp_nm(:)
    real, allocatable           :: sta_mod(:,:,:), stp_mod(:,:,:), stpr_mod(:)

    integer     time_st_dimid, src_st_dimid, src_len_dimid, src_num_dimid, stn_len_dimid, sta_num_dimid, stp_num_dimid
    integer     dimid_sta(3), dimid_stp(3), sta_dimid, stp_dimid
    integer     src_st_id, time_st_id, sta_id, stp_id, stpr_id
    integer     src_nm_id, sta_nm_id, stp_nm_id
    integer     var_st_id(2,MaxForm)

!-----------------------------------------------------------------
! work variables
    character(4), private   ::  YearNum
    character(2), private   ::  MonthNum, DayNum, PerNum
    character(800), private ::  fileName, fullName
    character(40), private  ::  tm_pfx(5)


contains


subroutine Atm_WriteNCF_Initial

    integer     aErr, s, itm, i
    character(40)   tmDir(5)
    character(120)  fp_pfx
    character(31)   tms

    if(AtmOutNCF6hourly .or. AtmOutNCFDaily .or. AtmOutNCFMonthly .or. AtmOutNCFYearly) then

    allocate(var3d(Imin:Imax,Jmin:Jmax,1:Atm_KMax,1:NumSrc,NumForm(Atm),1:2),&
            &var2ds(Imin:Imax,Jmin:Jmax,1:NumSurf,1:NumSrc,NumForm(Atm),3:5),&
            &var2dlc(Imin:Imax,Jmin:Jmax,1:NumSurf),&
            &src_ind(NumSrc), srf_ind(NumSurf), stat=aErr)
    if(aErr/=0) stop 'STOP: Memory allocation error in Atm_Output_Initial'

#if RTYPE==2
    allocate(sta_nm(1:NAir), stp_nm(1:NPrec), sta_mod(1:NAir,1:NumSrc,MaxForm), &
            &stp_mod(1:NPrec,1:NumSrc,MaxForm), stpr_mod(1:NPrec), stat=aErr)
    if(aErr/=0) stop 'STOP: Memory allocation error in Atm_WriteNCF_Initial'
#endif

! Set up matrix for output variables control
    fld_ncf_out(:,AC_3D) = (/.false.,.false., .true., .true., .true./)
    fld_ncf_out(:,MR_3D) = (/.false., .true.,.false., .true., .true./)
    fld_ncf_out(:,DD_2D) = (/.false.,.false., .true., .true., .true./)
    fld_ncf_out(:,WD_2D) = (/.false.,.false., .true., .true., .true./)
    fld_ncf_out(:,RE_2D) = (/.false.,.false., .true., .true., .true./)
    fld_ncf_out(:,AC_ST) = (/.false.,.false., .true., .true., .true./)
    fld_ncf_out(:,WD_ST) = (/.false.,.false., .true., .true., .true./)

! Fill in parts of output file names and headers
    tmDir(1) = trim(HourlyDir)
    tmDir(2) = trim(SixHourlyDir)
    tmDir(3) = trim(DailyDir)
    tmDir(4) = trim(MonthlyDir)
    tmDir(5) = trim(YearlyDir)
    tm_pfx(:)= (/'Hourly    ','Six hourly','Daily     ','Monthly   ','Annual    '/)

    do s = 1, NumSubsMedia(Atm)
        do itm=1, 5
            fp_pfx = trim(OutPath)//trim(AtmDir)//trim(FieldsNCFDir)//trim(tmDir(itm))//trim(SubsID(s))
            ncfPath(itm,AC_3D) = trim(fp_pfx)//trim(ConcName)
            ncfPath(itm,MR_3D) = trim(fp_pfx)//trim(MixRatName)
            ncfPath(itm,DD_2D) = trim(fp_pfx)//trim(DryName)
            ncfPath(itm,WD_2D) = trim(fp_pfx)//trim(WetName)
            ncfPath(itm,RE_2D) = trim(fp_pfx)//trim(RemName)
            ncfPath(itm,AC_ST) = trim(fp_pfx)//trim(ConcMonitName)
            ncfPath(itm,WD_ST) = trim(fp_pfx)//trim(WetMonitName)

            do i = 1, NCF_VARS_NO
                v_hdr(itm,i)= trim(tm_pfx(itm))//' '//trim(fld_nm(i))//' '//trim(SubsID(s))
            end do
        end do ! itm
    end do ! s

    do s = 1, NumSrc
        src_ind(s) = s
    end do
    do s = 1, NumSurf
        srf_ind(s) = s
    end do
    do s=1, Atm_Kmax
        levs(s) = Sigma(s)
    enddo

! Write supporting information for processing NetCDF model output
    write(tms, '(a12,i4,a15)') "years since ",BegDate(yr),"-01-01 00:00:00"
    call Atm_WriteMeshAreaNCF(trim(OutPath)//trim(AtmDir)//trim(FieldsNCFDir)//"MeshArea.nc","Mesh area",tms)   
    call Atm_WriteLandCoverNCF(trim(OutPath)//trim(AtmDir)//trim(FieldsNCFDir)//"LandCover.nc","Land Cover",tms) 

    end if ! If ncf output requested

end subroutine Atm_WriteNCF_Initial





subroutine Atm_WriteNCF_Final
    integer     dErr
    
    if(AtmOutNCF6hourly .or. AtmOutNCFDaily .or. AtmOutNCFMonthly .or. AtmOutNCFYearly) then
        deallocate(var3d, var2ds, var2dlc, src_ind, srf_ind, stat=dErr)
        if(dErr/=0) stop 'STOP: Memory deallocation error in Atm_WriteNCF_Final'

#if RTYPE==2
        deallocate(sta_nm, stp_nm, sta_mod, stp_mod, stpr_mod, stat=dErr)
        if(dErr/=0) stop 'STOP: Memory deallocation error in Atm_WriteNCF_Final'
#endif
    end if

end subroutine Atm_WriteNCF_Final


!*******************************************************************************************************
! Subroutine writing model data with contributions of emission sources and division by surface type
!*******************************************************************************************************
subroutine Atm_WriteFieldsNCF(prnType, nh)
    character(*), intent(in)    ::  prnType
    integer,intent(in),optional ::  nh                          ! current hour
    integer                     ::  s,tl, n, f, i, j, src
    character(31)               ::  tms
    character(40)               ::  fp
    real                        ::  tp

    write(YearNum,'(i4)') Year

    do s = 1, NumSubsMedia(Atm)

        select case(prnType)
            case('Hourly')
                tl = TM_HOURLY

            case('6hourly')
                tl = TM_6HOURLY
                tp = SecInDay
                write(MonthNum,'(i2.2)') Month
                write(DayNum,'(i2.2)') Day
                write(PerNum,'(i2.2)') Period
                fp = YearNum//MonthNum//DayNum//'_'//PerNum//'.nc'
                write(tms, '(a12,i4,a1,i2.2,a1,i2.2,a1,i2.2,a6)') "hours since ",Year,"-",Month,"-",Day," ",nh,":00:00"

!                var2d(Imin:Imax,Jmin:Jmax,1) = PxDay(Imin:Imax,Jmin:Jmax)               ! Pa
                var2d(Imin:Imax,Jmin:Jmax,1) = Px(Imin:Imax,Jmin:Jmax,Period,1)+Ptop     ! Surface pressure [Pa]
                var2d(Imin:Imax,Jmin:Jmax,2) = PrecipDay(Imin:Imax,Jmin:Jmax)            ! m

                do n = 1, sFormNum(Atm, s)
                    f = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    if(fld_ncf_out(tl,AC_3D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
#if RTYPE==1
                                var3d(i,j,1:Atm_KMax,1,n,AC_3D) = Atm_Conc(i,j,1:Atm_KMax,f)
#endif
#if RTYPE==2
                                do src = 1, NumSrc
                                    var3d(i,j,1:Atm_KMax,src,n,AC_3D) = Atm_Conc(i,j,1:Atm_KMax,f) * &
                                                    & Atm_Contrib(i,j,1:Atm_KMax,f,src)                     
                                end do ! src
#endif
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,MR_3D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
#if RTYPE==1
                                var3d(i,j,1:Atm_KMax,1,n,MR_3D) = Atm_MixRatio(i,j,1:Atm_KMax,f)     
#endif
#if RTYPE==2
                                do src = 1, NumSrc
                                    var3d(i,j,1:Atm_KMax,src,n,MR_3D) = Atm_MixRatio(i,j,1:Atm_KMax,f) * &
                                                    & Atm_Contrib(i,j,1:Atm_KMax,f,src)                     
                                end do ! src
#endif
                            end do ! i
                        end do ! j
                    end if
                end do ! n

            case('Daily')
                tl = TM_DAILY
                tp = SecInDay
                write(MonthNum,'(i2.2)') Month
                write(DayNum,'(i2.2)') Day
                fp = YearNum//MonthNum//DayNum//'.nc'
                write(tms, '(a12,i4,a1,i2.2,a1,i2.2,a9)') "days since  ",Year,"-",Month,"-",Day," 00:00:00"

                var2d(Imin:Imax,Jmin:Jmax,1) = PxDay(Imin:Imax,Jmin:Jmax)/tp + Ptop         ! Pa
                var2d(Imin:Imax,Jmin:Jmax,2) = PrecipDay(Imin:Imax,Jmin:Jmax)         ! m

                do n = 1, sFormNum(Atm, s)
                    f = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    if(fld_ncf_out(tl,AC_3D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var3d(i,j,1:Atm_KMax,1:NumSrc,n,AC_3D) = Atm_ConcDay(i,j,1:Atm_KMax,f,1:NumSrc)/tp
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,MR_3D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var3d(i,j,1:Atm_KMax,1:NumSrc,n,MR_3D) = Atm_MixDay(i,j,1:Atm_KMax,f,1:NumSrc)/tp
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,DD_2D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2ds(i,j,1:NumSurf,1:NumSrc,n,DD_2D) = DryDepDay(i,j,f,1:NumSurf,1:NumSrc)
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,WD_2D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2ds(i,j,1:NumSurf,1:NumSrc,n,WD_2D) = WetDepDay(i,j,f,1:NumSurf,1:NumSrc)
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,RE_2D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2ds(i,j,1:NumSurf,1:NumSrc,n,WD_2D) = DryRemDay(i,j,f,1:NumSurf,1:NumSrc)
                            end do ! i
                        end do ! j
                    end if
                end do ! n

            case('Monthly')
                tl = TM_MONTHLY
                tp = min(SecInDay*real(MonthDays(Month)), timeCalc)
                write(MonthNum,'(i2.2)') Month
                fp = YearNum//MonthNum//'.nc'
                write(tms, '(a12,i4,a1,i2.2,a12)') "months since",Year,"-",Month,"-01 00:00:00"

                var2d(Imin:Imax,Jmin:Jmax,1) = PxMonth(Imin:Imax,Jmin:Jmax)/tp + Ptop         ! Pa
                var2d(Imin:Imax,Jmin:Jmax,2) = PrecipMonth(Imin:Imax,Jmin:Jmax)         ! m

                do n = 1, sFormNum(Atm, s)
                    f = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    if(fld_ncf_out(tl,AC_3D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var3d(i,j,1:Atm_KMax,1:NumSrc,n,AC_3D) = Atm_ConcMonth(i,j,1:Atm_KMax,f,1:NumSrc)/tp
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,MR_3D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var3d(i,j,1:Atm_KMax,1:NumSrc,n,MR_3D) = Atm_MixMonth(i,j,1:Atm_KMax,f,1:NumSrc)/tp
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,DD_2D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2ds(i,j,1:NumSurf,1:NumSrc,n,DD_2D) = DryDepMonth(i,j,f,1:NumSurf,1:NumSrc)
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,WD_2D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2ds(i,j,1:NumSurf,1:NumSrc,n,WD_2D) = WetDepMonth(i,j,f,1:NumSurf,1:NumSrc)
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,RE_2D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2ds(i,j,1:NumSurf,1:NumSrc,n,RE_2D) = DryRemMonth(i,j,f,1:NumSurf,1:NumSrc)
                            end do ! i
                        end do ! j
                    end if
                end do ! n

            case('Yearly')
                tl = TM_YEARLY
                tp = timeCalc
                fp = YearNum//'.nc'
                write(tms, '(a12,i4,a15)') "years since ",Year,"-01-01 00:00:00"

                var2d(Imin:Imax,Jmin:Jmax,1) = PxYear(Imin:Imax,Jmin:Jmax)/tp + Ptop         ! Pa
                var2d(Imin:Imax,Jmin:Jmax,2) = PrecipYear(Imin:Imax,Jmin:Jmax)         ! m

                do n = 1, sFormNum(Atm, s)
                    f = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                    if(fld_ncf_out(tl,AC_3D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var3d(i,j,1:Atm_KMax,1:NumSrc,n,AC_3D) = Atm_ConcYear(i,j,1:Atm_KMax,f,1:NumSrc)/tp
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,MR_3D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var3d(i,j,1:Atm_KMax,1:NumSrc,n,MR_3D) = Atm_MixYear(i,j,1:Atm_KMax,f,1:NumSrc)/tp
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,DD_2D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2ds(i,j,1:NumSurf,1:NumSrc,n,DD_2D) = DryDepYear(i,j,f,1:NumSurf,1:NumSrc)
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,WD_2D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2ds(i,j,1:NumSurf,1:NumSrc,n,WD_2D) = WetDepYear(i,j,f,1:NumSurf,1:NumSrc)
                            end do ! i
                        end do ! j
                    end if
                    if(fld_ncf_out(tl,RE_2D)) then
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2ds(i,j,1:NumSurf,1:NumSrc,n,RE_2D) = DryRemYear(i,j,f,1:NumSurf,1:NumSrc)
                            end do ! i
                        end do ! j
                    end if
                end do ! n
        end select

        if(fld_ncf_out(tl,AC_3D)) call Atm_Write3dField(trim(ncfPath(tl,AC_3D))//'_'//trim(fp),v_hdr(tl,AC_3D),Atm,s,tms,AC_3D)
        if(fld_ncf_out(tl,MR_3D)) call Atm_Write3dField(trim(ncfPath(tl,MR_3D))//'_'//trim(fp),v_hdr(tl,MR_3D),Atm,s,tms,MR_3D)
        if(fld_ncf_out(tl,DD_2D)) call Atm_Write2dField(trim(ncfPath(tl,DD_2D))//'_'//trim(fp),v_hdr(tl,DD_2D),Atm,s,tms,DD_2D,0)
        if(fld_ncf_out(tl,WD_2D)) call Atm_Write2dField(trim(ncfPath(tl,WD_2D))//'_'//trim(fp),v_hdr(tl,WD_2D),Atm,s,tms,WD_2D,0)
        if(fld_ncf_out(tl,RE_2D)) call Atm_Write2dField(trim(ncfPath(tl,RE_2D))//'_'//trim(fp),v_hdr(tl,RE_2D),Atm,s,tms,RE_2D,0)

    end do ! s

end subroutine Atm_WriteFieldsNCF


!*******************************************************************************************************
! Writing 3d model dataset to netcdf file
!*******************************************************************************************************
subroutine Atm_Write3dField(fpath, v_hdr, med, sbs, tm_str, vr_id)
    character(*), intent(in)    :: fpath
    character(*), intent(in)    :: v_hdr
    integer, intent(in)         :: med
    integer, intent(in)         :: sbs
    character(*), intent(in)    :: tm_str
    integer, intent(in)         :: vr_id

    integer     nn, ncid
    integer     timeid, lonid, latid, levid, srcid, src_nm_id, psid, ptopid
    integer     time_dimid, lon_dimid, lat_dimid, lev_dimid, src_dimid, src_len_dimid, src_num_dimid
    integer     dimid3d(5), start3d(5), count3d(5)
    integer     dimid2d(3), start2d(3), count2d(3)
    integer     dimid1d(2), start1d(2), count1d(2)
    integer     dimidSrc(2)
    integer     k, src, f, varid(MaxForm), n, st                    

    character(40)   vr_nm

    nn = 0
! Open file
    st = nf90_create(trim(fpath), nf90_clobber, ncid)                               
    if(st /= nf90_noerr) then                                                                 
        print '(/,"STOP in Atm_Write3dField: Cannot open file ''",a,"''",/)', trim(fpath)   
        stop                                                                            
    endif                                                                               
! Define dimensions
    call check_NC(nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid), nn+2)
    call check_NC(nf90_def_dim(ncid, LAT_NAME, Jmax-Jmin+1, lat_dimid), nn+3)
    call check_NC(nf90_def_dim(ncid, LON_NAME, Imax-Imin+1, lon_dimid), nn+4)
    call check_NC(nf90_def_dim(ncid, LEV_NAME, Atm_KMax, lev_dimid), nn+5)
    call check_NC(nf90_def_dim(ncid, SRC_NAME, NumSrc, src_dimid), nn+6)
    call check_NC(nf90_def_dim(ncid, 'SrcStrLen', len(SourcID), src_len_dimid), nn+7)
    call check_NC(nf90_def_dim(ncid, 'NumSrcStr', NumSrc, src_num_dimid), nn+8)

! Define and put attributes
    call check_NC(nf90_def_var(ncid, TIME_NAME, NF90_REAL, time_dimid, timeid), nn+9)
    call check_NC(nf90_put_att(ncid, timeid, UNITS, tm_str), nn+10)
    call check_NC(nf90_put_att(ncid, timeid, CALENDAR, TIME_CLD), nn+11)
    call check_NC(nf90_put_att(ncid, timeid, AXIS, TIME_AXS), nn+12)
    call check_NC(nf90_put_att(ncid, timeid, LNG_NAME, TIME_LNG), nn+13)
    call check_NC(nf90_put_att(ncid, timeid, STD_NAME, TIME_STD), nn+14)

    call check_NC(nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, latid), nn+15)
    call check_NC(nf90_put_att(ncid, latid, UNITS, LAT_UNITS), nn+16)
    call check_NC(nf90_put_att(ncid, latid, AXIS, LAT_AXS), nn+17)
    call check_NC(nf90_put_att(ncid, latid, LNG_NAME, LAT_LNG), nn+18)
    call check_NC(nf90_put_att(ncid, latid, STD_NAME, LAT_STD), nn+19)

    call check_NC(nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lonid), nn+20)
    call check_NC(nf90_put_att(ncid, lonid, UNITS, LON_UNITS), nn+21)
    call check_NC(nf90_put_att(ncid, lonid, AXIS, LON_AXS), nn+22)
    call check_NC(nf90_put_att(ncid, lonid, LNG_NAME, LON_LNG), nn+23)
    call check_NC(nf90_put_att(ncid, lonid, STD_NAME, LON_STD), nn+24)

    call check_NC(nf90_def_var(ncid, LEV_NAME, NF90_REAL, lev_dimid, levid), nn+25)
    call check_NC(nf90_put_att(ncid, levid, UNITS, LEV_UNITS), nn+26)
    call check_NC(nf90_put_att(ncid, levid, "positive", "up"), nn+27)
    call check_NC(nf90_put_att(ncid, levid, AXIS, LEV_AXS), nn+28)
    call check_NC(nf90_put_att(ncid, levid, STD_NAME, LEV_STD), nn+29)
    call check_NC(nf90_put_att(ncid, levid, LNG_NAME, LEV_LNG), nn+30)
    call check_NC(nf90_put_att(ncid, levid, TERMS, LEV_TERM), nn+31)

    call check_NC(nf90_def_var(ncid, PTOP_NAME, NF90_REAL, ptopid), nn+32)
    call check_NC(nf90_put_att(ncid, ptopid, UNITS, PTOP_UNITS), nn+33)
    call check_NC(nf90_put_att(ncid, ptopid, STD_NAME, PTOP_STD), nn+34)
    call check_NC(nf90_put_att(ncid, ptopid, LNG_NAME, PTOP_LNG), nn+35)

    dimid2d = (/lon_dimid,lat_dimid,time_dimid/)
    call check_NC(nf90_def_var(ncid, PSUR_NAME, NF90_REAL, dimid2d, psid), nn+36)
    call check_NC(nf90_put_att(ncid, psid, UNITS, PSUR_UNITS), nn+37)
    call check_NC(nf90_put_att(ncid, psid, STD_NAME, PSUR_STD), nn+38)
    call check_NC(nf90_put_att(ncid, psid, LNG_NAME, PSUR_LNG), nn+39)

    dimidSrc = (/src_len_dimid,src_num_dimid/)
    call check_NC(nf90_def_var(ncid, SRC_ARR, NF90_CHAR, dimidSrc, src_nm_id), nn+40)
    call check_NC(nf90_def_var(ncid, SRC_NAME, NF90_INT, src_dimid, srcid), nn+41)
    call check_NC(nf90_put_att(ncid, srcid, UNITS, SRC_UNITS), nn+42)
    call check_NC(nf90_put_att(ncid, srcid, LNG_NAME, SRC_LNG), nn+43)
    call check_NC(nf90_put_att(ncid, srcid, STD_NAME, SRC_STD), nn+44)

    dimid3d = (/lon_dimid,lat_dimid,lev_dimid,src_dimid,time_dimid/)
    do n = 1, sFormNum(Atm, sbs)
        f = sFormInd(Atm, sbs, n)
        vr_nm = trim(v_nm(vr_id))//'_'//trim(SubsID(sbs))//trim(FormID(Atm,sbs,f))
        call check_NC(nf90_def_var(ncid, vr_nm, NF90_REAL, dimid3d, varid(f)), nn+45)
        call check_NC(nf90_put_att(ncid, varid(f), UNITS, trim(fld_unt(vr_id))), nn+46)
        call check_NC(nf90_put_att(ncid, varid(f), STD_NAME, vr_nm), nn+47)
        call check_NC(nf90_put_att(ncid, varid(f), LNG_NAME, v_hdr), nn+48)
    end do

    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "title", "GLEMOS model output"), nn+49)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Institution", "Meteorological Synthesizing Centre - East"), nn+50)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "institute_id", "MSC-E"), nn+51)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Contact", "MSC-E (msce@msceast.org)"), nn+52)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "project_id", "EMEP"), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model name", MODEL_NAME), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model long name", MODEL_LONG_NAME), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model version", MODEL_VERSION), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-2.0"), nn+54)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Grid code", GridCode), nn+55)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Substance", SubsID(sbs)), nn+56)

    call check_NC(nf90_enddef(ncid), nn+57)

! Write coordinates
    call check_NC(nf90_put_var(ncid, timeid, 0), nn+58)
    
    call check_NC(nf90_put_var(ncid, latid, LatMesh(Jmin:Jmax)/pi180), nn+59)
    call check_NC(nf90_put_var(ncid, lonid, LongMeshR(Imin:Imax)/pi180), nn+60)
    call check_NC(nf90_put_var(ncid, levid, levs), nn+61)
    call check_NC(nf90_put_var(ncid, srcid, src_ind), nn+62)
    call check_NC(nf90_put_var(ncid, src_nm_id, SourcID(1:NumSrc)), nn+63)
    call check_NC(nf90_put_var(ncid, ptopid, Ptop), nn+64)

! Prepare array for writing
    do n = 1, sFormNum(Atm, sbs)
        f = sFormInd(Atm, sbs, n)
        do k = 1, Atm_KMax
            do src = 1, NumSrc
                call GridDisAggreg(var3d(Imin:Imax,Jmin:Jmax,k,src,n,vr_id),1)
            end do ! src
        end do ! k
    end do ! n

! Write array to file
    start3d = (/1,1,1,1,1/)
    count3d = (/Imax-Imin+1,Jmax-Jmin+1,Atm_KMax,NumSrc,1/)
    do n = 1, sFormNum(Atm, sbs)
        f = sFormInd(Atm, sbs, n)
        call check_NC(nf90_put_var(ncid,varid(f),var3d(Imin:Imax,Jmin:Jmax,1:Atm_KMax,1:NumSrc,n,vr_id),start3d,count3d),nn+65)
    end do ! n

! Write supporting variables - Surface pressure
    call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,1),1)
    start2d = (/1,1,1/)
    count2d = (/Imax-Imin+1,Jmax-Jmin+1,1/)
    call check_NC(nf90_put_var(ncid, psid, var2d(Imin:Imax,Jmin:Jmax,1), start2d, count2d), nn+66)

! Close file
    call check_NC(nf90_close(ncid), nn+67)

end subroutine Atm_Write3dField


!*******************************************************************************************************
! Writing 2d model dataset to netcdf file
!*******************************************************************************************************
subroutine Atm_Write2dField(fpath, v_hdr, med, sbs, tm_str, vr_id, pr_id)
    character(*), intent(in)    :: fpath
    character(*), intent(in)    :: v_hdr
    integer, intent(in)         :: med
    integer, intent(in)         :: sbs
    character(*), intent(in)    :: tm_str
    integer, intent(in)         :: vr_id
    integer, intent(in)         :: pr_id

    integer     nn, ncid
    integer     timeid, lonid, latid, srcid, src_nm_id, srfid, srf_nm_id, prid
    integer     time_dimid, lon_dimid, lat_dimid, src_dimid, src_len_dimid, src_num_dimid, srf_dimid, srf_len_dimid, srf_num_dimid
    integer     dimid2ds(5), start2ds(5), count2ds(5)
    integer     dimid2d(3), start2d(3), count2d(3)
    integer     dimid1d(2), start1d(2), count1d(2)
    integer     dimidSrc(2), dimidSrf(2)
    integer     srf, src, varid(MaxForm), n, f, st              

    character(40)   vr_nm
    
    nn = 100
! Open file
    st= nf90_create(trim(fpath), nf90_clobber, ncid)                                
    if(st /= nf90_noerr) then                                                             
        print '(/,"STOP in Atm_Write2dField: Cannot open file ''",a,"''",/)', trim(fpath)   
        stop                                                                            
    endif                                                                               
! Define dimensions
    call check_NC(nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid), nn+2)
    call check_NC(nf90_def_dim(ncid, LAT_NAME, Jmax-Jmin+1, lat_dimid), nn+3)
    call check_NC(nf90_def_dim(ncid, LON_NAME, Imax-Imin+1, lon_dimid), nn+4)
    call check_NC(nf90_def_dim(ncid, SRC_NAME, NumSrc, src_dimid), nn+5)
    call check_NC(nf90_def_dim(ncid, 'SrcStrLen', len(SourcID), src_len_dimid), nn+6)
    call check_NC(nf90_def_dim(ncid, 'NumSrcStr', NumSrc, src_num_dimid), nn+7)
    call check_NC(nf90_def_dim(ncid, SRF_NAME, NumSurf, srf_dimid), nn+8)
    call check_NC(nf90_def_dim(ncid, 'SrfStrLen', len(SurfType), srf_len_dimid), nn+9)
    call check_NC(nf90_def_dim(ncid, 'NumSrfStr', NumSurf, srf_num_dimid), nn+10)

! Define and put attributes
    call check_NC(nf90_def_var(ncid, TIME_NAME, NF90_REAL, time_dimid, timeid), nn+11)
    call check_NC(nf90_put_att(ncid, timeid, UNITS, tm_str), nn+12)
    call check_NC(nf90_put_att(ncid, timeid, CALENDAR, TIME_CLD), nn+13)
    call check_NC(nf90_put_att(ncid, timeid, AXIS, TIME_AXS), nn+14)
    call check_NC(nf90_put_att(ncid, timeid, LNG_NAME, TIME_LNG), nn+15)
    call check_NC(nf90_put_att(ncid, timeid, STD_NAME, TIME_STD), nn+16)

    call check_NC(nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, latid), nn+17)
    call check_NC(nf90_put_att(ncid, latid, UNITS, LAT_UNITS), nn+18)
    call check_NC(nf90_put_att(ncid, latid, AXIS, LAT_AXS), nn+19)
    call check_NC(nf90_put_att(ncid, latid, LNG_NAME, LAT_LNG), nn+20)
    call check_NC(nf90_put_att(ncid, latid, STD_NAME, LAT_STD), nn+21)

    call check_NC(nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lonid), nn+22)
    call check_NC(nf90_put_att(ncid, lonid, UNITS, LON_UNITS), nn+23)
    call check_NC(nf90_put_att(ncid, lonid, AXIS, LON_AXS), nn+24)
    call check_NC(nf90_put_att(ncid, lonid, LNG_NAME, LON_LNG), nn+25)
    call check_NC(nf90_put_att(ncid, lonid, STD_NAME, LON_STD), nn+26)

    if(pr_id == 1) then
        dimid2d = (/lon_dimid,lat_dimid,time_dimid/)
        call check_NC( nf90_def_var(ncid, PREC_NAME, NF90_REAL, dimid2d, prid), nn+27)
        call check_NC( nf90_put_att(ncid, prid, UNITS, PREC_UNITS), nn+28)
        call check_NC( nf90_put_att(ncid, prid, STD_NAME, PREC_STD), nn+29)
        call check_NC( nf90_put_att(ncid, prid, LNG_NAME, PREC_LNG), nn+30)
    end if

    dimidSrc = (/src_len_dimid,src_num_dimid/)
    call check_NC(nf90_def_var(ncid, SRC_ARR, NF90_CHAR, dimidSrc, src_nm_id), nn+31)
    call check_NC(nf90_def_var(ncid, SRC_NAME, NF90_INT, src_dimid, srcid), nn+32)
    call check_NC(nf90_put_att(ncid, srcid, UNITS, SRC_UNITS), nn+33)
    call check_NC(nf90_put_att(ncid, srcid, LNG_NAME, SRC_LNG), nn+34)
    call check_NC(nf90_put_att(ncid, srcid, STD_NAME, SRC_STD), nn+35)

    dimidSrf = (/srf_len_dimid,srf_num_dimid/)
    call check_NC(nf90_def_var(ncid, SRF_ARR, NF90_CHAR, dimidSrf, srf_nm_id), nn+37)
    call check_NC(nf90_def_var(ncid, SRF_NAME, NF90_INT, srf_dimid, srfid), nn+38)
    call check_NC(nf90_put_att(ncid, srfid, UNITS, SRF_UNITS), nn+39)
    call check_NC(nf90_put_att(ncid, srfid, LNG_NAME, SRF_LNG), nn+40)
    call check_NC(nf90_put_att(ncid, srfid, STD_NAME, SRF_STD), nn+41)

    dimid2ds = (/lon_dimid,lat_dimid,srf_dimid,src_dimid,time_dimid/)
    do n = 1, sFormNum(Atm, sbs)
        f = sFormInd(Atm, sbs, n)
        vr_nm = trim(v_nm(vr_id))//'_'//trim(SubsID(sbs))//trim(FormID(Atm,sbs,f))
        call check_NC(nf90_def_var(ncid, vr_nm, NF90_REAL, dimid2ds, varid(f)), nn+42)
        call check_NC(nf90_put_att(ncid, varid(f), UNITS, trim(fld_unt(vr_id))), nn+43)
        call check_NC(nf90_put_att(ncid, varid(f), STD_NAME, vr_nm), nn+44)
        call check_NC(nf90_put_att(ncid, varid(f), LNG_NAME, v_hdr), nn+45)
    end do

    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "title", "GLEMOS model output"), nn+46)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Institution", "Meteorological Synthesizing Centre - East"), nn+47)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "institute_id", "MSC-E"), nn+48)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Contact", "MSC-E (msce@msceast.org)"), nn+49)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "project_id", "EMEP"), nn+50)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model name", MODEL_NAME), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model long name", MODEL_LONG_NAME), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model version", MODEL_VERSION), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-2.0"), nn+51)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Grid code", GridCode), nn+52)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Substance", SubsID(sbs)), nn+53)

    call check_NC(nf90_enddef(ncid), nn+54)

! Write coordinates
    call check_NC(nf90_put_var(ncid, timeid, 0), nn+55)
    call check_NC(nf90_put_var(ncid, latid, LatMesh(Jmin:Jmax)/pi180), nn+56)
    call check_NC(nf90_put_var(ncid, lonid, LongMeshR(Imin:Imax)/pi180), nn+57)
    call check_NC(nf90_put_var(ncid, srcid, src_ind), nn+58)
    call check_NC(nf90_put_var(ncid, src_nm_id, SourcID(1:NumSrc)), nn+59)
    call check_NC(nf90_put_var(ncid, srfid, srf_ind), nn+60)
    call check_NC(nf90_put_var(ncid, srf_nm_id, SurfType(1:NumSurf)), nn+61)

! Prepae array
    do n = 1, sFormNum(Atm, sbs)
        f = sFormInd(Atm, sbs, n)
        do srf = 1, NumSurf
            do src = 1, NumSrc
                call GridDisAggreg(var2ds(Imin:Imax,Jmin:Jmax,srf,src,n,vr_id),2)       ! Corrected disagr 0 - AG 08.02.17
            end do ! src
        end do ! srf
    end do ! n

! Write array
    start2ds = (/1,1,1,1,1/)
    count2ds = (/Imax-Imin+1,Jmax-Jmin+1,NumSurf,NumSrc,1/)
    do n = 1, sFormNum(Atm, sbs)
        f = sFormInd(Atm, sbs, n)
        call check_NC(nf90_put_var(ncid,varid(f),var2ds(Imin:Imax,Jmin:Jmax,1:NumSurf,1:NumSrc,n,vr_id),start2ds,count2ds),nn+62)
    end do

! Write supporting variables
    if(pr_id == 1) then
        call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,2),1)
        start2d = (/1,1,1/)
        count2d = (/Imax-Imin+1,Jmax-Jmin+1,1/)
        call check_NC(nf90_put_var(ncid, prid, var2d(Imin:Imax,Jmin:Jmax,2), start2d, count2d), nn+63)
    end if

! Close file
    call check_NC(nf90_close(ncid), nn+64)

end subroutine Atm_Write2dField

!*******************************************************************************************************
! Writing mesh area to netcdf file
!*******************************************************************************************************
subroutine Atm_WriteMeshAreaNCF(fpath,v_hdr,tm_str)
    character(*), intent(in)    :: fpath
    character(*), intent(in)    :: v_hdr
    character(*), intent(in)    :: tm_str

    integer     nn, ncid, st                                
    integer     timeid, lonid, latid, varid
    integer     i, j
    integer     time_dimid, lon_dimid, lat_dimid
    integer     dimid2d(3), start2d(3), count2d(3)

    nn = 200
! Open file
    st = nf90_create(trim(fpath), nf90_clobber, ncid)                                       
    if(st/= nf90_noerr) then                                                             
        print '(/,"STOP in Atm_WriteMeshAreaNCF: Cannot open file ''",a,"''",/)', trim(fpath)   
        stop                                                                            
    endif                                                                               
! Define dimensions
    call check_NC(nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid), nn+2)
    call check_NC(nf90_def_dim(ncid, LAT_NAME, Jmax-Jmin+1, lat_dimid), nn+3)
    call check_NC(nf90_def_dim(ncid, LON_NAME, Imax-Imin+1, lon_dimid), nn+4)

! Define and put attributes
    call check_NC(nf90_def_var(ncid, TIME_NAME, NF90_REAL, time_dimid, timeid), nn+5)
    call check_NC(nf90_put_att(ncid, timeid, UNITS, tm_str), nn+6)
    call check_NC(nf90_put_att(ncid, timeid, CALENDAR, TIME_CLD), nn+7)
    call check_NC(nf90_put_att(ncid, timeid, AXIS, TIME_AXS), nn+8)
    call check_NC(nf90_put_att(ncid, timeid, LNG_NAME, TIME_LNG), nn+9)
    call check_NC(nf90_put_att(ncid, timeid, STD_NAME, TIME_STD), nn+10)

    call check_NC(nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, latid), nn+11)
    call check_NC(nf90_put_att(ncid, latid, UNITS, LAT_UNITS), nn+12)
    call check_NC(nf90_put_att(ncid, latid, AXIS, LAT_AXS), nn+13)
    call check_NC(nf90_put_att(ncid, latid, LNG_NAME, LAT_LNG), nn+14)
    call check_NC(nf90_put_att(ncid, latid, STD_NAME, LAT_STD), nn+15)

    call check_NC(nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lonid), nn+16)
    call check_NC(nf90_put_att(ncid, lonid, UNITS, LON_UNITS), nn+17)
    call check_NC(nf90_put_att(ncid, lonid, AXIS, LON_AXS), nn+18)
    call check_NC(nf90_put_att(ncid, lonid, LNG_NAME, LON_LNG), nn+19)
    call check_NC(nf90_put_att(ncid, lonid, STD_NAME, LON_STD), nn+20)

    dimid2d = (/lon_dimid,lat_dimid,time_dimid/)
    call check_NC(nf90_def_var(ncid, "mesh_area", NF90_REAL, dimid2d, varid), nn+21)
    call check_NC(nf90_put_att(ncid, varid, UNITS, "m2"), nn+22)
    call check_NC(nf90_put_att(ncid, varid, STD_NAME, "mesh_area"), nn+23)
    call check_NC(nf90_put_att(ncid, varid, LNG_NAME, v_hdr), nn+24)

    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "title", "GLEMOS model output"), nn+25)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Institution", "Meteorological Synthesizing Centre - East"), nn+26)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "institute_id", "MSC-E"), nn+27)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Contact", "MSC-E (msce@msceast.org)"), nn+28)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "project_id", "EMEP"), nn+29)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model name", MODEL_NAME), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model long name", MODEL_LONG_NAME), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model version", MODEL_VERSION), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-2.0"), nn+30)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Grid code", GridCode), nn+31)

    call check_NC(nf90_enddef(ncid), nn+32)

! Write coordinates
    call check_NC(nf90_put_var(ncid, timeid, 0), nn+33)
    call check_NC(nf90_put_var(ncid, latid, LatMesh(Jmin:Jmax)/pi180), nn+34)
    call check_NC(nf90_put_var(ncid, lonid, LongMeshR(Imin:Imax)/pi180), nn+35)

! Disaggregate array - (added 08.02.17 - AG)
    do j = Jmin, Jmax
        do i = minI(j), maxI(j)
            var2d(i,j,1) = MeshArea(i,j)
        end do
    end do
    call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,1),2)

! Write array
    start2d = (/1,1,1/)
    count2d = (/Imax-Imin+1,Jmax-Jmin+1,1/)
    call check_NC(nf90_put_var(ncid, varid, var2d(Imin:Imax,Jmin:Jmax,1), start2d, count2d), nn+36)

! Close file
    call check_NC(nf90_close(ncid), nn+37)

end subroutine Atm_WriteMeshAreaNCF

!*******************************************************************************************************
! Writing land cover array to netcdf file
!*******************************************************************************************************
subroutine Atm_WriteLandCoverNCF(fpath,v_hdr,tm_str)
    character(*), intent(in)    :: fpath
    character(*), intent(in)    :: v_hdr
    character(*), intent(in)    :: tm_str

    integer     nn, ncid, st, ns, i, j
    integer     timeid, lonid, latid, varid
    integer     time_dimid, lon_dimid, lat_dimid
    integer     srf_dimid, srf_len_dimid, srf_num_dimid, dimidSrf(2), srfid, srf_nm_id
    integer     dimid2d(4), start2d(4), count2d(4)

    nn = 500

! Open file
    st = nf90_create(trim(fpath), nf90_clobber, ncid)                                       
    if(st /= nf90_noerr) then                                                             
        print '(/,"STOP in Atm_WriteLandCoverNCF: Cannot open file ''",a,"''",/)', trim(fpath)   
        stop                                                                            
    endif                                                                               
! Define dimensions
    call check_NC(nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid), nn+2)
    call check_NC(nf90_def_dim(ncid, LAT_NAME, Jmax-Jmin+1, lat_dimid), nn+3)
    call check_NC(nf90_def_dim(ncid, LON_NAME, Imax-Imin+1, lon_dimid), nn+4)
    call check_NC(nf90_def_dim(ncid, SRF_NAME, NumSurf, srf_dimid), nn+8)
    call check_NC(nf90_def_dim(ncid, 'SrfStrLen', len(SurfType), srf_len_dimid), nn+9)
    call check_NC(nf90_def_dim(ncid, 'NumSrfStr', NumSurf, srf_num_dimid), nn+10)

! Define and put attributes
    call check_NC(nf90_def_var(ncid, TIME_NAME, NF90_REAL, time_dimid, timeid), nn+5)
    call check_NC(nf90_put_att(ncid, timeid, UNITS, tm_str), nn+6)
    call check_NC(nf90_put_att(ncid, timeid, CALENDAR, TIME_CLD), nn+7)
    call check_NC(nf90_put_att(ncid, timeid, AXIS, TIME_AXS), nn+8)
    call check_NC(nf90_put_att(ncid, timeid, LNG_NAME, TIME_LNG), nn+9)
    call check_NC(nf90_put_att(ncid, timeid, STD_NAME, TIME_STD), nn+10)

    call check_NC(nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, latid), nn+11)
    call check_NC(nf90_put_att(ncid, latid, UNITS, LAT_UNITS), nn+12)
    call check_NC(nf90_put_att(ncid, latid, AXIS, LAT_AXS), nn+13)
    call check_NC(nf90_put_att(ncid, latid, LNG_NAME, LAT_LNG), nn+14)
    call check_NC(nf90_put_att(ncid, latid, STD_NAME, LAT_STD), nn+15)

    call check_NC(nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lonid), nn+16)
    call check_NC(nf90_put_att(ncid, lonid, UNITS, LON_UNITS), nn+17)
    call check_NC(nf90_put_att(ncid, lonid, AXIS, LON_AXS), nn+18)
    call check_NC(nf90_put_att(ncid, lonid, LNG_NAME, LON_LNG), nn+19)
    call check_NC(nf90_put_att(ncid, lonid, STD_NAME, LON_STD), nn+20)

    dimid2d = (/lon_dimid,lat_dimid,srf_dimid,time_dimid/)
    call check_NC(nf90_def_var(ncid, "land_cover", NF90_REAL, dimid2d, varid), nn+21)
    call check_NC(nf90_put_att(ncid, varid, UNITS, "fraction"), nn+22)
    call check_NC(nf90_put_att(ncid, varid, STD_NAME, "land_cover"), nn+23)
    call check_NC(nf90_put_att(ncid, varid, LNG_NAME, v_hdr), nn+24)

    dimidSrf = (/srf_len_dimid,srf_num_dimid/)
    call check_NC(nf90_def_var(ncid, SRF_ARR, NF90_CHAR, dimidSrf, srf_nm_id), nn+37)
    call check_NC(nf90_def_var(ncid, SRF_NAME, NF90_INT, srf_dimid, srfid), nn+38)
    call check_NC(nf90_put_att(ncid, srfid, UNITS, SRF_UNITS), nn+39)
    call check_NC(nf90_put_att(ncid, srfid, LNG_NAME, SRF_LNG), nn+40)
    call check_NC(nf90_put_att(ncid, srfid, STD_NAME, SRF_STD), nn+41)

    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "title", "GLEMOS model output"), nn+25)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Institution", "Meteorological Synthesizing Centre - East"), nn+26)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "institute_id", "MSC-E"), nn+27)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Contact", "MSC-E (msce@msceast.org)"), nn+28)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "project_id", "EMEP"), nn+29)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model name", MODEL_NAME), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model long name", MODEL_LONG_NAME), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Model version", MODEL_VERSION), nn+53)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-2.0"), nn+30)
    call check_NC( nf90_put_att(ncid, NF90_GLOBAL, "Grid code", GridCode), nn+31)

    call check_NC(nf90_enddef(ncid), nn+32)

! Write coordinates
    call check_NC(nf90_put_var(ncid, timeid, 0), nn+33)
    call check_NC(nf90_put_var(ncid, latid, LatMesh(Jmin:Jmax)/pi180), nn+34)
    call check_NC(nf90_put_var(ncid, lonid, LongMeshR(Imin:Imax)/pi180), nn+35)
    call check_NC(nf90_put_var(ncid, srfid, srf_ind), nn+60)
    call check_NC(nf90_put_var(ncid, srf_nm_id, SurfType(1:NumSurf)), nn+61)

! Disaggregate array - (added 31.10.18 - AG)
    do ns = 1, NumSurf
        do j = Jmin, Jmax
            do i = minI(j), maxI(j)
                var2dlc(i,j,ns) = LandCover(i,j,ns)
            end do
        end do
        call GridDisAggreg(var2dlc(Imin:Imax,Jmin:Jmax,ns),1)
    end do

! Write array
    start2d = (/1,1,1,1/)
    count2d = (/Imax-Imin+1,Jmax-Jmin+1,NumSurf,1/)
    call check_NC(nf90_put_var(ncid, varid, var2dlc(Imin:Imax,Jmin:Jmax,1:NumSurf), start2d, count2d), nn+36)

! Close file
    call check_NC(nf90_close(ncid), nn+37)

end subroutine Atm_WriteLandCoverNCF


!*******************************************************************************************************
! Subroutine writing model data at monitoring stations locations with contributions of sources
!*******************************************************************************************************

subroutine Atm_WriteMonitoringNCF(prnType)
    character(*), intent(in)    ::  prnType

    integer             s,tl, n, f, i, j, tm, src, nn
    real                tp, A11,A12,A21,A22,C11,C12,C21,C22
    integer             ncid, st, varid
    character(31)       tms
    character(15)       fp_suffix

    nn = 300
    write(YearNum,'(i4)') Year

    do s = 1, NumSubsMedia(Atm)

        select case(prnType)
            case('Hourly')
                tl = TM_HOURLY
            case('6hourly')
                tl = TM_6HOURLY

            case('Daily')
                tl = TM_DaiLY
                write(tms, '(a12,i4,a1,i2.2,a1,i2.2,a9)') "days since  ",Year,"-",Month,"-",Day," 00:00:00"
                write(DayNum,'(i2.2)') Day
                write(MonthNum,'(i2.2)') Month
                fp_suffix = YearNum//MonthNum//DayNum//'.nc'
                do src = 1, NumSrc
                ! air concentrations
                    do n = 1, sFormNum(Atm, s)
                        f = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                        var2d(Imin:Imax,Jmin:Jmax,1) = Atm_ConcDay(Imin:Imax,Jmin:Jmax,1,f,src)
                        call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,1),1)
                        do i = 1, NAir
                            A11 = AirStat(i)%Ast(1,1)
                            A12 = AirStat(i)%Ast(1,2)
                            A21 = AirStat(i)%Ast(2,1)
                            A22 = AirStat(i)%Ast(2,2)
                            C11 = var2d(AirStat(i)%iSt(1,1),AirStat(i)%jSt(1,1),1)
                            C12 = var2d(AirStat(i)%iSt(1,2),AirStat(i)%jSt(1,2),1)
                            C21 = var2d(AirStat(i)%iSt(2,1),AirStat(i)%jSt(2,1),1)
                            C22 = var2d(AirStat(i)%iSt(2,2),AirStat(i)%jSt(2,2),1)
                            sta_mod(i,src,f) = A11 * C11 + A12 * C12 + A21 * C21 + A22 * C22
                        end do ! i
                    end do ! n
                ! wet deposition
                    do n = 1, sFormNum(Atm, s)
                        f = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2d(i,j,1) = sum(WetDepDay(i,j,f,:,src))/MeshArea(i,j)
                            end do
                        end do
                        call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,1),1)
                        do i = 1, NPrec
                            A11 = PrecStat(i)%Ast(1,1)
                            A12 = PrecStat(i)%Ast(1,2)
                            A21 = PrecStat(i)%Ast(2,1)
                            A22 = PrecStat(i)%Ast(2,2)
                            C11 = var2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),1)
                            C12 = var2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),1)
                            C21 = var2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),1)
                            C22 = var2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),1)
                            stp_mod(i,src,f) = A11 * C11 + A12 * C12 + A21 * C21 + A22 * C22
                        end do ! i
                    end do ! n
                end do ! src

                var2d(Imin:Imax,Jmin:Jmax,2) = PrecipDay(Imin:Imax,Jmin:Jmax)
                call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,2),1)
                do i = 1, NPrec
                    A11 = PrecStat(i)%Ast(1,1)
                    A12 = PrecStat(i)%Ast(1,2)
                    A21 = PrecStat(i)%Ast(2,1)
                    A22 = PrecStat(i)%Ast(2,2)
                    C11 = var2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),2)
                    C12 = var2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),2)
                    C21 = var2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),2)
                    C22 = var2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),2)
                    stpr_mod(i) = A11 * C11 + A12 * C12 + A21 * C21 + A22 * C22
                end do

            case('Monthly')
                tl = TM_MONTHLY
                write(tms, '(a12,i4,a1,i2.2,a12)') "months since  ",Year,"-",Month,"-01 00:00:00"
                write(MonthNum,'(i2.2)') Month
                fp_suffix = YearNum//MonthNum//'.nc'
                do src = 1, NumSrc
                ! air concentrations
                    do n = 1, sFormNum(Atm, s)
                        f = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                        var2d(Imin:Imax,Jmin:Jmax,1) = Atm_ConcMonth(Imin:Imax,Jmin:Jmax,1,f,src)
                        call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,1),1)
                        do i = 1, NAir
                            A11 = AirStat(i)%Ast(1,1)
                            A12 = AirStat(i)%Ast(1,2)
                            A21 = AirStat(i)%Ast(2,1)
                            A22 = AirStat(i)%Ast(2,2)
                            C11 = var2d(AirStat(i)%iSt(1,1),AirStat(i)%jSt(1,1),1)
                            C12 = var2d(AirStat(i)%iSt(1,2),AirStat(i)%jSt(1,2),1)
                            C21 = var2d(AirStat(i)%iSt(2,1),AirStat(i)%jSt(2,1),1)
                            C22 = var2d(AirStat(i)%iSt(2,2),AirStat(i)%jSt(2,2),1)
                            sta_mod(i,src,f) = A11 * C11 + A12 * C12 + A21 * C21 + A22 * C22
                        end do ! i
                    end do ! n
                ! wet deposition
                    do n = 1, sFormNum(Atm, s)
                        f = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2d(i,j,1) = sum(WetDepMonth(i,j,f,:,src))/MeshArea(i,j)
                            end do
                        end do
                        call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,1),1)
                        do i = 1, NPrec
                            A11 = PrecStat(i)%Ast(1,1)
                            A12 = PrecStat(i)%Ast(1,2)
                            A21 = PrecStat(i)%Ast(2,1)
                            A22 = PrecStat(i)%Ast(2,2)
                            C11 = var2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),1)
                            C12 = var2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),1)
                            C21 = var2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),1)
                            C22 = var2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),1)
                            stp_mod(i,src,f) = A11 * C11 + A12 * C12 + A21 * C21 + A22 * C22
                        end do ! i
                    end do ! n
                end do

                var2d(Imin:Imax,Jmin:Jmax,2) = PrecipMonth(Imin:Imax,Jmin:Jmax)
                call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,2),1)
                do i = 1, NPrec
                    A11 = PrecStat(i)%Ast(1,1)
                    A12 = PrecStat(i)%Ast(1,2)
                    A21 = PrecStat(i)%Ast(2,1)
                    A22 = PrecStat(i)%Ast(2,2)
                    C11 = var2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),2)
                    C12 = var2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),2)
                    C21 = var2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),2)
                    C22 = var2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),2)
                    stpr_mod(i) = A11 * C11 + A12 * C12 + A21 * C21 + A22 * C22
                end do

            case('Yearly')
                tl = TM_YEARLY
                write(tms, '(a12,i4,a15)') "Years since ",Year,"-01-01 00:00:00"
                fp_suffix = YearNum//'.nc'
                do src = 1, NumSrc
                ! air concentrations
                    do n = 1, sFormNum(Atm, s)
                        f = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                        var2d(Imin:Imax,Jmin:Jmax,1) = Atm_ConcYear(Imin:Imax,Jmin:Jmax,1,f,src)
                        call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,1),1)
                        do i = 1, NAir
                            A11 = AirStat(i)%Ast(1,1)
                            A12 = AirStat(i)%Ast(1,2)
                            A21 = AirStat(i)%Ast(2,1)
                            A22 = AirStat(i)%Ast(2,2)
                            C11 = var2d(AirStat(i)%iSt(1,1),AirStat(i)%jSt(1,1),1)
                            C12 = var2d(AirStat(i)%iSt(1,2),AirStat(i)%jSt(1,2),1)
                            C21 = var2d(AirStat(i)%iSt(2,1),AirStat(i)%jSt(2,1),1)
                            C22 = var2d(AirStat(i)%iSt(2,2),AirStat(i)%jSt(2,2),1)
                            sta_mod(i,src,f) = A11 * C11 + A12 * C12 + A21 * C21 + A22 * C22
                        end do ! i
                    end do ! n
                ! wet deposition
                    do n = 1, sFormNum(Atm, s)
                        f = sFormInd(Atm, s, n)           ! Current form index in the arrays of forms
                        do j = Jmin, Jmax
                            do i = minI(j), maxI(j)
                                var2d(i,j,1) = sum(WetDepYear(i,j,f,:,src))/MeshArea(i,j)
                            end do
                        end do
                        call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,1),1)
                        do i = 1, NPrec
                            A11 = PrecStat(i)%Ast(1,1)
                            A12 = PrecStat(i)%Ast(1,2)
                            A21 = PrecStat(i)%Ast(2,1)
                            A22 = PrecStat(i)%Ast(2,2)
                            C11 = var2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),1)
                            C12 = var2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),1)
                            C21 = var2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),1)
                            C22 = var2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),1)
                            stp_mod(i,src,f) = A11 * C11 + A12 * C12 + A21 * C21 + A22 * C22
                        end do ! i
                    end do ! n
                end do

                var2d(Imin:Imax,Jmin:Jmax,2) = PrecipYear(Imin:Imax,Jmin:Jmax)
                call GridDisAggreg(var2d(Imin:Imax,Jmin:Jmax,2),1)
                do i = 1, NPrec
                    A11 = PrecStat(i)%Ast(1,1)
                    A12 = PrecStat(i)%Ast(1,2)
                    A21 = PrecStat(i)%Ast(2,1)
                    A22 = PrecStat(i)%Ast(2,2)
                    C11 = var2d(PrecStat(i)%iSt(1,1),PrecStat(i)%jSt(1,1),2)
                    C12 = var2d(PrecStat(i)%iSt(1,2),PrecStat(i)%jSt(1,2),2)
                    C21 = var2d(PrecStat(i)%iSt(2,1),PrecStat(i)%jSt(2,1),2)
                    C22 = var2d(PrecStat(i)%iSt(2,2),PrecStat(i)%jSt(2,2),2)
                    stpr_mod(i) = A11 * C11 + A12 * C12 + A21 * C21 + A22 * C22
                end do
        end select

! Write air conc stat data
        call Atm_WriteNCF_AirMonit(s, tms, tl, tm_pfx(tl), fp_suffix)
! Write wet dep stat data
        call Atm_WriteNCF_PrecMonit(s, tms, tl, tm_pfx(tl), fp_suffix)

    end do ! s

end subroutine Atm_WriteMonitoringNCF


subroutine Atm_WriteNCF_AirMonit(sbs, ts, itm, tprf, fpsfx)

    character(*), intent(in)  :: tprf
    integer, intent(in)       :: sbs
    integer, intent(in)       :: itm
    character(*), intent(in)  :: ts
    character(*), intent(in)  :: fpsfx

    integer                      nfid, k,i,n,f, nn, st
    character(80)                vnm, tmp_str
    integer                     start3(4), count3(4)

    nn = 400
    st= nf90_create(trim(ncfPath(itm,AC_ST))//trim(fpsfx), nf90_clobber, nfid)
    if(st /= nf90_noerr) then
        print '(/,"STOP in Atm_WriteNCF_CreateAirMonit: Cannot open file ''",a,"''",/)', &
                                    &trim(ncfPath(itm,AC_ST)//trim(fpsfx))
        stop
    endif

    call check_NC( nf90_def_dim(nfid, TIME_NAME, NF90_UNLIMITED, time_st_dimid), nn+2)
    call check_NC( nf90_def_dim(nfid, STA_NAME, NAir, sta_dimid), nn+3)
    call check_NC( nf90_def_dim(nfid, SRC_NAME, NumSrc, src_st_dimid), nn+4)

    call check_NC( nf90_def_dim(nfid, 'SrcStrLen', len(SourcID), src_len_dimid), nn+5)
    call check_NC( nf90_def_dim(nfid, 'NumSrcStr', NumSrc, src_num_dimid), nn+6)

    call check_NC( nf90_def_dim(nfid, 'StnStrLen', len(SourcID), stn_len_dimid), nn+7)
    call check_NC( nf90_def_dim(nfid, 'NumStnAirStr', NAir, sta_num_dimid), nn+8)

    call check_NC(nf90_def_var(nfid, SRC_ARR, NF90_CHAR, DimIDs=(/src_len_dimid,src_num_dimid/), VarID=src_nm_id), nn+9)
    call check_NC(nf90_def_var(nfid, TIME_NAME, NF90_REAL, time_st_dimid, time_st_id), nn+10)
    call check_NC(nf90_put_att(nfid, time_st_id, UNITS, ts), nn+11)
    call check_NC(nf90_put_att(nfid, time_st_id, CALENDAR, TIME_CLD), nn+12)
    call check_NC(nf90_put_att(nfid, time_st_id, AXIS, TIME_AXS), nn+13)
    call check_NC(nf90_put_att(nfid, time_st_id, LNG_NAME, TIME_LNG), nn+14)
    call check_NC(nf90_put_att(nfid, time_st_id, STD_NAME, TIME_STD), nn+15)

    call check_NC(nf90_def_var(nfid, SRC_NAME, NF90_INT, src_st_dimid, src_st_id), nn+16)
    call check_NC(nf90_put_att(nfid, src_st_id, UNITS, SRC_UNITS), nn+17)
    call check_NC(nf90_put_att(nfid, src_st_id, LNG_NAME, SRC_LNG), nn+18)
    call check_NC(nf90_put_att(nfid, src_st_id, STD_NAME, SRC_STD), nn+19)

    call check_NC(nf90_def_var(nfid, STA_ARR, NF90_CHAR, DimIDs=(/src_len_dimid,sta_num_dimid/), VarID=sta_nm_id), nn+20)
    call check_NC(nf90_def_var(nfid, STA_NAME, NF90_REAL, sta_dimid, sta_id), nn+21)
    call check_NC(nf90_put_att(nfid, sta_id, UNITS, "station"), nn+22)
    call check_NC(nf90_put_att(nfid, sta_id, LNG_NAME, STA_LNG), nn+23)
    call check_NC(nf90_put_att(nfid, sta_id, STD_NAME, STA_STD), nn+24)

    dimid_sta = (/sta_dimid,src_st_dimid,time_st_dimid/)

    do n = 1, sFormNum(Atm, sbs)
        f = sFormInd(Atm, sbs, n)           ! Current form index in the arrays of forms
    ! define and store substance specific attributes for 3d fields
        vnm = trim(v_nm(AC_ST))//'_'//trim(SubsID(sbs))//trim(FormID(Atm,sbs,f))
        call check_NC(nf90_def_var(nfid,vnm, NF90_REAL, dimid_sta, var_st_id(1,f)), nn+25)
        call check_NC(nf90_put_att(nfid,var_st_id(1,f),UNITS,trim(fld_unt(AC_ST))), nn+26)
        call check_NC(nf90_put_att(nfid,var_st_id(1,f),STD_NAME,vnm), nn+27)
        call check_NC(nf90_put_att(nfid,var_st_id(1,f),LNG_NAME,trim(tprf)//' '//trim(fld_nm(AC_ST))), nn+28)
    end do

    call check_NC( nf90_enddef(nfid), nn+29)

    do k = 1, NAir
        sta_nm(k) = AirStat(k)%indSt
    end do

    call check_NC(nf90_put_var(nfid, time_st_id, 0), nn+2)
    call check_NC(nf90_put_var(nfid, src_nm_id, SourcID(1:NumSrc)), nn+30)
    call check_NC(nf90_put_var(nfid, sta_nm_id, sta_nm), nn+31)

    start3 = (/1,1,1,1/)
    count3 = (/NAir,NumSrc,sFormNum(Atm, sbs),1/)
    do n = 1, sFormNum(Atm, sbs)
        f = sFormInd(Atm, sbs, n)
        call check_NC(nf90_put_var(nfid, var_st_id(1,f), sta_mod(1:NAir,1:NumSrc,f),start3,count3), nn+3)
    end do

    call check_NC(nf90_close(nfid), nn+32)

end subroutine Atm_WriteNCF_AirMonit



subroutine Atm_WriteNCF_PrecMonit(sbs, ts, itm, tprf, fpsfx)

    character(*), intent(in)  :: tprf
    integer, intent(in)       :: sbs
    integer, intent(in)       :: itm
    character(*), intent(in)  :: ts
    character(*), intent(in)  :: fpsfx

    integer                      nfid, k,i,n,f,nn, st
    character(80)                vnm, tmp_str
    integer                     start1(2), count1(2), start3(4), count3(4)

    nn = 500
    st= nf90_create(trim(ncfPath(itm,WD_ST))//trim(fpsfx), nf90_clobber, nfid)
    if(st /= nf90_noerr) then
        print '(/,"STOP in Atm_WriteNCF_CreatePrecMonit: Cannot open file ''",a,"''",/)', &
                                &trim(ncfPath(itm,WD_ST)//trim(fpsfx))
        stop
    endif

    call check_NC( nf90_def_dim(nfid, TIME_NAME, NF90_UNLIMITED, time_st_dimid), nn+2)
    call check_NC( nf90_def_dim(nfid, STP_NAME, NPrec, stp_dimid), nn+3)
    call check_NC( nf90_def_dim(nfid, SRC_NAME, NumSrc, src_st_dimid), nn+4)

    call check_NC( nf90_def_dim(nfid, 'SrcStrLen', len(SourcID), src_len_dimid), nn+5)
    call check_NC( nf90_def_dim(nfid, 'NumSrcStr', NumSrc, src_num_dimid), nn+6)

    call check_NC( nf90_def_dim(nfid, 'StnStrLen', len(SourcID), stn_len_dimid), nn+7)
    call check_NC( nf90_def_dim(nfid, 'NumStnPrecStr', NPrec, stp_num_dimid), nn+9)

    call check_NC(nf90_def_var(nfid, SRC_ARR, NF90_CHAR, DimIDs=(/src_len_dimid,src_num_dimid/), VarID=src_nm_id), nn+10)
    call check_NC(nf90_def_var(nfid, TIME_NAME, NF90_REAL, time_st_dimid, time_st_id), nn+11)
    call check_NC(nf90_put_att(nfid, time_st_id, UNITS, ts), nn+12)
    call check_NC(nf90_put_att(nfid, time_st_id, CALENDAR, TIME_CLD), nn+13)
    call check_NC(nf90_put_att(nfid, time_st_id, AXIS, TIME_AXS), nn+14)
    call check_NC(nf90_put_att(nfid, time_st_id, LNG_NAME, TIME_LNG), nn+15)
    call check_NC(nf90_put_att(nfid, time_st_id, STD_NAME, TIME_STD), nn+16)

    call check_NC(nf90_def_var(nfid, SRC_NAME, NF90_INT, src_st_dimid, src_st_id), nn+17)
    call check_NC(nf90_put_att(nfid, src_st_id, UNITS, SRC_UNITS), nn+18)
    call check_NC(nf90_put_att(nfid, src_st_id, LNG_NAME, SRC_LNG), nn+19)
    call check_NC(nf90_put_att(nfid, src_st_id, STD_NAME, SRC_STD), nn+20)

    call check_NC(nf90_def_var(nfid, STP_ARR, NF90_CHAR, DimIDs=(/src_len_dimid,stp_num_dimid/), VarID=stp_nm_id), nn+21)
    call check_NC(nf90_def_var(nfid, STP_NAME, NF90_REAL, stp_dimid, stp_id), nn+22)
    call check_NC(nf90_put_att(nfid, stp_id, UNITS, "station"), nn+23)
    call check_NC(nf90_put_att(nfid, stp_id, LNG_NAME, STP_LNG), nn+24)
    call check_NC(nf90_put_att(nfid, stp_id, STD_NAME, STP_STD), nn+25)

    call check_NC( nf90_def_var(nfid, STPR_NAME, NF90_REAL, stp_dimid, stpr_id), nn+26)
    call check_NC( nf90_put_att(nfid, stpr_id, UNITS, "station"), nn+27)
    call check_NC( nf90_put_att(nfid, stpr_id, LNG_NAME, STPR_LNG), nn+28)
    call check_NC( nf90_put_att(nfid, stpr_id, STD_NAME, STPR_STD), nn+29)

    dimid_stp = (/stp_dimid,src_st_dimid,time_st_dimid/)

    do n = 1, sFormNum(Atm, sbs)
        f = sFormInd(Atm, sbs, n)           ! Current form index in the arrays of forms
    ! define and store substance specific attributes for 3d fields
        vnm = trim(v_nm(WD_ST))//'_'//trim(SubsID(sbs))//trim(FormID(Atm,sbs,f))
        call check_NC(nf90_def_var(nfid, vnm, NF90_REAL, dimid_stp, var_st_id(2,f)), nn+30)
        call check_NC(nf90_put_att(nfid, var_st_id(2,f), UNITS, trim(fld_unt(WD_ST))), nn+31)
        call check_NC(nf90_put_att(nfid, var_st_id(2,f), STD_NAME, vnm), nn+32)
        call check_NC(nf90_put_att(nfid, var_st_id(2,f), LNG_NAME,trim(tprf)//' '//trim(fld_nm(WD_ST))), nn+33)
    end do

    call check_NC( nf90_enddef(nfid), nn+34)
    do k = 1, NPrec
        stp_nm(k) = PrecStat(k)%indSt
    end do

    call check_NC(nf90_put_var(nfid, time_st_id, 0), nn+6)
    call check_NC(nf90_put_var(nfid, src_nm_id, SourcID(1:NumSrc)), nn+35)
    call check_NC(nf90_put_var(nfid, stp_nm_id, stp_nm), nn+36)

    start3 = (/1,1,1,1/)
    count3 = (/NPrec,NumSrc,sFormNum(Atm, sbs),1/)
    do n = 1, sFormNum(Atm, sbs)
        f = sFormInd(Atm, sbs, n)           ! Current form index in the arrays of forms
        call check_NC(nf90_put_var(nfid, var_st_id(2,f), stp_mod(1:NPrec,1:NumSrc,f)), nn+7)
    end do

! Write precip stat data
    start1 = (/1,1/)
    count1 = (/NPrec,1/)
    call check_NC(nf90_put_var(nfid, stpr_id, stpr_mod(1:NPrec),start1,count1), nn+8)

    call check_NC(nf90_close(nfid), nn+37)

end subroutine Atm_WriteNCF_PrecMonit


!*******************************************************************************************************
! Subroutine for errors handling
!*******************************************************************************************************
subroutine check_NC(status,i)

	integer, intent (in) :: status,i
	character(80) str_er

	 if(status/=nf90_noerr) then
	   str_er=nf90_strerror(status)
	   print *, trim(str_er), i
	   stop 2
	 endif

end subroutine check_NC

end module Atm_NcfOutput
