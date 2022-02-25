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
! Module of the model output
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module TextOutputProc

    use GeneralParams
    use Geometry
    use Atm_Params
    
    implicit none
 
contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine writing 2d fields to text file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WriteFieldTxt(path, capt, fld, med, sbs, nforms, pr)

    real, intent(in)            :: fld(Imin:Imax,Jmin:Jmax,MaxForm)    ! Output field per forms
    character*(*), intent(in)   :: path                                ! Output path and filename
    character*(*), intent(in)   :: capt                                ! File header string
    integer, intent(in)         :: med                                 ! Media index
    integer, intent(in)         :: sbs                                 ! Substance index
    integer, intent(in), optional  :: nforms
    integer, intent(in), optional  :: pr

    character(3)    rep
    character(120)  frmt1, frmt2, frmt3
    integer         f, i, j, num, ind, FileStat                         


! Set number of forms
    if(present(nforms)) then
        num = nforms
    else
        num = sFormNum(med, sbs)
    end if

! Disaggregate fields
    do f = 1, num
        call GridDisAggreg(fld(:,:,sFormInd(med,sbs,f)), 1)
    end do

! write fields to file
    open(3, file=path, action='write', iostat=FileStat)             
    if(FileStat>0) then                                                                 
        print '(/,"STOP in WriteFieldTxt: Cannot open file ''",a,"''",/)', trim(path)   
        stop                                                                            
    endif                                                                               
    write(rep,'(i2)') num
    frmt1='("  i    j        long         lat  ",'//trim(adjustl(rep))//&
	    &'a15,/,"--------------------",'&
            &//trim(adjustl(rep))//'a15)'
    frmt2='(i3,2x,i3,2(1x,e11.4),'//trim(adjustl(rep))//'(2x,e13.6))'
    frmt3='("'//trim(adjustl(capt))//' '//trim(GridCode)//'",/)'
    write(3,frmt3)

    if(present(pr)) then
        write(3,frmt1) 'Precip',('---------------')
    elseif(sFormNum(med, sbs) == 1) then
        write(3,frmt1) trim(SubsID(sbs)),('---------------')
    else
        write(3,frmt1) (trim(SubsID(sbs))//trim(FormID(med,sbs,f)),f=1,num),('---------------',f=1,num)
    end if

    do j=Jmin, Jmax
        do i=Imin, Imax
            write(3,frmt2) i, j, LongMeshR(i)/pi180, LatMesh(j)/pi180, (fld(i,j,sFormInd(med,sbs,f)), f=1, num)
        enddo
    enddo
    close(3)

end subroutine WriteFieldTxt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine writing modelled source-receptor matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RTYPE==2

subroutine WriteMatrix(path, capt, matr)
    real(8), intent(in)            :: matr(1:NumSrc,1:NumRcp)    ! Output matrix array
    character*(*), intent(in)   :: path                       ! Output path and filename
    character*(*), intent(in)   :: capt                       ! File header string
    character(3)    rep
    character(120)  frmt1, frmt2, frmt3
    integer         i, j, Rcp, Src, FileStat                    

    open(13, file=path, action='write', iostat=FileStat)        
    if(FileStat>0) then                                           
        print '(/,"STOP in WriteMatrix: Cannot open file ''",a,"''",/)', trim(path)     
        stop                                                                            
    endif                                                                               
    write(rep,'(i2)') NumSrc
    frmt1='(8x,a8,'//trim(adjustl(rep))//'(6x,a8))'
    frmt2='(a6,2x,e12.6,'//trim(adjustl(rep))//'(2x,e12.6))'
    frmt3='("'//trim(adjustl(capt))//'",/)'
    write(13,frmt3)

    write(13,frmt1) 'RecArea',(trim(SourcID(Src)), Src=1, NumSrc)
    do Rcp=1, NumRcp
        write(13,frmt2) RecepID(Rcp), RecepArea(Rcp)*1.0e-6, &
    			&(matr(Src,Rcp), Src=1, NumSrc)		! Receptor area in km2
    enddo
    close(13)

end subroutine WriteMatrix

#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine writing modelled data at station locations to the text files
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine WriteMonitor(path, arr_stat, nst, tm_str, ind)

    character*(*), intent(in)       :: path                 ! Output path and filename
    type(station), intent(in)       :: arr_stat(NstatMax)   ! Array of station structures
    integer, intent(in)             :: nst                  ! Number of stations
    integer, intent(in)             :: ind                  ! Substance form index
    character*(*), intent(in)       :: tm_str               ! File header string

    character(120)  frmt1
    character(3)    rep
    integer         i, FileStat                             ! Corrected AG 19.01.17
    real            A11, A12, A21, A22, C11, C12, C21, C22, mod(NstatMax)

! Set format string
    write(rep,'(i3)') nst
    frmt1='(a17,'//trim(adjustl(rep))//'(2x,e10.4))'
! Write to file in the append mode
    open(4, file=trim(path), action='write',status='old',access='append', iostat=FileStat)  
    if(FileStat>0) then                                                                     
        print '(/,"STOP in WriteMonitor: Cannot open file ''",a,"''",/)', trim(path)        
        stop                                                                                
    endif                                                                                   

    do i = 1, nst
        A11 = arr_stat(i)%Ast(1,1)
        A12 = arr_stat(i)%Ast(1,2)
        A21 = arr_stat(i)%Ast(2,1)
        A22 = arr_stat(i)%Ast(2,2)
        C11 = arr_stat(i)%stMod(1,1,ind)
        C12 = arr_stat(i)%stMod(1,2,ind)
        C21 = arr_stat(i)%stMod(2,1,ind)
        C22 = arr_stat(i)%stMod(2,2,ind)
        mod(i) = A11 * C11 + A12 * C12 + A21 * C21 + A22 * C22
    end do

    write(4,frmt1) tm_str, (mod(i), i=1, nst)
    close(4)
    
end subroutine WriteMonitor


subroutine WriteMonitorHeader(path, capt, arr_stat, nst)

    character*(*), intent(in)   :: path                 ! Output path and filename
    character*(*), intent(in)   :: capt                 ! File header string
    type(station), intent(in)   :: arr_stat(NstatMax)   ! Array of station structures
    integer, intent(in)         :: nst                  ! Number of stations

    character(3)    rep
    character(120)  frmt1, frmt2
    integer         i, FileStat                         


! write fields to file
    open(5, file=path, action='write', iostat=FileStat)             
    if(FileStat>0) then                                             
        print '(/,"STOP in WriteMonitorHeader: Cannot open file ''",a,"''",/)', trim(path)  
        stop                                                                                
    endif                                                                                   
    frmt1='("'//trim(adjustl(capt))//'",/)'
    write(5,frmt1)

    write(rep,'(i3)') nst
    frmt2='(" Year Mm Dd Hr ",'//trim(adjustl(rep))//'a12,/,"---------------",'&
            &//trim(adjustl(rep))//'a12)'

    write(5,frmt2) (trim(arr_stat(i)%indSt),i=1,nst),('------------',i=1,nst)
    close(5)

end subroutine WriteMonitorHeader



end module TextOutputProc
