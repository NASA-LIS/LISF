!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_NESDISgfrac
! \label{read_NESDISgfrac}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!
! !INTERFACE:
subroutine read_NESDISgfrac(n, wt1, wt2, array1, array2)
! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use LIS_vegDataMod,     only : LIS_gfrac
  use LIS_timeMgrMod

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n
  real                :: wt1
  real                :: wt2
  real, intent(inout) :: array1(LIS_rc%ntiles(n))
  real, intent(inout) :: array2(LIS_rc%ntiles(n))

! !DESCRIPTION:
!  This subroutine retrieves the greenness fraction climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[mo]
!   time index (month or quarter)
!  \item[array]
!   output field with the retrieved greenness fraction
!  \end{description}
!
!EOP      

  character(len=LIS_CONST_PATH_LEN) :: filename1,filename2
  logical           :: gfracAlarmCheck
  logical           :: file_exists
  integer           :: i
  integer           :: yr,mo,da,hr,mn,ss,doy
  real              :: gmt
  real*8            :: ctime
  character*1       :: fyr(4)
  character(len=4)  :: fyr1
  character(len=3)  :: fdoy
  character*40      :: temp


! First check for the daily alarm. If alarm rings then 
! see if file exists for the current day. If exists, then read
! current day's file (time1) and the file 7 days later (time2). 
! 
! If file not read, then compute the weights based on current 
! time and the two bookends. 
! 
  if(LIS_gfrac(n)%firstInstance) then 
     gfracAlarmCheck = .true. 
  else
     gfracAlarmCheck = LIS_isAlarmRinging(LIS_rc,&
          "LIS gfrac read alarm",&           
       LIS_gfrac(n)%gfracIntervalType)
  endif

  if(gfracAlarmCheck) then

     if(LIS_gfrac(n)%firstInstance) then 
        LIS_gfrac(n)%firstInstance = .false. 

        array1 = LIS_rc%udef
        array2 = LIS_rc%udef

        yr = LIS_rc%yr
        mo = LIS_rc%mo
        da = LIS_rc%da
        hr = LIS_rc%hr
        mn = LIS_rc%mn
        ss = LIS_rc%ss 

        call LIS_tick(ctime,doy,gmt,yr,mo,da,hr,mn,ss,0.0)

        write(unit=temp,fmt='(I4)') yr
        read(unit=temp,fmt='(4a1)') (fyr(i),i=1,4)
        write(unit=fyr1,fmt='(i4.4)') yr
        write(unit=fdoy,fmt='(i3.3)') doy

        filename2 = trim(LIS_gfrac(n)%gfracfile)//"/"//trim(fyr1)&
             //'/NPR.VGWG.D'//trim(fyr(3))//trim(fyr(4))//trim(fdoy)//".GRIBF"

        inquire(file=trim(filename2), exist=file_exists)

        do while(.not.file_exists) 
           ! advance time until a file is found.            
           call LIS_tick(ctime,doy,gmt,yr,mo,da,hr,mn,ss,86400.0)

           write(unit=temp,fmt='(I4)') yr
           read(unit=temp,fmt='(4a1)') (fyr(i),i=1,4)
           write(unit=fyr1,fmt='(i4.4)') yr
           write(unit=fdoy,fmt='(i3.3)') doy

           filename2 = trim(LIS_gfrac(n)%gfracfile)//"/"//trim(fyr1)&
                //'/NPR.VGWG.D'//trim(fyr(3))//trim(fyr(4))//trim(fdoy)//".GRIBF"

           inquire(file=trim(filename2), exist=file_exists)
           if(file_exists) then       
              exit;
           endif
        enddo

        call LIS_tick(LIS_gfrac(n)%time2,doy,gmt,yr,mo,da,hr,mn,ss,0.0)
        write(unit=temp,fmt='(I4)') yr
        read(unit=temp,fmt='(4a1)') (fyr(i),i=1,4)
        write(unit=fyr1,fmt='(i4.4)') yr
        write(unit=fdoy,fmt='(i3.3)') doy
        filename2 = trim(LIS_gfrac(n)%gfracfile)//"/"//trim(fyr1)&
             //'/NPR.VGWG.D'//trim(fyr(3))//trim(fyr(4))//trim(fdoy)//".GRIBF"

        call LIS_tick(LIS_gfrac(n)%time1,doy,gmt,yr,mo,da,hr,mn,ss,-7*86400.0)
        write(unit=temp,fmt='(I4)') yr
        read(unit=temp,fmt='(4a1)') (fyr(i),i=1,4)
        write(unit=fyr1,fmt='(i4.4)') yr
        write(unit=fdoy,fmt='(i3.3)') doy
        filename1 = trim(LIS_gfrac(n)%gfracfile)//"/"//trim(fyr1)&
             //'/NPR.VGWG.D'//trim(fyr(3))//trim(fyr(4))//trim(fdoy)//".GRIBF"

        call get_NESDISgfrac(n,filename2,array2)
        call get_NESDISgfrac(n,filename1,array1)


     else
        write(unit=temp,fmt='(I4)') LIS_rc%yr
        read(unit=temp,fmt='(4a1)') (fyr(i),i=1,4)
        write(unit=fyr1,fmt='(i4.4)') LIS_rc%yr
        write(unit=fdoy,fmt='(i3.3)') LIS_rc%doy

        filename1 = trim(LIS_gfrac(n)%gfracfile)//"/"//trim(fyr1)&
             //'/NPR.VGWG.D'//trim(fyr(3))//trim(fyr(4))//trim(fdoy)//".GRIBF"

        inquire(file=trim(filename1), exist=file_exists)
        if(file_exists) then 
           array1 = LIS_rc%udef
           array2 = LIS_rc%udef

           call get_NESDISgfrac(n,filename1,array1)

           yr = LIS_rc%yr
           mo = LIS_rc%mo
           da = LIS_rc%da
           hr = LIS_rc%hr
           mn = LIS_rc%mn
           ss = LIS_rc%ss

           call LIS_tick(LIS_gfrac(n)%time1,doy,gmt,yr,mo,da,hr,mn,ss,0.0)
           call LIS_tick(LIS_gfrac(n)%time2,doy,gmt,yr,mo,da,hr,mn,ss,7.0*86400)

           write(unit=temp,fmt='(I4)') yr
           read(unit=temp,fmt='(4a1)') (fyr(i),i=1,4)
           write(unit=fyr1,fmt='(i4.4)') yr
           write(unit=fdoy,fmt='(i3.3)') doy
           filename2 = trim(LIS_gfrac(n)%gfracfile)//"/"//trim(fyr1)&
                //'/NPR.VGWG.D'//trim(fyr(3))//trim(fyr(4))//trim(fdoy)//".GRIBF"

           inquire(file=trim(filename2), exist=file_exists)
           if(.not.file_exists) then 
              write(LIS_logunit,*) 'GVF file ',trim(filename2), ' does not exist'
              call LIS_endrun()
           endif

           call get_NESDISgfrac(n,filename2,array2)

        endif
     endif
  end if
  
  !compute weights
  call LIS_tick(ctime,doy,gmt,LIS_rc%yr,LIS_rc%mo,LIS_rc%da,&
       LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,0.0)     
  
  wt2 = (ctime-LIS_gfrac(n)%time1)/(LIS_gfrac(n)%time2-LIS_gfrac(n)%time1)
  wt1 = (LIS_gfrac(n)%time2-ctime)/(LIS_gfrac(n)%time2-LIS_gfrac(n)%time1)

end subroutine read_NESDISgfrac

subroutine get_NESDISgfrac(n,filename, array)
  
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_vegDataMod,     only : LIS_gfrac
  use LIS_logMod

#if (defined USE_GRIBAPI)
  use grib_api
#endif


  implicit none

  integer, intent(in) :: n 
  character(len=*)    :: filename
  real                :: array(LIS_rc%ntiles(n))

  integer           :: iret
  integer           :: ftn 
  integer           :: nc,nr,npts
  integer           :: rc
  real, allocatable :: gvf(:)
  logical*1, allocatable :: li(:)
  real, allocatable :: gvfinterp(:)
  logical*1, allocatable :: lo(:)
  integer           :: i,j,t,c,r
  integer           :: str,enr,stc,enc
  real              :: gmt
  integer           :: yr,mo,da,hr,mn,ss,doy
  logical           :: file_exists
  real              :: missingValue
  real*8            :: t1
  integer           :: pds5_val, pds7_val
  integer           :: igrib
  

  allocate(gvfinterp(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

  ftn = 133
  nc = 2500
  nr = 1250
  npts = nc*nr 
  
#if (defined USE_GRIBAPI) 
  write(LIS_logunit,*) 'Reading GVF file ',trim(filename)

  call grib_open_file(ftn,trim(filename),'r',iret)
  call LIS_verify(iret, 'error grib_open_file in read_NESDISgfrac')
  
  call grib_new_from_file(ftn,igrib,iret)
  call LIS_verify(iret,'error in grib_new_from_file in read_NESDISgfrac')
  
  
  call grib_get(igrib,'indicatorOfParameter',pds5_val,iret)
  call LIS_verify(iret, 'error in grib_get: indicatorOfParameter in read_NESDISgfrac')
  
  call grib_get(igrib,'level',pds7_val,iret)
  call LIS_verify(iret, 'error in grib_get: level in read_NESDISgfrac')
  
  if(pds5_val.eq.87.and.pds7_val.eq.0) then 
     !Adding 8 because the internal file dimensions are inconsistent with the 
     ! overall domain dimensions. 
     allocate(gvf(npts+8))
     allocate(li(npts))
     
     call grib_get(igrib,'values',gvf,iret)
     call LIS_verify(iret, 'error in grib_get:values in read_NESDISgfrac')
     
     call grib_get(igrib,'missingValue',missingValue,iret)
     call LIS_verify(iret, 'error in grib_get:missingValue in read_NESDISgfrac')
     
     li = .false. 
     do i=1,npts
        if(gvf(i).ne.missingValue) li(i) = .true.
     enddo
     
     
     call bilinear_interp(LIS_rc%gridDesc(n,:), li, gvf, lo, &
          gvfinterp, npts, LIS_rc%lnc(n)*LIS_rc%lnr(n), &
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          LIS_gfrac(n)%w111, LIS_gfrac(n)%w121, &
          LIS_gfrac(n)%w211, LIS_gfrac(n)%w221, &
          LIS_gfrac(n)%n111, LIS_gfrac(n)%n121, &
          LIS_gfrac(n)%n211, LIS_gfrac(n)%n221, LIS_rc%udef, iret)
     
     
     do t=1,LIS_rc%ntiles(n)
        c = LIS_domain(n)%tile(t)%col
        r = LIS_domain(n)%tile(t)%row
        array(t) = gvfinterp(c+(r-1)*LIS_rc%lnc(n))
        if(array(t).lt.0) then 
           str = max(r-2,1)
           enr = min(r+2,LIS_rc%lnr(n))
           stc = max(c-2,1)
           enc = min(c+2,LIS_rc%lnc(n))
           do j=str,enr
              do i=stc,enc
                 if(gvfinterp(i+(j-1)*LIS_rc%lnc(n)).ne.-9999) then 
                    array(t) = gvfinterp(i+(j-1)*LIS_rc%lnc(n))
                    exit
                 endif
              enddo
           enddo
        endif
        ! if the neighbor search fails, then fill in a fixed value for now. 
        ! Need a better strategy. 
        if(array(t).lt.0) then 
           array(t) = 0.20
        endif
     enddo
     deallocate(gvf)
     deallocate(li)
     
     call grib_release(igrib,iret)
     call LIS_verify(iret,'error in grib_release in read_NESDISgfrac')
     
  else
     write(LIS_logunit,*) 'Could not retrieve entries in file: ',trim(filename)
     call LIS_endrun()
     
  endif
  
  call grib_close_file(ftn)
#endif
  
  deallocate(gvfinterp)
  deallocate(lo)

end subroutine get_NESDISgfrac
