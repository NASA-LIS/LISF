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
! !ROUTINE: read_VIIRSgfrac
! \label{read_VIIRSgfrac}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!  10 Aug 2010: Jonathan Case: Modified for daily SPoRT GVF composites.
!  30 Oct 2014: Jonathan Case: Modified to read gzipped SPoRT GVF composites.
!
! !INTERFACE:
subroutine read_VIIRSgfrac(n, wt1, wt2, array1, array2)
! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_logunit, LIS_getNextUnitNumber, &
                             LIS_releaseUnitNumber, LIS_endrun
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_vegDataMod, only : LIS_gfrac
  use LIS_fileIOMod,  only : LIS_readData
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

  character(len=LIS_CONST_PATH_LEN) :: filename1
  character(len=LIS_CONST_PATH_LEN) :: filename2
  logical       :: gfracAlarmCheck
  real          :: gmt
  real*8        :: ctime
  integer       :: t
  integer       :: doy
  integer       :: ftn
  integer       :: yr,mo,da,hr,mn,ss
  character*4   :: fyr
  character*2   :: fmo
  character*2   :: fda
  integer       :: rc
  integer       :: c,r,i,j

  if (LIS_gfrac(n)%firstInstance) then 
     gfracAlarmCheck = .true. 
  else
     gfracAlarmCheck = LIS_isAlarmRinging(LIS_rc,&
          "LIS gfrac read alarm",&           
          LIS_gfrac(n)%gfracIntervalType)
  endif

  if (gfracAlarmCheck) then
     
     if (LIS_gfrac(n)%firstInstance) LIS_gfrac(n)%firstInstance = .false.
     array1 = LIS_rc%udef
     array2 = LIS_rc%udef

     yr = LIS_rc%yr
     mo = LIS_rc%mo
     da = LIS_rc%da
     hr = LIS_rc%hr
     mn = LIS_rc%mn
     ss = LIS_rc%ss 

     if (LIS_gfrac(n)%realtimemode.eq.0) then
       call LIS_tick(LIS_gfrac(n)%time1,doy,gmt,yr,mo,da,hr,mn,ss,0.0)
     else
       call LIS_tick(LIS_gfrac(n)%time1,doy,gmt,yr,mo,da,hr,mn,ss,-2*86400.0)
     endif
     write(unit=fyr,fmt='(i4.4)') yr
     write(unit=fmo,fmt='(i2.2)') mo
     write(unit=fda,fmt='(i2.2)') da
     filename1 = trim(LIS_gfrac(n)%gfracfile)//"."//&
          trim(fyr)//trim(fmo)//trim(fda)//".1gd4r"

     call LIS_tick(LIS_gfrac(n)%time2,doy,gmt,yr,mo,da,hr,mn,ss,86400.0)
     write(unit=fyr,fmt='(i4.4)') yr
     write(unit=fmo,fmt='(i2.2)') mo
     write(unit=fda,fmt='(i2.2)') da
     filename2 = trim(LIS_gfrac(n)%gfracfile)//"."//&
          trim(fyr)//trim(fmo)//trim(fda)//".1gd4r"        

     call get_VIIRSgfrac(n,filename1,array1)
     call get_VIIRSgfrac(n,filename2,array2)

  endif

  !compute weights
  if (LIS_gfrac(n)%realtimemode.eq.0) then
    call LIS_tick(ctime,doy,gmt,LIS_rc%yr,LIS_rc%mo,LIS_rc%da,&
         LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,0.0)
  else
    call LIS_tick(ctime,doy,gmt,LIS_rc%yr,LIS_rc%mo,LIS_rc%da,&
         LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,-2*86400.0)
  endif
  
  wt1 = (LIS_gfrac(n)%time2-ctime)/(LIS_gfrac(n)%time2-LIS_gfrac(n)%time1)
  wt2 = (ctime-LIS_gfrac(n)%time1)/(LIS_gfrac(n)%time2-LIS_gfrac(n)%time1)

!  write (LIS_logunit,*) 'CTIME/GVFT1/GVFT2/WT1/WT2 = ',ctime,LIS_gfrac(n)%time1,LIS_gfrac(n)%time2,wt1,wt2

! reset ctime back to the current time, if realtimemode is used.
  if (LIS_gfrac(n)%realtimemode.eq.1) then
    call LIS_tick(ctime,doy,gmt,LIS_rc%yr,LIS_rc%mo,LIS_rc%da,&
           LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,2*86400.0)
  endif

end subroutine read_VIIRSgfrac


subroutine get_VIIRSgfrac(n,filename, array)
  
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use LIS_vegDataMod,     only : LIS_gfrac
  use LIS_logMod
  use VIIRSgfracMod

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! J.Case; Need to make gvf_nc/gvf_nr dynamic based on the input lat/lon and
! delta lat lon values in lis.config
!  integer, parameter     :: gvf_nc = 12000, gvf_nr = 5000
!  integer                :: gvf_nc, gvf_nr
  integer, intent(in)    :: n
  real                   :: array(LIS_rc%ntiles(n))
  integer                :: ftn
  character(len=*)       :: filename
  integer                :: c,r,t,i,j,stc,enc,str,enr
  integer                :: npts
  integer                :: iret
  real, allocatable      :: gvf(:)
  logical*1, allocatable :: li(:)
  real, allocatable      :: gvfinterp(:)
  logical*1, allocatable :: lo(:)
  logical                :: file_exists

! J.Case (10/29/2014) -- variables for reading gzipped GVF file.
  character(len=LIS_CONST_PATH_LEN) :: zname
  character*100 :: flag
  integer :: iretgz = 1   ! iretgz = 0 from readgvfviirs is "good"
  integer :: readgvfviirs

  external readgvfviirs

! J.Case (10/29/2014) Implemented reading of gzipped GVF files, as in LIS6.
  zname = trim(filename)//".gz"
  inquire(file=zname,exist=file_exists)

  if (file_exists) then

    write(LIS_logunit,*) 'Reading GVF file ',trim(zname)
    !gvf_nc = int(LIS_rc%gridDesc(n,32))
    !gvf_nr = int(LIS_rc%gridDesc(n,33))
    npts = gvf_nc*gvf_nr
    allocate(gvf(npts))
    gvf = 0.0
    flag = '1'
    write (LIS_logunit,*) "Before readgvfviirs; gvf_nc / gvf_nr / npts = ",gvf_nc,gvf_nr,npts
    iretgz = readgvfviirs( trim(zname)//char(0), gvf_nc, gvf_nr, flag, gvf )
    if ( iretgz > 0) then       ! FAILED TO READ
      write (LIS_logunit,*) "** Failed to read GVF gzipped file: ", trim(zname)
      deallocate(gvf)
      call LIS_endrun()
    endif
    write (LIS_logunit,*) "Successful readgvfviirs; gvf_nc / gvf_nr / npts = ",gvf_nc,gvf_nr,npts

  else ! read standard uncompressed real 4-byte binary data

    inquire(file=filename,exist=file_exists) 
    if(file_exists) then 
      write(LIS_logunit,*) 'Reading GVF file ',trim(filename)
      ftn = LIS_getNextUnitNumber()
      !gvf_nc = int(LIS_rc%gridDesc(n,32))
      !gvf_nr = int(LIS_rc%gridDesc(n,33))
      npts = gvf_nc*gvf_nr
      allocate(gvf(npts))
      gvf = 0.0
      open(ftn,file=filename,form='unformatted',access='direct',recl=npts*4)
      read(ftn,rec=1) gvf
      call LIS_releaseUnitNumber(ftn)
    else
      write(LIS_logunit,*) 'Could not open file: ',trim(filename)
      call LIS_endrun()
    endif

  endif  ! check for GVF zname or filename
! J.Case (10/29/2014) -- end of mods to read gzipped GVF file.

  allocate(li(npts))
  li = .true. 
  do i=1,npts
     if(gvf(i).lt.0) li(i) = .false.
  enddo
     
  allocate(gvfinterp(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

  call bilinear_interp(LIS_rc%gridDesc(n,:), li, gvf, lo, &
       gvfinterp, npts, LIS_rc%lnc(n)*LIS_rc%lnr(n), &
       LIS_domain(n)%lat, LIS_domain(n)%lon, &
       LIS_gfrac(n)%w111, LIS_gfrac(n)%w121, &
       LIS_gfrac(n)%w211, LIS_gfrac(n)%w221, &
       LIS_gfrac(n)%n111, LIS_gfrac(n)%n121, &
       LIS_gfrac(n)%n211, LIS_gfrac(n)%n221, LIS_rc%udef, iret)
     
  do t=1,LIS_rc%ntiles(n)
     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     array(t) = gvfinterp(c+(r-1)*LIS_rc%lnc(n))
! J.Case test (let's try just commenting out the neighbor search below)
! If the code below is allowed to run, data are all 0.20 everywhere in LIS domain.
!     if(array(t).lt.0) then 
!        str = max(r-2,1)
!        enr = min(r+2,LIS_rc%lnr(n))
!        stc = max(c-2,1)
!        enc = min(c+2,LIS_rc%lnc(n))
!        do j=str,enr
!           do i=stc,enc
!              if(gvfinterp(i+(j-1)*LIS_rc%lnc(n)).ne.-9999) then 
!                 array(t) = gvfinterp(i+(j-1)*LIS_rc%lnc(n))
!                 exit
!              endif
!           enddo
!        enddo
!     endif
!     ! if the neighbor search fails, then fill in a fixed value for now. 
!     ! Need a better strategy. 
!     if(array(t).lt.0) then 
!        array(t) = 0.20
!     endif
  enddo
     
  deallocate(gvf)
  deallocate(li)
  deallocate(gvfinterp)
  deallocate(lo)

end subroutine get_VIIRSgfrac

