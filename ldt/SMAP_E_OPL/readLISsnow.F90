!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: readLIS_snow
! \label{readLIS_snow}
! 
! !REVISION HISTORY: 
!  20 Feb 2022: Yonghwan Kwon, Initial Specification
! 
! !INTERFACE: 
subroutine readLIS_snow(n,yyyymmdd,hh,SnowDepth)
! !USES:
  use LDT_coreMod
  use LDT_logMod
  use LDT_smap_e_oplMod

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n 
  character*8         :: yyyymmdd
  character*2         :: hh
  real                :: SnowDepth(LDT_rc%lnc(n),LDT_rc%lnr(n))
  
!EOP 
  integer             :: c,r
  character*100       :: fname
  character*8         :: yyyymmdd_snow
  character*4         :: yyyy_snow
  character*2         :: mm_snow,dd_snow,hh_snow
  logical             :: file_exists
  integer             :: yr,mo,da,hr
  integer             :: nn
  integer             :: ierr

  SnowDepth = LDT_rc%udef

  yyyymmdd_snow = yyyymmdd
  yyyy_snow     = yyyymmdd(1:4)
  mm_snow       = yyyymmdd(5:6)
  dd_snow       = yyyymmdd(7:8)
  hh_snow       = hh

  read(yyyy_snow,*,iostat=ierr)  yr
  read(mm_snow,*,iostat=ierr) mo
  read(dd_snow,*,iostat=ierr) da
  read(hh_snow,*,iostat=ierr) hr

  nn = 0
  do while (nn.le.24)

     call create_LISsnowdepth_filename(SMAPeOPL%LISsnowdir, &
          yyyymmdd_snow, hh_snow, fname)

     inquire(file=trim(fname),exist=file_exists)  
     if(file_exists) then
        write(LDT_logunit,*) '[INFO] Reading snow outputs from ',trim(fname)
        call read_LISsnowdepth_data(n,fname,SnowDepth)
        write(LDT_logunit,*) '[INFO] Finished reading snow outputs from ',trim(fname)
 
        exit
     else
        hr = hr - 1
        if(hr.lt.0) then
           da = da - 1
           hr = 23 

           if(da.eq.0) then
              mo = mo - 1

              if(mo.eq.0) then
                 yr = yr -1
                 mo = 12
                 da = 31
              else
                 if(mo.eq.1.or.&
                  mo.eq.3.or.&
                  mo.eq.5.or.&
                  mo.eq.7.or.&
                  mo.eq.8.or.&
                  mo.eq.10.or.&
                  mo.eq.12) then
                    da = 31
                 elseif(mo.eq.2) then
                    if(mod(yr,4).eq.0) then
                       da = 29
                    else
                       da = 28
                    endif
                 else
                    da = 30
                 endif
              endif
           endif
        endif

        write(unit=yyyy_snow, fmt='(i4.4)') yr
        write(unit=mm_snow, fmt='(i2.2)') mo
        write(unit=dd_snow, fmt='(i2.2)') da
        write(unit=hh_snow, fmt='(i2.2)') hr
        yyyymmdd_snow = trim(yyyy_snow)//trim(mm_snow)//trim(dd_snow)

        nn = nn + 1
     endif
  enddo

end subroutine readLIS_snow

!BOP
! 
! !ROUTINE: read_LISsnowdepth_data
! \label{read_LISsnowdepth_data}
!
! !INTERFACE:
subroutine read_LISsnowdepth_data(n,fname,SnowDepth)
! 
! !USES:
  use LDT_logMod
  use LDT_coreMod
  use LDT_smap_e_oplMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer, intent(in)           :: n
  character (len=*)             :: fname
!EOP

  integer     :: ios, nid
  integer     :: snowdepthid             
  real        :: SnowDepth(LDT_rc%lnc(n),LDT_rc%lnr(n))

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LDT_verify(ios,'Error opening file '//trim(fname))

  ios = nf90_inq_varid(nid,'SnowDepth_tavg',snowdepthid)
  call LDT_verify(ios, 'Error nf90_inq_varid: SnowDepth_tavg')

  !values 

  ios = nf90_get_var(nid, snowdepthid, SnowDepth, &
        start=(/1,1/), &
        count=(/LDT_rc%lnc(n),LDT_rc%lnr(n)/))
  call LDT_verify(ios, 'Error nf90_get_var: SnowDepth')

  ios = nf90_close(ncid=nid)
  call LDT_verify(ios,'Error closing file '//trim(fname))

#endif

end subroutine read_LISsnowdepth_data

!BOP
! !ROUTINE: create_LISsnowdepth_filename
! \label{create_LISsnowdepth_filename}
! 
! !INTERFACE:
subroutine create_LISsnowdepth_filename(LISdir,yyyymmdd,hh,filename)
! !USES:

  implicit none
! !ARGUMENTS:
  character(len=*)  :: filename
  character (len=*) :: LISdir
  character*8       :: yyyymmdd
  character*6       :: yyyymm
  character*2       :: hh
!EOP

  yyyymm = trim(yyyymmdd(1:6))

  filename = trim(LISdir)//'/'//trim(yyyymm)//&
             '/LIS_HIST_'//trim(yyyymmdd)//trim(hh)//'00.d01.nc'

end subroutine create_LISsnowdepth_filename
