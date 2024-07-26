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
#include "LDT_NetCDF_inc.h"
!BOP
!
! !ROUTINE: readUSAFSI
! \label{readUSAFSI}
!
! !REVISION HISTORY:
!  13 Feb 2023: Eric Kemp, Initial Specification
!
! !INTERFACE:
subroutine readUSAFSI(n, yyyymmdd, hh, SnowDepth, rc)

! !USES:
  use LDT_coreMod
  use LDT_logMod
  use LDT_smap_e_oplMod

  implicit none

! !ARGUMENTS:
  integer, intent(in)     :: n
  character*8, intent(in) :: yyyymmdd
  character*2, intent(in) :: hh
  real, intent(out):: SnowDepth(LDT_rc%lnc(n),LDT_rc%lnr(n))
  integer, intent(out) :: rc

!EOP
  character*255       :: fname
  character*8         :: yyyymmdd_snow
  character*4         :: yyyy_snow
  character*2         :: mm_snow, dd_snow, hh_snow
  logical             :: file_exists
  integer             :: yr, mo, da, hr
  integer             :: nn
  integer             :: ierr
  integer :: rc1

  external :: create_USAFSI_filename
  external :: read_USAFSI_data

  rc = 1 ! Default to error, will update below if USAFSI file read in.
  SnowDepth = LDT_rc%udef

  yyyymmdd_snow = yyyymmdd
  yyyy_snow     = yyyymmdd(1:4)
  mm_snow       = yyyymmdd(5:6)
  dd_snow       = yyyymmdd(7:8)
  hh_snow       = hh

  read(yyyy_snow, *, iostat=ierr) yr
  read(mm_snow, *, iostat=ierr) mo
  read(dd_snow, *, iostat=ierr) da
  read(hh_snow, *, iostat=ierr) hr

  nn = 0
  do while (nn <= 24)

     call create_USAFSI_filename(SMAPeOPL%LISsnowdir, &
          yyyymmdd_snow, hh_snow, fname)
     inquire(file=trim(fname), exist=file_exists)

     if (file_exists) then
        write(LDT_logunit,*) '[INFO] Reading snow depth from ', trim(fname)
        call read_USAFSI_data(n, fname, SnowDepth, rc1)
        if (rc1 == 0) then
           write(LDT_logunit,*) '[INFO] Finished reading snow outputs from ', &
                trim(fname)
           rc = 0
           exit
        end if
     end if

     ! Go back one hour
     hr = hr - 1
     if (hr < 0) then
        da = da - 1
        hr = 23
        if (da == 0) then
           mo = mo - 1
           if (mo == 0) then
              yr = yr - 1
              mo = 12
              da = 31
           else
              if (mo == 1 .or. &
                   mo == 3 .or. &
                   mo == 5 .or. &
                   mo == 7 .or. &
                   mo == 8 .or. &
                   mo == 10 .or. &
                   mo == 12) then
                 da = 31
              elseif (mo == 2) then
                 if (mod(yr,4) == 0) then
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
     yyyymmdd_snow = trim(yyyy_snow) // trim(mm_snow) // trim(dd_snow)

     nn = nn + 1

  end do
end subroutine readUSAFSI

!BOP
!
! !ROUTINE: read_USAFSI_data
! \label{read_USAFSI_data}
!
! !INTERFACE:
subroutine read_USAFSI_data(n, fname, SnowDepth, rc)
!
! !USES:
  use LDT_logMod
  use LDT_coreMod
  use LDT_smap_e_oplMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS:
!
  integer, intent(in)      :: n
  character(*), intent(in) :: fname
  real, intent(inout)      :: SnowDepth(LDT_rc%lnc(n),LDT_rc%lnr(n))
  integer, intent(out)     :: rc
!EOP

  integer     :: ios, nid
  integer     :: snowdepthid

  rc = 1 ! Initialize as error, reset near bottom

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname), mode=NF90_NOWRITE, ncid=nid)
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Error opening file' // trim(fname)
     return
  end if

  ios = nf90_inq_varid(nid, 'snoanl', snowdepthid)
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Cannot find snoanl in ' // trim(fname)
     ios = nf90_close(ncid=nid)
     return
  end if

  ios = nf90_get_var(nid, snowdepthid, SnowDepth, &
        start=(/1, 1/), &
        count=(/LDT_rc%lnc(n), LDT_rc%lnr(n)/))
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Cannot read snoanl in ' // trim(fname)
     ios = nf90_close(ncid=nid)
     return
  end if

  ios = nf90_close(ncid=nid)
  if (ios .ne. 0) then
     write(LDT_logunit,*)'[WARN] Error closing ' // trim(fname)
     return
  end if

  rc = 0 ! No error detected!

#endif

end subroutine read_USAFSI_data

!BOP
! !ROUTINE: create_USAFSI_filename
! \label{create_USAFSI_filename}
!
! !INTERFACE:
subroutine create_USAFSI_filename(LISdir, yyyymmdd, hh, filename)
! !USES:

  implicit none

! !ARGUMENTS:
  character(*), intent(in)  :: LISdir
  character(8), intent(in)  :: yyyymmdd
  character(2), intent(in)  :: hh
  character(*), intent(out) :: filename
!EOP

  filename = trim(LISdir) // '/USAFSI_' // trim(yyyymmdd) &
       // trim(hh) // '.nc'

end subroutine create_USAFSI_filename

