!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: create_gfs_f0backup_filename
! \label{create_gfs_f0backup_filename}
!
! !INTERFACE:
subroutine create_gfs_f0backup_filename(option, name00, gfsdir, yr, mo, da, hr, status)
! !USES: 
  use LIS_timeMgrMod, only : LIS_tick

  implicit none
! !ARGUMENTS: 
  integer,          intent(in)    :: option
  character(len=*), intent(in)    :: gfsdir
  integer, intent(in)             :: yr, mo, da, hr
  character(len=*), intent (out)  :: name00
  integer                         :: status
! !DESCRIPTION:
!   This subroutine puts together the backup file for f00 forecast
! 
!  The arguments are:
!  \begin{description}
!  \item[option]
!    bookend option (1 or 2)
!   \item[name00]
!   name of the timestamped F00 backup file
!  \item[gfsdir]
!    Name of the GFS directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!  \item[status]
!   flag indicating if the filename was generated successfully
!  \end{description}
!
!EOP
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer :: uyr0, umo0, uda0, uhr0, umn0, uss0
  integer :: remainder
  character(len=2) :: initcode0, initcode1, fcstcode0, fcstcode1, fcstcode2
  character(len=6)  :: fdir
  character(len=10) :: ftime
  character(len=21) :: fsubs
  character(len=14), parameter :: fsubsprefix = '.gfs.sfluxgrbf'
  character(len=3),  parameter :: fsubsext = '.sg'
  real*8      :: time1,dumbtime
  integer     :: doy1,doy
  real        :: gmt1,gmt
!=== End Variable Definition ===============

!=== formats for filename segments
!-----------------------------------------------------------------
!  Make variables for the time used to create the file
!  We don't want these variables being passed out
!-----------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0

  uyr0 = yr
  umo0 = mo
  uda0 = da
  uhr0 = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn0 = 0
  uss0 = 0

  remainder = modulo(uhr,2)  !if even, then remainder equals zero
                             !if odd, then remainder equals one
!-----------------------------------------------------------------
!  hour 00 look for 18F6
!  hour 06 look for 00F6
!  hour 12 look for 06F6
!  hour 18 look for 12F6
!-----------------------------------------------------------------
  if(option.eq.1) then 
     if(uhr.eq.00.or.uhr.eq.06.or.uhr.eq.12.or.uhr.eq.18) then 
        if(uhr0.eq.0) then  !need to go back to the previous day 18z
           initcode1 = '18'
           call LIS_tick(dumbtime, doy, gmt, uyr0, umo0, uda0, uhr0, umn0, uss0, -6.0*60*60)
        endif
        
        if(uhr.eq.6)  initcode1 = '00'
        if(uhr.eq.12) initcode1 = '06'
        if(uhr.eq.18) initcode1 = '12'
        
        status = 0 
     else
        status = 1
     endif
  elseif(option.eq.2) then 
     if(uhr.eq.03.or.uhr.eq.09.or.uhr.eq.15.or.uhr.eq.21) then 
        if(uhr.eq.3)  initcode1 = '00'
        if(uhr.eq.9)  initcode1 = '06'
        if(uhr.eq.15) initcode1 = '12'
        if(uhr.eq.21) initcode1 = '18'
        
        status = 0 
     else
        status = 1
     endif
  endif

  fcstcode1 = '06'
   
  !name 00

  write(UNIT=fdir, fmt='(i4, i2.2)') uyr0, umo0

  write(UNIT=ftime, fmt='(i4, i2.2, i2.2, a2)') uyr0, umo0, uda0, initcode1

  fsubs = fsubsprefix // fcstcode1 // fsubsext

  name00 = trim(gfsdir) // '/' // fdir // '/' // ftime // fsubs

  return

end subroutine create_gfs_f0backup_filename

