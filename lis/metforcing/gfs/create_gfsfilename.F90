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
! !ROUTINE: create_gfsfilename
! \label{create_gfsfilename}
!
! !REVISION HISTORY:
!  16 Mar 2008: Sujay Kumar; Initial specification
!
! !INTERFACE:
subroutine create_gfsfilename(option, name00, name03, name06, &
     F06flag, gfsdir, yr, mo, da, hr )
! !USES: 
  use LIS_timeMgrMod, only : LIS_tick

  implicit none
! !ARGUMENTS: 
  integer,          intent(in)    :: option
  character(len=*), intent(in)    :: gfsdir
  integer, intent(in)             :: yr, mo, da, hr
  character(len=*), intent (out)  :: name00
  character(len=*), intent (out)  :: name03
  character(len=*), intent (out)  :: name06
  logical                         :: F06flag
! !DESCRIPTION:
!   This subroutine puts together GFS file names for 
!   different (0,3,6hr) forecast instances
! 
!  The arguments are:
!  \begin{description}
!  \item[option]
!    bookend option (1 or 2)
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
!   \item[name00]
!   name of the GFS file for reading instantaneous fields
!   \item[name03]
!   name of the 3hr GFS file for reading time averaged fields
!   \item[name06]
!   name of the 6hr GFS file for reading time averaged fields
!  \item[F06flag]
!    flag to indicate if 6hr forecast data is required for this interval
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
  character(len=3),  parameter :: fsubsext    = '.sg'
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
! bookend 2: 
!  hour 00 look for 00F3
!  hour 03 look for 00F3 and 00F6
!  hour 06 look for 06F3
!  hour 09 look for 06F3 and 06F6
!  hour 12 look for 12F3
!  hour 15 look for 12F3 and 12F6
!  hour 18 look for 18F3
!  hour 21 look for 18F3 and 18F6
! bookend 1: 
!  hour 00 look for 18F3 and 18F6
!  hour 03 look for 00F3
!  hour 06 look for 00F3 and 00F6
!  hour 09 look for 06F3
!  hour 12 look for 06F3 and 06F6
!  hour 15 look for 12F3
!  hour 18 look for 12F3 and 12F6
!  hour 21 look for 18F3
!-----------------------------------------------------------------
  F06flag = .false. 
  if(option.eq.2) then !bookend 2
     if(remainder.gt.0) F06flag = .true. ! need to read two files     
  else
     if(remainder.eq.0) F06flag = .true. ! need to read two files     
  endif

  if(option.eq.2) then !bookend 2 

     if(uhr.eq.0.or.uhr.eq.3) initcode1 = '00'
     if(uhr.eq.6.or.uhr.eq.9) initcode1 = '06'
     if(uhr.eq.12.or.uhr.eq.15) initcode1 = '12'
     if(uhr.eq.18.or.uhr.eq.21) initcode1 = '18'

     if(uhr.eq.21) then !need to go to the next day
        initcode0 = '00'
        call LIS_tick(dumbtime, doy, gmt, uyr, umo, uda, uhr, umn, uss, 24.0*60*60)
     endif
     if(uhr.eq.0) initcode0 = '00'
     if(uhr.eq.3.or.uhr.eq.6) initcode0 = '06'
     if(uhr.eq.9.or.uhr.eq.12) initcode0 = '12'
     if(uhr.eq.15.or.uhr.eq.18) initcode0 = '18'

  else !bookend 1
     if(uhr0.eq.0) then  !need to go back to the previous day 18z
        initcode1 = '18'
        call LIS_tick(dumbtime, doy, gmt, uyr0, umo0, uda0, uhr0, umn0, uss0, -24.0*60*60)
     endif

     if(uhr.eq.21) initcode1 = '18'
     if(uhr.eq.3.or.uhr.eq.6) initcode1 = '00'
     if(uhr.eq.9.or.uhr.eq.12) initcode1 = '06'
     if(uhr.eq.15.or.uhr.eq.18) initcode1 = '12'

     if(uhr.eq.0.or.uhr.eq.3) initcode0 = '00'
     if(uhr.eq.6.or.uhr.eq.9) initcode0 = '06'
     if(uhr.eq.12.or.uhr.eq.15) initcode0 = '12'
     if(uhr.eq.18.or.uhr.eq.21) initcode0 = '18'
  endif
 
  if(option.eq.2) then !bookend 2
     if(remainder.gt.0) then 
        fcstcode0 = '00'
     else
        fcstcode0 = '03'
     endif
  else
     if(remainder.gt.0) then 
        fcstcode0 = '03'
     else
        fcstcode0 = '00'
     endif
  endif

  fcstcode1 = '03'
  fcstcode2 = '06'
  
  !name00

  write(UNIT=fdir, fmt='(i4, i2.2)') uyr, umo

  write(UNIT=ftime, fmt='(i4, i2.2, i2.2, a2)') uyr, umo, uda, initcode0

  fsubs = fsubsprefix // fcstcode0 // fsubsext

  name00 = trim(gfsdir) // '/' // fdir // '/' // ftime // fsubs
  
  !name 03

  write(UNIT=fdir, fmt='(i4, i2.2)') uyr, umo

  write(UNIT=ftime, fmt='(i4, i2.2, i2.2, a2)') uyr0, umo0, uda0, initcode1

  fsubs = fsubsprefix // fcstcode1 // fsubsext

  name03 = trim(gfsdir) // '/' // fdir // '/' // ftime // fsubs

  !namef06

  write(UNIT=fdir, fmt='(i4, i2.2)') uyr, umo

  write(UNIT=ftime, fmt='(i4, i2.2, i2.2, a2)') uyr0, umo0, uda0, initcode1

  fsubs = fsubsprefix // fcstcode2 // fsubsext

  name06 = trim(gfsdir) // '/' // fdir // '/' // ftime // fsubs

  return

end subroutine create_gfsfilename

