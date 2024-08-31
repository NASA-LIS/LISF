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
! !ROUTINE: create_nam242filename
! \label{create_nam242filename}
!
! !INTERFACE:
subroutine create_nam242filename(option, name00, name03, name06, &
     F06flag, namdir, yr, mo, da, hr )
! !USES:
  ! Dagang add begin
  use LIS_timeMgrMod, only : LIS_tick, LIS_date2time
  use LIS_coreMod,    only : LIS_rc
  use LIS_logMod,     only : LIS_logunit
  ! Dagang add end

  implicit none
! !ARGUMENTS: 
  integer,          intent(in)    :: option
  character(len=*), intent(in)    :: namdir
  integer, intent(in)             :: yr, mo, da, hr
  character(len=*), intent (out)  :: name00
  character(len=*), intent (out)  :: name03
  character(len=*), intent (out)  :: name06
  logical                         :: F06flag
! !DESCRIPTION:
!   This subroutine puts together NAM file names for 
!   different (0,3,6hr) forecast instances
! 
!  The arguments are:
!  \begin{description}
!  \item[option]
!    bookend option (1 or 2)
!  \item[namdir]
!    Name of the NAM directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[name00]
!   name of the NAM file for reading instantaneous fields
!   \item[name03]
!   name of the 3hr NAM file for reading time averaged fields
!   \item[name06]
!   name of the 6hr NAM file for reading time averaged fields
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
  character(len=8) :: fdir
  character(len=26) :: fsubs
  character(len=5),  parameter :: fsubs_prefix = 'fh.00'
  character(len=19), parameter :: fsubs_suffix = '_tl.press_gr.awp242'
  real*8      :: time1,dumbtime
  integer     :: doy1,doy
  real        :: gmt1,gmt

  ! Dagang add begin
  ! for spin-up purpose only
  real*8      :: time_base1, time_base2, time_curr
  real        :: gmt_base, gmt_curr
  integer     :: doy_base, doy_curr
  ! Dagang add end
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

#if 0
  ! Dagang add begin
  ! For spin-up purpose only

  ! concept:
  ! The current run time is compared with base time 1, which is calculated based on
  ! LIS_rc%syr+1,LIS_rc%smo,LIS_rc%sda, and LIS_rc%shr
  ! (LIS_rc%sy: start year; LIS_rc%smo: start month; start day: LIS_rc%sda; LIS_rc%shr: start hour)
  ! if the current time is larger than the base time, uyr and uyr0 is back to LIS_rc%syr
  ! or LIS_rc%syr+1 depending on the comparison between current time and base time 2

  if (LIS_rc%forc_repeat == 1) then
     if (LIS_rc%startcode == 1) then
        call LIS_date2time(time_base1,doy_base,gmt_base,LIS_rc%spinup_yr+1,LIS_rc%spinup_mo, &
                                                        LIS_rc%spinup_da,LIS_rc%spinup_hr,0,0)
        call LIS_date2time(time_base2,doy_base,gmt_base,yr,LIS_rc%spinup_mo, &
                                                           LIS_rc%spinup_da,LIS_rc%spinup_hr,0,0)
        call LIS_date2time(time_curr,doy_curr,gmt_curr,yr,mo,da,hr,0,0)

        if (time_curr .ge. time_base1) then
           if (time_curr .ge. time_base2) then
              uyr = LIS_rc%spinup_yr
              uyr0= LIS_rc%spinup_yr
           else
              uyr = LIS_rc%spinup_yr+1
              uyr0= LIS_rc%spinup_yr+1
           endif
        endif
     else
        call LIS_date2time(time_base1,doy_base,gmt_base,LIS_rc%syr+1,LIS_rc%smo,LIS_rc%sda,LIS_rc%shr,0,0)
        call LIS_date2time(time_base2,doy_base,gmt_base,yr,LIS_rc%smo,LIS_rc%sda,LIS_rc%shr,0,0)
        call LIS_date2time(time_curr,doy_curr,gmt_curr,yr,mo,da,hr,0,0)
        
        if (time_curr .ge. time_base1) then
           if (time_curr .ge. time_base2) then
              uyr = LIS_rc%syr
              uyr0= LIS_rc%syr
           else
              uyr = LIS_rc%syr+1
              uyr0= LIS_rc%syr+1
           endif
        endif
     endif
  endif
  ! Dagang add end
#endif

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
        call LIS_tick(dumbtime, doy, gmt, uyr, umo, uda, uhr, umn, uss, 24*60*60.0)
     endif
     if(uhr.eq.0) initcode0 = '00'
     if(uhr.eq.3.or.uhr.eq.6) initcode0 = '06'
     if(uhr.eq.9.or.uhr.eq.12) initcode0 = '12'
     if(uhr.eq.15.or.uhr.eq.18) initcode0 = '18'

  else !bookend 1
     if(uhr0.eq.0) then  !need to go back to the previous day 18z
        initcode1 = '18'
        call LIS_tick(dumbtime, doy, gmt, uyr0, umo0, uda0, uhr0, umn0, uss0, -24*60*60.0)
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

  write(UNIT=fdir, fmt='(i4, i2.2, i2.2)') uyr, umo, uda

  write(UNIT=fsubs, fmt='(a5, a2, a19)') fsubs_prefix, fcstcode0, fsubs_suffix

  name00 = trim(namdir) // '/' // fdir // '/' // initcode0 // '/' // fsubs

  !name03

  write(UNIT=fdir, fmt='(i4, i2.2, i2.2)') uyr0, umo0, uda0

  write(UNIT=fsubs, fmt='(a5, a2, a19)') fsubs_prefix, fcstcode1, fsubs_suffix

  name03 = trim(namdir) // '/' // fdir // '/' // initcode1 // '/' // fsubs

  !name06

  write(UNIT=fdir, fmt='(i4, i2.2, i2.2)') uyr0, umo0, uda0

  write(UNIT=fsubs, fmt='(a5, a2, a19)') fsubs_prefix, fcstcode2, fsubs_suffix

  name06 = trim(namdir) // '/' // fdir // '/' // initcode1 // '/' // fsubs

  return

end subroutine create_nam242filename

