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
! !ROUTINE: create_nam242f9_filename
! \label{create_nam242f9_filename}
!
! !INTERFACE:
subroutine create_nam242f9_filename(option, name00, name03, &
     namdir, yr, mo, da, hr,status )
! !USES: 
  ! Dagang add begin
  use LDT_timeMgrMod, only : LDT_tick, LDT_date2time
  use LDT_coreMod,    only : LDT_rc
  use LDT_logMod,     only : LDT_logunit
  ! Dagang add end

  implicit none
! !ARGUMENTS: 
  integer,          intent(in)    :: option
  character(len=*), intent(in)    :: namdir
  integer, intent(in)             :: yr, mo, da, hr
  character(len=*), intent (out)  :: name00
  character(len=*), intent (out)  :: name03
  integer                         :: status
! !DESCRIPTION:
!   This subroutine puts together 9hr forecast NAM file name
! 
!  The arguments are:
!  \begin{description}
!  \item[option]
!    bookend option (1 or 2)
!   \item[name00]
!   name of the timestamped 9hr NAM file
!   \item[name03]
!   name of the timestamped 9hr NAM file
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
!  \item[status]
!   flag indicating if the 9hr forecast is available. 
!  \end{description}
!
!EOP
  integer   :: i, c
  integer   :: uyr, umo, uda, uhr, umn, uss, ts1
  integer   :: uyr0, umo0, uda0, uhr0, umn0, uss0
  integer   :: remainder
  character(len=2) :: initcode0, initcode1, fcstcode0, fcstcode1, fcstcode2
  character*1 :: fbase(80), fdir(13), ftime(10), fsubs(26)
  character(LEN=100) :: temp1,temp2
  real*8    :: time1,dumbtime
  integer   :: doy1,doy
  real      :: gmt1,gmt

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
  ! LDT_rc%syr+1,LDT_rc%smo,LDT_rc%sda, and LDT_rc%shr
  ! (LDT_rc%sy: start year; LDT_rc%smo: start month; start day: LDT_rc%sda; LDT_rc%shr: start hour)
  ! if the current time is larger than the base time, uyr and uyr0 is back to LDT_rc%syr
  ! or LDT_rc%syr+1 depending on the comparison between current time and base time 2

  if (LDT_rc%forc_repeat == 1) then
     if (LDT_rc%startcode == 1) then
        call LDT_date2time(time_base1,doy_base,gmt_base,LDT_rc%spinup_yr+1,LDT_rc%spinup_mo, &
                                                        LDT_rc%spinup_da,LDT_rc%spinup_hr,0,0)
        call LDT_date2time(time_base2,doy_base,gmt_base,yr,LDT_rc%spinup_mo, &
                                                           LDT_rc%spinup_da,LDT_rc%spinup_hr,0,0)
        call LDT_date2time(time_curr,doy_curr,gmt_curr,yr,mo,da,hr,0,0)

        if (time_curr .ge. time_base1) then
           if (time_curr .ge. time_base2) then
              uyr = LDT_rc%spinup_yr
              uyr0= LDT_rc%spinup_yr
           else
              uyr = LDT_rc%spinup_yr+1
              uyr0= LDT_rc%spinup_yr+1
           endif
        endif
     else
        call LDT_date2time(time_base1,doy_base,gmt_base,LDT_rc%syr+1,LDT_rc%smo,LDT_rc%sda,LDT_rc%shr,0,0)
        call LDT_date2time(time_base2,doy_base,gmt_base,yr,LDT_rc%smo,LDT_rc%sda,LDT_rc%shr,0,0)
        call LDT_date2time(time_curr,doy_curr,gmt_curr,yr,mo,da,hr,0,0)
        
        if (time_curr .ge. time_base1) then
           if (time_curr .ge. time_base2) then
              uyr = LDT_rc%syr
              uyr0= LDT_rc%syr
           else
              uyr = LDT_rc%syr+1
              uyr0= LDT_rc%syr+1
           endif
        endif
     endif 
  endif
  ! Dagang add end
#endif

  remainder = modulo(uhr,2)  !if even, then remainder equals zero
                             !if odd, then remainder equals one
!-----------------------------------------------------------------
!  hour 03 look for 18F9
!  hour 09 look for 00F9
!  hour 15 look for 06F9
!  hour 21 look for 12F9
!-----------------------------------------------------------------
  if(option.eq.1) then !bookend 1
     if(uhr.eq.03.or.uhr.eq.09.or.uhr.eq.15.or.uhr.eq.21) then 
        if(uhr0.eq.3) then  !need to go back to the previous day 18z
           initcode1 = '18'
           call LDT_tick(dumbtime, doy, gmt, uyr0, umo0, uda0, uhr0, umn0, uss0, -9*60*60.0)
        endif
        
        if(uhr.eq.9)  initcode1 = '00'
        if(uhr.eq.15) initcode1 = '06'
        if(uhr.eq.21) initcode1 = '12'
     else
        status = 1
     endif
  else
      if(uhr.eq.00.or.uhr.eq.06.or.uhr.eq.12.or.uhr.eq.18) then 
        if(uhr0.eq.0) then  !need to go back to the previous day 18z
           initcode1 = '18'
           call LDT_tick(dumbtime, doy, gmt, uyr0, umo0, uda0, uhr0, umn0, uss0, -6*60*60.0)
        endif
        
        if(uhr.eq.6)  initcode1 = '00'
        if(uhr.eq.12) initcode1 = '06'
        if(uhr.eq.18) initcode1 = '12'
     else
        status = 1
     endif
  endif

  fcstcode1 = '09'
   
  !name 00
  write(UNIT=temp1, fmt='(a40)') namdir  
  read(UNIT=temp1, fmt='(80a1)') (fbase(i), i=1,80)

  write(UNIT=temp1, fmt='(a1, i4, i2, i2, a1, a2, a1)') '/', uyr0, umo0, uda0, '/', initcode1, '/'
  read(UNIT=temp1, fmt='(13a1)') fdir
  do i = 1, 13
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do

  write(UNIT=temp1, fmt='(a5, a2, a19)') 'fh.00', fcstcode1, '_tl.press_gr.awp242'
  read (UNIT=temp1, fmt='(80a1)') (fsubs(i), i=1,26)

  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(UNIT=temp1, fmt='(80a1)') (fbase(i), i=1,c), (fdir(i), i=1,13), (fsubs(i), i=1,26)
  read(UNIT=temp1, fmt='(a80)') name00

  !namef03
  write(UNIT=temp2, fmt='(a40)') namdir  
  read(UNIT=temp2, fmt='(80a1)') (fbase(i), i=1,80)

  write(UNIT=temp2, fmt='(a1, i4, i2, i2, a1, a2, a1)') '/', uyr0, umo0, uda0, '/', initcode1, '/'
  read(UNIT=temp2, fmt='(13a1)') fdir
  do i = 1, 13
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do

  write(UNIT=temp2, fmt='(a5, a2, a19)') 'fh.00', fcstcode1, '_tl.press_gr.awp242'
  read (UNIT=temp2, fmt='(80a1)') (fsubs(i), i=1,26)

  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(UNIT=temp2, fmt='(80a1)') (fbase(i), i=1,c), (fdir(i), i=1,13), (fsubs(i), i=1,26)
  read(UNIT=temp2, fmt='(a80)') name03

  return

end subroutine create_nam242f9_filename

