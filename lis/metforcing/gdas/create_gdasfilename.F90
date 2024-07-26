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
! !ROUTINE: create_gdasfilename
! \label{create_gdasfilename}
!
! !INTERFACE:
subroutine create_gdasfilename(option, name00, name03, name06, &
     F06flag, gdasdir, yr, mo, da, hr )
! !USES: 
  use LIS_timeMgrMod, only : LIS_tick

  implicit none
! !ARGUMENTS: 
  integer,          intent(in)    :: option
  character(len=*), intent(in)    :: gdasdir
  integer, intent(in)             :: yr, mo, da, hr
  character(len=*), intent (out)  :: name00
  character(len=*), intent (out)  :: name03
  character(len=*), intent (out)  :: name06
  logical                         :: F06flag
! !DESCRIPTION:
!   This subroutine puts together GDAS file names for 
!   different (0,3,6hr) forecast instances
!
!   LIS uses a mix of instantaneous forcing values and averaged forcing
!   values from the GDAS forcing data.  LIS processes the GDAS forcing at
!   a 3-hourly period.  When it is time to read new GDAS forcing data,
!   LIS reads 3 hours ahead.  So when LIS' clock is at hour hr, LIS
!   processes data for the period ]hr, hr+3].  Here, LIS reads instantaneous
!   forcing valid at hour hr$+3$ and averaged forcing valid for [hr, hr+3]
!
!   GDAS files ``named'' 00.f00, 06.f00, 12.f00, and 18.f00 are the
!   analysis files, and they contain only instantaneous forcing values at
!   hours 00, 06, 12, and 18 respectively.  The files hr.f03, hr.f06,
!   and hr.f09 are, respectively, 3-hour, 6-hour, and 9-hour forecasts
!   from the analysis hour hr, where hr is either 00, 06, 12, or 18.
!
!   The hr.f03 files contain instantaneous forcing values at hour hr+3, and
!   they contain forcing values averaged over the 3-hour period [hr, hr+3].
!
!   The hr.f06 files contain instantaneous forcing values at hour hr+6, and
!   they contain forcing values averaged over the 6-hour period [hr, hr+6].
!
!   The hr.f09 files contain instantaneous forcing values at hour hr+9, and
!   they contain forcing values averaged over the 3-hour period [hr+6, hr+9].
! 
!   Again, LIS processes the GDAS forcing at a 3-hourly period.  When the
!   hour of LIS' clock, hr, is one of the analysis hours (00, 06, 12, 18)
!   then LIS reads ahead to extract instantaneous forcing values from
!   the hr.f03 file, and it reads ahead to extract averaged forcing values
!   from the hr.f03 file.
!
!   When the hour of LIS' clock, hr, is not one of the analysis hours
!   (in this case when hr is 03, 09, 15, or 21), then LIS reads instantaneous
!   forcing values from the from the next analysis hour.
!   E.g., when hr is 03, read from 06.f00; when hr is 09, read from 12.f00.
!   LIS computes averaged forcing values by reading the averaged values
!   from the 3-hour and 6-hour forecasts of the previous analysis hour.
!   E.g., when hr is 03, read from 00.f03 and 00.f06; when hr is 09, read
!   from 06.f03 and 06.f06.
!   For this case, the averaged forcing values are computed by subtracting
!   the 3-hour forecast from the 6-hour forecast both of the previous analysis
!   hour thusly:
!
!   Consider LIS' hour to be 15.  LIS needs averaged forcing valid for the
!   period [15, 18].  The previous analysis hour is then 12.
!   We want the averaged forcing values from 12.f03 and from 12.f06.
!
!   hr.f06 contains averaged forcing, denoted $f_6$, for the 
!   period [hr, hr+6] and hr.f03 contains averaged forcing, denoted $f_3$,
!   for the period [hr, hr+3].  Thus for forcing, $f_6$, $6*f_6$ is the
!   total $f$ at hour hr$+6$ and $3*f_3$ is the total $f$ at hour hr$+3$.
!   $( 6*f_6 - 3*f_3 )$ represents the total $f$ from hr$+3$ to hr$+6$,
!   so $( 6*f_6 - 3*f_3 ) / 3 \equiv 2*f_6 - f_3$ is the averaged $f$ for
!   the period [hr+3, hr+6].
!   
!   So to compute the averaged forcing for period [15, 15+3], LIS must
!   read 12.f03 and 12.f06.  $f_6$ represents averaged $f$ for [12, 18];
!   $f_3$ represents averaged $f$ for [12, 15]; $2*f_6 - f_3$ represents 
!   avergage $f$ for [15, 18], the desired period.
!
!   Note that the variable, F06flag, is used to indicate the above
!   situation when LIS must subtract $f_3$ values from $f_6$ values to
!   obtain the desired averaged forcing values.
!
!   Note that the logic for reading and creating names for bookend1, which
!   is only done at initialization, is the same as described above except
!   imagine that you are at start hour - 3 reading ahead to the start hour.
!
!  The arguments are:
!  \begin{description}
!  \item[option]
!    bookend option (1 or 2)
!  \item[gdasdir]
!    Name of the GDAS directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[name00]
!   name of the GDAS file for reading instantaneous fields
!   \item[name03]
!   name of the 3hr GDAS file for reading time averaged fields
!   \item[name06]
!   name of the 6hr GDAS file for reading time averaged fields
!  \item[F06flag]
!    flag to indicate if 6hr forecast data is required for this interval
!  \end{description}
!
!EOP
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss
  integer :: uyr0, umo0, uda0, uhr0, umn0, uss0
  logical :: is_analysis_hr
  character(len=2)  :: analysis_hour_inst, analysis_hour_avg, &
                       fcstcode0, fcstcode1, fcstcode2
  character(len=80) :: fbase
  character(len=8)  :: fdir
  character(len=10) :: ftime
  character(len=21) :: fsubs
  real*8      :: dumbtime
  integer     :: doy
  real        :: gmt
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

  if ( uhr == 0 .or. uhr == 6 .or. uhr == 12 .or. uhr == 18 ) then
     is_analysis_hr = .true.
  else
     !( uhr == 3 .or. uhr == 9 .or. uhr == 15 .or. uhr == 21 )
     is_analysis_hr = .false.
  endif

!-----------------------------------------------------------------
! Cheat sheet
!
! bookend 2: 
!  hour 00 look for 00F3
!  hour 03 look for 00F3 and 00F6
!  hour 06 look for 06F3
!  hour 09 look for 06F3 and 06F6
!  hour 12 look for 12F3
!  hour 15 look for 12F3 and 12F6
!  hour 18 look for 18F3
!  hour 21 look for 18F3 and 18F6
!
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

  if ( option == 2 ) then !bookend 2
     if ( .not. is_analysis_hr ) F06flag = .true. ! need to read two files     
  else
     if ( is_analysis_hr ) F06flag = .true. ! need to read two files     
  endif

  if ( option == 2 ) then !bookend 2 
     if ( uhr == 0 .or. uhr == 6 .or. uhr == 12 .or. uhr == 18 ) then
        ! read averged forcing from analysis_hour file ( 00, 06, 12, 18)
        write(unit=analysis_hour_avg,fmt='(i2.2)') uhr

        ! read instantaneous forcing from analysis_hour file ( 00, 06, 12, 18)
        write(unit=analysis_hour_inst,fmt='(i2.2)') uhr
        fcstcode0 = '03'
     else
        ! read averaged forcing from previous analysis_hour file
        ! (03 -> 00, 09 -> 06, 15 -> 12, 21 -> 18)
        write(unit=analysis_hour_avg,fmt='(i2.2)') uhr - 3

        ! read instantaneous forcing from next analysis_hour file
        ! (03 -> 06, 09 -> 12, 15 -> 18, 21 -> 00)
        write(unit=analysis_hour_inst,fmt='(i2.2)') uhr + 3
        fcstcode0 = '00'
     endif

     ! special case: need to go to the next day
     if ( uhr .eq. 21 ) then
        call LIS_tick(dumbtime,doy,gmt,uyr,umo,uda,uhr,umn,uss,24.0*60*60)
        analysis_hour_inst = '00'
     endif

  else !bookend 1
     if ( uhr == 0 .or. uhr == 6 .or. uhr == 12 .or. uhr == 18 ) then
        ! read averged forcing from previous analysis_hour file
        ! ( 00 -> -06, 06 -> 00, 12 -> 06, 18 -> 12)
        if(uhr.ne.0) &
             write(unit=analysis_hour_avg,fmt='(i2.2)') uhr - 6

        ! read instantaneous forcing from analysis_hour file ( 00, 06, 12, 18)
        write(unit=analysis_hour_inst,fmt='(i2.2)') uhr
        fcstcode0 = '00'
     else
        ! read averaged forcing from previous analysis_hour file
        ! (03 -> 00, 09 -> 06, 15 -> 12, 21 -> 18)
        write(unit=analysis_hour_avg,fmt='(i2.2)') uhr - 3

        ! read instantaneous forcing from previous analysis_hour file
        ! (03 -> 00, 09 -> 06, 15 -> 12, 21 -> 18)
        write(unit=analysis_hour_inst,fmt='(i2.2)') uhr - 3
        fcstcode0 = '03'
     endif

     ! special case: need to go to the previous day
     if ( uhr0 == 0 ) then
        call LIS_tick(dumbtime,doy,gmt,uyr0,umo0,uda0,uhr0,umn0,uss0,-24.0*60*60)
        analysis_hour_avg = '18'
     endif

  endif

  fcstcode1 = '03'
  fcstcode2 = '06'
  

  ! construct file names

  !name00
  write(unit=fdir, fmt='(a1, i4.4, i2.2, a1)') '/', uyr, umo, '/'
  write(unit=ftime, fmt='(i4.4, i2.2, i2.2, a2)') uyr, umo, uda, analysis_hour_inst
  write(unit=fsubs, fmt='(a16, a2, a3)') '.gdas1.sfluxgrbf', fcstcode0, '.sg'
  name00 = trim(gdasdir)//fdir//ftime//fsubs
  
  !name03
  write(unit=fdir, fmt='(a1, i4.4, i2.2, a1)') '/', uyr0, umo0, '/'
  write(unit=ftime, fmt='(i4.4, i2.2, i2.2, a2)') uyr0, umo0, uda0, analysis_hour_avg
  write(unit=fsubs, fmt='(a16, a2, a3)') '.gdas1.sfluxgrbf', fcstcode1, '.sg'
  name03 = trim(gdasdir)//fdir//ftime//fsubs

  !name06
  write(unit=fdir, fmt='(a1, i4.4, i2.2, a1)') '/', uyr0, umo0, '/'
  write(unit=ftime, fmt='(i4.4, i2.2, i2.2, a2)') uyr0, umo0, uda0, analysis_hour_avg
  write(unit=fsubs, fmt='(a16, a2, a3)') '.gdas1.sfluxgrbf', fcstcode2, '.sg'
  name06 = trim(gdasdir)//fdir//ftime//fsubs

end subroutine create_gdasfilename

