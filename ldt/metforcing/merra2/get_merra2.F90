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
!
! !ROUTINE: get_merra2
! \label{get_merra2}
!
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 08 Dec 2015: James Geiger, update timing logic
!
! !INTERFACE:
subroutine get_merra2(n, findex)
! !USES:
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_metforcingMod
  use merra2_forcingMod

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly MERRA2 forcing.
!
!  The MERRA2 forcing data are organized into daily files, where each
!  daily file contains 24 one-hourly records of forcing fields.  The data
!  are considered valid at the mid-point of the hourly interval.
!
!  In general, metforcing readers read the forcing data before the current
!  time, referred to as bookend1, and after the current time, referred to as
!  bookend2.  Then the readers temporally interpolate between bookend1 and
!  bookend2.  Here, each bookend contains 24 one-hourly records of forcing,
!  and, in general, the MERRA2 reader will be temporally interpolating
!  from one hour-interval to the next hour-interval, where both hour-intervals
!  are contained in bookend1.  Issues arise between the hours 23z of one day
!  and 1z of the next day, which are complicated by the size of LDT' running
!  time-step, say 15mn, 30mn, or 1hr.
!
!  Below are some examples to illustrate the timing logic of the
!  MERRA2 reader.
!
!  \begin{verbatim}
!          ---*---|---*---|---*---|---*---|---*---|---*---|---*---|---*---
!  hour          21      22      23       0       1       2       3
!  hr_int         <---22--X--23---X--24---X---1---X---2---X---3--->
!
!
!  where:
!  hour is the hour UTC
!  hr_int is the hour-interval
!  * marks the valid point for the interval <--- hr_int --->
!
!  For example, interval 2, is from 1z to 2z, valid at 01:30z.
!  \end{verbatim}
!
!
!  First, consider the situation where the start time is
!  2005-11-01T21:00:00 and the time-step is 15mn.  Here bookend1 contains
!  01 Nov data and bookend2 contains 02 Nov data.
!
!  \begin{tabular}{lll}
!  time-step & time     & intervals                         \cr
!  1         & 21:15:00 & 21 of bookend1 and 22 of bookend1 \cr
!  2         & 21:30:00 & 22 of bookend1 and 23 of bookend1 \cr
!  3         & 21:45:00 & 22 of bookend1 and 23 of bookend1 \cr
!  4         & 22:00:00 & 22 of bookend1 and 23 of bookend1 \cr
!  5         & 22:15:00 & 22 of bookend1 and 23 of bookend1 \cr
!  6         & 22:30:00 & 23 of bookend1 and 24 of bookend1 \cr
!  7         & 22:45:00 & 23 of bookend1 and 24 of bookend1 \cr
!  8         & 23:00:00 & 23 of bookend1 and 24 of bookend1 \cr
!  9         & 23:15:00 & 23 of bookend1 and 24 of bookend1 \cr
!  10        & 23:30:00 & 24 of bookend1 and  1 of bookend2 \cr
!  11        & 23:45:00 & 24 of bookend1 and  1 of bookend2 \cr
!  12        & 00:00:00 & 24 of bookend1 and  1 of bookend2 \cr
!  13        & 00:15:00 & 24 of bookend1 and  1 of bookend2
!  \end{tabular}
!
!  At 00:30:00, the MERRA2 reader moves bookend2 to bookend1, and reads
!  03 Nov as bookend2.
!
!  \begin{tabular}{lll}
!  14        & 00:30:00 & 1 of bookend1 and  2 of bookend1 \cr
!  15        & 00:45:00 & 1 of bookend1 and  2 of bookend1 \cr
!  16        & 01:00:00 & 1 of bookend1 and  2 of bookend1 \cr
!  17        & 01:15:00 & 1 of bookend1 and  2 of bookend1 \cr
!  18        & 01:30:00 & 2 of bookend1 and  3 of bookend1
!  \end{tabular}
!
!  Next, consider a similar situation where the start time is
!  2005-11-01T21:00:00 and the time-step is 30mn.  Here bookend1 contains
!  01 Nov data and bookend2 contains 02 Nov data.
!
!  \begin{tabular}{lll}
!  time-step & time     & intervals                         \cr
!  1         & 21:30:00 & 22 of bookend1 and 23 of bookend1 \cr
!  2         & 22:00:00 & 22 of bookend1 and 23 of bookend1 \cr
!  3         & 22:30:00 & 23 of bookend1 and 24 of bookend1 \cr
!  4         & 23:00:00 & 23 of bookend1 and 24 of bookend1 \cr
!  5         & 23:30:00 & 24 of bookend1 and  1 of bookend2 \cr
!  6         & 00:00:00 & 24 of bookend1 and  1 of bookend2
!  \end{tabular}
!
!  At 00:30:00, the MERRA2 reader moves bookend2 to bookend1, and reads
!  03 Nov as bookend2.
!
!  \begin{tabular}{lll}
!  7         & 00:30:00 &  1 of bookend1 and  2 of bookend1 \cr
!  8         & 01:00:00 &  1 of bookend1 and  2 of bookend1 \cr
!  9         & 01 30:00 &  2 of bookend1 and  3 of bookend1
!  \end{tabular}
!
!  Finally, consider the situation where the start time is
!  2005-11-01T21:00:00 and the time-step is 1hr.  Here bookend1 contains
!  01 Nov data and bookend2 contains 02 Nov data.
!
!  \begin{tabular}{lll}
!  time-step & time     & intervals                         \cr
!  1         & 22:00:00 & 22 of bookend1 and 23 of bookend1 \cr
!  2         & 23:00:00 & 23 of bookend1 and 24 of bookend1 \cr
!  3         & 00:00:00 & 24 of bookend1 and  1 of bookend2
!  \end{tabular}
!
!  At 01:00:00, the MERRA2 reader moves bookend2 to bookend1, and reads
!  03 Nov as bookend2.
!
!  \begin{tabular}{lll}
!  4         & 01:00:00 &  1 of bookend1 and  2 of bookend1 \cr
!  5         & 02:00:00 &  2 of bookend1 and  3 of bookend1
!  \end{tabular}
!
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    forcing dataset index
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    call to advance or retract time
!  \item[merra2files](\ref{merra2files}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_merra2](\ref{read_merra2}) \newline
!    call to read the MERRA2 data and perform spatial interpolation
!  \end{description}
!EOP
  integer           :: order
  integer           :: ferror
  character(len=LDT_CONST_PATH_LEN)     :: slvname, flxname, lfoname, radname
  integer           :: c, r
  integer           :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  real*8            :: time1, time2, timenow
  real              :: gmt1, gmt2
  real              :: ts1, ts2

  integer           :: hr_int1, hr_int2
  integer           :: movetime   ! Flag to move bookend2 files to bookend1

! _________________________________________________________

! Please note that the timing logic has been tested only for
! these scenarios:
!
! startime of 2005-11-01T00:30:00 with time-step of 15mn
! startime of 2005-11-01T00:00:00 with time-step of 30mn
! startime of 2005-11-01T00:00:00 with time-step of 1hr
! startime of 2005-11-02T00:00:00 with time-step of 15mn

  if( LDT_rc%nts(n).gt.3600 ) then   ! > 1-hr timestep
     write(LDT_logunit,*) '[ERR]  When running LDT with MERRA-2, the clock '
     write(LDT_logunit,*) '  should run with a timestep less than or '
     write(LDT_logunit,*) '  equal to one hour.'
     call LDT_endrun()
  endif

  merra2_struc(n)%findtime1 = 0
  merra2_struc(n)%findtime2 = 0
  movetime = 0

  ! Initialize ts1 and ts2 timepoints at beginning of run:
  if ( LDT_get_nstep(LDT_rc,n) == 1 .or. LDT_rc%rstflag(n) == 1 .or. &
       merra2_struc(n)%reset_flag ) then
     merra2_struc(n)%findtime1 = 1
     merra2_struc(n)%findtime2 = 1
     LDT_rc%rstflag(n) = 0
     merra2_struc(n)%reset_flag = .false.

     yr1=LDT_rc%yr
     mo1=LDT_rc%mo
     da1=LDT_rc%da
     hr1=0
     mn1=0
     ss1=0
     if ( LDT_rc%hr == 0 .and. LDT_rc%mn < 30 ) then
        ! initialize ringtime to today at 00:30z
        ts1=30*60
     else
        ! initialize ringtime to tomorrow at 00:30z
        ts1=86400 + 30*60 ! 1 day plus 30 minutes
     endif
     call LDT_tick(merra2_struc(n)%ringtime,doy1,gmt1, &
                   yr1,mo1,da1,hr1,mn1,ss1,ts1)
  endif

  !----------------------------------------------------------
  ! Determine current time
  !----------------------------------------------------------
  yr1=LDT_rc%yr
  mo1=LDT_rc%mo
  da1=LDT_rc%da
  hr1=LDT_rc%hr
  mn1=LDT_rc%mn
  ss1=0
  call LDT_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,0.0)

  ! If current time >= time for when to move bookend values, 2 --> 1.
  if ( timenow >= merra2_struc(n)%ringtime ) then
     merra2_struc(n)%findtime2 = 1
     if ( merra2_struc(n)%findtime1 == 0 ) then
        movetime = 1
     endif

     ! reset ringtime to tomorrow at 00:30z
     yr1=LDT_rc%yr
     mo1=LDT_rc%mo
     da1=LDT_rc%da
     hr1=0
     mn1=0
     ss1=0
     ts1=86400 + 30*60 ! 1 day plus 30 minutes
     call LDT_tick(merra2_struc(n)%ringtime,doy1,gmt1, &
                   yr1,mo1,da1,hr1,mn1,ss1,ts1)
  endif

  if ( merra2_struc(n)%findtime1 == 1 ) then
     !----------------------------------------------------------
     ! Determine MERRA-2 Forcing 1 Time
     !----------------------------------------------------------
     yr1 = LDT_rc%yr
     mo1 = LDT_rc%mo
     da1 = LDT_rc%da
     hr1 = LDT_rc%hr
     mn1 = LDT_rc%mn
     ss1 = 0

     if ( hr1 == 0 .and. mn1 < 30 ) then
        ! need yesterday for bookend1
        ts1 = -60*60
     else
        ! need today for bookend1
        ts1 = 0
     endif

     call LDT_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  endif

  if ( merra2_struc(n)%findtime2 == 1 ) then
     !----------------------------------------------------------
     ! Determine MERRA-2 Forcing 2 Time
     !----------------------------------------------------------
     yr2 = LDT_rc%yr
     mo2 = LDT_rc%mo
     da2 = LDT_rc%da
     hr2 = LDT_rc%hr
     mn2 = LDT_rc%mn
     ss2 = 0

     if ( merra2_struc(n)%findtime1 == 1 ) then
        if ( hr2 == 0 .and. mn2 < 30 ) then
           ! need today for bookend2
           ts2 = 0
        else
           ! need tomorrow for bookend2
           ts2 = 24*60*60
        endif
     else
        ! LDT_rc%hr == 0
        ! need tomorrow for bookend2
        ts2 = 24*60*60
     endif

     call LDT_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
  endif

  ! Read MERRA2 - Bookend 1 files:
  if ( merra2_struc(n)%findtime1 == 1 ) then
     order = 1
     call merra2files(n,findex,merra2_struc(n)%merra2dir, yr1, mo1, da1, &
                      slvname, flxname, lfoname, radname)
     call read_merra2(n, order, findex, slvname, flxname, lfoname, radname,&
                      merra2_struc(n)%merraforc1, ferror)
  endif

  ! Read MERRA2 - Bookend 2 files (next day - store values):
  if ( merra2_struc(n)%findtime2 == 1 ) then
     ! Move bookend 2 files to bookend 1 timeframe:
     if ( movetime == 1 ) then
        merra2_struc(n)%merraforc1 = merra2_struc(n)%merraforc2
        merra2_struc(n)%merraforc2 = LDT_rc%udef
     endif

     order = 2
     call merra2files(n,findex,merra2_struc(n)%merra2dir, yr2, mo2, da2, &
                      slvname, flxname, lfoname, radname)
     call read_merra2(n, order, findex, slvname, flxname, lfoname, radname, &
                      merra2_struc(n)%merraforc2, ferror)
  endif

  if ( timenow >= merra2_struc(n)%merra2time2 ) then
     yr1 = LDT_rc%yr
     mo1 = LDT_rc%mo
     da1 = LDT_rc%da
     hr1 = LDT_rc%hr
     mn1 = 0
     ss1 = 0

     yr2 = LDT_rc%yr
     mo2 = LDT_rc%mo
     da2 = LDT_rc%da
     hr2 = LDT_rc%hr
     mn2 = 0
     ss2 = 0

     if ( LDT_rc%mn < 30 ) then
        ts1 = -30*60
        ts2 =  30*60
     else
        ts1 = 30*60
        ts2 = 90*60
     endif
     call LDT_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
     call LDT_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

     if( LDT_rc%nts(n) == 3600 ) then   ! == 1-hr timestep
        if ( LDT_rc%hr == 23 ) then
           order = 1
           hr_int1 = 23
           hr_int2 = 24
        elseif ( LDT_rc%hr == 0 ) then
           order = 2
           hr_int1 = 24
           hr_int2 = 1
        else
           order = 1
           hr_int1 = hr1+1
           hr_int2 = hr2+1
        endif

     else  ! Timesteps < 1 hour
        if (LDT_rc%hr.eq.23) then
           if (LDT_rc%mn.ge.30) then
              order = 2
              hr_int1 = 24
              hr_int2 = 1
           else    ! If at hour 23 and LDT minute < 30:
              order = 1
              hr_int1 = 23
              hr_int2 = 24
           endif
        ! For all other hours (0-22Z):
        else
           ! If at hour=0Z and minute < 30:   ! Should this be done when hourly time step run??
           if ((LDT_rc%hr.eq.0).and.(LDT_rc%mn.lt.30)) then
              order = 2
              hr_int1 = 24
              hr_int2 = 1
           ! If at any other hour (unless 0 hr .and. >= 30 minutes):
           else
              order = 1
              hr_int1 = hr1+1
              hr_int2 = hr2+1
           endif
        endif
     endif

     ! Assign MERRA2 forcing fields to two LDT time-interp placeholders (metdata1,2):
     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if (LDT_domain(n)%gindex(c,r).ne.-1) then

              if ( order == 1 ) then
                   LDT_forc(n,findex)%metdata1(:,LDT_domain(n)%gindex(c,r)) = &
                         merra2_struc(n)%merraforc1(:,hr_int1,&  ! Store hour: Current hour (same day)
                         (c+(r-1)*LDT_rc%lnc(n)))
                   LDT_forc(n,findex)%metdata2(:,LDT_domain(n)%gindex(c,r)) = &
                         merra2_struc(n)%merraforc1(:,hr_int2,&  ! Store hour:  next hour (same day)
                         (c+(r-1)*LDT_rc%lnc(n)))
              else
                   LDT_forc(n,findex)%metdata1(:,LDT_domain(n)%gindex(c,r)) = &
                         merra2_struc(n)%merraforc1(:,hr_int1,&  ! Store hour: Current hour (same day)
                         (c+(r-1)*LDT_rc%lnc(n)))
                   LDT_forc(n,findex)%metdata2(:,LDT_domain(n)%gindex(c,r)) = &
                         merra2_struc(n)%merraforc2(:,hr_int2,&  ! Store hour:  next hour (same day)
                         (c+(r-1)*LDT_rc%lnc(n)))
              endif

            endif
         enddo
      enddo

      ! Assign the hourly times:
      merra2_struc(n)%merra2time2 = time2
      merra2_struc(n)%merra2time1 = time1

  endif

end subroutine get_merra2


!BOP
! !ROUTINE: merra2files
! \label{merra2files}
!
! !INTERFACE:
subroutine merra2files(n, findex, merra2dir, yr, mo, da, slvname, flxname, lfoname, radname)

! !USES:
  use LDT_coreMod
  use LDT_logMod
!  use LDT_forecastMod

  implicit none
! !ARGUMENTS:
  integer                       :: n 
  integer                       :: findex
  character(len=*), intent(in)  :: merra2dir
  integer, intent(in)           :: yr,mo,da
  character(len=*), intent(out) :: slvname
  character(len=*), intent(out) :: flxname
  character(len=*), intent(out) :: lfoname
  character(len=*), intent(out) :: radname

! !DESCRIPTION:
!   This subroutine puts together MERRA2 file names for
!   daily netcdf files
!
!  The arguments are:
!  \begin{description}
!  \item[merra2dir]
!    Name of the MERRA2 directory
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[slvname]
!   name of the timestamped single level file
!  \item[flxname]
!   name of the timestamped flux file
!  \item[lfoname]
!   name of the timestamped land surface forcings file
!  \item[radname]
!   name of the timestamped radiation forcings file
!  \end{description}
!
!EOP

  character*4  :: cyear
  character*2  :: cmonth
  character*8  :: cdate
  character*10 :: prefix
  character*17 :: slv_spec, flx_spec, lfo_spec, rad_spec
  character*20 :: dir
  integer      :: seed
  real         :: rand

!  if(LDT_rc%forecastMode.eq.0) then !hindcast run

     write(unit=cyear, fmt='(i4.4)') yr
     write(unit=cmonth,fmt='(i2.2)') mo
     write(unit=cdate, fmt='(i4.4,i2.2,i2.2)') yr,mo,da
     
     if (yr==1979 .and. mo>=2) then
        prefix = 'MERRA2_100'
     elseif (yr>1979 .and. yr<=1991) then
        prefix = 'MERRA2_100'
     elseif ( yr >= 1992 .and. yr <= 2000 ) then   ! Since 2000 is last full year
        prefix = 'MERRA2_200'
        !  elseif ( yr >= 2001 .and. yr <= 2009 ) then
     elseif ( yr >= 2001 .and. yr <= 2010 ) then   ! Corrected for longitudinal shift of data
        prefix = 'MERRA2_300'
        !  elseif ( yr >= 2010 ) then
     elseif ( yr >= 2011 ) then
        prefix = 'MERRA2_400'
     else
        write(LDT_logunit,*) '[ERR] merra2files: date out of range'
        write(LDT_logunit,*) '  Supported years are from 1979-2-1 through ...'
        call LDT_endrun()
     endif
     
     slv_spec = '.tavg1_2d_slv_Nx.'
     flx_spec = '.tavg1_2d_flx_Nx.'
     lfo_spec = '.tavg1_2d_lfo_Nx.'
     rad_spec = '.tavg1_2d_rad_Nx.'
     
     dir = prefix//'/Y'//cyear//'/M'//cmonth
     
     ! Single level fields:
     slvname = trim(merra2dir)//'/'//dir//'/'//prefix//slv_spec//cdate//'.nc4'
     
     ! Flux fields:
     flxname = trim(merra2dir)//'/'//dir//'/'//prefix//flx_spec//cdate//'.nc4'
     
     ! Land surface forcing level:
     lfoname = trim(merra2dir)//'/'//dir//'/'//prefix//lfo_spec//cdate//'.nc4'
     
     ! Radiation fields:
     radname = trim(merra2dir)//'/'//dir//'/'//prefix//rad_spec//cdate//'.nc4'

#if 0
  else !forecast mode
     !sample yr, mo, da

     call LDT_sample_forecastDate(n, findex, yr,mo,da)

     write(unit=cyear, fmt='(i4.4)') yr
     write(unit=cmonth,fmt='(i2.2)') mo
     write(unit=cdate, fmt='(i4.4,i2.2,i2.2)') yr,mo,da
     
     if (yr==1979 .and. mo>=2) then
        prefix = 'MERRA2_100'
     elseif (yr>1979 .and. yr<=1991) then
        prefix = 'MERRA2_100'
     elseif ( yr >= 1992 .and. yr <= 2000 ) then   ! Since 2000 is last full year
        prefix = 'MERRA2_200'
        !  elseif ( yr >= 2001 .and. yr <= 2009 ) then
     elseif ( yr >= 2001 .and. yr <= 2010 ) then   ! Corrected for longitudinal shift of data
        prefix = 'MERRA2_300'
        !  elseif ( yr >= 2010 ) then
     elseif ( yr >= 2011 ) then
        prefix = 'MERRA2_400'
     else
        write(LDT_logunit,*) '[ERR] merra2files: date out of range'
        write(LDT_logunit,*) '  Supported years are from 1979-2-1 through ...'
        call LDT_endrun()
     endif
     
     slv_spec = '.tavg1_2d_slv_Nx.'
     flx_spec = '.tavg1_2d_flx_Nx.'
     lfo_spec = '.tavg1_2d_lfo_Nx.'
     rad_spec = '.tavg1_2d_rad_Nx.'
     
     dir = prefix//'/Y'//cyear//'/M'//cmonth
     
     ! Single level fields:
     slvname = trim(merra2dir)//'/'//dir//'/'//prefix//slv_spec//cdate//'.nc4'
     
     ! Flux fields:
     flxname = trim(merra2dir)//'/'//dir//'/'//prefix//flx_spec//cdate//'.nc4'
     
     ! Land surface forcing level:
     lfoname = trim(merra2dir)//'/'//dir//'/'//prefix//lfo_spec//cdate//'.nc4'
     
     ! Radiation fields:
     radname = trim(merra2dir)//'/'//dir//'/'//prefix//rad_spec//cdate//'.nc4'
  endif
#endif

end subroutine merra2files

