!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_merraland
! \label{get_merraland}
!
!
! !REVISION HISTORY:
!  12 Oct 2009: Eric Kemp, initial version
!  27 May 2010: David Mocko, changed to hourly forcing for LIS6.0
!  25 Oct 2010: David Mocko, updated for LIS6.1 public release
!   5 Apr 2013: Sujay Kumar, updated for the direct use of files from GES-DISC
!  31 May 2013: David Mocko, changes to match MERRA-Land interpolation strategy
!  11 Nov 2015: KR Arsenault, fix to hourly time interval
!  10 Dec 2015: James Geiger, update timing logic
!
! !INTERFACE:
subroutine get_merraland(n,findex)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_metforcingMod
  use merraland_forcingMod

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly MERRA-Land forcing.
!
!  The MERRA-Land forcing data are organized into daily files, where each
!  daily file contains 24 one-hourly records of forcing fields.  The data
!  are considered valid at the mid-point of the hourly interval.
!
!  In general, metforcing readers read the forcing data before the current
!  time, referred to as bookend1, and after the current time, referred to as
!  bookend2.  Then the readers temporally interpolate between bookend1 and
!  bookend2.  Here, each bookend contains 24 one-hourly records of forcing,
!  and, in general, the MERRA-Land reader will be temporally interpolating
!  from one hour-interval to the next hour-interval, where both hour-intervals
!  are contained in bookend1.  Issues arise between the hours 23z of one day
!  and 1z of the next day, which are complicated by the size of LIS' running
!  time-step, say 15mn, 30mn, or 1hr.
!
!  Below are some examples to illustrate the timing logic of the
!  MERRA-Land reader.
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
!  At 00:30:00, the MERRA-Land reader moves bookend2 to bookend1, and reads
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
!  At 00:30:00, the MERRA-Land reader moves bookend2 to bookend1, and reads
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
!  At 01:00:00, the MERRA-Land reader moves bookend2 to bookend1, and reads
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
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    call to advance or retract time
!  \item[merralandfiles](\ref{merralandfiles}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_merraland](\ref{read_merraland}) \newline
!    call to read the MERRA-Land data and perform spatial interpolation
!  \end{description}
!EOP
  integer           :: c,r
  integer           :: order
  integer           :: ferror
  character*100     :: slvname, flxname, radname, mldname, lndname
  integer           :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  real*8            :: time1, time2, timenow
  real              :: gmt1,gmt2,ts1,ts2

  integer           :: hr_int1, hr_int2
  integer           :: movetime

! __________________________________________________________

! Please note that the timing logic has been tested only for
! these scenarios:
!
! startime of 2005-11-01T00:00:00 with time-step of 15mn
! startime of 2005-11-01T00:00:00 with time-step of 30mn
! startime of 2005-11-01T00:00:00 with time-step of 1hr

  if( LIS_rc%nts(n).gt.3600 ) then   ! > 1-hr timestep
     write(LIS_logunit,*) 'ERR: When running LIS with MERRA-Land, the clock '
     write(LIS_logunit,*) 'should run with a timestep less than or '
     write(LIS_logunit,*) 'equal to one hour.'
     call LIS_endrun()
  endif

  merraland_struc(n)%findtime1 = 0
  merraland_struc(n)%findtime2 = 0
  movetime = 0

  if ( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1 ) then
     merraland_struc(n)%findtime1 = 1
     merraland_struc(n)%findtime2 = 1
     LIS_rc%rstflag(n) = 0

     yr1=LIS_rc%yr
     mo1=LIS_rc%mo
     da1=LIS_rc%da
     hr1=0
     mn1=0
     ss1=0
     if ( LIS_rc%hr == 0 .and. LIS_rc%mn < 30 ) then
        ! initialize ringtime to today at 00:30z
        ts1=30*60
     else
        ! initialize ringtime to tomorrow at 00:30z
        ts1=86400 + 30*60 ! 1 day plus 30 minutes
     endif
     call LIS_tick(merraland_struc(n)%ringtime,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  endif

  !----------------------------------------------------------
  ! Determine current time
  !----------------------------------------------------------
  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,0.0)

  if ( timenow >= merraland_struc(n)%ringtime ) then
     merraland_struc(n)%findtime2 = 1
     if ( merraland_struc(n)%findtime1 == 0 ) then
        movetime = 1
     endif

     ! reset ringtime to tomorrow at 00:30z
     yr1=LIS_rc%yr
     mo1=LIS_rc%mo
     da1=LIS_rc%da
     hr1=0
     mn1=0
     ss1=0
     ts1=86400 + 30*60 ! 1 day plus 30 minutes
     call LIS_tick(merraland_struc(n)%ringtime,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  endif

  if ( merraland_struc(n)%findtime1 == 1 ) then
     !----------------------------------------------------------
     ! Determine MERRA-2 Forcing 1 Time
     !----------------------------------------------------------
     yr1 = LIS_rc%yr
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = LIS_rc%hr
     mn1 = LIS_rc%mn
     ss1 = 0

     if ( hr1 == 0 .and. mn1 < 30 ) then
        ! need yesterday for bookend1
        ts1 = -60*60
     else
        ! need today for bookend1
        ts1 = 0
     endif

     call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  endif

  if ( merraland_struc(n)%findtime2 == 1 ) then
     !----------------------------------------------------------
     ! Determine MERRA-2 Forcing 2 Time
     !----------------------------------------------------------
     yr2 = LIS_rc%yr
     mo2 = LIS_rc%mo
     da2 = LIS_rc%da
     hr2 = LIS_rc%hr
     mn2 = LIS_rc%mn
     ss2 = 0

     if ( merraland_struc(n)%findtime1 == 1 ) then
        if ( hr2 == 0 .and. mn2 < 30 ) then
           ! need today for bookend2
           ts2 = 0
        else
           ! need tomorrow for bookend2
           ts2 = 24*60*60
        endif
     else
        ! LIS_rc%hr == 0
        ! need tomorrow for bookend2
        ts2 = 24*60*60
     endif

     call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
  endif

  if ( merraland_struc(n)%findtime1 == 1 ) then
     order = 1
     call merralandfiles(merraland_struc(n)%merralanddir, yr1, mo1, da1, &
                         slvname, flxname, radname, mldname, lndname)
     call read_merraland(n,order,findex,slvname,flxname,&
                         radname,mldname,lndname,merraland_struc(n)%merraforc1,ferror)
  endif

  if ( merraland_struc(n)%findtime2 == 1 ) then
     if ( movetime == 1 ) then
        merraland_struc(n)%merraforc1 = merraland_struc(n)%merraforc2
        merraland_struc(n)%merraforc2 = LIS_rc%udef
     endif

     order = 2
     call merralandfiles(merraland_struc(n)%merralanddir,yr2,mo2,da2, &
                         slvname, flxname, radname, mldname, lndname)
     call read_merraland(n,order,findex,slvname,flxname,&
                         radname,mldname,lndname,merraland_struc(n)%merraforc2,ferror)
  endif

  if ( timenow >= merraland_struc(n)%merralandtime2 ) then

     yr1 = LIS_rc%yr
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = LIS_rc%hr
     mn1 = 0
     ss1 = 0

     yr2 = LIS_rc%yr
     mo2 = LIS_rc%mo
     da2 = LIS_rc%da
     hr2 = LIS_rc%hr
     mn2 = 0
     ss2 = 0

     if ( LIS_rc%mn < 30 ) then
        ts1 = -30*60
        ts2 =  30*60
     else
        ts1 = 30*60
        ts2 = 90*60
     endif

     call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
     call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

     if( LIS_rc%nts(n) == 3600 ) then   ! == 1-hr timestep
        if ( LIS_rc%hr == 23 ) then
           order = 1
           hr_int1 = 23
           hr_int2 = 24
        elseif ( LIS_rc%hr == 0 ) then
           order = 2
           hr_int1 = 24
           hr_int2 = 1
        else
           order = 1
           hr_int1 = hr1+1
           hr_int2 = hr2+1
        endif
     else  ! Timesteps < 1 hour
        if (LIS_rc%hr.eq.23) then
           if (LIS_rc%mn.ge.30) then
              order = 2
              hr_int1 = 24
              hr_int2 = 1
           else    ! If at hour 23 and LIS minute < 30:
              order = 1
              hr_int1 = 23
              hr_int2 = 24
           endif
        ! For all other hours (0-22Z):
        else
           ! If at hour=0Z and minute < 30:   ! Should this be done when hourly time step run??
           if ((LIS_rc%hr.eq.0).and.(LIS_rc%mn.lt.30)) then
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

  ! Assign MERRA-Land forcing fields to two LIS time-interp placeholders (metdata1,2):
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r).ne.-1) then

              if ( order == 1 ) then
                   merraland_struc(n)%metdata1(:,LIS_domain(n)%gindex(c,r)) = &
                         merraland_struc(n)%merraforc1(:,hr_int1,&  ! Store hour: Current hour (same day)
                         (c+(r-1)*LIS_rc%lnc(n)))
                   merraland_struc(n)%metdata2(:,LIS_domain(n)%gindex(c,r)) = &
                         merraland_struc(n)%merraforc1(:,hr_int2,&  ! Store hour:  next hour (same day)
                         (c+(r-1)*LIS_rc%lnc(n)))
              else
                   merraland_struc(n)%metdata1(:,LIS_domain(n)%gindex(c,r)) = &
                         merraland_struc(n)%merraforc1(:,hr_int1,&  ! Store hour: Current hour (same day)
                         (c+(r-1)*LIS_rc%lnc(n)))
                   merraland_struc(n)%metdata2(:,LIS_domain(n)%gindex(c,r)) = &
                         merraland_struc(n)%merraforc2(:,hr_int2,&  ! Store hour:  next hour (same day)
                         (c+(r-1)*LIS_rc%lnc(n)))
              endif

            endif
         enddo
      enddo

  ! Assign the hourly times:
    merraland_struc(n)%merralandtime1 = time1
    merraland_struc(n)%merralandtime2 = time2

  endif

end subroutine get_merraland


!BOP
! !ROUTINE: merralandfiles
! \label{merralandfiles}
!
! !INTERFACE:
subroutine merralandfiles(merralanddir,yr,mo,da, &
     slvname,flxname,radname,mldname,lndname)

! !USES:
  use LIS_logMod, only : LIS_endrun

  implicit none
! !ARGUMENTS:
  character(len=*), intent(in)  :: merralanddir
  integer, intent(in)           :: yr,mo,da
  character(len=*), intent(out) :: slvname
  character(len=*), intent(out) :: flxname
  character(len=*), intent(out) :: radname
  character(len=*), intent(out) :: mldname
  character(len=*), intent(out) :: lndname

! !DESCRIPTION:
!   This subroutine puts together MERRA-land file names for
!   daily netcdf files
!
!  The arguments are:
!  \begin{description}
!  \item[merralanddir]
!    Name of the MERRA-Land directory
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[slvname]
!   name of the timestamped met forcing file
!  \item[flxname]
!   name of the timestamped precip flux file
!  \item[radname]
!   name of the timestamped radidation file
!  \end{description}
!
!EOP

  character*6  :: ftime1
  character*8  :: ftime2
  character*10 :: prefix
  logical      :: file_exists

  write(unit=ftime1, fmt='(i4.4,i2.2)') yr,mo
  write(unit=ftime2, fmt='(i4.4,i2.2,i2.2)') yr,mo,da

  if(yr.ge.1993.and.yr.lt.2001) then
     prefix = 'MERRA20'
  elseif(yr.ge.2001) then
     prefix = 'MERRA30'
  else
     prefix = 'MERRA10'
  endif

! Single Level ("slv") input file:
  slvname = trim(merralanddir)//'/'//trim(ftime1)//&
       '/'//trim(prefix)//'0.prod.assim.tavg1_2d_slv_Nx.'//&
       trim(ftime2)//'.SUB.nc'

  inquire( file=slvname, exist=file_exists )
  if( .NOT. file_exists ) then
    slvname = trim(merralanddir)//'/'//trim(ftime1)//&
        '/'//trim(prefix)//'1.prod.assim.tavg1_2d_slv_Nx.'//&
        trim(ftime2)//'.SUB.nc'
    inquire( file=slvname, exist=file_exists )
    if( .NOT. file_exists ) then
       print *, " ## ERR MSG:  SLV file missing: ",trim(slvname)
       print *, " LIS endrun being called ..."
       call LIS_endrun
     endif
  endif

! Radiation fields ("rad") input file:
  radname = trim(merralanddir)//'/'//trim(ftime1)//&
       '/'//trim(prefix)//'1.prod.assim.tavg1_2d_rad_Nx.'//&
       trim(ftime2)//'.SUB.nc'

  inquire( file=radname, exist=file_exists )
  if( .NOT. file_exists ) then
    radname = trim(merralanddir)//'/'//trim(ftime1)//&
       '/'//trim(prefix)//'0.prod.assim.tavg1_2d_rad_Nx.'//&
       trim(ftime2)//'.SUB.nc'
    inquire( file=radname, exist=file_exists )
    if( .NOT. file_exists ) then
       print *, " ## ERR MSG:  RAD file missing: ",trim(radname)
       print *, " LIS endrun being called ..."
       call LIS_endrun
    endif
  endif

! Flux fields ("flx") input file:
  flxname = trim(merralanddir)//'/'//trim(ftime1)//&
       '/'//trim(prefix)//'0.prod.assim.tavg1_2d_flx_Nx.'//&
       trim(ftime2)//'.SUB.nc'

  inquire( file=flxname, exist=file_exists )
  if( .NOT. file_exists ) then
    flxname = trim(merralanddir)//'/'//trim(ftime1)//&
       '/'//trim(prefix)//'1.prod.assim.tavg1_2d_flx_Nx.'//&
       trim(ftime2)//'.SUB.nc'
    inquire( file=flxname, exist=file_exists )
    if( .NOT. file_exists ) then
       print *, " ## ERR MSG:  FLX file missing: ",trim(flxname)
       print *, " LIS endrun being called ..."
       call LIS_endrun
    endif
  endif

! MERRA-Land ("mld") fields (obs ppt) input file:
  mldname = trim(merralanddir)//'/'//trim(ftime1)//&
       '/'//trim(prefix)//'0.prod.simul.tavg1_2d_mld_Nx.'//&
       trim(ftime2)//'.SUB.nc'

  inquire( file=mldname, exist=file_exists )
  if( .NOT. file_exists ) then
    mldname = trim(merralanddir)//'/'//trim(ftime1)//&
         '/'//trim(prefix)//'1.prod.simul.tavg1_2d_mld_Nx.'//&
         trim(ftime2)//'.SUB.nc'
    inquire( file=mldname, exist=file_exists )
    if( .NOT. file_exists ) then
       print *, " ## ERR MSG:  MLD file missing: ",trim(mldname)
       print *, " LIS endrun being called ..."
       call LIS_endrun
    endif
  endif

! MERRA Land-based ("lnd") fields input file:
  lndname = trim(merralanddir)//'/'//trim(ftime1)//&
       '/'//trim(prefix)//'0.prod.assim.tavg1_2d_lnd_Nx.'//&
       trim(ftime2)//'.SUB.nc'

  inquire( file=lndname, exist=file_exists )
  if( .NOT. file_exists ) then
    lndname = trim(merralanddir)//'/'//trim(ftime1)//&
         '/'//trim(prefix)//'1.prod.assim.tavg1_2d_lnd_Nx.'//&
         trim(ftime2)//'.SUB.nc'
    inquire( file=lndname, exist=file_exists )
    if( .NOT. file_exists ) then
       print *, " ## ERR MSG:  LND file missing: ",trim(lndname)
       print *, " LIS endrun being called ..."
       call LIS_endrun
    endif
  endif

end subroutine merralandfiles

