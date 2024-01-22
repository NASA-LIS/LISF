!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_merra2_ac
! \label{get_merra2_ac}
!
!
! !REVISION HISTORY:
! 01 Jun 2022: Michel Bechtold, initial code (based on merra-2 data preprocessed
! to daily data)
! 17 Jan 2024: Louise Busschaert, AC71 implementation in NASA master
!
! !INTERFACE:
subroutine get_merra2_ac(n, findex)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_metforcingMod
  use merra2_ac_forcingMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates daily MERRA2 forcing.
!
!  The MERRA2 Aquacrop forcing data are organized into daily files, where each
!  daily file contains one daily record of forcing fields.  
!
!  In general, metforcing readers read the forcing data before the current
!  time, referred to as bookend1, and after the current time, referred to as
!  bookend2.  Then the readers temporally interpolate between bookend1 and
!  bookend2.  Here, each bookend contains one daily record of forcing,
!  and, in general, the MERRA2 reader will be temporally interpolating
!  from one daily-interval to the next daily-interval, where both daily-intervals
!  are contained in bookend1.  Issues arise between the hours 23z of one day
!  and 1z of the next day, which are complicated by the size of LIS' running
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
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    call to advance or retract time
!  \item[merra2files](\ref{merra2files}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_merra2_ac](\ref{read_merra2_ac}) \newline
!    call to read the MERRA2 data and perform spatial interpolation
!  \end{description}
!EOP
  integer           :: ferror
  character(len=LIS_CONST_PATH_LEN) :: ac71_daily_name
  integer           :: c, r,kk
  integer           :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  integer           :: yr_ac, mo_ac, da_ac, hr_ac, mn_ac, ss_ac, doy_ac
  real*8            :: time1, time2, timenow, time_merra2_ac
  real              :: gmt1, gmt2, gmt_ac
  real              :: ts1, ts2
    
  integer           :: hr_int1, hr_int2
  integer           :: movetime  ! Flag to move bookend2 files to bookend1

! _________________________________________________________

! Please note that the timing logic has been tested only for
! these scenarios:
!
! startime of 2005-11-01T00:30:00 with time-step of 15mn
! startime of 2005-11-01T00:00:00 with time-step of 30mn
! startime of 2005-11-01T00:00:00 with time-step of 1hr
! startime of 2005-11-02T00:00:00 with time-step of 15mn

  merra2_ac_struc(n)%findtime1 = 0
  merra2_ac_struc(n)%findtime2 = 0
  movetime = 0

  ! Initialize ts1 and ts2 timepoints at beginning of run:
  if ( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1 .or. &
       merra2_ac_struc(n)%reset_flag ) then
     merra2_ac_struc(n)%findtime1 = 1
     merra2_ac_struc(n)%findtime2 = 1
     LIS_rc%rstflag(n) = 0
     merra2_ac_struc(n)%reset_flag = .false.

     yr1=LIS_rc%yr
     mo1=LIS_rc%mo
     da1=LIS_rc%da
     hr1=0
     mn1=0
     ss1=0
     !if ( LIS_rc%hr == 0 .and. LIS_rc%mn < 30 ) then
     !   ! initialize ringtime to today at 00:30z
     !   ts1=30*60
     !else
     !   ! initialize ringtime to tomorrow at 00:30z
     !   ts1=86400 + 30*60 ! 1 day plus 30 minutes
     !endif
     ts1=86400 ! 1 day
     call LIS_tick(merra2_ac_struc(n)%ringtime,doy1,gmt1,&
                   yr1,mo1,da1,hr1,mn1,ss1,ts1)
  endif

  !----------------------------------------------------------
  ! Determine current time
  !----------------------------------------------------------
  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=0
  mn1=0
  ss1=0
  call LIS_tick(timenow,doy1,gmt1,&
                yr1,mo1,da1,hr1,mn1,ss1,0.0)

  ! MERRA2_AC daily files --> climate of the last 24 hours before 00:00
   ts1= -86400 ! 1 day 
  call LIS_tick(time_merra2_ac,doy1,gmt1,&
                yr1,mo1,da1,hr1,mn1,ss1,ts1)
  call LIS_time2date(time_merra2_ac,doy_ac,gmt_ac,yr_ac,mo_ac,da_ac,hr_ac,mn_ac)
  
  ! Read MERRA2 - Bookend 1 files:
     do kk= merra2_ac_struc(n)%st_iterid, merra2_ac_struc(n)%en_iterid
        call merra2files_ac(n,merra2_ac_struc(n)%merra2dir, yr_ac, mo_ac, da_ac, ac71_daily_name)
        call read_merra2_ac(n, mo1, &
             findex, ac71_daily_name,&
             merra2_ac_struc(n)%merraforc1(kk,:,:,:), ferror)
     enddo

 ! reset ringtime to tomorrow at 00:30z
 yr1=LIS_rc%yr
 mo1=LIS_rc%mo
 da1=LIS_rc%da
 hr1=0
 mn1=0
 ss1=0
 ts1=86400 ! 1 day 
 call LIS_tick(merra2_ac_struc(n)%ringtime,doy1,gmt1,&
               yr1,mo1,da1,hr1,mn1,ss1,ts1)


     call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

     ! Assign MERRA2 forcing fields to two LIS time-interp placeholders (metdata1,2):
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r).ne.-1) then

                   merra2_ac_struc(n)%metdata1(:,:,LIS_domain(n)%gindex(c,r)) = &
                         merra2_ac_struc(n)%merraforc1(:,:,1,&  ! Store hour: Current hour (same day)
                         (c+(r-1)*LIS_rc%lnc(n)))

            endif
         enddo
      enddo

      ! Assign the hourly times:
      merra2_ac_struc(n)%merra2time1 = time1
      merra2_ac_struc(n)%merra2time2 = time1


end subroutine get_merra2_ac


!BOP
! !ROUTINE: merra2files_ac
! \label{merra2files}
!
! !INTERFACE:
subroutine merra2files_ac(n, merra2dir, yr, mo, da, ac71_daily_name)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_forecastMod
  use LIS_timeMgrMod

  implicit none
! !ARGUMENTS:
  integer                       :: n 
  character(len=*), intent(in)  :: merra2dir
  integer, intent(in)           :: yr,mo,da
  character(len=*), intent(out) :: ac71_daily_name

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
!  \item[ac71_daily_name]
!   name of the timestamped single level file
!  \end{description}
!
!EOP

  character*4  :: cyear
  character*2  :: cmonth
  character*8  :: cdate
  character*10 :: prefix
  integer      :: seed
  real         :: rand
  integer      :: hr, mn, ss
  real*8       :: time
  integer      :: doy
  real         :: gmt

  hr = 0 
  mn = 0 
  ss = 0 

!hack for the synthetic tests
!  call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,&
!       -5*86400.0)

     write(unit=cyear, fmt='(i4.4)') yr
     write(unit=cmonth,fmt='(i2.2)') mo
     write(unit=cdate, fmt='(i4.4,i2.2,i2.2)') yr,mo,da
     
     if (yr<=1979 .and. mo<2) then
        write(LIS_logunit,*) '[ERR] merra2files: date out of range'
        write(LIS_logunit,*) '[ERR] Supported years are from 1979-2-1 through ...'
        call LIS_endrun()
     endif
     prefix= 'MERRA2_AC_'
     
     
     ! Daily fields:
     ac71_daily_name = trim(merra2dir)//'/'//prefix//cdate//'.nc'
     
end subroutine merra2files_ac

