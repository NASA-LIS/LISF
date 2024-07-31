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
! !ROUTINE: get_RFE2Daily
!  \label{get_RFE2Daily}
!
! !REVISION HISTORY:
!  30 May 2010; Soni Yatheendradas, Initial LDT version for FEWSNET
!  20 Mar 2013; KR Arsenault, Cleaned up code and fixed update time issue
!
! !INTERFACE:
subroutine get_RFE2Daily(n, findex)
! !USES:
  use ESMF
  use LDT_coreMod,           only : LDT_rc
  use LDT_metforcingMod,     only : LDT_forc
  use LDT_timeMgrMod,        only : LDT_calendar, LDT_get_nstep, &
                                    LDT_tick, LDT_date2time, LDT_time2date
  use LDT_logMod,            only : LDT_logunit, LDT_endrun, LDT_verify
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use RFE2Daily_forcingMod,  only : RFE2Daily_struc

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  
!  
! !DESCRIPTION:
!  Opens, reads, and upscales/interpolates daily RFE2.0 forcing.
!
!  At the beginning of a simulation, the code reads the relevant
!  daily data (nearest daily interval).  Since rain value is obtained 
!  from total daily accumulation, only one relevant file needs to be read. 
!  Current first-cut temporal interpolation implemented is
!  even distribution over sub-daily model run timesteps.
!  The strategy for missing data in spatial upscaling is to retain 
!  the corresponding base precip value/s (diurnal precip pattern 
!  then follows base precip and is not an even distribution). Upscaled 
!  values are missing only if ALL the contained RFE2 pixels are NODATA. 
!  The strategy for missing data in spatial interpolation is to set
!  such NODATA values in the RFE2Daily input to 0 (Ronald W Lietzow,
!  USGS, Personal E-mail Communication 05/27/2010). 
!
!  An upper daily 1000 mm cap is also implemented as per information 
!   from USGS staff before being input to the WRSI, MI, and BERM products.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the supplemental forcing source
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    determines the RFE2Daily data times
!  \item[RFE2Dailyfile](\ref{RFE2Dailyfile}) \newline
!    puts together the RFE2Daily file name
!  \item[readprecip\_RFE2Daily](\ref{readprecip_RFE2Daily}) \newline
!    reads the RFE2Daily data 
!  \end{description}
!EOP

!==== Local Variables=======================
  integer :: IsEndTime_RFE2Daily   ! 1=get a new file

! Time parameters for current LDAS time
  integer :: doyNow, yrNow, moNow, daNow, hrNow, mnNow, ssNow, tsNow

! Time parameters for RFE2Daily time start boundary
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1, ts1

! Time parameters for RFE2Daily time end boundary
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2

! Current LDAS time and end boundary time for precip data source
  real*8  :: ctime, EndTime_RFE2Daily

! GMT times for current LDAS time and begin and end boundary times
! for precip data sources
  real    :: gmtNow, gmt1, gmt2

  integer :: ferror_RFE2Daily ! Error flags for precip data sources

  character(len=LDT_CONST_PATH_LEN) :: filename    ! Filename variables for precip data sources

  integer      :: order
!=== End Variable Definition =======================

  IsEndTime_RFE2Daily = 0

!------------------------------------------------------------------------
! Determine required observed precip data times 
! (current, accumulation end, and filename timestamp start time)
! Model current time
!------------------------------------------------------------------------
  yrNow = LDT_rc%yr  !current time
  moNow = LDT_rc%mo
  daNow = LDT_rc%da
  hrNow = LDT_rc%hr
  mnNow = LDT_rc%mn
  ssNow = 0
  tsNow = 0
  call LDT_tick( ctime, doyNow, gmtNow, yrNow, moNow, daNow, hrNow, &
       mnNow, ssNow, real(tsNow))
  
  yr1 = LDT_rc%yr
  mo1 = LDT_rc%mo 
  da1 = LDT_rc%da
  hr1 = LDT_rc%hr
  mn1 = LDT_rc%mn
  ss1 = LDT_rc%ss

!------------------------------------------------------------------------
! Ensure that data is found during first time step
!------------------------------------------------------------------------
  if ( (LDT_get_nstep(LDT_rc,n).eq. 0) .OR. &
       (LDT_get_nstep(LDT_rc,n).eq. 1) .OR. &
       (LDT_rc%rstflag(n) .eq. 1) ) then
     IsEndTime_RFE2Daily = 1
     LDT_rc%rstflag(n) = 0

     ! Required for 1st time step reflecting data representative period
     ! of 0600 - 0600 GMT as per Ronald W Lietzow & James Rowland, USGS,
     ! E-mail Communication 01/07/2011)

     ! if the hour at initialization is less than
     ! the valid hour for the RFE2Daily data (specified by
     ! RFE2Daily_struc%hour_offset) then look back by one day.
     if ( hr1 < RFE2Daily_struc(n)%hour_offset ) then
        call LDT_tick(ctime, doy1, gmt1, yr1, mo1, da1, hr1, &
                      mn1, ss1, -24.0*60*60)
     endif
  endif

!------------------------------------------------------------------------
! Check for and get RFE2Daily Precipitation data
!------------------------------------------------------------------------
  filename=''
  if ( ctime > RFE2Daily_struc(n)%RFE2DailyEndTime ) then
     IsEndTime_RFE2Daily = 1
  endif

  if ( IsEndTime_RFE2Daily == 1 ) then  !get new time2 data

!------------------------------------------------------------------------
!   RFE2Daily product start time info to obtain relevant RFE2Daily file
!------------------------------------------------------------------------

    call RFE2Dailyfile( filename, RFE2Daily_struc(n)%RFE2DailyDir, &
                        yr1, mo1, da1 )

    write(LDT_logunit,*) "Locating RFE2Daily precip file: ",trim(filename)
    ferror_RFE2Daily = 0
    order = 2
    call readprecip_RFE2Daily( n, filename, findex, order, &
                               ferror_RFE2Daily )

!------------------------------------------------------------------------
!   RFE2Daily product end time till which acculumation value applies
!     SY: 0600 next day reflecting data representative
!         period of 0600 - 0600 GMT as per 
!         Ronald W Lietzow & James Rowland, 
!         USGS, E-mail Communication 01/07/2011)
!------------------------------------------------------------------------

    ! Set endtime to tomorrow at hour 0 plus "0600" offset.
    ! Subtract one minute to force reader to trigger on the hour
    ! instead of the timestep after the hour.
    yr2 = yr1
    mo2 = mo1
    da2 = da1
    hr2 = 0
    mn2 = 0
    ss2 = 0
    ts2 = 24*60*60 + RFE2Daily_struc(n)%hour_offset*60*60 - 60
    call LDT_tick( EndTime_RFE2Daily, doy2, gmt2, yr2, mo2, da2, hr2, &
                   mn2, ss2, real(ts2) )

    RFE2Daily_struc(n)%RFE2DailyEndTime = EndTime_RFE2Daily

  endif  !need new time2

end subroutine get_RFE2Daily

!BOP
! !ROUTINE: RFE2Dailyfile
! \label{RFE2Dailyfile}
!
! !INTERFACE:
subroutine RFE2Dailyfile( filename, RFE2DailyDir, yr, mo, da ) 

  implicit none

! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: RFE2DailyDir
  integer           :: yr, mo, da

! !DESCRIPTION:
!   This subroutine puts together the RFE2Daily file name
!
!  The arguments are:
!  \begin{description}
!  \item[RFE2DailyDir]
!    Name of the RFE2Daily directory
!  \item[yr]
!    year 
!  \item[mo]
!    month
!  \item[da]
!    day of month
!  \item[filename]
!    name of the timestamped RFE2Daily file
!  \end{description}
!
!EOP

!==== Local Variables=======================
   character*4  :: fyr
   character*2  :: fmo, fda
!=== End Variable Definition ===============

   write(unit=fyr, fmt='(i4.4)')  yr
   write(unit=fmo, fmt='(i2.2)')  mo
   write(unit=fda, fmt='(i2.2)')  da

  !=== Assemble CPC RFE2.0 filenames:

  ! [cpc_dir]:/YYYY/all_products.bin.YYYYMMDD

   filename = trim(RFE2DailyDir)//"/"//fyr//"/all_products.bin."//&
              fyr//fmo//fda

end subroutine RFE2Dailyfile


