!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_timeMgrMod
!BOP
! !MODULE: LIS_timeMgrMod
!
! !DESCRIPTION:
! 
! This module contains routines for time managment.The module provides
! routines for clock initialization, model timestepping, some basic
! alarm functions, and other useful time managment utilities. 
!
! !USES:
  use ESMF
  use LIS_PRIV_rcMod
  use LIS_logMod
!EOP
  implicit none
  
  PRIVATE 

  integer, parameter :: uninit_int = -999999999 

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LIS_timemgr_init     ! time manager initialization
  PUBLIC :: LIS_timemgr_print    ! print the details of the clock
  PUBLIC :: LIS_advance_timestep ! increment timestep number
  PUBLIC :: LIS_update_timestep  ! updates the LIS timestep 
  PUBLIC :: LIS_update_clock     ! updates the clock with the computed timestep
  PUBLIC :: LIS_timemgr_set      ! set the current time
  PUBLIC :: LIS_get_step_size    ! return step size in seconds
  PUBLIC :: LIS_get_nstep        ! return timestep number
  PUBLIC :: LIS_get_curr_date    ! return date components at end of current timestep
  PUBLIC :: LIS_get_curr_calday  ! return calendar day at end of current timestep
  PUBLIC :: LIS_get_julhr        ! returns time in julian hours
  PUBLIC :: LIS_get_timeoffset_sec        ! returns time offset in seconds
  PUBLIC :: LIS_is_last_step     ! return true on last timestep
  PUBLIC :: LIS_tick             ! clock advance/retract by a certain amount
  PUBLIC :: LIS_date2time        ! converts date to the LIS time format
  PUBLIC :: LIS_time2date        ! converts LIS time to the date format
  PUBLIC :: LIS_registerAlarm
  PUBLIC :: LIS_isAlarmRinging   ! checks if alarm is ringing 
  PUBLIC :: LIS_isDekadalAlarmRinging ! checks to see if the dekadal alarm is ringing
  PUBLIC :: LIS_getDekad         ! returns dekad of month
  PUBLIC :: LIS_getSecsInDekad   ! returns number of seconds in a given dekad
  PUBLIC :: LIS_finishDekadalAlarms ! completes the initialization of dekadal-based alarms
  PUBLIC :: LIS_computeTimeBookEnds ! compute the time interpolation weights 
  PUBLIC :: LIS_computeTemporalWeights !Compute monthly time interpolation weights
  PUBLIC :: LIS_julhr_date       ! converts julian hour to date format
  PUBLIC :: LIS_seconds2time     ! converts time in seconds to hours, minutes and seconds
  PUBLIC :: LIS_resetClock       ! resets the clock back to config file specifications
  PUBLIC :: LIS_tmjul4           ! converts julian time to zulu time. 
  PUBLIC :: LIS_doy2date         ! convert julian day to month and days
  PUBLIC :: LIS_localtime2gmt    ! computes GMT based on localtime
!  PUBLIC :: LIS_gmt2localtime
  PUBLIC :: LIS_parseTimeString
  PUBLIC :: LIS_resetClockForTimeWindow
  PUBLIC :: LIS_mon3char

  PUBLIC :: LIS_twStartTime
  PUBLIC :: LIS_twStopTime
  PUBLIC :: LIS_twMidTime
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: LIS_clock           ! LIS clock object
  PUBLIC :: LIS_calendar        ! LIS calendar object

  type(ESMF_Time),  save     :: LIS_twStartTime
  type(ESMF_Time),  save     :: LIS_twStopTime
  type(ESMF_Time),  save     :: LIS_twMidTime

  type, public :: lisalarmEntry
     character*50 :: name
     real         :: interval
     character*50 :: intervaltype
     real*8       :: alarmTime
     real         :: ts
     integer      :: prev_mo
     integer      :: monthCount
     integer      :: alarm_offset ! seconds to add to alarm time
                                  ! All alarms are based off 0z,
                                  ! so a 3-hourly alarm will ring
                                  ! at 0z, 3z, 6z, etc.
                                  ! To have a 3-hourly alarm ring
                                  ! at 1:30z, 4:30z, 7:30z, etc., set
                                  ! the alarm_offset to 90*60 seconds.

     ! These elements support dekadal alarms
     integer          :: ref_dekad
     integer          :: dek_offset ! seconds to add to to beginning 
                                    ! and ending of dekad periods 
     character(len=5) :: when   ! specifies when the alarm should ring;
                                ! a value of "begin" specifies that the alarm
                                ! should ring at the beginning of the dekad;
                                ! a value of "end" specifies that the alarm
                                ! should ring at the end of the dekad.
     logical           :: firstInstance
     type(lisalarmEntry), pointer :: next     
  end type lisalarmEntry

  type(lisalarmEntry), pointer :: LIS_alarmList

  type(ESMF_Clock), save     :: LIS_clock
  type(ESMF_Calendar), save  :: LIS_calendar  
!BOP
! !ROUTINE: LIS_isAlarmRinging
! \label{LIS_isAlarmRinging}
!
! !INTERFACE: 
   interface LIS_isAlarmRinging
! !PRIVATE MEMBER FUNCTIONS: 
      module procedure isMonthlyAlarmRinging
      module procedure isAlarmRinging

! !DESCRIPTION:  
!
!  checks if the monthly alarm is ringing. The private functions
!  have different arguments based on the input options specified. 
! 
!EOP
   end interface

contains
!BOP
! !ROUTINE: LIS_timemgr_init
! \label{LIS_timemgr_init}
!
! !INTERFACE:
  subroutine LIS_timemgr_init(LIS_rc)
! !USES: 
    use LIS_logMod,   only : LIS_verify

    implicit none
! !ARGUMENTS:
    type(lisrcdec) :: LIS_rc
! !DESCRIPTION:
! 
! Initialize the LIS time manager.
!
! NOTE - This assumes that the LIS time specific variables 
! pertaining to start and end times have been set before 
! this routine is called.  
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!  \end{description}
! 
! The calling sequence is:  
! \begin{description}
!  \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!   to convert current date to a floating point format
!  \item[LIS\_timemgr\_print](\ref{LIS_timemgr_print}) \newline
!   display the contents of the time manager
! \end{description}
!EOP
    type(ESMF_Time)         :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep, tw
    integer                 :: yr, mo
    integer                 :: status
    integer                 :: ts_frac_s, ts_frac_ms

    call LIS_date2time(LIS_rc%etime1,LIS_rc%edoy1,LIS_rc%egmt1, & 
         LIS_rc%eyr1,LIS_rc%emo1,LIS_rc%eda1,LIS_rc%ehr1,LIS_rc%emn1,LIS_rc%ess1)
    LIS_rc%tscount = 0
    LIS_calendar = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,&
         name="Gregorian",& 
         rc=status)

    ts_frac_s = nint(LIS_rc%ts)
    ts_frac_ms = (LIS_rc%ts - nint(LIS_rc%ts))*1000

    call ESMF_TimeIntervalSet(timeStep, s = ts_frac_s, &
         ms = ts_frac_ms, rc=status)
    call LIS_verify(status,'error in ESMF_TimeIntervalSet in LIS_timemgr_init')

    call ESMF_TimeSet(startTime, yy = LIS_rc%syr, &
         mm = LIS_rc%smo, &
         dd = LIS_rc%sda, &
         h  = LIS_rc%shr, &
         m  = LIS_rc%smn, & 
         s  = LIS_rc%sss, &
         calendar = LIS_calendar, & 
         rc = status)
    call LIS_verify(status,'error in ESMF_TimeSet:startTime in LIS_timemgr_init')

    call ESMF_TimeSet(stopTime, yy = LIS_rc%eyr, &
         mm = LIS_rc%emo, &
         dd = LIS_rc%eda, &
         h  = LIS_rc%ehr, &
         m  = LIS_rc%emn, & 
         s  = LIS_rc%ess, &
         calendar = LIS_calendar, & 
         rc = status)
    call LIS_verify(status,'error in ESMF_TimeSet:stopTime in LIS_timemgr_init')

    LIS_clock = ESMF_ClockCreate(name="LIS Clock", timestep=timestep, &
         startTime=startTime, &
         stopTime=stopTime, rc=status)

    call LIS_timemgr_print(LIS_rc)

    if(LIS_rc%twInterval.eq.2592000.0) then 
       
       LIS_twStartTime = startTime
       
       yr = LIS_rc%syr
       mo = LIS_rc%smo+ 1

       if(mo.gt.12) then 
          yr = LIS_rc%syr + 1
          mo = 1
       endif
       
       call ESMF_TimeSet(LIS_twStopTime, yy = yr, &
            mm = mo, &
            dd = 1, &
            h  = 0, &
            m  = 0,& 
            s  = 0, & 
            calendar = LIS_calendar, & 
            rc = status)
       call LIS_verify(status, &
            "ESMF_TimeSet failed in LIS_timemgr_init")    
       
       LIS_twMidTime = startTime
       
       if(LIS_twStopTime.gt.stopTime) then 
          LIS_twStopTime = stopTime
       endif
    else

       call ESMF_TimeIntervalSet(tw,s=nint(LIS_rc%twInterval),rc=status)

       LIS_twStartTime = startTime
       LIS_twMidTime   = startTime 
       LIS_twStopTime  = LIS_twStartTime + tw

       if(LIS_twStopTime.gt.stopTime) then 
          LIS_twStopTime = stopTime
       endif
    endif

    LIS_alarmList => null()
    
    LIS_rc%monthCount = 0 
    LIS_rc%alrm_prev_mo = -1

  end subroutine LIS_timemgr_init

!BOP
! !ROUTINE: LIS_resetClockForTimeWindow
! \label{LIS_resetClockForTimeWindow}
!
! !INTERFACE:
  subroutine LIS_resetClockForTimeWindow(LIS_rc)
! !USES: 
    implicit none
! !ARGUMENTS: 
    type(lisrcdec) :: LIS_rc
! !DESCRIPTION:
!
! resets the time manager clock based on the time window
! (currently assumed to be a month). This routine is 
! intended for the smoother DA run mode where at the end of
! each month the code sets the start time (of the time window)
! to be beginning of the previous month, the middle time to be the 
! beginning of the (current) month and the stop time to be the 
! the beginning of the next month. For e.g. if the routine is
! invoked on July 1st, then the twStartTime will be set to Jun 1, 
! twMiddleTime would be July 1 and twStopTime would be Aug 1.
! The clock is also reset to Jun 1.  
!
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!  \end{description}
!
!  The calling sequence is:  
! \begin{description}
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!    convert date to a flotating point format
! \end{description}
!   
!EOP
    type(ESMF_Time)         :: tTime,stopTime
    type(ESMF_TimeInterval) :: tw
    integer                 :: yr,mo,da,hr,mn,ss,ms
    integer                 :: status
    type(lisalarmEntry), pointer :: alrmEntry, current

    call ESMF_ClockGet(LIS_clock, stopTime=stopTime, rc=status)

    if(LIS_rc%twInterval.eq.2592000.0) then 
       mo = LIS_rc%mo - 1
       yr = LIS_rc%yr
       if(mo.lt.1) then    
          yr = LIS_rc%yr -1
          mo = 12
       endif
       
       ! This always going to back to the start of the previous month!
       
       call ESMF_TimeSet(tTime, yy = yr, &
            mm = mo, &
            dd = 1, &
            h  = 0, &
            m  = 0,& 
            s  = 0, & 
            calendar = LIS_calendar, & 
            rc = status)
       call LIS_verify(status,&
            "ESMF_TimeSet failed in LIS_resetClockForTimeWindow")    
       
       LIS_twStartTime = tTime
       
       write(LIS_logunit,*)'[INFO] Time Window reset month = ',mo
       
       mo = LIS_rc%mo 
       yr = LIS_rc%yr
       
       call ESMF_TimeSet(tTime, yy = yr, &
            mm = mo, &
            dd = 1, &
            h  = 0, &
            m  = 0,& 
            s  = 0, & 
            calendar = LIS_calendar, & 
            rc = status)
       call LIS_verify(status,&
            "ESMF_TimeSet failed in LIS_resetClockForTimeWindow")    
       LIS_twMidTime = tTime
       
       write(LIS_logunit,*)'[INFO] Time Window mid month = ',mo
       
       mo = LIS_rc%mo + 1
       yr = LIS_rc%yr
       if(mo.gt.12) then 
          mo = 1
          yr = LIS_rc%yr+1
       endif
       
       call ESMF_TimeSet(tTime, yy = yr, &
            mm = mo, &
            dd = 1, &
            h  = 0, &
            m  = 0,& 
            s  = 0, & 
            calendar = LIS_calendar, & 
            rc = status)
       call LIS_verify(status,&
            "ESMF_TimeSet failed in LIS_resetClockForTimeWindow")    
       LIS_twStopTime = tTime
       
       write(LIS_logunit,*)'[INFO] Time Window Stop month = ',mo
       
       
       mo = LIS_rc%mo - 1
       yr = LIS_rc%yr
       if(mo.lt.1) then    
          yr = yr -1
          mo = 12
       endif
       da = 1
       hr = 0 
       mn = 0 
       ss = 0 
       ms = 0 

       if(LIS_twStopTime.gt.stopTime) then 
          LIS_twStopTime = stopTime
       endif
       
    else
       
       yr = LIS_rc%yr
       mo = LIS_rc%mo
       da = LIS_rc%da
       hr = LIS_rc%hr
       mn = LIS_rc%mn
       ss = LIS_rc%ss
       
       call ESMF_TimeIntervalSet(tw,s=nint(LIS_rc%twInterval),rc=status)
       
       call ESMF_TimeSet(tTime, yy = yr, &
            mm = mo, &
            dd = da, &
            h  = 0, &
            m  = 0,& 
            s  = 0, & 
            calendar = LIS_calendar, & 
            rc = status)
       tTime = tTime - tw
       
       LIS_twStartTime = tTime       
       LIS_twMidTime = tTime + tw
       LIS_twStoptime = LIS_twMidTime + tw

       if(LIS_twStopTime.gt.stopTime) then 
          LIS_twStopTime = stopTime
       endif

       call ESMF_TimeGet(LIS_twStartTime, yy = yr, &
            mm = mo, &
            dd = da, &
            h  = hr, &
            m  = mn,& 
            s  = ss, & 
            calendar = LIS_calendar, & 
            rc = status)

    endif

    call LIS_timemgr_set(LIS_rc,yr,mo,da,hr,mn,ss,ms,0.0)
    
    LIS_rc%syr= yr
    LIS_rc%smo= mo
    LIS_rc%sda= da
    LIS_rc%shr= hr
    LIS_rc%smn= mn
    LIS_rc%sss= ss
    
    current => LIS_alarmList
    do while(associated(current%next)) 
       current%prev_mo    = -1
       current%monthCount = 0
       current%alarmTime  = 0.0
       if ( current%intervalType == "dekad" ) then
          if ( current%when == "begin" ) then
             current%ref_dekad = LIS_getDekad(LIS_rc,                    &
                                              offset=current%dek_offset, &
                                              which="now")
          else ! current%when == "end"
             current%ref_dekad = LIS_getDekad(LIS_rc,                    &
                                              offset=current%dek_offset, &
                                              which="next")
          endif
       endif
       current => current%next
    enddo

  end subroutine LIS_resetClockForTimeWindow


!BOP
! !ROUTINE: LIS_resetclock
! \label{LIS_resetclock}
!
! !INTERFACE:
  subroutine LIS_resetclock(LIS_rc)
! !USES:
    use LIS_logMod,  only : LIS_verify
    implicit none
! !ARGUMENTS:
    type(lisrcdec) :: LIS_rc
! !DESCRIPTION:
! 
! Resets the LIS time manager.
!
! NOTE - This assumes that the LIS time specific variables 
! pertaining to start and end times have been set before 
! this routine is called.  
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!  \end{description}
! 
! The calling sequence is:  
! \begin{description}
!  \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!   to convert current date to a floating point format
!  \item[LIS\_timemgr\_print](\ref{LIS_timemgr_print}) \newline
!   display the contents of the time manager
! \end{description}
!EOP
    type(ESMF_TimeInterval)      :: timeStep
    type(ESMF_Time)              :: startTime, stopTime
    integer                      :: status
    type(lisalarmEntry), pointer :: current
    integer                      :: ts_frac_s, ts_frac_ms

    call LIS_date2time(LIS_rc%etime1,LIS_rc%edoy1,LIS_rc%egmt1, & 
         LIS_rc%eyr1,LIS_rc%emo1,LIS_rc%eda1,LIS_rc%ehr1,LIS_rc%emn1,LIS_rc%ess1)
    LIS_rc%tscount = 0 
    call LIS_timemgr_print(LIS_rc)
    
    ts_frac_s = nint(LIS_rc%ts)
    ts_frac_ms = (LIS_rc%ts - nint(LIS_rc%ts))*1000

    call ESMF_TimeIntervalSet(timeStep, s = ts_frac_s, &
         ms = ts_frac_ms, & 
         rc=status)
    call LIS_verify(status,'error in ESMF_TimeIntervalSet in LIS_resetClock')

    call ESMF_TimeSet(startTime, yy = LIS_rc%syr, &
         mm = LIS_rc%smo, &
         dd = LIS_rc%sda, &
         h  = LIS_rc%shr, &
         m  = LIS_rc%smn, & 
         s  = LIS_rc%sss, &
         calendar = LIS_calendar, & 
         rc = status)
    call LIS_verify(status,&
         'error in ESMF_TimeSet:startTime in LIS_resetClock')

    call ESMF_TimeSet(stopTime, yy = LIS_rc%eyr, &
         mm = LIS_rc%emo, &
         dd = LIS_rc%eda, &
         h  = LIS_rc%ehr, &
         m  = LIS_rc%emn, & 
         s  = LIS_rc%ess, &
         calendar = LIS_calendar, & 
         rc = status)
    call LIS_verify(status,&
         'error in ESMF_TimeSet:stopTime in LIS_resetClock')

    LIS_clock = ESMF_ClockCreate(name="LIS Clock", timestep=timestep, &
         startTime=startTime, &
         stopTime=stopTime, rc=status)

    LIS_rc%yr= LIS_rc%syr
    LIS_rc%mo= LIS_rc%smo
    LIS_rc%da= LIS_rc%sda
    LIS_rc%hr= LIS_rc%shr
    LIS_rc%mn= LIS_rc%smn
    LIS_rc%ss= LIS_rc%sss
    LIS_rc%ms = 0 
    
    call LIS_date2time(LIS_rc%time,LIS_rc%doy,LIS_rc%gmt, &
         LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn,LIS_rc%ss)
    call LIS_timemgr_print(LIS_rc)
  
    LIS_rc%endtime = 0

    current => LIS_alarmList
    do while(associated(current%next)) 
       current%prev_mo    = -1
       current%monthCount = 0
       current%alarmTime  = 0.0
       if ( current%intervalType == "dekad" ) then
          if ( current%when == "begin" ) then
             current%ref_dekad = LIS_getDekad(LIS_rc,                    &
                                              offset=current%dek_offset, &
                                              which="now")
          else ! current%when == "end"
             current%ref_dekad = LIS_getDekad(LIS_rc,                    &
                                              offset=current%dek_offset, &
                                              which="next")
          endif
       endif
       current => current%next
    enddo
    
  end subroutine LIS_resetclock


!BOP
! !ROUTINE: LIS_timemgr_set
! \label{LIS_timemgr_set}
!
! !INTERFACE:
  subroutine LIS_timemgr_set(LIS_rc,yr,mo,da,hr,mn,ss,ms,ts)
! !USES: 

    implicit none
! !ARGUMENTS: 
    type(lisrcdec) :: LIS_rc
    integer :: yr,mo,da,hr,mn,ss,ms
    real    :: ts
! !DESCRIPTION:
!
! sets the time manager clock based on the input time 
! specification. This method is used to initialize the 
! time mangager when the clock is passed down to LIS from
! a parent component
!
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!   \item[yr]
!     year
!   \item[mn]
!     month
!   \item[da]
!     day of the month
!   \item[hr]
!     hour of the day 
!   \item[mn]
!     minute of the hour
!   \item[ss]
!     second
!  \end{description}
!
!  The calling sequence is:  
! \begin{description}
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!    convert date to a flotating point format
! \end{description}
!   
!EOP
    type(ESMF_Time) :: currTime
    integer         :: status
    real*8          :: time
    real            :: gmt
    integer         :: doy

    LIS_rc%yr = yr
    LIS_rc%mo = mo
    LIS_rc%da = da
    LIS_rc%hr = hr
    LIS_rc%mn = mn
    LIS_rc%ss = ss
    LIS_rc%ms = ms
    
    call LIS_date2time(LIS_rc%time,LIS_rc%doy,LIS_rc%gmt,& 
         LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
    
    call LIS_tick(time,doy,gmt,LIS_rc%yr,LIS_rc%mo,LIS_rc%da, &
        LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,ts)

    LIS_rc%tscount = LIS_rc%tscount + 1
    
    write(unit=LIS_logunit,fmt=24)'[INFO] GSFC-LIS time: ',LIS_rc%mo,'/',LIS_rc%da,'/', & 
         LIS_rc%yr,LIS_rc%hr,':',LIS_rc%mn,':',LIS_rc%ss
24  format(a23,i2.2,a1,i2.2,a1,i4,1x,i2.2,a1,i2.2,a1,i2.2)    ! MLF

    call ESMF_TimeSet(currTime, yy = yr, &
         mm = mo, &
         dd = da, &
         h  = hr, &
         m  = mn, & 
         s  = ss, &
         ms = ms, & 
         calendar = LIS_calendar, & 
         rc = status)
    call LIS_verify(status, &
         'error in ESMF_TimeSet:currTime in LIS_timemgr_set')

    call ESMF_ClockSet(LIS_clock, currTime = currTime, rc=status)
    call LIS_verify(status,&
         'error in ESMF_ClockSet in LIS_timemgr_set')

  end subroutine LIS_timemgr_set
!BOP
! !ROUTINE: LIS_timemgr_print
! \label{LIS_timemgr_print}
!
! !INTERFACE:
  subroutine LIS_timemgr_print(LIS_rc)
!EOP

    implicit none
! !ARGUMENTS: 
    type(lisrcdec) :: LIS_rc
!
! !DESCRIPTION:
! 
! Prints the details of the LIS time manager such as the start, stop
! times, current time, timestep, etc. 
!
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!  \end{description}
!
!EOP 

    write(unit=LIS_logunit,fmt=*)' [INFO] ************************************************'
    
    !write(unit=LIS_logunit,fmt=*)' [INFO] Timestep size (seconds):  ', LIS_rc%ts
    write(unit=LIS_logunit,fmt=*)' [INFO] Start date (ymd tod):     ', LIS_rc%syr, LIS_rc%smo, LIS_rc%sda
    write(unit=LIS_logunit,fmt=*)' [INFO] Stop date (ymd tod):      ', LIS_rc%eyr, LIS_rc%emo, LIS_rc%eda
    write(unit=LIS_logunit,fmt=*)' [INFO] Current date (ymd tod):   ', LIS_rc%yr, LIS_rc%mo, LIS_rc%da
    
    write(unit=LIS_logunit,fmt=*)' [INFO] ************************************************'
  end subroutine LIS_timemgr_print

!BOP
! !ROUTINE: LIS_advance_timestep
! \label{LIS_advance_timestep} 
!
! !INTERFACE:
  subroutine LIS_advance_timestep(LIS_rc)
! !USES: 

    implicit none
! !ARGUMENTS:    
    type(lisrcdec) :: LIS_rc
! !DESCRIPTION:
! 
! Increments time manager clock by the model timestep. In case of LIS 
! running multiple nests, each running at different timesteps, the 
! smallest timestep is chosen to advance the clock. 
!
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!  \end{description}
!
!  The calling sequence is:  
! \begin{description}
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!    convert date to a flotating point format
! \end{description}
!EOP

! Local variables
    type(ESMF_Time) :: currTime
    integer :: days(12),tda,n
    integer :: rc
    real    :: curr_time
    integer :: status
    data days /31,28,31,30,31,30,31,31,30,31,30,31/
    
24  format(a23,i2.2,a1,i2.2,a1,i4,1x,i2.2,a1,i2.2,a1,i2.2)  ! MLF

    call ESMF_ClockAdvance(LIS_clock,rc=status)

    call ESMF_ClockGet(LIS_clock, currTime=currTime, rc=rc)

    call ESMF_TimeGet(currTime,  yy=LIS_rc%yr, &
         mm = LIS_rc%mo, &
         dd = LIS_rc%da, &
         h =  LIS_rc%hr, &
         m =  LIS_rc%mn, &
         s =  LIS_rc%ss,&
         ms = LIS_rc%ms, &
         calendar = LIS_calendar, &
         rc=status)
    call LIS_verify(status, 'error in ESMF_TimeGet: LIS_timeMgrMod')

!    LIS_rc%ss = LIS_rc%ss + LIS_rc%ts
!    
!    do while(LIS_rc%ss .gt. 59) 
!       LIS_rc%ss = LIS_rc%ss - 60 
!       LIS_rc%mn = LIS_rc%mn + 1
!    enddo
!    
!    do while(LIS_rc%mn .gt.59)
!       LIS_rc%mn = LIS_rc%mn -60
!       LIS_rc%hr = LIS_rc%hr+1
!    enddo
!    
!    do while(LIS_rc%hr .ge.24) 
!       LIS_rc%hr = LIS_rc%hr -24
!       LIS_rc%da = LIS_rc%da +1
!    enddo
!  
!    if((mod(LIS_rc%yr,4) .eq. 0 .and. mod(LIS_rc%yr, 100).ne.0) &!leap year
!         .or.(mod(LIS_rc%yr,400) .eq.0)) then 
!       days(2) = 29
!    else 
!       days(2) = 28
!    endif
!    tda = days(LIS_rc%mo)
!    do while(LIS_rc%da.gt.tda)
!       LIS_rc%da = LIS_rc%da - days(LIS_rc%mo)
!       LIS_rc%mo = LIS_rc%mo + 1
!    enddo
!    
!    do while(LIS_rc%mo .gt. 12) 
!       LIS_rc%mo = LIS_rc%mo-12
!       LIS_rc%yr = LIS_rc%yr +1
!    enddo
    
    call LIS_date2time(LIS_rc%time,LIS_rc%doy,LIS_rc%gmt,& 
         LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
        
    curr_time = float(LIS_rc%hr)*3600+60*float(LIS_rc%mn)+float(LIS_rc%ss)
    do n=1,LIS_rc%nnest
       if(mod(curr_time,real(LIS_rc%nts(n))).eq.0) then 
          LIS_rc%tscount(n) = LIS_rc%tscount(n) + 1
       endif
    enddo
    
    if(LIS_rc%endcode.eq.0)then  !end at real-time date (tbd)
       write(*,*)'warning: do not know how to stop in real-time' 
    endif
    if(LIS_rc%endcode.eq.1)then  !end on date specified in lis configuration file
       call LIS_date2time(LIS_rc%etime,LIS_rc%edoy,LIS_rc%egmt, & 
            LIS_rc%eyr,LIS_rc%emo,LIS_rc%eda,LIS_rc%ehr,LIS_rc%emn,LIS_rc%ess)
       if(LIS_rc%time.ge.LIS_rc%etime)then
          LIS_rc%endtime=1
          write(unit=LIS_logunit,fmt=*) '[INFO] LIS cycle completed'
       endif
    endif
    write(unit=LIS_logunit,fmt=24)'[INFO] LIS cycle time: ',LIS_rc%mo,'/',LIS_rc%da,'/', & 
         LIS_rc%yr,LIS_rc%hr,':',LIS_rc%mn,':',LIS_rc%ss

  end subroutine LIS_advance_timestep

!BOP
! !ROUTINE: LIS_update_timestep
! \label{LIS_update_timestep}
! 
! !INTERFACE: 
  subroutine LIS_update_timestep(LIS_rc, n, ts)
    
    implicit none
! !ARGUMENTS: 
    type(lisrcdec) :: LIS_rc
    integer        :: n
    real           :: ts
!EOP

    if(ts.lt.LIS_rc%nts(n)) then 
       LIS_rc%nts(n) = ts
    endif
    LIS_rc%ts = min(LIS_rc%ts, LIS_rc%nts(n))

  end subroutine LIS_update_timestep

!BOP
! !ROUTINE: LIS_update_clock
! \label{LIS_update_clock}
! 
! !INTERFACE: 
  subroutine LIS_update_clock(ts)
! !USES: 
    use LIS_logMod,   only : LIS_verify
    
    implicit none
! !ARGUMENTS: 
    real       :: ts
!EOP

    type(ESMF_TimeInterval) :: timestep
    integer                 :: rc
    integer                 :: ts_frac_s, ts_frac_ms

    ts_frac_s = nint(ts)
    ts_frac_ms = (ts - nint(ts))*1000

    call ESMF_TimeIntervalSet(timeStep, s = ts_frac_s, &
         ms = ts_frac_ms, rc=rc)
    call LIS_verify(rc,'ESMF_TimeIntervalSet failed in LIS_timeMgrMod')

    call ESMF_ClockSet(clock=LIS_clock,timestep=timeStep,rc=rc)
    call LIS_verify(rc,'ESMF_ClockSet failed in LIS_timeMgrMod')
 
  end subroutine LIS_update_clock


!BOP
! !ROUTINE: LIS_get_step_size
! \label{LIS_get_step_size} 
!
! !INTERFACE:
  function LIS_get_step_size(LIS_rc)

    implicit none
! !ARGUMENTS: 
    type(lisrcdec) :: LIS_rc
    integer :: LIS_get_step_size
!
! !DESCRIPTION:
!
! Return the timestep used by the clock in the time manager. 
!
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!   \item[LIS\_get\_step\_size]
!     timestep value
!  \end{description}
!
!EOP

    LIS_get_step_size = LIS_rc%ts

  end function LIS_get_step_size

!BOP
! !ROUTINE: LIS_get_nstep
! \label{LIS_get_nstep} 
! 
! !INTERFACE:
  function LIS_get_nstep(LIS_rc,n)

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
    type(lisrcdec) :: LIS_rc
    integer :: LIS_get_nstep
!
! !DESCRIPTION: 
!
! Return the timestep number for each nest. 
!
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!   \item [n]
!     index of the nest
!   \item[LIS\_get\_nstep]
!     timestep number
!  \end{description}
!
!EOP
    LIS_get_nstep = LIS_rc%tscount(n)
    
  end function LIS_get_nstep

!BOP
! !ROUTINE: LIS_get_curr_day
! \label{LIS_get_curr_day}
!
! !INTERFACE:
  subroutine LIS_get_curr_date(LIS_rc, yr, mon, day, tod, offset)

    implicit none
   
! !ARGUMENTS: 
    type(lisrcdec) :: LIS_rc
    integer, intent(out) ::yr,mon,day,tod
    integer, optional, intent(in) :: offset
! 
! !DESCRIPTION: 
!
! Return date components valid at end of current timestep with an optional
! offset (positive or negative) in seconds.
!
!  The arguments are: 
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!   \item [yr]
!     year
!   \item [mon]
!     month
!   \item [day]
!     day of the month
!   \item[tod]
!     time of day (seconds past 0Z)
!   \item[offset]
!     offset from current time in seconds    
!     positive for future times, negative 
!     for previous times.
!  \end{description}
!EOP
    yr = LIS_rc%yr
    mon = LIS_rc%mo
    day = LIS_rc%da
    tod = LIS_rc%ss+LIS_rc%mn*60+LIS_rc%hr*3600
    
  end subroutine LIS_get_curr_date


!BOP
! !ROUTINE: LIS_get_julhr
! \label{LIS_get_julhr}
! 
! !INTERFACE:
  subroutine LIS_get_julhr(yr,mo,da,hr,mn,ss,julhr)

    implicit none
   
! !ARGUMENTS: 
    integer, intent(in)     :: yr
    integer, intent(in)     :: mo
    integer, intent(in)     :: da
    integer, intent(in)     :: hr
    integer, intent(in)     :: mn
    integer, intent(in)     :: ss
    integer                 :: julhr
!
! !DESCRIPTION: 
!
! Returns the julian hour. In this convention, julian day began at
! midnight at the beginning of May 24, 1968. Interestingly, this
! convention was introduced by NASA for the space program. 
! 
! The arguments are: 
!  \begin{description}
!   \item [yr]
!     year
!   \item [mo]
!     month
!   \item [da]
!     day of the month
!   \item [hr]
!     hour of day
!   \item [mn]
!     minute
!   \item[ss]
!     second 
!   \item [julhr]
!     julian hour
!  \end{description}
!EOP
    integer      :: years
    integer      :: flagly
    integer      :: days
    integer      :: lyears
    integer      :: month(12) 
    integer      :: temphr
    integer      :: thours
    integer      :: hours
    
    data month(1)  /0/
    data month(2)  /31/   
    data month(3)  /59/   
    data month(4)  /90/   
    data month(5)  /120/  
    data month(6)  /151/  
    data month(7)  /181/  
    data month(8)  /212/  
    data month(9)  /243/  
    data month(10) /273/  
    data month(11) /304/  
    data month(12) /334/ 
    
    julhr = 0 

!------------------------------------------------------------------
!   calculate number of 'complete' years past 1968.
!------------------------------------------------------------------
    
    years = yr - 1968 
    
!------------------------------------------------------------------
!        flag if this year is a leap year.
!------------------------------------------------------------------
    
    flagly = mod((yr - 1968), 4)  
         
!------------------------------------------------------------------
!        calculate number of days (only 'complete' years) since 1968.
! ------------------------------------------------------------------
    
    days = years * 365 
    
!------------------------------------------------------------------
!        calculate number of leap years                              
!        (not including 1968, if 68 is the passed yyyy)              
!        add 3 for rounding, since the operator '/' truncates        
!------------------------------------------------------------------
    
    years  = years + 3  
    lyears = years / 4 
    
!------------------------------------------------------------------
!       add one day for every 'complete' leap year to total days.
! ------------------------------------------------------------------
    
    days = days + lyears   
    
!------------------------------------------------------------------
!     add in this year's number of days.
! ------------------------------------------------------------------
    
    days = days + month(mo) + da
    
! ------------------------------------------------------------------
!        calculate the total number of hours                        
!------------------------------------------------------------------
    
    hours = days * 24  
    
!------------------------------------------------------------------
!        calculate today's number of hours. rounded to nearest hour. 
!        note :                                                      
!        since the operator '/' truncates for integer operations, we 
!        add 70 minutes to the passed hh, to facilitate being able to
!        properly round a time.  examples:                           
!           (1430 hrs  will 'round' properly to 1500 hrs if we add 70)
!             1430 + 70 = 1500 / 100 = 15                             
!           (1400 hrs will still 'round' properly to 1400 hrs)        
!             1400 + 70 = 1470 / 100 = 14                             
!------------------------------------------------------------------
    
!   temphr = hr + 70   
!   thours = temphr / 100  
    temphr = hr*3600 + mn*60 + ss
    thours = temphr/3600
    
!------------------------------------------------------------------
!       add today's number of hours to the total number of hours.
!------------------------------------------------------------------
    
    thours = thours + hours
    
!------------------------------------------------------------------
!       if this a leap year and past feb, add 24 hours.
!------------------------------------------------------------------
    
    if ((flagly .eq. 0) .and. (mo .gt. 2)) then
       thours = thours + 24
    endif
    
    julhr = thours

  end subroutine LIS_get_julhr


!BOP
! !ROUTINE: LIS_get_timeoffset_sec
! \label{LIS_get_timeoffset_sec}
! 
! !INTERFACE:
  subroutine LIS_get_timeoffset_sec(yr,mo,da,hr,mn,ss,julss)

    implicit none
   
! !ARGUMENTS: 
    integer, intent(in)     :: yr
    integer, intent(in)     :: mo
    integer, intent(in)     :: da
    integer, intent(in)     :: hr
    integer, intent(in)     :: mn
    integer, intent(in)     :: ss
    integer                 :: julss
!
! !DESCRIPTION: 
!
! Returns the julian hour. In this convention, julian day began at
! midnight at the beginning of May 24, 1968. Interestingly, this
! convention was introduced by NASA for the space program. 
! 
! The arguments are: 
!  \begin{description}
!   \item [yr]
!     year
!   \item [mo]
!     month
!   \item [da]
!     day of the month
!   \item [hr]
!     hour of day
!   \item [mn]
!     minute
!   \item[ss]
!     second 
!   \item [julss]
!     julian hour
!  \end{description}
!EOP
    integer      :: years
    integer      :: flagly
    integer      :: days
    integer      :: lyears
    integer      :: month(12) 
    integer      :: temphr
    integer      :: thours
    integer      :: hours
    
    data month(1)  /0/
    data month(2)  /31/   
    data month(3)  /59/   
    data month(4)  /90/   
    data month(5)  /120/  
    data month(6)  /151/  
    data month(7)  /181/  
    data month(8)  /212/  
    data month(9)  /243/  
    data month(10) /273/  
    data month(11) /304/  
    data month(12) /334/ 
    
    julss = 0 

!------------------------------------------------------------------
!   calculate number of 'complete' years past 1968.
!------------------------------------------------------------------
    
    years = yr - 1968 
    
!------------------------------------------------------------------
!        flag if this year is a leap year.
!------------------------------------------------------------------
    
    flagly = mod((yr - 1968), 4)  
         
!------------------------------------------------------------------
!        calculate number of days (only 'complete' years) since 1968.
! ------------------------------------------------------------------
    
    days = years * 365 
    
!------------------------------------------------------------------
!        calculate number of leap years                              
!        (not including 1968, if 68 is the passed yyyy)              
!        add 3 for rounding, since the operator '/' truncates        
!------------------------------------------------------------------
    
    years  = years + 3  
    lyears = years / 4 
    
!------------------------------------------------------------------
!       add one day for every 'complete' leap year to total days.
! ------------------------------------------------------------------
    
!    days = days + lyears   
!    Note : the routine is really used to compute the time offset
!    within a day. The added use of timespan in years is therefore
!    removed.     
    
!------------------------------------------------------------------
!     add in this year's number of days.
! ------------------------------------------------------------------
    
    days = days + month(mo) + da
    
! ------------------------------------------------------------------
!        calculate the total number of hours                        
!------------------------------------------------------------------
    
    hours = days * 24  
    
!------------------------------------------------------------------
!        calculate today's number of hours. rounded to nearest hour. 
!        note :                                                      
!        since the operator '/' truncates for integer operations, we 
!        add 70 minutes to the passed hh, to facilitate being able to
!        properly round a time.  examples:                           
!           (1430 hrs  will 'round' properly to 1500 hrs if we add 70)
!             1430 + 70 = 1500 / 100 = 15                             
!           (1400 hrs will still 'round' properly to 1400 hrs)        
!             1400 + 70 = 1470 / 100 = 14                             
!------------------------------------------------------------------
    
!   temphr = hr + 70   
!   thours = temphr / 100  
    temphr = hr*3600 + mn*60 + ss
    thours = temphr
    
!------------------------------------------------------------------
!       add today's number of hours to the total number of hours.
!------------------------------------------------------------------
    
    thours = thours + hours*3600
    
!------------------------------------------------------------------
!       if this a leap year and past feb, add 24 hours.
!------------------------------------------------------------------
    
    if ((flagly .eq. 0) .and. (mo .gt. 2)) then
       thours = thours + 24*3600
    endif
    
    julss = thours

  end subroutine LIS_get_timeoffset_sec


!BOP
! !ROUTINE: LIS_get_curr_calday
! \label{LIS_get_curr_calday}
! 
! !INTERFACE:
  function LIS_get_curr_calday(LIS_rc,offset)

   implicit none
   
! !ARGUMENTS: 
   type(lisrcdec) :: LIS_rc
   integer, optional, intent(in) :: offset
   real :: LIS_get_curr_calday
!
! !DESCRIPTION:
!
! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.
! 
! The arguments are: 
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!   \item [offset]
!     offset from current time in seconds, positive for future times, 
!     negative for previous times. 
!   \item [LIS\_get\_curr\_calday]
!     calendar day value
!  \end{description}
!
!EOP
   integer :: days(12)
   integer :: i
   data days /31,28,31,30,31,30,31,31,30,31,30,31/
   
   LIS_get_curr_calday = 0
   if((mod(LIS_rc%yr,4) .eq. 0 .and. mod(LIS_rc%yr, 100).ne.0) &!leap year
        .or.(mod(LIS_rc%yr,400) .eq.0)) then 
      days(2) = 29
   else 
      days(2) = 28
   endif
   if(LIS_rc%mo .ne. 1) then 
      do i=1,LIS_rc%mo-1
         LIS_get_curr_calday = LIS_get_curr_calday+days(i)
      enddo
   endif
   LIS_get_curr_calday = LIS_get_curr_calday + real(LIS_rc%da) + real(LIS_rc%hr)/24 + &
        real(LIS_rc%mn)/(24*60) + real(LIS_rc%ss)/(24*60*60)
   
   if (present(offset)) then
      if (offset > 0) then
         LIS_get_curr_calday = LIS_get_curr_calday + real(offset)/(24*60*60)
      else if (offset < 0) then
         LIS_get_curr_calday = LIS_get_curr_calday - real(offset)/(24*60*60)
      endif
   endif
   
 end function LIS_get_curr_calday


!BOP
! !ROUTINE: LIS_is_last_step
! \label{LIS_is_last_step}
!
! !INTERFACE:
function LIS_is_last_step(LIS_rc)
   implicit none
! !ARGUMENTS: 
   type(lisrcdec) :: LIS_rc
   logical :: LIS_is_last_step
! 
! !DESCRIPTION:
!
! Function returns true on last timestep.
!
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!   \item [LIS\_is\_last\_step]
!     result of the function 
!  \end{description}
!EOP

   LIS_is_last_step = .false.
   if(LIS_rc%time .ge. LIS_rc%etime1) then
      LIS_is_last_step = .true.
   endif
   
 end function LIS_is_last_step

!BOP
!
! !ROUTINE: isMonthlyAlarmRinging
! \label{isMonthlyAlarmRinging}
!
! !INTERFACE:
! !Private name: call using isAlarmRinging
   function isMonthlyAlarmRinging(LIS_rc, alarmTime, intervalType, &
       midmonth)

    implicit none
! !ARGUMENTS: 
    type(lisrcdec), intent(in)   :: LIS_rc
    real*8, intent(inout)        :: alarmTime
    character(len=*), intent(in) :: intervalType
    logical, intent(in)          :: midmonth
    logical                      :: isMonthlyAlarmRinging
!
! !DESCRIPTION:
!  checks if the monthly alarm is ringing. The function returns 
!  true when the elapsed alarm time is greater than the number of 
!  months specified in the alarm's interval. If the midmonth flag
!  is enabled, the elapsed time is counted from the middle of the
!  month to middle of another month, rather than from the 
!  beginning of the month to the end of another month. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[LIS\_rc]
!    instance of the {\tt lis\_module}
!   \item[alarmTime]
!    the elapsed alarm time
!   \item[intervalType]
!    the alarm's frequency
!   \item[midmonth]
!    flag to indicate if the elapsed time is to be counted
!    from the beginning or the middle of the month
!   \item[ringflag]
!    flag indicating the status of the call
!  \end{description}
!   
!  The calling sequence is:  
! \begin{description}
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!    convert date to a flotating point format
! \end{description}
!EOP
    integer :: yr2,mo2,yr
    integer :: numi, zeroi,doy1,doy2
    real    :: gmt2,gmt1
    real*8  :: time2,time
    real*8  :: jan31,apr30,jul31,oct31 ! dates of quarterly intervals
    integer :: janda,janmo          ! january 31 
    integer :: aprda,aprmo          ! april 30
    integer :: julda,julmo          ! july 31
    integer :: octda,octmo          ! october 31 
    real*8  :: quartTime 
    
    zeroi = 0 
    numi = 0 
    isMonthlyAlarmRinging = .false.
    if(intervalType.eq."none" ) then 
       if(alarmTime.ne.0) then 
          isMonthlyAlarmRinging = .false. 
       else 
          isMonthlyAlarmRinging = .true. 
          alarmTime = LIS_rc%udef
       endif
    elseif(intervalType.eq."monthly") then !monthly
       if(midmonth) then 
          if(LIS_rc%da.lt.16) then 
             mo2=LIS_rc%mo
             yr2=LIS_rc%yr
          else
             mo2=LIS_rc%mo+1
             yr2=LIS_rc%yr
             if(mo2.eq.13)then
                mo2=1
                yr2=LIS_rc%yr+1
             endif
          endif
          
          call LIS_date2time(time2,doy2,gmt2,yr2,mo2,&
               numi,zeroi,zeroi,zeroi)          
          if(time2.gt.alarmTime) then 
             isMonthlyAlarmRinging = .true.
             alarmTime = time2
          endif
       else ! check for the end of the month          
          mo2 = LIS_rc%mo + 1
          yr2 = LIS_rc%yr
          if(mo2.eq.13) then 
             mo2 = 1
             yr2 = LIS_rc%yr + 1
          endif
          call LIS_date2time(time2,doy2,gmt2,yr2,mo2,&
               numi,zeroi,zeroi,zeroi)
          if(time2.gt.alarmTime) then 
             isMonthlyAlarmRinging = .true. 
             alarmTime = time2
          endif          
       endif
    elseif(intervalType.eq."quarterly") then !quarterly
       zeroi = 0 
       time=LIS_rc%time
       yr=LIS_rc%yr
       janda=31
       janmo=01
       call LIS_date2time(jan31,doy1,gmt1,yr,janmo,&
            janda,zeroi,zeroi,zeroi)
       aprda=30
       aprmo=04
       call LIS_date2time(apr30,doy1,gmt1,yr,aprmo,&
            aprda,zeroi,zeroi,zeroi)
       julda=31
       julmo=07
       call LIS_date2time(jul31,doy1,gmt1,yr,julmo,&
            julda,zeroi,zeroi,zeroi)
       octda=31
       octmo=10
       call LIS_date2time(oct31,doy1,gmt1,yr,octmo,&
            octda,zeroi,zeroi,zeroi)
       if ( time.ge.jan31 .and. time.lt.apr30 ) then
          quartTime = 1
       elseif ( time.ge.apr30 .and. time.lt.jul31 ) then
          quartTime = 2
       elseif ( time.ge.jul31 .and. time.lt.oct31 ) then
          quartTime = 3
       elseif ( time.ge.oct31 ) then
          quartTime = 4
       elseif ( time.lt.jan31) then
          quartTime = 4
       endif
       
       if(alarmtime .ne. quartTime) then 
          isMonthlyAlarmRinging = .true.
          alarmTime = quartTime
       endif
    endif
  end function isMonthlyAlarmRinging

!BOP
! !ROUTINE: LIS_getDekad
! \label{LIS_getDekad}
!
! !INTERFACE:
function LIS_getDekad(LIS_rc, offset, which)

   implicit none

! !ARGUMENTS: 
   type(lisrcdec), intent(in)             :: LIS_rc
   integer, intent(in), optional          :: offset
   character(len=*), intent(in), optional :: which

! !DESCRIPTION:
!  This function returns the dekad of the month
!
! \begin{verbatim}
!  dekad 1 : 1  <= day <  11
!  dekad 2 : 11 <= day <  21
!  dekad 3 :       day >= 21
! \end{verbatim}
!
!  The arguments are: 
!  \begin{description}
!   \item[LIS\_rc]
!    instance of the {\tt lis\_module}
!   \item[offset]
!    offset in seconds to add to beginning and ending of dekad ranges 
!   \item[which]
!    specifies which to dekad to report.  A value of ``next'' will report \newline
!    the dekad corresponding to the next LIS time-step (LIS\_rc%ts). \newline
!    All other values will report the dekad corresponding \newline
!    to the current time-step.
!  \end{description}
!EOP

   integer :: LIS_getDekad

   type(ESMF_Time) :: refTime

   type(ESMF_Time) :: start_dekad1, end_dekad1, &
                      start_dekad2, end_dekad2, &
                      start_dekad3, end_dekad3

   type(ESMF_TimeInterval) :: offset_interval, &
                              dekad_interval,  &
                              month_interval,  &
                              timeStep

   integer :: yr, mo, da, hr, mn, ss
   integer :: rc

   integer, dimension(12) :: days
   data days /31,28,31,30,31,30,31,31,30,31,30,31/

   yr=LIS_rc%yr
   mo=LIS_rc%mo
   da=LIS_rc%da
   hr=LIS_rc%hr
   mn=LIS_rc%mn
   ss=LIS_rc%ss

!   write(*,'(a,i4,i3,i3,i3,i3,i3)') 'currTime     = ',yr,mo,da,hr,mn,ss

   if ( (mod(yr,4).eq.0.and.mod(yr,100).ne.0) .or. (mod(yr,400).eq.0) ) then
      ! leap year
      days(2) = 29
   else
      days(2) = 28
   endif

   call ESMF_TimeIntervalSet(dekad_interval, s=864000) ! 10 days in seconds
   call ESMF_TimeIntervalSet(month_interval, s=86400*days(mo))

   call ESMF_TimeSet(refTime, yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, rc=rc)

   call ESMF_TimeSet(start_dekad1, yy=yr, mm=mo, dd=1, h=0, m=0, s=0, rc=rc)
   !call ESMF_TimeSet(end_dekad1, yy=yr, mm=mo, dd=10, h=23, m=59, s=59, rc=rc)

   call ESMF_TimeSet(start_dekad2, yy=yr, mm=mo, dd=11, h=0, m=0, s=0, rc=rc)
   !call ESMF_TimeSet(end_dekad2, yy=yr, mm=mo, dd=20, h=23, m=59, s=59, rc=rc)

   call ESMF_TimeSet(start_dekad3, yy=yr, mm=mo, dd=21, h=0, m=0, s=0, rc=rc)

   if ( present(offset) ) then
      call ESMF_TimeIntervalSet(offset_interval, s=offset)
      start_dekad1 = start_dekad1 + offset_interval
      start_dekad2 = start_dekad2 + offset_interval
      start_dekad3 = start_dekad3 + offset_interval
   endif

   end_dekad1 = start_dekad1 + dekad_interval
   end_dekad2 = start_dekad2 + dekad_interval
   end_dekad3 = start_dekad1 + month_interval

   if ( present(which) ) then
      if ( which == "next" ) then
         call ESMF_TimeIntervalSet(timeStep, s=int(LIS_rc%ts))
         refTime = refTime + timeStep
      endif
   endif

   if ( refTime < start_dekad1 ) then
      LIS_getDekad = 3
   elseif ( start_dekad1 <= refTime .and. refTime < end_dekad1 ) then
      LIS_getDekad = 1
   elseif ( start_dekad2 <= refTime .and. refTime < end_dekad2 ) then
      LIS_getDekad = 2
   elseif ( start_dekad3 <= refTime .and. refTime < end_dekad3 ) then
      LIS_getDekad = 3
   elseif ( refTime >= end_dekad3 ) then
      LIS_getDekad = 1
   else
      LIS_getDekad = 0
   endif

#if 0 
   call ESMF_TimeGet(refTime, yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, rc=rc)
   write(*,'(a,i4,i3,i3,i3,i3,i3,i3)') 'refTime      = ',yr,mo,da,hr,mn,ss,LIS_getDekad
   call ESMF_TimeGet(start_dekad1, yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, rc=rc)
   write(*,'(a,i4,i3,i3,i3,i3,i3,i3)') 'start_dekad1 = ',yr,mo,da,hr,mn,ss,LIS_getDekad
   call ESMF_TimeGet(end_dekad1, yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, rc=rc)
   write(*,'(a,i4,i3,i3,i3,i3,i3,i3)') 'end_dekad1   = ',yr,mo,da,hr,mn,ss,LIS_getDekad
   call ESMF_TimeGet(start_dekad2, yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, rc=rc)
   write(*,'(a,i4,i3,i3,i3,i3,i3,i3)') 'start_dekad2 = ',yr,mo,da,hr,mn,ss,LIS_getDekad
   call ESMF_TimeGet(end_dekad2, yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, rc=rc)
   write(*,'(a,i4,i3,i3,i3,i3,i3,i3)') 'end_dekad2   = ',yr,mo,da,hr,mn,ss,LIS_getDekad
   call ESMF_TimeGet(start_dekad3, yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, rc=rc)
   write(*,'(a,i4,i3,i3,i3,i3,i3,i3)') 'start_dekad3 = ',yr,mo,da,hr,mn,ss,LIS_getDekad
   call ESMF_TimeGet(end_dekad3, yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, rc=rc)
   write(*,'(a,i4,i3,i3,i3,i3,i3,i3)') 'end_dekad3   = ',yr,mo,da,hr,mn,ss,LIS_getDekad
#endif

end function LIS_getDekad

!BOP
! !ROUTINE: LIS_isDekadalAlarmRinging
! \label{LIS_isDekadalAlarmRinging}
!
! !INTERFACE:
function LIS_isDekadalAlarmRinging(LIS_rc, current) 

   implicit none

! !ARGUMENTS: 
   type(lisrcdec), intent(in)               :: LIS_rc
   type(lisalarmEntry), pointer, intent(in) :: current

! !DESCRIPTION:
!  This routine checks whether the given dekadal alarm, represented by
!  current, is ringing.
!
!  The routine returns True when ref\_dekad stored in
!  current differs from the dekad returned by LIS\_getDekad.
!
!  The arguments are: 
!  \begin{description}
!   \item[LIS\_rc]
!    instance of the {\tt lis\_module}
!   \item[current]
!    representation of the current LIS dekadal alarm
!  \end{description}
!EOP

   logical :: LIS_isDekadalAlarmRinging

   integer          :: dekad
   character(len=4) :: which

   which = "now"
   if ( current%when == "end" ) then
      which = "next"
   endif

   dekad = LIS_getDekad(LIS_rc, offset=current%dek_offset, which=which)

   LIS_isDekadalAlarmRinging = .false.
   if ( dekad /= 0 ) then
      if ( dekad /= current%ref_dekad ) then
         current%ref_dekad = dekad
         LIS_isDekadalAlarmRinging = .true.
      endif
   endif

end function LIS_isDekadalAlarmRinging

!BOP
! !ROUTINE: LIS_getSecsInDekad
! \label{LIS_getSecsInDekad}
!
! !INTERFACE:
function LIS_getSecsInDekad(LIS_rc, dekad)

   implicit none

! !ARGUMENTS: 
   type(lisrcdec), intent(in) :: LIS_rc
   integer, intent(in) :: dekad

! !DESCRIPTION:
!  This function returns the number of seconds in the given dekad.
!
!  The arguments are: 
!  \begin{description}
!   \item[LIS\_rc]
!    instance of the {\tt lis\_module}
!   \item[dekad]
!    dekad number: 1, 2, or 3
!  \end{description}
!EOP

   integer :: LIS_getSecsInDekad

   type(ESMF_Time)         :: startTime, endTime
   type(ESMF_TimeInterval) :: seconds

   integer, dimension(12) :: days
   integer :: yr, mo, sda, eda

   data days /31,28,31,30,31,30,31,31,30,31,30,31/

   yr = LIS_rc%yr
   mo = LIS_rc%mo

   if ( (mod(yr,4).eq.0.and.mod(yr,100).ne.0) .or. (mod(yr,400).eq.0) ) then
      ! leap year
      days(2) = 29
   else
      days(2) = 28
   endif

   if ( dekad == 1 ) then
      sda = 1
      eda = 10
   elseif ( dekad == 2 ) then
      sda = 11
      eda = 20
   elseif (dekad == 3 ) then
      sda = 21
      eda = days(mo)
      LIS_getSecsInDekad = 86400*(eda-sda+1)
   else
      ! Add error checking
      ! this should not happen
      sda = 1
      eda = 1
   endif

   LIS_getSecsInDekad = 86400*(eda-sda+1)

#if 0
   call ESMF_TimeSet(startTime, yy=yr, mm=mo, dd=sda, h=0, m=0, s=0)
   call ESMF_TimeSet(endTime, yy=yr, mm=mo, dd=eda, h=23, m=59, s=59)

   seconds = endTime - startTime

   call ESMF_TimeIntervalGet(seconds, s=LIS_getSecsInDekad)

   LIS_getSecsInDekad = LIS_getSecsInDekad + 1
#endif

end function LIS_getSecsInDekad

!BOP
! !ROUTINE: LIS_finishDekadalAlarms
! \label{LIS_finishDekadalAlarms}
!
! !INTERFACE:
  subroutine LIS_finishDekadalAlarms(LIS_rc)
! !USES:
    use LIS_logMod,  only : LIS_verify
    implicit none
! !ARGUMENTS:
    type(lisrcdec) :: LIS_rc
! !DESCRIPTION:
! 
! Completes the initialization of dekadal-based alarms.
!
!  \begin{description}
!   \item [LIS\_rc]
!     instance of the {\tt lis\_module}
!  \end{description}
!EOP

    type(lisalarmEntry), pointer :: current

    current => LIS_alarmList
    do while(associated(current%next)) 
       if ( current%intervalType == "dekad" ) then
          if ( current%when == "begin" ) then
             current%ref_dekad = LIS_getDekad(LIS_rc,                    &
                                              offset=current%dek_offset, &
                                              which="now")
          else ! current%when == "end"
             current%ref_dekad = LIS_getDekad(LIS_rc,                    &
                                              offset=current%dek_offset, &
                                              which="next")
          endif
       endif
       current => current%next
    enddo
    
  end subroutine LIS_finishDekadalAlarms

!BOP
!
! !ROUTINE: LIS_computeTimeBookEnds
! \label{LIS_computeTimeBookEnds}
!
! !INTERFACE:
  subroutine LIS_computeTimeBookEnds(LIS_rc,interval,time1,time2)
    implicit none
! !ARGUMENTS: 
    type(lisrcdec), intent(in) :: LIS_rc
    integer, intent(in)   :: interval
    real*8, intent(inout) :: time1,time2
!
! !DESCRIPTION:
!  compute the time end points based on the interval. If the interval
!  is 3 hours, and the current time is 1:30Z, time1 will be set to 
!  1Z and time2 will be set to 4Z. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[LIS\_rc]
!    instance of the {\tt lis\_module}
!   \item[interval]
!    time interval in hours
!   \item[time1]
!    previous time
!   \item[time2]
!    next time
!  \end{description}
!
! The calling sequence is:  
! \begin{description}
!   \item[LIS\_tick] (\ref{LIS_tick}) \newline
!    advance the time by the specified increment
! \end{description}
!EOP
    integer :: yr1,mo1,da1,hr1,mn1,ss1,ts1
    integer :: yr2,mo2,da2,hr2,mn2,ss2,ts2
    real    :: gmt1,gmt2
    integer :: doy1,doy2

    yr1=LIS_rc%yr    !Previous Hour
    mo1=LIS_rc%mo
    da1=LIS_rc%da
    hr1=interval*((LIS_rc%hr)/interval)
    mn1=0
    ss1=0
    ts1=0
    call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,real(ts1))
    
    yr2=LIS_rc%yr    !Next Hour
    mo2=LIS_rc%mo
    da2=LIS_rc%da
    hr2=interval*((LIS_rc%hr)/interval)
    mn2=0
    ss2=0
    ts2=interval*60*60
    call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,real(ts2))

  end subroutine LIS_computeTimeBookEnds


!BOP
!
! !ROUTINE: LIS_computeTemporalWeights
! \label{LIS_computeTemporalWeights}
!
! !INTERFACE:
  subroutine LIS_computeTemporalWeights(LIS_rc,intervalType,time1,time2,&
       wt1,wt2, alarmName)

    implicit none
! !ARGUMENTS: 
    type(lisrcdec)   :: LIS_rc    
    character(len=*) :: intervalType
    integer          :: time1, time2
    real             :: wt1,wt2
    character(len=*), optional :: alarmName
! 
! !DESCRIPTION: 
! 
!  This routine computes the time interpolation weights based on the
!  climatology interval, current time, and the previous and next month. 
!  For example, a monthly climatology data to the current date can 
!  be interpolated as 
!      value  = value1 *wt1+valu2 *wt2
!   where value1 is the value from the previous month's climatology, 
!   and value2 is the value from the next month's climatology and 
!   value is the interpolated value. 
!   
!  The arguments are: 
!  \begin{description}
!   \item[LIS\_rc]
!    instance of the {\tt lis\_module}
!   \item[interval]
!    interval type for the climatology
!   \item[mo1]
!    previous month
!   \item[mo2]
!    next month
!   \item[wt1, wt2]
!    interpolation weights to be used on the climatology data to 
!    interpolate to the current day. 
!   \end{description}
! 
! The calling sequence is:  
! \begin{description}
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!    converts date to a floating point format
! \end{description}
!EOP
    integer   :: juld
    integer :: julm(13)
    real    :: day1
    real    :: day2
    real    :: rday
    integer :: mo1, mo2
    real    :: s1,s2
    type(lisalarmEntry), pointer :: current
    logical                      :: alarm_found
    real*8                       :: t1,t2,ctime
    real                         :: gmt
    integer                      :: rc
    type(ESMF_Time)              :: currTime
    integer                      :: yr,mo,da,hr,mn,ss,doy

    data julm /  0,  31,  59,  90, 120, 151, &
         181 , 212, 243, 273, 304, 334, 365/

    if(intervalType.eq."monthly") then  !monthly
       mo2 = LIS_rc%mo
       if(LIS_rc%da.gt.15) mo2 = mo2+1
       if(mo2.eq.1) mo2 = 13
       mo1 = mo2-1

       juld = julm(LIS_rc%mo)+LIS_rc%da
       if((LIS_rc%mo.eq.3).and.(mod(LIS_rc%yr,4).eq.0)) then 
          juld = juld+1
       endif
!     ------------------------------------------------------------------
!     to simplify time interpolation for days between dec 16 and jan 15,
!     jan 1 to jan 15 are defined as julian days 366 to 380 below,
!     i.e., a wrap around year.
!     ------------------------------------------------------------------
       
       if (juld .le. 15) juld = juld + 365

!     ------------------------------------------------------------------
!     determine the weightings that will be used to interpolate
!     bounding data to current day of year. account for leap years 
!     between feb 15 and march 15 (there are 29 days during this period
!     instead of 28).
!     ------------------------------------------------------------------
  
       if (juld .le. 15) juld = juld + 365
  
!     ------------------------------------------------------------------
!     determine the weightings that will be used to interpolate
!     bounding data to current day of year. account for leap years 
!     between feb 15 and march 15 (there are 29 days during this period
!     instead of 28).
!     ------------------------------------------------------------------
 
       day2 = float(julm(mo2) + 15)
      
       if ((mo2 .eq. 3) .and. (mod(LIS_rc%yr, 4) .eq. 0)) then
          day2 = day2 + 1
       end if
       day1  = float(julm(mo1) + 15)
       rday  = float(juld)
       
       wt1 = (day2 - rday) / (day2 - day1)
       wt2 = (rday - day1) / (day2 - day1)       
       
       if(mo2.eq.13) mo2 = 1
       if(mo1.eq.13) mo1 = 1

!Here everything else except the month is irrelevant
       time1 = mo1
       time2 = mo2

    elseif(intervalType.eq."dekad") then 

       time1 = LIS_getDekad(LIS_rc,offset=0,which="now") + 3*(LIS_rc%mo-1)
       ! There are only 36 dekads in a year.
       if ( time1 == 36 ) then
          time2 = 1
       else
          time2 = time1 + 1
       endif

       wt1 = 0.0
       wt2 = 1.0
    elseif(intervalType.eq."10-day") then 
       call computeDekadTimeIndex(secsFrmBegYr(LIS_rc%doy, &
            LIS_rc%hr, LIS_rc%mn, LIS_rc%ss),&
            LIS_rc%yr, 10, time1, time2, wt1, wt2)
      
    elseif(intervalType.eq."quarterly") then 
       juld = julm(LIS_rc%mo) + LIS_rc%da  

!     ------------------------------------------------------------------
!     case 1: 1 <= day-of-year < 31
!     interpolating from autumn to winter
!
!     This is a continuation of case 5.
!     Here the full from-autumn-to-winter period is
!     304 <= day-of-year < 365+31
!
!     Thus the left bookend is at day 304; the right bookend is at 396.
!     Since the day-of-year has wrapped around the calendar,
!     day-of-year should be mapped to day-of-year+365.
!
!     This gives:
!     s1 = 396 - (day-of-year+365)
!        = 31 - day-of-year
!
!     s2 = (day-of-year+365) - 304
!        = day-of-year + 61
!     ------------------------------------------------------------------

       if (juld .lt. 31) then
          
          s1   = float(31 - juld)
          s2   = float(juld + 61)
          wt1 = s1 / (396 - 304)
          wt2 = s2 / (396 - 304)
          mo1 = 4
          mo2 = 1

!     ------------------------------------------------------------------
!     case 2: 31 <= day-of-year < 120
!     interpolating from winter to spring
!     ------------------------------------------------------------------

       else if ( (juld .lt. 120) .and. (juld .ge. 31)) then
          
          s1   = float(120 - juld)
          s2   = float(juld - 31)
          wt1 = s1 / (120 - 31)
          wt2 = s2 / (120 - 31)
          mo1 = 1
          mo2 = 2
          
!     ------------------------------------------------------------------
!     case 3: 120 <= day-of-year < 212
!     interpolating from spring to summer
!     ------------------------------------------------------------------

       else if ( (juld .lt. 212) .and. (juld .ge. 120) ) then
          
          s1   = float(212 - juld)
          s2   = float(juld - 120)
          wt1 = s1 / (212 - 120)
          wt2 = s2 / (212 - 120)
          mo1 = 2
          mo2 = 3

!     ------------------------------------------------------------------
!     case 4: 212 <= day-of-year < 304
!     interpolating from summer to autumn
!     ------------------------------------------------------------------

       else if ( (juld .lt. 304) .and. (juld .ge. 212) ) then
          
          s1   =  float(304 - juld)
          s2   =  float(juld - 212)
          wt1 =  s1 / (304 - 212)
          wt2 =  s2 / (304 - 212)
          mo1 = 3
          mo2 = 4
          
!     ------------------------------------------------------------------
!     case 5: 304 <= day-of-year
!     interpolating from autumn to winter
!
!     Here the full from-autumn-to-winter period is
!     304 <= day-of-year < 365+31
!
!     Thus the left bookend is at day 304; the right bookend is at 396.
!     ------------------------------------------------------------------

       else
          
          s1 = float(396 - juld)
          s2 = float(juld - 304)
          wt1 = s1 / (396 - 304)
          wt2 = s2 / (396 - 304)
          mo1 = 4
          mo2 = 1          

       end if

       !Here everything else except the month entry is irrelevant,
       !where the month entry actually represents the quarter.
       time1 = mo1
       time2 = mo2

    elseif(intervalType.eq."weekly") then 
       print*, 'weekly alarm not supported'
       stop
    endif
  end subroutine LIS_computeTemporalWeights

  subroutine computeDekadTimeIndex(secsFrmBegYr, yr, &
       timeSchema_arg, time1, time2, wt1, wt2)

   integer*4, intent(in) :: secsFrmBegYr
   integer, intent(in) :: yr
   integer*4, intent(in) :: timeSchema_arg
   integer*4             :: time1
   integer*4             :: time2
   real                  :: wt1, wt2
    
   integer :: time1_val
   integer :: time2_val
   logical :: checkForDekadInterval
   logical :: firstInstance
!hardcoded
   integer*4, parameter :: gDEKAD_TIMESTEPS_PER_YEAR = 36
   integer*4, parameter :: gDEKAD_TIME_SCHEMA = 10
   
   integer*4 :: timeSchema 
   integer*4 :: ts
   integer*4, save :: dekOy = (0-1)
!!! Dekad time schema variables begin !!!
! Remove the guess work re. ts start hour:  We should by default always be using offset 0600 per corresp. Jul 2011
#ifdef NO_OFFSET_BY_0600
#undef NO_OFFSET_BY_0600
#endif
   ! Dekad time step seconds from beginning of the year non- leap year
   integer*4, parameter, dimension(gDEKAD_TIMESTEPS_PER_YEAR) :: dekTSsfbyNLY = (/&
#ifndef NO_OFFSET_BY_0600
        885600  ,&
        1749600 ,&
        2700000 ,&
        3564000 ,&
        4428000 ,&
        5119200 ,&
        5983200 ,&
        6847200 ,&
        7797600 ,&
        8661600 ,&
        9525600 ,&
        10389600,&
        11253600,&
        12117600,&
        13068000,&
        13932000,&
        14796000,&
        15660000,&
        16524000,&
        17388000,&
        18338400,&
        19202400,&
        20066400,&
        21016800,&
        21880800,&
        22744800,&
        23608800,&
        24472800,&
        25336800,&
        26287200,&
        27151200,&
        28015200,&
        28879200,&
        29743200,&
        30607200,&
        31557600 &
#else
        864000  ,&
        1728000 ,&
        2678400 ,&
        3542400 ,&
        4406400 ,&
        5097600 ,&
        5961600 ,&
        6825600 ,&
        7776000 ,&
        8640000 ,&
        9504000 ,&
        10368000,&
        11232000,&
        12096000,&
        13046400,&
        13910400,&
        14774400,&
        15638400,&
        16502400,&
        17366400,&
        18316800,&
        19180800,&
        20044800,&
        20995200,&
        21859200,&
        22723200,&
        23587200,&
        24451200,&
        25315200,&
        26265600,&
        27129600,&
        27993600,&
        28857600,&
        29721600,&
        30585600,&
        31536000 &
#endif
        /)
   ! Dekad time step seconds from beginning of the year leap year
   integer*4, parameter, dimension(gDEKAD_TIMESTEPS_PER_YEAR) :: dekTSsfbyLY = (/&
#ifndef NO_OFFSET_BY_0600
        885600  ,&
        1749600 ,&
        2700000 ,&
        3564000 ,&
        4428000 ,&
        5205600 ,&
        6069600 ,&
        6933600 ,&
        7884000 ,&
        8748000 ,&
        9612000 ,&
        10476000,&
        11340000,&
        12204000,&
        13154400,&
        14018400,&
        14882400,&
        15746400,&
        16610400,&
        17474400,&
        18424800,&
        19288800,&
        20152800,&
        21103200,&
        21967200,&
        22831200,&
        23695200,&
        24559200,&
        25423200,&
        26373600,&
        27237600,&
        28101600,&
        28965600,&
        29829600,&
        30693600,&
        31644000 &
#else
        864000  ,&
        1728000 ,&
        2678400 ,&
        3542400 ,&
        4406400 ,&
        5184000 ,&
        6048000 ,&
        6912000 ,&
        7862400 ,&
        8726400 ,&
        9590400 ,&
        10454400,&
        11318400,&
        12182400,&
        13132800,&
        13996800,&
        14860800,&
        15724800,&
        16588800,&
        17452800,&
        18403200,&
        19267200,&
        20131200,&
        21081600,&
        21945600,&
        22809600,&
        23673600,&
        24537600,&
        25401600,&
        26352000,&
        27216000,&
        28080000,&
        28944000,&
        29808000,&
        30672000,&
        31622400 &
#endif
        /)
!!! Dekad time schema variables end !!!
   
   
   checkForDekadInterval = .false.
   timeSchema = gDEKAD_TIME_SCHEMA ! assume dekad-time-schema* unless specified otherwise
   timeSchema = timeSchema_arg
   
   select case(timeSchema)
   case (gDEKAD_TIME_SCHEMA)

! From Tamuka Magadzire e-mail dated February 16, 2011:
!> For each month:
!> Dekad 1: the first 10 days of the month Dekad 2: the 11th to the 20th 
!> of the month Dekad 3: whatever days are left in the month after the 
!> 20th, whether its 8 days, 9, 10 or 11 days.

! Accordingly, see 'timestepTablesDekads2011.xls' (Excel 97-2003 workbook) 
! which was created to give the above table/s.
      
      if((mod(yr,4) /= 0) .or. ((mod(yr,100) == 0) &
           .and. (mod(yr,400) /= 0) )  ) then
         ! non- leap year
         
         checkForDekadInterval = .true. 
!find the closest dekad. 
         do ts=1, gDEKAD_TIMESTEPS_PER_YEAR, 1
            if( secsFrmBegYr .lt.dekTSsfbyNLY(ts) ) then
               time2 = ts
               time2_val = dekTSsfbyNLY(ts)
               time1 = max(1,time2 -1)
               time1_val = dekTSsfbyNLY(time1)
               exit;
            endif
         end do
      else
         ! is leap year
         do ts=1, gDEKAD_TIMESTEPS_PER_YEAR, 1
            if( secsFrmBegYr .lt.dekTSsfbyLY(ts) ) then
               time2 = ts
               time2_val = dekTSsfbyLY(ts)
               time1 = max(1,time2 -1)
               time1_val = dekTSsfbyLY(time1)
               exit;
            endif
         end do
      endif
! Seconds from beginning of year is used instead of other temporal representations because:
!   1) this was given the most thought and therefore closest to completion at time of this writing
!   2) this representation represents the minimal requirements for any code that may use a time loop*
!   3) this exacting approach should require less logic in LIS for discrete time steps
!  *including both LIS and standalone model codes

   case default
      ! unhandled time schema
   end select

!for linear weighting   
!   wt1 = (real(time2_val) - real(secsFrmBegYr))/&
!        (real(time2_val) - real(time1_val))
!   wt2 = (real(secsFrmBegYr)-real(time1_val))/&
!        (real(time2_val) - real(time1_val))
!constant weighting
   wt1 = 0
   wt2 = 1.0
 end subroutine computeDekadTimeIndex

!BOP
!
! !ROUTINE: LIS_date2time
! \label{LIS_date2time} 
! 
! !INTERFACE:
  subroutine LIS_date2time(time,doy,gmt,yr,mo,da,hr,mn,ss)

    implicit none
! !ARGUMENTS:
    integer,intent(in)   :: yr,mo,da,hr,mn,ss
    real*8, intent(out)  :: time
    integer,intent(out)  :: doy
    real,   intent(out)  :: gmt
!
! !DESCRIPTION:
! 
!  Determines the time, time in GMT, and the day of the year
!  based on the value of year, month, day of month, hour of 
!  the day, minute and second. This method is the inverse of 
!  time2date.
!
!   NOTE: This routine has been known to give round off error
!   problems when attempting to retrieving minutes and seconds
!   from the given number. Use at your own risk!
! 
!  The arguments are: 
!  \begin{description}
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of the month
!  \item[hr]
!   hour of day
!  \item[mn]
!   minute
!  \item[ss]
!   second
!  \item[time]
!   lis time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
    integer :: yrdays,days(13),k
    data days /31,28,31,30,31,30,31,31,30,31,30,31,30/

    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
         .or.(mod(yr,400).eq.0))then             !correct for y2k
       yrdays=366                  
    else
       yrdays=365
    endif
    
    doy=0
    do k=1,(mo-1)
       doy=doy+days(k)
    enddo
    doy=doy+da
    
    if(yrdays.eq.366.and.mo.gt.2)doy=doy+1
    
    time=(float(yr)+((((((float(ss)/60.d0)+float(mn))/60.d0)+ & 
         float(hr))/24.d0)+float(doy-1))/float(yrdays))
    
    gmt=( ( (float(ss)/60.0) +float(mn)) /60.0)+float(hr)
    return
  end subroutine LIS_date2time

!BOP
! 
! !ROUTINE: LIS_time2date
! \label{LIS_time2date}
! 
! !INTERFACE: 
  subroutine LIS_time2date(time,doy,gmt,yr,mo,da,hr,mn)

    implicit none
! !ARGUMENTS: 
    real*8  :: time
    integer :: yr,mo,da,hr,mn,ss,doy
    real    :: gmt
!
! !DESCRIPTION:
! 
!  Determines the value of the year, month, day of month, hour of 
!  the day, minute and second based on the specified time. This
!  method is the inverse of date2time.
! 
!  The arguments are: 
!  \begin{description}
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of the month
!  \item[hr]
!   hour of day
!  \item[mn]
!   minute
!  \item[ss]
!   second
!  \item[time]
!   lis time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
    real*8 :: tmp
    integer :: yrdays,days(13)
    data days /31,28,31,30,31,30,31,31,30,31,30,31,30/
    
    yr  = aint(time)
    tmp =     (time) 
    
    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
         .or.(mod(yr,400).eq.0))then               !correct for y2k
       yrdays=366                  
    else
       yrdays=365
    endif
    if (yrdays.eq.366) then
       days(2)=29
    else
       days(2)=28
    endif
    
    doy  = aint((tmp-yr)*float(yrdays))+1 
    tmp =      ((tmp-yr)*float(yrdays))+1 
    hr  = nint((tmp-doy)*24.0) 
    tmp =     ((tmp-doy)*24.0) 

    mn  = aint((tmp-hr)*60.0) 
    tmp =     ((tmp-hr)*60.0) 

    ss  = aint((tmp-mn)*60.0) 
    mo=1
    do while (doy.gt.0)
       doy=doy-days(mo)
       mo=mo+1
    enddo
    mo=mo-1
    da=doy+days(mo)
    
    gmt=(((float(ss)/60.0)+float(mn))/60.0)+float(hr)
    
    if(gmt.eq.24) then
       hr = 0 
       gmt=0
       da=da+1
       if (da.gt.days(mo)) then
          da=1
          mo=mo+1
          if (mo.gt.12) then
             mo=1
             yr=yr+1
          endif
       endif
    endif
    return
  end subroutine LIS_time2date
 
!BOP
! !ROUTINE: LIS_doy2date
! \label{LIS_doy2date}
! 
! !INTERFACE:
  subroutine LIS_doy2date(yr, doy, mo, da)

    implicit none
!
! !DESCRIPTION: 
!  This subroutine converts a given julian day to month and days
!EOP
    integer        :: yr
    integer        :: doy
    integer        :: mo
    integer        :: da

    integer :: dsum, i
    integer :: days(13), yrdays
    data days /31,28,31,30,31,30,31,31,30,31,30,31,30/

    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &    !correct for leap year
         .or.(mod(yr,400).eq.0))then               !correct for y2k
       yrdays=366                  
    else
       yrdays=365
    endif
    if (yrdays.eq.366) then
       days(2)=29
    else
       days(2)=28
    endif
    
    dsum = 0 
    mo = 1
    do while(dsum < doy) 
       dsum = dsum+days(mo)
       mo = mo+1
    enddo
    mo = mo-1
    dsum= 0
    do i=1,mo-1
       dsum = dsum + days(i)
    enddo
    da = doy-dsum

  end subroutine LIS_doy2date
    
! 
! !ROUTINE: LIS_seconds2time
! \label{LIS_seconds2time}
! 
! !INTERFACE: 
  subroutine LIS_seconds2time(secs,da, hr, mn, ss)
! !ARGUMENTS:     
    integer,     intent(in)   :: secs
    integer,     intent(out)  :: da
    integer,     intent(out)  :: hr
    integer,     intent(out)  :: mn
    integer,     intent(out)  :: ss
!
! !DESCRIPTION: 
!   Converts time in seconds to hours, minutes and seconds.
!EOP    

    integer                   :: tmp

    tmp = secs

    da = floor(float(tmp)/86400.0)
    if(da.gt.0) then 
       tmp = tmp - da*86400
    endif
    hr = floor(float(tmp)/3600.0)
    if(hr.gt.0) then 
       tmp = tmp - hr*3600 
    endif
    mn = floor(float(tmp)/60.0)
    if(mn.gt.0) then 
       tmp = tmp - mn*60
    endif
    ss = tmp
       
  end subroutine LIS_seconds2time

!BOP
! 
! !ROUTINE: LIS_tick
! \label{LIS_tick}
! 
! !INTERFACE:
  subroutine LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,ts)
    implicit none
! !ARGUMENTS:
    real*8  :: time
    integer :: yr,mo,da,hr,mn,ss,doy
    real    :: ts,gmt
! !DESCRIPTION:
! 
!  Method to advance or retract the time by the specified time
!  increment
! 
!  The arguments are: 
!  \begin{description}
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of the month
!  \item[hr]
!   hour of day
!  \item[mn]
!   minute
!  \item[ss]
!   second
!  \item[ts]
!   time increment in seconds
!  \item[time]
!   lis time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
    type(ESMF_Time)         :: time1
    type(ESMF_TimeInterval) :: timeStep
    integer days(13)
    integer prvmo   !previous month
    integer ts_frac_s, ts_frac_ms
    integer status
    
    data days/31,28,31,30,31,30,31,31,30,31,30,31,31/

    call ESMF_TimeSet(time1,  yy=yr, &
         mm = mo, &
         dd = da, &
         h =  hr, &
         m =  mn, &
         s =  ss,&
         calendar = LIS_calendar, &
         rc=status)
    call LIS_verify(status, 'error in ESMF_TimeSet: LIS_timeMgrMod')

    ts_frac_s = nint(ts)
    ts_frac_ms = (ts - nint(ts))*1000

    call ESMF_TimeIntervalSet(timeStep, s = ts_frac_s, &
         ms = ts_frac_ms, rc=status)
    call LIS_verify(status,'error in ESMF_TimeIntervalSet in LIS_timemgrMod')

    time1 = time1 +timeStep
    call ESMF_TimeGet(time1,  yy=yr, &
         mm = mo, &
         dd = da, &
         h =  hr, &
         m =  mn, &
         s =  ss,&
         calendar = LIS_calendar, &
         rc=status)
    call LIS_verify(status, 'error in ESMF_TimeGet: LIS_timeMgrMod')

    call LIS_date2time(time,doy,gmt,yr,mo,da,hr,mn,ss)
    return

  end subroutine LIS_tick


!BOP
!
! !ROUTINE: LIS_julhr_date
! \label{LIS_julhr_date}
!
! !REVISION HISTORY:
!     15 oct 1998  initial version.........................mr moore/dnxm
!     10 aug 1999  ported to ibm sp2.  added intent attributes to
!                  arguments...............................mr gayno/dnxm
!     29 oct 2005  Sujay Kumar, Adopted in LIS
!
! !INTERFACE:    
subroutine LIS_julhr_date( julhr, yyyy,mm,dd,hh)
! !USES:
  use LIS_logMod, only : lis_abort

  implicit none 
! !ARGUMENTS: 
  integer,          intent(in)   :: julhr   
  integer                        :: yyyy    
  integer                        :: mm       
  integer                        :: dd     
  integer                        :: hh     

!
! !DESCRIPTION:
!     to convert from a julian hour to a 10 digit date/time group
!     (yyyymmddhh)
!  
!     \subsubsection{Method}: 
!     call utility routine LIS\_tmjul4 to convert from julian hours
!     to year, month, day, hour. \newline
!     perform several checks to ensure date information passed
!     back from LIS\_tmjul4 is valid. \newline
!     if date is good, convert from integer data to a 10 digit
!     date/time group and pass back to calling routine. \newline
!
!   The arguments are: 
!   \begin{description}
!    \item[julhr]
!     input julian hour
!    \item[yyyy]
!     four digit year
!    \item[mm]
!     month of the year
!    \item[dd]
!     day of the month
!    \item[hh]
!     time of the day in hours
!    \end{description}
!
! The calling sequence is:  
! \begin{description}
!   \item[LIS\_tmjul4] (\ref{LIS_tmjul4}) \newline
!    convert julian hour to hour, day, month, and year
!   \item[LIS\_abort] (\ref{LIS_abort}) \newline
!    abort the code in case of error
! \end{description}
!
!EOP
  character*255                  :: message ( 20 )  

!     ------------------------------------------------------------------
!     executable code begins here... use LIS_tmjul4 to convert julhr to 
!     hour, day, month and year
!     ------------------------------------------------------------------

  call LIS_tmjul4( hh, dd, mm, yyyy, julhr )
 
!     ------------------------------------------------------------------
!     check for valid hour, day, month, and year
!     ------------------------------------------------------------------ 

  if( (  hh .lt.    0 .or.   hh .gt.   23) .or. &
       (  dd .lt.    1 .or.   dd .gt.   31) .or. &
       (  mm .lt.    1 .or.   mm .gt.   12) .or. &
       (yyyy .lt. 1968 .or. yyyy .gt. 9999) )then
     
     message(1) = 'program: unknown'
     message(2) = '  routine julhr_date10'
     message(3) = '  invalid julhr to date/time conversion'
     
     call lis_abort( message )
     return
  endif

!     ------------------------------------------------------------------ 
!     output errors:
!     flag invalid day/month combinations, for example...
!     if february check for correct number of days in month
!     ------------------------------------------------------------------ 
    
  if( mm .eq. 2 )then
     
     if( (mod((yyyy - 68), 4) .ne. 0) .and. (dd .gt. 28) )then
        
        message(1) = 'program: unknown'
        message(2) = '  routine julhr_date10'
        message(3) = '  oops, created 29 feb in non-leap year'
        
        call lis_abort( message )
     endif
     
  elseif (  (dd .gt. 30) .and. &
       ( (mm .eq. 4)  .or. &
       (mm .eq. 6)  .or. &
       (mm .eq. 9)  .or. &
       (mm .eq. 11) ) ) then  
     
     message(1) = 'program: unknown'
     message(2) = '  routine julhr_date10'
     message(3) = '  oops, created 31st day for month'
     message(4) = '  with only 30 days'
     
     call lis_abort( message )
  endif
  
  return
end subroutine LIS_julhr_date

!BOP
!
! !ROUTINE: LIS_tmjul4
! \label{LIS_tmjul4}
!
! !REVISION HISTORY:
!     18 mar 98 initial version.........................sra milburn/dnxm
!     08 jul 99 fixed error which initialzed done flag in a data
!               data statement.  variables must be initialized
!               using an assignment statement.  ported to ibm sp2.......
!               ...........................................mr gayno/dnxm
!     29 oct 2005  Sujay Kumar, Adopted in LIS
!
! !INTERFACE:    
subroutine LIS_tmjul4( hour, day, month, year, julhr ) 

  implicit none 
! !ARGUMENTS:   
  integer,  intent(out)       :: day  
  integer,  intent(out)       :: hour 
  integer,  intent(in)        :: julhr
  integer,  intent(out)       :: month
  integer,  intent(out)       :: year 
!
! !DESCRIPTION:
!     uses the julian hour to determine the 2-digit zulu time, day,  
!     month, and 4-digit year.    
!
!     
!    \subsubsection{Method}
!     - determine the current zulu hour   \newline
!     - determine the total number of elapsed days   \newline
!     - count forward to the current day/month/year   \newline
!
!    The arguments are: 
!    \begin{description}
!    \item[hour]
!     the zulu time of the julian hour   
!    \item[day]
!     day of the month (1..31)   
!    \item[month]
!     month of the year (1..12)  
!    \item[year]
!     four digit year
!    \item[julhr]
!     the julian hour being processed
!    \end{description}
!
!EOP  
  integer                     :: dypmon(12)   
  integer                     :: elapdy   
  integer                     :: i
  logical                     :: done 
  
  data dypmon   /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  
!     ------------------------------------------------------------------
!     initialize done flag to false.
!     ------------------------------------------------------------------

  done = .false.
      
!    ------------------------------------------------------------------
!     extract the zulu hour from the julian hour.   
!     ------------------------------------------------------------------
          
  hour = mod(julhr, 24) 

!    ------------------------------------------------------------------
!     determine the number of days that have elapsed since dec 31, 1967 
!     (julian hour 0).  
!     ------------------------------------------------------------------
    
  elapdy = julhr / 24   
          
!     ------------------------------------------------------------------
!     initialize the starting day, month, and year values.  
!     ------------------------------------------------------------------
    
  if (elapdy .gt. 0) then        
     year   = 1968   
     day    = 1  
     month  = 1  
     elapdy = elapdy - 1      
  else       
     year   = 1967   
     day    = 31 
     month  = 12 
  endif
      
!    ------------------------------------------------------------------
!     loop through the elapsed days to determine the year.  
!     ------------------------------------------------------------------
    
  do while(.not.done) 

     dypmon(2) = 28  

!    ------------------------------------------------------------------
!     determine if year value is a leap year.  leap years occur every 
!     4 years, with the exception of century years not evenly 
!     divisible by 400.   
!     ------------------------------------------------------------------
    
     if ((mod(year, 4)   .eq. 0) .and.   &
          ((mod(year, 100) .ne. 0) .or. &
          (mod(year, 400) .eq. 0))) then  
        
        dypmon(2) = 29
        
     endif

!     ------------------------------------------------------------------
!     if the elapsed number of days is more than a year's worth,  
!     subtract the appropriate number of days, and increment the year 
!     value.  
!     ------------------------------------------------------------------
    
     if (dypmon(2) .eq. 28) then      
        if (elapdy .ge. 365) then         
           year = year + 1 
           elapdy = elapdy - 365           
        else          
           done = .true.   
        endif
     else     
        if (elapdy .ge. 366) then         
           year = year + 1 
           elapdy = elapdy - 366           
        else          
           done = .true.           
        endif
     endif

!     ------------------------------------------------------------------
!     if the elapsed number of days is less than a year's worth, then   
!     exit loop.
!     ------------------------------------------------------------------    
  enddo

!     ------------------------------------------------------------------
!     count the days and months elapsed in the current year.
!     ------------------------------------------------------------------
    
  do i = 1, elapdy      
     day = day + 1        
     if (day .gt. dypmon(month)) then        
        day = 1   
        month = month + 1        
     endif     
  enddo
  return
  
end subroutine LIS_tmjul4
#if 0 
!BOP
! 
! !ROUTINE: LIS_gmt2localtime
!  \label{LIS_gmt2localtime}
! 
! !DESCRIPTION: 
! 
! Calculates GMT based on the local time 
! 
! !INTERFACE:
subroutine LIS_gmt2localtime (gmt,lon,lhour,zone)
! !ARGUMENTS:
  real:: gmt                 ! GMT time (0-23)
  real:: lon                 ! longitude in degrees
  real:: change              ! the change in number of hours between
  real:: lhour               ! local hour (0-23) 0= midnight, 23= 11:00 p.m.
       
  integer::  i              ! working integer
  integer:: zone            ! time zone (1-24)
!EOP
!----------------------------------------------------------------------
!  Determine into which time ZONE (15 degree interval) the
!  longitude falls.
!----------------------------------------------------------------------
  do i=1,25
     if (lon.lt.(-187.5+(15*i))) then
        zone=i
        if (zone.eq.25) zone=1
        go to 60
     end if
  end do
!----------------------------------------------------------------------
!     Calculate change (in number of hours) from GMT time to
!     local hour.  Change will be negative for zones < 13 and
!     positive for zones > 13.
!     There is also a correction for LHOUR < 0 and LHOUR > 23
!     to LHOUR between 0 and 23.
!----------------------------------------------------------------------
60 if (zone.lt.13) then
     change=zone-13
     lhour=gmt+change
  elseif (zone.eq.13) then
     lhour=gmt
  else
     change=zone-13
     lhour=gmt+change
  end if
  if (lhour.lt.0) lhour=lhour+24
  if (lhour.gt.23) lhour=lhour-24
  return
end subroutine LIS_gmt2localtime
#endif
!BOP
! 
! !ROUTINE: LIS_localtime2gmt
!  \label{LIS_localtime2gmt}
! 
! !DESCRIPTION: 
! 
! Calculates GMT based on the local time 
! 
! !INTERFACE:
subroutine LIS_localtime2gmt (gmt,lon,lhour,zone)
! !ARGUMENTS:
  real:: gmt                 ! GMT time (0-23)
  real:: lon                 ! longitude in degrees
  real:: change              ! the change in number of hours between
  real:: lhour               ! local hour (0-23) 0= midnight, 23= 11:00 p.m.
       
  integer::  i              ! working integer
  integer:: zone            ! time zone (1-24)
!EOP
!----------------------------------------------------------------------
!  Determine into which time ZONE (15 degree interval) the
!  longitude falls.
!----------------------------------------------------------------------
  do i=1,25
     if (lon.lt.(-187.5+(15*i))) then
        zone=i
        if (zone.eq.25) zone=1
        go to 60
     end if
  end do
!----------------------------------------------------------------------
!     Calculate change (in number of hours) from GMT time to
!     local hour.  Change will be negative for zones < 13 and
!     positive for zones > 13.
!     There is also a correction for LHOUR < 0 and LHOUR > 23
!     to LHOUR between 0 and 23.
!----------------------------------------------------------------------
60 if (zone.lt.13) then
     change=zone-13
!     lhour=gmt+change
     gmt = lhour-change
  elseif (zone.eq.13) then
!     lhour=gmt
     gmt = lhour
  else
     change=zone-13
!     lhour=gmt+change
     gmt = lhour-change
  end if
!  if (lhour.lt.0) lhour=lhour+24
!  if (lhour.gt.23) lhour=lhour-24
  return
end subroutine LIS_localtime2gmt

!BOP
! !ROUTINE: LIS_registerAlarm
! \label{LIS_registerAlarm}
!
! !INTERFACE:
subroutine LIS_registerAlarm(name, ts, interval, intervalType, alarm_offset, &
                             dek_offset, when)
  
! !USES:
  use LIS_logMod, only : LIS_logunit

  implicit none 
  
! !ARGUMENTS: 
  character(len=*)           :: name
  real                       :: ts
  real                       :: interval
  character(len=*), optional :: intervalType
  integer, optional          :: alarm_offset
  integer, optional          :: dek_offset
  character(len=*), optional :: when

  type(lisalarmEntry), pointer :: alrmEntry, current
  integer                      :: yr,mo,da,hr,mn,ss,rc
  type(ESMF_Time)              :: currTime
  type(ESMF_TimeInterval)      :: deltaT
  integer                      :: refjuld, sjuld

! !DESCRIPTION:
! \begin{description}
!   \item[name]
!     name of the alarm
!   \item[interval]
!     frequency of the alarm in seconds
!   \item[intervalType]
!     type of the alarm for non-interval-based alarms.
!     Acceptable values are:
!     ``monthly'', ``quarterly'', or ``dekad''
!   \item[alarm\_offset]
!     offset in seconds to add to beginning of alarm
!   \item[dek\_offset] (for dekadal alarms only)
!     offset in seconds to add to beginning and ending of dekad ranges 
!   \item[when] (for dekadal alarms only)
!     specifies when the alarm should ring; a value of ``begin'' specifies
!     that the alarm should ring at the beginning of the dekad;
!     a value of ``end'' specifies that the alarm should ring at the end
!     of the dekad.
! \end{description}
!EOP

  allocate(alrmEntry)

  alrmEntry%name         = name
  alrmEntry%interval     = interval
  alrmEntry%ts           = ts
  alrmEntry%prev_mo      = -1
  alrmEntry%monthCount   = 0
  alrmEntry%firstInstance = .true. 
  alrmEntry%alarm_offset = 0
  alrmEntry%dek_offset   = 0
  alrmEntry%ref_dekad    = 0
  alrmEntry%when         = ""

  if(present(intervalType)) then 
     alrmEntry%intervaltype = intervalType
  else
     alrmEntry%intervaltype = ""
  endif

  if (present(alarm_offset) ) then
     alrmEntry%alarm_offset = alarm_offset
  else
     alrmEntry%alarm_offset = 0
     !write(LIS_logunit,*) 'MSG: LIS_registerAlarm: '     // &
     !                     'alarm_offset for '//trim(alrmEntry%name) // &
     !                     ' is not defined.  Defaulting to 0.'
  endif

  if ( present(intervalType) ) then
    if ( intervalType == "dekad" ) then 
       if (present(dek_offset) ) then
          alrmEntry%dek_offset = dek_offset
       else
          alrmEntry%dek_offset = 0
          write(LIS_logunit,*) '[WARN]  LIS_registerAlarm: '     // &
                'dek_offset for '//trim(alrmEntry%name) // &
                ' is not defined.  Defaulting to 0.'
       endif

       if (present(when) ) then
          alrmEntry%when = when
       else
          alrmEntry%when = "end"
          write(LIS_logunit,*) '[WARN] LIS_registerAlarm: '     // &
                               'when for '//trim(alrmEntry%name) // &
                               ' is not defined.  Defaulting to "end".'
       endif
    endif
  endif

  alrmEntry%alarmTime = 0.0
  alrmEntry%next       => null()
       
  if(.not.associated(LIS_alarmlist)) then 
     LIS_alarmList => alrmEntry
  else
     current => LIS_alarmList
     do while(associated(current%next)) 
        current=> current%next
     enddo
     current%next => alrmEntry
  endif

end subroutine LIS_registerAlarm


!BOP
! !ROUTINE: isAlarmRinging
! \label{isAlarmRinging}
!
! !INTERFACE:
!  function isAlarmRinging(LIS_rc, model_TS, alarmInterval)
  function isAlarmRinging(LIS_rc, alarmName, intervaltype)
! !USES:

   implicit none
   ! !ARGUMENTS: 
   type(lisrcdec)             :: LIS_rc    
   character(len=*)           :: alarmName
   character(len=*), optional :: intervaltype
   logical                    :: isAlarmRinging
!
! !DESCRIPTION:
!  Checks if the alarm is ringing. Currently hourly, daily, 
!  and monthly interval alarms are supported. 
!
!  The arguments are: 
!  \begin{description}
!   \item[LIS\_rc]
!    instance of the {\tt lis\_module}
!   \item[alarmInterval]
!    alarmInterval in seconds
!   \item [model\_TS]
!    timestep of the model 
!   \item[isAlarmRinging]
!    flag indicating the status of the call
!  \end{description}
!   
!EOP
   logical                  :: alarm_found
   integer                  :: mfactor
   type(ESMF_Time)          :: currTime
   type(ESMF_TimeInterval)  :: ts
   integer                  :: yr, mo,da,hr,mn,ss
   integer                  :: status
   character*100            :: itype
   type(lisalarmEntry), pointer :: current
   integer :: yr2,mo2
   integer :: numi, zeroi,doy1,doy2
   real    :: gmt2,gmt1
   
   real*8  :: time2,time
   real*8  :: jan31,apr30,jul31,oct31 ! dates of quarterly intervals
   integer :: janda,janmo          ! january 31 
   integer :: aprda,aprmo          ! april 30
   integer :: julda,julmo          ! july 31
   integer :: octda,octmo          ! october 31 
   integer :: rc
   real*8  :: quartTime 
   real*8  :: t1
   real    :: gmt
   integer :: refjuld, sjuld,doy
   integer :: ts_frac_s, ts_frac_ms

   zeroi = 0 
   numi = 0 
   isAlarmRinging = .false. 
   alarm_found = .false. 
   numi = 0
   zeroi = 0
   if(present(intervalType)) then 
      itype = intervalType
   else
      itype = ""
   endif
   
   if(.not.associated(LIS_alarmList)) then 
      write(LIS_logunit,*) '[ERR] Alarm: '//alarmName//' is not registered'
      write(LIS_logunit,*) '[ERR] Program stopping..'
      call LIS_endrun()
   else
      current =>LIS_alarmList
      do while(current%name.ne.alarmName) 
         current => current%next
      enddo
      if(current%name.eq.alarmName) then 
         alarm_found = .true. 
      endif
   endif
   
   if(.not.alarm_found) then 
      write(LIS_logunit,*) '[ERR] Alarm: '//alarmName//' is not found in the '
      write(LIS_logunit,*) '[ERR] list of LIS alarms. Program stopping..'
      call LIS_endrun()
   endif
   
   ! "30-day" month interval based alarms:
   if(current%interval.eq.2592000.0.and.&
        itype.eq."monthly") then 
      if(LIS_rc%da.lt.16) then 
         mo2=LIS_rc%mo
         yr2=LIS_rc%yr
      else
         mo2=LIS_rc%mo+1
         yr2=LIS_rc%yr
         if(mo2.eq.13)then
            mo2=1
            yr2=LIS_rc%yr+1
         endif
      endif
      call LIS_date2time(time2,doy2,gmt2,yr2,mo2,&
           numi,zeroi,zeroi,zeroi)       

      if(time2.gt.current%alarmTime) then 
         isAlarmRinging = .true.
         current%alarmTime = time2
      endif

   ! Seasonal based alarms (valid: Jan31,Apr30,Jul31,and Oct31):
   elseif(itype.eq."quarterly") then 
      zeroi = 0 
      time=LIS_rc%time
      yr=LIS_rc%yr
      janda=31
      janmo=01
      call LIS_date2time(jan31,doy1,gmt1,yr,janmo,&
           janda,zeroi,zeroi,zeroi)
      aprda=30
      aprmo=04
      call LIS_date2time(apr30,doy1,gmt1,yr,aprmo,&
           aprda,zeroi,zeroi,zeroi)
      julda=31
      julmo=07
      call LIS_date2time(jul31,doy1,gmt1,yr,julmo,&
           julda,zeroi,zeroi,zeroi)
      octda=31
      octmo=10
      call LIS_date2time(oct31,doy1,gmt1,yr,octmo,&
           octda,zeroi,zeroi,zeroi)
      if ( time.ge.jan31 .and. time.lt.apr30 ) then
         quartTime = 1
      elseif ( time.ge.apr30 .and. time.lt.jul31 ) then
         quartTime = 2
      elseif ( time.ge.jul31 .and. time.lt.oct31 ) then
         quartTime = 3
      elseif ( time.ge.oct31 ) then
         quartTime = 4
      elseif ( time.lt.jan31) then
         quartTime = 4
      endif

      if(current%alarmtime .ne. quartTime) then 
         isAlarmRinging = .true.
         current%alarmTime = quartTime
      endif

   ! "30-day" month interval based alarms:
   elseif(current%interval.eq.2592000.0) then !monthly
      mfactor = current%interval/2592000

      ! for monthly timesteps, we check at the end of the month to trigger
      ! the alarm. 

      call ESMF_TimeSet(currTime,  yy=LIS_rc%yr, &
           mm = LIS_rc%mo, &
           dd = LIS_rc%da, &
           h = LIS_rc%hr, &
           m = LIS_rc%mn, &
           s = LIS_rc%ss,&
           calendar = LIS_calendar, &
           rc=status)
      call LIS_verify(status, 'error in ESMF_TimeSet: LIS_timeMgrMod')
      
      ts_frac_s = nint(current%ts)
      ts_frac_ms = (current%ts - nint(current%ts))*1000

      call ESMF_TimeIntervalSet(ts,s=ts_frac_s, ms = ts_frac_ms,rc=status)
      call LIS_verify(status, 'error in ESMF_TimeIntervalSet: LIS_timeMgrMod')

      currTime = currTime + ts
      call ESMF_TimeGet(currTime,  yy=yr, &
           mm = mo, &
           dd = da, &
           h = hr, &
           m = mn, &
           s = ss,&
           calendar = LIS_calendar, &
           rc=status)
      call LIS_verify(status, 'error in ESMF_TimeSet: LIS_timeMgrMod')

      if(mo.ne.current%prev_mo) then
         current%monthCount = current%monthCount + 1
         if(current%monthCount.eq.mfactor.and.&
              current%prev_mo.ne.-1) then 
            isAlarmRinging = .true.
         endif
         if(current%monthCount.eq.mfactor) then 
            current%monthCount = 0
         endif
         current%prev_mo = mo
      else
         isAlarmRinging = .false.
      endif

   ! Dekad-based (1st, 11th,& 21st of each month) alarm:
   !elseif(current%interval.eq.864000.0) then !dekadal
   elseif(itype == "dekad") then

      isAlarmRinging = LIS_isDekadalAlarmRinging(LIS_rc, current)

   ! Alarm for just 10-day interval:
   elseif(itype == "10-day") then 
      isAlarmRinging = checkForDekadInterval(&
           secsFrmBegYr(LIS_rc%doy, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss),&
           LIS_rc%yr, &
           timeSchema_arg=10, &
           firstInstance=current%firstInstance)
      current%firstInstance = .false. 

   ! Alarm for intervals at daily and less (e.g., hour, or minutes):
   elseif(mod(real(LIS_rc%hr)*3600+60*real(LIS_rc%mn)+    &
              real(LIS_rc%ss)+real(current%alarm_offset), &
          current%interval).eq.0.0) then 
      isAlarmRinging = .true. 
   endif

 end function IsAlarmRinging


  function checkForDekadInterval(secsFrmBegYr, yr, &
       timeSchema_arg, dekOy_arg, firstInstance)
    
   logical :: checkForDekadInterval
   logical :: firstInstance
!hardcoded
   integer*4, parameter :: gDEKAD_TIMESTEPS_PER_YEAR = 36
   integer*4, parameter :: gDEKAD_TIME_SCHEMA = 10
   
   integer*4, intent(in) :: secsFrmBegYr
   integer, intent(in) :: yr
   integer*4, intent(in), optional :: timeSchema_arg
   integer*4, intent(inout), optional :: dekOy_arg

   integer*4 :: timeSchema 
   integer*4 :: ts
   integer*4, save :: dekOy = (0-1)
!!! Dekad time schema variables begin !!!
! Remove the guess work re. ts start hour:  We should by default always be using offset 0600 per corresp. Jul 2011
   ! Dekad time step seconds from beginning of the year non- leap year
   integer*4, parameter, dimension(gDEKAD_TIMESTEPS_PER_YEAR) :: dekTSsfbyNLY = (/&
#if 0 
        885600  ,&
        1749600 ,&
        2700000 ,&
        3564000 ,&
        4428000 ,&
        5119200 ,&
        5983200 ,&
        6847200 ,&
        7797600 ,&
        8661600 ,&
        9525600 ,&
        10389600,&
        11253600,&
        12117600,&
        13068000,&
        13932000,&
        14796000,&
        15660000,&
        16524000,&
        17388000,&
        18338400,&
        19202400,&
        20066400,&
        21016800,&
        21880800,&
        22744800,&
        23608800,&
        24472800,&
        25336800,&
        26287200,&
        27151200,&
        28015200,&
        28879200,&
        29743200,&
        30607200,&
        31557600 &
        /)
#endif
        864000  ,&
        1728000 ,&
        2678400 ,&
        3542400 ,&
        4406400 ,&
        5097600 ,&
        5961600 ,&
        6825600 ,&
        7776000 ,&
        8640000 ,&
        9504000 ,&
        10368000,&
        11232000,&
        12096000,&
        13046400,&
        13910400,&
        14774400,&
        15638400,&
        16502400,&
        17366400,&
        18316800,&
        19180800,&
        20044800,&
        20995200,&
        21859200,&
        22723200,&
        23587200,&
        24451200,&
        25315200,&
        26265600,&
        27129600,&
        27993600,&
        28857600,&
        29721600,&
        30585600,&
        31536000 &
        /)

   ! Dekad time step seconds from beginning of the year leap year
   integer*4, parameter, dimension(gDEKAD_TIMESTEPS_PER_YEAR) :: dekTSsfbyLY = (/&
#if 0
        885600  ,&
        1749600 ,&
        2700000 ,&
        3564000 ,&
        4428000 ,&
        5205600 ,&
        6069600 ,&
        6933600 ,&
        7884000 ,&
        8748000 ,&
        9612000 ,&
        10476000,&
        11340000,&
        12204000,&
        13154400,&
        14018400,&
        14882400,&
        15746400,&
        16610400,&
        17474400,&
        18424800,&
        19288800,&
        20152800,&
        21103200,&
        21967200,&
        22831200,&
        23695200,&
        24559200,&
        25423200,&
        26373600,&
        27237600,&
        28101600,&
        28965600,&
        29829600,&
        30693600,&
        31644000 &
#endif
        864000  ,&
        1728000 ,&
        2678400 ,&
        3542400 ,&
        4406400 ,&
        5184000 ,&
        6048000 ,&
        6912000 ,&
        7862400 ,&
        8726400 ,&
        9590400 ,&
        10454400,&
        11318400,&
        12182400,&
        13132800,&
        13996800,&
        14860800,&
        15724800,&
        16588800,&
        17452800,&
        18403200,&
        19267200,&
        20131200,&
        21081600,&
        21945600,&
        22809600,&
        23673600,&
        24537600,&
        25401600,&
        26352000,&
        27216000,&
        28080000,&
        28944000,&
        29808000,&
        30672000,&
        31622400 &
        /)
!!! Dekad time schema variables end !!!
   
   
   checkForDekadInterval = .false.
   timeSchema = gDEKAD_TIME_SCHEMA ! assume dekad-time-schema* unless specified otherwise
   if(present(timeSchema_arg)) timeSchema = timeSchema_arg
   
   select case(timeSchema)
   case (gDEKAD_TIME_SCHEMA)

! From Tamuka Magadzire e-mail dated February 16, 2011:
!> For each month:
!> Dekad 1: the first 10 days of the month Dekad 2: the 11th to the 20th 
!> of the month Dekad 3: whatever days are left in the month after the 
!> 20th, whether its 8 days, 9, 10 or 11 days.

! Accordingly, see 'timestepTablesDekads2011.xls' (Excel 97-2003 workbook) 
! which was created to give the above table/s.
      
      if((mod(yr,4) /= 0) .or. ((mod(yr,100) == 0) &
           .and. (mod(yr,400) /= 0) )  ) then
         ! non- leap year
         
         if(firstInstance) then 
            checkForDekadInterval = .true. 
!find the closest dekad. 
            do ts=1, gDEKAD_TIMESTEPS_PER_YEAR, 1
               if( secsFrmBegYr .le.dekTSsfbyNLY(ts) ) then
                  dekOy = ts
                  exit;
               endif
            end do
         else
            do ts=1, gDEKAD_TIMESTEPS_PER_YEAR, 1
               if( secsFrmBegYr == dekTSsfbyNLY(ts) ) then
                  checkForDekadInterval = .true.
                  dekOy = ts
               endif
            end do
         endif
      else
         ! is leap year
         if(firstInstance) then 
            checkForDekadInterval = .true. 
!find the closest dekad. 
            do ts=1, gDEKAD_TIMESTEPS_PER_YEAR, 1
               if( secsFrmBegYr .le.dekTSsfbyLY(ts) ) then
                  dekOy = ts
                  exit;
               endif
            end do
         else
            do ts=1, gDEKAD_TIMESTEPS_PER_YEAR, 1
               if( secsFrmBegYr == dekTSsfbyLY(ts) ) then
                  checkForDekadInterval = .true.
                  dekOy = ts
               endif
            end do
         endif
      endif
! Seconds from beginning of year is used instead of other temporal representations because:
!   1) this was given the most thought and therefore closest to completion at time of this writing
!   2) this representation represents the minimal requirements for any code that may use a time loop*
!   3) this exacting approach should require less logic in LIS for discrete time steps
!  *including both LIS and standalone model codes

   case default
      ! unhandled time schema
   end select

   if(present(dekOy_arg)) dekOy_arg = dekOy
   
 end function checkForDekadInterval

 function secsFrmBegYr(doy, hr, mm, ss)
    
    integer*4 :: secsFrmBegYr
    
    integer, intent(in) :: doy
    integer, intent(in) :: hr
    integer, intent(in) :: mm
    integer, intent(in) :: ss
    
    secsFrmBegYr = (86400*doy) + (3600*hr) + (60*mm) + ss
    
  end function secsFrmBegYr
!BOP
! !ROUTINE: LIS_parseTimeString
! \label{LIS_parseTimeString}
!
! !INTERFACE: 
  subroutine LIS_parseTimeString(inputstr,time_value)
    
    implicit none
! !ARGUMENTS:    
    character(len=*)          :: inputstr
    real                      :: time_value
! 
! !DESCRIPTION: 
!  This routine parses the input time string and converts
!  it to an integer value, which represents the time in seconds. 
!  The input time string should have 2-character suffixes of 
!  ss (seconds), mn (minutes), hr (hours), da (days), mo (month)
!  or yr (years). 
!   
!EOP       
    character*2        :: suffix
    real               :: time_input
    integer            :: iloc

    iloc  = len(trim(inputstr))-1
    suffix = inputstr(iloc:len(inputstr))
    
    read(inputstr(1:iloc-1),*) time_input

    if(suffix.eq."ss") then 
       time_value = time_input
    elseif(suffix.eq."mn") then 
       time_value = time_input*60.0
    elseif(suffix.eq."hr") then
       time_value = time_input*3600.0
    elseif(suffix.eq."da") then 
       time_value = time_input*86400.0
    elseif(suffix.eq."mo") then 
       time_value = time_input*2592000.0
    elseif(suffix.eq."yr") then 
       time_value = time_input*31536000.0
    else
       print*, 'The time specification ',trim(inputstr),&
            'is incorrect. Program stopping....'
       stop
    endif
    
  end subroutine LIS_parseTimeString


  subroutine LIS_mon3char( mo2digchar, mo3char )

    implicit none
    character*2, intent(in)  :: mo2digchar
    character*3, intent(out) :: mo3char

    select case( mo2digchar )
     case( "01" )
       mo3char = "jan"
     case( "02" )
       mo3char = "feb"
     case( "03" )
       mo3char = "mar"
     case( "04" )
       mo3char = "apr"
     case( "05" )
       mo3char = "may"
     case( "06" )
       mo3char = "jun"
     case( "07" )
       mo3char = "jul"
     case( "08" )
       mo3char = "aug"
     case( "09" )
       mo3char = "sep"
     case( "10" )
       mo3char = "oct"
     case( "11" )
       mo3char = "nov"
     case( "12" )
       mo3char = "dec"
     case default
       write(*,*) " Don't recognize this 2-digit month ... stopping. "
       stop
    end select

  end subroutine LIS_mon3char

end module LIS_timeMgrMod
