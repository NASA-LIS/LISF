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
! !MODULE: LVT_timeMgrMod
! \label(LVT_timeMgrMod)
!
! !INTERFACE:
module LVT_timeMgrMod
! 
! !USES:
  use ESMF
  use LVT_PRIV_rcMod
  use LVT_logMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This module contains routines for time managment.The module provides
! routines for clock initialization, model timestepping, some basic
! alarm functions, and other useful time managment utilities. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  
  PRIVATE 

  integer, parameter :: uninit_int = -999999999 

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LVT_timemgr_init     ! time manager initialization
  PUBLIC :: LVT_timemgr_print    ! print the details of the clock
  PUBLIC :: LVT_advance_timestep ! increment timestep number
  PUBLIC :: LVT_update_clock
  PUBLIC :: LVT_update_timestep
  PUBLIC :: LVT_timemgr_set      ! set the current time
  PUBLIC :: LVT_get_step_size    ! return step size in seconds
  PUBLIC :: LVT_get_nstep        ! return timestep number
  PUBLIC :: LVT_get_curr_date  ! return date components at end of current timestep
  PUBLIC :: LVT_get_curr_calday ! return calendar day at end of current timestep
  PUBLIC :: LVT_get_julhr       ! returns time in julian hours
  PUBLIC :: LVT_get_julss       ! returns time in julian seconds
  PUBLIC :: LVT_is_last_step    ! return true on last timestep
  PUBLIC :: LVT_tick            ! clock advance/retract by a certain amount
  PUBLIC :: LVT_date2time       ! converts date to the LVT time format
  PUBLIC :: LVT_time2date       ! converts LVT time to the date format
  PUBLIC :: LVT_setAlarm        ! sets a monthly/weekly/hourly alarm
  PUBLIC :: LVT_setMonthlyAlarm ! sets a monthly alarm
  PUBLIC :: LVT_setHourlyAlarm  ! sets a hourly alamr
  PUBLIC :: LVT_isMonthlyAlarmRinging !checks to see if the monthly alarm is ringing
  PUBLIC :: LVT_isHourlyAlarmRinging !checks to see if the hourly alarm is ringing
  PUBLIC :: LVT_computeTimeBookEnds !Compute the time interpolation weights 
  PUBLIC :: LVT_computeTemporalWeights !Compute monthly time interpolation weights
  PUBLIC :: LVT_parseTimeString
  PUBLIC :: LVT_julhr_date  ! converts julian hour to date format
  PUBLIC :: LVT_seconds2time ! ?? 
  PUBLIC :: LVT_computeTimeSpan
  PUBLIC :: LVT_getTimeSpanIndex
  PUBLIC :: LVT_computeTimeSpanInYears
  PUBLIC :: LVT_getYearIndex
  PUBLIC :: LVT_doy2moda
  PUBLIC :: LVT_resetClock
  PUBLIC :: LVT_tmjul4
  PUBLIC :: LVT_localtime2gmt
  PUBLIC :: LVT_localtime
  PUBLIC :: LVT_clock
  PUBLIC :: LVT_calendar
  PUBLIC :: LVT_obsSmTwI
  PUBLIC :: LVT_obsSmTwL
!EOP  
  type(ESMF_Clock), save        :: LVT_clock
  type(ESMF_Calendar), save     :: LVT_calendar  
  type(ESMF_TimeInterval), save :: LVT_obsSmTwL
  type(ESMF_TimeInterval), save :: LVT_obsSmTwI


!BOP
! 
! !ROUTINE: isMonthlyAlarmRinging
! \label{isMonthlyAlarmRinging}
!
! !INTERFACE: 
   interface LVT_isMonthlyAlarmRinging
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  checks if the monthly alarm is ringing. The private functions
!  have different arguments based on the input options specified. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !PRIVATE MEMBER FUNCTIONS: 
      module procedure isMonthlyAlarmRinging1
      module procedure isMonthlyAlarmRinging2
!EOP
   end interface

!BOP
! 
! !ROUTINE: LVT_computeTemporalWeights
! \label{LVT_computeTemporalWeights}
!
! !INTERFACE:
   interface LVT_computeTemporalWeights
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  computes the temporal interpolation weights. The private functions
!  have different arguments based on the input options specified. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !PRIVATE MEMBER FUNCTIONS: 
      module procedure computeTemporalWeights1
      module procedure computeTemporalWeights2
!EOP
   end interface
contains

!BOP
! 
! !ROUTINE: LVT_timemgr_init
! \label{LVT_timemgr_init}
!
! !INTERFACE:
  subroutine LVT_timemgr_init(LVT_rc)
! 
! !USES: 

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! Initialize the LVT time manager.
!
! NOTE - This assumes that the LVT time specific variables 
! pertaining to start and end times have been set before 
! this routine is called.  
!  \begin{description}
!   \item [LVT\_rc]
!     instance of the {\tt lvt\_module}
!  \end{description}
! 
! The calling sequence is:  
! \begin{description}
!  \item[date2time](\ref{date2time}) \newline
!   to convert current date to a floating point format
!  \item[timemgr\_print]
!   display the contents of the time manager
! \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS:
    type(lvtrcdec) :: LVT_rc
!EOP
    type(ESMF_Time)         :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: status

    call LVT_date2time(LVT_rc%etime1,LVT_rc%edoy1,LVT_rc%egmt1, & 
         LVT_rc%eyr1,LVT_rc%emo1,LVT_rc%eda1,LVT_rc%ehr1,LVT_rc%emn1,LVT_rc%ess1)
    LVT_rc%tscount = 0
    LVT_calendar = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN, &
         name="Gregorian", rc=status)

    call ESMF_TimeIntervalSet(timeStep, s = LVT_rc%ts, &
         rc=status)
    call LVT_verify(status,'ESMF_TimeIntervalSet failed')

    call ESMF_TimeSet(startTime, yy = LVT_rc%syr, &
         mm = LVT_rc%smo, &
         dd = LVT_rc%sda, &
         h  = LVT_rc%shr, &
         m  = LVT_rc%smn, & 
         s  = LVT_rc%sss, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status, 'ESMF_TimeSet failed')

    call ESMF_TimeSet(stopTime, yy = LVT_rc%eyr, &
         mm = LVT_rc%emo, &
         dd = LVT_rc%eda, &
         h  = LVT_rc%ehr, &
         m  = LVT_rc%emn, & 
         s  = LVT_rc%ess, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status, 'ESMF_TimeSet failed')

    LVT_clock = ESMF_ClockCreate(name="LVT Clock", timestep=timestep, startTime=startTime, &
         stopTime=stopTime, rc=status)

    call LVT_timemgr_print(LVT_rc)
    
    LVT_rc%dayCount = 0 
    LVT_rc%monthCount = 0 

  end subroutine LVT_timemgr_init

!BOP
! 
! !ROUTINE: LVT_resetclock
! \label{LVT_resetclock}
!
! !INTERFACE:
  subroutine LVT_resetclock(LVT_rc)
! 
! !USES:

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! Resets the LVT time manager.
!
! NOTE - This assumes that the LVT time specific variables 
! pertaining to start and end times have been set before 
! this routine is called.  
!  \begin{description}
!   \item [LVT\_rc]
!     instance of the {\tt lvt\_module}
!  \end{description}
! 
! The calling sequence is:  
! \begin{description}
!  \item[date2time](\ref{date2time}) \newline
!   to convert current date to a floating point format
!  \item[timemgr\_print]
!   display the contents of the time manager
! \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS:
    type(lvtrcdec) :: LVT_rc
!EOP
    type(ESMF_TimeInterval)     :: timeStep
    type(ESMF_Time)             :: startTime, stopTime
    integer                     :: status

    call LVT_date2time(LVT_rc%etime1,LVT_rc%edoy1,LVT_rc%egmt1, & 
         LVT_rc%eyr1,LVT_rc%emo1,LVT_rc%eda1,LVT_rc%ehr1,LVT_rc%emn1,LVT_rc%ess1)
    LVT_rc%tscount = 0
    LVT_rc%endtime=0
    call LVT_timemgr_print(LVT_rc)
    
    call ESMF_TimeIntervalSet(timeStep, s = LVT_rc%ts, &
         rc=status)
    call LVT_verify(status, 'ESMF_TimeIntervalSet failed')

    call ESMF_TimeSet(startTime, yy = LVT_rc%syr, &
         mm = LVT_rc%smo, &
         dd = LVT_rc%sda, &
         h  = LVT_rc%shr, &
         m  = LVT_rc%smn, & 
         s  = LVT_rc%sss, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status, 'ESMF_TimeSet failed')

    call ESMF_TimeSet(stopTime, yy = LVT_rc%eyr, &
         mm = LVT_rc%emo, &
         dd = LVT_rc%eda, &
         h  = LVT_rc%ehr, &
         m  = LVT_rc%emn, & 
         s  = LVT_rc%ess, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status, 'ESMF_TimeSet failed')

    LVT_clock = ESMF_ClockCreate(name="LVT Clock", &
         timestep=timestep, startTime=startTime, &
         stopTime=stopTime, rc=status)

    LVT_rc%yr= LVT_rc%syr
    LVT_rc%mo= LVT_rc%smo
    LVT_rc%da= LVT_rc%sda
    LVT_rc%hr= LVT_rc%shr
    LVT_rc%mn= LVT_rc%smn
    LVT_rc%ss= LVT_rc%sss

    
    call LVT_date2time(LVT_rc%time,LVT_rc%doy,LVT_rc%gmt, &
         LVT_rc%yr,LVT_rc%mo,LVT_rc%da,LVT_rc%hr,LVT_rc%mn,LVT_rc%ss)

    LVT_rc%prev_mo_tavg = LVT_rc%mo
    LVT_rc%prev_yr_tavg = LVT_rc%yr
    LVT_rc%dayCount = 0 
    LVT_rc%monthCount = 0 

    LVT_rc%prev_mo_sout = LVT_rc%mo
    LVT_rc%prev_yr_sout = LVT_rc%yr
!    LVT_rc%prev_mo_sout = -1
!    LVT_rc%prev_yr_sout = -1
    LVT_rc%dayCount_sout = 0 
    LVT_rc%monthCount_sout = 0 
    LVT_rc%resetFlag = .true. 

  end subroutine LVT_resetclock

!BOP
! 
! !ROUTINE: timemgr_set
! \label{timemgr_set}
!
! !INTERFACE:
  subroutine LVT_timemgr_set(LVT_rc,yr,mo,da,hr,mn,ss)
! 
! !USES: 
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! sets the time manager clock based on the input time 
! specification. This method is used to initialize the 
! time mangager when the clock is passed down to LVT from
! a parent component
!
!  \begin{description}
!   \item [LVT\_rc]
!     instance of the {\tt lvt\_module}
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
!   \item[date2time](\ref{date2time}) \newline
!    convert date to a flotating point format
! \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    type(lvtrcdec) :: LVT_rc
    integer :: yr,mo,da,hr,mn,ss
!EOP
    type(ESMF_Time) :: currTime
    integer         :: status

    LVT_rc%yr = yr
    LVT_rc%mo = mo
    LVT_rc%da = da
    LVT_rc%hr = hr
    LVT_rc%mn = mn
    LVT_rc%ss = ss
    
    call LVT_date2time(LVT_rc%time,LVT_rc%doy,LVT_rc%gmt,& 
         LVT_rc%yr, LVT_rc%mo, LVT_rc%da, LVT_rc%hr, LVT_rc%mn, LVT_rc%ss)
    
    LVT_rc%tscount = LVT_rc%tscount + 1
    
    write(unit=LVT_logunit,fmt=24)' [INFO] LVT cycle time: ',LVT_rc%mo,'/',LVT_rc%da,'/', & 
         LVT_rc%yr,LVT_rc%hr,':',LVT_rc%mn,':',LVT_rc%ss
24  format(a24,i2.2,a1,i2.2,a1,i4,1x,i2.2,a1,i2.2,a1,i2.2)    ! MLF

    call ESMF_TimeSet(currTime, yy = yr, &
         mm = mo, &
         dd = da, &
         h  = hr, &
         m  = mn, & 
         s  = ss, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status, 'ESMF_TimeSet failed')

    call ESMF_ClockSet(LVT_clock, currTime = currTime, rc=status)
    call LVT_verify(status, 'ESMF_ClockSet failed')

  end subroutine LVT_timemgr_set

!BOP
! 
! !ROUTINE: timemgr_print
! \label{timemgr_print}
!
! !INTERFACE:
  subroutine LVT_timemgr_print(LVT_rc)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! Prints the details of the LVT time manager such as the start, stop
! times, current time, timestep, etc. 
!
!  \begin{description}
!   \item [LVT\_rc]
!     instance of the {\tt lvt\_module}
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
!EOP
! !ARGUMENTS: 
    type(lvtrcdec) :: LVT_rc
!EOP 

    write(unit=LVT_logunit,fmt=*)'[INFO] ************************************************'
    
    write(unit=LVT_logunit,fmt=*)'[INFO] Timestep size (seconds):  ', LVT_rc%ts
    write(unit=LVT_logunit,fmt=*)'[INFO] Start date (ymd tod):     ', LVT_rc%syr, LVT_rc%smo, LVT_rc%sda
    write(unit=LVT_logunit,fmt=*)'[INFO] Stop date (ymd tod):      ', LVT_rc%eyr, LVT_rc%emo, LVT_rc%eda
    write(unit=LVT_logunit,fmt=*)'[INFO] Current date (ymd tod):   ', LVT_rc%yr, LVT_rc%mo, LVT_rc%da
    
    write(unit=LVT_logunit,fmt=*)'[INFO] ************************************************'
  end subroutine LVT_timemgr_print

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_advance_timestep
! \label{LVT_advance_timestep} 
!
! !INTERFACE:
  subroutine LVT_advance_timestep(LVT_rc)
! !USES: 

    implicit none
! !ARGUMENTS:    
    type(lvtrcdec) :: LVT_rc
! !DESCRIPTION:
! 
! Increments time manager clock by the model timestep. In case of LVT 
! running multiple nests, each running at different timesteps, the 
! smallest timestep is chosen to advance the clock. 
!
!  \begin{description}
!   \item [LVT\_rc]
!     instance of the {\tt lvt\_module}
!  \end{description}
!
!  The calling sequence is:  
! \begin{description}
!   \item[date2time](\ref{date2time}) \newline
!    convert date to a flotating point format
! \end{description}
!EOP

! Local variables
    integer :: days(12),tda,n
    real    :: curr_time
    integer :: yr,mo,da,hr,mn,ss,status
    type(ESMF_Time) :: advTime,currTime
    type(ESMF_Time) :: start_dekad1,start_dekad2,start_dekad3
    type(ESMF_TimeInterval) :: dekad_interval, month_interval

    data days /31,28,31,30,31,30,31,31,30,31,30,31/
    
24  format(a24,i2.2,a1,i2.2,a1,i4,1x,i2.2,a1,i2.2,a1,i2.2)  ! MLF

    if(LVT_rc%tsconv.eq."regular") then 
       if(LVT_rc%ts.eq.2592000) then !special case of month
          
          LVT_rc%mo = LVT_rc%mo + 1
          
          do while(LVT_rc%mo .gt. 12) 
             LVT_rc%mo = LVT_rc%mo-12
             LVT_rc%yr = LVT_rc%yr +1
          enddo
          
          LVT_rc%tscount = LVT_rc%tscount + 1       
          
          call LVT_date2time(LVT_rc%time,LVT_rc%doy,LVT_rc%gmt,& 
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da, LVT_rc%hr, LVT_rc%mn, LVT_rc%ss)
          
       else
          call ESMF_ClockAdvance(LVT_clock,rc=status)
          
          LVT_rc%ss = LVT_rc%ss + LVT_rc%ts
          
          do while(LVT_rc%ss .gt. 59) 
             LVT_rc%ss = LVT_rc%ss - 60 
             LVT_rc%mn = LVT_rc%mn + 1
          enddo
          
          do while(LVT_rc%mn .gt.59)
             LVT_rc%mn = LVT_rc%mn -60
             LVT_rc%hr = LVT_rc%hr+1
          enddo
          
          do while(LVT_rc%hr .ge.24) 
             LVT_rc%hr = LVT_rc%hr -24
             LVT_rc%da = LVT_rc%da +1
          enddo
          
          if((mod(LVT_rc%yr,4) .eq. 0 .and. mod(LVT_rc%yr, 100).ne.0) &!leap year
               .or.(mod(LVT_rc%yr,400) .eq.0)) then 
             days(2) = 29
          else 
             days(2) = 28
          endif
          tda = days(LVT_rc%mo)
          do while(LVT_rc%da.gt.tda)
             LVT_rc%da = LVT_rc%da - days(LVT_rc%mo)
             LVT_rc%mo = LVT_rc%mo + 1
          enddo
          
          do while(LVT_rc%mo .gt. 12) 
             LVT_rc%mo = LVT_rc%mo-12
             LVT_rc%yr = LVT_rc%yr +1
          enddo
          
          call LVT_date2time(LVT_rc%time,LVT_rc%doy,LVT_rc%gmt,& 
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da, LVT_rc%hr, LVT_rc%mn, LVT_rc%ss)
          
          curr_time = float(LVT_rc%hr)*3600+60*float(LVT_rc%mn)+float(LVT_rc%ss)
          if(mod(curr_time,real(LVT_rc%ts)).eq.0) then 
             LVT_rc%tscount = LVT_rc%tscount + 1
          endif
       endif
    elseif(LVT_rc%tsconv.eq."dekad") then 
       yr=LVT_rc%yr
       mo=LVT_rc%mo
       da=LVT_rc%da
       hr=LVT_rc%hr
       mn=LVT_rc%mn
       ss=LVT_rc%ss
       
       if ( (mod(yr,4).eq.0.and.mod(yr,100).ne.0) .or. (mod(yr,400).eq.0) ) then
          ! leap year
          days(2) = 29
       else
          days(2) = 28
       endif
       
       call ESMF_TimeIntervalSet(dekad_interval, s=864000) ! 10 days in seconds
       call ESMF_TimeIntervalSet(month_interval, s=86400*days(mo))

       call ESMF_TimeSet(currTime, yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, rc=status)
       
       currTime = currTime + dekad_interval

       call ESMF_TimeGet(currTime, yy=yr, mm=mo, dd=da, h=hr, m=mn, s=ss, rc=status)

       call ESMF_TimeSet(start_dekad1, yy=yr, mm=mo, dd=1, h=0, &
            m=0, s=0, rc=status)       
       call ESMF_TimeSet(start_dekad2, yy=yr, mm=mo, dd=11, h=0, &
            m=0, s=0, rc=status)       
       call ESMF_TimeSet(start_dekad3, yy=yr, mm=mo, dd=21, h=0, &
            m=0, s=0, rc=status)
       
       if(currTime < start_dekad1) then 
          call ESMF_TimeSet(advTime, yy=yr, mm=mo, dd=1, &
               h=hr, m=mn, s=ss, rc=status)
       elseif(currTime .ge. start_dekad1 .and. currTime .lt. start_dekad2) then 
          call ESMF_TimeSet(advTime, yy=yr, mm=mo, dd=10, &
               h=hr, m=mn, s=ss, rc=status)
       elseif(currTime .ge. start_dekad2 .and. currTime .lt. start_dekad3) then 
          call ESMF_TimeSet(advTime, yy=yr, mm=mo, dd=20, &
               h=hr, m=mn, s=ss, rc=status)
       elseif(currTime.ge.start_dekad3) then 
          call ESMF_TimeSet(advTime, yy=yr, mm=mo, dd=days(mo), &
               h=hr, m=mn, s=ss, rc=status)          
       endif
       call ESMF_TimeGet(advTime, yy=LVT_rc%yr, mm=LVT_rc%mo, dd=LVT_rc%da, &
               h=LVT_rc%hr, m=LVT_rc%mn, s=LVT_rc%ss, rc=status)
       call LVT_date2time(LVT_rc%time,LVT_rc%doy,LVT_rc%gmt,& 
            LVT_rc%yr, LVT_rc%mo, LVT_rc%da, LVT_rc%hr, LVT_rc%mn, LVT_rc%ss)
       
       curr_time = float(LVT_rc%hr)*3600+60*float(LVT_rc%mn)+float(LVT_rc%ss)
       if(mod(curr_time,real(LVT_rc%ts)).eq.0) then 
          LVT_rc%tscount = LVT_rc%tscount + 1
       endif
    endif
    
    if(LVT_rc%endcode.eq.0)then  !end at real-time date (tbd)
       write(*,*)'warning: do not know how to stop in real-time' 
    endif
    if(LVT_rc%endcode.eq.1)then  !end on date specified in lvt configuration file
       call LVT_date2time(LVT_rc%etime,LVT_rc%edoy,LVT_rc%egmt, & 
            LVT_rc%eyr,LVT_rc%emo,LVT_rc%eda,LVT_rc%ehr,LVT_rc%emn,LVT_rc%ess)
       if(LVT_rc%time.ge.LVT_rc%etime)then
          LVT_rc%endtime=1
          write(unit=LVT_logunit,fmt=*) '[INFO] LVT cycle completed'
       endif
    endif
    
    write(unit=LVT_logunit,fmt=24)' [INFO] LVT cycle time: ',LVT_rc%mo,'/',LVT_rc%da,'/', & 
         LVT_rc%yr,LVT_rc%hr,':',LVT_rc%mn,':',LVT_rc%ss
    
  end subroutine LVT_advance_timestep

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_get_step_size
! \label{LVT_get_step_size} 
!
! !INTERFACE:
  function LVT_get_step_size(LVT_rc)

    implicit none
! !ARGUMENTS: 
    type(lvtrcdec) :: LVT_rc
    integer :: LVT_get_step_size
!
! !DESCRIPTION:
!
! Return the timestep used by the clock in the time manager. 
!
!  \begin{description}
!   \item [LVT\_rc]
!     instance of the {\tt lvt\_module}
!   \item[get\_step\_size]
!     timestep value
!  \end{description}
!
!EOP

    LVT_get_step_size = LVT_rc%ts

  end function LVT_get_step_size

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_get_nstep
! \label{LVT_get_nstep} 
! 
! !INTERFACE:
  function LVT_get_nstep(LVT_rc)

    implicit none
! !ARGUMENTS: 
    type(lvtrcdec) :: LVT_rc
    integer :: LVT_get_nstep
!
! !DESCRIPTION: 
!
! Return the timestep number for each nest. 
!
!  \begin{description}
!   \item [LVT\_rc]
!     instance of the {\tt lvt\_module}
!   \item [n]
!     index of the nest
!   \item[get\_nstep]
!     timestep number
!  \end{description}
!
!EOP
    LVT_get_nstep = LVT_rc%tscount
    
  end function LVT_get_nstep

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_get_curr_day
! \label{LVT_get_curr_day}
!
! !INTERFACE:
  subroutine LVT_get_curr_date(LVT_rc, yr, mon, day, tod, offset)

    implicit none
   
! !ARGUMENTS: 
    type(lvtrcdec) :: LVT_rc
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
!   \item [LVT\_rc]
!     instance of the {\tt lvt\_module}
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
    yr = LVT_rc%yr
    mon = LVT_rc%mo
    day = LVT_rc%da
    tod = LVT_rc%ss+LVT_rc%mn*60+LVT_rc%hr*3600
    
  end subroutine LVT_get_curr_date


!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: get_julhr
! \label{get_julhr}
! 
! !INTERFACE:
  subroutine LVT_get_julhr(yr,mo,da,hr,mn,ss,julhr)

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

  end subroutine LVT_get_julhr

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_get_julss
! \label{LVT_get_julss}
! 
! !INTERFACE:
  subroutine LVT_get_julss(yr,mo,da,hr,mn,ss,julss)

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

  end subroutine LVT_get_julss


!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_get_curr_calday
! \label{LVT_get_curr_calday}
! 
! !INTERFACE:
  function LVT_get_curr_calday(LVT_rc,offset)

   implicit none
   
! !ARGUMENTS: 
   type(lvtrcdec) :: LVT_rc
   integer, optional, intent(in) :: offset
   real :: LVT_get_curr_calday
!
! !DESCRIPTION:
!
! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.
! 
! The arguments are: 
!  \begin{description}
!   \item [LVT\_rc]
!     instance of the {\tt lvt\_module}
!   \item [offset]
!     offset from current time in seconds, positive for future times, 
!     negative for previous times. 
!   \item [get\_curr\_calday]
!     calendar day value
!  \end{description}
!
!EOP
   integer :: days(12)
   integer :: i
   data days /31,28,31,30,31,30,31,31,30,31,30,31/
   
   LVT_get_curr_calday = 0
   if((mod(LVT_rc%yr,4) .eq. 0 .and. mod(LVT_rc%yr, 100).ne.0) &!leap year
        .or.(mod(LVT_rc%yr,400) .eq.0)) then 
      days(2) = 29
   else 
      days(2) = 28
   endif
   if(LVT_rc%mo .ne. 1) then 
      do i=1,LVT_rc%mo-1
         LVT_get_curr_calday = LVT_get_curr_calday+days(i)
      enddo
   endif
   LVT_get_curr_calday = LVT_get_curr_calday + real(LVT_rc%da) + real(LVT_rc%hr)/24 + &
        real(LVT_rc%mn)/(24*60) + real(LVT_rc%ss)/(24*60*60)
   
   if (present(offset)) then
      if (offset > 0) then
         LVT_get_curr_calday = LVT_get_curr_calday + real(offset)/(24*60*60)
      else if (offset < 0) then
         LVT_get_curr_calday = LVT_get_curr_calday - real(offset)/(24*60*60)
      endif
   endif
   
 end function LVT_get_curr_calday


!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_is_last_step
! \label{LVT_is_last_step}
!
! !INTERFACE:
function LVT_is_last_step(LVT_rc)
   implicit none
! !ARGUMENTS: 
   type(lvtrcdec) :: LVT_rc
   logical :: LVT_is_last_step
! 
! !DESCRIPTION:
!
! Function returns true on last timestep.
!
!  \begin{description}
!   \item [LVT\_rc]
!     instance of the {\tt lvt\_module}
!   \item [is\_last\_step]
!     result of the function 
!  \end{description}
!EOP

   LVT_is_last_step = .false.
   if(LVT_rc%time .ge. LVT_rc%etime1) then
      LVT_is_last_step = .true.
   endif
   
 end function LVT_is_last_step


!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_setAlarm
! \label{LVT_setAlarm}
!
! !INTERFACE:
  subroutine LVT_setAlarm(LVT_rc, intervalType, alarm)
! !USES: 

! !ARGUMENTS: 
    type(lvtrcdec)    :: LVT_rc
    integer           :: intervalType
    type(ESMF_Alarm)  :: alarm
! 
! !DESCRIPTION: 
!  This routine initializes the ESMF alarms based on the 
!  specified intervaltype (weekly, monthly, quarterly). 
!EOP
    type(ESMF_Time)         :: alarmTime
    type(ESMF_Time)         :: refTime
    type(ESMF_Time)         :: currTime
    type(ESMF_Time)         :: qt1,qt2,qt3,qt4
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_TimeInterval) :: deltaT
    integer                 :: year, month
    integer                 :: refjuld, sjuld
    integer                 :: rc
        
    if(intervalType.eq.2) then !weekly
       call ESMF_AlarmGet(alarm,refTime=refTime,rc=rc)

       call ESMF_TimeGet(refTime,d=refjuld,rc=rc)
       call ESMF_ClockGet(LVT_clock,currTime=currTime, rc=rc)
       call ESMF_TimeGet(currTime,d=sjuld,rc=rc)
       call ESMF_TimeIntervalSet(deltaT,d=(7-mod((sjuld-refjuld),7)),rc=rc)
       alarmTime = currTime + deltaT
       call ESMF_TimeIntervalSet(alarmInterval,d=7,rc=rc)
       call ESMF_AlarmSet(alarm, ringTime=alarmTime, &
            ringInterval=alarmInterval,rc=rc)
       call LVT_verify(rc,'LVT_setalarm:alarmCreate')
    elseif(intervalType.eq.3) then 
!monthly alarm
       if(LVT_rc%sda.ge.16) then        
          if(LVT_rc%smo.eq.12) then 
             year = LVT_rc%syr + 1
             month = 1
          else
             year = LVT_rc%syr
             month = LVT_rc%smo +1
          endif
       else
          year = LVT_rc%syr
          month = LVT_rc%smo
       endif
       
       call ESMF_TimeSet(alarmTime, yy=year, mm=month,&
            dd=15, calendar = LVT_calendar, rc=rc)
       call ESMF_TimeIntervalSet(alarmInterval,mm=1,rc=rc)
       call LVT_verify(rc,'LVT_setalarm:timeintervalset')
       
       call ESMF_AlarmSet(alarm, ringTime=alarmTime, &
            ringInterval=alarmInterval,rc=rc)
       call LVT_verify(rc,'LVT_setalarm:alarmCreate')
    elseif(intervalType.eq.4) then
       call ESMF_ClockGet(LVT_clock,currTime=currTime, rc=rc)
       call ESMF_TimeGet(currTime,yy=year,mm=month,rc=rc)
             
       call ESMF_TimeSet(qt1,yy=year,mm=1,dd=30,&
            calendar=LVT_calendar,rc=rc)
       call ESMF_TimeSet(qt2,yy=year,mm=4,dd=30,&
            calendar=LVT_calendar,rc=rc)
       call ESMF_TimeSet(qt3,yy=year,mm=7,dd=30,&
            calendar=LVT_calendar,rc=rc)
       call ESMF_TimeSet(qt4,yy=year,mm=10,dd=30,&
            calendar=LVT_calendar,rc=rc)

       call ESMF_TimeIntervalSet(alarmInterval,mm=3,rc=rc)       
       if(currTime <= qt1 ) then 
          alarmTime = qt1 
       elseif(currTime > qt1 .and. currTime <= qt2) then 
          alarmTime = qt2
       elseif(currTime > qt2 .and. currTime <=qt3) then 
          alarmTime = qt3
       elseif(currTime > qt3 .and. currTime <=qt4) then 
          alarmTime = qt4
       else
          alarmTime = qt4 + alarmInterval
       endif

       call ESMF_AlarmSet(alarm, ringTime=alarmTime, &
            ringInterval=alarmInterval,rc=rc)
       call LVT_verify(rc,'LVT_setalarm:alarmCreate')
    endif

  end subroutine LVT_setAlarm
!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_setMonthlyAlarm
! \label{LVT_setMonthlyAlarm}
!
! !INTERFACE:
  subroutine LVT_setMonthlyAlarm(alarmTime)

    implicit none
! !ARGUMENTS: 
    real*8 :: alarmTime
!
! !DESCRIPTION:
!  Initializes the monthly alarm
! 
!  The arguments are: 
!  \begin{description}
!   \item[alarmTime]
!   time of the monthly alarm
!  \end{description}
! 
!EOP
    alarmTime = 0 

  end subroutine LVT_setMonthlyAlarm


!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! !ROUTINE: isMonthlyAlarmRinging1
! \label{isMonthlyAlarmRinging1}
!
! !INTERFACE:
! !Private name: call using isMonthlyAlarmRinging
  subroutine isMonthlyAlarmRinging1(LVT_rc, alarmTime, interval, &
       midmonth, ringFlag)

    implicit none
! !ARGUMENTS: 
    type(lvtrcdec), intent(in) :: LVT_rc
    real*8, intent(inout) :: alarmTime
    logical, intent(inout):: ringFlag
    integer, intent(in) :: interval
    logical, intent(in) :: midmonth
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
!   \item[LVT\_rc]
!    instance of the {\tt lvt\_module}
!   \item[alarmTime]
!    the elapsed alarm time
!   \item[interval]
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
!   \item[date2time](\ref{date2time}) \newline
!    convert date to a flotating point format
! \end{description}
!EOP
    integer :: yr2,mo2,yr
    integer :: numi, zeroi,doy1,doy2
    real :: gmt2,gmt1
    real*8 :: time2,time
    real*8 :: jan31,apr30,jul31,oct31 ! dates of quarterly intervals
    integer :: janda,janmo          ! january 31 
    integer :: aprda,aprmo          ! april 30
    integer :: julda,julmo          ! july 31
    integer :: octda,octmo          ! october 31 
    real*8 :: quartTime 
    
    zeroi = 0 
    numi = 0 
    ringFlag = .false.
    if(interval.eq.0 ) then 
       if(alarmTime.ne.0) then 
          ringFlag = .false. 
       else 
          ringFlag = .true. 
          alarmTime = LVT_rc%udef
       endif
    elseif(interval.eq.1) then !monthly
       if(midmonth) then 
          if(LVT_rc%da.lt.16) then 
             mo2=LVT_rc%mo
             yr2=LVT_rc%yr
          else
             mo2=LVT_rc%mo+1
             yr2=LVT_rc%yr
             if(mo2.eq.13)then
                mo2=1
                yr2=LVT_rc%yr+1
             endif
          endif
          
          call LVT_date2time(time2,doy2,gmt2,yr2,mo2,&
               numi,zeroi,zeroi,zeroi)          
          if(time2.gt.alarmTime) then 
             ringFlag = .true.
             alarmTime = time2
          endif
       else ! check for the end of the month          
          mo2 = LVT_rc%mo + 1
          yr2 = LVT_rc%yr
          if(mo2.eq.13) then 
             mo2 = 1
             yr2 = LVT_rc%yr + 1
          endif
          call LVT_date2time(time2,doy2,gmt2,yr2,mo2,&
               numi,zeroi,zeroi,zeroi)
          if(time2.gt.alarmTime) then 
             ringFlag = .true. 
             alarmTime = time2
          endif          
       endif
    elseif(interval.eq.3) then !quarterly
       zeroi = 0 
       time=LVT_rc%time
       yr=LVT_rc%yr
       janda=31
       janmo=01
       call LVT_date2time(jan31,doy1,gmt1,yr,janmo,&
            janda,zeroi,zeroi,zeroi)
       aprda=30
       aprmo=04
       call LVT_date2time(apr30,doy1,gmt1,yr,aprmo,&
            aprda,zeroi,zeroi,zeroi)
       julda=31
       julmo=07
       call LVT_date2time(jul31,doy1,gmt1,yr,julmo,&
            julda,zeroi,zeroi,zeroi)
       octda=31
       octmo=10
       call LVT_date2time(oct31,doy1,gmt1,yr,octmo,&
            octda,zeroi,zeroi,zeroi)
       if ( time.ge.jan31 .and. time.le.apr30 ) then
          quartTime = 1
       elseif ( time.ge.apr30 .and. time.le.jul31 ) then
          quartTime = 2
       elseif ( time.ge.jul31 .and. time.le.oct31 ) then
          quartTime = 3
       elseif ( time.ge.oct31 ) then
          quartTime = 4
       elseif ( time.lt.jan31) then
          quartTime = 5
       endif
       
       if(alarmtime .ne. quartTime) then 
          ringFlag = .true.
          alarmTime = quartTime
       endif
    endif
  end subroutine isMonthlyAlarmRinging1

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! !ROUTINE: isMonthlyAlarmRinging2
! \label{isMonthlyAlarmRinging2}
!
! !INTERFACE:
! !Private name: call using isMonthlyAlarmRinging
  subroutine isMonthlyAlarmRinging2(LVT_rc, alarmTime, interval, ringFlag)
!EOP
    implicit none
! !ARGUMENTS: 
    type(lvtrcdec) :: LVT_rc
    real*8 :: alarmTime
    integer :: interval
    logical :: ringFlag
!
! !DESCRIPTION:
!  checks if the monthly alarm is ringing. The function returns 
!  true when the elapsed alarm time is greater than the number of 
!  months specified in the alarm's interval.
!
!  The arguments are: 
!  \begin{description}
!   \item[LVT\_rc]
!    instance of the {\tt lvt\_module}
!   \item[alarmTime]
!    the elapsed alarm time
!   \item[interval]
!    the alarm's frequency
!   \item[ringflag]
!    flag indicating the status of the call
!  \end{description}
!
! The calling sequence is:  
! \begin{description}
!   \item[date2time](\ref{date2time}) \newline
!    convert date to a flotating point format
! \end{description}
!   
!EOP
    integer :: yr2,mo2,yr
    integer :: numi, zeroi,doy1,doy2
    real :: gmt2,gmt1

    real*8 :: time2,time
    real*8 :: jan31,apr30,jul31,oct31 ! dates of quarterly intervals
    integer :: janda,janmo          ! january 31 
    integer :: aprda,aprmo          ! april 30
    integer :: julda,julmo          ! july 31
    integer :: octda,octmo          ! october 31 

    real*8 :: quartTime 

    numi = 0 
    zeroi = 0   
    ringFlag = .false.
    if(interval.eq.0 ) then 
       if(alarmTime.ne.0) then 
          ringFlag = .false. 
       else 
          ringFlag = .true. 
          alarmTime = LVT_rc%udef
       endif
    elseif(interval.eq.1) then !monthly
       if(LVT_rc%da.lt.16) then 
          mo2=LVT_rc%mo
          yr2=LVT_rc%yr
       else
          mo2=LVT_rc%mo+1
          yr2=LVT_rc%yr
          if(mo2.eq.13)then
             mo2=1
             yr2=LVT_rc%yr+1
          endif
       endif
       call LVT_date2time(time2,doy2,gmt2,yr2,mo2,&
            numi,zeroi,zeroi,zeroi)          
       if(time2.gt.alarmTime) then 
          ringFlag = .true.
          alarmTime = time2
       endif
    elseif(interval.eq.3) then !quarterly
       zeroi = 0 
       time=LVT_rc%time
       yr=LVT_rc%yr
       janda=31
       janmo=01
       call LVT_date2time(jan31,doy1,gmt1,yr,janmo,&
            janda,zeroi,zeroi,zeroi)
       aprda=30
       aprmo=04
       call LVT_date2time(apr30,doy1,gmt1,yr,aprmo,&
            aprda,zeroi,zeroi,zeroi)
       julda=31
       julmo=07
       call LVT_date2time(jul31,doy1,gmt1,yr,julmo,&
            julda,zeroi,zeroi,zeroi)
       octda=31
       octmo=10
       call LVT_date2time(oct31,doy1,gmt1,yr,octmo,&
            octda,zeroi,zeroi,zeroi)
       if ( time.ge.jan31 .and. time.le.apr30 ) then
          quartTime = 1
       elseif ( time.ge.apr30 .and. time.le.jul31 ) then
          quartTime = 2
       elseif ( time.ge.jul31 .and. time.le.oct31 ) then
          quartTime = 3
       elseif ( time.ge.oct31 ) then
          quartTime = 4
       elseif ( time.lt.jan31) then
          quartTime = 5
       endif
       
       if(alarmtime .ne. quartTime) then 
          ringFlag = .true.
          alarmTime = quartTime
       endif
    endif

  end subroutine isMonthlyAlarmRinging2

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_setHourlyAlarm
! \label{LVT_setHourlyAlarm}
!
! !INTERFACE:
  subroutine LVT_setHourlyAlarm(alarmTime)

    implicit none
! !ARGUMENTS: 
    real*8 :: alarmTime
!
! !DESCRIPTION:
!  Initializes the hourly alarm
!
!  The arguments are: 
!  \begin{description}
!   \item[alarmTime]
!   time of the monthly alarm
!  \end{description}
!
!EOP
    alarmTime = 0 

  end subroutine LVT_setHourlyAlarm

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_isHourlyAlarmRinging
! \label{LVT_isHourlyAlarmRinging}
!
! !INTERFACE:
  subroutine LVT_isHourlyAlarmRinging(LVT_rc, alarmTime, interval, ringFlag)

    implicit none
! !ARGUMENTS: 
    type(lvtrcdec), intent(in) :: LVT_rc    
    real*8, intent(inout) :: alarmTime
    integer, intent(in)   :: interval
    logical, intent(inout) :: ringFlag
!
! !DESCRIPTION:
!  checks if the hourly alarm is ringing. The function returns 
!  true when the elapsed alarm time is greater than the number of 
!  hours specified in the alarm's interval.
!
!  The arguments are: 
!  \begin{description}
!   \item[LVT\_rc]
!    instance of the {\tt lvt\_module}
!   \item[alarmTime]
!    the elapsed alarm time
!   \item[interval]
!    the alarm's frequency in hours
!   \item[ringflag]
!    flag indicating the status of the call
!  \end{description}
!   
! The calling sequence is:  
! \begin{description}
!   \item[LVT\_tick] (\ref{LVT_tick}) \newline
!    advance the time by the specified increment
! \end{description}
!EOP

    integer :: yr1,mo1,da1,hr1,mn1,ss1,ts1
    integer :: yr2,mo2,da2,hr2,mn2,ss2,ts2
    real*8  :: time2,timenow
    integer :: doy1,doy2
    real    :: gmt1,gmt2
    
    ringFlag = .false. 

    yr1 = LVT_rc%yr  !current time
    mo1 = LVT_rc%mo
    da1 = LVT_rc%da
    hr1 = LVT_rc%hr
    mn1 = LVT_rc%mn
    ss1 = 0
    ts1 = 0
    call LVT_tick( timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
    
    yr2 = LVT_rc%yr  !next assimilation/forecast hour
    mo2 = LVT_rc%mo
    da2 = LVT_rc%da
    hr2 = interval*((LVT_rc%hr/interval))
    mn2 = 0
    ss2 = 0
    ts2 = interval*60*60
    call LVT_tick( time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

    if ( timenow .ge.alarmTime ) then
       ringFlag = .true. 
       alarmTime = time2
    endif
  end subroutine LVT_isHourlyAlarmRinging

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! !ROUTINE: computeTimeBookEnds
! \label{computeTimeBookEnds}
!
! !INTERFACE:
  subroutine LVT_computeTimeBookEnds(LVT_rc,interval,time1,time2)
    implicit none
! !ARGUMENTS: 
    type(lvtrcdec), intent(in) :: LVT_rc
    integer, intent(in) :: interval
    real*8, intent(inout) :: time1,time2
!
! !DESCRIPTION:
!  compute the time end points based on the interval. If the interval
!  is 3 hours, and the current time is 1:30Z, time1 will be set to 
!  1Z and time2 will be set to 4Z. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[LVT\_rc]
!    instance of the {\tt lvt\_module}
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
!   \item[LVT\_tick] (\ref{LVT_tick}) \newline
!    advance the time by the specified increment
! \end{description}
!EOP
    integer :: yr1,mo1,da1,hr1,mn1,ss1,ts1
    integer :: yr2,mo2,da2,hr2,mn2,ss2,ts2
    real    :: gmt1,gmt2
    integer :: doy1,doy2

    yr1=LVT_rc%yr    !Previous Hour
    mo1=LVT_rc%mo
    da1=LVT_rc%da
    hr1=interval*((LVT_rc%hr)/interval)
    mn1=0
    ss1=0
    ts1=0
    call LVT_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
    
    yr2=LVT_rc%yr    !Next Hour
    mo2=LVT_rc%mo
    da2=LVT_rc%da
    hr2=interval*((LVT_rc%hr)/interval)
    mn2=0
    ss2=0
    ts2=interval*60*60
    call LVT_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

  end subroutine LVT_computeTimeBookEnds

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! !ROUTINE: LVT_computeTemporalWeights1
! \label{LVT_computeTemporalWeights1}
!
! !INTERFACE:
  subroutine computeTemporalWeights1(LVT_rc,alarm,intervalType,&
       time1, time2,wt1,wt2)

    implicit none
! !ARGUMENTS: 
    type(lvtrcdec)   :: LVT_rc    
    type(ESMF_Alarm) :: alarm
    integer          :: intervalType
    type(ESMF_Time)  :: time1
    type(ESMF_Time)  :: time2
    real    :: wt1,wt2
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
!   \item[LVT\_rc]
!    instance of the {\tt lvt\_module}
!   \item[interval]
!    interval of the climatology
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
!   \item[date2time](\ref{date2time}) \newline
!    converts date to a floating point format
! \end{description}
!EOP
    type(ESMF_Time)         :: currTime, ringTime
    type(ESMF_TimeInterval) :: deltaT 
    integer                 :: yr,mo,da,hr
    integer                 :: rc

    call ESMF_ClockGet(LVT_clock, currTime=currTime, rc=rc)
    call ESMF_AlarmGet(alarm,ringTime=ringTime,rc=rc)

    if(intervalType.eq.2) then  !weekly
       call ESMF_TimeIntervalSet(deltaT,d=7,rc=rc)
    elseif(intervalType.eq.3) then  !monthly
       call ESMF_TimeIntervalSet(deltaT,mm=1,rc=rc)
    elseif(intervalType.eq.4) then 
       call ESMF_TimeIntervalSet(deltaT,mm=3,rc=rc)
    endif

    call ESMF_TimeGet(currTime,yy=yr,mm=mo,dd=da,h=hr)
    call ESMF_TimeGet(ringTime,yy=yr,mm=mo,dd=da,h=hr)

    if(currTime >=ringTime) then
       time1 = ringTime
       time2 = time1 + deltaT         
    else !initial
       time1 = ringTime -deltaT
       time2 = ringTime
    endif
    call ESMF_TimeGet(time1,yy=yr,mm=mo,dd=da,h=hr)
    call ESMF_TimeGet(time2,yy=yr,mm=mo,dd=da,h=hr)
    
    wt1 = (time2-currTime)/(time2-time1)
    wt2 = (currTime-time1)/(time2-time1)

  end subroutine computeTemporalWeights1


!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! !ROUTINE: computeTemporalWeights2
! \label{computeTemporalWeights2}
!
! !INTERFACE:
  subroutine computeTemporalWeights2(LVT_rc,intervalType,mo1,mo2,wt1,wt2)

    implicit none
! !ARGUMENTS: 
    type(lvtrcdec) :: LVT_rc    
    integer :: intervalType
    integer :: mo1,mo2
    real    :: wt1,wt2
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
!   \item[LVT\_rc]
!    instance of the {\tt lvt\_module}
!   \item[interval]
!    interval of the climatology
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
!   \item[date2time](\ref{date2time}) \newline
!    converts date to a floating point format
! \end{description}
!EOP
    integer   :: juld
    integer :: julm(13)
    real    :: day1
    real    :: day2
    real    :: rday
    real    :: s1,s2

    data julm /  0,  31,  59,  90, 120, 151, &
         181 , 212, 243, 273, 304, 334, 365/

    if(intervalType.eq.3) then  !monthly
       mo2 = LVT_rc%mo
       if(LVT_rc%da.gt.15) mo2 = mo2+1
       if(mo2.eq.1) mo2 = 13
       mo1 = mo2-1

       juld = julm(LVT_rc%mo)+LVT_rc%da
       if((LVT_rc%mo.eq.3).and.(mod(LVT_rc%yr,4).eq.0)) then 
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
      
       if ((mo2 .eq. 3) .and. (mod(LVT_rc%yr, 4) .eq. 0)) then
          day2 = day2 + 1
       end if
       day1  = float(julm(mo1) + 15)
       rday  = float(juld)
       
       wt1 = (day2 - rday) / (day2 - day1)
       wt2 = (rday - day1) / (day2 - day1)       
       
       if(mo2.eq.13) mo2 = 1
       if(mo1.eq.13) mo1 = 1

    elseif(intervalType.eq.4) then 
       juld = julm(LVT_rc%mo) + LVT_rc%da  

!     ------------------------------------------------------------------
!     between winter and autumn.
!     ------------------------------------------------------------------

       if (juld .le. 32) then
          
          s1   = float(32 - juld)
          s2   = float(juld + 30)
          wt1 = s2 / (s1+s2)
          wt2 = s1 / (s1+s2)
          mo1 = 4
          mo2 = 3

!     ------------------------------------------------------------------
!     between winter and spring.
!     ------------------------------------------------------------------

       else if ( (juld .le. 121) .and. (juld .gt. 32)) then
          
          s1   = float(121 - juld)
          s2   = float(juld - 32)
          wt1 = s2 / (s1+s2)
          wt2 = s1 / (s1+s2)
          mo1 = 1
          mo2 = 4
          
!     ------------------------------------------------------------------
!     between summer and spring.
!     ------------------------------------------------------------------

       else if ( (juld .le. 213) .and. (juld .gt. 121) ) then
          
          s1   = float(213 - juld)
          s2   = float(juld - 121)
          wt1 = s2 / (s1+s2)
          wt2 = s1 / (s1+s2)
          mo1 = 2
          mo2 = 1

!     ------------------------------------------------------------------
!     between autumn and summer.
!     ------------------------------------------------------------------

       else if ( (juld .le. 305) .and. (juld .gt. 213) ) then
          
          s1   =  float(305 - juld)
          s2   =  float(juld - 213)
          wt1 =  s2 / (s1+s2)
          wt2 =  s1 / (s1+s2)
          mo1 = 3
          mo2 = 2
!     ------------------------------------------------------------------
!     between winter and autumn..
!     ------------------------------------------------------------------

       else
          
          s1 = float(365 - juld + 32)
          s2 = float(juld - 305)
          wt1 = s2 / (s1+s2)
          wt2 = s1 / (s1+s2)
          mo1 = 4
          mo2 = 3          
       end if
    endif
  end subroutine ComputeTemporalWeights2

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! !ROUTINE: date2time
! \label{date2time} 
! 
! !INTERFACE:
  subroutine LVT_date2time(time,doy,gmt,yr,mo,da,hr,mn,ss)

    implicit none
! !ARGUMENTS:
    integer :: yr,mo,da,hr,mn,ss, doy
    real*8  :: time
    real    :: gmt
! !DESCRIPTION:
! 
!  determines the time, time in GMT, and the day of the year
!  based on the value of year, month, day of month, hour of 
!  the day, minute and second. This method is the inverse of 
!  time2date
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
!   lvt time
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
  end subroutine LVT_date2time

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ROUTINE: time2date
! \label{time2date}
! 
! !INTERFACE: 
  subroutine LVT_time2date(time,doy,gmt,yr,mo,da,hr,mn)

    implicit none
! !ARGUMENTS: 
    integer :: yr,mo,da,hr,mn,ss,doy
    real*8  :: time
    real    :: gmt
! !DESCRIPTION:
! 
!  determines the value of the year, month, day of month, hour of 
!  the day, minute and second based on the specified time. This
!  method is the inverse of date2time
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
!   lvt time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
    real*8 :: tmp
    integer :: yrdays,days(13)
    data days /31,28,31,30,31,30,31,31,30,31,30,31,30/
    
    yr  = dint(time)
    tmp =     (time) 
    
    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
         .or.(mod(yr,400).eq.0))then             !correct for y2k
       yrdays=366                  
    else
       yrdays=365
    endif
    if (yrdays.eq.366) then
       days(2)=29
    else
       days(2)=28
    endif
    
    doy  = dint((tmp-yr)*float(yrdays))+1 
    tmp =      ((tmp-yr)*float(yrdays))+1 
    hr  = nint((tmp-doy)*24.d0) 
    tmp =     ((tmp-doy)*24.d0) 
    
    mn  = dint((tmp-hr)*60.d0) 
    tmp =     ((tmp-hr)*60.d0) 
    
    ss  = dint((tmp-mn)*60.d0) 
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
  end subroutine LVT_time2date

  subroutine LVT_seconds2time(secs,da, hr, mn, ss)
    
    integer,     intent(in)   :: secs
    integer,     intent(out)  :: da
    integer,     intent(out)  :: hr
    integer,     intent(out)  :: mn
    integer,     intent(out)  :: ss
    
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
       
  end subroutine LVT_seconds2time

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ROUTINE: LVT_tick
! \label{LVT_tick}
! 
! !INTERFACE:
  subroutine LVT_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,ts)
    implicit none
! !ARGUMENTS:
    real*8  :: time
    integer :: yr,mo,da,hr,mn,ss,ts,doy
    real    :: gmt
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
!   lvt time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
    integer days(13)
    integer prvmo   !previous month
    
    data days/31,28,31,30,31,30,31,31,30,31,30,31,31/

143 format(a1,' yr',i6,' mo',i5,' dy',i5,' hr',i5, & 
         ' mn',i6,' ss',i8,' ts',i8)
    ss=ss+ts
    do while(ss.gt.59)
       ss=ss-60
       mn=mn+1
    enddo
    do while(ss.lt.0)
       ss=ss+60
       mn=mn-1
    enddo
    do while(mn.gt.59)
       mn=mn-60
       hr=hr+1
    enddo
    
    do while(mn.lt.0)
       mn=mn+60
       hr=hr-1
    enddo
    do while(hr.gt.23)
       hr=hr-24
       da=da+1
    enddo
    
    do while(hr.lt.0)
       hr=hr+24
       da=da-1
    enddo
    
    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &    !correct for leap year
         .or.(mod(yr,400).eq.0))then               !correct for y2k
       days(2)=29                  
    else
       days(2)=28
    endif
    
    do while(da.gt.days(mo))
       da=da-days(mo)
       mo=mo+1
    enddo
    
    do while(da.lt.1)
       
       prvmo=mo-1
       if(mo.eq.1) prvmo=12
       
       da=da+days(prvmo)
       
       if(prvmo.eq.12) then
          mo=prvmo
          yr=yr-1
       else
          mo=prvmo
       endif
    enddo
    do while(mo.gt.12)
       mo=mo-12
       yr=yr+1
    enddo

    do while(mo.lt.1)
       mo=mo+12
       yr=yr-1
    enddo
    call LVT_date2time(time,doy,gmt,yr,mo,da,hr,mn,ss)
    return

  end subroutine LVT_tick


!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! !ROUTINE: julhr_date
! \label{julhr_date}
!
! !REVISION HISTORY:
!     15 oct 1998  initial version.........................mr moore/dnxm
!     10 aug 1999  ported to ibm sp2.  added intent attributes to
!                  arguments...............................mr gayno/dnxm
!     29 oct 2005  Sujay Kumar, Adopted in LVT
!
! !INTERFACE:    
subroutine LVT_julhr_date( julhr, yyyy,mm,dd,hh)
! !USES:


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
!     call utility routine LVT\_tmjul4 to convert from julian hours
!     to year, month, day, hour. \newline
!     perform several checks to ensure date information passed
!     back from LVT\_tmjul4 is valid. \newline
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
!   \item[LVT\_tmjul4] (\ref{LVT_tmjul4}) \newline
!    convert julian hour to hour, day, month, and year
!   \item[lvt\_abort] (\ref{LVT_abort}) \newline
!    abort the code in case of error
! \end{description}
!
!EOP
  character*100                  :: message ( 20 )  

!     ------------------------------------------------------------------
!     executable code begins here... use LVT_tmjul4 to convert julhr to 
!     hour, day, month and year
!     ------------------------------------------------------------------

  call LVT_tmjul4( hh, dd, mm, yyyy, julhr )
 
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
     
     call lvt_abort( message )
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
        
        call lvt_abort( message )
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
     
     call lvt_abort( message )
  endif
  
  return
end subroutine LVT_julhr_date

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! !ROUTINE: LVT_tmjul4
! \label{LVT_tmjul4}
!
! !REVISION HISTORY:
!     18 mar 98 initial version.........................sra milburn/dnxm
!     08 jul 99 fixed error which initialzed done flag in a data
!               data statement.  variables must be initialized
!               using an assignment statement.  ported to ibm sp2.......
!               ...........................................mr gayno/dnxm
!     29 oct 2005  Sujay Kumar, Adopted in LVT
!
! !INTERFACE:    
subroutine LVT_tmjul4( hour, day, month, year, julhr ) 

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
!     - determine the current zulu hour  \newline
!     - determine the total number of elapsed days  \newline
!     - count forward to the current day/month/year  \newline
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
  
end subroutine LVT_tmjul4

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ROUTINE: LVT_localtime2gmt
!  \label{LVT_localtime2gmt}
! 
! !DESCRIPTION: 
! 
! Calculates the local time based on GMT
! 
! !INTERFACE:
subroutine LVT_localtime2gmt (gmt,lon,lhour,zone)
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
end subroutine LVT_localtime2gmt

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ROUTINE: LVT_localtime
!  \label{LVT_localtime}
! 
! !DESCRIPTION: 
! 
! Calculates the local time based on GMT
! 
! !INTERFACE:
subroutine LVT_localtime (gmt,lon,lhour,zone)
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
end subroutine LVT_localtime

subroutine LVT_doy2moda(yr, doy, mo, da)

  implicit none
! !ARGUMENTS: 
    integer :: yr,doy, mo, da

! !DESCRIPTION:
!   converts day of the year to month and day
! 
!  The arguments are: 
!  \begin{description}
!  \item[yr]
!    year
!  \item[doy]
!   day of the year
!  \item[mo]
!   month
!  \item[da]
!   day of the month
!  \end{description}
!EOP
    logical :: leapyear
    integer :: i
    integer :: days1(12), days2(12), cdays1(12), cdays2(12)

    data days1  /31,28,31,30,31,30,31,31,30,31,30,31/
    data days2  /31,29,31,30,31,30,31,31,30,31,30,31/
    data cdays1 /31,59,90,120,151,181,212,243,273,304,334,365/
    data cdays2 /31,60,91,121,152,182,213,244,274,305,335,366/

   
    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
         .or.(mod(yr,400).eq.0))then             !correct for y2k
       leapyear = .true.             
    else
       leapyear = .false. 
    endif
    
    do i=1,12
       if(leapyear) then 
          if(doy.le.cdays2(i)) then 
             mo = i
             if(i.ne.1) then 
                da = doy-cdays2(i-1)
             else
                da = doy
             endif
             exit;
          endif
       else
          if(doy.le.cdays1(i)) then 
             mo = i
             if(i.ne.1) then 
                da = doy-cdays1(i-1)
             else
                da = doy
             endif
             exit;
          endif
       endif
    enddo
       
end subroutine LVT_doy2moda

!BOP
! !ROUTINE: LVT_parseTimeString
! \label{LVT_parseTimeString}
!
! !INTERFACE: 
  subroutine LVT_parseTimeString(inputstr,time_value)
    
    implicit none
! !ARGUMENTS:     
    character(len=*)          :: inputstr
    integer                   :: time_value
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
    integer            :: iloc, ios

    iloc  = len(trim(inputstr))-1
    suffix = inputstr(iloc:len(inputstr))
    
    read(inputstr(1:iloc-1),*,iostat=ios) time_input
    if(ios.ne.0) then 
       write(LVT_logunit,*) '[ERR] '
       write(LVT_logunit,*) '[ERR] The time specification ',trim(inputstr), &
            ' is incorrect'
       write(LVT_logunit,*) '[ERR] Please use a two letter suffix for time specification'
       write(LVT_logunit,*) '[ERR] '
       call LVT_endrun()
    endif

    if(suffix.eq."ss") then 
       time_value = nint(time_input)
    elseif(suffix.eq."mn") then 
       time_value = nint(time_input*60.0)
    elseif(suffix.eq."hr") then
       time_value = nint(time_input*3600.0)
    elseif(suffix.eq."da") then 
       time_value = nint(time_input*86400.0)
    elseif(suffix.eq."mo") then 
       time_value = nint(time_input*2592000.0)
    elseif(suffix.eq."yr") then 
       time_value = nint(time_input*31536000.0)
    endif
    
  end subroutine LVT_parseTimeString

!BOP
! 
! !ROUTINE: LVT_computeTimeSpan
! \label{LVT_computeTimeSpan}
! 
! !INTERFACE: 
  subroutine LVT_computeTimeSpan(nsize, ts)
! !ARGUMENTS:     
    integer             :: nsize
    integer             :: ts
!
! !DESCRIPTION: 
!  This subroutine computes the time span of the LVT analysis
!  based on the start, stop times and time averaging interval. 
! 
!  The arguments are:
!  \begin{description}
!   \item[nyears] 
!     number of years
!  \end{description}
!EOP    
    
    type(ESMF_Time)         :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: status

 
    call ESMF_TimeIntervalSet(timeStep, s = ts, &
         rc=status)
    call LVT_verify(status, &
         'ESMF_TimeIntervalSet failed in LVT_computeTimeSpan')

    call ESMF_ClockGet(LVT_clock, startTime = startTime, &
         stopTime = stopTime, rc=status)
    call LVT_verify(status, 'ESMF_ClockGet failed in LVT_computeTimeSpan')

    nsize = nint((stopTime - startTime)/timestep) + 1

  end subroutine LVT_computeTimeSpan

!BOP
! 
! !ROUTINE: LVT_computeTimeSpanInYears
! \label{LVT_computeTimeSpanInYears}
! 
! !INTERFACE: 
  subroutine LVT_computeTimeSpanInYears(nyears)
! !ARGUMENTS:     
    integer             :: nyears

!
! !DESCRIPTION: 
!  This subroutine computes the time span of the LVT analysis
!  based on the start, stop times and time averaging interval. 
! 
!  The arguments are:
!  \begin{description}
!   \item[nsize] 
!     number of months computed by the routine
!   \item[ts] 
!     time step (in seconds) to be used for computing time span. 
!  \end{description}
!EOP    
    
    type(ESMF_Time)         :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: nsize
    integer                 :: status

 
    call ESMF_TimeIntervalSet(timeStep, s = 86400, &
         rc=status)
    call LVT_verify(status, &
         'ESMF_TimeIntervalSet failed in LVT_computeTimeSpanInYears')

    call ESMF_ClockGet(LVT_clock, startTime = startTime, &
         stopTime = stopTime, rc=status)
    call LVT_verify(status, 'ESMF_ClockGet failed in LVT_computeTimeSpanInYears')

    nsize = nint((stopTime - startTime)/timestep) + 1
    nyears = nint(real(nsize)/365.0) 

  end subroutine LVT_computeTimeSpanInYears


!BOP
! 
! !ROUTINE: LVT_getYearIndex
! \label{LVT_getYearIndex}
! 
! !INTERFACE: 
  subroutine LVT_getYearIndex(LVT_rc, yr_index, leap_year)
! !ARGUMENTS:     
    type(lvtrcdec)      :: LVT_rc
    integer             :: yr_index
    logical             :: leap_year

!
! !DESCRIPTION: 
!  This subroutine computes the time span of the LVT analysis
!  based on the start, stop times and time averaging interval. 
! 
!  The arguments are:
!  \begin{description}
!   \item[nsize] 
!     number of months computed by the routine
!   \item[ts] 
!     time step (in seconds) to be used for computing time span. 
!  \end{description}
!EOP    
    
    yr_index = LVT_rc%yr - LVT_rc%syr + 1
    
    if((mod(LVT_rc%yr,4) .eq. 0 .and. mod(LVT_rc%yr, 100).ne.0) &!leap year
         .or.(mod(LVT_rc%yr,400) .eq.0)) then 
       leap_year = .true. 
    else
       leap_year = .false.
    endif
       
  end subroutine LVT_getYearIndex


!BOP
! 
! !ROUTINE: LVT_getTimeSpanIndex
! \label{LVT_getTimeSpanIndex}
! 
! !INTERFACE: 
  subroutine LVT_getTimeSpanIndex(nsize, ts)
! !ARGUMENTS:     
    integer             :: nsize
    integer             :: ts
!
! !DESCRIPTION: 
!  This subroutine computes the time span of the LVT analysis
!  based on the start, stop times and time averaging interval. 
! 
!  The arguments are:
!  \begin{description}
!   \item[nsize] 
!     number of months computed by the routine
!   \item[ts] 
!     time step (in seconds) to be used for computing time span. 
!  \end{description}
!EOP    
    
    type(ESMF_Time)         :: startTime, currTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: status

 
    call ESMF_TimeIntervalSet(timeStep, s = ts, &
         rc=status)
    call LVT_verify(status, &
         'ESMF_TimeIntervalSet failed in LVT_getTimeSpanIndex')

    call ESMF_ClockGet(LVT_clock, startTime = startTime, &
         currTime = currTime, rc=status)
    call LVT_verify(status, 'ESMF_ClockGet failed in LVT_getTimeSpanIndex')

    nsize = nint((currTime - startTime)/timestep) + 1

  end subroutine LVT_getTimeSpanIndex

!BOP
! !ROUTINE: LVT_update_timestep
! \label{LVT_update_timestep}
! 
! !INTERFACE: 
  subroutine LVT_update_timestep(LVT_rc, ts)
    
    implicit none
! !ARGUMENTS: 
    type(lvtrcdec) :: LVT_rc
    integer        :: ts
!EOP

    if(ts.lt.LVT_rc%ts) then        
       LVT_rc%ts = ts

       if(LVT_rc%tavgInterval.lt.LVT_rc%ts) then 
          write(LVT_logunit,*) '[ERR] '
          write(LVT_logunit,*) '[ERR] Time averaging interval is less than the'
          write(LVT_logunit,*) '[ERR] data frequency of one of the datastreams'
          call LVT_endrun()
       endif
    endif

  end subroutine LVT_update_timestep

!BOP
! !ROUTINE: LVT_update_clock
! \label{LVT_update_clock}
! 
! !INTERFACE: 
  subroutine LVT_update_clock(ts)
! !USES: 
    use LVT_logMod,   only : LVT_verify
    
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
    call LVT_verify(rc,'ESMF_TimeIntervalSet failed in LVT_timeMgrMod')

    call ESMF_ClockSet(clock=LVT_clock,timestep=timeStep,rc=rc)
    call LVT_verify(rc,'ESMF_ClockSet failed in LVT_timeMgrMod')

  end subroutine LVT_update_clock

end module LVT_timeMgrMod
