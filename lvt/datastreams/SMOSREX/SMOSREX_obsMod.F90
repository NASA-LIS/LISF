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
! !MODULE: SMOSREX_obsMod
! \label(SMOSREX_obsMod)
!
! !INTERFACE:
module SMOSREX_obsMod
! 
! !USES: 
  use ESMF

  implicit none

  PRIVATE 
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
!  18 Apr 2009   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSREX_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSREXobs
!EOP
  type, public :: smosrexobsdec
     character*100        :: odir
     integer              :: nstns
     integer,     allocatable :: stnid(:)
     real,        allocatable :: stnlat(:)
     real,        allocatable :: stnlon(:)
     logical                 :: startFlag
     integer                 :: nts
     real                    :: udef
     type(ESMF_Time)         :: startTime, stopTime
     type(ESMF_TimeInterval) :: timeStep
     
     integer                 :: yr 

     real, allocatable           :: sm(:) 
     real, allocatable           :: rootsm(:,:) 
  end type smosrexobsdec

  type(smosrexobsdec), allocatable:: smosrexobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMOSREX_obsInit
! \label{SMOSREX_obsInit}
!
! !INTERFACE: 
  subroutine SMOSREX_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_timeMgrMod
    use LVT_histDataMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading SMOSREX data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: status, rc
    integer                 :: ftn, k
    integer                 :: ndata
    integer                 :: yr, mo, da, hr, mn, ss
    real*8                    :: jultime
    real                      :: smvalue
    integer                   :: offset
    type(ESMF_Time)           :: currtime


    if(.not.allocated(smosrexobs)) then 
       allocate(smosrexobs(LVT_rc%nDataStreams))
    endif

    smosrexobs(i)%startFlag = .true. 

    call ESMF_ConfigGetAttribute(LVT_Config, smosrexobs(i)%odir, &
         label='SMOSREX observation filename:', rc=status)
    call LVT_verify(status, 'SMOSREX observation filename: not defined')

!    call ESMF_TimeSet(smosrexobs(i)%starttime, yy=2002, mm=12,dd=19, h=0, m=0,s=0,&
!         calendar=LVT_calendar, rc=status)
!    call LVT_verify(status, 'smosrexobs(i)%starttime set failed in SMOSREX_obsMod')

    call ESMF_TimeSet(smosrexobs(i)%starttime, yy=2003, mm=1,dd=1, h=0, m=0,s=0,&
         calendar=LVT_calendar, rc=status)
    call LVT_verify(status, 'smosrexobs(i)%starttime set failed in SMOSREX_obsMod')

    call LVT_update_timestep(LVT_rc, 1800)

    call ESMF_TimeIntervalSet(smosrexobs(i)%timestep, s=1800,rc=status)
    call LVT_verify(status, 'timestep set failed in SMOSREX_obsMod')

    ndata = 87648 
    ftn=LVT_getNextUnitNumber()
    open(ftn,file=trim(smosrexobs(i)%odir),status='old')
    write(LVT_logunit,*) '[INFO] Reading SMOSREX obs file ',trim(smosrexobs(i)%odir)
    read(ftn,*) 

    allocate(smosrexobs(i)%sm(ndata))
    do k=1, ndata
       read(ftn,*) mo, da, yr, hr, mn, ss
       read(ftn,*) smvalue
       
!       print*, yr, mo, da, hr, mn, smvalue
!       write(LVT_logunit,*) jultime, smvalue
            
!       call convert_jdn2date(jultime, yr, mo, da, hr, mn, ss) 

       call ESMF_TimeSet(currtime, yy=yr, mm=mo,dd=da, h=hr, m=mn,s=ss,&
            calendar=LVT_calendar, rc=status)
       call LVT_verify(status, 'currtime set failed in SMOSREX_obsMod')

       offset = (currtime - smosrexobs(i)%starttime)/smosrexobs(i)%timestep + 1
!       print*, yr, mo, da, hr, mn, jultime

!       print*, offset, smvalue
       smosrexobs(i)%sm(offset) = smvalue
    enddo
    call LVT_releaseUnitNumber(ftn)
       
  end subroutine SMOSREX_obsinit


!BOP
! 
! !ROUTINE: convert_jdn2date
! \label(convert_jdn2date)
!
! !INTERFACE:
  subroutine convert_jdn2date(jultime, ryr,rmo, rda, rhr, rmn, rss)
! 
! !USES:   
    use LVT_timeMgrMod,   only : LVT_tick
    implicit none
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
    integer                 :: yr, mo, da, hr, mn, ss
    integer                 :: ryr, rmo, rda, rhr, rmn, rss, ts
    real*8                    :: jultime, jdn, jultime1
    real                      :: gmt
    integer                   :: doy
    real*8                    :: time
    real*8                    :: minval,diffval, prev_minval

    ryr = 2002
    rmo = 12
    rda = 19
    rhr = 0
    rmn = 0 
    rss = 0 
!assume 1 hr timestep 
    ts = 1800

    yr = ryr
    mo = rmo
    da = rda
    hr = rhr
    mn = rmn
    ss = rss
    minval = 0.0
    prev_minval = 1E8 

    jdn = 367*yr - (7*(yr+5001+(mo-9)/7))/4 + (275*mo)/9 + da + 1729777
    jultime1 = jdn + (hr-12)/24.0 + mn/1440.0 + ss/86400.0 

    diffval = abs(jultime - jultime1)

!    print*, 'diff ',abs(diffval-prev_minval)
!    do while(abs(diffval-prev_minval).gt.0.0) 
    do while(.true.)

       call LVT_tick(time, doy, gmt, yr, mo, da, hr, mn, ss, ts)

       !2000 jan 1 is JD = 2451545
       !wikipedia

       jdn = 367*yr - (7*(yr+5001+(mo-9)/7))/4 + (275*mo)/9 + da + 1729777
       jultime1 = jdn + (hr-12)/24.0 + mn/1440.0 + ss/86400.0 

       diffval = abs(jultime - jultime1)
!       print*, 'diff2 ',diffval, prev_minval
       if(diffval.lt.prev_minval) then 
          prev_minval = diffval
       endif
!       print*, 'diff3 ',diffval, prev_minval
       if(diffval.gt.prev_minval) then 
          ryr = yr
          rmo = mo
          rda = da
          rhr = hr
          rmn = mn
          rss = ss 
          
          call LVT_tick(time, doy, gmt, ryr, rmo, rda, rhr, rmn, rss, -ts)
!          print*, 'exiting ',diffval, prev_minval
          exit
       endif
    enddo

  end subroutine convert_jdn2date
end module SMOSREX_obsMod
