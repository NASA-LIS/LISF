!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module ESPboot_Mod
!BOP
!
! !MODULE: ESPboot_Mod
!
! !DESCRIPTION:
!    This module provides the routines to control the execution of 
!    the tau-omega model 
!
! !REVISION HISTORY:
! 28 Aug 2012: Sujay Kumar, initial specification based on 
!              the code from Wade Crow
!
! !USES:        

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_forecastMod
  use LIS_numerRecipesMod
  use LIS_ran2_gasdev

  implicit none
 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ESPboot_initialize
  public :: ESPboot_sampledate
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  type, public :: espboot_type_dec
     type(ESMF_TimeInterval), allocatable :: twInterval(:)
     integer                              :: numIters
     integer                              :: iterId
     integer, allocatable                 :: stime(:,:)
     integer, allocatable                 :: etime(:,:)
     integer, allocatable                 :: sampleCounter(:)
     integer, allocatable                 :: ref_da(:)
     type(ESMF_TimeInterval), allocatable :: offset(:)
     logical                              :: sampleFlag
  end type espboot_type_dec
  
  type(espboot_type_dec), allocatable :: Espboot_struc(:)

!EOP
  SAVE

contains
!BOP
! 
! !ROUTINE: ESPboot_initialize
! \label{ESPboot_initialize}
! 
! !INTERFACE:
  subroutine ESPboot_initialize()
! !USES:

! !DESCRIPTION:        
!
!EOP
    implicit none
    
    integer :: n,i,k,t,rc
    character*10             :: time_val
    type(ESMF_Time)          :: startTime, endTime
    type(ESMF_TimeInterval)  :: deltaT
    real                     :: twInterval
    integer, allocatable     :: stime(:,:),etime(:,:)
    integer, allocatable     :: niterations(:)
! _____________________________________________________

    allocate(ESPboot_struc(LIS_rc%nnest))

    allocate(LIS_forecast_struc(1)%seed(LIS_rc%nmetforc, NRANDSEED))
    allocate(ESPboot_struc(1)%twInterval(LIS_rc%nmetforc))
    allocate(ESPboot_struc(1)%sampleCounter(LIS_rc%nmetforc))
    allocate(ESPboot_struc(1)%stime(LIS_rc%nmetforc,3))
    allocate(ESPboot_struc(1)%etime(LIS_rc%nmetforc,3))
    allocate(ESPboot_struc(1)%offset(LIS_rc%nmetforc))
    allocate(ESPboot_struc(1)%ref_da(LIS_rc%nmetforc))
    allocate(stime(LIS_rc%nmetforc, 3))
    allocate(etime(LIS_rc%nmetforc, 3))
    allocate(niterations(LIS_rc%nmetforc))

    if(LIS_forecast_struc(1)%startMode.eq."coldstart") then 
       LIS_forecast_struc(1)%seed = -100
       call init_randseed( LIS_forecast_struc(1)%seed )
    else
       call LIS_forecast_readrestart()
    endif

! New code: KRA
    call ESMF_ConfigFindLabel(LIS_config,"ESP boot number of iterations:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,ESPboot_struc(n)%numIters,rc=rc)
       call LIS_verify(rc,"ESP boot number of iterations: not defined")
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,"ESP boot sampling time window interval:",rc=rc)
    do n=1,LIS_rc%nnest
       do k=1,LIS_rc%nmetforc
          call ESMF_ConfigGetAttribute(LIS_config,time_val,rc=rc)
          call LIS_verify(rc,"ESP boot sampling time window interval: not defined")
          
          call LIS_parseTimeString(time_val,twInterval)
          call ESMF_TimeIntervalSet(ESPboot_struc(n)%twInterval(k),&
               s=nint(twInterval),rc=rc)
       enddo
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,"ESP boot start time of the forcing archive:",rc=rc)
    do n=1,LIS_rc%nnest
       do k=1,LIS_rc%nmetforc
          do i=1,3
             call ESMF_ConfigGetAttribute(LIS_config,stime(k,i),rc=rc)
             call LIS_verify(rc,"ESP boot start time of the forcing archive: not defined")
          enddo
       enddo
    enddo
    call ESMF_ConfigFindLabel(LIS_config,"ESP boot end time of the forcing archive:",rc=rc)
    do n=1,LIS_rc%nnest
       do k=1,LIS_rc%nmetforc
          do i=1,3
             call ESMF_ConfigGetAttribute(LIS_config,etime(k,i),rc=rc)
             call LIS_verify(rc,"ESP boot end time of the forcing archive: not defined")
          enddo
       enddo
    enddo
    
    do n=1,LIS_rc%nnest
       do k=1,LIS_rc%nmetforc
          call ESMF_TimeSet(startTime, yy = stime(k,1), &
               mm = stime(k,2), &
               dd = stime(k,3), &
               h  = 0, &
               m  = 0, & 
               s  = 0, &
               calendar = LIS_calendar, & 
               rc = rc)
          call LIS_verify(rc,'error in ESMF_TimeSet:startTime in ESPboot_initialize')
          
          call ESMF_TimeSet(endTime, yy = etime(k,1), &
               mm = etime(k,2), &
               dd = etime(k,3), &
               h  = 0, &
               m  = 0, &
               s  = 0, &
               calendar = LIS_calendar, & 
               rc = rc)
          call LIS_verify(rc,'error in ESMF_TimeSet:endTime in ESPboot_initialize')
          
          call ESMF_TimeIntervalSet(deltaT,s=86400*365,rc=rc)
          call LIS_verify(rc,'error in ESMF_TimeIntervalSet in ESPboot_initialize')
          
!          nIterations(k) = 100   ! Former code
          ESPboot_struc(n)%stime(k,:) = stime(k,:)
          ESPboot_struc(n)%etime(k,:) = etime(k,:)
          ESPboot_struc(n)%sampleCounter(k) = 0          
       enddo
    enddo
    
    do n=1,LIS_rc%nnest
       LIS_forecast_struc(n)%nIterations = 0 
       ! New code (KRA):
       if( ESPboot_struc(n)%numIters .gt. LIS_forecast_struc(n)%nIterations) then
          LIS_forecast_struc(n)%nIterations = ESPboot_struc(n)%numIters 
       endif
! Former code
!       do k=1,LIS_rc%nmetforc
!          if(niterations(k).gt.LIS_forecast_struc(n)%nIterations) then 
!             LIS_forecast_struc(n)%nIterations = niterations(k)
!          endif
!       enddo
    enddo

  end subroutine ESPboot_initialize


  subroutine ESPboot_sampledate(n, k, yr,mo,da)
    
    integer                 :: n
    integer                 :: k
    integer                 :: yr
    integer                 :: mo
    integer                 :: da

    integer                 :: yr_in
    integer                 :: da1
    real                    :: rand
    integer                 :: twInterval
    type(ESMF_Time)         :: stTime, currTime, currTime1
    type(ESMF_TimeInterval) :: deltaT
    integer                 :: rc


    if(ESPboot_struc(n)%sampleCounter(k).eq.0) then 
       
       yr_in = yr
       
       do while(yr_in.eq.yr) 

          call nr_ran2(LIS_forecast_struc(n)%seed(k,:), rand)
          
          call ESMF_TimeSet(currTime1, yy = yr, &
               mm = mo, dd=da,h=0,m=0,s=0,&
               calendar = LIS_calendar, & 
               rc = rc)
          call LIS_verify(rc,'error in ESMF_TimeSet:currTime1 in ESPboot_sampledate')
       
          yr = ESPboot_struc(n)%stime(k,1) + & 
                nint(rand*(ESPboot_struc(n)%etime(k,1) - &
                ESPboot_struc(n)%stime(k,1) ))

          write(LIS_logunit,*) "[INFO] ESPboot member index and year:",&
                                  ESPboot_struc(n)%sampleCounter(k), yr
       enddo

       da1= da 
       ! If current date is 29 and sampled year is not a leap year       
       if(da.eq.29.and.&
            .not.((mod(yr,4) .eq. 0 .and. mod(yr, 100).ne.0) &
            .or.(mod(yr,400) .eq.0))) then 
          da1 = 28
       endif
       
       call ESMF_TimeSet(currTime, yy = yr, &
            mm = mo, dd=da1,h=0,m=0,s=0,&
            calendar = LIS_calendar, & 
            rc = rc)
       call LIS_verify(rc,'error in ESMF_TimeSet:currTime in ESPboot_sampledate')
       
       stTime = currTime - ESPboot_struc(n)%twInterval(k)
       
       !       call LIS_rand_func(ESPboot_struc(n)%seed(k),rand)  
       call nr_ran2(LIS_forecast_struc(n)%seed(k,:), rand)
       
       call ESMF_TimeIntervalGet(ESPboot_struc(n)%twInterval(k),&
            d=twInterval,rc=rc)
       twInterval = nint(rand*twInterval) 
       
       call ESMF_TimeIntervalSet(deltaT,d=twInterval,rc=rc)
       currTime = stTime + deltaT
       
       call ESMF_TimeGet(currTime, yy = yr, &
            mm = mo, dd=da,&
            calendar = LIS_calendar, & 
            rc = rc)
       call LIS_verify(rc,'error in ESMF_TimeGet:startTime in ESPboot_sampledate')
       
       ESPboot_struc(n)%offSet(k) = currTime1 - currTime
       
       ESPboot_struc(n)%ref_da(k) = da
       ESPboot_struc(n)%sampleCounter(k) = &
            ESPboot_struc(n)%sampleCounter(k) + 1
       
    else
       
       call ESMF_TimeSet(currTime, yy = yr, &
            mm = mo, dd=da,&
            calendar = LIS_calendar, & 
            rc = rc)
       call LIS_verify(rc,'error in ESMF_TimeSet:currTime in ESPboot_sampledate')
       
       currTime = currTime - ESPboot_struc(n)%offSet(k)
       
       call ESMF_TimeGet(currTime, yy = yr, &
            mm = mo, dd=da,&
            calendar = LIS_calendar, & 
            rc = rc)
       call LIS_verify(rc,'error in ESMF_TimeGet:currTime in ESPboot_sampledate')
       
       if(da.ne.ESPboot_struc(n)%ref_da(k)) then 
          ESPboot_struc(n)%sampleCounter(k) = ESPboot_struc(n)%sampleCounter(k) + 1
          ESPboot_struc(n)%ref_da(k) = da
       endif
       ! reset after 10 days. 
       if(ESPboot_struc(n)%sampleCounter(k).gt.10) then 
          ESPboot_struc(n)%sampleCounter(k) = 0 
       endif
    endif

  end subroutine ESPboot_sampledate
end module ESPboot_Mod

