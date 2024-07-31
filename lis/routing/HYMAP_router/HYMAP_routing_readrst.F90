!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine HYMAP_routing_readrst

  ! Augusto Getirana - 11/15/2011
  use ESMF
  use LIS_fileIOMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use HYMAP_routingMod, only : HYMAP_routing_struc

  implicit none

  integer       :: n 
  integer       :: ftn
  integer       :: i,j,k
  integer       :: ios,status
  character(len=LIS_CONST_PATH_LEN) :: filename
  logical       :: read_restart
  integer       :: yr,mo,da,hr,mn,ss,doy
  real*8        :: time
  real          :: gmt

  if(LIS_masterproc) then   
     do n=1, LIS_rc%nnest

        read_restart = .false. 

        if(HYMAP_routing_struc(n)%startmode.eq."restart") then !cold start
           read_restart = .true. 
        endif

        if(LIS_rc%runmode.eq."ensemble smoother") then 
           if(LIS_rc%iterationId(n).gt.1) then 
              read_restart = .true. 
              
              if(HYMAP_routing_struc(n)%rstInterval.eq.2592000) then 
                 !create the restart filename based on the timewindow start time
                 call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                      dd=da,calendar=LIS_calendar,rc=status)
                 hr = 0 
                 mn = 0 
                 ss = 0 
                 call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,&
                      (-1)*(HYMAP_routing_struc(n)%dt))
              else
                 call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                      dd=da,calendar=LIS_calendar,rc=status)
                 hr = 0 
                 mn = 0 
                 ss = 0 
              endif
              
              call LIS_create_restart_filename(n,filename,&
                   'ROUTING','HYMAP_router',&
                   yr, mo, da, hr, mn, ss, & 
                   wformat="binary")
              
              HYMAP_routing_struc(n)%rstfile = filename
           endif
        endif

        if(read_restart) then 
           write(LIS_logunit,*) '[INFO] HYMAP restart file used: ', &
                                 trim(HYMAP_routing_struc(n)%rstfile)
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=trim(HYMAP_routing_struc(n)%rstfile),&
                status='old',form='unformatted',&
                iostat = ios)
           if(ios.ne.0) then 
              write(LIS_logunit,*) '[ERR] File '//trim(HYMAP_routing_struc(n)%rstfile)//& 
                   'does not exist '
              call LIS_endrun()
           endif
           read(ftn) HYMAP_routing_struc(n)%rivsto
           read(ftn) HYMAP_routing_struc(n)%fldsto
           read(ftn) HYMAP_routing_struc(n)%rnfsto
           read(ftn) HYMAP_routing_struc(n)%bsfsto

           call LIS_releaseUnitNumber(ftn)     
        endif
     enddo
  end if
end subroutine HYMAP_routing_readrst
