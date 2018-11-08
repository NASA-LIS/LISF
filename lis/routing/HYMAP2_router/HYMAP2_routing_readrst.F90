!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine HYMAP2_routing_readrst

  ! Augusto Getirana - 11/15/2011
  use ESMF
  use LIS_fileIOMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use HYMAP2_routingMod, only : HYMAP2_routing_struc

  implicit none

  integer       :: n 
  integer       :: ftn
  integer       :: i,j,k
  integer       :: ios,status
  character*100 :: filename
  logical       :: read_restart
  integer           :: yr,mo,da,hr,mn,ss,doy
  real*8            :: time
  real              :: gmt

  if(LIS_masterproc) then   
     do n=1, LIS_rc%nnest

        read_restart = .false. 

        if(HYMAP2_routing_struc(n)%startmode.eq."restart") then !cold start
           read_restart = .true. 
        endif

        if(LIS_rc%runmode.eq."ensemble smoother") then 
           if(LIS_rc%iterationId(n).gt.1) then 
              read_restart = .true. 
              
              if(HYMAP2_routing_struc(n)%rstInterval.eq.2592000) then 
                 !create the restart filename based on the timewindow start time
                 call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                      dd=da,calendar=LIS_calendar,rc=status)
                 hr = 0 
                 mn = 0 
                 ss = 0 
                 call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,&
                      (-1)*(HYMAP2_routing_struc(n)%dt))
              else
                 call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                      dd=da,calendar=LIS_calendar,rc=status)
                 hr = 0 
                 mn = 0 
                 ss = 0 
              endif
              
              call LIS_create_restart_filename(n,filename,&
                   'ROUTING','HYMAP2_router',&
                   yr, mo, da, hr, mn, ss, & 
                   wformat="binary")
              
              HYMAP2_routing_struc(n)%rstfile = filename
           endif
        endif

        if(read_restart) then 
           write(LIS_logunit,*) 'HYMAP2 restart file used: ', &
                                 trim(HYMAP2_routing_struc(n)%rstfile)
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=trim(HYMAP2_routing_struc(n)%rstfile),&
                status='old',form='unformatted',&
                iostat = ios)
           if(ios.ne.0) then 
              write(LIS_logunit,*) 'File '//trim(HYMAP2_routing_struc(n)%rstfile)//& 
                   'does not exist '
              call LIS_endrun()
           endif
           read(ftn) HYMAP2_routing_struc(n)%rivsto
           read(ftn) HYMAP2_routing_struc(n)%fldsto
           read(ftn) HYMAP2_routing_struc(n)%rnfsto
           read(ftn) HYMAP2_routing_struc(n)%bsfsto
           if(HYMAP2_routing_struc(n)%flowtype==3)then
             read(ftn) HYMAP2_routing_struc(n)%rivout_pre
             read(ftn) HYMAP2_routing_struc(n)%rivdph_pre
             read(ftn) HYMAP2_routing_struc(n)%fldout_pre
             read(ftn) HYMAP2_routing_struc(n)%flddph_pre
           endif

           call LIS_releaseUnitNumber(ftn)     
        endif
     enddo
  end if
end subroutine HYMAP2_routing_readrst
