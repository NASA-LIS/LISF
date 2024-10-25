!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine NLDAS_routing_readrst

  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_logMod,  only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_endrun
  use LIS_timeMgrMod 
  use LIS_fileIOMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use NLDAS_routingMod, only : NLDAS_routing_struc

  implicit none

  integer       :: n 
  integer       :: ftn
  integer       :: i,j,k
  integer       :: ios,status
  character(len=LIS_CONST_PATH_LEN) :: filen
  integer           :: yr,mo,da,hr,mn,ss,doy
  real*8            :: time
  real              :: gmt
  logical       :: read_restart

  if(LIS_masterproc) then   
     do n=1, LIS_rc%nnest
        read_restart = .false. 
        if(trim(NLDAS_routing_struc(n)%startmode).eq."restart") then 
           read_restart = .true. 
        endif

        if(LIS_rc%runmode.eq."ensemble smoother") then 
           if(LIS_rc%iterationId(n).gt.1) then 
              read_restart = .true. 
              
              if(NLDAS_routing_struc(n)%rstInterval.ne.2592000) then 
                 write(LIS_logunit,*) 'restart interval must be set to a month' 
                 write(LIS_logunit,*) 'when running in the ensemble smoother mode'
                 call LIS_endrun()
              endif
              
              !create the restart filename based on the timewindow start time
              call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                   dd=da,calendar=LIS_calendar,rc=status)
              hr = 0 
              mn = 0 
              ss = 0 
              !call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,-86400.0)
              call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,-1*NLDAS_routing_struc(n)%dt)
              
              call LIS_create_restart_filename(n,filen,&
                   'ROUTING','NLDAS_router',yr,mo,da,hr,mn,ss,&
                   wformat="binary")
              
              NLDAS_routing_struc(n)%rstfile = filen
           endif
        endif

        if(.not. read_restart) then 
           write(LIS_logunit,*) 'Reading ', &
                trim(NLDAS_routing_struc(n)%initial_1)
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=trim(NLDAS_routing_struc(n)%initial_1),&
                status='old',form='unformatted',&
                iostat = ios)
           if(ios.ne.0) then 
              write(LIS_logunit,*) &
                   'File '//trim(NLDAS_routing_struc(n)%initial_1)//& 
                   'does not exist '
              call LIS_endrun()
           endif
           read(ftn) (((NLDAS_routing_struc(n)%runoff_intern(i,j,k), &
                i=1,NLDAS_routing_struc(n)%luh), j=1,LIS_rc%gnc(n)),&
                k=1,LIS_rc%gnr(n))
           
           call LIS_releaseUnitNumber(ftn)
           
           write(LIS_logunit,*) 'Reading ', &
                trim(NLDAS_routing_struc(n)%initial_2)
           
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=trim(NLDAS_routing_struc(n)%initial_2),&
                status='old',form='unformatted',&
                iostat = ios)
           if(ios.ne.0) then 
              write(LIS_logunit,*) 'File '//trim(NLDAS_routing_struc(n)%initial_2)//& 
                   'does not exist '
              call LIS_endrun()
           endif
           read(ftn) (((NLDAS_routing_struc(n)%runoff_trans(i,j,k), &
                i=1,NLDAS_routing_struc(n)%ltr), j=1,LIS_rc%gnc(n)),&
                k=1,LIS_rc%gnr(n))
           
           call LIS_releaseUnitNumber(ftn)        
        else
           
           write(LIS_logunit,*) 'Reading ', trim(NLDAS_routing_struc(n)%rstfile)
           ftn = LIS_getNextUnitNumber()
           open(ftn,file=trim(NLDAS_routing_struc(n)%rstfile),&
                status='old',form='unformatted',&
                iostat = ios)
           if(ios.ne.0) then 
              write(LIS_logunit,*) 'File '//trim(NLDAS_routing_struc(n)%rstfile)//& 
                   'does not exist '
              call LIS_endrun()
           endif
           read(ftn) NLDAS_routing_struc(n)%runoff_intern
           read(ftn) NLDAS_routing_struc(n)%runoff_trans
           
           call LIS_releaseUnitNumber(ftn)     
        endif
     enddo
  end if

end subroutine NLDAS_routing_readrst
