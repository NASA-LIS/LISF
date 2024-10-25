!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module GLSMod
!BOP
!
! !MODULE:
! 
! !DESCRIPTION:
! 
! !REVISION HISTORY: 
! 
  use ESMF

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GLS_init
  PUBLIC :: GLS_run
  PUBLIC :: GLS_output
  PUBLIC :: GLS_final
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------  
  public :: GLS_data
  public :: GLS_ctl

  type, public :: gls_varctl
     real, allocatable :: suscep_index(:)

     integer       :: out_ftn
     integer          :: outInterval
     type(ESMF_Alarm) :: outAlarm
     real             :: coeff1
     real             :: coeff2
     real, allocatable    :: ls_status(:)
     real, allocatable    :: ls_accum_status(:)   

     integer          :: rc !ring count
     integer       :: dt
     type(ESMF_Alarm) :: dt_alarm

     integer :: wform
     integer :: acc_interval1
     integer :: acc_interval2
     integer :: acc_interval3


     type(ESMF_Alarm) :: acc1_alarm
     type(ESMF_Alarm) :: acc2_alarm
     type(ESMF_Alarm) :: acc3_alarm

     logical         :: acc1_alarm_status
     logical         :: acc2_alarm_status
     logical         :: acc3_alarm_status

     real, allocatable :: rf_acc1(:)
     real, allocatable :: rf_acc2(:)
     real, allocatable :: rf_acc3(:)

     real, allocatable :: rf_acc1_sum(:,:)
     real, allocatable :: rf_acc2_sum(:,:)
     real, allocatable :: rf_acc3_sum(:,:)


  end type gls_varctl

  type, public :: gls_dec

     real :: rf_threshold1
     real :: rf_threshold2
     real :: rf_threshold3
     real :: suscep_index_threshold

  end type gls_dec
  
  type(gls_dec), allocatable     :: GLS_data(:)
  type(gls_varctl), save     :: GLS_ctl
!EOP

contains

!BOP
! 
! !ROUTINE: GLS_init
! \label{GLS_init}
! 
! !INTERFACE:
  subroutine GLS_init()
! !USES:   
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_domain
    use LIS_timeMgrMod, only : LIS_clock, LIS_calendar, LIS_seconds2time
    use LIS_fileIOMod, only : LIS_readDomainConfigSpecs, LIS_readData
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
!
! !DESCRIPTION:
!EOP
    character*100 :: suscep_map
    real          :: gls_gridDesc(1,6)
    integer       :: n,c,r,t
    integer       :: rc, status
    logical       :: file_exists
    integer       :: ftn
    integer       :: acc1, acc2, acc3
    real          :: threshold
    integer                 :: yr, mo, da, hr, mn, ss
    integer                 :: interval
    integer                 :: m 
    integer                 :: dt 
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: time2, alarmTime, currTime
    type(ESMF_TimeInterval) :: deltat
    real,     allocatable       :: suscep_index(:,:)
    integer   :: ts_frac_s, ts_frac_ms
    
    n = 1
    allocate(suscep_index(LIS_rc%lnc(n),LIS_rc%lnr(n)))
    allocate(GLS_ctl%suscep_index(LIS_rc%ntiles(n)))

    allocate(GLS_data(LIS_rc%nensem(n)))

    GLS_ctl%rc = 0
    
    GLS_ctl%acc1_alarm_status = .false. 
    GLS_ctl%acc2_alarm_status = .false. 
    GLS_ctl%acc3_alarm_status = .false. 

    call ESMF_ConfigGetAttribute(LIS_config,suscep_map,&
         label="GLS susceptibility index map:",rc=rc)
    call LIS_verify(rc,'GLS susceptibility index map: option not specified in the config file')
    
    call ESMF_ConfigGetAttribute(LIS_config,threshold,&
         label="GLS susceptibility index threshold:",rc=rc)
    call LIS_verify(rc,'GLS susceptibility index threshold: not specified')

    do m=1,LIS_rc%nensem(n)
       GLS_data(m)%suscep_index_threshold = threshold
    enddo

    call ESMF_ConfigGetAttribute(LIS_config,GLS_ctl%coeff1,&
         label="GLS slope value:",rc=rc)
    call LIS_verify(rc,'GLS slope value: not specified')

    call ESMF_ConfigGetAttribute(LIS_config,GLS_ctl%coeff2,&
         label="GLS y-intercept value:",rc=rc)
    call LIS_verify(rc,'GLS y-intercept value: not specified')


    call ESMF_ConfigGetAttribute(LIS_config,interval,&
         label="GLS rainfall accumulation interval1:",rc=rc)
    call LIS_verify(rc,'GLS rainfall accumulation interval1: not specified')

    GLS_ctl%acc_interval1 = interval
       
    call ESMF_ConfigGetAttribute(LIS_config,interval,&
         label="GLS rainfall accumulation interval2:",rc=rc)
    call LIS_verify(rc,'GLS rainfall accumulation interval2: not specified')

    GLS_ctl%acc_interval2 = interval


    call ESMF_ConfigGetAttribute(LIS_config,interval,&
         label="GLS rainfall accumulation interval3:",rc=rc)
    call LIS_verify(rc,'GLS rainfall accumulation interval3: not specified')

    GLS_ctl%acc_interval3 = interval

    call ESMF_ConfigGetAttribute(LIS_config,threshold,&
         label="GLS rainfall threshold for interval1:",rc=rc)
    call LIS_verify(rc,'GLS rainfall threshold for interval1: not specified')

    do m=1,LIS_rc%nensem(n)
       GLS_data(m)%rf_threshold1 = threshold
    enddo

    call ESMF_ConfigGetAttribute(LIS_config,threshold,&
         label="GLS rainfall threshold for interval2:",rc=rc)
    call LIS_verify(rc,'GLS rainfall threshold for interval2: not specified')

    do m=1,LIS_rc%nensem(n)
       GLS_data(m)%rf_threshold2 = threshold
    enddo

    call ESMF_ConfigGetAttribute(LIS_config,threshold,&
         label="GLS rainfall threshold for interval3:",rc=rc)
    call LIS_verify(rc,'GLS rainfall threshold for interval3: not specified')

    do m=1,LIS_rc%nensem(n)
       GLS_data(m)%rf_threshold3 = threshold
    enddo

    call ESMF_ConfigGetAttribute(LIS_config,dt,&
         label="GLS model timestep:",rc=rc)
    call LIS_verify(rc,'GLS model timestep: not specified')

    GLS_ctl%dt = dt

    call ESMF_ConfigGetAttribute(LIS_config, GLS_ctl%wform, &
         label="GLS output format:",rc=rc)
    call LIS_verify(rc,"GLS output format: not defined")
    
    call ESMF_ConfigGetAttribute(LIS_config, GLS_ctl%outInterval, &
         label="GLS output interval:",rc=rc)
    call LIS_verify(rc,"GLS output interval: not defined")

    acc1 = max(1,nint(GLS_ctl%acc_interval1/24.0))
    acc2 = max(1,nint(GLS_ctl%acc_interval2/24.0))
    acc3 = max(1,nint(GLS_ctl%acc_interval3/24.0))

    allocate(GLS_ctl%rf_acc1(LIS_rc%ntiles(n)))
    allocate(GLS_ctl%rf_acc2(LIS_rc%ntiles(n)))
    allocate(GLS_ctl%rf_acc3(LIS_rc%ntiles(n)))
       
    allocate(GLS_ctl%rf_acc1_sum(LIS_rc%ntiles(n),acc1))
    allocate(GLS_ctl%rf_acc2_sum(LIS_rc%ntiles(n),acc2))
    allocate(GLS_ctl%rf_acc3_sum(LIS_rc%ntiles(n),acc3))

    GLS_ctl%rf_acc1 = 0.0
    GLS_ctl%rf_acc2 = 0.0
    GLS_ctl%rf_acc3 = 0.0
       
    GLS_ctl%rf_acc1_sum = 0.0
    GLS_ctl%rf_acc2_sum = 0.0
    GLS_ctl%rf_acc3_sum = 0.0
    
    call LIS_readDomainConfigSpecs("GLS susceptibility map",gls_gridDesc)
    
    inquire(file=trim(suscep_map), exist=file_exists)
    if(.not.file_exists) then 
       write(LIS_logunit,*) 'susceptibility map ',trim(suscep_map),' not found'
       write(LIS_logunit,*) 'Program stopping ...'
       call LIS_endrun
    endif
    write(LIS_logunit,*) 'Reading susceptibilty map ',trim(suscep_map)

    ftn = LIS_getNextUnitNumber()

    open(ftn,file=trim(suscep_map),form='unformatted', &
         access='direct',recl=4,status='old')
    
    call LIS_readData(n,ftn,gls_gridDesc(n,:),suscep_index)
    
    call LIS_releaseUnitNumber(ftn)
    write(LIS_logunit,*) 'Done reading susceptibilty map ',trim(suscep_map)

    do t=1,LIS_rc%ntiles(n)
       c = LIS_domain(n)%tile(t)%col
       r = LIS_domain(n)%tile(t)%row
       GLS_ctl%suscep_index(t) = suscep_index(c,r)
    enddo
!    open(100,file='suscep.bin',form='unformatted')
!    write(100) GLS_data%suscep_index
!    close(100)
!    stop
! setting alarms
! model alarm - set it to the nearest 0z time
    call ESMF_TimeIntervalSet(alarmInterval, h=GLS_ctl%dt, rc=rc)
    call LIS_verify(rc, 'Time interval set failed in GLS_init')
    
    call ESMF_TimeSet(time2, yy = LIS_rc%yr, &
         mm = LIS_rc%mo, &
         dd = LIS_rc%da, &
         h  = 0, &
         m  = 0, & 
         s  = 0, &
         calendar = LIS_calendar,&
         rc = rc)
    call LIS_verify(rc, 'Error in time2 set in GLS_init')
    
    alarmTime = time2 + alarmInterval
    
    GLS_ctl%dt_alarm = ESMF_AlarmCreate(clock=LIS_clock, &
         ringTime = alarmTime, &
         ringInterval = alarmInterval, rc=rc)
    call LIS_verify(rc, 'Alarm Create failed in GLS_init')

!accumulation alarm 1
    call ESMF_TimeIntervalSet(alarmInterval, &
         h=GLS_ctl%acc_interval1, rc=rc)
    call LIS_verify(rc, 'Time interval set failed in GLS_init')

    call ESMF_TimeSet(time2, yy = LIS_rc%yr, &
         mm = LIS_rc%mo, &
         dd = LIS_rc%da, &
         h  = 0, &
         m  = 0, & 
         s  = 0, &
         calendar = LIS_calendar,&
         rc = rc)
    call LIS_verify(rc, 'Error in time2 set in GLS_init')

    alarmTime = time2 + alarmInterval
       
    GLS_ctl%acc1_alarm = ESMF_AlarmCreate(clock=LIS_clock, &
         ringTime = alarmTime, &
         ringInterval = alarmInterval, rc=rc)
    call LIS_verify(rc, 'Alarm Create failed in GLS_init')

!accumulation alarm 2
    call ESMF_TimeIntervalSet(alarmInterval, &
         h=GLS_ctl%acc_interval2, rc=rc)
    call LIS_verify(rc, 'Time interval set failed in GLS_init')
    
    call ESMF_TimeSet(time2, yy = LIS_rc%yr, &
         mm = LIS_rc%mo, &
         dd = LIS_rc%da, &
         h  = 0, &
         m  = 0, & 
         s  = 0, &
         calendar = LIS_calendar,&
         rc = rc)
    call LIS_verify(rc, 'Error in time2 set in GLS_init')

    alarmTime = time2 + alarmInterval
    
    call ESMF_TimeIntervalSet(alarmInterval, &
         h=GLS_ctl%acc_interval1, rc=rc)
    call LIS_verify(rc, 'Time interval set failed in GLS_init')
    
    GLS_ctl%acc2_alarm = ESMF_AlarmCreate(clock=LIS_clock, &
         ringTime = alarmTime, &
         ringInterval = alarmInterval, rc=rc)
    call LIS_verify(rc, 'Alarm Create failed in GLS_init')
    
    !accumulation alarm 3
    call ESMF_TimeIntervalSet(alarmInterval, &
         h=GLS_ctl%acc_interval3, rc=rc)
    call LIS_verify(rc, 'Time interval set failed in GLS_init')
    
    call ESMF_TimeSet(time2, yy = LIS_rc%yr, &
         mm = LIS_rc%mo, &
         dd = LIS_rc%da, &
         h  = 0, &
         m  = 0, & 
         s  = 0, &
         calendar = LIS_calendar,&
         rc = rc)
    call LIS_verify(rc, 'Error in time2 set in GLS_init')
    
    alarmTime = time2 + alarmInterval

    call ESMF_TimeIntervalSet(alarmInterval, &
         h=GLS_ctl%acc_interval1, rc=rc)
    call LIS_verify(rc, 'Time interval set failed in GLS_init')


    GLS_ctl%acc3_alarm = ESMF_AlarmCreate(clock=LIS_clock, &
         ringTime = alarmTime, &
         ringInterval = alarmInterval, rc=rc)
    call LIS_verify(rc, 'Alarm Create failed in GLS_init')
!output file 
    GLS_ctl%out_ftn = LIS_getNextUnitNumber()
    
    open(GLS_ctl%out_ftn, file='GLS_output.dat',form='formatted')

    allocate(GLS_ctl%ls_status(LIS_rc%ntiles(n)))
    allocate(GLS_ctl%ls_accum_status(LIS_rc%ntiles(n)))
    
    GLS_ctl%ls_status = 0 
    GLS_ctl%ls_accum_status = 0 

    call LIS_seconds2time(GLS_ctl%outInterval, da,hr,mn,ss)
       
    call ESMF_TimeIntervalSet(alarmInterval, &
         d = da, h = hr, m = mn, s = ss, &
         calendar = LIS_calendar, rc=status)
    
    call ESMF_ClockGet(LIS_clock,currTime=currTime,rc=status)       
    call ESMF_TimeGet(currTime,yy=yr,mm=mo,dd=da,h=hr,m=mn,s=ss,rc=status)
    
    if(GLS_ctl%outInterval.ge.3600) then 
       do while(mod(real(hr)*3600+60*real(mn)+float(ss),&
            real(GLS_ctl%outInterval)).ne.0)
          ts_frac_s = nint(LIS_rc%ts)
          ts_frac_ms = (LIS_rc%ts - nint(LIS_rc%ts))*1000
          
          call ESMF_TimeIntervalSet(deltaT,s=ts_frac_s, ms = ts_frac_ms,rc=status)
          currTime = currTime + deltaT
          call ESMF_TimeGet(currTime,yy=yr,mm=mo,dd=da,h=hr,m=mn,s=ss,rc=status)
       enddo
       call ESMF_TimeSet(alarmTime, yy = yr, &
            mm = mo, &
            dd = da, &
            h  = hr, &
            m  = mn, & 
            s  = ss, &
            calendar = LIS_calendar,&
            rc = status)
       call LIS_verify(status)
    else
       call ESMF_TimeSet(alarmTime, yy = LIS_rc%yr, &
            mm = LIS_rc%mo, &
            dd = LIS_rc%da, &
            h  = hr, &
            m  = mn, & 
            s  = ss, &
            calendar = LIS_calendar,&
            rc = status)
       call LIS_verify(status)
       
       alarmTime = alarmTime + alarmInterval
    endif

    GLS_ctl%outAlarm = ESMF_AlarmCreate(clock=LIS_Clock,&
         ringTime = alarmTime, &
         ringInterval = alarmInterval,&
         rc=status)
    call LIS_verify(status)

    deallocate(suscep_index)

  end subroutine GLS_init

!BOP
! 
! !ROUTINE: GLS_run
! \label{GLS_run}
! 
! !INTERFACE:
  subroutine GLS_run(n)
! !USES:
    use LIS_coreMod,  only : LIS_rc, LIS_domain
    use LIS_logMod,         only : LIS_verify, LIS_logunit
!
! !DESCRIPTION:
!   This routines computes the rainfall running sums, and 
!   writes out landslide forecasts based on susceptibility 
!   tests. 
!EOP
    implicit none

    integer, intent(IN) :: n 
    
    integer             :: c,r,t,g,m
    integer             :: gindex
    real                :: lat, lon
    integer             :: tindex
 !   real                :: p15
    logical             :: criteria1, criteria2, criteria3
    character*3         :: data_status
    integer             :: rc, status   

    if(ESMF_AlarmIsRinging(GLS_ctl%dt_alarm, rc=rc)) then 
       GLS_ctl%rc = GLS_ctl%rc + 1
    endif
       
! Do  accumulations
! Accumulation 1
    if(ESMF_AlarmIsRinging(GLS_ctl%acc1_alarm, rc=rc)) then 
       call ESMF_AlarmRingerOff(GLS_ctl%acc1_alarm, rc=rc)
 !resets the data for 1st accumulation interval.
       if(GLS_ctl%acc_interval1.le.24) then 
          tindex = 1
       else
          tindex = mod(GLS_ctl%rc,GLS_ctl%acc_interval1/24)+1
       endif

       GLS_ctl%rf_acc1 = GLS_ctl%rf_acc1_sum(:,tindex)
       GLS_ctl%rf_acc1_sum(:,tindex) = 0.0
       GLS_ctl%acc1_alarm_status = .true. 
    else
       call compute_rainf_accumulations (n, GLS_ctl%acc_interval1, &
            GLS_ctl%rc, &
            GLS_ctl%rf_acc1_sum) 
       
       GLS_ctl%rf_acc1 = LIS_rc%udef
    endif
    
! Accumulation 2
    if(ESMF_AlarmIsRinging(GLS_ctl%acc2_alarm, rc=rc)) then 
       call ESMF_AlarmRingerOff(GLS_ctl%acc2_alarm, rc=rc)       
       !resets the data for 2nd accumulation interval.       
       if(GLS_ctl%acc_interval2.le.24) then 
          tindex = 1
       else
          tindex = mod(GLS_ctl%rc,GLS_ctl%acc_interval2/24)+1
       endif
       
       GLS_ctl%rf_acc2 = GLS_ctl%rf_acc2_sum(:,tindex)
       GLS_ctl%rf_acc2_sum(:,tindex) = 0.0
       GLS_ctl%acc2_alarm_status = .true. 
    else
       call compute_rainf_accumulations (n, GLS_ctl%acc_interval2, &
            GLS_ctl%rc, &
            GLS_ctl%rf_acc2_sum) 
       
       GLS_ctl%rf_acc2 = LIS_rc%udef
    endif
    
! Accumulation 3
    if(ESMF_AlarmIsRinging(GLS_ctl%acc3_alarm, rc=rc)) then 
       call ESMF_AlarmRingerOff(GLS_ctl%acc3_alarm, rc=rc)       
!resets the data for 3rd accumulation interval.    
       if(GLS_ctl%acc_interval3.le.24) then 
          tindex = 1
       else   
          tindex = mod(GLS_ctl%rc,GLS_ctl%acc_interval3/24)+1
       endif
       
       GLS_ctl%rf_acc3 = GLS_ctl%rf_acc3_sum(:,tindex)
       GLS_ctl%rf_acc3_sum(:,tindex) = 0.0
       GLS_ctl%acc3_alarm_status = .true. 
    else
       call compute_rainf_accumulations (n, GLS_ctl%acc_interval3, &
            GLS_ctl%rc, &
            GLS_ctl%rf_acc3_sum) 
       
       GLS_ctl%rf_acc3 = LIS_rc%udef
    endif


! run the algorithm 
    if(ESMF_AlarmIsRinging(GLS_ctl%dt_alarm, rc=rc)) then 
          
       GLS_ctl%ls_status = 0 
       call ESMF_AlarmRingerOff(GLS_ctl%dt_alarm, rc=rc)
       do g=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
       
          do m=1,LIS_rc%nensem(n)      
             criteria1 = .false. 
             criteria2 = .false. 
             criteria3 = .false. 
                
             t = (g-1)*LIS_rc%nensem(n)+m
             
             c = LIS_domain(n)%tile(t)%col
             r = LIS_domain(n)%tile(t)%row
             
             if(GLS_ctl%suscep_index(t).ge.&
                  GLS_data(m)%suscep_index_threshold) then 
                if(GLS_ctl%rf_acc1(t).ne.LIS_rc%udef.and.&
                     GLS_ctl%rf_acc1(t).ge.GLS_data(m)%rf_threshold1.and.&
                     GLS_ctl%acc1_alarm_status) then 
                   criteria1 = .true. 
                endif
                if(GLS_ctl%rf_acc2(t).ne.LIS_rc%udef.and.&
                     GLS_ctl%rf_acc2(t).ge.GLS_data(m)%rf_threshold2.and.&
                     GLS_ctl%acc2_alarm_status) then                    
                   criteria2 = .true. 
                endif
                if(GLS_ctl%rf_acc3(t).ne.LIS_rc%udef.and.&
                     GLS_ctl%rf_acc3(t).ge.GLS_data(m)%rf_threshold3.and.&
                     GLS_ctl%acc3_alarm_status) then 
                   criteria3  = .true. 
                endif
                gindex = LIS_domain(n)%gindex(c,r)
                
                lat = LIS_domain(n)%grid(gindex)%lat
                lon = LIS_domain(n)%grid(gindex)%lon
                
                ! Scenario for printng output
                
                if(gindex.ne.-1) then 
                   if(criteria1) then  
                      !    if(criteria1.and.criteria2) then 
                      data_status = ' T1'
                      !     elseif(criteria1) then 
                      !       data_status = '  C'
                      !     elseif(criteria2) then 
                      !        data_status = ' T1'
                      !     endif
                      GLS_ctl%ls_status(t) = 1.0
                      GLS_ctl%ls_accum_status(t) = 1.0
                      if(GLS_ctl%wform.eq.1) then 
                         write(GLS_ctl%out_ftn,fmt='(I4.4,3I2.2, F7.3,A1,F8.3,A1,F3.1,A1,F9.3,A3)') &
                              LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, &
                              lat, ' ',lon, ' ',GLS_ctl%suscep_index(t), ' ',&
                              GLS_ctl%rf_acc1(t), data_status
                      endif
                   endif
                   
                   if(criteria2) then 
                      !  if(criteria1.and.criteria3) then 
                      !     data_status = ' A2'
                      !  elseif(criteria1) then 
                      !     data_status = '  C'
                      !  elseif(criteria3) then 
                      data_status = ' T2'
                      !  endif
                      GLS_ctl%ls_status(t) = 1.0
                      GLS_ctl%ls_accum_status(t) = 1.0
                      if(GLS_ctl%wform.eq.1) then 
                         write(GLS_ctl%out_ftn,fmt='(I4.4,3I2.2, F7.3,A1,F8.3,A1,F3.1,A1,F9.3,A3)') &
                              LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, &
                              lat, ' ',lon, ' ',GLS_ctl%suscep_index(t), ' ',&
                              GLS_ctl%rf_acc2(t),data_status
                      endif
                   endif
                   
                   if(criteria3) then 
                      !   if(criteria1.and.criteria4) then 
                      !      data_status = ' A3'
                      !   elseif(criteria1) then 
                      !      data_status = '  C'
                      !   elseif(criteria4) then 
                      data_status = ' T3'
                      !   endif
                      GLS_ctl%ls_status(t) = 1.0
                      GLS_ctl%ls_accum_status(t) = 1.0
                      if(GLS_ctl%wform.eq.1) then 
                         write(GLS_ctl%out_ftn,fmt='(I4.4,3I2.2, F7.3,A1,F8.3,A1,F3.1,A1,F9.3,A3)') &
                              LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, &
                              lat, ' ',lon, ' ',GLS_ctl%suscep_index(t), ' ',&
                              GLS_ctl%rf_acc3(t),data_status
                      endif
                   endif
                   
                endif
             endif
          enddo
       enddo
    endif
  end subroutine GLS_run

!BOP
! !ROUTINE: GLS_output
! 
! !INTERFACE: 
  subroutine GLS_output(n)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_masterproc
    use LIS_logMod, only : LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_fileIOMod,  only : LIS_create_output_directory
    use LIS_historyMod, only : LIS_writevar_bin
! !ARGUMENTS: 
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
! 
!EOP
    integer           :: ftn 
    character(len=12) :: cdate1
    logical        :: alarmCheck
    character*100  :: filename
    integer        :: c,r
    integer        :: iret

#if 0 
    alarmCheck = ESMF_AlarmIsRinging(GLS_ctl%outAlarm,rc=iret)
    
    if(alarmCheck) then 
       call ESMF_AlarmRingerOff(GLS_ctl%outAlarm,rc=iret)

       if(LIS_masterproc) then 
          ftn = LIS_getNextUnitNumber()        

          write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
               LIS_rc%yr, LIS_rc%mo, &
               LIS_rc%da, LIS_rc%hr,LIS_rc%mn
          
          filename = trim(LIS_rc%odir)//&
               '/GLS/GLS_'//cdate1//'.1gs4r'
          call LIS_create_output_directory('GLS')
          open(ftn,file=trim(filename), form='unformatted')
       endif
       
       call LIS_writevar_bin(ftn,n,GLS_ctl%ls_status)

       if(LIS_masterproc) then 
          call LIS_releaseUnitNumber(ftn)
       endif

!reset       
       GLS_ctl%ls_status = 0.0

    endif
#endif
    
    if(LIS_rc%endtime.eq.1) then 
       
       if(LIS_masterproc) then 
          ftn = LIS_getNextUnitNumber()        

          write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
               LIS_rc%yr, LIS_rc%mo, &
               LIS_rc%da, LIS_rc%hr,LIS_rc%mn
          
          filename = trim(LIS_rc%odir)//&
               '/GLS/GLS_final_'//cdate1//'.1gs4r'
          call LIS_create_output_directory('GLS')
          open(ftn,file=trim(filename), form='unformatted')
       endif
       
       call LIS_writevar_bin(ftn,n,GLS_ctl%ls_accum_status)

       if(LIS_masterproc) then 
          call LIS_releaseUnitNumber(ftn)
       endif
    endif

  end subroutine GLS_output
!BOP
! 
! !ROUTINE: compute_rainf_accumulations
!  
! !INTERFACE: 
  subroutine compute_rainf_accumulations(n, acc_interval, &
       rc, rf_acc_sum)
! !USES: 
    use LIS_coreMod,   only : LIS_rc, LIS_domain
    use LIS_FORC_AttributesMod 
    use LIS_metforcingMod, only : LIS_FORC_State 
    use LIS_logMod,         only : LIS_verify

    implicit none
    integer   :: n
    integer   :: acc_interval
    integer   :: rc
    real      :: rf_acc_sum(LIS_rc%ntiles(n),acc_interval)
! 
! !DESCRIPTION: 
!  This subroutine computing running sum of precip based on the 
!  accumulation interval and the timestep 
!
!EOP

    type(ESMF_Field)    :: pcpField
    real, pointer       :: pcp(:)
    integer             :: c, r, t
    integer             :: tindex
    integer             :: nts, i
    integer             :: status

    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),&
         pcpField,rc=status)   
    call LIS_verify(status,'GLS_run: error getting Rainf')
    call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
    call LIS_verify(status,'GLS_run: error retrieving pcp')
    
    nts = acc_interval/24.0

    if(rc.lt.nts) then !looking for the first ring
       do i=rc+1,1,-1 
          do t=1,LIS_rc%ntiles(n)
             c = LIS_domain(n)%tile(t)%col
             r = LIS_domain(n)%tile(t)%row
             rf_acc_sum(t,i) = rf_acc_sum(t,i) + &
                  pcp(LIS_domain(n)%gindex(c,r))*LIS_rc%ts !convert to mm
          enddo
       enddo
    else
!       tindex = rc - (rc/nts)* nts  !integer arithmetic
       do i=nts,1,-1             
          do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
             c = LIS_domain(n)%tile(t)%col
             r = LIS_domain(n)%tile(t)%row
             rf_acc_sum(t,i) = rf_acc_sum(t,i) + &
                  pcp(LIS_domain(n)%gindex(c,r))*LIS_rc%ts !convert to mm
          enddo
       enddo
    endif
  end subroutine compute_rainf_accumulations
!BOP
! 
! !ROUTINE: GLS_final
! \label{GLS_final}
! 
! !INTERFACE:
  subroutine GLS_final()
! !USES:
!
! !DESCRIPTION:
!EOP

  end subroutine GLS_final
end module GLSMod
