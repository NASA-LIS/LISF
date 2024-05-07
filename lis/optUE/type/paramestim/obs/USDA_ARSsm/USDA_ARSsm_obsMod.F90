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
! !MODULE: USDA_ARSsm_obsMod
! 
! !DESCRIPTION: 
!   
! !REVISION HISTORY: 
!  11 Jul 11    Ken Harrison;   Initial Specification
! 
module USDA_ARSsm_obsMod
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: USDA_ARSsm_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: USDA_ARSsm_obs_struc

  type, public :: soilmoisture_data_dec
     integer             :: numrecords
     integer             :: nobs_threshold    
     real, allocatable       :: soilmoisture(:) 
     integer, allocatable       :: nobs(:) 
     type(ESMF_Time)     :: start_time, end_time
     integer             :: rec_len_days
  end type soilmoisture_data_dec

  type(soilmoisture_data_dec), allocatable :: USDA_ARSsm_obs_struc(:)

contains
!BOP
! 
! !ROUTINE: USDA_ARSsm_obs_setup
! \label{USDA_ARSsm_obs_setup}
! 
! !INTERFACE: 
  subroutine USDA_ARSsm_obs_setup(Obs_State)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, LIS_vecGrid
    use LIS_timeMgrMod, only : LIS_calendar
    use LIS_logMod, only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obs_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!
!   The arguments are: 
!   \begin{description}
!    \item[Obj\_Space]   observation/Objective space object 
!   \end{description}
!EOP
    integer                   ::  n 
    integer                   ::  status
    integer                   ::  i 
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field)          ::  obsField
    character(len=LIS_CONST_PATH_LEN) ::  soilmoistureobsdir
    character(len=LIS_CONST_PATH_LEN) ::  soilmoistureobsmaskdir
    character*100             ::  vname
    character(len=LIS_CONST_PATH_LEN) ::  obsAttribFile(LIS_rc%nnest)
    integer                   ::  ftn
    integer                   ::  numrecords, ios
    integer                   ::  yr, jday
    integer                   ::  mo, da, hr, mn
    integer                   ::  yr_curr, mo_curr, da_curr, hr_curr, mn_curr
    type(ESMF_TimeInterval)   ::  dt,esmf_jday,  esmf_one_day, delta_t
    type(ESMF_Time)           ::  temp_time, t_temp, t_curr, t_curr_a, t_curr_d, t_prev_a, t_prev_d, t_next_a, t_next_d
    type(ESMF_TimeInterval)   ::  dt_curr_a, dt_curr_d, dt_prev_a, dt_prev_d, dt_next_a, dt_next_d,dt_a, dt_d, dt_min
    real                      ::  gridDesci(LIS_rc%nnest,50)
    real                      ::  sm_ave
    real                      ::  sm_std
    integer                   ::  day_index
    integer                   ::  k, p, f
    integer                   ::  rec_len_days
    integer                   ::  nobs

    n=1

    allocate(USDA_ARSsm_obs_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,& 
         rc=status)
    call LIS_verify(status)

    !Read lis.config entries
    call ESMF_ConfigFindLabel(LIS_config,'USDA ARS Soilmoisture Obs data directory:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,soilmoistureobsdir,&
            rc=status)
       call LIS_verify(status, 'Err: USDA ARS Soilmoisture Obs data directory: not defined')

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            soilmoistureobsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),'Data Update Status',&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'USDA ARS Soilmoisture observations attributes file:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsAttribFile(n),rc=status)
       call LIS_verify(status, 'Err: USDA ARS Soilmoisture observations attributes file: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'USDA ARS number of observations threshold:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,USDA_ARSsm_obs_struc(n)%nobs_threshold,&
            rc=status)
       call LIS_verify(status,'USDA ARS number of observations threshold: not defined')
    enddo

    !Determine indexing into arrays, 2 daily entries for asc- and desc-ending passes
    do n=1,LIS_rc%nnest
       call ESMF_TimeSet( &
            USDA_ARSsm_obs_struc(n)%start_time, &
            yy=2002, &
            mm=1,    &
            dd=1,    &
            h=0,     &
            m=0,     &
            calendar=LIS_calendar,rc=status)

       call ESMF_TimeSet( &
            USDA_ARSsm_obs_struc(n)%end_time, &
            yy=2011, &
            mm=12,    &
            dd=31,    &
            h=23,     &
            m=30,calendar=LIS_calendar,rc=status)

       dt=USDA_ARSsm_obs_struc(n)%end_time-USDA_ARSsm_obs_struc(n)%start_time

       !Retrieve number of days between start and stop
       call ESMF_TimeIntervalGet(dt,d=rec_len_days,rc=status)

       !initialize all prior to filling in
       allocate(USDA_ARSsm_obs_struc(n)%soilmoisture(rec_len_days))
       allocate(USDA_ARSsm_obs_struc(n)%nobs(rec_len_days))
       USDA_ARSsm_obs_struc(n)%soilmoisture = LIS_rc%udef
       USDA_ARSsm_obs_struc(n)%nobs = LIS_rc%udef
    enddo

    do n=1,LIS_rc%nnest
       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(obsAttribFile(n)),status='old')
       read(ftn,*)
       read(ftn,fmt='(a100)') vname
       obsField = ESMF_FieldCreate(arrayspec=realarrspec, &
            grid=LIS_vecGrid(n), &
            name=trim(vname), rc=status)
       call LIS_verify(status)

       call ESMF_StateAdd(Obs_State(n),(/obsField/),rc=status)
       call LIS_verify(status)

       call LIS_releaseUnitNumber(ftn)

       !Set interval of one day
       call ESMF_TimeIntervalSet(esmf_one_day,d=1,rc=status)

       !read to store the values
       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(soilmoistureobsdir),status='old')
       read(ftn,*)

       ios = 0 
       numrecords = 0 
       do while(ios.eq.0)
          numrecords = numrecords + 1
          read(ftn,*,iostat=ios) yr, jday, hr, mn, nobs, sm_ave, sm_std
          if(nobs.ge.USDA_ARSsm_obs_struc(n)%nobs_threshold) then
             if(hr.eq.0.and.mn.eq.0) then
                
                !Get current time into ESMF state
                !As given in 'julian' day, get interval as esmf interval
                call ESMF_TimeIntervalSet(esmf_jday,d=jday,rc=status)
                
                !Add interval to same year, hr, min but of day Jan 1
                call ESMF_TimeSet(temp_time, yy=yr, mm=1, dd=1,h=0,&
                     m=0,calendar=LIS_calendar,rc=status)
                !Put together, substracting by 1 day (consider Jan 1)
                t_curr = temp_time + esmf_jday - esmf_one_day
                
                
                
                !Finally get record time as esmf time
                call ESMF_TimeGet(t_curr, yy=yr_curr,mm=mo_curr,dd=da_curr,&
                     h=hr_curr,m=mn_curr,calendar=LIS_calendar,rc=status)
                call LIS_verify(status, 'error in timeget: readUSDA_ARSsm_obsMod.F90')


                !Get index for soilmoisture array
                delta_t = t_curr - USDA_ARSsm_obs_struc(n)%start_time
                call ESMF_TimeIntervalGet(delta_t, d=day_index, &
                     calendar=LIS_calendar,rc=status)
                call LIS_verify(status, 'error in timeget: USDA_ARSsm_obsMod.F90')


                !as fortran arrays are zero-based...
                day_index=day_index+1

                USDA_ARSsm_obs_struc(n)%soilmoisture(day_index)=sm_ave
                USDA_ARSsm_obs_struc(n)%nobs(day_index)=nobs
             endif
          endif
       enddo
!       Open(111, file="test.txt")
!       do i=1, rec_len_days
!          write(111, '(2I6, 9F14.4)') i, USDA_ARSsm_obs_struc(n)%nobs(i), &
!               USDA_ARSsm_obs_struc(n)%soilmoisture(i)
!       end do
!       close(111)
    enddo
    
    write(LIS_logunit,*) 'Created the States to hold the USDA ARS soilmoisture observations'
 
  end subroutine USDA_ARSsm_obs_setup

end module USDA_ARSsm_obsMod
