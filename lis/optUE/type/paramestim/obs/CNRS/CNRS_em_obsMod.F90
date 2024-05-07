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
! !MODULE: CNRS_em_obsMod
! 
! !DESCRIPTION: 
!   
! !REVISION HISTORY: 
!  11 Jul 11    Ken Harrison;   Initial Specification
! 
module CNRS_em_obsMod
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: CNRS_em_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: CNRS_em_obs_struc

  type, public :: emissivity_data_dec
     integer             :: numrecords
     !     integer             :: nobs_threshold    
     !     real, allocatable       :: lat(:)
     !     real, allocatable       :: lon(:)
     real, allocatable       :: emissivity(:,:,:,:) !numrecords, 7 channels, 2 passes, 3 platforms
     real, allocatable       :: temperature(:,:,:) ! numrecords, 2 passes, 3 platforms
     integer, allocatable       :: nobs(:,:,:) ! numrecords, 2 passes, 3 platforms
     real, allocatable       :: vap(:,:,:) ! numrecords, 2 passes, 3 platforms
     type(ESMF_Time)     :: start_time, end_time
     integer             :: overpass_hr_a
     integer             :: overpass_hr_d
     integer             :: mask_hr_a_lower  !pick closest to overpass hr, 0,3,...
     integer             :: mask_hr_a_upper  !also note that covers +/- 1.5hrs
     integer             :: mask_hr_d_lower  !also note that platforms have
     integer             :: mask_hr_d_upper  !slightly different overpass times
     integer             :: cloud_threshold
     integer             :: index_ascend
     integer             :: index_descend
     integer             :: rec_len_days
  end type emissivity_data_dec

  type(emissivity_data_dec), allocatable :: CNRS_em_obs_struc(:)

contains
!BOP
! 
! !ROUTINE: CNRS_em_obs_setup
! \label{CNRS_em_obs_setup}
! 
! !INTERFACE: 
  subroutine CNRS_em_obs_setup(Obs_State)
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
    integer, parameter        ::  numchannels=7
    integer, parameter        ::  numplatforms=3
    integer, parameter        ::  numdailypasses=2
    integer                   ::  status
    integer                   ::  i 
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field)          ::  obsField
    character(len=LIS_CONST_PATH_LEN) ::  emissivityobsdir
    character(len=LIS_CONST_PATH_LEN) ::  emissivityobsmaskdir
    character*100             ::  vname
    character(len=LIS_CONST_PATH_LEN) ::  obsAttribFile(LIS_rc%nnest)
    integer                   ::  ftn
    integer                   ::  numrecords, ios
    integer                   ::  jday, yr, lat, lon
    integer                   ::  id(numplatforms), igmt(numplatforms),nobs(numplatforms) ! for three platforms 13,14, 15
    real                      ::  tsurf(numplatforms), vap(numplatforms)
    integer                   ::  gmt_hr(numplatforms), gmt_mn(numplatforms) ! three platforms
    integer                   ::  mo, da, hr, mn
    integer                   ::  yr_curr, mo_curr, da_curr, hr_curr, mn_curr
    type(ESMF_TimeInterval)   ::  dt,esmf_jday,  esmf_one_day, delta_t
    type(ESMF_Time)           ::  temp_time, t_temp, t_curr, t_curr_a, t_curr_d, t_prev_a, t_prev_d, t_next_a, t_next_d
    type(ESMF_Time)           ::  masktime
    type(ESMF_TimeInterval)   ::  dt_curr_a, dt_curr_d, dt_prev_a, dt_prev_d, dt_next_a, dt_next_d,dt_a, dt_d, dt_min
    real                      ::  gridDesci(LIS_rc%nnest,50)
    integer, parameter        ::  i_undef=-999
!    integer                   ::  nobs_threshold
    real                      ::  em(numchannels,numplatforms)  ! for reading: 7 channels, 3 platforms
!    real, allocatable             ::  emiss(:,:,:,:) !for temp storage for averaging,day-channel-pltform-pass(a,d)
!    real                      ::  em_sum,em_avg
    logical                   ::  is_ascend_pass
    integer                   ::  day_index
    integer                   ::  k, p, f
    integer                   ::  rec_len_days
!    integer                   :: platform_count
    integer                   :: ipass
    logical                   :: apply_mask
    integer                   :: cloud
    allocate(CNRS_em_obs_struc(LIS_rc%nnest))


!!! 19V, 19H, 22V, [**NO 22H**], 37V, 37H,85V,85H
    do n=1,LIS_rc%nnest
       CNRS_em_obs_struc(n)%index_ascend  = 1
       CNRS_em_obs_struc(n)%index_descend = 2
    enddo
!!$    call ESMF_ArraySpecSet(realarrspec,rank=2,typekind=ESMF_TYPEKIND_R4,& !Note: rank=2
!!$         rc=status)
!!$    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,& !Note: rank=1
         rc=status)
    call LIS_verify(status)


    !Read lis.config entries
    call ESMF_ConfigFindLabel(LIS_config,'CNRS Emissivity Obs data directory:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,emissivityobsdir,&
            rc=status)
       call LIS_verify(status, 'Err: CNRS Emissivity Obs data directory: not defined')

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            emissivityobsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),'Data Update Status',&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'CNRS Emissivity Obs mask directory:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,emissivityobsmaskdir,&
            rc=status)
       call LIS_verify(status, 'Err: CNRS Emissivity Obs mask directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'CNRS Emissivity observations attributes file:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsAttribFile(n),rc=status)
       call LIS_verify(status, 'Err: CNRS Emissivity observations attributes file: not defined')
    enddo

!    call ESMF_ConfigFindLabel(LIS_config,'CNRS number of observations threshold:',&
!         rc=status)
!    do n=1,LIS_rc%nnest
!       call ESMF_ConfigGetAttribute(LIS_config,CNRS_em_obs_struc(n)%nobs_threshold,&
!            rc=status)
!       call LIS_verify(status,'CNRS number of observations threshold: not defined')
!    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Overpass hr descending:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,CNRS_em_obs_struc(n)%overpass_hr_d,&
            rc=status)
       call LIS_verify(status,'Overpass hr descending: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Overpass hr ascending:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,CNRS_em_obs_struc(n)%overpass_hr_a,&
            rc=status)
       call LIS_verify(status,'Overpass hr ascending: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Mask hr ascending lower:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,CNRS_em_obs_struc(n)%mask_hr_a_lower,&
            rc=status)
       call LIS_verify(status,'Mask hr ascending lower: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Mask hr ascending upper:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,CNRS_em_obs_struc(n)%mask_hr_a_upper,&
            rc=status)
       call LIS_verify(status,'Mask hr ascending upper: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Mask hr descending lower:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,CNRS_em_obs_struc(n)%mask_hr_d_lower,&
            rc=status)
       call LIS_verify(status,'Mask hr descending lower: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Mask hr descending upper:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,CNRS_em_obs_struc(n)%mask_hr_d_upper,&
            rc=status)
       call LIS_verify(status,'Mask hr descending upper: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Mask cloud threshold:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,CNRS_em_obs_struc(n)%cloud_threshold,&
            rc=status)
       call LIS_verify(status,'Mask cloud threshold: not defined')
    enddo

    !Determine indexing into arrays, 2 daily entries for asc- and desc-ending passes
    do n=1,LIS_rc%nnest
       call ESMF_TimeSet( &
            CNRS_em_obs_struc(n)%start_time, &
            yy=2004, &
            mm=7,    &
            dd=1,    &
            h=0,     &
            m=0,     &
            calendar=LIS_calendar,rc=status)

       call ESMF_TimeSet( &
            CNRS_em_obs_struc(n)%end_time, &
            yy=2007, &
            mm=7,    &
            dd=1,    &
            h=0,     &
            m=0,calendar=LIS_calendar,rc=status)

       dt=CNRS_em_obs_struc(n)%end_time-CNRS_em_obs_struc(n)%start_time

       !Retrieve number of days between start and stop
       call ESMF_TimeIntervalGet(dt,d=rec_len_days,rc=status)

       !Allocate dynamic arrays
       !allocate(CNRS_em_obs_struc(n)%time(numrecords))
       !allocate(CNRS_em_obs_struc(n)%lat(numrecords))
       !allocate(CNRS_em_obs_struc(n)%lon(numrecords))

!       allocate(emiss(rec_len_days,numchannels,numplatforms,numdailypasses))  !7 channels, 3 platforms, 2 passes

       !initialize all emissivities prior to filling in
       allocate(CNRS_em_obs_struc(n)%emissivity(rec_len_days,numchannels,numdailypasses, numplatforms))  !numrecords, numchannels, numpasses, numplatforms
       allocate(CNRS_em_obs_struc(n)%nobs(rec_len_days,numdailypasses, numplatforms))                  !numrecords, numpasses, numplatforms
       allocate(CNRS_em_obs_struc(n)%temperature(rec_len_days,numdailypasses, numplatforms))                  !numrecords, numpasses, numplatforms
       allocate(CNRS_em_obs_struc(n)%vap(rec_len_days,numdailypasses, numplatforms))                  !numrecords, numpasses, numplatforms
       CNRS_em_obs_struc(n)%emissivity = LIS_rc%udef
       CNRS_em_obs_struc(n)%nobs = LIS_rc%udef
       CNRS_em_obs_struc(n)%temperature = LIS_rc%udef
       CNRS_em_obs_struc(n)%vap = LIS_rc%udef
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

!!$       !initialize all emissivities to undefined
!!$       emiss = LIS_rc%udef

       !read to store the values
       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(emissivityobsdir),status='old')

       ios = 0 
       numrecords = 0 
       do while(ios.eq.0)
          numrecords = numrecords + 1
          read(ftn,*,iostat=ios) jday, yr, lat, lon, &
               id(1), igmt(1), nobs(1), em(:,1), tsurf(1), vap(1), &
               id(2), igmt(2), nobs(2), em(:,2), tsurf(2), vap(2), &
               id(3), igmt(3), nobs(3), em(:,3), tsurf(3), vap(3)
          ! Convert to 'normal' units and 'normal' time
          do i=1,numplatforms
             if (igmt(i).ne.i_undef) then
                tsurf(i)=tsurf(i)/10.0
                em(:,i)=em(:,i)/1000.0
                vap(i)=vap(i)/10.0
                gmt_hr(i)=int(igmt(i)/100)   !the CNRS time is, eg, 2252, meaning 22:52
                gmt_mn(i)=igmt(i)-gmt_hr(i)*100  
             else
                tsurf(i) =  LIS_rc%udef
                em(:,i)  =  LIS_rc%udef
                vap(i)   =  LIS_rc%udef
                gmt_hr(i) = LIS_rc%udef
                gmt_mn(i) = LIS_rc%udef
             end if
          end do

          !Set interval of one day
          call ESMF_TimeIntervalSet(esmf_one_day,d=1,rc=status)

          !additional complexity with the asc/desc diff's w/in record.  sometimes 2 ascends in one day
          !e.g., 00:03 followed by 23:59.  w/ current approach, assign to assumed overpass times
          !w/ current approach,  assign to the proper (assumed) day and pass time
          do i=1,numplatforms
             if(igmt(i).ne.i_undef) then

                !Get current time into ESMF state
                !As given in 'julian' day, get interval as esmf interval
                call ESMF_TimeIntervalSet(esmf_jday,d=jday,rc=status)

                !Add interval to same year, hr, min but of day Jan 1
                call ESMF_TimeSet(temp_time, yy=yr, mm=1, dd=1,h=gmt_hr(i),&
                     m=gmt_mn(i),calendar=LIS_calendar,rc=status)
                !Put together, substracting by 1 day (consider Jan 1)
                t_curr = temp_time + esmf_jday - esmf_one_day

                !Finally get record time as esmf time
                call ESMF_TimeGet(t_curr, yy=yr_curr,mm=mo_curr,dd=da_curr,&
                     h=hr_curr,m=mn_curr,calendar=LIS_calendar,rc=status)
                call LIS_verify(status, 'error in timeget: readCNRS_em_obsMod.F90')

                !Set time for same day but with overpass hr (for both passes)
                call ESMF_TimeSet(t_curr_a, yy=yr_curr, mm=mo_curr, dd=da_curr,&
                     h=CNRS_em_obs_struc(n)%overpass_hr_a, &
                     m=0,calendar=LIS_calendar,rc=status)

                call ESMF_TimeSet(t_curr_d, yy=yr_curr, mm=mo_curr, dd=da_curr,&
                     h=CNRS_em_obs_struc(n)%overpass_hr_d, &
                     m=0,calendar=LIS_calendar,rc=status)

                !Find overpass time for previous day and next day
                t_prev_a = t_curr_a - esmf_one_day
                t_prev_d = t_curr_d - esmf_one_day

                t_next_a = t_curr_a + esmf_one_day
                t_next_d = t_curr_d + esmf_one_day

                !Find which of curr, prev, or next overpass times and desc or ascend is closest to t_curr
                
                dt_curr_a = t_curr - t_curr_a
                dt_curr_a = ESMF_TimeIntervalAbsValue(dt_curr_a)

                dt_curr_d = t_curr - t_curr_d
                dt_curr_d = ESMF_TimeIntervalAbsValue(dt_curr_d)

                dt_prev_a = t_curr - t_prev_a
                dt_prev_a = ESMF_TimeIntervalAbsValue(dt_prev_a)

                dt_prev_d = t_curr - t_prev_d
                dt_prev_d = ESMF_TimeIntervalAbsValue(dt_prev_d)

                dt_next_a = t_curr - t_next_a
                dt_next_a = ESMF_TimeIntervalAbsValue(dt_next_a)

                dt_next_d = t_curr - t_next_d
                dt_next_d = ESMF_TimeIntervalAbsValue(dt_next_d)

                !initialize
                dt_min = dt_curr_a
                t_temp = t_curr_a
                is_ascend_pass = .true.

                if (dt_prev_a < dt_min) then
                   dt_min = dt_prev_a
                   t_temp = t_prev_a
                   is_ascend_pass = .true.
                endif

                if (dt_next_a < dt_min) then
                   dt_min = dt_next_a
                   t_temp = t_next_a
                   is_ascend_pass = .true.
                endif

                if (dt_curr_d < dt_min) then
                   dt_min = dt_curr_d
                   t_temp = t_curr_d
                   is_ascend_pass = .false.
                endif

                if (dt_prev_d < dt_min) then
                   dt_min = dt_prev_d
                   t_temp = t_prev_d
                   is_ascend_pass = .false.
                endif

                if (dt_next_d < dt_min) then
                   dt_min = dt_next_d
                   t_temp = t_next_d
                   is_ascend_pass = .false.
                endif

                !Get index for emissivity array
                delta_t = t_temp - CNRS_em_obs_struc(n)%start_time
                call ESMF_TimeIntervalGet(delta_t, d=day_index, &
                     calendar=LIS_calendar,rc=status)
                call LIS_verify(status, 'error in timeget: CNRS_em_obsMod.F90')

                !as fortran arrays are zero-based...
                day_index=day_index+1

                if (is_ascend_pass) then
                   ipass=CNRS_em_obs_struc(n)%index_ascend
!!$                   emiss(day_index,:,i,CNRS_em_obs_struc(n)%index_ascend)=em(:,i)
                else
!!$                   emiss(day_index,:,i,CNRS_em_obs_struc(n)%index_descend)=em(:,i)
                   ipass=CNRS_em_obs_struc(n)%index_descend
                endif
                CNRS_em_obs_struc(n)%emissivity(day_index,:,ipass,i)=em(:,i)
                CNRS_em_obs_struc(n)%nobs(day_index,ipass,i)=nobs(i)
                CNRS_em_obs_struc(n)%temperature(day_index,ipass,i)=tsurf(i)
                CNRS_em_obs_struc(n)%vap(day_index,ipass,i)=vap(i)
             endif
          enddo
       enddo



       ! Fill in with averages
!!$       call LIS_releaseUnitNumber(ftn)

!!$       do k=1,rec_len_days
!!$          do p=1,numdailypasses
!!$             em_sum=0
!!$             do f=1,numchannels
!!$                em_sum=0
!!$                platform_count=0
!!$                do i=1,numplatforms  !platform F13, 14, 15
!!$                   if (emiss(k,f,i,p).ne.LIS_rc%udef) then
!!$                      ! if(nobs(i).ge.nobs_threshold) then
!!$                      platform_count=platform_count+1
!!$                      em_sum=em_sum+emiss(k,f,i,p)
!!$                   endif
!!$                enddo
!!$                if(platform_count.ne.0) then
!!$                   em_avg=em_sum/dble(platform_count)
!!$                   CNRS_em_obs_struc(n)%emissivity(k,f,p) = em_avg
!!$                else
!!$                   CNRS_em_obs_struc(n)%emissivity(k,f,p) = LIS_rc%udef
!!$                endif
!!$             enddo
!!$          enddo
!!$       enddo
    enddo
    
! read mask file
! http://isccp.giss.nasa.gov/pub/documents/dn-022.pdf for definition
! of "three-hourly"; within 1.5 hr of stated time.
       !read to store the values

    n=1

    Open(111, file="F13-unfiltered-AM.txt")
    Open(112, file="F14-unfiltered-AM.txt")
    Open(113, file="F15-unfiltered-AM.txt")
    Open(211, file="F13-unfiltered-PM.txt")
    Open(212, file="F14-unfiltered-PM.txt")
    Open(213, file="F15-unfiltered-PM.txt")
    do i=1, rec_len_days
       write(111, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,1,1), &
            CNRS_em_obs_struc(n)%temperature(i, 1, 1), &
            CNRS_em_obs_struc(n)%vap(i, 1, 1), &
            CNRS_em_obs_struc(n)%emissivity(i,:,1,1)
       write(112, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,1,2), &
            CNRS_em_obs_struc(n)%temperature(i, 1, 2), &
            CNRS_em_obs_struc(n)%vap(i, 1, 2), &
            CNRS_em_obs_struc(n)%emissivity(i,:,1,2)
       write(113, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,1,3), &
            CNRS_em_obs_struc(n)%temperature(i, 1, 3), &
            CNRS_em_obs_struc(n)%vap(i, 1, 3), &
            CNRS_em_obs_struc(n)%emissivity(i,:,1,3)
       write(211, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,2,1), &
            CNRS_em_obs_struc(n)%temperature(i, 2, 1), &
            CNRS_em_obs_struc(n)%vap(i, 2, 1), &
            CNRS_em_obs_struc(n)%emissivity(i,:,2,1)
       write(212, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,2,2), &
            CNRS_em_obs_struc(n)%temperature(i, 2, 2), &
            CNRS_em_obs_struc(n)%vap(i, 2, 2), &
            CNRS_em_obs_struc(n)%emissivity(i,:,2,2)
       write(213, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,2,3), &
            CNRS_em_obs_struc(n)%temperature(i, 2, 3), &
            CNRS_em_obs_struc(n)%vap(i, 2, 3), &
            CNRS_em_obs_struc(n)%emissivity(i,:,2,3)
    end do
    close(111)
    close(112)
    close(113)
    close(211)
    close(212)
    close(213)

    ftn=LIS_getNextUnitNumber()
    open(ftn,file=trim(emissivityobsmaskdir),status='old')
    read(ftn,*)
    ios = 0 
    numrecords = 0 
    do while(ios.eq.0)
       numrecords = numrecords + 1
       read(ftn,*,iostat=ios) yr, mo, da, hr, cloud
       call ESMF_TimeSet(masktime, yy=yr,&
            mm=mo, dd=da,h=hr,m=0,&
            calendar = LIS_calendar, rc=status)
       call LIS_verify(status, 'ESMF_timeset error in CNRS_em_obsMod.F90')
       
       !determine if ascending/descending & if to test for clouds
       apply_mask=.false.
       if ( (hr.eq.CNRS_em_obs_struc(n)%mask_hr_a_lower) .or. (hr.eq.CNRS_em_obs_struc(n)%mask_hr_a_upper) ) then
          ipass=CNRS_em_obs_struc(n)%index_ascend
          apply_mask=.true.
       elseif ( (hr.eq.CNRS_em_obs_struc(n)%mask_hr_d_lower) .or. (hr.eq.CNRS_em_obs_struc(n)%mask_hr_d_upper) ) then
          ipass=CNRS_em_obs_struc(n)%index_descend
          apply_mask=.true.
       end if
       
!!!!!!!!!!changing strategy to focus on nobs!!!!!!!
!!!!!!!!!!                                  !!!!!!!
apply_mask=.false.
!!!!!!!!!!
!!!!!!!!!!

       if(apply_mask) then
          !Get day index for emissivity array
          delta_t = masktime - CNRS_em_obs_struc(n)%start_time
          call ESMF_TimeIntervalGet(delta_t, d=day_index, &
               calendar=LIS_calendar,rc=status)
          call LIS_verify(status, 'error in timeget: CNRS_em_obsMod.F90')
          
          !as fortran arrays are zero-based...
          day_index=day_index+1
          
          if(cloud.gt.CNRS_em_obs_struc(n)%cloud_threshold) then
             CNRS_em_obs_struc(n)%emissivity(day_index,:,ipass,:) = LIS_rc%udef
             CNRS_em_obs_struc(n)%nobs(day_index,ipass,:)         = LIS_rc%udef
             CNRS_em_obs_struc(n)%temperature(day_index,ipass,:)  = LIS_rc%udef
             CNRS_em_obs_struc(n)%vap(day_index,ipass,:)          = LIS_rc%udef
          endif
       endif
    enddo

    Open(111, file="F13-filtered-AM.txt")
    Open(112, file="F14-filtered-AM.txt")
    Open(113, file="F15-filtered-AM.txt")
    Open(211, file="F13-filtered-PM.txt")
    Open(212, file="F14-filtered-PM.txt")
    Open(213, file="F15-filtered-PM.txt")
    do i=1, rec_len_days
       write(111, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,1,1), &
            CNRS_em_obs_struc(n)%temperature(i, 1, 1), &
            CNRS_em_obs_struc(n)%vap(i, 1, 1), &
            CNRS_em_obs_struc(n)%emissivity(i,:,1,1)
       write(112, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,1,2), &
            CNRS_em_obs_struc(n)%temperature(i, 1, 2), &
            CNRS_em_obs_struc(n)%vap(i, 1, 2), &
            CNRS_em_obs_struc(n)%emissivity(i,:,1,2)
       write(113, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,1,3), &
            CNRS_em_obs_struc(n)%temperature(i, 1, 3), &
            CNRS_em_obs_struc(n)%vap(i, 1, 3), &
            CNRS_em_obs_struc(n)%emissivity(i,:,1,3)
       write(211, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,2,1), &
            CNRS_em_obs_struc(n)%temperature(i, 2, 1), &
            CNRS_em_obs_struc(n)%vap(i, 2, 1), &
            CNRS_em_obs_struc(n)%emissivity(i,:,2,1)
       write(212, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,2,2), &
            CNRS_em_obs_struc(n)%temperature(i, 2, 2), &
            CNRS_em_obs_struc(n)%vap(i, 2, 2), &
            CNRS_em_obs_struc(n)%emissivity(i,:,2,2)
       write(213, '(2I6, 9F14.4)') i, CNRS_em_obs_struc(n)%nobs(i,2,3), &
            CNRS_em_obs_struc(n)%temperature(i, 2, 3), &
            CNRS_em_obs_struc(n)%vap(i, 2, 3), &
            CNRS_em_obs_struc(n)%emissivity(i,:,2,3)
    end do
    close(111)
    close(112)
    close(113)
    close(211)
    close(212)
    close(213)
  

!do i=1,7
!   print*,
!enddo
    write(LIS_logunit,*) 'Created the States to hold the CNRS emissivity observations'
 
  end subroutine CNRS_em_obs_setup

end module CNRS_em_obsMod
