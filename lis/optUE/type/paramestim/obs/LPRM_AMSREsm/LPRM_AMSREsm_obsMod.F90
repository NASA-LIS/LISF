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
! !MODULE: LPRM_AMSREsm_obsMod
! 
! !DESCRIPTION: 
!   
! !REVISION HISTORY: 
!  11 Jul 11    Ken Harrison;   Initial Specification
! 
module LPRM_AMSREsm_obsMod
! !USES: 
  use ESMF
  use map_utils
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LPRM_AMSREsm_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: LPRM_AMSREsm_obs_struc

  type, public :: lprm_sm_data_dec
     logical                :: startMode
     integer                :: nc
     integer                :: nr
     real,     allocatable      :: smobs(:,:)
     real,     allocatable      :: smtime(:,:)

     integer                :: lprmnc, lprmnr
     type(proj_info)        :: lprmproj
     integer                :: rawdata
     integer, allocatable       :: n11(:)
     integer, allocatable       :: n12(:)
     integer, allocatable       :: n21(:)
     integer, allocatable       :: n22(:)
     real,    allocatable       :: w11(:)
     real,    allocatable       :: w12(:)
     real,    allocatable       :: w21(:)
     real,    allocatable       :: w22(:)

!!!     integer             :: numrecords
!!!     integer             :: nobs_threshold    
!!!     real, allocatable       :: lprm_sm(:,:,:,:,:,:) !row, col, numrecords, 6 freqs, 2 polarizations, 2 passes
!!!     real, allocatable       :: mpdi(:,:,:,:,:) !row, col, numrecords, 6 freqs, 2 passes
!!!     integer, allocatable    :: nobs(:,:,:,:) ! row, col, numrecords, 2 passes
!!!     type(ESMF_Time)     :: start_time, end_time
!!!     integer             :: overpass_hr_a
!!!     integer             :: overpass_hr_d
!!!     integer             :: mask_hr_a_lower
!!!     integer             :: mask_hr_d_lower
!!!     integer             :: mask_hr_a_upper
!!!     integer             :: mask_hr_d_upper
!!!     integer             :: index_ascend
!!!     integer             :: index_descend
!!!     integer             :: rec_len_days
  end type lprm_sm_data_dec

  type(lprm_sm_data_dec), allocatable :: LPRM_AMSREsm_obs_struc(:)

contains
!BOP
! 
! !ROUTINE: LPRM_AMSREsm_obs_setup
! \label{LPRM_AMSREsm_obs_setup}
! 
! !INTERFACE: 
  subroutine LPRM_AMSREsm_obs_setup(Obs_State)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod, only : LIS_logunit, LIS_verify, & 
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use map_utils, only: map_init, map_set, latlon_to_ij

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
!!!    integer                   ::  n 
!!!    integer, parameter        ::  numfreqs=6
!!!    integer, parameter        ::  numchannels=12
!!!    integer, parameter        ::  numpolarizations=2
!!!    integer, parameter        ::  numdailypasses=2
!!!    integer                   ::  status
!!!    type(ESMF_ArraySpec)      ::  realarrspec
!!!    type(ESMF_Field)          ::  obsField
!!!    character*100             ::  lprm_smobsdir
!!!    character*100             ::  lprm_smobsmaskdir
!!!    character*100             ::  vname
!!!    character*100             ::  obsAttribFile(LIS_rc%nnest)
!!!    integer                   ::  ftn
!!!    integer                   ::  numrecords, ios
!!!    integer                   ::  yr
!!!    real                      ::  lat, lon
!!!    real                      ::  LWP, TPW, tsurf
!!!    real                      ::  tb_temp(numchannels)
!!!    real                      ::  em_temp(numchannels)
!!!    real, allocatable             ::  tb(:,:,:)
!!!    real, allocatable             ::  em(:,:,:)
!!!    real, allocatable             ::  mpdi(:,:)
!!!    integer                   ::  mo, da, hr, mn
!!!    integer                   ::  yr_curr, mo_curr, da_curr, hr_curr, mn_curr
!!!    type(ESMF_TimeInterval)   ::  dt,esmf_jday,  esmf_one_day, delta_t
!!!    type(ESMF_Time)           ::  temp_time, t_temp, t_curr, t_curr_a, t_curr_d, t_prev_a, t_prev_d, t_next_a, t_next_d
!!!    type(ESMF_Time)           ::  masktime
!!!    type(ESMF_TimeInterval)   ::  dt_curr_a, dt_curr_d, dt_prev_a, dt_prev_d, dt_next_a, dt_next_d,dt_a, dt_d, dt_min
!!!!    real                      ::  gridDesci(LIS_rc%nnest,50)
!!!    integer, parameter        ::  i_undef=-999
!!!!    integer                   ::  nobs_threshold
!!!    logical                   ::  is_ascend_pass
!!!    integer                   ::  day_index
!!!    integer                   ::  k
!!!    integer                   ::  rec_len_days
!!!    integer                   ::  ipass
!!!    integer                   ::  col, row !lat lon indices
!!!    real                      ::  col_real, row_real
!!!    integer                   ::  c !channel counter
!!!    integer                   ::  f !freq counter
!!!    integer                   ::  p !polarization counter
!!!!    real                      ::  LPRM_AMSREsm_gridDesc(6)
!!!    real, allocatable       :: sum_em(:,:,:,:,:,:) !row, col, numrecords, 6 freqs, 2 polarizations, 2 passes
!!!    real, allocatable       :: sum_em_mpdi(:,:,:,:,:) !row, col, numrecords, 6 freqs, 2 passes
!!!    integer, allocatable    :: sum_obs(:,:,:,:) ! row, col, numrecords, 2 passes
!!!    
!!!    integer, parameter :: maxd = 4000  ! max number of days for daily overpass 
!!!    integer, parameter :: maxm = 108     ! max number of months                
!!!    
!!!    character*200 datfile(maxm)
!!!    
!!!    ! time management                                                          
!!!    character*6 cym                                                  
!!!    integer :: im, y0, m0                           
    
    integer                ::  n,i,t,kk
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  lprm_smobsdir
    character*100          ::  temp
    character(len=LIS_CONST_PATH_LEN) ::  obsAttribFile(LIS_rc%nnest)
    real,  allocatable         ::  obsstd(:)
    character*1            ::  vid(2)
    character*100          ::  vname
    real                   :: gridDesci(50)

    allocate(LPRM_AMSREsm_obs_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,& !Note: rank=1
         rc=status)
    call LIS_verify(status)


    !Read lis.config entries
    call ESMF_ConfigFindLabel(LIS_config,'LPRM AMSRE soil moisture data directory:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,lprm_smobsdir,&
            rc=status)
       call LIS_verify(status, 'Err: LPRM AMSRE soil moisture data directory: not defined')

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            lprm_smobsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),'Data Update Status',&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'LPRM AMSRE soil moisture observations attributes file:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsAttribFile(n),rc=status)
       call LIS_verify(status, 'Err: LPRM AMSRE soil moisture observations attributes file: not defined')
    enddo

!!!    call ESMF_ConfigFindLabel(LIS_config,"LPRM AMSRE soil moisture use raw data:",&
!!!         rc=status)
!!!    do n=1,LIS_rc%nnest
!!!       call ESMF_ConfigGetAttribute(LIS_config,LPRM_AMSREsm_obs_struc(n)%rawdata,&
!!!            rc=status)
!!!       call LIS_verify(status, 'LPRM AMSRE soil moisture use raw data: is missing')
!!!    enddo

!!!    call ESMF_ConfigFindLabel(LIS_config,'Overpass hr descending:',&
!!!         rc=status)
!!!    do n=1,LIS_rc%nnest
!!!       call ESMF_ConfigGetAttribute(LIS_config,LPRM_AMSREsm_obs_struc(n)%overpass_hr_d,&
!!!            rc=status)
!!!       call LIS_verify(status,'Overpass hr descending: not defined')
!!!    enddo
!!!
!!!    call ESMF_ConfigFindLabel(LIS_config,'Overpass hr ascending:', &
!!!         rc=status)
!!!    do n=1,LIS_rc%nnest
!!!       call ESMF_ConfigGetAttribute(LIS_config,LPRM_AMSREsm_obs_struc(n)%overpass_hr_a,&
!!!            rc=status)
!!!       call LIS_verify(status,'Overpass hr ascending: not defined')
!!!    enddo
!!!
!!!    call ESMF_ConfigFindLabel(LIS_config,'Mask hr ascending lower:', &
!!!         rc=status)
!!!    do n=1,LIS_rc%nnest
!!!       call ESMF_ConfigGetAttribute(LIS_config,LPRM_AMSREsm_obs_struc(n)%mask_hr_a_lower,&
!!!            rc=status)
!!!       call LIS_verify(status,'Mask hr ascending lower: not defined')
!!!    enddo
!!!
!!!    call ESMF_ConfigFindLabel(LIS_config,'Mask hr ascending upper:', &
!!!         rc=status)
!!!    do n=1,LIS_rc%nnest
!!!       call ESMF_ConfigGetAttribute(LIS_config,LPRM_AMSREsm_obs_struc(n)%mask_hr_a_upper,&
!!!            rc=status)
!!!       call LIS_verify(status,'Mask hr ascending upper: not defined')
!!!    enddo
!!!
!!!    call ESMF_ConfigFindLabel(LIS_config,'Mask hr descending lower:', &
!!!         rc=status)
!!!    do n=1,LIS_rc%nnest
!!!       call ESMF_ConfigGetAttribute(LIS_config,LPRM_AMSREsm_obs_struc(n)%mask_hr_d_lower,&
!!!            rc=status)
!!!       call LIS_verify(status,'Mask hr descending lower: not defined')
!!!    enddo
!!!
!!!    call ESMF_ConfigFindLabel(LIS_config,'Mask hr descending upper:', &
!!!         rc=status)
!!!    do n=1,LIS_rc%nnest
!!!       call ESMF_ConfigGetAttribute(LIS_config,LPRM_AMSREsm_obs_struc(n)%mask_hr_d_upper,&
!!!            rc=status)
!!!       call LIS_verify(status,'Mask hr descending upper: not defined')
!!!    enddo

    !Determine indexing into arrays, 2 daily entries for asc- and desc-ending passes
    do n=1,LIS_rc%nnest
       LPRM_AMSREsm_obs_struc(n)%nc = 1440
       LPRM_AMSREsm_obs_struc(n)%nr = 720
       allocate(LPRM_AMSREsm_obs_struc(n)%smobs(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(LPRM_AMSREsm_obs_struc(n)%smtime(LIS_rc%lnc(n),LIS_rc%lnr(n)))

!!!       call ESMF_TimeSet( &
!!!            LPRM_AMSREsm_obs_struc(n)%start_time, &
!!!            yy=y0, &
!!!            mm=m0,    &
!!!            dd=1,    &   !actual=8
!!!            h=0,     &   !actual=8
!!!            m=0,     &   !actual=44
!!!            calendar=LIS_calendar,rc=status)
!!!
!!!       call ESMF_TimeSet( &
!!!            LPRM_AMSREsm_obs_struc(n)%end_time, &
!!!            yy=2011, &
!!!            mm=7,    &
!!!            dd=31,   &
!!!            h=19,    &    !actual=19
!!!            m=7,     &    !actual=7
!!!            calendar=LIS_calendar,rc=status) 
!!!
!!!       dt=LPRM_AMSREsm_obs_struc(n)%end_time-LPRM_AMSREsm_obs_struc(n)%start_time
!!!
!!!       !Retrieve number of days between start and stop
!!!       call ESMF_TimeIntervalGet(dt,d=rec_len_days,rc=status)
!!!
!!!       !Allocate dynamic arrays
!!!       !allocate(LPRM_AMSREsm_obs_struc(n)%time(numrecords))
!!!       !allocate(LPRM_AMSREsm_obs_struc(n)%lat(numrecords))
!!!       !allocate(LPRM_AMSREsm_obs_struc(n)%lon(numrecords))
!!!
!!!!       allocate(emiss(rec_len_days,numfreqs,numpolarizations,numdailypasses))  !7 freqs, 3 polarizations, 2 passes
!!!
!!!
!!!       !initialize all emissivities prior to filling in
!!!       allocate(tb(numfreqs,numpolarizations,numdailypasses))
!!!       allocate(em(numfreqs,numpolarizations,numdailypasses))
!!!       allocate(mpdi(numfreqs,numdailypasses))
!!!
!!!       allocate(LPRM_AMSREsm_obs_struc(n)%lprm_sm(LIS_rc%gnc(n),LIS_rc%gnr(n),rec_len_days,numfreqs,numpolarizations, numdailypasses))
!!!       allocate(LPRM_AMSREsm_obs_struc(n)%mpdi(LIS_rc%gnc(n),LIS_rc%gnr(n),rec_len_days,numfreqs, numdailypasses))
!!!       allocate(LPRM_AMSREsm_obs_struc(n)%nobs(LIS_rc%gnc(n),LIS_rc%gnr(n),rec_len_days, numdailypasses))             
!!!       allocate(sum_em(LIS_rc%gnc(n),LIS_rc%gnr(n),rec_len_days,numfreqs,numpolarizations, numdailypasses))
!!!       allocate(sum_em_mpdi(LIS_rc%gnc(n),LIS_rc%gnr(n),rec_len_days,numfreqs, numdailypasses))
!!!       allocate(sum_obs(LIS_rc%gnc(n),LIS_rc%gnr(n),rec_len_days,numdailypasses))
    enddo

    do n=1,LIS_rc%nnest
       LPRM_AMSREsm_obs_struc(n)%lprmnc = 1440
       LPRM_AMSREsm_obs_struc(n)%lprmnr = 720

       call map_set(PROJ_LATLON, -89.875,-179.875,&
            0.0, 0.25,0.25, 0.0,&
            LPRM_AMSREsm_obs_struc(n)%lprmnc,LPRM_AMSREsm_obs_struc(n)%lprmnr,&
            LPRM_AMSREsm_obs_struc(n)%lprmproj)
       
       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = 1440
       gridDesci(3) = 720
       gridDesci(4) = -89.875
       gridDesci(5) = -179.875
       gridDesci(6) = 128
       gridDesci(7) = 89.875
       gridDesci(8) = 179.875
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64
       
       
       allocate(LPRM_AMSREsm_obs_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(LPRM_AMSREsm_obs_struc(n)%n12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(LPRM_AMSREsm_obs_struc(n)%n21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(LPRM_AMSREsm_obs_struc(n)%n22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       
       allocate(LPRM_AMSREsm_obs_struc(n)%w11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(LPRM_AMSREsm_obs_struc(n)%w12(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(LPRM_AMSREsm_obs_struc(n)%w21(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(LPRM_AMSREsm_obs_struc(n)%w22(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       
       call bilinear_interp_input(gridDesci,&
            LIS_rc%lnc(n), LIS_rc%lnr(n),&
            LIS_domain(n)%lat, LIS_domain(n)%lon,&
            LPRM_AMSREsm_obs_struc(n)%n11, &
            LPRM_AMSREsm_obs_struc(n)%n12, LPRM_AMSREsm_obs_struc(n)%n21, &
            LPRM_AMSREsm_obs_struc(n)%n22, LPRM_AMSREsm_obs_struc(n)%w11, &
            LPRM_AMSREsm_obs_struc(n)%w12, LPRM_AMSREsm_obs_struc(n)%w21, &
            LPRM_AMSREsm_obs_struc(n)%w22)
       
       call LIS_registerAlarm("LPRM AMSRE soil moisture read alarm",&
            86400.0, 86400.0)
       LPRM_AMSREsm_obs_struc(n)%startMode = .true. 

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
    enddo

!!!    do n=1,LIS_rc%nnest
!!!       !initializations
!!!       LPRM_AMSREsm_obs_struc(n)%lprm_sm = 0  !not LIS_rc%udef
!!!       LPRM_AMSREsm_obs_struc(n)%mpdi = 0        !not LIS_rc%udef
!!!       LPRM_AMSREsm_obs_struc(n)%nobs = 0        !LIS_rc%udef
!!!       
!!!       sum_em=0
!!!       sum_em_mpdi=0
!!!       sum_obs=0
!!!
!!!       Do  im=1, maxm                                                                
!!!          if (m0 .ne. 12) then                                                     
!!!             write(cym, '(I4, I2.2)') y0, m0                                         
!!!             m0=m0+1                                                                 
!!!          else                                                                     
!!!             write(cym, '(I4, I2.2)') y0, m0                                         
!!!             y0=y0+1                                                                 
!!!             m0=1                                                                    
!!!          end if
!!!          datfile(im) = trim(lprm_smobsdir)//"raw/emissivities_"//cym//"_AMSRE.output"                   
!!!          write(*, *)"Reading ", datfile(im)                                        
!!!
!!!          ftn=LIS_getNextUnitNumber()
!!!          open(ftn, file=datfile(im), form="formatted")                                
!!!          
!!!          ios = 0 
!!!          numrecords = 0 
!!!          
!!!          do while(ios.eq.0)
!!!             numrecords = numrecords + 1
!!!             read(ftn,'(i6,4i4,5f8.2,12f8.2,12f9.3)',iostat=ios) yr, mo, da, hr, mn, lat, lon, &
!!!                  LWP, TPW, tsurf, tb_temp(:), em_temp(:)
!!!             
!!!             !Set interval of one day
!!!             call ESMF_TimeIntervalSet(esmf_one_day,d=1,rc=status)
!!!             
!!!             
!!!             !additional complexity with the asc/desc diff's w/in record.  sometimes 2 ascends in one day
!!!             !e.g., 00:03 followed by 23:59.  w/ current approach, assign to assumed overpass times
!!!             !w/ current approach,  assign to the proper (assumed) day and pass time
!!!             
!!!             call ESMF_TimeSet(t_curr, yy=yr, mm=mo, dd=da,h=hr,&
!!!                  m=mn,calendar=LIS_calendar,rc=status)
!!!             
!!!             !Set time for same day but with overpass hr (for both passes)
!!!             call ESMF_TimeSet(t_curr_a, yy=yr, mm=mo, dd=da,&
!!!                  h=LPRM_AMSREsm_obs_struc(n)%overpass_hr_a, &
!!!                  m=0,calendar=LIS_calendar,rc=status)
!!!             
!!!             call ESMF_TimeSet(t_curr_d, yy=yr, mm=mo, dd=da,&
!!!                  h=LPRM_AMSREsm_obs_struc(n)%overpass_hr_d, &
!!!                  m=0,calendar=LIS_calendar,rc=status)
!!!             
!!!             !Find overpass time for previous day and next day
!!!             t_prev_a = t_curr_a - esmf_one_day
!!!             t_prev_d = t_curr_d - esmf_one_day
!!!             
!!!             t_next_a = t_curr_a + esmf_one_day
!!!             t_next_d = t_curr_d + esmf_one_day
!!!             
!!!             !Find which of curr, prev, or next overpass times and desc or ascend is closest to t_curr
!!!             
!!!             dt_curr_a = t_curr - t_curr_a
!!!             dt_curr_a = ESMF_TimeIntervalAbsValue(dt_curr_a)
!!!             
!!!             dt_curr_d = t_curr - t_curr_d
!!!             dt_curr_d = ESMF_TimeIntervalAbsValue(dt_curr_d)
!!!             
!!!             dt_prev_a = t_curr - t_prev_a
!!!             dt_prev_a = ESMF_TimeIntervalAbsValue(dt_prev_a)
!!!             
!!!             dt_prev_d = t_curr - t_prev_d
!!!             dt_prev_d = ESMF_TimeIntervalAbsValue(dt_prev_d)
!!!             
!!!             dt_next_a = t_curr - t_next_a
!!!             dt_next_a = ESMF_TimeIntervalAbsValue(dt_next_a)
!!!             
!!!             dt_next_d = t_curr - t_next_d
!!!             dt_next_d = ESMF_TimeIntervalAbsValue(dt_next_d)
!!!             
!!!             !initialize
!!!             dt_min = dt_curr_a
!!!             t_temp = t_curr_a
!!!             is_ascend_pass = .true.
!!!             
!!!             if (dt_prev_a < dt_min) then
!!!                dt_min = dt_prev_a
!!!                t_temp = t_prev_a
!!!                is_ascend_pass = .true.
!!!             endif
!!!             
!!!             if (dt_next_a < dt_min) then
!!!                dt_min = dt_next_a
!!!                t_temp = t_next_a
!!!                is_ascend_pass = .true.
!!!             endif
!!!             
!!!             if (dt_curr_d < dt_min) then
!!!                dt_min = dt_curr_d
!!!                t_temp = t_curr_d
!!!                is_ascend_pass = .false.
!!!             endif
!!!             
!!!             if (dt_prev_d < dt_min) then
!!!                dt_min = dt_prev_d
!!!                t_temp = t_prev_d
!!!                is_ascend_pass = .false.
!!!             endif
!!!             
!!!             if (dt_next_d < dt_min) then
!!!                dt_min = dt_next_d
!!!                t_temp = t_next_d
!!!                is_ascend_pass = .false.
!!!             endif
!!!             
!!!             !Get index for lprm_sm array
!!!             delta_t = t_temp - LPRM_AMSREsm_obs_struc(n)%start_time
!!!             call ESMF_TimeIntervalGet(delta_t, d=day_index, &
!!!                  calendar=LIS_calendar,rc=status)
!!!             call LIS_verify(status, 'error in timeget: LPRM_AMSREsm_obsMod.F90')
!!!             
!!!             !as fortran arrays are zero-based...
!!!             day_index=day_index+1
!!!             
!!!             if (is_ascend_pass) then
!!!                ipass=LPRM_AMSREsm_obs_struc(n)%index_ascend
!!!             else
!!!                ipass=LPRM_AMSREsm_obs_struc(n)%index_descend
!!!             endif
!!!             
!!!             !convert to more convenient data structure
!!!             do f=1,numfreqs
!!!                do p=1,numpolarizations
!!!                   c=(f-1)*2+p  !requires channels for both polarizations (V/H) for all freqs
!!!                   tb(f,p,ipass)=em_temp(c)
!!!                   em(f,p,ipass)=em_temp(c)
!!!                enddo
!!!             enddo
!!!             
!!!             !compute mpdi
!!!             do f=1,numfreqs
!!!                mpdi(f,ipass)=(em(f,1,ipass)-em(f,2,ipass))/(em(f,1,ipass)+em(f,2,ipass))
!!!             enddo
!!!
!!!             !compute i,j
!!!             call latlon_to_ij(LIS_domain(n)%lisproj, lat, lon, col_real, row_real)
!!!             col=nint(col_real)
!!!             row=nint(row_real)
!!!
!!!             if ((col.lt.1).or.(col.gt.LIS_rc%lnc(n)).or.(row.lt.1).or.(row.gt.LIS_rc%lnr(n))) then
!!!                !nothing--out of LIS run domain
!!!             else  ! in LIS run domain
!!!                sum_obs(col,row,day_index,ipass)= &
!!!                     sum_obs(col,row,day_index,ipass) + 1
!!!                
!!!                do f=1,numfreqs
!!!                   
!!!                   sum_em_mpdi(col,row,day_index,f,ipass)= &
!!!                        sum_em_mpdi(col,row,day_index,f,ipass)+ mpdi(f,ipass)
!!!                   
!!!                   do p=1,numpolarizations
!!!                      
!!!                      sum_em(col,row,day_index,f,p,ipass)= &
!!!                           sum_em(col,row,day_index,f,p,ipass)+ em(f,p,ipass)
!!!                      
!!!                   enddo
!!!                enddo
!!!             end if
!!!          end do
!!!          call LIS_releaseUnitNumber(ftn)
!!!       enddo
!!!       
!!!       LPRM_AMSREsm_obs_struc(n)%nobs= sum_obs
!!!       
!!!       do day_index=1,rec_len_days
!!!          do col=1, LIS_rc%lnc(n)
!!!             do row=1, LIS_rc%lnr(n)
!!!                do ipass=1,2
!!!                   do f=1,numfreqs
!!!                      if (sum_obs(col,row,day_index,ipass).lt.LPRM_AMSREsm_obs_struc(n)%nobs_threshold) then
!!!                         LPRM_AMSREsm_obs_struc(n)%mpdi(col,row,day_index,f,ipass)=LIS_rc%udef
!!!                      else
!!!                         LPRM_AMSREsm_obs_struc(n)%mpdi(col,row,day_index,f,ipass)= &
!!!                              sum_em_mpdi(col,row,day_index,f,ipass) &
!!!                              / sum_obs(col,row,day_index,ipass)
!!!                      end if
!!!                      do p=1,numpolarizations
!!!                         if (sum_obs(col,row,day_index,ipass).lt.LPRM_AMSREsm_obs_struc(n)%nobs_threshold) then
!!!                            LPRM_AMSREsm_obs_struc(n)%lprm_sm(col,row,day_index,f,p,ipass)= LIS_rc%udef
!!!                         else
!!!                            LPRM_AMSREsm_obs_struc(n)%lprm_sm(col,row,day_index,f,p,ipass)= &
!!!                                 sum_em(col,row,day_index,f,p,ipass) &
!!!                                 / sum_obs(col,row,day_index,ipass)
!!!                         end if
!!!                      enddo
!!!                   enddo
!!!                enddo
!!!             enddo
!!!          enddo
!!!       enddo
!!!    enddo

    write(LIS_logunit,*) 'Created the States to hold the AMSRE_SR lprm_sm observations'
 
  end subroutine LPRM_AMSREsm_obs_setup

end module LPRM_AMSREsm_obsMod
