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
! !MODULE: AMSRE_SR_em_obsMod
! 
! !DESCRIPTION: 
!   
! !REVISION HISTORY: 
!  11 Jul 11    Ken Harrison;   Initial Specification
! 
module AMSRE_SR_em_obsMod
! !USES: 
  use ESMF
  use map_utils
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: AMSRE_SR_em_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SRemobs

  type, public :: emissivity_data_dec
     character*2, allocatable      :: fnames(:)
     character*1, allocatable      :: pnames(:)
     integer                   :: numfreqs
     integer                   :: numchannels
     integer                   :: numpolarizations
     integer                   :: numdailypasses
     integer                   :: numrecords
     integer                   :: nobs_threshold    
     real, allocatable             :: emissivity(:,:,:,:,:,:) !col, row, numrecords, 6 freqs, 2 polarizations, 2 passes
     real, allocatable             :: mpdi(:,:,:,:,:)         !col, row, numrecords, 6 freqs, 2 passes
     integer, allocatable          :: nobs(:,:,:,:)           !col, row, numrecords, 2 passes
     type(ESMF_Time)           :: start_time, end_time
     integer                   :: overpass_hr_a
     integer                   :: overpass_hr_d
     integer                   :: mask_hr_a_lower
     integer                   :: mask_hr_d_lower
     integer                   :: mask_hr_a_upper
     integer                   :: mask_hr_d_upper
     integer                   :: index_ascend
     integer                   :: index_descend
     integer                   :: rec_len_days
!!!     type(proj_info)     :: projection
  end type emissivity_data_dec

  type(emissivity_data_dec), allocatable :: SRemobs(:)

contains
!BOP
! 
! !ROUTINE: AMSRE_SR_em_obs_setup
! \label{AMSRE_SR_em_obs_setup}
! 
! !INTERFACE: 
  subroutine AMSRE_SR_em_obs_setup(Obs_State)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_config, &
         LIS_vecGrid, LIS_domain, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_calendar
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
    integer                   ::  n
    integer                   ::  status
    type(ESMF_ArraySpec)      ::  realarrspec
    type(ESMF_Field)          ::  obsField
    character(len=LIS_CONST_PATH_LEN) ::  emissivityobsdir
    character(len=LIS_CONST_PATH_LEN) ::  emissivityobsmaskdir
    character*100             ::  vname
    character(len=LIS_CONST_PATH_LEN) ::  obsAttribFile(LIS_rc%nnest)
    integer                   ::  ftn
    integer                   ::  numrecords, ios
    integer                   ::  yr
    real                      ::  lat, lon
    real                      ::  LWP, TPW, tsurf
    real, allocatable             ::  tb_temp(:)
    real, allocatable             ::  em_temp(:)
    real, allocatable             ::  tb(:,:,:)
    real, allocatable             ::  em(:,:,:)
    real, allocatable             ::  mpdi(:,:)
    integer                   ::  mo, da, hr, mn
    integer                   ::  yr_curr, mo_curr, da_curr, hr_curr, mn_curr
    type(ESMF_TimeInterval)   ::  dt,esmf_jday,  esmf_one_day, delta_t
    type(ESMF_Time)           ::  temp_time, t_temp, &
         t_curr, t_curr_a, t_curr_d, t_prev_a, t_prev_d, t_next_a, t_next_d
    type(ESMF_Time)           ::  masktime
    type(ESMF_TimeInterval)   ::  dt_curr_a, dt_curr_d, &
         dt_prev_a, dt_prev_d, dt_next_a, dt_next_d,dt_a, dt_d, dt_min
    !    real                      ::  gridDesci(LIS_rc%nnest,50)
    integer, parameter        ::  i_undef=-999
    !    integer                   ::  nobs_threshold
    logical                   ::  is_ascend_pass
    integer                   ::  day_index
    integer                   ::  k,i
    integer                   ::  ipass
    integer                   ::  col, row !lat lon indices
    real                      ::  col_real, row_real
    integer                   ::  c !channel counter
    integer                   ::  f !freq counter
    integer                   ::  p !polarization counter
    !    real                      ::  amsre_sr_em_gridDesc(6)
    real, allocatable       :: sum_em(:,:,:,:,:,:) !row, col, numrecords, 6 freqs, 2 polarizations, 2 passes
    real, allocatable       :: sum_em_mpdi(:,:,:,:,:) !row, col, numrecords, 6 freqs, 2 passes
    integer, allocatable    :: sum_obs(:,:,:,:) ! row, col, numrecords, 2 passes

    integer, parameter :: maxd = 4000  ! max number of days for daily overpass 
    integer, parameter :: maxm = 108     ! max number of months    
    character(len=LIS_CONST_PATH_LEN) :: datfile(maxm)

    ! time management                                                          
    character*6 cym                                                  
    integer :: im, y0, m0                           

!!!  TEMP DEBUGGING --START
    type(ESMF_Time)         :: lis_time1
!!!  TEMP DEBUGGING -- FINISH

!!!    amsre_sr_em_gridDesc(1) = 34.125   !lower left lat
!!!    amsre_sr_em_gridDesc(2) = -99.875  !lower left lon
!!!    amsre_sr_em_gridDesc(3) = 38.875   !upper right lat
!!!    amsre_sr_em_gridDesc(4) = -95.125   !upper right lon
!!!    amsre_sr_em_gridDesc(5) = 0.25     !res
!!!    amsre_sr_em_gridDesc(6) = 0.25     !res


!!!!!!!!!!!   SEE HACK BELOW TO FIX OUTPUT OF LATLON_TO_IJ


    allocate(SRemobs(LIS_rc%nnest))
    do n=1,LIS_rc%nnest
       !       call map_set(PROJ_LATLON, 34.125, -99.875, &
       !            0.0, 0.25, 0.25, 0.0, 20, 20, SRemobs(n)%projection)
       SRemobs(n)%numfreqs=6        
       SRemobs(n)%numchannels=12    
       SRemobs(n)%numpolarizations=2
       SRemobs(n)%numdailypasses=2
       SRemobs(n)%index_ascend  = 1
       SRemobs(n)%index_descend = 2
       allocate(SRemobs(n)%fnames(SRemobs(n)%numfreqs))
       allocate(SRemobs(n)%pnames(SRemobs(n)%numpolarizations))
       SRemobs(n)%fnames=(/'07', '11', '19', '24', '37', '89'/)
       SRemobs(n)%pnames=(/'V','H'/)
       allocate(tb_temp(SRemobs(n)%numchannels))
       allocate(em_temp(SRemobs(n)%numchannels))
    enddo

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,& !Note: rank=1
         rc=status)
    call LIS_verify(status)


    !Read lis.config entries
    call ESMF_ConfigFindLabel(LIS_config,'AMSRE_SR Emissivity Obs data directory:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,emissivityobsdir,&
            rc=status)
       call LIS_verify(status, 'Err: AMSRE_SR Emissivity Obs data directory: not defined')

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            emissivityobsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),'Data Update Status',&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'AMSRE_SR Emissivity observations attributes file:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,obsAttribFile(n),rc=status)
       call LIS_verify(status, 'Err: AMSRE_SR Emissivity observations attributes file: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'AMSRE_SR number of observations threshold:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SRemobs(n)%nobs_threshold,&
            rc=status)
       call LIS_verify(status,'AMSRE_SR number of observations threshold: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Overpass hr descending:',&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SRemobs(n)%overpass_hr_d,&
            rc=status)
       call LIS_verify(status,'Overpass hr descending: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Overpass hr ascending:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SRemobs(n)%overpass_hr_a,&
            rc=status)
       call LIS_verify(status,'Overpass hr ascending: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Mask hr ascending lower:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SRemobs(n)%mask_hr_a_lower,&
            rc=status)
       call LIS_verify(status,'Mask hr ascending lower: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Mask hr ascending upper:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SRemobs(n)%mask_hr_a_upper,&
            rc=status)
       call LIS_verify(status,'Mask hr ascending upper: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Mask hr descending lower:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SRemobs(n)%mask_hr_d_lower,&
            rc=status)
       call LIS_verify(status,'Mask hr descending lower: not defined')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,'Mask hr descending upper:', &
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,SRemobs(n)%mask_hr_d_upper,&
            rc=status)
       call LIS_verify(status,'Mask hr descending upper: not defined')
    enddo

    ! data beginning 
    y0=2002                                                                    
    m0=8                                                                       

    !Determine indexing into arrays, 2 daily entries for asc- and desc-ending passes
    do n=1,LIS_rc%nnest
       call ESMF_TimeSet( &
            SRemobs(n)%start_time, &
            yy=y0, &
            mm=m0,    &
            dd=1,    &   !actual=8
            h=0,     &   !actual=8
            m=0,     &   !actual=44
            calendar=LIS_calendar,rc=status)

       call ESMF_TimeSet( &
            SRemobs(n)%end_time, &
            yy=2011, &
            mm=7,    &
            dd=31,   &
            h=19,    &    !actual=19
            m=7,     &    !actual=7
            calendar=LIS_calendar,rc=status) 

       dt=SRemobs(n)%end_time-SRemobs(n)%start_time

       !Retrieve number of days between start and stop
       call ESMF_TimeIntervalGet(dt,d=SRemobs(n)%rec_len_days,rc=status)

       !Allocate dynamic arrays
       !allocate(SRemobs(n)%time(numrecords))
       !allocate(SRemobs(n)%lat(numrecords))
       !allocate(SRemobs(n)%lon(numrecords))

       !       allocate(emiss(SRemobs(n)%rec_len_days,numfreqs,numpolarizations,numdailypasses))  !7 freqs, 3 polarizations, 2 passes


       !initialize all emissivities prior to filling in
       allocate(tb(SRemobs(n)%numfreqs,SRemobs(n)%numpolarizations,&
            SRemobs(n)%numdailypasses))
       allocate(em(SRemobs(n)%numfreqs,SRemobs(n)%numpolarizations,&
            SRemobs(n)%numdailypasses))
       allocate(mpdi(SRemobs(n)%numfreqs,SRemobs(n)%numdailypasses))

       allocate(SRemobs(n)%emissivity(LIS_rc%gnc(n),LIS_rc%gnr(n),&
            SRemobs(n)%rec_len_days+1,SRemobs(n)%numfreqs,&
            SRemobs(n)%numpolarizations, SRemobs(n)%numdailypasses))
       allocate(SRemobs(n)%mpdi(LIS_rc%gnc(n),LIS_rc%gnr(n),&
            SRemobs(n)%rec_len_days+1,SRemobs(n)%numfreqs, &
            SRemobs(n)%numdailypasses))
       allocate(SRemobs(n)%nobs(LIS_rc%gnc(n),LIS_rc%gnr(n),&
            SRemobs(n)%rec_len_days+1, SRemobs(n)%numdailypasses))             
       allocate(sum_em(LIS_rc%gnc(n),LIS_rc%gnr(n),&
            SRemobs(n)%rec_len_days+1,SRemobs(n)%numfreqs,&
            SRemobs(n)%numpolarizations, SRemobs(n)%numdailypasses))
       allocate(sum_em_mpdi(LIS_rc%gnc(n),LIS_rc%gnr(n),&
            SRemobs(n)%rec_len_days+1,SRemobs(n)%numfreqs, &
            SRemobs(n)%numdailypasses))
       allocate(sum_obs(LIS_rc%gnc(n),LIS_rc%gnr(n),&
            SRemobs(n)%rec_len_days+1,SRemobs(n)%numdailypasses))
    enddo

    do n=1,LIS_rc%nnest
       ftn=LIS_getNextUnitNumber()
       open(ftn,file=trim(obsAttribFile(n)),status='old')
       read(ftn,*)
       read(ftn,fmt='(a100)') vname

       do f=1,SRemobs(n)%numfreqs
          do p=1,SRemobs(n)%numpolarizations
             obsField = ESMF_FieldCreate(arrayspec=realarrspec, &
                  grid=LIS_vecGrid(n), &
                  name=trim(vname) // SRemobs(n)%fnames(f) // SRemobs(n)%pnames(p), rc=status)
             call LIS_verify(status)
             call ESMF_StateAdd(Obs_State(n),(/obsField/),rc=status)
             call LIS_verify(status)
          enddo
       enddo

       call LIS_releaseUnitNumber(ftn)
    enddo

    do n=1,LIS_rc%nnest
       !initializations
       SRemobs(n)%emissivity = 0  !not LIS_rc%udef
       SRemobs(n)%mpdi = 0        !not LIS_rc%udef
       SRemobs(n)%nobs = 0        !LIS_rc%udef

       sum_em=0
       sum_em_mpdi=0
       sum_obs=0

       Do  im=1, maxm                                                                
          if (m0 .ne. 12) then                                                     
             write(cym, '(I4, I2.2)') y0, m0                                         
             m0=m0+1                                                                 
          else                                                                     
             write(cym, '(I4, I2.2)') y0, m0                                         
             y0=y0+1                                                                 
             m0=1                                                                    
          end if
          datfile(im) = trim(emissivityobsdir)//"raw/emissivities_"//cym//"_AMSRE.output"                   
          write(*, *)"Reading ", datfile(im)                                        

          ftn=LIS_getNextUnitNumber()
          open(ftn, file=datfile(im), form="formatted")                                

          ios = 0 
          numrecords = 0 

          do while(ios.eq.0)
             numrecords = numrecords + 1
             read(ftn,'(i6,4i4,5f8.2,12f8.2,12f9.3)',iostat=ios) &
                  yr, mo, da, hr, mn, lat, lon, &
                  LWP, TPW, tsurf, tb_temp(:), em_temp(:)
!if (numrecords.le.50.and.LIS_masterproc) then
!print *, yr, mo, da, hr, mn, lat, lon, & 
!                  LWP, TPW, tsurf, tb_temp(:), em_temp(:)
!endif
             !Set interval of one day
             call ESMF_TimeIntervalSet(esmf_one_day,d=1,rc=status)


             !additional complexity with the asc/desc diff's w/in record.  sometimes 2 ascends in one day
             !e.g., 00:03 followed by 23:59.  w/ current approach, assign to assumed overpass times
             !w/ current approach,  assign to the proper (assumed) day and pass time

             call ESMF_TimeSet(t_curr, yy=yr, mm=mo, dd=da,h=hr,&
                  m=mn,calendar=LIS_calendar,rc=status)

             !Set time for same day but with overpass hr (for both passes)
             call ESMF_TimeSet(t_curr_a, yy=yr, mm=mo, dd=da,&
                  h=SRemobs(n)%overpass_hr_a, &
                  m=0,calendar=LIS_calendar,rc=status)

             call ESMF_TimeSet(t_curr_d, yy=yr, mm=mo, dd=da,&
                  h=SRemobs(n)%overpass_hr_d, &
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
             delta_t = t_temp - SRemobs(n)%start_time
             call ESMF_TimeIntervalGet(delta_t, d=day_index, &
                  calendar=LIS_calendar,rc=status)
             call LIS_verify(status, 'error in timeget: AMSRE_SR_em_obsMod.F90')

             !as fortran arrays are zero-based...
             day_index=day_index+1

             if (is_ascend_pass) then
                ipass=SRemobs(n)%index_ascend
             else
                ipass=SRemobs(n)%index_descend
             endif

             !convert to more convenient data structure
             do f=1,SRemobs(n)%numfreqs
                do p=1,SRemobs(n)%numpolarizations
                   c=(f-1)*2+p  !requires channels for both polarizations (V/H) for all freqs
                   tb(f,p,ipass)=tb_temp(c)
                   em(f,p,ipass)=em_temp(c)
                enddo
             enddo

             !compute mpdi
             do f=1,SRemobs(n)%numfreqs
                mpdi(f,ipass)=(em(f,1,ipass)-em(f,2,ipass))/(em(f,1,ipass)+em(f,2,ipass))
             enddo

             !compute i,j
             call latlon_to_ij(LIS_domain(n)%lisproj, lat, lon, col_real, row_real)
             col=nint(col_real)
             row=nint(row_real)


!!!  START HACK -- latlon_to_ij NOT SUITED FOR HANDLING CASE WHERE LON IS IN AND JUST LEFT OF LOWER LEFT CELL
             if (col.eq.1441) then
                col=1
             endif
!!!  FINISH HACK

             if ((col.lt.1).or.(col.gt.LIS_rc%lnc(n)).or.(row.lt.1).or.(row.gt.LIS_rc%lnr(n))) then
                !nothing--out of LIS run domain
!                print *, "Out of domain"
!                if (LIS_masterproc) then
!                   print *, LIS_domain(n)%lisproj,col,row,col_real,row_real,lon,lat
!                endif
             else  ! in LIS run domain
                sum_obs(col,row,day_index,ipass)= &
                     sum_obs(col,row,day_index,ipass) + 1

                do f=1,SRemobs(n)%numfreqs

                   sum_em_mpdi(col,row,day_index,f,ipass)= &
                        sum_em_mpdi(col,row,day_index,f,ipass)+ mpdi(f,ipass)

                   do p=1,SRemobs(n)%numpolarizations

                      sum_em(col,row,day_index,f,p,ipass)= &
                           sum_em(col,row,day_index,f,p,ipass)+ em(f,p,ipass)

                   enddo
                enddo
             end if
          end do
          call LIS_releaseUnitNumber(ftn)
       enddo

       SRemobs(n)%nobs= sum_obs

       do day_index=1,SRemobs(n)%rec_len_days+1
          do col=1, LIS_rc%lnc(n)
             do row=1, LIS_rc%lnr(n)
                do ipass=1,SRemobs(n)%numdailypasses
                   do f=1,SRemobs(n)%numfreqs
                      if (sum_obs(col,row,day_index,ipass).lt.&
                           SRemobs(n)%nobs_threshold) then
                         SRemobs(n)%mpdi(col,row,day_index,f,ipass)=LIS_rc%udef
                      else
                         SRemobs(n)%mpdi(col,row,day_index,f,ipass)= &
                              sum_em_mpdi(col,row,day_index,f,ipass) &
                              / sum_obs(col,row,day_index,ipass)
                      end if
                      do p=1,SRemobs(n)%numpolarizations
                         if (sum_obs(col,row,day_index,ipass).lt.&
                              SRemobs(n)%nobs_threshold) then
                            SRemobs(n)%emissivity(col,row,day_index,f,p,ipass)= LIS_rc%udef
                         else
                            SRemobs(n)%emissivity(col,row,day_index,f,p,ipass)= &
                                 sum_em(col,row,day_index,f,p,ipass) &
                                 / sum_obs(col,row,day_index,ipass)
                         end if
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo



!!!!    !TEMP TESTING
!!!!    if (LIS_masterproc) then
!!!!       n=1
!!!!       yr=2008
!!!!       mo=4
!!!!       da=1
!!!!       hr=0
!!!!       mn=0
!!!!       call ESMF_TimeSet(lis_time1, yy=yr, mm=mo, &
!!!!            dd=da,h=hr,&
!!!!            m=mn,calendar=LIS_calendar,rc=status)
!!!!       call LIS_verify(status, 'amsre sr obsmod --tempbugfix1')
!!!!       
!!!!       delta_t = lis_time1 - SRemobs(n)%start_time
!!!!       call ESMF_TimeIntervalGet(delta_t, d=day_index, &
!!!!            calendar=LIS_calendar,rc=status)
!!!!       call LIS_verify(status, 'amsre sr obsmod --tempbugfix2')
!!!!       
!!!!       !as fortran arrays are zero-based...
!!!!       day_index=day_index+1
!!!!       
!!!!       Open(213, file="AMSRE.txt")
!!!!       do i=day_index,SRemobs(n)%rec_len_days
!!!!          write(213, '(24F14.4)') i, SRemobs(n)%emissivity(1,1,i,2,2,1),SRemobs(n)%emissivity(1,1,i,2,2,2)
!!!!       end do
!!!!       close(213)
!!!!       
!!!!       write(LIS_logunit,*) 'Created the States to hold the AMSRE_SR emissivity observations'
!!!!    endif
  end subroutine AMSRE_SR_em_obs_setup
  
end module AMSRE_SR_em_obsMod
