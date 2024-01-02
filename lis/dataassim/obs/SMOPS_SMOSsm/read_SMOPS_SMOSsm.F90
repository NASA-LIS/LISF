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
!BOP
! !ROUTINE: read_SMOPS_SMOSsm
! \label{read_SMOPS_SMOSsm}
!
! !REVISION HISTORY:
!  17 Jun 2010: Sujay Kumar; Updated for use with LPRM AMSRE Version 5. 
!  20 Sep 2012: Sujay Kumar; Updated to the NETCDF version of the data. 
!  28 Sep 2017: Mahdi Navari; Updated to read SMOS from SMOPS V3
!
! !INTERFACE: 
subroutine read_SMOPS_SMOSsm(n, k, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use map_utils
  use LIS_pluginIndices
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use SMOPS_SMOSsm_Mod, only : SMOPS_SMOSsm_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the AMSRE soil moisture observations 
!  from NETCDF files and applies the spatial masking for dense
!  vegetation, rain and RFI. The data is then rescaled
!  to the land surface model's climatology using rescaling
!  algorithms. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  real, parameter        ::  minssdev = 0.05
  real, parameter        ::  maxssdev = 0.11
  real,  parameter       :: MAX_SM_VALUE=0.45, MIN_SM_VALUE=0.0001
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: smfield, pertField

  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  real                   :: obs_unsc(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: smobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: smtime(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: sm_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: dt
  integer                :: fnd
  real, allocatable      :: ssdev(:)
  integer                :: lis_julss
  real                   :: smvalue
  real                   :: model_delta(LIS_rc%obs_ngrid(k))
  real                   :: obs_delta(LIS_rc%obs_ngrid(k))
  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 
  obs_unsc = LIS_rc%udef
!-------------------------------------------------------------------------
!   Read both ascending and descending passes at 0Z and then store
!   the overpass time as 1.30AM for the descending pass and 1.30PM 
!   for the ascending pass. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SMOPS read alarm")
  
  if(alarmCheck.or.SMOPS_SMOSsm_struc(n)%startMode) then 
     SMOPS_SMOSsm_struc(n)%startMode = .false.

     SMOPS_SMOSsm_struc(n)%smobs = LIS_rc%udef
     SMOPS_SMOSsm_struc(n)%smtime = -1

     smobs = LIS_rc%udef
     smtime = LIS_rc%udef

     
     call create_SMOPS_SMOSsm_filename(smobsdir, &
          SMOPS_SMOSsm_struc(n)%useRealtime, &
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, LIS_rc%hr, SMOPS_SMOSsm_struc(n)%conv, fname)

     inquire(file=fname,exist=file_exists)

     if(file_exists) then 
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
        call read_RTSMOPS_SMOS_data(n,k,fname,smobs,smtime)
     else
        write(LIS_logunit,*) '[WARN] Missing SMOPS ',trim(fname)
     endif

     SMOPS_SMOSsm_struc(n)%smobs  = LIS_rc%udef
     SMOPS_SMOSsm_struc(n)%smtime = -1

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           grid_index = LIS_obs_domain(n,k)%gindex(c,r)
           if(grid_index.ne.-1) then 
              if(smobs(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then 
                 SMOPS_SMOSsm_struc(n)%smobs(c,r) = &
                      smobs(c+(r-1)*LIS_rc%obs_lnc(k))                 
                 SMOPS_SMOSsm_struc(n)%smtime(c,r) = &
                      smtime(c+(r-1)*LIS_rc%obs_lnc(k))
              endif
           endif
        enddo
     enddo

  endif
  
  call ESMF_StateGet(OBS_State,"Observation01",smfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')

  fnd = 0 
  sm_current = LIS_rc%udef
 
! dt is not defined as absolute value of the time difference to avoid
! double counting of the data in assimilation. 

  call LIS_get_timeoffset_sec(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
       LIS_rc%hr, LIS_rc%mn, LIS_rc%ss, lis_julss)

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           if(SMOPS_SMOSsm_struc(n)%smtime(c,r).ge.0) then 
              dt = (lis_julss-SMOPS_SMOSsm_struc(n)%smtime(c,r))
              if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
                 sm_current(c,r) = & 
                      SMOPS_SMOSsm_struc(n)%smobs(c,r)
                 fnd = 1
              endif           
           endif
        endif
     enddo
  enddo

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     

  if(fnd.ne.0) then        
     ! Store the unscaled obs (ie, before the rescaling)
     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = &
                   sm_current(c,r)
           end if
        end do
     end do
     
     call LIS_rescale_with_CDF_matching(    &
          n,k,                              & 
          SMOPS_SMOSsm_struc(n)%nbins,         & 
          SMOPS_SMOSsm_struc(n)%ntimes,         & 
          MAX_SM_VALUE,                        & 
          MIN_SM_VALUE,                        & 
          SMOPS_SMOSsm_struc(n)%model_xrange,  &
          SMOPS_SMOSsm_struc(n)%obs_xrange,    &
          SMOPS_SMOSsm_struc(n)%model_cdf,     &
          SMOPS_SMOSsm_struc(n)%obs_cdf,       &
          sm_current)
  endif

  obsl = LIS_rc%udef 
  if(fnd.ne.0) then 
     do r=1, LIS_rc%obs_lnr(k)
        do c=1, LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=sm_current(c,r)
           endif
        enddo
     enddo
  endif

  !-------------------------------------------------------------------------
  !  Apply LSM based QC and screening of observations
  !-------------------------------------------------------------------------  
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_SMOPS_SMOSsmobsId)//char(0),n,k,OBS_state)

  call LIS_checkForValidObs(n,k,obsl,fnd,sm_current)

  if(fnd.eq.0) then 
     data_upd_flag_local = .false. 
  else
     data_upd_flag_local = .true. 
  endif
        
#if (defined SPMD)
  call MPI_ALLGATHER(data_upd_flag_local,1, &
       MPI_LOGICAL, data_upd_flag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
  data_upd = .false.
  do p=1,LIS_npes
     data_upd = data_upd.or.data_upd_flag(p)
  enddo

!-------------------------------------------------------------------------
!  Depending on data update flag...
!-------------------------------------------------------------------------     
  
  if(data_upd) then 

     do t=1,LIS_rc%obs_ngrid(k)
        gid(t) = t
        if(obsl(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo
  
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true. , rc=status)
     call LIS_verify(status)

     if(LIS_rc%obs_ngrid(k).gt.0) then 
        call ESMF_AttributeSet(smField,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)
        
        call ESMF_AttributeSet(smField,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(smfield, "Unscaled Obs",&
             obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error in setting Unscaled Obs attribute')
     endif

     if(SMOPS_SMOSsm_struc(n)%useSsdevScal.eq.1) then
        call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')
        
        allocate(ssdev(LIS_rc%obs_ngrid(k)))
        ssdev = SMOPS_SMOSsm_struc(n)%ssdev_inp 

        if(SMOPS_SMOSsm_struc(n)%ntimes.eq.1) then 
           jj = 1
        else
           jj = LIS_rc%mo
        endif
        do t=1,LIS_rc%obs_ngrid(k)
           if(SMOPS_SMOSsm_struc(n)%obs_sigma(t,jj).gt.0) then 
              ssdev(t) = ssdev(t)*SMOPS_SMOSsm_struc(n)%model_sigma(t,jj)/&
                   SMOPS_SMOSsm_struc(n)%obs_sigma(t,jj)
              if(ssdev(t).lt.minssdev) then 
                 ssdev(t) = minssdev
              endif
           endif
        enddo
        
        if(LIS_rc%obs_ngrid(k).gt.0) then 
           call ESMF_AttributeSet(pertField,"Standard Deviation",&
                ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
        endif
        deallocate(ssdev)
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif

end subroutine read_SMOPS_SMOSsm

!BOP
! 
! !ROUTINE: read_RTSMOPS_data
! \label{read_RTSMOPS_data}
!
! !INTERFACE:
subroutine read_RTSMOPS_SMOS_data(n, k, fname, smobs_ip, smtime_ip)
! 
! !USES:   
#if(defined USE_GRIBAPI)
  use grib_api
#endif
  use LIS_coreMod,      only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use SMOPS_SMOSsm_Mod, only : SMOPS_SMOSsm_struc

  implicit none

  integer           :: n
  integer           :: k
  character (len=*) :: fname
  real              :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real              :: smtime_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

! !DESCRIPTION: 
!  This subroutine reads the SMOPS grib2 file and applies the data
!  quality flags to filter the data.
!
!  Quality flags are defined in:
!      NOAA NESDIS
!      CENTER FOR SATELLITE APPLICATIONS AND RESEARCH
!
!      SOIL MOISTURE OPERATIONAL PRODUCT SYSTEM (SMOPS)
!
!      ALGORITHM THEORETICAL BASIS DOCUMENT
!      Version 4.0
!
!  Found at http://www.ospo.noaa.gov/Products/land/smops/documents.html
!
!  The SMOPS QA flags are 16-bit (2-byte) integers, with the least
!  significant byte referred to as byte1 and the most significant byte
!  referred to as byte2.
!
!  bits: 15 | 14 | 13 | 12 | 11 | 10 | 9 | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 | 0
!                   byte2                    |           byte1
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTSMOPS AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP
  integer             :: param_SMOS, param_SMOS_qa
  integer             :: param_SMOS_hr, param_SMOS_mn
  integer*2,parameter :: SMOS_accept1 = b'0000000000000000'
  integer*2,parameter :: SMOS_accept2 = b'0000000000000001'
  integer*2,parameter :: SMOS_accept3 = b'0000000000001000'
  integer*2,parameter :: SMOS_accept4 = b'0000000000001001'
  integer*2,parameter :: SMOS_accept5 = b'0000000010000000'
  integer*2,parameter :: SMOS_accept6 = b'0000000010000001'
  integer*2,parameter :: SMOS_accept7 = b'0000000010001000'
  integer*2,parameter :: SMOS_accept8 = b'0000000010001001'

  real      :: sm_SMOS(SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr)
  real      :: sm_SMOS_t(SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr)
  real      :: sm_SMOS_hr(SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr)
  real      :: sm_SMOS_mn(SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr)

  real      :: sm_SMOS_qa(SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr)
  integer*2 :: sm_SMOS_qa_t(SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr)
  real      :: sm_data(SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr)
  real      :: sm_time(SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr)
  logical*1 :: sm_data_b(SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr)
  logical*1 :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer   :: hr_val, mn_val, julss
  integer   :: c,r,ios
  integer   :: ftn,iret,igrib,nvars
  integer   :: param_num
  logical   :: var_found
  logical   :: smDataNotAvailable
  integer*2 :: qavalue
  integer   :: updoy,yr1,mo1,da1,hr1,mn1,ss1
  real      :: upgmt
  real*8    :: timenow

  smDataNotAvailable = .false. 

#if(defined USE_GRIBAPI)
  yr1 = LIS_rc%yr
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  call LIS_date2time(timenow,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)
!<debug -- smops testing>
! Force version 3
write(LIS_logunit,*) '[Warning] Forcing LIS to read SMOPS version 3'
timenow = SMOPS_SMOSsm_struc(n)%version3_time
!</debug -- smops testing>
  if ( timenow < SMOPS_SMOSsm_struc(n)%version2_time ) then
     if ( SMOPS_SMOSsm_struc(n)%useRealtime == 1 ) then
        write(LIS_logunit,*) '[Warning] NRTSMOS is not availabe in '// &
                             'SMOPS version: 1.3'
        smDataNotAvailable = .true.   
        smobs_ip = LIS_rc%udef
        smtime_ip = LIS_rc%udef
     else     
        param_SMOS = 212; param_SMOS_qa = 233
        param_SMOS_hr = 221; param_SMOS_mn = 222
     endif
  elseif ( timenow >= SMOPS_SMOSsm_struc(n)%version2_time .and. &
           timenow <  SMOPS_SMOSsm_struc(n)%version3_time ) then
     if ( SMOPS_SMOSsm_struc(n)%useRealtime == 1 ) then
        param_SMOS = 211; param_SMOS_qa = 232
        param_SMOS_hr = 219; param_SMOS_mn = 220
     else
        param_SMOS = 212; param_SMOS_qa = 233
        param_SMOS_hr = 221; param_SMOS_mn = 222
     endif
  elseif ( timenow >= SMOPS_SMOSsm_struc(n)%version3_time ) then
     if ( SMOPS_SMOSsm_struc(n)%useRealtime == 1 ) then
        param_SMOS = 211; param_SMOS_qa = 241
        param_SMOS_hr = 222; param_SMOS_mn = 223
     else
        param_SMOS = 212; param_SMOS_qa = 242
        param_SMOS_hr = 224; param_SMOS_mn = 225
     endif
  else
     write(LIS_logunit,*) '[ERR] Invalid times for SMOPS versions'
     write(LIS_logunit,*) '      ', timenow
     write(LIS_logunit,*) '      ', SMOPS_SMOSsm_struc(n)%version2_time
     write(LIS_logunit,*) '      ', SMOPS_SMOSsm_struc(n)%version3_time
     call LIS_endrun()
  endif

  if ( smDataNotAvailable .eqv. .false. ) then
  call grib_open_file(ftn,trim(fname), 'r',iret)
  if(iret.ne.0) then 
     write(LIS_logunit,*) '[ERR] Could not open file: ',trim(fname)
     call LIS_endrun()
  endif
  call grib_multi_support_on

  do
     call grib_new_from_file(ftn,igrib,iret)

     if ( iret == GRIB_END_OF_FILE ) then
        exit
     endif

     call grib_get(igrib, 'parameterNumber',param_num, iret)
     call LIS_verify(iret, &
          'grib_get: parameterNumber failed in readSMOPS_SMOSsm_struc')

     var_found = .false. 
     if(SMOPS_SMOSsm_struc(n)%useSMOS.eq.1) then 
        if(param_num.eq.param_SMOS) then 
           var_found = .true.
        endif
     endif

     if(var_found) then
        call grib_get(igrib, 'values',sm_SMOS,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPS_SMOSsmObs')
        
        do r=1,SMOPS_SMOSsm_struc(n)%nr
           do c=1,SMOPS_SMOSsm_struc(n)%nc
              sm_SMOS_t(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc) = &
                 sm_SMOS(c+((SMOPS_SMOSsm_struc(n)%nr-r+1)-1)*&
                 SMOPS_SMOSsm_struc(n)%nc)
           enddo
        enddo     
        
     endif

     var_found = .false. 
     if(SMOPS_SMOSsm_struc(n)%useSMOS.eq.1) then 
        if(param_num.eq.param_SMOS_qa) then 
           var_found = .true.
        endif
     endif
     
     if(var_found) then
        call grib_get(igrib, 'values',sm_SMOS_qa,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPS_SMOSsmObs')
        
        do r=1,SMOPS_SMOSsm_struc(n)%nr
           do c=1,SMOPS_SMOSsm_struc(n)%nc
              sm_SMOS_qa_t(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc) = &
                 int(sm_SMOS_qa(c+((SMOPS_SMOSsm_struc(n)%nr-r+1)-1)*&
                 SMOPS_SMOSsm_struc(n)%nc))
           enddo
        enddo       
     endif

     var_found = .false. 
     if(SMOPS_SMOSsm_struc(n)%useSMOS.eq.1) then 
        if(param_num.eq.param_SMOS_hr) then 
           var_found = .true.
        endif
     endif
     
     if(var_found) then
        call grib_get(igrib, 'values',sm_SMOS_hr,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPS_SMOSsmObs')
     endif
     
     var_found = .false. 
     if(SMOPS_SMOSsm_struc(n)%useSMOS.eq.1) then 
        if(param_num.eq.param_SMOS_mn) then 
           var_found = .true.
        endif
     endif
     
     if(var_found) then
        call grib_get(igrib, 'values',sm_SMOS_mn,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPS_SMOSsmObs')
     endif

     call grib_release(igrib,iret)
     call LIS_verify(iret, 'error in grib_release in readRTSMOPS_SMOSsmObs')
  enddo

  call grib_close_file(ftn)

  sm_time = LIS_rc%udef

  ! Table 3.6.1 – SMOPS soil moisture product Quality Assessment (QA) bits.
  ! (c) SMOS Soil Moisture Product QA
  !
  ! Byte 1:
  !
  ! Bit |  Description
  ! --------------------------------------------------------
  ! 0   |  Spare bit
  ! 1   |  1 = RFI for H pol above threshold, 0 = otherwise
  ! 2   |  1 = RFI for V pol above threshold, 0 = otherwise
  ! 3   |  Spare bit
  ! 4   |  1 = No products are generated, 0 = otherwise
  ! 5   |  1 = Retrieval values outside range, 0 = otherwise
  ! 6   |  1 = High retrieval DQX, 0 = otherwise
  ! 7   |  1 = Poor fit quality, 0 = otherwise
  !
  ! Byte 2:
  !
  ! Bit |  Description
  ! -------------------------------------------------------------
  ! 0   |  1 = Presence of other than nominal soil; 0 = otherwise
  ! 1   |  1 = Rocks; 0 = not rocks
  ! 2   |  1 = Moderate or strong topography; 0 = otherwise
  ! 3   |  1 = Open water; 0 = not open water
  ! 4   |  1 = Snow; 0 = not snow
  ! 5   |  1 = Forest; 0 = not forest
  ! 6   |  1 = Flood risk; 0 = no flood risk
  ! 7   |  1 = Urban area; 0 = not urban area
  !
  !
  ! From bytes 1 and 2, we will reject an observation if
  !
  !     bit 1 is 1 or bit 2 is 1 or bit 4 is 1 or bit 5 is 1 or
  !     bit 6 is 1 or bit 7 is 1 or bit 8 is 1 or bit 9 is 1 or
  !     bit 10 is 1 or bit 11 is 1 or bit 12 is 1 or bit 13 is 1 or
  !     bit 14 is 1 or bit 15 is 1
  !
  ! Thus we will accept an observation only if
  !
  !     bit 0 is (0|1) and bit 1 is 0 and bit 2 is 0 and
  !     bit 3 is (0|1) and bit 1 is 0 and bit 2 is 0 and
  !     bit 4 is 0 and bit 5 is 0 and bit 6 is 0 and bit 4 is 0 and
  !     bit 5 is 0 and bit 6 is 0 and bit 7 is 0 and bit 8 is 0 and
  !     bit 9 is 0 and bit 10 is 0 and bit 11 is 0 and bit 12 is 0 and
  !     bit 13 is 0 and bit 14 is 0 and bit 15 is 0
  !
  ! I.e., accept when bytes 1 and 2 are either b'0000000000000000' or
  ! b'0000000000000001' or b'0000000000001000' or b'0000000000001001';
  ! otherwise reject.
  !
  do r=1, SMOPS_SMOSsm_struc(n)%nr
     do c=1, SMOPS_SMOSsm_struc(n)%nc
        qavalue = sm_SMOS_qa_t(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc)
        if ( qavalue .ne. 9999 ) then 
           if ( qavalue == SMOS_accept1 .or. &
                qavalue == SMOS_accept2 .or. &
                qavalue == SMOS_accept3 .or. &
                qavalue == SMOS_accept4 .or. &
                qavalue == SMOS_accept5 .or. &
                qavalue == SMOS_accept6 .or. &
                qavalue == SMOS_accept7 .or. &
                qavalue == SMOS_accept8 ) then
              hr_val = nint(sm_SMOS_hr(c+&
                            ((SMOPS_SMOSsm_struc(n)%nr-r+1)-1)*&
                             SMOPS_SMOSsm_struc(n)%nc))
              mn_val =  nint(sm_SMOS_mn(c+&
                             ((SMOPS_SMOSsm_struc(n)%nr-r+1)-1)*&
                              SMOPS_SMOSsm_struc(n)%nc))
              call LIS_get_timeoffset_sec(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
                                          hr_val, mn_val, 0, julss)
              sm_time(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc) = julss
              sm_data_b(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc) = .true. 
           else
              sm_SMOS_t(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc) = LIS_rc%udef
              sm_data_b(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc) = .false.
           endif
        else
           sm_SMOS_t(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc) = LIS_rc%udef
           sm_data_b(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc) = .false. 
        endif
        if(sm_SMOS_t(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc).lt.0.001) then 
           sm_SMOS_t(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc) = LIS_rc%udef
           sm_data_b(c+(r-1)*SMOPS_SMOSsm_struc(n)%nc) = .false.
        endif
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_SMOS_t, smobs_b_ip, smobs_ip, &
       SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       SMOPS_SMOSsm_struc(n)%rlat, SMOPS_SMOSsm_struc(n)%rlon, &
       SMOPS_SMOSsm_struc(n)%n11,  LIS_rc%udef, ios)

  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_time, smobs_b_ip, smtime_ip, &
       SMOPS_SMOSsm_struc(n)%nc*SMOPS_SMOSsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       SMOPS_SMOSsm_struc(n)%rlat, SMOPS_SMOSsm_struc(n)%rlon, &
       SMOPS_SMOSsm_struc(n)%n11,  LIS_rc%udef, ios)
 endif

#endif
  
end subroutine read_RTSMOPS_SMOS_data


!BOP
! !ROUTINE: create_SMOPS_SMOSsm_filename
! \label{create_SMOPS_SMOSsm_filename}
! 
! !INTERFACE: 
subroutine create_SMOPS_SMOSsm_filename(ndir, useRT, yr, mo,da, hr, conv, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: useRT
  integer           :: yr, mo, da,hr
  character (len=*) :: ndir
  character (len=*) :: conv
! 
! !DESCRIPTION: 
!  This subroutine creates the SMOPS filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the SMOPS soil moisture data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[conv] naming convention for the SMOPS data
!  \item[filename] Generated SMOPS filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
 
  if(useRT.eq.1) then 
     if ( conv == "LIS" ) then
        filename = trim(ndir)//'/'//trim(fyr)//'/NPR_SMOPS_CMAP_D' &
                   //trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//'.gr2'     
     else
        filename = trim(ndir)//'/smops_d' &
                   //trim(fyr)//trim(fmo)//trim(fda)//'_s'//trim(fhr)//'0000_cness.gr2'
     endif
  else
     filename = trim(ndir)//'/'//trim(fyr)//'/NPR_SMOPS_CMAP_D' &
          //trim(fyr)//trim(fmo)//trim(fda)//'.gr2'
  endif

end subroutine create_SMOPS_SMOSsm_filename




