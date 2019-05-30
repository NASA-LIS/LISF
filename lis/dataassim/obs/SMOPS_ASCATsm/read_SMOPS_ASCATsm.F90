!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_SMOPS_ASCATsm
! \label{read_SMOPS_ASCATsm}
!
! !REVISION HISTORY:
!  17 Jun 2010: Sujay Kumar; Updated for use with LPRM AMSRE Version 5. 
!  20 Sep 2012: Sujay Kumar; Updated to the NETCDF version of the data. 
!  27 Sep 2017: Mahdi Navari; Updated to read ASCAT from SMOPS V3
!  15 May 2017: Eric Kemp; LIS current date/time no longer directly passed
!               to SMOPS filename subroutine to avoid accidently changing
!               the time.
!  1  Apr 2019: Yonghwan Kwon: Upated for reading monthy CDF for the current month
!
! !INTERFACE: 
subroutine read_SMOPS_ASCATsm(n, k, OBS_State, OBS_Pert_State)
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
  use SMOPS_ASCATsm_Mod, only : SMOPS_ASCATsm_struc

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
  character*100          :: smobsdir
  character*100          :: fname
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
  integer :: yyyy,mm,dd,hh ! EMK

  
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
  
  if(alarmCheck.or.SMOPS_ASCATsm_struc(n)%startMode) then 
     SMOPS_ASCATsm_struc(n)%startMode = .false.

     SMOPS_ASCATsm_struc(n)%smobs = LIS_rc%udef
     SMOPS_ASCATsm_struc(n)%smtime = -1

     smobs = LIS_rc%udef
     smtime = LIS_rc%udef


     !EMK
     yyyy = LIS_rc%yr
     mm = LIS_rc%mo
     dd = LIS_rc%da
     hh = LIS_rc%hr
     !call create_SMOPS_ASCATsm_filename(smobsdir, &
     !     SMOPS_ASCATsm_struc(n)%useRealtime, &
     !     LIS_rc%yr, LIS_rc%mo, &
     !     LIS_rc%da, LIS_rc%hr, SMOPS_ASCATsm_struc(n)%conv, fname)
     call create_SMOPS_ASCATsm_filename(smobsdir, &
          SMOPS_ASCATsm_struc(n)%useRealtime, &
          yyyy,mm,dd,hh, &
          SMOPS_ASCATsm_struc(n)%conv, fname)

     inquire(file=fname,exist=file_exists)

     if(file_exists) then 
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
        call read_SMOPS_ASCAT_data(n,k,fname,smobs,smtime)
     else
        write(LIS_logunit,*) '[WARN] Missing SMOPS ',trim(fname)
     endif

     SMOPS_ASCATsm_struc(n)%smobs  = LIS_rc%udef
     SMOPS_ASCATsm_struc(n)%smtime = -1

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           grid_index = LIS_obs_domain(n,k)%gindex(c,r)
           if(grid_index.ne.-1) then 
              if(smobs(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then 
                 SMOPS_ASCATsm_struc(n)%smobs(c,r) = &
                      smobs(c+(r-1)*LIS_rc%obs_lnc(k))                 
                 SMOPS_ASCATsm_struc(n)%smtime(c,r) = &
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

! MN: temporally comment this bluck of code to reduce the fragmentation from temporal slicing so we can rule out 1) any interpolation errors 2) errors from CDF matching 
  call LIS_get_timeoffset_sec(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
       LIS_rc%hr, LIS_rc%mn, LIS_rc%ss, lis_julss)

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           if(SMOPS_ASCATsm_struc(n)%smtime(c,r).ge.0) then 
              dt = (lis_julss-SMOPS_ASCATsm_struc(n)%smtime(c,r))
              if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
                 sm_current(c,r) = & 
                      SMOPS_ASCATsm_struc(n)%smobs(c,r)
                 fnd = 1
              endif           
           endif
        endif
     enddo
  enddo

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     

  ! Read monthly CDF (only for the current month)
  if (SMOPS_ASCATsm_struc(n)%ntimes.gt.1.and.SMOPS_ASCATsm_struc(n)%cdf_read_opt.eq.1) then
     if (.not. SMOPS_ASCATsm_struc(n)%cdf_read_mon.or.LIS_rc%da.eq.1.and.LIS_rc%hr.eq.0.and.LIS_rc%mn.eq.0.and.LIS_rc%ss.eq.0) then
        call LIS_readMeanSigmaData(n,k,&
             SMOPS_ASCATsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             SMOPS_ASCATsm_struc(n)%modelcdffile, &
             "SoilMoist",&
             SMOPS_ASCATsm_struc(n)%model_mu,&
             SMOPS_ASCATsm_struc(n)%model_sigma,&
             LIS_rc%mo)

        call LIS_readMeanSigmaData(n,k,&
             SMOPS_ASCATsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             SMOPS_ASCATsm_struc(n)%obscdffile, &
             "SoilMoist",&
             SMOPS_ASCATsm_struc(n)%obs_mu,&
             SMOPS_ASCATsm_struc(n)%obs_sigma,&
             LIS_rc%mo)

        call LIS_readCDFdata(n,k,&
             SMOPS_ASCATsm_struc(n)%nbins,&
             SMOPS_ASCATsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             SMOPS_ASCATsm_struc(n)%modelcdffile, &
             "SoilMoist",&
             SMOPS_ASCATsm_struc(n)%model_xrange,&
             SMOPS_ASCATsm_struc(n)%model_cdf,&
             LIS_rc%mo)

        call LIS_readCDFdata(n,k,&
             SMOPS_ASCATsm_struc(n)%nbins,&
             SMOPS_ASCATsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             SMOPS_ASCATsm_struc(n)%obscdffile, &
             "SoilMoist",&
             SMOPS_ASCATsm_struc(n)%obs_xrange,&
             SMOPS_ASCATsm_struc(n)%obs_cdf,&
             LIS_rc%mo)

        SMOPS_ASCATsm_struc(n)%cdf_read_mon = .true.
     endif
  endif

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

     if (SMOPS_ASCATsm_struc(n)%ntimes.gt.1.and.SMOPS_ASCATsm_struc(n)%cdf_read_opt.eq.1) then
        call LIS_rescale_with_CDF_matching(     &
             n,k,                               &
             SMOPS_ASCATsm_struc(n)%nbins,         &
             1,                                 &
             MAX_SM_VALUE,                      &
             MIN_SM_VALUE,                      &
             SMOPS_ASCATsm_struc(n)%model_xrange,  &
             SMOPS_ASCATsm_struc(n)%obs_xrange,    &
             SMOPS_ASCATsm_struc(n)%model_cdf,     &
             SMOPS_ASCATsm_struc(n)%obs_cdf,       &
             sm_current)
     else
        call LIS_rescale_with_CDF_matching(    &
             n,k,                              & 
             SMOPS_ASCATsm_struc(n)%nbins,         & 
             SMOPS_ASCATsm_struc(n)%ntimes,         & 
             MAX_SM_VALUE,                        & 
             MIN_SM_VALUE,                        & 
             SMOPS_ASCATsm_struc(n)%model_xrange,  &
             SMOPS_ASCATsm_struc(n)%obs_xrange,    &
             SMOPS_ASCATsm_struc(n)%model_cdf,     &
             SMOPS_ASCATsm_struc(n)%obs_cdf,       &
             sm_current)                              
     endif  
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
       //trim(LIS_SMOPS_ASCATsmobsId)//char(0),n,k,OBS_state)

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

     if(SMOPS_ASCATsm_struc(n)%useSsdevScal.eq.1) then
        call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')
        
        allocate(ssdev(LIS_rc%obs_ngrid(k)))
        ssdev = SMOPS_ASCATsm_struc(n)%ssdev_inp 

        if (SMOPS_ASCATsm_struc(n)%cdf_read_opt .eq. 1) then
           jj = 1
        else if(SMOPS_ASCATsm_struc(n)%ntimes.eq.1) then 
           jj = 1
        else
           jj = LIS_rc%mo
        endif
        do t=1,LIS_rc%obs_ngrid(k)
           if(SMOPS_ASCATsm_struc(n)%obs_sigma(t,jj).gt.0) then 
              ssdev(t) = ssdev(t)*SMOPS_ASCATsm_struc(n)%model_sigma(t,jj)/&
                   SMOPS_ASCATsm_struc(n)%obs_sigma(t,jj)
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

end subroutine read_SMOPS_ASCATsm

!BOP
! 
! !ROUTINE: read_SMOPS_ASCAT_data
! \label{read_SMOPS_ASCAT_data}
!
! !INTERFACE:
subroutine read_SMOPS_ASCAT_data(n, k, fname, smobs_ip, smtime_ip)
! 
! !USES:   
#if(defined USE_GRIBAPI)
  use grib_api
#endif
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use SMOPS_ASCATsm_Mod, only : SMOPS_ASCATsm_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                          :: smtime_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the RTSMOPS grib2 file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  the estimated error is above a predefined threshold (the recommeded
!  value is 5%). 
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
  real, parameter :: err_threshold = 5 ! in percent
  integer         :: param_ASCAT_A, param_ASCAT_A_qa
  integer         :: param_ASCAT_A_hr, param_ASCAT_A_mn

  integer         :: param_ASCAT_B, param_ASCAT_B_qa
  integer         :: param_ASCAT_B_hr, param_ASCAT_B_mn

  real            :: sm_ASCAT_A(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  real            :: sm_ASCAT_A_t(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  real            :: sm_ASCAT_A_hr(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  real            :: sm_ASCAT_A_mn(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)

  real            :: sm_ASCAT_B(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  real            :: sm_ASCAT_B_t(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  real            :: sm_ASCAT_B_hr(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  real            :: sm_ASCAT_B_mn(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)

  real            :: sm_ASCAT_A_qa(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  integer*2       :: sm_ASCAT_A_qa_t(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)

  real            :: sm_ASCAT_B_qa(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  integer*2       :: sm_ASCAT_B_qa_t(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)

  real            :: sm_time_ASCAT_A(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  real            :: sm_time_ASCAT_B(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  real            :: sm_time(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  real            :: sm_data(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  logical*1       :: sm_data_b(SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr)
  logical*1       :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer         :: hr_val, mn_val, julss
  integer         :: c,r,ios
  integer         :: ftn,iret,igrib,nvars
  integer         :: param_num
  logical         :: var_found
  integer*2       :: qavalue
  integer*1       :: err, ql
  integer         :: updoy,yr1,mo1,da1,hr1,mn1,ss1,offset
  real            :: upgmt
  real*8          :: file_time


#if(defined USE_GRIBAPI)
  ! When we are reading the 6-hourly datasets, we read the file HR+6
  ! because it contains the previous 6 hours of retrievals.
  ! Thus we need to add 6 hours to the LIS time do determine the version
  ! of the SMOPS dataset.
  ! Otherwise, when reading the daily datasets, there is no need to
  ! adjust time.
  if ( SMOPS_ASCATsm_struc(n)%useRealtime == 1 ) then
     offset = 6
  else
     offset = 0
  endif

  if ( SMOPS_ASCATsm_struc(n)%version == '1.3' ) then
     file_time = SMOPS_ASCATsm_struc(n)%version2_time - 1.0
  elseif ( SMOPS_ASCATsm_struc(n)%version == '2.0' ) then
     file_time = SMOPS_ASCATsm_struc(n)%version2_time
  elseif ( SMOPS_ASCATsm_struc(n)%version == '3.0' ) then
     file_time = SMOPS_ASCATsm_struc(n)%version3_time
  else
     yr1 = LIS_rc%yr
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = LIS_rc%hr + offset
     mn1 = LIS_rc%mn
     ss1 = 0
     call LIS_date2time(file_time,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)
  endif

  if ( file_time < SMOPS_ASCATsm_struc(n)%version3_time ) then
     ! SMOPS version 1.3 or 2.0
     param_ASCAT_A = 213; param_ASCAT_A_qa = 234
     param_ASCAT_A_hr = 223; param_ASCAT_A_mn = 224
     param_ASCAT_B = 214; param_ASCAT_B_qa = 235
     param_ASCAT_B_hr = 225; param_ASCAT_B_mn = 226
     write(LIS_logunit,*) '[MSG] Reading SMOPS ASCAT dataset '//&
        'as SMOPS version 1.3/2.0'

  else ! ( file_time >= SMOPS_ASCATsm_struc(n)%version3_time ) then
     ! SMOPS version 3.0
     param_ASCAT_A = 213; param_ASCAT_A_qa = 243
     param_ASCAT_A_hr = 226; param_ASCAT_A_mn = 227
     param_ASCAT_B = 214; param_ASCAT_B_qa = 244
     param_ASCAT_B_hr = 228; param_ASCAT_B_mn = 229
     write(LIS_logunit,*) '[MSG] Reading SMOPS ASCAT dataset '//&
        'as SMOPS version 3.0'
  endif


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
          'grib_get: parameterNumber failed in readSMOPSsm_struc')

     var_found = .false.
     if(param_num.eq.param_ASCAT_A) then
        var_found = .true.
     endif

     if(var_found) then
        call grib_get(igrib, 'values',sm_ASCAT_A,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPSsmObs')

        do r=1,SMOPS_ASCATsm_struc(n)%nr
           do c=1,SMOPS_ASCATsm_struc(n)%nc
              sm_ASCAT_A_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = &
                   sm_ASCAT_A(c+((SMOPS_ASCATsm_struc(n)%nr-r+1)-1)*&
                   SMOPS_ASCATsm_struc(n)%nc)
           enddo
        enddo

     endif

     var_found = .false.
     if(param_num.eq.param_ASCAT_A_qa) then
        var_found = .true.
     endif

     if(var_found) then
        call grib_get(igrib, 'values',sm_ASCAT_A_qa,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPSsmObs')

        do r=1,SMOPS_ASCATsm_struc(n)%nr
           do c=1,SMOPS_ASCATsm_struc(n)%nc
              sm_ASCAT_A_qa_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = &
                   INT(sm_ASCAT_A_qa(c+((SMOPS_ASCATsm_struc(n)%nr-r+1)-1)*&
                   SMOPS_ASCATsm_struc(n)%nc))
           enddo
        enddo
     endif

     var_found = .false.
     if(param_num.eq.param_ASCAT_A_hr) then
        var_found = .true.
     endif

     if(var_found) then
        call grib_get(igrib, 'values',sm_ASCAT_A_hr,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPSsmObs')
     endif

     var_found = .false.
     if(param_num.eq.param_ASCAT_A_mn) then
        var_found = .true.
     endif

     if(var_found) then
        call grib_get(igrib, 'values',sm_ASCAT_A_mn,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPSsmObs')
     endif

     var_found = .false.
     if(param_num.eq.param_ASCAT_B) then
        var_found = .true.
     endif

     if(var_found) then
        call grib_get(igrib, 'values',sm_ASCAT_B,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPSsmObs')

        do r=1,SMOPS_ASCATsm_struc(n)%nr
           do c=1,SMOPS_ASCATsm_struc(n)%nc
              sm_ASCAT_B_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = &
                   sm_ASCAT_B(c+((SMOPS_ASCATsm_struc(n)%nr-r+1)-1)*&
                   SMOPS_ASCATsm_struc(n)%nc)
           enddo
        enddo

     endif

     var_found = .false.
     if(param_num.eq.param_ASCAT_B_qa) then
        var_found = .true.
     endif

     if(var_found) then
        call grib_get(igrib, 'values',sm_ASCAT_B_qa,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPSsmObs')

        do r=1,SMOPS_ASCATsm_struc(n)%nr
           do c=1,SMOPS_ASCATsm_struc(n)%nc
              sm_ASCAT_B_qa_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = &
                   INT(sm_ASCAT_B_qa(c+((SMOPS_ASCATsm_struc(n)%nr-r+1)-1)*&
                   SMOPS_ASCATsm_struc(n)%nc))
           enddo
        enddo
     endif

     var_found = .false.
     if(param_num.eq.param_ASCAT_B_hr) then
        var_found = .true.
     endif

     if(var_found) then
        call grib_get(igrib, 'values',sm_ASCAT_B_hr,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPSsmObs')
     endif

     var_found = .false.
     if(param_num.eq.param_ASCAT_B_mn) then
        var_found = .true.
     endif

     if(var_found) then
        call grib_get(igrib, 'values',sm_ASCAT_B_mn,iret)
        call LIS_warning(iret,'error in grib_get:values in readRTSMOPSsmObs')
     endif

     call grib_release(igrib,iret)
     call LIS_verify(iret, 'error in grib_release in readRTSMOPSsmObs')
  enddo

  call grib_close_file(ftn)

  ! Table 3.6.1 – SMOPS soil moisture product Quality Assessment (QA) bits.
  ! (d) ASCAT Soil Moisture Product QA
  !
  ! Byte |  Description
  ! -----------------------------------------------------------------------
  ! 0    |  Estimated Error in Soil Moisture. (Integer. Scale factor: 0.01)
  ! 1    |  Soil Moisture Quality (Integer, Scale factor: 0.01)
  !
  ! The retrievals are rejected when the estimated error is above
  ! a predefined threshold (the recommeded value is 5%).
  !
  ! Technically speaking, err_threshold should be 0.05 and
  ! err should be scaled by 0.01.  But below is consistent with the NESDIS
  ! documentation.
  !
  ! Note that I am assuming that Byte 0 above refers to the least
  ! significant byte and Byte 1 above refers to the most significant byte.
  ! Meaning estimated error is get_byte1, and quality is get_byte2.
  sm_time_ASCAT_A = LIS_rc%udef
  sm_time_ASCAT_B = LIS_rc%udef
  sm_data         = LIS_rc%udef
  sm_data_b       = .false.

  do r=1, SMOPS_ASCATsm_struc(n)%nr
     do c=1, SMOPS_ASCATsm_struc(n)%nc
        qavalue = sm_ASCAT_A_qa_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc)
        if ( qavalue .ne. 9999 ) then
           !estimated error
           err = get_byte1(qavalue)
           !quality flag - not used currently
           ql = get_byte2(qavalue)

           if(err.lt.err_threshold) then
              hr_val = nint(sm_ASCAT_A_hr(c+&
                   ((SMOPS_ASCATsm_struc(n)%nr-r+1)-1)*&
                   SMOPS_ASCATsm_struc(n)%nc))
              mn_val =  nint(sm_ASCAT_A_mn(c+&
                   ((SMOPS_ASCATsm_struc(n)%nr-r+1)-1)*&
                   SMOPS_ASCATsm_struc(n)%nc))
              call LIS_get_timeoffset_sec(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
                   hr_val, mn_val, 0, julss)
              sm_time_ASCAT_A(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = julss
              sm_data_b(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = .true.
           else
              sm_ASCAT_A_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = LIS_rc%udef
              sm_data_b(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = .false.
           endif
        else
           sm_ASCAT_A_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = LIS_rc%udef
           sm_data_b(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = .false.
        endif
        if(sm_ASCAT_A_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc).lt.0.001) then
           sm_ASCAT_A_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = LIS_rc%udef
           sm_data_b(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = .false.
        endif
     enddo
  enddo

  do r=1, SMOPS_ASCATsm_struc(n)%nr
     do c=1, SMOPS_ASCATsm_struc(n)%nc
        qavalue = sm_ASCAT_B_qa_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc)
        if ( qavalue .ne. 9999 ) then
           !estimated error
           err = get_byte1(qavalue)
           !quality flag - not used currently
           ql = get_byte2(qavalue)

           if(err.lt.err_threshold) then
              hr_val = nint(sm_ASCAT_B_hr(c+&
                   ((SMOPS_ASCATsm_struc(n)%nr-r+1)-1)*&
                   SMOPS_ASCATsm_struc(n)%nc))
              mn_val =  nint(sm_ASCAT_B_mn(c+&
                   ((SMOPS_ASCATsm_struc(n)%nr-r+1)-1)*&
                   SMOPS_ASCATsm_struc(n)%nc))
              call LIS_get_timeoffset_sec(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
                   hr_val, mn_val, 0, julss)
              sm_time_ASCAT_B(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = julss
           else
              sm_ASCAT_B_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = LIS_rc%udef
           endif
        else
           sm_ASCAT_B_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = LIS_rc%udef
        endif
        if(sm_ASCAT_B_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc).lt.0.001) then
           sm_ASCAT_B_t(c+(r-1)*SMOPS_ASCATsm_struc(n)%nc) = LIS_rc%udef
        endif
     enddo
  enddo

  sm_data = sm_ASCAT_A_t
  sm_time = sm_time_ASCAT_A
  where ( sm_ASCAT_B_t /= LIS_rc%udef )
     sm_data = sm_ASCAT_B_t
     sm_time = sm_time_ASCAT_B
     sm_data_b = .true.
  endwhere

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!--------------------------------------------------------------------------
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       SMOPS_ASCATsm_struc(n)%rlat, SMOPS_ASCATsm_struc(n)%rlon, &
       SMOPS_ASCATsm_struc(n)%n11,  LIS_rc%udef, ios)

  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_time, smobs_b_ip, smtime_ip, &
       SMOPS_ASCATsm_struc(n)%nc*SMOPS_ASCATsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       SMOPS_ASCATsm_struc(n)%rlat, SMOPS_ASCATsm_struc(n)%rlon, &
       SMOPS_ASCATsm_struc(n)%n11,  LIS_rc%udef, ios)

#endif

  contains

  integer*1 function get_byte1(i)
     implicit none
     integer*2, intent(in) :: i
     integer*2 :: j
     ! This function expects a 16-bit integer as input.  It returns
     ! the least significant byte (as an 8-bit integer), referred
     ! to as byte1 in the NESDIS documention.
     !
     ! For example,
     ! i = b'0000001000000001' <--> 00000010|00000001
     !                         <--> byte2|byte1
     !                         <--> 0x0201
     ! Here byte1 is b'00000001'; byte2 is b'00000010'
     j = iand(i, z'00ff')
     get_byte1 = j
  end function get_byte1

  integer*1 function get_byte2(i)
     implicit none
     integer*2, intent(in) :: i
     integer*2 :: j
     ! This function expects a 16-bit integer as input.  It returns
     ! the most significant byte (as an 8-bit integer), referred
     ! to as byte2 in the NESDIS documention.
     !
     ! For example,
     ! i = b'0000001000000001' <--> 00000010|00000001
     !                         <--> byte2|byte1
     !                         <--> 0x0201
     ! Here byte1 is b'00000001'; byte2 is b'00000010'
     j = ishft(i, -8)
     get_byte2 = j
  end function get_byte2
end subroutine read_SMOPS_ASCAT_data


!BOP
! !ROUTINE: create_SMOPS_ASCATsm_filename
! \label{create_SMOPS_ASCATsm_filename}
! 
! !INTERFACE: 
subroutine create_SMOPS_ASCATsm_filename(ndir, useRT, yr, mo,da, hr, conv, filename)
! !USES:   
  use ESMF
  use LIS_timeMgrMod, only : LIS_calendar

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: useRT
  integer           :: yr, mo, da, hr
  character (len=*) :: ndir
  character (len=*) :: conv
! 
! !DESCRIPTION: 
!  This subroutine creates the SMOPS filename based on the time and date 
!
!  On 2017-10-05, AGRMET ops switched from using the 
!  smops\_dYYYYMMDD\_sHH0000\_cness.gr2 naming convention
!  to using the NPR\_SMOPS\_CMAP\_DYYYYMMDDHH.gr2 naming convention.
!
!  The 6-hourly files contain retrievals for the past 6 hours.
!  For example, NPR_SMOPS_CMAP_D2017101218.gr2 contains retrievals
!  for 12 Oct 2017, hours 12, 13, 14, 15, 16, and 17.
!  So, when LIS is at hour 12, it must read the SMOPS hour 18 file to
!  process hours 12, 13, 14, 15, 16, and 17.
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

  logical, save           :: first_time=.true.
  type(ESMF_Time), save   :: naming3_time
  type(ESMF_Time)         :: file_time
  type(ESMF_TimeInterval) :: six_hours
  integer                 :: rc

  if ( first_time ) then
     call ESMF_TimeSet(naming3_time, &
        yy = 2017, &
        mm = 10,   &
        dd = 5,    &
        h  = 0,    &
        m  = 0,    &
        s  = 0,    &
        calendar = LIS_calendar, & 
        rc = rc)
     first_time = .false.
  endif

  if(useRT.eq.1) then 
     call ESMF_TimeIntervalSet(six_hours,h=6,rc=rc)
     call ESMF_TimeSet(file_time, &
        yy = yr, &
        mm = mo, &
        dd = da, &
        h  = hr, &
        m  = 0,  &
        s  = 0,  & 
        calendar = LIS_calendar, & 
        rc = rc)
     file_time = file_time + six_hours
     call ESMF_TimeGet(file_time, &
        yy = yr, &
        mm = mo, &
        dd = da, &
        h  = hr, &
        rc = rc)
  endif

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr

  if(useRT.eq.1) then 
     if ( conv == "LIS" ) then
        filename = trim(ndir)//'/'//trim(fyr)//'/NPR_SMOPS_CMAP_D' &
                   //trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//'.gr2'     
     else
        if ( file_time >= naming3_time) then
           filename = trim(ndir)//'/'//'/NPR_SMOPS_CMAP_D' &
              //trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//'.gr2'     
        else
           filename = trim(ndir)//'/smops_d' &
              //trim(fyr)//trim(fmo)//trim(fda)//'_s'//trim(fhr)//'0000_cness.gr2'
        endif
     endif
  else
     filename = trim(ndir)//'/'//trim(fyr)//'/NPR_SMOPS_CMAP_D' &
          //trim(fyr)//trim(fmo)//trim(fda)//'.gr2'
  endif

end subroutine create_SMOPS_ASCATsm_filename

