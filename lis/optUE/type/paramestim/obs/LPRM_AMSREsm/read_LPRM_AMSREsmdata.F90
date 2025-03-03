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
! !ROUTINE: read_LPRM_AMSREsm_obsdata
! \label{read_LPRM_AMSREsm_obsdata}
!
! !REVISION HISTORY:
!  31 Jan 2012: Ken Harrison; Initial Specification
!
! !INTERFACE: 
  subroutine read_LPRM_AMSREsm_obsdata(Obj_Space)
! !USES: 
    use ESMF
    use LIS_mpiMod
    use LIS_coreMod,  only : LIS_rc, LIS_domain, LIS_npes, LIS_localPet
    use LIS_timeMgrMod
    use LIS_logMod,     only : LIS_logunit, LIS_verify, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_fileIOMod,      only : LIS_readData
    use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
    use LPRM_AMSREsm_obsMod, only : LPRM_AMSREsm_obs_struc
    use map_utils

    implicit none
! !ARGUMENTS: 
    type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP

    real,    pointer    :: obse(:)
!!!    integer, parameter  :: numchannels=6
!!!    integer, parameter  :: numpolarizations=2
!!!    type(ESMF_Field)    :: smField
!!!    character*100       :: smobsdir 
!!!    logical             :: data_update
!!!    integer             :: status 
!!!    logical             :: found
!!!    integer             :: c,r
!!!    type(ESMF_TimeInterval) :: delta_t
!!!    type(ESMF_Time)         :: lis_time1
!!!    integer             :: t
!!!    integer             :: n 
!!!!    logical             :: is_ascend_pass_hr
!!!    logical             :: ob_in_curr_hr
!!!    integer             :: day_index
!!!    logical             :: is_overpass_hr
!!!    integer, parameter :: freq=1 ! start with 7
!!!    integer, parameter :: polar=2 ! start with 7H
!!!    integer            :: ipass
!!!    real                      :: ob
!!!    integer :: i, gid

  real,  parameter       :: MAX_SM_VALUE=0.45, MIN_SM_VALUE=0.0001
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname_A, fname_D
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: smfield

  integer, allocatable       :: gid(:)
  integer, allocatable       :: assimflag(:)
  real, allocatable          :: smobs_A(:)
  real, allocatable          :: smobs_D(:)
  real, allocatable          :: sm_current(:,:)
  real, allocatable          :: model_delta(:)
  real, allocatable          :: obs_delta(:)

  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd
  real                   :: dt
  real                   :: lon
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
  integer                :: fnd
  real                   :: smvalue
  integer                :: n  ! in DA reader, n is passed in

    n = 1
    allocate(gid(LIS_rc%ngrid(n)))
    allocate(assimflag(LIS_rc%ngrid(n)))
    allocate(smobs_A(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
    allocate(smobs_D(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
    allocate(sm_current(LIS_rc%lnc(n),LIS_rc%lnr(n)))
    allocate(model_delta(LIS_rc%ngrid(n)))
    allocate(obs_delta(LIS_rc%ngrid(n)))
!!!    !initialize some values
!!!    is_overpass_hr=.false.
!!!    found=.false.
    call ESMF_AttributeGet(Obj_Space,"Data Directory",&
         smobsdir, rc=status)
    call LIS_verify(status)
    call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
         data_upd, rc=status)
    call LIS_verify(status)
    
    data_upd = .false. 
    !-------------------------------------------------------------------------
    !   Read both ascending and descending passes at 0Z and then store
    !   the overpass time as 1.30AM for the descending pass and 1.30PM 
    !   for the ascending pass. 
    !-------------------------------------------------------------------------
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "LPRM AMSRE soil moisture read alarm")
    
    if(alarmCheck.or.LPRM_AMSREsm_obs_struc(n)%startMode) then 
       LPRM_AMSREsm_obs_struc(n)%startMode = .false.
       
       LPRM_AMSREsm_obs_struc(n)%smobs = LIS_rc%udef
       smobs_A = LIS_rc%udef
     smobs_D = LIS_rc%udef             
     LPRM_AMSREsm_obs_struc(n)%smtime = -1
     
     call create_LPRM_AMSREsm_obs_filename(smobsdir, 'A',&
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname_A)
     
     inquire(file=fname_A,exist=file_exists)

     if(file_exists) then 
        
        write(LIS_logunit,*) 'Reading ',trim(fname_A)
        call read_LPRM_obs_data(n,fname_A,smobs_A)
     endif

     call create_LPRM_AMSREsm_obs_filename(smobsdir, 'D',&
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname_D)

     inquire(file=fname_D,exist=file_exists)

     if(file_exists) then 
        
        write(LIS_logunit,*) 'Reading ',trim(fname_D)
        call read_LPRM_obs_data(n,fname_D,smobs_D)

     endif

     LPRM_AMSREsm_obs_struc(n)%smobs  = LIS_rc%udef
     LPRM_AMSREsm_obs_struc(n)%smtime = -1


     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           grid_index = LIS_domain(n)%gindex(c,r)
           if(grid_index.ne.-1) then 
              if(smobs_D(c+(r-1)*LIS_rc%lnc(n)).ne.-9999.0) then               
                 LPRM_AMSREsm_obs_struc(n)%smobs(c,r) = &
                      smobs_D(c+(r-1)*LIS_rc%lnc(n))
                 
                 lon = LIS_domain(n)%grid(grid_index)%lon
                 lhour = 1.5
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 LPRM_AMSREsm_obs_struc(n)%smtime(c,r) = gmt
              endif

              if(smobs_A(c+(r-1)*LIS_rc%lnc(n)).ne.-9999.0) then               
                 LPRM_AMSREsm_obs_struc(n)%smobs(c,r) = &
                      smobs_A(c+(r-1)*LIS_rc%lnc(n))                 
                 lon = LIS_domain(n)%grid(grid_index)%lon
                 lhour = 13.5 
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 LPRM_AMSREsm_obs_struc(n)%smtime(c,r) = gmt
              endif
           endif
        enddo
     enddo
  endif
  
  
  call ESMF_StateGet(Obj_Space,"LPRM AMSRE Surface Soil Moisture",smfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
  fnd = 0 
  sm_current = LIS_rc%udef
 
! dt is not defined as absolute value of the time difference to avoid
! double counting of the data in assimilation. 

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(LPRM_AMSREsm_obs_struc(n)%smtime(c,r).ge.0) then 
              dt = (LPRM_AMSREsm_obs_struc(n)%smtime(c,r) - LIS_rc%gmt)*3600.0
              if(dt.ge.0.and.dt.le.LIS_rc%ts) then 
                 sm_current(c,r) = & 
                      LPRM_AMSREsm_obs_struc(n)%smobs(c,r)
                 fnd = 1
              endif           
           endif
        endif
     enddo
  enddo

  obsl = LIS_rc%udef 
  if(fnd.ne.0) then 
     do r=1, LIS_rc%lnr(n)
        do c=1, LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              obsl(LIS_domain(n)%gindex(c,r))=sm_current(c,r)
           endif
        enddo
     enddo
  endif

!  NEED TO ADD CONCEPT TO PE OBS AS WELL.  OR HAVE GENERIC OBS HANDLER
!!!!!!lsm based qc
!!!!!  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
!!!!!       //trim(LIS_LPRM_AMSREsmobsId)//char(0),n, OBS_state)

!!!  call ESMF_StateGet(Obj_Space,"Observation01",smField,&
!!!       rc=status)
!!!  call LIS_verify(status)
!!!
!!!  call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
!!!  call LIS_verify(status)
!!!
!!!  do r =1,LIS_rc%lnr(n)
!!!     do c =1,LIS_rc%lnc(n)
!!!        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
!!!           sm_current(c,r) = obsl(LIS_domain(n)%gindex(c,r))
!!!        end if
!!!     end do
!!!  end do

!!!!-------------------------------------------------------------------------
!!!!  Transform data to the LSM climatology using a CDF-scaling approach
!!!!-------------------------------------------------------------------------     
!!!
!!!  if(LPRM_AMSREsm_obs_struc(n)%scal.ne.0.and.fnd.ne.0) then        
!!!     call LIS_rescale_with_CDF_matching(    &
!!!          n,                                   & 
!!!          LPRM_AMSREsm_obs_struc(n)%nbins,         & 
!!!          MAX_SM_VALUE,                        & 
!!!          MIN_SM_VALUE,                        & 
!!!          LPRM_AMSREsm_obs_struc(n)%model_xrange,  &
!!!          LPRM_AMSREsm_obs_struc(n)%obs_xrange,    &
!!!          LPRM_AMSREsm_obs_struc(n)%model_cdf,     &
!!!          LPRM_AMSREsm_obs_struc(n)%obs_cdf,       &
!!!          sm_current)
!!!
!!!  endif
!!!
!!!  fnd = 0 
!!!  data_upd_flag(LIS_localPet + 1) = .false.   
!!!  do r =1,LIS_rc%lnr(n)
!!!     do c =1,LIS_rc%lnc(n)
!!!        if(sm_current(c,r).ne.LIS_rc%udef) then
!!!           fnd = 1
!!!        endif
!!!     enddo
!!!  enddo
!!!
!!!  call ESMF_StateGet(Obj_Space,"Observation01",smField,&
!!!       rc=status)
!!!  call LIS_verify(status)
!!!
!!!  call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
!!!  call LIS_verify(status)
!!!  
!!!  obsl = LIS_rc%udef  
!!!  
!!!  if(fnd.eq.0) then
!!!     obsl = LIS_rc%udef
!!!  else
!!!     do r =1,LIS_rc%lnr(n)
!!!        do c =1,LIS_rc%lnc(n)
!!!           if (LIS_domain(n)%gindex(c,r) .ne. -1)then
!!!              obsl(LIS_domain(n)%gindex(c,r))=sm_current(c,r)
!!!           end if
!!!        end do
!!!     end do
!!!  endif

  if(fnd.eq.0) then 
     data_upd_flag(LIS_localPet+1) = .false. 
  else
     data_upd_flag(LIS_localPet+1) = .true. 
  endif
        
#if (defined SPMD)
  call MPI_ALLGATHER(data_upd_flag(LIS_localPet+1),1, &
       MPI_LOGICAL, data_upd_flag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
  data_upd = .false.
  do p=1,LIS_npes
     data_upd = data_upd.or.data_upd_flag(p)
  enddo
  
  if(data_upd) then 
!!!     do t=1,LIS_rc%ngrid(n)
!!!        gid(t) = t
!!!        if(obsl(t).ne.-9999.0) then 
!!!           assimflag(t) = 1
!!!        else
!!!           assimflag(t) = 0
!!!        endif
!!!     enddo
     
     call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
          .true. , rc=status)
     call LIS_verify(status)
     
!!!     if(LIS_rc%ngrid(n).gt.0) then 
!!!        call ESMF_AttributeSet(smField,"Grid Number",&
!!!             gid,itemCount=LIS_rc%ngrid(n),rc=status)
!!!        call LIS_verify(status)
!!!        
!!!        call ESMF_AttributeSet(smField,"Assimilation Flag",&
!!!             assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
!!!        call LIS_verify(status)
!!!     endif
  else
     call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif
deallocate(gid)
deallocate(assimflag)
deallocate(smobs_A)
deallocate(smobs_D)
deallocate(sm_current)
deallocate(model_delta)
deallocate(obs_delta)

end subroutine read_LPRM_AMSREsm_obsdata

!BOP
! 
! !ROUTINE: read_LPRM_obs_data
! \label{read_LPRM_obs_data}
!
! !INTERFACE:
subroutine read_LPRM_obs_data(n, fname, smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,   only : LIS_verify
  use map_utils,    only : latlon_to_ij
  use LPRM_AMSREsm_obsMod, only : LPRM_AMSREsm_obs_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname
  real                          :: smobs_ip(LIS_rc%lnc(n)*LIS_rc%lnr(n))


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the LPRM NETCDF file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  land surface temperature is below freezing, if rain is present, if 
!  RFI is present, if residual error is above 0.5 or if optical depth
!  is above 0.8. Finally the routine combines both the C-band and X-band
!  retrievals. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the LPRM AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  real                        :: rfi(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: tskin(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: rainf(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: optc(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: optx(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: smerrc(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: smerrx(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: smc(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: smx(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: sm_combined(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: sm_flags(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  real                        :: sm_cdf(LPRM_AMSREsm_obs_struc(n)%lprmnc,&
       LPRM_AMSREsm_obs_struc(n)%lprmnr)


  real                        :: sm_data(LPRM_AMSREsm_obs_struc(n)%lprmnc* & 
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  logical*1                   :: sm_data_b(LPRM_AMSREsm_obs_struc(n)%lprmnc* & 
       LPRM_AMSREsm_obs_struc(n)%lprmnr)
  logical*1                   :: smobs_b_ip(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  integer                     :: c,r,i,j
  real                        :: rlat,rlon,ri,rj
  integer                     :: nid, rfiid
  integer                     :: smcombId,smcdfId,smfId
  integer                     :: tskinId, rainfId, optCId, optXid
  integer                     :: smerrCid, smerrXid, smcId, smXId
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LIS_verify(ios,'Error opening file '//trim(fname))
  
  ios = nf90_inq_varid(nid, 'RFI_Daily',rfiid)
  call LIS_verify(ios, 'Error nf90_inq_varid: RFI_Daily')
  
  ios = nf90_inq_varid(nid, 'Land_Surface_Temperature',tskinid)
  call LIS_verify(ios, 'Error nf90_inq_varid: Land_Surface_Temperature')
  
  ios = nf90_inq_varid(nid, 'Rainfall_Diagnostic',rainfid)
  call LIS_verify(ios, 'Error nf90_inq_varid: Rainfall_Diagnostic')
  
  ios = nf90_inq_varid(nid, 'Optical_Depth_from_C_band',optCid)
  call LIS_verify(ios, 'Error nf90_inq_varid: Optical_Depth_from_C_band')
  
  ios = nf90_inq_varid(nid, 'Optical_Depth_from_X_band',optXid)
  call LIS_verify(ios, 'Error nf90_inq_varid: Optical_Depth_from_X_band')
  
  ios = nf90_inq_varid(nid, 'Soil_Moisture_Error_from_C_band',smerrCid)
  call LIS_verify(ios, &
       'Error nf90_inq_varid: Soil_Moisture_Error_from_C_band')
  
  ios = nf90_inq_varid(nid, 'Soil_Moisture_Error_from_X_band',smerrXid)
  call LIS_verify(ios, &
       'Error nf90_inq_varid: Soil_Moisture_Error_from_X_band')
  
  ios = nf90_inq_varid(nid, 'Soil_Moisture_from_C_band',smCid)
  call LIS_verify(ios, 'Error nf90_inq_varid: Soil_Moisture_from_C_band')
  
  ios = nf90_inq_varid(nid, 'Soil_Moisture_from_X_band',smXid)
  call LIS_verify(ios, 'Error nf90_inq_varid: Soil_Moisture_from_X_band')
  
  ios = nf90_inq_varid(nid, 'SM_Flags',smfid)
  call LIS_verify(ios, 'Error nf90_inq_varid: SM_Flags')

  ios = nf90_inq_varid(nid, 'SM_Combined',smcombid)
  call LIS_verify(ios, 'Error nf90_inq_varid: SM_Combined')

!  if(LPRM_AMSREsm_obs_struc(n)%rawdata.eq.0) then 
     ios = nf90_inq_varid(nid, 'SM_CDF',smcdfid)
     call LIS_verify(ios, 'Error nf90_inq_varid: SM_CDF')
!  endif

  !values
  ios = nf90_get_var(nid,rfiid, rfi)
  call LIS_verify(ios, 'Error nf90_get_var: RFI_Daily')
  
  ios = nf90_get_var(nid,tskinid, tskin)
  call LIS_verify(ios, 'Error nf90_get_var: Land_Surface_Temperature')
  
  ios = nf90_get_var(nid, rainfid, rainf)
  call LIS_verify(ios, 'Error nf90_get_var: Rainfall_Diagnostic')
  
  ios = nf90_get_var(nid, optcid, optc)
  call LIS_verify(ios, 'Error nf90_get_var: Optical_Depth_from_C_band')
  
  ios = nf90_get_var(nid, optxid, optx)
  call LIS_verify(ios, 'Error nf90_get_var: Optical_Depth_from_X_band')
  
  ios = nf90_get_var(nid, smerrcid, smerrc)
  call LIS_verify(ios, 'Error nf90_get_var: Soil_Moisture_Error_from_C_band')
  
  ios = nf90_get_var(nid, smerrxid,smerrx)
  call LIS_verify(ios, 'Error nf90_get_var: Soil_Moisture_Error_from_X_band')
  
  ios = nf90_get_var(nid, smfid, sm_flags)
  call LIS_verify(ios, 'Error nf90_get_var: SM_Flags')

  ios = nf90_get_var(nid, smcombid, sm_combined)
  call LIS_verify(ios, 'Error nf90_get_var: SM_Combined')

!  if(LPRM_AMSREsm_obs_struc(n)%rawdata.eq.0) then 
     ios = nf90_get_var(nid, smcdfid, sm_cdf)
     call LIS_verify(ios, 'Error nf90_get_var: SM_CDF')
!  endif
  

  ios = nf90_close(ncid=nid)
  call LIS_verify(ios,'Error closing file '//trim(fname))

  smc = 1.0
  smx = 1.0

  do r=1, LPRM_AMSREsm_obs_struc(n)%lprmnr
     do c=1, LPRM_AMSREsm_obs_struc(n)%lprmnc
!--------------------------------------------------------------------------
! Choose soil moisture retrievals only where the reported optical 
! depth values are between 0 and 0.8
!--------------------------------------------------------------------------
        if(optc(c,r)*0.01.gt.0.8.or.optc(c,r).lt.0) smc(c,r) = -1
        if(optx(c,r)*0.01.gt.0.8.or.optx(c,r).lt.0) smx(c,r) = -1
!--------------------------------------------------------------------------
! Reject retrievals when soil temperature is below freezing. 
!--------------------------------------------------------------------------
        if(tskin(c,r)*0.1.lt.273.15.or.tskin(c,r).lt.0) smc(c,r) = -1
        if(tskin(c,r)*0.1.lt.273.15.or.tskin(c,r).lt.0) smx(c,r) = -1
!--------------------------------------------------------------------------
! Reject retrievals when RFI is non zero
!--------------------------------------------------------------------------
        if(rfi(c,r).ne.0) smc(c,r) = -1
        if(rfi(c,r).ne.0.and.rfi(c,r).ne.9) smx(c,r) = -1
!--------------------------------------------------------------------------
! Reject retrievals when residual error is large (> 0.2)
!-------------------------------------------------------------------------- 
        if(smerrc(c,r)*0.01.gt.0.2.or.smerrc(c,r)*0.01.lt.0) smc(c,r) = -1
        if(smerrx(c,r)*0.01.gt.0.2.or.smerrx(c,r)*0.01.lt.0) smx(c,r) = -1
!--------------------------------------------------------------------------
! Reject retrievals when rain exists
!-------------------------------------------------------------------------- 
        if(rainf(c,r).ne.0) then 
           smc(c,r) = -1
           smx(c,r) = -1
        endif
!--------------------------------------------------------------------------
! Reject retrievals if flag is not set
! -99: Edge of the swath is masked (AMSR-E observations are disturbed along the edges)
!-10: Default (open water (oceans), no brightness temperature observations over land and land surface temperature below freezing)
!-9: No land surface temperature retrieval possible (Thomas Holmes developed some filtering of suspicious areas for retrieving LST)
!-5: Vegetation is too dense for regions where we use C-band
!-4: Vegetation is too dense for regions where we had to switch to X-band due to RFI in C-band
!-3: Both observations in C- and X-band frequencies are suspicious (happens sometimes in India)
!6: Areas where we used C-band (6.9 GHz) observations
!10: Areas where we used X-band (10.7 GHz) observations
!-------------------------------------------------------------------------- 
        if(sm_flags(c,r).lt.0) then 
           smc(c,r) = -1
           smx(c,r) = -1
        endif

!--------------------------------------------------------------------------
! Apply the QC flags
!-------------------------------------------------------------------------- 
        if(smc(c,r).gt.0.and.smx(c,r).gt.0) then 
!           if(LPRM_AMSREsm_obs_struc(n)%rawdata.eq.1) then 
!              sm_combined(c,LPRM_AMSREsm_obs_struc(n)%lprmnr-r+1)  = &
!                   sm_combined(c,r)/100.0
!           else
              if(sm_cdf(c,r).gt.0.and.sm_cdf(c,r).lt.50) then 
                 sm_combined(c,LPRM_AMSREsm_obs_struc(n)%lprmnr-r+1) =&
                      sm_cdf(c,r)/100.0
              else
                 sm_combined(c,LPRM_AMSREsm_obs_struc(n)%lprmnr-r+1) = LIS_rc%udef
              endif
!           endif
        else
           sm_combined(c,LPRM_AMSREsm_obs_struc(n)%lprmnr-r+1) = LIS_rc%udef
        endif
     enddo
  enddo
    
  do r=1, LPRM_AMSREsm_obs_struc(n)%lprmnr
     do c=1, LPRM_AMSREsm_obs_struc(n)%lprmnc
        sm_data(c+(r-1)*LPRM_AMSREsm_obs_struc(n)%lprmnc) = sm_combined(c,r)
!        sm_data(c+(r-1)*LPRM_AMSREsm_obs_struc(n)%lprmnc) = sm_cdf(c,r)
        if(sm_combined(c,r).ne.LIS_rc%udef) then 
           sm_data_b(c+(r-1)*LPRM_AMSREsm_obs_struc(n)%lprmnc) = .true. 
        else
           sm_data_b(c+(r-1)*LPRM_AMSREsm_obs_struc(n)%lprmnc) = .false.
        endif
     enddo
  enddo

!!!!!  open(100,file='sm.bin',form='unformatted')
!!!!!  write(100) sm_data
!!!!!  close(100)
!!!!!  stop

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call bilinear_interp(LIS_rc%gridDesc(n,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       LPRM_AMSREsm_obs_struc(n)%lprmnc*LPRM_AMSREsm_obs_struc(n)%lprmnr, &
       LIS_rc%lnc(n)*LIS_rc%lnr(n), &
       LIS_domain(n)%lat,LIS_domain(n)%lon,&
       LPRM_AMSREsm_obs_struc(n)%w11, LPRM_AMSREsm_obs_struc(n)%w12, &
       LPRM_AMSREsm_obs_struc(n)%w21, LPRM_AMSREsm_obs_struc(n)%w22, &
       LPRM_AMSREsm_obs_struc(n)%n11, LPRM_AMSREsm_obs_struc(n)%n12, &
       LPRM_AMSREsm_obs_struc(n)%n21, LPRM_AMSREsm_obs_struc(n)%n22, &
       LIS_rc%udef, ios)


#endif
  
end subroutine read_LPRM_obs_data


!BOP
! !ROUTINE: create_LPRM_AMSREsm_obs_filename
! \label{create_LPRM_AMSREsm_obs_filename}
! 
! !INTERFACE: 
subroutine create_LPRM_AMSREsm_obs_filename(ndir, path, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: path
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the LPRM AMSRE cmd based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the LPRM AMSRE soil moisture directory
!  \item[path] name of the sensor path (A-ascending, D-descending)
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated LPRM filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/'//trim(fyr)//'-'//trim(fmo)//'/AMSR_L3_LPRMv05_' &
       //trim(path)//'_' &
       //trim(fyr)//trim(fmo)//trim(fda)//'T000000_' &
       //trim(fyr)//trim(fmo)//trim(fda)//'T235959D.v05.nc'

  
end subroutine create_LPRM_AMSREsm_obs_filename

