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
! !ROUTINE: read_LPRM_AMSREsm
! \label{read_LPRM_AMSREsm}
!
! !REVISION HISTORY:
!  17 Jun 2010: Sujay Kumar; Updated for use with LPRM AMSRE Version 5. 
!  20 Sep 2012: Sujay Kumar; Updated to the NETCDF version of the data. 
!
! !INTERFACE: 
subroutine read_LPRM_AMSREsm(n, k, OBS_State, OBS_Pert_State)
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
  use LPRM_AMSREsm_Mod, only : LPRM_AMSREsm_struc

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
  !based on Liu et al. JHM 2011                   
  real, parameter        ::  minssdev = 0.01
  real, parameter        ::  maxssdev = 0.11
  real                   :: MAX_SM_VALUE, MIN_SM_VALUE
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: smobsdir, fname_A, fname_D
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
  real                   :: smobs_A(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: smobs_D(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: sm_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: dt
  real                   :: lon
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
  integer                :: fnd
  real                   :: smvalue
  real, allocatable      :: ssdev(:)
  real                   :: model_delta(LIS_rc%obs_ngrid(k))
  real                   :: obs_delta(LIS_rc%obs_ngrid(k))
  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  if(LPRM_AMSREsm_struc(n)%rawdata.eq.0) then 
     MAX_SM_VALUE=0.45
     MIN_SM_VALUE=0.0001
  else ! in degree of saturation
     MAX_SM_VALUE=1.0
     MIN_SM_VALUE=0.0001
  endif

  data_upd = .false. 
  obs_unsc = LIS_rc%udef
!-------------------------------------------------------------------------
!   Read both ascending and descending passes at 0Z and then store
!   the overpass time as 1.30AM for the descending pass and 1.30PM 
!   for the ascending pass. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "AMSR-E(LPRM) read alarm")
  
  if(alarmCheck.or.LPRM_AMSREsm_struc(n)%startMode) then 
     LPRM_AMSREsm_struc(n)%startMode = .false.

     LPRM_AMSREsm_struc(n)%smobs = LIS_rc%udef
     smobs_A = LIS_rc%udef
     smobs_D = LIS_rc%udef             
     LPRM_AMSREsm_struc(n)%smtime = -1

     call create_LPRM_AMSREsm_filename(smobsdir, 'A',&
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname_A)

     inquire(file=fname_A,exist=file_exists)

     if(file_exists) then 

        write(LIS_logunit,*) 'Reading ',trim(fname_A)
        call read_LPRM_data(n,k,fname_A,smobs_A)
     endif

     call create_LPRM_AMSREsm_filename(smobsdir, 'D',&
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname_D)

     inquire(file=fname_D,exist=file_exists)

     if(file_exists) then 

        write(LIS_logunit,*) 'Reading ',trim(fname_D)
        call read_LPRM_data(n,k,fname_D,smobs_D)

     endif

     LPRM_AMSREsm_struc(n)%smobs  = LIS_rc%udef
     LPRM_AMSREsm_struc(n)%smtime = -1


     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           grid_index = LIS_obs_domain(n,k)%gindex(c,r)
           if(grid_index.ne.-1) then 

              if(smobs_A(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then   
                 LPRM_AMSREsm_struc(n)%smobs(c,r) = &
                      smobs_A(c+(r-1)*LIS_rc%obs_lnc(k))                 
                 lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                 lhour = 13.5 
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 LPRM_AMSREsm_struc(n)%smtime(c,r) = gmt
              endif

              if(smobs_D(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then           
                 LPRM_AMSREsm_struc(n)%smobs(c,r) = &
                      smobs_D(c+(r-1)*LIS_rc%obs_lnc(k))

                 lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                 lhour = 1.5
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 LPRM_AMSREsm_struc(n)%smtime(c,r) = gmt
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


  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           grid_index = c+(r-1)*LIS_rc%obs_lnc(k)

           dt = (LIS_rc%gmt - LPRM_AMSREsm_struc(n)%smtime(c,r))*3600.0
           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
              sm_current(c,r) = & 
                   LPRM_AMSREsm_struc(n)%smobs(c,r)
              if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                 obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = &
                      sm_current(c,r)
              endif
              if(sm_current(c,r).ne.LIS_rc%udef) then 
                 fnd = 1
              endif
           endif
        endif
     enddo
  enddo

  !-------------------------------------------------------------------------
  !  Transform data to the LSM climatology using a CDF-scaling approach
  !-------------------------------------------------------------------------  
  if(LIS_rc%dascaloption(k).ne."none".and.fnd.ne.0) then
     call LIS_rescale_with_CDF_matching(    &
          n,k,                                 & 
          LPRM_AMSREsm_struc(n)%nbins,         & 
          LPRM_AMSREsm_struc(n)%ntimes,        & 
          MAX_SM_VALUE,                        & 
          MIN_SM_VALUE,                        & 
          LPRM_AMSREsm_struc(n)%model_xrange,  &
          LPRM_AMSREsm_struc(n)%obs_xrange,    &
          LPRM_AMSREsm_struc(n)%model_cdf,     &
          LPRM_AMSREsm_struc(n)%obs_cdf,       &
          sm_current)
  endif

  obsl = LIS_rc%udef 
  do r=1, LIS_rc%obs_lnr(k)
     do c=1, LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           obsl(LIS_obs_domain(n,k)%gindex(c,r))=sm_current(c,r)
        endif
     enddo
  enddo

  !-------------------------------------------------------------------------
  !  Apply LSM based QC and screening of observations
  !-------------------------------------------------------------------------  
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_LPRM_AMSREsmobsId)//char(0),n, k,OBS_state)

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

     if(LPRM_AMSREsm_struc(n)%useSsdevScal.eq.1.and.&
          LPRM_AMSREsm_struc(n)%ntimes.gt.1) then 

        call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')

        allocate(ssdev(LIS_rc%obs_ngrid(k)))
        ssdev = LPRM_AMSREsm_struc(n)%ssdev_inp
        if(LPRM_AMSREsm_struc(n)%ntimes.eq.1) then 
           jj = 1
        else
           jj = LIS_rc%mo
        endif

        do t=1,LIS_rc%obs_ngrid(k)
           if(LPRM_AMSREsm_struc(n)%obs_sigma(t,jj).gt.0) then 
              ssdev(t) = ssdev(t)*LPRM_AMSREsm_struc(n)%model_sigma(t,jj)/&
                   LPRM_AMSREsm_struc(n)%obs_sigma(t,jj)
              if(ssdev(t).gt.maxssdev) ssdev(t) = maxssdev
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
end subroutine read_LPRM_AMSREsm

!BOP
! 
! !ROUTINE: read_LPRM_data
! \label{read_LPRM_data}
!
! !INTERFACE:
subroutine read_LPRM_data(n, k, fname, smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,   only : LIS_verify
  use LPRM_AMSREsm_Mod, only : LPRM_AMSREsm_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


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
  real                        :: rfi(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: tskin(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: rainf(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: optc(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: optx(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: smerrc(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: smerrx(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: smc(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: smx(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: sm_combined(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: sm_raw(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: sm_flags(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)
  real                        :: sm_cdf(LPRM_AMSREsm_struc(n)%lprmnc,&
       LPRM_AMSREsm_struc(n)%lprmnr)


  real                        :: sm_data(LPRM_AMSREsm_struc(n)%lprmnc* & 
       LPRM_AMSREsm_struc(n)%lprmnr)
  logical*1                   :: sm_data_b(LPRM_AMSREsm_struc(n)%lprmnc* & 
       LPRM_AMSREsm_struc(n)%lprmnr)
  logical*1                   :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

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

  if(LPRM_AMSREsm_struc(n)%rawdata.eq.0) then 
     ios = nf90_inq_varid(nid, 'SM_CDF',smcdfid)
     call LIS_verify(ios, 'Error nf90_inq_varid: SM_CDF')
  endif

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

  if(LPRM_AMSREsm_struc(n)%rawdata.eq.0) then 
     ios = nf90_get_var(nid, smcdfid, sm_cdf)
     call LIS_verify(ios, 'Error nf90_get_var: SM_CDF')
  endif
  

  ios = nf90_close(ncid=nid)
  call LIS_verify(ios,'Error closing file '//trim(fname))

  smc = 1.0
  smx = 1.0
  sm_raw = -9999.0

  do r=1, LPRM_AMSREsm_struc(n)%lprmnr
     do c=1, LPRM_AMSREsm_struc(n)%lprmnc
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
           if(LPRM_AMSREsm_struc(n)%rawdata.eq.1) then 
              if(sm_combined(c,r).gt.0) then 
                 sm_raw(c,LPRM_AMSREsm_struc(n)%lprmnr-r+1)  = &
                      sm_combined(c,r)/100.0
              else
                 sm_raw(c,LPRM_AMSREsm_struc(n)%lprmnr-r+1)  = LIS_rc%udef
              endif
           else
              if(sm_cdf(c,r).gt.0) then 
                 sm_combined(c,LPRM_AMSREsm_struc(n)%lprmnr-r+1) =&
                      sm_cdf(c,r)/100.0
              else
                 sm_combined(c,LPRM_AMSREsm_struc(n)%lprmnr-r+1) = LIS_rc%udef
              endif
           endif
        else
           if(LPRM_AMSREsm_struc(n)%rawdata.eq.1) then 
              sm_raw(c,LPRM_AMSREsm_struc(n)%lprmnr-r+1) = LIS_rc%udef
           else
              sm_combined(c,LPRM_AMSREsm_struc(n)%lprmnr-r+1) = LIS_rc%udef
           endif
        endif
     enddo
  enddo
  
  if(LPRM_AMSREsm_struc(n)%rawdata.eq.1) then 
     sm_combined = sm_raw
  endif

  do r=1, LPRM_AMSREsm_struc(n)%lprmnr
     do c=1, LPRM_AMSREsm_struc(n)%lprmnc
        sm_data(c+(r-1)*LPRM_AMSREsm_struc(n)%lprmnc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LIS_rc%udef) then 
           sm_data_b(c+(r-1)*LPRM_AMSREsm_struc(n)%lprmnc) = .true. 
        else
           sm_data_b(c+(r-1)*LPRM_AMSREsm_struc(n)%lprmnc) = .false.
        endif
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the DA observation space
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       LPRM_AMSREsm_struc(n)%lprmnc*LPRM_AMSREsm_struc(n)%lprmnr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       LPRM_AMSREsm_struc(n)%rlat, LPRM_AMSREsm_struc(n)%rlon, &
       LPRM_AMSREsm_struc(n)%n11,  LIS_rc%udef, ios)

#endif
  
end subroutine read_LPRM_data


!BOP
! !ROUTINE: create_LPRM_AMSREsm_filename
! \label{create_LPRM_AMSREsm_filename}
! 
! !INTERFACE: 
subroutine create_LPRM_AMSREsm_filename(ndir, path, yr, mo,da, filename)
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

  
end subroutine create_LPRM_AMSREsm_filename




