!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_NASASMAPsm
! \label{read_NASASMAPsm}
!
! !REVISION HISTORY:
!  17 Jun 2010: Sujay Kumar; Updated for use with LPRM AMSRE Version 5. 
!  20 Sep 2012: Sujay Kumar; Updated to the NETCDF version of the data. 
!  15 Jul 2018: Mahdi Navari: Bug in the SMAP reader was fixed
!  31 Aug 2018: Mahdi Navari, Edited to read SPL3SMP.005 & SPL3SMP_E.002
!  1  Apr 2019: Yonghwan Kwon: Upated for reading monthy CDF for the current month
!
! !INTERFACE: 
subroutine read_NASASMAPsm(n, k, OBS_State, OBS_Pert_State)
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
  use NASASMAPsm_Mod, only : NASASMAPsm_struc

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
  real                   :: smobs_D(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: smobs_A(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: sm_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: dt
  real                   :: lon
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
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

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "NASASMAP read alarm")

  smobs_A = LIS_rc%udef
  smobs_D = LIS_rc%udef

  if(alarmCheck.or.NASASMAPsm_struc(n)%startMode) then 
     NASASMAPsm_struc(n)%startMode = .false.
!MN: bug fix: "SPL3SMP" reader was updated to read the new version of data 
!       if(NASASMAPsm_struc(n)%data_designation.eq."SPL3SMP_E") then 
     if ( (NASASMAPsm_struc(n)%data_designation.eq."SPL3SMP_E") .or. &
          (NASASMAPsm_struc(n)%data_designation.eq."SPL3SMP") ) then
        call create_NASASMAPsm_filename(smobsdir, &
             NASASMAPsm_struc(n)%data_designation,&
             LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname)
        
        inquire(file=fname,exist=file_exists)
        
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
           call read_NASASMAP_E_data(n,k,'D',fname,smobs_D)
        else
           write(LIS_logunit,*) '[WARN] Missing SMAP file: ',trim(fname)
        endif

        call create_NASASMAPsm_filename(smobsdir, &
             NASASMAPsm_struc(n)%data_designation,&
             LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname)
        
        inquire(file=fname,exist=file_exists)
        
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
           call read_NASASMAP_E_data(n,k,'A',fname,smobs_A)
        else
           write(LIS_logunit,*) '[WARN] Missing SMAP file: ',trim(fname)
        endif

        NASASMAPsm_struc(n)%smobs  = LIS_rc%udef
        NASASMAPsm_struc(n)%smtime = -1

!------------------------------------------------------------------------- 
!   Ascending pass assumed to be at 6pm localtime and the descending 
!   pass is assumed to be at 6am local time
!-------------------------------------------------------------------------
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              grid_index = LIS_obs_domain(n,k)%gindex(c,r)
              if(grid_index.ne.-1) then 
                 
                 if(smobs_D(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then   
                    NASASMAPsm_struc(n)%smobs(c,r) = &
                         smobs_D(c+(r-1)*LIS_rc%obs_lnc(k))                 
                    lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                    lhour = 6.0
                    call LIS_localtime2gmt (gmt,lon,lhour,zone)
                    NASASMAPsm_struc(n)%smtime(c,r) = gmt

                 endif
!-------------------------------------------------------------------------  
! The ascending data is used only over locations where descending data
! doesn't exist. 
!-------------------------------------------------------------------------
                 if(smobs_A(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0.and.&
                      NASASMAPsm_struc(n)%smobs(c,r).eq.-9999.0) then   
                    NASASMAPsm_struc(n)%smobs(c,r) = &
                         smobs_A(c+(r-1)*LIS_rc%obs_lnc(k))                 
                    lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                    lhour = 18.0
                    call LIS_localtime2gmt (gmt,lon,lhour,zone)
                    NASASMAPsm_struc(n)%smtime(c,r) = gmt
                 endif
              endif
           enddo
        enddo
     else
        NASASMAPsm_struc(n)%smobs = LIS_rc%udef
        smobs = LIS_rc%udef

        call create_NASASMAPsm_filename(smobsdir, &
             NASASMAPsm_struc(n)%data_designation,&
             LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname)
        
        inquire(file=fname,exist=file_exists)
        
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
           call read_NASASMAP_data(n,k,fname,smobs)
        else
           write(LIS_logunit,*) '[WARN] Missing SMAP file: ',trim(fname)
        endif
        
        NASASMAPsm_struc(n)%smobs  = LIS_rc%udef
        NASASMAPsm_struc(n)%smtime = -1

!-------------------------------------------------------------------------  
!  From the SMAP documentation: 
!  The current approach for the SPL3SMP product is to use the nearest 
!  6:00 a.m. LST criterion to perform Level-3 compositing for the 
!  descending passes. According to this criterion, for a given grid cell, 
!  an L2 data point acquired closest to 6:00 a.m. local solar time will 
!  make its way to the final Level-3 file; other late-coming L2 data 
!  points falling into the same grid cell will be ignored. For a given 
!  file whose time stamp (yyyy-mm-ddThh:mm:ss) is expressed in UTC, only 
!  the hh:mm:ss part is converted into local solar time. 
!  (O'Neill et al. 2012)
!-------------------------------------------------------------------------
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              grid_index = LIS_obs_domain(n,k)%gindex(c,r)
              if(grid_index.ne.-1) then 
                 
                 if(smobs(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then   
                    NASASMAPsm_struc(n)%smobs(c,r) = &
                         smobs(c+(r-1)*LIS_rc%obs_lnc(k))                 
                    lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                    lhour = 6.0 
                    call LIS_localtime2gmt (gmt,lon,lhour,zone)
                    NASASMAPsm_struc(n)%smtime(c,r) = gmt
                 endif
              endif
           enddo
        enddo

     endif

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

           dt = (LIS_rc%gmt - NASASMAPsm_struc(n)%smtime(c,r))*3600.0
           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
              sm_current(c,r) = & 
                   NASASMAPsm_struc(n)%smobs(c,r)
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

  ! Read monthly CDF (only for the current month)
  if (NASASMAPsm_struc(n)%ntimes.gt.1.and.NASASMAPsm_struc(n)%cdf_read_opt.eq.1) then
     if (.not. NASASMAPsm_struc(n)%cdf_read_mon.or.LIS_rc%da.eq.1.and.LIS_rc%hr.eq.0.and.LIS_rc%mn.eq.0.and.LIS_rc%ss.eq.0) then
        call LIS_readMeanSigmaData(n,k,&
             NASASMAPsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             NASASMAPsm_struc(n)%modelcdffile, &
             "SoilMoist",&
             NASASMAPsm_struc(n)%model_mu,&
             NASASMAPsm_struc(n)%model_sigma,&
             LIS_rc%mo)

        call LIS_readMeanSigmaData(n,k,&
             NASASMAPsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             NASASMAPsm_struc(n)%obscdffile, &
             "SoilMoist",&
             NASASMAPsm_struc(n)%obs_mu,&
             NASASMAPsm_struc(n)%obs_sigma,&
             LIS_rc%mo)

        call LIS_readCDFdata(n,k,&
             NASASMAPsm_struc(n)%nbins,&
             NASASMAPsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             NASASMAPsm_struc(n)%modelcdffile, &
             "SoilMoist",&
             NASASMAPsm_struc(n)%model_xrange,&
             NASASMAPsm_struc(n)%model_cdf,&
             LIS_rc%mo)

        call LIS_readCDFdata(n,k,&
             NASASMAPsm_struc(n)%nbins,&
             NASASMAPsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             NASASMAPsm_struc(n)%obscdffile, &
             "SoilMoist",&
             NASASMAPsm_struc(n)%obs_xrange,&
             NASASMAPsm_struc(n)%obs_cdf,&
             LIS_rc%mo)

        NASASMAPsm_struc(n)%cdf_read_mon = .true.
     endif
  endif

  if(LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then
     if (NASASMAPsm_struc(n)%ntimes.gt.1.and.NASASMAPsm_struc(n)%cdf_read_opt.eq.1) then   
        call LIS_rescale_with_CDF_matching(     &
             n,k,                               &
             NASASMAPsm_struc(n)%nbins,         &
             1,                                 &
             MAX_SM_VALUE,                      &
             MIN_SM_VALUE,                      &
             NASASMAPsm_struc(n)%model_xrange,  &
             NASASMAPsm_struc(n)%obs_xrange,    &
             NASASMAPsm_struc(n)%model_cdf,     &
             NASASMAPsm_struc(n)%obs_cdf,       &
             sm_current)                             
     else
        call LIS_rescale_with_CDF_matching(     &
             n,k,                               & 
             NASASMAPsm_struc(n)%nbins,         & 
             NASASMAPsm_struc(n)%ntimes,        & 
             MAX_SM_VALUE,                      & 
             MIN_SM_VALUE,                      & 
             NASASMAPsm_struc(n)%model_xrange,  &
             NASASMAPsm_struc(n)%obs_xrange,    &
             NASASMAPsm_struc(n)%model_cdf,     &
             NASASMAPsm_struc(n)%obs_cdf,       &
             sm_current)                      
     endif  
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
       //trim(LIS_NASASMAPsmobsId)//char(0),n,k,OBS_state)

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

     if(LIS_rc%dascaloption(k).eq."CDF matching") then 
        if(NASASMAPsm_struc(n)%useSsdevScal.eq.1) then
           call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                rc=status)
           call LIS_verify(status, 'Error: StateGet Observation01')
           
           allocate(ssdev(LIS_rc%obs_ngrid(k)))
           ssdev = NASASMAPsm_struc(n)%ssdev_inp 
           
           if(NASASMAPsm_struc(n)%ntimes.eq.1) then 
              jj = 1
           else
              jj = LIS_rc%mo
           endif
           do t=1,LIS_rc%obs_ngrid(k)
              if(NASASMAPsm_struc(n)%obs_sigma(t,jj).gt.0) then 
                 ssdev(t) = ssdev(t)*NASASMAPsm_struc(n)%model_sigma(t,jj)/&
                      NASASMAPsm_struc(n)%obs_sigma(t,jj)
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
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif

end subroutine read_NASASMAPsm

!BOP
! 
! !ROUTINE: read_NASASMAP_E_data
! \label{read_NASASMAP_E_data}
!
! !INTERFACE:
subroutine read_NASASMAP_E_data(n, k, pass, fname, smobs_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use NASASMAPsm_Mod, only : NASASMAPsm_struc
#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: pass
  character (len=*)             :: fname
  real                          :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the SMOS NESDIS binary file and applies the data
!  quality flags to filter the data. !\normalsize
!
!  tb_time_seconds
!  Arithmetic average of the same parameters found in the 
!  fore- and aft-looking groups in the input SPL1CTB granule. 
!  The resulting parameter thus describes the average of UTC 
!  acquisition times of SPL1BTB observations whose boresights 
!  fall within a 36 km EASE-Grid 2.0 cell. The result is then 
!  expressed in J2000 seconds (the number of seconds since 
!  11:58:55.816 on January 1, 2000 UT).
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTNASASMAP AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDF5)

  character*100,    parameter    :: sm_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: sm_field_name_D = "soil_moisture"
  character*100,    parameter    :: sm_qa_name_D = "retrieval_qual_flag"
  character*100,    parameter    :: sm_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: sm_field_name_A = "soil_moisture_pm"
  character*100,    parameter    :: sm_qa_name_A = "retrieval_qual_flag_pm"
! MN 
  character*100,    parameter    :: vwc_field_name_D = "vegetation_water_content"
  character*100,    parameter    :: vwc_field_name_A = "vegetation_water_content_pm"

  integer(hsize_t), allocatable  :: dims(:)
  integer(hsize_t), dimension(2) :: dimsm
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hid_t)                 :: memspace
  integer(hid_t)                 :: dataspace
  integer                        :: memrank = 2 ! scaler--> rank = 0 ; 2D array--> rank = 2
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer(hsize_t), dimension(2) :: offset_file = (/0,0/)
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: sm_gr_id_D,sm_field_id_D,sm_qa_id_D
  integer(hid_t)                 :: sm_gr_id_A,sm_field_id_A,sm_qa_id_A
  integer(hid_t)                 :: vwc_field_id_D ! MN 
  integer(hid_t)                 :: vwc_field_id_A ! MN
  real,             allocatable  :: sm_field(:,:)
  real,             allocatable  :: vwc_field(:,:)! MN
  integer,          allocatable  :: sm_qa(:,:)
  integer                        :: c,r,t
  logical*1                      :: sm_data_b(NASASMAPsm_struc(n)%nc*NASASMAPsm_struc(n)%nr)
  logical*1                      :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: sm_data(NASASMAPsm_struc(n)%nc*NASASMAPsm_struc(n)%nr)
  integer                        :: status,ios

  dimsm      = (/NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr/)
  count_file = (/NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr/)
  count_mem  = (/NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr/)
  
  allocate(sm_field(NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr))
  allocate(sm_qa(NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr))
  allocate(vwc_field(NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr))
  allocate(dims(2))

  dims(1) = NASASMAPsm_struc(n)%nc
  dims(2) = NASASMAPsm_struc(n)%nr

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LIS_verify(status, 'Error opening NASASMAP file ')
  
  if(pass.eq.'D') then 
     call h5gopen_f(file_id,sm_gr_name_D,sm_gr_id_D, status)
     call LIS_verify(status, 'Error opening SM group in NASASMAP file')
     
     call h5dopen_f(sm_gr_id_D,sm_field_name_D,sm_field_id_D, status)
     call LIS_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id_D, dataspace, status)
     call LIS_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LIS_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LIS_verify(status, 'Error in h5create_simple_f; read_NASASMAPsm')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_NASASMAPsm')
     
     call h5dread_f(sm_field_id_D, H5T_NATIVE_REAL,sm_field,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM field from NASASMAPfile')

     call h5dclose_f(sm_field_id_D,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

     call h5dopen_f(sm_gr_id_D,sm_qa_name_D,sm_qa_id_D, status)
     call LIS_verify(status, 'Error opening SM QA field in NASASMAP file')
     
     call h5dread_f(sm_qa_id_D, H5T_NATIVE_INTEGER,sm_qa,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM QA field from NASASMAPfile')
     
     call h5dclose_f(sm_qa_id_D,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

! MN get the vegetation water contnent 
     call h5dopen_f(sm_gr_id_D,vwc_field_name_D,vwc_field_id_D, status)
     call LIS_verify(status, 'Error opening Veg water content field in NASASMAP file')

     call h5dread_f(vwc_field_id_D, H5T_NATIVE_REAL,vwc_field,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting Veg water content (AM) field from NASASMAPfile')

     call h5dclose_f(vwc_field_id_D,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(sm_gr_id_D,status)
     call LIS_verify(status,'Error in H5GCLOSE call')

  else
     call h5gopen_f(file_id,sm_gr_name_A,sm_gr_id_A, status)
     call LIS_verify(status, 'Error opening SM group in NASASMAP file')
     
     call h5dopen_f(sm_gr_id_A,sm_field_name_A,sm_field_id_A, status)
     call LIS_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(sm_field_id_A, dataspace, status)
     call LIS_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LIS_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LIS_verify(status, 'Error in h5create_simple_f; read_NASASMAPsm')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_NASASMAPsm')
     
     call h5dread_f(sm_field_id_A, H5T_NATIVE_REAL,sm_field,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM field from NASASMAPfile')

     call h5dclose_f(sm_field_id_A,status)
     call LIS_verify(status,'Error in H5DCLOSE call')
     
     call h5dopen_f(sm_gr_id_A,sm_qa_name_A,sm_qa_id_A, status)
     call LIS_verify(status, 'Error opening SM QA field in NASASMAP file')
     
     call h5dread_f(sm_qa_id_A, H5T_NATIVE_INTEGER,sm_qa,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM QA field from NASASMAPfile')

     call h5dclose_f(sm_qa_id_A,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

! MN get the vegetation water contnent 
     call h5dopen_f(sm_gr_id_A,vwc_field_name_A,vwc_field_id_A, status)
     call LIS_verify(status, 'Error opening Veg water content field in NASASMAP file')

     call h5dread_f(vwc_field_id_A, H5T_NATIVE_REAL,vwc_field,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting Veg water content (AM) field from NASASMAPfile')

     call h5dclose_f(vwc_field_id_A,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(sm_gr_id_A,status)
     call LIS_verify(status,'Error in H5GCLOSE call')
     
  endif
  
  call h5fclose_f(file_id,status)
  call LIS_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LIS_verify(status,'Error in H5CLOSE call')

#if 0 
! =============================================================================
! MN the following section added 
!  1- to filter the densly vegetaded areas. 
!  2- combine the AM and PM data 
     
  sm_field = LIS_rc%udef
  do r=1,NASASMAPsm_struc(n)%nr
     do c=1,NASASMAPsm_struc(n)%nc
        if(vwc_field_D(c,r).gt. 5 ) then !MN Aply QC : if VWC > 5 kg/m2 
           sm_field_D(c,r) = LIS_rc%udef
	 else 
           if(sm_field_D(c,r).ne.LIS_rc%udef) then 
              sm_field(c,r) = sm_field_D(c,r)
           endif
        endif
     enddo
  enddo

! MN : Her we replace the PM observation with the AM opservations if 
!      there are both PM and AM observations. I think this is 
!      not a correct way to combine the PM and AM observations  ??????
  do r=1,NASASMAPsm_struc(n)%nr
     do c=1,NASASMAPsm_struc(n)%nc   
        if(vwc_field_A(c,r).gt. 5 ) then !MN Aply QC : if VWC > 5 kg/m2 
           sm_field_A(c,r) = LIS_rc%udef
	 else      
           if(sm_field_A(c,r).ne.LIS_rc%udef) then 
              sm_field(c,r) = sm_field_A(c,r)
           endif
        enddo
  enddo
! MN end added section 
! =============================================================================
#endif  
  sm_data_b = .false. 
  t = 1

! The retrieval_quality_field variable's binary representation consists of bits
! that indicate whether retrieval is performed or not at a given grid cell. 
! When retrieval is performed, it contains additional bits to further 
! indicate the exit status and quality of the retrieval. The first bit 
! indicates the recommended quality (0-means retrieval has recommended quality
!

  do r=1,NASASMAPsm_struc(n)%nr
     do c=1,NASASMAPsm_struc(n)%nc        
        sm_data(t) = sm_field(c,r)

        if(vwc_field(c,r).gt. 5 ) then !MN Aply QC : if VWC > 5 kg/m2 
           sm_data(t) = LIS_rc%udef
	 else 

           if(sm_data(t).ne.-9999.0) then 
              if(NASASMAPsm_struc(n)%qcFlag.eq.1) then 
                 if(ibits(sm_qa(c,r),0,1).eq.0) then 
                    sm_data_b(t) = .true.
                 else
                    sm_data(t) = -9999.0
                 endif
              else
                 sm_data_b(t) = .true.
              endif
           endif
        endif

        t = t+1
     enddo
  enddo

  deallocate(sm_qa)

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       NASASMAPsm_struc(n)%nc*NASASMAPsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       NASASMAPsm_struc(n)%rlat, NASASMAPsm_struc(n)%rlon,&
       NASASMAPsm_struc(n)%n11, LIS_rc%udef, ios)

  deallocate(sm_field)
  deallocate(dims)

#endif

end subroutine read_NASASMAP_E_data


!BOP
! 
! !ROUTINE: read_NASASMAP_data
! \label{read_NASASMAP_data}
!
! !INTERFACE:
subroutine read_NASASMAP_data(n, k, fname, smobs_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use NASASMAPsm_Mod, only : NASASMAPsm_struc
#if (defined USE_HDF5) 
  use hdf5
#endif

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
!
! !DESCRIPTION: 
!  This subroutine reads the SMOS NESDIS binary file and applies the data
!  quality flags to filter the data. !\normalsize
!
!  tb_time_seconds
!  Arithmetic average of the same parameters found in the 
!  fore- and aft-looking groups in the input SPL1CTB granule. 
!  The resulting parameter thus describes the average of UTC 
!  acquisition times of SPL1BTB observations whose boresights 
!  fall within a 36 km EASE-Grid 2.0 cell. The result is then 
!  expressed in J2000 seconds (the number of seconds since 
!  11:58:55.816 on January 1, 2000 UT).
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTNASASMAP AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDF5)
  character*100,    parameter    :: sm_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,    parameter    :: sm_field_name = "soil_moisture"
  character*100,    parameter    :: sm_qa_name = "retrieval_qual_flag"

  character*100,    parameter    :: sm_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: sm_field_name_D = "soil_moisture"
  character*100,    parameter    :: sm_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: sm_field_name_A = "soil_moisture_pm"
! MN 
  character*100,    parameter    :: vwc_field_name_D = "vegetation_water_content"
  character*100,    parameter    :: vwc_field_name_A = "vegetation_water_content_pm"

  integer(hid_t)                 :: file_id, sm_gr_id,sm_field_id
  integer(hid_t)                 :: sm_gr_id_D,sm_field_id_D
  integer(hid_t)                 :: sm_gr_id_A,sm_field_id_A
  integer(hid_t)                 :: vwc_field_id_D ! MN 
  integer(hid_t)                 :: vwc_field_id_A ! MN
  integer(hid_t)                 :: dataspace
  integer(hid_t)                 :: memspace
  integer                        :: memrank = 2
  integer(hsize_t), allocatable  :: dims(:)
  integer(hsize_t), dimension(2) :: dimsm
  integer(hsize_t), dimension(2) :: offset_file
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  real,             allocatable  :: sm_field(:,:)
  real,             allocatable  :: sm_field_D(:,:)
  real,             allocatable  :: sm_field_A(:,:)
  real,             allocatable  :: vwc_field_D(:,:)! MN
  real,             allocatable  :: vwc_field_A(:,:)! MN
  integer,          allocatable  :: sm_qa(:,:)
  integer                        :: c,r,t
  logical*1                      :: sm_data_b(NASASMAPsm_struc(n)%nc*NASASMAPsm_struc(n)%nr)
  logical*1                      :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: sm_data(NASASMAPsm_struc(n)%nc*NASASMAPsm_struc(n)%nr)
  integer                        :: status,ios

  dimsm      = (/NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr/)
  count_file = (/NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr/)
  count_mem  = (/NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr/)
  
  allocate(sm_field(NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr))
  allocate(sm_field_D(NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr))
  allocate(sm_field_A(NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr))
!  allocate(sm_qa(NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr))
  allocate(vwc_field_A(NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr))
  allocate(vwc_field_D(NASASMAPsm_struc(n)%nc, NASASMAPsm_struc(n)%nr))
  allocate(dims(2))

  dims(1) = NASASMAPsm_struc(n)%nc
  dims(2) = NASASMAPsm_struc(n)%nr

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LIS_verify(status, 'Error opening NASASMAP file ')
  
  call h5gopen_f(file_id,sm_gr_name_D,sm_gr_id_D, status)
  call LIS_verify(status, 'Error opening SM group in NASASMAP file')
  
  call h5dopen_f(sm_gr_id_D,sm_field_name_D,sm_field_id_D, status)
  call LIS_verify(status, 'Error opening SM field in NASASMAP file')
  
  call h5dget_space_f(sm_field_id_D, dataspace, status)
  call LIS_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
       start=offset_file, count=count_file, hdferr=status)
  call LIS_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
  
  call h5screate_simple_f(memrank,dimsm, memspace, status)
  call LIS_verify(status, 'Error in h5create_simple_f; read_NASASMAPsm')
  
  call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
       start=offset_mem, count=count_mem, hdferr=status)
  call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_NASASMAPsm')
  
  call h5dread_f(sm_field_id_D, H5T_NATIVE_REAL,sm_field_D,dims,status, &
       memspace, dataspace)
  call LIS_verify(status, 'Error extracting SM field from NASASMAPfile')

! MN get the vegetation water contnent DES (6 AM)
  call h5dopen_f(sm_gr_id_D,vwc_field_name_D,vwc_field_id_D, status)
  call LIS_verify(status, 'Error opening Veg water content field in NASASMAP file')

  call h5dread_f(vwc_field_id_D, H5T_NATIVE_REAL,vwc_field_D,dims,status, &
          memspace, dataspace)
  call LIS_verify(status, 'Error extracting Veg water content (AM) field from NASASMAPfile')

  call h5dclose_f(vwc_field_id_D,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

 
!Read the PM (ascending) data 
  call h5gopen_f(file_id,sm_gr_name_A,sm_gr_id_A, status)
  call LIS_verify(status, 'Error opening SM group in NASASMAP file')
  
  call h5dopen_f(sm_gr_id_A,sm_field_name_A,sm_field_id_A, status)
  call LIS_verify(status, 'Error opening SM field in NASASMAP file')
  
  call h5dget_space_f(sm_field_id_A, dataspace, status)
  call LIS_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
       start=offset_file, count=count_file, hdferr=status)
  call LIS_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
  
  call h5screate_simple_f(memrank,dimsm, memspace, status)
  call LIS_verify(status, 'Error in h5create_simple_f; read_NASASMAPsm')
  
  call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
       start=offset_mem, count=count_mem, hdferr=status)
  call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_NASASMAPsm')
  
  call h5dread_f(sm_field_id_A, H5T_NATIVE_REAL,sm_field_A,dims,status, &
       memspace, dataspace)
  call LIS_verify(status, 'Error extracting SM field from NASASMAPfile')

! MN get the vegetation water contnent ASC (6 PM)
  call h5dopen_f(sm_gr_id_A,vwc_field_name_A,vwc_field_id_A, status)
  call LIS_verify(status, 'Error opening Veg water content field in NASASMAP file')

  call h5dread_f(vwc_field_id_A, H5T_NATIVE_REAL,vwc_field_A,dims,status, &
          memspace, dataspace)
  call LIS_verify(status, 'Error extracting Veg water content (AM) field from NASASMAPfile')

  call h5dclose_f(vwc_field_id_A,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(sm_field_id_D,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(sm_field_id_A,status)
  call LIS_verify(status,'Error in H5DCLOSE call')


! MN why this part has been commented out    ???????

!  call h5dopen_f(sm_gr_id,sm_qa_name,sm_qa_id, status)
!  call LIS_verify(status, 'Error opening SM QA field in NASASMAP file')
  
!  call h5dread_f(sm_qa_id, H5T_NATIVE_INTEGER,sm_qa,dims,status, &
!       memspace, dataspace)
!  call LIS_verify(status, 'Error extracting SM QA field from NASASMAPfile')
!  
!  call h5dclose_f(sm_qa_id,status)
!  call LIS_verify(status,'Error in H5DCLOSE call')
  

  call h5gclose_f(sm_gr_id_D,status)
  call LIS_verify(status,'Error in H5GCLOSE call')

  call h5gclose_f(sm_gr_id_A,status)
  call LIS_verify(status,'Error in H5GCLOSE call')
  
  call h5fclose_f(file_id,status)
  call LIS_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LIS_verify(status,'Error in H5CLOSE call')



  sm_field = LIS_rc%udef
  do r=1,NASASMAPsm_struc(n)%nr
     do c=1,NASASMAPsm_struc(n)%nc
        if(vwc_field_D(c,r).gt. 5 ) then !MN Aply QC : if VWC > 5 kg/m2 
           sm_field_D(c,r) = LIS_rc%udef
	 else 
           if(sm_field_D(c,r).ne.LIS_rc%udef) then 
              sm_field(c,r) = sm_field_D(c,r)
           endif
        endif
     enddo
  enddo

! MN : Her we replace the PM observation with the AM opservations if 
!      there are both PM and AM observations. I think this is 
!      not a correct way to combine the PM and AM observations  ??????
  do r=1,NASASMAPsm_struc(n)%nr
     do c=1,NASASMAPsm_struc(n)%nc   
        if(vwc_field_A(c,r).gt. 5 ) then !MN Aply QC : if VWC > 5 kg/m2 
           sm_field_A(c,r) = LIS_rc%udef
	 else      
           if(sm_field_A(c,r).ne.LIS_rc%udef) then 
              sm_field(c,r) = sm_field_A(c,r)
           endif
        endif
     enddo
  enddo

  
  sm_data_b = .false. 
  t = 1

! Why qc flag has been disabled  ??????

  do r=1,NASASMAPsm_struc(n)%nr
     do c=1,NASASMAPsm_struc(n)%nc        
        sm_data(t) = sm_field(c,r)
        if(sm_data(t).ne.-9999.0) then 
!           if(NASASMAPsm_struc(n)%qcFlag.eq.1) then 
!              if(ibits(sm_qa(c,r),0,1).eq.0) then 
           sm_data_b(t) = .true. 
!              endif
!           else
!              sm_data_b(t) = .true.
!           endif
        endif
        t = t+1
     enddo
  enddo

!  deallocate(sm_qa)
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       NASASMAPsm_struc(n)%nc*NASASMAPsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       NASASMAPsm_struc(n)%rlat, NASASMAPsm_struc(n)%rlon,&
       NASASMAPsm_struc(n)%n11, LIS_rc%udef, ios)

  deallocate(sm_field)
  deallocate(dims)

#endif

end subroutine read_NASASMAP_data


!BOP
! !ROUTINE: create_NASASMAPsm_filename
! \label{create_NASASMAPsm_filename}
! 
! !INTERFACE: 
subroutine create_NASASMAPsm_filename(ndir, designation, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: designation
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the SMAP filename (from NSIDC)
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the SMAP soil moisture data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated SMAP filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  if(designation.eq."SPL3SMAP") then 
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_AP_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '_R12170_001.h5'
! MN:   
! The SMAP file names contain components that change in a way that
! is difficult to programatically generate. So after downloading
! a SMAP data file, a symbolic link was created to it which make it 
! easier to generate file name.
!   For example:
!   SMAP_L3_SM_P_20170902.h5 -> SMAP_L3_SM_P_20170902_R15152_001.h5
  elseif(designation.eq."SPL3SMP") then 
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
         '.h5'
! For example:
! SMAP_L3_SM_P_E_20180811.h5 -> SMAP_L3_SM_P_E_20180811_R16010_001.h5
! MN for sake of convenience instead of real name the symbolic link was used
  elseif(designation.eq."SPL3SMP_E") then 
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_E_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '.h5'

  endif

end subroutine create_NASASMAPsm_filename





