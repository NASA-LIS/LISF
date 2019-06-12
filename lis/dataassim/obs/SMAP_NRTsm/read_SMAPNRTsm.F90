!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_NASASMAPNRTsm
! \label{read_NASASMAPNRTsm}
!
! !REVISION HISTORY:
!  17 Jun 2010: Sujay Kumar; Updated for use with LPRM AMSRE Version 5. 
!  20 Sep 2012: Sujay Kumar; Updated to the NETCDF version of the data. 
!  1  Apr 2019: Yonghwan Kwon: Upated for reading monthy CDF for the current month
!
! !INTERFACE: 
subroutine read_SMAPNRTsm(n, k, OBS_State, OBS_Pert_State)
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
  use SMAPNRTsm_Mod, only : SMAPNRTsm_struc

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
  real                   :: sm_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: dt
  real                   :: lon
  real                   :: lhour
  integer                :: zone
  integer                :: fnd
  real, allocatable      :: ssdev(:)
  integer                :: lis_julss
  real                   :: smvalue
  real                   :: model_delta(LIS_rc%obs_ngrid(k))
  real                   :: obs_delta(LIS_rc%obs_ngrid(k))

  integer                :: mn_ind,orb_ind
  character*8            :: yyyymmdd
  character*4            :: yyyy
  character*2            :: mm,dd,hh
  integer                :: cyr, cmo, cda, chr, cmn, css
  integer                :: nyr, nmo, nda, nhr, nmn, nss
  integer                :: yr, mo, da, hr, mn, ss
  integer                :: doy
  integer                :: orbid(10),runid(10)
  real                   :: gmt
  real*8                 :: timenow, time1,time2,time3
  character*200          :: list_files
  character*100          :: temp1
  character*1            :: fproc(4)
  integer                :: ftn,ierr
  character*100          :: smap_filename(10),tstring(10)
  character(len=4) :: istring
  character(len=200) :: cmd
  integer :: rc

  integer, external :: create_filelist ! C function

  smap_filename = ""
  tstring       = ""
  orbid         = 0
  runid         = 0 

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 
  obs_unsc = LIS_rc%udef

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SMAP NRT read alarm")

  cyr = LIS_rc%yr
  cmo = LIS_rc%mo
  cda = LIS_rc%da
  chr = LIS_rc%hr
  cmn = LIS_rc%mn
  css = LIS_rc%ss

  call LIS_tick(time1,doy,gmt,cyr,cmo,cda,chr,cmn,css,0.0)

  nyr = LIS_rc%yr
  nmo = LIS_rc%mo
  nda = LIS_rc%da
  nhr = LIS_rc%hr
  nmn = LIS_rc%mn
  nss = LIS_rc%ss

  call LIS_tick(time2,doy,gmt,nyr,nmo,nda,nhr,nmn,nss,3600.0)
  
  nyr = LIS_rc%yr
  nmo = LIS_rc%mo
  nda = LIS_rc%da
  nhr = LIS_rc%hr
  nmn = LIS_rc%mn
  nss = LIS_rc%ss

  call LIS_tick(time3,doy,gmt,nyr,nmo,nda,nhr,nmn,nss,LIS_rc%ts)

  if(alarmCheck.or.SMAPNRTsm_struc(n)%startMode) then 
     SMAPNRTsm_struc(n)%startMode = .false.
     smobs = LIS_rc%udef
     SMAPNRTsm_struc(n)%smobs = LIS_rc%udef
     SMAPNRTsm_struc(n)%smtime = -1.0

     write(unit=temp1,fmt='(i4.4)') LIS_localPet
     read(unit=temp1,fmt='(4a1)')fproc
     
!  Time stamp in the SMAP file: 
!   Date/time in Universal Coordinated Time (UTC) of the first 
!   data element that appears in the product, where:
!    yyyymmdd	4-digit year, 2-digit month, 2-digit day
!   hhmmss	2-digit hour, 2-digit month, 2-digit second
! 
     write(yyyymmdd,'(i4.4,2i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
     write(yyyy,'(i4.4)') LIS_rc%yr
     write(mm,'(i2.2)') LIS_rc%mo
     write(dd,'(i2.2)') LIS_rc%da
     write(hh,'(i2.2)') LIS_rc%hr
     
     ! EMK...Make sure only one PET calls the file system to determine what
     ! SMAP files are available.  Then create a copy of the file list for
     ! every PET.
     if (LIS_masterproc) then
        ! EMK...Make sure we only process L2 files.
!        list_files = 'ls '//trim(smobsdir)//&
!             '/'//trim(yyyy)//'/'//trim(mm)//'/'//trim(dd)//&
!             '/*'//trim(yyyymmdd)//'T'//trim(hh)&
!             //"*.h5 > SMAP_filelist"//&
!             ".dat"
        ! EMK...Assume SMAP files are in common directory.
!        list_files = 'ls '//trim(smobsdir)//&
!             '/'//trim(yyyy)//'/'//trim(mm)//'/'//trim(dd)//&
!             '/SMAP_L2_*'//trim(yyyymmdd)//'T'//trim(hh)&
!             //"*.h5 > SMAP_filelist"//&
!             ".dat"
        !list_files = 'ls '//trim(smobsdir)//&
        !     '/SMAP_L2_*'//trim(yyyymmdd)//'T'//trim(hh)&
        !     //"*.h5 > SMAP_filelist"//&
        !     ".dat"
        !!fproc(1)//fproc(2)//fproc(3)//fproc(4)//".dat"
        !call system(trim(list_files))

        ! EMK...Avoid use of the 'ls' system command.
        list_files = trim(smobsdir)//&
             '/SMAP_L2_*'//trim(yyyymmdd)//'T'//trim(hh)&
             //"*.h5"
        write(LIS_logunit,*) &
             '[INFO] Searching for ',trim(list_files)

        rc = create_filelist(trim(list_files)//char(0), &
             "SMAP_filelist.dat"//char(0))
        if (rc .ne. 0) then
           write(LIS_logunit,*) &
                '[WARN] Problem encountered when searching for SMAP files'
           write(LIS_logunit,*) &
                'Was searching for ',trim(list_files)
           write(LIS_logunit,*) &
                'LIS will continue...'
        endif

        ! EMK TODO...Replace this with C code that uses the POSIX link
        ! function.
        do i = 0, LIS_npes-1
           write(istring,'(I4.4)') i
           cmd = 'cp SMAP_filelist.dat SMAP_filelist.'//istring//".dat"
           call system(trim(cmd))
        end do ! i
     end if
#if (defined SPMD)
     call mpi_barrier(lis_mpi_comm,ierr)
#endif

     i =1
     ftn = LIS_getNextUnitNumber()
     open(ftn,file="./SMAP_filelist."//&
          fproc(1)//fproc(2)//fproc(3)//fproc(4)//".dat",&
          status='old',iostat=ierr)

! if multiple files for the same time and orbits are present, the latest
! one will overwrite older ones, though multiple (redundant) reads occur. 
! This assumes that the 'ls command' will list the files in that order. 
 
     do while(ierr.eq.0) 
        read(ftn,'(a)',iostat=ierr) fname
        if(ierr.ne.0) then 
           exit
        endif
        ! From the filename, parse out minute, second
        mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh))
!        tstring(i) = fname(mn_ind:mn_ind+14)

!        read(fname(mn_ind+23:mn_ind+25),'(i3.3)') runid(i)

        mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh))+11        
        read(fname(mn_ind:mn_ind+1),'(i2.2)') mn
        ss=0
        call LIS_tick(timenow,doy,gmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
             LIS_rc%hr, mn, ss, 0.0)
        
!        orb_ind = index(fname,'SMAP_L2_SM_P_NRT_')+17
!        read(fname(orb_ind:orb_ind+4),'(i5.5)') orbid(i)

        smap_filename(i) = fname

        write(LIS_logunit,*) '[INFO] reading ',trim(smap_filename(i))
        call read_SMAPNRT_data(n,k,smap_filename(i),smobs,timenow)

        i = i +1 
     enddo
     
     call LIS_releaseUnitNumber(ftn)

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           SMAPNRTsm_struc(n)%smobs(c,r) = & 
                smobs(c+(r-1)*LIS_rc%obs_lnc(k))
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

           dt = (SMAPNRTsm_struc(n)%smtime(c,r)-time1)
           if(dt.ge.0.and.dt.lt.(time3-time1)) then 
              sm_current(c,r) = & 
                   SMAPNRTsm_struc(n)%smobs(c,r)
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
  if (SMAPNRTsm_struc(n)%ntimes.gt.1.and.SMAPNRTsm_struc(n)%cdf_read_opt.eq.1) then
     if (.not. SMAPNRTsm_struc(n)%cdf_read_mon.or.LIS_rc%da.eq.1.and.LIS_rc%hr.eq.0.and.LIS_rc%mn.eq.0.and.LIS_rc%ss.eq.0) then
        call LIS_readMeanSigmaData(n,k,&
             SMAPNRTsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             SMAPNRTsm_struc(n)%modelcdffile, &
             "SoilMoist",&
             SMAPNRTsm_struc(n)%model_mu,&
             SMAPNRTsm_struc(n)%model_sigma,&
             LIS_rc%mo)

        call LIS_readMeanSigmaData(n,k,&
             SMAPNRTsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             SMAPNRTsm_struc(n)%obscdffile, &
             "SoilMoist",&
             SMAPNRTsm_struc(n)%obs_mu,&
             SMAPNRTsm_struc(n)%obs_sigma,&
             LIS_rc%mo)

        call LIS_readCDFdata(n,k,&
             SMAPNRTsm_struc(n)%nbins,&
             SMAPNRTsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             SMAPNRTsm_struc(n)%modelcdffile, &
             "SoilMoist",&
             SMAPNRTsm_struc(n)%model_xrange,&
             SMAPNRTsm_struc(n)%model_cdf,&
             LIS_rc%mo)

        call LIS_readCDFdata(n,k,&
             SMAPNRTsm_struc(n)%nbins,&
             SMAPNRTsm_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             SMAPNRTsm_struc(n)%obscdffile, &
             "SoilMoist",&
             SMAPNRTsm_struc(n)%obs_xrange,&
             SMAPNRTsm_struc(n)%obs_cdf,&
             LIS_rc%mo)

        SMAPNRTsm_struc(n)%cdf_read_mon = .true.
     endif
  endif

  if(LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then
     if (SMAPNRTsm_struc(n)%ntimes.gt.1.and.SMAPNRTsm_struc(n)%cdf_read_opt.eq.1) then     
        call LIS_rescale_with_CDF_matching(     &
             n,k,                               &
             SMAPNRTsm_struc(n)%nbins,         &
             1,        &
             MAX_SM_VALUE,                      &
             MIN_SM_VALUE,                      &
             SMAPNRTsm_struc(n)%model_xrange,  &
             SMAPNRTsm_struc(n)%obs_xrange,    &
             SMAPNRTsm_struc(n)%model_cdf,     &
             SMAPNRTsm_struc(n)%obs_cdf,       &
             sm_current)
     else
        call LIS_rescale_with_CDF_matching(     &
             n,k,                               & 
             SMAPNRTsm_struc(n)%nbins,         & 
             SMAPNRTsm_struc(n)%ntimes,        & 
             MAX_SM_VALUE,                      & 
             MIN_SM_VALUE,                      & 
             SMAPNRTsm_struc(n)%model_xrange,  &
             SMAPNRTsm_struc(n)%obs_xrange,    &
             SMAPNRTsm_struc(n)%model_cdf,     &
             SMAPNRTsm_struc(n)%obs_cdf,       &
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
       //trim(LIS_SMAPNRTsmobsId)//char(0),n,k,OBS_state)

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
        if(SMAPNRTsm_struc(n)%useSsdevScal.eq.1) then
           call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                rc=status)
           call LIS_verify(status, 'Error: StateGet Observation01')
           
           allocate(ssdev(LIS_rc%obs_ngrid(k)))
           ssdev = SMAPNRTsm_struc(n)%ssdev_inp 
           
           if (SMAPNRTsm_struc(n)%cdf_read_opt .eq. 1) then
              jj = 1
           else if(SMAPNRTsm_struc(n)%ntimes.eq.1) then 
              jj = 1
           else
              jj = LIS_rc%mo
           endif
           do t=1,LIS_rc%obs_ngrid(k)
              if(SMAPNRTsm_struc(n)%obs_sigma(t,jj).gt.0) then 
                 ssdev(t) = ssdev(t)*SMAPNRTsm_struc(n)%model_sigma(t,jj)/&
                      SMAPNRTsm_struc(n)%obs_sigma(t,jj)
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

end subroutine read_SMAPNRTsm

!BOP
! 
! !ROUTINE: read_SMAPNRT_data
! \label{read_SMAPNRT_data}
!
! !INTERFACE:
subroutine read_SMAPNRT_data(n, k, fname, smobs_inp, time)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use SMAPNRTsm_Mod, only : SMAPNRTsm_struc
#if (defined USE_HDF5) 
  use hdf5
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n 
  integer                  :: k
  character (len=*)        :: fname
  real                     :: smobs_inp(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real*8                   :: time

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

  integer(hsize_t), dimension(1) :: dims
  integer(hsize_t), dimension(1) :: maxdims
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: dspace_id
  integer(hid_t)                 :: row_id, col_id
  integer(hid_t)                 :: sm_gr_id,sm_field_id, sm_qa_id
  integer(hid_t)                 :: sm_gr_id_A,sm_field_id_A
  real,             allocatable  :: sm_field(:)
  integer,          allocatable  :: sm_qa(:)
  integer,          allocatable  :: ease_row(:)
  integer,          allocatable  :: ease_col(:)
  integer                        :: c,r,t
  logical*1                      :: sm_data_b(SMAPNRTsm_struc(n)%nc*SMAPNRTsm_struc(n)%nr)
  logical*1                      :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: sm_data(SMAPNRTsm_struc(n)%nc*SMAPNRTsm_struc(n)%nr)
  real                           :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  integer                        :: status,ios

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LIS_verify(status, 'Error opening SMAP NRT file ')
  
  call h5gopen_f(file_id,sm_gr_name,sm_gr_id, status)
  call LIS_verify(status, 'Error opening SM group in SMAP NRT file')
  
  call h5dopen_f(sm_gr_id,sm_field_name,sm_field_id, status)
  call LIS_verify(status, 'Error opening SM field in SMAP NRT file')

  call h5dopen_f(sm_gr_id,"EASE_row_index",row_id, status)
  call LIS_verify(status, 'Error opening row index field in SMAP NRT file')

  call h5dopen_f(sm_gr_id,"EASE_column_index",col_id, status)
  call LIS_verify(status, 'Error opening column index field in SMAP NRT file')

  call h5dopen_f(sm_gr_id, sm_qa_name,sm_qa_id, status)
  call LIS_verify(status, 'Error opening QA field in SMAP NRT file')
  
  call h5dget_space_f(sm_field_id, dspace_id, status)
  call LIS_verify(status, 'Error in h5dget_space_f: readSMAP NRTObs')
  
! Size of the arrays
! This routine returns -1 on failure, rank on success. 
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, status) 
  if(status.eq.-1) then 
     call LIS_verify(status, 'Error in h5sget_simple_extent_dims_f: readSMAP NRTObs')
  endif
  
  allocate(sm_field(maxdims(1)))
  allocate(sm_qa(maxdims(1)))
  allocate(ease_row(maxdims(1)))
  allocate(ease_col(maxdims(1)))

  call h5dread_f(row_id, H5T_NATIVE_INTEGER,ease_row,dims,status)
  call LIS_verify(status, 'Error extracting row index from SMAP NRTfile')

  call h5dread_f(col_id, H5T_NATIVE_INTEGER,ease_col,dims,status)
  call LIS_verify(status, 'Error extracting col index from SMAP NRTfile')
  
  call h5dread_f(sm_field_id, H5T_NATIVE_REAL,sm_field,dims,status)
  call LIS_verify(status, 'Error extracting SM field from SMAP NRTfile')

  call h5dread_f(sm_qa_id, H5T_NATIVE_INTEGER,sm_qa,dims,status)
  call LIS_verify(status, 'Error extracting SM field from SMAP NRTfile')
  
  call h5dclose_f(sm_qa_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(row_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(col_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(sm_field_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')
  
  call h5gclose_f(sm_gr_id,status)
  call LIS_verify(status,'Error in H5GCLOSE call')
    
  call h5fclose_f(file_id,status)
  call LIS_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LIS_verify(status,'Error in H5CLOSE call')

  sm_data = LIS_rc%udef
  sm_data_b = .false. 

!grid the data in EASE projection
  do t=1,maxdims(1)
     if(ibits(sm_qa(t),0,1).eq.0) then 
        sm_data(ease_col(t) + (ease_row(t)-1)*SMAPNRTsm_struc(n)%nc) = sm_field(t) 
        if(sm_field(t).ne.-9999.0) then 
           sm_data_b(ease_col(t) + (ease_row(t)-1)*SMAPNRTsm_struc(n)%nc) = .true. 
        endif
     endif
  enddo

!  open(100,file='smobs.bin',form='unformatted')
!  write(100) sm_data
!  close(100)

  t = 1

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       SMAPNRTsm_struc(n)%nc*SMAPNRTsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       SMAPNRTsm_struc(n)%rlat, SMAPNRTsm_struc(n)%rlon,&
       SMAPNRTsm_struc(n)%n11, LIS_rc%udef, ios)

  deallocate(sm_field)
  deallocate(sm_qa)
  deallocate(ease_row)
  deallocate(ease_col)

!overwrite the input data 
  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(smobs_ip(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then 
           smobs_inp(c+(r-1)*LIS_rc%obs_lnc(k)) = &
                smobs_ip(c+(r-1)*LIS_rc%obs_lnc(k))

           SMAPNRTsm_struc(n)%smtime(c,r) = & 
                time
        endif
     enddo
  enddo

#endif

end subroutine read_SMAPNRT_data




