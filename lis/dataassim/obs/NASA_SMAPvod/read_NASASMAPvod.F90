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
! !ROUTINE: read_NASASMAPvod
! \label{read_NASASMAPvod}
!
! !REVISION HISTORY:
!  28 Mar 2019    Sujay Kumar; initial specification
!
! !INTERFACE: 
subroutine read_NASASMAPvod(n, k, OBS_State, OBS_Pert_State)
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
  use NASASMAPvod_Mod, only : NASASMAPvod_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the SMAP vegetation optical depth observations 
!  and rescales the data to a reference LAI data.
!  This routine essentially converts the vegetation 
!  optical depth datasets into the LAI space.  
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
  real,  parameter       :: MAX_LAI_VALUE=10.0, MIN_LAI_VALUE=0.0001
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: vodobsdir
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
  real                   :: vodobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: vodobs_D(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: vodobs_A(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: lai_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: dt
  real                   :: lon
  real                   :: lhour
  integer                :: zone
  integer                :: fnd
  real, allocatable      :: ssdev(:)
  integer                :: lis_julss
  real                   :: smvalue
  character(len=LIS_CONST_PATH_LEN) :: list_files
  character(len=LIS_CONST_PATH_LEN) :: smap_filename(10)
  real*8                 :: timenow, time1,time2,time3
  integer                :: doy
  real                   :: gmt
  character(len=4)       :: istring
  character(len=200)     :: cmd
  integer                :: ftn
  integer                :: ierr
  character*100          :: temp1
  character*1            :: fproc(4)
  integer                :: mn_ind
  integer                :: yr, mo, da, hr, mn, ss
  integer                :: cyr, cmo, cda, chr, cmn, css
  integer                :: nyr, nmo, nda, nhr, nmn, nss
  character*8            :: yyyymmdd
  character*4            :: yyyy
  character*2            :: mm,dd,hh
  real                   :: model_delta(LIS_rc%obs_ngrid(k))
  real                   :: obs_delta(LIS_rc%obs_ngrid(k))
  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       vodobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 
  obs_unsc = LIS_rc%udef

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "NASASMAP VOD read alarm")

  vodobs_A = LIS_rc%udef
  vodobs_D = LIS_rc%udef

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

  if(alarmCheck.or.NASASMAPvod_struc(n)%startMode) then 
     NASASMAPvod_struc(n)%startMode = .false.
     if ( (NASASMAPvod_struc(n)%data_designation.eq."SPL2SMP_E") .or. &
          (NASASMAPvod_struc(n)%data_designation.eq."SPL2SMP") ) then
       
        NASASMAPvod_struc(n)%vodobs = LIS_rc%udef
        NASASMAPvod_struc(n)%vodtime = -1.0

        write(temp1,fmt='(i4.4)') LIS_localPet
        read(temp1,fmt='(4a1)') fproc
        write(yyyymmdd,'(i4.4,2i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
        write(yyyy,'(i4.4)') LIS_rc%yr
        write(mm,'(i2.2)') LIS_rc%mo
        write(dd,'(i2.2)') LIS_rc%da
        write(hh,'(i2.2)') LIS_rc%hr
        
        if(LIS_masterproc) then 
           list_files = "ls "//trim((vodobsdir))//&
                "/"//trim(yyyy)//"."//trim(mm)//"."//dd//&
                "/SMAP_L2_*"//trim(yyyymmdd)//"T"//trim(hh)&
                //"*.h5 > SMAP_filelist.vod.dat"
        
           call system(trim(list_files))
           do i=0,LIS_npes-1
              write(istring,'(I4.4)') i
              cmd = 'cp SMAP_filelist.vod.dat SMAP_filelist.vod.'//istring//".dat"
              call system(trim(cmd))
           end do ! i
        end if
#if (defined SPMD)
        call mpi_barrier(lis_mpi_comm,ierr)
#endif

        i = 1
        ftn = LIS_getNextUnitNumber()
        open(ftn,file="./SMAP_filelist.vod."//&
             fproc(1)//fproc(2)//fproc(3)//fproc(4)//".dat",&
             status='old',iostat=ierr)

        do while(ierr.eq.0) 
           read(ftn,'(a)',iostat=ierr) fname
           if(ierr.ne.0) then 
              exit
           endif
           mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh))
           
           mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh))+11        
           read(fname(mn_ind:mn_ind+1),'(i2.2)') mn
           ss=0
           call LIS_tick(timenow,doy,gmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
                LIS_rc%hr, mn, ss, 0.0)
           
           smap_filename(i) = fname
           
           write(LIS_logunit,*) '[INFO] reading ',trim(smap_filename(i))
           
           call read_SMAPL2vod_data(n,k,smap_filename(i),&
                NASASMAPvod_struc(n)%vodobs,timenow)
           
           i = i+1
        enddo
        call LIS_releaseUnitNumber(ftn)

     elseif ( (NASASMAPvod_struc(n)%data_designation.eq."SPL3SMP_E") .or. &
          (NASASMAPvod_struc(n)%data_designation.eq."SPL3SMP") ) then
        call create_NASASMAPvod_filename(vodobsdir, &
             NASASMAPvod_struc(n)%data_designation,&
             LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname)
        
        inquire(file=fname,exist=file_exists)
        
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
           call read_NASASMAP_E_voddata(n,k,'D',fname,vodobs_D)
        else
           write(LIS_logunit,*) '[WARN] Missing SMAP file: ',trim(fname)
        endif

        call create_NASASMAPvod_filename(vodobsdir, &
             NASASMAPvod_struc(n)%data_designation,&
             LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname)
        
        inquire(file=fname,exist=file_exists)
        
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
           call read_NASASMAP_E_voddata(n,k,'A',fname,vodobs_A)
        else
           write(LIS_logunit,*) '[WARN] Missing SMAP file: ',trim(fname)
        endif

        NASASMAPvod_struc(n)%vodobs  = LIS_rc%udef
        NASASMAPvod_struc(n)%vodtime = -1

!------------------------------------------------------------------------- 
!   Ascending pass assumed to be at 6pm localtime and the descending 
!   pass is assumed to be at 6am local time
!-------------------------------------------------------------------------
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              grid_index = LIS_obs_domain(n,k)%gindex(c,r)
              if(grid_index.ne.-1) then 
                 
                 if(vodobs_D(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then   
                    NASASMAPvod_struc(n)%vodobs(c,r) = &
                         vodobs_D(c+(r-1)*LIS_rc%obs_lnc(k))                 
                    lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                    lhour = 6.0
                    call LIS_localtime2gmt (gmt,lon,lhour,zone)
                    NASASMAPvod_struc(n)%vodtime(c,r) = gmt

                 endif
!-------------------------------------------------------------------------  
! The ascending data is used only over locations where descending data
! doesn't exist. 
!-------------------------------------------------------------------------
                 if(vodobs_A(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0.and.&
                      NASASMAPvod_struc(n)%vodobs(c,r).eq.-9999.0) then   
                    NASASMAPvod_struc(n)%vodobs(c,r) = &
                         vodobs_A(c+(r-1)*LIS_rc%obs_lnc(k))                 
                    lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                    lhour = 18.0
                    call LIS_localtime2gmt (gmt,lon,lhour,zone)
                    NASASMAPvod_struc(n)%vodtime(c,r) = gmt
                 endif
              endif
           enddo
        enddo
     else
        NASASMAPvod_struc(n)%vodobs = LIS_rc%udef
        vodobs = LIS_rc%udef

        call create_NASASMAPvod_filename(vodobsdir, &
             NASASMAPvod_struc(n)%data_designation,&
             LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname)
        
        inquire(file=fname,exist=file_exists)
        
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
           call read_NASASMAP_voddata(n,k,fname,vodobs)
        else
           write(LIS_logunit,*) '[WARN] Missing SMAP file: ',trim(fname)
        endif
        
        NASASMAPvod_struc(n)%vodobs  = LIS_rc%udef
        NASASMAPvod_struc(n)%vodtime = -1

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
                 
                 if(vodobs(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then   
                    NASASMAPvod_struc(n)%vodobs(c,r) = &
                         vodobs(c+(r-1)*LIS_rc%obs_lnc(k))                 
                    lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                    lhour = 6.0 
                    call LIS_localtime2gmt (gmt,lon,lhour,zone)
                    NASASMAPvod_struc(n)%vodtime(c,r) = gmt
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
  lai_current = LIS_rc%udef
 
  ! dt is not defined as absolute value of the time difference to avoid
  ! double counting of the data in assimilation. 

  if ( (NASASMAPvod_struc(n)%data_designation.eq."SPL2SMP_E") .or. &
       (NASASMAPvod_struc(n)%data_designation.eq."SPL2SMP") ) then

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
              dt = (NASASMAPvod_struc(n)%vodtime(c,r)-time1)
              if(dt.ge.0.and.dt.lt.(time3-time1)) then 
                 lai_current(c,r) = & 
                      NASASMAPvod_struc(n)%vodobs(c,r)
                 if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                    obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = &
                         lai_current(c,r)
                 endif
                 if(lai_current(c,r).ne.LIS_rc%udef) then 
                    fnd = 1
                 endif
              endif
           endif
        enddo
     enddo

  else
     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
              
              dt = (LIS_rc%gmt - NASASMAPvod_struc(n)%vodtime(c,r))*3600.0
              if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
                 lai_current(c,r) = & 
                      NASASMAPvod_struc(n)%vodobs(c,r)
                 if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                    obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = &
                         lai_current(c,r)
                 endif
                 if(lai_current(c,r).ne.LIS_rc%udef) then 
                    fnd = 1
                 endif
              endif
           endif
        enddo
     enddo
  endif

!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     

  if(LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then
     
     call LIS_rescale_with_CDF_matching(     &
          n,k,                               & 
          NASASMAPvod_struc(n)%nbins,         & 
          NASASMAPvod_struc(n)%ntimes,        & 
          MAX_LAI_VALUE,                      & 
          MIN_LAI_VALUE,                      & 
          NASASMAPvod_struc(n)%model_xrange,  &
          NASASMAPvod_struc(n)%obs_xrange,    &
          NASASMAPvod_struc(n)%model_cdf,     &
          NASASMAPvod_struc(n)%obs_cdf,       &
          lai_current)
  endif
  
  obsl = LIS_rc%udef 
  do r=1, LIS_rc%obs_lnr(k)
     do c=1, LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           obsl(LIS_obs_domain(n,k)%gindex(c,r))=lai_current(c,r)
        endif
     enddo
  enddo
  !-------------------------------------------------------------------------
  !  Apply LSM based QC and screening of observations
  !-------------------------------------------------------------------------  
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_NASASMAPvodobsId)//char(0),n,k,OBS_state)

  call LIS_checkForValidObs(n,k,obsl,fnd,lai_current)

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
        if(NASASMAPvod_struc(n)%useSsdevScal.eq.1) then
           call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                rc=status)
           call LIS_verify(status, 'Error: StateGet Observation01')
           
           allocate(ssdev(LIS_rc%obs_ngrid(k)))
           ssdev = NASASMAPvod_struc(n)%ssdev_inp 
           
           if(NASASMAPvod_struc(n)%ntimes.eq.1) then 
              jj = 1
           else
              jj = LIS_rc%mo
           endif
           do t=1,LIS_rc%obs_ngrid(k)
              if(NASASMAPvod_struc(n)%obs_sigma(t,jj).gt.0) then 
                 ssdev(t) = ssdev(t)*NASASMAPvod_struc(n)%model_sigma(t,jj)/&
                      NASASMAPvod_struc(n)%obs_sigma(t,jj)
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

end subroutine read_NASASMAPvod

!BOP
! 
! !ROUTINE: read_SMAPL2vod_data
! \label{read_SMAPL2vod_data}
!
! !INTERFACE:
subroutine read_SMAPL2vod_data(n, k,fname, vodobs_inp, time)
! 
! !USES:   

  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use NASASMAPvod_Mod, only : NASASMAPvod_struc

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
  real                     :: vodobs_inp(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real*8                   :: time

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!
!
!EOP

#if (defined USE_HDF5)

  character*100,    parameter    :: vod_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,    parameter    :: vod_field_name = "vegetation_opacity_option3"
  character*100,    parameter    :: vod_qa_name = "retrieval_qual_flag"

  integer(hsize_t), dimension(1) :: dims
  integer(hsize_t), dimension(1) :: maxdims
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: dspace_id
  integer(hid_t)                 :: row_id, col_id
  integer(hid_t)                 :: vod_gr_id,vod_field_id, vod_qa_id
  integer(hid_t)                 :: vod_gr_id_A,vod_field_id_A
  real,             allocatable  :: vod_field(:)
  integer,          allocatable  :: vod_qa(:)
  integer,          allocatable  :: ease_row(:)
  integer,          allocatable  :: ease_col(:)
  integer                        :: c,r,t
  logical*1                      :: vod_data_b(NASASMAPvod_struc(n)%nc*NASASMAPvod_struc(n)%nr)
  logical*1                      :: vodobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: vod_data(NASASMAPvod_struc(n)%nc*NASASMAPvod_struc(n)%nr)
  real                           :: vodobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  integer                        :: status,ios,iret

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LIS_verify(status, 'Error opening SMAP L2 file ')
  
  call h5gopen_f(file_id,vod_gr_name,vod_gr_id, status)
  call LIS_verify(status, 'Error opening VOD group in SMAP L2 file')
  
  call h5dopen_f(vod_gr_id,vod_field_name,vod_field_id, status)
  call LIS_verify(status, 'Error opening VOD field in SMAP L2 file')

  call h5dopen_f(vod_gr_id,"EASE_row_index",row_id, status)
  call LIS_verify(status, 'Error opening row index field in SMAP L2 file')

  call h5dopen_f(vod_gr_id,"EASE_column_index",col_id, status)
  call LIS_verify(status, 'Error opening column index field in SMAP L2 file')

!  call h5dopen_f(sm_gr_id, sm_qa_name,sm_qa_id, status)
!  call LIS_verify(status, 'Error opening QA field in SMAP L2 file')
  
  call h5dget_space_f(vod_field_id, dspace_id, status)
  call LIS_verify(status, 'Error in h5dget_space_f: reaSMAP L2Obs')
  
! Size of the arrays
! This routine returns -1 on failure, rank on success. 
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, status) 
  if(status.eq.-1) then 
     call LIS_verify(status, 'Error in h5sget_simple_extent_dims_f: readSMAP L2Obs')
  endif
  
  allocate(vod_field(maxdims(1)))
!  allocate(sm_qa(maxdims(1)))
  allocate(ease_row(maxdims(1)))
  allocate(ease_col(maxdims(1)))

  call h5dread_f(row_id, H5T_NATIVE_INTEGER,ease_row,dims,status)
  call LIS_verify(status, 'Error extracting row index from SMAP L2 file')

  call h5dread_f(col_id, H5T_NATIVE_INTEGER,ease_col,dims,status)
  call LIS_verify(status, 'Error extracting col index from SMAP L2 file')
  
  call h5dread_f(vod_field_id, H5T_NATIVE_REAL,vod_field,dims,status)
  call LIS_verify(status, 'Error extracting VOD field from SMAP L2 file')

!  call h5dread_f(vod_qa_id, H5T_NATIVE_INTEGER,vod_qa,dims,status)
!  call LIS_verify(status, 'Error extracting VOD field from SMAP L2 file')
  
!  call h5dclose_f(vod_qa_id,status)
!  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(row_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(col_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(vod_field_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')
  
  call h5gclose_f(vod_gr_id,status)
  call LIS_verify(status,'Error in H5GCLOSE call')
    
  call h5fclose_f(file_id,status)
  call LIS_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LIS_verify(status,'Error in H5CLOSE call')

  vod_data = LIS_rc%udef
  vod_data_b = .false. 

!grid the data in EASE projection
  do t=1,maxdims(1)
     if(ease_col(t).gt.0.and.ease_row(t).gt.0) then 
        vod_data(ease_col(t) + &
             (ease_row(t)-1)*NASASMAPvod_struc(n)%nc) = vod_field(t) 
        if(vod_field(t).ne.-9999.0) then 
           vod_data_b(ease_col(t) + &
                (ease_row(t)-1)*NASASMAPvod_struc(n)%nc) = .true. 
        endif
     endif
  enddo
  
  t = 1
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:), vod_data_b, vod_data, &
       vodobs_b_ip, vodobs_ip, &
       NASASMAPvod_struc(n)%nc*NASASMAPvod_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       NASASMAPvod_struc(n)%rlat, NASASMAPvod_struc(n)%rlon,&
       NASASMAPvod_struc(n)%n11, LIS_rc%udef, ios)


  deallocate(vod_field)
!  deallocate(vod_qa)
  deallocate(ease_row)
  deallocate(ease_col)

!overwrite the input data 
  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(vodobs_ip(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then 
           vodobs_inp(c,r) = & 
                vodobs_ip(c+(r-1)*LIS_rc%obs_lnc(k))
           
           NASASMAPvod_struc(n)%vodtime(c,r) = & 
                time
        endif
     enddo
  enddo
#endif

end subroutine read_SMAPL2vod_data

!BOP
! 
! !ROUTINE: read_NASASMAP_E_voddata
! \label{read_NASASMAP_E_voddata}
!
! !INTERFACE:
subroutine read_NASASMAP_E_voddata(n, k, pass, fname, vodobs_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use NASASMAPvod_Mod, only : NASASMAPvod_struc
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
  real                          :: vodobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the SMAP (HDF5) file and applies the data
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
!  \item[fname]        name of the NASA SMAP file
!  \item[vodobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDF5)

  character*100,    parameter    :: vod_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: vod_field_name_D = "vegetation_opacity"
  character*100,    parameter    :: vod_qa_name_D = "retrieval_qual_flag"
  character*100,    parameter    :: vod_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: vod_field_name_A = "vegetation_opacity_pm"
  character*100,    parameter    :: vod_qa_name_A = "retrieval_qual_flag_pm"

  integer(hsize_t), allocatable  :: dims(:)
  integer(hsize_t), dimension(2) :: dimsm
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hid_t)                 :: memspace
  integer(hid_t)                 :: dataspace
  integer                        :: memrank = 2
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer(hsize_t), dimension(2) :: offset_file = (/0,0/)
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: vod_gr_id_D,vod_field_id_D,vod_qa_id_D
  integer(hid_t)                 :: vod_gr_id_A,vod_field_id_A,vod_qa_id_A
  real,             allocatable  :: vod_field(:,:)
  integer,          allocatable  :: vod_qa(:,:)
  integer                        :: c,r,t
  logical*1                      :: vod_data_b(NASASMAPvod_struc(n)%nc*NASASMAPvod_struc(n)%nr)
  logical*1                      :: vodobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: vod_data(NASASMAPvod_struc(n)%nc*NASASMAPvod_struc(n)%nr)
  integer                        :: status,ios

  dimsm      = (/NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr/)
  count_file = (/NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr/)
  count_mem  = (/NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr/)
  
  allocate(vod_field(NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr))
  allocate(vod_qa(NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr))
  allocate(dims(2))

  dims(1) = NASASMAPvod_struc(n)%nc
  dims(2) = NASASMAPvod_struc(n)%nr

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LIS_verify(status, 'Error opening NASASMAP file ')
  
  if(pass.eq.'D') then 
     call h5gopen_f(file_id,vod_gr_name_D,vod_gr_id_D, status)
     call LIS_verify(status, 'Error opening SM group in NASASMAP file')
     
     call h5dopen_f(vod_gr_id_D,vod_field_name_D,vod_field_id_D, status)
     call LIS_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(vod_field_id_D, dataspace, status)
     call LIS_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LIS_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LIS_verify(status, 'Error in h5create_simple_f; read_NASASMAPvod')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_NASASMAPvod')
     
     call h5dread_f(vod_field_id_D, H5T_NATIVE_REAL,vod_field,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM field from NASASMAPfile')

     call h5dclose_f(vod_field_id_D,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

     call h5dopen_f(vod_gr_id_D,vod_qa_name_D,vod_qa_id_D, status)
     call LIS_verify(status, 'Error opening SM QA field in NASASMAP file')
     
     call h5dread_f(vod_qa_id_D, H5T_NATIVE_INTEGER,vod_qa,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM QA field from NASASMAPfile')
     
     call h5dclose_f(vod_qa_id_D,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(vod_gr_id_D,status)
     call LIS_verify(status,'Error in H5GCLOSE call')

  else
     call h5gopen_f(file_id,vod_gr_name_A,vod_gr_id_A, status)
     call LIS_verify(status, 'Error opening SM group in NASASMAP file')
     
     call h5dopen_f(vod_gr_id_A,vod_field_name_A,vod_field_id_A, status)
     call LIS_verify(status, 'Error opening SM field in NASASMAP file')
     
     call h5dget_space_f(vod_field_id_A, dataspace, status)
     call LIS_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
     
     call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          start=offset_file, count=count_file, hdferr=status)
     call LIS_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
     
     call h5screate_simple_f(memrank,dimsm, memspace, status)
     call LIS_verify(status, 'Error in h5create_simple_f; read_NASASMAPvod')
     
     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_NASASMAPvod')
     
     call h5dread_f(vod_field_id_A, H5T_NATIVE_REAL,vod_field,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM field from NASASMAPfile')

     call h5dclose_f(vod_field_id_A,status)
     call LIS_verify(status,'Error in H5DCLOSE call')
     
     call h5dopen_f(vod_gr_id_A,vod_qa_name_A,vod_qa_id_A, status)
     call LIS_verify(status, 'Error opening SM QA field in NASASMAP file')
     
     call h5dread_f(vod_qa_id_A, H5T_NATIVE_INTEGER,vod_qa,dims,status, &
          memspace, dataspace)
     call LIS_verify(status, 'Error extracting SM QA field from NASASMAPfile')

     call h5dclose_f(vod_qa_id_A,status)
     call LIS_verify(status,'Error in H5DCLOSE call')

     call h5gclose_f(vod_gr_id_A,status)
     call LIS_verify(status,'Error in H5GCLOSE call')
     
  endif
  
  call h5fclose_f(file_id,status)
  call LIS_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LIS_verify(status,'Error in H5CLOSE call')
  
  vod_data_b = .false. 
  t = 1

! The retrieval_quality_field variable's binary representation consists of bits
! that indicate whether retrieval is performed or not at a given grid cell. 
! When retrieval is performed, it contains additional bits to further 
! indicate the exit status and quality of the retrieval. The first bit 
! indicates the recommended quality (0-means retrieval has recommended quality
!

  do r=1,NASASMAPvod_struc(n)%nr
     do c=1,NASASMAPvod_struc(n)%nc        
        vod_data(t) = vod_field(c,r)
        if(vod_data(t).ne.-9999.0) then 
           if(NASASMAPvod_struc(n)%qcFlag.eq.1) then 
              if(ibits(vod_qa(c,r),0,1).eq.0) then 
                 vod_data_b(t) = .true.
              else
                 vod_data(t) = -9999.0
              endif
           else
              vod_data_b(t) = .true.
           endif
        endif
        t = t+1
     enddo
  enddo

  deallocate(vod_qa)

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       vod_data_b, vod_data, vodobs_b_ip, vodobs_ip, &
       NASASMAPvod_struc(n)%nc*NASASMAPvod_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       NASASMAPvod_struc(n)%rlat, NASASMAPvod_struc(n)%rlon,&
       NASASMAPvod_struc(n)%n11, LIS_rc%udef, ios)

  deallocate(vod_field)
  deallocate(dims)

#endif

end subroutine read_NASASMAP_E_voddata


!BOP
! 
! !ROUTINE: read_NASASMAP_voddata
! \label{read_NASASMAP_voddata}
!
! !INTERFACE:
subroutine read_NASASMAP_voddata(n, k, fname, vodobs_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use NASASMAPvod_Mod, only : NASASMAPvod_struc
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
  real                          :: vodobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the SMAP (HDF5) file and applies the data
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
!  \item[fname]        name of the NASA SMAP file
!  \item[vodobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDF5)
  character*100,    parameter    :: vod_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,    parameter    :: vod_field_name = "vegetation_opacity"
  character*100,    parameter    :: vod_qa_name = "retrieval_qual_flag"

  character*100,    parameter    :: vod_gr_name_D = "Soil_Moisture_Retrieval_Data_AM"
  character*100,    parameter    :: vod_field_name_D = "vegetation_opacity"
  character*100,    parameter    :: vod_gr_name_A = "Soil_Moisture_Retrieval_Data_PM"
  character*100,    parameter    :: vod_field_name_A = "vegetation_opacity_pm"

  integer(hid_t)                 :: file_id, vod_gr_id,vod_field_id
  integer(hid_t)                 :: vod_gr_id_D,vod_field_id_D
  integer(hid_t)                 :: vod_gr_id_A,vod_field_id_A
  integer(hid_t)                 :: dataspace
  integer(hid_t)                 :: memspace
  integer                        :: memrank = 2
  integer(hsize_t), allocatable  :: dims(:)
  integer(hsize_t), dimension(2) :: dimsm
  integer(hsize_t), dimension(2) :: offset_file
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  real,             allocatable  :: vod_field(:,:)
  real,             allocatable  :: vod_field_D(:,:)
  real,             allocatable  :: vod_field_A(:,:)
  integer,          allocatable  :: vod_qa(:,:)
  integer                        :: c,r,t
  logical*1                      :: vod_data_b(NASASMAPvod_struc(n)%nc*NASASMAPvod_struc(n)%nr)
  logical*1                      :: vodobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: vod_data(NASASMAPvod_struc(n)%nc*NASASMAPvod_struc(n)%nr)
  integer                        :: status,ios

  dimsm      = (/NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr/)
  count_file = (/NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr/)
  count_mem  = (/NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr/)
  
  allocate(vod_field(NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr))
  allocate(vod_field_D(NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr))
  allocate(vod_field_A(NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr))
!  allocate(vod_qa(NASASMAPvod_struc(n)%nc, NASASMAPvod_struc(n)%nr))
  allocate(dims(2))

  dims(1) = NASASMAPvod_struc(n)%nc
  dims(2) = NASASMAPvod_struc(n)%nr

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')
  
  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status) 
  call LIS_verify(status, 'Error opening NASASMAP file ')
  
  call h5gopen_f(file_id,vod_gr_name_D,vod_gr_id_D, status)
  call LIS_verify(status, 'Error opening SM group in NASASMAP file')
  
  call h5dopen_f(vod_gr_id_D,vod_field_name_D,vod_field_id_D, status)
  call LIS_verify(status, 'Error opening SM field in NASASMAP file')
  
  call h5dget_space_f(vod_field_id_D, dataspace, status)
  call LIS_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
       start=offset_file, count=count_file, hdferr=status)
  call LIS_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
  
  call h5screate_simple_f(memrank,dimsm, memspace, status)
  call LIS_verify(status, 'Error in h5create_simple_f; read_NASASMAPvod')
  
  call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
       start=offset_mem, count=count_mem, hdferr=status)
  call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_NASASMAPvod')
  
  call h5dread_f(vod_field_id_D, H5T_NATIVE_REAL,vod_field_D,dims,status, &
       memspace, dataspace)
  call LIS_verify(status, 'Error extracting SM field from NASASMAPfile')
 
!Read the PM (ascending) data 
  call h5gopen_f(file_id,vod_gr_name_A,vod_gr_id_A, status)
  call LIS_verify(status, 'Error opening SM group in NASASMAP file')
  
  call h5dopen_f(vod_gr_id_A,vod_field_name_A,vod_field_id_A, status)
  call LIS_verify(status, 'Error opening SM field in NASASMAP file')
  
  call h5dget_space_f(vod_field_id_A, dataspace, status)
  call LIS_verify(status, 'Error in h5dget_space_f: readNASASMAPObs')
  
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
       start=offset_file, count=count_file, hdferr=status)
  call LIS_verify(status, 'Error setting hyperslab dataspace in readNASASMAPObs')
  
  call h5screate_simple_f(memrank,dimsm, memspace, status)
  call LIS_verify(status, 'Error in h5create_simple_f; read_NASASMAPvod')
  
  call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
       start=offset_mem, count=count_mem, hdferr=status)
  call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_NASASMAPvod')
  
  call h5dread_f(vod_field_id_A, H5T_NATIVE_REAL,vod_field_A,dims,status, &
       memspace, dataspace)
  call LIS_verify(status, 'Error extracting SM field from NASASMAPfile')


  call h5dclose_f(vod_field_id_D,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

  call h5dclose_f(vod_field_id_A,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

!  call h5dopen_f(vod_gr_id,vod_qa_name,vod_qa_id, status)
!  call LIS_verify(status, 'Error opening SM QA field in NASASMAP file')
  
!  call h5dread_f(vod_qa_id, H5T_NATIVE_INTEGER,vod_qa,dims,status, &
!       memspace, dataspace)
!  call LIS_verify(status, 'Error extracting SM QA field from NASASMAPfile')
!  
!  call h5dclose_f(vod_qa_id,status)
!  call LIS_verify(status,'Error in H5DCLOSE call')
  

  call h5gclose_f(vod_gr_id_D,status)
  call LIS_verify(status,'Error in H5GCLOSE call')

  call h5gclose_f(vod_gr_id_A,status)
  call LIS_verify(status,'Error in H5GCLOSE call')
  
  call h5fclose_f(file_id,status)
  call LIS_verify(status,'Error in H5FCLOSE call')
  
  call h5close_f(status)
  call LIS_verify(status,'Error in H5CLOSE call')

  vod_field = LIS_rc%udef
  do r=1,NASASMAPvod_struc(n)%nr
     do c=1,NASASMAPvod_struc(n)%nc        
        if(vod_field_D(c,r).ne.LIS_rc%udef) then 
           vod_field(c,r) = vod_field_D(c,r)
        endif
     enddo
  enddo
  do r=1,NASASMAPvod_struc(n)%nr
     do c=1,NASASMAPvod_struc(n)%nc        
        if(vod_field_A(c,r).ne.LIS_rc%udef) then 
           vod_field(c,r) = vod_field_A(c,r)
        endif
     enddo
  enddo

  
  vod_data_b = .false. 
  t = 1

  do r=1,NASASMAPvod_struc(n)%nr
     do c=1,NASASMAPvod_struc(n)%nc        
        vod_data(t) = vod_field(c,r)
        if(vod_data(t).ne.-9999.0) then 
!           if(NASASMAPvod_struc(n)%qcFlag.eq.1) then 
!              if(ibits(vod_qa(c,r),0,1).eq.0) then 
           vod_data_b(t) = .true. 
!              endif
!           else
!              vod_data_b(t) = .true.
!           endif
        endif
        t = t+1
     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       vod_data_b, vod_data, vodobs_b_ip, vodobs_ip, &
       NASASMAPvod_struc(n)%nc*NASASMAPvod_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       NASASMAPvod_struc(n)%rlat, NASASMAPvod_struc(n)%rlon,&
       NASASMAPvod_struc(n)%n11, LIS_rc%udef, ios)

  deallocate(vod_field)
  deallocate(dims)

#endif

end subroutine read_NASASMAP_voddata


!BOP
! !ROUTINE: create_NASASMAPvod_filename
! \label{create_NASASMAPvod_filename}
! 
! !INTERFACE: 
subroutine create_NASASMAPvod_filename(ndir, designation, yr, mo,da, filename)
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
  elseif(designation.eq."SPL3SMP_E") then 
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//'.'//&
          trim(fda)//'/SMAP_L3_SM_P_E_'&
          //trim(fyr)//trim(fmo)//trim(fda)//&
          '.h5'

  endif

end subroutine create_NASASMAPvod_filename





