!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_SMAPEOPLsm
! \label{read_SMAPEOPLsm}
!
! !REVISION HISTORY:
!  6 Jun 2022: Yonghwan Kwon; Updated for use with SMAP_E_OPL soil moisture
!  23 Feb 2023: Eric Kemp; Updated to read netCDF file.
!
! !INTERFACE:
subroutine read_SMAPEOPLsm(n, k, OBS_State, OBS_Pert_State)
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
   use SMAPEOPLsm_Mod, only: SMAPEOPLsm_struc

   implicit none
! !ARGUMENTS:
   integer, intent(in) :: n
   integer, intent(in) :: k
   type(ESMF_State)    :: OBS_State
   type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!
!  reads the SMAP_E_OPL soil moisture observations.
!  The data is then rescaled to the land surface model's 
!  climatology using rescaling algorithms.
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
   real, parameter       :: MAX_SM_VALUE = 0.45, MIN_SM_VALUE = 0.0001
   integer                :: status
   integer                :: grid_index
   character(len=LIS_CONST_PATH_LEN) :: smobsdir
   character(len=LIS_CONST_PATH_LEN) :: fname
   logical                :: alarmCheck, file_exists
   integer                :: t, c, r, jj
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
   real                   :: sm_current(LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k))
   real                   :: dt
   real                   :: lon
   real                   :: lhour
   real                   :: gmt
   integer                :: zone
   integer                :: fnd
   real, allocatable      :: ssdev(:)
   character*4            :: yyyy
   character*8            :: yyyymmdd
   character*2            :: mm, dd, hh
   integer                :: yr, mo, da, hr, mn, ss
   integer                :: cyr, cmo, cda, chr, cmn, css
   integer                :: nyr, nmo, nda, nhr, nmn, nss
   real*8                 :: timenow, time1,time2,time3
   integer                :: doy
   character(len=LIS_CONST_PATH_LEN) :: list_files
   integer                :: mn_ind
   integer                :: ftn, ierr
   integer                :: rc
   character(len=3)       :: CRID
   integer, external      :: create_filelist ! C function

   call ESMF_AttributeGet(OBS_State, "Data Directory", &
                          smobsdir, rc=status)
   call LIS_verify(status)
   call ESMF_AttributeGet(OBS_State, "Data Update Status", &
                          data_update, rc=status)
   call LIS_verify(status)

   data_upd = .false.
   obs_unsc = LIS_rc%udef

   alarmCheck = LIS_isAlarmRinging(LIS_rc, "SMAP_E_OPL read alarm")

   smobs_A = LIS_rc%udef
   smobs_D = LIS_rc%udef

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

   if (alarmCheck .or. SMAPEOPLsm_struc(n)%startMode) then
      SMAPEOPLsm_struc(n)%startMode = .false.
      SMAPEOPLsm_struc(n)%smobs = LIS_rc%udef
      SMAPEOPLsm_struc(n)%smtime = -1.0

      write(yyyymmdd,'(i4.4,2i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      write(hh,'(i2.2)') LIS_rc%hr

      if(LIS_masterproc) then
         !list_files = trim(smobsdir)//'/ARFS_SM_*' &
         !             //trim(yyyymmdd)//'T'//trim(hh)//'*.dat'
         !EMK...Use netCDF
         list_files = trim(smobsdir)//'/ARFS_SM_*' &
              //trim(yyyymmdd)//'T'//trim(hh)//'*.nc'
         write(LIS_logunit,*) &
               '[INFO] Searching for ',trim(list_files)
         rc = create_filelist(trim(list_files)//char(0), &
              "SMAP_filelist.sm.dat"//char(0))
         if (rc .ne. 0) then
            write(LIS_logunit,*) &
                 '[WARN] Problem encountered when searching for SMAP files'
            write(LIS_logunit,*) &
                 'Was searching for ',trim(list_files)
            write(LIS_logunit,*) &
                 'LIS will continue...'
         endif
      end if
#if (defined SPMD)
      call mpi_barrier(lis_mpi_comm,ierr)
#endif

      ftn = LIS_getNextUnitNumber()
      open(ftn,file="./SMAP_filelist.sm.dat",status='old',iostat=ierr)

      do while(ierr.eq.0)
         read(ftn,'(a)',iostat=ierr) fname
         if(ierr.ne.0) then
            exit
         endif

         mn_ind = index(fname,trim(yyyymmdd)//'T'//trim(hh))+11
         read(fname(mn_ind:mn_ind+1),'(i2.2)') mn
         ss=0
         call LIS_tick(timenow,doy,gmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
              LIS_rc%hr, mn, ss, 0.0)

         write(LIS_logunit,*) '[INFO] reading ',trim(fname)

         call read_SMAPEOPLsm_data(n,k,fname,&
              SMAPEOPLsm_struc(n)%smobs,timenow)
      enddo
      call LIS_releaseUnitNumber(ftn)

   endif ! alram

   call ESMF_StateGet(OBS_State, "Observation01", smfield, &
                      rc=status)
   call LIS_verify(status, 'Error: StateGet Observation01')

   call ESMF_FieldGet(smfield, localDE=0, farrayPtr=obsl, rc=status)
   call LIS_verify(status, 'Error: FieldGet')

   fnd = 0
   sm_current = LIS_rc%udef

   ! dt is not defined as absolute value of the time difference to avoid
   ! double counting of the data in assimilation.

   do r=1,LIS_rc%obs_lnr(k)
      do c=1,LIS_rc%obs_lnc(k)
         if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
            grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
            dt = (SMAPEOPLsm_struc(n)%smtime(c,r)-time1)
            if(dt.ge.0.and.dt.lt.(time3-time1)) then
               sm_current(c,r) = &
                    SMAPEOPLsm_struc(n)%smobs(c,r)

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
   if (SMAPEOPLsm_struc(n)%ntimes .gt. 1 .and. SMAPEOPLsm_struc(n)%cdf_read_opt .eq. 1) then
      if (.not. SMAPEOPLsm_struc(n)%cdf_read_mon .or. LIS_rc%da .eq. 1 .and. LIS_rc%hr .eq. 0 .and. &
          LIS_rc%mn .eq. 0 .and. LIS_rc%ss .eq. 0) then
         call LIS_readMeanSigmaData(n, k, &
                                    SMAPEOPLsm_struc(n)%ntimes, &
                                    LIS_rc%obs_ngrid(k), &
                                    SMAPEOPLsm_struc(n)%modelcdffile, &
                                    "SoilMoist", &
                                    SMAPEOPLsm_struc(n)%model_mu, &
                                    SMAPEOPLsm_struc(n)%model_sigma, &
                                    LIS_rc%mo)

         call LIS_readMeanSigmaData(n, k, &
                                    SMAPEOPLsm_struc(n)%ntimes, &
                                    LIS_rc%obs_ngrid(k), &
                                    SMAPEOPLsm_struc(n)%obscdffile, &
                                    "SoilMoist", &
                                    SMAPEOPLsm_struc(n)%obs_mu, &
                                    SMAPEOPLsm_struc(n)%obs_sigma, &
                                    LIS_rc%mo)

         call LIS_readCDFdata(n, k, &
                              SMAPEOPLsm_struc(n)%nbins, &
                              SMAPEOPLsm_struc(n)%ntimes, &
                              LIS_rc%obs_ngrid(k), &
                              SMAPEOPLsm_struc(n)%modelcdffile, &
                              "SoilMoist", &
                              SMAPEOPLsm_struc(n)%model_xrange, &
                              SMAPEOPLsm_struc(n)%model_cdf, &
                              LIS_rc%mo)

         call LIS_readCDFdata(n, k, &
                              SMAPEOPLsm_struc(n)%nbins, &
                              SMAPEOPLsm_struc(n)%ntimes, &
                              LIS_rc%obs_ngrid(k), &
                              SMAPEOPLsm_struc(n)%obscdffile, &
                              "SoilMoist", &
                              SMAPEOPLsm_struc(n)%obs_xrange, &
                              SMAPEOPLsm_struc(n)%obs_cdf, &
                              LIS_rc%mo)

         SMAPEOPLsm_struc(n)%cdf_read_mon = .true.
      endif
   endif

   if (LIS_rc%dascaloption(k) .eq. "CDF matching" .and. fnd .ne. 0) then
      if (SMAPEOPLsm_struc(n)%ntimes .gt. 1 .and. SMAPEOPLsm_struc(n)%cdf_read_opt .eq. 1) then
         call LIS_rescale_with_CDF_matching( &
            n, k, &
            SMAPEOPLsm_struc(n)%nbins, &
            1, &
            MAX_SM_VALUE, &
            MIN_SM_VALUE, &
            SMAPEOPLsm_struc(n)%model_xrange, &
            SMAPEOPLsm_struc(n)%obs_xrange, &
            SMAPEOPLsm_struc(n)%model_cdf, &
            SMAPEOPLsm_struc(n)%obs_cdf, &
            sm_current)
      else
         call LIS_rescale_with_CDF_matching( &
            n, k, &
            SMAPEOPLsm_struc(n)%nbins, &
            SMAPEOPLsm_struc(n)%ntimes, &
            MAX_SM_VALUE, &
            MIN_SM_VALUE, &
            SMAPEOPLsm_struc(n)%model_xrange, &
            SMAPEOPLsm_struc(n)%obs_xrange, &
            SMAPEOPLsm_struc(n)%model_cdf, &
            SMAPEOPLsm_struc(n)%obs_cdf, &
            sm_current)
      endif
   elseif(LIS_rc%dascaloption(k).eq."Linear scaling".and.fnd.ne.0) then
        call LIS_rescale_with_linear_scaling(    &
             n,                                   &
             k,                                   &
             SMAPEOPLsm_struc(n)%nbins,         &
             SMAPEOPLsm_struc(n)%ntimes,        &
             SMAPEOPLsm_struc(n)%obs_xrange,    &
             SMAPEOPLsm_struc(n)%obs_cdf,       &
             sm_current)
   elseif(LIS_rc%dascaloption(k).eq."Anomaly scaling".and.fnd.ne.0) then
        call LIS_rescale_with_anomaly(    &
             n,                                   &
             k,                                   &
             SMAPEOPLsm_struc(n)%nbins,         &
             SMAPEOPLsm_struc(n)%ntimes,        &
             SMAPEOPLsm_struc(n)%obs_mu,    &
             SMAPEOPLsm_struc(n)%model_mu,       &
             sm_current)
   endif

   obsl = LIS_rc%udef
   do r = 1, LIS_rc%obs_lnr(k)
      do c = 1, LIS_rc%obs_lnc(k)
         if (LIS_obs_domain(n, k)%gindex(c, r) .ne. -1) then
            obsl(LIS_obs_domain(n, k)%gindex(c, r)) = sm_current(c, r)
         endif
      enddo
   enddo
   !-------------------------------------------------------------------------
   !  Apply LSM based QC and screening of observations
   !-------------------------------------------------------------------------
   call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+" &
                        //trim(LIS_SMAPEOPLsmobsId)//char(0), n, k, OBS_state)

   call LIS_checkForValidObs(n, k, obsl, fnd, sm_current)

   if (fnd .eq. 0) then
      data_upd_flag_local = .false.
   else
      data_upd_flag_local = .true.
   endif

#if (defined SPMD)
   call MPI_ALLGATHER(data_upd_flag_local, 1, &
                      MPI_LOGICAL, data_upd_flag(:), &
                      1, MPI_LOGICAL, LIS_mpi_comm, status)
   data_upd = any(data_upd_flag)
#else
   data_upd = data_upd_flag_local
#endif

   if (data_upd) then

      do t = 1, LIS_rc%obs_ngrid(k)
         gid(t) = t
         if (obsl(t) .ne. -9999.0) then
            assimflag(t) = 1
         else
            assimflag(t) = 0
         endif
      enddo

      call ESMF_AttributeSet(OBS_State, "Data Update Status", &
                             .true., rc=status)
      call LIS_verify(status)

      if (LIS_rc%obs_ngrid(k) .gt. 0) then
         call ESMF_AttributeSet(smField, "Grid Number", &
                                gid, itemCount=LIS_rc%obs_ngrid(k), rc=status)
         call LIS_verify(status)

         call ESMF_AttributeSet(smField, "Assimilation Flag", &
                                assimflag, itemCount=LIS_rc%obs_ngrid(k), rc=status)
         call LIS_verify(status)

         call ESMF_AttributeSet(smfield, "Unscaled Obs", &
                                obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
         call LIS_verify(status, 'Error in setting Unscaled Obs attribute')

      endif

      if (LIS_rc%dascaloption(k) .eq. "CDF matching") then
         if (SMAPEOPLsm_struc(n)%useSsdevScal .eq. 1) then
            call ESMF_StateGet(OBS_Pert_State, "Observation01", pertfield, &
                               rc=status)
            call LIS_verify(status, 'Error: StateGet Observation01')

            allocate (ssdev(LIS_rc%obs_ngrid(k)))
            ssdev = SMAPEOPLsm_struc(n)%ssdev_inp

            if (SMAPEOPLsm_struc(n)%ntimes .eq. 1) then
               jj = 1
            else
               jj = LIS_rc%mo
            endif
            do t = 1, LIS_rc%obs_ngrid(k)
               if (SMAPEOPLsm_struc(n)%obs_sigma(t, jj) .gt. 0) then
                  ssdev(t) = ssdev(t)*SMAPEOPLsm_struc(n)%model_sigma(t, jj)/ &
                             SMAPEOPLsm_struc(n)%obs_sigma(t, jj)
                  if (ssdev(t) .lt. minssdev) then
                     ssdev(t) = minssdev
                  endif
               endif
            enddo

            if (LIS_rc%obs_ngrid(k) .gt. 0) then
               call ESMF_AttributeSet(pertField, "Standard Deviation", &
                                      ssdev, itemCount=LIS_rc%obs_ngrid(k), rc=status)
               call LIS_verify(status)
            endif
            deallocate (ssdev)
         endif
      endif
   else
      call ESMF_AttributeSet(OBS_State, "Data Update Status", &
                             .false., rc=status)
      call LIS_verify(status)
   endif

end subroutine read_SMAPEOPLsm

!BOP
! 
! !ROUTINE: read_SMAPEOPLsm_data
! \label{read_SMAPEOPLsm_data}
!
! !INTERFACE:
subroutine read_SMAPEOPLsm_data(n, k,fname, smobs_inp, time)
! 
! !USES: 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use SMAPEOPLsm_Mod, only : SMAPEOPLsm_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n
  integer                  :: k
  character (len=*)        :: fname
  real                     :: smobs_inp(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real*8                   :: time
!EOP
  integer,  parameter     :: nc=2560, nr=1920
  real*4                  :: sm_raw(SMAPEOPLsm_struc(n)%nc,SMAPEOPLsm_struc(n)%nr)
  real                    :: sm_in(SMAPEOPLsm_struc(n)%nc*SMAPEOPLsm_struc(n)%nr)
  real                    :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  logical*1               :: sm_data_b(SMAPEOPLsm_struc(n)%nc*SMAPEOPLsm_struc(n)%nr)
  logical*1               :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                 :: smid
  integer                 :: ios, nid
  integer                 :: c,r
  integer                 :: ftn1
  ! EMK
  logical :: file_exists
  character(255) :: map_projection
  integer :: ncid, dim_ids(3), var_id
  integer :: ntime, nlat, nlon
  real, allocatable :: tmp(:,:,:)
  integer :: rc

  ! Old code to use binary data
  !ftn1 = LIS_getNextUnitNumber()
  !open(unit=ftn1,file=fname,form='unformatted',access='direct',recl=4*nc*nr,status='old')

  ! read(ftn1, rec=1) sm_raw
  ! close(1)
  ! call LIS_releaseUnitNumber(ftn1)

  ! sm_data_b = .false.

  ! do r=1,SMAPEOPLsm_struc(n)%nr
  !    do c=1,SMAPEOPLsm_struc(n)%nc
  !       if (sm_raw(c,r)>=0.and.&
  !          sm_raw(c,r)<=100) then

  !          sm_in(c+(r-1)*SMAPEOPLsm_struc(n)%nc) = sm_raw(c,r)
  !          sm_data_b(c+(r-1)*SMAPEOPLsm_struc(n)%nc) = .true.
  !       else
  !          sm_in(c+(r-1)*SMAPEOPLsm_struc(n)%nc) = LIS_rc%udef
  !          sm_data_b(c+(r-1)*SMAPEOPLsm_struc(n)%nc) = .false.
  !       endif
  !    enddo
  ! enddo

  ! EMK...Read netCDF data
  sm_in = LIS_rc%udef
  sm_data_b = .false.

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ! See if the file exists.
  inquire(file=trim(fname), exist=file_exists)
  if (.not. file_exists) then
     write(LIS_logunit,*)'[ERR] Cannot find ', trim(fname)
     return
  end if

  ! Open the file
  rc = nf90_open(path=trim(fname), &
       mode=NF90_NOWRITE, &
       ncid=ncid)
  if (rc .ne. 0) then
     write(LIS_logunit,*)'[ERR] Cannot open ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     return
  end if

  ! Read the map projection
  rc = nf90_get_att(ncid=ncid, &
       varid=NF90_GLOBAL, &
       name='MAP_PROJECTION', &
       values=map_projection)
  if (rc .ne. 0) then
     write(LIS_logunit,*)'[ERR] Cannot read MAP_PROJECTION from ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Sanity check map projection
  ! TODO:  Support other map projections
  if (trim(map_projection) .ne. 'EQUIDISTANT CYLINDRICAL') then
     write(LIS_logunit,*) &
          '[ERR] Unrecognized map projection found in SMAP file!'
     write(LIS_logunit,*) '[ERR] Expected EQUIDISTANT CYLINDRICAL'
     write(LIS_logunit,*) '[ERR] Found ',trim(map_projection)
     write(LIS_logunit,*) '[ERR] LIS will skip file ', trim(fname)
     rc = nf90_close(ncid)
     return
  end if

  ! Get dimension IDs
  rc = nf90_inq_dimid(ncid=ncid, &
       name='time', &
       dimid=dim_ids(3))
  if (rc .ne. 0) then
     write(LIS_logunit,*)'[ERR] Cannot read time dimension from ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inq_dimid(ncid=ncid, &
       name='lat', &
       dimid=dim_ids(2))
  if (rc .ne. 0) then
     write(LIS_logunit,*)'[ERR] Cannot read lat dimension from ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inq_dimid(ncid=ncid, &
       name='lon', &
       dimid=dim_ids(1))
  if (rc .ne. 0) then
     write(LIS_logunit,*)'[ERR] Cannot read lon dimension from ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Get actual dimension sizes
  rc = nf90_inquire_dimension(ncid=ncid, &
       dimid=dim_ids(3), &
       len=ntime)
  if (rc .ne. 0) then
     write(LIS_logunit,*)'[ERR] Cannot read time dimension from ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inquire_dimension(ncid=ncid, &
       dimid=dim_ids(2), &
       len=nlat)
  if (rc .ne. 0) then
     write(LIS_logunit,*)'[ERR] Cannot read lat dimension from ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inquire_dimension(ncid=ncid, &
       dimid=dim_ids(1), &
       len=nlon)
  if (rc .ne. 0) then
     write(LIS_logunit,*)'[ERR] Cannot read lon dimension from ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Sanity check the dimensions
  if (ntime .ne. 1) then
     write(LIS_logunit,*)'[ERR] Expected time dimension to be 1'
     write(LIS_logunit,*)'[ERR] Found ', ntime, ' from ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if
  if (nlat .ne. SMAPEOPLsm_struc(n)%nr) then
     write(LIS_logunit,*)'[ERR] Expected lat dimension to be ', &
          SMAPEOPLsm_struc(n)%nr
     write(LIS_logunit,*)'[ERR] Found ', nlat, ' from ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if
  if (nlon .ne. SMAPEOPLsm_struc(n)%nc) then
     write(LIS_logunit,*)'[ERR] Expected lon dimension to be ', &
          SMAPEOPLsm_struc(n)%nc
     write(LIS_logunit,*)'[ERR] Found ', nlon, ' from ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Fetch the variable id
  rc = nf90_inq_varid(ncid=ncid, &
       name='arfs_sm', &
       varid=var_id)
  if (rc .ne. 0) then
     write(LIS_logunit,*)'[ERR] Cannot read arfs_sm ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Read the retrievals
  allocate(tmp(nlon, nlat, ntime))
  rc = nf90_get_var(ncid=ncid, &
       varid=var_id, &
       values=tmp)
  if (rc .ne. 0) then
     write(LIS_logunit,*)'[ERR] Cannot read arfs_sm ', trim(fname)
     write(LIS_logunit,*)'[ERR] LIS will continue...'
     rc = nf90_close(ncid)
     deallocate(tmp)
     return
  end if
  rc = nf90_close(ncid)

  do r = 1, nlat
     do c = 1, nlon
        if (tmp(c,r,1) >= 0 .and. &
             tmp(c,r,1) <= 1) then
           sm_in(c + (r-1)*nc) = tmp(c,r,1)*100
           sm_data_b(c + (r-1)*nc) = .true.
        end if
     end do
  end do
  deallocate(tmp)

#endif

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!--------------------------------------------------------------------------
  if(LIS_rc%obs_gridDesc(k,10).le.0.0937500) then
     call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
          sm_data_b, sm_in, smobs_b_ip, smobs_ip, &
          SMAPEOPLsm_struc(n)%nc*SMAPEOPLsm_struc(n)%nr, &
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          SMAPEOPLsm_struc(n)%rlat,SMAPEOPLsm_struc(n)%rlon,&
          SMAPEOPLsm_struc(n)%w11,SMAPEOPLsm_struc(n)%w12,&
          SMAPEOPLsm_struc(n)%w21,SMAPEOPLsm_struc(n)%w22,&
          SMAPEOPLsm_struc(n)%n11,SMAPEOPLsm_struc(n)%n12,&
          SMAPEOPLsm_struc(n)%n21,SMAPEOPLsm_struc(n)%n22,LIS_rc%udef,ios)
  else
     call upscaleByAveraging(SMAPEOPLsm_struc(n)%nc*SMAPEOPLsm_struc(n)%nr,&
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          LIS_rc%udef, SMAPEOPLsm_struc(n)%n11,&
          sm_data_b,sm_in, smobs_b_ip, smobs_ip)
  endif

!overwrite the input data
  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(smobs_ip(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then
           smobs_inp(c,r) = &
                smobs_ip(c+(r-1)*LIS_rc%obs_lnc(k))

           SMAPEOPLsm_struc(n)%smtime(c,r) = &
                time
        endif
     enddo
  enddo

end subroutine read_SMAPEOPLsm_data
