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
! !ROUTINE: read_SMOSNRTNNL2sm
! \label{read_SMOSNRTNNL2sm}
!
! !REVISION HISTORY:
!  8  May 2013: Sujay Kumar; initial specification
!  20 Jan 2021: Yonghwan Kwon; Updated to read SMOS NRT NN L2 soil moisture
!  20 Feb 2021: Mahdi Navari and Eric Kemp ; updated to read the DGG number
!               from the LDT output
!
! !INTERFACE:
subroutine read_SMOSNRTNNL2sm(n, k, OBS_State, OBS_Pert_State)
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
   use SMOSNRTNNL2sm_Mod, only: SMOSNRTNNL2sm_struc, SMOS_in_lis_gridbox
   use read_dgg_lookup_table

   implicit none
! !ARGUMENTS:
   integer, intent(in) :: n
   integer, intent(in) :: k
   type(ESMF_State)    :: OBS_State
   type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!
!  reads the SMOS NRT NN L2 soil moisture observations
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
   real, parameter        :: MAX_SM_VALUE = 0.45, MIN_SM_VALUE = 0.0001
   integer                :: status
   integer                :: grid_index
   character(len=LIS_CONST_PATH_LEN) :: smobsdir
   character(len=LIS_CONST_PATH_LEN) :: fname
   logical                :: alarmCheck, file_exists
   integer                :: t, c, r, i, j, p, jj
   real, pointer          :: obsl(:)
   type(ESMF_Field)       :: smfield, pertField
   integer                :: gid(LIS_rc%obs_ngrid(k))
   integer                :: assimflag(LIS_rc%obs_ngrid(k))
   real                   :: obs_unsc(LIS_rc%obs_ngrid(k))
   logical                :: data_update
   logical                :: data_upd_flag(LIS_npes)
   logical                :: data_upd_flag_local
   logical                :: data_upd
   !real                   :: smobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
   !real                   :: smtime(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
   real                   :: sm_current(LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k))
   real                   :: dt
   integer                :: fnd
   real                   :: gmt
   real, allocatable      :: ssdev(:)
   integer                :: lis_julss
   real                   :: smvalue
   real                   :: model_delta(LIS_rc%obs_ngrid(k))
   real                   :: obs_delta(LIS_rc%obs_ngrid(k))
   integer                :: yr, mo, da, hr, mn, ss
   integer                :: cyr, cmo, cda, chr, cmn, css
   integer                :: nyr, nmo, nda, nhr, nmn, nss
   character*4            :: yyyy
   character*8            :: yyyymmdd
   character*2            :: mm, dd, hh, hh0, hh1
   character(len=LIS_CONST_PATH_LEN) :: list_files
   real*8                 :: timenow, time1,time2,time3
   integer                :: doy
   integer                :: mn_ind
   integer                :: ftn, ierr
   character*1            :: fproc(4)
   character*100          :: temp1
   !real, allocatable    :: sm_dgg(:), smunct_dgg(:), obstime_dgg(:)
   integer, allocatable   :: num_dgg_glb(:,:)
   integer                :: Max_length
   !logical, allocatable   :: SMOS_assign_glb(:,:)
   integer, allocatable :: dgg_lookup_1d(:)
   type(SMOS_in_lis_gridbox), pointer :: SMOS_lookup_glb_1d(:)
   integer :: num_indices
   integer :: gindex
   integer :: ii , gc, gr
   integer :: leng
   logical, save :: dgg_file_already_read = .false. 
   integer, external :: create_filelist ! C function

   call ESMF_AttributeGet(OBS_State, "Data Directory", &
                          smobsdir, rc=status)
   call LIS_verify(status)
   call ESMF_AttributeGet(OBS_State, "Data Update Status", &
                          data_update, rc=status)
   call LIS_verify(status)

   data_upd = .false.
   obs_unsc = LIS_rc%udef

   if (SMOSNRTNNL2sm_struc(n)%start_day.ne.LIS_rc%da) then
      SMOSNRTNNL2sm_struc(n)%count_day = SMOSNRTNNL2sm_struc(n)%count_day + 1
      SMOSNRTNNL2sm_struc(n)%start_day = LIS_rc%da
   endif

   !read dgg lookup table form LDT output and store on local processes
   if (.not. dgg_file_already_read) then ! EMK 
      leng = 0
      If (LIS_masterproc) then
         call LIS_SMOS_DGG_lookup(n, dgg_lookup_1d)
         leng = size(dgg_lookup_1d)
      end if
#if (defined SPMD)
      ! For multiple processes, need to forward copy of array from the
      ! master process.
      call MPI_Bcast(leng, 1, MPI_INTEGER, 0, LIS_mpi_comm, ierr)
      if (.not. LIS_masterproc) then
         allocate(dgg_lookup_1d(leng))
      end if
      call MPI_Bcast(dgg_lookup_1d, leng, MPI_INTEGER, 0, LIS_mpi_comm, ierr)
#endif

      ! EMK Each process has a complete copy of the DGG array, but it is
      ! still run-length encoded.  It is convenient to convert to 1-d array
      ! of structures.
      allocate(SMOS_lookup_glb_1d(LIS_rc%gnc(n)*LIS_rc%gnr(n)))

      i = 1
      do r = 1, LIS_rc%gnr(n)
         do c = 1, LIS_rc%gnc(n)
            ii = c + (r-1)*LIS_rc%gnc(n)
            num_indices = dgg_lookup_1d(i)
            if (num_indices .eq. 0) then
               SMOS_lookup_glb_1d(ii)%dgg_assign = .false.
               i = i + 1
            else
               SMOS_lookup_glb_1d(ii)%dgg_assign = .true.
               allocate( &
                    SMOS_lookup_glb_1d(ii)%dgg_indices(num_indices))
               Do j = 1, num_indices
                  SMOS_lookup_glb_1d(ii)%dgg_indices(j) = &
                       dgg_lookup_1d(i+j)
               end do
               i = i + num_indices + 1
            end if
            if (i .gt. size(dgg_lookup_1d)) exit
         end do
         if (i .gt. size(dgg_lookup_1d)) exit
      end do
      deallocate(dgg_lookup_1d)

      ! EMK Each process now has a complete copy of the DGG array, with
      ! DGG indices easily accessible by GID.  We now pull out the values
      ! relevant to the local process.
      if (.not. associated(SMOSNRTNNL2sm_struc(n)%SMOS_lookup)) then
         allocate(SMOSNRTNNL2sm_struc(n)%SMOS_lookup( &
              LIS_rc%lnc(n), LIS_rc%lnr(n)))
      end if

      ! EMK: Now copy to local structure
      do r = 1, LIS_rc%lnr(n)
         do c = 1, LIS_rc%lnc(n)
            !gindex = LIS_domain(n)%gindex(c,r) ! EMK implementation
            gr = r + LIS_nss_ind(n,LIS_localPet+1) -1  
            gc = c + LIS_ews_ind(n,LIS_localPet+1) -1
            gindex = gc + (gr-1)*LIS_rc%gnc(n)

            if (gindex .eq. -1) cycle

            SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_assign = &
                 SMOS_lookup_glb_1d(gindex)%dgg_assign

            if (SMOS_lookup_glb_1d(gindex)%dgg_assign .eqv. .true.) then
               num_indices = size(SMOS_lookup_glb_1d(gindex)%dgg_indices)
               allocate(SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)% &
                    dgg_indices(num_indices))
               SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)% &
                    dgg_indices = &
                    SMOS_lookup_glb_1d(gindex)%dgg_indices
            end if
         end do
      end do

      ! EMK Clean up memory
      do i = 1, LIS_rc%gnc(n) * LIS_rc%gnr(n)
         if (SMOS_lookup_glb_1d(i)%dgg_assign .eqv. .true.) then
            deallocate(SMOS_lookup_glb_1d(i)%dgg_indices)
         end if
      end do
      deallocate(SMOS_lookup_glb_1d)

      ! EMK Ensure DGG file not read again
      dgg_file_already_read = .true.
   end if


   alarmCheck = LIS_isAlarmRinging(LIS_rc, "SMOSNRTNN read alarm")

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

   if (alarmCheck .or. SMOSNRTNNL2sm_struc(n)%startMode) then
      SMOSNRTNNL2sm_struc(n)%startMode = .false.

      SMOSNRTNNL2sm_struc(n)%smobs = LIS_rc%udef
      SMOSNRTNNL2sm_struc(n)%smtime = -1.0

      !smobs = LIS_rc%udef

      write(temp1,fmt='(i4.4)') LIS_localPet
      read(temp1,fmt='(4a1)') fproc
      write(yyyymmdd,'(i4.4,2i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
      write(yyyy,'(i4.4)') LIS_rc%yr
      write(mm,'(i2.2)') LIS_rc%mo
      write(dd,'(i2.2)') LIS_rc%da
      write(hh,'(i2.2)') LIS_rc%hr
      if (LIS_rc%hr-1 >= 0) then
         write(hh0,'(i2.2)') LIS_rc%hr-1
      endif
      write(hh1,'(i2.2)') LIS_rc%hr+1

      if(LIS_masterproc) then


         list_files = trim(smobsdir)// &
              '/'//trim(yyyy)//'.'//trim(mm)//'.'//dd// &
              '/W_XX-ESA,SMOS,NRTNN_C_LEMM_*_' &
              //trim(yyyymmdd)//trim(hh) &
              //'*.nc'

         status = create_filelist(trim(list_files)//char(0), &
                                  "file_01.txt"//char(0))
         if (status .ne. 0) then
            write(LIS_logunit,*) &
                 '[WARN] Problem encountered when searching for SMOS NRTNNL2 files'
            write(LIS_logunit,*) &
                 'Was searching for ',trim(list_files)
            write(LIS_logunit,*) &
                 'LIS will continue...'
         endif

         if (LIS_rc%hr-1 >= 0.and.LIS_rc%hr+1 <= 23) then
            list_files = trim(smobsdir)// &
                 '/'//trim(yyyy)//'.'//trim(mm)//'.'//dd// &
                 '/W_XX-ESA,SMOS,NRTNN_C_LEMM_*_' &
                 //trim(yyyymmdd)//trim(hh0)//'*_'//trim(yyyymmdd)//(hh1)  &
                 //'*.nc'

            status = create_filelist(trim(list_files)//char(0), &
                                     "file_02.txt"//char(0))
            if (status .ne. 0) then
               write(LIS_logunit,*) &
                    '[WARN] Problem encountered when searching for SMOS NRTNNL2 files'
               write(LIS_logunit,*) &
                    'Was searching for ',trim(list_files)
               write(LIS_logunit,*) &
                    'LIS will continue...'
            endif
         endif

         call system('cat file_01.txt file_02.txt >./SMOS_filelist_sm.dat 2>/dev/null')
         call system('rm -f file_01.txt file_02.txt')
      end if

#if (defined SPMD)
         call mpi_barrier(lis_mpi_comm,ierr)
#endif

         ftn = LIS_getNextUnitNumber()
         open(ftn,file="./SMOS_filelist_sm.dat", &
              action='read', status='old', iostat=ierr)

         do while(ierr.eq.0)
            read(ftn,'(a)',iostat=ierr) fname
            if(ierr.ne.0) then
               exit
            endif

            ss=0

            write(LIS_logunit,*) '[INFO] reading ',trim(fname)

            call read_SMOSNRTNNL2sm_data(n,k,fname,&
                 SMOSNRTNNL2sm_struc(n)%smobs,time1,chr)
         enddo
         call LIS_releaseUnitNumber(ftn)

   endif !alarm

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
            dt = (SMOSNRTNNL2sm_struc(n)%smtime(c,r)-time1)

            if(dt.ge.0.and.dt.lt.(time3-time1)) then

               sm_current(c,r) = &
                   SMOSNRTNNL2sm_struc(n)%smobs(c,r)
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
   if (SMOSNRTNNL2sm_struc(n)%ntimes .gt. 1 .and. SMOSNRTNNL2sm_struc(n)%cdf_read_opt .eq. 1) then
      if (.not. SMOSNRTNNL2sm_struc(n)%cdf_read_mon .or. LIS_rc%da .eq. 1 .and. LIS_rc%hr .eq. 0 .and. &
          LIS_rc%mn .eq. 0 .and. LIS_rc%ss .eq. 0) then
         call LIS_readMeanSigmaData(n, k, &
                                    SMOSNRTNNL2sm_struc(n)%ntimes, &
                                    LIS_rc%obs_ngrid(k), &
                                    SMOSNRTNNL2sm_struc(n)%modelcdffile, &
                                    "SoilMoist", &
                                    SMOSNRTNNL2sm_struc(n)%model_mu, &
                                    SMOSNRTNNL2sm_struc(n)%model_sigma, &
                                    LIS_rc%mo)

         call LIS_readMeanSigmaData(n, k, &
                                    SMOSNRTNNL2sm_struc(n)%ntimes, &
                                    LIS_rc%obs_ngrid(k), &
                                    SMOSNRTNNL2sm_struc(n)%obscdffile, &
                                    "SoilMoist", &
                                    SMOSNRTNNL2sm_struc(n)%obs_mu, &
                                    SMOSNRTNNL2sm_struc(n)%obs_sigma, &
                                    LIS_rc%mo)

         call LIS_readCDFdata(n, k, &
                              SMOSNRTNNL2sm_struc(n)%nbins, &
                              SMOSNRTNNL2sm_struc(n)%ntimes, &
                              LIS_rc%obs_ngrid(k), &
                              SMOSNRTNNL2sm_struc(n)%modelcdffile, &
                              "SoilMoist", &
                              SMOSNRTNNL2sm_struc(n)%model_xrange, &
                              SMOSNRTNNL2sm_struc(n)%model_cdf, &
                              LIS_rc%mo)

         call LIS_readCDFdata(n, k, &
                              SMOSNRTNNL2sm_struc(n)%nbins, &
                              SMOSNRTNNL2sm_struc(n)%ntimes, &
                              LIS_rc%obs_ngrid(k), &
                              SMOSNRTNNL2sm_struc(n)%obscdffile, &
                              "SoilMoist", &
                              SMOSNRTNNL2sm_struc(n)%obs_xrange, &
                              SMOSNRTNNL2sm_struc(n)%obs_cdf, &
                              LIS_rc%mo)

         SMOSNRTNNL2sm_struc(n)%cdf_read_mon = .true.
      endif
   endif

   if (LIS_rc%dascaloption(k) .eq. "CDF matching" .and. fnd .ne. 0) then
      if (SMOSNRTNNL2sm_struc(n)%ntimes .gt. 1 .and. SMOSNRTNNL2sm_struc(n)%cdf_read_opt .eq. 1) then
         call LIS_rescale_with_CDF_matching( &
            n, k, &
            SMOSNRTNNL2sm_struc(n)%nbins, &
            1, &
            MAX_SM_VALUE, &
            MIN_SM_VALUE, &
            SMOSNRTNNL2sm_struc(n)%model_xrange, &
            SMOSNRTNNL2sm_struc(n)%obs_xrange, &
            SMOSNRTNNL2sm_struc(n)%model_cdf, &
            SMOSNRTNNL2sm_struc(n)%obs_cdf, &
            sm_current)
      else
         call LIS_rescale_with_CDF_matching( &
            n, k, &
            SMOSNRTNNL2sm_struc(n)%nbins, &
            SMOSNRTNNL2sm_struc(n)%ntimes, &
            MAX_SM_VALUE, &
            MIN_SM_VALUE, &
            SMOSNRTNNL2sm_struc(n)%model_xrange, &
            SMOSNRTNNL2sm_struc(n)%obs_xrange, &
            SMOSNRTNNL2sm_struc(n)%model_cdf, &
            SMOSNRTNNL2sm_struc(n)%obs_cdf, &
            sm_current)
      endif
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
                        //trim(LIS_SMOSNRTNNL2smobsId)//char(0), n, k, OBS_state)

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
         if (SMOSNRTNNL2sm_struc(n)%useSsdevScal .eq. 1) then
            call ESMF_StateGet(OBS_Pert_State, "Observation01", pertfield, &
                               rc=status)
            call LIS_verify(status, 'Error: StateGet Observation01')

            allocate (ssdev(LIS_rc%obs_ngrid(k)))
            ssdev = SMOSNRTNNL2sm_struc(n)%ssdev_inp

            if (SMOSNRTNNL2sm_struc(n)%ntimes .eq. 1) then
               jj = 1
            else
               jj = LIS_rc%mo
            endif
            do t = 1, LIS_rc%obs_ngrid(k)
               if (SMOSNRTNNL2sm_struc(n)%obs_sigma(t, jj) .gt. 0) then
                  ssdev(t) = ssdev(t)*SMOSNRTNNL2sm_struc(n)%model_sigma(t, jj)/ &
                             SMOSNRTNNL2sm_struc(n)%obs_sigma(t, jj)
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

end subroutine read_SMOSNRTNNL2sm

!BOP
! 
! !ROUTINE: read_SMOSNRTNNL2sm_data
! \label{read_SMOSNRTNNL2sm_data}
!
! !INTERFACE:
subroutine read_SMOSNRTNNL2sm_data(n, k, fname, smobs_inp, time, chr)
! 
! !USES: 

  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_DAobservationsMod
  use SMOSNRTNNL2sm_Mod, only : SMOSNRTNNL2sm_struc
  !use read_dgg_lookup_table

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n
  integer                  :: k
  character (len=*)        :: fname
  real                     :: smobs_inp(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real*8                   :: time
  integer                  :: chr

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!
!
!EOP
   integer              :: ios
   integer              :: nid, smid, dimid, lonid, latid, smunctid, obstimeid, dggid
   integer              :: dgg_length
   real, allocatable    :: sm(:), smunct(:), obstime(:)
   real, allocatable    :: lon_dgg(:), lat_dgg(:)
   integer, allocatable :: DGG_id_number(:)
   integer              :: max_dgg_id_number, min_dgg_id_number
   real                 :: lat2d(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)), &
                           lon2d(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
   real                 :: smobs_sum(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
   integer              :: count_smobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
   real                 :: dx, dy
   real                 :: max_lon_dgg, min_lon_dgg, max_lat_dgg, min_lat_dgg
   character (len = 20) :: dimname
   integer              :: c,r,i,j
   integer, allocatable :: i_dgg(:)
   real, allocatable    :: sm_dgg(:), smunct_dgg(:), obstime_dgg(:)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

   ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
   call LIS_verify(ios,'Error opening file '//trim(fname))

   ios = nf90_inq_varid(nid, 'longitude',lonid)
   call LIS_verify(ios, 'Error nf90_inq_varid: longitude')

   ios = nf90_inq_varid(nid, 'latitude',latid)
   call LIS_verify(ios, 'Error nf90_inq_varid: latitude')

   ios = nf90_inq_varid(nid, 'soil_moisture',smid)
   call LIS_verify(ios, 'Error nf90_inq_varid: soil_moisture')

   ios = nf90_inq_varid(nid, 'soil_moisture_uncertainty',smunctid)
   call LIS_verify(ios, 'Error nf90_inq_varid: soil_moisture_uncertainty')

   ios = nf90_inq_varid(nid, 'seconds_since_midnight',obstimeid)
   call LIS_verify(ios, 'Error nf90_inq_varid: seconds_since_midnight')

   ios = nf90_inq_varid(nid, 'DGG_id_number',dggid)
   call LIS_verify(ios, 'Error nf90_inq_varid: DGG_id_number')

   !dimension
   ios = nf90_inq_dimid(nid, 'DGG_id_number', dimid)
   call LIS_verify(ios, 'Error nf90_inq_varid: soil_moisture dimension id')

   ios = nf90_inquire_dimension(nid, dimid, dimname, dgg_length)
   call LIS_verify(ios, 'Error nf90_inq_varid: soil_moisture dimension')

   !values
   if (dgg_length > 0) then
      allocate(sm(dgg_length))
      allocate(smunct(dgg_length))
      allocate(lon_dgg(dgg_length))
      allocate(lat_dgg(dgg_length))
      allocate(obstime(dgg_length))
      allocate(DGG_id_number(dgg_length))
      allocate(i_dgg(dgg_length))

      ios = nf90_get_var(nid, smid, sm)
      call LIS_verify(ios, 'Error nf90_get_var: soil_moisture')

      ios = nf90_get_var(nid, smunctid, smunct)
      call LIS_verify(ios, 'Error nf90_get_var: soil_moisture_uncertainty')

      ios = nf90_get_var(nid, lonid, lon_dgg)
      call LIS_verify(ios, 'Error nf90_get_var: longitude')

      ios = nf90_get_var(nid, latid, lat_dgg)
      call LIS_verify(ios, 'Error nf90_get_var: latitude')

      ios = nf90_get_var(nid, obstimeid, obstime)
      call LIS_verify(ios, 'Error nf90_get_var: seconds_since_midnight')

      ios = nf90_get_var(nid, dggid, DGG_id_number)
      call LIS_verify(ios, 'Error nf90_get_var: DGG_id_number')

      ! find max/min lon_dgg/lon_dgg
      max_lon_dgg = maxval(lon_dgg)
      min_lon_dgg = minval(lon_dgg)
      max_lat_dgg = maxval(lat_dgg)
      min_lat_dgg = minval(lat_dgg)

      ! find min/max DGG_id_number in a file
      max_dgg_id_number = maxval(DGG_id_number)
      min_dgg_id_number = minval(DGG_id_number)

      if (min_dgg_id_number > 0) then
         allocate(sm_dgg(max_dgg_id_number))
         allocate(smunct_dgg(max_dgg_id_number))
         allocate(obstime_dgg(max_dgg_id_number))

         sm_dgg = LIS_rc%udef
         smunct_dgg = LIS_rc%udef
         obstime_dgg = LIS_rc%udef

         ! put sm to corresponding DGG_id_number
         sm_dgg(DGG_id_number(:)) = sm(:)
         smunct_dgg(DGG_id_number(:)) = smunct(:)
         obstime_dgg(DGG_id_number(:)) = obstime(:)

         ! Put the SMOS soil moisture data in LIS input grid
         smobs_sum = 0 !initialize
         count_smobs = 0 !initialize
         i_dgg = (/(j, j = 1,dgg_length) /)

         do r=1,LIS_rc%obs_lnr(k)
            do c=1,LIS_rc%obs_lnc(k)
               lon2d(c,r) = LIS_obs_domain(n, k)%lon(c + (r - 1)*LIS_rc%obs_lnc(k))
               lat2d(c,r) = LIS_obs_domain(n, k)%lat(c + (r - 1)*LIS_rc%obs_lnc(k))

               !find dx and dy
               if (r==1.and.c==1) then
                  lon2d(c+1,r) = LIS_obs_domain(n, k)%lon((c + 1) + (r - 1)*LIS_rc%obs_lnc(k))
                  lat2d(c,r+1) = LIS_obs_domain(n, k)%lat(c + ((r + 1) - 1)*LIS_rc%obs_lnc(k))
                  dx = lon2d(c+1,r) - lon2d(c,r)
                  dy = lat2d(c,r+1) - lat2d(c,r)
               endif

               if (lon2d(c,r)+dx/2 >= min_lon_dgg.and.&
                   lon2d(c,r)-dx/2 <= max_lon_dgg.and.&
                   lat2d(c,r)+dy/2 >= min_lat_dgg.and.&
                   lat2d(c,r)-dy/2 <= max_lat_dgg) then
                  
                  if (SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_assign) then
                     if (size(SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices) > 0) then
                        do i=1,size(SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices)
                           if (SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices(i) <= max_dgg_id_number) then
                              obstime_dgg(SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices(i)) &
                                   = obstime_dgg(SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices(i))/3600


                              if (int(obstime_dgg(SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices(i))) &
                                                                                              == chr) then
                                 ! soil moisture uncertainty < 0.07 m3m-3
                                 if (smunct_dgg(SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices(i)) &
                                                                                            < 0.07) then
                                    if (sm_dgg(SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices(i)) >= 0) then
                                       smobs_sum(c,r) = smobs_sum(c,r) + &
                                                        sm_dgg(SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices(i))
                                       count_smobs(c,r) = count_smobs(c,r) + 1
                                    endif
                                 endif
                              endif
                           endif
                        enddo
                     endif
                  endif
               endif

            enddo
         enddo

         do r=1,LIS_rc%obs_lnr(k)
            do c=1,LIS_rc%obs_lnc(k)
               lon2d(c,r) = LIS_obs_domain(n, k)%lon(c + (r - 1)*LIS_rc%obs_lnc(k))
               lat2d(c,r) = LIS_obs_domain(n, k)%lat(c + (r - 1)*LIS_rc%obs_lnc(k))

               if (lon2d(c,r)+dx/2 >= min_lon_dgg.and.&
                   lon2d(c,r)-dx/2 <= max_lon_dgg.and.&
                   lat2d(c,r)+dy/2 >= min_lat_dgg.and.&
                   lat2d(c,r)-dy/2 <= max_lat_dgg) then

                  if (count_smobs(c,r) > 0) then
                     smobs_inp(c,r) = smobs_sum(c,r)/count_smobs(c,r)

                     SMOSNRTNNL2sm_struc(n)%smtime(c,r) = time


                  else ! count_smobs(c,r) == 0
                     ! find neighbor grids that have soil moisture values
                     if (r > 1) then
                        if (count_smobs(c,r-1) > 0) then
                           smobs_sum(c,r) = smobs_sum(c,r) + smobs_sum(c,r-1)/count_smobs(c,r-1)
                           count_smobs(c,r) = count_smobs(c,r) + 1
                        endif
                     endif
                     if (c > 1) then
                        if (count_smobs(c-1,r) > 0) then
                           smobs_sum(c,r) = smobs_sum(c,r) + smobs_sum(c-1,r)/count_smobs(c-1,r)
                           count_smobs(c,r) = count_smobs(c,r) + 1
                        endif
                     endif
                     if (c < LIS_rc%obs_lnc(k)) then
                        if (count_smobs(c+1,r) > 0) then
                           smobs_sum(c,r) = smobs_sum(c,r) + smobs_sum(c+1,r)/count_smobs(c+1,r)
                           count_smobs(c,r) = count_smobs(c,r) + 1
                        endif
                     endif
                     if (r < LIS_rc%obs_lnr(k)) then
                        if (count_smobs(c,r+1) > 0) then
                           smobs_sum(c,r) = smobs_sum(c,r) + smobs_sum(c,r+1)/count_smobs(c,r+1)
                           count_smobs(c,r) = count_smobs(c,r) + 1
                        endif
                     endif

                     ! get average smobs only if the number of neighbor grids having soil moisutre values
                     ! is 3 or greater
                     if (count_smobs(c,r) >= 3) then
                        smobs_inp(c,r) = smobs_sum(c,r)/count_smobs(c,r)

                        SMOSNRTNNL2sm_struc(n)%smtime(c,r) = time
                     endif
                  endif
               endif
            enddo
         enddo

         deallocate(sm_dgg)
         deallocate(smunct_dgg)
         deallocate(obstime_dgg)

      endif !min_dgg_id_number > 0

      deallocate(sm)
      deallocate(lon_dgg)
      deallocate(lat_dgg)
      deallocate(smunct)
      deallocate(obstime)
      deallocate(DGG_id_number)
      deallocate(i_dgg)
   endif

   ios = nf90_close(ncid=nid)
   call LIS_verify(ios,'Error closing file '//trim(fname))

#endif

end subroutine read_SMOSNRTNNL2sm_data


subroutine find_SMOS_Dgg_id_number(n,c,r, lon_dgg, lat_dgg, &
                                   dgg_length, i_dgg, DGG_id_number, &
                                   min_lon_dgg, max_lon_dgg, &
                                   min_lat_dgg, max_lat_dgg, &
                                   dx, dy, lon2d, lat2d)
! 
! !USES:
  use SMOSNRTNNL2sm_Mod, only : SMOSNRTNNL2sm_struc

  implicit none

! !INPUT PARAMETERS: 
! 
  integer              :: n,c,r
  integer              :: dgg_length
  real                 :: lon_dgg(dgg_length), lat_dgg(dgg_length)
  real                 :: max_lon_dgg, min_lon_dgg, max_lat_dgg, min_lat_dgg
  integer              :: DGG_id_number(dgg_length)
!
! !OUTPUT PARAMETERS:
!
!EOP
  integer              :: i_dgg(dgg_length)
  real                 :: lat2d, lon2d
  real                 :: dx, dy
  integer, allocatable :: indices(:)

  indices = pack([i_dgg], lon_dgg >= lon2d-dx/2.and.&
                          lon_dgg <= lon2d+dx/2.and.&
                          lat_dgg >= lat2d-dy/2.and.&
                          lat_dgg <= lat2d+dy/2)

  if (size(indices) > 0) then
     allocate(SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices(size(indices)))

     SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_indices = DGG_id_number(indices)
     SMOSNRTNNL2sm_struc(n)%SMOS_lookup(c,r)%dgg_assign = .true.

  endif

  deallocate(indices)

end subroutine find_SMOS_Dgg_id_number
