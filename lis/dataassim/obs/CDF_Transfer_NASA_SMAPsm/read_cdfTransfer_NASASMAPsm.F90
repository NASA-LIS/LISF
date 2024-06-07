!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_cdfTransfer_NASASMAPsm
! \label{read_cdfTransfer_NASASMAPsm}
!
! !REVISION HISTORY:
!  2 MAR 2022: Mahdi Navari: initial specification based on read_NASASMAPsm.F90
!
! !INTERFACE:
subroutine read_cdfTransfer_NASASMAPsm(n, k, OBS_State, OBS_Pert_State)
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
   use cdfTransfer_NASASMAPsm_Mod, only: cdfT_SMAPsm_struc
   !use cdfTransfer_NASASMAPsm_Mod, only: cdfT_SMAPsm_struc

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
   real, parameter       :: MAX_SM_VALUE = 0.45, MIN_SM_VALUE = 0.0001
   integer                :: status
   integer                :: grid_index
   character*100          :: smobsdir
   character*100          :: fname
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
   character*200          :: list_files
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

   alarmCheck = LIS_isAlarmRinging(LIS_rc, "NASASMAP read alarm")

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

   if (alarmCheck .or. cdfT_SMAPsm_struc(n)%startMode) then
      cdfT_SMAPsm_struc(n)%startMode = .false.
      if ( (cdfT_SMAPsm_struc(n)%data_designation.eq."SPL2SMP_E") .or. &
           (cdfT_SMAPsm_struc(n)%data_designation.eq."SPL2SMP") ) then

         cdfT_SMAPsm_struc(n)%smobs = LIS_rc%udef
         cdfT_SMAPsm_struc(n)%smtime = -1.0
 
         write(yyyymmdd,'(i4.4,2i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
         write(yyyy,'(i4.4)') LIS_rc%yr
         write(mm,'(i2.2)') LIS_rc%mo
         write(dd,'(i2.2)') LIS_rc%da
         write(hh,'(i2.2)') LIS_rc%hr

         if(LIS_masterproc) then
            list_files = trim(smobsdir)//'/'//trim(yyyy)//'.'//trim(mm)//'.'//&
                         trim(dd)//'/SMAP_L2_*' &
                         //trim(yyyy)//trim(mm)//trim(dd)//'T'//trim(hh)//'*.h5'
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

            call read_SMAPL2sm_data_cdfTransfer(n,k,fname,&
                 cdfT_SMAPsm_struc(n)%smobs,timenow)
         enddo
         call LIS_releaseUnitNumber(ftn)

      elseif (cdfT_SMAPsm_struc(n)%data_designation .eq. "SPL3SMP_E") then
!---------------------------------------------------------------------------
! MN: create filename for 9 km product
!---------------------------------------------------------------------------

         write (yyyy, '(i4.4)') LIS_rc%yr
         write (mm, '(i2.2)') LIS_rc%mo
         write (dd, '(i2.2)') LIS_rc%da
         write (CRID, '(a)') cdfT_SMAPsm_struc(n)%release_number

         ! EMK...Make sure only one PET calls the file system to determine what
         ! SMAP files are available.
         if (LIS_masterproc) then

            list_files = trim(smobsdir)//'/'//trim(yyyy)//'.'//trim(mm)//'.'// &
                         trim(dd)//'/SMAP_L3_SM_P_E_' &
                         //trim(yyyy)//trim(mm)//trim(dd)//'_'// &
                         trim(CRID)//'*.h5'
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
         end if
#if (defined SPMD)
         call mpi_barrier(lis_mpi_comm, ierr)
#endif

         ftn = LIS_getNextUnitNumber()
         open (ftn, file="./SMAP_filelist.dat", &
               action='read', status='old', iostat=ierr)

! if multiple files for the same time and orbits are present, the latest
! one will overwrite older ones, though multiple (redundant) reads occur.
! This assumes that the 'ls command' will list the files in that order.

         do while (ierr .eq. 0)
            read (ftn, '(a)', iostat=ierr) fname
            if (ierr .ne. 0) then
               exit
            endif

            write (LIS_logunit, *) '[INFO] Reading descending pass ', trim(fname)
            call read_NASASMAP_E_data_cdfTransfer(n, k, 'D', fname, smobs_D)

            write (LIS_logunit, *) '[INFO] Reading ascending pass ', trim(fname)
            call read_NASASMAP_E_data_cdfTransfer(n, k, 'A', fname, smobs_A)
         enddo

         cdfT_SMAPsm_struc(n)%smobs = LIS_rc%udef
         cdfT_SMAPsm_struc(n)%smtime = -1
         call LIS_releaseUnitNumber(ftn)
!-------------------------------------------------------------------------
!   Ascending pass assumed to be at 6pm localtime and the descending
!   pass is assumed to be at 6am local time
!-------------------------------------------------------------------------
         do r = 1, LIS_rc%obs_lnr(k)
            do c = 1, LIS_rc%obs_lnc(k)
               grid_index = LIS_obs_domain(n, k)%gindex(c, r)
               if (grid_index .ne. -1) then

                  if (smobs_D(c + (r - 1)*LIS_rc%obs_lnc(k)) .ne. -9999.0) then
                     cdfT_SMAPsm_struc(n)%smobs(c, r) = &
                        smobs_D(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lon = LIS_obs_domain(n, k)%lon(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lhour = 6.0
                     call LIS_localtime2gmt(gmt, lon, lhour, zone)
                     cdfT_SMAPsm_struc(n)%smtime(c, r) = gmt

                  endif
!-------------------------------------------------------------------------
! The ascending data is used only over locations where descending data
! doesn't exist.
!-------------------------------------------------------------------------
                  if (smobs_A(c + (r - 1)*LIS_rc%obs_lnc(k)) .ne. -9999.0 .and. &
                      cdfT_SMAPsm_struc(n)%smobs(c, r) .eq. -9999.0) then
                     cdfT_SMAPsm_struc(n)%smobs(c, r) = &
                        smobs_A(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lon = LIS_obs_domain(n, k)%lon(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lhour = 18.0
                     call LIS_localtime2gmt(gmt, lon, lhour, zone)
                     cdfT_SMAPsm_struc(n)%smtime(c, r) = gmt
                  endif
               endif
            enddo
         enddo

      elseif (cdfT_SMAPsm_struc(n)%data_designation .eq. "SPL3SMP") then
!-------------------------------------------------------------------------
! MN: create filename for 36km product  (SMAP_L3_SM_P_)
!-------------------------------------------------------------------------

         write (yyyy, '(i4.4)') LIS_rc%yr
         write (mm, '(i2.2)') LIS_rc%mo
         write (dd, '(i2.2)') LIS_rc%da
         write (CRID, '(a)') cdfT_SMAPsm_struc(n)%release_number

         ! EMK...Make sure only one PET calls the file system to determine what
         ! SMAP files are available.
         if (LIS_masterproc) then
            list_files = trim(smobsdir)//'/'//trim(yyyy)//'.'//trim(mm)//'.'// &
                         trim(dd)//'/SMAP_L3_SM_P_' &
                         //trim(yyyy)//trim(mm)//trim(dd)//'_'// &
                         trim(CRID)//'*.h5'
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
         end if
#if (defined SPMD)
         call mpi_barrier(lis_mpi_comm, ierr)
#endif

         ftn = LIS_getNextUnitNumber()
         open (ftn, file="./SMAP_filelist.dat", &
               action='read', status='old', iostat=ierr)

! if multiple files for the same time and orbits are present, the latest
! one will overwrite older ones, though multiple (redundant) reads occur.
! This assumes that the 'ls command' will list the files in that order.

         do while (ierr .eq. 0)
            read (ftn, '(a)', iostat=ierr) fname
            if (ierr .ne. 0) then
               exit
            endif

            write (LIS_logunit, *) '[INFO] Reading descending pass ', trim(fname)
            call read_NASASMAP_data_cdfTransfer(n, k, 'D', fname, smobs_D)

            write (LIS_logunit, *) '[INFO] Reading ascending pass ', trim(fname)
            call read_NASASMAP_data_cdfTransfer(n, k, 'A', fname, smobs_A)
         enddo

         cdfT_SMAPsm_struc(n)%smobs = LIS_rc%udef
         cdfT_SMAPsm_struc(n)%smtime = -1
         call LIS_releaseUnitNumber(ftn)
!-------------------------------------------------------------------------
!   Ascending pass assumed to be at 6pm localtime and the descending
!   pass is assumed to be at 6am local time
!-------------------------------------------------------------------------
         do r = 1, LIS_rc%obs_lnr(k)
            do c = 1, LIS_rc%obs_lnc(k)
               grid_index = LIS_obs_domain(n, k)%gindex(c, r)
               if (grid_index .ne. -1) then

                  if (smobs_D(c + (r - 1)*LIS_rc%obs_lnc(k)) .ne. -9999.0) then
                     cdfT_SMAPsm_struc(n)%smobs(c, r) = &
                        smobs_D(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lon = LIS_obs_domain(n, k)%lon(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lhour = 6.0
                     call LIS_localtime2gmt(gmt, lon, lhour, zone)
                     cdfT_SMAPsm_struc(n)%smtime(c, r) = gmt

                  endif
!-------------------------------------------------------------------------
! The ascending data is used only over locations where descending data
! doesn't exist.
!-------------------------------------------------------------------------
                  if (smobs_A(c + (r - 1)*LIS_rc%obs_lnc(k)) .ne. -9999.0 .and. &
                      cdfT_SMAPsm_struc(n)%smobs(c, r) .eq. -9999.0) then
                     cdfT_SMAPsm_struc(n)%smobs(c, r) = &
                        smobs_A(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lon = LIS_obs_domain(n, k)%lon(c + (r - 1)*LIS_rc%obs_lnc(k))
                     lhour = 18.0
                     call LIS_localtime2gmt(gmt, lon, lhour, zone)
                     cdfT_SMAPsm_struc(n)%smtime(c, r) = gmt
                  endif
               endif
            enddo
         enddo

      endif ! sensor
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

   if ( (cdfT_SMAPsm_struc(n)%data_designation.eq."SPL2SMP_E") .or. &
        (cdfT_SMAPsm_struc(n)%data_designation.eq."SPL2SMP") ) then

      do r=1,LIS_rc%obs_lnr(k)
         do c=1,LIS_rc%obs_lnc(k)
            if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
               grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
               dt = (cdfT_SMAPsm_struc(n)%smtime(c,r)-time1)
               if(dt.ge.0.and.dt.lt.(time3-time1)) then
                  sm_current(c,r) = &
                       cdfT_SMAPsm_struc(n)%smobs(c,r)
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

   else

      do r = 1, LIS_rc%obs_lnr(k)
         do c = 1, LIS_rc%obs_lnc(k)
            if (LIS_obs_domain(n, k)%gindex(c, r) .ne. -1) then
               grid_index = c + (r - 1)*LIS_rc%obs_lnc(k)

               dt = (LIS_rc%gmt - cdfT_SMAPsm_struc(n)%smtime(c, r))*3600.0
               if (dt .ge. 0 .and. dt .lt. LIS_rc%ts) then
                  sm_current(c, r) = &
                        cdfT_SMAPsm_struc(n)%smobs(c, r)
                  if (LIS_obs_domain(n, k)%gindex(c, r) .ne. -1) then
                     obs_unsc(LIS_obs_domain(n, k)%gindex(c, r)) = &
                        sm_current(c, r)
                  endif
                  if (sm_current(c, r) .ne. LIS_rc%udef) then
                     fnd = 1
                  endif
               endif
            endif
         enddo
      enddo
   endif


!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a geolocation independent
!  CDF-scaling approach
!-------------------------------------------------------------------------
   if (cdfT_SMAPsm_struc(n)%useCDFtransfer .gt. 0) then
      if (LIS_rc%da .eq. 1 .and. LIS_rc%hr .eq. 0 .and. &
          LIS_rc%mn .eq. 0 .and. LIS_rc%ss .eq. 0) then
         call read_CDFtransferdata_all(n,k,&
              cdfT_SMAPsm_struc(n)%nbins,&
              cdfT_SMAPsm_struc(n)%ntimes,&
              cdfT_SMAPsm_struc(n)%n_strat_bins, &
              cdfT_SMAPsm_struc(n)%modelcdffile, &
              "SoilMoist",&
              cdfT_SMAPsm_struc(n)%model_xrange,&
              cdfT_SMAPsm_struc(n)%model_cdf)

         call read_CDFtransferdata_all(n,k,&
              cdfT_SMAPsm_struc(n)%nbins,&
              cdfT_SMAPsm_struc(n)%ntimes,&
              cdfT_SMAPsm_struc(n)%n_strat_bins, &
              cdfT_SMAPsm_struc(n)%obscdffile, &
              "SoilMoist",&
              cdfT_SMAPsm_struc(n)%obs_xrange,&
              cdfT_SMAPsm_struc(n)%obs_cdf)
      endif
      if (fnd .ne. 0) then

         call LIS_rescale_with_stratified_CDF(&
              n,             &
              k,             &
              cdfT_SMAPsm_struc(n)%nbins,         &
              cdfT_SMAPsm_struc(n)%ntimes,        &
              MAX_SM_VALUE, &
              MIN_SM_VALUE, &
              cdfT_SMAPsm_struc(n)%model_xrange,  &
              cdfT_SMAPsm_struc(n)%obs_xrange,    &
              cdfT_SMAPsm_struc(n)%model_cdf,     &
              cdfT_SMAPsm_struc(n)%obs_cdf,       &
              cdfT_SMAPsm_struc(n)%ref_p_climo_maxval, &
              cdfT_SMAPsm_struc(n)%target_p_climo,&
              cdfT_SMAPsm_struc(n)%n_strat_bins,&
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
                        //trim(LIS_CDFTRANSFERNASASMAPsmobsId)//char(0), n, k, OBS_state)

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

      if(cdfT_SMAPsm_struc(n)%useSsdevScal.eq.1) then
         write(LIS_logunit,*) '[ERR] "use scaled standard deviation model" does not work '
         write(LIS_logunit,*) '[ERR] for CDF transfer method'
         call LIS_endrun()
      endif
   else
      call ESMF_AttributeSet(OBS_State, "Data Update Status", &
                             .false., rc=status)
      call LIS_verify(status)
   endif

end subroutine read_cdfTransfer_NASASMAPsm

!BOP
!
! !ROUTINE: read_SMAPL2sm_data_cdfTransfer
! \label{read_SMAPL2sm_data_cdfTransfer}
!
! !INTERFACE:
subroutine read_SMAPL2sm_data_cdfTransfer(n, k,fname, smobs_inp, time)
!
! !USES:

  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use cdfTransfer_NASASMAPsm_Mod, only : cdfT_SMAPsm_struc

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
  real                     :: smobs_inp(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real*8                   :: time

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION:
!
!
!EOP

#if (defined USE_HDF5)

  character*100,    parameter    :: sm_gr_name = "Soil_Moisture_Retrieval_Data"
  character*100,    parameter    :: sm_field_name = "soil_moisture"
  character*100,    parameter    :: sm_qa_name = "retrieval_qual_flag"
!YK
  character*100,    parameter    :: vwc_field_name = "vegetation_water_content"

  integer(hsize_t), dimension(1) :: dims
  integer(hsize_t), dimension(1) :: maxdims
  integer(hid_t)                 :: file_id
  integer(hid_t)                 :: dspace_id
  integer(hid_t)                 :: row_id, col_id
  integer(hid_t)                 :: sm_gr_id,sm_field_id, sm_qa_id
  integer(hid_t)                 :: sm_gr_id_A,sm_field_id_A
  integer(hid_t)                 :: vwc_field_id ! YK
  real,             allocatable  :: sm_field(:)
  real,             allocatable  :: vwc_field(:)! YK
  integer,          allocatable  :: sm_qa(:)
  integer,          allocatable  :: ease_row(:)
  integer,          allocatable  :: ease_col(:)
  integer                        :: c,r,t
  logical*1                      :: sm_data_b(cdfT_SMAPsm_struc(n)%nc*cdfT_SMAPsm_struc(n)%nr)
  logical*1                      :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: sm_data(cdfT_SMAPsm_struc(n)%nc*cdfT_SMAPsm_struc(n)%nr)
  real                           :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  integer                        :: status,ios,iret

  call h5open_f(status)
  call LIS_verify(status, 'Error opening HDF fortran interface')

  call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F, file_id, status)
  call LIS_verify(status, 'Error opening SMAP L2 file ')

  call h5gopen_f(file_id,sm_gr_name,sm_gr_id, status)
  call LIS_verify(status, 'Error opening SM group in SMAP L2 file')

  call h5dopen_f(sm_gr_id,sm_field_name,sm_field_id, status)
  call LIS_verify(status, 'Error opening SM field in SMAP L2 file')

  call h5dopen_f(sm_gr_id,"EASE_row_index",row_id, status)
  call LIS_verify(status, 'Error opening row index field in SMAP L2 file')

  call h5dopen_f(sm_gr_id,"EASE_column_index",col_id, status)
  call LIS_verify(status, 'Error opening column index field in SMAP L2 file')

!YK
  call h5dopen_f(sm_gr_id, sm_qa_name,sm_qa_id, status)
  call LIS_verify(status, 'Error opening QA field in SMAP L2 file')

!YK
  call h5dopen_f(sm_gr_id, vwc_field_name,vwc_field_id, status)
  call LIS_verify(status, 'Error opening Veg water content field in SMAP L2 file')

  call h5dget_space_f(sm_field_id, dspace_id, status)
  call LIS_verify(status, 'Error in h5dget_space_f: reaSMAP L2Obs')

! Size of the arrays
! This routine returns -1 on failure, rank on success.
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, status)
  if(status.eq.-1) then
     call LIS_verify(status, 'Error in h5sget_simple_extent_dims_f: readSMAP L2Obs')
  endif

  allocate(sm_field(maxdims(1)))
  allocate(sm_qa(maxdims(1)))    !YK
  allocate(vwc_field(maxdims(1)))    !YK
  allocate(ease_row(maxdims(1)))
  allocate(ease_col(maxdims(1)))

  call h5dread_f(row_id, H5T_NATIVE_INTEGER,ease_row,dims,status)
  call LIS_verify(status, 'Error extracting row index from SMAP L2 file')

  call h5dread_f(col_id, H5T_NATIVE_INTEGER,ease_col,dims,status)
  call LIS_verify(status, 'Error extracting col index from SMAP L2 file')

  call h5dread_f(sm_field_id, H5T_NATIVE_REAL,sm_field,dims,status)
  call LIS_verify(status, 'Error extracting SM field from SMAP L2 file')

!YK
  call h5dread_f(sm_qa_id, H5T_NATIVE_INTEGER,sm_qa,dims,status)
  call LIS_verify(status, 'Error extracting SM field from SMAP L2 file')

!YK get the vegetation water content
  call h5dread_f(vwc_field_id, H5T_NATIVE_REAL,vwc_field,dims,status)
  call LIS_verify(status, 'Error extracting Veg water content (AM) field from SMAP L2 file')

!YK
  call h5dclose_f(sm_qa_id,status)
  call LIS_verify(status,'Error in H5DCLOSE call')

!YK
  call h5dclose_f(vwc_field_id,status)
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

!--------------------------------------------------------------------YK
!grid the data in EASE projection
! The retrieval_quality_field variable's binary representation consists of bits
! that indicate whether retrieval is performed or not at a given grid cell. 
! When retrieval is performed, it contains additional bits to further 
! indicate the exit status and quality of the retrieval. The first bit 
! indicates the recommended quality (0-means retrieval has recommended quality).
  do t=1,maxdims(1)
     if(ease_col(t).gt.0.and.ease_row(t).gt.0) then
        sm_data(ease_col(t) + &
             (ease_row(t)-1)*cdfT_SMAPsm_struc(n)%nc) = sm_field(t)

        if(vwc_field(t).gt.5) then !YK Aply QC : if VWC > 5 kg/m2
           sm_data(ease_col(t) + &
                (ease_row(t)-1)*cdfT_SMAPsm_struc(n)%nc) = LIS_rc%udef
        else
           if(sm_data(ease_col(t) + &
                   (ease_row(t)-1)*cdfT_SMAPsm_struc(n)%nc).ne.-9999.0) then
              if(ibits(sm_qa(t),0,1).eq.0) then
                 sm_data_b(ease_col(t) + &
                    (ease_row(t)-1)*cdfT_SMAPsm_struc(n)%nc) = .true.
              else
                 sm_data(ease_col(t) + &
                      (ease_row(t)-1)*cdfT_SMAPsm_struc(n)%nc) = LIS_rc%udef
              endif
           endif
        endif
     endif
  enddo
!-----------------------------------------------------------------------

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!--------------------------------------------------------------------------
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:), sm_data_b, sm_data, &
       smobs_b_ip, smobs_ip, &
       cdfT_SMAPsm_struc(n)%nc*cdfT_SMAPsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       cdfT_SMAPsm_struc(n)%rlat, cdfT_SMAPsm_struc(n)%rlon,&
       cdfT_SMAPsm_struc(n)%n11, LIS_rc%udef, ios)


  deallocate(sm_field)
  deallocate(ease_row)
  deallocate(ease_col)

!overwrite the input data
  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(smobs_ip(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then
           smobs_inp(c,r) = &
                smobs_ip(c+(r-1)*LIS_rc%obs_lnc(k))

           cdfT_SMAPsm_struc(n)%smtime(c,r) = &
                time
        endif
     enddo
  enddo
#endif

end subroutine read_SMAPL2sm_data_cdfTransfer

!BOP
!
! !ROUTINE: read_NASASMAP_E_data_cdfTransfer
! \label{read_NASASMAP_E_data_cdfTransfer}
!
! !INTERFACE:
subroutine read_NASASMAP_E_data_cdfTransfer(n, k, pass, fname, smobs_ip)
!
! !USES:

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use cdfTransfer_NASASMAPsm_Mod, only : cdfT_SMAPsm_struc
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
  logical*1                      :: sm_data_b(cdfT_SMAPsm_struc(n)%nc*cdfT_SMAPsm_struc(n)%nr)
  logical*1                      :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: sm_data(cdfT_SMAPsm_struc(n)%nc*cdfT_SMAPsm_struc(n)%nr)
  integer                        :: status,ios

  dimsm      = (/cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr/)
  count_file = (/cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr/)
  count_mem  = (/cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr/)

  allocate(sm_field(cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr))
  allocate(sm_qa(cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr))
  allocate(vwc_field(cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr))
  allocate(dims(2))

  dims(1) = cdfT_SMAPsm_struc(n)%nc
  dims(2) = cdfT_SMAPsm_struc(n)%nr

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
     call LIS_verify(status, 'Error in h5create_simple_f; read_cdfTransfer_NASASMAPsm')

     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_cdfTransfer_NASASMAPsm')

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
     call LIS_verify(status, 'Error in h5create_simple_f; read_cdfTransfer_NASASMAPsm')

     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_cdfTransfer_NASASMAPsm')

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

  sm_data_b = .false.
  t = 1

! The retrieval_quality_field variable's binary representation consists of bits
! that indicate whether retrieval is performed or not at a given grid cell.
! When retrieval is performed, it contains additional bits to further
! indicate the exit status and quality of the retrieval. The first bit
! indicates the recommended quality (0-means retrieval has recommended quality).

  do r=1,cdfT_SMAPsm_struc(n)%nr
     do c=1,cdfT_SMAPsm_struc(n)%nc
        sm_data(t) = sm_field(c,r)

        if(vwc_field(c,r).gt. 5 ) then !MN Aply QC : if VWC > 5 kg/m2
           sm_data(t) = LIS_rc%udef
	 else

           if(sm_data(t).ne.-9999.0) then
              if(ibits(sm_qa(c,r),0,1).eq.0) then
                 sm_data_b(t) = .true.
              else
                 sm_data(t) = -9999.0
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
       cdfT_SMAPsm_struc(n)%nc*cdfT_SMAPsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       cdfT_SMAPsm_struc(n)%rlat, cdfT_SMAPsm_struc(n)%rlon,&
       cdfT_SMAPsm_struc(n)%n11, LIS_rc%udef, ios)

  deallocate(sm_field)
  deallocate(dims)

#endif

end subroutine read_NASASMAP_E_data_cdfTransfer

! MN: the data structure in both 36 km and 9 km products is the same therefore
!         read_NASASMAP_E_data_cdfTransfer is similar to read_NASASMAP_data_cdfTransfer

!BOP
!
! !ROUTINE: read_NASASMAP_data_cdfTransfer
! \label{read_NASASMAP_data_cdfTransfer}
!
! !INTERFACE:
subroutine read_NASASMAP_data_cdfTransfer(n, k, pass, fname, smobs_ip)
!
! !USES:

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use cdfTransfer_NASASMAPsm_Mod, only : cdfT_SMAPsm_struc
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
  logical*1                      :: sm_data_b(cdfT_SMAPsm_struc(n)%nc*cdfT_SMAPsm_struc(n)%nr)
  logical*1                      :: smobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                           :: sm_data(cdfT_SMAPsm_struc(n)%nc*cdfT_SMAPsm_struc(n)%nr)
  integer                        :: status,ios

  dimsm      = (/cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr/)
  count_file = (/cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr/)
  count_mem  = (/cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr/)

  allocate(sm_field(cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr))
  allocate(sm_qa(cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr))
  allocate(vwc_field(cdfT_SMAPsm_struc(n)%nc, cdfT_SMAPsm_struc(n)%nr))
  allocate(dims(2))

  dims(1) = cdfT_SMAPsm_struc(n)%nc
  dims(2) = cdfT_SMAPsm_struc(n)%nr

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
     call LIS_verify(status, 'Error in h5create_simple_f; read_cdfTransfer_NASASMAPsm')

     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_cdfTransfer_NASASMAPsm')

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
     call LIS_verify(status, 'Error in h5create_simple_f; read_cdfTransfer_NASASMAPsm')

     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
          start=offset_mem, count=count_mem, hdferr=status)
     call LIS_verify(status, 'Error in h5sselect_hyperslab_f: read_cdfTransfer_NASASMAPsm')

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

  sm_data_b = .false.
  t = 1

! The retrieval_quality_field variable's binary representation consists of bits
! that indicate whether retrieval is performed or not at a given grid cell.
! When retrieval is performed, it contains additional bits to further
! indicate the exit status and quality of the retrieval. The first bit
! indicates the recommended quality (0-means retrieval has recommended quality).

  do r=1,cdfT_SMAPsm_struc(n)%nr
     do c=1,cdfT_SMAPsm_struc(n)%nc
        sm_data(t) = sm_field(c,r)

        if(vwc_field(c,r).gt. 5 ) then !MN Aply QC : if VWC > 5 kg/m2
           sm_data(t) = LIS_rc%udef
	 else

           if(sm_data(t).ne.-9999.0) then
              if(ibits(sm_qa(c,r),0,1).eq.0) then
                 sm_data_b(t) = .true.
              else
                 sm_data(t) = -9999.0
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
       cdfT_SMAPsm_struc(n)%nc*cdfT_SMAPsm_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       cdfT_SMAPsm_struc(n)%rlat, cdfT_SMAPsm_struc(n)%rlon,&
       cdfT_SMAPsm_struc(n)%n11, LIS_rc%udef, ios)

  deallocate(sm_field)
  deallocate(dims)

#endif


end subroutine read_NASASMAP_data_cdfTransfer


