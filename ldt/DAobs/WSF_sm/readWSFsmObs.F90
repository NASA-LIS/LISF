!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: readWSFsmObs
! \label{readWSFsmObs}
!
! !REVISION HISTORY:
!   2025: Initial Specification
!
! !INTERFACE:
subroutine readWSFsmObs(n)
! !USES:
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_DAobsDataMod
  use WSFsm_obsMod

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
!
! !DESCRIPTION:
!
!  This subroutine provides the data reader for
!  WSF soil moisture data for CDF generation.
!  It searches for all available hourly NetCDF files
!  for the current day and composites them.
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c, r
  character*200     :: fname
  character*200     :: list_files
  character*8       :: yyyymmdd
  character*2       :: hh
  real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer           :: ftn, ierr, rc
  integer, external :: create_filelist  ! C function

!-----------------------------------------------------------------------
!  It is assumed that CDF is computed using daily observations.
!  Read all available hourly files for the current day.
!-----------------------------------------------------------------------
  WSFsmobs(n)%smobs = LDT_rc%udef
  smobs = LDT_rc%udef

  write(yyyymmdd,'(i4.4,2i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da

  ! Search for all hourly files for this day
  ! File pattern: ARFS_SM_*YYYYMMDDT*.nc
  list_files = trim(WSFsmobs(n)%odir)//'/ARFS_SM_*' &
       //trim(yyyymmdd)//'T*.nc'

  write(LDT_logunit,*) '[INFO] Searching for ',trim(list_files)
  rc = create_filelist(trim(list_files)//char(0), &
       "WSF_filelist_ldt.sm.dat"//char(0))

  if (rc .ne. 0) then
     write(LDT_logunit,*) &
          '[WARN] Problem encountered when searching for WSF files'
     write(LDT_logunit,*) &
          'Was searching for ',trim(list_files)
     write(LDT_logunit,*) &
          'LDT will continue...'
  endif

  ftn = LDT_getNextUnitNumber()
  open(ftn, file="./WSF_filelist_ldt.sm.dat", &
       status='old', iostat=ierr)

  do while(ierr.eq.0)
     read(ftn,'(a)',iostat=ierr) fname
     if(ierr.ne.0) then
        exit
     endif

     inquire(file=trim(fname), exist=file_exists)
     if(file_exists) then
        write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
        call read_WSFsm_data_ldt(n, fname, smobs)
     endif
  enddo
  call LDT_releaseUnitNumber(ftn)

  ! Transfer interpolated 1D data to 2D observation array
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(smobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then
           WSFsmobs(n)%smobs(c,r) = smobs(c+(r-1)*LDT_rc%lnc(n))
        endif
     enddo
  enddo

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%soilmoist_obs,&
     WSFsmobs(n)%smobs,vlevel=1)

end subroutine readWSFsmObs


!BOP
!
! !ROUTINE: read_WSFsm_data_ldt
! \label{read_WSFsm_data_ldt}
!
! !INTERFACE:
subroutine read_WSFsm_data_ldt(n, fname, smobs_ip)
!
! !USES:
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_coreMod
  use LDT_logMod
  use WSFsm_obsMod

  implicit none
!
! !INPUT PARAMETERS:
!
  integer                  :: n
  character (len=*)        :: fname
  real                     :: smobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
!
! !DESCRIPTION:
!   This subroutine reads the WSF soil moisture NetCDF file and
!   interpolates to the LDT running domain.
!
!   The arguments are:
!   \begin{description}
!    \item[n]    index of the nest
!    \item[fname] name of the WSF SM NetCDF file
!    \item[smobs\_ip] SM data processed to the LDT domain
!   \end{description}
!
!EOP

  real                    :: sm_in(WSFsmobs(n)%nc*WSFsmobs(n)%nr)
  logical*1               :: sm_data_b(WSFsmobs(n)%nc*WSFsmobs(n)%nr)
  logical*1               :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer                 :: c, r
  integer                 :: ios
  ! NetCDF variables
  logical                 :: file_exists
  character(255)          :: map_projection
  integer                 :: ncid, dim_ids(3), var_id
  integer                 :: ntime, nlat, nlon
  real, allocatable       :: tmp(:,:,:)
  integer                 :: rc

  sm_in = LDT_rc%udef
  sm_data_b = .false.

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  ! Check if file exists
  inquire(file=trim(fname), exist=file_exists)
  if (.not. file_exists) then
     write(LDT_logunit,*) '[ERR] Cannot find ', trim(fname)
     return
  end if

  ! Open the NetCDF file
  rc = nf90_open(path=trim(fname), &
       mode=NF90_NOWRITE, &
       ncid=ncid)
  if (rc .ne. 0) then
     write(LDT_logunit,*) '[ERR] Cannot open ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     return
  end if

  ! Read the map projection
  rc = nf90_get_att(ncid=ncid, &
       varid=NF90_GLOBAL, &
       name='MAP_PROJECTION', &
       values=map_projection)
  if (rc .ne. 0) then
     write(LDT_logunit,*) &
          '[ERR] Cannot read MAP_PROJECTION from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Sanity check map projection
  if (trim(map_projection) .ne. 'EQUIDISTANT CYLINDRICAL') then
     write(LDT_logunit,*) &
          '[ERR] Unrecognized map projection found in WSF file!'
     write(LDT_logunit,*) '[ERR] Expected EQUIDISTANT CYLINDRICAL'
     write(LDT_logunit,*) '[ERR] Found ',trim(map_projection)
     write(LDT_logunit,*) '[ERR] LDT will skip file ', trim(fname)
     rc = nf90_close(ncid)
     return
  end if

  ! Get dimension IDs
  rc = nf90_inq_dimid(ncid=ncid, name='time', dimid=dim_ids(3))
  if (rc .ne. 0) then
     write(LDT_logunit,*) &
          '[ERR] Cannot read time dimension from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inq_dimid(ncid=ncid, name='lat', dimid=dim_ids(2))
  if (rc .ne. 0) then
     write(LDT_logunit,*) &
          '[ERR] Cannot read lat dimension from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inq_dimid(ncid=ncid, name='lon', dimid=dim_ids(1))
  if (rc .ne. 0) then
     write(LDT_logunit,*) &
          '[ERR] Cannot read lon dimension from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Get actual dimension sizes
  rc = nf90_inquire_dimension(ncid=ncid, &
       dimid=dim_ids(3), len=ntime)
  if (rc .ne. 0) then
     write(LDT_logunit,*) &
          '[ERR] Cannot read time dimension from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inquire_dimension(ncid=ncid, &
       dimid=dim_ids(2), len=nlat)
  if (rc .ne. 0) then
     write(LDT_logunit,*) &
          '[ERR] Cannot read lat dimension from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inquire_dimension(ncid=ncid, &
       dimid=dim_ids(1), len=nlon)
  if (rc .ne. 0) then
     write(LDT_logunit,*) &
          '[ERR] Cannot read lon dimension from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Sanity check the dimensions
  if (ntime .ne. 1) then
     write(LDT_logunit,*) '[ERR] Expected time dimension to be 1'
     write(LDT_logunit,*) '[ERR] Found ', ntime, ' from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  if (nlat .ne. WSFsmobs(n)%nr) then
     write(LDT_logunit,*) '[ERR] Expected lat dimension to be ', &
          WSFsmobs(n)%nr
     write(LDT_logunit,*) '[ERR] Found ', nlat, ' from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  if (nlon .ne. WSFsmobs(n)%nc) then
     write(LDT_logunit,*) '[ERR] Expected lon dimension to be ', &
          WSFsmobs(n)%nc
     write(LDT_logunit,*) '[ERR] Found ', nlon, ' from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Fetch the variable id
  rc = nf90_inq_varid(ncid=ncid, name='arfs_sm', varid=var_id)
  if (rc .ne. 0) then
     write(LDT_logunit,*) &
          '[ERR] Cannot read arfs_sm from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Read the retrievals
  allocate(tmp(nlon, nlat, ntime))
  rc = nf90_get_var(ncid=ncid, varid=var_id, values=tmp)
  if (rc .ne. 0) then
     write(LDT_logunit,*) &
          '[ERR] Cannot read arfs_sm from ', trim(fname)
     write(LDT_logunit,*) '[ERR] LDT will continue...'
     rc = nf90_close(ncid)
     deallocate(tmp)
     return
  end if
  rc = nf90_close(ncid)

! Apply QC: valid range 0.0 to 1.0 (guard against NaN)
  do r = 1, nlat
     do c = 1, nlon
        if (tmp(c,r,1) .eq. tmp(c,r,1)) then  ! NaN /= NaN, so this filters NaNs
           if (tmp(c,r,1) >= 0.0 .and. &
                tmp(c,r,1) <= 1.0) then
              sm_in(c + (r-1)*WSFsmobs(n)%nc) = tmp(c,r,1)
              sm_data_b(c + (r-1)*WSFsmobs(n)%nc) = .true.
           end if
        end if
     end do
  end do
  deallocate(tmp)

#endif

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!--------------------------------------------------------------------------
  if(LDT_isLDTatAfinerResolution(n,0.093701)) then
     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          sm_data_b, sm_in, smobs_b_ip, smobs_ip, &
          WSFsmobs(n)%nc*WSFsmobs(n)%nr, &
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          WSFsmobs(n)%w11, WSFsmobs(n)%w12,&
          WSFsmobs(n)%w21, WSFsmobs(n)%w22,&
          WSFsmobs(n)%n11, WSFsmobs(n)%n12,&
          WSFsmobs(n)%n21, WSFsmobs(n)%n22,LDT_rc%udef,ios)
  else
     call upscaleByAveraging(WSFsmobs(n)%nc*WSFsmobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_rc%udef, WSFsmobs(n)%n11,&
          sm_data_b, sm_in, smobs_b_ip, smobs_ip)
  endif

end subroutine read_WSFsm_data_ldt
