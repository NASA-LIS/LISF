!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readSMAPEOPL_SMObs
! \label{readSMAPEOPL_SMObs}
! 
! !REVISION HISTORY: 
! 26 Apr 2023: Mahdi Navari, Initial Specification
! 
! !INTERFACE: 
subroutine readSMAPEOPL_SMObs(n)
! !USES: 
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod,   only : LVT_get_julss, LVT_tick
  use SMAPEOPLSMobsMod,  only : SMAPEOPLsmobs 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for 
! SMAP_E_OPL soil moisture data
!
!EOP

  real*8            :: timenow
  logical           :: alarmCheck
  real              :: smobs(LVT_rc%lnc,LVT_rc%lnr)
  character*200     :: fname
  integer           :: mn_ind
  integer           :: mn, ss
  integer           :: doy
  character*8       :: yyyymmdd
  character*2       :: hh
  character*200     :: list_files
  character*200     :: smap_filename(10)
  integer           :: i
  integer           :: ftn, ierr
  real              :: gmt
  integer           :: rc

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations.
!-----------------------------------------------------------------------
  SMAPEOPLsmobs(n)%smobs = LVT_rc%udef
  smobs = LVT_rc%udef

  if (LVT_rc%ts .gt. 3600) then
     write (LVT_logunit, *) '[ERR] Please set the LVT timestep to 1hr or less'
     write (LVT_logunit, *) '[ERR] This is required for SMAP_E_OPL data processing'
     call LVT_endrun()
  endif

  timenow = float(LVT_rc%dhr(n))*3600 + 60*LVT_rc%dmn(n) + &
       LVT_rc%dss(n)
  alarmcheck = (mod(timenow, 3600.0).eq.0)


  if(alarmCheck) then

  write (yyyymmdd, '(i4.4,2i2.2)') LVT_rc%dyr(n), LVT_rc%dmo(n), LVT_rc%dda(n)
  write (hh, '(i2.2)') LVT_rc%dhr(n)


  list_files = 'ls '//trim(SMAPEOPLsmobs(n)%odir)// &
               '/ARFS_SM_*'//trim(yyyymmdd)//'T'//trim(hh) &
               //"*.nc > SMAPEOPL_filelist.dat"

  call system(trim(list_files))

  i = 1
  ftn = LVT_getNextUnitNumber()
  open (ftn, file="./SMAPEOPL_filelist.dat", &
        status='old', iostat=ierr)

  do while (ierr .eq. 0)
     read (ftn, '(a)', iostat=ierr) fname
     if (ierr .ne. 0) then
        exit
     endif
     mn_ind = index(fname, trim(yyyymmdd)//'T'//trim(hh))

     mn_ind = index(fname, trim(yyyymmdd)//'T'//trim(hh)) + 11
     read (fname(mn_ind:mn_ind + 1), '(i2.2)') mn
     ss = 0
     call LVT_tick(timenow, doy, gmt, LVT_rc%dyr(n), LVT_rc%dmo(n), LVT_rc%dda(n), &
                   LVT_rc%dhr(n), mn, ss, 0)

     smap_filename(i) = fname

     write (LVT_logunit, *) '[INFO] reading ', trim(smap_filename(i))

     call read_SMAPEOPLsm_data(n, smap_filename(i), smobs, timenow)

     i = i + 1
  enddo
  call LVT_releaseUnitNumber(ftn)
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, n,&
       smobs,vlevel=1,units="m3/m3")


end subroutine readSMAPEOPL_SMObs

!BOP
! 
! !ROUTINE: read_SMAPEOPLsm_data
! \label{read_SMAPEOPLsm_data}
!
! !INTERFACE:
subroutine read_SMAPEOPLsm_data(n, fname, smobs_inp, time)
! 
! !USES:

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SMAPEOPLSMobsMod,  only : SMAPEOPLsmobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                  :: n
  character (len=*)        :: fname
  real                     :: smobs_inp(LVT_rc%lnc,LVT_rc%lnr)
  real*8                   :: time
!EOP
  integer,  parameter     :: nc=2560, nr=1920
  real*4                  :: sm_raw(SMAPEOPLsmobs(n)%nc,SMAPEOPLsmobs(n)%nr)
  real                    :: sm_in(SMAPEOPLsmobs(n)%nc*SMAPEOPLsmobs(n)%nr)
  real                    :: smobs_ip(LVT_rc%lnc*LVT_rc%lnr)
  logical*1               :: sm_data_b(SMAPEOPLsmobs(n)%nc*SMAPEOPLsmobs(n)%nr)
  logical*1               :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)
  integer                 :: smid
  integer                 :: ios, nid
  integer                 :: c,r
  integer                 :: ftn1
  logical                 :: file_exists
  character(255)          :: map_projection
  integer                 :: ncid, dim_ids(3), var_id
  integer                 :: ntime, nlat, nlon
  real, allocatable       :: tmp(:,:,:)
  integer                 :: rc

  sm_in = LVT_rc%udef
  sm_data_b = .false.

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ! See if the file exists.
  inquire(file=trim(fname), exist=file_exists)
  if (.not. file_exists) then
     write(LVT_logunit,*)'[ERR] Cannot find ', trim(fname)
     return
  end if

  ! Open the file
  rc = nf90_open(path=trim(fname), &
       mode=NF90_NOWRITE, &
       ncid=ncid)
  if (rc .ne. 0) then
     write(LVT_logunit,*)'[ERR] Cannot open ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     return
  end if

  ! Read the map projection
  rc = nf90_get_att(ncid=ncid, &
       varid=NF90_GLOBAL, &
       name='MAP_PROJECTION', &
       values=map_projection)
  if (rc .ne. 0) then
     write(LVT_logunit,*)'[ERR] Cannot read MAP_PROJECTION from ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Sanity check map projection
  ! TODO:  Support other map projections
  if (trim(map_projection) .ne. 'EQUIDISTANT CYLINDRICAL') then
     write(LVT_logunit,*) &
          '[ERR] Unrecognized map projection found in SMAP file!'
     write(LVT_logunit,*) '[ERR] Expected EQUIDISTANT CYLINDRICAL'
     write(LVT_logunit,*) '[ERR] Found ',trim(map_projection)
     write(LVT_logunit,*) '[ERR] LVT will skip file ', trim(fname)
     rc = nf90_close(ncid)
     return
  end if

  ! Get dimension IDs
  rc = nf90_inq_dimid(ncid=ncid, &
       name='time', &
       dimid=dim_ids(3))
  if (rc .ne. 0) then
     write(LVT_logunit,*)'[ERR] Cannot read time dimension from ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inq_dimid(ncid=ncid, &
       name='lat', &
       dimid=dim_ids(2))
  if (rc .ne. 0) then
     write(LVT_logunit,*)'[ERR] Cannot read lat dimension from ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inq_dimid(ncid=ncid, &
       name='lon', &
       dimid=dim_ids(1))
  if (rc .ne. 0) then
     write(LVT_logunit,*)'[ERR] Cannot read lon dimension from ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  ! Get actual dimension sizes
  rc = nf90_inquire_dimension(ncid=ncid, &
       dimid=dim_ids(3), &
       len=ntime)
  if (rc .ne. 0) then
     write(LVT_logunit,*)'[ERR] Cannot read time dimension from ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inquire_dimension(ncid=ncid, &
       dimid=dim_ids(2), &
       len=nlat)
  if (rc .ne. 0) then
     write(LVT_logunit,*)'[ERR] Cannot read lat dimension from ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  rc = nf90_inquire_dimension(ncid=ncid, &
       dimid=dim_ids(1), &
       len=nlon)
  if (rc .ne. 0) then
     write(LVT_logunit,*)'[ERR] Cannot read lon dimension from ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Sanity check the dimensions
  if (ntime .ne. 1) then
     write(LVT_logunit,*)'[ERR] Expected time dimension to be 1'
     write(LVT_logunit,*)'[ERR] Found ', ntime, ' from ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  if (nlat .ne. SMAPEOPLsmobs(n)%nr) then
     write(LVT_logunit,*)'[ERR] Expected lat dimension to be ', &
          SMAPEOPLsmobs(n)%nr
     write(LVT_logunit,*)'[ERR] Found ', nlat, ' from ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if
  if (nlon .ne. SMAPEOPLsmobs(n)%nc) then
     write(LVT_logunit,*)'[ERR] Expected lon dimension to be ', &
          SMAPEOPLsmobs(n)%nc
     write(LVT_logunit,*)'[ERR] Found ', nlon, ' from ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Fetch the variable id
  rc = nf90_inq_varid(ncid=ncid, &
       name='arfs_sm', &
       varid=var_id)
  if (rc .ne. 0) then
     write(LVT_logunit,*)'[ERR] Cannot read arfs_sm ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     return
  end if

  ! Read the retrievals
  allocate(tmp(nlon, nlat, ntime))
  rc = nf90_get_var(ncid=ncid, &
       varid=var_id, &
       values=tmp)
  if (rc .ne. 0) then
     write(LVT_logunit,*)'[ERR] Cannot read arfs_sm ', trim(fname)
     write(LVT_logunit,*)'[ERR] LVT will continue...'
     rc = nf90_close(ncid)
     deallocate(tmp)
     return
  end if
  rc = nf90_close(ncid)

  do r = 1, nlat
     do c = 1, nlon
        if (tmp(c,r,1) >= 0 .and. &
             tmp(c,r,1) <= 1) then
           sm_in(c + (r-1)*nc) = tmp(c,r,1) !*100  multiplying by 100 was a bug  
           sm_data_b(c + (r-1)*nc) = .true.
        end if
     end do
  end do
  deallocate(tmp)

#endif



  if(LVT_isAtAfinerResolution(0.0937500)) then
     call bilinear_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_in, smobs_b_ip, smobs_ip, &
          SMAPEOPLsmobs(n)%nc*SMAPEOPLsmobs(n)%nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          SMAPEOPLsmobs(n)%rlat, &
          SMAPEOPLsmobs(n)%rlon, &
          SMAPEOPLsmobs(n)%w11,SMAPEOPLsmobs(n)%w12,&
          SMAPEOPLsmobs(n)%w21,SMAPEOPLsmobs(n)%w22,&
          SMAPEOPLsmobs(n)%n11,SMAPEOPLsmobs(n)%n12,&
          SMAPEOPLsmobs(n)%n21,SMAPEOPLsmobs(n)%n22,LVT_rc%udef,ios)
  else
     call upscaleByAveraging(SMAPEOPLsmobs(n)%nc*SMAPEOPLsmobs(n)%nr,&
          LVT_rc%lnc*LVT_rc%lnr, &
          LVT_rc%udef, SMAPEOPLsmobs(n)%n11,&
          sm_data_b,sm_in, smobs_b_ip, smobs_ip)
  endif

!overwrite the input data
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(smobs_ip(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then
           smobs_inp(c,r) = &
                smobs_ip(c+(r-1)*LVT_rc%lnc)
        endif
     enddo
  enddo

end subroutine read_SMAPEOPLsm_data


