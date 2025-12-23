!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: USAFSI_espcdMod
!
! REVISION HISTORY:
! 17 Jul 2024  Eric Kemp      First version.  (Based on USAF_gofsMod.F90)
! 16 Dec 2024  Yeosang Yoon   Updated ESPC-D SST file format (depth dimensions)
!
! DESCRIPTION:
! Source code for reading US Navy ESPC-D data.
!-------------------------------------------------------------------------

#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module USAFSI_espcdMod

  ! Defaults
  implicit none
  private

  ! Public routines
  public :: process_espcd_sst
  public :: process_espcd_cice

contains

  ! Find ESPCD CICE file on file system
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  subroutine find_espcd_cice_file(rootdir, region, &
       yyyy, mm, dd, hh, filename, &
       aice, nlon, nlat)

    ! Imports
    use netcdf
    use LDT_constantsMod, only: LDT_CONST_PATH_LEN
    use LDT_logMod, only: LDT_logunit, LDT_endrun
    use LDT_timeMgrMod, only: LDT_get_julhr, LDT_julhr_date
    use USAFSI_paramsMod, only: program_name, msglns
    use USAFSI_utilMod, only: error_message

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: rootdir
    character*3, intent(in) :: region
    integer, intent(in) :: yyyy
    integer, intent(in) :: mm
    integer, intent(in) :: dd
    integer, intent(in) :: hh
    character(len=LDT_CONST_PATH_LEN), intent(out) :: filename
    real, allocatable, intent(inout) :: aice(:,:,:)
    integer, intent(out) :: nlon
    integer, intent(out) :: nlat

    ! Locals
    integer :: julhr, julhr_orig
    logical :: file_exists
    integer :: yyyy_local, mm_local, dd_local, hh_local
    integer :: fh_local
    character(len=LDT_CONST_PATH_LEN) :: message (msglns)
    character*20 :: routine_name
    character*10 :: yyyymmddhh
    integer, parameter :: nlat_arc = 2501
    integer, parameter :: nlat_ant = 1549
    integer :: ncid, aice_varid
    integer, allocatable :: dimids(:), lens(:)
    integer :: ndims
    logical :: good
    logical :: first_time
    integer :: ierr
    integer :: i

    nlon = 9000

    message = ''
    routine_name = 'find_espcd_cice_file'

    ! Build the file name.  Note that all ESPC-D CICE runs start at 12Z.
    ! NOTE:  CICE output is 12-hrly.
    call LDT_get_julhr(yyyy, mm, dd, hh, 0, 0, julhr)
    write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') yyyy, mm, dd, hh
    if (hh == 12) then
       fh_local = 0
    else if (hh == 18) then
       fh_local = 0
       julhr = julhr - 6
    else if (hh == 00) then
       fh_local = 12
       julhr = julhr - 12
    else if (hh == 06) then
       fh_local = 12
       julhr = julhr - 18
    else
       write(LDT_logunit,*)'[ERR] Bad USAFSI hour ', hh
       write(LDT_logunit,*)'[ERR] Must be 00, 06, 12, or 18'
       write(LDT_logunit,*)'[ERR] LDT will exit...'
       call LDT_endrun()
    end if

    julhr_orig = julhr

    ! Start looping for earlier files
    first_time = .true.
    do

       if (.not. first_time) then
          fh_local = fh_local + 24
          julhr = julhr - 24 ! Roll back to previous 12Z cycle
       end if
       first_time = .false.

       ! Give up after 5 days
       if ( (julhr_orig - julhr) > 24*5) then
          write(LDT_logunit,*)&
               '[WARN] *** GIVING UP ON ESPC-D CICE FOR ', &
               trim(region),' ***'
          write(LDT_logunit,*) &
               '[WARN] *** NO ESPC-D CICE DATA FOR ',trim(region), &
               ' AVAILABLE!!! ***'
          filename = 'NONE'
          return
       end if
       call LDT_julhr_date(julhr, yyyy_local, mm_local, dd_local, &
            hh_local)

       call construct_espcd_cice_filename(rootdir, region, &
            yyyy_local, mm_local, dd_local, hh_local, fh_local, filename)
       inquire(file=trim(filename), exist=file_exists)
       if (.not. file_exists) then
          message(1) = '[WARN] CANNOT FIND FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          cycle
       end if

       ! Try opening the file.
       ierr = nf90_open(path=trim(filename), &
            mode=nf90_nowrite, &
            ncid=ncid)
       if (ierr .ne. nf90_noerr) then
          message(1) = '[WARN] CANNOT OPEN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          cycle
       end if

       ! See if aice is in the file.
       ierr = nf90_inq_varid(ncid, "aice", aice_varid)
       if (ierr .ne. nf90_noerr) then
          message(1) = '[WARN] CANNOT FIND aice IN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          ierr = nf90_close(ncid)
          cycle
       end if

       ! Check the dimensions
       ierr = nf90_inquire_variable(ncid, aice_varid, ndims=ndims)
       if (ierr .ne. nf90_noerr) then
          message(1) = '[WARN] CANNOT GET DIMENSIONS FOR aice IN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          ierr = nf90_close(ncid)
          cycle
       end if
       allocate(dimids(ndims))
       ierr = nf90_inquire_variable(ncid, aice_varid, dimids=dimids)
       if (ierr .ne. nf90_noerr) then
          message(1) = '[WARN] CANNOT GET DIMENSIONS FOR aice IN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          deallocate(dimids)
          ierr = nf90_close(ncid)
          cycle
       end if
       allocate(lens(ndims))
       good = .true.
       do i = 1, ndims
          ierr = nf90_inquire_dimension(ncid, dimids(i), len=lens(i))
          if (ierr .ne. nf90_noerr) then
             message(1) = '[WARN] CANNOT GET DIMENSIONS FOR aice IN FILE'
             message(2) = '[WARN] PATH = ' // trim(filename)
             call error_message(program_name, routine_name, yyyymmddhh, &
                  message)
             deallocate(dimids)
             deallocate(lens)
             ierr = nf90_close(ncid)
             good = .false.
             exit
          end if
       end do
       if (.not. good) cycle

       deallocate(dimids)

       ! Sanity check the dimensions
       if (region .eq. 'arc') then
          if (lens(3) .ne. 1 .or. &
               lens(2) .ne. nlat_arc .or. &
               lens(1) .ne. nlon) then
             message(1) = '[WARN] BAD DIMENSIONS FOR aice IN FILE'
             message(2) = '[WARN] PATH = ' // trim(filename)
             call error_message(program_name, routine_name, yyyymmddhh, &
                  message)
             deallocate(lens)
             ierr = nf90_close(ncid)
             cycle
          end if
          ! Good dimensions
          nlat = nlat_arc
       else if (region .eq. 'ant') then
          if (lens(3) .ne. 1 .or. &
               lens(2) .ne. nlat_ant .or. &
               lens(1) .ne. nlon) then
             message(1) = '[WARN] BAD DIMENSIONS FOR aice IN FILE'
             message(2) = '[WARN] PATH = ' // trim(filename)
             call error_message(program_name, routine_name, yyyymmddhh, &
                  message)
             deallocate(lens)
             ierr = nf90_close(ncid)
             cycle
          end if
          ! Good dimensions
          nlat = nlat_ant
       end if

       ! Allocate the aice array
       deallocate(lens) ! Don't need this anymore
       allocate(aice(nlon, nlat, 1))

       ! Try reading the file
       ierr = nf90_get_var(ncid, aice_varid, aice)
       if (ierr .ne. nf90_noerr) then
          message(1) = '[WARN] CANNOT READ aice IN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          deallocate(aice)
          ierr = nf90_close(ncid)
          cycle
       end if

       ! We have a winner.
       ierr = nf90_close(ncid)
       write(LDT_logunit,*)'[INFO] Will use ', trim(filename)
       return

    end do

  end subroutine find_espcd_cice_file

#else

  ! Dummy version w/o netCDF
  subroutine find_espcd_cice_file(rootdir, region, &
       yyyy, mm, dd, hh, filename, &
       aice, nlon, nlat)

    ! Imports
    use LDT_constantsMod, only: LDT_CONST_PATH_LEN
    use LDT_logMod, only: LDT_logunit, LDT_endrun

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: rootdir
    character*3, intent(in) :: region
    integer, intent(in) :: yyyy
    integer, intent(in) :: mm
    integer, intent(in) :: dd
    integer, intent(in) :: hh
    character(len=LDT_CONST_PATH_LEN), intent(out) :: filename
    real, allocatable, intent(inout) :: aice(:,:,:)
    integer, intent(out) :: nlon
    integer, intent(out) :: nlat

    write(LDT_logunit,*) &
         '[ERR] LDT was compiled without netCDF support!'
    write(LDT_logunit,*) "[ERR] Recompile and try again!"
    call LDT_endrun()

  end subroutine find_espcd_cice_file
#endif

  ! Builds path to ESPC-D CICE netcdf file
  subroutine construct_espcd_cice_filename(rootdir, region, &
       yyyy, mm, dd, hh, fh, filename)

    ! Imports
    use LDT_constantsMod, only: LDT_CONST_PATH_LEN

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: rootdir
    character*3, intent(in) :: region
    integer, intent(in) :: yyyy
    integer, intent(in) :: mm
    integer, intent(in) :: dd
    integer, intent(in) :: hh
    integer, intent(in) :: fh
    character(len=LDT_CONST_PATH_LEN), intent(out) :: filename

    ! Local variables
    character*10 :: yyyymmddhh
    character*5  :: thhhh

    write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') yyyy, mm, dd, hh
    write(thhhh,'(a1,i4.4)') 't', fh

    filename = trim(rootdir) // '/espc-d_031_cice_' // trim(region) &
         // 'u0.04_' // yyyymmddhh // '_' // thhhh // '.nc'

  end subroutine construct_espcd_cice_filename

  ! Find ESPC-D SST file on file system
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  subroutine find_espcd_sst_file(rootdir, yyyy, mm, dd, hh, &
       filename, water_temp, nlat, nlon)

    ! Imports
    use netcdf
    use LDT_constantsMod, only: LDT_CONST_PATH_LEN
    use LDT_logMod, only: LDT_logunit, LDT_endrun
    use LDT_timeMgrMod, only: LDT_get_julhr, LDT_julhr_date
    use USAFSI_paramsMod, only: program_name, msglns
    use USAFSI_utilMod, only: error_message

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: rootdir
    integer, intent(in) :: yyyy
    integer, intent(in) :: mm
    integer, intent(in) :: dd
    integer, intent(in) :: hh
    character(len=LDT_CONST_PATH_LEN), intent(inout) :: filename
    real, allocatable, intent(inout) :: water_temp(:,:,:,:)
    integer, intent(out) :: nlat
    integer, intent(out) :: nlon

    ! Locals
    integer :: julhr, julhr_orig
    integer :: yyyy_local, mm_local, dd_local, hh_local, fh_local
    logical :: file_exists
    character(len=LDT_CONST_PATH_LEN) :: message (msglns)
    character*20 :: routine_name
    character*10 :: yyyymmddhh
    integer :: ncid, water_temp_varid
    integer :: ndims
    integer, allocatable :: dimids(:), lens(:)
    logical :: good
    logical :: first_time
    integer :: ierr
    integer :: i

    nlat = 8001
    nlon = 9000

    message = ''
    routine_name = 'find_espcd_sst_file'

    ! Build the file name.  Note that all ESPC-D SST runs start at 12Z.
    call LDT_get_julhr(yyyy, mm, dd, hh, 0, 0, &
         julhr)
    write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') yyyy, mm, dd, hh
    if (hh == 12) then
       fh_local = 0
    else if (hh == 18) then
       fh_local = 6
    else if (hh == 00) then
       fh_local = 12
    else if (hh == 06) then
       fh_local = 18
    else
       write(LDT_logunit,*)'[ERR] Bad USAFSI hour ', hh
       write(LDT_logunit,*)'[ERR] Must be 00, 06, 12, or 18'
       write(LDT_logunit,*)'[ERR] LDT will exit...'
       call LDT_endrun()
    end if
    julhr = julhr - fh_local
    julhr_orig = julhr

    ! Loop through possible files
    first_time = .true.
    do
       if (.not. first_time) then
          fh_local = fh_local + 24
          julhr = julhr - 24 ! Roll back to previous 12Z cycle
       end if
       first_time = .false.

       ! Give up after 5 days
       if ( (julhr_orig - julhr) > 24*5) then
          write(LDT_logunit,*)"[WARN] *** GIVING UP ON ESPC-D SST! ***"
          write(LDT_logunit,*) &
               "[WARN] *** NO ESPC-D SST AVAILABLE!!! ***"
          filename = "NONE"
          return
       end if

       call LDT_julhr_date(julhr, yyyy_local, mm_local, dd_local, &
            hh_local)
       call construct_espcd_sst_filename(rootdir, &
            yyyy_local, mm_local, dd_local, hh_local, fh_local, filename)

       ! See if file exists
       inquire(file=trim(filename), exist=file_exists)
       if (.not. file_exists) then
          message(1) = '[WARN] CANNOT FIND FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          cycle
       end if

       ! Try opening the file
       ierr = nf90_open(path=trim(filename), &
            mode=nf90_nowrite, &
            ncid=ncid)
       if (ierr .ne. nf90_noerr) then
          message(1) = '[WARN] CANNOT OPEN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          cycle
       end if

       ! See if water_temp is in the file.
       ierr = nf90_inq_varid(ncid, "water_temp", &
            water_temp_varid)
       if (ierr .ne. nf90_noerr) then
          message(1) = '[WARN] CANNOT FIND water_temp IN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          ierr = nf90_close(ncid)
          cycle
       end if

       ! Check the dimensions
       ierr = nf90_inquire_variable(ncid, water_temp_varid, ndims=ndims)
       if (ierr .ne. nf90_noerr) then
          message(1) = &
               '[WARN] CANNOT GET DIMENSIONS FOR water_temp IN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          ierr = nf90_close(ncid)
          cycle
       end if
       if (ndims .ne. 4) then
          message(1) = &
               '[WARN] BAD DIMENSIONS FOR water_temp IN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          ierr = nf90_close(ncid)
          cycle
       end if
       allocate(dimids(ndims))
       ierr = nf90_inquire_variable(ncid, water_temp_varid, dimids=dimids)
       if (ierr .ne. nf90_noerr) then
          message(1) = &
               '[WARN] CANNOT GET DIMENSIONS FOR water_temp IN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          deallocate(dimids)
          ierr = nf90_close(ncid)
          cycle
       end if
       good = .true.
       allocate(lens(ndims))
       do i = 1, ndims
          ierr = nf90_inquire_dimension(ncid, dimids(i), &
               len=lens(i))
          if (ierr .ne. nf90_noerr) then
             message(1) = &
                  '[WARN] CANNOT GET DIMENSIONS FOR water_temp IN FILE'
             message(2) = '[WARN] PATH = ' // trim(filename)
             call error_message(program_name, routine_name, yyyymmddhh, &
                  message)
             deallocate(dimids)
             deallocate(lens)
             ierr = nf90_close(ncid)
             good = .false.
             exit
          end if
       end do
       if (.not. good) cycle
       deallocate(dimids)

       if (lens(1) .ne. nlon .or. &
            lens(2) .ne. nlat .or. &
            lens(3) .ne. 2 .or. &          !depth=2, updated 2024/10/31 
            lens(4) .ne. 1) then
          message(1) = &
               '[WARN] CANNOT GET DIMENSIONS FOR water_temp IN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          deallocate(lens)
          ierr = nf90_close(ncid)
          cycle
       end if
       deallocate(lens)

       ! Allocate a subset of water temp
       allocate(water_temp(nlon, nlat, 1, 1))
       water_temp = 0

       ierr = nf90_get_var(ncid, water_temp_varid, water_temp)
       if (ierr .ne. nf90_noerr) then
          message(1) = &
               '[WARN] CANNOT READ water_temp IN FILE'
          message(2) = '[WARN] PATH = ' // trim(filename)
          call error_message(program_name, routine_name, yyyymmddhh, &
               message)
          deallocate(water_temp)
          ierr = nf90_close(ncid)
          cycle
       end if

       ! We have a winner
       write(LDT_logunit,*)'[INFO] Will use ',trim(filename)
       return

    end do

  end subroutine find_espcd_sst_file

#else
  ! Dummy version w/o netCDF
  subroutine find_espcd_sst_file(rootdir, yyyy, mm, dd, hh, &
       filename, water_temp, nlat, nlon)

    ! Imports
    use LDT_constantsMod, only: LDT_CONST_PATH_LEN
    use LDT_logMod, only: LDT_logunit, LDT_endrun

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: rootdir
    integer, intent(in) :: yyyy
    integer, intent(in) :: mm
    integer, intent(in) :: dd
    integer, intent(in) :: hh
    character(len=LDT_CONST_PATH_LEN), intent(inout) :: filename
    real, allocatable, intent(inout) :: water_temp(:,:,:,:)
    integer, intent(out) :: nlat
    integer, intent(out) :: nlon

    write(LDT_logunit,*) &
         '[ERR] LDT was compiled without netCDF support!'
    write(LDT_logunit,*) "[ERR] Recompile and try again!"
    call LDT_endrun()

  end subroutine find_espcd_sst_file
#endif

  ! Builds path to ESPC-D SST netcdf file
  subroutine construct_espcd_sst_filename(rootdir, &
       yyyy, mm, dd, hh, fh, filename)

    ! Imports
    use LDT_constantsMod, only: LDT_CONST_PATH_LEN

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: rootdir
    integer, intent(in) :: yyyy
    integer, intent(in) :: mm
    integer, intent(in) :: dd
    integer, intent(in) :: hh
    integer, intent(in) :: fh
    character(len=LDT_CONST_PATH_LEN), intent(out) :: filename

    ! Locals
    character*10 :: yyyymmddhh
    character*5  :: thhhh

    write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') yyyy, mm, dd, hh
    write(thhhh,'(a1,i4.4)') 't', fh

    filename = trim(rootdir) // "/espc-d_hycom_sfc_u_" // yyyymmddhh // &
         "_" // thhhh // ".nc"

  end subroutine construct_espcd_sst_filename


#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  ! Read ESPC-D sea surface temperature and reproject to LDT grid
  subroutine process_espcd_sst(rootdir, nc, nr, landmask, sst, &
       yyyy, mm, dd, hh, ierr)

    ! Imports
    use LDT_constantsMod, only: LDT_CONST_PATH_LEN
    use LDT_coreMod, only: LDT_rc, LDT_domain
    use LDT_logMod, only: LDT_verify, LDT_endrun
    use netcdf

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: rootdir
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    real, intent(in) :: landmask(nc,nr)
    real, intent(inout) :: sst(nc,nr)
    integer, intent(inout) :: yyyy
    integer, intent(inout) :: mm
    integer, intent(inout) :: dd
    integer, intent(inout) :: hh
    integer, intent(out) :: ierr

    ! Locals
    integer :: nlat
    integer :: nlon
    character(LDT_CONST_PATH_LEN) :: filename
    real, allocatable :: water_temp(:,:,:,:)
    real, allocatable :: water_temp_1d(:)
    real, allocatable :: sst_1d(:)
    integer :: c, r, c1, r1
    logical*1, allocatable :: lb(:)
    logical*1, allocatable :: lo(:)
    real :: griddesci(50)
    real, allocatable :: n11(:)
    integer :: gindex
    real :: rlat

    ! External subroutines
    external :: upscaleByAveraging_input
    external :: upscaleByAveraging

    ! Find a valid file on the file system
    call find_espcd_sst_file(rootdir, yyyy, mm, dd, hh, filename, &
         water_temp, nlat, nlon)
    if (trim(filename) == "NONE") then
       ierr = 1
       return
    end if

    ! We need to interpolate to the LDT grid.  First, copy to 1D array
    allocate(water_temp_1d(nlon*nlat*1*1))
    water_temp_1d = -9999.0
    allocate(lb(nlon*nlat*1*1))
    lb = .false.
    do r = 1, nlat
       do c = 1, nlon
          if (water_temp(c,r,1,1) .eq. -30000) cycle ! Missing value
          ! Convert from Celsius to Kelvin, taking into account the scale
          ! factor and offset.  Also, wrap the data so it starts at 180W
          if (c .gt. 4500) then
             c1 = c - 4500
             r1 = r
          else
             c1 = c + 4500
             r1 = r
          end if
          water_temp_1d(c1 + (r1-1)*nlon) = &
               (0.001*water_temp(c,r,1,1)) + 20.0 + 273.15
          lb(c1 + (r1-1)*nlon) = .true.
       end do ! c
    end do ! r
    deallocate(water_temp)

    ! Set up interpolation weights
    gridDesci = 0
    gridDesci(1) = 0          ! Lat/lon projection
    gridDesci(2) = nlon       ! Number of columns
    gridDesci(3) = nlat       ! Number of rows
    gridDesci(4) =  -80.0     ! Lower-left latitude (deg)
    gridDesci(5) = -180.0     ! Lower-left longitude (deg)
    gridDesci(6) =  128       ! Not used
    gridDesci(7) =   80.0     ! Upper-right latitude (deg)
    gridDesci(8) =  180.0     ! Upper-right longitude(deg)
    gridDesci(9) =    0.0400390625        ! Delta longitude (deg)
    gridDesci(10) =   0.01999664306640625 ! Delta latitude (deg)
    gridDesci(20) =  64       ! E-W ordering
    allocate(n11(nlon*nlat))
    call upscaleByAveraging_input(gridDesci, LDT_rc%gridDesc, &
         nlon*nlat, nc*nr, n11)

    ! Now interpolate
    allocate(sst_1d(nc*nr))
    sst_1d = -9999.
    allocate(lo(nc*nr))
    lo = .false.
    call upscaleByAveraging(nlon*nlat, nc*nr, -9999., &
         n11, lb, water_temp_1d, lo, sst_1d)

    ! Since SST is missing north of 80N, we need to set water points in
    ! this region to a reasonable value.  We follow the typical
    ! UKMET SURF value of 271.35K.
    sst = -1
    do r = 1, nr
       do c = 1, nc
          ! Skip land points
          if (landmask(c,r) >= 0.5) cycle

          gindex = c + (r-1)*nc
          rlat = LDT_domain(1)%lat(gindex)
          if (rlat >= 80.) then
             sst(c,r) = 271.35
          else
             if (sst_1d(gindex) > 0) then
                sst(c,r) = sst_1d(gindex)
             end if
          end if
       end do ! c
    end do ! r

    ! Clean up
    deallocate(water_temp_1d)
    deallocate(lb)
    deallocate(lo)
    deallocate(sst_1d)
    deallocate(n11)

    ! The end
    ierr = 0
  end subroutine process_espcd_sst

#else

   ! Dummy version with no netCDF support
  subroutine process_espcd_sst(rootdir, nc, nr, landmask, sst, &
       yyyy, mm, dd, hh, fh, vierr)

    ! Imports
    use LDT_logMod, only: LDT_logunit, LDT_endrun

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: rootdir
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    real, intent(in) :: landmask(nc,nr)
    real, intent(inout) :: sst(nc,nr)
    integer, intent(inout) :: yyyy
    integer, intent(inout) :: mm
    integer, intent(inout) :: dd
    integer, intent(inout) :: hh
    integer, intent(inout) :: fh
    integer, intent(out) :: ierr

    write(LDT_logunit,*) &
         '[ERR] LDT was compiled without netCDF support!'
    write(LDT_logunit,*) "[ERR] Recompile and try again!"
    ierr = 1
    call LDT_endrun()
  end subroutine process_espcd_sst

#endif

  ! Read ESPC-D sea ice and reproject to LDT grid
  subroutine process_espcd_cice(rootdir, nc, nr, landmask, icecon, &
       yyyy, mm, dd, hh, ierr)

    ! Imports
    use LDT_coreMod, only: LDT_domain

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: rootdir
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    real, intent(in) :: landmask(nc,nr)
    real, intent(inout) :: icecon(nc,nr)
    integer, intent(inout) :: yyyy
    integer, intent(inout) :: mm
    integer, intent(inout) :: dd
    integer, intent(inout) :: hh
    integer, intent(out) :: ierr

    ! Locals
    real, allocatable :: icecon_arc(:,:)
    real, allocatable :: icecon_ant(:,:)
    integer :: c, r
    integer :: gindex
    real :: rlat

    ! First handle Arctic region
    call process_espcd_cice_region('arc', rootdir, nc, nr, landmask, &
         yyyy, mm, dd, hh, icecon_arc, ierr)
    if (ierr .ne. 0) then
       if (allocated(icecon_arc)) deallocate(icecon_arc)
       return
    end if

    ! Next handle Antarctic region
    call process_espcd_cice_region('ant', rootdir, nc, nr, landmask, &
         yyyy, mm, dd, hh, icecon_ant, ierr)
    if (ierr .ne. 0) then
       if (allocated(icecon_arc)) deallocate(icecon_arc)
       if (allocated(icecon_ant)) deallocate(icecon_ant)
       return
    end if

    ! Merge the two regions together
    icecon = -1
    do r = 1, nr
       do c = 1, nc
          ! Skip land points
          if (landmask(c,r) > 0.5) cycle

          gindex = c + (r-1)*nc
          rlat = LDT_domain(1)%lat(gindex)
          if (rlat >= 0) then
             icecon(c,r) = icecon_arc(c,r)
          else
             icecon(c,r) = icecon_ant(c,r)
          end if

       end do ! c
    end do ! r

    ! Clean up
    deallocate(icecon_arc)
    deallocate(icecon_ant)
    ierr = 0
  end subroutine process_espcd_cice

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  ! Process a particular region of ESPC-D CICE data (Arctic or Antarctic
  subroutine process_espcd_cice_region(region, rootdir, nc, nr, &
       landmask, yyyy, mm, dd, hh, icecon, ierr)

    ! Imports
    use LDT_constantsMod, only: LDT_CONST_PATH_LEN
    use LDT_coreMod, only: LDT_rc
    use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify
    use netcdf

    ! Arguments
    character*3, intent(in) :: region
    character(len=*), intent(in) :: rootdir
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    real, intent(in) :: landmask(nc,nr)
    integer, intent(inout) :: yyyy
    integer, intent(inout) :: mm
    integer, intent(inout) :: dd
    integer, intent(inout) :: hh
    real, allocatable, intent(out) :: icecon(:,:)
    integer, intent(out) :: ierr

    ! Locals
    integer, parameter :: nlat_arc = 2501
    integer, parameter :: nlat_ant = 1549
    integer :: nlon = 9000
    character(LDT_CONST_PATH_LEN) :: filename
    real, allocatable :: aice(:,:,:)
    real, allocatable :: aice_1d(:)
    real, allocatable :: icecon_1d(:)
    integer :: c, r
    logical*1, allocatable :: lb(:)
    logical*1, allocatable :: lo(:)
    real :: griddesci(50)
    real, allocatable :: n11(:)
    integer :: gindex, nlat

    ! External subroutines
    external :: upscaleByAveraging_input
    external :: upscaleByAveraging

    ! Sanity check the region
    if (region .ne. 'arc' .and. &
         region .ne. 'ant') then
       write(LDT_logunit,*)'[ERR] Invalid ESPC-D region for cice: ' &
            // region
       write(LDT_logunit,*)'[ERR] Must be either arc or ant'
       ierr = 1
       call LDT_endrun()
    end if

    ! Find a valid file on the file system
    call find_espcd_cice_file(rootdir, region, yyyy, mm, dd, hh, &
         filename, &
         aice, nlon, nlat)
    if (trim(filename) == "NONE") then
       ierr = 1
       return
    end if

    allocate(aice_1d(nlon*nlat*1))
    aice_1d = -9999
    allocate(lb(nlon*nlat*1))
    lb = .false.
    do r = 1, nlat
       do c = 1, nlon
          if (aice(c,r,1) .eq. -30000) cycle

          ! Take into account the scale factor, and convert to %
          aice_1d(c + (r-1)*nlon) = &
               aice(c,r,1)*0.0001*100
          lb(c + (r-1)*nlon) = .true.
       end do ! c
    end do ! r
    deallocate(aice)

    ! Set up interpolation weights
    if (region .eq. 'arc') then
       gridDesci = 0
       gridDesci(1) = 0 ! Lat/lon projection
       gridDesci(2) = nlon    ! Number of columns
       gridDesci(3) = nlat    ! Number of rows
       gridDesci(4) =   40.   ! Lower-left latitude (deg N)
       gridDesci(5) = -180.0  ! Lower-left longitude (deg E)
       gridDesci(6) = 128     ! Not used
       gridDesci(7) =   90.0          ! Upper-right latitude (deg N)
       gridDesci(8) =  179.9599609375 ! Upper-right longitude (deg E)
       gridDesci(9) =    0.039978027344005795 ! delta-lon (deg)
       gridDesci(10) =   0.0200004577637 !     delta-lat (deg)
       gridDesci(20) = 64  ! East-west ordering
    else if (region .eq. 'ant') then
       gridDesci = 0
       gridDesci(1) = 0 ! Lat/lon projection
       gridDesci(2) = nlon  ! Number of columns
       gridDesci(3) = nlat  ! Number of rows
       gridDesci(4) =  -80.4800033569336 ! Lower-left latitude (deg N)
       gridDesci(5) = -180.0  ! Lower-left longitude (deg E)
       gridDesci(6) =  128    ! Not used
       gridDesci(7) =  -49.5200004577637 ! Upper-right latitude (deg N)
       gridDesci(8) =  179.9599609375    ! Upper-right longitude (deg E)
       gridDesci(9) =    0.039978027344005795 ! Delta-lon (deg)
       gridDesci(10) =   0.020004272460894867 ! Delta-lat (deg)
       gridDesci(20) = 64  ! East-west ordering
    end if
    allocate(n11(nlon*nlat))

    call upscaleByAveraging_input(gridDesci, LDT_rc%gridDesc, &
         nlon*nlat, nc*nr, n11)

    ! Now interpolate
    allocate(icecon_1d(nc*nr))
    icecon_1d = -9999
    allocate(lo(nc*nr))
    lo = .false.
    call upscaleByAveraging(nlon*nlat, nc*nr, -9999., &
         n11, lb, aice_1d, lo, icecon_1d)

    ! Just copy the non-missing values to the output array.  This should
    ! prevent overwriting of data outside of the ESPC-D polar region.
    allocate(icecon(nc,nr))
    icecon = 0.0
    do r = 1, nr
       do c = 1, nc
          ! Skip land points
          if (landmask(c,r) >= 0.5) cycle

          gindex = c + (r-1)*nc
          if (icecon_1d(gindex) .ne. -9999) then
             icecon(c,r) = icecon_1d(gindex)
          end if
       end do ! c
    end do ! r

    ! Clean up
    deallocate(aice_1d)
    deallocate(lb)
    deallocate(lo)
    deallocate(icecon_1d)
    deallocate(n11)

    ! The end
    ierr = 0
  end subroutine process_espcd_cice_region

#else
  ! Dummy version
  subroutine process_espcd_cice_region(region, rootdir, nc, nr, &
       landmask, icecon, yyyy, mm, dd, hh, fh, ierr)

    ! Imports
    use LDT_logMod, only: LDT_logunit, LDT_endrun

    ! Arguments
    character*3, intent(in) :: region
    character(len=*), intent(in) :: rootdir
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    real, intent(in) :: landmask(nc,nr)
    real, intent(inout) :: icecon(nc,nr)
    integer, intent(inout) :: yyyy
    integer, intent(inout) :: mm
    integer, intent(inout) :: dd
    integer, intent(inout) :: hh
    integer, intent(inout) :: fh
    integer, intent(out) :: ierr

    write(LDT_logunit,*) &
         '[ERR] LDT was compiled without netCDF support!'
    write(LDT_logunit,*) "[ERR] Recompile and try again!"
    ierr = 1
    call LDT_endrun()

  end subroutine process_espcd_cice_region

#endif

end module USAFSI_espcdMod
