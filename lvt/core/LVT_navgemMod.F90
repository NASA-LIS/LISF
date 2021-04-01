!-----------------------BEGIN NOTICE -- DO NOT EDIT----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT----------------------------

#include "LVT_misc.h"
#include "LVT_NetCDF_inc.h"

!------------------------------------------------------------------------------
!
! MODULE: LVT_navgemMod
!
! DESCRIPTION:
! Contains routines for reading skin temperature and sea ice from
! NAVGEM HDF5 restart files on thinned Gaussian grids, and interpolate to
! LVT grid.  Intended to run as part of 557post mode for Air Force operations.
!
! REVISION HISTORY:
! 01 Apr 2021: Eric Kemp (SSAI), Initial implementation.  Read codes borrow
!              heavily from sample Python code provided by FNMOC.
!------------------------------------------------------------------------------

module LVT_navgemMod

  ! Defaults
  implicit none
  private

  ! Public routines

contains

  subroutine construct_navgem_filename(rootdir, year, month, day, hour, &
       fcst_hr, filename)

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in)  :: rootdir
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in) :: hour
    integer, intent(in) :: fcst_hr
    character(len=*), intent(out) :: filename

    ! Local variables
    character(len=10) :: yyyymmddhh
    character(len=6)  :: hhhhhh

    write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') year, month, day, hour
    write(hhhhhh,'(i6.6)') fcst_hr

    filename = trim(rootdir) // '/navgem_restart_T0681L060_slthin_quad_' &
         // yyyymmddhh // '_' // hhhhhh // '.h5'

  end subroutine construct_navgem_filename

  subroutine get_navgem_filename(filename, &
       year, month, day, hour, fcst_hr)

    ! Modules
    use LVT_coreMod, only: LVT_rc
    use LVT_logMod, only: LVT_logunit
    use LVT_timeMgrMod, only: LVT_get_julhr, LVT_julhr_date

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(inout) :: filename
    integer, intent(out) :: year
    integer, intent(out) :: month
    integer, intent(out) :: day
    integer, intent(out) :: hour
    integer, intent(out) :: fcst_hr

    ! Locals
    integer :: navgem_julhr, lvt_julhr
    logical :: file_exists

    ! FIXME...Add dynamic search for nearest NAVGEM file
    year = 2021
    month = 03
    day = 31
    hour = 18
    fcst_hr = 0
    call construct_navgem_filename('./navgem', &
         year, month, day, hour, fcst_hr, filename)

    write(LVT_logunit,*)'[INFO] *** Searching for NAVGEM file ', &
         trim(filename)
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
       write(LVT_logunit,*)'[INFO] Will use ', trim(filename)
       return
    end if

    ! FIXME...Add dynamic search for NAVGEM file
    write(LVT_logunit,*)'[ERR] Cannot find NAVGEM file!'
    stop

  end subroutine get_navgem_filename

  subroutine fetch_navgem_fields(sst, cice, lat, lon)

    ! Modules
#if (defined USE_HDF5)
    use HDF5
#endif

    ! Defaults
    implicit none

    ! Arguments
    real, allocatable, intent(out) :: sst(:)
    real, allocatable, intent(out) :: cice(:)
    real, allocatable, intent(out) :: lat(:)
    real, allocatable, intent(out) :: lon(:)

    ! Locals
    character(len=250) :: filename
    integer :: year, month, day, hour, fcst_hr
    logical :: fail
    integer :: hdferr
#if (defined USE_HDF5)
    integer(HID_T) :: file_id, dataset_id, datatype_id
#endif

    ! Get NAVGEM filename
    call get_navgem_filename(filename, year, month, day, hour, fcst_hr)

#if (defined USE_HDF5)
    ! Initialize IDs.  Useful later for error handling.
    file_id = -1
    dataset_id = -1
    datatype_id = -1

    ! Initialize HDF5 Fortran interface
    call open_hdf5_f_interface(fail)
    if (fail) goto 100

    ! Open the file
    call open_navgem_file(filename, file_id, fail)
    if (fail) goto 100

    ! Cleanup before returning
100 continue
    if (file_id .gt. -1) call close_navgem_file(filename, file_id, fail)
    call close_hdf5_f_interface(fail)

#endif
  end subroutine fetch_navgem_fields

#if (defined USE_HDF5)
  subroutine open_hdf5_f_interface(fail)
    use HDF5
    use LVT_logMod, only: LVT_logunit
    implicit none
    logical,intent(out) :: fail
    integer :: hdferr
    fail = .false.
    call h5open_f(hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot initialize HDF5 ', &
            'Fortran interface!'
       fail = .true.
    end if
  end subroutine open_hdf5_f_interface
#endif

#if (defined USE_HDF5)
  subroutine close_navgem_file(filename, file_id, fail)
    use HDF5
    use LVT_logMod, only: LVT_logunit
    implicit none
    character(len=*), intent(in) :: filename
    integer(HID_T), intent(inout) :: file_id
    logical, intent(out) :: fail
    integer :: hdferr
    fail = .false.
    call h5fclose_f(file_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot close file ', trim(filename)
       fail = .true.
    else
       write(LVT_logunit,*) &
            '[INFO] Closed NAVGEM file ',trim(filename)
    end if
    file_id = -1
  end subroutine close_navgem_file
#endif

#if (defined USE_HDF5)
  subroutine open_navgem_file(filename, file_id, fail)
    use HDF5
    use LVT_logMod, only: LVT_logunit
    implicit none
    character(len=*), intent(in) :: filename
    integer(HID_T), intent(out) :: file_id
    logical, intent(out) :: fail
    integer :: hdferr
    fail = .false.
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*)&
            '[ERR] Cannot open file ', trim(filename)
       fail = .true.
    else
       write(LVT_logunit,*) &
            '[INFO] Opened NAVGEM file ', trim(filename)
    end if
  end subroutine open_navgem_file
#endif

#if (defined USE_HDF5)
  subroutine close_hdf5_f_interface(fail)
    use HDF5
    use LVT_logMod, only: LVT_logunit
    implicit none
    logical, intent(out) :: fail
    integer :: hdferr
    fail = .false.
    call h5close_f(hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot close HDF5 Fortran interface!'
       fail = .true.
    end if
  end subroutine close_hdf5_f_interface
#endif
end module LVT_navgemMod
