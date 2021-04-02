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
! 01 Apr 2021: Eric Kemp (SSAI), Initial implementation.  Basic logic for
!              pulling fields and calculating latitudes and longitudes is
!              based on sample Python code provided by FNMOC.  HDF5 logic
!              borrows from IMERG reader in LIS.
!------------------------------------------------------------------------------

module LVT_navgemMod

  ! Defaults
  implicit none
  private

  ! Public routines
  public :: LVT_upscaleByAveraging_input_navgem

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
    use LVT_logMod, only: LVT_logunit

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
    integer(HSIZE_T), allocatable :: dims(:)
#endif
    real, allocatable :: tmp_gt(:,:)
    real, allocatable :: tmp_conice(:,:)
    real, allocatable :: tmp_latitudes(:,:)
    integer, allocatable :: tmp_points_per_lat(:,:)
    integer :: rank
    integer :: im, itmp, t_number

    ! Handle case where LVT was not compiled with HDF5 support
#if (!defined USE_HDF5)
    write(LVT_logunit,*)'[ERR] Cannot read NAVGEM HDF5 file!'
    write(LVT_logunit,*) &
         '[ERR] Reconfigure with HDF5, recompile, and try again!'
    stop
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

    ! Get the gt ("ground temperature") field
    call open_navgem_dataset(file_id, "/Grid/gt", dataset_id, fail)
    if (fail) goto 100
    call get_navgem_datatype(dataset_id, datatype_id, fail)
    if (fail) goto 100
    call check_navgem_type(datatype_id, H5T_IEEE_F32LE, fail)
    if (fail) goto 100
    call check_navgem_units(dataset_id, "K", fail)
    if (fail) goto 100
    call get_navgem_dims(dataset_id, rank, dims, fail)
    if (fail) goto 100
    if (rank .ne. 2) then
       write(LVT_logunit,*)'[ERR] HDF5 dataset /Grid/gt has wrong rank!'
       write(LVT_logunit,*)'Expected 2, found ', rank
       goto 100
    end if
    if (dims(2) .ne. 1) then
       write(LVT_logunit,*) &
            '[ERR] Unexpected first dimension for HDF5 dataset /Grid/gt!'
       write(LVT_logunit,*) 'Expected 1, found ', dims(2)
       goto 100
    end if
    allocate(tmp_gt(dims(1), dims(2)))
    tmp_gt = 0
    call h5dread_f(dataset_id, H5T_IEEE_F32LE, tmp_gt, dims, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*)'[ERR] Cannot read HDF5 dataset /Grid/gt!'
       goto 100
    end if

    ! Close the /Grid/gt types
    if (datatype_id .gt. -1) call close_navgem_datatype(datatype_id, fail)
    if (dataset_id .gt. -1) call close_navgem_dataset(dataset_id, fail)

    ! Save the data into the sst array.
    allocate(sst(dims(2)))
    sst = tmp_gt(:,1)
    deallocate(tmp_gt)
    deallocate(dims)

    ! Get the conice (sea ice area fraction) field
    call open_navgem_dataset(file_id, "/Grid/conice", dataset_id, fail)
    if (fail) goto 100
    call get_navgem_datatype(dataset_id, datatype_id, fail)
    if (fail) goto 100
    call check_navgem_type(datatype_id, H5T_IEEE_F32LE, fail)
    if (fail) goto 100
    call get_navgem_dims(dataset_id, rank, dims, fail)
    if (fail) goto 100
    if (rank .ne. 2) then
       write(LVT_logunit,*)'[ERR] HDF5 dataset /Grid/conice has wrong rank!'
       write(LVT_logunit,*)'Expected 2, found ', rank
       goto 100
    end if
    if (dims(2) .ne. 1) then
       write(LVT_logunit,*) &
            '[ERR] Unexpected first dimension for HDF5 dataset /Grid/conice!'
       write(LVT_logunit,*) 'Expected 1, found ', dims(2)
       goto 100
    end if
    allocate(tmp_conice(dims(1), dims(2)))
    tmp_conice = 0
    call h5dread_f(dataset_id, H5T_IEEE_F32LE, tmp_conice, dims, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*)'[ERR] Cannot read HDF5 dataset /Grid/conice!'
       goto 100
    end if

    ! Close the /Grid/conice types
    if (datatype_id .gt. -1) call close_navgem_datatype(datatype_id, fail)
    if (dataset_id .gt. -1) call close_navgem_dataset(dataset_id, fail)

    ! Save the data into the cice array.
    allocate(cice(dims(2)))
    cice = tmp_conice(:,1)
    deallocate(tmp_conice)
    deallocate(dims)

    ! Get the Gaussian latitudes
    call open_navgem_dataset(file_id, "/Geometry/Latitudes", dataset_id, fail)
    if (fail) goto 100
    call get_navgem_datatype(dataset_id, datatype_id, fail)
    if (fail) goto 100
    call check_navgem_type(datatype_id, H5T_IEEE_F32LE, fail)
    if (fail) goto 100
    call get_navgem_dims(dataset_id, rank, dims, fail)
    if (fail) goto 100
    if (rank .ne. 2) then
       write(LVT_logunit,*) &
            '[ERR] HDF5 dataset /Geometry/Latitudes has wrong rank!'
       write(LVT_logunit,*)'Expected 2, found ', rank
       goto 100
    end if
    if (dims(2) .ne. 1) then
       write(LVT_logunit,*) &
            '[ERR] Unexpected first dimension for HDF5 dataset ', &
            '/Geometry/Latitudes!'
       write(LVT_logunit,*) 'Expected 1, found ', dims(2)
       goto 100
    end if
    allocate(tmp_latitudes(dims(1), dims(2)))
    tmp_latitudes = 0
    call h5dread_f(dataset_id, H5T_IEEE_F32LE, tmp_latitudes, dims, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot read HDF5 dataset /Geometry/Latitudes!'
       goto 100
    end if

    ! Close the /Geometry/Latitudes types
    if (datatype_id .gt. -1) call close_navgem_datatype(datatype_id, fail)
    if (dataset_id .gt. -1) call close_navgem_dataset(dataset_id, fail)
    deallocate(dims)

    ! Get the points per latitudes
    call open_navgem_dataset(file_id, "/Geometry/Points_per_lat", dataset_id, &
         fail)
    if (fail) goto 100
    call get_navgem_datatype(dataset_id, datatype_id, fail)
    if (fail) goto 100
    call check_navgem_type(datatype_id, H5T_STD_I32LE, fail)
    if (fail) goto 100
    call get_navgem_dims(dataset_id, rank, dims, fail)
    if (fail) goto 100
    if (rank .ne. 2) then
       write(LVT_logunit,*) &
            '[ERR] HDF5 dataset /Geometry/Points_per_lat has wrong rank!'
       write(LVT_logunit,*)'Expected 2, found ', rank
       goto 100
    end if
    if (dims(2) .ne. 1) then
       write(LVT_logunit,*) &
            '[ERR] Unexpected first dimension for HDF5 dataset ', &
            '/Geometry/Points_per_lat!'
       write(LVT_logunit,*) 'Expected 1, found ', dims(2)
       goto 100
    end if
    allocate(tmp_points_per_lat(dims(1), dims(2)))
    tmp_latitudes = 0
    call h5dread_f(dataset_id, H5T_STD_I32LE, tmp_points_per_lat, dims, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot read HDF5 dataset /Geometry/Points_per_lat!'
       goto 100
    end if

    ! Close the /Geometry/Points_per_lat types
    if (datatype_id .gt. -1) call close_navgem_datatype(datatype_id, fail)
    if (dataset_id .gt. -1) call close_navgem_dataset(dataset_id, fail)

    ! Calculate dimension im
    itmp = floor(dims(1) / 2.) + 1
    im = tmp_points_per_lat(itmp,1)
    deallocate(dims)

    ! Calculate t_number
    call get_navgem_truncation(im, t_number)
    if (t_number .ne. 681) then
       write(LVT_logunit,*)'[ERR] Unexpected T-number for NAVGEM!'
       write(LVT_logunit,*)'Expected T681, found T', t_number
       goto 100
    end if

    ! Now we need to calculate the lat and lon of each sst and cice point.
    call calc_navgem_latlons(tmp_latitudes, tmp_points_per_lat, size(sst), &
       lat, lon)

    ! Clean up temporary arrays
    deallocate(tmp_latitudes)
    deallocate(tmp_points_per_lat)

    write(LVT_logunit,*)'[INFO] Read data from NAVGEM file'

    ! Cleanup before returning
100 continue
    if (allocated(dims)) deallocate(dims)
    if (allocated(tmp_gt)) deallocate(tmp_gt)
    if (allocated(tmp_conice)) deallocate(tmp_conice)
    if (allocated(tmp_latitudes)) deallocate(tmp_latitudes)
    if (allocated(tmp_points_per_lat)) deallocate(tmp_points_per_lat)
    if (datatype_id .gt. -1) call close_navgem_datatype(datatype_id, fail)
    if (dataset_id .gt. -1) call close_navgem_dataset(dataset_id, fail)
    if (file_id .gt. -1) call close_navgem_file(filename, file_id, fail)
    call close_hdf5_f_interface(fail)

#endif
  end subroutine fetch_navgem_fields

#if (defined USE_HDF5)
  subroutine open_hdf5_f_interface(fail)
    use HDF5
    use LVT_logMod, only: LVT_logunit
    implicit none
    logical, intent(out) :: fail
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
            '[ERR] Cannot open NAVGEM file ', trim(filename)
       fail = .true.
    else
       write(LVT_logunit,*) &
            '[INFO] Opened NAVGEM file ', trim(filename)
    end if
  end subroutine open_navgem_file
#endif

#if (defined USE_HDF5)
  subroutine open_navgem_dataset(file_id, dataset_name, dataset_id, fail)
    use HDF5
    use LVT_logMod, only: LVT_logunit
    implicit none
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name
    integer(HID_T), intent(out) :: dataset_id
    integer :: hdferr
    logical, intent(out) :: fail
    fail = .false.
    call h5dopen_f(file_id, trim(dataset_name), dataset_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*)&
            '[ERR] Cannot open HDF5 dataset ', trim(dataset_name)
       fail = .true.
    end if
  end subroutine open_navgem_dataset
#endif

#if (defined USE_HDF5)
  subroutine get_navgem_datatype(dataset_id, datatype_id, fail)
    use HDF5
    use LVT_logMod, only: LVT_logunit
    implicit none
    integer(HID_T), intent(in) :: dataset_id
    integer(HID_T), intent(out) :: datatype_id
    logical, intent(out) :: fail
    integer :: hdferr
    fail = .false.
    call h5dget_type_f(dataset_id, datatype_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*)&
            '[ERR] Cannot determine HDF5 datatype'
       fail = .true.
    end if
  end subroutine get_navgem_datatype
#endif

#if (defined USE_HDF5)
  subroutine check_navgem_type(datatype_id, datatype, fail)
    use HDF5
    use LVT_logMod, only: LVT_logunit
    integer(HID_T), intent(in) :: datatype_id
    integer(HID_T), intent(in) :: datatype
    logical, intent(out) :: fail
    logical :: flag
    integer :: hdferr
    fail = .false.
    call h5tequal_f(datatype_id, datatype, flag, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot confirm HDF5 datatype!'
       fail = .true.
       return
    end if
    if (.not. flag) then
       write(LVT_logunit,*)&
            '[ERR] HDF5 datatype is wrong type!'
       fail = .true.
       return
    end if
  end subroutine check_navgem_type
#endif

#if (defined USE_HDF5)
  subroutine check_navgem_units(dataset_id, units, fail)

    ! Modules
    use HDF5
    use ISO_C_BINDING
    use LVT_logMod, only: LVT_logunit

    ! Defaults
    implicit none

    ! Arguments
    integer(HID_T), intent(in) :: dataset_id
    character(len=*), intent(in) :: units
    logical, intent(out) :: fail

    ! Local variables
    integer(HID_T) :: attr_id, type_id, space_id, memtype_id
    integer :: hdferr
    integer(size_t) :: size
    integer(SIZE_T), parameter :: sdim = 5
    integer(HSIZE_T), dimension(1:1) :: dims = (/1/)
    integer(HSIZE_T), dimension(1:1) :: maxdims
    character(len=sdim), dimension(:), allocatable, target :: rdata
    type(C_PTR) :: f_ptr
    integer :: i

    fail = .false.

    ! Open the attribute
    call h5aopen_f(dataset_id, 'units', attr_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot open HDF5 attribute'
       fail = .true.
       return
    end if

    ! Get the attribute datatype
    call h5aget_type_f(attr_id, type_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot get HDF5 attribute datatype'
       call h5aclose_f(attr_id, hdferr)
       fail = .true.
       return
    end if

    ! Get the size of the attribute datatype, and sanity check.
    call h5tget_size_f(type_id, size, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot get HDF5 attribute ', &
            'datatype size'
       call h5tclose_f(type_id, hdferr)
       call h5aclose_f(attr_id, hdferr)
       fail = .true.
       return
    end if
    if (size .gt. sdim+1) then
       write(LVT_logunit,*) &
            '[ERR] Expected smaller HDF5 attribute',&
            'datatype size'
       write(LVT_logunit,*)'Expected ',sdim+1
       write(LVT_logunit,*)'Found ',size
       call h5tclose_f(type_id, hdferr)
       call h5aclose_f(attr_id, hdferr)
       fail = .true.
       return
    end if

    ! Get the attribute dataspace
    call h5aget_space_f(attr_id, space_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot get HDF5 attribute', &
            'dataspace'
       call h5tclose_f(type_id, hdferr)
       call h5aclose_f(attr_id, hdferr)
       fail = .true.
       return
    end if

    ! Get the dimensions of the dataspace
    call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot get HDF5 attribute ', &
            'dataspace dimensions'
       call h5sclose_f(space_id, hdferr)
       call h5tclose_f(type_id, hdferr)
       call h5aclose_f(attr_id, hdferr)
       fail = .true.
       return
    end if

    ! Create the memory datatype
    call h5tcopy_f(H5T_FORTRAN_S1, memtype_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot copy HDF5 attribute ', &
            'memory datatype.'
       call h5sclose_f(space_id, hdferr)
       call h5tclose_f(type_id, hdferr)
       call h5aclose_f(attr_id, hdferr)
       fail = .true.
       return
    end if
    call h5tset_size_f(memtype_id, sdim, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot set HDF5 attribute ', &
            'memory datatype size.'
       call h5tclose_f(memtype_id, hdferr)
       call h5sclose_f(space_id, hdferr)
       call h5tclose_f(type_id, hdferr)
       call h5aclose_f(attr_id, hdferr)
       fail = .true.
       return
    end if

    ! Read the attribute
    allocate(rdata(1:dims(1)))
    f_ptr = C_LOC(rdata(1)(1:1))
    call h5aread_f(attr_id, memtype_id, f_ptr, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot read HDF5 attribute.'
       deallocate(rdata)
       call h5tclose_f(memtype_id, hdferr)
       call h5sclose_f(space_id, hdferr)
       call h5tclose_f(type_id, hdferr)
       call h5aclose_f(attr_id, hdferr)
       fail = .true.
       return
    end if

    ! Check the units
    if (trim(rdata(1)) .ne. trim(units)) then
       write(LVT_logunit,*) &
            '[ERR] Found wrong HDF5 data', &
                 'units'
       write(LVT_logunit,*) 'Expected ', trim(units)
       write(LVT_logunit,*) 'Found ',trim(rdata(1))
       deallocate(rdata)
       call h5tclose_f(memtype_id, hdferr)
       call h5sclose_f(space_id, hdferr)
       call h5tclose_f(type_id, hdferr)
       call h5aclose_f(attr_id, hdferr)
       fail = .true.
       return
    end if

    ! Clean up
    deallocate(rdata)
    call h5tclose_f(memtype_id, hdferr)
    call h5sclose_f(space_id, hdferr)
    call h5tclose_f(type_id, hdferr)
    call h5aclose_f(attr_id, hdferr)

  end subroutine check_navgem_units
#endif

#if (defined USE_HDF5)
  subroutine get_navgem_dims(dataset_id, rank, dims, fail)

    ! Modules
    use HDF5
    use LVT_logMod, only: LVT_logunit

    ! Defaults
    implicit none

    ! Arguments
    integer(HID_T), intent(in) :: dataset_id
    integer, intent(out) :: rank
    integer(HSIZE_T), allocatable, intent(out) :: dims(:)
    logical, intent(out) :: fail

    ! Local variables
    integer(HID_T) :: dataspace_id
    integer(HSIZE_T), allocatable :: dataspace_maxdims(:)
    integer :: hdferr
    logical :: flag
    integer :: i

    ! First, get the dataspace for the dataset
    call h5dget_space_f(dataset_id, dataspace_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*)&
            '[ERR] Could not get HDF5 dataspace'
       fail = .true.
       return
    end if

    ! Sanity check:  Make sure this dataspace is "simple"
    call h5sis_simple_f(dataspace_id, flag, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot determine if ', &
            'HDF5 dataspace is simple'
       fail = .true.
       return
    end if
    if (.not. flag) then
       write(LVT_logunit,*) &
            '[ERR] HDF5 dataspace is not simple'
       fail = .true.
       return
    end if

    ! Get the rank (number of dimensions)
    call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*)&
            '[ERR] Cannot get rank of HDF5 dataspace '
       fail = .true.
       return
    end if

    ! Get the dimensions
    allocate(dims(rank))
    allocate(dataspace_maxdims(rank))
    call h5sget_simple_extent_dims_f(dataspace_id, dims, &
         dataspace_maxdims, hdferr)
    if (hdferr .ne. rank) then
       write(LVT_logunit,*) &
            '[ERR] Cannot get dims for HDF5 dataspace'
       deallocate(dims)
       deallocate(dataspace_maxdims)
       fail = .true.
       return
    end if

    ! Clean up
    deallocate(dataspace_maxdims)
    call h5sclose_f(dataspace_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot close HDF5 dataspace'
       fail = .true.
       return
    end if

  end subroutine get_navgem_dims
#endif

#if (defined USE_HDF5)
  subroutine close_navgem_datatype(datatype_id, fail)
    use HDF5
    use LVT_logMod, only: LVT_logunit
    integer(HID_T), intent(inout) :: datatype_id
    logical, intent(out) :: fail
    integer :: hdferr
    fail = .false.
    call h5tclose_f(datatype_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot close HDF5 datatype '
       fail = .true.
    end if
    datatype_id = -1
  end subroutine close_navgem_datatype
#endif

#if (defined USE_HDF5)
  subroutine close_navgem_dataset(dataset_id, fail)
    use HDF5
    use LVT_logMod, only: LVT_logunit
    integer(HID_T), intent(inout) :: dataset_id
    logical, intent(out) :: fail
    integer :: hdferr
    fail = .false.
    call h5dclose_f(dataset_id, hdferr)
    if (hdferr .ne. 0) then
       write(LVT_logunit,*) &
            '[ERR] Cannot close HDF5 dataset '
       fail = .true.
    end if
    dataset_id = -1
  end subroutine close_navgem_dataset
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
            '[ERR] Cannot close NAVGEM file ', trim(filename)
       fail = .true.
    else
       write(LVT_logunit,*) &
            '[INFO] Closed NAVGEM file ',trim(filename)
    end if
    file_id = -1
  end subroutine close_navgem_file
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

  subroutine get_navgem_truncation(im, t_number)
    implicit none
    integer, intent(in) :: im
    integer, intent(out) :: t_number
    t_number = int( 2 * int( int ( int( int(im-1)/3) + 1) / 2 ) )
    t_number = t_number - 1
  end subroutine get_navgem_truncation

  ! Calculate the latitude and longitude of each data point on the reduced
  ! Gaussian grid.  This logic borrows heavily from Python code provided by
  ! FNMOC.
  subroutine calc_navgem_latlons(tmp_latitudes, tmp_points_per_lat, dim, &
       lats, lons)

    ! Defaults
    implicit none

    ! Arguments
    real, intent(in) :: tmp_latitudes(:,:) ! The latitude of each parallel
    integer, intent(in) :: tmp_points_per_lat(:,:) ! Points per parallel
    integer, intent(in) :: dim ! Total number of points
    real, allocatable, intent(out) :: lats(:)
    real, allocatable, intent(out) :: lons(:)

    ! Locals
    integer :: num_lons
    real :: d_lon
    integer :: r, c, i, jm

    jm = size(tmp_latitudes, 1) ! Number of parallels

    ! First calculate the latitudes at each point
    allocate(lats(dim))
    lats = 0
    i = 0
    do r = 1, jm
       num_lons = tmp_points_per_lat(r,1)
       do c = 1, num_lons
          i = i + 1
          lats(i) = tmp_latitudes(r,1)
       end do ! c
    end do ! r

    ! Next, calculate the longitudes at each point
    allocate(lons(dim))
    lons = 0
    i = 0
    do r = 1, jm
       num_lons = tmp_points_per_lat(r,1)
       d_lon = 360. / num_lons
       do c = 1, num_lons
          i = i + 1
          lons(i) = 0.0 + (c-1)*d_lon
          if (lons(i) .gt. 180.) then
             lons(i) = lons(i) - 360.
          end if
       end do ! c
    end do ! r
  end subroutine calc_navgem_latlons

  ! Special version of upscaleByAveraging_input.  This variant skips the
  ! internal calculation of latitudes and longitudes by compute_earth_coord,
  ! since that subroutine doesn't support thinned (quasi-regular) Gaussian
  ! grids.  Instead, the latitudes and longitudes are passed in as additional
  ! arguments.
  subroutine LVT_upscaleByAveraging_input_navgem(gridDesco, mi, rlat, rlon, &
       mo, n11)

    ! Defaults
    implicit none

    ! Arguments
    real, intent(in) :: gridDesco(50)
    integer, intent(in) :: mi
    real, intent(in) :: rlat(mi)
    real, intent(in) :: rlon(mi)
    integer, intent(in) :: mo
    integer, intent(out) :: n11(mi)

    ! Locals
    integer             :: n
    integer             :: i, j
    real                :: xi, yi
    real                :: xpts(mi), ypts(mi)
    real                :: xpts1(mo), ypts1(mo)
    real                :: rlat1(mo), rlon1(mo)
    integer             :: nv
    real, parameter     :: fill = -9999.0

    ! External functions
    integer, external   :: get_fieldpos

    ! Find the x,y coordinates of the input points on the output grid
    call compute_grid_coord(gridDesco, mi,fill, xpts, ypts, rlon, rlat, nv)

    ! Determine which grid box in the outer grid each inner point resides in.
    do n = 1, mi
       xi = xpts(n)
       yi = ypts(n)
       if (xi.ne.fill .and. yi.ne.fill) then
          i = nint(xi)
          j = nint(yi)
          n11(n) = get_fieldpos(i, j, gridDesco)
       else
          n11(n) = 0
       endif
    end do ! n
  end subroutine LVT_upscaleByAveraging_input_navgem

end module LVT_navgemMod
