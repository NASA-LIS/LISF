!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LDT_misc.h"

module LDT_CHELSAV21_climpptMod

  ! DESCRIPTION:
  ! This module implements routines to read, average, and interpolate monthly
  ! precipitation values from the CHELSAV21 dataset.
  !
  ! REFERENCES:
  ! Beck, H, E F Wood, T McVicar, M Zambrano-Bigiarini, C Alvarez-Garreton,
  !   O Baez-Villanueva, J Sheffield, and D N Karger, 2020: Bias correction of
  !   global high-resolution precipitation climatologies using streamflow
  !   observations from 9372 catchments. Journal of Climate, 33, 1299-1315,
  !   https://doi.org/10.1175/JCLI-D-19-0332.1.
  ! Karger, D N, and N E Zimmerman, 2021: CHELSA V2.1: Technical specification.
  !   https://chelsa-climate.org/wp-admin/download-page/\
  !   CHELSA_tech_specification_V2.pdf
  ! Karger, D N, O Conrad, J Bohner, T Kawohl, H Kreft, R W Soria-Auza,
  !   N E Zimmermann, H P Linder, and M Kessler, 2017: Climatologies at high
  !   resolution for the Earth land surface areas. Scientific Data, 4,
  !   170122, https://doi.org/10.1038/sdata.2017.122.
  ! Karger, D N, O Conrad, J Bohner, T Kawohl, H Kreft, R W Soria-Auza,
  !   N E Zimmermann, H P Linder, and M Kessler, 2021: Climatologies at high
  !   resolution for the Earth's land surface areas. EnviDat,
  !   https://doi.org/10.16904/envidat.228.v2.1.
  !
  ! REVISION HISTORY:
  ! 07 Feb 2022: Eric Kemp/SSAI: First version.

  ! Defaults
  implicit none
  private

  ! Class type for reading CHELSAV21 data from GeoTIFF file, interpolating
  ! to LIS grid, summing, and averaging each month.
  type, public :: LDT_CHELSAV21_climppt_t
     private
     character(500) :: topdir_native
     integer :: nlat_native
     integer :: nlon_native
     real :: gridDesc_native(20)
     real, allocatable :: pcp_native(:,:)
     integer :: nrows_out
     integer :: ncols_out
     real :: gridDesc_out(20)
     real, allocatable :: pcp_out(:,:)
     integer :: count
     real :: missing_value_native
   contains
     ! Object methods
     procedure :: new => LDT_CHELSAV21_climppt_new
     procedure :: delete => LDT_CHELSAV21_climppt_delete
     procedure :: process => LDT_CHELSAV21_climppt_process
  end type LDT_CHELSAV21_climppt_t

  ! Public reader that uses the object under the hood
  public :: LDT_read_CHELSAV21_climppt

contains

  ! Public reader
  subroutine LDT_read_CHELSAV21_climppt(nest, ncols_out, nrows_out, &
       gridDesc_out, pcp_out)

    ! Imports
    use LDT_climateParmsMod, only: LDT_climate_struc

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: nest
    integer, intent(in) :: ncols_out
    integer, intent(in) :: nrows_out
    real, intent(in) :: gridDesc_out(20)
    real, intent(out) :: pcp_out(ncols_out, nrows_out)

    ! Locals
    type(LDT_CHELSAV21_climppt_t) :: chelsav21
    integer :: imonth

    call chelsav21%new(nest, ncols_out, nrows_out, gridDesc_out)
    imonth = LDT_climate_struc(nest)%climpptimonth
    call chelsav21%process(nest, imonth, ncols_out, nrows_out, pcp_out)
    call chelsav21%delete()

  end subroutine LDT_read_CHELSAV21_climppt

  ! Constructor method
  subroutine LDT_CHELSAV21_climppt_new(this, nest, ncols_out, nrows_out, &
       gridDesc_out)

    ! Imports
    use LDT_climateParmsMod, only: LDT_climate_struc

    ! Defaults
    implicit none

    ! Arguments
    class(LDT_CHELSAV21_climppt_t), intent(inout) :: this
    integer, intent(in) :: nest
    integer, intent(in) :: ncols_out
    integer, intent(in) :: nrows_out
    real, intent(in) :: gridDesc_out(20)

    this%topdir_native = trim(LDT_climate_struc(nest)%climpptdir)

    ! Set up CHELSA grid.  Note:  The native grid has origin in the upper-left.
    ! The reader will flip to a more conventional grid before the
    ! interpolation; the gridDesc_native reflects the 'flipped' grid.
    this%nlat_native = 20880
    this%nlon_native = 43200
    this%gridDesc_native(:) = 0
    this%gridDesc_native(1) = 0 ! Lat/lon projection
    this%gridDesc_native(2) = this%nlon_native ! Number of columns
    this%gridDesc_native(3) = this%nlat_native ! Number of rows
    this%gridDesc_native(4) =  -89.99597222215 ! Lower-left latitude
    this%gridDesc_native(5) = -179.99597222215 ! Lower-left longitude
    this%gridDesc_native(6) =  128             ! Not used
    this%gridDesc_native(7) =   83.99569444445001  ! Upper-right latitude
    this%gridDesc_native(8) =  179.99569444444998  ! Upper-right longitude
    this%gridDesc_native(9) =    0.008333333300000 ! Delta longitude
    this%gridDesc_native(10) =   0.008333333300000 ! Delta latitude
    this%gridDesc_native(11) =   64 ! Not used
    this%gridDesc_native(20) =  255 ! Indicates E-W ordering of data
    allocate(this%pcp_native(this%nlon_native, this%nlat_native))
    this%pcp_native = 0

    ! Set up LDT grid
    this%ncols_out = ncols_out
    this%nrows_out = nrows_out
    this%gridDesc_out = gridDesc_out
    allocate(this%pcp_out(this%ncols_out, this%nrows_out))
    this%pcp_out = 0

    ! Other initializations.  These will change as data are read, interpolated,
    ! and averaged.
    this%count = 0
    this%missing_value_native = 0

  end subroutine LDT_CHELSAV21_climppt_new

  ! Destructor method
  subroutine LDT_CHELSAV21_climppt_delete(this)

    ! Defaults
    implicit none

    ! Arguments
    class(LDT_CHELSAV21_climppt_t), intent(inout) :: this

    ! Clean up data structure
    this%topdir_native = 'NULL'
    this%nlat_native = 0
    this%nlon_native = 0
    this%gridDesc_native = 0
    deallocate(this%pcp_native)
    this%nrows_out = 0
    this%ncols_out = 0
    this%gridDesc_out = 0
    deallocate(this%pcp_out)
    this%count = 0
    this%missing_value_native = 0
  end subroutine LDT_CHELSAV21_climppt_delete

#if (defined USE_GDAL)

  ! Handle processing of one month of CHELSAV21 precipitation data
  subroutine LDT_CHELSAV21_climppt_process(this, nest, imonth, &
       ncols_out, nrows_out, pcp_out)

    ! Imports
    use LDT_climateParmsMod, only: LDT_climate_struc
    use LDT_coreMod, only: LDT_rc
    use LDT_logMod, only: LDT_logunit, LDT_endrun
    use, intrinsic :: iso_c_binding
    use fortranc
    use gdal

    ! Defaults
    implicit none

    ! Arguments
    class(LDT_CHELSAV21_climppt_t), intent(inout) :: this
    integer, intent(in) :: nest
    integer, intent(in) :: imonth
    integer, intent(in) :: ncols_out
    integer, intent(in) :: nrows_out
    real, intent(out) :: pcp_out(ncols_out, nrows_out)

    ! Locals
    type(gdaldriverh) :: driver
    type(gdaldataseth) :: ds
    type(gdalrasterbandh) :: band
    real(kind=c_double) :: gt(6)
    integer(kind=c_int) :: ierr, xsize, ysize
    logical :: found_inq
    real :: dres, nodata
    integer, allocatable :: row_strip(:,:)
    integer, allocatable :: n11(:)
    real, allocatable :: gi1(:)
    logical*1, allocatable :: li1(:)
    real, allocatable :: go1(:)
    logical*1, allocatable :: lo1(:)
    character(500) :: filename
    integer :: mi, mo
    integer :: pb_success
    integer :: i, j, ij, iyear, ipass

    external :: upscaleByAveraging_input, upscaleByAveraging

    ! Loop through each year for the selected month.  Note: The last year
    ! varies for different months, and this is accounted for below.
    ipass = 0
    do iyear = 1979, 2019 ! Years in climatology
       call create_filename(this, imonth, iyear, filename)
       inquire(file=trim(filename), exist=found_inq)
       if (.not. found_inq) then
          if (iyear .gt. 2018 .and. imonth .gt. 6) exit
          write(ldt_logunit,*)'[ERR] Cannot find ', trim(filename)
          write(ldt_logunit,*)'[ERR] Stopping...'
          call LDT_endrun()
       end if
       write(ldt_logunit,*)'[INFO] Reading ', trim(filename)
       ipass = ipass + 1

       call GDALAllRegister()
       driver = gdalgetdriverbyname('Tif'//char(0))
       ds = gdalopen(trim(filename)//char(0), GA_READONLY)
       if (.not. gdalassociated(ds)) then
          write(ldt_logunit,*)'[ERR] Failed to open ', trim(filename)
          write(ldt_logunit,*)'[ERR] Stopping...'
          call LDT_endrun()
       end if
       ierr = gdalgetgeotransform(ds, gt)
       dres = abs(gt(6))
       xsize = gdalgetrasterxsize(ds)
       ysize = gdalgetrasterysize(ds)

       write(ldt_logunit,*)'[INFO] CHELSAV21 Dimensions: ', xsize, ysize
       write(ldt_logunit,*)'[INFO] CHELSAV21 Resolution (deg): ', dres

       band = gdalgetrasterband(ds, 1)
       if (.not. gdalassociated(band)) then
          write(ldt_logunit,*)'[ERR] Failed to get band from ', trim(filename)
          write(ldt_logunit,*)'[ERR] Stopping...'
          call LDT_endrun()
       end if
       nodata = gdalgetrasternodatavalue(band, pb_success)
       write(ldt_logunit,*)'[INFO] Missing data flag is ', nodata

       ! Read each row of the TIFF band, convert to kg m^-2,
       ! and flip the y-axis
       allocate(row_strip(this%nlon_native,1))
       do j = 1, this%nlat_native
          ierr = gdalrasterio_f(band, GF_READ, 0, j-1, row_strip)
          if (ierr .ne. 0) then
             write(ldt_logunit,*)'[ERR] Failed to read data from ', &
                  trim(filename)
             write(ldt_logunit,*)'[ERR] Stopping...'
             call LDT_endrun()
          end if
          do i = 1, this%nlon_native
             if (row_strip(i,1) == nodata) then
                this%pcp_native(i,this%nlat_native - j + 1) = LDT_rc%udef
             else
                this%pcp_native(i,this%nlat_native - j + 1) = &
                     real(row_strip(i,1)) * 0.01 ! Convert to kg m^-2
             end if
          end do
       end do
       deallocate(row_strip)
       call gdalclose(ds)

       ! Screen out bad region in Antarctica
       do j = 1, 3
          do i = 16949, 17232
             this%pcp_native(i,j) = LDT_rc%udef
          end do
       end do

       ! Interpolate to LDT grid
       mi = this%nlon_native * this%nlat_native
       mo = this%ncols_out * this%nrows_out
       select case (LDT_climate_struc(nest)%clim_gridtransform)
       case ("average")
          if (ipass .eq. 1) then
             allocate(n11(mi)) ; n11 = 0
             write(ldt_logunit,*)'[INFO] Averaging CHELSA data to LDT grid'
             call upscaleByAveraging_input(this%gridDesc_native, &
                  this%gridDesc_out, mi, mo, n11)
             allocate(gi1(mi))
             allocate(li1(mi))
             allocate(go1(mo))
             allocate(lo1(mo))
          end if

          gi1 = LDT_rc%udef
          li1 = .false.
          go1 = 0
          lo1 = .false.
          do j = 1, this%nlat_native
             do i = 1, this%nlon_native
                ij = i + (j-1)*this%nlon_native
                gi1(ij) = this%pcp_native(i,j)
                if (gi1(ij) .ne. LDT_rc%udef) then
                   li1(ij) = .true.
                end if
             end do
          end do

          call upscaleByAveraging(mi, mo, LDT_rc%udef, n11, &
               li1, gi1, lo1, go1)

          do j = 1, this%nrows_out
             do i = 1, this%ncols_out
                ij = i + (j-1)*this%ncols_out
                if (go1(ij) < 0.) then
                   this%pcp_out(i,j) = LDT_rc%udef
                else
                   this%pcp_out(i,j) = this%pcp_out(i,j) + go1(ij)
                end if
             end do
          end do

       case default
          write(ldt_logunit,*)"[ERR] Invalid spatial transform for CHELSA21"
          write(ldt_logunit,*)"[ERR] Please select 'average' for now"
          write(ldt_logunit,*)"[ERR] Found ", &
               trim(LDT_climate_struc(nest)%clim_gridtransform)
          call LDT_endrun()
       end select

    end do ! iyear

    ! Clean up
    if (allocated(n11)) deallocate(n11)
    if (allocated(li1)) deallocate(li1)
    if (allocated(gi1)) deallocate(gi1)
    if (allocated(lo1)) deallocate(lo1)
    if (allocated(go1)) deallocate(go1)

    ! Calculate average value
    do j = 1, this%nrows_out
       do i = 1, this%ncols_out
          if (this%pcp_out(i,j) .ne. LDT_rc%udef) then
             this%pcp_out(i,j) = &
                  this%pcp_out(i,j) / real(ipass)
          end if
       end do
    end do

    do j = 1, this%nrows_out
       do i = 1, this%ncols_out
          pcp_out(i,j) = this%pcp_out(i,j)
       end do
    end do
  end subroutine LDT_CHELSAV21_climppt_process

#else

  ! Dummy subroutine if not compiled with GDAL
  subroutine LDT_CHELSAV21_climppt_process(this, nest, imonth, &
       ncols_out, nrows_out, pcp_out)
    use LDT_logMod, only: ldt_logunit, LDT_endrun
    implicit none
    class(LDT_CHELSAV21_climppt_t), intent(inout) :: this
    integer, intent(in) :: nest
    integer, intent(in) :: imonth
    integer, intent(in) :: ncols_out
    integer, intent(in) :: nrows_out
    real, intent(out) :: pcp_out(ncols_out,nrows_out)
    write(ldt_logunit,*)'[ERR] LDT was compiled without GDAL support'
    write(ldt_logunit,*)'[ERR] Cannot process CHELSA21 GeoTIFF files!'
    write(ldt_logunit,*)'[ERR] Recompile LDT with GDAL support and try again.'
    write(ldt_logunit,*)'[ERR] Stopping...'
    call LDT_endrun()
  end subroutine LDT_CHELSAV21_climppt_process
#endif

  subroutine create_filename(this, imonth, iyear, filename)
    implicit none
    class(LDT_CHELSAV21_climppt_t), intent(in) :: this
    integer, intent(in) :: imonth
    integer, intent(in) :: iyear
    character(500), intent(out) :: filename
    write(filename,'(A,A,I2.2,A,I4.4,A)') &
         trim(this%topdir_native), &
         '/CHELSA_pr_', imonth, '_', iyear, '_V.2.1.tif'
  end subroutine create_filename

end module LDT_CHELSAV21_climpptMod
