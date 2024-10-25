!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LDT_misc.h"

module LDT_NAFPA_back_climpptMod

  ! DESCRIPTION:
  ! This module implements routines to read, average, and interpolate monthly
  ! precipitation values from NAFPA background fields (GFS and GALWEM),
  ! already processed by LVT.
  !
  ! REVISION HISTORY:
  ! 13 May 2022: Eric Kemp/SSAI: First version.

  ! Defaults
  implicit none
  private

  ! Class type for reading NAFPA background fields from LVT files,
  ! interpolating to LIS grid, summing, and averaging each month
  type, public :: LDT_NAFPA_back_climppt_t
     private
     character(255) :: topdir_input
     integer :: nlat_input
     integer :: nlon_input
     real :: gridDesc_input(20)
     real, allocatable :: pcp_input(:,:)
     integer :: nrows_out
     integer :: ncols_out
     real :: gridDesc_out(20)
     real, allocatable :: pcp_out(:,:)
     integer :: count
     real :: missing_value_input
     character(6) :: source ! GFS or GALWEM
   contains
     ! Object methods
     procedure :: new => LDT_NAFPA_back_climppt_new
     procedure :: delete => LDT_NAFPA_back_climppt_delete
     procedure :: process => LDT_NAFPA_back_climppt_process
  end type LDT_NAFPA_back_climppt_t

  ! Public readers that use the object under the hood
  public :: LDT_read_NAFPA_back_gfs_climppt
  public :: LDT_read_NAFPA_back_galwem_climppt

contains

  ! Public reader
  subroutine LDT_read_NAFPA_back_gfs_climppt(nest, ncols_out, nrows_out, &
       gridDesc_out, pcp_out)

    ! Imports
    use LDT_climateParmsMod, only: LDT_climate_struc
    use LDT_logMod, only: LDT_logunit

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: nest
    integer, intent(in) :: ncols_out
    integer, intent(in) :: nrows_out
    real, intent(in) :: gridDesc_out(20)
    real, intent(out) :: pcp_out(ncols_out, nrows_out)

    ! Locals
    type(LDT_nafpa_back_climppt_t) :: nafpa_back
    integer :: imonth

    call nafpa_back%new(nest, ncols_out, nrows_out, gridDesc_out, "GFS")
    imonth = LDT_climate_struc(nest)%climpptimonth
    call nafpa_back%process(nest, imonth, ncols_out, nrows_out, pcp_out, &
         "GFS")
    call nafpa_back%delete()

  end subroutine LDT_read_NAFPA_back_gfs_climppt

    ! Public reader
  subroutine LDT_read_NAFPA_back_galwem_climppt(nest, ncols_out, nrows_out, &
       gridDesc_out, pcp_out)

    ! Imports
    use LDT_climateParmsMod, only: LDT_climate_struc
    use LDT_logMod, only: LDT_logunit

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: nest
    integer, intent(in) :: ncols_out
    integer, intent(in) :: nrows_out
    real, intent(in) :: gridDesc_out(20)
    real, intent(out) :: pcp_out(ncols_out, nrows_out)

    ! Locals
    type(LDT_nafpa_back_climppt_t) :: nafpa_back
    integer :: imonth

    call nafpa_back%new(nest, ncols_out, nrows_out, gridDesc_out, "GALWEM")
    imonth = LDT_climate_struc(nest)%climpptimonth
    call nafpa_back%process(nest, imonth, ncols_out, nrows_out, pcp_out, &
         "GALWEM")
    call nafpa_back%delete()

  end subroutine LDT_read_NAFPA_back_galwem_climppt

  ! Constructor method
  subroutine LDT_NAFPA_back_climppt_new(this, nest, ncols_out, nrows_out, &
       gridDesc_out, source)

    ! Imports
    use LDT_climateParmsMod, only: LDT_climate_struc
    use LDT_logMod, only: LDT_logunit, LDT_endrun

    ! Defaults
    implicit none

    ! Arguments
    class(LDT_NAFPA_back_climppt_t), intent(inout) :: this
    integer, intent(in) :: nest
    integer, intent(in) :: ncols_out
    integer, intent(in) :: nrows_out
    real, intent(in) :: gridDesc_out(20)
    character(*), intent(in) :: source

    if (trim(source) .eq. "GALWEM") then
       this%topdir_input = trim(LDT_climate_struc(nest)%climpptdir)
    else if (trim(source) .eq. "GFS") then
       this%topdir_input = trim(LDT_climate_struc(nest)%climpptdir2)
    else
       write(LDT_logunit,*)'[ERR] Invalid source for NAFPA background!'
       call LDT_endrun()
    end if

    ! Set up NAFPA n1280 grid.
    this%nlat_input = 1920
    this%nlon_input = 2560
    this%gridDesc_input = 0
    this%gridDesc_input(1) = 0 ! Lat/lon projection
    this%gridDesc_input(2) = this%nlon_input ! Number of columns
    this%gridDesc_input(3) = this%nlat_input ! Number of rows
    this%gridDesc_input(4) =  -89.9531250 ! Lower-left latitude
    this%gridDesc_input(5) = -179.9296875 ! Lower-left longitude
    this%gridDesc_input(6) =  128             ! Not used
    this%gridDesc_input(7) =   89.9531250  ! Upper-right latitude
    this%gridDesc_input(8) =  179.9296875  ! Upper-right longitude
    this%gridDesc_input(9) =    0.1406250 ! Delta longitude
    this%gridDesc_input(10) =   0.0937500 ! Delta latitude
    this%gridDesc_input(11) =   64 ! Not used
    this%gridDesc_input(20) =  255 ! Indicates E-W ordering of data
    allocate(this%pcp_input(this%nlon_input, this%nlat_input))
    this%pcp_input = 0

    ! Set up LDT grid
    this%ncols_out = ncols_out
    this%nrows_out = nrows_out
    this%gridDesc_out = gridDesc_out
    allocate(this%pcp_out(this%ncols_out, this%nrows_out))
    this%pcp_out = 0

    ! Other initializations.  These will change as data are read, interpolated,
    ! and averaged.
    this%count = 0
    this%missing_value_input = 0
  end subroutine LDT_NAFPA_back_climppt_new

  ! Destructor method
  subroutine LDT_NAFPA_back_climppt_delete(this)

    ! Defaults
    implicit none

    ! Arguments
    class(LDT_NAFPA_back_climppt_t), intent(inout) :: this

    ! Clean up data structure
    this%topdir_input = 'NULL'
    this%nlat_input = 0
    this%nlon_input = 0
    this%gridDesc_input = 0
    deallocate(this%pcp_input)
    this%nrows_out = 0
    this%ncols_out = 0
    this%gridDesc_out = 0
    deallocate(this%pcp_out)
    this%count = 0
    this%missing_value_input = 0

  end subroutine LDT_NAFPA_back_climppt_delete

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  ! Handle processing of one month of NAFPA background precipitation data
  subroutine LDT_NAFPA_back_climppt_process(this, nest, imonth, &
       ncols_out, nrows_out, pcp_out, source)

    ! Imports
    use ESMF
    use LDT_climateParmsMod, only: LDT_climate_struc
    use LDT_coreMod, only: LDT_rc, LDT_domain
    use LDT_logMod, only: LDT_logunit, LDT_endrun
    use netcdf

    ! Defaults
    implicit none

    ! Arguments
    class(LDT_NAFPA_back_climppt_t), intent(inout) :: this
    integer, intent(in) :: nest
    integer, intent(in) :: imonth
    integer, intent(in) :: ncols_out
    integer, intent(in) :: nrows_out
    real, intent(out) :: pcp_out(ncols_out, nrows_out)
    character(*), intent(in) :: source

    ! Locals
    integer :: ipass
    integer :: iyear, iyear_start, iyear_end
    character(255) :: filename
    logical :: found_inq
    integer :: iret
    integer :: ftn, varid
    integer, allocatable :: n11(:)
    real, allocatable :: gi1(:)
    logical*1, allocatable :: li1(:)
    real, allocatable :: go1(:)
    logical*1, allocatable :: lo1(:)
    integer, allocatable :: n112(:,:)  ! Map array for conservative interp
    integer, allocatable :: n122(:,:)
    integer, allocatable :: n212(:,:)
    integer, allocatable :: n222(:,:)
    real, allocatable :: w112(:,:), w122(:,:)
    real, allocatable :: w212(:,:), w222(:,:)
    integer :: mi, mo
    integer :: i, j, ij

    external :: upscaleByAveraging_input, upscaleByAveraging
    external :: conserv_interp_input, conserv_interp

    ! Loop through each year for the selected month.  We have data for
    ! Oct 2017 through May 2022.
    ipass = 0
    if (imonth .ge. 1 .and. imonth .le. 5) then
       iyear_start = 2018
       iyear_end = 2022
    else if (imonth .ge. 6 .and. imonth .le. 9) then
       iyear_start = 2018
       iyear_end = 2021
    else
       iyear_start = 2017
       iyear_end = 2021
    end if

    do iyear = iyear_start, iyear_end

       ! Read data from LVT netCDF file.
       call create_filename(this, imonth, iyear, filename)
       inquire(file=trim(filename), exist=found_inq)
       if (.not. found_inq) then
          write(LDT_logunit,*)'[ERR] Cannot find ', trim(filename)
          write(LDT_logunit,*)'[ERR] Stopping...'
          call LDT_endrun()
       end if
       write(LDT_logunit,*)'[INFO] Reading ', trim(filename)
       ipass = ipass + 1

       iret = nf90_open(path=trim(filename), mode=NF90_NOWRITE, ncid=ftn)
       if (iret .ne. 0) then
          write(LDT_logunit,*)'[ERR] nf90_open failed for ',trim(filename)
          call LDT_endrun()
       end if

       iret = nf90_inq_varid(ftn, 'TotalPrecip', varid)
       if (iret .ne. 0) then
          write(LDT_logunit,*)'[ERR] nf90_inq_varid failed for TotalPrecip'
          call LDT_endrun()
       end if

       iret = nf90_get_var(ftn, varid, this%pcp_input)
       if (iret .ne. 0) then
          write(LDT_logunit,*)'[ERR] nf90_get_var failed for TotalPrecip'
          call LDT_endrun()
       end if

       iret = nf90_close(ftn)
       if (iret .ne. 0) then
          write(LDT_logunit,*)'[ERR] nf90_close failed'
          call LDT_endrun()
       end if

       ! Interpolate to LDT grid
       mi = this%nlon_input * this%nlat_input
       mo = this%ncols_out * this%nrows_out
       select case (LDT_climate_struc(nest)%clim_gridtransform)
       case ("average")
          !write(ldt_logunit,*) &
          !     '[INFO] Averaging NAFPA background data to LDT grid'
          if (ipass .eq. 1) then
             allocate(n11(mi)) ; n11 = 0
             call upscaleByAveraging_input(this%gridDesc_input, &
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
          do j = 1, this%nlat_input
             do i = 1, this%nlon_input
                ij = i + (j-1)*this%nlon_input
                gi1(ij) = this%pcp_input(i,j)
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

       case ("budget-bilinear")
          !write(ldt_logunit,*) &
          !     '[INFO] Budget-Bilinear Interpolating NAFPA background data to LDT grid'
          if (ipass .eq. 1) then
             allocate(gi1(mi))
             allocate(li1(mi))
             allocate(go1(mo))
             allocate(lo1(mo))
             allocate( n112(mo,25), n122(mo,25), n212(mo,25), n222(mo,25) )
             n112 = 0; n122 = 0; n212 = 0; n222 = 0
             allocate( w112(mo,25), w122(mo,25), w212(mo,25), w222(mo,25) )
             w112 = 0; w122 = 0; w212 = 0; w222 = 0
             call conserv_interp_input(nest, this%gridDesc_input, &
                  n112, n122, n212, n222, w112, w122, w212, w222 )
          end if

          gi1 = LDT_rc%udef
          li1 = .false.
          go1 = 0
          lo1 = .false.
          do j = 1, this%nlat_input
             do i = 1, this%nlon_input
                ij = i + (j-1)*this%nlon_input
                gi1(ij) = this%pcp_input(i,j)
                if (gi1(ij) .ne. LDT_rc%udef) then
                   li1(ij) = .true.
                end if
             end do
          end do

          call conserv_interp( this%gridDesc_out, li1, gi1, lo1(:), go1(:), &
               mi, mo, LDT_domain(nest)%lat, LDT_domain(nest)%lon, &
               w112, w122, w212, w222, n112, n122, n212, n222, &
               LDT_rc%udef, iret)

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

       case ("none")

          do j = 1, this%nrows_out
             do i = 1, this%ncols_out
                if (this%pcp_input(i,j) < 0.) then
                   this%pcp_out(i,j) = LDT_rc%udef
                else
                   this%pcp_out(i,j) = this%pcp_out(i,j) + this%pcp_input(i,j)
                end if
             end do
          end do

       case default
          write(ldt_logunit,*)"[ERR] Invalid spatial transform for NAFPA"
          write(ldt_logunit,*) &
               "[ERR] Please select average, budget-bilinear, or none for now"
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
    if (allocated(n112)) deallocate(n112)
    if (allocated(n122)) deallocate(n122)
    if (allocated(n212)) deallocate(n212)
    if (allocated(n222)) deallocate(n222)
    if (allocated(w112)) deallocate(w112)
    if (allocated(w122)) deallocate(w122)
    if (allocated(w212)) deallocate(w212)
    if (allocated(w222)) deallocate(w222)

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
  end subroutine LDT_NAFPA_back_climppt_process

#else

  ! Dummy subroutine if not compiled with netcdf
  subroutine LDT_NAFPA_back_climppt_process(this, nest, imonth, &
       ncols_out, nrows_out, pcp_out, gfs_or_galwem)
    use LDT_logMod, only: ldt_logunit, LDT_endrun
    implicit none
    class(LDT_NAFPA_back_climppt_t), intent(inout) :: this
    integer, intent(in) :: nest
    integer, intent(in) :: imonth
    integer, intent(in) :: ncols_out
    integer, intent(in) :: nrows_out
    real, intent(out) :: pcp_out(ncols_out, nrows_out)
    character(*), intent(in) :: gfs_or_galwem
    write(LDT_logunit,*)'[ERR] LDT was compiled without netCDF support!'
    write(LDT_logunit,*)'[ERR] Cannot process NAFPA background files!'
    write(LDT_logunit,*) &
         '[ERR] Recompile LDT with netCDF support and try again.'
    write(LDT_logunit,*)'[ERR] Stopping...'
    call LDT_endrun()
  end subroutine LDT_NAFPA_back_climppt_process
#endif

  subroutine create_filename(this, imonth, iyear, filename)

    ! Defaults
    implicit none

    ! Arguments
    class(LDT_NAFPA_back_climppt_t), intent(in) :: this
    integer, intent(in) :: imonth
    integer, intent(in) :: iyear
    character(255), intent(out) :: filename

    ! Locals
    integer :: imonth_file, iyear_file

    ! Each NAFPA file has data for the previous month.  So, for Sep 2010,
    ! we look in SUM_TS.201010010000.d01.nc.  Apply that logic here.
    iyear_file = iyear
    imonth_file = imonth + 1
    if (imonth .ge. 12) then
       imonth_file = 1
       iyear_file = iyear + 1
    end if

    write(filename,'(A,A,I4.4,I2.2,A)') &
         trim(this%topdir_input), &
         '/SUM_TS.', iyear_file, imonth_file, '010000.d01.nc'
  end subroutine create_filename
end module LDT_NAFPA_back_climpptMod
