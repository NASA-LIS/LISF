!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: LVT_USAFSIpostMod
!
! REVISION HISTORY:
! 13 May 2019  Eric Kemp  Initial version
! 13 Dec 2019  Eric Kemp  Changed to USAFSI.
! 09 Oct 2020  Eric Kemp  Added legacy SNODEP files.
! 14 Jul 2021  Eric Kemp  Fixed bug in creating GRIB output directory.
! 26 Jul 2022  Eric Kemp  Corrected GRIB2 output name.
!
! DESCRIPTION:
! Source code for reading USAFSI netCDF file, writing back out in GRIB,
! interpolating to predefined Air Force grids, and outputting those grids in
! GRIB1 files.
!------------------------------------------------------------------------------

#include "LVT_misc.h"
#include "LVT_NetCDF_inc.h"

module LVT_USAFSIpostMod

   implicit none
   private

   type, public :: LVT_USAFSIpost_t
      private
      character(len=255) :: input_nc_file
      integer :: nc
      integer :: nr
      real, allocatable :: snoanl(:,:)
      real, allocatable :: snoage(:,:)
      real, allocatable :: icecon(:,:)
      real, allocatable :: icemask(:,:)
      real, allocatable :: iceage(:,:)
      character(len=255) :: map_projection
      real :: sw_corner_lat
      real :: sw_corner_lon
      real :: ne_corner_lat
      real :: ne_corner_lon
      real :: dx
      real :: dy

   contains

      ! Object methods
      procedure new
      procedure delete
      procedure read_usafsi_ncfile
      procedure output_grib2
      procedure interp_and_output_grib1

   end type LVT_USAFSIpost_t

   character(len=15), parameter, public :: GLOBAL_LL0P25 = 'global_ll0p25'
   character(len=15), parameter, public :: NH_PS16 = 'nh_ps16'
   character(len=15), parameter, public :: SH_PS16 = 'sh_ps16'

contains

   ! Constructor for LVT_USAFSIpost_t object
   subroutine new(this)

      ! Imports
      use LVT_coreMod, only: LVT_rc
      use LVT_logMod, only: LVT_logunit, LVT_endrun

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_USAFSIpost_t), intent(inout) :: this

      ! Local variables
      logical :: file_exists

      ! Construct the input netCDF filename
      this%input_nc_file = &
           trim(LVT_rc%input_dir) // &
           '/' // &
           trim(LVT_rc%input_prefix) // &
           trim(LVT_rc%yyyymmddhh) // &
           ".nc"

      ! See if the file exists
      inquire(file=trim(this%input_nc_file), exist=file_exists)
      if (.not. file_exists) then
         write(LVT_logunit,*)'[ERR] Cannot find file ',trim(this%input_nc_file)
         write(LVT_logunit,*)'[ERR] LVT will exit gracefully'
         call LVT_endrun()
      end if

      ! Set dummy values for scalars.  These will be updated later.
      this%nc = 0
      this%nr = 0
      this%map_projection = "NULL"
      this%sw_corner_lat = 0
      this%sw_corner_lon = 0
      this%ne_corner_lat = 0
      this%ne_corner_lon = 0
      this%dx = 0
      this%dy = 0

   end subroutine new

   ! Destructor of LVT_USAFSIpost_t object
   subroutine delete(this)

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_USAFSIpost_t), intent(inout) :: this

      this%input_nc_file = "NULL"
      this%nc = 0
      this%nr = 0
      if (allocated(this%snoanl)) deallocate(this%snoanl)
      if (allocated(this%snoage)) deallocate(this%snoage)
      if (allocated(this%icecon)) deallocate(this%icecon)
      if (allocated(this%icemask)) deallocate(this%icemask)
      if (allocated(this%iceage)) deallocate(this%iceage)
      this%map_projection = "NULL"
      this%sw_corner_lat = 0
      this%sw_corner_lon = 0
      this%ne_corner_lat = 0
      this%ne_corner_lon = 0
      this%dx = 0
      this%dy = 0

   end subroutine delete

   ! Method for reading USAFSI netCDF file.
   subroutine read_usafsi_ncfile(this)

      ! Imports
      use LVT_coreMod, only: LVT_rc
      use LVT_logMod, only: LVT_logunit, LVT_endrun, LVT_verify
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
      use netcdf
#endif

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_USAFSIpost_t), intent(inout) :: this

      ! Local variables
      integer :: ncid
      integer :: dim_ids(3)
      integer :: ntime, nlat, nlon
      real, allocatable :: tmp(:,:,:)
      integer :: snoanl_varid, snoage_varid, icecon_varid, icemask_varid, &
           iceage_varid
      integer :: c,r

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
      write(LVT_logunit,*)'[INFO] Reading USAFSI file ', &
           trim(this%input_nc_file)

      ! Open the file for reading
      call LVT_verify(nf90_open(path=trim(this%input_nc_file), &
           mode=NF90_NOWRITE, &
           ncid=ncid), &
           '[ERR] Error in nf90_open for ' // trim(this%input_nc_file))

      ! Read the map projection
      call LVT_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name='MAP_PROJECTION', &
           values=this%map_projection), &
           '[ERR] nf90_get_att failed for MAP_PROJECTION')

      ! Sanity check map projection
      ! TODO:  Support other map projections
      if (trim(this%map_projection) .ne. 'EQUIDISTANT CYLINDRICAL') then
         write(LVT_logunit,*) &
              '[ERR] Unrecognized map projection found in USAFSI file!'
         write(LVT_logunit,*) '[ERR] Expected EQUIDISTANT CYLINDRICAL'
         write(LVT_logunit,*) '[ERR] Found ',trim(this%map_projection)
         write(LVT_logunit,*) '[ERR] LVT will exit gracefully'
         call LVT_endrun()
      end if

      ! Get the dimension IDs
      ! TODO:  Support projections other than lat/lon
      call LVT_verify(nf90_inq_dimid(ncid=ncid, &
           name='time', &
           dimid=dim_ids(3)), &
           '[ERR] Error in nf90_inq_dimid for dimension time')
      call LVT_verify(nf90_inq_dimid(ncid=ncid, &
           name='lat', &
           dimid=dim_ids(2)), &
           '[ERR] Error in nf90_inq_dimid for dimension lat')
      call LVT_verify(nf90_inq_dimid(ncid=ncid, &
           name='lon', &
           dimid=dim_ids(1)), &
           '[ERR] Error in nf90_inq_dimid for dimension lon')

      ! Get the actual dimension sizes
      call LVT_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(3), &
           len=ntime), &
           '[ERR] Error in nf90_inquire_dimension for dimension time')
      call LVT_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(2), &
           len=nlat), &
           '[ERR] Error in nf90_inquire_dimension for dimension lat')
      call LVT_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(1), &
           len=nlon), &
           '[ERR] Error in nf90_inquire_dimension for dimension lon')

      ! Sanity check the ntime length -- should be one
      if (ntime .ne. 1) then
         write(LVT_logunit,*)'[ERR] Expected time dimension to be 1'
         write(LVT_logunit,*)'[ERR] Found ',ntime
         write(LVT_logunit,*)'[ERR] LVT will exit gracefully'
         call LVT_endrun()
      end if

      this%nc = nlon
      this%nr = nlat

      ! Fetch the USAFSI analysis variable IDs
      call LVT_verify(nf90_inq_varid(ncid=ncid, &
           name='snoanl', &
           varid=snoanl_varid), &
           '[ERR] Error in nf90_inq_varid for snoanl')
      call LVT_verify(nf90_inq_varid(ncid=ncid, &
           name='snoage', &
           varid=snoage_varid), &
           '[ERR] Error in nf90_inq_varid for snoanl')
      call LVT_verify(nf90_inq_varid(ncid=ncid, &
           name='icecon', &
           varid=icecon_varid), &
           '[ERR] Error in nf90_inq_varid for snoanl')
      call LVT_verify(nf90_inq_varid(ncid=ncid, &
           name='icemask', &
           varid=icemask_varid), &
           '[ERR] Error in nf90_inq_varid for snoanl')
      call LVT_verify(nf90_inq_varid(ncid=ncid, &
           name='iceage', &
           varid=iceage_varid), &
           '[ERR] Error in nf90_inq_varid for snoanl')

      ! Read the USAFSI variables
      allocate(tmp(nlon, nlat, ntime)) ! Need 3D array
      call LVT_verify(nf90_get_var(ncid=ncid, &
           varid=snoanl_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for snoanl')
      allocate(this%snoanl(nlon,nlat))
      do r = 1, nlat
         do c = 1, nlon
            if (tmp(c,r,1) < 0) then
               this%snoanl(c,r) = LVT_rc%udef
            else
               this%snoanl(c,r) = tmp(c,r,1)
            end if
         end do
      end do

      call LVT_verify(nf90_get_var(ncid=ncid, &
           varid=snoage_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for snoage')
      allocate(this%snoage(nlon,nlat))
      do r = 1, nlat
         do c = 1, nlon
            if (tmp(c,r,1) < 0) then
               this%snoage(c,r) = LVT_rc%udef
            else
               this%snoage(c,r) = tmp(c,r,1)
            end if
         end do
      end do

      call LVT_verify(nf90_get_var(ncid=ncid, &
           varid=icecon_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for icecon')
      allocate(this%icecon(nlon,nlat))
      do r = 1, nlat
         do c = 1, nlon
            if (tmp(c,r,1) < 0) then
               this%icecon(c,r) = LVT_rc%udef
            else
               this%icecon(c,r) = tmp(c,r,1)
            end if
         end do
      end do

      call LVT_verify(nf90_get_var(ncid=ncid, &
           varid=icemask_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for icemask')
      allocate(this%icemask(nlon,nlat))
      do r = 1, nlat
         do c = 1, nlon
            if (tmp(c,r,1) < 0) then
               this%icemask(c,r) = LVT_rc%udef
            else
               this%icemask(c,r) = tmp(c,r,1)
            end if
         end do
      end do

      call LVT_verify(nf90_get_var(ncid=ncid, &
           varid=iceage_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for icemask')
      allocate(this%iceage(nlon,nlat))
      do r = 1, nlat
         do c = 1, nlon
            if (tmp(c,r,1) < 0) then
               this%iceage(c,r) = LVT_rc%udef
            else
               this%iceage(c,r) = tmp(c,r,1)
            end if
         end do
      end do

      ! Clean up
      deallocate(tmp)

      ! Read additional global attributes
      call LVT_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name='SOUTH_WEST_CORNER_LAT', &
           values=this%sw_corner_lat), &
           '[ERR] nf90_get_att failed for SOUTH_WEST_CORNER_LAT')
      call LVT_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name='SOUTH_WEST_CORNER_LON', &
           values=this%sw_corner_lon), &
           '[ERR] nf90_get_att failed for SOUTH_WEST_CORNER_LON')
      call LVT_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name='NORTH_EAST_CORNER_LAT', &
           values=this%ne_corner_lat), &
           '[ERR] nf90_get_att failed for NORTH_EAST_CORNER_LAT')
      call LVT_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name='NORTH_EAST_CORNER_LON', &
           values=this%ne_corner_lon), &
           '[ERR] nf90_get_att failed for NORTH_EAST_CORNER_LON')
      call LVT_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name='DX', &
           values=this%dx), &
           '[ERR] nf90_get_att failed for DX')
      call LVT_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name='DY', &
           values=this%dy), &
           '[ERR] nf90_get_att failed for DY')

      ! Close the file
      call LVT_verify(nf90_close(ncid), &
           '[ERR] Error in nf90_close for ' // trim(this%input_nc_file))

#else
      write(LVT_logunit,*) &
           '[ERR] Must compile LVT with netCDF support for USAFSIpost mode!'
      write(LVT_logunit,*) &
           '[ERR] LVT will exit gracefully.'
      call LVT_endrun()
#endif

   end subroutine read_usafsi_ncfile

   ! Output USAFSI data in GRIB2 format.
   subroutine output_grib2(this)

      ! Imports
      use grib_api
      use LVT_coreMod, only: LVT_rc
      use LVT_logMod, only: LVT_logunit, LVT_endrun

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_USAFSIpost_t), intent(inout) :: this

      ! Local variables
      real :: griddesci(50)
      character(len=255) :: fname
      integer :: ftn, rc, status2
      character(len=255) :: msg
      real, allocatable :: go(:)
      logical :: found
      integer :: c, r

      integer, external :: LVT_create_subdirs

      ! Make sure output directory exists
      inquire(file=trim(LVT_rc%output_dir), &
           exist=found)
      if (.not. found) then
         rc = LVT_create_subdirs(len_trim(LVT_rc%output_dir), &
              trim(LVT_rc%output_dir))
         if (rc .ne. 0) then
            write(LVT_logunit,*)'[ERR] Cannot create directory ', &
                 trim(LVT_rc%output_dir)
            write(LVT_logunit,*)'[ERR] Program will stop'
            call LVT_endrun()
         end if
      end if

      ! Set the USAFSI grid description
      call set_griddesci(this, griddesci)

      ! Construct the GRIB2 filename
      call build_filename_g2(LVT_rc%output_dir, &
           LVT_rc%yyyymmddhh, fname)

      ! Open the GRIB2 file
      call grib_open_file(ftn, fname, 'w', rc)
      if (rc .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_open_file'
         call grib_get_error_string(rc, msg, status2)
         write(LVT_logunit,*)'[ERR] ', trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      else
         write(LVT_logunit,*)'[INFO] Writing to ', trim(fname)
      end if

      allocate(go(this%nc * this%nr))

      ! Handle snoanl
      do r = 1, this%nr
         do c = 1, this%nc
            go(c + (r-1)*this%nc) = this%snoanl(c,r)
         end do ! c
      end do ! r
      call write_grib2(ftn, griddesci, this%nc, this%nr, go, &
           dspln=0, cat=1, num=11, typegenproc=12, fcsttime=0)

      ! Handle snoage
      do r = 1, this%nr
         do c = 1, this%nc
            go(c + (r-1)*this%nc) = this%snoage(c,r)
         end do ! c
      end do ! r
      call write_grib2(ftn, griddesci, this%nc, this%nr, go, &
           dspln=0, cat=1, num=17, typegenproc=12, fcsttime=0)

      ! Handle icecon
      do r = 1, this%nr
         do c = 1, this%nc
            go(c + (r-1)*this%nc) = this%icecon(c,r)
         end do ! c
      end do ! r
      call write_grib2(ftn, griddesci, this%nc, this%nr, go, &
           dspln=10, cat=2, num=0, typegenproc=12, fcsttime=0)

      ! Handle icemask
      do r = 1, this%nr
         do c = 1, this%nc
            go(c + (r-1)*this%nc) = this%icemask(c,r)
         end do ! c
      end do ! r
      ! FIXME:  Find parameter number for ice mask
      call write_grib2(ftn, griddesci, this%nc, this%nr, go, &
           dspln=10, cat=2, num=192, typegenproc=12, fcsttime=0)

      ! Handle iceage
      do r = 1, this%nr
         do c = 1, this%nc
            go(c + (r-1)*this%nc) = this%iceage(c,r)
         end do ! c
      end do ! r
      ! FIXME:  Find parameter number for ice age
      call write_grib2(ftn, griddesci, this%nc, this%nr, go, &
           dspln=10, cat=2, num=193, typegenproc=12, fcsttime=0)

      ! Close the GRIB2 file
      call grib_close_file(ftn, rc)
      if (rc .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_close_file'
         call grib_get_error_string(rc, msg, status2)
         write(LVT_logunit,*)'[ERR] ', trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if

      ! Clean up
      deallocate(go)
   end subroutine output_grib2

   ! Interpolate and output USAFSI data in GRIB1 format.
   subroutine interp_and_output_grib1(this, gridID)

      ! Imports
      use grib_api
      use LVT_coreMod, only: LVT_rc
      use LVT_logMod, only: LVT_logunit, LVT_endrun

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_USAFSIpost_t), intent(inout) :: this
      character(len=*), intent(in) :: gridID

      ! Local variables
      real :: griddesci(50), griddesco(50)
      real :: xmesh, xpnmcaf, ypnmcaf, orient, xj, xi, alat, alon
      integer :: nc_out, nr_out
      integer, allocatable :: n11(:)
      real, allocatable :: rlat_bin(:)
      real, allocatable :: rlon_bin(:)
      integer, allocatable :: n11_bin(:)
      integer, allocatable :: n12_bin(:)
      integer, allocatable :: n21_bin(:)
      integer, allocatable :: n22_bin(:)
      real, allocatable :: w11_bin(:)
      real, allocatable :: w12_bin(:)
      real, allocatable :: w21_bin(:)
      real, allocatable :: w22_bin(:)
      integer, allocatable :: n11_neighbor(:)
      real, allocatable :: rlat_neighbor(:)
      real, allocatable :: rlon_neighbor(:)
      logical*1, allocatable :: li(:), lo(:), lo_bin(:), lo_neighbor(:)
      real, allocatable :: gi(:), go(:), go_bin(:), go_neighbor(:)
      real, allocatable :: go2d(:,:)
      character(len=255) :: fname
      integer :: ftn, rc, status2, iret
      integer :: c,r
      integer :: igrib
      character(len=8) :: yyyymmdd
      character(len=4) :: hhmm
      integer :: idate8, idate4
      integer :: gridDefinitionTemplateNumber
      character(len=255) :: msg
      integer :: grid_definition
      logical :: write_fullgrib_file
      logical :: write_snodep_file
      logical :: found

      integer, external :: LVT_create_subdirs

      ! Make sure output directory exists
      inquire(file=trim(LVT_rc%output_dir), &
           exist=found)
      if (.not. found) then
         rc = LVT_create_subdirs(len_trim(LVT_rc%output_dir), &
              trim(LVT_rc%output_dir))
         if (rc .ne. 0) then
            write(LVT_logunit,*)'[ERR] Cannot create directory ', &
                 trim(LVT_rc%output_dir)
            write(LVT_logunit,*)'[ERR] Program will stop'
            call LVT_endrun()
         end if
      end if

      ! Sanity check the gridID.
      call check_gridID(gridID)

      ! Set the grid definition number
      grid_definition = set_grid_definition(gridID)

      ! Set the USAFSI grid description
      call set_griddesci(this, griddesci)

      ! Set the output grid description
      if (trim(gridID) .eq. trim(GLOBAL_LL0P25)) then
         call set_griddesco_global_ll0p25(this, griddesco)
      else if (trim(gridID) .eq. trim(NH_PS16)) then
         call set_griddesco_nh_ps16(this, griddesco)
      else if (trim(gridID) .eq. trim(SH_PS16)) then
         call set_griddesco_sh_ps16(this, griddesco)
      end if

      ! Useful grid variables for later
      nc_out = griddesco(2)
      nr_out = griddesco(3)
      if (griddesco(1) == 0) then
         gridDefinitionTemplateNumber = 0
      else if (griddesco(1) == 5) then
         gridDefinitionTemplateNumber = 20
      end if

      ! Calculate weights for upscaling (averaging or mode)
      allocate(n11(this%nc*this%nr))
      call upscaleByAveraging_input(griddesci, griddesco, &
           (this%nc*this%nr), (nc_out*nr_out), n11)

      ! Calculate weights for bilinear interpolation
      if (griddesco(1) == 5) then
         allocate(rlat_bin(nc_out*nr_out))
         allocate(rlon_bin(nc_out*nr_out))
         if (trim(gridID) .eq. NH_PS16) then
            do r = 1, nr_out
               do c = 1, nc_out
                  call pstoll(1, 1, float(c), float(nr_out - r + 1), 16, &
                       alat, alon)
                  if (alon > 180.) then
                     alon = alon - 360.
                  end if
                  rlat_bin(c + (r-1)*nc_out) = alat
                  rlon_bin(c + (r-1)*nc_out) = alon
               end do ! c
            end do ! r
         else if (trim(gridID) .eq. SH_PS16) then
            do r = 1, nr_out
               do c = 1, nc_out
                  call pstoll(2, 1, float(c), float(nr_out - r + 1), 16, &
                       alat, alon)
                  if (alon > 180.) then
                     alon = alon - 360.
                  end if
                  rlat_bin(c + (r-1)*nc_out) = alat
                  rlon_bin(c + (r-1)*nc_out) = alon
               end do ! c
            end do ! r
         end if
         allocate(n11_bin(nc_out*nr_out))
         allocate(n12_bin(nc_out*nr_out))
         allocate(n21_bin(nc_out*nr_out))
         allocate(n22_bin(nc_out*nr_out))
         allocate(w11_bin(nc_out*nr_out))
         allocate(w12_bin(nc_out*nr_out))
         allocate(w21_bin(nc_out*nr_out))
         allocate(w22_bin(nc_out*nr_out))
         call bilinear_interp_input_usaf(griddesci, griddesco, &
              (nc_out*nr_out), &
              rlat_bin, rlon_bin, &
              n11_bin, n12_bin, n21_bin, n22_bin, &
              w11_bin, w12_bin, w21_bin, w22_bin, gridID)
      end if

      ! Calculate weights for neighbor interpolation
      if (griddesco(1) == 5) then
         allocate(rlat_neighbor(nc_out*nr_out))
         allocate(rlon_neighbor(nc_out*nr_out))
         if (trim(gridID) .eq. NH_PS16) then
            do r = 1, nr_out
               do c = 1, nc_out
                  call pstoll(1, 1, float(c), float(nr_out - r + 1), 16, &
                       alat, alon)
                  if (alon > 180.) then
                     alon = alon - 360.
                  end if
                  rlat_neighbor(c + (r-1)*nc_out) = alat
                  rlon_neighbor(c + (r-1)*nc_out) = alon
               end do ! c
            end do ! r
         else if (trim(gridID) .eq. SH_PS16) then
            do r = 1, nr_out
               do c = 1, nc_out
                  call pstoll(2, 1, float(c), float(nr_out - r + 1), 16, &
                       alat, alon)
                  if (alon > 180.) then
                     alon = alon - 360.
                  end if
                  rlat_neighbor(c + (r-1)*nc_out) = alat
                  rlon_neighbor(c + (r-1)*nc_out) = alon
               end do ! c
            end do ! r
         end if
         allocate(n11_neighbor(nc_out*nr_out))
         call neighbor_interp_input_usaf(griddesci, griddesco, &
              (nc_out*nr_out), &
              rlat_neighbor, rlon_neighbor, &
              n11_neighbor, gridID)
      end if

      ! Allocate memory for interpolation
      allocate(li(this%nc*this%nr))
      allocate(gi(this%nc*this%nr))
      allocate(lo(nc_out*nr_out))
      allocate(go(nc_out*nr_out))
      if (griddesco(1) == 5) then
         allocate(go2d(nc_out, nr_out))
         allocate(go_bin(nc_out*nr_out))
         allocate(lo_bin(nc_out*nr_out))
         allocate(go_neighbor(nc_out*nr_out))
         allocate(lo_neighbor(nc_out*nr_out))
      end if

      ! Only create full grib1 file in certain cases.
      write_fullgrib_file = .false.
      if (trim(gridID) .eq. GLOBAL_LL0P25 .or. &
           (trim(gridID) .eq. NH_PS16 .and. LVT_rc%output_nh_ps16) .or. &
           (trim(gridID) .eq. SH_PS16 .and. LVT_rc%output_sh_ps16)) then
         write_fullgrib_file = .true.
      end if

      if (write_fullgrib_file) then

         ! Construct the GRIB1 filename
         call build_filename_g1(gridID, LVT_rc%output_dir, &
              LVT_rc%yyyymmddhh, fname)

         ! Open the GRIB1 file
         call grib_open_file(ftn, fname, 'w', rc)
         if (rc .ne. GRIB_SUCCESS) then
            write(LVT_logunit,*)'[ERR] Error from grib_open_file'
            call grib_get_error_string(rc, msg, status2)
            write(LVT_logunit,*)'[ERR] ', trim(msg)
            write(LVT_logunit,*)'[ERR] LVT will stop'
            call LVT_endrun()
         else
            write(LVT_logunit,*)'[INFO] Writing to ', trim(fname)
         end if

         ! Interpolate snoanl
         do r = 1, this%nr
            do c = 1, this%nc
               if (this%snoanl(c,r) < 0) then
                  li(c + (r-1)*this%nc) = .false.
                  gi(c + (r-1)*this%nc) = LVT_rc%udef
               else
                  li(c + (r-1)*this%nc) = .true.
                  gi(c + (r-1)*this%nc) = this%snoanl(c,r)
               end if
            end do ! c
         end do ! r
         if (griddesco(1) == 0) then
            call upscaleByAveraging((this%nc*this%nr), &
                 (nc_out*nr_out), LVT_rc%udef, n11, li, gi, lo, go)
         else if (griddesco(1) == 5) then
            call bilinear_interp(griddesco, li, gi, lo_bin, go_bin, &
                 (this%nc*this%nr), (nc_out*nr_out), &
                 rlat_bin, rlon_bin, &
                 w11_bin, w12_bin, w21_bin, w22_bin, &
                 n11_bin, n12_bin, n21_bin, n22_bin, &
                 LVT_rc%udef, iret)
            do r = 1, nr_out
               do c = 1, nc_out
                  go(c + (r-1)*nc_out) = go_bin(c + (r-1)*nc_out)
               end do
            end do
            ! Filter out points that are outside of the hemisphere.
            if (trim(gridID) .eq. NH_PS16) then
               do r = 1, nr_out
                  do c = 1, nc_out
                     if (rlat_bin(c + (r-1)*nc_out) < 0) then
                        go(c + (r-1)*nc_out) = LVT_rc%udef
                     end if
                  end do ! c
               end do ! r
            else if (trim(gridID) .eq. SH_PS16) then
               do r = 1, nr_out
                  do c = 1, nc_out
                     if (rlat_bin(c + (r-1)*nc_out) > 0) then
                        go(c + (r-1)*nc_out) = LVT_rc%udef
                     end if
                  end do ! c
               end do ! r
            end if
            ! If using Air Force polar stereographic, we must flip the grid so
            ! the origin is in the upper-left corner instead of lower-left
            do r = 1, nr_out
               do c = 1, nc_out
                  go2d(c,nr_out - r + 1) = go(c + (r-1)*nc_out)
               end do ! c
            end do ! r
            do r = 1, nr_out
               do c = 1, nc_out
                  go(c + (r-1)*nc_out) = go2d(c,r)
               end do ! c
            end do ! r
         end if
         call write_grib1(ftn, griddesco, nc_out, nr_out, go, param=66, &
              decimal_scale_factor=2, bits_per_value=8, &
              grid_definition=grid_definition)

         ! Interpolate snoage
         do r = 1, this%nr
            do c = 1, this%nc
               if (this%snoage(c,r) < 0) then
                  li(c + (r-1)*this%nc) = .false.
                  gi(c + (r-1)*this%nc) = LVT_rc%udef
               else
                  li(c + (r-1)*this%nc) = .true.
                  gi(c + (r-1)*this%nc) = this%snoage(c,r)
               end if
            end do ! c
         end do ! r
         if (griddesco(1) == 0) then
            call upscaleByMode((this%nc*this%nr), &
                 (nc_out*nr_out), LVT_rc%udef, n11, li, gi, lo, go)
         else if (griddesco(1) == 5) then
            call neighbor_interp(griddesco, li, gi, lo_neighbor, go_neighbor, &
                 (this%nc*this%nr), (nc_out*nr_out), &
                 rlat_neighbor, rlon_neighbor, &
                 n11_neighbor, &
                 LVT_rc%udef, iret)
            do r = 1, nr_out
               do c = 1, nc_out
                  go(c + (r-1)*nc_out) = go_neighbor(c + (r-1)*nc_out)
               end do
            end do
            ! Filter out points that are outside of the hemisphere.
            if (trim(gridID) .eq. NH_PS16) then
               do r = 1, nr_out
                  do c = 1, nc_out
                     if (rlat_bin(c + (r-1)*nc_out) < 0) then
                        go(c + (r-1)*nc_out) = LVT_rc%udef
                     end if
                  end do ! c
               end do ! r
            else if (trim(gridID) .eq. SH_PS16) then
               do r = 1, nr_out
                  do c = 1, nc_out
                     if (rlat_bin(c + (r-1)*nc_out) > 0) then
                        go(c + (r-1)*nc_out) = LVT_rc%udef
                     end if
                  end do ! c
               end do ! r
            end if
            ! If using Air Force polar stereographic, we must flip the grid so
            ! the origin is in the upper-left corner instead of lower-left
            do r = 1, nr_out
               do c = 1, nc_out
                  go2d(c,nr_out - r + 1) = go(c + (r-1)*nc_out)
               end do ! c
            end do ! r
            do r = 1, nr_out
               do c = 1, nc_out
                  go(c + (r-1)*nc_out) = go2d(c,r)
               end do ! c
            end do ! r
         end if
         call write_grib1(ftn, griddesco, nc_out, nr_out, go, param=175, &
              decimal_scale_factor=0, bits_per_value=7, &
              grid_definition=grid_definition)

         ! Handle icecon
         do r = 1, this%nr
            do c = 1, this%nc
               if (this%icecon(c,r) < 0) then
                  li(c + (r-1)*this%nc) = .false.
                  gi(c + (r-1)*this%nc) = LVT_rc%udef
               else
                  li(c + (r-1)*this%nc) = .true.
                  gi(c + (r-1)*this%nc) = this%icecon(c,r)*100 ! GRIB1 is in %
               end if
            end do ! c
         end do ! r
         if (griddesco(1) == 0) then
            call upscaleByAveraging((this%nc*this%nr), &
                 (nc_out*nr_out), LVT_rc%udef, n11, li, gi, lo, go)
         else if (griddesco(1) == 5) then
            call bilinear_interp(griddesco, li, gi, lo_bin, go_bin, &
                 (this%nc*this%nr), (nc_out*nr_out), &
                 rlat_bin, rlon_bin, &
                 w11_bin, w12_bin, w21_bin, w22_bin, &
                 n11_bin, n12_bin, n21_bin, n22_bin, &
                 LVT_rc%udef, iret)
            do r = 1, nr_out
               do c = 1, nc_out
                  go(c + (r-1)*nc_out) = go_bin(c + (r-1)*nc_out)
               end do
            end do
            ! Filter out points that are outside of the hemisphere.
            if (trim(gridID) .eq. NH_PS16) then
               do r = 1, nr_out
                  do c = 1, nc_out
                     if (rlat_bin(c + (r-1)*nc_out) < 0) then
                        go(c + (r-1)*nc_out) = LVT_rc%udef
                     end if
                  end do ! c
               end do ! r
            else if (trim(gridID) .eq. SH_PS16) then
               do r = 1, nr_out
                  do c = 1, nc_out
                     if (rlat_bin(c + (r-1)*nc_out) > 0) then
                        go(c + (r-1)*nc_out) = LVT_rc%udef
                     end if
                  end do ! c
               end do ! r
            end if
            ! If using Air Force polar stereographic, we must flip the grid so
            ! the origin is in the upper-left corner instead of lower-left
            do r = 1, nr_out
               do c = 1, nc_out
                  go2d(c,nr_out - r + 1) = go(c + (r-1)*nc_out)
               end do ! c
            end do ! r
            do r = 1, nr_out
               do c = 1, nc_out
                  go(c + (r-1)*nc_out) = go2d(c,r)
               end do ! c
            end do ! r
         end if
         call write_grib1(ftn, griddesco, nc_out, nr_out, go, param=128, &
              decimal_scale_factor=0, bits_per_value=7, &
              grid_definition=grid_definition)

         ! Handle icemask
         do r = 1, this%nr
            do c = 1, this%nc
               if (this%icemask(c,r) < 0) then
                  li(c + (r-1)*this%nc) = .false.
                  gi(c + (r-1)*this%nc) = LVT_rc%udef
               else
                  li(c + (r-1)*this%nc) = .true.
                  gi(c + (r-1)*this%nc) = this%icemask(c,r)
               end if
            end do ! c
         end do ! r
         if (griddesco(1) == 0) then
            call upscaleByMode((this%nc*this%nr), &
                 (nc_out*nr_out), LVT_rc%udef, n11, li, gi, lo, go)
         else if (griddesco(1) == 5) then
            call neighbor_interp(griddesco, li, gi, lo_neighbor, go_neighbor, &
                 (this%nc*this%nr), (nc_out*nr_out), &
                 rlat_neighbor, rlon_neighbor, &
                 n11_neighbor, &
                 LVT_rc%udef, iret)
            do r = 1, nr_out
               do c = 1, nc_out
                  go(c + (r-1)*nc_out) = go_neighbor(c + (r-1)*nc_out)
               end do
            end do
            ! Filter out points that are outside of the hemisphere.
            if (trim(gridID) .eq. NH_PS16) then
               do r = 1, nr_out
                  do c = 1, nc_out
                     if (rlat_bin(c + (r-1)*nc_out) < 0) then
                        go(c + (r-1)*nc_out) = LVT_rc%udef
                     end if
                  end do ! c
               end do ! r
            else if (trim(gridID) .eq. SH_PS16) then
               do r = 1, nr_out
                  do c = 1, nc_out
                     if (rlat_bin(c + (r-1)*nc_out) > 0) then
                        go(c + (r-1)*nc_out) = LVT_rc%udef
                     end if
                  end do ! c
               end do ! r
            end if
            ! If using Air Force polar stereographic, we must flip the grid so
            ! the origin is in the upper-left corner instead of lower-left
            do r = 1, nr_out
               do c = 1, nc_out
                  go2d(c,nr_out - r + 1) = go(c + (r-1)*nc_out)
               end do ! c
            end do ! r
            do r = 1, nr_out
               do c = 1, nc_out
                  go(c + (r-1)*nc_out) = go2d(c,r)
               end do ! c
            end do ! r
         end if
         call write_grib1(ftn, griddesco, nc_out, nr_out, go, param=91, &
              decimal_scale_factor=0, bits_per_value=1, &
              grid_definition=grid_definition)

         ! Handle iceage
         do r = 1, this%nr
            do c = 1, this%nc
               if (this%icemask(c,r) < 0) then
                  li(c + (r-1)*this%nc) = .false.
                  gi(c + (r-1)*this%nc) = LVT_rc%udef
               else
                  li(c + (r-1)*this%nc) = .true.
                  gi(c + (r-1)*this%nc) = this%iceage(c,r)
               end if
            end do ! c
         end do ! r
         if (griddesco(1) == 0) then
            call upscaleByMode((this%nc*this%nr), &
                 (nc_out*nr_out), LVT_rc%udef, n11, li, gi, lo, go)
         else if (griddesco(1) == 5) then
            call neighbor_interp(griddesco, li, gi, lo_neighbor, go_neighbor, &
                 (this%nc*this%nr), (nc_out*nr_out), &
                 rlat_neighbor, rlon_neighbor, &
                 n11_neighbor, &
                 LVT_rc%udef, iret)
            do r = 1, nr_out
               do c = 1, nc_out
                  go(c + (r-1)*nc_out) = go_neighbor(c + (r-1)*nc_out)
               end do
            end do
            ! Filter out points that are outside of the hemisphere.
            if (trim(gridID) .eq. NH_PS16) then
               do r = 1, nr_out
                  do c = 1, nc_out
                     if (rlat_bin(c + (r-1)*nc_out) < 0) then
                        go(c + (r-1)*nc_out) = LVT_rc%udef
                     end if
                  end do ! c
               end do ! r
            else if (trim(gridID) .eq. SH_PS16) then
               do r = 1, nr_out
                  do c = 1, nc_out
                     if (rlat_bin(c + (r-1)*nc_out) > 0) then
                        go(c + (r-1)*nc_out) = LVT_rc%udef
                     end if
                  end do ! c
               end do ! r
            end if
            ! If using Air Force polar stereographic, we must flip the grid so
            ! the origin is in the upper-left corner instead of lower-left
            do r = 1, nr_out
               do c = 1, nc_out
                  go2d(c,nr_out - r + 1) = go(c + (r-1)*nc_out)
               end do ! c
            end do ! r
            do r = 1, nr_out
               do c = 1, nc_out
                  go(c + (r-1)*nc_out) = go2d(c,r)
               end do ! c
            end do ! r
         end if
         call write_grib1(ftn, griddesco, nc_out, nr_out, go, param=129, &
           decimal_scale_factor=0, bits_per_value=9, &
           grid_definition=grid_definition)

         ! Close the GRIB1 file
         call grib_close_file(ftn, rc)
         if (rc .ne. GRIB_SUCCESS) then
            write(LVT_logunit,*)'[ERR] Error from grib_close_file'
            call grib_get_error_string(rc, msg, status2)
            write(LVT_logunit,*)'[ERR] ', trim(msg)
            write(LVT_logunit,*)'[ERR] LVT will stop'
            call LVT_endrun()
         end if

      end if ! if write_fullgrib_file


      ! EMK...Additional SNODEP file for 16th mesh only
      write_snodep_file = .false.
      if (trim(gridID) .eq. trim(NH_PS16) .and. &
           LVT_rc%output_nh_ps16_snodep) then
         write_snodep_file = .true.
      end if
      if (trim(gridID) .eq. trim(SH_PS16) .and. &
           LVT_rc%output_sh_ps16_snodep) then
         write_snodep_file = .true.
      end if
      if (write_snodep_file) then
         call build_filename_g1_snodep(gridID, LVT_rc%output_dir, &
              LVT_rc%yyyymmddhh, fname)

         ! Open the GRIB1 file
         call grib_open_file(ftn, fname, 'w', rc)
         if (rc .ne. GRIB_SUCCESS) then
            write(LVT_logunit,*)'[ERR] Error from grib_open_file'
            call grib_get_error_string(rc, msg, status2)
            write(LVT_logunit,*)'[ERR] ', trim(msg)
            write(LVT_logunit,*)'[ERR] LVT will stop'
            call LVT_endrun()
         else
            write(LVT_logunit,*)'[INFO] Writing to ', trim(fname)
         end if

         ! Interpolate snoanl
         do r = 1, this%nr
            do c = 1, this%nc
               if (this%snoanl(c,r) < 0) then
                  li(c + (r-1)*this%nc) = .false.
                  gi(c + (r-1)*this%nc) = LVT_rc%udef
               else
                  li(c + (r-1)*this%nc) = .true.
                  gi(c + (r-1)*this%nc) = this%snoanl(c,r)
               end if
            end do ! c
         end do ! r
         call bilinear_interp(griddesco, li, gi, lo_bin, go_bin, &
              (this%nc*this%nr), (nc_out*nr_out), &
              rlat_bin, rlon_bin, &
              w11_bin, w12_bin, w21_bin, w22_bin, &
              n11_bin, n12_bin, n21_bin, n22_bin, &
              LVT_rc%udef, iret)
         do r = 1, nr_out
            do c = 1, nc_out
               go(c + (r-1)*nc_out) = go_bin(c + (r-1)*nc_out)
            end do
         end do

         ! Filter out points that are outside of the hemisphere.
         if (trim(gridID) .eq. NH_PS16) then
            do r = 1, nr_out
               do c = 1, nc_out
                  if (rlat_bin(c + (r-1)*nc_out) < 0) then
                     go(c + (r-1)*nc_out) = LVT_rc%udef
                  end if
               end do ! c
            end do ! r
         else if (trim(gridID) .eq. SH_PS16) then
            do r = 1, nr_out
               do c = 1, nc_out
                  if (rlat_bin(c + (r-1)*nc_out) > 0) then
                     go(c + (r-1)*nc_out) = LVT_rc%udef
                  end if
               end do ! c
            end do ! r
         end if

         ! If using Air Force polar stereographic, we must flip the grid so
         ! the origin is in the upper-left corner instead of lower-left
         do r = 1, nr_out
            do c = 1, nc_out
               go2d(c,nr_out - r + 1) = go(c + (r-1)*nc_out)
            end do ! c
         end do ! r
         do r = 1, nr_out
            do c = 1, nc_out
               go(c + (r-1)*nc_out) = go2d(c,r)
            end do ! c
         end do ! r

         ! Write the message
         call write_grib1(ftn, griddesco, nc_out, nr_out, go, param=66, &
              decimal_scale_factor=2, bits_per_value=8, &
              grid_definition=grid_definition)

         ! Close the GRIB1 file
         call grib_close_file(ftn, rc)
         if (rc .ne. GRIB_SUCCESS) then
            write(LVT_logunit,*)'[ERR] Error from grib_close_file'
            call grib_get_error_string(rc, msg, status2)
            write(LVT_logunit,*)'[ERR] ', trim(msg)
            write(LVT_logunit,*)'[ERR] LVT will stop'
            call LVT_endrun()
         end if

      end if ! if write_snodep_file

      ! Clean up
      deallocate(n11)
      if (allocated(rlat_bin)) deallocate(rlat_bin)
      if (allocated(rlon_bin)) deallocate(rlon_bin)
      if (allocated(n11_bin)) deallocate(n11_bin)
      if (allocated(n12_bin)) deallocate(n12_bin)
      if (allocated(n21_bin)) deallocate(n21_bin)
      if (allocated(n22_bin)) deallocate(n22_bin)
      if (allocated(w11_bin)) deallocate(w11_bin)
      if (allocated(w12_bin)) deallocate(w12_bin)
      if (allocated(w21_bin)) deallocate(w21_bin)
      if (allocated(w22_bin)) deallocate(w22_bin)
      if (allocated(rlat_neighbor)) deallocate(rlat_neighbor)
      if (allocated(rlon_neighbor)) deallocate(rlon_neighbor)
      if (allocated(n11_neighbor)) deallocate(n11_neighbor)
      deallocate(li)
      deallocate(gi)
      deallocate(lo)
      deallocate(go)
      if (allocated(go2d)) deallocate(go2d)
      if (allocated(lo_bin)) deallocate(lo_bin)
      if (allocated(go_bin)) deallocate(go_bin)
      if (allocated(lo_neighbor)) deallocate(lo_neighbor)
      if (allocated(go_neighbor)) deallocate(go_neighbor)
   end subroutine interp_and_output_grib1

   ! Internal subroutine for checking gridID
   subroutine check_gridID(gridID)

      ! Imports
      use LVT_logMod, only: LVT_logunit, LVT_endrun

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: gridID

      if (trim(gridID) .ne. trim(GLOBAL_LL0P25) .and. &
           trim(gridID) .ne. trim(NH_PS16) .and. &
           trim(gridID) .ne. trim(SH_PS16)) then
         write(LVT_logunit,*) &
              '[ERR] Unrecognized output gridID!'
         write(LVT_logunit,*) &
              '[ERR] Expected one of: '
         write(LVT_logunit,*) trim(GLOBAL_LL0P25)
         write(LVT_logunit,*) trim(NH_PS16)
         write(LVT_logunit,*) trim(SH_PS16)
         write(LVT_logunit,*) &
              '[ERR] Found ',trim(gridID)
         write(LVT_logunit,*) &
              '[ERR] LVT will exit gracefully'
         call LVT_endrun()
      end if
   end subroutine check_gridID

   ! Build the grib2 filename
   subroutine build_filename_g2(output_dir, yyyymmddhh, filename)

      ! Imports
      use LVT_coreMod, only: LVT_rc

      ! Defaults
      implicit none

      ! Arguments
      character(len=255), intent(in) :: output_dir
      character(len=10), intent(in) :: yyyymmddhh
      character(len=255), intent(out) :: filename

      ! filename = trim(output_dir)  &
      !      // '/PS.557WW_SC.' &
      !      // trim(LVT_rc%security_class)//'_DI.' &
      !      // trim(LVT_rc%data_category)//'_GP.' &
      !      // 'LIS-SNOWICE_GR.C0P09DEG_AR.' &
      !      // trim(LVT_rc%area_of_data)//'_PA.' &
      !      //'USAFSI_DD.' &
      !      // yyyymmddhh(1:8)//'_DT.' &
      !      // yyyymmddhh(9:10)//'00_DF.GR2'
      filename = trim(output_dir)  &
           // '/PS.557WW_SC.' &
           // trim(LVT_rc%security_class)//'_DI.' &
           // trim(LVT_rc%data_category)//'_GP.' &
           // 'USAFSI_GR.C0P09DEG_AR.' &
           // trim(LVT_rc%area_of_data)//'_PA.' &
           // 'SNOW-ICE_DD.' &
           // yyyymmddhh(1:8)//'_DT.' &
           // yyyymmddhh(9:10)//'00_DF.GR2'

   end subroutine build_filename_g2

   ! Build the grib1 filename
   subroutine build_filename_g1(gridID, output_dir, yyyymmddhh, filename)

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: gridID
      character(len=255), intent(in) :: output_dir
      character(len=10), intent(in) :: yyyymmddhh
      character(len=255), intent(out) :: filename

      ! Local variables
      character(len=10) :: grid
      character(len=10) :: area
      character(len=8)  :: yyyymmdd
      character(len=2)  :: hh

      if (trim(gridID) .eq. trim(GLOBAL_LL0P25)) then
         grid = 'C0P25DEG'
         area = 'GLOBAL'
      else if (trim(gridID) .eq. trim(NH_PS16)) then
         grid = 'P24KM'
         area = 'N-HEM'
      else if (trim(gridID) .eq. trim(SH_PS16)) then
         grid = 'P24KM'
         area = 'S-HEM'
      end if
      yyyymmdd = yyyymmddhh(1:8)
      hh = yyyymmddhh(9:10)

      filename = trim(output_dir) // '/' // &
           'PS.AFWA_SC.U_DI.D_GP.SNODEP-U_GR.' // &
           trim(grid) // '_AR.' // &
           trim(area) // '_PA.SNODEP_DD.' // &
           trim(yyyymmdd) // '_DT.' // &
           trim(hh) // '00_DF.GR1'

   end subroutine build_filename_g1

   ! Build the grib1 filename just for snodep
   subroutine build_filename_g1_snodep(gridID, output_dir, yyyymmddhh, &
        filename)

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: gridID
      character(len=255), intent(in) :: output_dir
      character(len=10), intent(in) :: yyyymmddhh
      character(len=255), intent(out) :: filename

      ! Local variables
      character(len=10) :: area

      if (trim(gridID) .eq. trim(NH_PS16)) then
         area = 'NH'
      else if (trim(gridID) .eq. trim(SH_PS16)) then
         area = 'SH'
      end if

      filename = trim(output_dir) // '/' // &
           'SNODEP_16_' // trim(area) // '_' // &
           trim(yyyymmddhh) // '.GR1'

   end subroutine build_filename_g1_snodep

   ! Internal subroutine for setting griddesci
   ! FIXME:  Add support for non-lat/lon projections
   subroutine set_griddesci(this, griddesci)

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_USAFSIpost_t), intent(inout) :: this
      real, intent(inout) ::  griddesci(50)

      griddesci(:) = 0

      griddesci(1) = 0 ! lat/lon
      griddesci(2) = this%nc
      griddesci(3) = this%nr
      griddesci(4) = this%sw_corner_lat
      griddesci(5) = this%sw_corner_lon
      griddesci(6) = 128
      griddesci(7) = this%ne_corner_lat
      griddesci(8) = this%ne_corner_lon
      griddesci(9) = this%dx
      griddesci(10) = this%dy
      griddesci(11) = 64
      griddesci(20) = 64
      griddesci(30) = 0
      griddesci(32) = this%nc
      griddesci(33) = this%nr
      griddesci(34) = this%sw_corner_lat
      griddesci(35) = this%sw_corner_lon
      griddesci(36) = 128
      griddesci(37) = this%ne_corner_lat
      griddesci(38) = this%ne_corner_lon
      griddesci(39) = this%dx
      griddesci(40) = this%dy

   end subroutine set_griddesci

   ! Internal subroutine for setting griddesco for global lat/lon 0.25 deg grid
   subroutine set_griddesco_global_ll0p25(this, griddesco)

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_USAFSIpost_t), intent(inout) :: this
      real, intent(inout) ::  griddesco(50)

      griddesco(:) = 0
      griddesco(1) = 0 ! lat/lon
      griddesco(2) = 1440
      griddesco(3) =  720
      griddesco(4) =  -89.875000
      griddesco(5) = -179.875000
      griddesco(6) = 128
      griddesco(7) =   89.875000
      griddesco(8) =  179.875000
      griddesco(9) =  0.250000
      griddesco(10) = 0.250000
      griddesco(11) = 64
      griddesco(20) = 64
      griddesco(30) = 0
      griddesco(32) = 1440
      griddesco(33) =  720
      griddesco(34) =  -89.875000
      griddesco(35) = -179.875000
      griddesco(36) = 128
      griddesco(37) =   89.875000
      griddesco(38) =  179.875000
      griddesco(39) = 0.250000
      griddesco(40) = 0.250000

   end subroutine set_griddesco_global_ll0p25

   ! Internal subroutine for setting griddesco for Northern Hemisphere 16th
   ! Mesh polar stereographic grid, in GRIB1.
   subroutine set_griddesco_nh_ps16(this, griddesco)

      ! Imports
      use LVT_logMod, only: LVT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_USAFSIpost_t), intent(inout) :: this
      real, intent(inout) ::  griddesco(50)

      ! Local variables
      real :: xmesh, xpnmcaf, ypnmcaf, orient, xj, xi, alat, alon

      xmesh = 23.813 ! Per 557WW GRIB1 manual
      xpnmcaf = 513
      ypnmcaf = 513
      orient = 100.0

      ! We need to use the USAF code to calculate the lat/lon.  However,
      ! the Air Force grid specifies the origin in the upper-left corner,
      ! while the NCEP interpolation code specifies in the origin in the
      ! lower-left corner.  Since we must first interpolate, we will use
      ! the NCEP convention at this step
      call pstoll(1, 1, float(1), float(1024), 16, alat, alon)
      if (alon > 180.) then
         alon = alon - 360.
      end if

      griddesco(:) = 0
      griddesco(1) = 5
      griddesco(2) = 1024
      griddesco(3) = 1024
      griddesco(4) = alat
      griddesco(5) = alon
      griddesco(6) = 8
      griddesco(7) = orient
      griddesco(8) = xmesh
      griddesco(9) = xmesh
      griddesco(10) = 60.0
      griddesco(11) = orient
      griddesco(20) = 64

      ! Stash away the upper-left lat/lon for later use
      call pstoll(1, 1, float(1), float(1), 16, alat, alon)
      griddesco(30) = alat
      griddesco(31) = alon
   end subroutine set_griddesco_nh_ps16

   ! Internal subroutine for setting griddesco for Southern Hemisphere 16th
   ! Mesh polar stereographic grid, in GRIB1.
   subroutine set_griddesco_sh_ps16(this, griddesco)

      ! Imports
      use LVT_logMod, only: LVT_logunit

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_USAFSIpost_t), intent(inout) :: this
      real, intent(inout) ::  griddesco(50)

      ! Local variables
      real :: xmesh, xpnmcaf, ypnmcaf, orient, xj, xi, alat, alon

      xmesh = 23.813 ! Per 557WW GRIB1 manual
      xpnmcaf = 513
      ypnmcaf = 513
      orient = 100.

      ! We need to use the USAF code to calculate the lat/lon.  However,
      ! the Air Force grid specifies the origin in the upper-left corner,
      ! while the NCEP interpolation code  specifies in the origin in the
      ! lower-left corner.  Since we must first interpolate, we will use
      ! the NCEP convention at this step
      call pstoll(2, 1, float(1), float(1024), 16, alat, alon)

      griddesco(:) = 0
      griddesco(1) = 5
      griddesco(2) = 1024
      griddesco(3) = 1024
      griddesco(4) = alat
      griddesco(5) = alon
      griddesco(6) = 8
      griddesco(7) = orient
      griddesco(8) = xmesh
      griddesco(9) = xmesh
      griddesco(10) = -60.0
      griddesco(11) = orient
      griddesco(20) = 64

      ! Stash away the upper-left lat/lon for later use
      call pstoll(2, 1, float(1), float(1), 16, alat, alon)
      griddesco(30) = alat
      griddesco(31) = alon

   end subroutine set_griddesco_sh_ps16

   ! Internal subroutine for writing grib2 message
   subroutine write_grib2(ftn, griddesco, &
        nc_out, nr_out, go, dspln, cat, num, typegenproc, fcsttime)

      ! Imports
      use grib_api
      use LVT_coreMod, only: LVT_rc
      use LVT_gribWrapperMod, only: LVT_grib_set
      use LVT_logMod, only: LVT_logunit, LVT_endrun, LVT_verify

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: ftn
      real, intent(in) :: griddesco(50)
      integer, intent(in) :: nc_out
      integer, intent(in) :: nr_out
      real, intent(in) :: go(nc_out*nr_out)
      integer, intent(in) :: dspln
      integer, intent(in) :: cat
      integer, intent(in) :: num
      integer, intent(in) :: typegenproc
      integer, intent(in) :: fcsttime

      ! Local variables
      integer :: igrib, rc, status2
      character(len=8) :: yyyymmdd
      character(len=4) :: hhmm
      integer :: idate8, idate4
      character(len=255) :: msg

#if (defined USE_ECCODES)
      call grib_new_from_samples(igrib, "GRIB2", rc)
      if (rc .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_new_from_samples'
         call grib_get_error_string(rc, msg, status2)
         write(LVT_logunit,*)'[ERR] ',trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if
#else
      call grib_new_from_template(igrib, "GRIB2", rc)
      if (rc .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_new_from_template'
         call grib_get_error_string(rc, msg, status2)
         write(LVT_logunit,*)'[ERR] ',trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if
#endif

      ! Section 0: Indicator
      ! Octet 7
      call LVT_grib_set(igrib, 'discipline', dspln)

      ! Section 1: Identification
      ! Octets 6-7
      call LVT_grib_set(igrib, 'centre', LVT_rc%grib_center_id)
      ! Octets 8-9
      call LVT_grib_set(igrib, 'subCentre', LVT_rc%grib_subcenter_id)
      ! Octet 10
      call LVT_grib_set(igrib, 'tablesVersion', LVT_rc%grib_table)
      ! Octet 11
      call LVT_grib_set(igrib, 'localTablesVersion', 1)
      ! Octet 12
      call LVT_grib_set(igrib, 'significanceOfReferenceTime', 0)
      ! Octet 13-16
      yyyymmdd = LVT_rc%yyyymmddhh(1:8)
      read(yyyymmdd,'(i8)') idate8
      call LVT_grib_set(igrib, 'dataDate', idate8)
      ! Octet 17-19
      hhmm = LVT_rc%yyyymmddhh(9:10) // '00'
      read(hhmm,'(i4)') idate4
      call LVT_grib_set(igrib, 'dataTime', idate4)
      ! Octet 20
      call LVT_grib_set(igrib, 'productionStatusOfProcessedData', 0)
      ! Octet 21
      call LVT_grib_set(igrib, 'typeOfProcessedData', 0)

      ! Section 2:  Local Use Section (Optional) -- none for now

      ! Section 3: Grid
      if (griddesco(1) == 0) then
         call LVT_grib_set(igrib, 'gridDefinitionTemplateNumber', 0)
      else if (griddesco(1) == 5) then
         call LVT_grib_set(igrib, 'gridDefinitionTemplateNumber', 20)
      end if
      ! Octet 15
      if (griddesco(1) == 0) then
         call LVT_grib_set(igrib, 'shapeOfTheEarth', 0)
      else if (griddesco(1) == 5) then
         call LVT_grib_set(igrib, 'shapeOfTheEarth', 1)
         call LVT_grib_set(igrib, 'scaleFactorOfRadiusOfSphericalEarth', 10.)
         call LVT_grib_set(igrib, 'scaledValueOfRadiusOfSphericalEarth', &
              63712213)
      end if

      ! Set dimensions
      if (griddesco(1) == 0) then
         ! Latitude/longitude
         call LVT_grib_set(igrib, 'Ni', nc_out)
         call LVT_grib_set(igrib, 'Nj', nr_out)
         call LVT_grib_set(igrib, 'basicAngleOfTheInitialProductionDomain', &
              0)
         call LVT_grib_set(igrib, 'latitudeOfFirstGridPointInDegrees', &
              griddesco(4))
         call LVT_grib_set(igrib, 'longitudeOfFirstGridPointInDegrees', &
              griddesco(5))
         call LVT_grib_set(igrib, 'latitudeOfLastGridPointInDegrees', &
              griddesco(7))
         call LVT_grib_set(igrib, 'longitudeOfLastGridPointInDegrees', &
              griddesco(8))
         call LVT_grib_set(igrib, 'iDirectionIncrementInDegrees', &
              griddesco(9))
         call LVT_grib_set(igrib, 'jDirectionIncrementInDegrees', &
              griddesco(10))
         call LVT_grib_set(igrib, 'iScansNegatively', 0)
         call LVT_grib_set(igrib, 'jScansPositively', 1)
         call LVT_grib_set(igrib, 'jPointsAreConsecutive', 0)

      else if (griddesco(1) == 5) then
         ! Polar stereographic
         call LVT_grib_set(igrib, 'Nx', nc_out)
         call LVT_grib_set(igrib, 'Ny', nr_out)
         call LVT_grib_set(igrib, 'latitudeOfFirstGridPointInDegrees', &
              griddesco(30))
         ! NOTE:  ECCODES will not accept a value of -125, so we use the
         ! equivalent of 235 in the southern hemisphere
         call LVT_grib_set(igrib, 'longitudeOfFirstGridPointInDegrees', &
              griddesco(31))
         call LVT_grib_set(igrib, 'resolutionAndComponentFlags', 0)
         call LVT_grib_set(igrib, 'LaD', 1e6*griddesco(10))
         call LVT_grib_set(igrib, 'orientationOfTheGrid', &
              1e6*griddesco(7))
         call LVT_grib_set(igrib, 'DxInMetres', &
              abs(1000*griddesco(8)))
         call LVT_grib_set(igrib, 'DyInMetres', &
              abs(1000*griddesco(9)))
         if (griddesco(4) < 0) then
            call LVT_grib_set(igrib, 'projectionCentreFlag', 0)
         else
            call LVT_grib_set(igrib, 'projectionCentreFlag', 1)
         end if
         call LVT_grib_set(igrib, 'scanningMode', 0)
      end if

      ! Section 4: Product Definition Section
      ! Octets 8-9
      call LVT_grib_set(igrib, 'productDefinitionTemplateNumber', 0)
      ! Octet 10
      call LVT_grib_set(igrib, 'parameterCategory', cat)
      ! Octet 11
      call LVT_grib_set(igrib, 'parameterNumber', num)
      ! Octet 12
      call LVT_grib_set(igrib, 'typeOfGeneratingProcess', typegenproc)
      ! Octet 13
      call LVT_grib_set(igrib, 'backgroundGeneratingProcessIdentifier', &
           LVT_rc%grib_process_id)
      ! Octet 14
      call LVT_grib_set(igrib, 'generatingProcessIdentifier', &
           LVT_rc%grib_process_id)
      ! Octet 15-17 is skipped
      ! Octet 18...Use hours
      call LVT_grib_set(igrib, 'indicatorOfUnitOfTimeRange', 1)
      ! Octets 19-22
      call LVT_grib_set(igrib, 'forecastTime', fcsttime)
      ! Octets 23-34...Ground or water surface
      call LVT_grib_set(igrib, 'typeOfFirstFixedSurface', 1)

      ! Section 5: Data Representation
      call LVT_grib_set(igrib, 'packingType', LVT_rc%grib_packing_type)
      call LVT_grib_set(igrib, 'missingValue', LVT_rc%udef)

      ! Section 6: Bit-Map
      call LVT_grib_set(igrib, 'bitmapPresent', 1)

      ! Section 7
      call LVT_grib_set(igrib, 'values', go)

      call grib_write(igrib, ftn, rc)
      if (rc .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_write'
         call grib_get_error_string(rc, msg, status2)
         write(LVT_logunit,*)'[ERR] ', trim(msg)
         write(LVT_Logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if

      call grib_release(igrib, rc)
      if (rc .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_release'
         call grib_get_error_string(rc, msg, status2)
         write(LVT_logunit,*)'[ERR] ', trim(msg)
         write(LVT_Logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if

   end subroutine write_grib2

   ! Internal subroutine for writing grib1 message
   subroutine write_grib1(ftn, griddesco, &
        nc_out, nr_out, go, param, decimal_scale_factor, &
        bits_per_value, grid_definition)

      ! Imports
      use grib_api
      use LVT_coreMod, only: LVT_rc
      use LVT_gribWrapperMod, only: LVT_grib_set
      use LVT_logMod, only: LVT_logunit, LVT_endrun, LVT_verify

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: ftn
      real, intent(in) :: griddesco(50)
      integer, intent(in) :: nc_out
      integer, intent(in) :: nr_out
      real, intent(in) :: go(nc_out*nr_out)
      integer, intent(in) :: param
      integer, intent(in) :: decimal_scale_factor
      integer, intent(in) :: bits_per_value
      integer, intent(in) :: grid_definition

      ! Local variables
      integer :: igrib, rc, status2
      character(len=255) :: msg
      character(len=4) :: cyyyy
      character(len=2) :: cmm, cdd, chh
      integer :: iyyyy, imm, idd, ihh, iyear, iyoc, ic

#if (defined USE_ECCODES)
      call grib_new_from_samples(igrib, "GRIB1", rc)
      if (rc .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_new_from_samples'
         call grib_get_error_string(rc, msg, status2)
         write(LVT_logunit,*)'[ERR] ',trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if
#else
      call grib_new_from_template(igrib, "GRIB1", rc)
      if (rc .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_new_from_template'
         call grib_get_error_string(rc, msg, status2)
         write(LVT_logunit,*)'[ERR] ',trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if
#endif

      ! Section 1: Product Definition Section
      ! Octet 4
      call LVT_grib_set(igrib, 'table2Version', 128)
      ! Octet 5
      call LVT_grib_set(igrib, 'centre', 57)
      ! Octet 6
      call LVT_grib_set(igrib, 'generatingProcessIdentifier', 35)
      ! Octet 7
      call LVT_grib_set(igrib, 'gridDefinition', grid_definition)
      ! Octet 9
      call LVT_grib_set(igrib, 'indicatorOfParameter', param)
      ! Octet 10
      call LVT_grib_set(igrib, 'indicatorOfTypeOfLevel', 1)
      ! Octet 11-12
      call LVT_grib_set(igrib, 'level', 0)
      ! Octet 13
      iyear = 1000 * (ichar(LVT_rc%yyyymmddhh(1:1)) - 48) +             &
               100 * (ichar(LVT_rc%yyyymmddhh(2:2)) - 48) +             &
                10 * (ichar(LVT_rc%yyyymmddhh(3:3)) - 48) +             &
                     (ichar(LVT_rc%yyyymmddhh(4:4)) - 48)
      iyoc = mod(iyear,1000)
      if (iyoc .eq. 0) iyoc = 100
      call LVT_grib_set(igrib, 'yearOfCentury', iyoc)

      ! Octet 14
      cmm = LVT_rc%yyyymmddhh(5:6)
      read(cmm, '(i2)') imm
      call LVT_grib_set(igrib, 'month', imm)
      ! Octet 15
      cdd = LVT_rc%yyyymmddhh(7:8)
      read(cdd, '(i2)') idd
      call LVT_grib_set(igrib, 'day', idd)
      ! Octet 16
      chh = LVT_rc%yyyymmddhh(9:10)
      read(chh, '(i2)') ihh
      call LVT_grib_set(igrib, 'hour', ihh)
      ! Octet 17
      call LVT_grib_set(igrib, 'minute', 0)
      ! Octet 18
      call LVT_grib_set(igrib, 'unitOfTimeRange', 1)
      ! Octet 19-20
      call LVT_grib_set(igrib, 'P1', 0)
      call LVT_grib_set(igrib, 'P2', 0)
      ! Octet 21
      call LVT_grib_set(igrib, 'timeRangeIndicator', 0)
      ! Octet 22-23
      call LVT_grib_set(igrib, 'numberIncludedInAverage', 0)
      ! Octet 24
      call LVT_grib_set(igrib, 'numberMissingFromAveragesOrAccumulations', 0)
      ! Octet 25
      if (iyoc == 100) then
         ic = iyear / 100
      else
         ic = (iyear / 100) + 1
      end if
      call LVT_grib_set(igrib, 'centuryOfReferenceTimeOfData', ic)
      ! Octet 26
      call LVT_grib_set(igrib, 'subCentre', 2)
      ! Octet 27-28
      call LVT_grib_set(igrib, 'decimalScaleFactor', decimal_scale_factor)

      ! Section 2: Grid Description Section
      ! Octet 4
      call LVT_grib_set(igrib, 'numberOfVerticalCoordinateValues', 0)
      ! Octet 5
      call LVT_grib_set(igrib, 'pvlLocation', 255)
      ! Octet 6
      if (griddesco(1) == 0) then
         call LVT_grib_set(igrib, 'dataRepresentationType', 0)
      else if (griddesco(1) == 5) then
         call LVT_grib_set(igrib, 'dataRepresentationType', 5)
      end if
      ! Octets 7-32
      if (griddesco(1) == 0) then
         call LVT_grib_set(igrib, 'Ni', nc_out)
         call LVT_grib_set(igrib, 'Nj', nr_out)
         call LVT_grib_set(igrib, 'latitudeOfFirstGridPointInDegrees', &
              griddesco(4))
         call LVT_grib_set(igrib, 'longitudeOfFirstGridPointInDegrees', &
              griddesco(5))
         call LVT_grib_set(igrib, 'latitudeOfLastGridPointInDegrees', &
              griddesco(7))
         call LVT_grib_set(igrib, 'longitudeOfLastGridPointInDegrees', &
              griddesco(8))
         call LVT_grib_set(igrib, 'iDirectionIncrementInDegrees', &
              griddesco(9))
         call LVT_grib_set(igrib, 'jDirectionIncrementInDegrees', &
              griddesco(10))
         call LVT_grib_set(igrib, 'iScansNegatively', 0)
         call LVT_grib_set(igrib, 'jScansPositively', 1)
         call LVT_grib_set(igrib, 'jPointsAreConsecutive', 0)

      else if (griddesco(1) == 5) then
         ! Polar stereographic
         call LVT_grib_set(igrib, 'Nx', nc_out)
         call LVT_grib_set(igrib, 'Ny', nr_out)
         call LVT_grib_set(igrib, 'latitudeOfFirstGridPointInDegrees', &
              griddesco(30))
         call LVT_grib_set(igrib, 'longitudeOfFirstGridPointInDegrees', &
              griddesco(31))
         call LVT_grib_set(igrib, 'resolutionAndComponentFlags', 136)
         call LVT_grib_set(igrib, 'orientationOfTheGrid', &
              1e3*griddesco(7))
         call LVT_grib_set(igrib, 'DxInMetres', &
              abs(1000*griddesco(8)))
         call LVT_grib_set(igrib, 'DyInMetres', &
              abs(1000*griddesco(9)))
         if (griddesco(4) < 0) then
            call LVT_grib_set(igrib, 'projectionCentreFlag', 0)
         else
            call LVT_grib_set(igrib, 'projectionCentreFlag', 128)
         end if
         call LVT_grib_set(igrib, 'scanningMode', 0)
      end if

      ! Section 3: Bit-map section
      ! This is handled implicitly by ECCODES based on the following
      call LVT_grib_set(igrib, 'missingValue', LVT_rc%udef)
      call LVT_grib_set(igrib, 'bitmapPresent', 1)
      call LVT_grib_set(igrib, 'bitsPerValue', bits_per_value)
      call LVT_grib_set(igrib, 'values', go)

      ! Section 4:  Binary data section -- handled implicitly

      ! Write the message
      call grib_write(igrib, ftn, rc)
      if (rc .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_write'
         call grib_get_error_string(rc, msg, status2)
         write(LVT_logunit,*)'[ERR] ', trim(msg)
         write(LVT_Logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if

      call grib_release(igrib, rc)
      if (rc .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_release'
         call grib_get_error_string(rc, msg, status2)
         write(LVT_logunit,*)'[ERR] ', trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if

   end subroutine write_grib1

   ! Internal subroutine for writing global lat/lon output in netCDF.
   ! For testing purposes.
   subroutine write_netcdf_latlon(griddesco, nc_out, nr_out, go)

      ! Imports
      use LVT_coreMod, only: LVT_rc
      use LVT_logMod, only: LVT_logunit, LVT_verify, LVT_endrun
      use netcdf

      ! Defaults
      implicit none

      ! Arguments
      real, intent(in) :: griddesco(50)
      integer, intent(in) :: nc_out
      integer, intent(in) :: nr_out
      real, intent(in) :: go(nc_out*nr_out)

      ! Local variables
      character(len=255) :: outfilename
      integer :: shuffle, deflate, deflate_level
      integer :: iret, ncid
      integer :: dim_ids(3)
      real :: dlon, dlat, swlat, swlon, nelat, nelon
      integer :: snoanl_varid, snoanl2_varid
      integer :: lon_varid, lat_varid, time_varid
      character*4 :: cyyyy
      character*2 :: cmm,cdd,chh
      character*120 :: time_units
      real, allocatable :: lats(:), lons(:)
      integer :: i,j,c,r
      real, allocatable :: snoanl(:,:)

      outfilename = "foo.nc"
      write(LVT_logunit,*)'[INFO] Creating netCDF file ', trim(outfilename)

      ! Copy the netcdf compression settings
      shuffle = NETCDF_shuffle
      deflate = NETCDF_deflate
      deflate_level = NETCDF_deflate_level

      ! Create the output file
      iret=nf90_create(path=trim(outfilename), &
           cmode=nf90_netcdf4, ncid=ncid)
      call LVT_verify(iret, &
           '[ERR] nf90_create failed')

      ! Write out dimensions headers
      call LVT_verify(nf90_def_dim(ncid,'time',1,dim_ids(3)), &
           '[ERR] nf90_def_dim failed')
      call LVT_verify(nf90_def_dim(ncid,'lat',nr_out,dim_ids(2)), &
           '[ERR] nf90_def_dim failed')
      call LVT_verify(nf90_def_dim(ncid,'lon',nc_out,dim_ids(1)), &
           '[ERR] nf90_def_dim failed')

      ! Map projection
      select case ("latlon")
      case ("latlon")
         dlon  = gridDesco(9)
         dlat  = gridDesco(10)
         swlat = gridDesco(4)
         swlon = gridDesco(5)
         nelat = gridDesco(7)
         nelon = gridDesco(8)

         call LVT_verify(nf90_put_att(ncid,nf90_global,&
              "MAP_PROJECTION", "EQUIDISTANT CYLINDRICAL"), &
              '[ERR] nf90_put_att failed')
         call LVT_verify(nf90_put_att(ncid,nf90_global,&
              "SOUTH_WEST_CORNER_LAT", swlat), &
              '[ERR] nf90_put_att failed')
         call LVT_verify(nf90_put_att(ncid,nf90_global,&
              "SOUTH_WEST_CORNER_LON", swlon), &
              '[ERR] nf90_put_att failed')
         call LVT_verify(nf90_put_att(ncid,nf90_global, &
              "NORTH_EAST_CORNER_LAT", nelat), &
              '[ERR] nf90_put_att failed')
         call LVT_verify(nf90_put_att(ncid,nf90_global, &
              "NORTH_EAST_CORNER_LON", nelon), &
              '[ERR] nf90_put_att failed')
         call LVT_verify(nf90_put_att(ncid,nf90_global, &
              "DX",dlon),&
              '[ERR] nf90_put_att failed')
         call LVT_verify(nf90_put_att(ncid,nf90_global, &
              "DY",dlat), &
              '[ERR] nf90_put_att failed')

      case default
         write(LVT_logunit,*) &
              '[ERR] Only latlon map projection supported!'
         call LVT_endrun()
      end select

      ! Include the water points
      call LVT_verify(nf90_put_att(ncid,nf90_global, &
           "INC_WATER_PTS","true"), &
           '[ERR] nf90_put_att failed')

      ! Construct the longitudes
      ! FIXME:  Add support for other map projections
      call LVT_verify(nf90_def_var(ncid,"lon",nf90_float,dim_ids(1), &
           lon_varid),'[ERR] nf90_def_var failed')
      call LVT_verify(nf90_def_var_deflate(ncid,&
           lon_varid, shuffle, deflate, deflate_level),&
           '[ERR] nf90_def_var_deflate')
      call LVT_verify(nf90_put_att(ncid,lon_varid, &
           "units","degrees_east"), &
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,lon_varid, &
           "long_name","longitude"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,lon_varid, &
           "standard_name","longitude"),&
           '[ERR] nf90_put_att failed')

      ! Construct the latitudes
      ! FIXME:  Add support for other map projections
      call LVT_verify(nf90_def_var(ncid,"lat",nf90_float,dim_ids(2), &
           lat_varid),'[ERR] nf90_def_var failed')
      call LVT_verify(nf90_def_var_deflate(ncid,&
           lat_varid, shuffle, deflate, deflate_level),&
           '[ERR] nf90_def_var_deflate')
      call LVT_verify(nf90_put_att(ncid,lat_varid, &
           "units","degrees_north"), &
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,lat_varid, &
           "long_name","latitude"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,lat_varid, &
           "standard_name","latitude"),&
           '[ERR] nf90_put_att failed')

      ! Define the time array.  The valid time will be written as an
      ! attribute.
      call LVT_verify(nf90_def_var(ncid,'time',nf90_double,&
           dimids=dim_ids(3), varid=time_varid), &
           '[ERR] nf90_def_var failed')
      cyyyy = '2018'
      cmm = '12'
      cdd = '08'
      chh = '00'
      write(time_units,'(A)') &
           "seconds since "//trim(cyyyy)//"-"//trim(cmm)//"-"//trim(cdd)//&
           " "//trim(chh)//":00:00"
      call LVT_verify(nf90_put_att(ncid,time_varid, &
           "units",trim(time_units)),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,time_varid, &
           "long_name","time"),&
           '[ERR] nf90_put_att failed')

      ! Define the snow depth analysis
      call LVT_verify(nf90_def_var(ncid,"snoanl",nf90_float, &
           dimids=dim_ids, &
           varid=snoanl_varid),'[ERR] nf90_def_var failed')
      call LVT_verify(nf90_def_var_deflate(ncid,&
           snoanl_varid, shuffle, deflate, deflate_level), &
           '[ERR] nf90_def_var_deflate failed')
      call LVT_verify(nf90_put_att(ncid,snoanl_varid, &
           "units","m"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,snoanl_varid, &
           "long_name","depth of surface snow over land"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,snoanl_varid, &
           '_FillValue',LVT_rc%udef), &
           '[ERR] nf90_put_att failed for SNOANL')

      ! Define the snow depth analysis
      call LVT_verify(nf90_def_var(ncid,"snoanl2",nf90_float, &
           dimids=dim_ids, &
           varid=snoanl2_varid),'[ERR] nf90_def_var failed')
      call LVT_verify(nf90_def_var_deflate(ncid,&
           snoanl2_varid, shuffle, deflate, deflate_level), &
           '[ERR] nf90_def_var_deflate failed')
      call LVT_verify(nf90_put_att(ncid,snoanl2_varid, &
           "units","m"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,snoanl2_varid, &
           "long_name","depth of surface snow over land"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,snoanl2_varid, &
           '_FillValue', LVT_rc%udef), &
           '[ERR] nf90_put_att failed for SNOANL2')

      ! Miscellaneous header information
      call LVT_verify(nf90_put_att(ncid,nf90_global,"Conventions", &
           "CF-1.7"), &
           '[ERR] nf90_put_att failed')

      ! We are ready to write the actual data.  This requires taking NETCDF
      ! out of define mode.
      call LVT_verify(nf90_enddef(ncid), &
           '[ERR] ncf90_enddef failed')

      ! Write the latitude data
      allocate(lats(nr_out))
      do j = 1, nr_out
         lats(j) = swlat + (j-1)*(dlat)
      end do
      call LVT_verify(nf90_put_var(ncid,lat_varid, &
           lats(:),(/1/),(/nr_out/)), &
           '[ERR] nf90_put_var failed for lats')
      deallocate(lats)

      ! Write the longitude data
      allocate(lons(nc_out))
      do i = 1, nc_out
         lons(i) = swlon + (i-1)*(dlon)
      end do
      call LVT_verify(nf90_put_var(ncid,lon_varid,lons(:),&
           (/1/),(/nc_out/)), &
           '[ERR] nf90_put_var failed for lon')
      deallocate(lons)

      ! Write the time data
      call LVT_verify(nf90_put_var(ncid,time_varid,0.0), &
           '[ERR] nf90_put_var failed for time')

      ! Write the USAFSI fields
      allocate(snoanl(nc_out,nr_out))
      do r = 1, nr_out
         do c = 1, nc_out
            snoanl(c,r) = go(c + (r-1)*nc_out)
         end do
      end do
      call LVT_verify(nf90_put_var(ncid,snoanl_varid,&
           snoanl(:,:), &
           (/1,1/),(/nc_out,nr_out/)), &
           '[ERR] nf90_put_var failed for snoanl')
      call LVT_verify(nf90_put_var(ncid,snoanl2_varid,&
           snoanl(:,:), &
           (/1,1/),(/nc_out,nr_out/)), &
           '[ERR] nf90_put_var failed for snoanl2')
      deallocate(snoanl)

      call LVT_verify(nf90_close(ncid), &
           '[ERR] nf90_close failed!')

   end subroutine write_netcdf_latlon

   ! Internal subroutine for writing polar stereographic output in netCDF.
   ! For testing purposes.
   subroutine write_netcdf_ps(griddesco, nc_out, nr_out, go)

      ! Imports
      use LVT_coreMod, only: LVT_rc
      use LVT_logMod, only: LVT_logunit, LVT_verify, LVT_endrun
      use netcdf

      ! Defaults
      implicit none

      ! Arguments
      real, intent(in) :: griddesco(50)
      integer, intent(in) :: nc_out
      integer, intent(in) :: nr_out
      real, intent(in) :: go(nc_out*nr_out)

      ! Local variables
      character(len=255) :: outfilename
      integer :: shuffle, deflate, deflate_level
      integer :: iret, ncid
      integer :: dim_ids(2)
      integer :: snoanl_varid
      integer :: lon_varid, lat_varid, xc_varid, yc_varid
      character*4 :: cyyyy
      character*2 :: cmm,cdd,chh
      character*120 :: time_units
      real, allocatable :: lats(:,:), lons(:,:)
      real, allocatable :: xc(:), yc(:)
      integer :: i,j,c,r
      real, allocatable :: snoanl(:,:)
      real :: alat, alon

      outfilename = "foo.nc"
      write(LVT_logunit,*)'[INFO] Creating netCDF file ', trim(outfilename)

      ! Copy the netcdf compression settings
      shuffle = NETCDF_shuffle
      deflate = NETCDF_deflate
      deflate_level = NETCDF_deflate_level

      ! Create the output file
      iret=nf90_create(path=trim(outfilename), &
           cmode=nf90_netcdf4, ncid=ncid)
      call LVT_verify(iret, &
           '[ERR] nf90_create failed')

      ! Write out dimensions headers
      call LVT_verify(nf90_def_dim(ncid,'yc',nr_out,dim_ids(2)), &
           '[ERR] nf90_def_dim failed')
      call LVT_verify(nf90_def_dim(ncid,'xc',nc_out,dim_ids(1)), &
           '[ERR] nf90_def_dim failed')

      call LVT_verify(nf90_def_var(ncid, "xc", nf90_float, dim_ids(1), &
           xc_varid),'[ERR] nf90_def_var failed')
      call LVT_verify(nf90_def_var_deflate(ncid, &
           xc_varid, shuffle, deflate, deflate_level), &
           '[ERR] nf90_def_var_deflate')
       call LVT_verify(nf90_put_att(ncid,xc_varid, &
           "axis","X"), &
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,xc_varid, &
           "long_name","x-coordinate in Cartesian system"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,xc_varid, &
           "units","m"),&
           '[ERR] nf90_put_att failed')

      call LVT_verify(nf90_def_var(ncid, "yc", nf90_float, dim_ids(2), &
           yc_varid),'[ERR] nf90_def_var failed')
      call LVT_verify(nf90_def_var_deflate(ncid, &
           yc_varid, shuffle, deflate, deflate_level), &
           '[ERR] nf90_def_var_deflate')
       call LVT_verify(nf90_put_att(ncid,yc_varid, &
           "axis","Y"), &
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,yc_varid, &
           "long_name","y-coordinate in Cartesian system"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,yc_varid, &
           "units","m"),&
           '[ERR] nf90_put_att failed')

      call LVT_verify(nf90_def_var(ncid, "lon", nf90_float, dim_ids, &
           lon_varid),'[ERR] nf90_def_var failed')
      call LVT_verify(nf90_def_var_deflate(ncid, &
           lon_varid, shuffle, deflate, deflate_level), &
           '[ERR] nf90_def_var_deflate')
      call LVT_verify(nf90_put_att(ncid, lon_varid, &
           "long_name","longitude"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid, lon_varid, &
           "units","degrees_east"),&
           '[ERR] nf90_put_att failed')

      call LVT_verify(nf90_def_var(ncid, "lat", nf90_float, dim_ids, &
           lat_varid),'[ERR] nf90_def_var failed')
      call LVT_verify(nf90_def_var_deflate(ncid, &
           lat_varid, shuffle, deflate, deflate_level), &
           '[ERR] nf90_def_var_deflate')
      call LVT_verify(nf90_put_att(ncid, lat_varid, &
           "long_name","latitude"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid, lat_varid, &
           "units","degrees_north"),&
           '[ERR] nf90_put_att failed')

      call LVT_verify(nf90_def_var(ncid,"snoanl",nf90_float, &
           dimids=dim_ids, &
           varid=snoanl_varid),'[ERR] nf90_def_var failed')
      call LVT_verify(nf90_def_var_deflate(ncid,&
           snoanl_varid, shuffle, deflate, deflate_level), &
           '[ERR] nf90_def_var_deflate failed')
      call LVT_verify(nf90_put_att(ncid,snoanl_varid, &
           "units","m"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,snoanl_varid, &
           "long_name","depth of surface snow over land"),&
           '[ERR] nf90_put_att failed')
      call LVT_verify(nf90_put_att(ncid,snoanl_varid, &
           '_FillValue', LVT_rc%udef), &
           '[ERR] nf90_put_att failed for SNOANL')

      ! Miscellaneous header information
      call LVT_verify(nf90_put_att(ncid,nf90_global,"Conventions", &
           "CF-1.7"), &
           '[ERR] nf90_put_att failed')

      ! We are ready to write the actual data.  This requires taking NETCDF
      ! out of define mode.
      call LVT_verify(nf90_enddef(ncid), &
           '[ERR] ncf90_enddef failed')

      allocate(xc(nc_out))
      do c = 1, nc_out
         xc(c) = c
      end do
      allocate(yc(nr_out))
      do r = 1, nr_out
         yc(r) = r
      end do
      allocate(lats(nc_out, nr_out))
      allocate(lons(nc_out, nr_out))
      do r = 1, nr_out
         do c = 1, nc_out
            call pstoll(2, 1, float(c), float(r), 16, alat, alon)
            lats(c,r) = alat
            lons(c,r) = alon
         end do
      end do
      call LVT_verify(nf90_put_var(ncid,lat_varid, &
           lats,(/1,1/),(/nc_out,nr_out/)), &
           '[ERR] nf90_put_var failed for lats')
      call LVT_verify(nf90_put_var(ncid,lon_varid, &
           lons,(/1,1/),(/nc_out,nr_out/)), &
           '[ERR] nf90_put_var failed for lons')
      call LVT_verify(nf90_put_var(ncid,xc_varid, &
           xc,(/1,1/),(/nc_out/)), &
           '[ERR] nf90_put_var failed for xc')
      call LVT_verify(nf90_put_var(ncid,yc_varid, &
           yc,(/1/),(/nr_out/)), &
           '[ERR] nf90_put_var failed for yc')
      deallocate(lats)
      deallocate(lons)
      deallocate(xc)
      deallocate(yc)
      allocate(snoanl(nc_out,nr_out))
      do r = 1, nr_out
         do c = 1, nc_out
            snoanl(c,r) = go(c + (r-1)*nc_out)
         end do
      end do
      call LVT_verify(nf90_put_var(ncid,snoanl_varid,&
           snoanl(:,:), &
           (/1,1/),(/nc_out,nr_out/)), &
           '[ERR] nf90_put_var failed for snoanl')
      deallocate(snoanl)

      call LVT_verify(nf90_close(ncid), &
           '[ERR] nf90_close failed!')

   end subroutine write_netcdf_ps

   ! Internal subroutine to set up weights for bilinear interpolation.
   ! Based on bilinear_interp_input, but handles Air Force 16th
   ! mesh polar stereographic grids.
   subroutine bilinear_interp_input_usaf(gridDesci, gridDesco, npts, &
        rlat, rlon, n11, n12, n21, n22, w11, w12, w21, w22, afwa_grid)

      ! Defaults
      implicit none

      ! Arguments
      real, intent(in)    :: gridDesci(50)
      real, intent(in)    :: gridDesco(50)
      integer, intent(in) :: npts
      real, intent(inout) :: rlat(npts)
      real, intent(inout) :: rlon(npts)
      integer, intent(inout) :: n11(npts)
      integer, intent(inout) :: n12(npts)
      integer, intent(inout) :: n21(npts)
      integer, intent(inout) :: n22(npts)
      real, intent(inout) :: w11(npts)
      real, intent(inout) :: w12(npts)
      real, intent(inout) :: w21(npts)
      real, intent(inout) :: w22(npts)
      character(len=*), intent(in) :: afwa_grid

      ! Local variables
      integer :: mo
      integer :: nv
      real :: xpts(npts), ypts(npts)
      integer :: n
      real :: xi, xf, yi, yf
      integer :: i1, i2, j1, j2
      integer, external :: get_fieldpos
      real, parameter     :: fill = -9999.0

      mo = npts
      if (trim(afwa_grid) .ne. NH_PS16 .and. &
           trim(afwa_grid) .ne. SH_PS16) then
         if (gridDesco(1).ge.0) then
            call compute_earth_coord(gridDesco, mo, fill, &
                 xpts, ypts, rlon, rlat,nv)
         endif
      end if
      call compute_grid_coord(gridDesci, mo, fill, xpts, ypts, rlon, rlat, nv)
      do n=1,mo
         xi=xpts(n)
         yi=ypts(n)
         if(xi.ne.fill.and.yi.ne.fill) then
            i1=xi
            i2=i1+1
            j1=yi
            j2=j1+1
            xf=xi-i1
            yf=yi-j1
            n11(n)=get_fieldpos(i1, j1, gridDesci)
            n21(n)=get_fieldpos(i2, j1, gridDesci)
            n12(n)=get_fieldpos(i1, j2, gridDesci)
            n22(n)=get_fieldpos(i2, j2, gridDesci)
            if(min(n11(n),n21(n),n12(n),n22(n)).gt.0) then
               w11(n)=(1-xf)*(1-yf)
               w21(n)=xf*(1-yf)
               w12(n)=(1-xf)*yf
               w22(n)=xf*yf
            else
               n11(n)=0
               n21(n)=0
               n12(n)=0
               n22(n)=0
            endif
         else
            n11(n)=0
            n21(n)=0
            n12(n)=0
            n22(n)=0
         endif
      enddo

   end subroutine bilinear_interp_input_usaf

   ! Internal subroutine to set up weights for neighbor interpolation.
   ! Based on neighbor_interp_input, but handles Air Force 16th
   ! mesh polar stereographic grids.
   subroutine neighbor_interp_input_usaf(griddesci, griddesco, npts, &
        rlat2, rlon2, n112, afwa_grid)

      ! Defaults
      implicit none

      ! Arguments
      real, intent(in) :: griddesci(50)
      real, intent(in) :: griddesco(50)
      integer, intent(in) :: npts
      real, intent(inout) :: rlat2(npts)
      real, intent(inout) :: rlon2(npts)
      integer, intent(inout) :: n112(npts)
      character(len=*), intent(in) :: afwa_grid

      ! Local variables
      integer             :: n
      integer             :: mo, nv
      real                :: xpts(npts), ypts(npts)
      integer             :: i1, j1
      real                :: xi, yi
      integer, external   :: get_fieldpos
      real, parameter     :: fill = -9999.0

      mo = npts
      if (trim(afwa_grid) .ne. NH_PS16 .and. &
           trim(afwa_grid) .ne. SH_PS16) then
         if (gridDesco(1) .ge. 0) then
            call compute_earth_coord(griddesco, mo, fill, xpts, ypts, &
                 rlon2, rlat2, nv)
         end if
      end if
      call compute_grid_coord(griddesci, mo, fill, xpts, ypts, &
           rlon2, rlat2, nv)
      do n=1, mo
         xi = xpts(n)
         yi = ypts(n)
         if (xi .ne. fill .and. yi .ne. fill) then
            i1 = nint(xi)
            j1 = nint(yi)
            n112(n) = get_fieldpos(i1, j1, griddesci)
         else
            n112(n) = 0
         endif
      enddo

   end subroutine neighbor_interp_input_usaf

   ! Internal function for setting GRIB1 grid definition
   function set_grid_definition(gridID) result (grid_definition)
      character(len=*), intent(in) :: gridID
      integer :: grid_definition
      grid_definition = 255
      if (trim(gridID) .eq. NH_PS16) then
         grid_definition = 212
      else if (trim(gridID) .eq. SH_PS16) then
         grid_definition = 213
      end if
   end function set_grid_definition

end module LVT_USAFSIpostMod
