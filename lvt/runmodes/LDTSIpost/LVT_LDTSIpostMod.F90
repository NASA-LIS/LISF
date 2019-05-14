!-----------------------BEGIN NOTICE -- DO NOT EDIT----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT----------------------------
!
! MODULE: LVT_LDTSIpostMod
!
! REVISION HISTORY:
! 13 May 2019  Eric Kemp  Initial version
!
! DESCRIPTION:
! Source code for reading LDTSI netCDF files, interpolating to predefined
! Air Force grids, and outputting GRIB2 files.
!------------------------------------------------------------------------------

module LVT_LDTSIpostMod

   implicit none
   private

   type, public :: LVT_LDTSIpost_t
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
      procedure read_ldtsi_ncfile
      procedure interp_and_output_grib2

   end type LVT_LDTSIpost_t

   character(len=50), parameter :: GLOBAL_LL0P25 = 'global_ll0p25'
   character(len=50), parameter :: NH_PS16 = 'nh_ps16'
   character(len=50), parameter :: SH_PS16 = 'sh_ps16'

contains
   
   ! Constructor for LVT_LDTSIpost_t object
   subroutine new(this)
      
      ! Imports
      use LVT_coreMod, only: LVT_rc
      use LVT_logMod, only: LVT_logunit, LVT_endrun

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_LDTSIpost_t), intent(inout) :: this

      ! Local variables
      logical :: file_exists

      ! Construct the input netCDF filename
      this%input_nc_file = &
           trim(LVT_rc%input_dir) // &
           "/ldtsi_" // &
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

   ! Destructor of LVT_LDTSIpost_t object
   subroutine delete(this)

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_LDTSIpost_t), intent(inout) :: this

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

   ! Method for reading LDTSI netCDF file.
   subroutine read_ldtsi_ncfile(this)

      ! Imports
      use LVT_logMod, only: LVT_logunit, LVT_endrun
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
      use netcdf
#endif

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_LDTSIpost_t), intent(inout) :: this

      ! Local variables
      integer :: ncid
      integer :: dims_ids(3)
      integer :: ntime, nlat, nlon
      real, allocatable :: tmp(:,:,:)
      integer :: snoanl_varid, snoage_varid, icecon_varid, icemask_varid, &
           iceage_varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
      write(LVT_logunit,*)'[INFO] Reading LDTSI file ',trim(this%input_nc_file)
      
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
              '[ERR] Unrecognized map projection found in LDTSI file!'
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

      ! Fetch the LDTSI analysis variable IDs
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

      ! Read the LDTSI variables
      allocate(tmp(nlon, nlat, ntime)) ! Need 3D array
      call LVT_verify(nf90_get_var(ncid=ncid, &
           varid=snoanl_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for snoanl')
      this%snoanl(:,:) = tmp(:,:,1)

      call LVT_verify(nf90_get_var(ncid=ncid, &
           varid=snoage_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for snoage')
      this%snoage(:,:) = tmp(:,:,1)

      call LVT_verify(nf90_get_var(ncid=ncid, &
           varid=icecon_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for icecon')
      ! A unit transform is needed here
      do r = 1, nlat
         do c = 1, nlon
            if (tmp(c,r,1) < 0) then
               this%icecon(c,r) = -1
            else
               this%icecon(c,r) = 100*tmp(c,r,1)
            end if
         end do
      end do

      call LVT_verify(nf90_get_var(ncid=ncid, &
           varid=icemask_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for icemask')
      this%icemask(:,:) = tmp(:,:,1)

      call LVT_verify(nf90_get_var(ncid=ncid, &
           varid=iceage_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for icemask')
      this%iceage(:,:) = tmp(:,:,1)

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
           '[ERR] Error in nf90_close for ' // trim(this%input_nc_file)

#else
      write(LVT_logunit,*) &
           '[ERR] Must compile LVT with netCDF support for LDTSIpost mode!'
      write(LVT_logunit,*) &
           '[ERR] LVT will exit gracefully.'
      call LVT_endrun()
#endif

   end subroutine read_ldtsi_ncfile

   ! Interpolate and output LDTSI data in GRIB2 format.
   subroutine interp_and_output_grib2(this, gridID)

      ! Imports
      use grib_api
      use LVT_coreMod, only: LVT_rc
      use LVT_logMod, only: LVT_logunit, LVT_endrun

      ! Defaults
      implicit none

      ! Arguments
      class(LVT_LDTSIpost_t), intent(inout) :: this
      character(len=*), intent(in) :: gridID

      ! Local variables
      real :: griddesci(50), griddesco(50)
      real :: xmesh, xpnmcaf, ypnmcaf, orient, xj, xi, alat, alon
      integer :: nc_out, nr_out
      integer, allocatable :: n11(:)
      logical*1, allocatable :: li(:), lo(:)
      real, allocatable :: gi(:), go(:)
      character(len=255) :: fname
      integer :: ftn, rc
      integer :: c,r
      integer :: igrib
      character(len=8) :: yyyymmdd
      character(len=4) :: hhmm
      integer :: idate8, idate4
      integer :: gridDefinitionTemplateNumber

      ! Sanity check the gridID. 
      call check_gridID(gridID)

      ! Set the LDTSI grid description
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
      griddesci(20) = 255
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

      ! Set the output grid description
      if (trim(gridID) .eq. trim(GLOBAL_LL0P25)) then
         gridDefinitionTemplateNumber = 0
         nc_out = 1440
         nr_out =  720
         griddesco(:) = 0
         griddesco(1) = 0 ! lat/lon
         griddesco(2) = nc_out
         griddesco(3) = nr_out
         griddesco(4) =  -89.875000
         griddesco(5) = -179.875000
         griddesco(6) = 128
         griddesco(7) =   89.875000
         griddesco(8) =  179.875000
         griddesco(9) =  0.250000
         griddesco(10) = 0.250000
         griddesco(11) = 64
         griddesco(20) = 255
         griddesco(30) = 0
         griddesco(32) = nc_out
         griddesco(33) = nr_out
         griddesco(34) =  -89.875000
         griddesco(35) = -179.875000
         griddesco(36) = 128
         griddesco(37) =   89.875000
         griddesco(38) =  179.875000
         griddesco(39) = 0.250000
         griddesco(40) = 0.250000
         
      else if (trim(gridID) .eq. trim(NH_PS16)) then

         gridDefinitionTemplateNumber = 20
         xmesh = 47.625/2
         xpnmcaf = 513
         ypnmcaf = 513
         orient = 100.0
         xj = float(1) - ypnmcaf
         xi = float(1) - xpnmcaf
         call polarToLatLon(xi,xj,xmesh,orient,alat,alon)
         
         nc_out = 1024
         nr_out = 1024
         griddesco(:) = 0
         griddesco(1) = 5
         griddesco(2) = nc_out
         griddesco(3) = nr_out
         griddesco(4) = alat
         griddesco(5) = alon
         griddesco(6) = 8
         griddesco(7) = orient
         griddesco(8) = xmesh
         griddesco(9) = xmesh
         griddesco(10) = 0.0
         griddesco(11) = 128
         griddesco(13) = 1
         griddesco(20) = 128
         
      else if (trim(gridID) .eq. trim(SH_PS16)) then

         gridDefinitionTemplateNumber = 20

         xmesh = -1*47.625/2
         xpnmcaf = 513
         ypnmcaf = 513
         orient = 280
         xj = float(1) - ypnmcaf
         xi = float(1) - xpnmcaf
         call polarToLatLon(xi,xj,xmesh,orient,alat,alon)

         nc_out = 1024
         nr_out = 1024         
         griddesco(:) = 0
         griddesco(1) = 5
         griddesco(2) = nc_out
         griddesco(3) = nr_out
         griddesco(4) = alat
         griddesco(5) = alon
         griddesco(6) = 8
         griddesco(7) = orient
         griddesco(8) = xmesh
         griddesco(9) = xmesh
         griddesco(10) = 0.0
         griddesco(11) = 128
         griddesco(13) = 1
         griddesco(20) = 128

      end if

      ! Calculate neighbor weights for upscaling
      allocate(n11(this%nc*this%nr))
      call upscaleByAveraging_input(griddesci, griddesco, &
           (this%nc*this%nr), (nc_out*nr_out), n11)

      ! Construct the GRIB2 filename
      call build_filename_g2(gridID, LVT_rc%yyyymmddhh, fname)
      
      ! Open the GRIB2 file
      call grib_open_file(ftn, fname, 'w', rc)

      ! Allocate memory for interpolation
      allocate(li(this%nc*this%nr))
      allocate(gi(this%nr*this%nr))
      allocate(lo(nc_out*nr_out))
      allocate(go(nc_out*nr_out))

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
      call upscaleByAveraging((this%nc*this%nr), &
           (nc_out*nr_out), LVT_rc%udef, li, gi, lo, go)

#if (defined USE_GRIBAPI)
      call grib_new_from_template(igrib, "GRIB2", rc)
#else
      call grib_new_from_samples(igrib, "GRIB2", rc)
#endif

      ! Section 0: Indicator
      ! Octet 7
      call grib_set(igrib, 'discipline', 0, rc)

      ! Section 1: Identification
      ! Octets 6-7
      call grib_set(igrib, 'centre', LVT_rc%grib_center_id, rc)
      ! Octets 8-9
      call grib_set(igrib, 'subCentre', LVT_rc%grib_subcenter_id, rc)
      ! Octet 10
      call grib_set(igrib, 'tablesVersion', LVT_rc%grib_table, rc)
      ! Octet 11
      call grib_set(igrib, 'localTabelsVersion', 1, rc)
      ! Octet 12
      call grib_set(igrib, 'significanceOfReferenceTime', 0, rc)
      ! Octet 13-16
      yyyymmdd = LVT_rc%yyyymmddhh(1:8)
      read(idate8,'(i8)') yyyymmdd
      call grib_set(igrib, 'dataDate', idate8, rc)
      ! Octet 17-19
      hhmm = LVT_rc%yyyymmddhh(9:10) // '00'
      read(idate4,'(i4)') hhmm
      call grib_set(igrib, 'dataTime', idate4, rc)
      ! Octet 20
      call grib_set(igrib, 'productionStatusOfProcessedData', 0, rc)
      ! Octet 21
      call grib_set(igrib, 'typeOfProcessedData', 0, rc)

      ! Section 2:  Local Use Section (Optional) -- none for now

      ! Section 3: Grid
      call grib_set(igrib, 'gridDefinitionTemplateNumber', &
           gridDefinitionTemplateNumber, rc)
      ! Shape of the earth.  (Spherical earth, radius = 6,367,470.0 m)
      call grib_set(igrib, 'shapeOfTheEarth', 0, rc)
      ! Change data order
      call grib_set(igrib, 'swapScanningLat', 1, rc)
      ! Set dimensions
      if (gridDefinitionTemplateNumber == 0) then
         call grib_set(igrib, 'Ni', nc_out, rc)
         call grib_set(igrib, 'Nj', nr_out, rc)
         call grib_set(igrib, 'latitudeOfFirstGridPointInDegrees', &
              -89.875000, rc)
         call grib_set(igrib, 'longitudeOfFirstGridPointInDegrees', &
              -179.875000, rc)
         call grib_set(igrib, 'latitudeOfLastGridPointInDegrees', &
              89.875000, rc)
         call grib_set(igrib, 'longitudeOfLastGridPointInDegrees', &
              179.875000, rc)
         call grib_set(igrib, 'gridType', 'regular_ll', rc)
         call grib_set(igrib, 'iDirectionIncrementInDegrees', &
              griddesco(9), rc)
         call grib_set(igrib, 'jDirectionIncrementInDegrees', &
              griddesco(10), rc)
      else if (gridDefinitionTemplateNumber == 0) then
         call grib_set(igrib, 'Nx', nc_out, rc)
         call grib_set(igrib, 'Ny', nr_out, rc)
      end if

      ! Section 4: Product Definition Section
      ! Octets 8-9
      call grib_set(igrib, 'productionDefinitionTemplateNumber', 0, rc)
      ! Octet 10
      call grib_set(igrib, 'parameterCategory', 1, rc)
      ! Octet 11
      call grib_set(igrib, 'parameterNumber', 11, rc)
      ! Octet 12
      call grib_set(igrib, 'typeOfGeneratingProcess', 0, rc)
      ! Octet 13...Mark as SNODEP (unmodified)
      call grib_set(igrib, 'backgroundGeneratingProcessIdentifier', &
           35, rc)
      ! Octet 14...Mark as SNODEP (unmodified)
      call grib_set(igrib, 'generatingProcessIdentifier', 35, rc)
      ! Octet 15-17 is skipped
      ! Octet 18...Use hours
      call grib_set(igrib, 'indicatorOfUnitOfTimeRange', 1, rc)
      ! Octets 19-22 
      call grib_set(igrib, 'forecastTime', 0, rc)
      ! Octets 23-34
      call grib_set(igrib, 'typeOfFirstFixedSurface', 1, rc)
      call grib_set(igrib, 'typeOfSecondFixedSurface', 255, rc)
      call grib_set(igrib, 'scaleFactorOfFirstFixedSurface', 0, rc)
      call grib_set(igrib, 'scaledValueOfFirstFixedSurface', 0, rc)
      call grib_set(igrib, 'scaleFactorOfSecondFixedSurface', 255, rc)
      call grib_set(igrib, 'scaledValueOfSecondFixedSurface', 255, rc)

      ! Section 5: Data Representation
      call grib_set(igrib, 'packingType', LVT_rc%grib_packing_type, rc)
      call grib_set(igrib, 'missingValue', LVT_rc%udef, rc)

      ! Section 6: Bit-Map
      call grib_set(igrib, 'bitmapPresent', 1, rc)

      ! Section 7
      call grib_set(igrib, 'values', go, rc)
      call grib_write(igrib, ftn, rc)
      call grib_release(igrib, rc)

      ! Close the GRIB2 file
      call grib_close_file(ftn, rc)

      ! Clean up
      deallocate(n11)
      deallocate(li)
      deallocate(gi)
      deallocate(go)


   end subroutine interp_and_output_grib2

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
   subroutine build_filename_g2(gridID, yyyymmddhh, filename)

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: gridID
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

      filename = &
           'PS.AFWA_SC.U_DI.D_GP.SNODEP_GR.' // &
           trim(grid) // '_AR.' // &
           trim(area) // '_PA.SNODEP_DD' // &
           trim(yyyymmdd) // '_DT.' // &
           trim(hh) // '00_DF.GR2'
           
   end subroutine build_filename_g2
end module LVT_LDTSIpostMod
