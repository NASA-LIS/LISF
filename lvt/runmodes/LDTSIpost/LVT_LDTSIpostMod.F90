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
      
   end type LVT_LDTSIpost_t

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

end module LVT_LDTSIpostMod
