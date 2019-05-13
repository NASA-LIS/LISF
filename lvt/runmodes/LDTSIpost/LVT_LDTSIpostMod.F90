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


end module LVT_LDTSIpostMod
