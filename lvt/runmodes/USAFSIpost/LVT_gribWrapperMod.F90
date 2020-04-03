!-----------------------BEGIN NOTICE -- DO NOT EDIT----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT----------------------------
!
! MODULE: LVT_gribWrapperMod
!
! REVISION HISTORY:
! 16 May 2019  Eric Kemp  Initial version
!
! DESCRIPTION:
! Provides some error-checking wrapper routines for GRIB_API/ECCODES.
!------------------------------------------------------------------------------

module LVT_gribWrapperMod
   
   implicit none
   private

   public :: LVT_grib_set

   interface LVT_grib_set
      module procedure LVT_grib_set_int
      module procedure LVT_grib_set_real4
      module procedure LVT_grib_set_real4_array
      module procedure LVT_grib_set_string
   end interface
contains

   subroutine LVT_grib_set_int(gribid, key, value)
      use grib_api
      use LVT_logMod, only: LVT_logunit, LVT_endrun
      implicit none
      integer, intent(in) :: gribid
      character(len=*), intent(in) :: key
      integer, intent(in) :: value
      integer :: status1, status2
      character(len=255) :: msg
      call grib_set(gribid, key, value, status1)
      if (status1 .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_set for key ',trim(key), &
              ', value ',value
         call grib_get_error_string(status1, msg, status2)
         write(LVT_logunit,*)'[ERR] ',trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if
   end subroutine LVT_grib_set_int

   subroutine LVT_grib_set_real4(gribid, key, value)
      use grib_api
      use LVT_logMod, only: LVT_logunit, LVT_endrun
      implicit none
      integer, intent(in) :: gribid
      character(len=*), intent(in) :: key
      real, intent(in) :: value
      integer :: status1, status2
      character(len=255) :: msg
      call grib_set(gribid, key, value, status1)
      if (status1 .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_set for key ',trim(key), &
              ', value ',value
         call grib_get_error_string(status1, msg, status2)
         write(LVT_logunit,*)'[ERR] ',trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if
   end subroutine LVT_grib_set_real4

   subroutine LVT_grib_set_real4_array(gribid, key, value)
      use grib_api
      use LVT_logMod, only: LVT_logunit, LVT_endrun
      implicit none
      integer, intent(in) :: gribid
      character(len=*), intent(in) :: key
      real, dimension(:), intent(in) :: value
      integer :: status1, status2
      character(len=255) :: msg
      call grib_set(gribid, key, value, status1)
      if (status1 .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_set for key ',trim(key), &
              ', value ',value
         call grib_get_error_string(status1, msg, status2)
         write(LVT_logunit,*)'[ERR] ',trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if
   end subroutine LVT_grib_set_real4_array

   subroutine LVT_grib_set_string(gribid, key, value)
      use grib_api
      use LVT_logMod, only: LVT_logunit, LVT_endrun
      implicit none
      integer, intent(in) :: gribid
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: value
      integer :: status1, status2
      character(len=255) :: msg
      call grib_set(gribid, key, value, status1)
      if (status1 .ne. GRIB_SUCCESS) then
         write(LVT_logunit,*)'[ERR] Error from grib_set for key ',trim(key), &
              ', value ',value
         call grib_get_error_string(status1, msg, status2)
         write(LVT_logunit,*)'[ERR] ',trim(msg)
         write(LVT_logunit,*)'[ERR] LVT will stop'
         call LVT_endrun()
      end if
   end subroutine LVT_grib_set_string

end module LVT_gribWrapperMod
