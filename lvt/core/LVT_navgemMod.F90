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

  ! Modules
  use LVT_logMod, only: LVT_logunit

  ! Defaults
  implicit none
  private

  ! Public routines
  public :: LVT_get_navgem_filename

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

  subroutine LVT_get_navgem_filename(filename, &
       year, month, day, hour, fcst_hr)

    ! Modules
    use LVT_coreMod, only: LVT_rc
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

  end subroutine LVT_get_navgem_filename
end module LVT_navgemMod
