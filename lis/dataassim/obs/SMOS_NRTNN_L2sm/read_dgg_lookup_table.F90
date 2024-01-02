!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!#include "LIS_NetCDF_inc.h"
module read_dgg_lookup_table
!BOP
! !MODULE: read_dgg_lookup_table
!
! !DESCRIPTION:
!  The code in this file read the SMOS DGG lookup table from a netCDF file
!
! !REVISION HISTORY:
!  15 Feb 2021: Mahdi Navari; Initial version
!
   use LIS_constantsMod, only : LIS_CONST_PATH_LEN
   use LIS_logMod
!#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
!  use netcdf
!#endif
   use SMOSNRTNNL2sm_Mod

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_SMOS_DGG_lookup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

contains

!BOP
!
! !ROUTINE: LIS_SMOS_DGG_lookup
! \label{LIS_SMOS_DGG_lookup}
!
 subroutine LIS_SMOS_DGG_lookup(n, dgg_lookup_1d)

! USES:
   use LIS_logMod,        only : LIS_endrun
   use LIS_coreMod,       only : LIS_rc
   use SMOSNRTNNL2sm_Mod

   implicit none

   ! Arguments
    integer, intent(in)       :: n
    integer, allocatable, intent(inout) :: dgg_lookup_1d(:)

    ! Locals
    character(len=LIS_CONST_PATH_LEN) :: fname, fname_tmp
    integer       :: ncid, leng
    logical       :: file_exists

   fname_tmp = SMOSNRTNNL2sm_struc(n)%obscdffile
   leng = LEN_TRIM(fname_tmp)
   fname = fname_tmp(1:leng-3)//'_DGG_lookupTable.nc'

   inquire(file=fname,exist=file_exists)
   if(file_exists) then

      call read_dgg_netcdf(fname, dgg_lookup_1d)

   else
      write(LIS_logunit,*) '[ERR] input file '//trim(fname)
      write(LIS_logunit,*) '[ERR] does not exist '
      write(LIS_logunit,*) '[ERR] failed in LIS_SMOS_DGG_lookup'
      call LIS_endrun()
   endif

  end subroutine LIS_SMOS_DGG_lookup


  subroutine read_dgg_netcdf(filename, dgg_lookup_1d)

      ! Imports
      use LIS_logMod, only: LIS_verify , LIS_logunit
      use netcdf

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in)       :: filename
      integer, allocatable, intent(inout) :: dgg_lookup_1d(:)

      ! Locals
      integer                            :: leng
      integer                            :: ncid, varid, dimid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

       call LIS_verify(nf90_open(filename, NF90_NOWRITE, ncid=ncid), &
            'nf90_open failed in read_dgg_netcdf')

       call LIS_verify(nf90_inq_dimid(ncid, "x",dimid),&
            'nf90_inq_dimid failed in read_dgg_netcdf')

       call LIS_verify(nf90_inquire_dimension(ncid, dimid, len=leng),&
           'nf90_inquire_dimension failed in read_dgg_netcdf')

       allocate(dgg_lookup_1d(leng))

       call LIS_verify(nf90_inq_varid(ncid,"dgg_lookup" ,varid), &
            'nf90_inq_varid failed in read_dgg_netcdf')

       call LIS_verify(nf90_get_var(ncid, varid, dgg_lookup_1d), &
            'nf90_get_var failed in read_dgg_netcdf')

       call LIS_verify(nf90_close(ncid),&
            'nf90_close failed in read_dgg_netcdf')

       write(LIS_logunit,*) 'Successfully read dgg_lookup file ',trim(filename)
#endif

  end subroutine read_dgg_netcdf

end module read_dgg_lookup_table
