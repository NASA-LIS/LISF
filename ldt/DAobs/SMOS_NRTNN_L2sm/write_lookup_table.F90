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
#include "LDT_NetCDF_inc.h"
module write_lookup_table
!BOP
! !MODULE: write_lookup_table
!
! !DESCRIPTION:
!  The code in this file wirte the SMOS DGG lookup table into a netCDF file
!
! !REVISION HISTORY:
!  11 Feb 2021: Mahdi Navari; Initial version
!
  !use ESMF
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_SMOS_DGG_lookup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

contains

!BOP
! 
! !ROUTINE: LDT_SMOS_DGG_lookup
! \label{LDT_SMOS_DGG_lookup}
!
 subroutine LDT_SMOS_DGG_lookup(n)

! USES:
   use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify
   use LDT_coreMod,       only : LDT_rc, LDT_domain   
   use SMOSNRTNNL2sm_obsMod
   implicit none

   ! Arguments
    character(len=LDT_CONST_PATH_LEN) :: fname
    integer       :: n
    integer       :: iret, ncid, length
    integer, allocatable :: data(:) 


    length = size (SMOSNRTNNsmobs(n)%dgg_lookup_1d)
    allocate (data(length))
    data = SMOSNRTNNsmobs(n)%dgg_lookup_1d


    call system('mkdir -p '//(LDT_rc%odir))
         write(LDT_logunit,*) "Writing to LDT output directory: ",&
         trim(LDT_rc%odir)

         fname = trim(LDT_rc%odir)//'/'//&
               trim(LDT_rc%dapreprocfile)//'_DGG_lookupTable.nc'


    call write_netcdf(fname, length, data)
    deallocate(data)
    
  end subroutine LDT_SMOS_DGG_lookup


  subroutine write_netcdf(filename, leng, dgg_lookup)

      ! Imports
      use LDT_logMod, only: LDT_verify
      use netcdf

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in)       :: filename
      integer, parameter                 :: NDIMS = 1
      integer                            :: NX, leng
      integer                            :: dgg_lookup(1:leng)
      integer                            :: iret, ncid, varid, dimids(NDIMS)
      integer                            :: x_dimid

      NX = leng 

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
#if (defined USE_NETCDF3)
      write(LDT_logunit,*) 'Writing DGG lookup table file ',trim(filename)
      iret=nf90_create(path=trim(filename),cmode=nf90_clobber,&
           ncid=ncid)
#endif
#if (defined USE_NETCDF4)
      write(LDT_logunit,*) 'Writing DGG lookup table file ',trim(filename)
      iret=nf90_create(path=trim(filename),cmode=nf90_netcdf4,&
           ncid=ncid)
#endif

      call LDT_verify(nf90_def_dim(ncid, "x", NX, x_dimid)) 
      dimids =  (/x_dimid /)
      call LDT_verify(nf90_def_var(ncid, "dgg_lookup", NF90_INT, dimids, varid ),&
              'nf90_def_var failed for dgg_lookup')

       call LDT_verify(nf90_enddef(ncid))
       call LDT_verify(nf90_put_var(ncid, varid, dgg_lookup),&
         'nf90_put_att failed for dgg_lookup')
       iret=nf90_close(ncid)
       write(LDT_logunit,*) 'Successfully wrote dgg_lookup file ',trim(filename)
#endif

  end subroutine write_netcdf

end module write_lookup_table
