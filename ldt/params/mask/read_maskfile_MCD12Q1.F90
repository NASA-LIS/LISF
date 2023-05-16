!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! SUBROUTINE: MCD12Q1
!
! REVISION HISTORY:
!  03 June 2022: Sujay Kumar; Initial Specification
!
! DESCRIPTION:
! Source code for reading and interpolating land mask data in netCDF format
! produced from MCD12Q1 land use data.  The input maskfile is expected to be
! at the same exact grid and resolution of the LDT target grid.
! No spatial transformation is performed. 
!------------------------------------------------------------------------------

#include "LDT_misc.h"
 subroutine read_maskfile_MCD12Q1(n,vegtype, fgrd, localmask)

   use LDT_coreMod
   use LDT_logMod
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

   implicit none

   integer, intent(in) :: n
   real,    intent(in)  :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real,    intent(in)  :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
   real,    intent(out) :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))

   integer              :: ftn, ierr,maskId

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
   ierr = nf90_open(path=trim(LDT_rc%mfile(n)),mode=NF90_NOWRITE,ncid=ftn)
   call LDT_verify(ierr,'error opening MCD12Q1 landmask data')
   
   ierr = nf90_inq_varid(ftn,'LANDMASK',maskId)
   call LDT_verify(ierr, 'nf90_inq_varid failed for LANDMASK in read_maskfile_MCD12Q1')

   ierr = nf90_get_var(ftn,maskId, localmask)
   call LDT_verify(ierr, 'nf90_get_var failed for LANDMASK in read_maskfile_MCD12Q1')

   ierr = nf90_close(ftn)
   call LDT_verify(ierr)
#endif


 end subroutine read_maskfile_MCD12Q1
