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
!BOP
!
! !ROUTINE: read_ALMIPII_dsoil
!  \label{read_ALMIPII_dsoil}
!
! !REVISION HISTORY:
!  09 Jul 2012: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_ALMIPII_dsoil(n,dsoil)
! !USES:
  use LDT_coreMod,        only : LDT_rc
  use LDT_logMod,         only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_pluginIndices
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer :: n  
  real    :: dsoil(LDT_rc%lnc(n),LDT_rc%lnr(n))

!
! !DESCRIPTION:
!  This subroutine reads the ALMIPII soil depth data
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[dsoil]
!    soil depth values
!   \end{description}
!EOP      

  real :: isum
  integer :: ftn
  integer :: nc,nr
  integer :: latid, lonid, dsoilid
  integer :: t, c,r,ierr, ios1
  logical :: file_exists
  real    :: dsoil_t(LDT_rc%lnc(n),LDT_rc%lnr(n))


  inquire(file=trim(LDT_rc%dsoilfile(n)),exist=file_exists) 
  if(.not.file_exists) then 
     write(LDT_logunit,*) 'soil depth map: ',trim(LDT_rc%dsoilfile(n)),&
          ' does not exist'
     write(LDT_logunit,*) 'program stopping ...'
     call LDT_endrun
  endif

  write(LDT_logunit,*) 'Reading ALMIPII soil depth data ',&
       trim(LDT_rc%dsoilfile(n))

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  ierr = nf90_open(path=trim(LDT_rc%dsoilfile(n)),mode=NF90_NOWRITE,ncid=ftn)
  call LDT_verify(ierr,'error opening ALMIPII soil depth data')

  ierr = nf90_inq_varid(ftn,'dsoil',dsoilid)
  call LDT_verify(ierr, 'nf90_inq_varid failed for dsoil in read_ALMIPII_dsoil')
  
  ierr = nf90_inq_dimid(ftn,'longitude',lonid)
  call LDT_verify(ierr,'nf90_inq_dimid failed for longitude in read_ALMIPII_dsoil')

  ierr = nf90_inq_dimid(ftn,'latitude',latid)
  call LDT_verify(ierr,'nf90_inq_dimid failed for latitude in read_ALMIPII_dsoil')

  ierr = nf90_inquire_dimension(ftn,lonid,len=nc)
  call LDT_verify(ierr,'nf90_inquire_dimension for longitude')

  ierr = nf90_inquire_dimension(ftn,latid,len=nr)
  call LDT_verify(ierr,'nf90_inquire_dimension for latitude')


  if(nc.ne.LDT_rc%lnc(n).or.nr.ne.LDT_rc%lnr(n)) then 
     write(LDT_logunit,*) 'ALMIPII domain dimensions do not match the'
     write(LDT_logunit,*) 'LIS domain dimensions '
     write(LDT_logunit,*) 'Program stopping .... '
     call LDT_endrun()
  endif

  ierr = nf90_get_var(ftn,dsoilid, dsoil_t)
  call LDT_verify(ierr, 'nf90_get_var failed for dsoil')

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
!        dsoil(c,r) = dsoil_t(c,LDT_rc%lnr(n)-r+1)
        dsoil(c,r) = dsoil_t(c,r)
     enddo
  enddo
  ierr = nf90_close(ftn)
  call LDT_verify(ierr, 'nf90_close failed in read_ALMIPII_dsoil')

#endif
end subroutine read_ALMIPII_dsoil
