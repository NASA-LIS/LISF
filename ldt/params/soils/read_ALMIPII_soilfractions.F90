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
! !ROUTINE: read_ALMIPII_soilfractions
!  \label{read_ALMIPII_soilfractions}
!
! !REVISION HISTORY:
!  28 Jun 2012: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_ALMIPII_soilfractions(n,num_bins, soilsfgrd, &
     sandave, clayave, siltave)
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
  integer, intent(in)   :: n          ! nest index
  integer, intent(in)   :: num_bins   ! number of bins for tiling
  real,    intent(out)  :: soilsfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: sandave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: clayave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: siltave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

!
! !DESCRIPTION:
!  This subroutine reads the ALMIPII landcover data and returns the 
!  distribution of vegetation in each grid cell, in a lat/lon
!  projection. The data has 12 vegetation types.   
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[sand]
!    sand fraction values
!   \end{description}
!EOP      

  real :: isum
  integer :: ftn
  integer :: nc,nr
  integer :: latid, lonid, sandid,clayid
  integer :: t, c,r,ierr, ios1
  logical :: file_exists
  real    :: sand_f(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: clay_f(LDT_rc%lnc(n),LDT_rc%lnr(n))


  inquire(file=(LDT_rc%safile(n)),exist=file_exists) 
  if(.not.file_exists) then 
     write(LDT_logunit,*) 'sand fraction map: ',(trim(LDT_rc%safile(n))),&
          ' does not exist'
     write(LDT_logunit,*) 'program stopping ...'
     call LDT_endrun
  endif

  write(LDT_logunit,*) 'Reading ALMIPII sand fraction data ',&
       (LDT_rc%safile(n))

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  ierr = nf90_open(path=(LDT_rc%safile(n)),mode=NF90_NOWRITE,ncid=ftn)
  call LDT_verify(ierr,'error opening ALMIPII landcover data')

  ierr = nf90_inq_varid(ftn,'sand',sandid)
  call LDT_verify(ierr, 'nf90_inq_varid failed for sand in read_ALMIPII_soilfractions')
  
  ierr = nf90_inq_dimid(ftn,'longitude',lonid)
  call LDT_verify(ierr,'nf90_inq_dimid failed for longitude in read_ALMIPII_soilfractions')

  ierr = nf90_inq_dimid(ftn,'latitude',latid)
  call LDT_verify(ierr,'nf90_inq_dimid failed for latitude in read_ALMIPII_soilfractions')

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

  ierr = nf90_get_var(ftn,sandid, sand_f)
  call LDT_verify(ierr, 'nf90_get_var failed for sand')
  
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
!        sandave(c,r,1) = sand_f(c,LDT_rc%lnr(n)-r+1)
        sandave(c,r,1) = sand_f(c,r)
     enddo
  enddo
  ierr = nf90_close(ftn)
  call LDT_verify(ierr, 'nf90_close failed in read_ALMIPII')


  inquire(file=(LDT_rc%clfile(n)),exist=file_exists) 
  if(.not.file_exists) then 
     write(LDT_logunit,*) 'clay fraction map: ',(trim(LDT_rc%clfile(n))),&
          ' does not exist'
     write(LDT_logunit,*) 'program stopping ...'
     call LDT_endrun
  endif

  write(LDT_logunit,*) 'Reading ALMIPII sand fraction data ',&
       (LDT_rc%clfile(n))


  ierr = nf90_open(path=(LDT_rc%clfile(n)),mode=NF90_NOWRITE,ncid=ftn)
  call LDT_verify(ierr,'error opening ALMIPII landcover data')

  ierr = nf90_inq_varid(ftn,'clay',clayid)
  call LDT_verify(ierr, 'nf90_inq_varid failed for clay in read_ALMIPII_soilfractions')
  
  ierr = nf90_inq_dimid(ftn,'longitude',lonid)
  call LDT_verify(ierr,'nf90_inq_dimid failed for longitude in read_ALMIPII_soilfractions')

  ierr = nf90_inq_dimid(ftn,'latitude',latid)
  call LDT_verify(ierr,'nf90_inq_dimid failed for latitude in read_ALMIPII_soilfractions')

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

  ierr = nf90_get_var(ftn,clayid, clay_f)
  call LDT_verify(ierr, 'nf90_get_var failed for clay')
  
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
!        clayave(c,r,1) = clay_f(c,LDT_rc%lnr(n)-r+1)
        clayave(c,r,1) = clay_f(c,r)
     enddo
  enddo

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        siltave(c,r,1) = 1.0 - sandave(c,r,1) - clayave(c,r,1)
     enddo
  enddo
  ierr = nf90_close(ftn)
  call LDT_verify(ierr, 'nf90_close failed in read_ALMIPII')

  soilsfgrd = 0.0
  soilsfgrd(:,:,1) = 1.0
 
#endif
end subroutine read_ALMIPII_soilfractions
