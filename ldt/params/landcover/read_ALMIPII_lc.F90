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
! !ROUTINE: read_ALMIPII_lc
!  \label{read_ALMIPII_lc}
!
! !REVISION HISTORY:
!  28 Jun 2012: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_ALMIPII_lc(n, num_types, fgrd, maskarray)
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
  integer, intent(in) :: num_types
  real    :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_types)
  real, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))
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
!   \item[fgrd]
!    fraction of grid covered by each vegetation type
!   \item[maskarray]
!     landmask for the region of interest
!   \end{description}
!EOP      

  real :: isum
  integer :: ftn
  integer :: nc,nr
  integer :: latid, lonid, lcid
  integer :: t, c,r,ierr, ios1
  logical :: file_exists
  real    :: patch(LDT_rc%lnc(n),LDT_rc%lnr(n),num_types)

  fgrd(:,:,:) = 0.0

  if( LDT_rc%lc_type(n) == "ECOCLIMAP2" ) then
     LDT_rc%bareclass    = 1
     LDT_rc%urbanclass   = 0
     LDT_rc%snowclass    = 0
     LDT_rc%waterclass   = 0
     LDT_rc%wetlandclass = 0
     LDT_rc%glacierclass = 0
  else  ! not supported options
     write(LDT_logunit,*) ' The land classification: ',&
          trim(LDT_rc%lc_type(n)),' does not exist for ALMIPII source.'
      write(LDT_logunit,*) " -- Please select:  ECOCLIMAP2 "
     write(LDT_logunit,*) 'Program stopping ...'
     call LDT_endrun
  endif
  
  inquire(file=trim(LDT_rc%vfile(n)),exist=file_exists) 
  if(.not.file_exists) then 
     write(LDT_logunit,*) 'landcover map: ',trim(LDT_rc%vfile(n)),&
          ' does not exist'
     write(LDT_logunit,*) 'program stopping ...'
     call LDT_endrun
  endif

  write(LDT_logunit,*) 'Reading ALMIPII landcover data ',trim(LDT_rc%vfile(n))

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  ierr = nf90_open(path=trim(LDT_rc%vfile(n)),mode=NF90_NOWRITE,ncid=ftn)
  call LDT_verify(ierr,'error opening ALMIPII landcover data')

  ierr = nf90_inq_varid(ftn,'patch',lcId)
  call LDT_verify(ierr, 'nf90_inq_varid failed for patch in read_ALMIPII_lc')
  
  ierr = nf90_inq_dimid(ftn,'longitude',lonid)
  call LDT_verify(ierr,'nf90_inq_dimid failed for longitude in read_ALMIPII_lc')

  ierr = nf90_inq_dimid(ftn,'latitude',latid)
  call LDT_verify(ierr,'nf90_inq_dimid failed for latitude in read_ALMIPII_lc')

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

  ierr = nf90_get_var(ftn,lcId, patch)
  call LDT_verify(ierr, 'nf90_get_var failed for patch')

  maskarray = LDT_rc%udef

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
!        fgrd(c,r,:) = patch(c,LDT_rc%lnr(n)-r+1,:)
        fgrd(c,r,:) = patch(c,r,:)
        if(sum(fgrd(c,r,:)).ne.0) then 
           maskarray(c,r) = 1.0
!           do t=1,num_types
!              if(fgrd(c,r,t).gt.0) then 
!                 sfctype(c,r,t) = 1
!              endif
!           enddo
        endif
     enddo
  enddo
  ierr = nf90_close(ftn)
  call LDT_verify(ierr, 'nf90_close failed in read_ALMIPII')

#endif
end subroutine read_ALMIPII_lc
