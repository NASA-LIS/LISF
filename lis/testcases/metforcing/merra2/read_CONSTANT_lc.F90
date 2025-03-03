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
! !ROUTINE: read_CONSTANT_lc
!  \label{read_CONSTANT_lc}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  23  May 2012: KR Arsenault; Implemented new features to read different
!                 resolutions and generate landmask from landcover
!
! !INTERFACE:
subroutine read_CONSTANT_lc(n, num_types, fgrd, maskarray)

! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_verify, LDT_endrun

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_types
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
  real, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine sets a constant landcover type for basic tests
!   or hypothetical test cases.  Also, the landmask is simply assigned.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[fgrd]
!     fraction of grid covered by each vegetation type
!   \item[maskarray]
!     landmask for the region of interest
!   \end{description}
!EOP      

   real   :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
!__________________________________________________________________

  write(LDT_logunit, *) "[INFO] Setting a CONSTANT land cover value ..."

!- Assign land cover types universal water class value: 
   select case ( LDT_rc%lc_type(n) )
    case ( "UMD" )
      LDT_rc%waterclass   = 14
    case ( "IGBP","IGBPNCEP" )
      LDT_rc%waterclass   = 17
    case ( "USGS" )
      LDT_rc%waterclass   = 16
    case ( "MOSAIC" )
      LDT_rc%waterclass   = 7
    case ( "ISA" )
      LDT_rc%waterclass   = 13   ! Originally 0
    case default  ! Non-supported options
      write(LDT_logunit,*) "The land classification: ",trim(LDT_rc%lc_type(n)),&
                           " does not exist."
      write(LDT_logunit,*) " -- Please select either: UMD, IGBP, IGBPNCEP, USGS, MOSAIC, ISA "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select

!- Initialize local variables:
   fgrd = 0.
   vegtype = 0.
   maskarray = 0.

!- LAND COVER:
!- Assign the fraction of land class for the water class layer:
! - All other fractional areas are set to 0.
!<debug -- merra-land/merr2 testing>
!   fgrd(:,:,LDT_rc%waterclass) = 1.
!</debug -- merra-land/merr2 testing>
   fgrd(:,:,1) = 1.


!- LANDMASK:
!- "READ-IN" land mask file, if user-specified:
   if( LDT_rc%mask_type(n) == "readin" ) then
     call read_maskfile( n, vegtype, fgrd, maskarray )

!- "CREATE" land mask and surface type fields (user-specified):
   elseif( LDT_rc%mask_type(n) == "create" ) then
      maskarray = 1.0
   end if


end subroutine read_CONSTANT_lc
