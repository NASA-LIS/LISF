!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
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
!  30  Nov 2018: David Mocko; Added Bondville landcover classification
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
    case ( "Bondville" )
      LDT_rc%waterclass   = 17
      LDT_rc%urbanclass   = 13
      LDT_rc%cropclass1   = 12
      LDT_rc%snowclass    = 15
      LDT_rc%bareclass    = 16
      LDT_rc%wetlandclass = 11
      LDT_rc%glacierclass = 15
    case default  ! Non-supported options
      write(LDT_logunit,*) "The land classification: ",trim(LDT_rc%lc_type(n)),&
                           " does not exist."
      write(LDT_logunit,*) " -- Please select either: UMD, IGBP, "
      write(LDT_logunit,*) " IGBPNCEP, USGS, MOSAIC, ISA, or Bondville."
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
   fgrd(:,:,LDT_rc%waterclass) = 1.


!- LANDMASK:
!- "READ-IN" land mask file, if user-specified:
   if( LDT_rc%mask_type(n) == "readin" ) then
     call read_maskfile( n, vegtype, fgrd, maskarray )

!- "CREATE" land mask and surface type fields (user-specified):
   elseif( LDT_rc%mask_type(n) == "create" ) then
      maskarray = 1.0
   end if


end subroutine read_CONSTANT_lc
