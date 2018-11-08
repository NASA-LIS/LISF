!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT)
!
! See RELEASE_NOTES.txt for more information.
!
! The LDT source code and documentation are not in the public domain
! and may not be freely distributed.  Only qualified entities may receive 
! the source code and documentation. 
!
! Qualified entities must be covered by a Software Usage Agreement. 
! The Software Usage Agreement contains all the terms and conditions
! regarding the release of the LDT software.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See the Software Usage Agreement for the full disclaimer of warranty.
!
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: geowrsi2_set_gcoords
!  \label{geowrsi2_set_gcoords}
!
! !REVISION HISTORY:
!  28 Jun 2013: James Geiger; Initial implementation
!  25 Oct 2013: KR Arsenault;  Added GeoWRSI2.0 model to LDT
!
! !INTERFACE:
subroutine geowrsi2_set_gcoords(n)

! !USES:
   use LDT_coreMod, only : LDT_rc
   use fbil_module

   implicit none

! !ARGUMENTS: 
   integer :: n

! !DESCRIPTION:
!  This routine sets the gCoords data structure.  This data
!  structure contains the running domain specifications that are
!  used by the BIL routines.
!
!  Note that the BIL library expects grid-cells to be defined
!  by their upper left corner.  LDT defines grid-cells on their
!  centers.  When setting gCoords, nudge the values to conform
!  to the BIL library.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

 ! Attention: gCoords must have been allocated by the calling routine.
   gCoords%minLat = LDT_rc%gridDesc(n, 4) + LDT_rc%gridDesc(n,10)/2.0
   gCoords%maxLat = LDT_rc%gridDesc(n, 7) + LDT_rc%gridDesc(n,10)/2.0
   gCoords%minLon = LDT_rc%gridDesc(n, 5) - LDT_rc%gridDesc(n,9)/2.0
   gCoords%maxLon = LDT_rc%gridDesc(n, 8) - LDT_rc%gridDesc(n,9)/2.0
   gCoords%pixLat = LDT_rc%gridDesc(n, 10)
   gCoords%pixLon = LDT_rc%gridDesc(n, 9)

end subroutine geowrsi2_set_gcoords
