!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: geowrsi2_set_gcoords
!  \label{geowrsi2_set_gcoords}
!
! !REVISION HISTORY:
!  28 Jun 2013: James Geiger; Initial implementation
!  25 Oct 2013: KR Arsenault;  Added GeoWRSI2.0 model to LIS-7
!
! !INTERFACE:
subroutine geowrsi2_set_gcoords(n)

! !USES:
   use LIS_coreMod, only : LIS_rc
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
!  by their upper left corner.  LIS defines grid-cells on their
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
   gCoords%minLat = LIS_rc%gridDesc(n, 4) + LIS_rc%gridDesc(n,10)/2.0
   gCoords%maxLat = LIS_rc%gridDesc(n, 7) + LIS_rc%gridDesc(n,10)/2.0
   gCoords%minLon = LIS_rc%gridDesc(n, 5) - LIS_rc%gridDesc(n,9)/2.0
   gCoords%maxLon = LIS_rc%gridDesc(n, 8) - LIS_rc%gridDesc(n,9)/2.0
   gCoords%pixLat = LIS_rc%gridDesc(n, 10)
   gCoords%pixLon = LIS_rc%gridDesc(n, 9)

end subroutine geowrsi2_set_gcoords
