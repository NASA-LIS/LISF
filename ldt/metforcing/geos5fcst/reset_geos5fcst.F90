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
! !MODULE: reset_geos5fcst
!  \label{reset_geos5fcst}
!
! !REVISION HISTORY: 
! 7 Mar 2013: Sujay Kumar, initial specification
! 
! !INTERFACE:
subroutine reset_geos5fcst()
! !USES:
  use LDT_coreMod,       only : LDT_rc
  use geos5fcst_forcingMod, only : geos5fcst_struc
!
! !DESCRIPTION:
!  Routine to reset GEOS5 forecast forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LDT_rc%nnest
     geos5fcst_struc(n)%fcsttime1 = 3000.0
     geos5fcst_struc(n)%fcsttime2 = 0.0
  enddo

end subroutine reset_geos5fcst
