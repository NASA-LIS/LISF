!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: finalize_PALSmetdata
!  \label{finalize_PALSmetdata}
!
! !REVISION HISTORY: 
! 7 Mar 2013: Sujay Kumar, initial specification
! 
! !INTERFACE:
subroutine finalize_PALSmetdata(findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use PALSmetdata_forcingMod, only : PALSmetdata_struc
!
! !DESCRIPTION:
!  Routine to finalize PALS forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex
  
 deallocate(PALSmetdata_struc)

end subroutine finalize_PALSmetdata
