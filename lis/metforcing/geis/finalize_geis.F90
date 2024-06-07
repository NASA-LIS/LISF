!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: finalize_geis
!  \label{finalize_geis}
!
! !REVISION HISTORY: 
! 15 Nov 2023; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_geis(findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use geis_forcingMod, only : geis_struc
!
! !DESCRIPTION:
!  Routine to cleanup geis forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex
  
 deallocate(geis_struc)

end subroutine finalize_geis
