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
! !MODULE: reset_mrms_grib
!  \label{reset_mrms_grib}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 12 Feb 2015; Jonathan Case, revised for MRMS QPE
! 07 Sep 2017; Jessica Erlingis, revised for operational MRMS 
!
! !INTERFACE:
subroutine reset_mrms_grib()

! !USES:
  use LIS_coreMod, only : LIS_rc
  use mrms_grib_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup MRMS forcing related memory allocations.   
! 
!EOP
  implicit none
  integer   :: n

  do n=1,LIS_rc%nnest  
     mrms_grib_struc(n)%mrms_grib_time = 0.0
  enddo

end subroutine reset_mrms_grib
