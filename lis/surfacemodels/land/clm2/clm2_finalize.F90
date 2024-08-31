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
! !ROUTINE: clm2_finalize
! \label{clm2_finalize}
!
! !REVISION HISTORY:
!
! 26 Oct 2005: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine clm2_finalize()
! !USES:
  use LIS_coreMod, only : LIS_rc
  use clm2_lsmMod
!
! !DESCRIPTION:
!  
!  This routine cleans up the allocated memory structures in CLM
!  
!EOP
  implicit none

  deallocate(clm2_struc)

end subroutine clm2_finalize
