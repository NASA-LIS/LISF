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
! !ROUTINE: mos_finalize
! \label{mos_finalize}
!
! !REVISION HISTORY:
!
! 04 Oct 2005: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine mos_finalize()
! !USES:
  use LIS_coreMod, only : LIS_rc
  use mos_lsmMod

!
! !DESCRIPTION:
!  
!  This routine cleans up the allocated memory structures in Mosaic
!  
!EOP
  implicit none

  deallocate(mos_struc)

end subroutine mos_finalize
