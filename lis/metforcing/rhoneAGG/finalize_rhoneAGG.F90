!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: finalize_rhoneAGG
! \label{finalize_rhoneAGG}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_rhoneAGG
! !USES:
  use LIS_coreMod, only : LIS_rc
  use rhoneAGG_forcingMod, only : rhoneAGG_struc
!
! !DESCRIPTION:
!  Routine to cleanup RHONEAGG forcing related memory allocations.   
! 
!EOP
  implicit none

 deallocate(rhoneAGG_struc)
end subroutine finalize_rhoneAGG
