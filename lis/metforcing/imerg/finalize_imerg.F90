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
! !ROUTINE: finalize_imerg
! \label{finalize_imerg}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 09Mar2015: Jon Case;  Added IMERG precipitation reader
! 
! !INTERFACE:
subroutine finalize_imerg(findex)

! !USES:
  use imerg_forcingMod, only : imerg_struc
! !DESCRIPTION: 
!  Routine to cleanup IMERG forcing related memory allocations.   
!
!EOP
  implicit none
  integer, intent(IN) :: findex  

  deallocate(imerg_struc)

end subroutine finalize_imerg
