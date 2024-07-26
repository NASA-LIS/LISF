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
! !MODULE: finalize_gswp1
! \label{finalize_gswp1}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_gswp1
! !USES:
  use gswp1_forcingMod, only : gswp1_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GSWP1 forcing. 
!
!EOP  
  implicit none
  
  deallocate(gswp1_struc)

end subroutine finalize_gswp1
