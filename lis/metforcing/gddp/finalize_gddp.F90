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
! !MODULE: finalize_gddp
! \label{finalize_gddp}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_gddp
! !USES:
  use gddp_forcingMod, only : gddp_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GDDP forcing. 
!
!EOP  
  implicit none
  
  deallocate(gddp_struc)

end subroutine finalize_gddp
