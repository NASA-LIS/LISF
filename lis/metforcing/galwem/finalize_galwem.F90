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
! !MODULE: finalize_galwem
! \label{finalize_galwem}
!
! !REVISION HISTORY:
! 04 Apr 2022; Yeosang Yoon, Initial Code
!
! !INTERFACE:
subroutine finalize_galwem
! !USES:
  use galwem_forcingMod, only : galwem_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GALWEM forcing.
!
!EOP
  implicit none

  deallocate(galwem_struc)

end subroutine finalize_galwem
