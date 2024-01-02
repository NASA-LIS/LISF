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
! !MODULE: finalize_galwemge
! \label{finalize_galwemge}
!
! !REVISION HISTORY:
! 09 May 2022; Yeosang Yoon, Initial Code
!
! !INTERFACE:
subroutine finalize_galwemge
! !USES:
  use galwemge_forcingMod, only : galwemge_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GALWEM-GE forcing.
!
!EOP
  implicit none

  deallocate(galwemge_struc)

end subroutine finalize_galwemge
