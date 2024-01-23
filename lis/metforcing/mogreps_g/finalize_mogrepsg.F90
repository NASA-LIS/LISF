!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: finalize_mogrepsg
! \label{finalize_mogrepsg}
!
! !REVISION HISTORY:
! 26 Jan 2023; Yeosang Yoon, Initial Code
!
! !INTERFACE:
subroutine finalize_mogrepsg
! !USES:
  use mogrepsg_forcingMod, only : mogrepsg_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for MOGREPS-G forcing.
!
!EOP
  implicit none

  deallocate(mogrepsg_struc)

end subroutine finalize_mogrepsg
