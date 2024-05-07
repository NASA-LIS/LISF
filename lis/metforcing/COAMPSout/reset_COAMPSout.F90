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
! !MODULE: reset_COAMPSout
!  \label{reset_COAMPSout}
!
! !REVISION HISTORY: 
! 
! !INTERFACE:
subroutine reset_COAMPSout()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use COAMPSout_forcingMod, only : COAMPSout_struc
!
! !DESCRIPTION:
!  Routine to reset COAMPSout forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n

  do n=1,LIS_rc%nnest
     COAMPSout_struc(n)%COAMPSouttime1 = 3000.0
     COAMPSout_struc(n)%COAMPSouttime2 = 0.0
  enddo

end subroutine reset_COAMPSout
