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
! !MODULE: reset_WRFout
!  \label{reset_WRFout}
!
! !REVISION HISTORY: 
! 
! !INTERFACE:
subroutine reset_WRFout()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use WRFout_forcingMod, only : WRFout_struc
!
! !DESCRIPTION:
!  Routine to reset WRFout forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n

  do n=1,LIS_rc%nnest
     WRFout_struc(n)%WRFouttime1 = 3000.0
     WRFout_struc(n)%WRFouttime2 = 0.0
  enddo

end subroutine reset_WRFout
