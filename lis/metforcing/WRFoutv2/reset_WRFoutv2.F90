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
! !MODULE: reset_WRFoutv2
!  \label{reset_WRFoutv2}
!
! !REVISION HISTORY: 
! 20 Nov 2020; K.R. Arsenault, Updated for different WRF output files
! 
! !INTERFACE:
subroutine reset_WRFoutv2()
! !USES:
  use LIS_coreMod,         only : LIS_rc
  use WRFoutv2_forcingMod, only : WRFoutv2_struc
!
! !DESCRIPTION:
!  Routine to reset WRFoutv2 forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n

  do n=1,LIS_rc%nnest
     WRFoutv2_struc(n)%WRFouttime1 = 3000.0
     WRFoutv2_struc(n)%WRFouttime2 = 0.0
  enddo

end subroutine reset_WRFoutv2
