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
! !MODULE: reset_WRF_AKdom
!  \label{reset_WRF_AKdom}
!
! !REVISION HISTORY: 
!  21 Jun 2021: K.R. Arsenault; Updated for different WRF AK files
! 
! !INTERFACE:
subroutine reset_WRF_AKdom()
! !USES:
  use LIS_coreMod,          only : LIS_rc
  use WRF_AKdom_forcingMod, only : WRFAK_struc
!
! !DESCRIPTION:
!  Routine to reset WRF-AK forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n

  do n=1,LIS_rc%nnest
     WRFAK_struc(n)%WRFouttime1 = 3000.0
     WRFAK_struc(n)%WRFouttime2 = 0.0
  enddo

end subroutine reset_WRF_AKdom
