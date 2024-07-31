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
! !MODULE: finalize_WRF_AKdom
!  \label{finalize_WRF_AKdom}
!
! !REVISION HISTORY: 
!  21 Jun 2021: K.R. Arsenault; Updated for different WRF output files
! 
! !INTERFACE:
subroutine finalize_WRF_AKdom(findex)
! !USES:
  use LIS_coreMod,          only : LIS_rc
  use WRF_AKdom_forcingMod, only : WRFAK_struc
!
! !DESCRIPTION:
!  Routine to cleanup WRF Alaska (AK) forcing related memory allocations.   
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing
!  \end{description}
!EOP
  implicit none
  
  integer   :: findex

  integer   :: n
  
  deallocate(WRFAK_struc)

end subroutine finalize_WRF_AKdom
