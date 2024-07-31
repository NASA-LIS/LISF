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
! !MODULE: finalize_WRFoutv2
!  \label{finalize_WRFoutv2}
!
! !REVISION HISTORY: 
! 20 Nov 2020; K.R. Arsenault, Updated for different WRF output files
! 
! !INTERFACE:
subroutine finalize_WRFoutv2(findex)
! !USES:
  use LIS_coreMod,         only : LIS_rc
  use WRFoutv2_forcingMod, only : WRFoutv2_struc
!
! !DESCRIPTION:
!  Routine to cleanup WRFoutv2 forcing related memory allocations.   
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing
!  \end{description}
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex
  
  deallocate(WRFoutv2_struc)

end subroutine finalize_WRFoutv2
