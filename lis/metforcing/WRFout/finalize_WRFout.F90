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
! !MODULE: finalize_WRFout
!  \label{finalize_WRFout}
!
! !REVISION HISTORY: 
! 
! !INTERFACE:
subroutine finalize_WRFout(findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use WRFout_forcingMod, only : WRFout_struc
!
! !DESCRIPTION:
!  Routine to cleanup WRFout forcing related memory allocations.   
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
  
  deallocate(WRFout_struc)

end subroutine finalize_WRFout
