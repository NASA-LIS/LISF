!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: set_AVHRR_lai_attribs
!  \label{set_AVHRR_lai_attribs}
!
! !REVISION HISTORY:
!  19 Aug 2014: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine set_AVHRR_lai_attribs()

! !USES:
  use LDT_laisaiMod

  implicit none

! !ARGUMENTS: 

! !DESCRIPTION:
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[num_times]
!     number of climatology time points (e.g., months)
!   \end{description}
!EOP      
!
  LDT_laisai_struc(:)%lai%num_bins = 1
  LDT_laisai_struc(:)%lai%num_times = 12

end subroutine set_AVHRR_lai_attribs

