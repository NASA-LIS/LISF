!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: set_ZOBLER_GFS_texture_attribs
!  \label{set_ZOBLER_GFS_texture_attribs}
!
! !REVISION HISTORY:
!  19 Aug 2014: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine set_ZOBLER_GFS_texture_attribs()

! !USES:
  use LDT_paramDataMod

  implicit none

! !ARGUMENTS: 

! !DESCRIPTION:
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[fgrd]
!     fraction of grid covered by each vegetation type
!   \item[maskarray]
!     landmask for the region of interest
!   \end{description}
!EOP      
!
  LDT_LSMparam_struc(:)%texture%num_bins = 9
  LDT_LSMparam_struc(:)%texture%vlevels = 1

end subroutine set_ZOBLER_GFS_texture_attribs
