!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: set_CONSTANT_texture_attribs
!  \label{set_CONSTANT_texture_attribs}
!
! !REVISION HISTORY:
!  19 Aug 2014: Sujay Kumar; Initial Specification
!   7 Oct 2016: KR Arsenault; Added back CONSTANT soil texture option
!
! !INTERFACE:
subroutine set_CONSTANT_texture_attribs()

! !USES:
  use LDT_paramDataMod

  implicit none

! !ARGUMENTS: 

! !DESCRIPTION:
!
!  The arguments are:
!  \begin{description}
!   \item[num_bins]
!    Number of categories
!   \item[vlevels]
!    Number of vertical levels
!   \end{description}
!EOP      
!
  LDT_LSMparam_struc(:)%texture%num_bins = 16
  LDT_LSMparam_struc(:)%texture%vlevels = 1

end subroutine set_CONSTANT_texture_attribs
