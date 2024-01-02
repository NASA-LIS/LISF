!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: finalize_AWRAL
!  \label{finalize_AWRAL}
!
! !REVISION HISTORY: 
! 30 Jan 2017; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_AWRAL(findex)

! !USES:
  use LIS_coreMod, only : LIS_rc
  use AWRAL_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup Stage IV forcing related memory allocations.   
! 
!EOP
  implicit none
  integer, intent(IN) :: findex
  integer   :: n
  
  do n=1,LIS_rc%nnest
    if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

       deallocate(AWRAL_struc(n)%n111)
       deallocate(AWRAL_struc(n)%n121)
       deallocate(AWRAL_struc(n)%n211)
       deallocate(AWRAL_struc(n)%n221)
       deallocate(AWRAL_struc(n)%w111)
       deallocate(AWRAL_struc(n)%w121)
       deallocate(AWRAL_struc(n)%w211)
       deallocate(AWRAL_struc(n)%w221)

    elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

       deallocate(AWRAL_struc(n)%n112)
       deallocate(AWRAL_struc(n)%n122)
       deallocate(AWRAL_struc(n)%n212)
       deallocate(AWRAL_struc(n)%n222)
       deallocate(AWRAL_struc(n)%w112)
       deallocate(AWRAL_struc(n)%w122)
       deallocate(AWRAL_struc(n)%w212)
       deallocate(AWRAL_struc(n)%w222)
    endif
 enddo
 deallocate(AWRAL_struc)

end subroutine finalize_AWRAL
