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
! !MODULE: finalize_stg4
!  \label{finalize_stg4}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 17Jul2006; K. Arsenault, Added Stage IV
! 
! !INTERFACE:
subroutine finalize_stg4(findex)

! !USES:
  use LIS_coreMod, only : LIS_rc
  use stg4_forcingMod
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

       deallocate(stg4_struc(n)%n111)
       deallocate(stg4_struc(n)%n121)
       deallocate(stg4_struc(n)%n211)
       deallocate(stg4_struc(n)%n221)
       deallocate(stg4_struc(n)%w111)
       deallocate(stg4_struc(n)%w121)
       deallocate(stg4_struc(n)%w211)
       deallocate(stg4_struc(n)%w221)

    elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

       deallocate(stg4_struc(n)%n112)
       deallocate(stg4_struc(n)%n122)
       deallocate(stg4_struc(n)%n212)
       deallocate(stg4_struc(n)%n222)
       deallocate(stg4_struc(n)%w112)
       deallocate(stg4_struc(n)%w122)
       deallocate(stg4_struc(n)%w212)
       deallocate(stg4_struc(n)%w222)
    endif
 enddo
 deallocate(stg4_struc)

end subroutine finalize_stg4
