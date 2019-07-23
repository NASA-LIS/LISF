!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
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
  use LDT_coreMod, only : LDT_rc
  use stg4_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup Stage IV forcing related memory allocations.   
! 
!EOP
  implicit none
  integer, intent(IN) :: findex
  integer   :: n
  
  do n=1,LDT_rc%nnest
    if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 

       deallocate(stg4_struc(n)%n111)
       deallocate(stg4_struc(n)%n121)
       deallocate(stg4_struc(n)%n211)
       deallocate(stg4_struc(n)%n221)
       deallocate(stg4_struc(n)%w111)
       deallocate(stg4_struc(n)%w121)
       deallocate(stg4_struc(n)%w211)
       deallocate(stg4_struc(n)%w221)

    elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 

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
