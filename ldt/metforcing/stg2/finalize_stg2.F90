!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: finalize_stg2
!  \label{finalize_stg2}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 17Jul2006; K. Arsenault, Added Stage II
! 
! !INTERFACE:
subroutine finalize_stg2(findex)

! !USES:
  use LDT_coreMod, only : LDT_rc
  use stg2_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup Stage II forcing related memory allocations.   
! 
!EOP
  implicit none
  integer, intent(IN) :: findex
  
  integer   :: n
  
  do n=1,LDT_rc%nnest

    select case( LDT_rc%met_gridtransform(findex) )

     case( "bilinear" )
       deallocate(stg2_struc(n)%n111)
       deallocate(stg2_struc(n)%n121)
       deallocate(stg2_struc(n)%n211)
       deallocate(stg2_struc(n)%n221)
       deallocate(stg2_struc(n)%w111)
       deallocate(stg2_struc(n)%w121)
       deallocate(stg2_struc(n)%w211)
       deallocate(stg2_struc(n)%w221)

     case( "budget-bilinear" )

       deallocate(stg2_struc(n)%n112)
       deallocate(stg2_struc(n)%n122)
       deallocate(stg2_struc(n)%n212)
       deallocate(stg2_struc(n)%n222)
       deallocate(stg2_struc(n)%w112)
       deallocate(stg2_struc(n)%w122)
       deallocate(stg2_struc(n)%w212)
       deallocate(stg2_struc(n)%w222)
   end select
 enddo
 deallocate(stg2_struc)

end subroutine finalize_stg2
