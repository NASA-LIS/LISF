!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: finalize_rdhm356
!  \label{finalize_rdhm356}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 17Jul2006; K. Arsenault, Added dmip II
! 19Dec2013; Shugong Wang, RDHM356 
! !INTERFACE:
subroutine finalize_rdhm356(findex)

! !USES:
  use LDT_coreMod, only : LDT_rc
  use rdhm356_forcingMod, only : rdhm356_struc_precip, &
                                 rdhm356_struc_temper, &
                                 const_wind 
!
! !DESCRIPTION:
!  Routine to cleanup dmip II forcing related memory allocations.   
! 
!EOP
  implicit none
  integer :: findex
  
  integer   :: n
  
  do n=1,LDT_rc%nnest
    if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 

       deallocate(rdhm356_struc_precip(n)%n111)
       deallocate(rdhm356_struc_precip(n)%n121)
       deallocate(rdhm356_struc_precip(n)%n211)
       deallocate(rdhm356_struc_precip(n)%n221)
       deallocate(rdhm356_struc_precip(n)%w111)
       deallocate(rdhm356_struc_precip(n)%w121)
       deallocate(rdhm356_struc_precip(n)%w211)
       deallocate(rdhm356_struc_precip(n)%w221)

       deallocate(rdhm356_struc_temper(n)%n111)
       deallocate(rdhm356_struc_temper(n)%n121)
       deallocate(rdhm356_struc_temper(n)%n211)
       deallocate(rdhm356_struc_temper(n)%n221)
       deallocate(rdhm356_struc_temper(n)%w111)
       deallocate(rdhm356_struc_temper(n)%w121)
       deallocate(rdhm356_struc_temper(n)%w211)
       deallocate(rdhm356_struc_temper(n)%w221)

    elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 

       deallocate(rdhm356_struc_precip(n)%n112)
       deallocate(rdhm356_struc_precip(n)%n122)
       deallocate(rdhm356_struc_precip(n)%n212)
       deallocate(rdhm356_struc_precip(n)%n222)
       deallocate(rdhm356_struc_precip(n)%w112)
       deallocate(rdhm356_struc_precip(n)%w122)
       deallocate(rdhm356_struc_precip(n)%w212)
       deallocate(rdhm356_struc_precip(n)%w222)

       deallocate(rdhm356_struc_temper(n)%n112)
       deallocate(rdhm356_struc_temper(n)%n122)
       deallocate(rdhm356_struc_temper(n)%n212)
       deallocate(rdhm356_struc_temper(n)%n222)
       deallocate(rdhm356_struc_temper(n)%w112)
       deallocate(rdhm356_struc_temper(n)%w122)
       deallocate(rdhm356_struc_temper(n)%w212)
       deallocate(rdhm356_struc_temper(n)%w222)

    elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then

       deallocate(rdhm356_struc_precip(n)%n113)

       deallocate(rdhm356_struc_temper(n)%n113)

    endif
 enddo
 deallocate(rdhm356_struc_precip)
 deallocate(rdhm356_struc_temper)
 deallocate(const_wind)
end subroutine finalize_rdhm356
