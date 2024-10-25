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
! !MODULE: finalize_nldas2
!  \label{finalize_nldas2}
!
! !REVISION HISTORY: 
! 25 Oct 2005; Sujay Kumar, Initial Code
! 24 Aug 2007: Chuck Alonge; Modified for use with NLDAS2 data
! 
! !INTERFACE:
subroutine finalize_nldas2(findex)

! !USES:
  use LDT_coreMod,       only : LDT_rc
  use nldas2_forcingMod, only : nldas2_struc
!
! !DESCRIPTION:
!  Routine to cleanup nldas2 forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex
  
  do n=1,LDT_rc%nnest
    if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 

       deallocate(nldas2_struc(n)%n111)
       deallocate(nldas2_struc(n)%n121)
       deallocate(nldas2_struc(n)%n211)
       deallocate(nldas2_struc(n)%n221)
       deallocate(nldas2_struc(n)%w111)
       deallocate(nldas2_struc(n)%w121)
       deallocate(nldas2_struc(n)%w211)
       deallocate(nldas2_struc(n)%w221)
    elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 

       deallocate(nldas2_struc(n)%n111)
       deallocate(nldas2_struc(n)%n121)
       deallocate(nldas2_struc(n)%n211)
       deallocate(nldas2_struc(n)%n221)
       deallocate(nldas2_struc(n)%w111)
       deallocate(nldas2_struc(n)%w121)
       deallocate(nldas2_struc(n)%w211)
       deallocate(nldas2_struc(n)%w221)

       deallocate(nldas2_struc(n)%n112)
       deallocate(nldas2_struc(n)%n122)
       deallocate(nldas2_struc(n)%n212)
       deallocate(nldas2_struc(n)%n222)
       deallocate(nldas2_struc(n)%w112)
       deallocate(nldas2_struc(n)%w122)
       deallocate(nldas2_struc(n)%w212)
       deallocate(nldas2_struc(n)%w222)
    elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then

       deallocate(nldas2_struc(n)%n113)
    endif

#if 0
    if(LDT_rc%met_ecor(findex).ne."none") then 
       deallocate(nldas2_struc(n)%orig_ediff)
    endif
#endif
 enddo
 deallocate(nldas2_struc)

end subroutine finalize_nldas2
