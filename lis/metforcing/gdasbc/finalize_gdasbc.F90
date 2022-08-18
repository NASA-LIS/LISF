!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: finalize_gdasbc
!  \label{finalize_gdasbc}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 24 Aug 2007: Chuck Alonge; Modified for use with GDASBC data
! 
! !INTERFACE:
subroutine finalize_gdasbc(findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use gdasbc_forcingMod, only : gdasbc_struc
!
! !DESCRIPTION:
!  Routine to cleanup gdasbc forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex
  
  do n=1,LIS_rc%nnest
    if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

       deallocate(gdasbc_struc(n)%n111)
       deallocate(gdasbc_struc(n)%n121)
       deallocate(gdasbc_struc(n)%n211)
       deallocate(gdasbc_struc(n)%n221)
       deallocate(gdasbc_struc(n)%w111)
       deallocate(gdasbc_struc(n)%w121)
       deallocate(gdasbc_struc(n)%w211)
       deallocate(gdasbc_struc(n)%w221)
    elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

       deallocate(gdasbc_struc(n)%n111)
       deallocate(gdasbc_struc(n)%n121)
       deallocate(gdasbc_struc(n)%n211)
       deallocate(gdasbc_struc(n)%n221)
       deallocate(gdasbc_struc(n)%w111)
       deallocate(gdasbc_struc(n)%w121)
       deallocate(gdasbc_struc(n)%w211)
       deallocate(gdasbc_struc(n)%w221)

       deallocate(gdasbc_struc(n)%n112)
       deallocate(gdasbc_struc(n)%n122)
       deallocate(gdasbc_struc(n)%n212)
       deallocate(gdasbc_struc(n)%n222)
       deallocate(gdasbc_struc(n)%w112)
       deallocate(gdasbc_struc(n)%w122)
       deallocate(gdasbc_struc(n)%w212)
       deallocate(gdasbc_struc(n)%w222)
    elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then

       deallocate(gdasbc_struc(n)%n113)
    endif

 enddo
 deallocate(gdasbc_struc)

end subroutine finalize_gdasbc
