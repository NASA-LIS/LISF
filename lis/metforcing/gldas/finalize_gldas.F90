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
! !MODULE: finalize_gldas
! \label{finalize_gldas}
! 
! !REVISION HISTORY: 
!  19 Sept 2008: Sujay Kumar: Initial Implementation
! 
! !INTERFACE:
subroutine finalize_gldas(findex)
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use gldas_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GDAS forcing. 
!
!EOP  
  implicit none

  integer :: n 
  integer :: findex
  
  do n=1,LIS_rc%nnest  
     if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
        deallocate(gldas_struc(n)%n111)
        deallocate(gldas_struc(n)%n121)
        deallocate(gldas_struc(n)%n211)
        deallocate(gldas_struc(n)%n221)
        deallocate(gldas_struc(n)%w111)
        deallocate(gldas_struc(n)%w121)
        deallocate(gldas_struc(n)%w211)
        deallocate(gldas_struc(n)%w221)
     elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
        deallocate(gldas_struc(n)%n111)
        deallocate(gldas_struc(n)%n121)
        deallocate(gldas_struc(n)%n211)
        deallocate(gldas_struc(n)%n221)
        deallocate(gldas_struc(n)%w111)
        deallocate(gldas_struc(n)%w121)
        deallocate(gldas_struc(n)%w211)
        deallocate(gldas_struc(n)%w221)
        deallocate(gldas_struc(n)%n112)
        deallocate(gldas_struc(n)%n122)
        deallocate(gldas_struc(n)%n212)
        deallocate(gldas_struc(n)%n222)
        deallocate(gldas_struc(n)%w112)
        deallocate(gldas_struc(n)%w122)
        deallocate(gldas_struc(n)%w212)
        deallocate(gldas_struc(n)%w222)
     endif
  enddo
  deallocate(gldas_struc)
end subroutine finalize_gldas
