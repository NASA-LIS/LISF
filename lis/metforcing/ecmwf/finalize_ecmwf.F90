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
! !ROUTINE: finalize_ecmwf
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_ecmwf(findex)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use ecmwf_forcingMod, only : ecmwf_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for ECMWF forcing. 
!
!EOP
  implicit none
  integer :: n 
  integer :: findex

  do n=1,LIS_rc%nnest
     if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
        deallocate(ecmwf_struc(n)%n111)
        deallocate(ecmwf_struc(n)%n121)
        deallocate(ecmwf_struc(n)%n211)
        deallocate(ecmwf_struc(n)%n221)
        deallocate(ecmwf_struc(n)%w111)
        deallocate(ecmwf_struc(n)%w121)
        deallocate(ecmwf_struc(n)%w211)
        deallocate(ecmwf_struc(n)%w221)
     elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
        deallocate(ecmwf_struc(n)%n111)
        deallocate(ecmwf_struc(n)%n121)
        deallocate(ecmwf_struc(n)%n211)
        deallocate(ecmwf_struc(n)%n221)
        deallocate(ecmwf_struc(n)%w111)
        deallocate(ecmwf_struc(n)%w121)
        deallocate(ecmwf_struc(n)%w211)
        deallocate(ecmwf_struc(n)%w221)
        deallocate(ecmwf_struc(n)%n112)
        deallocate(ecmwf_struc(n)%n122)
        deallocate(ecmwf_struc(n)%n212)
        deallocate(ecmwf_struc(n)%n222)
        deallocate(ecmwf_struc(n)%w112)
        deallocate(ecmwf_struc(n)%w122)
        deallocate(ecmwf_struc(n)%w212)
        deallocate(ecmwf_struc(n)%w222)
     endif
  enddo
  deallocate(ecmwf_struc)

end subroutine finalize_ecmwf
