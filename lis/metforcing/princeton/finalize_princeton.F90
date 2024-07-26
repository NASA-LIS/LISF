!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <LIS_misc.h>
!BOP
! !ROUTINE: finalize_princeton
! \label{finalize_princeton}
! 
! !REVISION HISTORY: 
! 26 Jan 2007; Hiroko Kato, Initial Code
! 
! !INTERFACE:
subroutine finalize_princeton(findex)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use princeton_forcingMod, only : princeton_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for princeton forcing. 
!
!EOP
  implicit none
  integer :: findex
  integer :: n 
  
  do n=1,LIS_rc%nnest
     if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
        deallocate(princeton_struc(n)%n111)
        deallocate(princeton_struc(n)%n121)
        deallocate(princeton_struc(n)%n211)
        deallocate(princeton_struc(n)%n221)
        deallocate(princeton_struc(n)%w111)
        deallocate(princeton_struc(n)%w121)
        deallocate(princeton_struc(n)%w211)
        deallocate(princeton_struc(n)%w221)
     elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
        deallocate(princeton_struc(n)%n111)
        deallocate(princeton_struc(n)%n121)
        deallocate(princeton_struc(n)%n211)
        deallocate(princeton_struc(n)%n221)
        deallocate(princeton_struc(n)%w111)
        deallocate(princeton_struc(n)%w121)
        deallocate(princeton_struc(n)%w211)
        deallocate(princeton_struc(n)%w221)

        deallocate(princeton_struc(n)%n112)
        deallocate(princeton_struc(n)%n122)
        deallocate(princeton_struc(n)%n212)
        deallocate(princeton_struc(n)%n222)
        deallocate(princeton_struc(n)%w112)
        deallocate(princeton_struc(n)%w122)
        deallocate(princeton_struc(n)%w212)
        deallocate(princeton_struc(n)%w222)
     endif
  enddo
  deallocate(princeton_struc)
end subroutine finalize_princeton
