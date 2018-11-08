!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: finalize_pildas
! \label{finalize_pildas}
!
! !REVISION HISTORY: 
! 25 Apr 2013, Sujay Kumar, initial specification
! 14 Jul 2016: Mahdi Navari - Modified for PILDAS
! 
! !INTERFACE:

subroutine finalize_pildas(findex)

! !USES:
  use LIS_coreMod,     only : LIS_rc
  use pildas_forcingMod, only : pildas_struc
!
! !DESCRIPTION:
!  Routine to cleanup PILDAS forcing related memory allocations.   
! 
!EOP
  implicit none

  integer :: findex
  integer :: n

  do n=1,LIS_rc%nnest
    if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
       deallocate(pildas_struc(n)%n111)
       deallocate(pildas_struc(n)%n121)
       deallocate(pildas_struc(n)%n211)
       deallocate(pildas_struc(n)%n221)
       deallocate(pildas_struc(n)%w111)
       deallocate(pildas_struc(n)%w121)
       deallocate(pildas_struc(n)%w211)
       deallocate(pildas_struc(n)%w221)
    elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
       deallocate(pildas_struc(n)%n111)
       deallocate(pildas_struc(n)%n121)
       deallocate(pildas_struc(n)%n211)
       deallocate(pildas_struc(n)%n221)
       deallocate(pildas_struc(n)%w111)
       deallocate(pildas_struc(n)%w121)
       deallocate(pildas_struc(n)%w211)
       deallocate(pildas_struc(n)%w221)
       deallocate(pildas_struc(n)%n112)
       deallocate(pildas_struc(n)%n122)
       deallocate(pildas_struc(n)%n212)
       deallocate(pildas_struc(n)%n222)
       deallocate(pildas_struc(n)%w112)
       deallocate(pildas_struc(n)%w122)
       deallocate(pildas_struc(n)%w212)
       deallocate(pildas_struc(n)%w222)
    endif
 enddo
 deallocate(pildas_struc)
end subroutine finalize_pildas
