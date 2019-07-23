!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: finalize_geos
! \label{finalize_geos}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_geos(findex)
! !USES:
  use LIS_coreMod,     only : LIS_rc
  use geos_forcingMod, only : geos_struc
!
! !DESCRIPTION:
!  Routine to cleanup GEOS forcing related memory allocations.   
! 
!EOP
  implicit none
  integer :: findex
  integer :: n

  do n=1,LIS_rc%nnest
    if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

       deallocate(geos_struc(n)%n111)
       deallocate(geos_struc(n)%n121)
       deallocate(geos_struc(n)%n211)
       deallocate(geos_struc(n)%n221)
       deallocate(geos_struc(n)%w111)
       deallocate(geos_struc(n)%w121)
       deallocate(geos_struc(n)%w211)
       deallocate(geos_struc(n)%w221)
    elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

       deallocate(geos_struc(n)%n111)
       deallocate(geos_struc(n)%n121)
       deallocate(geos_struc(n)%n211)
       deallocate(geos_struc(n)%n221)
       deallocate(geos_struc(n)%w111)
       deallocate(geos_struc(n)%w121)
       deallocate(geos_struc(n)%w211)
       deallocate(geos_struc(n)%w221)

       deallocate(geos_struc(n)%n112)
       deallocate(geos_struc(n)%n122)
       deallocate(geos_struc(n)%n212)
       deallocate(geos_struc(n)%n222)
       deallocate(geos_struc(n)%w112)
       deallocate(geos_struc(n)%w122)
       deallocate(geos_struc(n)%w212)
       deallocate(geos_struc(n)%w222)
    endif
 enddo
 deallocate(geos_struc)
end subroutine finalize_geos
