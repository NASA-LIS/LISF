!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: finalize_merraland
! \label{finalize_merraland}
!
! !REVISION HISTORY: 
! 12 Oct 2009, Eric Kemp
! 
! !INTERFACE:

subroutine finalize_merraland(findex)

! !USES:
  use LIS_coreMod,     only : LIS_rc
  use merraland_forcingMod, only : merraland_struc
!
! !DESCRIPTION:
!  Routine to cleanup MERRA-Land forcing related memory allocations.   
! 
!EOP
  implicit none

  integer :: findex
  integer :: n

  do n=1,LIS_rc%nnest
    if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

       deallocate(merraland_struc(n)%n111)
       deallocate(merraland_struc(n)%n121)
       deallocate(merraland_struc(n)%n211)
       deallocate(merraland_struc(n)%n221)
       deallocate(merraland_struc(n)%w111)
       deallocate(merraland_struc(n)%w121)
       deallocate(merraland_struc(n)%w211)
       deallocate(merraland_struc(n)%w221)
    elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

       deallocate(merraland_struc(n)%n111)
       deallocate(merraland_struc(n)%n121)
       deallocate(merraland_struc(n)%n211)
       deallocate(merraland_struc(n)%n221)
       deallocate(merraland_struc(n)%w111)
       deallocate(merraland_struc(n)%w121)
       deallocate(merraland_struc(n)%w211)
       deallocate(merraland_struc(n)%w221)

       deallocate(merraland_struc(n)%n112)
       deallocate(merraland_struc(n)%n122)
       deallocate(merraland_struc(n)%n212)
       deallocate(merraland_struc(n)%n222)
       deallocate(merraland_struc(n)%w112)
       deallocate(merraland_struc(n)%w122)
       deallocate(merraland_struc(n)%w212)
       deallocate(merraland_struc(n)%w222)
    elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 

       deallocate(merraland_struc(n)%n113)
    endif
 enddo
 deallocate(merraland_struc)
end subroutine finalize_merraland
