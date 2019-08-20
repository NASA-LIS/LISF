!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: finalize_nldas1
!  \label{finalize_nldas1}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_nldas1(findex)
! !USES:
  use LIS_coreMod,      only : LIS_rc
  use nldas1_forcingMod, only : nldas1_struc
!
! !DESCRIPTION:
!  Routine to cleanup nldas1 forcing related memory allocations.   
! 
!EOP
  implicit none
  integer, intent(IN) :: findex
  
  integer   :: n
  
  do n=1,LIS_rc%nnest
    if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

       deallocate(nldas1_struc(n)%n111)
       deallocate(nldas1_struc(n)%n121)
       deallocate(nldas1_struc(n)%n211)
       deallocate(nldas1_struc(n)%n221)
       deallocate(nldas1_struc(n)%w111)
       deallocate(nldas1_struc(n)%w121)
       deallocate(nldas1_struc(n)%w211)
       deallocate(nldas1_struc(n)%w221)
    elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

       deallocate(nldas1_struc(n)%n111)
       deallocate(nldas1_struc(n)%n121)
       deallocate(nldas1_struc(n)%n211)
       deallocate(nldas1_struc(n)%n221)
       deallocate(nldas1_struc(n)%w111)
       deallocate(nldas1_struc(n)%w121)
       deallocate(nldas1_struc(n)%w211)
       deallocate(nldas1_struc(n)%w221)

       deallocate(nldas1_struc(n)%n112)
       deallocate(nldas1_struc(n)%n122)
       deallocate(nldas1_struc(n)%n212)
       deallocate(nldas1_struc(n)%n222)
       deallocate(nldas1_struc(n)%w112)
       deallocate(nldas1_struc(n)%w122)
       deallocate(nldas1_struc(n)%w212)
       deallocate(nldas1_struc(n)%w222)
    elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
       deallocate(nldas1_struc(n)%n113)
    endif
 enddo
 deallocate(nldas1_struc)

end subroutine finalize_nldas1
