!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: finalize_ecmwfreanal
! \label{finalize_ecmwfreanal}
! 
! !REVISION HISTORY: 
! 25 Oct 2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_ecmwfreanal(findex)
! !USES:
  use LDT_coreMod, only : LDT_rc
  use ecmwfreanal_forcingMod, only : ecmwfreanal_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for ecmwf reanalysis forcing. 
!
!EOP
  implicit none
  integer :: findex
  integer :: n 
  
  do n=1,LDT_rc%nnest
     if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 

        deallocate(ecmwfreanal_struc(n)%n111)
        deallocate(ecmwfreanal_struc(n)%n121)
        deallocate(ecmwfreanal_struc(n)%n211)
        deallocate(ecmwfreanal_struc(n)%n221)
        deallocate(ecmwfreanal_struc(n)%w111)
        deallocate(ecmwfreanal_struc(n)%w121)
        deallocate(ecmwfreanal_struc(n)%w211)
        deallocate(ecmwfreanal_struc(n)%w221)
     elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 

        deallocate(ecmwfreanal_struc(n)%n111)
        deallocate(ecmwfreanal_struc(n)%n121)
        deallocate(ecmwfreanal_struc(n)%n211)
        deallocate(ecmwfreanal_struc(n)%n221)
        deallocate(ecmwfreanal_struc(n)%w111)
        deallocate(ecmwfreanal_struc(n)%w121)
        deallocate(ecmwfreanal_struc(n)%w211)
        deallocate(ecmwfreanal_struc(n)%w221)

        deallocate(ecmwfreanal_struc(n)%n112)
        deallocate(ecmwfreanal_struc(n)%n122)
        deallocate(ecmwfreanal_struc(n)%n212)
        deallocate(ecmwfreanal_struc(n)%n222)
        deallocate(ecmwfreanal_struc(n)%w112)
        deallocate(ecmwfreanal_struc(n)%w122)
        deallocate(ecmwfreanal_struc(n)%w212)
        deallocate(ecmwfreanal_struc(n)%w222)
     endif
  enddo
  deallocate(ecmwfreanal_struc)
end subroutine finalize_ecmwfreanal
