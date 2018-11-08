!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
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
  use LIS_coreMod, only : LIS_rc
  use ecmwfreanal_forcingMod, only : ecmwfreanal_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for ecmwf reanalysis forcing. 
!
!EOP
  implicit none
  integer :: findex
  integer :: n 
  
  do n=1,LIS_rc%nnest
     if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

        deallocate(ecmwfreanal_struc(n)%n111)
        deallocate(ecmwfreanal_struc(n)%n121)
        deallocate(ecmwfreanal_struc(n)%n211)
        deallocate(ecmwfreanal_struc(n)%n221)
        deallocate(ecmwfreanal_struc(n)%w111)
        deallocate(ecmwfreanal_struc(n)%w121)
        deallocate(ecmwfreanal_struc(n)%w211)
        deallocate(ecmwfreanal_struc(n)%w221)
     elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

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
