!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: finalize_nldas3sw
!  \label{finalize_nldas3sw}
!
! !REVISION HISTORY:
! 27 Dec 2024: David Mocko, Initial Specification
!                           (derived from finalize_nldas20.F90)
!
! !INTERFACE:
subroutine finalize_nldas3sw(findex)
! !USES:
  use LIS_coreMod,         only : LIS_rc
  use nldas3sw_forcingMod, only : nldas3sw_struc

  implicit none
! !ARGUMENTS:
  integer   :: findex
!
! !DESCRIPTION:
!  Routine to cleanup NLDAS-3 SWdown forcing related memory allocations.
!
!EOP
  integer   :: n

  do n = 1,LIS_rc%nnest
     if (trim(LIS_rc%met_interp(findex)).eq."bilinear") then
        deallocate(nldas3sw_struc(n)%n111)
        deallocate(nldas3sw_struc(n)%n121)
        deallocate(nldas3sw_struc(n)%n211)
        deallocate(nldas3sw_struc(n)%n221)
        deallocate(nldas3sw_struc(n)%w111)
        deallocate(nldas3sw_struc(n)%w121)
        deallocate(nldas3sw_struc(n)%w211)
        deallocate(nldas3sw_struc(n)%w221)
     elseif (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear")     &
          then
        deallocate(nldas3sw_struc(n)%n111)
        deallocate(nldas3sw_struc(n)%n121)
        deallocate(nldas3sw_struc(n)%n211)
        deallocate(nldas3sw_struc(n)%n221)
        deallocate(nldas3sw_struc(n)%w111)
        deallocate(nldas3sw_struc(n)%w121)
        deallocate(nldas3sw_struc(n)%w211)
        deallocate(nldas3sw_struc(n)%w221)
        deallocate(nldas3sw_struc(n)%n112)
        deallocate(nldas3sw_struc(n)%n122)
        deallocate(nldas3sw_struc(n)%n212)
        deallocate(nldas3sw_struc(n)%n222)
        deallocate(nldas3sw_struc(n)%w112)
        deallocate(nldas3sw_struc(n)%w122)
        deallocate(nldas3sw_struc(n)%w212)
        deallocate(nldas3sw_struc(n)%w222)
     elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
        deallocate(nldas3sw_struc(n)%n113)
     endif

  enddo
  deallocate(nldas3sw_struc)

end subroutine finalize_nldas3sw

