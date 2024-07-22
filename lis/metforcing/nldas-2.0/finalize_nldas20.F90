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
! !ROUTINE: finalize_nldas20
!  \label{finalize_nldas20}
!
! !REVISION HISTORY:
! 11 Jul 2024: David Mocko, Initial Specification
!                           (derived from finalize_nldas2.F90)
!
! !INTERFACE:
subroutine finalize_nldas20(findex)
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use nldas20_forcingMod, only : nldas20_struc

  implicit none
! !ARGUMENTS:
  integer   :: findex
!
! !DESCRIPTION:
!  Routine to cleanup nldas20 forcing related memory allocations.
!
!EOP
  integer   :: n

  do n = 1,LIS_rc%nnest
     if (trim(LIS_rc%met_interp(findex)).eq."bilinear") then
        deallocate(nldas20_struc(n)%n111)
        deallocate(nldas20_struc(n)%n121)
        deallocate(nldas20_struc(n)%n211)
        deallocate(nldas20_struc(n)%n221)
        deallocate(nldas20_struc(n)%w111)
        deallocate(nldas20_struc(n)%w121)
        deallocate(nldas20_struc(n)%w211)
        deallocate(nldas20_struc(n)%w221)
     elseif (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") &
          then
        deallocate(nldas20_struc(n)%n111)
        deallocate(nldas20_struc(n)%n121)
        deallocate(nldas20_struc(n)%n211)
        deallocate(nldas20_struc(n)%n221)
        deallocate(nldas20_struc(n)%w111)
        deallocate(nldas20_struc(n)%w121)
        deallocate(nldas20_struc(n)%w211)
        deallocate(nldas20_struc(n)%w221)
        deallocate(nldas20_struc(n)%n112)
        deallocate(nldas20_struc(n)%n122)
        deallocate(nldas20_struc(n)%n212)
        deallocate(nldas20_struc(n)%n222)
        deallocate(nldas20_struc(n)%w112)
        deallocate(nldas20_struc(n)%w122)
        deallocate(nldas20_struc(n)%w212)
        deallocate(nldas20_struc(n)%w222)
     elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
        deallocate(nldas20_struc(n)%n113)
     endif

     if (LIS_rc%met_ecor(findex).ne."none") then
        deallocate(nldas20_struc(n)%orig_ediff)
     endif
  enddo
  deallocate(nldas20_struc)

end subroutine finalize_nldas20

