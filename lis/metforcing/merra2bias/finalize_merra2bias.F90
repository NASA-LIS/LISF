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
! !ROUTINE: finalize_merra2bias
! \label{finalize_merra2bias}
!
! !REVISION HISTORY:
! 29 Jun 2026: Kristen Whitney, initial code (based on geos-itbias)
!
! !INTERFACE:
subroutine finalize_merra2bias(findex)

! !USES:
  use LIS_coreMod,     only   : LIS_rc
  use merra2bias_forcingMod, only : merra2bias_struc
!
! !DESCRIPTION:
!  Routine to cleanup MERRA2bias forcing related memory allocations.
!
!EOP
  implicit none

  integer :: findex
  integer :: n

  do n = 1,LIS_rc%nnest
     select case( LIS_rc%met_interp(findex) )

     case( "bilinear" )
        deallocate(merra2bias_struc(n)%n111)
        deallocate(merra2bias_struc(n)%n121)
        deallocate(merra2bias_struc(n)%n211)
        deallocate(merra2bias_struc(n)%n221)
        deallocate(merra2bias_struc(n)%w111)
        deallocate(merra2bias_struc(n)%w121)
        deallocate(merra2bias_struc(n)%w211)
        deallocate(merra2bias_struc(n)%w221)

     case( "budget-bilinear" )
        deallocate(merra2bias_struc(n)%n111)
        deallocate(merra2bias_struc(n)%n121)
        deallocate(merra2bias_struc(n)%n211)
        deallocate(merra2bias_struc(n)%n221)
        deallocate(merra2bias_struc(n)%w111)
        deallocate(merra2bias_struc(n)%w121)
        deallocate(merra2bias_struc(n)%w211)
        deallocate(merra2bias_struc(n)%w221)
        deallocate(merra2bias_struc(n)%n112)
        deallocate(merra2bias_struc(n)%n122)
        deallocate(merra2bias_struc(n)%n212)
        deallocate(merra2bias_struc(n)%n222)
        deallocate(merra2bias_struc(n)%w112)
        deallocate(merra2bias_struc(n)%w122)
        deallocate(merra2bias_struc(n)%w212)
        deallocate(merra2bias_struc(n)%w222)

     case( "neighbor" )
        deallocate(merra2bias_struc(n)%n113)
     end select
  enddo
  deallocate(merra2bias_struc)

end subroutine finalize_merra2bias

