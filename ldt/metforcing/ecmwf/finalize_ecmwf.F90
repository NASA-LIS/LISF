!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !ROUTINE: finalize_ecmwf
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_ecmwf(findex)
! !USES:
  use LDT_coreMod, only : LDT_rc
  use ecmwf_forcingMod, only : ecmwf_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for ECMWF forcing. 
!
!EOP
  implicit none
  integer :: n 
  integer :: findex

  do n=1,LDT_rc%nnest
    select case( LDT_rc%met_gridtransform(findex) )

     case( "bilinear" ) 
        deallocate(ecmwf_struc(n)%n111)
        deallocate(ecmwf_struc(n)%n121)
        deallocate(ecmwf_struc(n)%n211)
        deallocate(ecmwf_struc(n)%n221)
        deallocate(ecmwf_struc(n)%w111)
        deallocate(ecmwf_struc(n)%w121)
        deallocate(ecmwf_struc(n)%w211)
        deallocate(ecmwf_struc(n)%w221)

     case( "budget-bilinear" ) 
        deallocate(ecmwf_struc(n)%n111)
        deallocate(ecmwf_struc(n)%n121)
        deallocate(ecmwf_struc(n)%n211)
        deallocate(ecmwf_struc(n)%n221)
        deallocate(ecmwf_struc(n)%w111)
        deallocate(ecmwf_struc(n)%w121)
        deallocate(ecmwf_struc(n)%w211)
        deallocate(ecmwf_struc(n)%w221)
        deallocate(ecmwf_struc(n)%n112)
        deallocate(ecmwf_struc(n)%n122)
        deallocate(ecmwf_struc(n)%n212)
        deallocate(ecmwf_struc(n)%n222)
        deallocate(ecmwf_struc(n)%w112)
        deallocate(ecmwf_struc(n)%w122)
        deallocate(ecmwf_struc(n)%w212)
        deallocate(ecmwf_struc(n)%w222)
     end select

  enddo
  deallocate(ecmwf_struc)

end subroutine finalize_ecmwf
