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
! !ROUTINE: finalize_merra2
! \label{finalize_merra2}
!
! !REVISION HISTORY: 
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 12 Nov 2015: KR Arsenault, added to LDT
! 
! !INTERFACE:

subroutine finalize_merra2(findex)

! !USES:
  use LDT_coreMod,       only : LDT_rc
  use merra2_forcingMod, only : merra2_struc
!
! !DESCRIPTION:
!  Routine to cleanup MERRA2 forcing related memory allocations.   
! 
!EOP
  implicit none

  integer :: findex
  integer :: n

  do n=1,LDT_rc%nnest
    select case( LDT_rc%met_gridtransform(findex) )

     case( "bilinear" )
       deallocate(merra2_struc(n)%n111)
       deallocate(merra2_struc(n)%n121)
       deallocate(merra2_struc(n)%n211)
       deallocate(merra2_struc(n)%n221)
       deallocate(merra2_struc(n)%w111)
       deallocate(merra2_struc(n)%w121)
       deallocate(merra2_struc(n)%w211)
       deallocate(merra2_struc(n)%w221)

     case( "budget-bilinear" )
       deallocate(merra2_struc(n)%n111)
       deallocate(merra2_struc(n)%n121)
       deallocate(merra2_struc(n)%n211)
       deallocate(merra2_struc(n)%n221)
       deallocate(merra2_struc(n)%w111)
       deallocate(merra2_struc(n)%w121)
       deallocate(merra2_struc(n)%w211)
       deallocate(merra2_struc(n)%w221)
       deallocate(merra2_struc(n)%n112)
       deallocate(merra2_struc(n)%n122)
       deallocate(merra2_struc(n)%n212)
       deallocate(merra2_struc(n)%n222)
       deallocate(merra2_struc(n)%w112)
       deallocate(merra2_struc(n)%w122)
       deallocate(merra2_struc(n)%w212)
       deallocate(merra2_struc(n)%w222)

     case( "neighbor" )
       deallocate(merra2_struc(n)%n113)
    end select

 enddo
 deallocate(merra2_struc)

end subroutine finalize_merra2
