!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: finalize_merraland
! \label{finalize_merraland}
!
! !REVISION HISTORY: 
! 12 Oct 2009, Eric Kemp
! 22 Jan 2015: KR Arsenault, added to LDT
! 
! !INTERFACE:

subroutine finalize_merraland(findex)

! !USES:
  use LDT_coreMod,     only : LDT_rc
  use merraland_forcingMod, only : merraland_struc
!
! !DESCRIPTION:
!  Routine to cleanup MERRA-Land forcing related memory allocations.   
! 
!EOP
  implicit none

  integer :: findex
  integer :: n

  do n=1,LDT_rc%nnest
   
    select case( LDT_rc%met_gridtransform(findex) )

     case( "bilinear" )
       deallocate(merraland_struc(n)%n111)
       deallocate(merraland_struc(n)%n121)
       deallocate(merraland_struc(n)%n211)
       deallocate(merraland_struc(n)%n221)
       deallocate(merraland_struc(n)%w111)
       deallocate(merraland_struc(n)%w121)
       deallocate(merraland_struc(n)%w211)
       deallocate(merraland_struc(n)%w221)

     case( "budget-bilinear" )

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

     case( "neighbor" )
       deallocate(merraland_struc(n)%n113)

    end select
 enddo
 deallocate(merraland_struc)

end subroutine finalize_merraland
