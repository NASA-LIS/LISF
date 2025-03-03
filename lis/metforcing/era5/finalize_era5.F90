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
! !ROUTINE: finalize_era5
! \label{finalize_era5}
!
! !REVISION HISTORY: 
! 23 Dec 2019: Sujay Kumar, initial code
! 
! !INTERFACE:

subroutine finalize_era5(findex)

! !USES:
  use LIS_coreMod,     only : LIS_rc
  use era5_forcingMod, only : era5_struc
!
! !DESCRIPTION:
!  Routine to cleanup ERA5 forcing related memory allocations.   
! 
!EOP
  implicit none

  integer :: findex
  integer :: n

  do n=1,LIS_rc%nnest
    select case( LIS_rc%met_interp(findex) )
     case( "bilinear" )
       deallocate(era5_struc(n)%n111)
       deallocate(era5_struc(n)%n121)
       deallocate(era5_struc(n)%n211)
       deallocate(era5_struc(n)%n221)
       deallocate(era5_struc(n)%w111)
       deallocate(era5_struc(n)%w121)
       deallocate(era5_struc(n)%w211)
       deallocate(era5_struc(n)%w221)

     case( "budget-bilinear" )
       deallocate(era5_struc(n)%n111)
       deallocate(era5_struc(n)%n121)
       deallocate(era5_struc(n)%n211)
       deallocate(era5_struc(n)%n221)
       deallocate(era5_struc(n)%w111)
       deallocate(era5_struc(n)%w121)
       deallocate(era5_struc(n)%w211)
       deallocate(era5_struc(n)%w221)
       deallocate(era5_struc(n)%n112)
       deallocate(era5_struc(n)%n122)
       deallocate(era5_struc(n)%n212)
       deallocate(era5_struc(n)%n222)
       deallocate(era5_struc(n)%w112)
       deallocate(era5_struc(n)%w122)
       deallocate(era5_struc(n)%w212)
       deallocate(era5_struc(n)%w222)

     case( "neighbor" )
       deallocate(era5_struc(n)%n113)
    end select
 enddo
 deallocate(era5_struc)

end subroutine finalize_era5
