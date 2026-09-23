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
! !ROUTINE: finalize_era5cds
! \label{finalize_era5cds}
!
! !REVISION HISTORY: 
! 23 Dec 2019: Sujay Kumar, initial code
! 17 Apr 2025: Hiroko Beudoing, adopted ERA5 routines for the public CDS
!                               data format
! 
! !INTERFACE:

subroutine finalize_era5cds(findex)

! !USES:
  use LDT_coreMod,     only : LDT_rc
  use era5cds_forcingMod, only : era5cds_struc
!
! !DESCRIPTION:
!  Routine to cleanup ERA5 forcing related memory allocations.   
! 
!EOP
  implicit none

  integer :: findex
  integer :: n

  do n=1,LDT_rc%nnest
    select case( LDT_rc%met_gridtransform(findex) )
     case( "bilinear" )
       deallocate(era5cds_struc(n)%n111)
       deallocate(era5cds_struc(n)%n121)
       deallocate(era5cds_struc(n)%n211)
       deallocate(era5cds_struc(n)%n221)
       deallocate(era5cds_struc(n)%w111)
       deallocate(era5cds_struc(n)%w121)
       deallocate(era5cds_struc(n)%w211)
       deallocate(era5cds_struc(n)%w221)

     case( "budget-bilinear" )
       deallocate(era5cds_struc(n)%n111)
       deallocate(era5cds_struc(n)%n121)
       deallocate(era5cds_struc(n)%n211)
       deallocate(era5cds_struc(n)%n221)
       deallocate(era5cds_struc(n)%w111)
       deallocate(era5cds_struc(n)%w121)
       deallocate(era5cds_struc(n)%w211)
       deallocate(era5cds_struc(n)%w221)
       deallocate(era5cds_struc(n)%n112)
       deallocate(era5cds_struc(n)%n122)
       deallocate(era5cds_struc(n)%n212)
       deallocate(era5cds_struc(n)%n222)
       deallocate(era5cds_struc(n)%w112)
       deallocate(era5cds_struc(n)%w122)
       deallocate(era5cds_struc(n)%w212)
       deallocate(era5cds_struc(n)%w222)

     case( "neighbor" )
       deallocate(era5cds_struc(n)%n113)

     case( "average" )
       deallocate(era5cds_struc(n)%n111)
    end select
 enddo
 deallocate(era5cds_struc)

end subroutine finalize_era5cds
