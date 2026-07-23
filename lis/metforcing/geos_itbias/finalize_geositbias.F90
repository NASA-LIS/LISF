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
! !ROUTINE: finalize_geositbias
! \label{finalize_geositbias}
!
! !REVISION HISTORY:
! 02 Oct 2025: Fadji Maina, initial code (based on geos-it)
! 09 Jan 2026: Eric Kemp, reformatted.
!
! !INTERFACE:
subroutine finalize_geositbias(findex)

! !USES:
  use LIS_coreMod,     only   : LIS_rc
  use geositbias_forcingMod, only : geositbias_struc
!
! !DESCRIPTION:
!  Routine to cleanup GEOS-ITbias forcing related memory allocations.
!
!EOP
  implicit none

  integer :: findex
  integer :: n

  do n = 1,LIS_rc%nnest
     select case( LIS_rc%met_interp(findex) )

     case( "bilinear" )
        deallocate(geositbias_struc(n)%n111)
        deallocate(geositbias_struc(n)%n121)
        deallocate(geositbias_struc(n)%n211)
        deallocate(geositbias_struc(n)%n221)
        deallocate(geositbias_struc(n)%w111)
        deallocate(geositbias_struc(n)%w121)
        deallocate(geositbias_struc(n)%w211)
        deallocate(geositbias_struc(n)%w221)

     case( "budget-bilinear" )
        deallocate(geositbias_struc(n)%n111)
        deallocate(geositbias_struc(n)%n121)
        deallocate(geositbias_struc(n)%n211)
        deallocate(geositbias_struc(n)%n221)
        deallocate(geositbias_struc(n)%w111)
        deallocate(geositbias_struc(n)%w121)
        deallocate(geositbias_struc(n)%w211)
        deallocate(geositbias_struc(n)%w221)
        deallocate(geositbias_struc(n)%n112)
        deallocate(geositbias_struc(n)%n122)
        deallocate(geositbias_struc(n)%n212)
        deallocate(geositbias_struc(n)%n222)
        deallocate(geositbias_struc(n)%w112)
        deallocate(geositbias_struc(n)%w122)
        deallocate(geositbias_struc(n)%w212)
        deallocate(geositbias_struc(n)%w222)

     case( "neighbor" )
        deallocate(geositbias_struc(n)%n113)
     end select
  enddo
  deallocate(geositbias_struc)

end subroutine finalize_geositbias

