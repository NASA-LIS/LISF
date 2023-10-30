!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: finalize_geosit
! \label{finalize_geosit}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 20 Apr 2023: David Mocko,  initial code (based on merra2)
!
! !INTERFACE:
      subroutine finalize_geosit(findex)

! !USES:
      use LIS_coreMod,     only   : LIS_rc
      use geosit_forcingMod, only : geosit_struc
!
! !DESCRIPTION:
!  Routine to cleanup GEOS-IT forcing related memory allocations.
!
!EOP
      implicit none

      integer :: findex
      integer :: n

      do n = 1,LIS_rc%nnest
         select case( LIS_rc%met_interp(findex) )

         case( "bilinear" )
            deallocate(geosit_struc(n)%n111)
            deallocate(geosit_struc(n)%n121)
            deallocate(geosit_struc(n)%n211)
            deallocate(geosit_struc(n)%n221)
            deallocate(geosit_struc(n)%w111)
            deallocate(geosit_struc(n)%w121)
            deallocate(geosit_struc(n)%w211)
            deallocate(geosit_struc(n)%w221)

         case( "budget-bilinear" )
            deallocate(geosit_struc(n)%n111)
            deallocate(geosit_struc(n)%n121)
            deallocate(geosit_struc(n)%n211)
            deallocate(geosit_struc(n)%n221)
            deallocate(geosit_struc(n)%w111)
            deallocate(geosit_struc(n)%w121)
            deallocate(geosit_struc(n)%w211)
            deallocate(geosit_struc(n)%w221)
            deallocate(geosit_struc(n)%n112)
            deallocate(geosit_struc(n)%n122)
            deallocate(geosit_struc(n)%n212)
            deallocate(geosit_struc(n)%n222)
            deallocate(geosit_struc(n)%w112)
            deallocate(geosit_struc(n)%w122)
            deallocate(geosit_struc(n)%w212)
            deallocate(geosit_struc(n)%w222)

         case( "neighbor" )
            deallocate(geosit_struc(n)%n113)
         end select
      enddo
      deallocate(geosit_struc)

      end subroutine finalize_geosit

