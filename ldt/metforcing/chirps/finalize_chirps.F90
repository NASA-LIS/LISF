!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: finalize_chirps
!  \label{finalize_chirps}
!
! !REVISION HISTORY: 
! 16 July 2015: K. Arsenault Hall; Initial LDT version
! 29 Oct  2025; J Nattala; Updated to support both CHIRPS v2.0 and v3.0
!
! !INTERFACE:
subroutine finalize_chirps(findex)
! !USES:
  use LDT_coreMod,          only : LDT_rc
  use chirps_forcingMod, only : chirps_struc
!
! !DESCRIPTION:
!  Routine to cleanup CHIRPS forcing related memory allocations.   
! 
!EOP

  implicit none

  integer   :: n
  integer   :: findex

  do n=1,LDT_rc%nnest

    select case( LDT_rc%met_gridtransform(findex) )

     case( "average" )   ! Upscaling 
        deallocate(chirps_struc(n)%n111)

     case( "budget-bilinear" )
        deallocate(chirps_struc(n)%n112)
        deallocate(chirps_struc(n)%n122)
        deallocate(chirps_struc(n)%n212)
        deallocate(chirps_struc(n)%n222)
        deallocate(chirps_struc(n)%w112)
        deallocate(chirps_struc(n)%w122)
        deallocate(chirps_struc(n)%w212)
        deallocate(chirps_struc(n)%w222)

      case( "neighbor" )
        deallocate(chirps_struc(n)%n113)

    end select

  enddo

  deallocate(chirps_struc)

end subroutine finalize_chirps
