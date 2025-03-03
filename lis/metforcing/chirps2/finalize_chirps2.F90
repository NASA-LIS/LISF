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
! !MODULE: finalize_chirps2
!  \label{finalize_chirps2}
!
! !REVISION HISTORY: 
! 16 July 2015: K. Arsenault Hall; Initial LIS version
!
! !INTERFACE:
subroutine finalize_chirps2(findex)
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use chirps2_forcingMod, only : chirps2_struc
!
! !DESCRIPTION:
!  Routine to cleanup CHIRPS 2 forcing related memory allocations.   
! 
!EOP

  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LIS_rc%nnest

     select case( LIS_rc%met_interp(findex) )

     case( "average" )   ! Upscaling 
        deallocate(chirps2_struc(n)%n111)

     case( "budget-bilinear" )
        deallocate(chirps2_struc(n)%n112)
        deallocate(chirps2_struc(n)%n122)
        deallocate(chirps2_struc(n)%n212)
        deallocate(chirps2_struc(n)%n222)
        deallocate(chirps2_struc(n)%w112)
        deallocate(chirps2_struc(n)%w122)
        deallocate(chirps2_struc(n)%w212)
        deallocate(chirps2_struc(n)%w222)

      case( "neighbor" ) 
        deallocate(chirps2_struc(n)%n113)

    end select

  enddo

  deallocate(chirps2_struc)

end subroutine finalize_chirps2


