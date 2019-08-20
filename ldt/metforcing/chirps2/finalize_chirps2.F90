!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.1
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: finalize_chirps2
!  \label{finalize_chirps2}
!
! !REVISION HISTORY: 
! 16 July 2015: K. Arsenault Hall; Initial LDT version
!
! !INTERFACE:
subroutine finalize_chirps2(findex)
! !USES:
  use LDT_coreMod,          only : LDT_rc
  use chirps2_forcingMod, only : chirps2_struc
!
! !DESCRIPTION:
!  Routine to cleanup CHIRPS 2 forcing related memory allocations.   
! 
!EOP

  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LDT_rc%nnest

    select case( LDT_rc%met_gridtransform(findex) ) 

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


