!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_capa
! \label{finalize_capa}
! 
! !REVISION HISTORY: 
! 23 Nov 2009: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine finalize_capa(findex)
! !USES:

  use LIS_coreMod,     only : LIS_rc
  use capa_forcingMod, only : capa_struc

  implicit none
! !ARGUMENTS:
  integer, intent(IN) :: findex

! !DESCRIPTION: 
!  Routine to cleanup CAPA forcing related memory allocations.   
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!EOP

  integer :: n
  do n=1,LIS_rc%nnest  
     if ( LIS_rc%met_interp(findex) == "budget-bilinear" ) then
        deallocate(capa_struc(n)%n112)
        deallocate(capa_struc(n)%n122)
        deallocate(capa_struc(n)%n212)
        deallocate(capa_struc(n)%n222)
        deallocate(capa_struc(n)%w112)
        deallocate(capa_struc(n)%w122)
        deallocate(capa_struc(n)%w212)
        deallocate(capa_struc(n)%w222)
     endif
  enddo
  deallocate(capa_struc)

end subroutine finalize_capa
