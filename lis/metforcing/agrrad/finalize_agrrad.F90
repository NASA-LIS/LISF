!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_agrrad
! \label{finalize_agrrad}
!
! !REVISION HISTORY:
! 29 Jul 2005; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine finalize_agrrad(findex)
! !USES:

  use LIS_coreMod,       only : LIS_rc
  use agrrad_forcingMod, only : agrrad_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  This routine deallocates and cleans up allocated memory 
!  structures for AGRMET radiation forcing implementation.
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing
!  \end{description}
!EOP

  integer :: n
  
  do n=1,LIS_rc%nnest  
     if(LIS_rc%met_interp(findex).eq."bilinear") then 
        deallocate(agrrad_struc(n)%n111)
        deallocate(agrrad_struc(n)%n121)
        deallocate(agrrad_struc(n)%n211)
        deallocate(agrrad_struc(n)%n221)
        deallocate(agrrad_struc(n)%w111)
        deallocate(agrrad_struc(n)%w121)
        deallocate(agrrad_struc(n)%w211)
        deallocate(agrrad_struc(n)%w221)
     endif
  enddo
  deallocate(agrrad_struc)

end subroutine finalize_agrrad
