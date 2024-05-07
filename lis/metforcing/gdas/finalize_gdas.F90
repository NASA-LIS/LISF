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
! !MODULE: finalize_gdas
! \label{finalize_gdas}
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine finalize_gdas(findex)
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use gdas_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GDAS forcing. 
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing
!  \end{description}
!EOP  
  implicit none

  integer, intent(in) :: findex

  integer :: n 
  integer :: rc

  do n=1,LIS_rc%nnest
     deallocate(gdas_struc(n)%n111, stat=rc)
     deallocate(gdas_struc(n)%n121, stat=rc)
     deallocate(gdas_struc(n)%n211, stat=rc)
     deallocate(gdas_struc(n)%n221, stat=rc)
     deallocate(gdas_struc(n)%w111, stat=rc)
     deallocate(gdas_struc(n)%w121, stat=rc)
     deallocate(gdas_struc(n)%w211, stat=rc)
     deallocate(gdas_struc(n)%w221, stat=rc)

     deallocate(gdas_struc(n)%n112, stat=rc)
     deallocate(gdas_struc(n)%n122, stat=rc)
     deallocate(gdas_struc(n)%n212, stat=rc)
     deallocate(gdas_struc(n)%n222, stat=rc)
     deallocate(gdas_struc(n)%w112, stat=rc)
     deallocate(gdas_struc(n)%w122, stat=rc)
     deallocate(gdas_struc(n)%w212, stat=rc)
     deallocate(gdas_struc(n)%w222, stat=rc)
  enddo

  deallocate(gdas_struc)

end subroutine finalize_gdas
