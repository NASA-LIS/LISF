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
! !ROUTINE: noah32_write_soilm
! \label{noah32_write_soilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine noah32_write_soilm(ftn,n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noah32_lsmMod
  use LIS_historyMod, only : LIS_writevar_restart
  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: ftn
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  integer :: t
  real    :: tmp(LIS_rc%ntiles(n))
  
  do t=1,LIS_rc%ntiles(n)
     tmp(t) = noah32_struc(n)%noah(t)%smc(1)
  enddo  
  call LIS_writevar_restart(ftn,n,tmp)

  do t=1,LIS_rc%ntiles(n)
     tmp(t) = noah32_struc(n)%noah(t)%smc(2)
  enddo
  call LIS_writevar_restart(ftn,n,tmp)

  do t=1,LIS_rc%ntiles(n)
     tmp(t) = noah32_struc(n)%noah(t)%smc(3)
  enddo
  call LIS_writevar_restart(ftn,n,tmp)

  do t=1,LIS_rc%ntiles(n)
     tmp(t) = noah32_struc(n)%noah(t)%smc(4)
  enddo
  call LIS_writevar_restart(ftn,n,tmp)

end subroutine noah32_write_soilm

