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
! !ROUTINE: noahmp401_write_veg
! \label{noahmp401_write_veg}
!
! !REVISION HISTORY:
! 13 Feb 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp401_write_veg(ftn,n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noahmp401_lsmMod
  use LIS_historyMod, only : LIS_writevar_restart
  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: ftn
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the esioisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  integer :: t
  real, allocatable    :: tmp(:)

  allocate(tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = NOAHMP401_struc(n)%noahmp401(t)%smc(1)
  enddo  
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = NOAHMP401_struc(n)%noahmp401(t)%smc(2)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = NOAHMP401_struc(n)%noahmp401(t)%smc(3)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = NOAHMP401_struc(n)%noahmp401(t)%smc(4)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)
  deallocate(tmp)

end subroutine noahmp401_write_veg

