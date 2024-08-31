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
! !ROUTINE: jules5x_write_soilm
! \label{jules5x_write_soilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 21 Dec 2018: Mahdi Navari; Modified for JULES 5.3
! 12 Feb 2020: Shugong Wang; Modifed for JULES 5.5 and later version
!
! !INTERFACE:
subroutine jules5x_write_soilm(ftn,n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use jules5x_lsmMod
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
  real, allocatable    :: tmp(:)

  allocate(tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = jules5x_struc(n)%jules5x(t)%smcl_soilt(1)
  enddo  
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = jules5x_struc(n)%jules5x(t)%smcl_soilt(2)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = jules5x_struc(n)%jules5x(t)%smcl_soilt(3)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = jules5x_struc(n)%jules5x(t)%smcl_soilt(4)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)
  deallocate(tmp)

end subroutine jules5x_write_soilm

