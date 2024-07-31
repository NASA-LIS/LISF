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
! !ROUTINE: clsmf25_write_tws
! \label{clsmf25_write_tws}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 29Sep2011: Ben Zaitchik: Applied to GRACE
!
! !INTERFACE:
subroutine clsmf25_write_tws(ftn,n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use clsmf25_lsmMod
  use LIS_historyMod, only : LIS_writevar_restart
  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: ftn
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the TWS related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  integer :: t
  real, allocatable :: tmp(:)
  
  allocate(tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = clsmf25_struc(n)%cat_progn(t)%catdef
  enddo  
  call LIS_writevar_restart(ftn,n,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = clsmf25_struc(n)%cat_progn(t)%srfexc
  enddo
  call LIS_writevar_restart(ftn,n,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = clsmf25_struc(n)%cat_progn(t)%wesn(1)
  enddo
  call LIS_writevar_restart(ftn,n,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = clsmf25_struc(n)%cat_progn(t)%wesn(2)
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = clsmf25_struc(n)%cat_progn(t)%wesn(3)
  enddo
    
  call LIS_writevar_restart(ftn,n,tmp)
  deallocate(tmp)

end subroutine clsmf25_write_tws

