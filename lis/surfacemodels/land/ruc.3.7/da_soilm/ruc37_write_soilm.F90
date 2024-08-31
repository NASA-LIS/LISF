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
! !ROUTINE: RUC37_write_soilm
! \label{RUC37_write_soilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine RUC37_write_soilm(ftn, n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc
  use LIS_historyMod, only : LIS_writevar_restart
  use RUC37_lsmMod,   only : RUC37_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: ftn
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
  real,allocatable   :: soilm(:)
  integer :: m, t
 
  allocate(soilm(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  do m=1,RUC37_struc(n)%nsoil
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        soilm(t)= RUC37_struc(n)%ruc37(t)%smc(m)
     enddo
     call LIS_writevar_restart(ftn,n,1,soilm)
  enddo

  deallocate(soilm)

end subroutine RUC37_write_soilm

