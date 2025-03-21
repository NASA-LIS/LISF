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
! !ROUTINE: NoahMP401_getsoilm
! \label{NoahMP401_getsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari; Modified for NoahMP401
! 20 Mar 2025: Eric Kemp; Changed to use top soil layer only.  Based on
!  code from Sujay Kumar and Wanshu Nie used for, e.g., HydroGlobe.
! !INTERFACE:
subroutine NoahMP401_getsoilm(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only: LIS_rc
  use LIS_logMod,  only: LIS_verify
  use NoahMP401_lsmMod, only: noahmp401_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State), intent(in) :: LSM_State
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
  type(ESMF_Field)       :: sm1Field
  integer                :: status
  real, pointer          :: soilm1(:)
  integer                :: t

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in NoahMP401_getsoilm')
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm1 in NoahMP401_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
  enddo

end subroutine NoahMP401_getsoilm
