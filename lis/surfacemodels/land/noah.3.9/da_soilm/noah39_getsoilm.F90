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
! !ROUTINE: Noah39_getsoilm
! \label{Noah39_getsoilm}
!
! !REVISION HISTORY:
! 21Oct2018: Mahdi Navari; Sujay Kumar ; Initial Specification
! 20 Mar 2025: Eric Kemp; Changed to use top soil layer only.  Based on
!  Noah39 code from Sujay Kumar and Wanshu Nie used for, e.g.,
!  HydroGlobe.
! !INTERFACE:
subroutine Noah39_getsoilm(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only: LIS_rc
  use LIS_logMod,  only: LIS_verify
  use Noah39_lsmMod, only: noah39_struc

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
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in Noah39_getsoilm')
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm1 in Noah39_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = noah39_struc(n)%noah(t)%smc(1)
  enddo

end subroutine Noah39_getsoilm
