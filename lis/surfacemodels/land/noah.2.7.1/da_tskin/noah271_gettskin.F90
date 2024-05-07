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
! !ROUTINE: noah271_gettskin
! \label{noah271_gettskin}
!
! !REVISION HISTORY:
! 13Nov2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah271_gettskin(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the tskin related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
!  type(ESMF_Field)       :: avgsurftField
  type(ESMF_Field)       :: soiltField
  integer                :: t
  integer                :: status
!  real, allocatable          :: avgsurft(:)
  real, allocatable          :: soilt(:)

!  call ESMF_StateGet(LSM_State,"Skin Temperature",avgsurftField,rc=status)
!  call LIS_verify(status)

!  call ESMF_FieldGet(avgsurftField,localDE=0,farrayPtr=avgsurft,rc=status)

  call ESMF_StateGet(LSM_State,"Soil Temperature",soiltField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(soiltField,localDE=0,farrayPtr=soilt,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
!     avgsurft(t) = noah271_struc(n)%noah(t)%t1
     soilt(t) = noah271_struc(n)%noah(t)%stc(1)
  enddo

end subroutine noah271_gettskin

