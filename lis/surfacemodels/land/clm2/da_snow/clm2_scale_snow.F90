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
! !ROUTINE: clm2_scale_snow
! \label{clm2_scale_snow}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clm2_scale_snow(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use clm2_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the snow related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP

  type(ESMF_Field)       :: h2osnoField
  type(ESMF_Field)       :: snowdpField

  integer                :: t
  integer                :: status

  real, allocatable          :: h2osno(:)  ! clm2.0 total snow water (kg m^-2)
  real, allocatable          :: snowdp(:)  ! clm2.0 total snow depth

! ----------------------------------------------------------------------------
 
#if 0 
  call ESMF_StateGet(LSM_State,"Total Snow Water",h2osnoField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Total Snow Depth",snowdpField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(h2osnoField,localDE=0,farrayPtr=h2osno,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snowdpField,localDE=0,farrayPtr=snowdp,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     snowdp(t) = clm2_struc(n)%clm(t)%snowdp*1000.0
  enddo
#endif
end subroutine clm2_scale_snow

