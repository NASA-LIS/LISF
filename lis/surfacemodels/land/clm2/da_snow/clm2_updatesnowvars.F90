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
! !ROUTINE: clm2_updatesnowvars
!  \label{clm2_updatesnowvars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 06Oct2008: Gabrielle De Lannoy: Update the state vars for scf/a DA in CLM
! 04Nov2008: Gabrielle De Lannoy: only do update for 2 total state vars in CLM
!
!
! !INTERFACE:
subroutine clm2_updatesnowvars(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use clm2_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
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
!  \item[LSM\_Incr\_State] ESMF State container for LSM state increments \newline
!  \end{description}
! 
!EOP

  type(ESMF_Field)       :: h2osnoField, h2osnoIncrField
  type(ESMF_Field)       :: snowdpField, snowdpIncrField

  integer                :: t
  integer                :: status

  real, pointer          :: h2osno(:), h2osnoincr(:)  ! total water in snow (kg m^-2)
  real, pointer          :: snowdp(:), snowdpincr(:)  ! total snow depth

! ----------------------------------------------------------------------------------

  call ESMF_StateGet(LSM_State,"Total Snow Water",h2osnoField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Total Snow Depth",snowdpField,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"Total Snow Water",h2osnoIncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Total Snow Depth",snowdpIncrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(h2osnoField,localDE=0,farrayPtr=h2osno,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snowdpField,localDE=0,farrayPtr=snowdp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(h2osnoIncrField,localDE=0,farrayPtr=h2osnoIncr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snowdpIncrField,localDE=0,farrayPtr=snowdpIncr,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
!     print*, 'upd ',t, h2osno(t), h2osnoincr(t), snowdp(t), snowdpincr(t)
     h2osno(t) = h2osno(t) + h2osnoincr(t)
     snowdp(t) = snowdp(t) + snowdpincr(t)
     
  enddo
end subroutine clm2_updatesnowvars

