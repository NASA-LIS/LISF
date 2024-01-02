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
! !ROUTINE: clm2_getsnowvars
! \label{clm2_getsnowvars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 06Oct2008: Gabrielle De Lannoy: state variables for scf/a EnKF in CLM2.0
! 04Nov2008: Gabrielle De Lannoy: reduce state var to h2osno and snowdp
!
! !INTERFACE:
subroutine clm2_getsnowvars(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use clm2_lsmMod
  use clm2_varcon,    only : denh2o

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
!  NOW: swe (h2osno) and snod (snowdp) are diagnostics, also
!       dz and z (snow layer thickness, level) are 'diagnostic'!
!  ==> NOTE: possibly snow layer temp could be included as prognostic
!      (cfr. KA)
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

  real, pointer          :: h2osno(:)  ! clm2.0 total snow water (kg m^-2)
  real, pointer          :: snowdp(:)  ! clm2.0 total snow depth

! ----------------------------------------------------------------------------
 
  call ESMF_StateGet(LSM_State,"Total Snow Water",h2osnoField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Total Snow Depth",snowdpField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(h2osnoField,localDE=0,farrayPtr=h2osno,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snowdpField,localDE=0,farrayPtr=snowdp,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     h2osno(t) = clm2_struc(n)%clm(t)%h2osno
     snowdp(t) = clm2_struc(n)%clm(t)%snowdp
  enddo

end subroutine clm2_getsnowvars

