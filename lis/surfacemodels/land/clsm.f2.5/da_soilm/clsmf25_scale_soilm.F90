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
! !ROUTINE: clsmf25_scale_soilm
! \label{clsmf25_scale_soilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clsmf25_scale_soilm(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use clsmf25_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Scales soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP

  type(ESMF_Field)       :: catdefField
  type(ESMF_Field)       :: rzexcField
  type(ESMF_Field)       :: srfexcField
  integer                :: t
  integer                :: status
  real, pointer          :: catdef(:)
  real, pointer          :: rzexc(:)
  real, pointer          :: srfexc(:)

  character*100          :: lsm_state_objs(4)

  call ESMF_StateGet(LSM_State,"Catchment deficit",catdefField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for catdef in clsmf25_getsoilm')
  call ESMF_StateGet(LSM_State,"Root zone excess",rzexcField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for rzexc in clsmf25_getsoilm')
  call ESMF_StateGet(LSM_State,"Surface excess",srfexcField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for srfexc in clsmf25_getsoilm')

  call ESMF_FieldGet(catdefField,localDE=0,farrayPtr=catdef,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for catdef in clsmf25_getsoilm')
  call ESMF_FieldGet(rzexcField,localDE=0,farrayPtr=rzexc,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for rzexc in clsmf25_getsoilm')
  call ESMF_FieldGet(srfexcField,localDE=0,farrayPtr=srfexc,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for srfexc in clsmf25_getsoilm')

!  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     catdef(t)   = catdef(t)/500.0
!     rzexc(t)    = rzexc(t)/50.0
!     srfexc(t)   = srfexc(t)/5.0
!  enddo

end subroutine clsmf25_scale_soilm

