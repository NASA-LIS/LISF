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
! !ROUTINE: clsmf25_qcsoilm
! \label{clsmf25_qcsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clsmf25_qcsoilm(n, LSM_State)

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
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: catdefField
  integer                :: t
  real                   :: smmax1,smmin1
  integer                :: status
  real, pointer          :: catdef(:)

  call ESMF_StateGet(LSM_State,"Catchment deficit",catdefField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for catdef in clsmf25_getsoilm')

  call ESMF_FieldGet(catdefField,localDE=0,farrayPtr=catdef,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for catdef in clsmf25_getsoilm')

  call ESMF_AttributeGet(catdefField,"Max Value",smmax1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Max Value failed in clsmf25_qcsoilm")
  call ESMF_AttributeGet(catdefField,"Min Value",smmin1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Min Value failed in clsmf25_qcsoilm")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(catdef(t).gt.smmax1) catdef(t) = smmax1
     if(catdef(t).lt.smmin1) catdef(t) = smmin1
  enddo

end subroutine clsmf25_qcsoilm

