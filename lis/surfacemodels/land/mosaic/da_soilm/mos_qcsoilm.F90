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
! !ROUTINE: mos_qcsoilm
! \label{mos_qcsoilm}
!
! !REVISION HISTORY:
! 06 Oct 2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine mos_qcsoilm(n, LSM_State)

! !USES:
  use ESMF
  use mos_lsmMod
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Qc's the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field

  integer                :: t
  integer                :: status
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)

  real                   :: smmax1,smmax2,smmax3
  real                   :: smmin1,smmin2,smmin3
 
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(sm1Field,"Max Value",smmax1,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm1Field,"Min Value",smmin1,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm2Field,"Max Value",smmax2,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm2Field,"Min Value",smmin2,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm3Field,"Max Value",smmax3,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sm3Field,"Min Value",smmin3,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(soilm1(t).gt.smmax1) &
          soilm1(t) = smmax1
     if(soilm2(t).gt.smmax2) &
          soilm2(t) = smmax2
     if(soilm3(t).gt.smmax3) &
          soilm3(t) = smmax3

     if(soilm1(t).lt.smmin1) soilm1(t) = smmin1
     if(soilm2(t).lt.smmin2) soilm2(t) = smmin2
     if(soilm3(t).lt.smmin3) soilm3(t) = smmin3
  enddo

end subroutine mos_qcsoilm

