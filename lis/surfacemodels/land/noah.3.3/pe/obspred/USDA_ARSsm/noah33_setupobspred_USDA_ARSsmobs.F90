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
! !ROUTINE: noah33_setupobspred_USDA_ARSsmObs
!  \label{noah33_setupobspred_USDA_ARSsmObs}
!
! !REVISION HISTORY:
! 16 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah33_setupobspred_USDA_ARSsmObs(OBSPred)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_vecPatch
  use LIS_logMod,       only : LIS_verify

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: OBSPred
!
! !DESCRIPTION:
!  
!  This routine creates an entry in the Obs pred object used for 
!  parameter estimation
! 
!EOP
  integer                :: n
  type(ESMF_ArraySpec)   :: realarrspec
  type(ESMF_Field)       :: sfsmField
  integer                :: status

  n = 1
  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
       rc=status)
  call LIS_verify(status)

  sfsmField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecPatch(n,LIS_rc%lsm_index), &
       name="USDA ARS Surface Soil Moisture", rc=status)
  call LIS_verify(status)
  
  call ESMF_StateAdd(OBSPred,(/sfsmField/),rc=status)
  call LIS_verify(status)

end subroutine noah33_setupobspred_USDA_ARSsmObs

