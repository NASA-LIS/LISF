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
! !ROUTINE: noah271_setupobspred_wgPBMRsm
!  \label{noah271_setupobspred_wgPBMRsm}
!
! !REVISION HISTORY:
! 16 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah271_setupobspred_wgPBMRsm(OBSPred)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_vecTile
  use LIS_logMod,       only : LIS_verify

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: OBSPred
!
! !DESCRIPTION:
!  
!  This routine generates a placeholder entry in the OBSPRED object
!  for creating the obspred from the LSM for the walnut
!  gulch PBMR soil moisture observations.  
! 
!EOP
  integer                :: n
  type(ESMF_ArraySpec)   :: realarrspec
  type(ESMF_Field)       :: smField
  integer                :: status

  n = 1
  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
       rc=status)
  call LIS_verify(status)

  smField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecTile(n), &
       name="PBMR soil moisture", rc=status)
  call LIS_verify(status)
  
  call ESMF_StateAdd(OBSPred,(/smField/),rc=status)
  call LIS_verify(status)

end subroutine noah271_setupobspred_wgPBMRsm

