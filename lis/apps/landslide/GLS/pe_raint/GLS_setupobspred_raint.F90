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
! !ROUTINE: GLS_setupobspred_raint
!  \label{GLS_setupobspred_raint}
!
! !REVISION HISTORY:
! 16 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine GLS_setupobspred_raint(OBSPred)
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
! 
!EOP
  integer                :: n
  type(ESMF_ArraySpec)   :: realarrspec
  type(ESMF_Field)       :: lsField
  integer                :: status

  n = 1
  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
       rc=status)
  call LIS_verify(status)

  lsField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecTile(n), &
       name="Landslide estimate", rc=status)
  call LIS_verify(status)
  
  call ESMF_StateAdd(OBSPred,(/lsField/),rc=status)
  call LIS_verify(status)

end subroutine GLS_setupobspred_raint

