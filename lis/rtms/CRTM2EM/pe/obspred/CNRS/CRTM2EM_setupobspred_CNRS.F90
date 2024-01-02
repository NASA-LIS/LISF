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
! !ROUTINE: CRTM2EM_setupobspred_CNRS
!  \label{CRTM2EM_setupobspred_CNRS}
!
! !REVISION HISTORY:
! 16 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine CRTM2EM_setupobspred_CNRS(OBSPred)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_vecTile
  use LIS_logMod,       only : LIS_verify

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: OBSPred
!
! !DESCRIPTION:
!  This routine creates an entry in the Obs pred object used for 
!  parameter estimation!  
! 
!EOP
  integer                :: n
  type(ESMF_ArraySpec)   :: realarrspec
  type(ESMF_Field)       :: emField
  integer                :: status

  n = 1
  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
       rc=status)
  call LIS_verify(status)

!   for now, one channel at a time; if multiple, then rank=2
!   call ESMF_ArraySpecSet(realarrspec, rank=2,typekind=ESMF_TYPEKIND_R4,&
!        rc=status)
!   call LIS_verify(status)

  emField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecTile(n), &
       name="Emissivity", rc=status)
  call LIS_verify(status)
  
  call ESMF_StateAdd(OBSPred,(/emField/),rc=status)
  call LIS_verify(status)

end subroutine CRTM2EM_setupobspred_CNRS

