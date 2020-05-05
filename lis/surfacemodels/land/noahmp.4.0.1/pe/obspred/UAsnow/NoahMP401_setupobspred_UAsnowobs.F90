!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP401_setupobspred_UAsnowobs
!  \label{NoahMP401_setupobspred_UAsnowobs}
!
! !REVISION HISTORY:
! 2 May 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine NoahMP401_setupobspred_UAsnowobs(OBSPred)
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
!  This routine creates an entry in the Obs pred object from 
!  NoahMP4.0.1 used for parameter estimation
! 
!EOP
  integer                :: n
  type(ESMF_ArraySpec)   :: realarrspec
  type(ESMF_Field)       :: snowField
  integer                :: status

  n = 1
  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
       rc=status)
  call LIS_verify(status)

  snowField = ESMF_FieldCreate(arrayspec=realarrspec, &
       grid=LIS_vecPatch(n,LIS_rc%lsm_index), &
       name="UA_snow", rc=status)
  call LIS_verify(status)
  
  call ESMF_StateAdd(OBSPred,(/snowField/),rc=status)
  call LIS_verify(status)

end subroutine NoahMP401_setupobspred_UAsnowobs

