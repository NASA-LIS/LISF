!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LMobjFunc_Mod
!BOP
! !MODULE: LMobjFunc_Mod
! 
! !DESCRIPTION: 
!  This module provides the objects and methods to compute a residuals array
!  as required  by Levenberg-Marquadt algorithm (LM)
!  for use in parameter estimation. 
! 
! !REVISION HISTORY: 
!  
! 3 Aug 2009: Ken Harrison; Initial implementation
! 
  use ESMF

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public ::  initializeLMObjFunc
!EOP
contains

!BOP
! !ROUTINE: intializeLMObjFunc
! \label{initializeLMObjFunc}
! 
! !INTERFACE: 
  subroutine initializeLMObjFunc()
! !USES: 
    use LIS_coreMod,         only : LIS_vecTile
    use LIS_optUEMod,        only : LIS_ObjectiveFunc
    use LIS_logMod,          only : LIS_verify
! 
! !DESCRIPTION:
!  This method initializes the objects to be used in the LM residuals array 
!  computations
!EOP    
    implicit none

    type(ESMF_ArraySpec) :: arrspec1
    type(ESMF_Field)     :: numobsField
    integer,   pointer   :: numobs(:)
    type(ESMF_Field)     :: errField
    real,   pointer      :: err(:,:)  !index first by tile, then ob
    integer              :: n 
    integer              :: status

    n = 1

!Requires matrix of rank 2 (1+ that of grid) as the residuals array is needed
    call ESMF_ArraySpecSet(arrspec1,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
!Hardcoded upper bound of obs dimension; correct to maximum after first pass
    errField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Error field",rc=status, &
         ungriddedLBound=(/1/), ungriddedUBound=(/33/))
    call LIS_verify(status)

    call ESMF_FieldGet(errField, localDE=0, farrayPtr=err, rc=status)
    call LIS_verify(status)
    err = 0.0

    call ESMF_StateAdd(LIS_ObjectiveFunc, (/errField/), rc=status)
    call LIS_verify(status)  

!Store count of obs
    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)
    numobsField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Number of observations",rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(numobsField, localDE=0, farrayPtr=numobs, rc=status)
    call LIS_verify(status)
    numobs = 0

    call ESMF_StateAdd(LIS_ObjectiveFunc, (/numobsField/), rc=status)
    call LIS_verify(status)  

  end subroutine initializeLMObjFunc
  

end module LMobjFunc_Mod
