!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LLobjFunc_Mod
!BOP
! !MODULE: LLobjFunc_Mod
! 
! !DESCRIPTION: 
!  This module provides the objects and methods to compute a Least Square (LL) metric
!  for use in parameter estimation. 
! 
! !REVISION HISTORY: 
!  
! 15 Jul 2009: Ken Harrison and Soni Yatheendradas; Initial implementation
! 
! !USES:
  use ESMF

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public ::  initializeLLObjFunc
!EOP
contains

!BOP
! !ROUTINE: intializeLLObjFunc
! \label{initializeLLObjFunc}
! 
! !INTERFACE: 
  subroutine initializeLLObjFunc()
! !USES: 
    use LIS_coreMod,         only : LIS_vecTile
    use LIS_optUEMod,        only : LIS_ObjectiveFunc
    use LIS_logMod,          only : LIS_verify
! 
! !DESCRIPTION:
!  This method initializes the objects to be used in the least square 
!  computations
!EOP    
    implicit none

    type(ESMF_ArraySpec) :: arrspec1
    type(ESMF_Field)     :: objfuncField
    type(ESMF_Field)     :: minField
    type(ESMF_Field)     :: maxField
    real,   pointer      :: objfunc(:)
    integer              :: n 
    integer              :: status

    n = 1

    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    objfuncField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Objective Function Value",rc=status)
    call LIS_verify(status)

    minField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Min Criteria Value",rc=status)
    call LIS_verify(status)

    maxField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Max Criteria Value",rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(objfuncField, localDE=0, farrayPtr=objfunc, rc=status)
    call LIS_verify(status)
    objfunc = 0.0

    call ESMF_StateAdd(LIS_ObjectiveFunc, (/objfuncField/), rc=status)
    call LIS_verify(status)  

    call ESMF_StateAdd(LIS_ObjectiveFunc, (/minField/), rc=status)
    call LIS_verify(status)  

    call ESMF_StateAdd(LIS_ObjectiveFunc, (/maxField/), rc=status)
    call LIS_verify(status)  
    
  end subroutine initializeLLObjFunc
  

end module LLobjFunc_Mod
