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
! !ROUTINE: updateLMestimate
!  \label{updateLMestimate}
! 
! !REVISION HISTORY: 
!  
! 15 Jul 2009: Sujay Kumar; Initial implementation
!
! !INTERFACE: 
subroutine updateLMestimate()
! !USES: 
  use ESMF
  use LIS_coreMod,         only : LIS_rc, LIS_domain
  use LIS_optUEMod,        only : LIS_ObjectiveFunc
  use LIS_logMod,          only : LIS_verify
  use LIS_PE_HandlerMod,       only : LIS_PEOBS_State, LIS_PEOBSPred_State

! 
! !DESCRIPTION:
!  This method updates the running sum of the squared error. Here the 
!  squared error is computed using the observations used for parameter
!  estimation and the model simulated values (obspred). 
! 
!EOP
  implicit none

  character*100, allocatable         :: peobsname(:), peobspredname(:)
  type(ESMF_Field)               :: numobsField,errField, peobsField, peobspredField
  real, pointer                  :: err(:,:), peobs(:), obspred(:)
  integer, pointer               :: numobs(:)
  integer                        :: t, index1
  integer                        :: n
  integer                        :: status

  n = 1

  allocate(peobsname(1))
  allocate(peobspredname(1))

  call ESMF_StateGet(LIS_PEOBS_State, itemNameList=peobsname,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_PEOBSPred_State, itemNameList=peobspredname, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_PEOBS_State, trim(peobsname(1)), peobsField, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(peobsField, localDE=0, farrayPtr=peobs, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_PEOBSPred_State, trim(peobspredname(1)), peobspredField, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(peobspredField, localDE=0, farrayPtr=obspred, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_ObjectiveFunc,"Number of observations",numobsField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(numobsField, localDE=0, farrayPtr=numobs, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_ObjectiveFunc,"Error field",errField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(errField, localDE=0, farrayPtr=err, rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     numobs(t)=numobs(t)+1
     err(t,numobs(t)) = peobs(index1) - obspred(t)
  enddo

  deallocate(peobsname)
  deallocate(peobspredname)
  
 end subroutine updateLMestimate

