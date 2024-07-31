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
! 
! !ROUTINE: computeLMestimate
! \label{computeLMestimate}
! 
! !REVISION HISTORY: 
!  
! 03 Aug 2009: Ken Harrison; Initial implementation
! 
! !INTERFACE: 
 subroutine computeLMestimate()
! !USES: 
   use ESMF
   use LIS_coreMod,         only : LIS_rc
   use LIS_optUEMod, only : LIS_ObjectiveFunc
   use LIS_logMod,          only : LIS_logUnit, LIS_verify
! 
! !DESCRIPTION: 
!  This routine computes the objective function values to be used in the 
!  optimization algorithm.
! 
!EOP
   implicit none
   
   type(ESMF_Field)      :: numobsField,errField
   integer, pointer      :: numobs(:)
   real, pointer         :: err(:,:)
   integer               :: t
   integer               :: m
   integer               :: n 
   integer               :: status
   
   n = 1
   call ESMF_StateGet(LIS_ObjectiveFunc,"Number of observations",numobsField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(numobsField, localDE=0, farrayPtr=numobs, rc=status)
   call LIS_verify(status)

   call ESMF_StateGet(LIS_ObjectiveFunc,"Error field",errField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(errField, localDE=0, farrayPtr=err, rc=status)
   call LIS_verify(status)

! Nothing to actually do....da da da da da.....

 end subroutine computeLMestimate
