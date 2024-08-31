!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine resetLMestimate()

  use ESMF
  use LIS_coreMod,         only : LIS_rc
  use LIS_optUEMod, only : LIS_ObjectiveFunc
  use LIS_logMod,          only : LIS_verify

  implicit none

  type(ESMF_Field)               :: errField
  real, pointer                  :: err(:,:)
  type(ESMF_Field)               :: numobsField
  integer,   pointer             :: numobs(:)
  integer                        :: n
  integer                        :: status

  n = 1

  call ESMF_StateGet(LIS_ObjectiveFunc,"Error field",errField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(errField, localDE=0, farrayPtr=err, rc=status)
  call LIS_verify(status)
  
  err = 0.0

  call ESMF_StateGet(LIS_ObjectiveFunc,"Number of observations",numobsField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(numobsField, localDE=0, farrayPtr=numobs, rc=status)
  call LIS_verify(status)

  numobs=0
 end subroutine resetLMestimate

