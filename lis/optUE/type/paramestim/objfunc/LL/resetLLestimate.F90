!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine resetLLestimate()

  use ESMF
  use LIS_coreMod,         only : LIS_rc
  use LIS_optUEMod, only : LIS_ObjectiveFunc
  use LIS_logMod,          only : LIS_verify

  implicit none

  type(ESMF_Field)               :: lnLLField
  real, pointer                  :: lnLL(:)
  integer                        :: n
  integer                        :: status

  n = 1

  call ESMF_StateGet(LIS_ObjectiveFunc,"Objective Function Value",lnLLField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(lnLLField, localDE=0, farrayPtr=lnLL, rc=status)
  call LIS_verify(status)
  
  lnLL = 0.0
 end subroutine resetLLestimate

