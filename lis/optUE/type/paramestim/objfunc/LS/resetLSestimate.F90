!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine resetLSestimate()

  use ESMF
  use LIS_coreMod,         only : LIS_rc
  use LSobjFunc_Mod !,       only : ls_ctl
  use LIS_optUEMod,        only : LIS_ObjectiveFunc
  use LIS_logMod,          only : LIS_verify

  implicit none

  type(ESMF_Field)               :: sqerrField
  real, pointer                  :: sqerr(:)
  type(ESMF_Field)               :: numobsField
  real, pointer                  :: numobs(:)

  type(ESMF_Field)               :: modelvField
  real, pointer                  :: modelv(:)
  type(ESMF_Field)               :: nummodelvField
  real, pointer                  :: nummodelv(:)
  
  integer                        :: n
  integer                        :: status

  n = 1

  call ESMF_StateGet(LIS_ObjectiveFunc,"Objective Function Value",sqerrField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sqerrField, localDE=0, farrayPtr=sqerr, rc=status)
  call LIS_verify(status)
  
  sqerr = 0

  call ESMF_StateGet(LIS_ObjectiveFunc,"Number of observations",numobsField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(numobsField, localDE=0, farrayPtr=numobs, rc=status)
  call LIS_verify(status)

  numobs=0

 if(ls_ctl%LSobjfunc_mode.eq.3) then 
     call ESMF_StateGet(LIS_ObjectiveFunc,"Model obspred",modelvField,rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(modelvField, localDE=0, farrayPtr=modelv, rc=status)
     call LIS_verify(status)
     modelv = 0
     
     call ESMF_StateGet(LIS_ObjectiveFunc,"Count Model obspred",nummodelvField,rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(nummodelvField, localDE=0, farrayPtr=nummodelv, rc=status)
     call LIS_verify(status)
     nummodelv = 0 

  endif

 end subroutine resetLSestimate

