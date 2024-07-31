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
! !ROUTINE: GLS_getpeobspred_raint
!  \label{GLS_getpeobspred_raint}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine GLS_getpeobspred_raint(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_soilsMod,  only : LIS_soils
  use LIS_logMod,       only : LIS_verify
  use GLSMod, only           : GLS_ctl

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: Obj_Func
!
! !DESCRIPTION:
!  
! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: lsField
  real, pointer          :: ls_data(:)
  integer                :: t,col,row
  integer                :: i
  integer                :: status

  n = 1
  call ESMF_StateGet(Obj_Func,"Landslide estimate",lsField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lsField,localDE=0,farrayPtr=ls_data,rc=status)
  call LIS_verify(status)

#if 0 
! this is for the final landslide prediction
  do t=1,LIS_rc%ntiles(n)
     ls_data(t) = GLS_ctl%ls_accum_status(t)
  enddo
#endif

  do t=1,LIS_rc%ntiles(n)
     ls_data(t) = GLS_ctl%ls_status(t)
  enddo
end subroutine GLS_getpeobspred_raint



