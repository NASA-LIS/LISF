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
! !ROUTINE: NoahMP36_getpeobspred_ISMNsmobs
!  \label{NoahMP36_getpeobspred_ISMNsmobs}
!
! !REVISION HISTORY:
! 02 Feb 2018: Soni Yatheendradas; Initial Specification
!
! !INTERFACE:
subroutine NoahMP36_getpeobspred_ISMNsmobs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_soilsMod,  only : LIS_soils
  use NoahMP36_lsmMod, only : NoahMP36_struc
  use LIS_logMod,       only : LIS_verify, LIS_logunit

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: Obj_Func
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP3.6 model variables. 
! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: smcField
  real, pointer          :: smc(:)
!  type(ESMF_Field)       :: smstdField
!  real, pointer          :: smstd(:)
  integer                :: t
  integer                :: i
  integer                :: status

!  write(LIS_logunit,*) '[INFO] Here 1 in NoahMP36_getpeobspred_ISMNsmobs '

  n = 1

  call ESMF_StateGet(Obj_Func,"ISMN_sm",smcField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=smc,rc=status)
  call LIS_verify(status)

!  write(LIS_logunit,*) '[INFO] Here 2 in NoahMP36_getpeobspred_ISMNsmobs '

!  call ESMF_StateGet(Obj_Func,"ISMNsm standard deviation of soil moisture",smstdField,rc=status)
!  call LIS_verify(status)
!
!  call ESMF_FieldGet(smstdField,localDE=0,farrayPtr=smstd,rc=status)
!  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     smc(t) = NoahMP36_struc(n)%noahmp36(t)%smc(1)
!     smstd(t) = NoahMP36_struc(n)%noahmp36(t)%smc_std
  enddo

!  write(LIS_logunit,*) '[INFO] Finished NoahMP36_getpeobspred_ISMNsmobs, smc = ', smc

end subroutine NoahMP36_getpeobspred_ISMNsmobs



