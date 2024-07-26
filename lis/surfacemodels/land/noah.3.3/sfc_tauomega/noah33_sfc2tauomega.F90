!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: noah33_sfc2tauomega
! \label{noah33_sfc2tauomega}
!
! !REVISION HISTORY:
!  28 Aug 2012: Sujay Kumar; Initial Code
!
! !INTERFACE:
subroutine noah33_sfc2tauomega(n, sfcState)
! !USES:      
  use ESMF
  use LIS_coreMod
  use LIS_logMod,    only : LIS_verify
  use LIS_constantsMod,  only : LIS_CONST_RHOFW

  use noah33_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  type(ESMF_State)    :: sfcState
  
! FUNCTIONS

! 
! !DESCRIPTION: 
! This subroutine assigns the noah33 specific surface variables
! to TAUOMEGA. 
!
!EOP
  type(ESMF_Field)    :: smcField, stcField
  real, pointer       :: soil_moisture_content(:), soil_temperature(:)
  integer             :: t,status

  call ESMF_StateGet(sfcState,"Soil Moisture Content",smcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=soil_moisture_content,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Soil Temperature",stcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(stcField,localDE=0,farrayPtr=soil_temperature, rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soil_moisture_content(t) = noah33_struc(n)%noah(t)%smc(1)
     soil_temperature(t) = noah33_struc(n)%noah(t)%stc(1)
  enddo

end subroutine noah33_sfc2tauomega


