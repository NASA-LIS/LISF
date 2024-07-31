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
! !ROUTINE: CRTM2EM_getpeobspred_CNRS
!  \label{CRTM2EM_getpeobspred_CNRS}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine CRTM2EM_getpeobspred_CNRS(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_soilsMod,  only : LIS_soils
  use LIS_logMod,       only : LIS_verify
#if (defined RTMS) 
  use CRTM2_EMMod, only : crtm_struc
#endif

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: Obj_Func
!
! !DESCRIPTION:
!  
! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: emField
!  real, pointer          :: em_data(:,:)
  real, pointer          :: em_data(:) ! just focusing on one channel now
  integer                :: t,col,row
  integer                :: i
  integer                :: status
  integer                :: k
  integer, parameter     :: cnrs_numchannels=7  !19V,19H,22V,37V,37H,85V,85H

  n = 1
  call ESMF_StateGet(Obj_Func,"Emissivity",emField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(emField,localDE=0,farrayPtr=em_data,rc=status)
  call LIS_verify(status)

#if (defined RTMS) 
  do t=1,LIS_rc%ntiles(n)  ! ln is channel
     em_data(t) = crtm_struc(n)%SfcOptics(1, t)%Emissivity(2, 1 ) ! channel is first index of SfcOptics
  enddo
#endif

end subroutine CRTM2EM_getpeobspred_CNRS



