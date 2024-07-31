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
! !ROUTINE: CRTM2EM_getpeobspred_CNRS_MPDI
!  \label{CRTM2EM_getpeobspred_CNRS_MPDI}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine CRTM2EM_getpeobspred_CNRS_MPDI(Obj_Func)
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
  real :: em_v
  real :: em_h

  n = 1
  call ESMF_StateGet(Obj_Func,"Emissivity",emField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(emField,localDE=0,farrayPtr=em_data,rc=status)
  call LIS_verify(status)

#if (defined RTMS) 
!   allocate(em_data(LIS_rc%ntiles(n),cnrs_numchannels))
!   allocate(em_data(LIS_rc%ntiles(n)))

!   do t=1,LIS_rc%ntiles(n)
!      em_data(t,1) = crtm_struc(n)%emV(t,2)  
!      em_data(t,2) = crtm_struc(n)%emV(t,1)  
!      em_data(t,3) = crtm_struc(n)%emV(t,4)  
!      em_data(t,4) = crtm_struc(n)%emV(t,6)  
!      em_data(t,5) = crtm_struc(n)%emV(t,5)  
!      em_data(t,6) = crtm_struc(n)%emV(t,8)  
!      em_data(t,7) = crtm_struc(n)%emV(t,7)  
!   enddo
  do t=1,LIS_rc%ntiles(n)
     em_v= crtm_struc(n)%SfcOptics(1, t)%Emissivity(1, 1 )
     em_h= crtm_struc(n)%SfcOptics(2, t)%Emissivity(1, 1 )
     em_data(t) = (em_v-em_h)/(em_v+em_h)
!      em_data(t,2) = crtm_struc(n)%emV(t,1)  
!      em_data(t,3) = crtm_struc(n)%emV(t,4)  
!      em_data(t,4) = crtm_struc(n)%emV(t,6)  
!      em_data(t,5) = crtm_struc(n)%emV(t,5)  
!      em_data(t,6) = crtm_struc(n)%emV(t,8)  
!      em_data(t,7) = crtm_struc(n)%emV(t,7)  
  enddo
#endif

end subroutine CRTM2EM_getpeobspred_CNRS_MPDI



