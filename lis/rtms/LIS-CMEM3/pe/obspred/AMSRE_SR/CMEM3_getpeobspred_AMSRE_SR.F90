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
! !ROUTINE: CMEM3_getpeobspred_AMSRE_SR
!  \label{CMEM3_getpeobspred_AMSRE_SR}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine CMEM3_getpeobspred_AMSRE_SR(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_soilsMod,  only : LIS_soils
  use LIS_logMod,       only : LIS_verify
#if (defined RTMS) 
  use CMEM3_Mod  , only : CMEM3_struc
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
  real, pointer          :: em_data(:) ! just focusing on one channel now
  integer                :: t,col,row
  integer                :: i
  integer                :: status
  integer                :: k
  real :: em_v
  real :: em_h
  character*2, parameter :: fnames(6)=(/'07', '11', '19', '24', '37', '89'/)
  character*1, parameter :: pnames(2)=(/'V','H'/)
  integer                   ::  f !freq counter
  integer                   ::  p !polarization counter
  integer                   :: numfreqs=6        
  integer                   :: numpolarizations=2

  n = 1

  do f=1, numfreqs
     do p=1, numpolarizations 
        call ESMF_StateGet(Obj_Func,"Emissivity" // fnames(f) // pnames(p),emField,rc=status)
        call LIS_verify(status)
        
        call ESMF_FieldGet(emField,localDE=0,farrayPtr=em_data,rc=status)
        call LIS_verify(status)

#if (defined RTMS) 
        do t=1,LIS_rc%ntiles(n)  ! ln is channel
           if(pnames(p).eq.'H') then
              em_data(t) = cmem3_struc(n)%emH(t,f)
           else !V
              em_data(t) = cmem3_struc(n)%emV(t,f)
           endif
        enddo
#endif
     enddo
  enddo

end subroutine CMEM3_getpeobspred_AMSRE_SR
