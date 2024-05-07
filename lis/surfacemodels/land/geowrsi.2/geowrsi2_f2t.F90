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
!
! !ROUTINE: geowrsi2_f2t
! \label{geowrsi2_f2t}
!
! !REVISION HISTORY: 
! 31Jul2011: Brad Wind;   Initial Definition
! 28Feb2013: KR Arsenault/JV Geiger; Removed original FEWSNET forcing calls 
! 25Oct2013: KR Arsenault;  Added GeoWRSI2.0 model to LIS-7
! 
! !INTERFACE:
subroutine geowrsi2_f2t(n)

! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_surface
  use LIS_logMod,         only : LIS_verify
  use LIS_metforcingMod,  only : LIS_FORC_State
  use LIS_FORC_AttributesMod 
  use geowrsi2_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
!  This routine transfers the LIS provided forcing onto the WRSI
!  model tiles
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

  integer            :: t, status
  integer            :: tid
  type(ESMF_Field)   :: tpcpField, rETfield
  real,pointer       :: tpcp(:),rET(:)
  real               :: total_precip, total_pet
  integer            :: rate_to_total
! __________________________________________________________

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),tpcpField,&
       rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_PET%varname(1)),rETfield,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(tpcpField,localDE=0,farrayPtr=tpcp,rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(rETfield,localDE=0,farrayPtr=rET,rc=status)
  call LIS_verify(status)

! Assign number seconds per model timestep to local rate:
  rate_to_total = LIS_rc%nts(n)

  do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)

   !transform t to the patch
    tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id

    if(tpcp(tid).ne.LIS_rc%udef) then
     ! Convert rate to total precipitation (mm) applicable for this LIS time-step
       total_precip = tpcp(tid) * rate_to_total
    else
       total_precip = 0.0
    endif
    geowrsi2_struc(n)%wrsi(t)%ACC_PPT = geowrsi2_struc(n)%wrsi(t)%ACC_PPT + &
                                    total_precip

    if( rET(tid).ne.LIS_rc%udef ) then
     ! Convert rate to total potential evapotranspiration (mm) applicable 
     ! for this LIS time-step
       total_pet = rET(tid) * rate_to_total
    else
       total_pet = 0.0
    endif
    geowrsi2_struc(n)%wrsi(t)%ACC_PET = geowrsi2_struc(n)%wrsi(t)%ACC_PET + &
                                    total_pet

  enddo

end subroutine geowrsi2_f2t

