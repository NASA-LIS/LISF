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
! !ROUTINE: ac71_getirrigationstates
! \label{ac71_getirrigationstates}
! 
! !INTERFACE:
subroutine ac71_getirrigationstates(n,irrigState)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use ac71_lsmMod


! !DESCRIPTION:

! Gets the irrigation from the AquaCrop structure and transfers
! it to the irrigation structure to be written out as output.
! All irrigaiton calculations are performed in AquaCrop,
! the obtained values in mm/day are converted to kg m-2 s-1
! to remain consistent with other irrigation outputs.
! 
! 
!
! REVISION HISTORY:
!
! 18 Apr 2024: Louise Busschaert; Initial implementation
!
!EOP
  implicit none

  integer              :: n,rc
  integer              :: t
  type(ESMF_State)     :: irrigState
  type(ESMF_Field)     :: irrigRateField

  real,  pointer       :: irrigRate(:)
  real				         :: gthresh

  call ESMF_StateGet(irrigState, "Irrigation rate",irrigRateField,rc=rc)
  call LIS_verify(rc,'ESMF_StateGet failed for Irrigation rate')    
  call ESMF_FieldGet(irrigRateField, localDE=0,farrayPtr=irrigRate,rc=rc)
  call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation rate')


!----------------------------------------------------------------------
! Extract irrigation from LSM and convert to the right units
!----------------------------------------------------------------------

    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      ! Extract irrigation from ac71 structure
        irrigRate(t) = AC71_struc(n)%ac71(t)%Irrigation/86400

        !---------------------------------------------------------------------
	      ! For next time step: dynamic irrigation (if option turned on)
	      !----------------------------------------------------------------------
      if (LIS_rc%irrigation_dveg .eq. 1) then
        ! Calculate threshold
        gthresh = AC71_struc(n)%ac71(t)%Crop%CCini &
          + (LIS_rc%irrigation_GVFparam1 + LIS_rc%irrigation_GVFparam2*&
          (AC71_struc(n)%ac71(t)%Crop%CCx - AC71_struc(n)%ac71(t)%Crop%CCini)) &
          * (AC71_struc(n)%ac71(t)%Crop%CCx - AC71_struc(n)%ac71(t)%Crop%CCini)
        if (AC71_struc(n)%ac71(t)%CCiActual .ge. gthresh) then
            ! Irrigation allowed
            AC71_struc(n)%ac71(t)%IrriInfoRecord1%TimeInfo = &
                  int(LIS_rc%irrigation_thresh)
        else ! very large threshold to block irrigation
            AC71_struc(n)%ac71(t)%IrriInfoRecord1%TimeInfo = 400 
        endif
      else ! irrigation always allowed
        AC71_struc(n)%ac71(t)%IrriInfoRecord1%TimeInfo = &
              int(LIS_rc%irrigation_thresh)
      endif
    end do

  end subroutine ac71_getirrigationstates