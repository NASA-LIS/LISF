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
module templateOpenWaterMod
!BOP
!
! !MODULE: templateOpenWaterMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! !USES:        

  use templateOpenWater_module

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: template_openwater_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: templateopenwater_struc
!EOP
  type, public :: templateopenwater_type_dec
     real                       :: ts     
     type(templateopenwaterdec), allocatable :: templateopenwater(:)
  end type templateopenwater_type_dec

  type(templateopenwater_type_dec), allocatable :: templateopenwater_struc(:)

  SAVE
contains
!BOP
! 
! !ROUTINE: template_openwater_ini
! \label{template_openwater_ini}
! 
! !INTERFACE:
  subroutine template_openwater_ini()
! !USES:
   use ESMF
   use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
   use LIS_coreMod, only : LIS_rc
   use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
        LIS_update_timestep, LIS_registerAlarm
   use LIS_logMod,       only : LIS_verify
! !DESCRIPTION:        
!
!EOP
   implicit none
   integer :: n
   integer                 :: yr, mo, da, hr, mn, ss
   integer                 :: status

   allocate(templateopenwater_struc(LIS_rc%nnest))

   call readtemplateopenwatercrd()
   do n = 1, LIS_rc%nnest
      allocate(templateopenwater_struc(n)%templateopenwater(LIS_rc%npatch(n,&
           LIS_rc%openwater_index)))

      call LIS_update_timestep(LIS_rc, n, templateopenwater_struc(n)%ts)

!      LIS_sfmodel_struc(n)%nsm_layers = 1
!      LIS_sfmodel_struc(n)%nst_layers = 1
!      allocate(LIS_sfmodel_struc(n)%lyrthk(1))
!      LIS_sfmodel_struc(n)%lyrthk(1) = 1
      LIS_sfmodel_struc(n)%ts = templateopenwater_struc(n)%ts
   enddo
  
 end subroutine template_openwater_ini

end module templateOpenWaterMod

