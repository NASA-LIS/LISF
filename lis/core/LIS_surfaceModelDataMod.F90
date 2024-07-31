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
module LIS_surfaceModelDataMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: LIS_sfmodel_struc

  type, public :: sfmodeldec
     real          :: outInterval 
     character*50  :: outIntervalType
     integer       :: nsm_layers
     integer       :: nst_layers
     real          :: ts
     logical       :: stats_file_open
     real, allocatable :: lyrthk(:)
     real, allocatable :: lyrthk2(:)  !for CLSM nsm_layers/=nst_layers
     character*100 :: models_used
  end type sfmodeldec

  type(sfmodeldec), allocatable :: LIS_sfmodel_struc(:)

end module LIS_surfaceModelDataMod
