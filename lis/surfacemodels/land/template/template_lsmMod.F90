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
module template_lsmMod
!BOP
!
! !MODULE: template_lsmMod
!
! !DESCRIPTION:
!  Module for 1-D land model driver variable initialization
!  
! \begin{description}
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[numout]
!    number of output times 
!  \item[outInterval]
!    output writing interval
!  \item[templateopen]
!    variable to keep track of opened files
!  \item[template]
!   Template LSM specific variables
! \end{description} 
!
! !REVISION HISTORY:
!    Apr 2003; Sujay Kumar, Initial Code
! 23 Oct 2007; Kristi Arsenault, Updated for V5.0
!
! !USES:        

  use template_module

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: template_lsm_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: template_struc
!EOP
  type, public :: template_type_dec
     integer                    :: templateopen
     integer                    :: numout
     real                       :: ts
     type(templatedec), allocatable :: template(:)
  end type template_type_dec
  type(template_type_dec), allocatable :: template_struc(:)

  SAVE
contains
!BOP
! 
! !ROUTINE: template_lsm_ini
! \label{template_lsm_ini}
! 
! !INTERFACE:
  subroutine template_lsm_ini()
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

   allocate(template_struc(LIS_rc%nnest))

   call template_readcrd()
   do n = 1, LIS_rc%nnest
      allocate(template_struc(n)%template(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      template_struc(n)%numout = 0

      call LIS_update_timestep(LIS_rc, n, template_struc(n)%ts)

      LIS_sfmodel_struc(n)%nsm_layers = 1
      LIS_sfmodel_struc(n)%nst_layers = 1
      allocate(LIS_sfmodel_struc(n)%lyrthk(1))
      LIS_sfmodel_struc(n)%lyrthk(1) = 1
      LIS_sfmodel_struc(n)%ts = template_struc(n)%ts
   enddo
  
  end subroutine template_lsm_ini

end module template_lsmMod

