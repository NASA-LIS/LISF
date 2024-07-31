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
module vic411_lsmMod
!BOP
!
! !MODULE: vic411_lsmMod
!
! !DESCRIPTION:
! This module provides the definition of derived data type used to 
! control the operation of the VIC LSM. It also provides the entry method
! for the initialization of VIC-specific variables. The derived
! data type {\tt vic411\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!  
! \begin{description}
!  \item[vicopen]
!    variable to keep track of opened files
!  \item[numout]
!    number of output times 
!  \item[rstInterval]
!    restart writing interval
!  \item[ts]
!    model time-step
!  \item[snowstep]
!    snow-step for water-balance mode
!  \item[tscount]
!    time-step count, used to know which snow-step is being processed when
!    running in water-balance mode.
!  \item[veg\_tiling\_scheme]
!    specifies whether sub-grid tiling is determined by LIS or by VIC. \newline
!       0 = VIC \newline
!       1 = LIS
!  \item[NT] total number of vegetation types for the vegetation parameter file
!  \item[debugging\_convert\_units]
!    debugging flag -- remove
!    specifies whether to convert the units of temperature from Celcius
!    to Kelvin.         \newline
!    0 = do not convert \newline
!    1 = convert        \newline
!    Final version of forcing reader will always perform this conversion.
!  \item[MAX\_SNOW\_TEMP]
!    maximum snow temperature
!  \item[global\_param]
!    VIC's global parameter file/configuration file
!    When running in energy-balance mode, VIC should run every time-step.
!    However, when running in water-balance mode, VIC runs daily, but it
!    expects to be given an array of forcing values (sized by the number
!    of snow-steps).  This alarm is necessary to allow LIS to run at a
!    sub-daily time-step to gather all the forcing needed to run VIC daily.
!  \item[vic]
!   VIC LSM specific variables
! \end{description} 
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 14 Aug 2013; Shugong Wang, Implementation of VIC 4.1.1 into LIS-7. 
! !USES:        
  use vic411_module

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: vic411_lsm_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: vic411_struc
!EOP
  type, public :: vic411_type_dec
     integer                    :: vicopen
     integer                    :: numout
     real                       :: rstInterval
     real                       :: ts
     integer                    :: snowstep
     integer                    :: tscount
     integer                    :: veg_tiling_scheme
     integer                    :: NT
     integer                    :: debugging_convert_units
     real*8                     :: MAX_SNOW_TEMP
     character*100              :: global_param
     type(vicdec), allocatable :: vic(:)
  end type vic411_type_dec
  type(vic411_type_dec), allocatable :: vic411_struc(:)

  SAVE
contains
!BOP
! 
! !ROUTINE: vic411_lsm_ini
! \label{vic411_lsm_ini}
! 
! !INTERFACE:
  subroutine vic411_lsm_ini()
! !USES:
   use LIS_coreMod, only : LIS_rc
   use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
                                LIS_update_timestep, LIS_registerAlarm
   use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
   use LIS_logMod,       only : LIS_verify
! !DESCRIPTION:        
!
!EOP
   implicit none
   integer :: n, t, i
   integer :: status
   character*3 :: fnest

   allocate(vic411_struc(LIS_rc%nnest))

   call vic411_readcard()

   do n = 1, LIS_rc%nnest
      allocate(vic411_struc(n)%vic(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      vic411_struc(n)%numout = 0
      vic411_struc(n)%snowstep = 1
      vic411_struc(n)%tscount = 1

!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
       call LIS_update_timestep(LIS_rc, n, vic411_struc(n)%ts)

       write(fnest,'(i3.3)') n
  
       call LIS_registerAlarm("VIC411 model alarm "//trim(fnest), &
            vic411_struc(n)%ts,&
            vic411_struc(n)%ts)

       call LIS_registerAlarm("VIC411 restart alarm "//trim(fnest), &
            vic411_struc(n)%ts,&
            vic411_struc(n)%rstInterval)

       ! Added by Shugong Wang 07/28/2014 to resolve parallel run
       ! issue
       ! VIC has 3 layers. This is hard coded. 
       LIS_sfmodel_struc(n)%nsm_layers = 3
       LIS_sfmodel_struc(n)%nst_layers = 3
       allocate(LIS_sfmodel_struc(n)%lyrthk(3))
       do i = 1,3
          !!! 100 is a dummy. VIC manages its own layer thickness 
          LIS_sfmodel_struc(n)%lyrthk(i) = 100 
       enddo
       LIS_sfmodel_struc(n)%ts = vic411_struc(n)%ts
   enddo
  
  end subroutine vic411_lsm_ini

end module vic411_lsmMod

