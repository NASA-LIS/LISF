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
module vic412_lsmMod
!BOP
!
! !MODULE: vic412_lsmMod
!
! !DESCRIPTION:
! This module provides the definition of derived data type used to 
! control the operation of the VIC LSM. It also provides the entry method
! for the initialization of VIC-specific variables. The derived
! data type {\tt vic412\_struc} includes the variables that specify
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
  use vic412_module
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: vic412_lsm_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: vic412_struc
!EOP
  type, public :: vic412_type_dec
     integer                    :: vicopen
     integer                    :: numout
     real                       :: rstInterval
     real                       :: ts
     integer                    :: snowstep
     integer                    :: tscount
     integer                    :: veg_tiling_scheme
     integer                    :: NT
     integer                    :: debugging_convert_units
     integer                    :: state_chunk_size
     real*8                     :: MAX_SNOW_TEMP
     character*100              :: global_param
     character(len=LIS_CONST_PATH_LEN) :: rfile 
     character*32               :: rfile_format 
     type(vicdec), allocatable :: vic(:)
  end type vic412_type_dec
  type(vic412_type_dec), allocatable :: vic412_struc(:)

  SAVE
contains
!BOP
! 
! !ROUTINE: vic412_lsm_ini
! \label{vic412_lsm_ini}
! 
! !INTERFACE:
  subroutine vic412_lsm_ini()
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

   allocate(vic412_struc(LIS_rc%nnest))

   call vic412_readcard()

   do n = 1, LIS_rc%nnest
      allocate(vic412_struc(n)%vic(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      ! allocate memory for state chunk of each tile
      !do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
      !  allocate(vic412_struc(n)%vic(t)%state_chunk(vic412_struc(n)%state_chunk_size))
      !enddo

      vic412_struc(n)%numout = 0
      vic412_struc(n)%snowstep = 1
      vic412_struc(n)%tscount = 1

!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
       call LIS_update_timestep(LIS_rc, n, vic412_struc(n)%ts)

       write(fnest,'(i3.3)') n
       call LIS_registerAlarm("VIC412 model alarm "//trim(fnest),&
            vic412_struc(n)%ts,&
            vic412_struc(n)%ts)

       call LIS_registerAlarm("VIC412 restart alarm "//trim(fnest),&
            vic412_struc(n)%ts,&
            vic412_struc(n)%rstInterval)
       
! Added by Shugong Wang 06/18/2014 to resolve parallel run issue
! VIC has 3 layers. This is hard coded. 
       LIS_sfmodel_struc(n)%nsm_layers = 3
       LIS_sfmodel_struc(n)%nst_layers = 3
       allocate(LIS_sfmodel_struc(n)%lyrthk(3))
! Layer thicknesses given below are _ONLY_ for LIS_historyMod,
! especially for specification of layer depths in GRIB output.
! Note that these thicknesses are currently set to a constant
! 1.0cm to represent the first, second, and third soil layers.
! The layer thicknesses in the LIS-VIC GRIB and netCDF files
! are not representative of the actual VIC layer thicknesses,
! which vary spatially.  20 Mar 2015 - DMM
! VIC manages its own layer thicknesses elsewhere.
       do i = 1,3
          LIS_sfmodel_struc(n)%lyrthk(i) = 1.0
       enddo
       LIS_sfmodel_struc(n)%ts = vic412_struc(n)%ts
   enddo
  
  end subroutine vic412_lsm_ini

end module vic412_lsmMod

