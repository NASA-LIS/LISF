!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
module FLake1_Mod
!BOP
!
! !MODULE: FLake1_Mod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of FLake1 model. It also provides the entry method
! for the initialization of FLake1-specific variables. The derived
! data type {\tt FLake1\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the FLake1 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for FLake1
! \item[ts]
!   FLake1 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[flake1]
!  FLake1 model specific variables
! \item[forc\_count]
!   counter of forcing data ??????
! \item[height\_wind]
!   height where wind speed is measured
! \item[height\_tq]
!   height where temperature and humidity are measured
! \item[flake\_dt]
!   model time step
! \item[lon]
!   longitude of lake center
! \item[lat]
!   latitude of lake center
! \item[depth\_w]
!   lake depth
! \item[fetch]
!   typical wind fetch
! \item[depth\_bs]
!   depth of the thermally active layer of the bottom sediments
! \item[T\_bs]
!   temperature at the outer edge of the thermally active layer of the bottom sediments
! \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  6/4/13: Shugong Wang Initial implementation for LIS 7 and FLake1
!
! !USES:
  use FLake1_module
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none
  
  PRIVATE
!-------------------------------------------------------------------------
! PUBLIC MEMBER FUNCTIONS
!-------------------------------------------------------------------------
  public :: FLake1_ini
!-------------------------------------------------------------------------
! PUBLIC TYPES
!-------------------------------------------------------------------------
  public :: FLake1_struc
!EOP
  type, public :: FLake1_type_dec
     logical            :: startFlag
     character(len=LIS_CONST_PATH_LEN) :: rfile
     character*256      :: rformat
!-------------------------------------------------------------------------
! ts, Count, rstInterval, outInterval
!-------------------------------------------------------------------------
     real               :: ts
     integer            :: count
     real               :: rstInterval
     real               :: outInterval
     integer            :: forc_count
!-------------------------------------------------------------------------
! Model State
!-------------------------------------------------------------------------
     real               :: init_T_snow
     real               :: init_T_ice
     real               :: init_T_mnw
     real               :: init_T_wML
     real               :: init_T_bot
     real               :: init_T_b1
     real               :: init_C_T
     real               :: init_H_snow
     real               :: init_H_ice
     real               :: init_H_ML
     real               :: init_H_B1
     real               :: init_T_sfc
     real               :: init_albedo_water
     real               :: init_albedo_ice
     real               :: init_albedo_snow
!-------------------------------------------------------------------------
! Uniform Parameter
!-------------------------------------------------------------------------
     real               :: height_wind
     real               :: height_tq
     real               :: flake_dt
!-------------------------------------------------------------------------
! Spatial Parameter
!-------------------------------------------------------------------------
     real, allocatable      :: lon(:,:)
     real, allocatable      :: lat(:,:)
     real, allocatable      :: depth_w(:,:)
     real, allocatable      :: fetch(:,:)
     real, allocatable      :: depth_bs(:,:)
     real, allocatable      :: T_bs(:,:)
     type(FLake1dec), allocatable :: flake1(:)
  end type FLake1_type_dec
  
  type(FLake1_type_dec), allocatable :: FLAKE1_struc(:)
  
contains 

!BOP
!
! !ROUTINE: FLake1_ini
! \label{FLake1_ini}
!
! !INTERFACE:
  subroutine FLake1_ini()
! !USES:
    use LIS_coreMod, only             : LIS_rc
    use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
    use LIS_timeMgrMod, only          : LIS_clock, LIS_calendar, &
         LIS_update_timestep,     &
         LIS_registerAlarm
    use LIS_logMod, only              : LIS_verify
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for FLake1-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for FLake1 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[FLake1\_readcrd](\ref{FLake1_readcrd}) \newline
!    reads the runtime options for FLake1 model
!  \end{description}
!EOP
    implicit none
    integer                 :: i, n,t
    integer                 :: status
    character*3             :: fnest
    
    ! allocate memory for nest 
    allocate(FLAKE1_struc(LIS_rc%nnest))
    
    ! read configuation information from lis.config file
    call FLake1_readcrd()
    
    do n=1, LIS_rc%nnest
       ! allocate memory for all tiles in current nest 
       allocate(FLAKE1_struc(n)%flake1(LIS_rc%npatch(n, LIS_rc%lake_index)))
!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
       FLAKE1_struc(n)%forc_count = 0
       FLAKE1_struc(n)%startFlag = .true. 
       
       call LIS_update_timestep(LIS_rc, n, FLAKE1_struc(n)%ts)
  
       write(fnest,'(i3.3)') n
       call LIS_registerAlarm("FLAKE1 model alarm"//trim(fnest),&
            FLAKE1_struc(n)%ts, &
            FLAKE1_struc(n)%ts)
       
       call LIS_registerAlarm("FLAKE1 restart alarm"//trim(fnest), &
            FLAKE1_struc(n)%ts,&
            FLAKE1_struc(n)%rstInterval)
!------------------------------------------------------------------------
! TODO: setup number of soil moisture/temperature layers and depth here  
!------------------------------------------------------------------------
       ! TODO: set number of soil moisture layers in surface model
       LIS_sfmodel_struc(n)%nsm_layers = 0
       ! TODO: set number of soil temperature layers in surface model
       LIS_sfmodel_struc(n)%nst_layers = 0
       
       if(.not.allocated(LIS_sfmodel_struc(n)%lyrthk)) then 
          allocate(LIS_sfmodel_struc(n)%lyrthk(0))
       endif

       LIS_sfmodel_struc(n)%ts = FLAKE1_struc(n)%ts
       
       do t=1, LIS_rc%npatch(n, LIS_rc%lake_index)
          FLAKE1_struc(n)%flake1(t)%tair =0 
          FLAKE1_struc(n)%flake1(t)%qair =0
          FLAKE1_struc(n)%flake1(t)%swdown =0
          FLAKE1_struc(n)%flake1(t)%lwdown = 0
          FLAKE1_struc(n)%flake1(t)%wind_e =0 
          FLAKE1_struc(n)%flake1(t)%wind_n = 0
          FLAKE1_struc(n)%flake1(t)%psurf =0
       enddo
    enddo
  end subroutine FLake1_ini
end module FLake1_Mod
