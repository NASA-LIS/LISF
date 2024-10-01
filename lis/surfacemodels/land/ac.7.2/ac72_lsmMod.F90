!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
module AC72_lsmMod
!BOP
!
! !MODULE: AC72_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of AC72 model. It also provides the entry method
! for the initialization of AC72-specific variables. The derived
! data type {\tt AC72\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the AC72 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for AC72
! \item[LDT\_ncvar\_soiltype]
!   LDT NetCDF variable name for soil type index
! \item[ts]
!   AC72 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[ac72]
!  AC72 model specific variables
! \item[forc\_count]
!   counter of forcing data
! \item[soil\_tbl\_name]
!   AquaCrop model soil parameter table
! \item[soil\_scheme\_name]
!   soil classification scheme
! \item[dt]
!   time step in seconds
! \item[PathNameSimul]
!   path of simulation folder
! \item[PathCropFiles]
!   path of crop files folder
! \item[CO2_Filename]
!   CO2 filename (should be in simul folder)
! \item[Management_Filename]
!   Management filename (should be in simul folder)
! \item[Irrigation_Filename]
!   Irrigation filename (should be in simul folder)
! \item[Sim_AnnualStartDay]
!   annual start day of simulation period
! \item[Sim_AnnualStartMonth]
!   annual start month of simulation period
! \item[Crop_AnnualStartDay]
!   annual start day of cropping period
! \item[Crop_AnnualStartMonth]
!   annual start month of cropping period
! \item[NrSoilLayers]
!   number of soil layers !!different from compartments!!
! \item[max_No_compartments]
!   maximum number of compartments (=12)
! \item[Thickness]
!   thickness of soil layers
! \item[refz_tq]
!   reference height of forcings T and q
! \item[refz_uv]
!   reference height of forcings u and v

!
! !REVISION HISTORY:
!  18 JAN 2024, Louise Busschaert; initial implementation for AC72
!
! !USES:
    use AC72_module

    implicit none

    PRIVATE
    !-------------------------------------------------------------------------
    ! PUBLIC MEMBER FUNCTIONS
    !-------------------------------------------------------------------------
    public :: AC72_ini
    !-------------------------------------------------------------------------
    ! PUBLIC TYPES
    !-------------------------------------------------------------------------
    public :: AC72_struc
!EOP
    type, public :: AC72_type_dec
        character*256      :: rfile
        character*256      :: rformat
        !-------------------------------------------------------------------------
        ! Parameter file names
        !-------------------------------------------------------------------------
        character*128      :: LDT_ncvar_soiltype
        !-------------------------------------------------------------------------
        ! ts, Count, rstInterval, outInterval
        !-------------------------------------------------------------------------
        real               :: ts
        integer            :: count
        real               :: rstInterval
        integer            :: outInterval
        integer            :: forc_count
        !-------------------------------------------------------------------------
        ! Initial Model State for cold start
        !-------------------------------------------------------------------------
        real, pointer      :: init_smc(:)

        !-------------------------------------------------------------------------
        ! Constant Parameter
        !-------------------------------------------------------------------------
        character(len=256) :: soil_tbl_name
        character(len=256) :: soil_scheme_name
        real               :: dt

        character(len=256) :: PathNameOutp
        character(len=256) :: PathNameSimul
        character(len=256) :: PathNameList
        character(len=256) :: PathNameParam
        character(len=256) :: PathCropFiles
        character(len=256) :: CO2_Filename
        character(len=256) :: Management_Filename
        character(len=256) :: Irrigation_Filename
        integer            :: Sim_AnnualStartDay
        integer            :: Sim_AnnualStartMonth
        integer            :: Crop_AnnualStartDay
        integer            :: Crop_AnnualStartMonth
        integer            :: NrSoilLayers
        integer            :: max_No_compartments
        real, pointer      :: Thickness(:)
        real               :: refz_tq
        real               :: refz_uv
        type(AC72dec), pointer :: ac72(:)
    end type AC72_type_dec

    type(AC72_type_dec), pointer :: AC72_struc(:)
 
contains 

!BOP
!
! !ROUTINE: AC72_ini
! \label{AC72_ini}
!
! !INTERFACE:
    subroutine AC72_ini()
! !USES:
        use LIS_coreMod, only : LIS_rc
        use LIS_logMod, only : LIS_verify
        use LIS_timeMgrMod, only : LIS_clock,  LIS_calendar, &
            LIS_update_timestep, LIS_registerAlarm
        use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for AC72-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for AC72 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[AC72\_readcrd](\ref{AC72_readcrd}) \newline
!    reads the runtime options for AC72 model
!  \end{description}
!EOP
        implicit none        
        integer  :: n, t, i    
        integer  :: status   

        ! allocate memory for nest 
        allocate(AC72_struc(LIS_rc%nnest))
 
        ! read configuation information from lis.config file
        call AC72_readcrd()

        do n=1, LIS_rc%nnest
            ! allocate memory for all tiles in current nest 
            allocate(AC72_struc(n)%ac72(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            !------------------------------------------------------------------------
            ! allocate memory for vector variables passed to model interfaces        
            ! TODO: check the following allocation statements carefully!
            !------------------------------------------------------------------------
            AC72_struc(n)%max_No_compartments = 12 ! hard coded
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(AC72_struc(n)%ac72(t)%smc(AC72_struc(n)%max_No_compartments))
            enddo

            ! allocate memory for Trecord arrays
            do t=1,LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(AC72_struc(n)%ac72(t)%Tmax_record(366))
                allocate(AC72_struc(n)%ac72(t)%Tmin_record(366))
            enddo
            
            ! initialize forcing variables to zeros
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC72_struc(n)%ac72(t)%tair = 0.0
                AC72_struc(n)%ac72(t)%tmax = 0.0
                AC72_struc(n)%ac72(t)%tmin = 0.0
                AC72_struc(n)%ac72(t)%tdew = 0.0
                AC72_struc(n)%ac72(t)%wndspd = 0.0            
                AC72_struc(n)%ac72(t)%psurf = 0.0
                AC72_struc(n)%ac72(t)%prcp = 0.0
                AC72_struc(n)%ac72(t)%eto = 0.0

                !LB: Initialize HarvestNow (new in AC7.1)
                AC72_struc(n)%ac72(t)%HarvestNow = .false.
                
            enddo ! end of tile (t) loop
!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
            AC72_struc(n)%forc_count = 0

            call LIS_update_timestep(LIS_rc, n, AC72_struc(n)%ts)

            call LIS_registerAlarm("AC72 model alarm",&
                                   AC72_struc(n)%ts, &
                                   AC72_struc(n)%ts)

            call LIS_registerAlarm("AC72 restart alarm", &
                                   AC72_struc(n)%ts,&
                                   AC72_struc(n)%rstInterval)
            ! Set number of soil moisture layers in surface model
            LIS_sfmodel_struc(n)%nsm_layers = AC72_struc(n)%max_No_compartments
            allocate(LIS_sfmodel_struc(n)%lyrthk(AC72_struc(n)%max_No_compartments))
            ! Note, the default compartment size (0.10) is stored in the 
            ! SURFACEMODEL output file, but these are spatially variable and included in the
            ! input file produced by LDT
            LIS_sfmodel_struc(n)%lyrthk(:) = 0.10
            LIS_sfmodel_struc(n)%ts = AC72_struc(n)%ts
        enddo
    end subroutine AC72_ini
end module AC72_lsmMod
