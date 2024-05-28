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
module Ac71_lsmMod
!BOP
!
! !MODULE: Ac71_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of Ac71 model. It also provides the entry method
! for the initialization of Ac71-specific variables. The derived
! data type {\tt Ac71\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the Ac71 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for Ac71
! \item[LDT\_ncvar\_soiltype]
!   LDT NetCDF variable name for soil type index
! \item[ts]
!   Ac71 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[ac71]
!  Ac71 model specific variables
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
! \item[Sim_AnnualEndDay]
!   annual end day of simulation period
! \item[Sim_AnnualEndMonth]
!   annual end month of simulation period
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
! \item[refz_forc]
!   reference height of forcings

!
! !REVISION HISTORY:
!  18 JAN 2024, Louise Busschaert; initial implementation for AC71
!
! !USES:
    use Ac71_module

    implicit none

    PRIVATE
    !-------------------------------------------------------------------------
    ! PUBLIC MEMBER FUNCTIONS
    !-------------------------------------------------------------------------
    public :: Ac71_ini
    !-------------------------------------------------------------------------
    ! PUBLIC TYPES
    !-------------------------------------------------------------------------
    public :: Ac71_struc
!EOP
    type, public :: Ac71_type_dec
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

        character(len=256) :: PathNameSimul
        character(len=256) :: PathCropFiles
        character(len=256) :: CO2_Filename
        character(len=256) :: Management_Filename
        character(len=256) :: Irrigation_Filename
        integer            :: Sim_AnnualStartDay
        integer            :: Sim_AnnualEndDay
        integer            :: Sim_AnnualStartMonth
        integer            :: Sim_AnnualEndMonth
        integer            :: Crop_AnnualStartDay
        integer            :: Crop_AnnualStartMonth
        integer            :: GDD_Mode
        integer            :: NrSoilLayers
        integer            :: max_No_compartments
        real, pointer      :: Thickness(:)
        real               :: refz_forc
        type(Ac71dec), pointer :: ac71(:)
        type(ac71_trecord), pointer :: Trecord(:)
    end type Ac71_type_dec

    type(Ac71_type_dec), pointer :: AC71_struc(:)
 
contains 

!BOP
!
! !ROUTINE: Ac71_ini
! \label{Ac71_ini}
!
! !INTERFACE:
    subroutine Ac71_ini()
! !USES:
        use LIS_coreMod, only : LIS_rc
        use LIS_logMod, only : LIS_verify
        use LIS_timeMgrMod, only : LIS_clock,  LIS_calendar, &
            LIS_update_timestep, LIS_registerAlarm
        use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for Ac71-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for Ac71 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[Ac71\_readcrd](\ref{Ac71_readcrd}) \newline
!    reads the runtime options for Ac71 model
!  \end{description}
!EOP
        implicit none        
        integer  :: n, t, p    
        integer  :: status   

        ! allocate memory for nest 
        allocate(AC71_struc(LIS_rc%nnest))
 
        ! read configuation information from lis.config file
        call Ac71_readcrd()

        do n=1, LIS_rc%nnest
            ! allocate memory for all tiles in current nest 
            allocate(AC71_struc(n)%ac71(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            !------------------------------------------------------------------------
            ! allocate memory for vector variables passed to model interfaces        
            ! TODO: check the following allocation statements carefully!
            !------------------------------------------------------------------------
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(AC71_struc(n)%ac71(t)%smc(AC71_struc(n)%max_No_compartments))
            enddo

            ! allocate memory for Trecord arrays
            allocate(AC71_struc(n)%Trecord(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
            do p=1,LIS_rc%gnc(n)*LIS_rc%gnr(n)
                allocate(AC71_struc(n)%Trecord(p)%Tmax_record(366))
                allocate(AC71_struc(n)%Trecord(p)%Tmin_record(366))
            enddo

            
            ! initialize forcing variables to zeros
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%tair = 0.0
                AC71_struc(n)%ac71(t)%tmax = 0.0
                AC71_struc(n)%ac71(t)%tmin = 0.0
                AC71_struc(n)%ac71(t)%tdew = 0.0
                AC71_struc(n)%ac71(t)%wndspd = 0.0            
                AC71_struc(n)%ac71(t)%psurf = 0.0
                AC71_struc(n)%ac71(t)%prcp = 0.0
                AC71_struc(n)%ac71(t)%eto = 0.0

                !LB: Initialize HarvestNow (new in AC7.1)
                AC71_struc(n)%ac71(t)%HarvestNow = .false.
                
            enddo ! end of tile (t) loop
!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
            AC71_struc(n)%forc_count = 0

            call LIS_update_timestep(LIS_rc, n, AC71_struc(n)%ts)

            call LIS_registerAlarm("Ac71 model alarm",&
                                   AC71_struc(n)%ts, &
                                   AC71_struc(n)%ts)

            call LIS_registerAlarm("Ac71 restart alarm", &
                                   AC71_struc(n)%ts,&
                                   AC71_struc(n)%rstInterval)
            ! Set number of soil moisture layers in surface model
            LIS_sfmodel_struc(n)%nsm_layers = AC71_struc(n)%max_No_compartments
            allocate(LIS_sfmodel_struc(n)%lyrthk(AC71_struc(n)%max_No_compartments))
            LIS_sfmodel_struc(n)%lyrthk(:) = 0.1
            LIS_sfmodel_struc(n)%ts = AC71_struc(n)%ts
        enddo
    end subroutine Ac71_ini
end module Ac71_lsmMod
