!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
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
! \item[LDT\_ncvar\_vegetype]
!   LDT NetCDF variable name for land cover type index
! \item[LDT\_ncvar\_soiltype]
!   LDT NetCDF variable name for soil type index
! \item[LDT\_ncvar\_slopetype]
!   LDT NetCDF variable name for slope type for Noah baseflow
! \item[LDT\_ncvar\_tbot]
!   LDT NetCDF variable name for deep-layer soil temperature
! \item[LDT\_ncvar\_pblh]
!   LDT NetCDF variable name for planetary boundary layer height
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
! \item[landuse\_tbl\_name]
!   Noah model landuse parameter table
! \item[soil\_tbl\_name]
!   Noah model soil parameter table
! \item[gen\_tbl\_name]
!   Noah model general parameter table
! \item[ac\_tbl\_name]
!   Ac parameter table
! \item[landuse\_scheme\_name]
!   landuse classification scheme
! \item[soil\_scheme\_name]
!   soil classification scheme
! \item[dveg\_opt]
!   vegetation model
! \item[crs\_opt]
!   canopy stomatal resistance
! \item[btr\_opt]
!   soil moisture factor for stomatal resistance
! \item[run\_opt]
!   runoff and groundwater
! \item[sfc\_opt]
!   surface layer drag coefficients (CH \& CM)
! \item[frz\_opt]
!   supercooled liquid water
! \item[inf\_opt]
!   frozen soil permeability
! \item[rad\_opt]
!   radiation transfer
! \item[alb\_opt]
!   snow surface albedo
! \item[snf\_opt]
!   rainfall \& snowfall
! \item[tbot\_opt]
!   lower boundary of soil temperature
! \item[stc\_opt]
!   snow/soil temperature time scheme
! \item[nslcats]
!   the number of total soil types in parameter table
! \item[nlucats]
!   the number of total land cover types in parameter table
! \item[nslpcats]
!   the number of total slope category for Noah baseflow
! \item[dt]
!   time step in seconds
! \item[nsoil]
!   number of soil layers
! \item[sldpth]
!   thickness of soil layers
! \item[nsnow]
!   maximum number of snow layers
! \item[urban\_vegetype]
!   urban land cover type index
! \item[ice\_flag]
!   ice flag: 0 = no ice, 1 = ice
! \item[st\_flag]
!   surface type 1=soil, 2=lake
! \item[sc\_idx]
!   soil color type
! \item[iz0tlnd]
!   option of Chen adjustment of Czil
! \item[zlvl]
!   reference height of temperature and humidity
! \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  18 JAN 2024, Louise Busschaert; initial implementation for LIS 7 and AC71
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
        character*128      :: LDT_ncvar_shdfac_monthly
        character*128      :: LDT_ncvar_vegetype
        character*128      :: LDT_ncvar_soiltype
        character*128      :: LDT_ncvar_slopetype
        character*128      :: LDT_ncvar_smceq
        character*128      :: LDT_ncvar_tbot
        character*128      :: LDT_ncvar_pblh
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
        !!! MB: AC71
        integer            :: init_daynri

        !-------------------------------------------------------------------------
        ! Constant Parameter
        !-------------------------------------------------------------------------
        character(len=256) :: landuse_tbl_name
        character(len=256) :: soil_tbl_name
        character(len=256) :: gen_tbl_name
        character(len=256) :: ac_tbl_name
        character(len=256) :: landuse_scheme_name
        integer            :: nlucats
        character(len=256) :: soil_scheme_name
        real               :: dt
        integer            :: nsoil
        !!! MB: AC71
        integer            :: daynrinextclimaterecord
        character(len=256) :: PathNameOutp
        character(len=256) :: PathNameSimul
        character(len=256) :: PathNameList
        character(len=256) :: PathNameParam
        character(len=256) :: Climate_Filename
        character(len=256) :: ETo_Filename
        character(len=256) :: Rain_Filename
        character(len=256) :: CO2_Filename
        character(len=256) :: Crop_Filename
        character(len=256) :: Management_Filename
        character(len=256) :: Irrigation_Filename
        integer            :: Crop_AnnualStartDay
        integer            :: Crop_AnnualEndDay
        integer            :: Crop_AnnualStartMonth
        integer            :: Crop_AnnualEndMonth
        integer            :: NrSoilLayers
        integer            :: max_No_compartments
        integer            :: Tmin_windowsize
        real, pointer      :: Thickness(:)
        !!! MB: AC71
        type(Ac71dec), pointer :: ac71(:)
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
        integer  :: n, t     
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
            !!! MB: AC71
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(AC71_struc(n)%ac71(t)%smc(AC71_struc(n)%max_No_compartments))
            enddo
            
            ! initialize forcing variables to zeros
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%PREC_ac = 0.0
                AC71_struc(n)%ac71(t)%TMIN_ac = 0.0
                AC71_struc(n)%ac71(t)%TMAX_ac = 0.0
                AC71_struc(n)%ac71(t)%ETo_ac = 0.0

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
!------------------------------------------------------------------------
! TODO: setup number of soil moisture/temperature layers and depth here  
!------------------------------------------------------------------------
            ! TODO: set number of soil moisture layers in surface model
            LIS_sfmodel_struc(n)%nsm_layers = AC71_struc(n)%max_No_compartments
            ! TODO: set number of soil temperature layers in surface model
            LIS_sfmodel_struc(n)%nst_layers = AC71_struc(n)%max_No_compartments
            allocate(LIS_sfmodel_struc(n)%lyrthk(AC71_struc(n)%max_No_compartments))
            LIS_sfmodel_struc(n)%lyrthk(:) = 0.1
            LIS_sfmodel_struc(n)%ts = AC71_struc(n)%ts
        enddo
    end subroutine Ac71_ini
end module Ac71_lsmMod
