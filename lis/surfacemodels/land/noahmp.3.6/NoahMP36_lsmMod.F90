!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
module NoahMP36_lsmMod
!BOP
!
! !MODULE: NoahMP36_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of NoahMP36 model. It also provides the entry method
! for the initialization of NoahMP36-specific variables. The derived
! data type {\tt NoahMP36\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the NoahMP36 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for NoahMP36
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
!   NoahMP36 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[noahmp36]
!  NoahMP36 model specific variables
! \item[forc\_count]
!   counter of forcing data
! \item[landuse\_tbl\_name]
!   Noah model landuse parameter table
! \item[soil\_tbl\_name]
!   Noah model soil parameter table
! \item[gen\_tbl\_name]
!   Noah model general parameter table
! \item[noahmp\_tbl\_name]
!   NoahMP parameter table
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
!  9/4/14: Shugong Wang Initial implementation for LIS 7 and NoahMP36
!
! !USES:
    use NoahMP36_module

    implicit none

    PRIVATE
    !-------------------------------------------------------------------------
    ! PUBLIC MEMBER FUNCTIONS
    !-------------------------------------------------------------------------
    public :: NoahMP36_ini
    !-------------------------------------------------------------------------
    ! PUBLIC TYPES
    !-------------------------------------------------------------------------
    public :: NoahMP36_struc
!EOP
    type, public :: NoahMP36_type_dec
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
        real               :: init_albold
        real               :: init_sneqvo
        real, pointer      :: init_stc(:)
        real, pointer      :: init_sh2o(:)
        real, pointer      :: init_smc(:)
        real               :: init_tah
        real               :: init_eah
        real               :: init_fwet
        real               :: init_canliq
        real               :: init_canice
        real               :: init_tv
        real               :: init_tg
        real               :: init_qsnow
        !integer            :: init_isnow
        !real, pointer      :: init_zss(:)
        real               :: init_snowh
        real               :: init_sneqv
        !real, pointer      :: init_snowice(:)
        !real, pointer      :: init_snowliq(:)
        real               :: init_zwt
        real               :: init_wa
        real               :: init_wt
        real               :: init_wslake
        real               :: init_lfmass
        real               :: init_rtmass
        real               :: init_stmass
        real               :: init_wood
        real               :: init_stblcp
        real               :: init_fastcp
        real               :: init_lai
        real               :: init_sai
        real               :: init_cm
        real               :: init_ch
        real               :: init_tauss
        real               :: init_smcwtd
        real               :: init_deeprech
        real               :: init_rech
        real               :: init_zlvl 
        !-------------------------------------------------------------------------
        ! Constant Parameter
        !-------------------------------------------------------------------------
        character(len=256) :: landuse_tbl_name
        character(len=256) :: soil_tbl_name
        character(len=256) :: gen_tbl_name
        character(len=256) :: noahmp_tbl_name
        character(len=256) :: landuse_scheme_name
        character(len=256) :: soil_scheme_name
        integer            :: dveg_opt
        integer            :: crs_opt
        integer            :: btr_opt
        integer            :: run_opt
        integer            :: sfc_opt
        integer            :: frz_opt
        integer            :: inf_opt
        integer            :: rad_opt
        integer            :: alb_opt
        integer            :: snf_opt
        integer            :: tbot_opt
        integer            :: stc_opt
        integer            :: nslcats
        integer            :: nlucats
        integer            :: nslpcats
        real               :: dt
        integer            :: nsoil
        real, pointer      :: sldpth(:)
        integer            :: nsnow
        integer            :: urban_vegetype
        integer            :: ice_flag
        integer            :: st_flag
        integer            :: sc_idx
        integer            :: iz0tlnd
        !real               :: zlvl
        type(NoahMP36dec), pointer :: noahmp36(:)
    end type NoahMP36_type_dec

    type(NoahMP36_type_dec), pointer :: NOAHMP36_struc(:)
 
contains 

!BOP
!
! !ROUTINE: NoahMP36_ini
! \label{NoahMP36_ini}
!
! !INTERFACE:
    subroutine NoahMP36_ini()
! !USES:
        use LIS_coreMod, only : LIS_rc
        use LIS_logMod, only : LIS_verify
        use LIS_timeMgrMod, only : LIS_clock,  LIS_calendar, &
            LIS_update_timestep, LIS_registerAlarm
        use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for NoahMP36-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for NoahMP36 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[NoahMP36\_readcrd](\ref{NoahMP36_readcrd}) \newline
!    reads the runtime options for NoahMP36 model
!  \end{description}
!EOP
        implicit none        
        integer  :: n, t     
        integer  :: status   

        ! allocate memory for nest 
        allocate(NOAHMP36_struc(LIS_rc%nnest))
 
        ! read configuation information from lis.config file
        call NoahMP36_readcrd()

        do n=1, LIS_rc%nnest
            ! allocate memory for all tiles in current nest 
            allocate(NOAHMP36_struc(n)%noahmp36(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            !------------------------------------------------------------------------
            ! allocate memory for vector variables passed to model interfaces        
            ! TODO: check the following allocation statements carefully!
            !------------------------------------------------------------------------
            ! allocate memory for multilevel spatial parameter
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(NOAHMP36_struc(n)%noahmp36(t)%shdfac_monthly(12))
                allocate(NOAHMP36_struc(n)%noahmp36(t)%smceq(NOAHMP36_struc(n)%nsoil))
            enddo
            ! allocate memory for state variables
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(NOAHMP36_struc(n)%noahmp36(t)%sstc( NOAHMP36_struc(n)%nsoil + NOAHMP36_struc(n)%nsnow))
                allocate(NOAHMP36_struc(n)%noahmp36(t)%sh2o(NOAHMP36_struc(n)%nsoil))
                allocate(NOAHMP36_struc(n)%noahmp36(t)%smc(NOAHMP36_struc(n)%nsoil))
                allocate(NOAHMP36_struc(n)%noahmp36(t)%zss( NOAHMP36_struc(n)%nsoil + NOAHMP36_struc(n)%nsnow))
                allocate(NOAHMP36_struc(n)%noahmp36(t)%snowice(NOAHMP36_struc(n)%nsnow))
                allocate(NOAHMP36_struc(n)%noahmp36(t)%snowliq(NOAHMP36_struc(n)%nsnow))
            enddo
!            ! allocate memory for intiali state variables
!            allocate(NOAHMP36_struc(n)%init_stc( NOAHMP36_struc(n)%nsoil))
!            allocate(NOAHMP36_struc(n)%init_sh2o(NOAHMP36_struc(n)%nsoil))
!            allocate(NOAHMP36_struc(n)%init_smc(NOAHMP36_struc(n)%nsoil))
            !allocate(NOAHMP36_struc(n)%init_zss( NOAHMP36_struc(n)%nsoil + NOAHMP36_struc(n)%nsnow))
            !allocate(NOAHMP36_struc(n)%init_snowice(NOAHMP36_struc(n)%nsnow))
            !allocate(NOAHMP36_struc(n)%init_snowliq(NOAHMP36_struc(n)%nsnow))
            
            ! initialize forcing variables to zeros
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                NOAHMP36_struc(n)%noahmp36(t)%lwdown = 0.0
                NOAHMP36_struc(n)%noahmp36(t)%swdown = 0.0
                NOAHMP36_struc(n)%noahmp36(t)%psurf = 0.0
                NOAHMP36_struc(n)%noahmp36(t)%prcp = 0.0
                NOAHMP36_struc(n)%noahmp36(t)%tair = 0.0
                NOAHMP36_struc(n)%noahmp36(t)%qair = 0.0
                NOAHMP36_struc(n)%noahmp36(t)%wind_e = 0.0
                NOAHMP36_struc(n)%noahmp36(t)%wind_n = 0.0
                
                NOAHMP36_struc(n)%noahmp36(t)%albd = -9999.0
                NOAHMP36_struc(n)%noahmp36(t)%albi = -9999.0
                NOAHMP36_struc(n)%noahmp36(t)%alb_upd_flag = .false.
            enddo ! end of tile (t) loop
!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
            NOAHMP36_struc(n)%forc_count = 0

            call LIS_update_timestep(LIS_rc, n, NOAHMP36_struc(n)%ts)

            call LIS_registerAlarm("NoahMP36 model alarm",&
                                   NOAHMP36_struc(n)%ts, &
                                   NOAHMP36_struc(n)%ts)

            call LIS_registerAlarm("NoahMP36 restart alarm", &
                                   NOAHMP36_struc(n)%ts,&
                                   NOAHMP36_struc(n)%rstInterval)
!------------------------------------------------------------------------
! TODO: setup number of soil moisture/temperature layers and depth here  
!------------------------------------------------------------------------
            ! TODO: set number of soil moisture layers in surface model
            LIS_sfmodel_struc(n)%nsm_layers = NOAHMP36_struc(n)%nsoil
            ! TODO: set number of soil temperature layers in surface model
            LIS_sfmodel_struc(n)%nst_layers = NOAHMP36_struc(n)%nsoil
            allocate(LIS_sfmodel_struc(n)%lyrthk(NOAHMP36_struc(n)%nsoil))
            LIS_sfmodel_struc(n)%lyrthk(:) = NOAHMP36_struc(n)%sldpth(:) 
            LIS_sfmodel_struc(n)%ts = NOAHMP36_struc(n)%ts
        enddo
    end subroutine NoahMP36_ini
end module NoahMP36_lsmMod
