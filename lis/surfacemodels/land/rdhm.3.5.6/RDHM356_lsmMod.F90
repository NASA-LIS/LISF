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
module RDHM356_lsmMod
!BOP
!
! !MODULE: RDHM356_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of RDHM356 model. It also provides the entry method
! for the initialization of RDHM356-specific variables. The derived
! data type {\tt RDHM356\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the RDHM356 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for RDHM356
! \item[LDT\_ncvar\_SoilAlb]
!   LDT NetCDF variable name for snow free ALBEDO (default value 0.15)
! \item[LDT\_ncvar\_SnowAlb]
!   LDT NetCDF variable name for snow ALBEDO (default value 0.7)
! \item[LDT\_ncvar\_SOILTYP]
!   LDT NetCDF variable name for Soil type
! \item[LDT\_ncvar\_VEGETYP]
!   LDT NetCDF variable name for Vegetation type
! \item[LDT\_ncvar\_UZTWM]
!   LDT NetCDF variable name for upper zone tension water maximum storage
! \item[LDT\_ncvar\_UZFWM]
!   LDT NetCDF variable name for upper zone free water maximum storage
! \item[LDT\_ncvar\_UZK]
!   LDT NetCDF variable name for upper zone free water latent depletion rate
! \item[LDT\_ncvar\_PCTIM]
!   LDT NetCDF variable name for impervious fraction of the watershad area
! \item[LDT\_ncvar\_ADIMP]
!   LDT NetCDF variable name for additional impervious area
! \item[LDT\_ncvar\_RIVA]
!   LDT NetCDF variable name for riparian vegetation area
! \item[LDT\_ncvar\_ZPERC]
!   LDT NetCDF variable name for maximum percolation rate
! \item[LDT\_ncvar\_REXP]
!   LDT NetCDF variable name for exponent of the percolation equation (percolation parameter)
! \item[LDT\_ncvar\_LZTWM]
!   LDT NetCDF variable name for lower zone tension water maximum storage
! \item[LDT\_ncvar\_LZFSM]
!   LDT NetCDF variable name for lower zone supplemental free water (fast) maximum storage
! \item[LDT\_ncvar\_LZFPM]
!   LDT NetCDF variable name for lower zone primary free water (slow) maximum storage
! \item[LDT\_ncvar\_LZSK]
!   LDT NetCDF variable name for lower zone supplemental free water depletion rate
! \item[LDT\_ncvar\_LZPK]
!   LDT NetCDF variable name for lower zone primary free water depletion rate
! \item[LDT\_ncvar\_PFREE]
!   LDT NetCDF variable name for fraction percolation from upper to lower free water storage
! \item[LDT\_ncvar\_SIDE]
!   LDT NetCDF variable name for ratio of deep recharge to channel base flow
! \item[LDT\_ncvar\_RSERV]
!   LDT NetCDF variable name for fraction of lower zone free water not transferable to tension water
! \item[LDT\_ncvar\_EFC]
!   LDT NetCDF variable name for fraction of forest cover
! \item[LDT\_ncvar\_TBOT]
!   LDT NetCDF variable name for bottom boundary soil temperature
! \item[LDT\_ncvar\_RSMAX]
!   LDT NetCDF variable name for maximum residual porosity (usually = 0.58)
! \item[LDT\_ncvar\_CKSL]
!   LDT NetCDF variable name for ratio of frozen to non-frozen surface (increase in frozen ground contact, usually = 8 s/m)
! \item[LDT\_ncvar\_ZBOT]
!   LDT NetCDF variable name for lower boundary depth (negative value, usually = -2.5 m)
! \item[LDT\_ncvar\_ALON]
!   LDT NetCDF variable name for logitude
! \item[LDT\_ncvar\_ALAT]
!   LDT NetCDF variable name for latitude
! \item[LDT\_ncvar\_SCF]
!   LDT NetCDF variable name for snow fall correction factor
! \item[LDT\_ncvar\_MFMAX]
!   LDT NetCDF variable name for maximum melt factor
! \item[LDT\_ncvar\_MFMIN]
!   LDT NetCDF variable name for minimum melt factor
! \item[LDT\_ncvar\_NMF]
!   LDT NetCDF variable name for maximum negative melt factor
! \item[LDT\_ncvar\_UADJ]
!   LDT NetCDF variable name for the average wind function during rain-on-snow periods
! \item[LDT\_ncvar\_SI]
!   LDT NetCDF variable name for areal water-equivalent above which 100 percent areal snow cover
! \item[LDT\_ncvar\_MBASE]
!   LDT NetCDF variable name for base temperature for non-rain melt factor
! \item[LDT\_ncvar\_PXTEMP]
!   LDT NetCDF variable name for temperature which spereates rain from snow
! \item[LDT\_ncvar\_PLWHC]
!   LDT NetCDF variable name for maximum amount of liquid-water held against gravity drainage
! \item[LDT\_ncvar\_TIPM]
!   LDT NetCDF variable name for antecedent snow temperature index parameter
! \item[LDT\_ncvar\_GM]
!   LDT NetCDF variable name for daily ground melt
! \item[LDT\_ncvar\_ELEV]
!   LDT NetCDF variable name for elevation
! \item[LDT\_ncvar\_LAEC]
!   LDT NetCDF variable name for snow-rain split temperature
! \item[ts]
!   RDHM356 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[rdhm356]
!  RDHM356 model specific variables
! \item[forc\_count]
!   counter of forcing data
! \item[TempHeight]
!   observation height of temperature of humidity
! \item[WindHeight]
!   observation height of wind
! \item[DT\_SAC\_SNOW17]
!   simulation time interval of SAC model and Snow-17
! \item[DT\_FRZ]
!   simulation time interval of frozen soil model
! \item[FRZ\_VER\_OPT]
!   version number of frozen soil model. 1: old version, 2: new version
! \item[SNOW17\_OPT]
!   option for snow-17. If SNOW17\_OPT=1, use SNOW-17, otherwise, don't use
! \item[NSTYP]
!   number of soil types
! \item[NVTYP]
!   number of vegetation types
! \item[NDINTW]
!   number of desired soil layers for total and liquid soil moisture
! \item[NDSINT]
!   number of desired soil layers for soil temperature
! \item[NORMALIZE]
!   normalization flag for total and liquid soil moisture output (1-normalized, 0-not)
! \item[DSINTW]
!   thickness of desired soil layers for liquid and total soil moisture
! \item[DSINT]
!   thickness of desired soil layers for soil temperature
! \item[PETADJ\_MON]
!   adjustment of PET for 12 months
! \item[vegRCMIN]
!   minimal stomatal resistance table for SACHTET, 14 values
! \item[climRCMIN]
!   climate dependent miminal stomatal resistance for SACHTET, 14 values
! \item[RGL]
!   solar radiation threshold table for SACHTET, 14 values
! \item[HS]
!   vapor pressure resistance factor table for SACHTET, 14 values
! \item[LAI]
!   leaf area index table for SACHTET, 14 values
! \item[D50]
!   the depth (cm) table at which 50\% roots are allocated for SACHTET, 14 values
! \item[CROOT]
!   root distribution parameter table for SACHTET, 14 values
! \item[Z0]
!   roughness coefficient of surface
! \item[CLAY]
!   clay content for SACHTET, 12 values
! \item[SAND]
!   sand content for sACHTET, 12 values
! \item[SATDK]
!   saturated hydraulic conductivityfor SACHTET, 12 values
! \item[CZIL]
!   default=0.12 Zilitinkevich
! \item[FXEXP]
!   FXEXP(fxexp),(default=2.0) bare soil
! \item[vegRCMAX]
!   RCMAX,(default=5000s/m) maximum stomatal resistance
! \item[TOPT]
!   TOPT,(default=298K)optimum air
! \item[PC]
!   plant coef. default pc = -1, 0.6 - 0.8
! \item[PET\_OPT]
!   if PET\_OPT = 0, use non Penmann-based ETP;if penpt > 0 empirical Penmann equation; if penpt < 0, use energy based Pennman
! \item[RDST]
!   default=1 means noah option,this constant allows selection of tension water redistribution option, if rdst = 0 (ohd), use OHD version of SRT subroutine this SRT uses reference gradient instead an actual. if rdst = 1 ( noah), use Noah version of SRT subroutine
! \item[thresholdRCMIN]
!   this constant allows change of RCMIN (0.5)
! \item[SFCREF]
!   reference wind speed for PET adjustment (4 m s-1)
! \item[BAREADJ]
!   Ek-Chen evaporation threshold switch. Bare soil evaporation option changes according to greenness.
! \item[SNOW17\_SWITCH]
!   switch variable change liquid water freezing version, 0: Victor's version, 1: Eric's version
! \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  11/5/13: Shugong Wang Initial implementation for LIS 7 and RDHM356
!
! !USES:
    use RDHM356_module
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none

    PRIVATE
    !-------------------------------------------------------------------------
    ! PUBLIC MEMBER FUNCTIONS
    !-------------------------------------------------------------------------
    public :: RDHM356_ini
    !-------------------------------------------------------------------------
    ! PUBLIC TYPES
    !-------------------------------------------------------------------------
    public :: RDHM356_struc
!EOP
    type, public :: RDHM356_type_dec
        character(len=LIS_CONST_PATH_LEN) :: rfile
        character*256      :: rformat
        character(len=LIS_CONST_PATH_LEN) :: tmxmn_dir
        !-------------------------------------------------------------------------
        ! Parameter file names
        !-------------------------------------------------------------------------
        character*128      :: LDT_ncvar_Tair_min
        character*128      :: LDT_ncvar_Tair_max
        character*128      :: LDT_ncvar_PET_MON
        character*128      :: LDT_ncvar_GRN_MON
        character*128      :: LDT_ncvar_SoilAlb
        character*128      :: LDT_ncvar_SnowAlb
        character*128      :: LDT_ncvar_SOILTYP
        character*128      :: LDT_ncvar_VEGETYP
        character*128      :: LDT_ncvar_UZTWM
        character*128      :: LDT_ncvar_UZFWM
        character*128      :: LDT_ncvar_UZK
        character*128      :: LDT_ncvar_PCTIM
        character*128      :: LDT_ncvar_ADIMP
        character*128      :: LDT_ncvar_RIVA
        character*128      :: LDT_ncvar_ZPERC
        character*128      :: LDT_ncvar_REXP
        character*128      :: LDT_ncvar_LZTWM
        character*128      :: LDT_ncvar_LZFSM
        character*128      :: LDT_ncvar_LZFPM
        character*128      :: LDT_ncvar_LZSK
        character*128      :: LDT_ncvar_LZPK
        character*128      :: LDT_ncvar_PFREE
        character*128      :: LDT_ncvar_SIDE
        character*128      :: LDT_ncvar_RSERV
        character*128      :: LDT_ncvar_EFC
        character*128      :: LDT_ncvar_TBOT
        character*128      :: LDT_ncvar_RSMAX
        character*128      :: LDT_ncvar_CKSL
        character*128      :: LDT_ncvar_ZBOT
        character*128      :: LDT_ncvar_vegRCMIN
        character*128      :: LDT_ncvar_climRCMIN
        character*128      :: LDT_ncvar_RGL
        character*128      :: LDT_ncvar_HS
        character*128      :: LDT_ncvar_LAI
        character*128      :: LDT_ncvar_D50
        character*128      :: LDT_ncvar_CROOT
        character*128      :: LDT_ncvar_Z0
        character*128      :: LDT_ncvar_CLAY
        character*128      :: LDT_ncvar_SAND
        character*128      :: LDT_ncvar_SATDK
        character*128      :: LDT_ncvar_ALON
        character*128      :: LDT_ncvar_ALAT
        character*128      :: LDT_ncvar_SCF
        character*128      :: LDT_ncvar_MFMAX
        character*128      :: LDT_ncvar_MFMIN
        character*128      :: LDT_ncvar_NMF
        character*128      :: LDT_ncvar_UADJ
        character*128      :: LDT_ncvar_SI
        character*128      :: LDT_ncvar_MBASE
        character*128      :: LDT_ncvar_PXTEMP
        character*128      :: LDT_ncvar_PLWHC
        character*128      :: LDT_ncvar_TIPM
        character*128      :: LDT_ncvar_GM
        character*128      :: LDT_ncvar_ELEV
        character*128      :: LDT_ncvar_LAEC
        character*128      :: LDT_ncvar_ADC
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
        real               :: init_UZTWC_ratio
        real               :: init_UZFWC_ratio
        real               :: init_LZTWC_ratio
        real               :: init_LZFPC_ratio
        real               :: init_LZFSC_ratio
        real               :: init_ADIMC_ratio
        real               :: init_TS0
        real               :: init_TS1
        real               :: init_TS2
        real               :: init_TS3
        real               :: init_TS4
        real               :: init_UZTWH_ratio
        real               :: init_UZFWH_ratio
        real               :: init_LZTWH_ratio
        real               :: init_LZFSH_ratio
        real               :: init_LZFPH_ratio
        real               :: init_SMC(6)
        real               :: init_SH2O(6)
        real               :: init_WE
        real               :: init_LIQW
        real               :: init_NEGHS
        real               :: init_TINDEX
        real               :: init_ACCMAX
        real               :: init_SNDPT
        real               :: init_SNTMP
        real               :: init_SB
        real               :: init_SBAESC
        real               :: init_SBWS
        real               :: init_STORAGE
        real               :: init_AEADJ
        real               :: init_EXLAG(7)
        integer            :: init_NEXLAG
        real               :: init_TA_PREV
        !-------------------------------------------------------------------------
        ! Constant Parameter
        !-------------------------------------------------------------------------
        real               :: TempHeight
        real               :: WindHeight
        real               :: DT_SAC_SNOW17
        real               :: DT_FRZ
        integer            :: FRZ_VER_OPT
        integer            :: SNOW17_OPT
        integer            :: SACHTET_OPT
        integer            :: NSTYP
        integer            :: NVTYP
        integer            :: NDINTW
        integer            :: NDSINT
        integer            :: NORMALIZE
        real, allocatable      :: DSINTW(:)
        real, allocatable      :: DSINT(:)
        real               :: PETADJ_MON(12)
        real               :: CZIL
        real               :: FXEXP
        real               :: vegRCMAX
        real               :: TOPT
        real               :: PC
        integer            :: PET_OPT
        integer            :: RDST
        real               :: thresholdRCMIN
        real               :: SFCREF
        real               :: BAREADJ
        integer            :: SNOW17_SWITCH
        !-------------------------------------------------------------------------
        ! Lookup Parameter
        !-------------------------------------------------------------------------
!        real, allocatable      :: vegRCMIN(:)
!        real, allocatable      :: climRCMIN(:)
!        real, allocatable      :: RGL(:)
!        real, allocatable      :: HS(:)
!        real, allocatable      :: LAI(:)
!        real, allocatable      :: D50(:)
!        real, allocatable      :: CROOT(:)
!        real, allocatable      :: Z0(:)
!        real, allocatable      :: CLAY(:)
!        real, allocatable      :: SAND(:)
!        real, allocatable      :: SATDK(:)
         integer            :: day
        type(RDHM356dec), allocatable :: rdhm356(:)
    end type RDHM356_type_dec

    type(RDHM356_type_dec), allocatable :: RDHM356_struc(:)
 
contains 

!BOP
!
! !ROUTINE: RDHM356_ini
! \label{RDHM356_ini}
!
! !INTERFACE:
    subroutine RDHM356_ini()
! !USES:
        use LIS_coreMod, only             : LIS_rc
        use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
        use LIS_timeMgrMod, only          : LIS_clock, LIS_calendar, &
                                            LIS_update_timestep,     &
                                            LIS_registerAlarm
        use LIS_logMod, only              : LIS_verify
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for RDHM356-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for RDHM356 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[RDHM356\_readcrd](\ref{RDHM356_readcrd}) \newline
!    reads the runtime options for RDHM356 model
!  \end{description}
!EOP
        implicit none
        integer                 :: n, t 
        integer                 :: status

        ! allocate memory for nest 
        allocate(RDHM356_struc(LIS_rc%nnest))
 
        ! read configuation information from lis.config file
        call RDHM356_readcrd()

        do n=1, LIS_rc%nnest
            ! allocate memory for all tiles in current nest 
            allocate(RDHM356_struc(n)%rdhm356(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            !------------------------------------------------------------------------
            ! allocate memory for vector variables passed to model interfaces        
            ! TODO: check the following allocation statements carefully!
            !------------------------------------------------------------------------
            ! allocate memory for lookup parameter
!            allocate(RDHM356_struc(n)%vegRCMIN(RDHM356_struc(n)%NVTYP))
!            allocate(RDHM356_struc(n)%climRCMIN(RDHM356_struc(n)%NVTYP))
!            allocate(RDHM356_struc(n)%RGL(RDHM356_struc(n)%NVTYP))
!            allocate(RDHM356_struc(n)%HS(RDHM356_struc(n)%NVTYP))
!            allocate(RDHM356_struc(n)%LAI(RDHM356_struc(n)%NVTYP))
!            allocate(RDHM356_struc(n)%D50(RDHM356_struc(n)%NVTYP))
!            allocate(RDHM356_struc(n)%CROOT(RDHM356_struc(n)%NVTYP))
!            allocate(RDHM356_struc(n)%Z0(RDHM356_struc(n)%NVTYP))
!            allocate(RDHM356_struc(n)%CLAY(RDHM356_struc(n)%NSTYP))
!            allocate(RDHM356_struc(n)%SAND(RDHM356_struc(n)%NSTYP))
!            allocate(RDHM356_struc(n)%SATDK(RDHM356_struc(n)%NSTYP))
            ! allocate memory for constant parameter
!            allocate(RDHM356_struc(n)%DSINTW(RDHM356_struc(n)%NDINTW))
!            allocate(RDHM356_struc(n)%DSINT(RDHM356_struc(n)%NDSINT))
            ! allocate memory for output variables
            ! allocate memory for output variables
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(RDHM356_struc(n)%rdhm356(t)%SWINT(RDHM356_struc(n)%NDINTW))
                allocate(RDHM356_struc(n)%rdhm356(t)%SWHINT(RDHM356_struc(n)%NDINTW))
                allocate(RDHM356_struc(n)%rdhm356(t)%TSINT(RDHM356_struc(n)%NDSINT))
            enddo
            !------------------------------------------------------------------------
            ! Model timestep Alarm
            !------------------------------------------------------------------------
            RDHM356_struc(n)%forc_count = 0

            do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
               RDHM356_struc(n)%rdhm356(t)%tair = 0 
               RDHM356_struc(n)%rdhm356(t)%qair = 0 
               RDHM356_struc(n)%rdhm356(t)%swdown = 0 
               RDHM356_struc(n)%rdhm356(t)%lwdown = 0
               RDHM356_struc(n)%rdhm356(t)%Wind_E = 0
               RDHM356_struc(n)%rdhm356(t)%Wind_N = 0 
               RDHM356_struc(n)%rdhm356(t)%psurf = 0 
               RDHM356_struc(n)%rdhm356(t)%rainf = 0 
               RDHM356_struc(n)%rdhm356(t)%snowf = 0 
            enddo

            call LIS_update_timestep(LIS_rc, n, RDHM356_struc(n)%ts)

            call LIS_registerAlarm("RDHM356 model alarm",&
                                   RDHM356_struc(n)%ts, &
                                   RDHM356_struc(n)%ts)

            call LIS_registerAlarm("RDHM356 restart alarm", &
                                   RDHM356_struc(n)%ts,&
                                   RDHM356_struc(n)%rstInterval)
            !------------------------------------------------------------------------
            ! TODO: setup number of soil moisture/temperature layers and depth here  
            !------------------------------------------------------------------------
            ! TODO: set number of soil moisture layers in surface model
            LIS_sfmodel_struc(n)%nsm_layers = 0
            ! TODO: set number of soil temperature layers in surface model
            LIS_sfmodel_struc(n)%nst_layers = 0
            allocate(LIS_sfmodel_struc(n)%lyrthk(0))
            LIS_sfmodel_struc(n)%ts = RDHM356_struc(n)%ts
        enddo
    end subroutine RDHM356_ini
end module RDHM356_lsmMod
