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
module AWRAL600_lsmMod
!BOP
!
! !MODULE: AWRAL600_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of AWRAL600 model. It also provides the entry method
! for the initialization of AWRAL600-specific variables. The derived
! data type {\tt AWRAL600\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the AWRAL600 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for AWRAL600
! \item[LDT\_ncvar\_k\_rout]
!   LDT NetCDF variable name for rate coefficient controlling discharge to stream
! \item[LDT\_ncvar\_kssat]
!   LDT NetCDF variable name for saturated hydraulic conductivity of shallow soil layer
! \item[LDT\_ncvar\_prefr]
!   LDT NetCDF variable name for reference value for precipitation
! \item[LDT\_ncvar\_s0max]
!   LDT NetCDF variable name for maximum storage of the surface soil layer
! \item[LDT\_ncvar\_slope]
!   LDT NetCDF variable name for slope of the land surface
! \item[LDT\_ncvar\_ssmax]
!   LDT NetCDF variable name for maximum storage of the shallow soil layer
! \item[LDT\_ncvar\_k\_gw]
!   LDT NetCDF variable name for groundwater drainage coefficient
! \item[LDT\_ncvar\_kr\_sd]
!   LDT NetCDF variable name for routing delay factor for the deep layer
! \item[LDT\_ncvar\_kr\_0s]
!   LDT NetCDF variable name for routing delay factor for the surface layer
! \item[LDT\_ncvar\_k0sat]
!   LDT NetCDF variable name for saturated hydraulic conductivity of surface soil layer
! \item[LDT\_ncvar\_sdmax]
!   LDT NetCDF variable name for maximum storage of the deep soil layer
! \item[LDT\_ncvar\_kdsat]
!   LDT NetCDF variable name for saturated hydraulic conductivity of shallow soil layer
! \item[LDT\_ncvar\_ne]
!   LDT NetCDF variable name for effective porosity
! \item[ts]
!   AWRAL600 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[awral600]
!  AWRAL600 model specific variables
! \item[forc\_count]
!   counter of forcing data
! \item[slope\_coeff]
!   scaling factor for slope
! \item[pair]
!   air pressure
! \item[kr\_coeff]
!   scaling factor for ratio of saturated hydraulic conductivity
! \item[nhru]
!   number of hydrologic response units
! \item[nhypsbins]
!   number of hypsometric curve distribution percentile bins
! \item[hypsperc]
!   hypsometric curve distribution percentile bins
! \item[alb\_dry]
!   dry soil albedo for each hru
! \item[alb\_wet]
!   wet soil albedo for each hru
! \item[cgsmax]
!   coefficient relating vegetation photosynthetic capacity to maximum stomatal conductance for each hru
! \item[er\_frac\_ref]
!   specific ratio of the mean evaporation rate and the mean rainfall intensity during storms for each hru
! \item[fsoilemax]
!   soil evaporation scaling factor corresponding to unlimited soil water supply for each hru
! \item[lairef]
!   reference leaf area index (at which fv = 0.63) for each hru
! \item[rd]
!   rooting depth for each hru
! \item[s\_sls]
!   specific canopy rainfall storage per unit leaf area for each hru
! \item[sla]
!   specific leaf area for each hru
! \item[tgrow]
!   characteristic time scale for vegetation growth towards equilibrium for each hru
! \item[tsenc]
!   characteristic time scale for vegetation senescence towards equilibrium for each hru
! \item[ud0]
!   maximum possible root water uptake from the deep soil store for each hru
! \item[us0]
!   maximum possible root water uptake from the shallow soil store for each hru
! \item[vc]
!   vegetation photosynthetic capacity index per unit canopy cover for each hru
! \item[w0lime]
!   limiting the value of the relative soil moisture content of the top soil layer at which evaporation is reduced for each hru
! \item[w0ref\_alb]
!   Reference value of w0 that determines the rate of albedo decrease with wetness for each hru
! \item[wdlimu]
!   water-limiting relative water content of the deep soil store for each hru
! \item[wslimu]
!   water-limiting relative water content of the shallow soil store for each hruWRAL600_struc
! \item[timesteps]
!   number of daily timesteps
! \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  12/18/18: Wendy Sharples, Shugong Wang Initial implementation for LIS 7 and AWRAL600
!
! !USES:
    use AWRAL600_module
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none

    PRIVATE
    !-------------------------------------------------------------------------
    ! PUBLIC MEMBER FUNCTIONS
    !-------------------------------------------------------------------------
    public :: AWRAL600_lsm_ini
    !-------------------------------------------------------------------------
    ! PUBLIC TYPES
    !-------------------------------------------------------------------------
    public :: AWRAL600_struc
!EOP
    type, public :: AWRAL600_type_dec
        character(len=LIS_CONST_PATH_LEN) :: rfile
        character*256      :: rformat
        !-------------------------------------------------------------------------
        ! Parameter file names
        !-------------------------------------------------------------------------
        character*128      :: LDT_ncvar_k_rout
        character*128      :: LDT_ncvar_kssat
        character*128      :: LDT_ncvar_prefr
        character*128      :: LDT_ncvar_s0max
        character*128      :: LDT_ncvar_slope
        character*128      :: LDT_ncvar_ssmax
        character*128      :: LDT_ncvar_k_gw
        character*128      :: LDT_ncvar_kr_sd
        character*128      :: LDT_ncvar_kr_0s
        character*128      :: LDT_ncvar_k0sat
        character*128      :: LDT_ncvar_sdmax
        character*128      :: LDT_ncvar_kdsat
        character*128      :: LDT_ncvar_ne
        character*128      :: LDT_ncvar_height
        character*128      :: LDT_ncvar_fhru
        character*128      :: LDT_ncvar_hveg
        character*128      :: LDT_ncvar_laimax
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
        real               :: init_sr
        real               :: init_sg
        real, allocatable  :: init_s0(:)
        real, allocatable  :: init_ss(:)
        real, allocatable  :: init_sd(:)
        real, allocatable  :: init_mleaf(:)
        !-------------------------------------------------------------------------
        ! Constant Parameter
        !-------------------------------------------------------------------------
        real               :: slope_coeff
        real               :: pair
        real               :: kr_coeff
        integer            :: nhru
        integer            :: nhypsbins
        real, allocatable  :: hypsperc(:)
        real, allocatable  :: alb_dry(:)
        real, allocatable  :: alb_wet(:)
        real, allocatable  :: cgsmax(:)
        real, allocatable  :: er_frac_ref(:)
        real, allocatable  :: fsoilemax(:)
        real, allocatable  :: lairef(:)
        real, allocatable  :: rd(:)
        real, allocatable  :: s_sls(:)
        real, allocatable  :: sla(:)
        real, allocatable  :: tgrow(:)
        real, allocatable  :: tsenc(:)
        real, allocatable  :: ud0(:)
        real, allocatable  :: us0(:)
        real, allocatable  :: vc(:)
        real, allocatable  :: w0lime(:)
        real, allocatable  :: w0ref_alb(:)
        real, allocatable  :: wdlimu(:)
        real, allocatable  :: wslimu(:)
        real, allocatable  :: hy
        integer            :: timesteps
        type(AWRAL600dec), allocatable :: awral600(:)
    end type AWRAL600_type_dec

    type(AWRAL600_type_dec), allocatable :: AWRAL600_struc(:)
 
contains 

!BOP
!
! !ROUTINE: AWRAL600_lsm_ini
! \label{AWRAL600_lsm_ini}
!
! !INTERFACE:
    subroutine AWRAL600_lsm_ini()
! !USES:
        use LIS_coreMod, only : LIS_rc
        use LIS_logMod, only : LIS_verify
        use LIS_timeMgrMod, only : LIS_clock,  LIS_calendar, &
            LIS_update_timestep, LIS_registerAlarm
        use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for AWRAL600-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for AWRAL600 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[AWRAL600\_readcrd](\ref{AWRAL600_readcrd})\\
!    reads the runtime options for AWRAL600 model
!  \end{description}
!EOP
        implicit none        
        integer  :: n, t     
        integer  :: status   

        ! allocate memory for nest 
        allocate(AWRAL600_struc(LIS_rc%nnest))
 
        ! read configuation information from lis.config file
        call AWRAL600_readcrd()

        do n=1, LIS_rc%nnest
            ! allocate memory for all tiles in current nest 
            allocate(AWRAL600_struc(n)%awral600(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            !------------------------------------------------------------------------
            ! allocate memory for vector variables passed to model interfaces        
            !------------------------------------------------------------------------
            ! allocate memory for multilevel spatial parameter
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(AWRAL600_struc(n)%awral600(t)%height(AWRAL600_struc(n)%nhypsbins))
                allocate(AWRAL600_struc(n)%awral600(t)%fhru(AWRAL600_struc(n)%nhru))
                allocate(AWRAL600_struc(n)%awral600(t)%hveg(AWRAL600_struc(n)%nhru))
                allocate(AWRAL600_struc(n)%awral600(t)%laimax(AWRAL600_struc(n)%nhru))
            enddo
            ! allocate memory for state variables
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(AWRAL600_struc(n)%awral600(t)%s0(AWRAL600_struc(n)%nhru))
                allocate(AWRAL600_struc(n)%awral600(t)%ss(AWRAL600_struc(n)%nhru))
                allocate(AWRAL600_struc(n)%awral600(t)%sd(AWRAL600_struc(n)%nhru))
                allocate(AWRAL600_struc(n)%awral600(t)%mleaf(AWRAL600_struc(n)%nhru))
            enddo
            !------------------------------------------------------------------------
            ! Model timestep Alarm
            !------------------------------------------------------------------------
            AWRAL600_struc(n)%forc_count = 0

            call LIS_update_timestep(LIS_rc, n, AWRAL600_struc(n)%ts)

            call LIS_registerAlarm("AWRAL600 model alarm",&
                                   AWRAL600_struc(n)%ts, &
                                   AWRAL600_struc(n)%ts)

            call LIS_registerAlarm("AWRAL600 restart alarm", &
                                   AWRAL600_struc(n)%ts,&
                                   AWRAL600_struc(n)%rstInterval)
            !------------------------------------------------------------------------
            ! setup number of soil moisture/temperature layers and depth here  
            !------------------------------------------------------------------------
            ! set number of soil moisture layers in surface model
            LIS_sfmodel_struc(n)%nsm_layers = 1
            ! set number of soil temperature layers in surface model
            LIS_sfmodel_struc(n)%nst_layers = 1
            allocate(LIS_sfmodel_struc(n)%lyrthk(1))
            LIS_sfmodel_struc(n)%ts = AWRAL600_struc(n)%ts
        enddo
    end subroutine AWRAL600_lsm_ini
end module AWRAL600_lsmMod
