!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#define ESMF_STDERRORCHECK(rc) ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)
#define FILENAME "LIS_NUOPC_Gluecode"
#define MODNAME "LIS_NUOPC_Gluecode"

!-------------------------------------------------------------------------------
! Define ESMF real kind to match Appplications single/double precision
!-------------------------------------------------------------------------------
#if defined(REAL4)
#define ESMF_KIND_RX ESMF_KIND_R4
#define ESMF_TYPEKIND_RX ESMF_TYPEKIND_R4
#else
#define ESMF_KIND_RX ESMF_KIND_R8
#define ESMF_TYPEKIND_RX ESMF_TYPEKIND_R8
#endif

#define VERBOSITY_MIN 0
#define VERBOSITY_MAX 255
#define VERBOSITY_DBG 1023
#define UNMASKED 0
#define UNINITIALIZED -999

module LIS_NUOPC_Gluecode
!BOP
!
! !MODULE: LIS_NUOPC_Gluecode
!
! !DESCRIPTION:
!   This module connects NUOPC initialize, advance,
!   and finalize to LIS.
!
! !REVISION HISTORY:
!  13Oct15    Dan Rosen  Initial Specification
!
! !USES:
  use ESMF
  use NUOPC
  use LIS_coreMod, only: &
    LIS_config_init, &
    LIS_core_init, &
    LIS_finalize, &
    LIS_rc, &
    LIS_domain, &
    LIS_initialized, &
    LIS_vm, &
    LIS_masterproc, &
    LIS_localPet, &
    LIS_npes, &
    LIS_ews_ind, &
    LIS_ewe_ind, &
    LIS_nss_ind, &
    LIS_nse_ind, &
    LIS_ews_halo_ind, &
    LIS_ewe_halo_ind, &
    LIS_nss_halo_ind, &
    LIS_nse_halo_ind
  use LIS_domainMod, only: &
    LIS_domain_init         !initialize specified domains
  use LIS_surfaceModelMod, only: &
    LIS_surfaceModel_init, &
    LIS_surfaceModel_setup, &
    LIS_surfaceModel_readrestart, &
    LIS_surfaceModel_run, &
    LIS_surfaceModel_f2t, &
    LIS_surfaceModel_output, &
    LIS_surfaceModel_writerestart, &
    LIS_surfaceModel_perturb_states, &
    LIS_surfaceModel_setexport
  use LIS_metforcingMod, only: &
    LIS_metforcing_init, &  ! initialize met forcing setup
    LIS_get_met_forcing, &  ! retrieve data and interpolate spatially and temporally
    LIS_perturb_forcing, &  ! perturbs the met forcing variables
    LIS_FORC_State
  use LIS_DAobservationsMod, only: &
    LIS_initDAobservations, &
    LIS_readDAobservations, &
    LIS_perturb_DAobservations
  use LIS_perturbMod, only: &
    LIS_perturb_init, &     ! initialization for perturbation routines
    LIS_perturb_readrestart ! read perturbation restart files
  use LIS_dataAssimMod, only: &
    LIS_dataassim_init, &   !initialize data assimilation routines
    LIS_dataassim_run, &    !invoke data assimilation algorithms
    LIS_dataassim_output    !invoke data assimilation output
  use LIS_paramsMod, only: &
    LIS_param_init, &  !initializes structures, read static data
    LIS_setDynParams !read time dependent data
  use LISWRFGridCompMod, only: &
    LISWRF_alloc_states, &
    LISWRF_export
  use LIS_tbotAdjustMod, only: &
    LIS_createTmnUpdate
  use LIS_timeMgrMod, only: &
    LIS_timemgr_set      ! set the current time
  use LIS_FORC_AttributesMod, only: &
    forc_attrib_type, &     ! forcing arrribute type
    LIS_FORC_Tair, &
    LIS_FORC_Qair, &
    LIS_FORC_SWdown, &
    LIS_FORC_SWdirect, &
    LIS_FORC_SWdiffuse, &
    LIS_FORC_LWdown, &
    LIS_FORC_Wind_E, &
    LIS_FORC_Wind_N, &
    LIS_FORC_Psurf, &
    LIS_FORC_Rainf, &
    LIS_FORC_Snowf, &
    LIS_FORC_CRainf, &
    LIS_FORC_Forc_Hgt, &
    LIS_FORC_Ch, &
    LIS_FORC_Cm, &
    LIS_FORC_Emiss, &
    LIS_FORC_Q2sat, &
    LIS_FORC_Cosz, &
    LIS_FORC_Alb, &
    LIS_FORC_XICE, &
    LIS_FORC_QSFC, &
    LIS_FORC_CHS2, &
    LIS_FORC_CQS2, &
    LIS_FORC_T2, &
    LIS_FORC_Q2, &
    LIS_FORC_TH2, &
    LIS_FORC_TMN, &
    LIS_FORC_LPressure, &
    LIS_FORC_O3, &    !absorber
    LIS_FORC_PET, & ! SY for FEWSNET
    LIS_FORC_RefET, & ! SY for FEWSNET
    LIS_FORC_CMFORC, & ! dmm
    LIS_FORC_CHFORC, & ! dmm for NLDAS-2
    LIS_FORC_CAPE, & ! dmm for NLDAS-2
    LIS_FORC_PARDR, &
    LIS_FORC_PARDF, &
    LIS_FORC_SWNET, &
    LIS_FORC_SNOWFLAG, &
    LIS_FORC_DENSITY, &
    LIS_FORC_VAPORPRESS, &
    LIS_FORC_VAPORPRESSDEFICIT, &
    LIS_FORC_WIND, &
    LIS_FORC_Z0, &
    LIS_FORC_GVF, &
    LIS_FORC_CO2
  use NUOPC_LogUtility
  use NUOPC_FileReadUtility
  use NUOPC_FileWriteUtility

  IMPLICIT NONE

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_NUOPC_Init      ! init method for nuopc cpl mode
  public :: LIS_NUOPC_Run       ! run method for nuopc cpl mode
  public :: LIS_NUOPC_Final     ! finalize method for nuopc cpl mode
  public :: LIS_GridCreate
  public :: LIS_TimestepGet
  public :: LIS_NestCntGet
  public :: LIS_RunModeGet
  public :: LIS_Unknown
  public :: LIS_Offline
  public :: LIS_Coupled
  public :: LIS_Hybrid
  public :: LIS_FieldDictionaryAdd
  public :: LIS_Log
  public :: LIS_FieldList
  public :: LIS_FieldListLog

!-----------------------------------------------------------------------------
! Interface definitions for copy from LIS 1D tiled data to 2D
!-----------------------------------------------------------------------------

  interface LIS_CopyToLIS
    module procedure LIS_FieldCopyToLisField
    module procedure LIS_FieldCopyToLisFarray
    module procedure LIS_ArrayCopyToLisArray
    module procedure LIS_ArrayCopyToLisFarray
    module procedure LIS_FarrayI4CopyToLisFarrayI4
    module procedure LIS_FarrayI8CopyToLisFarrayI8
    module procedure LIS_FarrayR8CopyToLisFarrayR4
    module procedure LIS_FarrayR4CopyToLisFarrayR4
    module procedure LIS_FarrayR8CopyToLisFarrayR8
  end interface

  interface LIS_CopyFromLIS
    module procedure LIS_FieldCopyFromLisField
    module procedure LIS_FieldCopyFromLisFarray
    module procedure LIS_ArrayCopyFromLisArray
    module procedure LIS_ArrayCopyFromLisFarray
    module procedure LIS_FarrayI4CopyFromLisFarrayI4
    module procedure LIS_FarrayI8CopyFromLisFarrayI8
    module procedure LIS_FarrayR8CopyFromLisFarrayR4
    module procedure LIS_FarrayR4CopyFromLisFarrayR4
    module procedure LIS_FarrayR8CopyFromLisFarrayR8
  end interface

!-----------------------------------------------------------------------------
! !LOCAL VARIABLES:
!-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: LIS_Unknown = -1
  INTEGER, PARAMETER :: LIS_Offline =  0
  INTEGER, PARAMETER :: LIS_Coupled =  1
  INTEGER, PARAMETER :: LIS_Hybrid  =  2

  type LIS_FieldHookup
    real,pointer,dimension(:,:)     :: exportArray => null()
    real,pointer,dimension(:)       :: exportArray_t => null()
  end type

  type LIS_Field
    character(len=64)                  :: stdname = " "
    character(len=10)                  :: units = " "
    character(len=20)                  :: transferOffer = " "
    logical                            :: lisForc = .FALSE.
    logical                            :: lisForcSelect = .FALSE.
    logical                            :: lisExport = .FALSE.
    logical                            :: lisExportSelect = .FALSE.
    character(len=100)                 :: lisForcVarname = " "
    type(LIS_FieldHookup), allocatable :: hookup(:) ! Individual hookup for each nest
  end type

  type(LIS_Field),dimension(54)  :: LIS_FieldList = (/ &
    LIS_Field(stdname='aerodynamic_roughness_length', &
      units='m',transferOffer='will provide'), &
    LIS_Field(stdname='canopy_moisture_storage', &
      units='kg m-2',transferOffer='will provide'), & 
    LIS_Field(stdname='carbon_dioxide', &
      units='mol?',transferOffer='will provide'), &
    LIS_Field(stdname='cosine_zenith_angle', &
      units='?',transferOffer='will provide'), &
    LIS_Field(stdname='exchange_coefficient_heat', &
      units='?',transferOffer='will provide'), &
    LIS_Field(stdname='exchange_coefficient_heat_height2m', &
      units='?',transferOffer='will provide'), &
    LIS_Field(stdname='exchange_coefficient_moisture_height2m', &
      units='?',transferOffer='will provide'), &
    LIS_Field(stdname='ice_mask', &
      units='1',transferOffer='will provide'), &
    LIS_Field(stdname='inst_down_lw_flx', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='inst_down_sw_flx', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='inst_height_lowest', &
      units='m',transferOffer='will provide'), &
    LIS_Field(stdname='inst_merid_wind_height_lowest', &
      units='m s-1',transferOffer='will provide'), &
    LIS_Field(stdname='inst_pres_height_lowest', &
      units='Pa',transferOffer='will provide'), &
    LIS_Field(stdname='inst_pres_height_surface', &
      units='Pa',transferOffer='will provide'), &
    LIS_Field(stdname='inst_spec_humid_height_lowest', &
      units='kg kg-1',transferOffer='will provide'), &
    LIS_Field(stdname='inst_temp_height_lowest', &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdname='mean_temp_height_surface', &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdname='inst_wind_speed_height_lowest', &
      units='m s-1',transferOffer='will provide'), &
    LIS_Field(stdname='inst_zonal_wind_height_lowest', &
      units='m s-1',transferOffer='will provide'), &
    LIS_Field(stdname='liquid_water_content_of_soil_layer_1', &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdname='liquid_water_content_of_soil_layer_2', &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdname='liquid_water_content_of_soil_layer_3', &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdname='liquid_water_content_of_soil_layer_4', &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdname='liquid_water_content_of_surface_snow', &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_cprec_rate', &
      units='kg s-1 m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_down_lw_flx', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_down_sw_flx', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_fprec_rate', &
      units='kg s-1 m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_grnd_sensi_heat_flx', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_laten_heat_flx', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_laten_heat_flx_atm_into_lnd', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_laten_heat_flx_kinematic', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_net_lw_flx', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_prec_rate', &
      units='kg s-1 m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_sensi_heat_flx', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_sensi_heat_flx_atm_into_lnd', &
      units='W m-2',transferOffer='will provide'), &
    LIS_Field(stdname='mean_surface_albedo', &
      units='lm lm-1',transferOffer='will provide'), &
    LIS_Field(stdname='moisture_content_of_soil_layer_1', &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdname='moisture_content_of_soil_layer_2', &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdname='moisture_content_of_soil_layer_3', &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdname='moisture_content_of_soil_layer_4', &
      units='kg m-2',transferOffer='will provide'), &
    LIS_Field(stdname='saturated_mixing_ratio', &
      units='kg kg-1',transferOffer='will provide'), &
    LIS_Field(stdname='effective_mixing_ratio', &
      units='kg kg-1',transferOffer='will provide'), &
    LIS_Field(stdname='soil_temperature_bottom', &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdname='subsurface_runoff_flux', &
      units='kg s-1 m-2',transferOffer='will provide'), &
    LIS_Field(stdname='surface_microwave_emissivity', &
      units='?',transferOffer='will provide'), &
    LIS_Field(stdname='surface_runoff_flux', &
      units='kg s-1 m-2',transferOffer='will provide'), &
    LIS_Field(stdname='surface_snow_thickness', &
      units='m',transferOffer='will provide'), &
    LIS_Field(stdname='fractional_snow_cover', &
      units='1',transferOffer='will provide'), &
    LIS_Field(stdname='temperature_of_soil_layer_1', &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdname='temperature_of_soil_layer_2', &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdname='temperature_of_soil_layer_3', &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdname='temperature_of_soil_layer_4', &
      units='K',transferOffer='will provide'), &
    LIS_Field(stdname='volume_fraction_of_total_water_in_soil', &
      units='m3 m-3',transferOffer='will provide')/)

  
!EOP

contains

  !-----------------------------------------------------------------------------
  ! LIS_NUOPC_Init: Allocate memory and initialize domain and start values
  !-----------------------------------------------------------------------------

!BOP
! !ROUTINE: LIS_NUOPC_Init
!
! !INTERFACE:
  subroutine LIS_NUOPC_Init(vm,rc)
! !ARGUMENTS:
    type(ESMF_VM),intent(in)                :: vm
    integer,intent(out)                     :: rc
!
! !DESCRIPTION:
!  This routine defines the set of steps required from LIS during the
!  initialization of the LIS-NUOPC system.
! \begin{description}
!  \item[LIS\_config\_init] (\ref{LIS_config_init}) \\
!    initialize the LIS configuration
!  \item[LIS\_domain\_init] (\ref{LIS_domain_init}) \\
!    initialize the LIS domains
!  \item[LIS\_createTmnUpdate] (\ref{LIS_createTmnUpdate}) \\
!    no description
!  \item[LIS\_param\_init] (\ref{LIS_param_init}) \\
!    initialize parameters
!  \item[LIS\_perturb\_init] (\ref{LIS_perturb_init}) \\
!    no description
!  \item[LIS\_surfaceModel\_init] (\ref{LIS_surfaceModel_init}) \\
!    initialize the surface model.
!  \item[LIS\_surfaceModel\_init] (\ref{LIS_surfaceModel_init}) \\
!    no description
!  \item[LIS\_metforcing\_init] (\ref{LIS_metforcing_init}) \\
!    initialize the met forcing
!  \item[LIS\_initDAObservations] (\ref{LIS_initDAobservations}) \\
!    initialize structures needed to read observations for
!    data assimilation
!  \item[LIS\_dataassim\_init] (\ref{LIS_dataassim_init}) \\
!    no description
!  \item[LIS\_surfaceModel\_setup] (\ref{LIS_surfaceModel_setup}) \\
!    complete the LSM setups
!  \item[LIS\_surfaceModel\_readrestart] (\ref{LIS_surfaceModel_readrestart}) \\
!    read the restart files
!  \item[LIS\_perturb\_readrestart] (\ref{LIS_perturb_readrestart}) \\
!    no description
!  \item[LIS\_core\_init] (\ref{LIS_core_init}) \\
!    no description
! \end{description}
!
!EOP
!
! ! LOCAL VARIABLES
    integer                                 :: coupled

    rc = ESMF_SUCCESS

    if(.not.LIS_initialized) then
       call LIS_config_init(vm) ! LIS runtime configuration for offline (uncoupled simulation)
       call LIS_domain_init
       call LIS_createTmnUpdate
       call LIS_param_init
       call LIS_perturb_init
       call LIS_surfaceModel_init

       LIS_rc%met_nf(:) = 9 ! WRFout sets to 17

       ! call LIS_metforcing_init(coupled) ! Initialize coupling metforcing subroutines runmode: <coupled>
       call LIS_metforcing_init ! Initialize offline metforcing subroutines runmode: retrospective
       call LIS_initDAObservations
       call LIS_dataassim_init
       call LIS_surfaceModel_setup
       call LIS_surfaceModel_readrestart
       call LIS_perturb_readrestart

       call LISWRF_alloc_states
       call LIS_core_init

       call LIS_HookupInit(rc)
       if(ESMF_STDERRORCHECK(rc)) return ! bail out

       LIS_initialized = .true.
    endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! LIS_NUOPC_Run: Advance nest calculation based on clock
  !-----------------------------------------------------------------------------

!BOP
! !ROUTINE: LIS_NUOPC_Run
!
! !INTERFACE:
  subroutine LIS_NUOPC_Run(nest,mode,slice,importState,exportState,clock,rc)
! !ARGUMENTS:
    integer,intent(in)                     :: nest
    integer,intent(in)                     :: mode
    integer,intent(in)                     :: slice
    type(ESMF_State),intent(inout)         :: importState
    type(ESMF_State),intent(inout)         :: exportState
    type(ESMF_Clock),intent(in)            :: clock
    integer,intent(out)                    :: rc
!
! !DESCRIPTION:
!  This routine defines the set of steps required from LIS during a coupled
!  NUOPC run. The routine translates the NUOPC import fields to NUOPC forcing
!  variables, and then calls the LIS LSM physics routines. Finally the
!  export state from LIS to NUOPC is generated.
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=10)           :: sliceStr
    type(ESMF_Time)             :: currTime
    type(ESMF_Time)             :: stopTime
    type(ESMF_TimeInterval)     :: timeStep
    integer                     :: yy, mm, dd, h, m, s

    rc = ESMF_SUCCESS

    if (slice > 999999999) then
      sliceStr = '999999999+'
    else
      write (sliceStr,"(I0)") slice
    endif

    ! use incoming clock
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    ! LIS time manager expects end of the interval
    stopTime = currTime + timeStep

    ! Confirm if the timemgr should receive current time or stop time
    call ESMF_TimeGet(stopTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    call LIS_timemgr_set(LIS_rc, yy, mm, dd, h, m, s, 0, 0.0)

    select case (mode)
      case (LIS_Offline)
        ! Read in met forcing data from Met forcing sources listed in lis.config
        call LIS_get_met_forcing(nest)
      case (LIS_Coupled)
        ! Copy data from import state
        call LIS_ImportFieldsCopy(nest,importState,rc)
        if(ESMF_STDERRORCHECK(rc)) return ! bail out
      case (LIS_Hybrid)
        ! Read in met forcing data from Met forcing sources listed in lis.config
        call LIS_get_met_forcing(nest)
        ! Copy data from import state
        call LIS_ImportFieldsCopy(nest,importState,rc)
        if(ESMF_STDERRORCHECK(rc)) return ! bail out        
      case default
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg="Running mode is unknown.", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return  ! bail out
    end select

    call LIS_setDynparams(nest)
    call LIS_perturb_forcing(nest)
    call LIS_surfaceModel_f2t(nest)
    call LIS_surfaceModel_run(nest)
    call LIS_surfaceModel_perturb_states(nest)
    call LIS_readDAobservations(nest)
    call LIS_perturb_DAobservations(nest)
    call LIS_dataassim_run(nest)
    call LIS_dataassim_output(nest)

    call LIS_surfaceModel_output(nest)
    call LIS_surfaceModel_writerestart(nest)

    ! =========================================================
    ! Write LIS output data to export state
    ! =========================================================

    ! See LIS_lsmcpl_pluginMod for subroutines registered to: "retrospective", "WRF coupling", "NUOPC coupling"
    call LIS_surfaceModel_setexport(nest)

    call LIS_ExportFieldsCopy(nest,exportState,rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

  end subroutine

  !-----------------------------------------------------------------------------
  ! LIS_NUOPC_Final: Cleanup allocated memory and set endtime.
  !-----------------------------------------------------------------------------

!BOP
! !ROUTINE: LIS_NUOPC_Final
!
! !INTERFACE:
  subroutine LIS_NUOPC_Final(nest,clock,rc)
! !ARGUMENTS:
    integer,intent(in)                      :: nest
    type(ESMF_Clock), intent(in)            :: clock
    integer,intent(out)                     :: rc
!
! !DESCRIPTION:
!  This routine defines the set of steps required from LIS to finalize
!  a coupled LIS-NUOPC run. The routine writes LIS output and restart files
!  valid at the end of the run.
!
!EOP
!
! !LOCAL VARIABLES:
    integer                    :: stat
    type(ESMF_Time)            :: currTime
    integer                    :: yy, mm, dd, h, m, s
    rc = ESMF_SUCCESS

    ! use incoming clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    ! Confirm if the timemgr should receive current time or stop time
    call ESMF_TimeGet(currTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    call LIS_timemgr_set(LIS_rc, yy, mm, dd, h, m, s, 0, 0.0)

    LIS_rc%endtime = 1

  end subroutine

  !-----------------------------------------------------------------------------
  ! LIS_HookupInit: Initialize the field list hookups
  !-----------------------------------------------------------------------------

  subroutine LIS_HookupInit(rc)
    ! ARGUMENTS
    integer,intent(out) :: rc
    ! LOCAL VARIABLES
    integer :: stat
    integer :: fIndex
    integer :: nIndex
    character(len=10)  :: nStr

    rc = ESMF_SUCCESS

    if (LIS_rc%nnest > 999999999) then
      nStr = '999999999+'
    else
      write (nStr,"(I0)") LIS_rc%nnest
    endif

    if (.NOT.allocated(LISWRF_export)) then
      call ESMF_LogSetError(ESMF_RC_OBJ_INIT, &
        msg="LISWRF_export is not initialized.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif

    if (size(LISWRF_export) < LIS_rc%nnest) then
      call ESMF_LogSetError(ESMF_RC_OBJ_INIT, &
        msg="LISWRF_export is not initialized for each nest.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif

    do fIndex=1, size(LIS_FieldList)
      allocate(LIS_FieldList(fIndex)%hookup(LIS_rc%nnest), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Allocation of field hookup memory failed.", &
        line=__LINE__,file=__FILE__)) &
      return  ! bail out
    enddo

    do nIndex=1, LIS_rc%nnest
    do fIndex=1, size(LIS_FieldList)
      select case (trim(LIS_FieldList(fIndex)%stdname))
        case ('aerodynamic_roughness_length')              ! (01)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Z0%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Z0%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Z0%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%znt
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%znt_t
        case ('canopy_moisture_storage')                   ! (02)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%cmc
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%cmc_t
        case ('carbon_dioxide')                            ! (03)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('cosine_zenith_angle')                       ! (04)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Cosz%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Cosz%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Cosz%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('exchange_coefficient_heat')                 ! (05)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Ch%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Ch%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Ch%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('exchange_coefficient_heat_height2m')        ! (06)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Chs2%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Chs2%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Chs2%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%chs2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%chs2_t
        case ('exchange_coefficient_moisture_height2m')    ! (07)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Cqs2%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Cqs2%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Cqs2%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%cqs2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%cqs2_t
        case ('ice_mask')                                  ! (08)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_XICE%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_XICE%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_XICE%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%xice
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%xice_t
        case ('inst_down_lw_flx')                          ! (09)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_LWdown%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_LWdown%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_LWdown%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('inst_down_sw_flx')                          ! (10)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_SWdown%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_SWdown%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_SWdown%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('inst_height_lowest')                        ! (11)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Forc_Hgt%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Forc_Hgt%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Forc_Hgt%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('inst_merid_wind_height_lowest')             ! (12)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Wind_N%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Wind_N%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Wind_N%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('inst_pres_height_lowest')                   ! (13)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('inst_pres_height_surface')                  ! (14)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Psurf%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Psurf%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Psurf%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('inst_spec_humid_height_lowest')             ! (15)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Qair%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Qair%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Qair%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('inst_temp_height_lowest')                   ! (16)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Tair%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Tair%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Tair%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('mean_temp_height_surface')                  ! (17)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%avgsurft
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%avgsurft_t
        case ('inst_wind_speed_height_lowest')             ! (18)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_WIND%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_WIND%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_WIND%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('inst_zonal_wind_height_lowest')             ! (19)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Wind_E%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Wind_E%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Wind_E%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('liquid_water_content_of_soil_layer_1')      ! (20)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%sh2o1
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%sh2o1_t
        case ('liquid_water_content_of_soil_layer_2')      ! (21)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%sh2o2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%sh2o2_t
        case ('liquid_water_content_of_soil_layer_3')      ! (22)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%sh2o3
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%sh2o3_t
        case ('liquid_water_content_of_soil_layer_4')      ! (23)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%sh2o4
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%sh2o4_t
        case ('liquid_water_content_of_surface_snow')      ! (24)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%snow
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%snow_t
        case ('mean_cprec_rate')                           ! (25)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_CRainf%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_CRainf%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_CRainf%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('mean_down_lw_flx')                          ! (26)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('mean_down_sw_flx')                          ! (27)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('mean_fprec_rate')                           ! (28)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Snowf%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Snowf%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Snowf%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('mean_grnd_sensi_heat_flx')                  ! (29)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qg
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qg_t
        case ('mean_laten_heat_flx')                       ! (30)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qle
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qle_t
        case ('mean_laten_heat_flx_atm_into_lnd')          ! (31)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qle
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qle_t
        case ('mean_laten_heat_flx_kinematic')             ! (32)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%eta_kinematic
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%eta_kinematic_t
        case ('mean_net_lw_flx')                           ! (33)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('mean_prec_rate')                            ! (34)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Rainf%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Rainf%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Rainf%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('mean_sensi_heat_flx')                       ! (35)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qh
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qh_t
        case ('mean_sensi_heat_flx_atm_into_lnd')          ! (36)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qh
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qh_t
        case ('mean_surface_albedo')                       ! (37)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Alb%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Alb%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Alb%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%albedo
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%albedo_t
        case ('moisture_content_of_soil_layer_1')          ! (38)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%smc1
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%smc1_t
        case ('moisture_content_of_soil_layer_2')          ! (39)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%smc2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%smc2_t
        case ('moisture_content_of_soil_layer_3')          ! (40)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%smc3
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%smc3_t
        case ('moisture_content_of_soil_layer_4')          ! (41)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%smc4
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%smc4_t
        case ('saturated_mixing_ratio')                    ! (42)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Q2sat%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Q2sat%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Q2sat%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('effective_mixing_ratio')                    ! (43)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%q1
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%q1_t
        case ('soil_temperature_bottom')                   ! (44)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_TMN%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_TMN%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_TMN%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case ('subsurface_runoff_flux')                    ! (45)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qsb
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qsb_t
        case ('surface_microwave_emissivity')              ! (46)
          LIS_FieldList(fIndex)%lisForc=.TRUE.
          if (allocated(LIS_FORC_Emiss%varname)) &
            LIS_FieldList(fIndex)%lisForcVarname=LIS_FORC_Emiss%varname(1)
          LIS_FieldList(fIndex)%lisForcSelect=(LIS_FORC_Emiss%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%emiss
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%emiss_t
        case ('surface_runoff_flux')                       ! (47)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%qs
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%qs_t
        case ('surface_snow_thickness')                    ! (48)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%snowh
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%snowh_t
        case ('fractional_snow_cover')                     ! (49)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%snocvr
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%snocvr_t
        case ('temperature_of_soil_layer_1')               ! (50)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!             LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%stc1
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%stc1_t
        case ('temperature_of_soil_layer_2')               ! (51)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%stc2
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%stc2_t
        case ('temperature_of_soil_layer_3')               ! (52)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%stc3
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%stc3_t
        case ('temperature_of_soil_layer_4')               ! (53)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
          LIS_FieldList(fIndex)%lisExport=.TRUE.
          LIS_FieldList(fIndex)%lisExportSelect=.TRUE.
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%stc4
          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%stc4_t
        case ('volume_fraction_of_total_water_in_soil')    ! (54)
!          LIS_FieldList(fIndex)%lisForc=.FALSE.
!          if (allocated(UNKNOWN%varname)) &
!            LIS_FieldList(fIndex)%lisForcVarname=UNKNOWN%varname(1)
!          LIS_FieldList(fIndex)%lisForcSelect=(UNKNOWN%selectOpt == 1)
!          LIS_FieldList(fIndex)%lisExport=.FALSE.
!          LIS_FieldList(fIndex)%lisExportSelect=.FALSE.
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray=>LISWRF_export(nIndex)%UNKNOWN
!          LIS_FieldList(fIndex)%hookup(nIndex)%exportArray_t=>LISWRF_export(nIndex)%UNKNOWN_T
        case default
          call ESMF_LogWrite("LIS: Field hookup information missing. " //&
            "Skipping hookup: "//trim(LIS_FieldList(fIndex)%stdName), &
            ESMF_LOGMSG_WARNING)
      end select
    enddo
    enddo

  end subroutine

  !-----------------------------------------------------------------------------
  ! Copy NUOPC Import Fields to LIS Forcing State
  !-----------------------------------------------------------------------------

  subroutine LIS_ImportFieldsCopy(nest,importState,rc)
    ! ARGUMENTS
    integer,intent(in)                      :: nest
    type(ESMF_State),intent(inout)          :: importState
    integer,intent(out)                     :: rc
    ! LOCAL VARIABLES
    type(ESMF_StateItem_Flag)  :: itemType
    type(ESMF_StateItem_Flag)  :: lisItemType
    type(ESMF_Field)           :: importField
    type(ESMF_Field)           :: lisImportField
    integer                    :: fIndex
    
    rc = ESMF_SUCCESS

    do fIndex = 1,size(LIS_FieldList)
      if (LIS_FieldList(fIndex)%lisForc) then
        ! Check itemType to see if field exists in import state
        call ESMF_StateGet(importState, &
          itemName=trim(LIS_FieldList(fIndex)%stdname), &
          itemType=itemType, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call LIS_ForcFieldGet(LIS_FieldList(fIndex)%lisForcVarname, &
          nest=nest,itemType=lisItemType,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        if (itemType == ESMF_STATEITEM_FIELD) then
          if (lisItemType == ESMF_STATEITEM_FIELD ) then
!            call ESMF_LogWrite( "LIS: Copying from import state: "// &
!              trim(LIS_FieldList(fIndex)%stdname),ESMF_LOGMSG_INFO)
            call ESMF_StateGet(importState, &
              itemName=trim(LIS_FieldList(fIndex)%stdname), &
              field=importField,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return ! bail out
            call LIS_ForcFieldGet(LIS_FieldList(fIndex)%lisForcVarname, &
              nest=nest,field=lisImportField,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return ! bail out
            call LIS_CopyToLIS(field=importField, &
              fieldLIS=lisImportField, &
              nest=nest,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return ! bail out
          else
            call ESMF_LogWrite( &
              "LIS: Field is not present in LIS_Forc state. Skipping copy: "// &
              trim(LIS_FieldList(fIndex)%stdname),ESMF_LOGMSG_WARNING)
          endif
        else
!          call ESMF_LogWrite( &
!            "LIS: Field is not present in import state. Skipping copy: "// &
!            trim(LIS_FieldList(fIndex)%stdname),ESMF_LOGMSG_WARNING) 
          cycle 
       endif
      endif
    enddo

  end subroutine

  !-----------------------------------------------------------------------------
  ! Copy LIS Export Fields to NUOPC Export State
  !-----------------------------------------------------------------------------

  subroutine LIS_ExportFieldsCopy(nest,exportState,rc)
! !ARGUMENTS:
    integer,intent(in)                     :: nest
    type(ESMF_State),intent(inout)         :: exportState
    integer,intent(out)                    :: rc
!EOP
! !LOCAL VARIABLES:
    type(ESMF_Field)           :: exportField
    integer                    :: fIndex
    type(ESMF_StateItem_Flag)  :: itemType

    rc = ESMF_SUCCESS

    do fIndex = 1,size(LIS_FieldList)
      if (LIS_FieldList(fIndex)%lisExport) then
        ! Check itemType to see if field exists in import state
        call ESMF_StateGet(exportState, &
          itemName=trim(LIS_FieldList(fIndex)%stdname), &
          itemType=itemType, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        if (itemType == ESMF_STATEITEM_FIELD) then
          if (.NOT. allocated(LIS_FieldList(fIndex)%hookup)) then
            call ESMF_LogWrite( "LIS: Export field hookups have not been allocated. Skipping copy: "// &
              trim(LIS_FieldList(fIndex)%stdname),ESMF_LOGMSG_WARNING)
            cycle
          endif
          if (.NOT. associated(LIS_FieldList(fIndex)%hookup(nest)%exportArray_t)) then
            call ESMF_LogWrite( "LIS: Export array is missing. Skipping copy: "// &
              trim(LIS_FieldList(fIndex)%stdname),ESMF_LOGMSG_WARNING)
            cycle
          endif
!          call ESMF_LogWrite( "LIS: Copying to export state: "// &
!            trim(LIS_FieldList(fIndex)%stdname),ESMF_LOGMSG_INFO)
          call ESMF_StateGet(exportState, &
            itemName=trim(LIS_FieldList(fIndex)%stdname), &
            field=exportField,rc=rc)
          if(ESMF_STDERRORCHECK(rc)) return ! bail out
          call LIS_CopyFromLIS( &
            farrayLIS=LIS_FieldList(fIndex)%hookup(nest)%exportArray_t, &
            field=exportField,nest=nest,rc=rc)
          if(ESMF_STDERRORCHECK(rc)) return ! bail out
        else
 !         call ESMF_LogWrite( "LIS: Field is not present in export state. Skipping copy: "// &
 !           trim(LIS_FieldList(fIndex)%stdname),ESMF_LOGMSG_WARNING)
          cycle
        endif
      endif
    enddo
  end subroutine

  !-----------------------------------------------------------------------------

  function LIS_GridCreate(nest,rc)
    !RETURN VALUE
    type(ESMF_Grid) :: LIS_GridCreate
    !ARGUEMENTS
    integer,intent(in)                      :: nest
    integer,intent(out)                     :: rc

    ! LOCAL VARIABLES
    integer                                    :: stat
    integer,allocatable                        :: deBlockList(:,:,:)
    integer,allocatable                        :: petMap(:)
    type(ESMF_DELayout)                        :: deLayout
!   type(ESMF_DistGridConnection), allocatable :: connectionList(:)
    type(ESMF_DistGrid)                        :: distGrid
    character(len=10)   :: did
    integer             :: petID, deID
    integer             :: istart, jstart
    integer             :: id, col, row
    integer             :: col_halo, row_halo, gindex
    real                :: lat_cor, lon_cor
    real                :: lat_cen, lon_cen
    integer             :: lbnd(2),ubnd(2)
    real(ESMF_KIND_RX), pointer    :: coordXcenter(:,:)
    real(ESMF_KIND_RX), pointer    :: coordYcenter(:,:)
    real(ESMF_KIND_RX), pointer    :: coordXcorner(:,:)
    real(ESMF_KIND_RX), pointer    :: coordYcorner(:,:)
    integer(ESMF_KIND_I4), pointer :: gridmask(:,:)
    real(ESMF_KIND_R8), pointer    :: gridarea(:,:)
    integer,dimension(2)           :: halowidth_x(2), halowidth_y(2)

    rc = ESMF_SUCCESS

    ! Setup list of 2-D decompositions
    ! One block per PET
    allocate(deBlockList(2,2,LIS_npes),petMap(LIS_npes),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='Allocation of deBlockList and petMap memory failed.', &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) return ! bail out
    do petID=0, LIS_npes-1
      deID = petID + 1
      deBlockList(1,1,deID) = LIS_ews_ind(nest,deID) ! East-West Start
      deBlockList(1,2,deID) = LIS_ewe_ind(nest,deID) ! East-West End
      deBlockList(2,1,deID) = LIS_nss_ind(nest,deID) ! North-South Start
      deBlockList(2,2,deID) = LIS_nse_ind(nest,deID) ! North-South End
      petMap(deID)=petID
    enddo

    deLayout = ESMF_DELayoutCreate(petMap, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

!    allocate(connectionList(1),stat=stat)
!    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
!      msg='Allocation of connectionList memory failed.', &
!      line=__LINE__, file=FILENAME, rcToReturn=rc)) return ! bail out
!    call ESMF_DistGridConnectionSet(LIS_ConnectionList(1), tileIndexA=1, &
!      tileIndexB=1, positionVector=(/LIS_rc%gnc, 0/), rc=rc)
!    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    distGrid = ESMF_DistGridCreate( &
      minIndex=(/1,1/), maxIndex=(/LIS_rc%gnc(nest),LIS_rc%gnr(nest)/), &
      indexflag = ESMF_INDEX_DELOCAL, &
      deBlockList=deBlockList, &
      delayout=deLayout, &
!     connectionList=LIS_ConnectionList, &
      rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call LIS_DecompGet(distgrid,istart=istart,jstart=jstart,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

!    deallocate(LIS_ConnectionList,stat=stat)
!    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
!      msg='Deallocation of connection list memory failed.', &
!      line=__LINE__,file=FILENAME,rcToReturn=rc)) return ! bail out

    deallocate(deBlockList,petMap,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='Deallocation of deBlockList and petMap memory failed.', &
      line=__LINE__,file=FILENAME,rcToReturn=rc)) return ! bail out

    write(did,"(I0)") nest
    LIS_GridCreate = ESMF_GridCreate(name='LIS_Grid_'//trim(did), distgrid=distgrid, &
      coordSys = ESMF_COORDSYS_SPH_DEG, &
      coordTypeKind = ESMF_TYPEKIND_RX, &
      gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
      rc = rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! Add Center Latitude & Longitude Coordinates
    call ESMF_GridAddCoord(LIS_GridCreate, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! Get pointer to Center Latitude array
    call ESMF_GridGetCoord(LIS_GridCreate, coordDim=1, &
      staggerloc=ESMF_STAGGERLOC_CENTER, localDE=0, &
      computationalLBound=lbnd, computationalUBound=ubnd, &
      farrayPtr=coordXcenter, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    ! Get pointer to Center Longitude array
    call ESMF_GridGetCoord(LIS_GridCreate, coordDim=2, &
      staggerloc=ESMF_STAGGERLOC_CENTER, localDE=0, &
      farrayPtr=coordYcenter, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! error checking
    ! lnc_red and lnr_red are "reduced" to not include the halo region
    if ((ubnd(1)-lbnd(1)+1 /= LIS_rc%lnc_red(nest)) .OR. &
    (ubnd(2)-lbnd(2)+1 /= LIS_rc%lnr_red(nest))) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg="Size of local coordinate arrays not equal to size of local DE", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif

    call NUOPC_NetcdfReadIXJX("lon",trim(LIS_rc%paramfile(nest)),(/istart,jstart/),coordXcenter,rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call NUOPC_NetcdfReadIXJX("lat",trim(LIS_rc%paramfile(nest)),(/istart,jstart/),coordYcenter,rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

!    ! determine halo widths
!    halowidth_x(1) = LIS_ews_ind(nest,LIS_localPet+1) - LIS_ews_halo_ind(nest,LIS_localPet+1)
!    halowidth_x(2) = LIS_ewe_halo_ind(nest,LIS_localPet+1) - LIS_ewe_ind(nest,LIS_localPet+1)
!    halowidth_y(1) = LIS_nss_ind(nest,LIS_localPet+1) - LIS_nss_halo_ind(nest,LIS_localPet+1)
!    halowidth_y(2) = LIS_nse_halo_ind(nest,LIS_localPet+1) - LIS_nse_ind(nest,LIS_localPet+1)
!
!    ! coordinate indices are local, starting at 1
!    ! exlude halo region
!    do row=lbnd(2), ubnd(2)
!    do col=lbnd(1), ubnd(1)
!      col_halo = halowidth_x(1) + col ! Map column to grid with halo
!      row_halo = halowidth_y(1) + row ! Map row to grid with halo
!      id = col_halo+(row_halo-1)*LIS_rc%lnc(nest)! Map to local 1-D grid array index
!
!      coordXCenter(col,row) = LIS_domain(nest)%lon(id)
!      coordYCenter(col,row) = LIS_domain(nest)%lat(id)
!    enddo
!    enddo

    ! Add Grid Mask
    call ESMF_GridAddItem(LIS_GridCreate, itemFlag=ESMF_GRIDITEM_MASK, itemTypeKind=ESMF_TYPEKIND_I4, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    ! Get pointer to Grid Mask array
    call ESMF_GridGetItem(LIS_GridCreate, itemflag=ESMF_GRIDITEM_MASK, localDE=0, &
      staggerloc=ESMF_STAGGERLOC_CENTER, &
      farrayPtr=gridmask, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call NUOPC_NetcdfReadIXJX("LANDMASK",trim(LIS_rc%paramfile(nest)),(/istart,jstart/),gridmask,rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

!    ! Add Grid Area
!    call ESMF_GridAddItem(LIS_GridCreate, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, &
!       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
!    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
!    ! Get pointer to Grid Area array
!    call ESMF_GridGetItem(LIS_GridCreate, itemflag=ESMF_GRIDITEM_AREA, localDE=0, &
!      staggerloc=ESMF_STAGGERLOC_CENTER, &
!      farrayPtr=gridarea, rc=rc)
!    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

  end function

  !-----------------------------------------------------------------------------

  subroutine LIS_DecompGet(distgrid,istart,iend,jstart,jend,rc)
    type(ESMF_DistGrid),intent(in)          :: distGrid
    integer,intent(out),optional            :: istart,iend,jstart,jend
    integer,intent(out)                     :: rc

    ! LOCAL VARIABLES
    integer                :: stat
    integer,allocatable    :: indexCountPDe(:,:)
    integer,allocatable    :: iIndexList(:), jIndexList(:)

    rc = ESMF_SUCCESS

    !! Get the grid distribution for this pet
    allocate(indexCountPDe(2, 0:(LIS_npes - 1)),stat=stat) ! (dimCount, deCount)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='Allocation of indexCountPDE memory failed.', &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) return ! bail out
    call ESMF_DistGridGet(distgrid, indexCountPDe=indexCountPDe, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    allocate(iIndexList(indexCountPDe(1, LIS_localPet)),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='Allocation of iIndexList memory failed.', &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) return ! bail out
    call ESMF_DistGridGet(distgrid, localDe=0, dim=1, &
      indexList=iIndexList, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    allocate(jIndexList(indexCountPDe(2, LIS_localPet)),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='Allocation of jIndexList memory failed.', &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) return ! bail out
    call ESMF_DistGridGet(distgrid, localDe=0, dim=2, &
      indexList=jIndexList, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (present(istart)) istart = minVal(iIndexList)
    if (present(iend)) iend = maxVal(iIndexList)
    if (present(jstart)) jstart = minVal(jIndexList)
    if (present(jend)) jend = maxVal(jIndexList)

    deallocate(indexCountPDe,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='Deallocation of indexCountPDe memory failed.', &
      line=__LINE__,file=FILENAME,rcToReturn=rc)) return ! bail out
    deallocate(iIndexList,jIndexList,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='Deallocation of iIndexList and jIndexList memory failed.', &
      line=__LINE__,file=FILENAME,rcToReturn=rc)) return ! bail out

  end subroutine

  !-----------------------------------------------------------------------------
  ! Copy Data To/From LIS 1D Array and 2D Array
  !-----------------------------------------------------------------------------

  subroutine LIS_FieldCopyToLisField(field,fieldLIS,nest,rc)
! !ARGUMENTS:
    type(ESMF_Field),intent(in)            :: field
    type(ESMF_Field),intent(inout)         :: fieldLIS
    integer,intent(in)                     :: nest
    integer,intent(out)                    :: rc
! !ARGUMENTS:
    type(ESMF_Array)                        :: array
    type(ESMF_Array)                        :: arrayLIS

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field=field,array=array,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call ESMF_FieldGet(field=fieldLIS,array=arrayLIS,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call LIS_CopyToLIS(array=array,arrayLIS=arrayLIS,nest=nest,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FieldCopyToLisFarray(field,farrayLIS,nest,rc)
! !ARGUMENTS:
    type(ESMF_Field),intent(in)            :: field
    real,intent(inout),pointer             :: farrayLIS(:)
    integer,intent(in)                     :: nest
    integer,intent(out)                    :: rc
! !ARGUMENTS:
    type(ESMF_Array)                        :: array

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field=field,array=array,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call LIS_CopyToLIS(array=array,farrayLIS=farrayLIS,nest=nest,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FieldCopyFromLisField(fieldLIS,field,nest,rc)
! !ARGUMENTS:
    type(ESMF_Field),intent(in)            :: fieldLIS
    type(ESMF_Field),intent(inout)         :: field
    integer,intent(in)                     :: nest
    integer,intent(out)                    :: rc
! !ARGUMENTS:
    type(ESMF_Array)                       :: arrayLIS
    type(ESMF_Array)                       :: array

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field=fieldLIS,array=arrayLIS,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call ESMF_FieldGet(field=field,array=array,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call LIS_CopyFromLIS(arrayLIS=arrayLIS,array=array,nest=nest,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FieldCopyFromLisFarray(farrayLIS,field,nest,rc)
! !ARGUMENTS:
    real,intent(in),pointer                :: farrayLIS(:)
    type(ESMF_Field),intent(inout)         :: field
    integer,intent(in)                     :: nest
    integer,intent(out)                    :: rc
! !ARGUMENTS:
    type(ESMF_Array)                       :: array

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field=field,array=array,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call LIS_CopyFromLIS(farrayLIS=farrayLIS,array=array,nest=nest,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_ArrayCopyToLisArray(array,arrayLIS,nest,rc)
! !ARGUMENTS:
    type(ESMF_Array),intent(in)             :: array
    type(ESMF_Array),intent(inout)          :: arrayLIS
    integer,intent(in)                      :: nest
    integer,intent(out)                     :: rc
! !LOCAL VARIABLES:
    integer                         :: localDeCount, localDeCountLIS
    type(ESMF_TypeKind_Flag)        :: typekind, typekindLIS
    integer                         :: rank,rankLIS
    integer                         :: deIndex
    integer(ESMF_KIND_I4),pointer   :: farrayLIS_I4(:)
    integer(ESMF_KIND_I8),pointer   :: farrayLIS_I8(:)
    real(ESMF_KIND_R4),pointer      :: farrayLIS_R4(:)
    real(ESMF_KIND_R8),pointer      :: farrayLIS_R8(:)
    integer(ESMF_KIND_I4),pointer   :: farray_I4(:,:)
    integer(ESMF_KIND_I8),pointer   :: farray_I8(:,:)
    real(ESMF_KIND_R4),pointer      :: farray_R4(:,:)
    real(ESMF_KIND_R8),pointer      :: farray_R8(:,:)
!
! !DESCRIPTION:
!  This routine copies from the cap to LIS or from LIS to the cap
!
!EOP
    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(arrayLIS,typekind=typekindLIS,rank=rankLIS, &
      localDeCount=localDeCountLIS,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_ArrayGet(array,typekind=typekind,rank=rank, &
      localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (rank /= 2) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. Array is not a 2D array.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif
    if (rankLIS /= 1) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. LIS array is not a 1D tile array.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif
    if (localDeCount /= localDeCountLIS) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. LIS array does not match array decomposition.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif
    if (typekind /= typekindLIS) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. LIS array typekind does not match array typekind.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif

    if(typekind==ESMF_TYPEKIND_I4) then
      call ESMF_ArrayGet(array,farrayPtr=farray_I4,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call ESMF_ArrayGet(arrayLIS,farrayPtr=farrayLIS_I4,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyToLIS(farray=farray_I4,farrayLIS=farrayLIS_I4,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    elseif(typekind==ESMF_TYPEKIND_I8) then
      call ESMF_ArrayGet(array,farrayPtr=farray_I8,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call ESMF_ArrayGet(arrayLIS,farrayPtr=farrayLIS_I8,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyToLIS(farray=farray_I8,farrayLIS=farrayLIS_I8,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    elseif(typekind==ESMF_TYPEKIND_R4) then
      call ESMF_ArrayGet(array,farrayPtr=farray_R4,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call ESMF_ArrayGet(arrayLIS,farrayPtr=farrayLIS_R4,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyToLIS(farray=farray_R4,farrayLIS=farrayLIS_R4,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    elseif(typekind==ESMF_TYPEKIND_R8) then
      call ESMF_ArrayGet(array,farrayPtr=farray_R8,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call ESMF_ArrayGet(arrayLIS,farrayPtr=farrayLIS_R8,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyToLIS(farray=farray_R8,farrayLIS=farrayLIS_R8,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    else
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="Typekind copy not implemented.",rcToReturn=rc)
      return
    endif
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_ArrayCopyToLisFarray(array,farrayLIS,nest,rc)
! !ARGUMENTS:
    type(ESMF_Array),intent(in)             :: array
    real,intent(inout),pointer              :: farrayLIS(:)
    integer,intent(in)                      :: nest
    integer,intent(out)                     :: rc
! !LOCAL VARIABLES:
    integer                         :: localDeCount
    type(ESMF_TypeKind_Flag)        :: typekind
    integer                         :: rank
    real(ESMF_KIND_R4),pointer      :: farray_R4(:,:)
    real(ESMF_KIND_R8),pointer      :: farray_R8(:,:)
!
! !DESCRIPTION:
!  This routine copies from the cap to LIS or from LIS to the cap
!
!EOP
    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank, &
      localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (rank /= 2) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. Array is not a 2D array.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif
    if (localDeCount /= 1) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. Array local decomposition count must be 1.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif
    if (typekind /= ESMF_TYPEKIND_R4 .AND. typekind /= ESMF_TYPEKIND_R8) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. LIS array typekind does not match array typekind.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif

    if(typekind==ESMF_TYPEKIND_R4) then
      call ESMF_ArrayGet(array,farrayPtr=farray_R4,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyToLIS(farray=farray_R4,farrayLIS=farrayLIS,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    elseif(typekind==ESMF_TYPEKIND_R8) then
      call ESMF_ArrayGet(array,farrayPtr=farray_R8,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyToLIS(farray=farray_R8,farrayLIS=farrayLIS,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    else
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="Typekind copy not implemented.",rcToReturn=rc)
      return
    endif
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_ArrayCopyFromLisArray(arrayLIS,array,nest,rc)
! !ARGUMENTS:
    type(ESMF_Array),intent(inout)          :: arrayLIS
    type(ESMF_Array),intent(in)             :: array
    integer,intent(in)                      :: nest
    integer,intent(out)                     :: rc
! !LOCAL VARIABLES:
    integer                         :: localDeCount, localDeCountLIS
    type(ESMF_TypeKind_Flag)        :: typekind, typekindLIS
    integer                         :: rank,rankLIS
    integer                         :: deIndex
    integer(ESMF_KIND_I4),pointer   :: farrayLIS_I4(:)
    integer(ESMF_KIND_I8),pointer   :: farrayLIS_I8(:)
    real(ESMF_KIND_R4),pointer      :: farrayLIS_R4(:)
    real(ESMF_KIND_R8),pointer      :: farrayLIS_R8(:)
    integer(ESMF_KIND_I4),pointer   :: farray_I4(:,:)
    integer(ESMF_KIND_I8),pointer   :: farray_I8(:,:)
    real(ESMF_KIND_R4),pointer      :: farray_R4(:,:)
    real(ESMF_KIND_R8),pointer      :: farray_R8(:,:)
!
! !DESCRIPTION:
!  This routine copies from the cap to LIS or from LIS to the cap
!
!EOP
    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(arrayLIS,typekind=typekindLIS,rank=rankLIS, &
      localDeCount=localDeCountLIS,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_ArrayGet(array,typekind=typekind,rank=rank, &
      localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (rank /= 2) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. Array is not a 2D array.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif
    if (rankLIS /= 1) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. LIS array is not a 1D tile array.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif
    if (localDeCount /= localDeCountLIS) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. LIS array does not match array decomposition.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif
    if (typekind /= typekindLIS) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. LIS array typekind does not match array typekind.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif

    if(typekind==ESMF_TYPEKIND_I4) then
      call ESMF_ArrayGet(array,farrayPtr=farray_I4,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call ESMF_ArrayGet(arrayLIS,farrayPtr=farrayLIS_I4,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyFromLIS(farrayLIS=farrayLIS_I4,farray=farray_I4,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    elseif(typekind==ESMF_TYPEKIND_I8) then
      call ESMF_ArrayGet(array,farrayPtr=farray_I8,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call ESMF_ArrayGet(arrayLIS,farrayPtr=farrayLIS_I8,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyFromLIS(farrayLIS=farrayLIS_I8,farray=farray_I8,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    elseif(typekind==ESMF_TYPEKIND_R4) then
      call ESMF_ArrayGet(array,farrayPtr=farray_R4,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call ESMF_ArrayGet(arrayLIS,farrayPtr=farrayLIS_R4,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyFromLIS(farrayLIS=farrayLIS_R4,farray=farray_R4,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    elseif(typekind==ESMF_TYPEKIND_R8) then
      call ESMF_ArrayGet(array,farrayPtr=farray_R8,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call ESMF_ArrayGet(arrayLIS,farrayPtr=farrayLIS_R8,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyFromLIS(farrayLIS=farrayLIS_R8,farray=farray_R8,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    else
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="Typekind copy not implemented.",rcToReturn=rc)
      return
    endif
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_ArrayCopyFromLisFarray(farrayLIS,array,nest,rc)
! !ARGUMENTS:
    real,intent(in),pointer                 :: farrayLIS(:)
    type(ESMF_Array),intent(inout)          :: array
    integer,intent(in)                      :: nest
    integer,intent(out)                     :: rc
! !LOCAL VARIABLES:
    integer                         :: localDeCount
    type(ESMF_TypeKind_Flag)        :: typekind
    integer                         :: rank
    real(ESMF_KIND_R4),pointer      :: farray_R4(:,:)
    real(ESMF_KIND_R8),pointer      :: farray_R8(:,:)
!
! !DESCRIPTION:
!  This routine copies from the cap to LIS or from LIS to the cap
!
!EOP
    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank, &
      localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (rank /= 2) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. Array is not a 2D array.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif
    if (localDeCount /= 1) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Local array decomposition count must be 1.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif
    if (typekind /= ESMF_TYPEKIND_R4 .AND. typekind /= ESMF_TYPEKIND_R8) then
      call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
        msg="Cannot copy. Array typekind does not match array typekind.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif

    if(typekind==ESMF_TYPEKIND_R4) then
      call ESMF_ArrayGet(array,farrayPtr=farray_R4,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyFromLIS(farrayLIS=farrayLIS,farray=farray_R4,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    elseif(typekind==ESMF_TYPEKIND_R8) then
      call ESMF_ArrayGet(array,farrayPtr=farray_R8,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
      call LIS_CopyFromLIS(farrayLIS=farrayLIS,farray=farray_R8,nest=nest,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    else
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="Typekind copy not implemented.",rcToReturn=rc)
      return
    endif
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FarrayI4CopyToLisFarrayI4(farray,farrayLIS,nest,rc)
! !ARGUMENTS:
    integer(ESMF_KIND_I4),intent(in),pointer    :: farray(:,:)
    integer(ESMF_KIND_I4),intent(inout),pointer :: farrayLIS(:)
    integer,intent(in)                          :: nest
    integer,intent(out)                         :: rc
! !LOCAL VARIABLES:
    integer                         :: tile, col, row
! !DESCRIPTION:
!  This routine copies from a 2D array to an LIS 1D array
!EOP
    rc = ESMF_SUCCESS
    do tile=1,LIS_rc%ntiles(nest)
      col = LIS_domain(nest)%tile(tile)%col
      row = LIS_domain(nest)%tile(tile)%row
      farrayLIS(tile) = farray(col,row)
    enddo
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FarrayI8CopyToLisFarrayI8(farray,farrayLIS,nest,rc)
! !ARGUMENTS:
    integer(ESMF_KIND_I8),intent(in),pointer    :: farray(:,:)
    integer(ESMF_KIND_I8),intent(inout),pointer :: farrayLIS(:)
    integer,intent(in)                          :: nest
    integer,intent(out)                         :: rc
! !LOCAL VARIABLES:
    integer                         :: tile, col, row
! !DESCRIPTION:
!  This routine copies from a 2D array to an LIS 1D array
!EOP
    rc = ESMF_SUCCESS
    do tile=1,LIS_rc%ntiles(nest)
      col = LIS_domain(nest)%tile(tile)%col
      row = LIS_domain(nest)%tile(tile)%row
      farrayLIS(tile) = farray(col,row)
    enddo
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FarrayR8CopyToLisFarrayR4(farray,farrayLIS,nest,rc)
! !ARGUMENTS:
    real(ESMF_KIND_R8),intent(in),pointer    :: farray(:,:)
    real(ESMF_KIND_R4),intent(inout),pointer :: farrayLIS(:)
    integer,intent(in)                       :: nest
    integer,intent(out)                      :: rc
! !LOCAL VARIABLES:
    integer                         :: tile, col, row
! !DESCRIPTION:
!  This routine copies from a 2D array to an LIS 1D array
!EOP
    rc = ESMF_SUCCESS
    do tile=1,LIS_rc%ntiles(nest)
      col = LIS_domain(nest)%tile(tile)%col
      row = LIS_domain(nest)%tile(tile)%row
      farrayLIS(tile) = farray(col,row)
    enddo
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FarrayR4CopyToLisFarrayR4(farray,farrayLIS,nest,rc)
! !ARGUMENTS:
    real(ESMF_KIND_R4),intent(in),pointer    :: farray(:,:)
    real(ESMF_KIND_R4),intent(inout),pointer :: farrayLIS(:)
    integer,intent(in)                       :: nest
    integer,intent(out)                      :: rc
! !LOCAL VARIABLES:
    integer                         :: tile, col, row
! !DESCRIPTION:
!  This routine copies from a 2D array to an LIS 1D array
!EOP
    rc = ESMF_SUCCESS
    do tile=1,LIS_rc%ntiles(nest)
      col = LIS_domain(nest)%tile(tile)%col
      row = LIS_domain(nest)%tile(tile)%row
      farrayLIS(tile) = farray(col,row)
    enddo
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FarrayR8CopyToLisFarrayR8(farray,farrayLIS,nest,rc)
! !ARGUMENTS:
    real(ESMF_KIND_R8),intent(in),pointer    :: farray(:,:)
    real(ESMF_KIND_R8),intent(inout),pointer :: farrayLIS(:)
    integer,intent(in)                       :: nest
    integer,intent(out)                      :: rc
! !LOCAL VARIABLES:
    integer                         :: tile, col, row
! !DESCRIPTION:
!  This routine copies from a 2D array to an LIS 1D array
!EOP
    rc = ESMF_SUCCESS
    do tile=1,LIS_rc%ntiles(nest)
      col = LIS_domain(nest)%tile(tile)%col
      row = LIS_domain(nest)%tile(tile)%row
      farrayLIS(tile) = farray(col,row)
    enddo
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FarrayI4CopyFromLisFarrayI4(farrayLIS,farray,nest,rc)
! !ARGUMENTS:
    integer(ESMF_KIND_I4),intent(in),pointer    :: farrayLIS(:)
    integer(ESMF_KIND_I4),intent(inout),pointer :: farray(:,:)
    integer,intent(in)                          :: nest
    integer,intent(out)                         :: rc
! !LOCAL VARIABLES:
    integer                         :: tile, col, row
! !DESCRIPTION:
!  This routine copies from an LIS 1D array to a 2D array
!EOP
    rc = ESMF_SUCCESS
    do row=1,LIS_rc%lnr(nest)
    do col=1,LIS_rc%lnc(nest)
      tile = LIS_domain(nest)%gindex(col,row)
      if(LIS_domain(nest)%gindex(col,row).ne.-1) then
        farray(col,row) = farrayLIS(tile)
      else
        farray(col,row) = UNMASKED
      endif
    enddo
    enddo
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FarrayI8CopyFromLisFarrayI8(farrayLIS,farray,nest,rc)
! !ARGUMENTS:
    integer(ESMF_KIND_I8),intent(in),pointer    :: farrayLIS(:)
    integer(ESMF_KIND_I8),intent(inout),pointer :: farray(:,:)
    integer,intent(in)                          :: nest
    integer,intent(out)                         :: rc
! !LOCAL VARIABLES:
    integer                         :: tile, col, row
! !DESCRIPTION:
!  This routine copies from an LIS 1D array to a 2D array
!EOP
    rc = ESMF_SUCCESS
    do row=1,LIS_rc%lnr(nest)
    do col=1,LIS_rc%lnc(nest)
      tile = LIS_domain(nest)%gindex(col,row)
      if(LIS_domain(nest)%gindex(col,row).ne.-1) then
        farray(col,row) = farrayLIS(tile)
      else
        farray(col,row) = UNMASKED
      endif
    enddo
    enddo
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FarrayR8CopyFromLisFarrayR4(farrayLIS,farray,nest,rc)
! !ARGUMENTS:
    real(ESMF_KIND_R4),intent(in),pointer       :: farrayLIS(:)
    real(ESMF_KIND_R8),intent(inout),pointer    :: farray(:,:)
    integer,intent(in)                          :: nest
    integer,intent(out)                         :: rc
! !LOCAL VARIABLES:
    integer                         :: tile, col, row
! !DESCRIPTION:
!  This routine copies from an LIS 1D array to a 2D array
!EOP
    rc = ESMF_SUCCESS
    do row=1,LIS_rc%lnr(nest)
    do col=1,LIS_rc%lnc(nest)
      tile = LIS_domain(nest)%gindex(col,row)
      if(LIS_domain(nest)%gindex(col,row).ne.-1) then
        farray(col,row) = farrayLIS(tile)
      else
        farray(col,row) = UNMASKED
      endif
    enddo
    enddo
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FarrayR4CopyFromLisFarrayR4(farrayLIS,farray,nest,rc)
! !ARGUMENTS:
    real(ESMF_KIND_R4),intent(in),pointer       :: farrayLIS(:)
    real(ESMF_KIND_R4),intent(inout),pointer    :: farray(:,:)
    integer,intent(in)                          :: nest
    integer,intent(out)                         :: rc
! !LOCAL VARIABLES:
    integer                         :: tile, col, row
! !DESCRIPTION:
!  This routine copies from an LIS 1D array to a 2D array
!EOP
    rc = ESMF_SUCCESS
    do row=1,LIS_rc%lnr(nest)
    do col=1,LIS_rc%lnc(nest)
      tile = LIS_domain(nest)%gindex(col,row)
      if(LIS_domain(nest)%gindex(col,row).ne.-1) then
        farray(col,row) = farrayLIS(tile)
      else
        farray(col,row) = UNMASKED
      endif
    enddo
    enddo
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FarrayR8CopyFromLisFarrayR8(farrayLIS,farray,nest,rc)
! !ARGUMENTS:
    real(ESMF_KIND_R8),intent(in),pointer       :: farrayLIS(:)
    real(ESMF_KIND_R8),intent(inout),pointer    :: farray(:,:)
    integer,intent(in)                          :: nest
    integer,intent(out)                         :: rc
! !LOCAL VARIABLES:
    integer                         :: tile, col, row
! !DESCRIPTION:
!  This routine copies from an LIS 1D array to a 2D array
!EOP
    rc = ESMF_SUCCESS
    do row=1,LIS_rc%lnr(nest)
    do col=1,LIS_rc%lnc(nest)
      tile = LIS_domain(nest)%gindex(col,row)
      if(LIS_domain(nest)%gindex(col,row).ne.-1) then
        farray(col,row) = farrayLIS(tile)
      else
        farray(col,row) = UNMASKED
      endif
    enddo
    enddo
  end subroutine

  !-----------------------------------------------------------------------------
  ! Retrieve ESMF Field from LIS_FORC_State
  !-----------------------------------------------------------------------------

  !BOP
! !ROUTINE: LIS_ForcFieldGet(varname,nest,field,itemType,rc)
!
! !INTERFACE:
  subroutine LIS_ForcFieldGet(varname,nest,field,itemType,rc)
! !ARGUMENTS:
    character(len=*),intent(in)                     :: varname
    integer,intent(in)                              :: nest
    type(ESMF_Field),intent(out),optional           :: field
    type(ESMF_StateItem_Flag),intent(out),optional  :: itemType
    integer,intent(out)                             :: rc
! !LOCAL VARIABLES:
    type(ESMF_StateItem_Flag)                       :: l_itemType
! !DESCRIPTION:
!   Return ESMF_Field for forcing variable
!
!EOP

    rc = ESMF_SUCCESS

    if (nest > size(LIS_FORC_State)) then
      if (present(itemType)) itemType = ESMF_STATEITEM_NOTFOUND
      return
    endif

    call ESMF_StateGet(LIS_FORC_State(nest), &
      itemName=trim(varname), &
      itemType=l_itemType, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    if (present(itemType)) itemType = l_itemType

    if (l_itemType == ESMF_STATEITEM_FIELD) then
      if (present(field)) then
        call ESMF_StateGet(LIS_FORC_State(nest), &
          itemName=trim(varname), &
          field=field,rc=rc)
        if(ESMF_STDERRORCHECK(rc)) return ! bail out
      endif 
    endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! Retrieve timestep of nest from LIS_rc
  !-----------------------------------------------------------------------------

!BOP
! !FUNCTION: LIS_TimestepGet(nest, rc)
! !INTERFACE:
  function LIS_TimestepGet(nest,rc)
! !RETURN VALUE:
    real               :: LIS_TimestepGet
! !ARGUMENTS:
    integer,intent(in)                     :: nest
    integer,intent(out)                    :: rc
! !DESCRIPTION:
!   Return timestep for LIS nest
!
!EOP
!
! !LOCAL VARIABLES:

    rc = ESMF_SUCCESS
    LIS_TimestepGet = LIS_rc%nts(nest)

  end function

  !-----------------------------------------------------------------------------
  ! Retrieve nest count from LIS_rc
  !-----------------------------------------------------------------------------

!BOP
! !FUNCTION: LIS_NestCntGet(rc)
! !INTERFACE:
  function LIS_NestCntGet(rc)
! !RETURN VALUE:
    integer               :: LIS_NestCntGet
! !ARGUMENTS:
    integer,intent(out)           :: rc
! !DESCRIPTION:
!   Return timestep for LIS nest
!
!EOP
!
! !LOCAL VARIABLES:

    rc = ESMF_SUCCESS
    LIS_NestCntGet = LIS_rc%nnest

  end function

  !-----------------------------------------------------------------------------
  ! Retrieve NUOPC coupling mode
  !-----------------------------------------------------------------------------

!BOP
! !FUNCTION: LIS_RunModeGet(fieldList,rc)
! !INTERFACE:
  function LIS_RunModeGet(fieldList,state,rc)
! !RETURN VALUE:
    integer :: LIS_RunModeGet
! !ARGUMENTS:
    type(LIS_Field),intent(in)      :: fieldList(:)
    type(ESMF_State),intent(in)     :: state
    integer,intent(out)             :: rc
! !DESCRIPTION:
!   Return NUOPC coupling for LIS nest
!
!EOP
!
! !LOCAL VARIABLES:
    integer                   :: fIndex
    integer                   :: selectCount
    integer                   :: connectedCount
    type(ESMF_StateItem_Flag) :: itemType

    rc = ESMF_SUCCESS

    LIS_RunModeGet = LIS_Unknown
    selectCount = 0
    connectedCount = 0

    do fIndex=1, size(fieldList)
      if(fieldList(fIndex)%lisForcSelect) then
        selectCount = selectCount + 1
        ! Check itemType to see if field exists in state
        call ESMF_StateGet(state, &
          itemName=trim(fieldList(fIndex)%stdname), &
          itemType=itemType, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        if (itemType == ESMF_STATEITEM_FIELD) then
          if (NUOPC_IsConnected(state, fieldName=trim(fieldList(fIndex)%stdname))) then
            connectedCount = connectedCount + 1
          endif
        endif
      endif
    enddo

    if( connectedCount == 0 ) then
      LIS_RunModeGet = LIS_Offline
    elseif ( connectedCount == selectCount ) then
      LIS_RunModeGet = LIS_Coupled
    else
      LIS_RunModeGet = LIS_Hybrid
    endif

  end function

  !-----------------------------------------------------------------------------
  ! Dictionary Utility
  !-----------------------------------------------------------------------------

  subroutine LIS_FieldDictionaryAdd(rc)
    ! ARGUMENTS
    integer,intent(out)                     :: rc
    ! LOCAL VARIABLES
    integer                    :: fIndex
    logical                    :: isPresent

    rc = ESMF_SUCCESS

    do fIndex=1,size(LIS_FieldList)
      isPresent = NUOPC_FieldDictionaryHasEntry( &
        trim(LIS_FieldList(fIndex)%stdname), &
        rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      if (.not.isPresent) then
        call NUOPC_FieldDictionaryAddEntry( &
          trim(LIS_FieldList(fIndex)%stdname), &
          trim(LIS_FieldList(fIndex)%units), &
          rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
    enddo

  end subroutine

  !-----------------------------------------------------------------------------
  ! Log Utilities
  !-----------------------------------------------------------------------------

  subroutine LIS_Log(label,rc)
! !INTERFACE:
    character(*),intent(in),optional      :: label
    integer,intent(out)                   :: rc
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    integer                    :: nIndex
    integer                    :: mIndex
    character(len=4)           :: nestStr
    character(ESMF_MAXSTR)     :: logMsg

    rc = ESMF_SUCCESS

    if(present(label)) then
      l_label = trim(label)
    else
      l_label = 'LIS_Log'
    endif

    rc = ESMF_SUCCESS

    ! LIS State
    write (logMsg,"(A,A,L1)") trim(l_label), &
      " Initialized: ",LIS_initialized
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,L1)") trim(l_label), &
      " Master process: ",LIS_masterproc
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Local PET: ",LIS_localPet
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Number of PETs: ",LIS_npes
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    ! Run Config
    write (logMsg,"(A,A,A)") trim(l_label), &
      " Running mode: ",LIS_rc%runmode
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,A)") trim(l_label), &
      " Start code: ",LIS_rc%startcode
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Number of nests: ",LIS_rc%nnest
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,2(I0,A))") trim(l_label), &
      " Number of processors (East-West,North-South): (", &
      LIS_rc%npesx,",",LIS_rc%npesy,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    ! Run time
    write (logMsg,"(A,A,6(I0,A))") trim(l_label), &
      " Current time LIS (yr-mo-da-hr:mn:ss): (", &
      LIS_rc%yr,"-",LIS_rc%mo,"-",LIS_rc%da,"_", &
      LIS_rc%hr,":",LIS_rc%mn,":",LIS_rc%ss,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Number of forcing variables: ", &
      LIS_rc%nf
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Number of meteorological forcing datasets: ", &
      LIS_rc%nmetforc
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)

    if (allocated(LIS_rc%met_nf)) then
      do mIndex=1,size(LIS_rc%met_nf)
        write (logMsg,"(A,A,I0,A,I0)") trim(l_label), &
          " Number of forcing variables in met forcing dataset (", &
          mIndex,"): ",LIS_rc%met_nf(mIndex)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
      enddo    
    endif

    if (allocated(LIS_rc%metforc)) then
      do mIndex=1,size(LIS_rc%metforc)
        write (logMsg,"(A,A,I0,A,A)") trim(l_label), &
          " Met forcing (", &
          mIndex,"): ",trim(LIS_rc%metforc(mIndex))
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
      enddo
    endif

    do nIndex=1,LIS_rc%nnest
      if (nIndex > 999999999) then
        write (nestStr,"(A)") '999999999+'
      else
        write (nestStr,"(I0)") nIndex
      endif
      call LIS_LogNest(nIndex,trim(l_label)//" Nest: "//trim(nestStr),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    enddo

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_LogNest(nest,label,rc)
! !INTERFACE:
    integer,intent(in)                    :: nest
    character(*),intent(in),optional      :: label
    integer,intent(out)                   :: rc
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    character(len=10)          :: nestStr
    integer                    :: lbounds(3), ubounds(3)
    character(ESMF_MAXSTR)     :: logMsg

    rc = ESMF_SUCCESS

    if(present(label)) then
      l_label = trim(label)
    else
      if (nest > 999999999) then
        write (nestStr,"(A)") '999999999+'
      else
        write (nestStr,"(I0)") nest
      endif
      l_label = 'LIS_LogNest_'//trim(nestStr)
    endif

    write (logMsg,"(A,A,2(I0,A))") trim(l_label), &
      " Grid dimension (East-West,North-South)=(", &
      LIS_rc%gnc(nest),",",LIS_rc%gnr(nest),")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,I0)") trim(l_label), &
      " Size of gridspace: ",LIS_rc%glbngrid(nest)
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,F0.5)") trim(l_label), &
      " Timestep for nest: ",LIS_rc%nts(nest)
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_LogDomainForNest(nest,label,rc)
! !INTERFACE:
    integer,intent(in)                    :: nest
    character(*),intent(in),optional      :: label
    integer,intent(out)                   :: rc
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    character(len=10)          :: nestStr
    integer                    :: lbounds(3), ubounds(3)
    character(ESMF_MAXSTR)     :: logMsg

    rc = ESMF_SUCCESS

    if(present(label)) then
      l_label = trim(label)
    else
      if (nest > 999999999) then
        write (nestStr,"(A)") '999999999+'
      else
        write (nestStr,"(I0)") nest
      endif
      l_label = 'LIS_LogDomainForNest_'//trim(nestStr)
    endif

    ! LIS Domain
    lbounds(1) = lbound(LIS_domain(nest)%lat,1)
    ubounds(1) = ubound(LIS_domain(nest)%lat,1)
    write (logMsg,"(A,A,2(F0.5,A))") trim(l_label), &
      " Latitude array values=(", &
      LIS_domain(nest)%lat(lbounds(1)), ":", &
      LIS_domain(nest)%lat(ubounds(1)), ")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    lbounds(1) = lbound(LIS_domain(nest)%lon,1)
    ubounds(1) = ubound(LIS_domain(nest)%lon,1)
    write (logMsg,"(A,A,2(F0.5,A))") trim(l_label), &
      " Longitude array values=(", &
      LIS_domain(nest)%lon(lbounds(1)), ":", &
      LIS_domain(nest)%lon(ubounds(1)), ")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,2(I0,A))") trim(l_label), &
      " Grid index (min,max)=(", &
      MINVAL(LIS_domain(nest)%gindex),":",MAXVAL(LIS_domain(nest)%gindex),")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,2(F0.5,A))") trim(l_label), &
      " Latitude (min,max)=(", &
      LIS_domain(nest)%minLat,",",LIS_domain(nest)%maxLat,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,2(F0.5,A))") trim(l_label), &
      " Longitude (min,max)=(", &
      LIS_domain(nest)%minLon,",",LIS_domain(nest)%maxLon,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,2(F0.3,A))") trim(l_label), &
      " Domain Resolution (DX,DY)=(", &
      LIS_domain(nest)%dx,",",LIS_domain(nest)%dy,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,F0.3)") trim(l_label), &
      " Lat increment for lat/lon: ",LIS_domain(nest)%lisproj%dlat
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,F0.3)") trim(l_label), &
      " Lon increment for lat/lon: ",LIS_domain(nest)%lisproj%dlon
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_LogGridDescForNest(nest,label,rc)
! !INTERFACE:
    integer,intent(in)                    :: nest
    character(*),intent(in),optional      :: label
    integer,intent(out)                   :: rc
! !LOCAL VARIABLES:
    character(len=64)          :: l_label
    character(len=10)          :: nestStr
    character(ESMF_MAXSTR)     :: logMsg

    rc = ESMF_SUCCESS

    if(present(label)) then
      l_label = trim(label)
    else
      if (nest > 999999999) then
        write (nestStr,"(A)") '999999999+'
      else
        write (nestStr,"(I0)") nest
      endif
      l_label = 'LIS_LogGridDescForNest'//trim(nestStr)
    endif

    ! Grid Description
    select case ( INT(LIS_rc%gridDesc(nest,1)) )
      case (0) ! 0 = indicateds lat/lon projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(0) lat/lon", ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.0)") trim(l_label), &
          " Number of columns in the domain: ",LIS_rc%gridDesc(nest,2)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.0)") trim(l_label), &
           " Number of rows in the domain: ",LIS_rc%gridDesc(nest,3)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Lower left corner latitude: ",LIS_rc%gridDesc(nest,4)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Lower left corner longitude: ",LIS_rc%gridDesc(nest,5)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Upper right corner latitude: ",LIS_rc%gridDesc(nest,7)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Upper right corner longitude: ",LIS_rc%gridDesc(nest,8)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.3)") trim(l_label), &
          " East-West spatial resolution: ",LIS_rc%gridDesc(nest,9)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.3)") trim(l_label), &
          " North-South spatial resolution: ",LIS_rc%gridDesc(nest,10)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
      case (1) ! 1 = indicates mercator projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(1) mercator", ESMF_LOGMSG_INFO)
      case (3) ! 3 = indicates lambert projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(3) lambert", ESMF_LOGMSG_INFO)
      case (4) ! 4 = indicates gaussian projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(4) gaussian", ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.0)") trim(l_label), &
          " Number of columns in the domain: ",LIS_rc%gridDesc(nest,2)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.0)") trim(l_label), &
          " Number of rows in the domain: ",LIS_rc%gridDesc(nest,3)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Lower left corner latitude: ",LIS_rc%gridDesc(nest,4)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Lower left corner longitude: ",LIS_rc%gridDesc(nest,5)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Upper right corner latitude: ",LIS_rc%gridDesc(nest,7)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.5)") trim(l_label), &
          " Upper right corner longitude: ",LIS_rc%gridDesc(nest,8)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.3)") trim(l_label), &
          " East-West spatial resolution: ",LIS_rc%gridDesc(nest,9)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,A,F0.3)") trim(l_label), &
          " Number of latitude circles: ",LIS_rc%gridDesc(nest,10)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
      case (5) ! 5 = indicates polar stereographic projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(5) polar stereographic", ESMF_LOGMSG_INFO)
      case (7) ! 7 = indicates UTM projection
        call ESMF_LogWrite(trim(l_label)//" Grid projection=(7) UTM", ESMF_LOGMSG_INFO)
      case default
        write (logMsg,"(A,A,I0,A)") trim(l_label), &
          " Grid projection=(",LIS_rc%gridDesc(nest,1),") UNKNOWN"
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    end select

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FieldListLog(label,values,rc)
    ! ARGUMENTS
    character(len=*),intent(in),optional    :: label
    logical,intent(in),optional             :: values
    integer,intent(out)                     :: rc
    ! LOCAL VARIABLES
    character(len=64)          :: l_label
    integer                    :: fIndex

    rc = ESMF_SUCCESS

    if(present(label)) then
      l_label = trim(label)
    else
      l_label = 'LIS_FieldListLog'
    endif

    do fIndex=1, size(LIS_FieldList)
      call LIS_FieldLog(LIS_FieldList(fIndex), &
        label=trim(l_label)//" "//trim(LIS_FieldList(fIndex)%stdname), &
        values=values,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out
    enddo

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine LIS_FieldLog(field,label,values,rc)
    ! ARGUMENTS
    type(LIS_Field),intent(in)              :: field
    character(len=*),intent(in),optional    :: label
    logical,intent(in),optional             :: values
    integer,intent(out)                     :: rc
    ! LOCAL VARIABLES
    character(len=64)          :: l_label
    logical                    :: l_values
    integer                    :: nIndex
    character(len=10)          :: nStr
    character(ESMF_MAXSTR)     :: logMsg
    logical                    :: iAssoc
    logical                    :: eAssoc
    logical                    :: tAssoc
    type(ESMF_StateItem_Flag)  :: itemType
    type(ESMF_Field)           :: importField

    rc = ESMF_SUCCESS

    if(present(label)) then
      l_label = trim(label)
    else
      l_label = 'LIS_FieldLog '//trim(field%stdname)
    endif
    if(present(values)) then
      l_values = values
    else
      l_values = .FALSE.
    endif

    write (logMsg, "(A,(A,A))") trim(l_label), &
      ' Transfer offer: ',trim(field%transferOffer)
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(l_label), &
      ' Forcing available: ',field%lisForc
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(l_label), &
      ' Forcing selected: ',field%lisForcSelect
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(l_label), &
      ' Export available: ',field%lisExport
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(l_label), &
      ' Export selected: ',field%lisExportSelect
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)

    if (allocated(field%hookup)) then
      do nIndex=1, size(field%hookup)
        if (nIndex > 999999999) then
          write (nStr,"(A)") '999999999+'
        else
          write (nStr,"(I0)") nIndex
        endif
        call LIS_ForcFieldGet(field%lisForcVarname, &
          nest=nIndex,itemType=itemType,rc=rc)
        if(ESMF_STDERRORCHECK(rc)) return ! bail out
        if (itemType == ESMF_STATEITEM_FIELD) then
          iAssoc = .TRUE.
        else
          iAssoc = .FALSE.
        endif
        eAssoc = associated(field%hookup(nIndex)%exportArray)
        tAssoc = associated(field%hookup(nIndex)%exportArray_t)
        write (logMsg,"(A,(A,A),(A,L1))") trim(l_label), &
          " Nest: ",trim(nStr), &
          " Import hookup: ",iAssoc
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,(A,A),(A,L1))") trim(l_label), &
          " Nest: ",trim(nStr), &
          " Export hookup_2D: ",eAssoc
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        write (logMsg,"(A,(A,A),(A,L1))") trim(l_label), &
          " Nest: ",trim(nStr), &
          " Export hookup_1D: ",tAssoc
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
        if (l_values) then
          if (iAssoc) then
            call LIS_ForcFieldGet(field%lisForcVarname, &
              nest=nIndex,field=importField,rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return ! bail out
            call NUOPC_LogFieldValue(importField, &
              label=trim(l_label)//" Nest: "//trim(nStr)//" import",rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return ! bail out
          endif
          if (eAssoc) then
            call NUOPC_LogFarrayValue(field%hookup(nIndex)%exportArray, &
              label=trim(l_label)//" Nest: "//trim(nStr)//" export2D",rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return ! bail out
          endif
          if (tAssoc) then
            call NUOPC_LogFarrayValue(field%hookup(nIndex)%exportArray_t, &
              label=trim(l_label)//" Nest: "//trim(nStr)//" export1D",rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return ! bail out
          endif
        endif
      enddo
    else
      write (logMsg,"(A,A)") trim(l_label)//" ", &
        field%stdname//" does not have field hookups."
      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO, &
          line=__LINE__, file=FILENAME)
    endif

  end subroutine

end module
