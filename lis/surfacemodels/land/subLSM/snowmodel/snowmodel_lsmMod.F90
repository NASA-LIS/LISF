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
module snowmodel_lsmMod
!BOP
!
! !MODULE: snowmodel_lsmMod
!
! !DESCRIPTION:
!  
! This module provides the definition of derived data type used to
! control the operation of SnowModel. It also provides the entry
! method for the initialization of SnowModel-specific variables.
! The derived data type {\tt snowmodel\_struc} includes the variables
! that specify the runtime options and other control variables as
! described below:
!
! \begin{description}
!  \item[rfile]
!    name of the SnowModel restart file
!  \item[rformat]
!    format of restart file (binary or netcdf) for SnowModel
!  \item[vfile]
!    name of the static vegetation parameter table
!  \item[parfile]
!    name of the main input parameter for SnowModel
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[nsnow]
!    number of snow layers
!  \item[nvegp]
!    number of static vegetation parameters in the table
!  \item[inittemp]
!    initial soil temperatures for a cold start run
!  \item[rstInterval]
!   restart writing interval
!  \item[ht_rhobs]
!   reference height of T and q forcing
!  \item[ht_windobs]
!   reference height of u and v forcing
!  \item[lyrthk]
!   thickness of soil layers
!  \item[snowmodel]
!   SnowModel specific variables
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Added G. Liston's SnowModel
!
! !USES:        
  use snowmodel_module
  use snowmodel_vars
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: snowmodel_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: snowmodel_struc
!EOP
  type, public ::  snowmodel_type_dec

     ! LIS related parameters:
     real                                        :: rstInterval
     character(len=LIS_CONST_PATH_LEN)           :: rfile
     character(20)                               :: rformat

     ! LIS-SnowModel related parameters:
     character(len=LIS_CONST_PATH_LEN)           :: parfile
     character(len=LIS_CONST_PATH_LEN)           :: vfile
     real                                        :: ht_windobs
     real                                        :: ht_rhobs
     integer                                     :: call_sm_preproc
     integer                                     :: write_sm_metfields
     character(10)                               :: sm_params_opt
     character(10)                               :: sm_micromet_opt

     integer                                     :: nsnow
     real, allocatable                           :: lyrthk(:)
     real, allocatable                           :: inittemp(:)
     real                                        :: initsnowdepth
     real                                        :: initsnowequiv

     integer                                     :: forc_count

     real                                        :: ts
     integer                                     :: iter

     ! SnowModel parameter - snowmodel.h
!     integer                  :: nx_max
!     integer                  :: ny_max
!     integer                  :: nstns_max
!     integer                  :: nvegtypes
!     integer                  :: max_time_steps
!     integer                  :: max_obs_dates
!     integer                  :: nz_max
!     integer                  :: n_print_vars

     ! SnowModel parameter - snowmodel_vars.h
     integer                  :: max_iter
     integer                  :: nx
     integer                  :: ny
     real                     :: deltax
     real                     :: deltay

     real                     :: run_micromet
     real                     :: run_enbal
     real                     :: run_snowpack
     real                     :: run_snowtran

     real                     :: usum_glb
     real                     :: vsum_glb
     real                     :: windspdflg_glb
     real                     :: wslopemax_glb
     real                     :: curvemax_glb
     real                     :: bsflag_glb

     type(snowmodeldec), allocatable :: sm(:)

  end type snowmodel_type_dec

  type(snowmodel_type_dec), allocatable :: snowmodel_struc(:)
  SAVE

contains
!BOP
!
! !ROUTINE: snowmodel_init
! \label{snowmodel_init}
!
! !INTERFACE:
  subroutine snowmodel_init(eks)
! !USES:
    use ESMF
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
         LIS_update_timestep, LIS_registerAlarm
    use LIS_surfaceModelDataMod
    use LIS_lsmMod
    use snowmodel_inc
    use snowmodel_vars
!
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for
!  SnowModel-specific variables. It also invokes the routine to
!  read the runtime specific options for SnowModel from the
!  configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[snowmodel\_readconfig](\ref{snowmodel_readconfig}) \newline
!    reads the runtime options for SnowModel 
!  \end{description}
!EOP
    implicit none
    integer, intent(in)  :: eks

    integer              :: n,t
    character(3)         :: fnest
    integer              :: status

    type(ESMF_ArraySpec) :: arrspec1
    type(ESMF_Field)     :: sweField, snwdField

! _____________________________________________________

    allocate(snowmodel_struc(LIS_rc%nnest))

    ! Read in LIS config entries for SnowModel:
    call snowmodel_readconfig()

    do n=1,LIS_rc%nnest

      allocate(snowmodel_struc(n)%sm(LIS_rc%npatch(n,LIS_rc%lsm_index)))

      ! Number of snow layers (e.g., 2, 25)
      snowmodel_struc(n)%nsnow = 1   ! Set to 1 for now, can be high as 25

      ! Call Snowmodel read parameter subroutine to read in
      !  primary parameter and input settings to run the model:

      snowmodel_dot_par_fname = trim(snowmodel_struc(n)%parfile)

      if( LIS_localPet==0 ) then
         snowmodel_masterproc = .true.
      else
         snowmodel_masterproc = .false.
      endif

      call READPARAM_CODE(&
       dt,deltax,deltay,Utau_t_flag,&
       subgrid_flag,twolayer_flag,snowmodel_dot_par_fname,&
       bc_flag,curve_len_scale,slopewt,curvewt,ht_windobs,&
       ht_rhobs,ro_snow,snow_d_init_const,const_veg_flag,&
       vegsnowdepth,LIS_rc%lnc(n),LIS_rc%lnr(n),max_iter,met_input_fname,xmn,ymn,&  ! KRA: nx,ny
       iyear_init,imonth_init,iday_init,xhour_init,undef,ifill,&
       iobsint,dn,xlat,i_tair_flag,i_rh_flag,i_wind_flag,&
       i_solar_flag,i_prec_flag,isingle_stn_flag,igrads_metfile,&
       windspd_min,icond_flag,run_micromet,run_enbal,run_snowpack,&
       run_snowtran,topoflag,topoveg_fname,snowtran_output_fname,&
       micromet_output_fname,enbal_output_fname,Utau_t_const,&
       snowpack_output_fname,print_micromet,print_enbal,&
       print_snowpack,print_snowtran,i_longwave_flag,print_user,&
       ascii_topoveg,topo_ascii_fname,veg_ascii_fname,&
       irun_data_assim,lapse_rate_user_flag,&
       iprecip_lapse_rate_user_flag,use_shortwave_obs,&
       use_longwave_obs,use_sfc_pressure_obs,calc_subcanopy_met,&
       sfc_sublim_flag,gap_frac,cloud_frac_factor,&
       albedo_snow_forest,albedo_snow_clearing,albedo_glacier,&
       barnes_lg_domain,n_stns_used,tabler_dir,slope_adjust,&
       lat_solar_flag,UTC_flag,iveg_ht_flag,ihrestart_flag,&
       ihrestart_inc,i_dataassim_loop,tsls_threshold,dz_snow_min,&
       print_multilayer,multilayer_snowpack,max_layers,&
       multilayer_output_fname,izero_snow_date,curve_lg_scale_flag,&
       check_met_data,seaice_run,snowmodel_line_flag,wind_lapse_rate,&
       iprecip_scheme,cf_precip_flag,snowfall_frac,print_inc,&
       output_path_wo_assim,output_path_wi_assim,Tabler_1_flag,&
       Tabler_2_flag,tabler_sfc_path_name,print_var)

      ! Check SnowModel parameter source origins:
      if( snowmodel_struc(n)%sm_params_opt == "LDT" ) then
        write(LIS_logunit,*) "[INFO] Reading in SnowModel LSM parameters from LDT"
        ascii_topoveg = 2.0    ! No file read in by SnowModel; use LDT input
      else
        write(LIS_logunit,*) "[INFO] Reading in SnowModel LSM parameters from snowmodel.par file "
      endif

      ! Check SnowModel Micromet source calls:
      if( snowmodel_struc(n)%sm_micromet_opt == "SnowModel" ) then
         write(LIS_logunit,*) "[INFO] Calling MicroMet routines within SnowModel "
      elseif( snowmodel_struc(n)%sm_micromet_opt == "LIS" ) then
         write(LIS_logunit,*) "[INFO] Calling MicroMet routines from LIS metforcing layer "
      else
         write(LIS_logunit,*) "[ERR] Incorrect option set for "
         write(LIS_logunit,*) "[ERR]  SnowModel MicroMet input source: "
         write(LIS_logunit,*) "[ERR]  See documentation in configs/lis.config.adoc "
         write(LIS_logunit,*) "[ERR]  for option details. "
         call LIS_endrun
      endif

       ! Check inputs to driving SnowModel:
       write(LIS_logunit,*) "[INFO] Reading in SnowModel input options"
       write(LIS_logunit,*) ' Run Micromet   : ',run_micromet
       write(LIS_logunit,*) ' Run EnBal      : ',run_enbal
       write(LIS_logunit,*) ' Run Snowpack   : ',run_snowpack
       write(LIS_logunit,*) ' Run Snowtran   : ',run_snowtran
       write(LIS_logunit,*) ' Run data assimilation: ',irun_data_assim
       write(LIS_logunit,*) ' Data assimilation loop: ',i_dataassim_loop
       write(LIS_logunit,*) ' Undefined value: ',undef
       write(LIS_logunit,*) ''
       write(LIS_logunit,*) "[INFO] SnowModel grid information: "
       write(LIS_logunit,*) "[INFO] Note: SnowModel typically runs on an"
       write(LIS_logunit,*) "  equal-area grid, e.g., UTM or Albers Equal Area"
       write(LIS_logunit,*) ' Number of x cells in grid: ',nx
       write(LIS_logunit,*) ' Number of y cells in grid: ',ny
       write(LIS_logunit,*) ' Increment in x direction(m): ',deltax
       write(LIS_logunit,*) ' Increment in y direction(m): ',deltay
       write(LIS_logunit,*) ' LL_X Location in meters (UTM): ',xmn
       write(LIS_logunit,*) ' LL_Y Location in meters (UTM): ',ymn
       write(LIS_logunit,*) ' '
       write(LIS_logunit,*) "[INFO] SnowModel date/time information: "
       write(LIS_logunit,*) ' Model time step: ',dt
       write(LIS_logunit,*) ' Start year of input data: ',iyear_init
       write(LIS_logunit,*) ' Start month of input data: ',imonth_init
       write(LIS_logunit,*) ' Start day of input data: ',iday_init
       write(LIS_logunit,*) ' Start hour of input data: ',xhour_init
       write(LIS_logunit,*) ' Maximum iterations: ',max_iter
       write(LIS_logunit,*) ' '
       write(LIS_logunit,*) "[INFO] SnowModel veg and topo information: "
       write(LIS_logunit,*) ' Topo and veg input file type: ',ascii_topoveg
       write(LIS_logunit,*) ' Topo and veg file name: ',trim(topoveg_fname)
       write(LIS_logunit,*) ' Constant veg type: ',const_veg_flag
       write(LIS_logunit,*) ' Veg height flag: ',iveg_ht_flag
       write(LIS_logunit,*) ' '
       write(LIS_logunit,*) '[INFO] SnowModel solar radiation calculations: '  
       write(LIS_logunit,*) ' UTC flag(local=0; UTC: GRADS_long=-1,TXT_long=1): ',UTC_flag
       write(LIS_logunit,*) ' Latitude solar flag (same as UTC flag options): ',lat_solar_flag
       write(LIS_logunit,*) ' Domain center latitude (used when lat_solar_flag=0): ',xlat
       ! Input height of wind and rh observations (ENBAL and SNOWTRAN):
       write(LIS_logunit,*) '[INFO] Other SnowModel Inputs (ENBAL AND SNOWTRAN)'
       write(LIS_logunit,*) ' Height of wind obs: ',ht_windobs
       write(LIS_logunit,*) ' Height of rel hum obs: ',ht_rhobs
       write(LIS_logunit,*) ' '
       ! Micromet Model setup::
       write(LIS_logunit,*) '[INFO] MICROMET ENTRIES ... '
       write(LIS_logunit,*) ' Tair output : ',i_tair_flag 
       write(LIS_logunit,*) ' RH output   : ',i_rh_flag 
       write(LIS_logunit,*) ' Wind output : ',i_wind_flag 
       write(LIS_logunit,*) ' SWrad output: ',i_solar_flag 
       write(LIS_logunit,*) ' LWdown output: ',i_longwave_flag 
       write(LIS_logunit,*) ' Prec output : ',i_prec_flag 
       write(LIS_logunit,*) ' '
       write(LIS_logunit,*) ' Met Input Filename: ',trim(met_input_fname)
       write(LIS_logunit,*) ' Run met forcing check: ',check_met_data
       write(LIS_logunit,*) ' MicroMet station flag: ',isingle_stn_flag
       write(LIS_logunit,*) ' MicroMet GrADS Metfile: ',igrads_metfile
       write(LIS_logunit,*) ' Barnes station interp: ',barnes_lg_domain
       write(LIS_logunit,*) ' Barnes nearest stations number: ',n_stns_used
       write(LIS_logunit,*) ' Radius of influence: ',iobsint
       write(LIS_logunit,*) ' Radius of influence - obs interval: ',dn

       write(LIS_logunit,*) ' '
       ! These may be useful and considered:
       write(LIS_logunit,*) ' Run grid as 1d vector (not 2d): ', snowmodel_line_flag

       ! Preprocess, MicroMet, SnowTran3D 
       write(LIS_logunit,*) '[INFO] Wind model curvature entries (MicroMet, SNOWTRAN)'
       write(LIS_logunit,*) ' Wind model curvature (len scale): ',curve_len_scale
       write(LIS_logunit,*) ' Wind model curve weight: ',curvewt
       write(LIS_logunit,*) ' Wind model slope weight: ',slopewt
       write(LIS_logunit,*) ' Wind 2nd len scale curve affect: ',curve_lg_scale_flag

       ! MicroMet option for mininum windspeed threshold
       write(LIS_logunit,*) ' Min. windspeed threshold: ',windspd_min

       ! Default is 0, but if user wants to change lapse-rate values,
       !  they need to do it within micromet_code.f:
       write(LIS_logunit,*) ' Monthly lapse rate values: ',lapse_rate_user_flag

       ! Default is 0, but if users want to change precip LR scaling values,
       !  they need to do it within micromet_code.f:
       write(LIS_logunit,*) ' Precip adjustment factor (precip LR scaling): ',iprecip_lapse_rate_user_flag

       ! Two precip-increase-with-elev schemes; 1-orig, 2-van Pelt scheme:
       write(LIS_logunit,*) ' Precip-increase-elev scheme: ',iprecip_scheme

       ! Define rain-snow fraction of precip water equiv (forcing) calculation scheme:
       write(LIS_logunit,*) ' Rain-snow fraction WE calculation scheme: ',snowfall_frac

       ! Turn on wind lapse-rate factor:
       write(LIS_logunit,*) ' Wind lapse-rate factor: ',wind_lapse_rate

       ! Turn on sub-forest-canopy estimates for wind-speed, solar rad, and LW rad:
       write(LIS_logunit,*) ' Calc Sub-forest-canopy estimates: ',calc_subcanopy_met

       ! Define canopy gap fraction (0-1) for solar radiation reaching snow below canopy:
       write(LIS_logunit,*) ' Canopy gap fraction: ',gap_frac

       ! Cloud fraction factor can be used to reduce cloud-cover fraction, 
       !  like in valley locations.
       write(LIS_logunit,*) ' Cloud fraction factor: ',cloud_frac_factor

       ! Select whether SW radiation obs will be assimilated:
       write(LIS_logunit,*) ' SW radiation obs DA: ',use_shortwave_obs

       ! Select whether LW radiation obs will be assimilated:
       write(LIS_logunit,*) ' LW radiation obs DA: ',use_longwave_obs

       ! Select whether sfc pressure obs will be assimilated:
       write(LIS_logunit,*) ' Surface pressure obs DA: ',use_sfc_pressure_obs

       ! User-defined precip correction factor, either 2d field applied, or constant factor:
       !  default=0, turned off
       write(LIS_logunit,*) ' Precipitation correction factor: ',cf_precip_flag

       write(LIS_logunit,*) '** ENDING MICROMET ENTRIES ... '

       ! EnBAL Model setup::
       write(LIS_logunit,*) ' '
       write(LIS_logunit,*) '[INFO] ENBAL ENTRIES ... '

       ! Surface energy balance calculation accounting for non-0 conduction term
       ! 0-single layer; 1-multilayer snowpack physics
       write(LIS_logunit,*) ' Energy balance calc with non-0 cond. term: ', icond_flag

       ! Define albedo for melting snow cover under forest canopy
       write(LIS_logunit,*) ' ENBAL: Forest snow melt albedo value: ',albedo_snow_forest

       ! Define albedo for melting snow cover in non-forested areas:
       write(LIS_logunit,*) ' ENBAL: Nonforest snow melt albedo value: ',albedo_snow_clearing

       ! Define albedo for glacier surface:
       write(LIS_logunit,*) ' ENBAL: Glacier surface albedo value: ',albedo_glacier
   
       write(LIS_logunit,*) '** ENDING ENBAL ENTRIES ... '

       write(LIS_logunit,*) ' '
       ! SnowPack Model setup::
       write(LIS_logunit,*) '[INFO] SNOWPACK ENTRIES ... '
       
       ! Set static-surface (non-blowing snow) sublimation effect on:
       !  If on, Qle from ENBAL is used to add/remove snow via sublimation:
       write(LIS_logunit,*) ' SNOWPACK: non-blowing snow sublimation effect: ',sfc_sublim_flag    

       ! Turn on multilayer snowpack option:
       write(LIS_logunit,*) ' SNOWPACK: multilayer option: ',multilayer_snowpack

       ! Define time since last snowfall to determine if new snow layer gets created:
       write(LIS_logunit,*) ' SNOWPACK: time since snowfall for creating new snow layer: ',tsls_threshold 

       ! Define the min snow layer thickness (m) for multi-layer snowpack model.
       write(LIS_logunit,*) ' SNOWPACK: min snow layer thickness(m): ',dz_snow_min

       ! Set date to reset snowpack to 0 to not allow glacier formation in multi-year runs:
       !  To disable setting = 999999
       write(LIS_logunit,*) ' SNOWPACK: date for zeroing out snowpack: ',izero_snow_date 

       write(LIS_logunit,*) '** ENDING SNOWPACK ENTRIES ... '

       write(LIS_logunit,*) ' '
       ! SnowTran-3d Model setup::
       write(LIS_logunit,*) '[INFO] SNOWTRAN3D ENTRIES ... '

       ! Threshold surface shear velocity as function of AirTemp and Wspd:
       write(LIS_logunit,*) ' SNOWTRAN: Surface shear velocity flag (T+Wspd): ',Utau_t_flag

       ! Surface shear velocity threshold (when Utau_t_flag=0):
       write(LIS_logunit,*) ' SNOWTRAN: Surface shear velocity threshold: ',Utau_t_const

       ! Flag to turn on the Tabler subgrid algorithm (0-no snow redistribution)
       !  Not appropriate for grid increments > 30 m. If turned on, topoflag = 1.
       write(LIS_logunit,*) ' SNOWTRAN: Tabler subsgrid snow redistribution option: ',subgrid_flag

       ! When Tabler flag on, user can define dominant wind direction (of 8 directions).
       !  Value is set constant over domain area:
       write(LIS_logunit,*) ' SNOWTRAN: Tabler constant wind directory: ',tabler_dir

       ! Adjust slopes of Tabler surfaces using 'slope-adjust' parameter, which
       !  is multiplied by calc. Tabler equilibrium drift sufrace slope.
       write(LIS_logunit,*) ' SNOWTRAN: Tabler slope surface adjustment: ',slope_adjust 

       ! Account for two snow-layer movement:
       write(LIS_logunit,*) ' SNOWTRAN: Snow-layer movement: ',twolayer_flag

       ! Define upwind boundary to have 0 incoming transport flux:
       write(LIS_logunit,*) ' SNOWTRAN: Upwind boundary incoming transport flux: ',bc_flag

       ! Snowtran wo Snowpack -- Provide Snow density:
       write(LIS_logunit,*) ' SNOWTRAN: snow density (when Snowpack off): ',ro_snow

       ! Define initial snow-depth distributions:
       write(LIS_logunit,*) ' SNOWTRAN: Initial snow-depth distribution: ',snow_d_init_const

       ! Flag to set surface topography to ground topography:
       write(LIS_logunit,*) ' SNOWTRAN: Set surface topography to ground topography: ',topoflag

       write(LIS_logunit,*) '** END SNOWTRAN3D ENTRIES ... WILL CHECK WITH OTHER SUBMODELS'
       write(LIS_logunit,*) ' '

       ! SeaIce Model:
       write(LIS_logunit,*) '[INFO] SEAICE ENTRIES ... '
       ! Define whether to simulate snow on sea ice:
       write(LIS_logunit,*) 'SEAICE: Turn on sea-ice option: ',seaice_run
       write(LIS_logunit,*) '** END SEAICE ENTRIES ... '

       write(LIS_logunit,*) ' '
       write(LIS_logunit,*) '[INFO] Output directory information '
       write(LIS_logunit,*) ' OUTPUT: Output directory: ',trim(output_path_wo_assim)
       write(LIS_logunit,*) ' OUTPUT: DA Output directory: ',trim(output_path_wi_assim)
       write(LIS_logunit,*) ' OUTPUT: Write individual files(1=yes): ',print_user
       ! This option below defines the write output increment, etc., 8 = 3 hourly
       write(LIS_logunit,*) ' OUTPUT: Print increment, (times)X(day): ',print_inc
       ! Output print options:
       write(LIS_logunit,*) '[INFO] SnowModel values averaged over the period'
       write(LIS_logunit,*) ' OUTPUT: Tair (degC)                     : ',print_var(1)
       write(LIS_logunit,*) ' OUTPUT: RH (%)                          : ',print_var(2)
       write(LIS_logunit,*) ' OUTPUT: Wind spd (m/s)                  : ',print_var(3)
       write(LIS_logunit,*) ' OUTPUT: Inc solar radiation (W/m2)      : ',print_var(4)
       write(LIS_logunit,*) ' OUTPUT: Inc LW radiation (W/m2)         : ',print_var(5)
       write(LIS_logunit,*) ' OUTPUT: Emitted LW rad (W/m2)           : ',print_var(6)
       write(LIS_logunit,*) ' OUTPUT: albedo (-)                      : ',print_var(7)
       write(LIS_logunit,*) ' OUTPUT: wind direction (deg)            : ',print_var(8)
       write(LIS_logunit,*) ' OUTPUT: Cloud fraction (ave over period): ',print_var(21)
 
       write(LIS_logunit,*) '[INFO] SnowModel values summed over period'
       write(LIS_logunit,*) ' OUTPUT: Wat eqv precip (m/ts)           : ',print_var(9)
       write(LIS_logunit,*) ' OUTPUT: Liq precip (m/ts)               : ',print_var(10)
       write(LIS_logunit,*) ' OUTPUT: Solid precip (m/ts)             : ',print_var(11)
       write(LIS_logunit,*) ' OUTPUT: SWE melt (m)                    : ',print_var(12)
       write(LIS_logunit,*) ' OUTPUT: Surface sublmiation (m)         : ',print_var(13)
       write(LIS_logunit,*) ' OUTPUT: Runoff from snowpack base (m/ts): ',print_var(14)
       write(LIS_logunit,*) ' OUTPUT: SWE melt from glacier ice (m)   : ',print_var(15)
 
       write(LIS_logunit,*) '[INFO] SnowModel values saved at the end of period'
       write(LIS_logunit,*) ' OUTPUT: Snow depth (m)                  : ',print_var(16)
       write(LIS_logunit,*) ' OUTPUT: Snow density (kg/m3)            : ',print_var(17)
       write(LIS_logunit,*) ' OUTPUT: SWE (m)                         : ',print_var(18)
       write(LIS_logunit,*) ' OUTPUT: Summed snow precip year (m)     : ',print_var(19)
       write(LIS_logunit,*) ' OUTPUT: Summed SWE melt for year (m)    : ',print_var(20)

       ! Output print options:
!       write(LIS_logunit,*) ' OUTPUT: Micromet variables: ',print_micromet
!       write(LIS_logunit,*) ' OUTPUT: EnBal variables: ',print_enbal
!       write(LIS_logunit,*) ' OUTPUT: Snowpack variables: ',print_snowpack
!       write(LIS_logunit,*) ' OUTPUT: Snowpack multilayer variables: ',print_multilayer
!       write(LIS_logunit,*) ' OUTPUT: SnowTran variables: ',print_snowtran
!       write(LIS_logunit,*) ' OUTPUT: SnowTran Tabler1 fields: ',Tabler_1_flag
!       write(LIS_logunit,*) ' OUTPUT: SnowTran Tabler2 fields: ',Tabler_2_flag

       ! ------------------------------

       ! SnowModel Main code calls and options for reading in and 
       !  preprocessing necessary inputs to the submodel components:
       ! All below from the original snomodel_main.f ...

       ! This loop runs the correction/data assimilation adjustment
       !  iterations.
       if (ihrestart_flag.ge.0) then
         if (i_dataassim_loop.lt.0.0) then
           i_corr_start = 2
         else
           i_corr_start = 1
         endif
       else
         i_corr_start = 1
       endif

       do icorr_factor_loop=i_corr_start,irun_data_assim+1

       ! Perform the correction (precipitation and melt) factor
       ! calculations.
        if (irun_data_assim.eq.1 .and. icorr_factor_loop.eq.2) then
          write(LIS_logunit,*) "[WARN] No call to 'DATAASSIM_USER' at this time"    
          CALL DATAASSIM_USER(LIS_rc%lnc(n),LIS_rc%lnr(n),icorr_factor_index,&  ! KRA: nx,ny
            corr_factor,max_iter,deltax,deltay,xmn,ymn,nobs_dates,&
            print_inc,iday_init,imonth_init,iyear_init,dt,&
            output_path_wo_assim,xhour_init)
          if (ihrestart_flag.ge.-1) then
          write(LIS_logunit,*) "[WARN] No call to 'HRESTART_SAVE_DA' at this time"    
            CALL HRESTART_SAVE_DA(LIS_rc%lnc(n),LIS_rc%lnr(n),max_iter,corr_factor,& ! KRA: nx,ny
              icorr_factor_index,nobs_dates)
          endif
        endif

       end do

       ! ------------------------------

       ! Passing the required SnowModel fields to LIS for what LIS requires ...
       snowmodel_struc(n)%iter = 0

       ! Initialize forcing fields:
       snowmodel_struc(n)%forc_count = 0
       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          snowmodel_struc(n)%sm(t)%tair = 0
          snowmodel_struc(n)%sm(t)%qair = 0
          snowmodel_struc(n)%sm(t)%swdown = 0
          snowmodel_struc(n)%sm(t)%lwdown = 0
          snowmodel_struc(n)%sm(t)%uwind = 0
          snowmodel_struc(n)%sm(t)%vwind = 0
          snowmodel_struc(n)%sm(t)%psurf = 0
          snowmodel_struc(n)%sm(t)%rainf = 0
          snowmodel_struc(n)%sm(t)%snowf = 0
       enddo

!------------------------------------------------------------------------
!      Model timestep Alarm
!------------------------------------------------------------------------
       call LIS_update_timestep(LIS_rc, n, snowmodel_struc(n)%ts)

       write(fnest,'(i3.3)') n
       call LIS_registerAlarm("SnowModel model alarm "//trim(fnest),&
            snowmodel_struc(n)%ts,&
            snowmodel_struc(n)%ts)

       call LIS_registerAlarm("SnowModel restart alarm "//trim(fnest),&
            snowmodel_struc(n)%ts,&
            snowmodel_struc(n)%rstInterval)

       LIS_sfmodel_struc(n)%ts = snowmodel_struc(n)%ts


!------------------------------------------------------------------------
!      Create fields for LSM2SUBLSM exchanges
!------------------------------------------------------------------------
       call ESMF_ArraySpecSet(arrspec1,&
            rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status, &
            'ESMF_ArraySpecSet failed in snowmodel_init')

       sweField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1, &
            name="Total SWE",&
            rc=status)
       call LIS_verify(status,&
            'ESMF_FieldCreate failed for SWE in snowmodel_init')

       snwdField = ESMF_FieldCreate(&
            grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
            arrayspec=arrspec1, &
            name="Total snowdepth",&
            rc=status)
       call LIS_verify(status,&
            'ESMF_FieldCreate failed for snowdepth in snowmodel_init')

       call ESMF_StateAdd(LIS_SUBLSM2LSM_State(n,eks),&
            (/sweField/),rc=status)
       call LIS_verify(status,&
            'ESMF_StateAdd failed for swe in snowmodel_init')
       call ESMF_StateAdd(LIS_SUBLSM2LSM_State(n,eks),&
            (/snwdField/),rc=status)
       call LIS_verify(status,&
            'ESMF_StateAdd failed for snwd in snowmodel_init')

    enddo

  end subroutine snowmodel_init

end module snowmodel_lsmMod

