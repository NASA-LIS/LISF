!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module jules51_module
!BOP
!
! !MODULE: jules51_module.F90
!
! !DESCRIPTION:
!  Declare forcing-only option (jules51) variables
!
!  \begin{description}
!   \item[forcing]
!     Array of meteorological forcing
!   \item{tair}
!     2m air temperature forcing
!   \item{qair}
!     2m specific humidity forcing
!   \item{swdown}
!     downward shortwave forcing
!   \item{lwdown}
!     downward longwave forcing
!   \item{uwind}
!     u-wind component forcing
!   \item{vwind}
!     v-wind component forcing
!   \item{psurf}
!     surface pressure forcing
!   \item{rainf}
!     total rainfall forcing
!   \item{rainf\_c}
!     convective rainfall forcing
!   \item{snowf}
!     total snowfall forcing
!   \end{description}
!
!
! !REVISION HISTORY:
! 16 May 2016; Shugong Wang; initial implementation for JULES 4.3
! 01 Feb 2018; Shugong Wang; updated for JULES.5.1 
! 
!EOP
  implicit none

  type jules51dec
     integer :: pft 
! state variable 
!     real, allocatable :: canopy(:)          ! Surface/canopy water for snow-free land tiles (kg m-2)  
!     real, allocatable :: cs(:)              ! Soil carbon (kg C/m2)
!     real :: gs                 ! "Stomatal" conductance to evaporation (m s-1)
!     real, allocatable :: snow_tile(:)       ! Lying snow on tiles (kg m-2)
!     real, allocatable :: sthuf(:)           ! sthu is used as storage for total wetness 
!     real, allocatable :: t_soil(:)          ! Sub-surface temperatures (K) 
!     real, allocatable :: tstar_tile(:)      ! Tile surface temperatures (K) 
!     real, allocatable :: rho_snow(:)        ! Snowpack bulk density (kg/m3)
!     real, allocatable :: snow_depth(:)      ! Snow depth on ground on tiles (m) 
     real, allocatable :: frac(:)            ! fractional cover of each surface type
     real :: frac_agr           ! fraction of agriculture 
!     real, allocatable :: co2_mmr(:)         ! CO2 Mass Mixing Ratio (parameter?)







!-------------------------------------------------------------------------
! spatial parameters
!-------------------------------------------------------------------------
     real :: b
     real :: sathh
     real :: satcon
     real :: sm_sat
     real :: sm_crit
     real :: sm_wilt
     real :: hcap
     real :: hcon
     real :: albsoil
!-------------------------------------------------------------------------
! Forcing Variables
!-------------------------------------------------------------------------
     real :: tair
     real :: qair
     real :: swdown
     real :: lwdown
     real :: wind_e
     real :: wind_n
     real :: psurf
     real :: rainf
     real :: rainf_c
     real :: snowf

!-------------------------------------------------------------------------
! State Variables
!-------------------------------------------------------------------------

! prognostics 
    integer, allocatable :: nsnow(:)       ! ntiles  number of snow layers on ground on tiles
    real, allocatable :: tsoil_deep(:)     ! ns_deep  deep soil temperatures (k)
    real, allocatable :: sice(:,:)         ! ntiles, nsmax, Snow layer ice mass on tiles (Kg/m2)
    real, allocatable :: sliq(:,:)         ! ntiles, nsmax, Snow layer liquid mass on tiles (Kg/m2)
    real, allocatable :: snowdepth(:)      ! ntiles, Snow depth on ground on tiles (m)
    real, allocatable :: tsnow(:,:)        ! ntiles, nsmax,Snow layer temperature (K)
    real, allocatable :: rgrainl(:,:)      ! ntiles, nsmax, Snow layer grain size on tiles (microns)
    real, allocatable :: rho_snow_grnd(:)  ! ntiles, Snowpack bulk density (kg/m3)
    real, allocatable :: rho_snow(:,:)     ! ntiles, nsmax, Snow layer densities (m)
    real, allocatable :: ds(:,:)           ! ntiles, nsmax, Snow layer thickness (m)
    real              :: wood_prod_fast    ! Fast-turnover wood product C pool.
    real              :: wood_prod_med     ! Medium-turnover wood product C pool. 
    real              :: wood_prod_slow    ! Slow-turnover wood product C pool.
    real              :: frac_agr_prev     ! Agricultural fraction from previous TRIFFID call
    real              :: frac_past_prev    ! Pasture fraction from previous TRIFFID call
    real, allocatable :: n_inorg_soilt_lyrs(:) ! Gridbox Inorganic N pool on soil levels (kg N/m2)
    real, allocatable :: n_inorg_avail_pft(:,:) ! Availabile inorganic N for PFTs (depends on roots)
    real              :: n_inorg           ! Gridbox Inorganic N pool (kg N m-2)
    real, allocatable :: ns(:,:)           ! dim_cslayer, dim_cs1, Soil Organic Nitrogen (kg N m-2)
    real              :: triffid_co2_gb    !  Atmospheric CO2 fluxes from TRIFFID (kgC/m2/yr)
    real, allocatable :: canht_ft(:)       ! npft, Canopy height (m)
    real, allocatable :: canopy(:)         ! ntiles, Surface/canopy water for snow-free land tiles (kg m-2)
    real              :: canopy_gb         ! Gridbox canopy water content (kg m-2)
    real, allocatable :: cs_pool_soilt(:,:)             !cs_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)  dim_cs1, Soil carbon (kg C/m2)
    real              :: di                !  "Equivalent thickness" of sea-ice GBM agregate (m)
    real, allocatable :: di_ncat(:)        ! nice, "Equivalent thickness" of sea-ice catagories (m)
    real, allocatable :: k_sice(:)         ! nice, Sea ice effective conductivity (2*kappai/de)
    real, allocatable :: gc_surft(:)             ! ntiles, Stomatal" conductance to evaporation for land tiles(m s-1)
    real              :: gs_gb                !  "Stomatal" conductance to evaporation (m s-1)
    real, allocatable :: lai(:)            ! npft LAI of plant functional types
    real, allocatable :: rgrain(:)         ! ntiles, Snow surface grain size on tiles (microns)
    real              :: smc_soilt               !  Soil moisture in a layer at the surface (kg m-2).
    real, allocatable :: smcl_soilt(:)           ! sm_levels, Soil moisture content of layers (kg m-2)
    real, allocatable :: snow_tile(:)      ! ntiles, Lying snow on tiles (kg m-2)
    real, allocatable :: snow_grnd(:)      ! ntiles, Snow on the ground (kg m-2)
    real              :: snow_mass_ij        !  Gridbox snowmass (kg m-2)
    real, allocatable :: snow_mass_sea_sicat(:) !nice_use,  Snow on category sea-ice (Kg/m2)
    real              :: soot_ij              !  Snow soot content (kg kg-1)
    real, allocatable :: t_soil(:)         ! sm_levels, Sub-surface temperatures (K)
    real, allocatable :: ti_sicat(:)                !  Sea-ice surface layer
    real, allocatable :: tstar_tile(:)     ! ntiles, Tile surface temperatures (K)
    real, allocatable :: tsurf_elev_surft(:)  ! Tiled land-ice bedrock subsurface temperatures (K)
    real              :: z0msea            !  Sea-surface roughness length for momentum (m).
!    real, allocatable :: routestore(:,:)                ! channel storage (kg). This is defined at all
!    integer, allocatable ::seed_rain(:)                 ! Seeding number for subdaily rainfall for use                ! in IMOGEN or when l_daily_disagg = T

!-------------------------------------------------------------------------
! Plant and Soil Variables and Parameters
!-------------------------------------------------------------------------

    real, allocatable :: p_s_b(:)                  !  Exponent for soil moisture characteristic function Clapp-Hornberger model: b is the Clapp-Hornberger exponent,  van Genuchten model: b=1/(n-1)  (metres)
    real, allocatable :: p_s_sathh(:)              !  Parameter for soil moisture characteristic functions Clapp-Hornberger model: sathh is the saturated soil water pressure (m), van Genuchten model: sathh=1/alpha
    real, allocatable :: p_s_hcap(:)               !  Soil heat capacity (J/K/m3)
    real, allocatable :: p_s_hcon(:)               !  Soil thermal conductivity (W/m/K)
    real, allocatable :: p_s_satcon(:)             !  Saturated hydraulic conductivity (kg m-2/s)
    real, allocatable :: p_s_smvccl(:)             !  Critical volumetric SMC (cubic m per cubic m of soil)    
    real, allocatable :: p_s_smvcst(:)             !  Volumetric saturation point (m^3 m-3 of soil)
    real, allocatable :: p_s_smvcwt(:)             !  Volumetric wilting point (cubic m per cubic m of soil)
    real              :: p_s_albsoil               !  Soil albedo
    real              :: p_s_albobs_sw             !  Obs SW albedo    
    real              :: p_s_albobs_vis            !  Obs VIS albedo    
    real              :: p_s_albobs_nir            !  Obs NIR albedo
    real, allocatable :: p_s_catch(:)              !  Surface/canopy water capacity of snow-free land tiles (kg m-2)
    real, allocatable :: p_s_catch_snow(:)         !  Snow interception capacity (kg m-2)
    real              :: p_s_cosz                  !  Cosine of the zenith angle    
    real, allocatable :: p_s_infil_tile(:)         !  Maximum possible surface infiltration for tiles (kg m-2/s)
    real, allocatable :: p_s_z0_tile(:)            !  Surface roughness on tiles (m).
    real, allocatable :: p_s_z0h_tile_bare(:)      !  Surface thermal roughness on tiles before allowance for snow cover (m).
    real, allocatable :: p_s_sthu(:)               !  Unfrozen soil moisture content of the layers as a fraction of saturation.
    real, allocatable :: p_s_sthf(:)               !  Frozen soil moisture content of the layers as a fraction of saturation.
    real, allocatable :: p_s_clay_frac_soilt(:)    !  Soil clay fraction
    real, allocatable :: p_s_sthu_min(:)           !  Minimum unfrozen water content for each layer. Used to normalise thaw depth calculation based on unfrozen water content fraction.
    real, allocatable :: p_s_v_close_pft(:,:)      !  Volumetric soil moisture concentration below which stomata are fully closed. (m3 H2O/m3 soil).  If l_use_pft_psi=F, this variable is not used, and the soil ancil variable sm_wilt is used instead. 
    real, allocatable :: p_s_v_open_pft(:,:)       !  concentration above which stomatal aperture is not limited by soil water (m3 H2O/m3 soil). If l_use_pft_psi=F, this variable is not used, and the soil ancillary variable sm_crit and the pft variable fsmc_p0 are used instead
    real              :: p_s_clay_soilt            !  Fraction of clay
    real, allocatable :: p_s_soil_ph_soilt(:)      !  Soil pH, defined on soil layers.
    real              :: p_s_z0m_soil_gb               !  Bare soil roughness, for momentum (m).
!-------------------------------------------------------------------------
!   Ancilary data/parameters 
!-------------------------------------------------------------------------
    integer :: ssi_pts                !  Number of sea or sea-ice points
    integer :: sea_pts                !  Number of sea points
    integer :: sice_pts               !  Number of sea-ice points

    integer :: ssi_index           ! index of sea and sea-ice points
    integer :: sea_index           ! index of sea points
    integer :: sice_index          ! index of sea-ice points
    integer, allocatable :: sice_pts_ncat(:)       ! number of points for each sea-ice category
    integer, allocatable :: sice_index_ncat(:)     ! index of points for each sea-ice category
    integer :: soilt_pts           ! !  Number of points for each soil tile
    logical :: l_soil_point        !  TRUE if a soil point  FALSE otherwise
    logical :: l_lice_point        !  TRUE if a land ice point  FALSE otherwise
    logical, allocatable :: l_lice_surft(:)        !   TRUE if a land ice (surface) tile, FALSE otherwise
    real ::  fssi              !  Fraction of gridbox covered by sea or sea-ice
    real ::  sea_frac            !  Fraction of gridbox covered by sea (converted to single vector array)
    real ::  sice_frac           !  Fraction of gridbox covered by sea-ice converted to single vector array)
    real, allocatable ::  sice_frac_ncat(:)    !  Fraction of gridbox covered by each sea-ice category (converted to single vector array)
    real :: frac_soilt           ! Fraction of gridbox for each soil tile 
    integer ::  halo_i                 !  Size of halo in i direction
    integer ::  halo_j                 !  Size of halo in j direction
    integer ::  n_rows                 !  Number of rows in a v field
    integer ::  off_x                  !  Size of small halo in i
    integer ::  off_y                  !  Size of small halo in j
    integer ::  row_length             !  Number of points on a row
    integer ::  rows                   !  Number of rows in a theta field
    
    integer ::  co2_dim_len            !  Length of a CO2 field row
    integer ::  co2_dim_row            !  Number of CO2 field rows
    integer ::  land_pts               !  No. of land points
    integer ::  land_pts_trif          !  For dimensioning land fields in TRIFFID
    integer ::  lice_pts               !  Number of land ice points
    integer ::  npft_trif              !  For dimensioning pft fields in TRIFFID =npft when TRIFFID on, otherwise =1
    integer ::  nsurft                 !  Number of surface tiles
    integer ::  soil_pts               !  Number of soil points
    integer ::  dim_cs1                !  size of second dimension in soil carbon (cs) and related respiration variables
    integer ::  dim_cs2                !  size used for some variables that are only used with TRIFFID. If not using TRIFFID these variables are set to be smaller to save some space
    integer ::  dim_cslayer            ! size of depth dimension in soil carbon (cs)  and related respiration variables
    integer :: land_index          !  index of land points
    integer :: soilt_index          !  index of land points
    integer, allocatable :: tile_index(:)        !  indices of land points which include the nth surface type
    integer :: soil_index          !  index of soil points (i.e. land point number for each soil point)
    integer :: lice_index          !  index of land ice points (i.e. land point number for each land ice point)
    integer, allocatable :: tile_pts(:)            !  Number of land points which include the nth surface type

    ! real, allocatable :: frac(:)              !  fractional cover of each surface type
    real ::  z1_tq             !  height of temperature data
    real ::  z1_uv             !  height of wind data
    real ::  ice_fract         !  fraction of gridbox covered by sea-ice  (decimal fraction)
    real, allocatable ::  ice_fract_ncat(:)  !  fraction of gridbox covered by sea-ice on catagories
    real, allocatable ::  ti_cat(:)          ! sea ice surface temperature on categories
    real, allocatable ::  pond_frac_cat(:)   ! Meltpond fraction on sea ice categories 
    real, allocatable ::  pond_depth_cat(:)  ! Meltpond depth on sea ice categories (m)  
    real ::  sstfrz            ! Salinity-dependent sea surface freezing temperature (K)
 
    logical  :: land_mask       !  t if land, f elsewhere

    !!! BVOC diagnostics
    real :: isoprene         ! Gridbox mean isoprene emission flux (kgC/m2/s) 
    real, allocatable :: isoprene_ft(:)   ! Isoprene emission flux on PFTs (kgC/m2/s) 
    real :: terpene          ! Gridbox mean (mono-)terpene emission flux (kgC/m2/s) 
    real, allocatable :: terpene_ft(:)  ! (Mono-)Terpene emission flux on PFTs (kgC/m2/s) 
    real :: methanol      ! Gridbox mean methanol emission flux (kgC/m2/s) 
    real, allocatable :: methanol_ft(:) ! Methanol emission flux on PFTs (kgC/m2/s) 
    real :: acetone       ! Gridbox mean acetone emission flux (kgC/m2/s) 
    real, allocatable :: acetone_ft(:)  ! Acetone emission flux on PFTs (kgC/m2/s) 
    !!! the variables for TRIFFID and plant phenology
    integer :: asteps_since_triffid       ! Number of atmospheric timesteps since last call to TRIFFID
    real, dimension(:), allocatable :: g_leaf_acc ! Accumulated leaf turnover rate
    real, dimension(:), allocatable :: npp_ft_acc ! Accumulated NPP_FT
    real, dimension(:), allocatable :: g_leaf_phen_acc ! Accumulated leaf turnover rate including phenology
    real, dimension(:), allocatable :: resp_w_ft_acc   ! Accum RESP_W_FT
    real, dimension(:,:), allocatable :: resp_s_acc_soilt     ! Accumulated RESP_S
    real :: gpp ! Gross primary productivity (kg C/m2/s)
    real :: npp ! Net primary productivity (kg C/m2/s)
    real  :: resp_p ! Plant respiration (kg C/m2/s)
    real, dimension(:), allocatable :: g_leaf ! Leaf turnover rate (/360days)
    real, dimension(:), allocatable :: g_leaf_phen ! Mean leaf turnover rate over phenology period(/360days)
    real, dimension(:), allocatable :: gpp_ft ! Gross primary productivity on PFTs (kg C/m2/s)
    real, dimension(:), allocatable :: npp_ft ! Net primary productivity on PFTs (kg C/m2/s)
    real, dimension(:), allocatable :: resp_p_ft ! Plant respiration on PFTs (kg C/m2/s)
    real, dimension(:,:), allocatable :: resp_s_soilt ! Soil respiration (kg C/m2/s)
    real, dimension(:), allocatable :: resp_w_ft ! Wood maintenance respiration (kg C/m2/s)
    real, dimension(:), allocatable :: lai_phen ! LAI of PFTs after phenology. Required as separate variable for top-level argument list matching with VEG_IC2A
    real, dimension(:), allocatable :: c_veg ! Total carbon content of the vegetation (kg C/m2)
    real  :: cv ! Gridbox mean vegetation carbon (kg C/m2)
    real, dimension(:), allocatable :: g_leaf_day ! Mean leaf turnover rate for input to PHENOL (/360days)
    real, dimension(:), allocatable :: g_leaf_dr_out ! Mean leaf turnover rate for driving TRIFFID (/360days)
    real, dimension(:), allocatable :: lit_c ! Carbon Litter (kg C/m2/360days)
    real  :: lit_c_mn ! Gridbox mean carbon litter (kg C/m2/360days)
    real, dimension(:), allocatable :: npp_dr_out ! Mean NPP for driving TRIFFID (kg C/m2/360days)
    real, dimension(:), allocatable :: resp_w_dr_out ! Mean wood respiration for driving TRIFFID (kg C/m2/360days)
    real, dimension(:,:), allocatable :: resp_s_dr_out_gb ! Mean soil respiration for driving TRIFFID (kg C/m2/360days)
    !real  :: frac_agr ! Fraction of agriculture
    !!! variables for Topmodel and PDM
    real  :: fexp ! Decay factor in Sat. Conductivity in water table layer
    real  :: gamtot ! integrated complete gamma function dbc gamtot doesn't need to be in a module in this version, but left there for now for compatability.
    real  :: ti_mean ! Mean topographic index
    real  :: ti_sig ! Standard dev. of topographic index
    real  :: fsat ! Surface saturation fraction
    real  :: fwetl ! Wetland fraction
    real  :: zw ! Water table depth (m)
    real  :: drain ! Drainage out of bottom (nshyd) soil layer (kg m-2/s)
    real  :: dun_roff ! Dunne part of sfc runoff (kg m-2/s)
    real  :: qbase ! Base flow (kg m-2/s)
    real  :: qbase_zw ! Base flow from ZW layer (kg m-2/s)
    real  :: fch4_wetl ! Scaled wetland methane flux (10^-9 kg C/m2/s)
    real  :: fch4_wetl_cs_soilt  ! Scaled wetland methane flux using soil carbon as substrate (kg C/m2/s).
    real  :: fch4_wetl_npp_soilt ! Scaled wetland methane flux using NPP as substrate (kg C/m2/s).
    real  :: fch4_wetl_resps_soilt !Scaled wetland methane flux using soil respiration as substrate (kg C/m2/s).  
    real :: inlandout_atm  ! TRIP inland basin outflow (for land points only)(kg m-2/s)
    real :: sthzw          ! soil moist fraction in deep (water table) layer.
    real :: a_fsat         ! Fitting parameter for Fsat in LSH model
    real :: c_fsat         ! Fitting parameter for Fsat in LSH model
    real :: a_fwet         ! Fitting parameter for Fwet in LSH model
    real :: c_fwet         ! Fitting parameter for Fwet in LSH model
    
    !!! JULES internal
    real, allocatable :: unload_backgrnd(:)

    !!! fluxes 
    real, allocatable :: anthrop_heat_surft(:) ! Additional heat source on surface tiles used for anthropgenic urban heat source (W/m2)
    real, allocatable :: surf_ht_store_surft(:)! Diagnostic to store values of C*(dT/dt) during calculation of energy balance 
    real :: sub_surf_roff                    !  Sub-surface runoff (kg m-2/s)
    real :: surf_roff                        !  Surface runoff (kg m-2/s) 
    real, allocatable :: alb_tile(:,:)       !  Albedo for surface tiles, 1: direct beam visible, 2: diffuse visible, 3: direct beam near-IR, 4: diffuse near-IR   
    real :: tstar                            !  GBM surface temperature (K)
    real :: e_sea                            !  Evaporation from sea times leads fraction. Zero over land (kg m-2/s)
    real :: fqw_1                            !  Moisture flux between layers, FQW(,1) is total water flux from surface, 'E' (kg m-2/s)
    real, allocatable :: fsmc(:)             !  Moisture availability factor
    real :: ftl_1                            !  FTL(,K) contains net turbulent sensible heat flux into layer K from below, so FTL(,1) is the surface sensible heat, H.(W m-2)
    real, allocatable :: ftl_tile(:)         !  Surface FTL for land tiles (W m-2)
    real, allocatable :: le_tile(:)          !  Surface latent heat flux for land tiles (W m-2) 
    real :: h_sea                            !  Surface sensible heat flux over sea times leads fraction (W m-2)
    real :: taux_1                           !  W'ly component of surface wind stress (N/sq m)
    real :: tauy_1                           !  S'ly component of surface wind stress (N/sq m)
    real, allocatable :: fqw_tile(:)         !  Surface FQW (Moisture flux between layers) for land tiles, (kg m-2/s)
    real, allocatable :: fqw_ice(:)          !  Surface FQW (Moisture flux between layers) for sea-ice, (kg m-2/s)
    real, allocatable :: ftl_ice(:)          !  Surface FTL for sea-ice (W m-2) 
    real :: ecan                             !  Gridbox mean evaporation from canopy/surface store (kg m-2/s)
    real, allocatable :: esoil_tile(:)       !  ESOIL for snow-free land tiles
    real, allocatable :: sea_ice_htf(:)      !  Heat flux through sea-ice (W m-2, positive downwards)
    real :: surf_ht_flux                     !  Net downward heat flux at surface over land and sea-ice fraction of gridbox (W m-2)
    real, allocatable :: surf_htf_tile(:)    !  Surface heat flux on land tiles (W m-2)
    real, allocatable :: land_albedo(:)      !  GBM albedo, 1-direct beam visible, 2-diffuse visible, 3-direct beam near-IR, 4-diffuse near-IR
    real :: ei                               !  Sublimation from lying snow or sea-ice (kg m-2/s)
    real, allocatable :: ei_tile(:)          !  EI for land tiles (kg m-2/s)
    real, allocatable :: ecan_tile(:)        !  Canopy evaporation from for snow-free land tiles (kg m-2/s)
    real :: esoil                            !  Surface evapotranspiration from soil moisture store (kg m-2/s)
    real, allocatable :: ext(:)              !  Extraction of water from each soil layer (kg m-2/s)
    real :: snow_melt_gb                     !  Snowmelt on land points (kg m-2/s)
    real :: snow_melt_ij                     !  Snowmelt (kg/m2/s)
    real, allocatable :: melt_tile(:)        !  Snowmelt on land tiles (kg m-2/s)
    real :: hf_snow_melt                     !  Gridbox snowmelt heat flux (W m-2)
    real, allocatable :: radnet_tile(:)      !  Surface net radiation on tiles ( W m-2)
    real, allocatable :: sw_tile(:)          !  Surface net shortwave on tiles (W m-2)
    real, allocatable :: emis_tile(:)        !  Tile emissivity
    real              :: rflow               !  River outflow on model grid (kg m-2/s)
    real              :: rrun                !  Runoff after river routing on model grid (kg m-2/s)
    real              :: snomlt_sub_htf      !  Sub-canopy snowmelt heat flux (W m-2)
    real              :: tot_tfall           !  Total throughfall (kg m-2/s)
    real, allocatable :: snow_soil_htf(:)    !  Heat flux under snow to subsurface on tiles (W/m2)
    real              :: u10m
    real              :: v10m
    real, allocatable :: surft_frac(:)
    integer           :: end_k = 0 

    ! EMK...Special output variables for Air Force
    real :: tair_agl_min
    real :: rhmin
end type jules51dec

end module jules51_module
