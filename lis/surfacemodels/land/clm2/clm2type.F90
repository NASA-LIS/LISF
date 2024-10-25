!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module clm2type
!
!BOP
! !MODULE:  clm2type
!
! !DESCRIPTION:
!  The code in this file provides a description of the 
!  data structure containing the CLM 1-d variables. The 
!  variables specified in the data structure include: 
!
!  \begin{description}
!   \item[nstep]
!    time step number
!   \item[kpatch]
!    tile index
!   \item[itypveg]
!    vegetation type
!   \item[itypwat]
!    water type
!   \item[itypprc]
!    precipitation type
!   \item[isoicol]
!    color class for soil albedos
!   \item[snl]
!    number of snow layers
!   \item[frac\_veg\_nosno]
!    fraction of vegetation not covered by snow 
!   \item[frac\_veg\_nosno\_alb]
!    fraction of vegetation not covered by snow
!   \item[imelt]
!    flag for meltiing (=1), freezing (=2), not=0(new)
!   \item[lakpoi]
!    flag for lakpoint (true=lake point)
!   \item[do\_capsnow]
!    flag to indicate snow capping (true=do sno caping)
!   \item[present]
!    whether PFT is present in the current patch 
!   \item[lat]
!    latitude of the patch
!   \item[lon]
!    longitude of the patch
!   \item[dtime]
!    model timestep
!   \item[zi]
!    interface level below a "z" level (m)
!   \item[dz]
!    layer depth (m)
!   \item[z]
!    layer thickness (m)
!   \item[bsw]
!    Clapp and Hornberger "b"
!   \item[watsat]
!    volumetric soil water at saturation (porosity)
!   \item[hksat]
!    hydraulic conductivity at saturation (mm H2O /s)
!   \item[sucsat]
!    minimum soil suction (mm)
!   \item[csol]
!    heat capacity, soil solids (J/m**3/Kelvin)
!   \item[tkmg]
!    thermal conductivity, soil minerals  [W/m-K]  (new)
!   \item[tkdry]
!    thermal conductivity, dry soil       (W/m/Kelvin)
!   \item[tksatu]
!    thermal conductivity, saturated soil [W/m-K]  (new)
!   \item[rootfr]
!    fraction of roots in each soil layer
!   \item[rootr]
!    effective fraction of roots in each layer
!   \item[begwb]
!    water mass at the beginning of the time step
!   \item[endwb]
!    water mass at the end of the time step
!   \item[forc\_t]
!    atmospheric temperature (K)
!   \item[forc\_u]
!    atmospheric wind speed in east direction (m s-1)
!   \item[forc\_v]
!    atmospheric wind speed in north direction (m s-1)
!   \item[forc\_q]
!    atmospheric specific humidity (kg kg-1)
!   \item[forc\_hgt]
!    atmospheric reference height (m) 
!   \item[forc\_hgt\_u]
!    observational height of wind [m] 
!   \item[forc\_hgt\_v]
!    observational height of temperature [m] 
!   \item[forc\_hgt\_q]
!    observational height of humidity [m] (new)
!   \item[forc\_pbot]
!    atmospheric pressure (Pa)
!   \item[forc\_th]
!    atmospheric potential temperature (Kelvin)
!   \item[forc\_vp]
!    atmospheric vapor pressure (Pa)
!   \item[forc\_rho]
!    density (kg/m**3)
!   \item[forc\_lwrad]
!    downward infrared (longwave) radiation (W/m**2)
!   \item[forc\_solad]
!    direct beam radiation (vis=forc\_sols , nir=forc\_soll )
!   \item[forc\_solai]
!    diffuse radiation     (vis=forc\_solsd, nir=forc\_solld)
!   \item[forc\_ch]
!    heat/moisture exchange coefficient
!   \item[forc\_rain]
!    rain rate [mm s-1]
!   \item[forc\_snow]
!    snow rate [mm s-1]
!   \item[rssun]
!    sunlit stomatal resistance (s/m)
!   \item[rssha]
!    shaded stomatal resistance (s/m)
!   \item[psnsun]
!    sunlit leaf photosynthesis (umol CO2 /m**2/ s) 
!   \item[psnsha]
!    shaded leaf photosynthesis (umol CO2 /m**2/ s)
!   \item[laisun]
!    sunlit leaf area
!   \item[laisha]
!    shaded leaf area
!   \item[sabg]
!    solar radiation absorbed by ground (W/m**2)
!   \item[sabv]
!    solar radiation absorbed by vegetation (W/m**2)
!   \item[fsa]
!    solar radiation absorbed (total) (W/m**2)
!   \item[taux]
!    wind stress: e-w (kg/m s-1**2)
!   \item[tauy]
!    wind stress: n-s (kg/m s-1**2)
!   \item[eflx\_lwrad\_out]
!    emitted infrared (longwave) radiation (W/m**2) 
!   \item[eflx\_lwrad\_net]
!    net infrared (longwave) rad (W/m**2) [+ = to atm]
!   \item[eflx\_sh\_tot]
!    total sensible heat flux (W/m**2) [+ to atm]
!   \item[eflx\_sh\_veg]
!    sensible heat flux from leaves (W/m**2) [+ to atm]
!   \item[eflx\_sh\_grnd]
!    sensible heat flux from ground (W/m**2) [+ to atm]
!   \item[eflx\_lh\_tot]
!    total latent heat flux (W/m8*2)  [+ to atm] 
!   \item[eflx\_soil\_grnd]
!    soil heat flux (W/m**2) [+ = into soil]
!   \item[t\_veg]
!    vegetation temperature (Kelvin)
!   \item[t\_grnd]
!    ground temperature (Kelvin)
!   \item[t\_rad]
!    radiative temperature (Kelvin)
!   \item[t\_ref2m]
!    2 m height surface air temperature (Kelvin)
!   \item[t\_soisno]
!    soil temperature (Kelvin)
!   \item[qflx\_infl]
!    infiltration (mm H2O /s) 
!   \item[qflx\_surf]
!    surface runoff (mm H2O /s) 
!   \item[qflx\_drain]
!    sub-surface runoff (mm H2O /s) 
!   \item[qflx\_top\_soil]
!    net water input into soil from top (mm s-1)
!   \item[qflx\_evap\_soi]
!    soil evaporation (mm H2O/s) (+ = to atm)
!   \item[qflx\_evap\_veg]
!    vegetation evaporation (mm H2O/s) (+ = to atm)
!   \item[qflx\_tran\_veg]
!    vegetation transpiration (mm H2O/s) (+ = to atm)
!   \item[qflx\_snomelt]
!    snow melt (mm H2O /s)
!   \item[qflx\_evap\_tot]
!    qflx\_evap\_soi + qflx\_evap\_veg + qflx\_tran\_veg
!   \item[qflx\_rain\_grnd]
!    rain on ground after interception (mm H2O/s) [+]
!   \item[qflx\_evap\_grnd]
!    ground surface evaporation rate (mm H2O/s) [+]
!   \item[qflx\_dew\_grnd]
!    ground surface dew formation (mm H2O /s) [+]
!   \item[qflx\_sub\_snow]
!    sublimation rate from snow pack (mm H2O /s) [+]
!   \item[qflx\_dew\_snow]
!    surface dew added to snow pack (mm H2O /s) [+]
!   \item[qflx\_snowcap]
!    excess precipitation due to snow capping (mm H2O /s) [+]
!   \item[qflx\_qrgwl]
!    qflx\_surf at glaciers, wetlands, lakes
!   \item[h2osno]
!    snow water (mm H2O)
!   \item[h2ocan]
!    canopy water (mm H2O)
!   \item[h2osoi\_liq]
!    liquid water (kg m-2) 
!   \item[h2osoi\_ice]
!    ice lens (kg m-2) 
!   \item[h2osoi\_vol]
!    volumetric soil water (0<=h2osoi\_vol<=watsat) [m^3 m-3]
!   \item[snowdp]
!    snow height (m) 
!   \item[snowage]
!    non dimensional snow age [-] 
!   \item[h2osno\_old]
!    snow mass for previous time step (kg m-2) (new)
!   \item[frac\_sno]
!    fraction of ground covered by snow (0 to 1) 
!   \item[frac\_iceold]
!    fraction of ice relative to the total water (new)
!   \item[eff\_porosity]
!    effective porosity = porosity - vol\_ice
!   \item[parsun]
!    average absorbed PAR for sunlit leaves (W/m**2)
!   \item[albgrd]
!    ground albedo (direct)                      
!   \item[albgri]
!    ground albedo (diffuse)                     
!   \item[fabd]
!    flux absorbed by veg per unit direct flux   
!   \item[fabi]
!    flux absorbed by veg per unit diffuse flux  
!   \item[ftdd]
!    down direct flux below veg per unit dir flx 
!   \item[ftid]
!    down diffuse flux below veg per unit dir flx
!   \item[ftii]
!    down diffuse flux below veg per unit dif flx
!   \item[fsun]
!    sunlit fraction of canopy
!   \item[surfalb]
!    instantaneous all-wave surface albedo
!   \item[snoalb]
!    instantaneous all-wave snow albedo
!   \item[hbot]
!    canopy bottom (m)
!   \item[htop]
!    canopy top (m)
!   \item[tlai]
!    one-sided leaf area index, no burying by snow
!   \item[tsai]
!    one-sided stem area index, no burying by snow
!   \item[elai]
!    one-sided leaf area index with burying by snow
!   \item[esai]
!    one-sided stem area index with burying by snow
!   \item[fwet]
!    fraction of canopy that is wet (0 to 1)
!   \item[fdry]
!    fraction of foliage that is green and dry [-] (new)
!   \item[annpsn]
!    annual photosynthesis (umol CO2 /m**2)
!   \item[annpsnpot]
!    annual potential photosynthesis (same units)
!   \item[wf]
!    soil water as frac. of whc for top 0.5 m
!   \item[z0mr]
!    ratio of momentum roughness length to canopy top height [-]
!   \item[z0m]
!    momentum roughness length [m]
!   \item[displar]
!    ratio of displacement height to canopy top height [-]
!   \item[displa]
!    displacement height [m]
!   \item[dleaf]
!    leaf dimension [m]
!   \item[xl]
!    pft\_varcon leaf/stem orientation index
!   \item[rhol]
!     pft\_varcon leaf reflectance  : 1=vis, 2=nir 
!   \item[rhos]
!    pft\_varcon stem reflectance  : 1=vis, 2=nir 
!   \item[taul]
!    pft\_varcon leaf transmittance: 1=vis, 2=nir 
!   \item[taus]
!    pft\_varcon stem transmittance: 1=vis, 2=nir 
!   \item[qe25]
!    quantum efficiency at 25c (umol co2 / umol photon)
!   \item[vcmx25]
!    maximum rate of carboxylation at 25c (umol co2/m**2/s)
!   \item[mp]
!    slope for conductance-to-photosynthesis relationship
!   \item[c3psn]
!    photosynthetic pathway: 0. = c4, 1. = c3
!   \item[totfsa]
!     solar absorbed solar radiation [W m-2]
!   \item[toteflx\_lwrad\_net]
!     net longwave radiation [W m-2]
!   \item[toteflx\_lh\_tot]
!    total latent heat flux [W m-2]
!   \item[toteflx\_sh\_tot]
!    total sensible heat flux [W m-2]
!   \item[toteflx\_soil\_grnd]
!    ground heat flux [W m-2]
!   \item[toqflx\_snomelt]
!    snowmelt heat flux [W m-2]
!   \item[totrain]
!    accumulation of rain [mm]
!   \item[totsnow]
!    accumulation of snow [mm]
!   \item[totqflx\_evap]
!    total evaporation [mm]
!   \item[totqflx\_surf]
!    surface runoff [mm]
!   \item[totqflx\_drain]
!    subsurface runoff [mm]
!   \item[totqflx\_ecanop]
!     interception evaporation [W m-2]
!   \item[totqflx\_tran\_veg]
!     Total vegetation transpiration
!   \item[totqflx\_evap\_grnd]
!     Total ground surface evaporation 
!   \item[totqflx\_sub\_snow]
!     Total sublimation rate from snow pack 
!   \item[acond]
!     aerodyamic conductance
!   \item[soilmtc\_prev]
!     Total column soil moisture for the prev.timestep
!   \item[h2osno\_prev]
!     Total column snow water equivalent for the prev.timestep
!  \end{description}
!
! !USES: 
  use LIS_precisionMod, only : r8
!EOP
  implicit none

  PRIVATE 

  type, public :: clm1d

     integer  :: nstep          
     integer  :: kpatch         
     integer  :: itypveg        
     integer  :: itypwat        
     integer  :: itypprc        
     integer  :: isoicol        
     integer  :: snl             
     integer  :: frac_veg_nosno 
     integer  :: frac_veg_nosno_alb
     integer  :: imelt(-4:10)      
     logical  :: lakpoi       
     logical  :: do_capsnow   
     logical  :: present   
     real(r8) :: lat       
     real(r8) :: lon       
     real(r8) :: dtime
     real(r8) :: zi(-5:10)
     real(r8) :: dz(-4:10)
     real(r8) :: z (-4:10)
     real(r8) :: bsw   (10)  
     real(r8) :: watsat(10)
     real(r8) :: hksat (10)
     real(r8) :: sucsat(10)      
     real(r8) :: csol  (10)      
     real(r8) :: tkmg  (10)       
     real(r8) :: tkdry (10)       
     real(r8) :: tksatu(10)       
     real(r8) :: rootfr(10)       
     real(r8) :: rootr(10)        
     real(r8) :: begwb
     real(r8) :: endwb
     real(r8) :: forc_t
     real(r8) :: forc_u
     real(r8) :: forc_v
     real(r8) :: forc_q       
     real(r8) :: forc_hgt     
     real(r8) :: forc_hgt_u   
     real(r8) :: forc_hgt_t   
     real(r8) :: forc_hgt_q   
     real(r8) :: forc_pbot    
     real(r8) :: forc_th      
     real(r8) :: forc_rho        
     real(r8) :: forc_lwrad      
     real(r8) :: forc_solad(2)  
     real(r8) :: forc_solai(2)  
     real(r8) :: forc_ch
     real(r8) :: forc_chs2
     real(r8) :: forc_cqs2
     real(r8) :: forc_q2sat
     real(r8) :: forc_rain           !rain rate [mm s-1]
     real(r8) :: forc_snow           !snow rate [mm s-1]
     real(r8) :: forc_cosz

!*************************************************************************
! biogeophys
!*************************************************************************

! Surface solar radiation 

     real(r8) :: rssun
     real(r8) :: rssha          
     real(r8) :: psnsun         
     real(r8) :: psnsha         
     real(r8) :: laisun         
     real(r8) :: laisha         
     real(r8) :: sabg        
     real(r8) :: sabv        
     real(r8) :: fsa           

! Surface energy fluxes

     real(r8) :: taux           
     real(r8) :: tauy           
     real(r8) :: eflx_lwrad_out 
     real(r8) :: eflx_lwrad_net 
     real(r8) :: eflx_sh_tot   
     real(r8) :: eflx_sh_veg    
     real(r8) :: eflx_sh_grnd   
     real(r8) :: eflx_lh_tot   
     real(r8) :: eflx_soil_grnd

! Temperatures

     real(r8) :: t_veg                        
     real(r8) :: t_grnd                       
     real(r8) :: t_rad                        
     real(r8) :: t_ref2m                      
     real(r8) :: t_soisno(-4:10) 


!*************************************************************************
! hydrology
!*************************************************************************

     real(r8) :: qflx_infl                       
     real(r8) :: qflx_surf                       
     real(r8) :: qflx_drain                      
     real(r8) :: qflx_top_soil                   
     real(r8) :: qflx_evap_soi                   
     real(r8) :: qflx_evap_veg                   
     real(r8) :: qflx_tran_veg                   
     real(r8) :: qflx_snomelt                    
     real(r8) :: qflx_evap_tot                   
     real(r8) :: qflx_rain_grnd                  
     real(r8) :: qflx_evap_grnd                  
     real(r8) :: qflx_dew_grnd                   
     real(r8) :: qflx_sub_snow
     real(r8) :: qflx_dew_snow
     real(r8) :: qflx_snowcap 
     real(r8) :: qflx_qrgwl   
     real(r8) :: h2osno                          
     real(r8) :: h2ocan                          
     real(r8) :: h2osoi_liq(-4:10)  
     real(r8) :: h2osoi_ice(-4:10)  
     real(r8) :: h2osoi_vol(10)             
     real(r8) :: snowdp                          
     real(r8) :: snowage          
     real(r8) :: h2osno_old                      
     real(r8) :: frac_sno                        
     real(r8) :: frac_iceold(-4:10) 
     real(r8) :: eff_porosity(10)           



!*************************************************************************
! surfacealbedo (for next time step)
!*************************************************************************

     real(r8) :: parsun         
     real(r8) :: albgrd(2) 
     real(r8) :: albgri(2) 
     real(r8) :: fabd(2)   
     real(r8) :: fabi(2)   
     real(r8) :: ftdd(2)   
     real(r8) :: ftid(2)   
     real(r8) :: ftii(2)   
     real(r8) :: fsun      
     real(r8) :: surfalb   
     real(r8) :: snoalb    

!*************************************************************************
! ecosysdynamics
!*************************************************************************

     real(r8) :: hbot           
     real(r8) :: htop           
     real(r8) :: tlai           
     real(r8) :: tsai           
     real(r8) :: elai           
     real(r8) :: esai                
     real(r8) :: fwet           
     real(r8) :: fdry           
     real(r8) :: annpsn         
     real(r8) :: annpsnpot      
     real(r8) :: wf             

!*************************************************************************
! terms from pft_varcon - to avoid indirect indexing
!*************************************************************************

     real(r8) :: z0mr           
     real(r8) :: z0m            
     real(r8) :: displar        
     real(r8) :: displa         
     real(r8) :: dleaf           
     real(r8) :: xl             
     real(r8) :: rhol(2)   
     real(r8) :: rhos(2)   
     real(r8) :: taul(2)   
     real(r8) :: taus(2)   
     real(r8) :: qe25      
     real(r8) :: vcmx25    
     real(r8) :: mp        
     real(r8) :: c3psn     

! -----------------------------------------------------------------

!     real(r8) ::  totfsa             
!     real(r8) ::  toteflx_lwrad_net  
!     real(r8) ::  toteflx_lh_tot     
!     real(r8) ::  toteflx_sh_tot     
!     real(r8) ::  toteflx_soil_grnd  
!     real(r8) ::  totqflx_snomelt    
!     real(r8) ::  totrain            
!     real(r8) ::  totsnow            
!     real(r8) ::  totqflx_evap       
!     real(r8) ::  totqflx_surf       
!     real(r8) ::  totqflx_drain      
!     real(r8) ::  totqflx_ecanop      
!     real(r8) ::  totqflx_tran_veg
!     real(r8) ::  totqflx_evap_grnd
!     real(r8) ::  totqflx_sub_snow
     real(r8) ::  acond
     real(r8) ::  ch
     real(r8) ::  chs2
     real(r8) ::  cqs2
     real(r8) ::  canopint
     real(r8) ::  soilmtc_prev      
     real(r8) ::  h2osno_prev       

!Variables removed from the module: 
!     real(r8) :: litterag       !above ground litter
!     real(r8) :: litterbg       !above ground litter
!     real(r8) :: cpool_fast
!     real(r8) :: cpool_slow
!     real(r8) :: k_fast_ave
!     real(r8) :: k_slow_ave
!     real(r8) :: litter_decom_ave

!     real(r8) ::  totsolisbd          ! total downward surface shortwave radiation [W m-2]
!     real(r8) ::  totforc_lwrad       ! atmospheric infrared (longwave) radiation [W m-2]

!    real(r8) :: dewmx                 !Maximum allowed dew [mm]
!    real(r8) :: wtfact                !Fraction of model area with high water table
!     real(r8) :: trsmx0                !Max transpiration for moist soil+100% veg. [mm s-1]
!    real(r8) :: zlnd                  !Roughness length for soil [m] (new)             
!    real(r8) :: zsno                  !Roughness length for snow [m] (new)             
!    real(r8) :: csoilc                !Drag coefficient for soil under canopy [-] (new)
!    real(r8) :: cnfac                 !Crank Nicholson factor (between 0 and 1) (new)
!    real(r8) :: capr                  !Tuning factor to turn first layer T into surface T (new)  
!    real(r8) :: ssi                   !Irreducible water saturation of snow (new)
!    real(r8) :: wimp                  !Water impremeable if porosity less than wimp (new)
!    real(r8) :: pondmx                !Ponding depth (mm) (new)
!    real(r8) :: smpmax                !wilting point potential in mm (new)
!    real(r8) :: smpmin                !restriction for min of soil potential (mm) (new)
!     real(r8) :: errh2o                !water conservation error (mm H2O)
!     real(r8) :: errsoi                !soil/lake energy conservation error (W/m**2)
!     real(r8) :: errseb                !surface energy conservation error (W/m**2)
!     real(r8) :: errsol                !solar radiation conservation error (W/m**2)
!     real(r8) :: errlon                !longwave radiation conservation error (W/m**2)
!     real(r8) :: acc_errh2o            !accumulation of water balance error
!     real(r8) :: acc_errseb            !accumulation of surface energy balance error
!     real(r8) :: acc_errsoi            !accumulation of energy balance error
!     real(r8) :: acc_errsol            !accumulation of energy balance error
!    real(r8) :: forc_co2            !atmospheric CO2 concentration (Pa)
!    real(r8) :: forc_o2             !atmospheric O2 concentration (Pa)
!     real(r8) :: forc_swc1        ! Layer 1 (0-10 cm) soil water content (% capacity for GEOS)
!     real(r8) :: forc_sdepth      ! Model liquid equivalent snow depth (kg m-2) (same for GEOS)
!     real(r8) :: ndvi           !Normalized Difference Vegetation Index
!     real(r8) :: fsr            !solar radiation reflected (W/m**2)
     real(r8) :: eflx_lh_vege   !veg evaporation heat flux (W/m**2) [+ to atm]
     real(r8) :: eflx_lh_vegt   !veg transpiration heat flux (W/m**2) [+ to atm]
     real(r8) :: eflx_lh_grnd   !ground evaporation heat flux (W/m**2) [+ to atm]   
     real(r8) :: eflx_snomelt   !snow melt heat flux (W/m**2)
!     real(r8) :: eflx_impsoil   !implicit evaporation for soil temperature equation (W/m**2)
     real(r8) :: t_lake(1:10)            !lake temperature (Kelvin)
     real(r8) :: t_snow                       !vertically averaged snow temperature
!     real(r8) :: dt_veg                       !change in t_veg, last iteration (Kelvin)
!     real(r8) :: dt_grnd                      !change in t_grnd, last iteration (Kelvin)
     real(r8) :: btran          !transpiration wetness factor (0 to 1) 

!*************************************************************************
! biogeochem
!*************************************************************************

!     real(r8) :: fpsn           !photosynthesis (umol CO2 /m**2 /s)
!     real(r8) :: frm            !total maintenance respiration (umol CO2 /m**2/s)
!     real(r8) :: frmf           !leaf maintenance respiration  (umol CO2 /m**2 /s)
!     real(r8) :: frms           !stem maintenance respiration  (umol CO2 /m**2 /s)
!     real(r8) :: frmr           !root maintenance respiration  (umol CO2 /m**2 /s)
!     real(r8) :: frg            !growth respiration (umol CO2 /m**2 /s)
!     real(r8) :: fmicr          !microbial respiration (umol CO2 /m**2 /s)
!     real(r8) :: fco2           !net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]
!     real(r8) :: dmi            !total dry matter production (ug /m**2 /s)

! dust model

!     real(r8) :: flx_mss_vrt_dst(4) !surface dust emission (kg/m**2/s) [ + = to atm]
!     real(r8) :: vwc_thr               !threshold soil moisture based on clay content
!     real(r8) :: mss_frc_cly_vld       ![frc] Mass fraction clay limited to 0.20
!     real(r8) :: mbl_bsn_fct

! voc model

!     real(r8) :: vocflx(4)          !VOC flux [ug C m-2 h-1]
     real(r8) :: qflx_prec_intr                  !interception of precipitation [mm s-1]
     real(r8) :: qflx_prec_grnd                  !water onto ground including canopy runoff [kg/(m2 s)]
!     real(r8) :: qflx_snow_grnd                  !snow on ground after interception (mm H2O/s) [+]
! velocities

!     real(r8) :: u10                       !10-m wind (m s-1)
!     real(r8) :: fv                        !friction velocity (m s-1)
     real(r8) :: snowice                         !average snow ice lens
     real(r8) :: snowliq                         !average snow liquid water
!     real(r8) :: rsw                             !soil water content for root zone
!     real(r8) :: sfact                           !term for implicit correction to evaporation
!     real(r8) :: sfactmax                        !maximim of "sfact"
!     real(r8) :: parsha         !average absorbed PAR for shaded leaves (W/m**2)
!     real(r8) :: albd(2)   !surface albedo (direct)                     
!     real(r8) :: albi(2)   !surface albedo (diffuse)                    
!     real(r8) :: top_t1_f
!     real(r8) :: top_t2_f
!     real(r8) :: bot_t1_f
!     real(r8) :: bot_t2_f
!     real(r8) :: surfalb_old    !previous instantaneous all-wave surface albedo
!     real(r8) :: snoalb_old     !previous instantaneous all-wave snow albedo

!*************************************************************************
! terms due to splitting the code into Biogeophys1 and Biogeophys2
!*************************************************************************

!     real(r8) cgrnd  ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
!     real(r8) cgrndl ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
!     real(r8) cgrnds ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
!     real(r8) tg     ! ground surface temperature [K]
!     real(r8) tssbef(-5:10)  ! soil/snow temperature before update
!     real(r8) qg     ! ground specific humidity [kg kg-1]
!     real(r8) dqgdT  ! d(qg)/dT
!     real(r8) emg    ! ground emissivity
!     real(r8) emv    ! vegetation emissivity
!     real(r8) htvp   ! latent heat of vapor of water (or sublimation) [j/kg]
!     real(r8) z0mg   ! roughness length over ground, momentum [m]
!     real(r8) z0hg   ! roughness length over ground, sensible heat [m]
!     real(r8) z0qg   ! roughness length over ground, latent heat [m]
!     real(r8) z0mv   ! roughness length over vegetation, momentum [m]
!     real(r8) z0hv   ! roughness length over vegetation, sensible heat [m]
!     real(r8) z0qv   ! roughness length over vegetation, latent heat [m]
!     real(r8) beta   ! coefficient of convective velocity [-]
!     real(r8) zii    ! convective boundary height [m]
!     real(r8) thm    ! intermediate variable (forc_t+0.0098*forc_hgt_t)
!     real(r8) thv    ! virtual potential temperature (kelvin)
!     real(r8) ur     ! wind speed at reference height [m s-1] ***DO WE NEED THIS???
!     real(r8) dlrad  ! downward longwave radiation below the canopy [W m-2]
!     real(r8) ulrad  ! upward longwave radiation above the canopy [W m-2]
     real(r8) qmelt  ! snow melt [mm s-1]
!     real(r8) :: agdd0          !accumulated growing degree days above 0 deg C
!     real(r8) :: agdd5          !accumulated growing degree days above -5
!     real(r8) :: agddtw         !accumulated growing degree days above twmax
!     real(r8) :: agdd           !accumulated growing degree days above 5
!     real(r8) :: t10            !10-day running mean of the 2 m temperature (K)
!     real(r8) :: t_mo           !30-day average temperature (Kelvin)
!     real(r8) :: t_mo_min       !annual min of t_mo (Kelvin)
!     real(r8) :: fnpsn10        !10-day running mean net photosynthesis
!     real(r8) :: prec365        !365-day running mean of tot. precipitation
!     real(r8) :: agdd20         !20-yr running mean of agdd
!     real(r8) :: tmomin20       !20-yr running mean of tmomin
!     real(r8) :: t10min         !annual minimum of 10-day running mean (K)
!     real(r8) :: tsoi25         !soil temperature to 0.25 m (Kelvin)
!     real(r8) :: dphen          !phenology [0 to 1]
!     real(r8) :: leafon         !leafon days
!     real(r8) :: leafof         !leafoff days
!     real(r8) :: nind           !number of individuals
!     real(r8) :: lm_ind         !individual leaf mass
!     real(r8) :: sm_ind         !individual stem mass
!     real(r8) :: hm_ind         !individual heartwood mass
!     real(r8) :: rm_ind         !individual root mass
!     real(r8) :: lai_ind        !LAI per individual
!     real(r8) :: fpcinc         !
!     real(r8) :: fpcgrid        !
!     real(r8) :: crownarea      !
!     real(r8) :: bm_inc         !
!     real(r8) :: afmicr         !
!     real(r8) :: firelength     !fire season in days
  end type clm1d

  SAVE

end module clm2type



