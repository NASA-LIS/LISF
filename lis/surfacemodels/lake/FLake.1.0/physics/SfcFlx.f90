! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE SfcFlx

!------------------------------------------------------------------------------
!
! Description:
!
!  The main program unit of 
!  the atmospheric surface-layer parameterization scheme "SfcFlx".
!  "SfcFlx" is used to compute fluxes 
!  of momentum and of sensible and latent heat over lakes.
!  The surface-layer scheme developed by Mironov (1991) was used as the starting point.
!  It was modified and further developed to incorporate new results as to 
!  the roughness lenghts for scalar quantities,
!  heat and mass transfer in free convection,
!  and the effect of limited fetch on the momentum transfer.
!  Apart from the momentum flux and sensible and latent heat fluxes,
!  the long-wave radiation flux from the water surface and
!  the long-wave radiation flux from the atmosphere can also be computed.
!  The atmospheric long-wave radiation flux is computed with simple empirical formulae,
!  where the atmospheric emissivity is taken to be dependent on 
!  the water vapour pressure and cloud fraction.
!
!  A description of SfcFlx is available from the author.
!  Dmitrii Mironov 
!  German Weather Service, Kaiserleistr. 29/35, D-63067 Offenbach am Main, Germany. 
!  dmitrii.mironov@dwd.de 
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!  Lines embraced/marked with "!_dbg" are used
!  for debugging purposes only.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE data_parameters  , ONLY :   &
  ireals                      , & ! KIND-type parameter for real variables
  iintegers                       ! KIND-type parameter for "normal" integer variables

USE flake_parameters , ONLY :   &
  tpl_grav                    , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_f                     , & ! Fresh water freezing point [K]
  tpl_rho_w_r                 , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_c_w                     , & ! Specific heat of water [J kg^{-1} K^{-1}]
  tpl_L_f                     , & ! Latent heat of fusion [J kg^{-1}]
  h_Ice_min_flk                   ! Minimum ice thickness [m]

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Dimensionless constants in the Monin-Obukhov surface-layer 
!  similarity relations and in the expressions for the roughness lengths.
REAL (KIND = ireals), PARAMETER ::   &
  c_Karman      = 0.40_ireals      , & ! The von Karman constant 
  Pr_neutral    = 1.0_ireals       , & ! Turbulent Prandtl number at neutral static stability
  Sc_neutral    = 1.0_ireals       , & ! Turbulent Schmidt number at neutral static stability
  c_MO_u_stab   = 5.0_ireals       , & ! Constant of the MO theory (wind, stable stratification)
  c_MO_t_stab   = 5.0_ireals       , & ! Constant of the MO theory (temperature, stable stratification)
  c_MO_q_stab   = 5.0_ireals       , & ! Constant of the MO theory (humidity, stable stratification)
  c_MO_u_conv   = 15.0_ireals      , & ! Constant of the MO theory (wind, convection)
  c_MO_t_conv   = 15.0_ireals      , & ! Constant of the MO theory (temperature, convection)
  c_MO_q_conv   = 15.0_ireals      , & ! Constant of the MO theory (humidity, convection)
  c_MO_u_exp    = 0.25_ireals      , & ! Constant of the MO theory (wind, exponent)
  c_MO_t_exp    = 0.5_ireals       , & ! Constant of the MO theory (temperature, exponent)
  c_MO_q_exp    = 0.5_ireals       , & ! Constant of the MO theory (humidity, exponent)
  z0u_ice_rough = 1.0E-03_ireals   , & ! Aerodynamic roughness of the ice surface [m] (rough flow)
  c_z0u_smooth  = 0.1_ireals       , & ! Constant in the expression for z0u (smooth flow) 
  c_z0u_rough   = 1.23E-02_ireals  , & ! The Charnock constant in the expression for z0u (rough flow)
  c_z0u_rough_L = 1.00E-01_ireals  , & ! An increased Charnock constant (used as the upper limit)
  c_z0u_ftch_f  = 0.70_ireals      , & ! Factor in the expression for fetch-dependent Charnock parameter
  c_z0u_ftch_ex = 0.3333333_ireals , & ! Exponent in the expression for fetch-dependent Charnock parameter
  c_z0t_rough_1 = 4.0_ireals       , & ! Constant in the expression for z0t (factor) 
  c_z0t_rough_2 = 3.2_ireals       , & ! Constant in the expression for z0t (factor)
  c_z0t_rough_3 = 0.5_ireals       , & ! Constant in the expression for z0t (exponent) 
  c_z0q_rough_1 = 4.0_ireals       , & ! Constant in the expression for z0q (factor)
  c_z0q_rough_2 = 4.2_ireals       , & ! Constant in the expression for z0q (factor)
  c_z0q_rough_3 = 0.5_ireals       , & ! Constant in the expression for z0q (exponent)
  c_z0t_ice_b0s = 1.250_ireals     , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b0t = 0.149_ireals     , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b1t = -0.550_ireals    , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b0r = 0.317_ireals     , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b1r = -0.565_ireals    , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b2r = -0.183_ireals    , & ! Constant in the expression for z0t over ice
  c_z0q_ice_b0s = 1.610_ireals     , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b0t = 0.351_ireals     , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b1t = -0.628_ireals    , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b0r = 0.396_ireals     , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b1r = -0.512_ireals    , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b2r = -0.180_ireals    , & ! Constant in the expression for z0q over ice
  Re_z0s_ice_t  = 2.5_ireals       , & ! Threshold value of the surface Reynolds number 
                                       ! used to compute z0t and z0q over ice (Andreas 2002)
  Re_z0u_thresh = 0.1_ireals           ! Threshold value of the roughness Reynolds number 
                                       ! [value from Zilitinkevich, Grachev, and Fairall (200),
                                       ! currently not used] 

!  Dimensionless constants 
REAL (KIND = ireals), PARAMETER ::   &
  c_free_conv   = 0.14_ireals          ! Constant in the expressions for fluxes in free convection

!  Dimensionless constants 
REAL (KIND = ireals), PARAMETER ::   &
  c_lwrad_emis  = 0.99_ireals          ! Surface emissivity with respect to the long-wave radiation

!  Thermodynamic parameters
REAL (KIND = ireals), PARAMETER ::        &
  tpsf_C_StefBoltz = 5.67E-08_ireals    , & ! The Stefan-Boltzmann constant [W m^{-2} K^{-4}]
  tpsf_R_dryair    = 2.8705E+02_ireals  , & ! Gas constant for dry air [J kg^{-1} K^{-1}]
  tpsf_R_watvap    = 4.6151E+02_ireals  , & ! Gas constant for water vapour [J kg^{-1} K^{-1}]
  tpsf_c_a_p       = 1.005E+03_ireals   , & ! Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
  tpsf_L_evap      = 2.501E+06_ireals   , & ! Specific heat of evaporation [J kg^{-1}]
  tpsf_nu_u_a      = 1.50E-05_ireals    , & ! Kinematic molecular viscosity of air [m^{2} s^{-1}]
  tpsf_kappa_t_a   = 2.20E-05_ireals    , & ! Molecular temperature conductivity of air [m^{2} s^{-1}]
  tpsf_kappa_q_a   = 2.40E-05_ireals        ! Molecular diffusivity of air for water vapour [m^{2} s^{-1}]

!  Derived thermodynamic parameters
REAL (KIND = ireals), PARAMETER ::                        &
  tpsf_Rd_o_Rv  = tpsf_R_dryair/tpsf_R_watvap           , & ! Ratio of gas constants (Rd/Rv)
  tpsf_alpha_q  = (1._ireals-tpsf_Rd_o_Rv)/tpsf_Rd_o_Rv     ! Diemsnionless ratio 

!  Thermodynamic parameters
REAL (KIND = ireals), PARAMETER ::     &
  P_a_ref             = 1.0E+05_ireals   ! Reference pressure [N m^{-2} = kg m^{-1} s^{-2}]


!  The variables declared below
!  are accessible to all program units of the MODULE "SfcFlx"
!  and to the driving routines that use "SfcFlx".
!  These are basically the quantities computed by SfcFlx.
!  Apart from these quantities, there a few local scalars 
!  used by SfcFlx routines mainly for security reasons.
!  All variables declared below have a suffix "sf".

!  SfcFlx variables of type REAL

!  Roughness lengths
REAL (KIND = ireals) ::    &
  z0u_sf                 , & ! Roughness length with respect to wind velocity [m]
  z0t_sf                 , & ! Roughness length with respect to potential temperature [m]
  z0q_sf                     ! Roughness length with respect to specific humidity [m]

!  Fluxes in the surface air layer
REAL (KIND = ireals) ::    &
  u_star_a_sf            , & ! Friction velocity [m s^{-1}]
  Q_mom_a_sf             , & ! Momentum flux [N m^{-2}]
  Q_sens_a_sf            , & ! Sensible heat flux [W m^{-2}]
  Q_lat_a_sf             , & ! Laten heat flux [W m^{-2}]
  Q_watvap_a_sf              ! Flux of water vapout [kg m^{-2} s^{-1}]

!  Security constants
REAL (KIND = ireals), PARAMETER ::   &
  u_wind_min_sf  = 1.0E-02_ireals  , & ! Minimum wind speed [m s^{-1}]
  u_star_min_sf  = 1.0E-04_ireals  , & ! Minimum value of friction velocity [m s^{-1}]
  c_accur_sf     = 1.0E-07_ireals  , & ! A small number (accuracy)
  c_small_sf     = 1.0E-04_ireals      ! A small number (used to compute fluxes)

!  Useful constants
REAL (KIND = ireals), PARAMETER ::     &
  num_1o3_sf = 1._ireals/3._ireals       ! 1/3

!==============================================================================
! Procedures 
!==============================================================================

CONTAINS

!==============================================================================
!  The codes of the SfcFlx procedures are stored in separate "*.incf" files
!  and are included below.
!------------------------------------------------------------------------------

!==============================================================================
include 'SfcFlx_lwradatm.incf'

!==============================================================================

!==============================================================================
include 'SfcFlx_lwradwsfc.incf'

!==============================================================================

!==============================================================================
include 'SfcFlx_momsenlat.incf'

!==============================================================================

!==============================================================================
include 'SfcFlx_rhoair.incf'

!==============================================================================

!==============================================================================
include 'SfcFlx_roughness.incf'

!==============================================================================

!==============================================================================
include 'SfcFlx_satwvpres.incf'

!==============================================================================

!==============================================================================
include 'SfcFlx_spechum.incf'

!==============================================================================

!==============================================================================
include 'SfcFlx_wvpreswetbulb.incf'

!==============================================================================

END MODULE SfcFlx

