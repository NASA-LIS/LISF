! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flake_parameters

!------------------------------------------------------------------------------
!
! Description:
!
!  Values of empirical constants of the lake model FLake 
!  and of several thermodynamic parameters are set.
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

USE data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables 
  iintegers                    ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Dimensionless constants 
!  in the equations for the mixed-layer depth 
!  and for the shape factor with respect to the temperature profile in the thermocline
REAL (KIND = ireals), PARAMETER ::         &
  c_cbl_1       = 0.17_ireals            , & ! Constant in the CBL entrainment equation
  c_cbl_2       = 1._ireals              , & ! Constant in the CBL entrainment equation
  c_sbl_ZM_n    = 0.5_ireals             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
  c_sbl_ZM_s    = 10._ireals             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
  c_sbl_ZM_i    = 20._ireals             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
  c_relax_h     = 0.030_ireals           , & ! Constant in the relaxation equation for the SBL depth
  c_relax_C     = 0.0030_ireals              ! Constant in the relaxation equation for the shape factor
                                             ! with respect to the temperature profile in the thermocline

!  Parameters of the shape functions 
!  Indices refer to T - thermocline, S - snow, I - ice,
!  B1 - upper layer of the bottom sediments, B2 - lower layer of the bottom sediments.
!  "pr0" and "pr1" denote zeta derivatives of the corresponding shape function 
!  at "zeta=0" ad "zeta=1", respectively.
REAL (KIND = ireals), PARAMETER ::         &
  C_T_min       = 0.5_ireals             , & ! Minimum value of the shape factor C_T (thermocline)
  C_T_max       = 0.8_ireals             , & ! Maximum value of the shape factor C_T (thermocline)
  Phi_T_pr0_1   = 40._ireals/3._ireals   , & ! Constant in the expression for the T shape-function derivative 
  Phi_T_pr0_2   = 20._ireals/3._ireals   , & ! Constant in the expression for the T shape-function derivative 
  C_TT_1        = 11._ireals/18._ireals  , & ! Constant in the expression for C_TT (thermocline)
  C_TT_2        = 7._ireals/45._ireals   , & ! Constant in the expression for C_TT (thermocline)
  C_B1          = 2._ireals/3._ireals    , & ! Shape factor (upper layer of bottom sediments)
  C_B2          = 3._ireals/5._ireals    , & ! Shape factor (lower layer of bottom sediments)
  Phi_B1_pr0    = 2._ireals              , & ! B1 shape-function derivative 
  C_S_lin       = 0.5_ireals             , & ! Shape factor (linear temperature profile in the snow layer)
  Phi_S_pr0_lin = 1._ireals              , & ! S shape-function derivative (linear profile) 
  C_I_lin       = 0.5_ireals             , & ! Shape factor (linear temperature profile in the ice layer)
  Phi_I_pr0_lin = 1._ireals              , & ! I shape-function derivative (linear profile) 
  Phi_I_pr1_lin = 1._ireals              , & ! I shape-function derivative (linear profile) 
  Phi_I_ast_MR  = 2._ireals              , & ! Constant in the MR2004 expression for I shape factor
  C_I_MR        = 1._ireals/12._ireals   , & ! Constant in the MR2004 expression for I shape factor
  H_Ice_max     = 3._ireals                  ! Maximum ice tickness in 
                                             ! the Mironov and Ritter (2004, MR2004) ice model [m] 

!  Security constants
REAL (KIND = ireals), PARAMETER ::         &
  h_Snow_min_flk = 1.0E-5_ireals         , & ! Minimum snow thickness [m]
  h_Ice_min_flk  = 1.0E-9_ireals         , & ! Minimum ice thickness [m]
  h_ML_min_flk   = 1.0E-2_ireals         , & ! Minimum mixed-layer depth [m]
  h_ML_max_flk   = 1.0E+3_ireals         , & ! Maximum mixed-layer depth [m]
  H_B1_min_flk   = 1.0E-3_ireals         , & ! Minimum thickness of the upper layer of bottom sediments [m]
  u_star_min_flk = 1.0E-6_ireals             ! Minimum value of the surface friction velocity [m s^{-1}]

!  Security constant(s)
REAL (KIND = ireals), PARAMETER ::         &
  c_small_flk    = 1.0E-10_ireals            ! A small number

!  Thermodynamic parameters
REAL (KIND = ireals), PARAMETER ::        &
  tpl_grav          = 9.81_ireals       , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_r           = 277.13_ireals     , & ! Temperature of maximum density of fresh water [K]
  tpl_T_f           = 273.15_ireals     , & ! Fresh water freezing point [K]
  tpl_a_T           = 1.6509E-05_ireals , & ! Constant in the fresh-water equation of state [K^{-2}]
  tpl_rho_w_r       = 1.0E+03_ireals    , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_rho_I         = 9.1E+02_ireals    , & ! Density of ice [kg m^{-3}]
  tpl_rho_S_min     = 1.0E+02_ireals    , & ! Minimum snow density [kg m^{-3}]
  tpl_rho_S_max     = 4.0E+02_ireals    , & ! Maximum snow density [kg m^{-3}]
  tpl_Gamma_rho_S   = 2.0E+02_ireals    , & ! Empirical parameter [kg m^{-4}]  
                                            ! in the expression for the snow density 
  tpl_L_f           = 3.3E+05_ireals    , & ! Latent heat of fusion [J kg^{-1}]
  tpl_c_w           = 4.2E+03_ireals    , & ! Specific heat of water [J kg^{-1} K^{-1}]
  tpl_c_I           = 2.1E+03_ireals    , & ! Specific heat of ice [J kg^{-1} K^{-1}]
  tpl_c_S           = 2.1E+03_ireals    , & ! Specific heat of snow [J kg^{-1} K^{-1}]
  tpl_kappa_w       = 5.46E-01_ireals   , & ! Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_I       = 2.29_ireals       , & ! Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_min   = 0.2_ireals        , & ! Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_max   = 1.5_ireals        , & ! Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_Gamma_kappa_S = 1.3_ireals            ! Empirical parameter [J m^{-2} s^{-1} K^{-1}] 
                                            ! in the expression for the snow heat conductivity 

!==============================================================================

END MODULE flake_parameters

