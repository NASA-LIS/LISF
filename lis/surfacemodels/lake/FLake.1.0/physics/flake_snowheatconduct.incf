! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION flake_snowheatconduct (h_snow)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the snow heat conductivity,
!  using an empirical approximation from Heise et al. (2003).
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

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_parameters , ONLY : &
  tpl_rho_w_r               , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_kappa_S_min           , & ! Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_max           , & ! Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_Gamma_kappa_S             ! Empirical parameter [J m^{-2} s^{-1} K^{-1}] in the expression for kappa_S

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = ireals), INTENT(IN) :: &
  h_snow                              ! Snow thickness [m]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Snow heat conductivity [J m^{-1} s^{-1} K^{-1} = kg m s^{-3} K^{-1}]

  flake_snowheatconduct = flake_snowdensity( h_snow )   ! Compute snow density
  flake_snowheatconduct = MIN( tpl_kappa_S_max, tpl_kappa_S_min                      &
                        + h_snow*tpl_Gamma_kappa_S*flake_snowheatconduct/tpl_rho_w_r )

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_snowheatconduct

