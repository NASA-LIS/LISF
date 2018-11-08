! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION flake_snowdensity (h_snow)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the snow density,
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
  tpl_rho_S_min             , & ! Minimum snow density [kg m^{-3}]
  tpl_rho_S_max             , & ! Maximum snow density [kg m^{-3}]
  tpl_Gamma_rho_S           , & ! Empirical parameter [kg m^{-4}] in the expression for the snow density
  c_small_flk                   ! A small number

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

! Snow density [kg m^{-3}]

!  Security. Ensure that the expression in () does not become negative at a very large h_snow.
  flake_snowdensity = MAX( c_small_flk, (1._ireals - h_snow*tpl_Gamma_rho_S/tpl_rho_w_r) )
  flake_snowdensity = MIN( tpl_rho_S_max, tpl_rho_S_min/flake_snowdensity )

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_snowdensity 

