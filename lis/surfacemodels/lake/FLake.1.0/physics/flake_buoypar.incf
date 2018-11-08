! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION flake_buoypar (T_water)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the buoyancy parameter,
!  using a quadratic equation of state for the fresh-water.
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
  tpl_grav                  , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_r                   , & ! Temperature of maximum density of fresh water [K]
  tpl_a_T                       ! Constant in the fresh-water equation of state [K^{-2}]

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = ireals), INTENT(IN) :: &
  T_water                             ! Water temperature [K]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Buoyancy parameter [m s^{-2} K^{-1}]

  flake_buoypar = tpl_grav*tpl_a_T*(T_water-tpl_T_r)

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_buoypar

