! cable_albedo.f90
!
! Source file containing albedo routine for CABLE
!
! Development by Ying-Ping Wang, Eva Kowalczyk, Ray Leuning
! Gab Abramowitz, Martin Dix, Harvey Davies, Mike Raupach, 
!
! bugs to bernard.pak@csiro.au
!
! This file contains modules:
!   albedo_module
! The subroutines included are:
!   surface_albedo
!
! Most user-defined types (e.g. met%tk) are defined in cable_types module
! in cable_variables.f90

MODULE cable_albedo
  USE cable_dimensions, ONLY: r_1,i_d,mp_patch,nrb
  USE cable_other_constants, ONLY: refl,taul
  USE cable_types
  IMPLICIT NONE
  ! This module contains the following subroutines:
  PRIVATE
  PUBLIC surface_albedo
CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE surface_albedo(ssoil, met, rad, canopy)
    TYPE (soil_snow_type),INTENT(INOUT)   :: ssoil
    TYPE (met_type),INTENT(INOUT)         :: met
    TYPE (radiation_type),INTENT(INOUT)   :: rad
    TYPE (canopy_type),INTENT(INOUT)      :: canopy
    REAL(r_1), DIMENSION(nrb) :: c1 ! sqrt(1. - taul - refl)                                       
    LOGICAL, DIMENSION(mp_patch)   :: mask ! select points for calculation               
    REAL(r_1), DIMENSION(mp_patch,nrb) :: rhocbm ! modified canopy beam reflectance (6.21) 
    INTEGER(i_d) :: b ! rad. band 1=visible, 2=near-infrared, 3=long-wave
    REAL(r_1), DIMENSION(nrb)    :: rhoch ! canopy reflection black horizontal leaves(6.19)

    ! Initialise effective conopy beam reflectance:
    rad%reffbm = ssoil%albsoilsn
    rad%reffdf = ssoil%albsoilsn
    
    ! Define vegetation mask:
    mask = canopy%vlaiw > 1e-2 .AND. met%fsd > 1.0e-2
    
    c1 = SQRT(1. - taul - refl)
    ! Define canopy reflection black horizontal leaves(6.19)
    rhoch = (1.0 - c1) / (1.0 + c1)
    ! Update extinction coefficients and fractional transmittance for 
    ! leaf transmittance and reflection (ie. NOT black leaves):
    DO b = 1, 2 ! 1 = visible, 2 = nir radiaition
       
       rad%extkdm(:,b) = rad%extkd * c1(b)
       ! Define canopy diffuse transmittance (fraction):
       rad%cexpkdm(:,b) = EXP(-rad%extkdm(:,b) * canopy%vlaiw)
       ! Calculate effective diffuse reflectance (fraction):
       rad%reffdf(:,b) = rad%rhocdf(:,b) + (ssoil%albsoilsn(:,b) &
            - rad%rhocdf(:,b)) * rad%cexpkdm(:,b)**2     
       WHERE (mask) ! i.e. vegetation and sunlight are present
          rad%extkbm(:,b) = rad%extkb * c1(b)
          ! Canopy reflection (6.21) beam:
          rhocbm(:,b) = 2.*rad%extkb/(rad%extkb+rad%extkd)*rhoch(b)
          ! Canopy beam transmittance (fraction):
          rad%cexpkbm(:,b) = EXP(-rad%extkbm(:,b)*canopy%vlaiw)
          ! Calculate effective beam reflectance (fraction):
          rad%reffbm(:,b) = rhocbm(:,b) + (ssoil%albsoilsn(:,b) &
               - rhocbm(:,b))*rad%cexpkbm(:,b)*rad%cexpkbm(:,b)
       END WHERE
       ! Define albedo:
       rad%albedo(:,b) = (1.0-rad%fbeam)*rad%reffdf(:,b) &
            + rad%fbeam*rad%reffbm(:,b)
    END DO
    ! Define IR albedo - CURRENTLY NOT USED elsewhere
    rad%albedo(:,3) = 0.05
    
  END SUBROUTINE surface_albedo

END MODULE cable_albedo
