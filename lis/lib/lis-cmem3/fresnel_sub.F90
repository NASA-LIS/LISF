!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!!!! Insert banner here !!!!
!SUBROUTINE FRESNEL 

! Purpose :
!  Fresnel Law to computes the reflectivities of flat surfaces

!  Reference:
!  Njoku and Kong, 1977: Theory for passive microwave remote sensing
!    of near-surface soil moisture.
!    Journal of Geophysical Research, Vol. 82, No. 20, 3108-3118.

! Input variables: 
! theta		incidence angle (degrees) 	
! eps_surf  	dielectric constant of the surface

! Output variable: 
! r_s(2) : reflectivity of smooth soil surface (-), index 1 = h-pol., 2 = v-pol.

!---------------------------------------------------------------------------
SUBROUTINE FRESNEL (theta, eps_surf, r_s )

IMPLICIT NONE
real :: theta, r_s(2) 
COMPLEX :: eps_surf

real :: sintheta, costheta, pi
COMPLEX :: g
!---------------------------------------------------------------------------

pi = acos(-1.0)
costheta = cos(theta*pi/180.0) 
sintheta = sin(theta*pi/180.0) 

g = sqrt( eps_surf - sintheta*sintheta )

r_s(1) = abs((costheta-g)/(costheta+g)) ** 2.
r_s(2) = abs((costheta*eps_surf-g)/(costheta*eps_surf+g)) ** 2.

END SUBROUTINE FRESNEL
