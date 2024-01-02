!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!SUBROUTINE cmem_desert 

! Purpose :
!   Calculate desert emissivity using Grody and Weng, 2008 
! Interface :
! ---------
! Reference : Grody and Weng, 2008, IEEE TGARS, 46, 361-375. 
! ---------

! Author :
! ------
!   21-Mar-2011 Yudong Tian

! Modifications :
! -------------
! End Modifications
!
! Input Variables (all are real scalars):
!------------------- sensors -----------------------------
!   fghz         frequency (GHz)
!   theta        zenith angle (degree)
!   tsoil        desert surfce temperature (K)
!   sr	         effective surface roughness height (m)
!   vfrac        volume fraction of desert sand (0-1) 
!   rmm          sand particle size (mm) 

! Output Variables:
!   soil_em      Soil surface emissivity (H and V) 
!   tb_soil      Tb from soil (H and V) 

! Internal Variables: 
!   r_s          reflectivity of smooth surface (h-pol and v-pol)
!   r_r          reflectivity of rough surface (h-pol and v-pol)
!   yR, yI       y parameters (Grody and Weng, 2008, Eq. A15) 
!   omega 	 single-particle albedo (ibid, Eq. A16) 
!   g     	 asymmetry parameter (ibid) 
!   a     	 asymmetry parameter (ibid, Eq. 3b) 
!------------------------------------------------------------------------------

subroutine cmem_desert(fghz, theta, tsoil, sr, vfrac, rmm, soil_em, tb_soil)

implicit none

real :: fghz, theta, tsoil, sr, vfrac, rmm, soil_em(2), tb_soil(2)
real :: pi, k, kr

! values used in (Grody and Weng 2008)
!real, parameter :: vfrac = 0.6 
!real, parameter :: rmm = 0.5 

integer :: i 
complex :: eps      ! diel. const of desert 
real   :: r_r(2), r_s(2)
real   :: yR, yI, omega, g, a, ftmp

!------------------------------------------------------------------------------
! 0.0 Initialize
!------------------------------------------------------------------------------
pi = acos(-1.0)
k = 2.0 * pi * fghz / 300.0    !  1/mm
kr = k * rmm

soil_em = 0.0 
tb_soil = 0.0 


! dielectric constant of desert, Matzler (1998)
 CALL diel_desert (fghz, eps)

! or Grody and Weng (2008)
!eps = (4.0, 0.08)

#if (defined DEBUG)
write(*, *) "cmem_desert.F90: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"                fghz         diele. c. of desert  " 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(A10, F10.3, A, F7.3, A, F7.3, A)') "fghz:e_drt", fghz, & 
          '     (', real(eps) , ' + i', imag(eps),  ') '
write(*, *)
#endif

! 2B.3. Compute Smooth Surface Reflectivity
 CALL FRESNEL (theta, eps, r_s)            
   
!--- skip rough surface correction, assuming desert is smooth and emissivity is already high
! 4. Compute Rough Surface Reflectivity
! CALL RGHWEGM(fghz, theta, sr, r_s, r_r)

! Try Choudahury-Wang, to compare with CRTM:
! call RGHdesert(fghz, theta, 0.0, 0.0, 0.0, r_s, r_r)

r_r = r_s

yR = ( Real(eps) - 1.0 ) / ( Real(eps) + 2.0)
yI = 3.0 * Aimag(eps)/ ( Real(eps) + 2.0)**2

ftmp = ( 1.0 - vfrac )**4
ftmp = ftmp * kr**3 * yR * yR 
omega = ftmp / (ftmp  + 1.5 * ( 1.0 + 2.0 * vfrac )**2 * yI )   ! Eq. A16
g = 0.23 * kr * kr
a = Sqrt ( ( 1.0 - omega) / (1.0 - omega * g) )    ! Eq. 3b

soil_em = ( 1.0 - r_r ) * ( 2.0 * a / ( ( 1.0 + a) - ( 1.0 -a )*r_r )  )

tb_soil = tsoil  * soil_em

return 
 
END SUBROUTINE cmem_desert

