!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!SUBROUTINE cmem_water 

! Purpose :
!   Calculate inland water  emissivity 
! Interface :
! ---------
! Reference : CMEM
! ---------

! Modifications :
! -------------
! End Modifications
!   alpha : fitting parameter (Wang and Schmugge, 1980)
! z : depth of center of layer [m]
! dz : layer thickness [m]
!   medium : choose pure water (0), sea water (1), soil water (2)
!   isal : choose model for diel_wat, stogryn (1), Klein&Swift (2)
!
! Input Variables (all are real scalars):
!------------------- sensors -----------------------------
!   fghz         frequency (GHz)
!   theta        zenith angle (degree)
!   tsurf        surface temperature (K)

! Output Variables:
!   water_em      water surface emissivity (H and V) 
!   tb_water      Tb from water (H and V) 

! Internal Variables: 
!   r_s          reflectivity of smooth surface (h-pol and v-pol)
!   r_r          reflectivity of rough surface (h-pol and v-pol)
!   sal_water     soil water salinity (psu = ppt(weight) ) 
!------------------------------------------------------------------------------

subroutine cmem_water(fghz, theta, tsurf,  &
                     water_em, tb_water)

implicit none

real :: fghz, theta, tsurf, water_em(2), tb_water(2)

real, parameter :: tfreeze = 273.15
real, parameter :: sal_water = 0.0 
complex, parameter :: ef = (5.0, 0.5)

integer :: i 
integer :: medium, isal
real   :: tc, NrH, NrV, hrmodel, sand, clay, wc 
real   :: t_eff(2)      !  Effective temperature of soil medium
complex :: ew     ! diel. const of water 
real   :: r_s(2)

!------------------------------------------------------------------------------
! 0.0 Initialize
!------------------------------------------------------------------------------

water_em = 0.0 
tb_water = 0.0 

! 1. roughness parameters as functions of veg type 

NrH = 0.0
NrV = 0.0 
hrmodel =  0.0 

! 2B. Soil tiles
!------------------------------------------------------------------------------
! 2B.1 Effective temperature of soil medium
 t_eff(:) = tsurf 
 tc = t_eff(1) - tfreeze

! dielectric constant of surface medium  
  
 if (tc < -0.5 ) then 
     CALL DIEL_ICE (fghz, tc+tfreeze, ew)
 else 
     isal = 2
     medium = 0     ! fresh water
     CALL DIEL_WAT (fghz, tc, sal_water, wc, sand, clay, medium, isal, ew)

#if (defined DEBUG)
write(*, *) "cmem_water.F90: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"        tc      diele. const of fresh water"
write(*, *)"---------------------------------------------------------------------------"
write(*, '(F10.3, 1(A, F7.1, A, F7.1, A))') tc,  &
          '      (', real(ew) , ' + i', imag(ew),  ') '
write(*, *)
#endif

 endif

! 2B.3. Compute Smooth Surface Reflectivity
 CALL FRESNEL (theta, ew, r_s)            
   
! 5. Compute Surface Emissivity 
DO i = 1,2    ! h- and v-polarization
  water_em(i) = 1.0 - r_s(i)
ENDDO

! 6. Compute surface brightness temperatures
!     --------------------------------------
DO i = 1,2    ! h- and v-polarization
  tb_water(i) = t_eff(i) * water_em(i)
ENDDO

return 
 
END SUBROUTINE CMEM_water 
