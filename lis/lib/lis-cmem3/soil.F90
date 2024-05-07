!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!SUBROUTINE cmem_SOIL

! Purpose :
!   Calculate soil layer  emissivity using Dobson et al. 1985. 
! Interface :
! ---------
! Reference : CMEM
! ---------

! Author :
! ------
!   23-Aug-2006 Thomas Holmes   *ECMWF*
! Revised:      
!     Patricia de Rosnay  August 2007
!                         Nov 2008 (Tskin for water bodies: lakes and sea)

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
!   vtype        land cover type (UMD, 0.0-13.0) 
!   tsoil        top layer soil temperature (K)
!   wc           volumetric soil water content (cm3/cm3)
!   sr	         effective surface roughness height (m)
!   sand         sand fraction (% of dry weight of soil)
!   clay         clay fraction (% of dry weight of soil)

! Output Variables:
!   soil_em      Soil surface emissivity (H and V) 
!   tb_soil      Tb from soil (H and V) 

! Internal Variables: 
!   ef 	         diel. const. of frozen soil (Hallikainen et al. 1984, 1985)
!   r_s          reflectivity of smooth surface (h-pol and v-pol)
!   r_r          reflectivity of rough surface (h-pol and v-pol)
!   sal_soil     soil water salinity (psu = ppt(weight) ) 
!------------------------------------------------------------------------------

subroutine cmem_soil(fghz, theta, &
                     vtype, &
                     tsoil, wc, sr, sand, clay, &
                     soil_em, tb_soil)

implicit none

integer :: vtype
real :: fghz, theta, tsoil, wc, sr, sand, clay, soil_em(2), tb_soil(2)

real, parameter :: tfreeze = 273.15
real, parameter :: sal_soil = 0.0 
complex, parameter :: ef = (5.0, 0.5)

integer :: i 
integer :: medium, isal
real   :: frostfrac, tc
real   :: fsr           ! functions 
real   :: t_eff(2)      !  Effective temperature of soil medium
complex :: ew, eps      ! diel. const of water and soil
real   :: r_r(2), r_s(2)

!------------------------------------------------------------------------------
! 0.0 Initialize
!------------------------------------------------------------------------------

soil_em = 0.0 
tb_soil = 0.0 

! 1. roughness height (m) as empirical function of land cover type 

!sr = fsr (vtype) 

! 2B. Soil tiles
!------------------------------------------------------------------------------
! 2B.1 Effective temperature of soil medium
 t_eff(:) = tsoil  
 tc = t_eff(1) - tfreeze

! dielectric constant of surface medium  
  
 if (tc < -0.5 ) then 
     CALL DIEL_ICE (fghz, tc+tfreeze, ew)
 else 
     isal = 2   ! 1 Stogryn model 2  Klein and Swift model 
     medium = 2  ! 0 pure water, 1 sea water 2 soil water 
     CALL DIEL_WAT (fghz, tc, sal_soil, wc, sand, clay, medium, isal, ew)

 endif
 CALL DIELDOBSON (ew, wc, sand, clay, eps)

#if (defined DEBUG)
write(*, *) "cmem_soil.F90: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"                fghz        tc      diele. c. of soil water   diele. c. of soil" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(A10, 2F10.3, 2(A, F7.1, A, F7.1, A))')"fghz:tc:e", fghz, tc, & 
          '     (', real(ew) , ' + i', imag(ew),  ') ', &
          '     (', real(eps) , ' + i', imag(eps),  ') '
write(*, *)
#endif

 
 ! mix dielectric constant of frozen and non-frozen soil
 IF (tc < -5.) THEN ! frozen soil
   frostfrac = 1.
 ELSEIF (tc < -0.5) THEN
   frostfrac = 0.5
 ELSE 
   frostfrac = 0.
 ENDIF
 eps = eps * (1.-frostfrac) + ef * frostfrac 
 

! 2B.3. Compute Smooth Surface Reflectivity
 CALL FRESNEL (theta, eps, r_s)            
   
! 4. Compute Rough Surface Reflectivity
 CALL RGHWEGM(fghz, theta, sr, r_s, r_r)

! Try Choudahury-Wang, to compare with CRTM: 
! call RGHCHOU(fghz, theta, 0.0, 0.0, 0.0, r_s, r_r)

! 5. Compute Surface Emissivity 
DO i = 1,2    ! h- and v-polarization
  soil_em(i) = 1.0 - r_r(i)
ENDDO

! 6. Compute surface brightness temperatures
!     --------------------------------------
DO i = 1,2    ! h- and v-polarization
  tb_soil(i) = t_eff(i) * soil_em(i)
ENDDO

return 
 
END SUBROUTINE CMEM_SOIL

! Empirical function to get soil surface roughness height (m) 
! Usually a few cm

function fsr(vtype)
 implicit none
 integer vtype
 real :: fsr, nsr(0:13)   ! UMD land cover
 ! tentative values
 nsr = (/ 0.0,          &    ! 0     WATER
          0.05,         &    ! 1     EVERGREEN NEEDLELEAF FOREST
          0.05,         &    ! 2     EVERGREEN BROADLEAF FOREST
          0.05,         &    ! 3     DECIDUOUS NEEDLELEAF FOREST
          0.05,         &    ! 4     DECIDUOUS BROADLEAF FOREST
          0.05,         &    ! 5     MIXED FOREST
          0.05,         &    ! 6     WOODLAND
          0.05,         &    ! 7     WOODED GRASSLAND
          0.03,         &    ! 8     CLOSED SHRUBLAND
          0.03,         &    ! 9     OPEN SHRUBLAND
          0.02,         &    ! 10     GRASSLAND
          0.02,         &    ! 11     CROPLAND
          0.01,         &    ! 12     BARE GROUND
          0.03  /)           ! 13     URBAN AND BUILT-UP

  select case (vtype)
    case (0:13)
      fsr = nsr(vtype)
    case default
      fsr =  0.0
  end select

 return
end


! other functions to support Choudahury and Wang's roughness model.

function fNrH(vtype)
 implicit none
 integer vtype 
 real :: fNrH, nrh(0:13)   ! UMD land cover
 ! tentative values
 nrh = (/ 0.0,         &    ! 0     WATER 
          1.75,        &    ! 1     EVERGREEN NEEDLELEAF FOREST                
          1.0,         &    ! 2     EVERGREEN BROADLEAF FOREST                 
          1.0,         &    ! 3     DECIDUOUS NEEDLELEAF FOREST                
          1.0,         &    ! 4     DECIDUOUS BROADLEAF FOREST                 
          1.0,         &    ! 5     MIXED FOREST
          1.0,         &    ! 6     WOODLAND
          1.0,         &    ! 7     WOODED GRASSLAND                              
          1.0,         &    ! 8     CLOSED SHRUBLAND                          
          1.0,         &    ! 9     OPEN SHRUBLAND                           
          1.0,         &    ! 10     GRASSLAND                                
          0.0,         &    ! 11     CROPLAND                        
          0.0,         &    ! 12     BARE GROUND
          1.0  /)           ! 13     URBAN AND BUILT-UP            
  
  select case (vtype)  
    case (0:13)
      fNrH = nrh(vtype)
    case default 
      fNrH =  0.0
  end select 

 return 
end 
 
  
function fNrV(vtype)
 implicit none
 integer vtype
 real :: fNrV, nrv(0:13)   ! UMD land cover
 ! tentative values
 nrv = (/ 0.0,         &    ! 0     WATER
          0.0,         &    ! 1     EVERGREEN NEEDLELEAF FOREST
          0.0,         &    ! 2     EVERGREEN BROADLEAF FOREST
          2.0,         &    ! 3     DECIDUOUS NEEDLELEAF FOREST
          2.0,         &    ! 4     DECIDUOUS BROADLEAF FOREST
          0.0,         &    ! 5     MIXED FOREST
          0.0,         &    ! 6     WOODLAND
          0.0,         &    ! 7     WOODED GRASSLAND
          0.0,         &    ! 8     CLOSED SHRUBLAND
          0.0,         &    ! 9     OPEN SHRUBLAND
          0.0,         &    ! 10     GRASSLAND
         -1.0,         &    ! 11     CROPLAND
         -1.0,         &    ! 12     BARE GROUND
         -1.0  /)           ! 13     URBAN AND BUILT-UP

  select case (vtype)
    case (0:13)
      fNrV = nrv(vtype)
    case default
      fNrV =  0.0
  end select

 return
end

function fhrmodel(vtype)
 implicit none
 integer vtype
 real :: fhrmodel, hrm(0:13)   ! UMD land cover
 ! tentative values
 hrm = (/ 1.0,         &    ! 0     WATER
          1.0,         &    ! 1     EVERGREEN NEEDLELEAF FOREST
          1.0,         &    ! 2     EVERGREEN BROADLEAF FOREST
          1.0,         &    ! 3     DECIDUOUS NEEDLELEAF FOREST
          1.0,         &    ! 4     DECIDUOUS BROADLEAF FOREST
          1.0,         &    ! 5     MIXED FOREST
          1.0,         &    ! 6     WOODLAND
          1.0,         &    ! 7     WOODED GRASSLAND
          1.0,         &    ! 8     CLOSED SHRUBLAND
          1.0,         &    ! 9     OPEN SHRUBLAND
          1.0,         &    ! 10     GRASSLAND
          1.0,         &    ! 11     CROPLAND
          1.0,         &    ! 12     BARE GROUND
          1.0  /)           ! 13     URBAN AND BUILT-UP

  select case (vtype)
    case (0:13)
      fhrmodel = hrm(vtype)
    case default
      fhrmodel =  1.0
  end select

 return
end

