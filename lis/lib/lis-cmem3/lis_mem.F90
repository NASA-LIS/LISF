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

! main subroutine to compute land surface emissivity and TOA Tb. 

! Input Variables (all are real scalars except vtype): 
!------------------- sensors -----------------------------
!   fghz         frequency (GHz)                         
!   theta        zenith angle (degree)
!------------------- atmospheric contribution ------------
!   tsurf        air temperature at land surface (K)
!   tau_atm0     atmospheric zenith optical depth
!   tb_ad        atmospheric downwelling Tb 
!   tb_au        atmospheric upwelling Tb 
!------------------- soil properties ---------------------
!   tsoil        top layer soil temperature (K)
!   wc           volumetric soil water content (cm3/cm3)
!   sr           effective surface roughness height (m)
!   sand         sand fraction (% of dry weight of soil)
!   clay         clay fraction (% of dry weight of soil)
!------------------- veg properties ----------------------
!   vtype        vegetation cover type, real (UMD class numbers, 0-13) 
!   t_veg        vegetation temperature (K) 
!   h_veg        vegetation height (m)
!   lai          leaf area index
!   wc_veg       vegetation water content (kg/m2)
!   rho_veg      vegetation density (kg/m3)
!   m_d          dry mass fraction of vegetation
!   d_leaf       leaf thickness (mm)
!------------------- snow properties ---------------------
!   sn_t         snow temperature (K)
!   sn_moist     snow moisture [cm3/cm3]  (CMEM default : 0.1 )
!   sn_density   snow density [g/cm^3]
!   sn_depth     snow depth [m]
!   sn_gsize     snow grain sizze [mm]
!---------------------------------------------------------

! Output Variables (all are real scalars):
!   emH          Land surface (bottom of atmosphere) emissivity, H-POL
!   emV          Land surface (bottom of atmosphere) emissivity, V-POL
!   toaH         Top of atmosphere Tb, H-POL 
!   toaV         Top of atmosphere Tb, V-POL 
!
! Internal Variables: 
!   tau_veg      effective vegetation opacity (index 1.h-pol, 2.v-pol)
!   tb_tov       top of land surface (veg) Tb ( H and V) 
!---------------------------------------------------------

subroutine lis_mem (                                           &
   fghz, theta,                                                & ! sensor
   tsurf, tau_atm0, tb_ad, tb_au,                              & ! atmosphere
   tsoil, wc, sr, sand, clay,                                  & ! soil 
   vtype, t_veg, h_veg, lai, wc_veg,                           & ! veg
   rho_veg, m_d, d_leaf,                                       & 
   sn_t, sn_moist, sn_density, sn_depth, sn_gsize,             & ! snow 
   emH, emV, toaH, toaV)                                         ! output

implicit none 

real :: fghz, theta
real :: tsurf, tau_atm0, tb_ad, tb_au
real :: tsoil, wc, sr, sand, clay
integer :: vtype
real :: t_veg, h_veg, lai, wc_veg, rho_veg, m_d, d_leaf
real :: sn_t, sn_moist, sn_density, sn_depth, sn_gsize
real :: emH, emV, toaH, toaV, t_boa       ! bottom of atmosphere temperature 

real :: pi, costheta
real :: tb_soil(2), tb_veg(2), tb_snow(2), tb_water(2), tb_toa(2), tb_tov(2)
real :: soil_em(2), snow_em(2), tau_veg(2), water_em(2) 
real :: tau_atm
real :: vfrac, rmm      ! desert params: volume fraction and grain size (mm)

pi = acos(-1.0)
costheta = cos(theta*pi/180.0)

tb_soil = 0.0
tb_veg  = 0.0
tb_snow = 0.0
tb_water = 0.0
tb_toa = 0.0
tb_tov = 0.0
soil_em = 0.0
snow_em = 0.0
tau_veg = 0.0     ! effective vegetation opacity (index 1.h-pol, 2.v-pol)

tau_atm = tau_atm0/costheta

!#  safeguarding the valid temperature range for calculation of dielectric constant 
!   of water  (0-40C) 

tsurf=min(tsurf, 273.15+40)
tsoil=min(tsoil, 273.15+40)
t_veg=min(t_veg, 273.15+40)


#if (defined DEBUG)
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *) "======================== Start of debug info =============================" 
write(*, *) &
"    vtype     tsoil        wc        sr      sand      clay     t_veg    wc_veg" 
write(*, '(I10.0, 7F10.3)')vtype, tsoil, wc, sr, sand, clay, t_veg, wc_veg
write(*, *) "lis_mem atm inputs: "  
write(*, *)  & 
"     fghz     theta  tau_atm0     tb_ad     tb_au " 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(5F10.3)')fghz, theta, tau_atm0, tb_ad, tb_au 
write(*, *)
#endif

If (vtype .GT. 0 ) then     ! not water.  
  ! soil layer
  call cmem_soil(fghz, theta, &
                     vtype, &
                     tsoil, wc, sr, sand, clay, &
                     soil_em, tb_soil)

#if (defined DEBUG)
write(*, *) "cmem_soil outputs: "  
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)   "    vtype     tsoil        wc        sr      sand      clay" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(I10.0, 5F10.3)') vtype, tsoil, wc, sr, sand, clay
write(*, *)
write(*, *)"  soil_emH  soil_emV  tb_soilH  tb_soilV sn_depth" 
write(*, *)"---------------------------------------------------------------------------"

write(*, '(5F10.3)') soil_em(1), soil_em(2), tb_soil(1), tb_soil(2), sn_depth 
write(*, *)
#endif

  
  if (sn_depth .GT. 0.001 ) then   ! snow-modified terrain emissivity and Tb
    call cmem_snow(fghz, theta, tsoil, &
                   sn_t, sn_moist, sn_density, sn_depth, sn_gsize, &    ! snow data
                   (1.0 - soil_em),  &            ! reflectivity between snow and soil
                   soil_em, tb_soil)
#if (defined DEBUG)
write(*, *) "cmem_snow inputs: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
" sn_depth  sn_gsize sn_moist   sn_t     sn_density"
write(*, *)"---------------------------------------------------------------------------"
write(*, '(5F10.3)') sn_depth, sn_gsize, sn_moist, sn_t, sn_density 
write(*, *)
write(*, *) "cmem_snow outputs: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
" snow_emH  snow_emV  tb_soilH  tb_soilV"
write(*, *)"---------------------------------------------------------------------------"
write(*, '(4F10.3)') soil_em(1), soil_em(2), tb_soil(1), tb_soil(2)
write(*, *)
#endif

  end if 

  t_boa = tsoil       ! bottom of atmosphere temp

  if (h_veg .GE. sn_depth) then   
    ! veg layer 
    call cmem_veg(fghz, theta, vtype, t_veg, h_veg-sn_depth, lai, &
                    wc_veg, &    ! veg params
                    rho_veg, m_d, d_leaf,  &
                    tau_veg, tb_veg)

    !YDT t_boa = t_veg       ! bottom of atmosphere temp

#if (defined DEBUG)
write(*, *) "cmem_veg outputs: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"    vtype     t_veg       lai  tau_vegH  tau_vegV   tb_vegH   tb_vegV" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(I10, 6F10.3)') vtype, t_veg, lai, tau_veg(1), tau_veg(2), tb_veg(1), tb_veg(2) 
write(*, *) 
#endif

  end if
 
  tb_tov = tb_soil * exp(-tau_veg)                              &     ! soil 
         + tb_veg + tb_veg * (1.0 - soil_em) * exp(-tau_veg)    &     ! veg
         + tb_ad  * (1.0 - soil_em) * exp(-2.0*tau_veg)               ! atmosphere
  
#if (defined DEBUG)
write(*, *) "lis_mem.F90: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"   tb_tovH   tb_tovV  tau_vegH  tau_vegV" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(4F10.3)') tb_tov(1), tb_tov(2), tau_veg(1), tau_veg(2) 
write(*, *)
write(*, *) "lis_mem TB components: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *) "                                   H-pol     V-pol"
write(*, '(A, 2F10.3)') & 
"Soil contributed TB:          ", tb_soil(1) * exp(-tau_veg(1)), &
                                  tb_soil(2) * exp(-tau_veg(2))
write(*, '(A, 2F10.3)') "Veg  contributed TB:          ", & 
     tb_veg(1) + tb_veg(1) * (1.0 - soil_em(1) ) * exp(-tau_veg(1)),  &
     tb_veg(2) + tb_veg(2) * (1.0 - soil_em(2) ) * exp(-tau_veg(2))
write(*, '(A, 2F10.3)') "Atm  contributed TB:          ", & 
        tb_ad * (1.0 - soil_em(1)) * exp(-2.0*tau_veg(1)), &
        tb_ad * (1.0 - soil_em(2)) * exp(-2.0*tau_veg(2))
write(*, *)"---------------------------------------------------------------------------"
write(*, '(A, 2F10.3)') "Total TOV TB:                 ", tb_tov(1), tb_tov(2) 
write(*, *)
#endif


Else    ! water. inland  fresh water body, smooth surface  

  tb_veg = 0.0
  tb_soil = 0.0
  tau_veg = 0.0
  call cmem_water(fghz, theta, tsurf, water_em, tb_water) 

  t_boa = tsurf       ! bottom of atmosphere temp

#if (defined DEBUG)
write(*, *) "cmem_water outputs: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"    tsurf water_emH water_emV tb_waterH tb_waterV" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(5F10.3)') tsurf, water_em(1), water_em(2), tb_water(1), tb_water(2) 
write(*, *)
#endif

  tb_tov = tb_water                                             &     ! soil 
         + tb_ad  * (1.0 - water_em)                                  ! atmosphere

End if 

! Output
toaH = tb_tov(1) * exp(-tau_atm) + tb_au
toaV = tb_tov(2) * exp(-tau_atm) + tb_au
emH = tb_tov(1) / t_boa 
emV = tb_tov(2) / t_boa

#if (defined DEBUG)
write(*, *) "lis-mem outputs: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"       emH       emV   tb_tovH   tb_tovV     t_boa" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(5F10.3)') emH, emV, tb_tov(1), tb_tov(2), t_boa
write(*, *)
#endif

return 

End subroutine lis_mem 

