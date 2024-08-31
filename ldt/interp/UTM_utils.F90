!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: UTM_utils
! 
! !DESCRIPTION: Performs conversions for UTM coordinates
!
! !SOURCE:  Defense Mapping Agency (DMA) Technical Manual 8358.2
!             "The Universal Grids:  Universal Transverse Mercator (UTM)
!                 and Universal Polar Stereographic (UPS)"
!             Defense Mapping Agency, Fairfax VA, September 1989
!           obtained on the web at http://earth-info.nga.mil/GandG/pubs.html
!
! !PURPOSE:  primarily for LIS-ARMS on Windows platforms
!
! !REVISION HISTORY:
!  28Oct2004  Matthew Garcia  Initial specification of UTM2Geo subroutine
!  29Oct2004  Matthew Garcia  Initial specification of Geo2UTM subroutine
!  04Dec2004  Matthew Garcia  Edit/update for incorporation in LISv3.2
!
!
module UTM_utils
!
!
  implicit none

  PRIVATE

  public    :: UTM2geo
!EOP 
! global module parameters
  real(4), parameter :: k_0 = 0.9996 ! UTM central scale factor
  real(8), parameter :: pi = 3.14159265359 ! if I gotta tell you...
  real(4), parameter :: FNNH = 0.0 ! False Northing for Northern Hemisphere [m]
  real(8), parameter :: FNSH = 10000000.0 ! False Northing for Southern Hemisphere [m]
  real(8), parameter :: FE = 500000.0 ! False Easting [m]
!
! WGS-84 geoid parameters
  real(8), parameter :: a = 6378137.00000 ! semi-major axis length [m]
  real(8), parameter :: invf = 298.257223563 ! inverse flattening (see below)
!
!
contains
!
!
real(8) function rad(deg)
!
! argument variables
  real(8), intent(IN) :: deg
!
  rad =(deg / 180.0) * pi
!
  return
end function rad
!
!
real(8) function deg(rad)
!
! argument variables
  real(8), intent(IN) :: rad
!
  deg = (rad / pi) * 180.0
!
  return
end function deg
!
!
real(8) function sinsq(phi)
!
! argument variables
  real(8), intent(IN) :: phi
!
  sinsq = sin(phi) * sin(phi)
!
  return
end function sinsq
!
!
real(8) function cossq(phi)
!
! argument variables
  real(8), intent(IN) :: phi
!
  cossq = cos(phi) * cos(phi)
!
  return
end function cossq
!
!
real(8) function tansq(phi)
!
! argument variables
  real(8), intent(INOUT) :: phi
!
  tansq = sinsq(phi) / cossq(phi)
!
  return
end function tansq
!
!
subroutine prime_func(f,A_prime,B_prime,C_prime,D_prime,E_prime)
!
! finds values of "prime" coefficients for calculation of S
!
! argument variables
  real(8), intent(IN) :: f
  real(8), intent(OUT) :: A_prime,B_prime,C_prime,D_prime,E_prime
!
! local variables
  real(8) :: n
!
  n = f / (2.0 - f)
  A_prime = a * (1 - n + (5.0/4.0)*(n**2 - n**3) + &
                 (81.0/64.0)*(n**4 - n**5))
  B_prime = (3.0/2.0) * a * (n - n**2 + (7.0/8.0)*(n**3 - n**4) + &
                             (55.0/64.0)*n**5)
  C_prime = (15.0/16.0) * a * (n**2 - n**3 + (3.0/4.0)*(n**4 - n**5))
  D_prime = (35.0/48.0) * a * (n**3 - n**4 + (11.0/16.0)*n**5)
  E_prime = (315.0/512.0) * a * (n**4 - n**5)
!
  return
end subroutine prime_func
!
!
real(8) function S_func(f,lat)
!
! !PURPOSE:  finds meridional arc (S) at latitude (lat)
!
! argument variables
  real(8), intent(INOUT) :: f
  real(8), intent(IN) :: lat
!
! local variables
  real(8) :: A_prime,B_prime,C_prime,D_prime,E_prime
  real(8) :: S_term1,S_term2,S_term3,S_term4,S_term5
!
  call prime_func(f,A_prime,B_prime,C_prime,D_prime,E_prime)
  S_term1 = A_prime * lat
  S_term2 = B_prime * sin(2*lat)
  S_term3 = C_prime * sin(4*lat)
  S_term4 = D_prime * sin(6*lat)
  S_term5 = E_prime * sin(8*lat)
  S_func = S_term1 - S_term2 + S_term3 - S_term4 + S_term5
!
  return
end function S_func
!
!
real(8) function T1_func(f,lat)
!
! !PURPOSE:  finds calculation term (T1)
!
! argument variables
  real(8), intent(INOUT) :: f,lat
!
  T1_func = S_func(f,lat) * k_0
!
  return
end function T1_func
!
!
real(8) function rho_func(fesq,lat)
!
! !PURPOSE:  finds radius of curvature in the meridian (rho)
!
! argument variables
  real(8), intent(IN) :: fesq
  real(8), intent(INOUT) :: lat
!
! local variables
  real(8) :: rho_num,rho_den
!
  rho_num = a * (1.0 - fesq)
  rho_den = (1.0 - fesq * sinsq(lat))**(3.0/2.0)
  rho_func = rho_num / rho_den 
!
  return
end function rho_func
!
!
real(8) function nu_func(fesq,sesq,lat)
!
! !PURPOSE:  finds radius of curvature in the prime vertical (nu)
!
! argument variables
  real(8), intent(IN) :: sesq
  real(8), intent(INOUT) :: fesq,lat
!
  nu_func = rho_func(fesq,lat) * (1.0 + sesq * cossq(lat))
!
  return
end function nu_func
!
!
real(8) function nor_kern_func(lat,nu)
!
! !PURPOSE:  finds northing calculation kernel (nor_kern)
!
! argument variables
  real(8), intent(IN) :: lat,nu
!
  nor_kern_func = nu * sin(lat) * cos(lat) * k_0 / 2
!
  return
end function nor_kern_func
!
!
real(8) function T2_func(lat,nu)
!
! !PURPOSE:  finds northing calculation term (T2)
!
! argument variables
  real(8), intent(INOUT) :: lat,nu
!
  T2_func = nor_kern_func(lat,nu)
!
  return
end function T2_func
!
!
real(8) function T3_func(lat,sesq,nu)
!
! !PURPOSE:  finds northing calculation term (T3)
!
! argument variables
  real(8), intent(INOUT) :: lat,nu
  real(8), intent(IN) :: sesq
!
! local variables
  real(8) :: T3_coeff,T3_term1,T3_term2,T3_term3,T3_term4
!
  T3_coeff = T2_func(lat,nu) * cossq(lat) / 12
  T3_term1 = 5
  T3_term2 = tansq(lat)
  T3_term3 = 9 * sesq * cossq(lat)
  T3_term4 = 4 * sesq**2 * cossq(lat)**2
  T3_func = T3_coeff * (T3_term1 - T3_term2 + T3_term3 + T3_term4)
!
  return
end function T3_func
!
!
real(8) function T4_func(lat,sesq,nu)
!
! !PURPOSE:  finds northing calculation term (T4)
!
! argument variables
  real(8), intent(INOUT) :: lat,nu
  real(8), intent(IN) :: sesq
!
! local variables
  real(8) :: T4_coeff,T4_term1,T4_term2,T4_term3
  real(8) :: T4_term4,T4_term5,T4_term6,T4_term7
  real(8) :: T4_term8,T4_term9,T4_term10,T4_term11
!
  T4_coeff = T2_func(lat,nu) * cossq(lat)**2 / 360
  T4_term1 = 61
  T4_term2 = 58 * tansq(lat)
  T4_term3 = tansq(lat)**2
  T4_term4 = 270 * sesq * cossq(lat)
  T4_term5 = 330 * tansq(lat) * sesq * cossq(lat)
  T4_term6 = 445 * sesq**2 * cossq(lat)**2
  T4_term7 = 324 * sesq**3 * cossq(lat)**3
  T4_term8 = 680 * tansq(lat) * sesq**2 * cossq(lat)**2
  T4_term9 = 88 * sesq**4 * cossq(lat)**4
  T4_term10 = 600 * tansq(lat) * sesq**3 * cossq(lat)**3
  T4_term11 = 192 * tansq(lat) * sesq**4 * cossq(lat)**4
  T4_func = T4_coeff * (T4_term1 - T4_term2 + T4_term3 + &
                        T4_term4 - T4_term5 + T4_term6 + &
						T4_term7 - T4_term8 + T4_term9 - &
						T4_term10 - T4_term11)
!
  return
end function T4_func
!
!
real(8) function T5_func(lat,nu)
!
! !PURPOSE:  finds northing calculation term (T5)
!
! argument variables
  real(8), intent(INOUT) :: lat,nu
!
! local variables
  real(8) :: T5_coeff,T5_term1,T5_term2,T5_term3,T5_term4
!
  T5_coeff = T2_func(lat,nu) * cossq(lat)**3 / 20160
  T5_term1 = 1385
  T5_term2 = 3111 * tansq(lat)
  T5_term3 = 543 * tansq(lat)**2
  T5_term4 = tansq(lat)**3
  T5_func = T5_coeff * (T5_term1 - T5_term2 + T5_term3 - T5_term4)
!
  return
end function T5_func
!
!
real(8) function eas_kern_func(lat,nu)
!
! !PURPOSE:  finds easting calculation kernel (eas_kern)
!
! argument variables
  real(8), intent(IN) :: lat,nu
!
  eas_kern_func = nu * cos(lat) * k_0
!
  return
end function eas_kern_func
!
!
real(8) function T6_func(lat,nu)
!
! !PURPOSE:  finds easting calculation term (T6)
!
! argument variables
  real(8), intent(INOUT) :: lat,nu
!
  T6_func = eas_kern_func(lat,nu)
!
  return
end function T6_func
!
!
real(8) function T7_func(lat,sesq,nu)
!
! !PURPOSE:  finds easting calculation term (T7)
!
! argument variables
  real(8), intent(INOUT) :: lat,nu
  real(8), intent(IN) :: sesq
!
! local variables
  real(8) :: T7_coeff,T7_term1,T7_term2,T7_term3
!
  T7_coeff = T6_func(lat,nu) * cossq(lat) / 6
  T7_term1 = 1.0
  T7_term2 = tansq(lat)
  T7_term3 = sesq * cossq(lat)
  T7_func = T7_coeff * (T7_term1 - T7_term2 + T7_term3)
!
  return
end function T7_func
!
!
real(8) function T8_func(lat,sesq,nu)
!
! !PURPOSE:  finds easting calculation term (T8)
!
! argument variables
  real(8), intent(INOUT) :: lat,nu
  real(8), intent(IN) :: sesq
!
! local variables
  real(8) :: T8_coeff,T8_term1,T8_term2,T8_term3
  real(8) :: T8_term4,T8_term5,T8_term6,T8_term7
  real(8) :: T8_term8,T8_term9
!
  T8_coeff = T6_func(lat,nu) * cossq(lat)**2 / 120
  T8_term1 = 5.0
  T8_term2 = 18 * tansq(lat)
  T8_term3 = tansq(lat)**2
  T8_term4 = 14 * sesq * cossq(lat)
  T8_term5 = 58 * tansq(lat) * sesq * cossq(lat)
  T8_term6 = 13 * sesq**2 * cossq(lat)**2
  T8_term7 = 4 * sesq**3 * cossq(lat)**3
  T8_term8 = 64 * tansq(lat) * sesq**2 * cossq(lat)**2
  T8_term9 = 24 * tansq(lat) * sesq**3 * cossq(lat)**3
  T8_func = T8_coeff * (T8_term1 - T8_term2 + T8_term3 + &
                        T8_term4 - T8_term5 + T8_term6 + &
						T8_term7 - T8_term8 - T8_term9)
!
  return
end function T8_func
!
!
real(8) function T9_func(lat,nu)
!
! !PURPOSE:  finds easting calculation term (T9)
!
! argument variables
  real(8), intent(INOUT) :: lat,nu
!
! local variables
  real(8) :: T9_coeff,T9_term1,T9_term2,T9_term3,T9_term4
!
  T9_coeff = T6_func(lat,nu) * cossq(lat)**3 / 5040
  T9_term1 = 61.0
  T9_term2 = 479 * tansq(lat)
  T9_term3 = 179 * tansq(lat)**2
  T9_term4 = tansq(lat)**3
  T9_func = T9_coeff * (T9_term1 - T9_term2 + T9_term3 - T9_term4)
!
  return
end function T9_func
!
!
real(8) function lat_kern_func(rho,nu)
!
! !PURPOSE:  finds lat calculation kernel (lat_kern)
!
! argument variables
  real(8), intent(IN) :: rho,nu
!
  lat_kern_func = 2 * rho * nu * k_0**2
!
  return
end function lat_kern_func
!
!
real(8) function T10_func(lat_prime,rho,nu)
!
! !PURPOSE:  finds lat calculation term (T10)
!
! argument variables
  real(8), intent(IN) :: lat_prime
  real(8), intent(INOUT) :: rho,nu
!
  T10_func = tan(lat_prime) / lat_kern_func(rho,nu)
!
  return
end function T10_func
!
!
real(8) function T11_func(lat_prime,sesq,rho,nu)
!
! !PURPOSE:  finds lat calculation term (T11)
!
! argument variables
  real(8), intent(INOUT) :: lat_prime,rho,nu
  real(8), intent(IN) :: sesq
!
! local variables
  real(8) :: T11_coeff,T11_term1,T11_term2,T11_term3
  real(8) :: T11_term4,T11_term5
!
  T11_coeff = T10_func(lat_prime,rho,nu) / (12 * nu**2 * k_0**2)
  T11_term1 = 5.0
  T11_term2 = 3 * tansq(lat_prime)
  T11_term3 = sesq * cossq(lat_prime)
  T11_term4 = 4 * sesq**2 * cossq(lat_prime)**2
  T11_term5 = 9 * tansq(lat_prime) * sesq * cossq(lat_prime)
  T11_func = T11_coeff * (T11_term1 + T11_term2 + T11_term3 - &
                          T11_term4 - T11_term5)
!
  return
end function T11_func
!
!
real(8) function T12_func(lat_prime,sesq,rho,nu)
!
! !PURPOSE:  finds lat calculation term (T12)
!
! argument variables
  real(8), intent(INOUT) :: lat_prime,rho,nu
  real(8), intent(IN) :: sesq
!
! local variables
  real(8) :: T12_coeff,T12_term1,T12_term2,T12_term3,T12_term4
  real(8) :: T12_term5,T12_term6,T12_term7,T12_term8,T12_term9
  real(8) :: T12_term10,T12_term11,T12_term12,T12_term13
!
  T12_coeff = T10_func(lat_prime,rho,nu) / (360 * nu**4 * k_0**4)
  T12_term1 = 61.0
  T12_term2 = 90 * tansq(lat_prime)
  T12_term3 = 46 * sesq * cossq(lat_prime)
  T12_term4 = 45 * tansq(lat_prime)**2
  T12_term5 = 252 * tansq(lat_prime) * sesq * cossq(lat_prime)
  T12_term6 = 3 * sesq**2 * cossq(lat_prime)**2
  T12_term7 = 100 * sesq**3 * cossq(lat_prime)**3
  T12_term8 = 66 * tansq(lat_prime) * sesq**2 * cossq(lat_prime)**2
  T12_term9 = 90 * tansq(lat_prime)**2 * sesq * cossq(lat_prime)
  T12_term10 = 88 * sesq**4 * cossq(lat_prime)**4
  T12_term11 = 225 * tansq(lat_prime)**2 * sesq**2 * cossq(lat_prime)**2
  T12_term12 = 84 * tansq(lat_prime) * sesq**3 * cossq(lat_prime)**3
  T12_term13 = 192 * tansq(lat_prime) * sesq**4 * cossq(lat_prime)**4
  T12_func = T12_coeff * (T12_term1 + T12_term2 + T12_term3 + &
                          T12_term4 - T12_term5 - T12_term6 + &
						  T12_term7 - T12_term8 - T12_term9 + &
						  T12_term10 + T12_term11 + T12_term12 - &
						  T12_term13)
!
  return
end function T12_func
!
!
real(8) function T13_func(lat_prime,rho,nu)
!
! !PURPOSE:  finds lat calculation term (T13)
!
! argument variables
  real(8), intent(INOUT) :: lat_prime,rho,nu
!
! local variables
  real(8) :: T13_coeff
  real(8) :: T13_term1,T13_term2,T13_term3,T13_term4
!
  T13_coeff = T10_func(lat_prime,rho,nu) / (20160 * nu**6 * k_0**6)
  T13_term1 = 1385.0
  T13_term2 = 3633 * tansq(lat_prime)
  T13_term3 = 4095 * tansq(lat_prime)**2
  T13_term4 = 1575 * tansq(lat_prime)**3
  T13_func = T13_coeff * (T13_term1 + T13_term2 + T13_term3 + T13_term4)
!
  return
end function T13_func
!
!
real(8) function lon_kern_func(lat_prime,nu)
!
! !PURPOSE:  finds lon calculation kernel (lon_kern)
!
! argument variables
  real(8), intent(IN) :: lat_prime,nu
!
  lon_kern_func = nu * cos(lat_prime) * k_0
!
  return
end function lon_kern_func
!
!
real(8) function T14_func(lat_prime,nu)
!
! !PURPOSE:  finds lon calculation term (T14)
!
! argument variables
  real(8), intent(INOUT) :: lat_prime,nu
!
  T14_func = 1 / lon_kern_func(lat_prime,nu)
!
  return
end function T14_func
!
!
real(8) function T15_func(lat_prime,sesq,nu)
!
! !PURPOSE:  finds lon calculation term (T15)
!
! argument variables
  real(8), intent(INOUT) :: lat_prime,nu
  real(8), intent(IN) :: sesq
!
! local variables
  real(8) :: T15_coeff,T15_term1,T15_term2,T15_term3
!
  T15_coeff = T14_func(lat_prime,nu) / (6 * nu**2 * k_0**2)
  T15_term1 = 1.0
  T15_term2 = 2 * tansq(lat_prime)
  T15_term3 = sesq * cossq(lat_prime)
  T15_func = T15_coeff * (T15_term1 + T15_term2 + T15_term3)
!
  return
end function T15_func
!
!
real(8) function T16_func(lat_prime,sesq,nu)
!
! !PURPOSE:  finds lon calculation term (T16)
!
! argument variables
  real(8), intent(INOUT) :: lat_prime,nu
  real(8), intent(IN) :: sesq
!
! local variables
  real(8) :: T16_coeff,T16_term1,T16_term2,T16_term3,T16_term4
  real(8) :: T16_term5,T16_term6,T16_term7,T16_term8,T16_term9
!
  T16_coeff = T14_func(lat_prime,nu) / (120 * nu**4 * k_0**4)
  T16_term1 = 5.0
  T16_term2 = 6 * sesq * cossq(lat_prime)
  T16_term3 = 28 * tansq(lat_prime)
  T16_term4 = 3 * sesq**2 * cossq(lat_prime)**2
  T16_term5 = 8 * tansq(lat_prime) * sesq * cossq(lat_prime)
  T16_term6 = 24 * tansq(lat_prime)**2
  T16_term7 = 4 * sesq**3 * cossq(lat_prime)**3
  T16_term8 = 4 * tansq(lat_prime) * sesq**2 * cossq(lat_prime)**2
  T16_term9 = 24 * tansq(lat_prime) * sesq**3 * cossq(lat_prime)**3
  T16_func = T16_coeff * (T16_term1 + T16_term2 + T16_term3 - &
                          T16_term4 + T16_term5 + T16_term6 - &
				  	      T16_term7 + T16_term8 + T16_term9)
!
  return
end function T16_func
!
!
real(8) function T17_func(lat_prime,nu)
!
! !PURPOSE:  finds lon calculation term (T17)
!
! argument variables
  real(8), intent(INOUT) :: lat_prime,nu
!
! local variables
  real(8) :: T17_coeff,T17_term1,T17_term2,T17_term3,T17_term4
!
  T17_coeff = T14_func(lat_prime,nu) / (5040 * nu**6 * k_0**6)
  T17_term1 = 61.0
  T17_term2 = 662 * tansq(lat_prime)
  T17_term3 = 1320 * tansq(lat_prime)**2
  T17_term4 = 720 * tansq(lat_prime)**3
  T17_func = T17_coeff * (T17_term1 + T17_term2 + &
                          T17_term3 + T17_term4)
!
  return
end function T17_func
!
!
subroutine UTM2Geo(zone,northing,easting,latdeg,londeg)
! 
! !DESCRIPTION: Converts UTM grid coordinates to geographic coordinates
!
  implicit none
!
! argument variables
  integer, intent(IN) :: zone
  real, intent(IN) :: northing,easting
  real, intent(OUT) :: latdeg,londeg
!
! general parameters
  real :: convg = 1.0 ! northing convergence criterion [m]
  integer :: maxiter = 100 ! maximum iterations for convergence
!
! local variables
  integer :: counter
  real(8) :: f,b,fesq,sesq,rho,nu
  real(8) :: norest,nordiff,norfrac
  real(8) :: FN,DeltaE,zone_lon
  real(8) :: convd,latdiff,latfrac,lat,lat2,lon
  real(8) :: lat_prime,lat_term2,lat_term3,lat_term4,lat_term5
  real(8) :: lon_nought,lon_term2,lon_term3,lon_term4,lon_term5
!
! calculate location-independent ellipsoid parameters
  f = 1 / invf ! flattening or ellipticity for Earth (a prolate spheroid)
  b = a * (1 - f) ! WGS84 geoid semi-minor axis length [m]
  fesq =  (a*a - b*b) / (a*a) ! first eccentricity, squared
  sesq =  fesq / (1.0 - fesq) ! second eccentricity, squared 
!
! find the footpoint latitude (lat_prime) by successive approximation
  counter = 0
  lat = 0.0 ! initial guess (start at Equator)
  do 
!
! determine appropriate value for False Northing
    if (lat.ge.0.0) then ! Northern Hemisphere
	  FN = FNNH
	else ! Southern Hemisphere
	  FN = FNSH
	end if
!
! calculate estimated Northing value
    norest = FN + T1_func(f,lat) ! proper formulation?
!
! calculate difference and assess convergence
    nordiff = norest - northing
	norfrac = nordiff / northing
    if (abs(nordiff).le.convg) then ! convergence accomplished
	  lat_prime = lat ! establish (lat_prime) for (lat) calculations
	  exit ! from do loop
	else ! convergence not yet accomplished
      lat = lat - (norfrac / 2) ! adjust lat estimate by bisection method
	end if ! convergence assessment
!
! determine if too many iterations
	counter = counter + 1
	if (counter.gt.maxiter) then
      print*,'MSG: UTM2Geo - convergence on lat_prime failed'
	  stop
	end if
  end do ! iteration to find (lat_prime)
!
! calculate DeltaE
  DeltaE = easting - FE
!
! find the actual latitude (lat) by successive approximation
  counter = 0
  lat = lat_prime ! initial guess
  convd = convg / (pi * (a + b)) ! latitude convergence criterion [deg]
  do
!
! calculate location-dependent ellipsoid parameters 
    rho = rho_func(fesq,lat)
    nu =  nu_func(fesq,sesq,lat)
!
! calculate latitude as in DMA TM 8358.2, section 2-6
    lat_term2 = DeltaE**2 * T10_func(lat_prime,rho,nu)
    lat_term3 = DeltaE**4 * T11_func(lat_prime,sesq,rho,nu)
	lat_term4 = DeltaE**6 * T12_func(lat_prime,sesq,rho,nu)
    lat_term5 = DeltaE**8 * T13_func(lat_prime,rho,nu)
    lat2 = lat_prime - lat_term2 + lat_term3 - lat_term4 + lat_term5
!
! calculate difference and assess convergence
    latdiff = lat - lat2
	latfrac = latdiff / lat
    if (abs(latdiff).le.convd) then ! convergence accomplished
	  exit ! from do loop
	else ! convergence not yet accomplished
	  lat = lat - (latdiff / 2) ! adjust lat estimate by bisection method
	end if ! convergence assessment
!
! determine if too many iterations
	counter = counter + 1
	if (counter.gt.maxiter) then
      print*,'MSG: UTM2Geo - convergence on lat failed'
	  stop
	end if
  end do ! iteration to find (lat)
!
! calculate longitude of projection origin (i.e. UTM zone center)
  zone_lon = -177.0 + 6.0 * (zone - 1)
  lon_nought = rad(zone_lon) ! centers of 6deg zones starting at 180degW
!
! calculate longitude as in DMA TM 8358.2, section 2-6
  lon_term2 = DeltaE * T14_func(lat_prime,nu)
  lon_term3 = DeltaE**3 * T15_func(lat_prime,sesq,nu)
  lon_term4 = DeltaE**5 * T16_func(lat_prime,sesq,nu)
  lon_term5 = DeltaE**7 * T17_func(lat_prime,nu)
  lon = lon_nought + lon_term2 - lon_term3 + lon_term4 - lon_term5
!
! convert lat and lon from radians to degrees
  latdeg = deg(lat)
  londeg = deg(lon)
!
  return
end subroutine UTM2Geo
!
!
subroutine Geo2UTM(latdeg,londeg,zone,northing,easting)
! 
! !DESCRIPTION: Converts geographic coordinates to UTM grid coordinates
!
  implicit none
!
! argument variables
  real(4), intent(IN) :: latdeg,londeg
  integer, intent(OUT) :: zone
  real(4), intent(OUT) :: northing,easting
!
! local variables
  real(8) :: f,b,fesq,sesq,nu
  real(8) :: FN,DeltaL,lon_nought
  real(8) :: lat,lon,zone_lon
  real(8) :: nor_term1,nor_term2,nor_term3,nor_term4,nor_term5
  real(8) :: eas_term1,eas_term2,eas_term3,eas_term4
!
! calculate location-independent ellipsoid parameters
  f = 1 / invf ! flattening or ellipticity for Earth (a prolate spheroid)
  b = a * (1 - f) ! WGS84 geoid semi-minor axis length [m]
  fesq =  (a*a - b*b) / (a*a) ! first eccentricity, squared
  sesq =  fesq / (1.0 - fesq) ! second eccentricity, squared 
!
! convert lat and lon from degrees to radians
  lat = rad(dble(latdeg))
  lon = rad(dble(londeg))
!
! determine appropriate value for False Northing
  if (latdeg.ge.0.0) then ! Northern Hemisphere
	FN = FNNH
  else ! Southern Hemisphere
	FN = FNSH
  end if
!
! find zone number (1..60)
  zone = 1
  do
    zone_lon = -180.0 + 6.0 * (zone - 1)
    if (londeg.ge.zone_lon.and. &
	    londeg.lt.(zone_lon + 6.0)) then
	  lon_nought = zone_lon + 3.0 ! centers of 6deg zones starting at 180degW
	  exit
	else
	  zone = zone + 1
	end if
  end do ! find zone number
!
! calculate DeltaL
  DeltaL = lon - rad(lon_nought)
!
! calculate location-dependent ellipsoid parameter 
  nu =  nu_func(fesq,sesq,lat)
!
! calculate northing as in DMA TM 8358.2, section 2-5
  nor_term1 = T1_func(f,lat)
  nor_term2 = DeltaL**2 * T2_func(lat,nu)
  nor_term3 = DeltaL**4 * T3_func(lat,sesq,nu)
  nor_term4 = DeltaL**6 * T4_func(lat,sesq,nu)
  nor_term5 = DeltaL**8 * T5_func(lat,nu)
  northing = FN + nor_term1 + nor_term2 + nor_term3 + nor_term4 + nor_term5
!
! calculate easting as in DMA TM 8358.2, section 2-5
  eas_term1 = DeltaL * T6_func(lat,nu)
  eas_term2 = DeltaL**3 * T7_func(lat,sesq,nu)
  eas_term3 = DeltaL**5 * T8_func(lat,sesq,nu)
  eas_term4 = DeltaL**7 * T9_func(lat,nu)
  easting = FE + eas_term1 + eas_term2 + eas_term3 + eas_term4
!
  return
end subroutine Geo2UTM
!
!
end module UTM_utils
