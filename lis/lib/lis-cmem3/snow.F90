!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!!!! insert banner here !!!!

!subroutine cmem_snow
! purpose :
! -------
!   hut-snow emission model for single snow layer
! reference :
! ---------
! jp/8.12.1997, last modification 14.4.1998
! jpw : traduction fortran f77 ( 02/01)
! author :
! ------
!   2-oct-2006 thomas holmes   *ecmwf*
!   january 2008 patricia de rosnay, ecmwf coding standards
! modifications :
! -------------
! end modifications

! internal variables: 
! not used :  tsky  = sky brightness temperature (k)
! ppq2mod : empirical parameter of modified hut model
! eps_ice : dielectric constant of ice
! eir : real part of ice epsilon 
! eii : imaginary part of ice epsilon 
! sal_snow : salinity of snow (ppm)

! input variables:
!   fghz         frequency (ghz)
!   theta        zenith angle (degree)
!   tsoil        top layer soil temperature (k) 
!   sn_t 	 snow temperature (k) 
!   sn_moist 	 snow moisture [cm3/cm3]  (cmem default : 0.1 ) 
!   sn_density   snow density [g/cm^3]
!   sn_depth     snow depth [m] 
!   sn_gsize     snow grain sizze [mm] 
!   rsn          reflectivity between the snow and ground at (1, h-pol. 2, v.)

! output variables: 
! esn            emissivity of snow covered terrain at 1: h-polarization 2: v-polarization
! tb_tos         snow covered terrain brightness temperature at h and v-polarization
!------------------------------------------------------------------------------

subroutine cmem_snow(fghz, theta, tsoil, & 
                       sn_t, sn_moist, sn_density, sn_depth, sn_gsize, &    ! snow data
                       rsn,  &            ! reflectivity between snow and soil
                       esn,  tb_tos)

implicit none
real, parameter :: clight = 2.998e8         ! light speed,  m/s
real, parameter :: tfreeze = 273.15 
real, parameter :: rhowat = 1000.0     ! water density (kg/m^3) 
real, parameter :: eps_winf = 4.9      ! diel. constant at infinite frequency (stogryn 1971)
real, parameter :: eps_0 = 8.854e-12   ! permittivity of free space 
                                       !(klein and swift 1977) [farads/meter]
real :: pi, fghz, theta, f, omega, k, sintheta, costheta
real :: rsn(2), tb_tos(2), esn(2) 
real :: tsoil, sn_t, sn_moist, sn_density, sn_depth, sn_gsize
real ::       tlumi, xmv     , pdsw      , dlumi,    ppdiam
complex :: cim, epsa, epsb, epsc, epsmarka, cn2, rf

integer :: i

real :: pice, y0
real, parameter :: ppq2mod = 0.96
real, parameter :: ppsal_snow = 0.
real :: ee(2), ysuora(2)
real :: pds, weq, reds, xn1
real :: a, b, c, ap, bp, cp, aa, ab, ac
real :: epssw, ff0

real :: eir, eii, eiis, eiip, deltaeii
real :: vi, xieds, f0a, f0b, f0c, epsinfa, epsinfb, epsinfc
real :: epssa, epssb, epssc, rews, xiews, alfa, beta
real :: k2x, ppp, qq, xjakaja, o2, tcoh, tsa, trans, sec,ked,ke
real :: kabsd, kabs, xinta, xinte, suhde, l2, la, xr1
real :: ks, l2apu, xml, emiss1, emiss2, sulapu, sul

complex :: eps_ice
!------------------------------------------------------------------------------
! mapping variable names   
xmv = sn_moist
pdsw = sn_density
dlumi = sn_depth
ppdiam = sn_gsize 
tlumi = sn_t-tfreeze

 pi = acos(-1.0)
 f = fghz * 1.0e9
 k = 2.0 * pi * f / clight
 omega = 2.0 * pi * f
 costheta = cos(theta*pi/180.0)
 sintheta = sin(theta*pi/180.0)


! snow temperature [c] set equal to the surface soil temperature
! constants:
 cim=(0.,1.)
 pice = 0.916 ! density of ice
 y0 = 4.*pi*1.e-7! myy nolla

pds = (pdsw - xmv)/(1.0 - xmv) ! density of dry snow
weq = rhowat * dlumi * pdsw   ! snow water equivalent
! transmissivity from ground to snow (e = 1 - [s0v s0h])
ee(1) = 1. - rsn(1) ! h
ee(2) = 1. - rsn(2) ! v

! -------- imaginary part of ice epsilon
call diel_ice (fghz, tlumi+tfreeze,eps_ice)
eir = real(eps_ice)
eii = aimag(eps_ice)

#if (defined DEBUG)
write(*, *) "snow.f90: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"                fghz  tsnow(c)      diele. c. of ice" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(a10, 4f10.3)')"fghz:tn:e", fghz, tlumi, eir, eii 
write(*, *)
#endif

 ! -------- impure ice -5 c (matzler)
  a=0.0026
  b=0.00023
  c=0.87
  eiis=a/fghz+b*(fghz**c)
 ! -------- pure ice -5 c (matzler)
  ap=6.e-4               
  bp=6.5e-5
  cp=1.07
  eiip=ap/fghz+bp*(fghz**cp)
 ! -------- difference between pure and impure ice at -5 c
  deltaeii = eiis - eiip 
 ! -------- effect of salinity (using of the difference at -5 c):
 ! -------- imaginary part of ice epsilon
  eii = eii + deltaeii*ppsal_snow/13. 

! -------- real part of dry snow epsilon
reds = 1.+ 1.58*pds/(1.-0.365*pds)
xn1 = sqrt(y0/eps_0)

! -------- imaginary part of dry snow epsilon (pvs model)
vi=pds/pice
xieds = 3.*vi*eii*(reds**2.)*(2.*reds+1.) /((eir+2.*reds)*(eir+2.*reds**2.))

! -------- dielectric constant of wet snow (matzler 1997)
if (xmv > 0.) then
    aa = 0.005
    ab = 0.4975
    ac = 0.4975
    epssw = 88.
    ff0= 9.

    f0a = ff0 *(1.+ aa*(epssw-eps_winf)/(reds+aa*(eps_winf-reds)) )
    f0b = ff0 *(1.+ ab*(epssw-eps_winf)/(reds+ab*(eps_winf-reds)) )
    f0c = ff0 *(1.+ ac*(epssw-eps_winf)/(reds+ac*(eps_winf-reds)) )

    epsinfa = (xmv/3.)*(eps_winf-reds)/(1.+aa*((eps_winf/reds)-1.))
    epsinfb = (xmv/3.)*(eps_winf-reds)/(1.+ab*((eps_winf/reds)-1.))
    epsinfc = (xmv/3.)*(eps_winf-reds)/(1.+ac*((eps_winf/reds)-1.))

    epssa = (xmv/3.) * (epssw-reds)/(1.+aa*((epssw/reds)-1.))
    epssb = (xmv/3.) * (epssw-reds)/(1.+ab*((epssw/reds)-1.))
    epssc = (xmv/3.) * (epssw-reds)/(1.+ac*((epssw/reds)-1.))

    epsa = epsinfa + (epssa-epsinfa)/(1.+cim*fghz/f0a)
    epsb = epsinfb + (epssb-epsinfb)/(1.+cim*fghz/f0b)
    epsc = epsinfc + (epssc-epsinfc)/(1.+cim*fghz/f0c)
       
    epsmarka = epsa + epsb + epsc + (reds-cim*xieds)
    rews = real(epsmarka)
    xiews = -1. * aimag(epsmarka)
else
    rews = reds
    xiews = xieds
endif
  
alfa = k*abs(aimag(sqrt(rews-cim*xiews)))
beta = k*real(sqrt(rews-cim*xiews))

k2x = k * sintheta
ppp = 2. * alfa * beta
qq = beta**2. - alfa**2. - (k**2.)*(sintheta)**2.
xjakaja = (1./(sqrt(2.)))*sqrt(sqrt(ppp**2.+qq**2.)+qq)
! -------- propagation angle in snow
o2 = atan(k2x/xjakaja)                             
!etkullumi=o2*180./pi;      % non useful

#if (defined DEBUG)
write(*, *) "snow.f90: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"                fghz  tsnow(c)    die.c.  dry snow     die.c. wet snow"
write(*, *)"---------------------------------------------------------------------------"
write(*, '(a10, 6f10.3)')"fghz:tn:s", fghz, tlumi, reds, xieds, rews, xiews
write(*, *)
#endif


! -------- wave imbedance in snow
 cn2 = sqrt((y0/eps_0)/(rews-cim*xiews))                  

! polarization h, v
do i=1,2 
    
  ! -------- fresnel reflection coefficients between snow and air
  if (i == 1) then
    rf = (cn2*costheta-xn1*cos(o2))/(cn2*costheta+xn1*cos(o2))
  else
    rf = (xn1*costheta-cn2*cos(o2))/(xn1*costheta+cn2*cos(o2))
  endif

  ! -------- power transmission coefficient between air and snow:
  tcoh = (1. - (abs(rf))**2.)

  tsa=tcoh
  trans=tsa

  ! -------- dry snow extinction coefficient (db/m => np/m)
  sec = 1./cos(o2)
  ! -------- hallikainen 1986, if fghz<60  ked=f(fghz,ppdiam)
  ked = 0.0018*(fghz**2.8)*(ppdiam**2.0)
  ked = ked/4.3429  ! -------- (np/m)

  ! -------- dry snow absorption coefficient (np/m)
  kabsd = 2.*omega*sqrt(y0*eps_0*reds) * sqrt(0.5*(sqrt(1.+(xieds/reds)**2.)-1.))
  ! YDT: should be this? 
  !kabsd = 2.* k*sqrt(y0*eps_0*reds) * sqrt(0.5*(sqrt(1.+(xieds/reds)**2.)-1.))
  ! ----- correction to extinction coefficient for small frequencies (f<18 ghz):
  if (ked < kabsd)  ked = kabsd

  ! -------- absorption coefficient of wet snow:
  kabs = 2.*omega* sqrt(y0*eps_0*rews) * sqrt(0.5*(sqrt(1.+(xiews/rews)**2.)-1.))
  ! YDT: should be this?
  !kabs = 2.* k * sqrt(y0*eps_0*rews) * sqrt(0.5*(sqrt(1.+(xiews/rews)**2.)-1.))

  ! -------- attenuation due to absorption in wet snow (np/m):
  xinta=0.001*kabs/pdsw*weq*sec 

  ! -------- total extinction:             ! (assuming that scattering 
  !        coeff. (kappas) is the same as in the case of dry snow!!!
   ke = (ked - kabsd) + kabs
  ! -------- total attenuation (np/m)
  xinte = 0.001*ke/pdsw*weq*sec           
  suhde=xinte/xinta

#if (defined DEBUG)
write(*, *) "snow.f90 outputs: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"         k        y0     eps_0      reds" 
write(*, '(4f10.3)') k, y0, eps_0, reds 
write(*, *)  &
"      pds      pdsw       xmv     dlumi      reds     tlumi       eir       eii" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(8f10.3)')pds, pdsw, xmv, dlumi, reds, tlumi, eir, eii
write(*, *)  &
"      ked     kabsd      kabs        ke      pdsw       weq       sec" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(7f10.3)') ked, kabsd, kabs, ke, pdsw, weq, sec
write(*, *) &
"    xinte     xinta" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(2f10.3)') xinte,  xinta
write(*, *)
#endif

  ! -------- total attenuation
!YDT  l2=exp(xinte)       
  ! -------- attenuation due to absorption
!YDT  la=exp(xinta)       
  
  ! -------- reflectivity of snow-air boundary:
  xr1 = 1. - trans

  ! -------- ground emission contribution
  ks = ke - kabs
  l2apu = exp((ke-ppq2mod*ks)*sec*dlumi)

#if (defined DEBUG)
write(*, *) "snow.f90 outputs2: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"     ee(i)       xr1        ke        ks     l2apu   ppq2mod     exp()" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(7f10.3)')ee(i), xr1, ke, ks, l2apu, ppq2mod,(ke-ppq2mod*ks)*sec*dlumi
write(*, *)
#endif

  xml = ee(i) *tsoil/l2apu * (1.-xr1)/(1.-(1.-rsn(i) )*(xr1)/l2apu**2.)
  ysuora(1)=xml
  emiss1 = ysuora(1)/tsoil

  ! -------- snow emission contribution
  sulapu = (1.-xr1) * (tlumi+tfreeze) * (kabs/(ke-ppq2mod*ks)) * (1. - 1./l2apu)
  sul = (1. + rsn(i)/l2apu) * sulapu/(1.-rsn(i)*(xr1)/l2apu**2.)
  ysuora(2)=sul 
  emiss2 = ysuora(2)/(tlumi+tfreeze)

  ! -------- brightness temperature of snow covered terrain:
  tb_tos(i) = ysuora(1)+ysuora(2)

  ! -------- emissivity of snow covered terrain
  esn(i) = emiss1 + emiss2
  
enddo  
 
end subroutine cmem_snow
