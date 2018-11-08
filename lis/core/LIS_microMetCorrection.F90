!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: LIS_microMetCorrection
! \label{LIS_microMetCorrection}
!
! !REVISION HISTORY:
!
!  21 Dec 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine LIS_microMetCorrection(nest)
! !USES:
  use LIS_logMod, only : LIS_logunit

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: nest
! !DESCRIPTION:
!  Corrects Temperature, Pressure, Humidity and Longwave Radiation
!  values for topography (based on Liston et al; 2005) 
!
!  WARNING: This routine is not currently implemented
!
! The arguments are:  
! \begin{description}
!  \item[nest]
!   index of the nest 
! \end{description}
!EOP
#if 0 
  real, parameter ::a=611.21,b=22.452, c=240.97, pi = 3.14
  real, parameter :: grav = 9.81, rdry = 287., Sstar=1367.0
  real, parameter :: sigma = 5.67E-8
  integer :: t
  real :: lapse(12),tdlapse(12)
  real :: force_tmp,force_hum,force_prs,elevdiff
  real :: force_pcp,force_swd, hcforce
  real :: pcforce,Td,tcforce,tbar,e
  real :: theta, omegas, thetad, wind, Ww,force_u, force_v
  real :: sigma_c, zenith, dec,tau,mu,cosi, psi_dir,psi_dif
  real :: lhour,czenith
  real :: force_lwd, ea, epsilona,epsilonb
  integer :: zone
  real :: elev700, T700, Td700,Rh700,es
  data lapse/0.0044,0.0059,0.0071,0.0078,0.0081,0.0082,0.0081,&
       0.0081,0.0077,0.0068,0.0055,0.0047/
  data tdlapse/0.00041,0.00042,0.00040,0.00039,0.00038,0.00036,&
       0.00033,0.00033,0.00036,0.00037,0.00040,0.00040/

  do t=1,lis%ngrid(nest)
     force_tmp = lisdom(nest)%grid(t)%forcing(1)
     force_hum = lisdom(nest)%grid(t)%forcing(2)
     force_v = lisdom(nest)%grid(t)%forcing(6)
     force_u = lisdom(nest)%grid(t)%forcing(5)
     force_swd = lisdom(nest)%grid(t)%forcing(3)
     force_prs = lisdom(nest)%grid(t)%forcing(7)
     force_pcp = lisdom(nest)%grid(t)%forcing(8)
     elevdiff = lisdom(nest)%grid(t)%elev-modelelev(t)
     tcforce = force_tmp-lapse(lis%mo)*elevdiff
     tbar=(force_tmp+tcforce)/2.

     pcforce =force_prs/(exp((grav*elevdiff)/(rdry*tbar)))
     e = force_hum*force_prs/(0.622+force_hum*0.378)
     Td = c*log((e/a))/(b-log((e/a)))
     Td = Td -tdlapse(lis%mo)*c*elevdiff/b
     e = a*exp((b*Td)/(c+Td))
     hcforce = 0.622*e/(pcforce - 0.378*e)
     theta = 3*pi/2-atan(force_v/force_u)
     omegas = grid(t)%slope*cos(theta-grid(t)%aspect)
     Ww = 1+0.5*omegas+0.5*grid(t)%curv
     wind = Ww*sqrt(force_u*force_u + force_v*force_v)
     thetad = -0.5*omegas*sin(2*(grid(t)%aspect-theta))
     force_u = -wind*sin(theta+thetad)
     force_v = -wind*cos(theta+thetad)
     force_pcp = force_pcp*((1+0.00035*elevdiff)/(1-0.00035*elevdiff))

! Longwave similar to agrmet implementation.
     ea = e * 0.01
     epsilona = 0.70 + (5.95e-5 * ea * exp(1500 / tcforce))
     epsilonb = -0.792 + (3.161 * epsilona) - (1.573 * epsilona * epsilona) 
     force_lwd = epsilonb*5.67E-8*(tcforce)**4

!shortwave     
     elev700 = modelelev(t)+(70000-force_prs)/(1.29*9.8)
!     elev700 = rdry*force_tmp*log(force_prs/70000)/grav
!     elev700 = 5537.0986-modelelev(t)
     elev700 = elev700-modelelev(t)
     T700 =  force_tmp-lapse(lis%mo)*elev700-273.16
     e = force_hum*force_prs/(0.622)
     Td = c*log((e/a))/(b-log((e/a)))
     Td700 = Td -tdlapse(lis%mo)*c*elev700/b
     e = a*exp(b*Td700/(c+Td700))
     es = a*exp(b*T700/(c+T700))
     Rh700 = 100* e/es
     sigma_c = 0.832*exp((Rh700-100)/41.6)
     if(sigma_c > 1) sigma_c = 1.0
     call localtime(lis%gmt,lisdom(nest)%grid(t)%lon,lhour,zone)
     call coszenith(lisdom(nest)%grid(t)%lon,lisdom(nest)%grid(t)%lat,lhour,zone,lis%doy,&
          czenith,dec,tau)
     zenith = acos(czenith)
     mu = asin(cos(dec)*sin(tau)/sin(zenith))
     cosi = cos(lisdom(nest)%grid(t)%slope)*cos(zenith)+sin(lisdom(nest)%grid(t)%slope)*sin(zenith)*&
          cos(mu-lisdom(nest)%grid(t)%aspect)
!     psi_dir = (1-sigma_c)*(0.6-0.2*cos(zenith))
!     psi_dif = sigma_c*(0.3-0.1*cos(zenith))
     psi_dir = 0.85-0.65*sigma_c
     psi_dif = 0.15+0.65*sigma_c
!     force_swd = Sstar*(psi_dir*cosi+psi_dif*cos(zenith))
     force_swd = force_swd*(psi_dir*cos(grid(t)%slope)+psi_dir*&
          sin(grid(t)%slope)*tan(zenith)*cos(mu-grid(t)%aspect)+&
          psi_dif)

     lisdom(nest)%grid(t)%forcing(1) = tcforce
     lisdom(nest)%grid(t)%forcing(2) = hcforce
     lisdom(nest)%grid(t)%forcing(3) = force_swd
     lisdom(nest)%grid(t)%forcing(4) = force_lwd
     lisdom(nest)%grid(t)%forcing(7) = pcforce
     lisdom(nest)%grid(t)%forcing(5) = force_u
     lisdom(nest)%grid(t)%forcing(6) = force_v
     lisdom(nest)%grid(t)%forcing(8) = force_pcp
  enddo
#endif
  write(LIS_logunit,*) 'This routine is not currently supported'

end subroutine LIS_microMetCorrection
