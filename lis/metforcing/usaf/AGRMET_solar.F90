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
! !ROUTINE: AGRMET_solar
!  \label{AGRMET_solar}
!
! !REVISION HISTORY:
!
!     10 feb 88  initial version..........................capt rice/sddc
!     21 apr 89  testing, updating, error corrections....rice&moore/sddc   
!     07 sep 99  ported to ibm-sp2.  updated prolog. removed all
!                8th mesh "box" logic.  replaced call to grd2ll with
!                new utility pstoll........................mr gayno/dnxm
!     10 jun 02  removed all references to rtneph..........mr gayno/dnxm  
!     03 aug 05; Sujay Kumar, Adopted in LIS
!
! !INTERFACE:    
subroutine AGRMET_solar ( coszen, albdo, &
     r1,r2,r3,t1,t2,t3,rsolin,yr1,mo1,da1,hr1)   
! !USES:
  use LIS_constantsMod, only : LIS_CONST_PI, LIS_CONST_SOLAR

  implicit none
! !ARGUMENTS: 
  real                         :: coszen
  real,        intent(in)      :: albdo
  real,        intent(out)     :: rsolin   
  real                         :: r1   
  real                         :: r2   
  real                         :: r3   
  real                         :: t1   
  real                         :: t2   
  real                         :: t3   
  integer                      :: yr1,mo1,da1,hr1
!
! !DESCRIPTION:
!  
!     to compute the net downward shortwave ( solar )  
!     radiation value at the earth's surface.  
!  
!
!     \textbf{Method} \newline
!     
!     shapiro's method is based on a 3-layer plane-parallel 
!     atmosphere.  each layer ( high, middle, and low clouds )   
!     transmits and reflects some of the solar radiation incident   
!     on it from above and from below.  each layer's trans- 
!     missivity and reflectivity value depends on the layer's   
!     cloud type and amount and on the solar zenith angle.  
!     empirically derived reflectivity and transmissivity values
!     for each of the three layers and the fraction of the  
!     solar radiation entering the top of the troposphere (for  
!     the given location are used to calculate the surface  
!     insolation for the point. 
!
!  The arguments and variables are: 
!  \begin{description}
!  \item[alat]
!    latitude of point.  
!  \item[albdo]
!    surface albedo of the point  
!  \item[alon]
!    longitude of the point.  
!  \item[asolcn]
!   amount of solar radiation entering the top of
!   the atmosphere for a particular time of the year.
!  \item[cdel]
!   cosine of the sun's declination angle. 
!  \item[clat]
!   cosine of the sun's elevation angle.
!  \item[coszen]
!   solar zenith angle.
!  \item[d2]
!   double back scatter from air and clouds. 
!  \item[deltim]
!   longitudnal zulu time difference.
!  \item[fog]
!   present weather that indicates fog is present
!  \item[frac]
!   intermediate step to fracsq. 
!  \item[fracsq]
!   fraction of sunlight based on time of year.  
!  \item[hemi]
!   hemisphere (1=nh, 2=sh)  
!  \item[hrangl]
!   hour angle of the sun.   
!  \item[iclamt]
!   low, middle, and high cloud amounts. 
!  \item[icltyp]
!   low, middle, and high cloud types.
!  \item[julday]
!   julian day of the calendar year. 
!  \item[pi]
!   geometric pi = 3.14. 
!  \item[pid12]
!   coefficient for hour angle.  
!  \item[pid180]
!   coefficient for sun' angle.
!  \item[r1]
!   reflectivity coefficient for high clouds.
!  \item[r2]
!   reflectivity coefficient for mid clouds. 
!  \item[r3]
!   reflectivity coefficient for low clouds. 
!  \item[rsolin]
!   solar radiation (w m-2)  
!  \item[rtop]
!   shortwave radiation at the top of the atmosphere.
!  \item[sdec]
!   solar declination angle
!  \item[sdel]
!   solar zenith angle
!  \item[slat]
!   sine of sun's elevation angle.   
!  \item[solcon]
!    solar constant (W m-2)  
!  \item[t1]
!   transmissivity coefficient for high clouds.  
!  \item[t2]
!   transmissivity coefficient for mid clouds.   
!  \item[t3]
!   transmissivity coefficient for low clouds.   
!  \item[ztime]
!   current zulu time.   
!  \end{description}
!
!
!  The routines invoked are: 
!  \begin{description}
!  \item[agrmet\_bakfac](\ref{AGRMET_bakfac}) \newline
!    computes the backscatter factor 
!  \end{description}
!EOP
  integer                      :: julday
  real                         :: asolcn   
  real,        external        :: AGRMET_bakfac
  real                         :: d2
  real                         :: frac 
  real                         :: fracsq   
!  real                         :: pid12
!  real                         :: pid180   
  real                         :: rtop 
  integer                      :: days   ( 12 ) 


  data days   / 0,31,59,90,120,151,181,212,243,273,304,334 /
!  data pid12  / 0.2617993878 /  
!  data pid180 / 0.0174532925 /  

!     ------------------------------------------------------------------
!     compute a julian day using month and day fm array datime.
!     ------------------------------------------------------------------

  julday = days(mo1) + da1

!     ------------------------------------------------------------------
!     adjust the julian day for leap years.  
!     ------------------------------------------------------------------
!
  if ( (mo1 .gt. 2) .and. (mod(yr1,4) .eq. 0) ) then  
     julday = julday + 1
  endif
!     
!     ------------------------------------------------------------------
!     adjust the solar constant based on julian day.
!     ------------------------------------------------------------------
  
  frac   = (1.00014 + &
       (0.016726 * cos(2.0 * LIS_CONST_PI * (julday - 2) / 365.0)))
  fracsq = frac * frac  
  asolcn = LIS_CONST_SOLAR * fracsq  

  if ( coszen .ge. 0.01 ) then  
!     ------------------------------------------------------------------
!       compute the top-of-atmosphere insolation then retrieve all of   
!       the transmission and reflection coefficients for the given  
!       cld types and amounts.  
!     ------------------------------------------------------------------

     rtop = asolcn * coszen  

!     ------------------------------------------------------------------
!       determine the effect of backscattering on incoming solar
!       radiation.  
!     ------------------------------------------------------------------

     d2 = AGRMET_bakfac( t2, t3, r1, r2, r3, albdo ) 

!     ------------------------------------------------------------------
!       calc the dwnwrd solar radiation at the earth's surface based
!       on the insolation at the top of the atmosphere and upon the 
!       atmosphere's transmissivity, reflection, and scattering.
!     ------------------------------------------------------------------
     rsolin = ( (t1 * t2 * t3 * rtop) / d2 ) 
     
  else
     rsolin = 0.0
  endif

  return

end subroutine AGRMET_solar
