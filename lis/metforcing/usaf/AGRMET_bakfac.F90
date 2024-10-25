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
! !ROUTINE: AGRMET_bakfac
!  \label{AGRMET_bakfac}
!
! !REVISION HISTORY:
!
!     15 feb 88  initial version...............................rice/sddc
!     07 sep 99  ported to ibm sp-2.  updated prolog.  added intent
!                attributes to arguments...................mr gayno/dnxm
!    03 aug 05; Sujay Kumar, Adopted in LIS
!     
!
! !INTERFACE:    
function AGRMET_bakfac( t2, t3, r1, r2, r3, albdo ) 

  implicit none
! !ARGUMENTS
  real,        intent(in)      :: albdo
  real,        intent(in)      :: r1   
  real,        intent(in)      :: r2   
  real,        intent(in)      :: r3   
  real,        intent(in)      :: t2   
  real,        intent(in)      :: t3   
  real                         :: AGRMET_bakfac   

! !DESCRIPTION:
!
!     to compute the double backscattering effect that the air and   
!     clouds have on solar radiation.
! 
!    \textbf{Method} \newline
!     
!    calculate the double backscatter effect (card2) based on the   
!    transmissivity and reflectivity coefficients of the sky and
!    the surface albedo.
!  
!    The arguments and variables are: 
!    \begin{description}
!     \item[albdo]
!      surface albedo.  
!     \item[agrmet\_bakfac]
!      overall backscattering effect.   
!     \item[capd1]
!      intermediate step for backscatter equation.  
!     \item[capd2]
!       the backscatter equation as described by shapiro.
!     \item[capdo]
!       intermediate step for backscatter equation.  
!     \item[d1]
!       intermediate step for backscatter equation as a  
!       function of reflectivity.
!     \item[d2]
!       intermediate step for backscatter equation as a  
!       function of reflectivity.
!     \item[d3]
!       intermediate step for backscatter equation as a  
!       function of reflectivity and albedo. 
!     \item[r1]
!       reflectivity coefficient of low clouds.  
!     \item[r2]
!       reflectivity coefficient of mid clouds.  
!     \item[r3]
!       reflectivity coefficient of high clouds. 
!     \item[t2]
!       transmissivity coefficient of mid clouds.
!     \item[t3]
!       transmissivity coefficient of high clouds.   
!    \end{description}
! 
!EOP
  real                         :: capd1
  real                         :: capd2
  real                         :: capdo
  real                         :: d1   
  real                         :: d2   
  real                         :: d3   

!     ------------------------------------------------------------------
!     execution starts here...
!     calculate the intermediate steps to shapiro's backscatter eqn.
!     the backscattering effect is a function of reflectivity (r1 to
!     r3), albedo, and transmissivity (t2 and t3).  
!     ------------------------------------------------------------------

  d1 = 1.0 - (r1 * r2)   
  d2 = 1.0 - (r2 * r3)   
  d3 = 1.0 - (r3 * albdo)
  
  capdo = (d1 * r2 * albdo * t3 * t3)
  
  capd1 = (d2 * d1) - (r1 * r3 * t2 * t2)

!     ------------------------------------------------------------------
!     combine above steps to calculate backscatter (capd2 and bakfac).  
!     ------------------------------------------------------------------

  capd2  = (d3 * capd1) - (capdo)
  AGRMET_bakfac = capd2 - (r1 * albdo * ((t2 * t3) * (t2 * t3)))   

  return
  
end function AGRMET_bakfac
