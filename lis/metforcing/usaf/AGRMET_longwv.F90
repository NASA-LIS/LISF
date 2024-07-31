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
! !ROUTINE: AGRMET_longwv
! \label{AGRMET_longwv}
!
! !REVISION HISTORY:
!     15 may 88  initial version..........................capt rice/sddc  
!     07 sep 99  ported to ibm sp-2.  added intent attributes to
!                arguments.  updated prolog................mr gayno/dnxm
!
!     29 jul 05  Sujay Kumar, Adopted in LIS
!
! !INTERFACE:    
subroutine AGRMET_longwv( sfctmp, e,iclamt, rldown )  

  implicit none

! !ARGUMENTS:   
  real,  intent(in)            :: iclamt  ( 3 ) 
  real,     intent(in)         :: e  
  real,     intent(out)        :: rldown   
  real,     intent(in)         :: sfctmp   

  
! !DESCRIPTION:
!     to compute the net downward longwave radiation at the
!     earth's surface.
!  
!     \textbf{Method} \newline
!     
!     - calculate the emissivity of the clear sky. \newline
!     - calculate downwelling longwave radiation of the clear sky. \newline
!     - add contribution of low, middle and high clouds to the 
!       clear sky portion. \newline
!
!     References: \newline
!    
!     dr idso's paper in the j. of geophys. research,   
!     no 74, pp 5397-5403.  \newline
!  
!     dr r.f.wachtmann's paper in the digest of preprints,  
!     topical meeting on remote sensing of the atmosphere,  
!     anaheim,ca, optical society of america, entitled, 
!     "expansion of atmospheric temperature-moisture
!     profiles in empirical orthogonal functions for remote 
!     sensing applications", 1975  \newline
!
!  The arguments and variables are: 
!  \begin{description}
!  \item[e]
!   sfc vapor pressure (pascals) 
!  \item[iclamt]
!   low, mid, and high cloud amounts in percent
!  \item[rldown]
!   net downward longwave irradiance.
!  \item[sfctmp]
!   sfc temperature ( deg k )
!  \item[cldfrt]
!   fraction of low, mid, and high cloud amounts
!  \item[clrsky]
!   irradiance contribution from the clear sky  
!  \item[emb]
!   sfc vapor pressure (millibars)   
!  \item[emissa]
!   idso clr sky emissivity (all wavelengths)
!  \item[emissb]
!   sasc clr sky emissivity (adjusted emissa value)
!  \item[hcterm]
!   high cloud emissivity term 
!  \item[lcterm]
!   low cloud emissivity term  
!  \item[mcterm]
!   mid cloud emissivity term  
!  \item[sigma]
!   stefan boltzman constant
!  \item[zh]
!   high cloud height coefficient   
!  \item[zl]
!   low cloud height coefficient
!  \item[zm]
!   mid cloud height coefficient
!  \end{description}
!EOP
  real                         :: cldfrt  ( 3 )
  real                         :: clrsky   

  real                         :: emb  
  real                         :: emissa   
  real                         :: emissb   
  real                         :: hcterm   
  real                         :: lcterm   
  real                         :: mcterm   
  real,     parameter          :: sigma = 5.67e-08
  real,     parameter          :: zh    = 8.0
  real,     parameter          :: zl    = 1.3
  real,     parameter          :: zm    = 3.1  

!     ------------------------------------------------------------------
!     executable code starts here...compute the cloud amount 
!     in fraction of overcast (.1 to 1.0). 
!     ------------------------------------------------------------------

  cldfrt(1) = iclamt(1)  / 100.0
  cldfrt(2) = iclamt(2)  / 100.0
  cldfrt(3) = iclamt(3)  / 100.0
  
!     ------------------------------------------------------------------
!     convert vapor pressure units from pascals to millibars for use
!     in determining emissivities.  
!     ------------------------------------------------------------------

  emb = e * 0.01

!     ------------------------------------------------------------------
!     compute the effective clr sky emissivity for all wavelengths  
!     (emissa) using idso's equation.   
!     ------------------------------------------------------------------

  emissa = 0.700 + (5.95e-5 * emb * exp(1500 / sfctmp))

!     ------------------------------------------------------------------
!     use emissa in wachtmann's model for sky irradiance to calc a  
!     resultant longwave downward radiation value.  first calc a sasc   
!     emmisivity (emissb), which is an adjusted idso emmisivity.
!     then use emissb to calculate the blackbody irradiance of the  
!     clear sky (the 1st term of wachtmann's equation).
!     ------------------------------------------------------------------

  emissb = -0.792 + (3.161 * emissa) - (1.573 * emissa * emissa) 
  clrsky =  emissb * sigma * ( sfctmp * sfctmp * sfctmp * sfctmp )

!     ------------------------------------------------------------------
!     now compute the irradiance contribution from the low, middle, 
!     and hi cloud layers (the 2nd thru 4th terms in wachtmann' eqn).
!     ------------------------------------------------------------------

  lcterm = (80.0 - (5.0 * zl)) * cldfrt(1) 
  mcterm = (80.0 - (5.0 * zm)) * (1.0 - cldfrt(1)) * cldfrt(2)  
  hcterm = (80.0 - (5.0 * zh)) * (1.0 - cldfrt(1)) * &
       (1.0 - cldfrt(2)) * cldfrt(3)   

!     ------------------------------------------------------------------
!     put it all together to get a resultant downwrd longwave irrad.
!     ------------------------------------------------------------------

  rldown = clrsky + hcterm + mcterm + lcterm
  
  return

end subroutine AGRMET_longwv
