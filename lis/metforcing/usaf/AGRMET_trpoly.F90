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
! !ROUTINE: AGRMET_trpoly
! \label{AGRMET_trpoly}
!
! !REVISION HISTORY:
!
!     01 mar 88   initial version.........................capt rice/sddc 
!     25 apr 89   testing and logic improvements.......rice & moore/sddc
!     07 sep 99   ported to ibm sp-2.  updated prolog.  added intent
!                 attributes to arguments..................mr gayno/dnxm
!  
!    03 aug 05; Sujay Kumar, Adopted in LIS
!     
!
! !INTERFACE:    
subroutine AGRMET_trpoly( cldtyp, clrtyp, coszen, level, sky, trn, ref )

  implicit none
! !ARGUMENTS:
  integer,     intent(inout)   :: cldtyp     ( 3 ) 
  integer,     intent(in)      :: clrtyp     ( 3 ) 
  integer,     intent(in)      :: level
  logical,     intent(in)      :: sky   
  real,        intent(in)      :: coszen     
  real,        intent(out)     :: trn  
  real,        intent(out)     :: ref  
!
! !DESCRIPTION:
!     to compute the transmissivity and reflectivity coefficients
!     using shapiro's method. 
!  
!     \textbf{Method} \newline
!
!     - compute the cosine function for reflectivity and transmissivity,   
!       then the clear sky coefficient, first for transmissivity, then 
!       for reflectivity. \newline
!     - repeat these steps for cloudy skies. \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[a0]        reflectivity allocatable for cloudy sky polynomial
!   \item[a1]        reflectivity allocatable for cloudy sky polynomial  
!   \item[a2]        reflectivity allocatable for cloudy sky polynomial  
!   \item[a3]        reflectivity allocatable for cloudy sky polynomial  
!   \item[b0]        transmissivity allocatable for cloudy sky polynomial
!   \item[b1]        transmissivity allocatable for cloudy sky polynomial
!   \item[b2]        transmissivity allocatable for cloudy sky polynomial
!   \item[b3]        transmissivity allocatable for cloudy sky polynomial
!   \item[c0]        trans/reflect argument for function fac
!   \item[c1]        trans/reflect argument for function fac
!   \item[c2]        trans/reflect argument for function fac
!   \item[c3]        trans/reflect argument for function fac
!   \item[cldtyp]    low,mid, and high cloud types
!   \item[clrtyp]    low,mid, and high clear cloud types
!   \item[cosz]    cosine of the solar zenith angle
!   \item[coszen]    cosine of the solar zenith angle
!   \item[level]     allocatable for cloud level
!   \item[rc]        clear sky reflectivity polynomial values
!   \item[ref]       the reflectivity value sent back to trcalc
!   \item[rhoc]      cloudy sky reflectivity polynomial values   
!   \item[sky]       flag for clear or cloudy sky
!   \item[tauc]      cloudy sky transmissivity polynomial values
!   \item[tc]        clear sky transmissivity polynomial values  
!   \item[trn]       final transmissivity value sent back to trcalc
!  \end{description}
!EOP

  
  real                         :: a0   
  real                         :: a1   
  real                         :: a2   
  real                         :: a3   
  real                         :: b0   
  real                         :: b1   
  real                         :: b2   
  real                         :: b3   
  real                         :: c0
  real                         :: c1
  real                         :: c2
  real                         :: c3
  real                         :: cosz
  real                         :: fac
  real                         :: rc         ( 4, 4 )  
  real                         :: rhoc       ( 4, 4 )  
  real                         :: tauc       ( 4, 4 )  
  real                         :: tc         ( 4, 4 )  

!     ------------------------------------------------------------------
!     clear-sky polynomial coefficients for transmissivity: 
!     ------------------------------------------------------------------

  data  tc / 0.76977 , 0.49407 , -0.44647 ,  0.11558 ,  &
       0.69318 , 0.68227 , -0.64289 ,  0.17910 ,  &
       0.68679 , 0.71012 , -0.71463 ,  0.22339 , & 
       0.55336 , 0.61511 , -0.29816 , -0.06663 /  

!     ------------------------------------------------------------------
!     cloudy-sky polynomial coefficients for transmissivity:
!     ------------------------------------------------------------------

  data  tauc / 0.63547 , 0.35229 ,  0.08709 , -0.22902 , &
       0.43562 , 0.26094 ,  0.36428 , -0.38556 , &
       0.23865 , 0.20143 , -0.01183 , -0.07892 , &
       0.15785 , 0.32410 , -0.14458 ,  0.01457 /

!     ------------------------------------------------------------------
!     clear-sky polynomial coefficients for reflectivity:   
!     ------------------------------------------------------------------

  data  rc / 0.12395 , -0.34765 , 0.39478 , -0.14627 ,  &
       0.15325 , -0.39620 , 0.42095 , -0.14200 ,  &
       0.15946 , -0.42185 , 0.48800 , -0.18493 ,  &
       0.27436 , -0.43132 , 0.26920 , -0.00447 /  

!     ------------------------------------------------------------------
!     cloudy-sky polynomial coefficients for reflectivity:  
!     ------------------------------------------------------------------

  data  rhoc / 0.25674 , -0.18077 , -0.21961 , 0.25272 , &
       0.42111 , -0.04002 , -0.51833 , 0.40540 , &
       0.61394 , -0.01469 , -0.17400 , 0.14215 , &
       0.69143 , -0.14419 ,  0.05100 , 0.06682 /
  
!     ------------------------------------------------------------------
!     function fac, defined here, is for computing both 
!     transmissivity and reflectivity:  
!                                   2            3  
!       fac= c  + c (cos z) + c (cos  z) + c (cos  z)   
!             0    1           2            3   
!   
!     as seen in afgl tr 82-0039, pp. 38-39.  fac is computed using 
!     horner's method, in order to speed computation:   
!     ------------------------------------------------------------------

  fac(c0, c1, c2, c3, cosz) = &
       (((((c3 * cosz) + c2) * cosz) + c1) * cosz) + c0
  
!     ------------------------------------------------------------------
!     if sky condition flag for level is false (i.e., clear ) then  
!     ------------------------------------------------------------------

  if ( .not. sky ) then 

!     ------------------------------------------------------------------
!       clear sky coefficients will be computed.
!       first transmissivity,...
!     ------------------------------------------------------------------

     b0 = tc( 1, clrtyp(level) ) 
     b1 = tc( 2, clrtyp(level) ) 
     b2 = tc( 3, clrtyp(level) ) 
     b3 = tc( 4, clrtyp(level) ) 

     trn = fac ( b0, b1, b2, b3, coszen )

!     ------------------------------------------------------------------
!       then reflectivity   
!     ------------------------------------------------------------------

     a0 = rc( 1, clrtyp(level) ) 
     a1 = rc( 2, clrtyp(level) ) 
     a2 = rc( 3, clrtyp(level) ) 
     a3 = rc( 4, clrtyp(level) ) 
     
     ref = fac ( a0, a1, a2, a3, coszen )

!     ------------------------------------------------------------------
!     if sky condition flag for level is true (i.e., cloudy) then   
!     ------------------------------------------------------------------

  else  

!     ------------------------------------------------------------------
!       reduce range of cldtyp by chging all fives to fours.
!     ------------------------------------------------------------------

     if ( cldtyp(level) .eq. 5 ) cldtyp(level) = 4   

!     ------------------------------------------------------------------
!       then compute transmissivity values...   
!     ------------------------------------------------------------------

     b0 = tauc( 1, cldtyp(level) )   
     b1 = tauc( 2, cldtyp(level) )   
     b2 = tauc( 3, cldtyp(level) )   
     b3 = tauc( 4, cldtyp(level) )   

     trn = fac ( b0, b1, b2, b3, coszen )

!     ------------------------------------------------------------------
!       and reflectivity values.
!     ------------------------------------------------------------------

     a0 = rhoc ( 1, cldtyp(level) )  
     a1 = rhoc ( 2, cldtyp(level) )  
     a2 = rhoc ( 3, cldtyp(level) )  
     a3 = rhoc ( 4, cldtyp(level) )  
     
     ref = fac ( a0, a1, a2, a3, coszen )

  endif

  return

end subroutine AGRMET_trpoly
