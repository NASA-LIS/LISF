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
! !ROUTINE: AGRMET_w
! \label{AGRMET_w}
!
! !REVISION HISTORY:
!
!    01 mar 88  initial version................................rice/sddc  
!    25 apr 89  improve coeff and testing................rice&moore/sddc 
!    07 sep 99  ported to ibm sp-2.  added intent attributes to
!               arguments..................................mr gayno/dnxm
!    03 aug 05; Sujay Kumar, Adopted in LIS
!     
!
! !INTERFACE:
function AGRMET_w( cldtyp, cldamt, level, coszen ) 

  implicit none
! !ARGUMENTS: 
  integer, intent(inout)      :: cldtyp ( 3 )  
  integer, intent(in)         :: level 
  real,    intent(in)         :: cldamt    ( 3 )  
  real,    intent(in)         :: coszen   
!
! !DESCRIPTION:
!
!     to calculate a weighting coefficient to be used in   
!     determining transmissivity and reflectivity values   
!     for partly cloudy conditions from clear and cloudy   
!     transmissivity and reflectivity values.  
!
!     \textbf{Method} \newline
!     - a weight is assigned based on the level of cloud types and 
!       amounts plus the effect of the solar zenith angle using  
!       shapiro's polynomials. \newline
!     - this weight is then applied to the cloud amount at
!       a particular level. \newline
!
!  The arguments and variables are: 
!  \begin{description}  
!   \item[c0]         allocatables for weighting polynomial.   
!   \item[c1]         allocatables for weighting polynomial.   
!   \item[c2]        allocatables for weighting polynomial.   
!   \item[c3]         allocatables for weighting polynomial.   
!   \item[c4]        allocatables for weighting polynomial.   
!   \item[c5]         allocatables for weighting polynomial.   
!   \item[ca]         dummy cld amt variable   
!   \item[cldamt]     fraction of low, mid, and high cloud amounts.
!   \item[cldtyp]     low, mid, and high cloud types.  
!   \item[coef]       shapiro's polynomial coefficients.   
!   \item[coszen]     solar zenith angle.  
!   \item[level]      allocatable for cloud level. 
!   \item[agrmet\_w]          final weighting coefficient passed back to trcalc.   
!   \item[w1]         intermediate weighting coefficient.  
!   \item[w2]         intermediate weighting coefficient.  
!  \end{description}
!EOP
  real                        :: c0   
  real                        :: c1   
  real                        :: c2   
  real                        :: c3   
  real                        :: c4   
  real                        :: c5   
  real                        :: ca   
  real                        :: coef      ( 6, 4 )   
  real                        :: AGRMET_w
  real                        :: w1   
  real                        :: w2   
  
  data  coef / 0.675 , -3.432 ,  1.929 , 0.842 ,  2.693 , -1.354 ,  &
       1.552 , -1.957 , -1.762 , 2.067 ,  0.448 ,  0.932 , & 
       1.429 , -1.207 , -2.008 , 0.853 ,  0.324 ,  1.582 , & 
       1.512 , -1.176 , -2.160 , 1.420 , -0.032 ,  1.422 /  
  
!     ------------------------------------------------------------------
!     execution starts here...reduce range of cldtyp by replacing
!     fives with fours.
!     ------------------------------------------------------------------

  if ( cldtyp(level) .eq. 5 )  cldtyp(level) = 4

!     ------------------------------------------------------------------
!     retrieve polynomial coefficients based on cld type at level.
!     ------------------------------------------------------------------

  c0 = coef ( 1, cldtyp(level) )
  c1 = coef ( 2, cldtyp(level) )
  c2 = coef ( 3, cldtyp(level) )
  c3 = coef ( 4, cldtyp(level) )
  c4 = coef ( 5, cldtyp(level) )
  c5 = coef ( 6, cldtyp(level) )

!     ------------------------------------------------------------------
!     compute the weight (w) based on the coefficients above, the   
!     cloud amount and the solar zenith angle.  
!     ------------------------------------------------------------------

  ca = cldamt(level)
  
  w1 = c0 + (c1 * coszen) + (c2 * ca)   
  w2 = ( (c4 * coszen) + (c3 * ca) ) * coszen + (c5 * ca * ca)  
  
  AGRMET_w  = w1 + w2   
  
  return
end function AGRMET_w
