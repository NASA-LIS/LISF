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
! !ROUTINE: AGRMET_trnref
! \label{AGRMET_trnref}
!
! !REVISION HISTORY:
!
!    15 feb 88  initial version..........................rice,moore/sddc  
!    21 apr 89  testing, logic improvements..............rice&moore/sddc  
!    07 sep 99  ported to ibm sp-2.  updated prolog.  added intent
!               attributes to arguments....................mr gayno/dnxm
!    03 aug 05; Sujay Kumar, Adopted in LIS
!
! !INTERFACE:    
subroutine AGRMET_trnref( cldtyp, cldamt, clrtyp, coszen, &
     t1, t2, t3, r1, r2, r3 )  

  implicit none
! !ARGUMENTS:   
  integer,      intent(in)     :: cldtyp  ( 3 ) 
  integer,      intent(in)     :: clrtyp  ( 3 )   
  real,         intent(in)     :: cldamt     ( 3 ) 
  real,         intent(in)     :: coszen   
  real,         intent(out)    :: r1   
  real,         intent(out)    :: r2   
  real,         intent(out)    :: r3   
  real,         intent(out)    :: t1   
  real,         intent(out)    :: t2   
  real,         intent(out)    :: t3   

!
! !DESCRIPTION:
!     to compute all three (low, middle, and high cloud)   
!     transmissivities, and relflectivities.   
!  
!     \textbf{Method} \newline
!
!     - decide if radiation incident on the layer in question
!       is direct or diffuse (n/a for top layer). \newline
!     - call trcalc routine to calculate the required values. \newline
!
!  The arguments and variables are: 
!  \begin{description}
!  \item[cldamt]    fraction of low, mid, and high cloud amounts.
!  \item[cldtyp]    low, mid, and high cloud types.  
!  \item[clrtyp]    shapiro low, mid, and high cloud types.  
!  \item[coszen]    solar zenith angle.  
!  \item[difuse]    logical operator for diffuse radiation.  
!  \item[level]     identifies cloud level (low, mid, or high).  
!  \item[r1]        high level reflectivity coefficient. 
!  \item[r2]        mid level reflectivity coefficient.  
!  \item[r3]        low level reflectivity coefficient.  
!  \item[t1]        high level transmissivity coefficient.   
!  \item[t2]        mid level transmissivity coefficient.
!  \item[t3]        low level transmissivity coefficient.
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[agrmet\_trcalc](\ref{AGRMET_trcalc}) \newline
!   computes tranmissivity and reflectivities
!  \end{description}
!EOP
  integer                      :: level 
  logical                      :: difuse   

!     ------------------------------------------------------------------
!     execution starts here...upper level processing:   
!     ------------------------------------------------------------------

  level  = 1  
  difuse = .false.   

!     ------------------------------------------------------------------
!     call trcalc to calculate reflectivity and transmissivity  
!     coefficients for first level ( t1 and r1 ).
!     ------------------------------------------------------------------

  call AGRMET_trcalc( cldamt, cldtyp, clrtyp, coszen, &
       level, difuse, t1, r1)   

!     ------------------------------------------------------------------
!     middle level processing:  
!     ------------------------------------------------------------------

  level = 2  

!     ------------------------------------------------------------------
!     first, determine if the radiation incident on middle layer from   
!     above is diffuse (i.e., if there are clds above the layer)
!     set the 'difuse' flag accordingly.
!     ------------------------------------------------------------------

  if( (cldamt(1) .ge. 0.875) .and. (cldtyp(1) .eq. 2) ) &
       difuse = .true. 

!     ------------------------------------------------------------------
!     retrieve coefficients for level 2 ( t2 and r2 ) by calling trcalc.
!     ------------------------------------------------------------------

  call AGRMET_trcalc ( cldamt, cldtyp, clrtyp, coszen, &
       level, difuse, t2, r2 )   

!     ------------------------------------------------------------------
!     lower level processing:   
!     ------------------------------------------------------------------

  level = 3  

!     ------------------------------------------------------------------
!     determine if radiation incident on lowest layer from above is 
!     diffuse  (i.e., if hi or mid lvl clds are above this layer)   
!     set 'difuse' flag accordingly.
!     ------------------------------------------------------------------

  if (cldamt(2) .ge. 0.875)  difuse = .true. 

!     ------------------------------------------------------------------
!     retrieve coefficients for level 3 ( t3 and r3 ) by calling trcalc.
!     ------------------------------------------------------------------

  call AGRMET_trcalc( cldamt, cldtyp, clrtyp, coszen, level, &
       difuse, t3, r3 )   
  
  return

end subroutine AGRMET_trnref
