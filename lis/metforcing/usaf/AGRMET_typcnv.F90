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
! !ROUTINE: AGRMET_typcnv
! \label{AGRMET_typcnv}
!
!  
! !REVISION HISTORY:
!
!     10 feb 88  initial version..........................capt rice/sddc
!     21 apr 89  testing, logic improvements.............rice&moore/sddc   
!     07 sep 99  ported to ibm sp-2. added intent attributes to
!                arguments.  updated prolog................mr gayno/dnxm
!     10 jun 02  change all references from rtneph to cdfs2.............
!                ..........................................mr gayno/dnxm
!     03 aug 05  Sujay Kumar, Adopted in LIS
!
! !INTERFACE:    
subroutine AGRMET_typcnv( icltyp, iclamt, fog, cldtyp, clrtyp, cldamt )
!EOP  
  implicit none

  integer,      intent(out)    :: cldtyp  ( 3 ) 
  integer,      intent(out)    :: clrtyp  ( 3 ) 
  integer                      :: hi
  integer                      :: hityp   ( 0:9 )  
  integer,      intent(in)     :: iclamt  ( 3 ) 
  integer,      intent(in)     :: icltyp  ( 3 ) 
  integer                      :: low   
  integer                      :: lowtyp  ( 0:4 )   
  integer                      :: mid   
  integer                      :: midtyp  ( 0:7 )   
  
  logical,      intent(in)     :: fog   
  
  real,         intent(out)    :: cldamt  ( 3 ) 
!  
! !DESCRIPTION:

!     to convert the cdfs2 cloud type and amts to shapiro 
!     types and amts.  
!   
!     \textbf{Method} \newline
!     -cld types are used as indices to shapiro cld types. \newline
!     -then cdfs2 cld amts are converted from whole percentage  
!      values to tenths values. \newline
!     -fog below any cloud layer is accounted for in
!      solar radiation calculation. \newline
! 
!  The arguments and variables are: 
!  \begin{description}
!   \item[cldamt]     fraction of low, mid, and high cloud amounts (shapiro)
!   \item[cldtyp]     low, mid, and high cloud types (shapiro)
!   \item[clrtyp]     low, mid, and high shapiro clear cloud types
!   \item[fog]        identifies fog below any layer
!   \item[hi]         high cloud allocatable
!   \item[hityp]      shapiro high cloud type numbers  
!   \item[iclamt]     cdfs2 low, mid, and high cloud amounts. 
!   \item[icltyp]     cdfs2 low, mid, and high cloud types.   
!   \item[low]        low cloud allocatable.   
!   \item[lowtyp]     shapiro low cloud type numbers   
!   \item[mid]        mid cloud allocatable
!   \item[midtyp]     shapiro mid cloud types numbers  
!  \end{description}
!EOP
  data hityp  / 0, 0, 0, 0, 0, 0, 0, 0, 2, 1 /   
  data lowtyp / 0, 5, 4, 4, 5 / 
  data midtyp / 0, 0, 0, 0, 0, 3, 4, 3 /
  
!     ------------------------------------------------------------------
!     execution starts here...convert whole percentage
!     cdfs2 cloud amounts to tenths.
!     ------------------------------------------------------------------
  
  cldamt(3) = float( iclamt(1) ) / 100.0
  cldamt(2) = float( iclamt(2) ) / 100.0
  cldamt(1) = float( iclamt(3) ) / 100.0
  
!     ------------------------------------------------------------------
!     convert low level cloud types from cdfs2 to shapiro.
!     ------------------------------------------------------------------

  low       = icltyp(1)   
  cldtyp(3) = lowtyp( low )  
  
!     ------------------------------------------------------------------
!     then mid-level clouds.  
!     ------------------------------------------------------------------

  mid       = icltyp(2)   
  cldtyp(2) = midtyp( mid )  

!     ------------------------------------------------------------------
!     then high-level clouds. 
!     ------------------------------------------------------------------

  hi        = icltyp(3)
  cldtyp(1) = hityp( hi )
  
  if ( ( hi .eq. 8 ) .or. ( hi .eq. 9 ) ) then  
     if ( cldamt(1) .lt. 0.875 )  cldtyp(1) = 1  
  endif
  
!     ------------------------------------------------------------------
!     if there is no cld type for a level, zero-out the level's cld amt. 
!     ------------------------------------------------------------------

  if (cldtyp(3) .eq. 0)  cldamt(3) = 0.0
  if (cldtyp(2) .eq. 0)  cldamt(2) = 0.0 
  if (cldtyp(1) .eq. 0)  cldamt(1) = 0.0 

!     ------------------------------------------------------------------
!     set each level's clear type.
!     ------------------------------------------------------------------

  clrtyp(1) = 1 
  clrtyp(2) = 2 
  clrtyp(3) = 3 

!     ------------------------------------------------------------------
!     set level 3 clear type to 4 if fog flag is .true. 
!     ------------------------------------------------------------------

  if ( fog ) clrtyp(3) = 4 
  
  return
end subroutine AGRMET_typcnv
      
