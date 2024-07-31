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
! !ROUTINE: AGRMET_tr_coeffs
!  \label{AGRMET_tr_coeffs}
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
!     31 MAR 2010 added handling of 16th native grid in pstoll call.....
!                 ...................................Michael Shaw/WXE
!
! !INTERFACE:    
subroutine AGRMET_tr_coeffs( i,j,hemi, fog, icltyp, iclamt,   &
     r1,r2,r3,t1,t2,t3,coszen,yr1,mo1,da1,hr1,n)   
! !USES:
  use AGRMET_forcingMod,only : agrmet_struc
  implicit none
! !ARGUMENTS: 

  integer,     intent(in)      :: iclamt ( 3 )  
  integer,     intent(in)      :: icltyp ( 3 )  
  logical,     intent(in)      :: fog   
  integer,     intent(in)      :: hemi
  integer,     intent(in)      :: i
  integer,     intent(in)      :: j
  integer,     intent(in)      :: n 
  real                         :: r1   
  real                         :: r2   
  real                         :: r3   
  real                         :: t1   
  real                         :: t2   
  real                         :: t3   
  real                         :: coszen
  integer                      :: yr1,mo1,da1,hr1
!
! !DESCRIPTION:
!  
!     to compute the transmissivity and reflectivity coefficients
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
!  \item[albdo]
!    surface albedo of the point  
!  \item[cldamt]
!   fraction of low, middle, and high cloud amounts.
!  \item[cldtyp]
!   low, middle, and high cloud types. 
!  \item[clrtyp]
!   shapiro low, middle, and high clear cloud types. 
!  \item[coszen]
!   solar zenith angle.
!  \item[d2]
!   double back scatter from air and clouds. 
!  \item[fog]
!   present weather that indicates fog is present
!  \item[hemi]
!   hemisphere (1=nh, 2=sh)  
!  \item[iclamt]
!   low, middle, and high cloud amounts. 
!  \item[icltyp]
!   low, middle, and high cloud types.
!  \item[pi]
!   geometric pi = 3.14. 
!  \item[r1]
!   reflectivity coefficient for high clouds.
!  \item[r2]
!   reflectivity coefficient for mid clouds. 
!  \item[r3]
!   reflectivity coefficient for low clouds. 
!  \item[rtop]
!   shortwave radiation at the top of the atmosphere.
!  \item[solcon]
!    solar constant (W m-2)  
!  \item[t1]
!   transmissivity coefficient for high clouds.  
!  \item[t2]
!   transmissivity coefficient for mid clouds.   
!  \item[t3]
!   transmissivity coefficient for low clouds.   
!  \end{description}
!
!
!  The routines invoked are: 
!  \begin{description}
!  \item[agrmet\_typcnv](\ref{AGRMET_typcnv}) \newline
!    cloud type convertion routine
!  \item[agrmet\_trnref](\ref{AGRMET_trnref}) \newline
!    computes transmissivities
!  \item[agrmet\_bakfac](\ref{AGRMET_bakfac}) \newline
!    computes the backscatter factor 
!  \end{description}
!EOP
  integer                      :: cldtyp ( 3 )  
  integer                      :: clrtyp ( 3 )  
  real,        external        :: AGRMET_bakfac
  real                         :: cldamt    ( 3 )  


!     ------------------------------------------------------------------
!     execution starts here...get 'shapiro' cloud types and
!     cloud amounts.  
!     ------------------------------------------------------------------
  call AGRMET_typcnv ( icltyp, iclamt, fog, cldtyp, clrtyp, cldamt)  

  call AGRMET_coszen(yr1, mo1, da1, hr1, hemi, i, j, &
                     int(agrmet_struc(n)%imax/64), coszen)

!     ------------------------------------------------------------------
!     if the sun is above the horizon, it is daytime and there will 
!     be some incoming solar radiation.  thus, continue processing...   
!     ------------------------------------------------------------------

  if ( coszen .ge. 0.01 ) then  

     call AGRMET_trnref ( cldtyp, cldamt, clrtyp, coszen, &
          t1, t2, t3, r1, r2, r3 ) 
  else
     r1 = 0.0 
     r2 = 0.0 
     r3 = 0.0

     t1 = 0.0 
     t2 = 0.0 
     t3 = 0.0
  endif
  
  return

end subroutine AGRMET_tr_coeffs
