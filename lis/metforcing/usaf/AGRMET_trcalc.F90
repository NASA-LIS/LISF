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
! !ROUTINE: AGRMET_trcalc
! \label{AGRMET_trcalc}
!
! !REVISION HISTORY:
!
!     01 feb 88  initial version.........................rice,moore/sddc
!     21 apr 89  testing and logic improvements..........rice&moore/sddc
!     07 sep 99  ported to ibm sp-2.  updated prolog.  added intent
!                attributes to arguments...................mr gayno/dnxm
!    03 aug 05; Sujay Kumar, Adopted in LIS
!     
!
! !INTERFACE:    
subroutine AGRMET_trcalc( cldamt, cldtyp, clrtyp, coszen, level, difuse,  &
     trans, reflec)  

  implicit none
  
  integer,   intent(in)        :: cldtyp  ( 3 )
  integer,   intent(in)        :: clrtyp  ( 3 )
  integer,   intent(in)        :: level
  logical,   intent(in)        :: difuse   
  real,      intent(in)        :: cldamt  ( 3 )   
  real,      intent(in)        :: coszen   
  real,      intent(out)       :: reflec   
  real,      intent(out)       :: trans

!
! !DESCRIPTION:
!
!     to calculate the transmissivity and reflectivity 
!     coefficients for low, mid, and high clouds at the level in   
!     question with a direct or difuse insolation fm above.
!   
!     \textbf{Method} \newline
!
!     - compute clear sky diffuse transmissivity and reflectivity   
!       coefficients. \newline
!     - if overcast, computes cloudy layer coefficients.   \newline
!     - if partly cloudy, compute the weighting factor for
!       transmissivity and reflectivity. weight = 1.0 for diffuse sky.
!       the non-diffuse sky coefficients are computed in trpoly.  \newline
!  
!  The arguments and variables are: 
!  \begin{description}
!   \item[cldamt]    fraction of low, mid, and high cloud amounts.
!   \item[cldtyp]    low, mid, and high cloud types.  
!   \item[clrtyp]    shapiro low, mid, and high clear cloud types.
!   \item[coszen]    cosize of solar zenith angle.  
!   \item[difuse]    logical identifier for diffuse sky.  
!   \item[level]     allocatable for low, mid, and high cloud level.  
!   \item[phi]       equal to cloud amount times weighting coefficient.   
!   \item[rcldy]     reflectivity coefficients for cloudy sky.
!   \item[rclr]      reflectivity coefficients for clear sky. 
!   \item[rdfs]      reflectivity coefficients for clear diffuse sky. 
!   \item[reflec]    final reflectivity coefficient to trnref.
!   \item[rhodfs]    reflectivity coefficients for diffuse cloudy sky.
!   \item[taudfs]    transmissivity coefficients for diffuse cloudy sky   
!   \item[tcldy]     transmissivity coefficients for cloudy sky.  
!   \item[tclr]      transmissivity coefficients for clear sky.   
!   \item[tdfs]      transmissivity coefficients for diffuse clear sky.   
!   \item[trans]     final transmissivity coefficient to trnref.  
!   \item[weight]    weight for non-diffuse sky.  
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[agrmet\_trpoly](\ref{AGRMET_trcalc}) \newline
!   transmissivity and reflectivity polynomial routine  
!  \item[agrmet\_w](\ref{AGRMET_w}) \newline
!   computes weighting coefficient 
!  \end{description}
!EOP
  real                         :: phi  
  real                         :: rcldy
  real                         :: rclr 
  real                         :: rdfs    ( 4 )   
  real                         :: rhodfs  ( 5 )   
  real                         :: taudfs  ( 5 )   
  real                         :: tcldy
  real                         :: tclr 
  real                         :: tdfs    ( 4 )   

  real,      external          :: AGRMET_w
  real                         :: weight   

  data rdfs   / 0.0, 0.040 , 0.045 , 0.116 / 
  data rhodfs / 0.0, 0.0, 0.560 , 0.609 , 0.520 / 
  data taudfs / 0.0, 0.0, 0.361 , 0.311 , 0.400 / 
  data tdfs   / 0.0, 0.905 , 0.900 , 0.788 / 

!     ------------------------------------------------------------------
!     executable code starts here...if there are no clds at this level,   
!     ------------------------------------------------------------------
  
  if ( cldamt(level) .eq. 0. ) then 
     
!     ------------------------------------------------------------------
!       and if there is diffuse insolation from above,  
!     ------------------------------------------------------------------

     if ( difuse ) then  

!     ------------------------------------------------------------------
!         then use the data'd constant values for t and r   
!     ------------------------------------------------------------------

        trans  = tdfs ( clrtyp ( level ) ) 
        reflec = rdfs ( clrtyp ( level ) )

!     ------------------------------------------------------------------
!       otherwise if the insolation from above is direct,   
!     ------------------------------------------------------------------

     else

!     ------------------------------------------------------------------
!         call trpoly for the clr case to obtain the output t and r.
!     ------------------------------------------------------------------

        call AGRMET_trpoly( cldtyp, clrtyp, coszen, level, .false.,&
             trans, reflec) 

     endif

!     ------------------------------------------------------------------
!     if there are clouds at this level,
!     ------------------------------------------------------------------

  else  

!     ------------------------------------------------------------------
!       call trpoly for the cldy case t and r.
!     ------------------------------------------------------------------

     call AGRMET_trpoly( cldtyp, clrtyp, coszen, level, .true.,&
          tcldy, rcldy ) 

!     ------------------------------------------------------------------
!       if the level is not completely overcast...  
!     ------------------------------------------------------------------

     if ( cldamt(level) .ne. 1 ) then

!     ------------------------------------------------------------------
!         call trpoly for the clr case t and r.
!     ------------------------------------------------------------------

        call AGRMET_trpoly( cldtyp, clrtyp, coszen, level, .false.,&
             tclr, rclr)

!     ------------------------------------------------------------------
!         and if there is diffuse insolation from above.
!     ------------------------------------------------------------------

        if ( difuse ) then

!     ------------------------------------------------------------------
!           set 'weight' to 1.  
!     ------------------------------------------------------------------

           weight = 1.0

!     ------------------------------------------------------------------
!         otherwise, if the insolation from above is direct.
!     ------------------------------------------------------------------

        else  

!     ------------------------------------------------------------------
!           use the function 'w' to obtain 'weight.' 
!     ------------------------------------------------------------------

           weight = AGRMET_w( cldtyp, cldamt, level, coszen )   

        endif

!     ------------------------------------------------------------------
!         use 'weight' to determine partly cloudy t and r.
!     ------------------------------------------------------------------

        phi    =  cldamt(level) * weight  
        trans  = ( phi * tcldy ) + ( 1 - phi ) * tclr  
        reflec = ( phi * rcldy ) + ( 1 - phi ) * rclr 

!     ------------------------------------------------------------------
!       if this level is completely overcast... 
!     ------------------------------------------------------------------

     else

!     ------------------------------------------------------------------
!         and if there is diffuse insolation from above,
!     ------------------------------------------------------------------

        if ( difuse ) then

!     ------------------------------------------------------------------
!           then use the data'd constant values for t and r. 
!     ------------------------------------------------------------------

           trans  = taudfs ( cldtyp(level) )
           reflec = rhodfs ( cldtyp(level) )   

!     ------------------------------------------------------------------
!         otherwise, if the insolation from above is direct.
!     ------------------------------------------------------------------

        else  

!     ------------------------------------------------------------------
!           use the cldy case t and r as the output t and r.
!     ------------------------------------------------------------------

           trans  = tcldy   
           reflec = rcldy  

        endif
        
     endif
     
  endif
  
  return

end subroutine AGRMET_trcalc
