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
! !ROUTINE: agrmet_calc_albedo
!  \label{AGRMET_calc_albedo}
!
! !REVISION HISTORY:
!
!     20 oct 99  initial version based on ncep routines.................
!                ..........................lt col mitchell/mr gayno/dnxm
!     12 apr 02  pass in variables salp and snup instead of hard-wiring
!                them.  modified for new vegetation type database (USGS)
!                ..........................................mr gayno/dnxm
!     14 may 02  initialized albedo array to zero..........mr gayno/dnxm
!
!      8 aug 2005: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine agrmet_calc_albedo ( alb, albedo, shdfac, snoalb, &
  sneqv, snup, salp, bare)

  implicit none
! !ARGUMENTS:   
  real,           intent(in)     :: alb    
  real,           intent(out)    :: albedo 
  real,           intent(in)     :: salp   
  real,           intent(in)     :: shdfac 
  real,           intent(in)     :: sneqv  
  real,           intent(in)     :: snoalb 
  real,           intent(in)     :: snup   
  logical,        intent(in)     :: bare
!
! !DESCRIPTION:
!
!     calculate the snow-free albedo for the current day and
!     then modify for current snow cover.
!
!     \textbf{Method}
!
!     - determine the current julian day. \newline
!     - determine the seasonal snow-free albedo files that bound
!       the current day. \newline
!     - time interpolate snow-free albedo to the current day. \newline
!     - calculate albedo based on snow-free value and the current
!       snow cover. \newline
!     - snow-free and snow-adjusted albedo are passed back to
!       driver for use in ncep land-sfc model. \newline
!
!
!  The arguments and variables are: 
!  \begin{description}
!  \item[alb]    
!     snow-free albedo on the agrmet grid for the current
!     day
!  \item[albedo]
!    albedo modified for snow cover (if any) for the
!    agrmet grid for the current day
!   \item[rsnow]
!     ratio of sneqv to snup
!    \item[shdfac]
!     greenness fraction
!    \item[sneqv]
!     snow liquid equivalent
!    \item[snoalb]
!     maximum snow albedo
!    \item[snofac]
!     snow-depth factor
!    \item[snup]
!      threshold snow depth (in water equivalent m) that
!      implies 100% snow cover
!   \end{description}
!EOP
  real                           :: snofac
  real                           :: rsnow

!     ------------------------------------------------------------------
!     now that we have the snow-free albedo, modify it for snowcover.
!     note: this modified albedo is also calculated in the land
!     surface model as one of its first steps.  however, the radiation
!     schemes in flux3, which are run before the land-surface model,
!     need the modified albedo as well.  therefore, calculate it here.
!
!     the algorithm here was lifted from the land-sfc model,
!     routines sflx and redprm.  therefore, if this calculation 
!     ever changes in a new version of the model it MUST be changed
!     here.
!     ------------------------------------------------------------------

  albedo = 0.0
    
  if ( sneqv .gt. 0.001) then
     
     if (sneqv .lt. snup) then
        rsnow  = sneqv / snup
        snofac = 1. - ( exp(-salp*rsnow) - rsnow*exp(-salp))
     else
        snofac = 1.0
     end if
           
!     ------------------------------------------------------------------
!           to be consistent with mm5 (see routine prmveg) the 
!           greenness is zeroed out over vegetation types 
!           where you would expect no plants.  take this into account
!           in the albedo calculation.  
!     ------------------------------------------------------------------
           
     if(bare) then 
        albedo = alb + snofac * &
             (snoalb - alb)
     else
        albedo = alb + (1.0 - shdfac) * &
             snofac * (snoalb - alb)
     end if
     
!     ------------------------------------------------------------------
!         no snow, set albedo to snow-free value.
!     ------------------------------------------------------------------

  else
     
     albedo = alb
     
  end if
  
  return
end subroutine agrmet_calc_albedo
 
