!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! 
!BOP
! 
! !ROUTINE: AGRMET_getpwe
! \label{AGRMET_getpwe}
!
! !REVISION HISTORY: 
!     20 oct 89  initial version..................mr moore/sddc(agromet)
!     20 may 91  added check for octal 7 value indicating missing rpt.  
!                updtd prolog.....................mr moore/sddc(agromet)
!     15 jun 96  added new variables (arguements) newd and oldd to be   
!                passed thru to setprc routine. updated prolog and  
!                brought up to stds...............mr moore/sysm(agromet)
!     03 apr 97  added the variables wmonum and ztime to the called 
!                routine setprc so that they may be used by the 
!                routine................................ssgt miller/sysm
!     04 dec 97  revised routine to call rainbo for bogus values instead
!                of retrieving the values from database. brought up to
!                standards and updated prolog..capt andrus/dnxm(agromet)
!      7 oct 99  ported to ibm sp-2, updated prolog, incorporated
!                FORTRAN 90 features.................capt hidalgo/agrmet
!     21 feb 01  removed call to routine setprc (which placed estimates
!                on the grid) and moved its logic within this routine.
!                moved observation loop and data checks to this
!                routine because routines onhr1, offhr1 and offhr2
!                were replaced with new precip decoder.....mr gayno/dnxm
!      9 sep 10  Modified to process data from JMOBS.Chris Franks/16WS/WXE/SEMS
!
! !INTERFACE: 
subroutine AGRMET_getpwe (n, nsize, isize, network, plat_id, ilat, ilon,&
     month, pastwx, preswx, wmoblk, oldd, &
     prcpwe, imax, jmax) 
! !USES: 
  use LIS_coreMod,  only : LIS_domain
  use map_utils, only: latlon_to_ij

  implicit none

! !ARGUMENTS: 
  integer,       intent(in)    :: n 
  integer,       intent(in)    :: isize  
  character*10,  intent(in)    :: network(isize)
  character*10,  intent(in)    :: plat_id(isize)
  integer,       intent(in)    :: ilat(isize)
  integer,       intent(in)    :: ilon(isize)
  integer,       intent(in)    :: imax
  integer,       intent(in)    :: jmax
  integer,       intent(in)    :: month
  integer,       intent(in)    :: nsize
  integer,       intent(in)    :: pastwx(isize)
  integer,       intent(in)    :: preswx(isize)
  integer,       intent(in)    :: wmoblk(isize)
  real,          intent(inout) :: oldd(imax,jmax)
  real,          intent(inout) :: prcpwe(imax,jmax)

! 
! !DESCRIPTION: 
!     to determine precipitation amounts using past and
!     present weather codes.
! 
!     \textbf{Method}
!
!     - loop through all observations for this time period. \newline
!     - check for valid lat/lon, wmo station id, present/past
!       weather code. \newline
!       - if observation is valid, call routine makpwe to 
!         calculate an estimate. \newline
!       - place this estimate at the nearest grid point.  if
!         there is more than one estimate within one grid box,
!         use the nearest one. \newline
!
!  The arguments and variables are: 
!  \begin{description}
!  \item[a,b,c,d]   holds i/j coordinates of estimate/grid point in function
!                  sumsqr
!  \item[network]   JMOBS network
!  \item[plat\_id]   JMOBS platform id
!  \item[hemi]      hemisphere flag (1 - northern, 2 - southern)
!  \item[ilat]      observation latitude
!  \item[ilon]      observation longitude
!  \item[imax]      grid dimension, i-direction
!  \item[irecord]   loop index
!  \item[isize]     maximum size of observation arrays
!  \item[jmax]      grid dimension, j-direction
!  \item[month]     current month
!  \item[newd]      distance (in grid lengths) between estimate and
!                    the nearest grid point 
!  \item[nsize]     number of observations in the arrays
!  \item[oldd]      holds distances between estimates and the
!                   nearest grid point
!  \item[pastwx]    past weather code
!  \item[preswx]    present weather code
!  \item[ri/rj]     i/j coordinate of an estimate in grid point space
!  \item[rlat/rlon] lat/lon of an estimate
!  \item[sumsqr]    function to calculate the distance between an estimate
!                   and the nearest grid point        
!  \item[wmoblk]    First two digits of WMO station id, created for non-WMO obs
!  \item[wxest]     precipitation estimate passed back from makpwe
!   \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[AGRMET\_makpwe](\ref{AGRMET_makpwe}) \newline
!   create a present weather estimate 
!  \item[lltops](\ref{lltops}) \newline
!   converts a lat lon values to a corresponding point on the AGRMET
!   grid
!   \end{description}
!EOP

  integer                      :: irecord
  integer                      :: wxest
  real                         :: a
  real                         :: b
  real                         :: c
  real                         :: d
  real                         :: newd
  real                         :: ri
  real                         :: rj
  real                         :: rlat
  real                         :: rlon
  real                         :: sumsqr

! 
!     ------------------------------------------------------------------
!     internal function used to calculate the distance between an
!     observation and a grid point.
!     ------------------------------------------------------------------
  
  sumsqr(a,b,c,d) = ((a-b)**2) + ((c-d)**2)

!     ------------------------------------------------------------------
!     loop through the observations for this hour.  if an observation
!     has a valid lat/lon/block station number and pres/past weather
!     report, then call makpwe to calculate an estimate.
!     ------------------------------------------------------------------
  
  PRES_PAST_WX : do irecord = 1, nsize
     
     if ( (ilat(irecord) > -99999998) .and.  &
          (ilon(irecord) > -99999998) ) then
        rlat = float(ilat(irecord)) * 0.01
        rlon = float(ilon(irecord)) * 0.01
     else
        cycle PRES_PAST_WX
     end if
     
     if ( preswx(irecord) > -99999998 .and. &
          pastwx(irecord) > -99999998 ) then
        
        call AGRMET_makpwe (preswx(irecord), pastwx(irecord), &
             rlat, month, wmoblk(irecord), wxest)
        
!     ------------------------------------------------------------------
!         call lltops to convert from lat/lon space to gridpoint
!         space.  note: dumhemi is passed in because routine lltops
!         can modify the hemisphere (for observations on the equator).
!         this can cause an infinite loop in the program as variable
!         hemi is a loop index in the driver routine.
!     ------------------------------------------------------------------

        call latlon_to_ij(LIS_domain(n)%lisproj, rlat, rlon, ri, rj)
        
        newd = sumsqr(ri, float(nint(ri)), rj, float(nint(rj))) &
             * 100.0
        
!     ------------------------------------------------------------------
!         if there is more than one estimate within one grid box,
!         pick the closest one.
!     ------------------------------------------------------------------
        if(nint(ri).ge.1.and.nint(ri).le.imax.and.&
             nint(rj).ge.1.and.nint(rj).le.jmax) then 
           if ( newd < oldd(nint(ri),nint(rj)) ) then
              
              prcpwe(nint(ri),nint(rj)) = float(wxest)
              
              oldd(nint(ri),nint(rj)) = newd
              
           end if
        endif
     end if

  enddo PRES_PAST_WX
  return
  
end subroutine AGRMET_getpwe
