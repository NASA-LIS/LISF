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
! !ROUTINE: AGRMET_loadcloud
!  \label{AGRMET_loadcloud}
!
! !REVISION HISTORY:
!     11 feb 91  initial version..................mr moore/sddc(agromet)   
!     05 may 91  changed estwt arguments to thres, made thres an array   
!                deleted the assignment statements using estwt and the   
!                nint function, delted the estwt array, updated the  
!                prolog.......................capt bertone/sddc(agromet)   
!     01 jun 92  chgd range of scalar latbnd to match new size of array  
!                pcoef elsewhere in pgm. updtd prolog...................
!                .................................mr moore/sysm(agromet)   
!     02 sep 99  ported to ibm sp-2.  removed "box" logic.  updated
!                prolog.  added intent attributes to arguments..........
!                ..........................................mr gayno/dnxm
!     10 jun 02  modified to use cdfs2 data................mr gayno/dnxm
!
!     29Jul2005; Sujay Kumar, Adopted in LIS
!     31 MAR 2010 Add handling of 16th native grid..see code comments
!                 for where..............................Michael Shaw/WXE
!
! !INTERFACE:    
subroutine AGRMET_loadcloud(hemi, land,thres, times, amounts, tops, types, &
     cldtyp, cldamt, fog, julhr, imax,jmax)

  implicit none

! !ARGUMENTS: 
  integer, intent(in)          :: imax
  integer, intent(in)          :: jmax
  byte, intent(in)             :: amounts   ( 4, 1024, 1024)
  byte, intent(in)             :: types     ( 4, 1024, 1024)
  integer, intent(in)          :: land      (imax,jmax)
  integer, intent(out)         :: cldamt    ( 3, imax,jmax)
  integer, intent(out)         :: cldtyp    ( 3, imax,jmax)
  logical, intent(out)         :: fog       (imax,jmax)
  integer, intent(in)          :: hemi  
  integer, intent(in)          :: julhr
  integer, intent(in)          :: thres     ( 5 ) 
  integer*4, intent(in)        :: times     (1024, 1024)
  integer*2, intent(in)        :: tops      (4, 1024, 1024)

!  
! 
! !DESCRIPTION:
!
!     to convert cdfs2 data for use by the flux3 radiation routines.
!  
!     \textbf{Method} \newline
!
!     - check if the point is land. \newline
!     - loop over all four cdfs2 layers.
!       - check time of cdfs2 data, if within acceptable time range
!         then \newline
!         - convert from cdfs2 cloud type to flux3 cloud type. \newline
!           - if cdfs2 says type is cb, see if cloud top is
!             below user-specified threshold.  if it is, set
!             flux3 cloud type to cumulus instead. \newline
!         - set flux3 cloud amount for its three layers to the
!           maximum amount of low/mid/high clouds from cdfs2. \newline
!       - if data is old, set flux3 cloud arrays to zero (clear skies). \newline
!
!  The arguments and variables are: 
!  \begin{description}  
!   \item[amounts]   cdfs2 cloud amounts
!   \item[cldamt]    output cloud amounts for use in radiation calcs
!   \item[cldtim]    cdfs2 pixel times
!   \item[cldtyp]    output cloud types for use in radiation calcs
!   \item[fog]       fog present flag (logical)  
!   \item[hemi]      current hemisphere (n = 1, s = 2)   
!   \item[i]         loop counter
!   \item[ii]        corresponding agrmet i-coordinate on the cdfs2 grid
!   \item[ixx]       agrmet grid dimension, i-direction
!   \item[j]         loop counter
!   \item[jj]        corresponding agrmet j-coordinate on the cdfs2 grid
!   \item[julhr]     current julian hour
!   \item[jyy]       agrmet grid dimension, j-direction
!   \item[land]      point processing switches (land-sea mask)
!   \item[latbnd]    latitude band corresponding to cb thresholds
!   \item[lonslc]    longitude slice   
!   \item[lvl]       loop counter (cdfs2 cloud level)  
!   \item[nefamt]    cdfs2 cloud amount
!   \item[neftop]    cdfs2 cloud top height 
!   \item[neftyp]    cdfs2 cloud type
!   \item[thres]     cb cloud top thresholds from flux3 control file
!   \item[times]     cdfs2 pixel times
!   \item[tops]      cdfs2 cloud tops
!   \item[types]     cdfs2 cloud types
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[AGRMET\_bndslc](\ref{AGRMET_bndslc}) \newline
!   determine latitude band and longitude slice
!  \end{description}
!EOP
  integer                      :: i,j
  integer                      :: ii
  integer                      :: jj
  integer                      :: cldtim
  integer                      :: latbnd
  integer                      :: lonslc
  integer                      :: lvl   
  integer                      :: nefamt
  integer                      :: neftop
  integer                      :: neftyp


!     ------------------------------------------------------------------
!     executable code begins here...initialize the output arrays.
!     note: the rtneph had the capability to analyze for "fog."  
!     however, cdfs2 does not analyze for "fog."  therefore,
!     initialize the array to false.  keep the "fog" logic throughout
!     AGRMET for now.  someday we may want to use it.
!     ------------------------------------------------------------------

  cldtyp = 0
  cldamt = 0
  fog    = .false.

!     ------------------------------------------------------------------
!     loop thru the four cdfs2 cloud layers.
!     ------------------------------------------------------------------

  do lvl = 1, 4  
     do j = 1, jmax
        do i = 1, imax

!     ------------------------------------------------------------------
!           only process over land areas.
!     ------------------------------------------------------------------
           if ( land(i,j) .gt. 0) then

!     ------------------------------------------------------------------
!             Michael Shaw - extract cloud type, amount, and top ht from database.
!             cloud data is 16th mesh, so every other grid point
!             corresponds to the agrmet grid if the agrmet native grid
!             is 8th mesh.
!     ------------------------------------------------------------------
!             (i*2)-1 = 1,3,5,...,1023 (i = 1,2,3,...512)
!             i       = 1,2,3,...,1024 (i = 1,2,3,...1024)
!     ------------------------------------------------------------------

              if(imax.eq.512)then
                ii = (i*2) - 1
                jj = (j*2) - 1
              elseif(imax.eq.1024)then
                ii = i 
                jj = j
              endif 
              neftyp = types   (lvl,ii,jj)
              nefamt = amounts (lvl,ii,jj)
              neftop = tops    (lvl,ii,jj)
              cldtim = times   (ii,jj)              
!     ------------------------------------------------------------------
!             if cdfs2 cloud type is of known type, proceed... 
!     ------------------------------------------------------------------

              if ( (neftyp.gt.0) .and. (neftyp.le.9) ) then

!     ------------------------------------------------------------------
!               retrieve this pt's lat band and lon slice.
!     ------------------------------------------------------------------

                 call AGRMET_bndslc(i,j,hemi, latbnd, lonslc )  

!     ------------------------------------------------------------------
!               if this pt's cld type is cirrus or cirrostratus...
!     ------------------------------------------------------------------

                 if ( (neftyp .eq. 9) .or. (neftyp .eq. 8) ) then

!     ------------------------------------------------------------------
!                 and it's valid time is newer than 6 hours old...
!     ------------------------------------------------------------------

                    if ( ((julhr-cldtim) .lt.  6) .and.&
                         ((julhr-cldtim) .gt. -3) )  then

!     ------------------------------------------------------------------
!                   save it's cld type and amount in the '3' arrays.
!     ------------------------------------------------------------------

                       cldtyp(3,i,j) = neftyp
                       cldamt(3,i,j) = max( cldamt(3,i,j), nefamt )
                       
                    endif

!     ------------------------------------------------------------------
!               if the pt's cloud type is a mid layer type...   
!     ------------------------------------------------------------------

                 elseif ( neftyp .ge. 5 ) then

!     ------------------------------------------------------------------
!                 and it's valid time is newer than 6 hours old...
!     ------------------------------------------------------------------
                    if ( ((julhr-cldtim) .lt.  6) .and.&
                         ((julhr-cldtim) .gt. -3) )  then

!     ------------------------------------------------------------------
!                   save it's cld type and amount in the '2' arrays 
!     ------------------------------------------------------------------

                       cldtyp(2,i,j) = neftyp
                       cldamt(2,i,j) = max( cldamt(2,i,j), nefamt )
                       
                    endif

!     ------------------------------------------------------------------
!               if the pt's cld type is a low layer type... 
!     ------------------------------------------------------------------

                 elseif ( neftyp .ge. 1 ) then 
                    
!     ------------------------------------------------------------------
!                 and it's valid time is newer than 3 hours old.
!     ------------------------------------------------------------------

                    if ( ((julhr-cldtim) .lt.  3) .and. &
                         ((julhr-cldtim) .gt. -3) )  then

!     ------------------------------------------------------------------
!                   and it's cld type is not cumulonimbus...
!     ------------------------------------------------------------------

                       if ( neftyp .ne. 1) then 

!     ------------------------------------------------------------------
!                     save it's cld type and amount in the '1' arrays.
!     ------------------------------------------------------------------

                          cldtyp(1,i,j) = neftyp
                          cldamt(1,i,j) = max( cldamt(1,i,j), nefamt )  

!     ------------------------------------------------------------------
!                   cld type is cumulonimbus.  however, if the 
!                   cloud top is below the latitudinal threshold
!                   set in the flux3 control file, set it to
!                   cumulus instead.
!     ------------------------------------------------------------------

                       else

                          if ( latbnd .eq. 5 ) then 
                             
                             if ( neftop .gt. thres(1) ) then
                                cldtyp(1,i,j) = 1 
                             else
                                cldtyp(1,i,j) = 4 
                             endif
                             cldamt(1,i,j) = max( cldamt(1,i,j), nefamt )
                             
                          elseif ( latbnd .eq. 4 ) then 
                             
                             if ( neftop .gt. thres(2) ) then
                                cldtyp(1,i,j) = 1 
                             else
                                cldtyp(1,i,j) = 4 
                             endif
                             cldamt(1,i,j) = max( cldamt(1,i,j), nefamt )
                             
                          elseif ( latbnd .eq. 3 ) then 
                             
                             if ( neftop .gt. thres(3) ) then
                                cldtyp(1,i,j) = 1 
                             else
                                cldtyp(1,i,j) = 4 
                             endif
                             cldamt(1,i,j) = max( cldamt(1,i,j), nefamt )
                             
                          elseif ( latbnd .eq. 2 ) then 
                             
                             if ( neftop .gt. thres(4) ) then
                                cldtyp(1,i,j) = 1 
                             else
                                cldtyp(1,i,j) = 4 
                             endif
                             cldamt(1,i,j) = max( cldamt(1,i,j), nefamt )
                             
                          else  
                             
                             if ( neftop .gt. thres(5) ) then
                                cldtyp(1,i,j) = 1 
                             else
                                cldtyp(1,i,j) = 4 
                             endif
                             cldamt(1,i,j) = max( cldamt(1,i,j), nefamt )
                             
                          endif
                       endif
                    endif
                 endif
              endif
           endif
        enddo
     enddo
  enddo
  
  return

end subroutine AGRMET_loadcloud
