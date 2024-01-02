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
! !ROUTINE: AGRMET_calest
!  \label{AGRMET_calest}
!
! !REVISION HISTORY:
!
!     1 jun 89  initial version...................mr moore/sddc(agromet)
!    11 feb 91  added use of precip probabilities in adjustment of
!               precip coef prior to calculation of final precip 
!               estimate.  eliminated use of cld amts in same process.
!               resized pcoef and pcout arrays.  isolated debug stmt....
!               ..................................mr moore/sddc(agromet)
!     4 may 91  chgd if stmts to allow ra, prcpwe, and rpcoef precip 
!               rates of 0.0 to be used. allowed cldtyps of 0 to be 
!               processed, changed prob array to allow for cldtyp of
!               0, updated the prolog.......capt bertone/sddc(agromet)  
!     3 mar 92  made chgs necessary to include cldamt as an index of
!               pcoef array. eliminated array pcout, use pcoef instead..
!               updtd prolog.......................mr moore/sys(agromet)
!     1 jun 92  expanded pcoef array to include rh class and stability  
!               class as indices. added cloud amount class' to condense 
!               10 cld amts into 5 keeping pcoef array smaller.  added  
!               use of pwgt........................mr moore sys(agromet)
!     1 feb 94  resized the pcoef array, deleting the stability class   
!               rh class, cloud type and lat band indices and adding
!               cld top ht class and tropopause ht class indices. chgd  
!               use of precip probabilities, making them a function of  
!               cld amt only.  started using pwgte rather than pwgt.
!               grouped comments and updtd prolog.......................
!               ...................mr moore, capt bertone, sysm(agromet)
!    25 apr 94  added 'files accessed' section as 'id'ed during the 
!               94-01 qaa....................capt bertone, sysm(agromet)
!    15 jun 96  reduced 'e' array size from (64,64,12) to (64,64,4).
!               replaced chk of ca with chk of rpcoef and eliminated
!               calc of pprob which is now done in parent routine.  
!               updated prolog and brought up to stds...................
!               ..................................mr moore sysm(agromet)
!    07 mar 97  removed precip coefficient, probability, and tropopause 
!               cloud height calculations.  added calculations for climo
!               estimate.  altered code to run more efficiently.
!               brought code up to current sysm and afgwc standards.
!               ........................capt andrus, ssgt mccormick/sysm
!    24 mar 97  added 'source' array which sets each type of input  
!               precip estimate with a specified integer source flag.   
!               array values are set according to determined precedence.
!               ........................................ssgt miller/sysm
!    19 may 97  redimensioned source array to write each 3-hour sourced 
!               estimate into its own file..............ssgt miller/sysm
!    04 dec 97  added ability to determine source of bogus, rtneph, 
!               and climo estimates. modified to allow only singular
!               estimate type per grid point rather than a combined
!               estimate. updated prolog and brought up to standards....
!               ...............................capt andrus/dnxm(agromet)
!    08 oct 99  ported to ibm sp-2, updated prolog, incorporated
!               FORTRAN 90 features, included the geo_precip estimate
!               into the heirarchy...................capt hidalgo/agrmet
!    01 jun 00  incorporated grnk values for geo_precip and used it to 
!               sort the order in which estimates will be used. Removed
!               estimated weights from the sort.....ssgt campbell/agrmet
!    21 feb 01  reduced third dimensions of all estimate arrays
!               (number of three hourly time periods) from 4 to 2
!               so module precip can run in 6 hourly cycles.  changed
!               variable name e to estpcp..................mr gayno/dnxm
!    10 jun 02  modified references and variable names to reflect
!               change from rtneph to cdfs2................mr gayno/dnxm
!
!     3 nov 05 Sujay Kumar, incorporated into LIS
!     1 aug 13 Ryan Ruhge, added CMORPH and GFS, removed cdfs2
!
! !INTERFACE:
subroutine AGRMET_calest( n, cest, fg_data, estpcp, gest, k, prcpwe,&
     ra, cmorphdata, cmorphpixel, source, imax, jmax, grnk )

! !USES:
  use LIS_coreMod,       only : LIS_domain
  use LIS_LMLCMod,       only : LIS_LMLC
  use AGRMET_forcingMod, only : agrmet_struc
  use LIS_snowMod,       only : LIS_snow_struc
  use LIS_logMod,        only : LIS_logunit, LIS_endrun

  implicit none
! !ARGUMENTS:
  integer, intent(in)               :: n
  integer, intent(in)               :: imax
  integer, intent(in)               :: jmax
  integer, intent(in)               :: grnk(imax,jmax)
  integer, intent(in)               :: k
  integer, intent(inout)            :: source(imax,jmax,4)
  real, intent(in)                  :: cest(imax,jmax)
  real, intent(inout)               :: estpcp(imax,jmax,4)
  real, intent(in)                  :: gest(imax,jmax)
  real, intent(in)                  :: prcpwe(imax,jmax,4)
  real, intent(in)                  :: ra(imax,jmax)
  real, intent(in)                  :: cmorphdata(imax, jmax)
  real, intent(in)                  :: fg_data(imax,jmax)
  logical, intent(out)              :: cmorphpixel(imax,jmax,4)

!
! !DESCRIPTION:
!
!    calculate a 'final' 3-hour estimated precipitation amount based
!    upon 3-hour precipitation estimate inputs, using a hierarchy
!    also to identify the source of the estimate placed in the 'final'
!    3-hour estimated precipitation array.
! 
!    \textbf{Method} \newline
!   
!    - determine the precedence given by the hierarchical values.  \newline
!    - for each gridpoint, check which estimates are available and
!      place the highest ranked estimate in the estimation field.  \newline
!    - identify the source of the estimate. \newline
! 
!  The arguments and variables are:
!  \begin{description}
!   \item[cest]       precip est based purely on climatology (mm/3hr)
!   \item[cmorphdata] cmorph precip estimate
!   \item[cmorphpixel] indicates if the current pixel comes from cmorph data
!   \item[ehit]       used to determine valid points
!   \item[estpcp]     final 3-hour estimated precip amt (mm/3hr)
!   \item[gest]       precip est from geo\_precip (mm/3hr)
!   \item[fg\_data]   3 hour precipitation forecast from first guess
!                     (GFS or GALWEM)
!   \item[grnk]       geo\_precip ranking value (1-5) \newline
!                               1 = better than ssmi \newline
!                               2 = better than present weather \newline
!                               3 = better than gfs \newline
!                               4 = default value for geo\_precip \newline
!                               5 = worse than climatological \newline
!   \item[hourindex]  indicates which index in sfctmp hourly array to use
!   \item[hemi]       hemisphere ( 1 = nh, 2 = sh)
!   \item[i]          loop counter, i-coordinate of the grid
!   \item[imax]       number of gridpoints in the east/west direction
!   \item[j]          loop counter, j-coordinate of the grid
!   \item[jmax]       number of gridpoints in the north/south direction
!   \item[k]          3-hour time counter
!   \item[l]          loop counter
!   \item[land]       land/sea mask (0 = water, 1 = land)
!   \item[latfield]   latitude field of entire grid
!   \item[latpcp]     precip value after interpolation by latitude
!   \item[latwt]      weight of latitude for interpolation
!   \item[nc1]        number of columns in nh
!   \item[nc2]        number of columns in sh
!   \item[nr1]        number of rows in nh
!   \item[nr2]        number of rows in sh
!   \item[nhlatfield] latitude field of northern hemisphere
!   \item[order]      array to sort the desired hierarchical rank
!   \item[prcpwe]     present weather estimated precip rate (mm/3hr)
!   \item[ra]         ssmi rain rate (mm/3hr)
!   \item[shlatfield] latitude field of southern hemisphere
!   \item[source]     array that contains estimate source information
!               with specified source flags \newline
!                              2 = ssmi est. \newline
!                              3 = present weather est. \newline
!                              4 = gfs est. \newline
!                              5 = geo\_precip/cmorph est. \newline
!                              6 = climatological est. \newline
!   \item[temp]       dummy scalar
!   \item[temp\_order] array used to hold shuffled values of order array
!   \item[threshdiff] difference in min and max temp thresholds for computing weights
!   \item[tmppcp]     precip value after interploation by temperature
!   \item[tmpwt]      weight of temperature for interpolation
!  \end{description}
!EOP
  integer                           :: hourindex
  integer                           :: i
  integer                           :: j
  integer                           :: l
  integer                           :: order(5)
  logical                           :: ehit
  integer                           :: temp_order(5)
  integer                           :: gid
  real                              :: latwt
  real                              :: tempwt
  real                              :: threshdiff
  real                              :: latpcp
  real                              :: tmppcp
  real, allocatable                 :: nhlatfield(:,:)
  real, allocatable                 :: shlatfield(:,:)
  real, allocatable                 :: latfield(:,:)
  integer                           :: nc1
  integer                           :: nc2
  integer                           :: nr1
  integer                           :: nr2
  real                              :: snowdepth

!     ------------------------------------------------------------------
!     executable code starts here ...
!     initialize order array
!     such that the first element 'order(1)' contains the
!     source number referring to the most important estimate.
!     ------------------------------------------------------------------

  if ( agrmet_struc(n)%gridspan == 1 ) then ! NH only
     nc1=agrmet_struc(n)%hemi_nc(1)
     nr1=agrmet_struc(n)%hemi_nr(1)
     nc2 = 0
     nr2 = 0
     allocate(nhlatfield(nc1, nr1))
     allocate(latfield(nc1, nr1))
     nhlatfield=RESHAPE(agrmet_struc(n)%rlat1_nh, (/nc1,nr1/))
     latfield(1:nc1,1:nr1)=nhlatfield
     deallocate(nhlatfield)
  elseif ( agrmet_struc(n)%gridspan == 2 ) then ! SH only
     nc1 = 0
     nr1 = 0
     nc2=agrmet_struc(n)%hemi_nc(2)
     nr2=agrmet_struc(n)%hemi_nr(2)
     allocate(shlatfield(nc2, nr2))
     allocate(latfield(nc2, nr2))
     shlatfield=RESHAPE(agrmet_struc(n)%rlat1_sh, (/nc2,nr2/))
     latfield(1:nc2,1:nr2)=shlatfield
     deallocate(shlatfield)
  elseif ( agrmet_struc(n)%gridspan == 3 ) then ! NH+SH
     nc1=agrmet_struc(n)%hemi_nc(1)
     nr1=agrmet_struc(n)%hemi_nr(1)
     nc2=agrmet_struc(n)%hemi_nc(2)
     nr2=agrmet_struc(n)%hemi_nr(2)
     allocate(nhlatfield(nc1, nr1))
     allocate(shlatfield(nc2, nr2))
     allocate(latfield(nc1, nr1+nr2))
     nhlatfield=RESHAPE(agrmet_struc(n)%rlat1_nh, (/nc1,nr1/))
     shlatfield=RESHAPE(agrmet_struc(n)%rlat1_sh, (/nc2,nr2/))
     latfield(1:nc2,1:nr2)=shlatfield
     latfield(1:nc1,nr2+1:nr1+nr2)=nhlatfield
     deallocate(nhlatfield)
     deallocate(shlatfield)
  else
     write(LIS_logunit,*) 'agrmet_struc%gridspan is not properly defined.'
     write(LIS_logunit,*) 'Program stopping.'
     call LIS_endrun()
  endif

  do i = 1,5
     order(i) = i
  enddo

!     ------------------------------------------------------------------
!     loop thru points in the hemi
!     if land pt, initialize hit flag and proceed ...
!     ------------------------------------------------------------------

  do j = 1, jmax
     do i = 1, imax
        cmorphpixel(i,j,k) = .false.
        if( LIS_LMLC(n)%landmask(i,j) .gt. 0 )then
           ehit = .false.

!     ------------------------------------------------------------------
!           if geo_precip and geo ranks are being used, initialize
!           temp_order array and change order of it based
!           upon the grnk of the geo_precip estimate.
!                       1 = ssmi
!                       2 = present weather
!                       3 = gfs
!                       4 = geo_precip
!                       5 = climatological
!           For grnk values of 1-3 the geo_precip estimate is moved to
!           that position in the array, meaning it is better than the
!           estimate that was there before, and there is a legit geo
!           estimate.  If grnk is a 5 then climo will be used before the
!           geo_precip estimate, since there is always climo, no need to
!           put 4 in position 5 of the temp_order array.
!           In the case of CMORPH being used instead of GEO_PRECIP/SSMI,
!           temp_order(4) is alway set to 5 and CMORPH is given highest
!           preference, followed by present weather, gfs, and climo
!     __________________________________________________________________

           temp_order = order

           if(grnk(i,j).ne.-9999.0) then 
              if( grnk(i,j) .eq. 5 ) then
                 temp_order(4) = 5
              else
                 temp_order(grnk(i,j))=4
                 temp_order(4)=grnk(i,j)
              endif
           endif
!     ------------------------------------------------------------------
!           loop thru 5 estimate types from the highest rank
!           to the lowest rank until an estimate is found
!     ------------------------------------------------------------------

            do l = 1, 5
               if( .not. ehit )then
!     ------------------------------------------------------------------
!               estimate not selected for this point yet.
!               proceed to select an estimate and identify a source.
!     ------------------------------------------------------------------


                  gid = LIS_domain(n)%gindex(i,j)
                  if(gid.eq.-1) then 
                     snowdepth = 0.0
                  else
                     snowdepth = LIS_snow_struc(n)%snowdepth(gid)
                  endif

                  if ( (k .eq. 1) .or. (k .eq. 3) ) then
                     hourindex = 2
                  else if ( (k .eq. 2) .or. (k .eq. 4) ) then
                     hourindex=5
                  endif

!     ------------------------------------------------------------------
!                 Only use CMORPH where no snow is on the ground and
!                 surface temperature is above minimum temperature
!                 threshold.  For areas outside the CMORPH domain
!                 (60S-60N) missing values (-9999) are used in the
!                 cmorphdata array.
!     ------------------------------------------------------------------!
                  if( (temp_order(l).eq.1) .and. (cmorphdata(i,j).gt.-9990.0) &
                     .and. (agrmet_struc(n)%sfctmp(hourindex,i,j) .ge. &
                     (agrmet_struc(n)%cmorminthresh)) &
                     .and. (snowdepth .le. 0) ) then


!     ------------------------------------------------------------------
!                 use CMORPH estimate
!     ------------------------------------------------------------------



!     ------------------------------------------------------------------
!                       CMORPH goes to 60 N/S, so do a linear
!                       interpolation between CMORPH and GFS in the
!                       50-60 latitude band to smooth the transition.
!     -------------------------------------------------------------------
                     if ((agrmet_struc(n)%sfctmp(hourindex,i,j) .lt. &
                        agrmet_struc(n)%cmormaxthresh)) then
                        threshdiff=agrmet_struc(n)%cmormaxthresh-agrmet_struc(n)%cmorminthresh
                        tempwt=(agrmet_struc(n)%sfctmp(hourindex,i,j)-agrmet_struc(n)%cmorminthresh)/threshdiff
                        if ((abs(latfield(i,j)) .gt. 50) .and. (abs(latfield(i,j)) .le. 60) &
                           .and. (fg_data(i,j) .gt. -9990.0)) then
!     -------------------------------------------------------------------
!                          If surface temperature is between the minimum
!                          and maximum threshold do an interpolation
!                          between CMORPH and GFS based on the temp and
!                          latitude.
!     -------------------------------------------------------------------
                           latwt=(abs(latfield(i,j))-50.0)/10.0
                           latpcp=((1-latwt)*cmorphdata(i,j))+(latwt*fg_data(i,j))
                           tmppcp=(tempwt*cmorphdata(i,j))+((1-tempwt)*fg_data(i,j))
                           estpcp(i,j,k)=(0.5*latpcp)+(0.5*tmppcp)
                           source(i,j,k)=5
                           ehit = .true.
                        else if (fg_data(i,j) .gt. -9990.0) then
! ------------------------------------------------------------------------
!                          Latitude is equatorward of 50 N/S so just
!                          interpolate between CMORPH and GFS based on
!                          temperature.
! ------------------------------------------------------------------------
                           estpcp(i,j,k)=(tempwt*cmorphdata(i,j))+((1-tempwt)*fg_data(i,j))
                           source(i,j,k)=5
                           ehit = .true.
                        endif
                     elseif ((abs(latfield(i,j)) .gt. 50) .and. (abs(latfield(i,j)) .le. 60) &
                        .and. (fg_data(i,j) .gt. -9990.0)) then
!     ------------------------------------------------------------------
!                       CMORPH goes to 60 N/S, so do a linear
!                       interpolation between CMORPH and GFS in the
!                       50-60 latitude band to smooth the transition.
!     -------------------------------------------------------------------
                        latwt=(abs(latfield(i,j))-50.0)/10.0
                        estpcp(i,j,k)=((1-latwt)*cmorphdata(i,j))+(latwt*fg_data(i,j))
                        source(i,j,k)=5
                        ehit = .true.
                     else
! -----------------------------------------------------------------------
!                       Inside all thresholds so just insert cmorph
!                       data.
! -----------------------------------------------------------------------
                        estpcp(i,j,k) = cmorphdata(i,j)
                        source(i,j,k) = 5
                        ehit = .true.
                        cmorphpixel(i,j,k) = .true.
                     endif
                  elseif( (temp_order(l).eq.1) .and. (ra(i,j).gt.-9990.0) )then

!     ------------------------------------------------------------------
!                 use ssmi estimate
!     ------------------------------------------------------------------

                     estpcp(i,j,k) = ra(i,j)
                     source(i,j,k) = 2
                     ehit = .true.
                  elseif( (temp_order(l).eq.2).and.&
                     (prcpwe(i,j,k).gt.-9990.0) )then

!     ------------------------------------------------------------------
!                 use present wx based estimate
!     ------------------------------------------------------------------

                     estpcp(i,j,k) = prcpwe(i,j,k)
                     source(i,j,k) = 3
                     ehit = .true.
                  elseif( (temp_order(l) .eq. 3) .and. &
                     (fg_data(i,j).gt.-9990.0) )then

!     ------------------------------------------------------------------
!                 use gfs precipitation based estimate
!     ------------------------------------------------------------------
                     estpcp(i,j,k) = fg_data(i,j)
                     source(i,j,k) = 4
                     ehit = .true.

!     ------------------------------------------------------------------
!                 Only use geoprecip where no snow is on the ground and
!                 surface temperature is above minimum temperature
!                 threshold.  For areas outside the geoprecip domain
!                 (50S-50N) missing values (-9999) are used in the
!                 gest array.
!     ------------------------------------------------------------------!
                  else if( (temp_order(l).eq.4) .and. (gest(i,j).gt.-9990.0) &
                     .and. (agrmet_struc(n)%sfctmp(hourindex,i,j) .ge. &
                     (agrmet_struc(n)%geominthresh)) &
                     .and. (snowdepth .le. 0) ) then


!     ------------------------------------------------------------------
!                 use geo_precip estimate
!     ------------------------------------------------------------------



!     ------------------------------------------------------------------
!                       geo_precip goes to 50 N/S, so do a linear
!                       interpolation between geo_precip and GFS in the
!                       40-50 latitude band to smooth the transition.
!     -------------------------------------------------------------------
                     if ((agrmet_struc(n)%sfctmp(hourindex,i,j) .lt. &
                        agrmet_struc(n)%geomaxthresh)) then
                        threshdiff=agrmet_struc(n)%geomaxthresh-agrmet_struc(n)%geominthresh
                        tempwt=abs(agrmet_struc(n)%sfctmp(hourindex,i,j)-agrmet_struc(n)%geominthresh)/threshdiff
                        if ((abs(latfield(i,j)) .gt. 40) .and. (abs(latfield(i,j)) .le. 50) &
                           .and. (fg_data(i,j) .gt. -9990.0)) then
!     -------------------------------------------------------------------
!                          If surface temperature is between the minimum
!                          and maximum threshold do an interpolation
!                          between GEORECIP and GFS based on the temp and
!                          latitude.
!     -------------------------------------------------------------------
                           latwt=(abs(latfield(i,j))-40.0)/10.0
                           latpcp=((1-latwt)*gest(i,j))+(latwt*fg_data(i,j))
                           tmppcp=(tempwt*gest(i,j))+((1-tempwt)*fg_data(i,j))
                           estpcp(i,j,k)=(0.5*latpcp)+(0.5*tmppcp)
                           source(i,j,k)=5
                           ehit = .true.
                        else if (fg_data(i,j) .gt. -9990.0) then
! ------------------------------------------------------------------------
!                          Latitude is equatorward of 40 N/S so just
!                          interpolate between GEOPRECIP and GFS based on
!                          temperature.
! ------------------------------------------------------------------------
                           estpcp(i,j,k)=(tempwt*gest(i,j))+((1-tempwt)*fg_data(i,j))
                           source(i,j,k)=5
                           ehit = .true.
                        endif
                     elseif ((abs(latfield(i,j)) .gt. 40) .and. (abs(latfield(i,j)) .le. 50) &
                        .and. (fg_data(i,j) .gt. -9990.0)) then
!     ------------------------------------------------------------------
!                       geo_precip goes to 50 N/S, so do a linear
!                       interpolation between geo_precip and GFS in the
!                       40-50 latitude band to smooth the transition.
!     -------------------------------------------------------------------
                        latwt=(abs(latfield(i,j))-40.0)/10.0
                        estpcp(i,j,k)=((1-latwt)*gest(i,j))+(latwt*fg_data(i,j))
                        source(i,j,k)=5
                        ehit = .true.
                     else
! -----------------------------------------------------------------------
!                       Inside all thresholds so just insert GEOPRECIP
!                       data.
! -----------------------------------------------------------------------
                        estpcp(i,j,k) = gest(i,j)
                        source(i,j,k) = 5
                        ehit = .true.
                     endif
                  elseif( (temp_order(l).eq.5).and.&
                     (cest(i,j).gt.-9990.0) )then

!     ------------------------------------------------------------------
!                 use pure climatological estimate
!     ------------------------------------------------------------------

                     estpcp(i,j,k) = cest(i,j)
                     source(i,j,k) = 6
                     ehit = .true.
                  endif
               endif
            enddo
         endif

      enddo
  enddo
  deallocate(latfield)
  return

end subroutine AGRMET_calest
