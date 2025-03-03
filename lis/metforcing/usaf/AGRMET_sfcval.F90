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
! !ROUTINE: AGRMET_sfcval
! \label{AGRMET_sfcval}
! 
! !REVISION HISTORY: 
!     10 Apr 1989  Initial version.......................MSgt Neill/SDDC   
!     11 Feb 1991  added relative humidity value ceiling (100%) and
!                  floor (10%). added declaration for i and j ..........
!                  ........................................Mr Moore/SDDC
!     08 Oct 1997  variables jul, adr, data, fldbld, z6, and chour6
!                  added to argument list so retrieval of 9 parameters 
!                  could be moved here from driver.  variables tmp10,
!                  tmp85, tmp70, dval10, dval85, dval70, relh10, and
!                  relh85 deleted from argument list since equivalencing
!                  of data array to these variables was moved here from 
!                  driver. These variables also became local variables. 
!                  Brought code up to current afgwc and sysm standards..
!                  ...............Mr Moore(AGROMET), SSgt McCormick/SYSM
!     04 Mar 1999  changed to supergrid structure. corrected linear
!                  interpolation error. removed julian hour check for 
!                  each box. Added vertical interpolation of relative 
!                  humidity. Set consistant flagging of water values for 
!                  proper use by barnes. modified code for opening and
!                  reading fldbld files.......capt andrus, mr moore/dnxm
!     11 Aug 1999  minor bug fixes during port to IBM SP-2.  added
!                  intent attributes to arguments.  made grid specific
!                  variables dynamically allocatable. removed 
!                  nogaps conversion of dewpoint depression to
!                  relative humidity as this is now done in module
!                  fldbld.  pass in field build directory path
!                  instead of hardwiring it here...........mr gayno/dnxm
!     01 Aug 2002  modified routine to be more generic so it could
!                  handle any number of isobaric levels....mr gayno/dnxm
! 
! !INTERFACE: 
subroutine AGRMET_sfcval( alt, jul, &!land, 
     minwnd, sfcprs, sfcrlh, &
     sfcspd, sfctmp, imax, jmax, &
     kprs, prslvls, kprs1,tmp, hgt, rlh, wndspd, pres)

  use LIS_LMLCMod, only : LIS_LMLC
  implicit none
! !ARGUMENTS:   

  integer,  intent(in)       :: kprs
  integer,  intent(in)       :: kprs1
  integer,  intent(in)       :: imax
  integer,  intent(in)       :: jmax
  real,     intent(in)       :: alt  ( imax, jmax )
  real,     intent(in)       :: hgt   ( imax, jmax, kprs1 )
  real,     intent(in)       :: rlh    ( imax, jmax, kprs1 )
  real,     intent(in)       :: tmp   ( imax, jmax, kprs1 )
  real,     intent(in)       :: wndspd  ( imax,jmax)
  real,     intent(in)       :: pres  ( imax,jmax)
  integer,  intent(in)       :: jul
!  real,     intent(in)       :: land      ( imax, jmax)
  integer,  intent(in)       :: prslvls   ( 30 )
  real,     intent(in)       :: minwnd
  real,     intent(out)      :: sfcprs    ( imax,jmax )
  real,     intent(out)      :: sfcrlh    ( imax,jmax )
  real,     intent(out)      :: sfcspd    ( imax,jmax )
  real,     intent(out)      :: sfctmp    ( imax,jmax )

! !DESCRIPTION: 
!  
!     vertically interpolate input fldbld fields to the model surface
!     
!     \textbf{Method} \newline
!     - allocate grid specific arrays. \newline
!     - set file names and retrieve all input fldbld parameters \newline
!     - over land, vertically interpolate isobaric fldbld data to the
!       model surface.  nogaps dewpoint depression was converted
!       to relative humidity in module fldbld. \newline
!       - if model terrain height is below the highest pressure level,
!         set surface temperature, pressure and relative humidity to
!         the values at that highest level. \newline
!       - if model terrain height is above the lowest pressure level,
!         set surface temperature, pressure and relative humidity to
!         the values at that lowest pressure \newline
!       - if model terrain height is between the pressure levels,
!         vertically interpolate temperature and relative humidity.
!         surface pressure is estimated using the hypsometric equation. \newline
!     - surface winds are not vertically interpolated, they are
!       simply set to the values calculated in module fldbld. \newline
!     - range check all data. \newline
!     - deallocate grid specific arrays. \newline
!    
!  The arguments and variables are: 
!  \begin{description}
!  \item[alt]
!   terrain elevation for the agrmet grid
!  \item[jul]
!   current julian hour
!  \item[land]
!   array of land mass points for agrmet grid
!  \item[minwnd]
!    minimum allowable windspeed on the agrmet grid   
!  \item[sfcprs]
!    first guess pressure fields on the AGRMET grid interpolated to 
!    the current time
!  \item[sfcrlh]
!    first guess rh fields on the AGRMET grid interpolated to 
!    the current time
!  \item[sfcspd]
!    first guess wind speed fields on the AGRMET grid interpolated to 
!    the current time   
!  \item[sfctmp]
!    first guess temperature fields on the AGRMET grid interpolated to 
!    the current time
!  \item[imax]
!    east/west dimension of agrmet grid
!  \item[jmax]
!    north/south dimension of agrmet grid
!  \item[kprs1]
!    number of isobaric levels
!  \item[prslvls]
!   isobaric levels in first guess data
!  \item[tmp]
!   temperature on the agrmet grid 
!  \item[hgt]
!   geopotential heights on the agrmet grid
!  \item[rlh]
!   relative humidity on the agrmet grid
!  \item[wndspd]
!   wind speeds on the agrmet grid
!  \item[bashgt]
!   height of bottom of interpolation layer
!  \item[basp]
!   pressure of bottom of interpolation layer
!  \item[basrh]
!   relative humidity of bottom of interpolation layer
!  \item[bastmp]
!   temperature at bottom of interpolation layer
!  \item[calprs]
!    function for calculating surface pressure on the
!    agrmet grid
!  \item[h1]
!    function that calculates (R * T) / g which is
!    used to hypsometrically calculate surface
!    pressure on the agrmet grid
!  \item[linint]
!    function for vertically interpolating isobaric
!    data to the agrmet grid
!  \item[thick]
!    thickness of interpolation layer
!  \item[toprh]
!    relative humidity of top of interpolation layer
!  \item[toptmp]
!    temperature of top of interpolation layer
!  \item[ttop]
!   temperature of top of interpolation layer 
!  \end{description}
!
!EOP

  integer                       :: i
  integer                       :: j 
  integer                       :: kk   
  real                          :: aalt
  real                          :: abas
  real                          :: athick
  real                          :: bashgt
  real                          :: basp  
  real                          :: basrh
  real                          :: bastmp
  real                          :: calprs
  real                          :: h1
  real                          :: linint
  real                          :: p 
  real                          :: t1
  real                          :: t2
  real                          :: tbas
  real                          :: thick 
  real                          :: toprh
  real                          :: toptmp
  real                          :: ttop
  real                          :: za
  real                          :: zb
  real       :: alt_t  ( imax, jmax )
!  real :: ttmp(512,512)

!     ------------------------------------------------------------------
!     define internal functions h1 and calprs for calculation of
!     final surface pressure using the hypsometric equation.  also 
!     define internal function linint for vertically interpolating
!     from isobaric levels to the model surface. 
!     ------------------------------------------------------------------

  h1(t1, t2) = 287.04 * ((t1 + t2) * 0.5) / 9.80665
  
  calprs(t1, t2, za, zb, p) = max (min  &
       (p * exp(-((za - zb) /  &
       (h1 (t1, t2)))) * 100.0, &
       110000.0), 40000.0)
  
  linint(ttop, tbas, athick, aalt, abas) = tbas + ((ttop - tbas) *   &
       (aalt - abas) / athick)
    
!     ------------------------------------------------------------------
!     loop through each model grid point.
!     ------------------------------------------------------------------
  alt_t = alt
  JGRID : do j = 1, jmax
     IGRID : do i = 1, imax
!     ------------------------------------------------------------------
!         if grid point is a land point, proceed ...
!     ------------------------------------------------------------------
!        LAND_TEST : if( land(i,j) .gt. 0 )then
           if(alt_t(i,j).eq.-9999.0) alt_t(i,j) = 0.0
           VERT_TEST : if( alt_t(i,j) <= hgt(i,j,1) )then
!     ------------------------------------------------------------------
!             elevation is less than or equal to the height of the
!             highest mb level. set all outputs to this mb's values.
!     ------------------------------------------------------------------
              
              sfctmp(i,j) = tmp(i,j,1)
              sfctmp(i,j) = max( min( sfctmp(i,j), 350.0 ), 200.0 )
!              sfcprs(i,j) = float(prslvls(1)) * 100.0
!              sfcprs(i,j) = (float(prslvls(1)) - &
!                   float(prslvls(1))-float(prslvls(2))*(hgt(i,j,1)/&
!                   (hgt(i,j,1)-hgt(i,j,2))))*100.0
              sfcprs(i,j) = pres(i,j)
              sfcrlh(i,j) = rlh(i,j,1)
              sfcrlh(i,j) = max( min( sfcrlh(i,j), 1.0 ), 0.1 )
           else if ( alt_t(i,j) >= hgt(i,j,kprs) ) then
!     ------------------------------------------------------------------
!             elevation is above the height of the lowest mb level.
!             set all outputs to this mb's values.
!     ------------------------------------------------------------------
                 
              sfctmp(i,j) = tmp(i,j,kprs)
              sfctmp(i,j) = max( min( sfctmp(i,j), 350.0 ), 200.0 )
              sfcprs(i,j) = float(prslvls(kprs)) * 100.
              sfcrlh(i,j) = rlh(i,j,kprs)
              sfcrlh(i,j) = max( min( sfcrlh(i,j), 1.0 ), 0.1 )
           else

!     ------------------------------------------------------------------
!             elevation is between two mb levels.  search to find the
!             two levels that surround the terrain height.  then, 
!             verically interpolate.
!     ------------------------------------------------------------------
                 
              SEARCH : do kk = 1, kprs - 1

                 if ( (alt_t(i,j) >  hgt(i,j,kk)) .and. &
                      (alt_t(i,j) <= hgt(i,j,kk+1))  ) then
                    
                    basp   = float(prslvls(kk))
                    bastmp = tmp(i,j,kk)
                    toptmp = tmp(i,j,kk+1)
                    bashgt = hgt(i,j,kk)
                    thick  = hgt(i,j,kk+1) - hgt(i,j,kk)

!     ------------------------------------------------------------------
!                 using the variables calculated above and the
!                 functions defined earlier, calculate surface temp and
!                 pressure.  also limit the values.  note:  pressure
!                 is limited in the function which computes it.
!     ------------------------------------------------------------------
                    sfctmp(i,j) = linint( toptmp, bastmp, thick, &
                         alt_t(i,j), bashgt )
                    sfctmp(i,j) = max( min( sfctmp(i,j), 350.0 ), 200.0 )
                    
                    sfcprs(i,j) = calprs( bastmp, toptmp, alt_t(i,j),  &
                         bashgt, basp )

!     ------------------------------------------------------------------
!                 calculate the surface relative humidity from the 
!                 values at 1000 and 850 mb. the rh must be between 
!                 10% and 100%.
!     ------------------------------------------------------------------                
                    
                    basrh       = rlh(i,j,kk)
                    toprh       = rlh(i,j,kk+1)
                    sfcrlh(i,j) = linint( toprh, basrh, thick,  &
                         alt_t(i,j), bashgt )
                    sfcrlh(i,j) = max( min( sfcrlh(i,j), 1.0 ), 0.1 )
                    
                 endif
                 
              enddo SEARCH
           end if VERT_TEST

!     ------------------------------------------------------------------
!           set the output surface wind speed to the input resultant
!           wind speed.  limit the values between a max and min.
!     ------------------------------------------------------------------
           
           sfcspd(i,j) = wndspd(i,j)
           sfcspd(i,j) = max( min( sfcspd(i,j), 75.0 ), minwnd )           
!         else ! NOT LAND
             
! !     ------------------------------------------------------------------
! !           grid pt is not a land pt. flag all output values.
! !           flags must all be the same for use in the barnes routine.
! !     ------------------------------------------------------------------
              
!            sfcrlh(i,j) = -1.0
!            sfcprs(i,j) = -1.0
!            sfctmp(i,j) = -1.0
!            sfcspd(i,j) = -1.0
           
!         endif LAND_TEST
        
     enddo IGRID
  enddo JGRID
!     ------------------------------------------------------------------
!     deallocate grid specific variables
!     ------------------------------------------------------------------

end subroutine AGRMET_sfcval
