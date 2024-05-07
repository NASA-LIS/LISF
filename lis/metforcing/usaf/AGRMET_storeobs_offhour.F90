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
! !ROUTINE: AGRMET_storeobs_offhour
! \label{AGRMET_storeobs_offhour}
! 
! !REVISION HISTORY: 
!    11 may 11   Adapted from AGRMET_storeobs...Chris Franks/16WS/WXE/SEMS
!
! !INTERFACE: 
subroutine AGRMET_storeobs_offhour(nsize, isize, obs, ilat, ilon,  &
     mscprc, sixprc, twfprc, network, plat_id, cdms_flag, bsn, &
     duration, stncnt)

  implicit none
  
  integer,    intent(in)         :: isize
  character*10, intent(in)       :: network(isize)
  character*10, intent(in)       :: plat_id(isize)
  logical,    intent(in)         :: cdms_flag
  integer,    intent(in)         :: bsn(isize)
  integer,    intent(in)         :: duration(isize)
  integer,    intent(in)         :: ilat(isize)
  integer,    intent(in)         :: ilon(isize)
  integer,    intent(in)         :: mscprc(isize)
  integer,    intent(in)         :: nsize
  integer,    intent(in)         :: sixprc(isize)
  integer,    intent(inout)      :: stncnt
  integer,    intent(in)         :: twfprc(isize)   

!
! !DESCRIPTION: 
!    performs some preprocessing on raw 3-hourly observations for 
!    India and Sri Lanka and stores valid data in a storage array.
!    \textbf{Method} \newline
!    - ensure the observation is for India or Sri Lanka
!    - ensure observation has a valid latitude, longitude, 
!      and at least one valid precip amount. \newline
!    - if there is a valid 24 hourly amount... \newline
!      - if amount is zero, then we know the 12 and 6 hourly amounts
!        must be zero as well. \newline
!      - if amount is less than 1 mm, as indicated by a flag
!        of -91 thru -98, set amount to 1 mm. \newline
!      - store amount in 24 hourly part of obs data structure.
!    - if there is a valid 12 hourly amount... \newline
!      - if amount is zero, then we know the 6 hourly amount
!        must be zero as well. \newline
!      - if amount is less than 1 mm, as indicated by a flag
!        of -91 thru -98, set amount to 1 mm. \newline
!      - store amount in 12 hourly part of obs data structure
!    - if the observation is from india... \newline
!      - if amount is less than 1 mm, as indicated by a flag
!        of -91 thru -98, set amount to 1 mm. \newline
!      - store amount in miscellaneous and 24-hourly part of 
!        obs data structure
!    - if the observation is from Sri Lanka \& the duration is 3 hour... \newline
!      - if amount is less than 1 mm, as indicated by a flag
!        of -91 thru -98, set amount to 1 mm. \newline
!      - store amount in miscellaneous part of 
!        obs data structure
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[bsn]       wmo block station number array
!   \item[duration]  valid period of precipitation
!   \item[ilat]      array of observation latitudes
!   \item[ilon]      array of observation longitudes
!   \item[irecord]   loop index
!   \item[isize]     maximum size of observation arrays
!   \item[mscprc]    miscellaneous precip amounts
!   \item[network]   JMOBS network, i.e. WMO, ICAO, FAA, etc.
!   \item[nsize]     number of obs returned from database
!   \item[obs]       array of processed observations 
!   \item[net]       JMOBS network of the observation
!   \item[platform]  JMOBS platform ID of the observation
!   \item[amt6]      six hourly precip amount
!   \item[amt12]     twelve hourly precip amount
!   \item[amt24]     twenty-four hourly precip amount
!   \item[lat]       latitude
!   \item[lon]       longitude
!   \item[wmonum]    wmo block station number
!   \item[plat\_id]   JMOBS platform ID
!   \item[rlat]      temporary holding variable for a latitude
!   \item[rlon]      temporary holding variable for a latitude
!   \item[sixprc]    array of six hourly precip amounts
!   \item[stncnt]    number of observations with a precip report
!   \item[temp6]     temporary holding variable for a 6 hourly amount
!   \item[temp12]    temporary holding variable for a 12 hourly amount
!   \item[temp24]    temporary holding variable for a 24 hourly amount
!   \item[tempmsc]   temporary holding variable for a misc amount
!   \item[twfprc]    array of 24 hourly precip amounts
!   \item[india\_lowlimit]      First 5-digit WMO or 6-digit BSN for India
!   \item[india\_highlimit]     Last 5-digit WMO or 6-digit BSN for India
!   \item[srilanka\_lowlimit]   First 5-digit WMO or 6-digit BSN for Sri Lanka
!   \item[srilanka\_highlimit]  Last 5-digit WMO or 6-digit BSN for Sri Lanka
!   \end{description}  
!EOP
  integer                        :: india_lowlimit
  integer                        :: india_highlimit
  integer                        :: srilanka_lowlimit
  integer                        :: srilanka_highlimit
  integer                        :: irecord
  integer                        :: temp6
  integer                        :: temp12
  integer                        :: temp24
  integer                        :: tempmsc
  real                           :: rlat
  real                           :: rlon
  
  type rain_obs
     sequence
     character*10                 :: net
     character*10                 :: platform
     integer                      :: wmonum
     real                         :: lat
     real                         :: lon
     integer                      :: amt24
     integer                      :: amt12
     integer                      :: amt6
     integer                      :: amtmsc
  end type rain_obs
  
  type(rain_obs), intent(inout)  :: obs(isize)

!     ------------------------------------------------------------------
!     If observations were retrieved from CDMS use 6-digit BSN limits 
!     for the specially processed regions
!     Else the observations were retrieved from JMOBS use 5-digit WMO 
!     numbers
!     ------------------------------------------------------------------

  if (cdms_flag) then
     india_lowlimit = 420010
     india_highlimit = 433999
     srilanka_lowlimit = 434000
     srilanka_highlimit = 434970
  else
     india_lowlimit = 42001
     india_highlimit = 43399
     srilanka_lowlimit = 43400
     srilanka_highlimit = 43497
  end if

!     ------------------------------------------------------------------
!     loop through all obs in the arrays retrieved from database
!     ------------------------------------------------------------------

  RECORD : do irecord = 1, nsize, 1

!     ------------------------------------------------------------------
!     Only process obs from India or Sri Lanka
!     ------------------------------------------------------------------

     if ((bsn(irecord) < india_lowlimit) .or. &
         (bsn(irecord) > srilanka_highlimit)) then
        cycle RECORD
     end if

!     ------------------------------------------------------------------
!       check for valid latitude and longitude.  if invalid skip
!       to next ob.
!     ------------------------------------------------------------------
     
     if ( (ilat(irecord) > -99999998) .and.  &
          (ilon(irecord) > -99999998) ) then
        rlat = float(ilat(irecord)) * 0.01
        rlon = float(ilon(irecord)) * 0.01
     else
        cycle RECORD
     end if
     
!     ------------------------------------------------------------------
!       check for valid wmo block station number (bsn) and valid
!       precipitation totals.  if they exist, process this ob.
!     ------------------------------------------------------------------

     if ( (twfprc(irecord) <= -99999998) .and. &
          (mscprc(irecord) <= -99999998) .and. &
          (sixprc(irecord) <= -99999998) ) cycle RECORD
     
        temp24  = -99999999
        temp12  = -99999999
        temp6   = -99999999
        tempmsc = -99999999

!     ------------------------------------------------------------------
!           in our database, negative numbers correspond to amounts
!           less than one mm.  for example, 0.1 mm is stored as -91.
!           set these small amounts to 0.01 inches.  0.8 mm is
!           stored as -98.  there are no numbers less than -98
!     ------------------------------------------------------------------

        temp24 = twfprc(irecord)
        
        if ( temp24 < 0 .and. temp24 > -99 ) temp24 = 1

!     ------------------------------------------------------------------
!           if there was no rain during the previous 24 hours, then
!           the 12 and 6 hourly amounts must be zero as well.
!     ------------------------------------------------------------------

        if (temp24 == 0) then
           
           temp12 = 0
           temp6  = 0
           
        else
           
           if (duration(irecord) == 12) then
              
              temp12 = mscprc(irecord)
              if (temp12 < 0 .and. temp12 > -99) temp12 = 1
              tempmsc = mscprc(irecord)
              if (tempmsc < 0 .and. tempmsc > -99) tempmsc = 1

!     -----------------------------------------------------------------
!             India report precipitation as accumulations
!             from 03z, not 24, 12 or 6 hourly amounts.  place it 
!             in the miscellaneous part of the obs array.
!             Save the miscellanous as 24-hour precip in order to 
!             separate 0-3Z from 3-6Z accumulation in processobs.
!     -----------------------------------------------------------------

           elseif ((bsn(irecord) >= india_lowlimit) .and. &
                   (bsn(irecord) <= india_highlimit)) then

              tempmsc = mscprc(irecord)
              if (tempmsc < 0 .and. tempmsc > -99) tempmsc = 1
              if (temp24 == -99999999) temp24 = tempmsc

              
!     -----------------------------------------------------------------
!             Sri Lanka reports 3-hourly precipitation 
!             from 03z, not 24, 12 or 6 hourly amounts.  place it 
!             in the miscellaneous part of the obs array.
!     -----------------------------------------------------------------

           elseif ((duration(irecord) == 3) .and. &
                ((bsn(irecord) >= srilanka_lowlimit) .and. &
                (bsn(irecord) <= srilanka_highlimit))) then

              tempmsc = mscprc(irecord)
              if (tempmsc < 0 .and. tempmsc > -99) tempmsc = 1

           end if

!     ------------------------------------------------------------------
!             if there was no rain during the past 12 hours, then
!             the 6-hourly amount must also be zero.
!     ------------------------------------------------------------------

           if (temp12 == 0) then
              
              temp6   = 0
              
           else
              
              temp6 = sixprc(irecord)

!     ------------------------------------------------------------------
!               set small amounts less than 1 mm to 0.01 inches. 
!     ------------------------------------------------------------------

              if ( temp6 < 0 .and. &
                   temp6 > -99) temp6 = 1
              
           end if
           
        end if

!     ------------------------------------------------------------------
!           if the observation has at least one reported precip 
!           amount, store it in the obs array.
!     ------------------------------------------------------------------
        
        if (temp24 >=0 .or. temp12  >=0 .or. &
             temp6  >=0 .or. tempmsc >= 0) then
           
           stncnt = stncnt + 1
           
           obs(stncnt)%net    = network(irecord)
           obs(stncnt)%platform = plat_id(irecord)
           obs(stncnt)%wmonum = bsn(irecord)
           obs(stncnt)%lat    = rlat
           obs(stncnt)%lon    = rlon
           
           obs(stncnt)%amt24  = temp24
           obs(stncnt)%amt12  = temp12
           obs(stncnt)%amt6   = temp6
           obs(stncnt)%amtmsc = tempmsc
           
        end if
        
  enddo RECORD

  return

end subroutine AGRMET_storeobs_offhour
