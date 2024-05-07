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
! !ROUTINE: AGRMET_storeobs
! \label{AGRMET_storeobs}
! 
! !REVISION HISTORY: 
!    21 feb 01   initial version...........................mr gayno/dnxm
!    25 feb 02   store misc amounts for russian obs........mr gayno/dnxm
!     9 sep 10   modified for JMOBS...........Chris Franks/16WS/WXE/SEMS
!     7 feb 11   enable use of either JMOBS platform id or CDMS
!                block stn number.............Chris Franks/16WS/WXE/SEMS
!    11 may 11   improve processing of obs from India and Sri Lanka
!                .............................Chris Franks/16WS/WXE/SEMS
!
! !INTERFACE: 
subroutine AGRMET_storeobs(nsize, nsize3, isize, obs, obs3, ilat, ilon,  &
     mscprc, sixprc, twfprc, network, plat_id, cdms_flag, bsn, &
     duration, julhr, stncnt)

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
  integer,    intent(in)         :: nsize3
  integer,    intent(in)         :: sixprc(isize)
  integer,    intent(in)         :: julhr
  integer,    intent(inout)      :: stncnt
  integer,    intent(in)         :: twfprc(isize)   

!
! !DESCRIPTION: 
!    performs some preprocessing on the raw observations and stores
!    valid data in a storage array.
!    \textbf{Method} \newline
!    - ensure observation has a valid latitude, longitude \newline
!    - if processing the first WMO ob set a flag \newline
!    - if done processing 6-hourly WMO obs and 3-hourly obs remain, 
!        insert those 3-hourly obs \newline
!    - if a 3-hourly ob needs to be inserted here, do so \newline
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
!      - if the obs is for 6Z and there is a 3Z ob for this station
!        sum the miscellanous amounts from both and store it in 
!        the miscellaneous part of the obs data structure \newline
!    - else store the miscellanous amount in the miscellaneous 
!        part of the obs data structure \newline
!        - if the miscellaneous amount is less than 1 mm,  
!          set amount to 1 mm. \newline
!    - if the duration is 3 hours and this is a Sri Lankan ob... \newline
!        - if there is an ob for this station for 3 hours ago 
!          sum the miscellanous amounts from both and store it in 
!          the miscellaneous part of the obs data structure \newline
!        - else store the miscellanous amount in the miscellaneous 
!          part of the obs data structure \newline
!        - if the miscellaneous amount is less than 1 mm,  
!          set amount to 1 mm. \newline
!    - if there is a valid 6 hourly amount... \newline
!      - if amount is less than 1 mm, as indicated by a flag
!        of -91 thru -98, set amount to 1 mm. \newline
!      - store amount in 6 hourly part of obs data structure
!    - if the observation is from russia (or any part of the
!        former soviet union) store the misc amount.  russia
!        reports every 12 hours but does not use a valid
!        duration flag.  therefore, the misc amount will be
!        used as a 12 hourly amount in later routines. \newline
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
!   \item[obs]       array of processed observations \newline
!       net      JMOBS network of the observation
!       platform JMOBS platform ID of the observation
!       amt6     six hourly precip amount \newline
!       amt12    twelve hourly precip amount \newline
!       amt24    twenty-four hourly precip amount \newline
!       lat      latitude \newline
!       lon      longitude \newline
!       wmonum   wmo block station number \newline
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
!   \item[obs3\_ptr]  allocatable to next 3-hour obervation to process
!   \item[begin\_wmo] flag indicating that processing of CDMS or WMO obs has begun
!   \item[india\_lowlimit]      First 5-digit WMO num or 6-digit BSN for India
!   \item[india\_highlimit]     Last 5-digit WMO num or 6-digit BSN for India
!   \item[srilanka\_lowlimit]   First 5-digit WMO num or 6-digit BSN for Sri Lanka
!   \item[srilanka\_highlimit]  Last 5-digit WMO num or 6-digit BSN for Sri Lanka
!   \item[russia\_lowlimit]     First 5-digit WMO num or 6-digit BSN for Russia
!   \item[russia\_highlimit]    Last 5-digit WMO num or 6-digit BSN for Russia
!   \end{description}  
!EOP
  character*10                   :: date10
  integer                        :: india_lowlimit
  integer                        :: india_highlimit
  integer                        :: srilanka_lowlimit
  integer                        :: srilanka_highlimit
  integer                        :: irecord
  integer                        :: obs3_ptr
  integer                        :: temp6
  integer                        :: temp12
  integer                        :: temp24
  integer                        :: tempmsc
  real                           :: rlat
  real                           :: rlon
  integer                        :: russia_lowlimit
  integer                        :: russia_highlimit
  logical                        :: begin_wmo
  
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
  type(rain_obs), intent(in)     :: obs3(isize)

!     ------------------------------------------------------------------
!     If observations were retrieved from CDMS use 6-digit BSN limits 
!     for the specially processed regions
!     Else the observations were retrieved from JMOBS use 5-digit WMO 
!     numbers
!     ------------------------------------------------------------------

  begin_wmo = .FALSE.
  if (cdms_flag) then
     russia_lowlimit = 200000
     russia_highlimit = 390000
     india_lowlimit = 420010
     india_highlimit = 433990
     srilanka_lowlimit = 434000
     srilanka_highlimit = 434970
  else
     russia_lowlimit = 20000
     russia_highlimit = 39000
     india_lowlimit = 42001
     india_highlimit = 43399
     srilanka_lowlimit = 43400
     srilanka_highlimit = 43497
  end if
  obs3_ptr = 1

!     ------------------------------------------------------------------
!     Convert date/time to julian hour
!     ------------------------------------------------------------------

  call AGRMET_julhr_date10(julhr, date10)

!     ------------------------------------------------------------------
!     loop through all obs in the arrays retrieved from database
!     ------------------------------------------------------------------

  RECORD : do irecord = 1, nsize, 1

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
!     If bsn is not zero either the observations came from CDMS or they
!     came from JMOBS and this is a WMO observation.  Set the flag to 
!     indicate that we've begun processing these
!     Else if bsn is zero and we did process some WMO observations, 
!     insert any remaining 3-hourly observations in the array (this 
!     should never be executed when using CDMS)
!     ------------------------------------------------------------------

      if ((bsn(irecord) > 0) .and. (begin_wmo .eqv. .FALSE.)) then
         begin_wmo = .TRUE.
      else if ((begin_wmo .eqv. .TRUE.) .and. &
               (bsn(irecord) .eq. 0)   .and. &
               (obs3_ptr < nsize3)) then

         do while (obs3_ptr <= nsize3) 

           stncnt = stncnt + 1
           obs(stncnt)%net    = obs3(obs3_ptr)%net
           obs(stncnt)%platform = obs3(obs3_ptr)%platform
           obs(stncnt)%wmonum = obs3(obs3_ptr)%wmonum
           obs(stncnt)%lat    = obs3(obs3_ptr)%lat
           obs(stncnt)%lon    = obs3(obs3_ptr)%lon
           obs(stncnt)%amt24  = obs3(obs3_ptr)%amt24
           obs(stncnt)%amt12  = obs3(obs3_ptr)%amt12
           obs(stncnt)%amt6   = obs3(obs3_ptr)%amt6
           obs(stncnt)%amtmsc = obs3(obs3_ptr)%amtmsc

           obs3_ptr = obs3_ptr + 1 

         end do
      end if

!     ------------------------------------------------------------------
!       If there are obs for a station in India or Sri Lanka from 3 hours
!       ago which are not in this set of obs, store it in the obs array.
!     ------------------------------------------------------------------

      if ((nsize3 > 0) .and. (obs3_ptr <= nsize3) .and. &
          (bsn(irecord) > obs3(obs3_ptr)%wmonum)) then

         do while ((bsn(irecord) > obs3(obs3_ptr)%wmonum) .and. &
                (obs3_ptr <= nsize3)) 

           stncnt = stncnt + 1
           obs(stncnt)%net    = obs3(obs3_ptr)%net
           obs(stncnt)%platform = obs3(obs3_ptr)%platform
           obs(stncnt)%wmonum = obs3(obs3_ptr)%wmonum
           obs(stncnt)%lat    = obs3(obs3_ptr)%lat
           obs(stncnt)%lon    = obs3(obs3_ptr)%lon
           obs(stncnt)%amt24  = obs3(obs3_ptr)%amt24
           obs(stncnt)%amt12  = obs3(obs3_ptr)%amt12
           obs(stncnt)%amt6   = obs3(obs3_ptr)%amt6
           obs(stncnt)%amtmsc = obs3(obs3_ptr)%amtmsc

           obs3_ptr = obs3_ptr + 1 

         end do
      end if

!     ------------------------------------------------------------------
!       check for valid precipitation totals.  if they exist, process this ob.
!       Else check for an Indian or Sri Lankan ob from 3 hours ago for 
!       this bsn.  If one exists, copy it to the output array.
!     ------------------------------------------------------------------

     if ( (twfprc(irecord) <= -99999998) .and. &
          (mscprc(irecord) <= -99999998) .and. &
          (sixprc(irecord) <= -99999998) ) then
        if ((bsn(irecord) == obs3(obs3_ptr)%wmonum) .and. &
            (obs3_ptr <= nsize3)) then

           stncnt = stncnt + 1
           obs(stncnt)%net    = obs3(obs3_ptr)%net
           obs(stncnt)%platform = obs3(obs3_ptr)%platform
           obs(stncnt)%wmonum = obs3(obs3_ptr)%wmonum
           obs(stncnt)%lat    = obs3(obs3_ptr)%lat
           obs(stncnt)%lon    = obs3(obs3_ptr)%lon
           obs(stncnt)%amt24  = obs3(obs3_ptr)%amt24
           obs(stncnt)%amt12  = obs3(obs3_ptr)%amt12
           obs(stncnt)%amt6   = obs3(obs3_ptr)%amt6
           obs(stncnt)%amtmsc = obs3(obs3_ptr)%amtmsc

           obs3_ptr = obs3_ptr + 1 

        end if
        cycle RECORD
     end if
     
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
              if ((bsn(irecord) >= india_lowlimit) .and. &
                  (bsn(irecord) <= india_highlimit)) then
                 tempmsc = mscprc(irecord)
              end if

!     -----------------------------------------------------------------
!             India reports precipitation as accumulations
!             from 03z, not 24, 12 or 6 hourly amounts.  place it 
!             in the miscellaneous part of the obs array. 
!             If we're processing 6Z and there's a matching ob from 3Z
!             set the miscellaneous value to the sum of the two.
!             else place the value from the 6-hourly ob in the 
!             miscellaneous part of the obs array. 
!     -----------------------------------------------------------------

           elseif ((bsn(irecord) >= india_lowlimit) .and. &
                   (bsn(irecord) <= india_highlimit)) then

              if (( date10(9:10) == "06" ) .and. &
                  (bsn(irecord) == obs3(obs3_ptr)%wmonum) .and. &
                  (obs3_ptr <= nsize3)) then
                 tempmsc = mscprc(irecord) + obs3(obs3_ptr)%amtmsc
                 temp24 = obs3(obs3_ptr)%amt24
                 obs3_ptr = obs3_ptr + 1
              else
                 tempmsc = mscprc(irecord)
              end if
              if (tempmsc < 0 .and. tempmsc > -99) tempmsc = 1

!     -----------------------------------------------------------------
!             Sri Lanka reports 3 hour precipitation accumulations
!             not 24, 12 or 6 hourly amounts.  place it 
!             in the miscellaneous part of the obs array and add the 
!             amount from 3 hours ago.
!     -----------------------------------------------------------------

           elseif ((duration(irecord) == 3) .and. &
                ((bsn(irecord) >= srilanka_lowlimit) .and. &
                (bsn(irecord) <= srilanka_highlimit))) then

              if ((bsn(irecord) == obs3(obs3_ptr)%wmonum) .and. &
                  (obs3_ptr <= nsize3)) then
                 tempmsc = mscprc(irecord) + obs3(obs3_ptr)%amtmsc
                 obs3_ptr = obs3_ptr + 1
              else
                 tempmsc = mscprc(irecord)
              end if
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
!           many russian stations report every 12 hours, but do not 
!           report a valid duration.  therefore, we can't trust the
!           above logic to get us a valid amount.  however, most
!           russian obs report a misc amount, which we will use as
!           a 12 hourly amount in later routines.
!     ------------------------------------------------------------------

        if (bsn(irecord) >= russia_lowlimit .and. &
             bsn(irecord) <  russia_highlimit ) then

           tempmsc = mscprc(irecord)
           if (tempmsc < 0 .and. tempmsc > -99) tempmsc = 1
           
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

end subroutine AGRMET_storeobs
