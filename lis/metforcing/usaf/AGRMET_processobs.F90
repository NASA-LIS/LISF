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
! !ROUTINE: AGRMET_processobs
! \label{AGRMET_processobs}
!
! !REVISION HISTORY: 
!     21 Feb 01  Initial version...........................Mr Gayno/DNXM
!     25 Feb 02  Modified to use misc array over Russia....Mr Gayno/DNXM
!     07 Jan 03  Expanded use of Russian obs at 00/12 utc.  Use all
!                available obs with block station numbers between 
!                20000 and 38000...........................Mr Gayno/DNXM
!     17 Mar 04  Removed duplicate "/" from filenames...Mr Lewiston/DNXM
!     14 Aug 08  Added blacklist code from AGRMET_precip................
!                ............................Chris Franks/2WXG/WEA(SEMS)
!     27 Jan 09  Removed lines setting p6 & p12 to missing since these are
!                now the same grid in both hemispheres.Chris Franks/2WG/WEA
!     23 Jun 09  Moved call of AGRMET_julhr_date10 to the beginning in
!                so that julhr is populated prior to being used.......
!                ............................Chris Franks/2WXG/WEA(SEMS)
!     11 Mar 10  Changed program names in messages to LIS.  Left argument
!                to LIS_alert with old program name to keep alert numbers
!                unique........................Chris Franks/16WS/WXE/SEMS
!      9 Sep 10  Modified for JMOBS............Chris Franks/16WS/WXE/SEMS
!     11 Mar 11  Enable use of JMOBS or CDMS...Chris Franks/16WS/WXE/SEMS
!     11 May 11  Improve processing for India and Sri Lanka
!                ..............................Chris Franks/16WS/WXE/SEMS
!      5 Sep 13  Add checks for negative indicies in grids where they
!                weren't included..............Chris Franks/16WS/WXE/SEMS
!
! !INTERFACE: 
subroutine AGRMET_processobs(n, obs, isize, stncnt, hemi, julhr, &
     imax, jmax, p6, pathpcp, p12,&
     use_twelve, cdms_flag, quad9r, alert_number, &
     precip6, precip12)

! !USES: 
  use AGRMET_forcingMod, only: agrmet_struc
  use LIS_coreMod, only : LIS_masterproc, LIS_domain
  use LIS_logMod, only : LIS_alert, LIS_logunit
  use LIS_mpiMod
  use USAF_bratsethMod, only: USAF_ObsData, USAF_assignObsData
  use map_utils, only: latlon_to_ij

  implicit none
! !ARGUMENTS:   
  integer,       intent(in)     :: n 
  character(len=*), intent(in)     :: pathpcp
  integer,  intent(in)          :: imax
  integer,  intent(in)          :: jmax
  real,  intent(inout)          :: p6(imax,jmax)
  real,  intent(inout)          :: p12(imax,jmax)
  real,  intent(in)             :: quad9r  
  integer,  intent(inout)       :: alert_number
  integer,  intent(in)          :: hemi
  integer,  intent(in)          :: isize
  integer,  intent(in)          :: julhr
  integer,  intent(in)          :: stncnt
  logical, intent(in)           :: use_twelve
  logical, intent(in)           :: cdms_flag
  type(USAF_ObsData), intent(inout) :: precip6
  type(USAF_ObsData), intent(inout) :: precip12
!
! !DESCRIPTION:
!    Processes rain gauge data then maps it to the model grid.
!    The idea behind the processing is to make as many 6 and 12 hourly
!    rain gauge amounts as possible using the information that
!    each station provides.
!
!    \textbf{Method} \newline
!    
!    - Convert date/time to julian hour
!    - Filter out duplicate obs (those with the same WMO number. \newline
!    - Read in observations from 6 hours ago (if available). \newline
!    - Read in observations from 12 hours ago (if available). \newline
!    - If the current 12 hourly amount is missing, try to calculate it
!      from current 6 hourly amts and 6 hourly amts from 6 six
!      hours ago. \newline
!    - If the above is not possible, then try to create the current
!      12 hourly amt from a current 24 hourly amt and a 12
!      hourly amt from 12 hours ago. \newline
!    - If the current 6 hourly amt is missing, then try to calculate it
!      using a current 12 hourly amt and a 6 hourly amt from six
!      hours ago. \newline
!    - If the station is from India/Bangledesh, then calculate 
!      current 12 and 6 hourly observations according to their own
!      renegade reporting practices. \newline
!    - Write the final array of current observations to an ascii
!      file for use the next time precip runs. \newline
!    - If this is an 06 or 18 utc run, put the current 6 hourly amts
!      on the model grid. \newline
!    - If this is an 00 or 12 utc run, put the current 12 hourly amts
!      on the model grid. \newline
!      - Data over S America and Russia require special handling. \newline
! 
! The arguments are variables are: 
! \begin{description}
!   \item[a,b,c,d]         Holds i/j coordinates of estimate/grid point in
!                     function sumsqr
!   \item[alert\_number]    Number of alert messages
!   \item[chemi]           Character representation of hemisphere 
!   \item[count]           Position of a filtered observation in the obs\_cur
!                     array
!   \item[count6]          Number of observations from 6 hours ago
!   \item[count6obs]       Number of 6 hourly rain gauge observations at
!                     the current time
!   \item[count12]         Number of observations from 12 hours ago
!   \item[count12obs]      Number of 12 hourly rain gauge observations at
!                     the current time
!   \item[count\_dup]       Number of duplicate observations
!   \item[date10]          Date/time group (yyyymmddhh) of current time
!   \item[date10\_min6]     Date/time group of current time minus six hours
!   \item[date10\_min12]    Date/time group of current time minus twelve hours
!   \item[dumhemi]         Dummy variable to hold hemisphere flag
!   \item[filename]        Name, including path, of the current observation 
!                     file
!   \item[filename\_min6]   Name, including path, of the observation file
!                     from six hours ago
!   \item[filename\_min12]  Name, including path, of the observation file
!                     from twelve hours ago
!   \item[firsts6]         Index to the first instance of each network in obs\_6
!   \item[firsts12]        Index to the first instance of each network in obs\_12
!   \item[fltrcnt]         Number of observations after duplicates are 
!                     filtered out
!   \item[hemi]            Hemisphere flag (1 - N, 2 - S)
!   \item[i]               Loop index
!   \item[index6]          Array index of an observation (from 6 hours ago)
!   \item[index12]         Array index of an observation (from 12 hours ago)
!   \item[imax]            Grid dimension - E/W direction
!   \item[iofunc]          I/O function - read or write. Used for diagnostic
!                         prints
!   \item[isize]           Maximum array size of observation arrays
!   \item[istat]           i/o status
!   \item[jmax]            Grid dimension - n/s direction
!   \item[julhr]           Julian hour of the current time
!   \item[julhr\_min6]      Julian hour of the current time minus 6 hrs
!   \item[julhr\_min12]     Julian hour of the current time minus 12 hrs
!   \item[lasts6]          Index to the last instance of each network in obs\_6
!   \item[lasts12]         Index to the last instance of each network in obs\_12
!   \item[message]         Alert message
!   \item[net\_count6]      The number of networks in obs\_6
!   \item[net\_count12]     The number of networks in obs\_12
!   \item[networks6]       List of networks in the obs\_6 list
!   \item[networks12]      List of networks in the obs\_12 list
!   \item[newd]            Distance in grid lengths between an observation
!                     and the nearest model grid point
!   \item[obs]             Array of unsorted observations (by WMO number)
!   \item[obs\_cur]         Array of sorted observations (by WMO number)
!                     for the current time
!   \item[obs\_6]           Array of observations from 6 hours ago
!   \item[obs\_12]          Array of observations from 12 hours ago
!   \item[oldd]            Holds distances between an observation and the
!                     nearest grid point
!   \item[pathpcp]         Directory path of precip data
!   \item[p6]              Six hourly precip amts on the model grid
!   \item[p12]             Twelve hourly precip amts on the model grid
!   \item[prior\_net]       The JMOBS of the previous ob in the obs\_6 or obs\_12 list
!   \item[quad9r]          The missing indicator (9999.)
!   \item[rain\_obs]        Abstract data type used to hold rain gauge obs.
!                     contains the following fields: \newline
!                     amt6   - six hourly precip amount \newline
!                     amt12  - twelve hourly precip amount \newline
!                     amt24  - twenty-four hourly precip amount \newline
!                     amtmsc - miscellaneous precip amount \newline
!                     lat    - latitude \newline
!                     lon    - longitude \newline
!                     wmonum - WMO block station number \newline
!   \item[ri]              I coordinate of an observation on the model grid
!   \item[rj]              J coordinate of an observation on the model grid
!   \item[sixyes]          Logical flag - does obs data exist from 
!                     six hours ago
!   \item[stncnt]          Number of observations with a precip report
!   \item[sumsqr]          Function to calculate the distance between an
!                     Estimate and the nearest grid point        
!   \item[twelveyes]       Logical flag - does obs data exist from 
!                     twelve hours ago
!   \item[use\_twelve]      Logical flag to put 12 hourly precip amts on 
!                     the model grid (set to true for 00 and 12 utc
!                     runs).  when false (at 06 and 18 utc), put 6
!                     hourly precip amts on the model grid 
!   \end{description}
!      
!     \textbf{Remarks} \newline
!     This routine has special logic to handle some countries' renegade
!     reporting practices (for example, India).  Therefore, the
!     maintainer should remain alert to changes in reporting practices
!     and modify the code as necessary.
!
!
!  The routines invoked are: 
!  \begin{description}
!  \item[AGRMET\_julhr\_date10](\ref{AGRMET_julhr_date10}) \newline
!    converts julian hour to a 10 character date string
!  \item[lltops](\ref{lltops}) \newline
!    converts lat lon values to points on the AGRMET grid
!  \item[pcpobs\_search](\ref{AGRMET_pcpobs_search}) \newline
!    finds a specific observation in an array of observations sorted
!    by WMO number.
!  \item[LIS\_alert](\ref{LIS_alert}) \newline
!   prints out an alert message
!  \end{description}
!EOP

  character*4                   :: chemi(2)
  character*10                  :: prior_net
  character*10                  :: networks6(25)
  character*10                  :: networks12(25)
  character*10                  :: date10
  character*10                  :: date10_min6
  character*10                  :: date10_min12
  character*120                 :: filename
  character*120                 :: filename_min6
  character*120                 :: filename_min12
  character*5                   :: iofunc
  character*100                 :: message(20)
  integer                       :: count
  integer                       :: count6
  integer                       :: count6obs
  integer                       :: count12
  integer                       :: count12obs
  integer                       :: count_dup
  integer                       :: dumhemi
  integer                       :: firsts6(25)
  integer                       :: firsts12(25)
  integer                       :: fltrcnt

  integer                       :: i

  integer                       :: index6
  integer                       :: index12

  integer                       :: india_lowlimit
  integer                       :: india_highlimit
  integer                       :: srilanka_lowlimit
  integer                       :: srilanka_highlimit

  integer                       :: istat
  integer                       :: julhr_min6
  integer                       :: julhr_min12
  integer                       :: lasts6(25)
  integer                       :: lasts12(25)
  integer                       :: net_count6
  integer                       :: net_count12

  integer                       :: russia_limit1
  integer                       :: russia_limit2
  integer                       :: russia_limit3
  integer                       :: russia_limit4
  integer                       :: s_amer_lowlimit
  integer                       :: s_amer_highlimit

  real                          :: a
  real                          :: b
  real                          :: c
  real                          :: d
  real                          :: newd
  real                          :: oldd(imax,jmax)
  real                          :: ri
  real                          :: rj
  real                          :: sumsqr
  real                          :: rlon
  
  logical                       :: sixyes      
  logical                       :: twelveyes

  integer                       :: ierr
  type rain_obs
     sequence
     character*10                :: net
     character*10                :: platform
     integer                     :: wmonum
     real                        :: lat
     real                        :: lon
     integer                     :: amt24
     integer                     :: amt12
     integer                     :: amt6
     integer                     :: amtmsc
  end type rain_obs

  type(rain_obs), intent(inout) :: obs(isize)
  type(rain_obs), allocatable   :: obs_cur(:)
  type(rain_obs), allocatable   :: obs_6(:)
  type(rain_obs), allocatable   :: obs_12(:)

  character*32 :: net32, platform32

  data chemi / '_nh.', '_sh.' /
  
  sumsqr(a,b,c,d) = ((a-b)**2) + ((c-d)**2)

#if (defined SPMD)
      call MPI_Barrier(LIS_mpi_comm, ierr)
#endif

  if (cdms_flag) then
     russia_limit1 = 200000
     russia_limit2 = 290009
     russia_limit3 = 320000
     russia_limit4 = 390009
     india_lowlimit = 420010
     india_highlimit = 433999
     srilanka_lowlimit = 434000
     srilanka_highlimit = 434979
     s_amer_lowlimit = 800000
     s_amer_highlimit = 889999
  else
     russia_limit1 = 20000
     russia_limit2 = 29000
     russia_limit3 = 32000
     russia_limit4 = 39000
     india_lowlimit = 42000
     india_highlimit = 43399
     srilanka_lowlimit = 43400
     srilanka_highlimit = 43497
     s_amer_lowlimit = 80000
     s_amer_highlimit = 889999
  end if

!     ------------------------------------------------------------------
!     Convert date/time to julian hour
!     ------------------------------------------------------------------

  call AGRMET_julhr_date10(julhr, date10)

!     ------------------------------------------------------------------
!     sort observation array based on WMO number using a selection
!     sort.
!     ------------------------------------------------------------------
  
  count_dup = 0
  
  SORT : do i = 1, stncnt-1
     
!     ------------------------------------------------------------------
!       keep track of duplicate observations - those with the same
!       WMO number. these are filtered from the observation array below.
!     ------------------------------------------------------------------
     
     if ( (i > 2) ) then 
        if(((obs(i)%net == obs(i-1)%net)) .and. &
           ((obs(i)%platform == obs(i-1)%platform))) then
           count_dup = count_dup + 1
        end if
     endif
     
  enddo SORT

!     ------------------------------------------------------------------
!     filter out duplicate observations and store in obs_cur array.
!     ------------------------------------------------------------------

  fltrcnt = stncnt - count_dup
  
  allocate ( obs_cur(fltrcnt) )
  
  obs_cur(1)%net    = obs(1)%net
  obs_cur(1)%platform = obs(1)%platform
  obs_cur(1)%wmonum = obs(1)%wmonum
  obs_cur(1)%lat    = obs(1)%lat
  obs_cur(1)%lon    = obs(1)%lon
  obs_cur(1)%amt24  = obs(1)%amt24
  obs_cur(1)%amt12  = obs(1)%amt12
  obs_cur(1)%amt6   = obs(1)%amt6
  obs_cur(1)%amtmsc = obs(1)%amtmsc
  
  count = 1
  
  FILTER : do i = 2, stncnt
     

     DUP_CHECK : if ((obs(i)%net == obs(i-1)%net) .and. &
                     (obs(i)%platform == obs(i-1)%platform)) then
        
!     ------------------------------------------------------------------
!         remove duplicate obs.  i think most duplicate obs are 
!         specials that are reported sometime during the past hour.
!         since there is no precip group associated with these 
!         asynoptic obs, database puts a zero in the six hourly
!         slot.  so, duplicate obs tend to look like the following,
!         where the amts are inches times 100:
!      
!                     24hr        12hr         6hr
!         1st dup     48          30          -99999999
!         2nd dup    -99999999   -99999999     0
!
!         the logic is written to select 48 (for 24hr), 30 (for 12hr)
!         and -99999999 (for 6hr).  because the zero is associated 
!         with an asynoptic ob, it is not selected.  i looked at
!         hundreds of examples over a 24 hr time period and selecting
!         the zero is wrong because the precip amts won't add up.
!         it is correct to leave the 6hrly as missing and let
!         the logic below calculate the 6 hourly amt using
!         data from the current and previous reports.
!         confused?  
!     ------------------------------------------------------------------
 
        if (obs(i)%amt24 > obs(i-1)%amt24) then
           obs_cur(count)%amt24 = obs(i)%amt24
        end if
        
        if (obs(i)%amt12 > obs(i-1)%amt12) then
           obs_cur(count)%amt12 = obs(i)%amt12
        end if
        
        if (obs(i)%amtmsc > obs(i-1)%amtmsc) then
           obs_cur(count)%amtmsc = obs(i)%amtmsc
        end if
        
        if ( (obs_cur(count)%amt24 > 0 .or. &
             obs_cur(count)%amt12 > 0) .and. &
             obs(i-1)%amt6 < -9999 .and. &
             obs(i)%amt6 == 0) then
           
           obs_cur(count)%amt6 = -99999999
           
        else if ( (obs_cur(count)%amt24 > 0 .or. &
             obs_cur(count)%amt12 > 0) .and. &
             obs(i)%amt6 < -9999 .and. &
             obs(i-1)%amt6 == 0) then

           obs_cur(count)%amt6 = -99999999
           
        else
           
           if (obs(i)%amt6 > obs(i-1)%amt6) then
              obs_cur(count)%amt6 = obs(i)%amt6
           end if
           
        end if
        
     else DUP_CHECK 

        count = count + 1
        
        obs_cur(count)%net    = obs(i)%net
        obs_cur(count)%platform = obs(i)%platform
        obs_cur(count)%wmonum = obs(i)%wmonum
        obs_cur(count)%lat    = obs(i)%lat
        obs_cur(count)%lon    = obs(i)%lon
        obs_cur(count)%amt24  = obs(i)%amt24
        obs_cur(count)%amt12  = obs(i)%amt12
        obs_cur(count)%amt6   = obs(i)%amt6
        obs_cur(count)%amtmsc = obs(i)%amtmsc
        
     end if DUP_CHECK
     
  enddo FILTER

!-----------------------------------------------------------------------
!     Adjust fltrcnt to reflect blacklisted obs.
!-----------------------------------------------------------------------

      fltrcnt = count

!-----------------------------------------------------------------------
!     now read in precip observations from 6 hours ago.  use this
!     data to calculate any missing 12 or 6 hourly amounts at the
!     current time.  if the file does not exist (if we are cold 
!     starting precip, for example) or is not read in properly, 
!     then set flag to false to turn off the calculations.  by the way,
!     if precip is ever cold started it will take a few cycles to
!     get the maximum number of rain gauge reports.
!-----------------------------------------------------------------------

  julhr_min6 = julhr - 6
  
  call AGRMET_julhr_date10(julhr_min6, date10_min6)
  
  filename_min6= trim(pathpcp) // "/presav" // chemi(hemi) // &
       "06hr." // date10_min6
  
  sixyes = .false.
  istat  = 0
  
  inquire (file=trim(filename_min6), exist=sixyes)
  
!-----------------------------------------------------------------------
!     read in precip data from 6 hours ago.
!-----------------------------------------------------------------------

  if (sixyes) then
     
     write(LIS_logunit,*)'[INFO] READING ', trim(filename_min6)
     
     iofunc = "OPEN "
     open(8, file=trim(filename_min6), iostat=istat, err=100)
     
     iofunc = "READ "
     read(8, *, iostat=istat, err=100, end=100) count6
     allocate(obs_6(count6))
     
     prior_net = "NULL"
     net_count6 = 0

     do i = 1, count6
        
        read(8, 6000, iostat=istat, err=100, end=100) &
             obs_6(i)%net, obs_6(i)%platform, obs_6(i)%lat, obs_6(i)%lon, &
             obs_6(i)%amt24, obs_6(i)%amt12, obs_6(i)%amt6,&
             obs_6(i)%amtmsc

        if ((obs_6(i)%net .eq. "WMO") .or. (obs_6(i)%net .eq. "CDMS")) then
           read(obs_6(i)%platform, FMT='(I10)') obs_6(i)%wmonum
        else
           obs_6(i)%wmonum = 0
        end if

        if (prior_net .eq. "NULL") then
           firsts6(1) = i
           networks6(1) = obs_6(i)%net
           prior_net = obs_6(i)%net
           net_count6 = 1
        else if (obs_6(i)%net .ne. prior_net) then
           lasts6(net_count6) = i - 1
           net_count6 = net_count6 + 1 
           firsts6(net_count6) = i
           networks6(net_count6) = obs_6(i)%net
           prior_net = obs_6(i)%net
        end if

     enddo
     lasts6(net_count6) = count6
     
     close (8)
     
  else
     
     write(LIS_logunit,*)' '
     write(LIS_logunit,*)"******************************************************"
     write(LIS_logunit,*)"[WARN] PRECIP SAVE FILE FROM 6 HOURS AGO DOES NOT EXIST."
     write(LIS_logunit,*)"[WARN] FILE NAME IS ", trim(filename_min6)
     write(LIS_logunit,*)"[WARN] OBSERVATION COUNT WILL BE REDUCED."
     write(LIS_logunit,*)"******************************************************"
     write(LIS_logunit,*)' '
     
  end if
  
100 continue

!-----------------------------------------------------------------------
!     a bad read on the 6 hour old data means we can't do some
!     calculations below.  precip can continue to run, but we 
!     won't have as many observations.  don't abort, just print
!     out a message.
!-----------------------------------------------------------------------

  if (istat /= 0) then
     
     write(LIS_logunit,*)' '
     write(LIS_logunit,*)"*******************************************************"
     write(LIS_logunit,*)"[WARN] BAD ", trim(iofunc), " ON FILE ",trim(filename_min6)
     write(LIS_logunit,*)"[WARN] ISTAT = ", istat
     write(LIS_logunit,*)"[WARN] OBSERVATION COUNT WILL BE REDUCED."
     write(LIS_logunit,*)"*******************************************************"
     write(LIS_logunit,*)' '
     
     message=' '
     message(1)='program:  LIS'
     message(2)='  routine:  AGRMET_processobs'
     message(3)='  bad ' // trim(iofunc) // ' on file'
        message(4)='  presav' // chemi(hemi) // '06hr.' // &
             date10_min6 // '.'
        message(5)='  observation count will be reduced.'
        alert_number = alert_number + 1
        if(LIS_masterproc) then 
           call lis_alert('precip              ',alert_number,message)
        endif
        sixyes = .false.

     end if
 
!----------------------------------------------------------------------------
!     now read in precip observations from 12 hours ago.  use this
!     data to calculate any missing 12 or 6 hourly amounts at the
!     current time - that is if the file exists and is read in properly.
!----------------------------------------------------------------------------

     julhr_min12 = julhr - 12
     
     call AGRMET_julhr_date10(julhr_min12, date10_min12)
     
     filename_min12= trim(pathpcp) // "/presav" // chemi(hemi) // &
          "06hr." // date10_min12
     
     twelveyes = .false.
     istat     = 0
     
     inquire (file=trim(filename_min12), exist=twelveyes)
     
     if (twelveyes) then
        
        iofunc = "OPEN "
        open(9, file=trim(filename_min12), iostat=istat, err=200)
        
        write(LIS_logunit,*)'[INFO] READING ', trim(filename_min12) 
        
        iofunc = "READ "
        read(9, *, iostat=istat, err=200, end=200) count12
        
        allocate(obs_12(count12))
        
        prior_net = "NULL"
        net_count12 = 0
 
        do i = 1, count12
           
           read (9, 6000, iostat=istat, err=200, end=200)  &
                obs_12(i)%net, obs_12(i)%platform, obs_12(i)%lat, &
                obs_12(i)%lon, obs_12(i)%amt24, &
                obs_12(i)%amt12, obs_12(i)%amt6, obs_12(i)%amtmsc

           if ((obs_12(i)%net .eq. "WMO") .or. (obs_12(i)%net .eq. "CDMS")) then
              read(obs_12(i)%platform, FMT='(I10)') obs_12(i)%wmonum
           else
              obs_12(i)%wmonum = 0
           end if

           if (prior_net .eq. "NULL") then
              firsts12(1) = i
              networks12(1) = obs_12(i)%net
              prior_net = obs_12(i)%net
              net_count12 = 1
           else if (obs_12(i)%net .ne. prior_net) then
              lasts12(net_count12) = i - 1
              net_count12 = net_count12 + 1 
              firsts12(net_count12) = i
              networks12(net_count12) = obs_12(i)%net
              prior_net = obs_12(i)%net
           end if
           
        enddo
        lasts12(net_count12) = count12
        
        close (9)
        
     else
        
        write(LIS_logunit,*)' '
        write(LIS_logunit,*)"*******************************************************"
        write(LIS_logunit,*)"[WARN] PRECIP SAVE FILE FROM 12 HOURS AGO DOES NOT EXIST."
        write(LIS_logunit,*)"[WARN] FILE NAME IS ", trim(filename_min12)
        write(LIS_logunit,*)"[WARN] OBSERVATION COUNT WILL BE REDUCED."
        write(LIS_logunit,*)"*******************************************************"
        write(LIS_logunit,*)' '
        
     end if

!-----------------------------------------------------------------------
!     a bad read on the 12 hour old data means we can't do some
!     calculations below.  precip can continue to run, but we 
!     won't have as many observations.  don't abort, just print
!     out a message.
!-----------------------------------------------------------------------

200  continue

     if (istat /= 0) then
        
        write(LIS_logunit,*)' '
        write(LIS_logunit,*)"*******************************************************"
        write(LIS_logunit,*)"[WARN] BAD ", trim(iofunc), " ON FILE ",trim(filename_min12)
        write(LIS_logunit,*)"[WARN] ISTAT = ", istat
        write(LIS_logunit,*)"[WARN] OBSERVATION COUNT WILL BE REDUCED."
        write(LIS_logunit,*)"*******************************************************"
        write(LIS_logunit,*)' '

        message=' '
        message(1)='program:  LIS'
        message(2)='  routine:  AGRMET_processobs'
        message(3)='  bad ' // trim(iofunc) // ' on file'
        message(4)='  presav' // chemi(hemi) // '06hr.' // &
             date10_min12 // '.'
        message(5)='  observation count will be reduced.'
        alert_number = alert_number + 1
        if(LIS_masterproc) then 
           call lis_alert('precip              ',alert_number,message)
        endif
        twelveyes = .false.

      end if

!-----------------------------------------------------------------------
!     if there is no 12 hourly amount for the current time, 
!     try to calculate it using a current 6 hourly amount and
!     a 6-hourly amount from 6 hours ago.  
!-----------------------------------------------------------------------

#if (defined SPMD)
      call MPI_Barrier(LIS_mpi_comm, ierr)
#endif
      
      if (sixyes) then
         
         write(LIS_logunit,*)' '
         write(LIS_logunit,*)'[INFO] MAKING 12 HRLY AMTS FROM PREVIOUS 6 HRLY AMTS'
         
         MAKE12_FROM_6 : do i = 1, fltrcnt
            
!     ------------------------------------------------------------------
!         India and Sri Lanka are a special case and are handled
!         in a separate loop below. 
!     ------------------------------------------------------------------

            if (obs_cur(i)%wmonum >= india_lowlimit .and.  &
                 obs_cur(i)%wmonum <= srilanka_highlimit) cycle MAKE12_FROM_6

            if (obs_cur(i)%amt12 <= -99999998) then
               
               if (obs_cur(i)%amt6 > -99999998) then
                  
                  call AGRMET_pcpobs_search(obs_cur(i)%net, obs_cur(i)%platform, &
                                            networks6, firsts6, lasts6, net_count6, &
                                            count6, obs_6, index6)

!     ------------------------------------------------------------------
!             an index of -9999 means a station did not report 
!             6 hours ago.
!     ------------------------------------------------------------------

                  if (index6 /= -9999) then
                     
                     if (obs_6(index6)%amt6 > -99999998) then
                        
                        obs_cur(i)%amt12 = obs_cur(i)%amt6 + &
                             obs_6(index6)%amt6
                        if (obs_cur(i)%amt12 < 0) then
                           obs_cur(i)%amt12 = 0
                        end if
                        
                     end if
                     
                  end if
                  
               end if
               
            end if
            
         enddo MAKE12_FROM_6
         
      end if

!-----------------------------------------------------------------------
!     if there is no 12 hourly amount for the current time, 
!     try to calculate it using a current 24 hourly amount and
!     a 12-hourly amount from 12 hours ago.  
!-----------------------------------------------------------------------

#if (defined SPMD)
      call MPI_Barrier(LIS_mpi_comm, ierr)
#endif

      if (twelveyes) then
         
         write(LIS_logunit,*)'[INFO] MAKING 12 HRLY AMTS FROM CURRENT 24 AND OLD 12 HRLY AM  TS'

         MAKE12_FROM_24 : do i = 1, fltrcnt

!     ------------------------------------------------------------------
!         India and Sri Lanka are a special case and are handled
!         in a separate loop below. 
!     ------------------------------------------------------------------

            if (obs_cur(i)%wmonum >= india_lowlimit .and.  &
                 obs_cur(i)%wmonum <= srilanka_highlimit) cycle MAKE12_FROM_24

            if (obs_cur(i)%amt12 <= -99999998) then
               
               if (obs_cur(i)%amt24 > -99999998) then
                  
                  call AGRMET_pcpobs_search(obs_cur(i)%net, obs_cur(i)%platform, &
                                            networks12, firsts12, lasts12, net_count12, &
                                            count12, obs_12, index12)

!     ------------------------------------------------------------------
!             an index of -9999 means a station did not report 
!             12 hours ago.
!     ------------------------------------------------------------------

                  if (index12 /= -9999) then
                     
                     if (obs_12(index12)%amt12 > -99999998) then
                        
                        obs_cur(i)%amt12 = obs_cur(i)%amt24 - &
                             obs_12(index12)%amt12
                        if (obs_cur(i)%amt12 < 0) then
                           obs_cur(i)%amt12 = 0
                        end if
                        
                     end if
                     
                  end if
                  
               end if
               
            end if
            
         enddo MAKE12_FROM_24
         
      end if

!-----------------------------------------------------------------------
!     if there is no 6 hourly amount for the current time, 
!     try to calculate it using a current 12 hourly amount and
!     a 6-hourly amount from 6 hours ago.  
!-----------------------------------------------------------------------

#if (defined SPMD)
      call MPI_Barrier(LIS_mpi_comm, ierr)
#endif

      if (sixyes) then
         
         MAKE6_FROM_12 : do i = 1, fltrcnt

!     ------------------------------------------------------------------
!         India and Sri Lanka are a special case and are handled
!         in a separate loop below. 
!     ------------------------------------------------------------------

            if (obs_cur(i)%wmonum >= india_lowlimit .and.  &
                 obs_cur(i)%wmonum <= srilanka_highlimit) cycle MAKE6_FROM_12
            
            if (obs_cur(i)%amt6   <= -99999998 .and. &
                 obs_cur(i)%amt12  >  -99999998) then

               call AGRMET_pcpobs_search(obs_cur(i)%net, obs_cur(i)%platform, &
                                         networks6, firsts6, lasts6, net_count6, &
                                         count6, obs_6, index6)

!     ------------------------------------------------------------------
!           an index of -9999 means a station did not report 
!           6 hours ago.
!     ------------------------------------------------------------------

               if (index6 == -9999) cycle MAKE6_FROM_12
               
               if (obs_6(index6)%amt6 > -99999998) then
                  
                  obs_cur(i)%amt6  = obs_cur(i)%amt12 - &
                       obs_6(index6)%amt6
                  if (obs_cur(i)%amt6 < 0) then
                     obs_cur(i)%amt6 = 0
                  end if

               end if
               
            end if
    
         enddo MAKE6_FROM_12
         
      end if

!     ------------------------------------------------------------------
!     handle India and Sri Lanka, which have their own goofy 
!     reporting standards.
!     ------------------------------------------------------------------

#if (defined SPMD)
      call MPI_Barrier(LIS_mpi_comm, ierr)
#endif

      INDIA : do i = 1, fltrcnt

!     ------------------------------------------------------------------
!       process obs with block station number from India
!       Find any ob from 6 hours ago for this station
!     ------------------------------------------------------------------

         if (obs_cur(i)%wmonum >= india_lowlimit .and.  &
              obs_cur(i)%wmonum <= india_highlimit) then

            if (sixyes) then
               call AGRMET_pcpobs_search(obs_cur(i)%net, obs_cur(i)%platform, &
                                         networks6, firsts6, lasts6, net_count6, &
                                         count6, obs_6, index6)
            end if

!         -----------------------------------------------------------------
!         India reports precip by accumulating it from 
!         03z.  Any accumulation will be placed in the amtmsc slot.
!         If amtmsc is missing set it to 0.  Set amt12 to the the amt6 from
!         6 hours ago if present, else set amt12 to 0.
!         -----------------------------------------------------------------

            if (obs_cur(i)%amtmsc <= -99999998) then

               obs_cur(i)%amtmsc = 0

               if ((sixyes) .and. (index6 /= -9999)) then
                  obs_cur(i)%amt12 = obs_6(index6)%amtmsc
               else
                  obs_cur(i)%amt12 = 0
               end if

!         ---------------------------------------------------------------
!         If there's a valid amtmsc and amt6 is missing
!              if there's an ob from 6 hours ago with a valid amtmsc the 
!              current amt6 is the current amtmsc less the amtmsc 
!              from 6 hours ago
!                 if there's a valid amt6 from 6 hours ago set amt12 to 
!                 the sum of the current & 6-hour old amt6
!         ---------------------------------------------------------------
               
            elseif ( obs_cur(i)%amtmsc >  -99999998 .and. &
                 obs_cur(i)%amt6   <= -99999998 .and. &
                 sixyes )                              then
               
               if (index6 /= -9999) then
                  
                  if ((obs_6(index6)%amtmsc > -99999998) .and. &
                      (obs_6(index6)%amtmsc <= obs_cur(i)%amtmsc)) then
                     
                     obs_cur(i)%amt6 = &
                          obs_cur(i)%amtmsc - obs_6(index6)%amtmsc

                  else

                     obs_cur(i)%amt6 = 0

                  end if

                  if (obs_6(index6)%amt6 > -99999998) then
                        
                     obs_cur(i)%amt12 = obs_cur(i)%amt6 +  &
                          obs_6(index6)%amt6
                     if (obs_cur(i)%amt12 < 0) obs_cur(i)%amt12 = 0
                        
                  else

                     obs_cur(i)%amt12 = obs_cur(i)%amt6

                  end if
                     
!-----------------------------------------------------------------------
!              Presume first precipitation of the 03 to 03Z day
!-----------------------------------------------------------------------

               else

                  obs_cur(i)%amt6 = obs_cur(i)%amtmsc 
                  obs_cur(i)%amt12 = obs_cur(i)%amtmsc 
                  
               end if  ! index6 is valid

!-----------------------------------------------------------------------
!           If amt6 has been set before this module calculate amt12
!-----------------------------------------------------------------------

            elseif ( obs_cur(i)%amt6   > -99999998 ) then
               
               if ((sixyes) .and. (index6 /= -9999) .and. &
                   (obs_6(index6)%amt6 > -99999998)) then

                   obs_cur(i)%amt12 = obs_cur(i)%amt6 +  &
                          obs_6(index6)%amt6

               else

                   obs_cur(i)%amt12 = obs_cur(i)%amt6

               end if

            end if

!-----------------------------------------------------------------------
!           At 6Z we need to carry only the 3-6Z accumultion forward
!           as amtmsc.  The 24-hour accumulation from 3z the previous
!           day to 3z today was saved in amt24, subtract this from amtmsc
!-----------------------------------------------------------------------

            if ( date10(9:10)      == "06"      .and. &
                 obs_cur(i)%amtmsc >  0 .and. &
                 obs_cur(i)%amt24  >  0 ) then
               obs_cur(i)%amtmsc = obs_cur(i)%amtmsc - obs_cur(i)%amt24
               if (obs_cur(i)%amtmsc < 0) obs_cur(i)%amtmsc = 0
            end if
            
         end if  ! WMO number is India
         
!     ------------------------------------------------------------------
!       Process obs with block station number from Sri Lanka
!       Search for an ob for this station 6 hours ago
!     ------------------------------------------------------------------

         if (obs_cur(i)%wmonum >= srilanka_lowlimit .and.  &
              obs_cur(i)%wmonum <= srilanka_highlimit) then
            if (sixyes) then
               call AGRMET_pcpobs_search(obs_cur(i)%net, obs_cur(i)%platform, &
                                         networks6, firsts6, lasts6, net_count6, &
                                         count6, obs_6, index6)
            end if

!     ------------------------------------------------------------------
!           If there's no precip set amtmsc to 0
!              If there's an amount from 6 hours ago set amt12 to that
!     ------------------------------------------------------------------

            if (obs_cur(i)%amtmsc <= -99999998 .and. &
                 obs_cur(i)%amt6  == 0) then
               obs_cur(i)%amtmsc = 0

               if ((sixyes) .and. (index6 /= -9999)) then
                  obs_cur(i)%amt12 = obs_6(index6)%amtmsc
               else
                  obs_cur(i)%amt12 = 0
               end if

!     ------------------------------------------------------------------
!           Else there's a valid amtmsc, set amt6 to that
!              If there's an ob from 6 hours ago set amt12 to current
!              amt6 plus that from 6 hours ago
!              Else amt12 equals amt6
!     ------------------------------------------------------------------

            else
               obs_cur(i)%amt6 = obs_cur(i)%amtmsc 
               if ((sixyes) .and. (index6 /= -9999)) then
                  obs_cur(i)%amt12 = obs_cur(i)%amt6 +  obs_6(index6)%amt6
               else
                  obs_cur(i)%amt12 = obs_cur(i)%amt6
               end if
            end if
         end if
      enddo INDIA

!-----------------------------------------------------------------------
!     print out new and improved data, only by the master processor
!-----------------------------------------------------------------------

      if(LIS_masterproc) then 
         count6obs  = 0
         count12obs = 0
      
         filename = trim(pathpcp) // "/presav" // chemi(hemi) // &
              "06hr." // date10
         
         istat = 0
         
         iofunc = "OPEN"
         open (11, file=trim(filename), iostat=istat, err=300)
         
         write(LIS_logunit,*)' '
         write(LIS_logunit,*)"[INFO] WRITING ",trim(filename)
         
         iofunc = "WRITE"
         write(11,*, iostat=istat, err=300) fltrcnt
         
         do i = 1, fltrcnt
            
            write(11, 6000, iostat=istat, err=300) &
                 obs_cur(i)%net, obs_cur(i)%platform, &
                 obs_cur(i)%lat, obs_cur(i)%lon,&
                 obs_cur(i)%amt24, &
                 obs_cur(i)%amt12, obs_cur(i)%amt6, &
                 obs_cur(i)%amtmsc
            
            if (obs_cur(i)%amt12 > -99999998) count12obs = count12obs + 1
            if (obs_cur(i)%amt6  > -99999998) count6obs  = count6obs + 1
         
         enddo
         
         close (11)
         
         write(LIS_logunit,*)' '
         write(LIS_logunit,*)'[INFO] NUMBER OF 12 AND 6 HOURLY OBS IS ',count12obs, count6obs
         write(LIS_logunit,*)' '
      
300      continue
!-----------------------------------------------------------------------
!     send alert message on bad write.
!-----------------------------------------------------------------------

         if (istat /= 0) then
            
            write(LIS_logunit,*)' '
            write(LIS_logunit,*)"*******************************************************"
            write(LIS_logunit,*)"[WARN] BAD ", trim(iofunc), " ON FILE ",trim(filename)
            write(LIS_logunit,*)"[WARN] ISTAT = ", istat
            write(LIS_logunit,*)"*******************************************************"
            write(LIS_logunit,*)' '
            
            message=' '
            message(1)='program:  LIS'
            message(2)='  routine:  AGRMET_processobs'
            message(3)='  bad ' // trim(iofunc) // ' on file'
            message(4)='  presav' // chemi(hemi) // '06hr.' // date10 // '.'
            message(5)='  precip analysis will be degraded next cycle.'
            alert_number = alert_number + 1
            call lis_alert('precip              ',alert_number,message)
            
         end if
      endif
      
#if (defined SPMD)
      call MPI_Barrier(LIS_mpi_comm, ierr)
#endif

!-----------------------------------------------------------------------
!     if the run until time is 06 or 18 utc, will use all 6 hrly 
!     reports. place the 6-hourly precip observations on the agrmet grid
!     using a nearest neighbor approach.  if more than one observation
!     is within one grid length of a particular grid point, use the
!     closest one.  this is achieved by calculating the distance 
!     between the observation and grid point using function sumsqr
!     and storing it in array oldd.  first, initialize oldd to a
!     very large number so that a "first" observation will always
!     be placed at a grid point.  note: russian obs only report
!     every 12 hours, so no need to process them here.
!-----------------------------------------------------------------------

      USE_6 : if ( .not. use_twelve ) then
         
         write(LIS_logunit,*)"[INFO] PUTTING 6 HOURLY RAIN GAUGE OBSERVATIONS ON GRID."
         
         oldd = 999.0
         
         OBS_2_GRID6 : do i = 1, fltrcnt
            
            if ( obs_cur(i)%wmonum >= russia_limit1 .and. &
                 obs_cur(i)%wmonum  < russia_limit4 ) cycle OBS_2_GRID6
            
            if (obs_cur(i)%amt6 >= 0 .and. obs_cur(i)%amt6 <= 1016) then
               
!-----------------------------------------------------------------------
!           in rare cases (like when an observation is located at the
!           equator) the hemisphere argument passed to lltops will
!           come out of the routine as the other hemisphere.  this can
!           cause problems because variable hemi is used as a loop
!           counter.  to avoid this potential problem, send in a
!           dummy argument for hemisphere.
!-----------------------------------------------------------------------

               dumhemi = hemi

               call latlon_to_ij(LIS_domain(n)%lisproj, &
                    obs_cur(i)%lat, obs_cur(i)%lon,ri,rj)

               ! EMK...Add observation
               net32 = obs_cur(i)%net
               platform32 = obs_cur(i)%platform
               call USAF_assignObsData(precip6, &
                    net32, platform32, &
                    float(obs_cur(i)%amt6) * 0.1, &
                    obs_cur(i)%lat, obs_cur(i)%lon,&
                    agrmet_struc(n)%bratseth_precip_gauge_sigma_o_sqr, &
                    0.)

!-----------------------------------------------------------------------
!           calculate distance between the current observation and 
!           grid point.  if it is closer than a previous observation, 
!           use the current observation.  then convert units from
!           mm times 10 to mm.  then store the new distance in
!           oldd.
!-----------------------------------------------------------------------

            newd = sumsqr(ri, float(nint(ri)), rj, float(nint(rj)))  &
                 * 100.0
            
            if(nint(ri).ge.1.and. nint(ri).le.imax.and.&
                 nint(rj).ge.1.and.nint(rj).le.jmax) then                

               if ( newd < oldd(nint(ri),nint(rj)) ) then
                  
                  p6(nint(ri),nint(rj))   = float(obs_cur(i)%amt6) * 0.1
                  
                  oldd(nint(ri),nint(rj)) = newd
                  
               end if
            endif
         end if
         
      enddo OBS_2_GRID6
      
   end if USE_6
   
!-----------------------------------------------------------------------
!     if the run until time is 00 or 12 utc, will use all 12 hrly 
!     reports. place the 12-hourly precip observations on the agrmet grid
!     using a nearest neighbor approach.  if more than one observation
!     is within one grid length of a particular grid point, use the
!     closest one.  this is achieved by calculating the distance 
!     between the observation and grid point using function sumsqr
!     and storing it in array oldd.  first, initialize oldd to a
!     very large number so that a "first" observation will always
!     be placed at a grid point. note: russia is a special case
!     that is handled below.
!-----------------------------------------------------------------------

   USE_12 : if (use_twelve) then

      write(LIS_logunit,*)"[INFO] PUTTING 12 HOURLY RAIN GAUGE OBSERVATIONS ON GRID."
      
      oldd = 999.0
      
      OBS_2_GRID12 : do i = 1, fltrcnt
         
         if ( obs_cur(i)%wmonum >= russia_limit1 .and. &
              obs_cur(i)%wmonum  < russia_limit4 ) cycle OBS_2_GRID12
         
         if (obs_cur(i)%amt12 >= 0 .and. obs_cur(i)%amt12 <= 1016) then
            
            dumhemi = hemi

            call latlon_to_ij(LIS_domain(n)%lisproj, obs_cur(i)%lat, &
                 obs_cur(i)%lon,ri,rj)

            ! EMK...Add observation
            net32 = obs_cur(i)%net
            platform32 = obs_cur(i)%platform
            call USAF_assignObsData(precip12, &
                 net32, platform32, &
                 float(obs_cur(i)%amt12) * 0.1, &
                 obs_cur(i)%lat, obs_cur(i)%lon, &
                 agrmet_struc(n)%bratseth_precip_gauge_sigma_o_sqr, &
                 0.)
            
            newd = sumsqr(ri, float(nint(ri)), rj, float(nint(rj)))  &
                 * 100.0
            
            if(nint(ri).ge.1.and.nint(ri).le.imax.and.&
                 nint(rj).ge.1.and.nint(rj).le.jmax) then 

               if ( newd < oldd(nint(ri),nint(rj)) ) then
  
!-----------------------------------------------------------------------
!             convert from mm times 10 to mm.
!-----------------------------------------------------------------------

                  p12(nint(ri),nint(rj))  = float(obs_cur(i)%amt12) * 0.1
               
                  oldd(nint(ri),nint(rj)) = newd
               
               end if
            endif
         end if
         
      enddo OBS_2_GRID12

!-----------------------------------------------------------------------
!       special case: south america at 12 utc.  s america reports are
!       as follows: 6 hrly amounts at 00, 06 and 18 utc; 24 hourly
!       amounts at 12 utc.  because the database stores "zeros" in
!       the six hourly slot (i think this is because when an location
!       has no precip, it does not report a precip group in the ob.
!       as a result, you know there was no precip, but you don't know
!       the duration of that "zero."  conversely, when there is precip,
!       the duration must be reported in the ob.  therefore, database
!       shop is able to put a "non-zero" in the proper database slot.)
!       because there are relatively few obs at 06z (guess they are
!       sleeping), the logic above will make very few 12 hrly zeroes
!       at 12z.  therefore, i assume the "zeros" in the 6 hrly slot are
!       are really 12 hrly amounts.  i.e, the duration of the zero
!       is the last time they reported, which was 12 hours ago.
!       i spent a week looking at s amer reporting practices and this
!       is a pretty good assumption.  you don't have this problem with
!       obs that reported precip at 12z (a 24 hrly amt) because it
!       is easy to make a 12 hrly non-zero amount from the 6 hourly
!       data at 18 (the previous day) and 00 utc.  confused?
!-----------------------------------------------------------------------
      
      if ( date10(9:10) == "12") then
         
         SAMER_12Z : do i = 1, fltrcnt
            
            if (obs_cur(i)%amt12 < -99999 .and.  &
                 obs_cur(i)%amt6 == 0 .and. &
                 obs_cur(i)%wmonum >= 80000 .and. &
                 obs_cur(i)%wmonum <= 88999) then

              dumhemi = hemi
              
              call latlon_to_ij(LIS_domain(n)%lisproj, obs_cur(i)%lat, &
                   obs_cur(i)%lon,ri,rj)

              ! EMK...Add observation
              net32 = obs_cur(i)%net
              platform32 = obs_cur(i)%platform
              call USAF_assignObsData(precip12, &
                   net32, platform32, &
                   0.0, &
                   obs_cur(i)%lat, obs_cur(i)%lon, &
                   agrmet_struc(n)%bratseth_precip_gauge_sigma_o_sqr, &
                   0.)

              newd = sumsqr(ri, float(nint(ri)), rj, float(nint(rj)))  &
                   * 100.0

              if(nint(ri).ge.1.and.nint(ri).le.imax.and.&
                 nint(rj).ge.1.and.nint(rj).le.jmax) then 
                 
                 if ( newd < oldd(nint(ri),nint(rj)) ) then
                    
                    p12(nint(ri),nint(rj))   = 0.0
                    
                    oldd(nint(ri),nint(rj)) = newd
                    
                 end if
              endif
           end if
           
        enddo SAMER_12Z
        
     end if
     
!-----------------------------------------------------------------------
!       russia reports every 12 hours and does not typically report
!       a valid duration flag.  therefore, use what was in the misc
!       array.
!-----------------------------------------------------------------------

     if ( (date10(9:10) == "12" .or.  &
          date10(9:10) == "00") .and. hemi == 1 ) then

        RUSSIA_0012: do i = 1, fltrcnt
           
           if ( obs_cur(i)%wmonum >= russia_limit1  .and. &
                obs_cur(i)%wmonum  < russia_limit4 ) then

              if (obs_cur(i)%amtmsc >= 0 .and. &
                   obs_cur(i)%amtmsc <  1016) then

                 dumhemi = hemi
 
                 call latlon_to_ij(LIS_domain(n)%lisproj, &
                      obs_cur(i)%lat, obs_cur(i)%lon,&
                      ri,rj)

                 ! EMK...Add observation
                 net32 = obs_cur(i)%net
                 platform32 = obs_cur(i)%platform
                 call USAF_assignObsData(precip12, &
                      net32, platform32, &
                      float(obs_cur(i)%amtmsc) * 0.1, &
                      obs_cur(i)%lat, obs_cur(i)%lon, &
                      agrmet_struc(n)%bratseth_precip_gauge_sigma_o_sqr, &
                      0.)

                newd = sumsqr(ri, float(nint(ri)), rj, float(nint(rj)))  &
                     * 100.0

                if(nint(ri).ge.1.and.nint(ri).le.imax.and.&
                   nint(rj).ge.1.and.nint(rj).le.jmax) then 

                   if ( newd < oldd(nint(ri),nint(rj)) ) then 
                      
                      p12(nint(ri),nint(rj))  =  &
                           float(obs_cur(i)%amtmsc) * 0.1
                      
                      oldd(nint(ri),nint(rj)) = newd
                      
                   end if
                endif
             end if
             
          end if
          
       enddo RUSSIA_0012

!----------------------------------------------------------------------
!         large parts of russia only report 12 hourly amts at 06 and
!         18z.  therefore, you can't make a 12 hourly amount at 00
!         or 12z.  because there are not as many estimates at this 
!         latitude as there are in the tropics, just assume that
!         a 12 hourly report ending at 18/06z really ended at 00/12z.
!         the previous agrmet decoder made this same assumption.
!----------------------------------------------------------------------

       if (sixyes) then
                    
          RUSSIA_0618 : do i = 1, count6
             
             if ( (obs_6(i)%wmonum >= russia_limit1 .and. &
                  obs_6(i)%wmonum  < russia_limit2) .or. &
                  (obs_6(i)%wmonum >= russia_limit3 .and. &
                  obs_6(i)%wmonum  < russia_limit4) ) then

                if ( obs_6(i)%amtmsc >= 0 .and. &
                     obs_6(i)%amtmsc <= 1016 ) then 
                   
                   dumhemi = hemi
                   
                   call latlon_to_ij(LIS_domain(n)%lisproj, &
                        obs_6(i)%lat, obs_6(i)%lon,&
                        ri,rj)

                   ! EMK...Add observation
                   net32 = obs_6(i)%net
                   platform32 = obs_6(i)%platform
                   call USAF_assignObsData(precip12, &
                        net32, platform32, &
                        float(obs_6(i)%amtmsc) * 0.1, &
                        obs_6(i)%lat, obs_6(i)%lon, &
                        agrmet_struc(n)%bratseth_precip_gauge_sigma_o_sqr, &
                        0.)
                   
                   newd = sumsqr(ri, float(nint(ri)), &
                        rj, float(nint(rj))) * 100.0

                   if(nint(ri).ge.1.and.nint(ri).le.imax.and.&
                      nint(rj).ge.1.and.nint(rj).le.jmax) then 

                      if ( newd < oldd(nint(ri),nint(rj)) ) then
                      
!-----------------------------------------------------------------------
!                   convert from mm times 10 to mm.
!-----------------------------------------------------------------------

                         p12(nint(ri),nint(rj)) =  &
                              float(obs_6(i)%amtmsc) * 0.1
                               
                         oldd(nint(ri),nint(rj)) = newd
                      
                      end if
                   endif
                end if
                
             end if
             
          enddo RUSSIA_0618
          
       end if ! if sixyes
       
    end if ! if (n hemis) and (00 or 12)
    
 end if USE_12
 
 
 deallocate (obs_cur)
 
 if ( allocated(obs_6) )  deallocate(obs_6)
 if ( allocated(obs_12) ) deallocate(obs_12)
 
 return 

!-----------------------------------------------------------------------
!     format statements.
!-----------------------------------------------------------------------

! EMK 15 Apr 2022...Cray compiler does not like the original format statement.
!6000 format(1x, a9, 1x, a9, 1x, f7.2, 1x, f7.2, 1x, ' 24hr ', 1x, & 
!            i9, 1x, '12hr ', i9, ' 6hr ', i9, ' misc ', i9)
6000 format(1x, a9, 1x, a9, 1x, f7.2, 1x, f7.2, 8x, i9, 6x, i9, 5x, i9, &
          6x, i9)
 
end subroutine AGRMET_processobs

