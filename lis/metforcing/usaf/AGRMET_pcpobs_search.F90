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
! !ROUTINE: AGRMET_pcpobs_search
! \label{AGRMET_pcpobs_search}
!
! !REVISION HISTORY: 
!    21 feb 01   initial version...........................mr gayno/dnxm
!     9 Sep 10   modified for JMOBS...........Chris Franks/16WS/WXE/SEMS
!
! !INTERFACE: 
subroutine AGRMET_pcpobs_search(network, plat_id, networks, firsts, lasts, &
                                net_count, isize, obs, index)

  implicit none
! !ARGUMENTS:   

  character*10, intent(in)     :: network
  character*10, intent(in)     :: networks(25)
  character*10, intent(in)     :: plat_id
  integer,  intent(in)         :: firsts(25)
  integer,  intent(out)        :: index
  integer,  intent(in)         :: isize
  integer,  intent(in)         :: lasts(25)
  integer,  intent(in)         :: net_count
  
! !DESCRIPTION: 
!    finds a specific observation in an array of observations sorted
!    by network and platform ID.
!    
!    \textbf{Method} \newline
!    - search for an observation using a binary search. \newline
!    - if found, pass back its array index. \newline
!    - if not found, pass back an array index of -9999 (missing). \newline
!    
! The arguments are variables are: 
! \begin{description}
!   \item[network]    JMOBS network of the obs, i.e. WMO, ICAO, FAA, etc.
!   \item[plat\_id]    JMOBS platform ID, WMO number, call letters, etc.
!   \item[net\_count]  The number of networks represented in the obs array
!   \item[networks]   List of networks represented in the obs array
!   \item[firsts]     Array of the first instance of each network in the obs array
!   \item[lasts]      Array of the last instance of each network in the obs array
!   \item[index]      array index of the observation we were searching for.
!               if not found, set to -9999
!   \item[isize]      number of elements in the observation array
!   \item[obs]        array of observations sorted by network \& plat\_id \newline
!      net      JMOBS network 
!      platform JMOBS platform ID
!      amt6     six hourly precip amount \newline
!      amt12    twelve hourly precip amount \newline
!      amt24    twenty-four hourly precip amount \newline
!      lat      latitude \newline
!      lon      longitude \newline
!      wmonum   wmo block station number \newline
! Other variables
!   \item[first]      array index of the beginning of the search region
!   \item[found]      logical flag - true/false - observation found
!   \item[last]       array index of the end of the search region
!   \item[middle]     array index of the middle of the search region
!  \end{description}
!EOP
  integer                      :: last
  integer                      :: middle
  integer                      :: net_index
  integer                      :: first
  logical                      :: found
  
  type rain_obs
     sequence
     character*10               :: net
     character*10               :: platform
     integer                    :: wmonum
     real                       :: lat
     real                       :: lon
     integer                    :: amt24
     integer                    :: amt12
     integer                    :: amt6
     integer                    :: amtmsc
  end type rain_obs
  
  type(rain_obs), intent(in)   :: obs(isize)
  
!-----------------------------------------------------------------------
!     begin binary search here.  initialize array index to -9999 -
!     the "not found" indicator.
!     First find the range where the observations from this network 
!     reside
!-----------------------------------------------------------------------

  index  =  -9999
  first  =  1
  last   =  isize
  found  =  .false.

  net_index = 1

  do while ((net_index <= net_count) .and. (.not. found))

     if (network .eq. networks(net_index)) then
        first = firsts(net_index)
        last = lasts(net_index)
        found = .true.
     end if
     net_index = net_index + 1

  end do
  
  if (found) then

     found = .false.

     do
     
        if ( (first > last) .or. found) return

!-----------------------------------------------------------------------
!       cut the search area in half every time through loop - the
!       definition of a binary search.
!-----------------------------------------------------------------------

        middle = (first + last) / 2
     
        if (plat_id < obs(middle)%platform) then
        
           last = middle  - 1
        
        elseif (plat_id > obs(middle)%platform) then
        
           first = middle + 1
        
        else
           
           found = .true.
           index = middle
        
        end if
        
     enddo

  end if
     
  return
     
end subroutine AGRMET_pcpobs_search

