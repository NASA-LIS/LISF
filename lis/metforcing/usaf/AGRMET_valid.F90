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
! !ROUTINE: AGRMET_valid
!  \label{AGRMET_valid}
!
! !REVISION HISTORY:
!
!
!     15 jun 96  initial version..................mr moore/sysm(agromet)
!     10 apr 97  brought up to software standards.  modified to create
!                supergrid of new source files.  also now begins summing
!                12-hr estimate and real quad9r values while in box
!                arrays for greater cpu and memory efficiency...........
!                .......................................ssgt miller/sysm
!     02 may 97  changed array indices to prevent 'thrashing' of memory.
!                added check of source value to ensure validity of est.
!                updated prolog and brought up to standards.............
!                .......................................capt andrus/sysm
!      7 oct 99  ported to ibm sp-2, updated prolog, incorporated 
!                FORTRAN 90 features.................capt hidalgo/agrmet
!     05 jan 01  renamed variables to reflect 6-hourly change to
!                6-hourly cycles...........................mr gayno/dnxm
!     10 jun 02  modified to reflect change from rtneph to cdfsii.......
!                ..........................................mr gayno/dnxm
!      3 nov 05 Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine AGRMET_valid( n, pcap, mrg, est, est6, src,&
     rel, rel6, cdfsii3, cdfsii6, srcwts, julhr,&
     pathpcp)
! !USES: 
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_logunit,LIS_endrun, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_fileIOMod, only : LIS_readData
  use LIS_timeMgrMod,only : LIS_julhr_date

  implicit none

! !ARGUMENTS:       

  integer,       intent(in)     :: n
  integer                       :: src(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real,          intent(inout)  :: cdfsii3(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real,          intent(out)    :: cdfsii6(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real                          :: est(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real,          intent(out)    :: est6(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real,          intent(out)    :: mrg(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real                          :: rel(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real,          intent(out)    :: rel6(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real,          intent(in)     :: pcap
  real,          intent(in)     :: srcwts(8)
  integer,       intent(in)     :: julhr
  character(len=*), intent(in)    ::  pathpcp
!  
! !DESCRIPTION:
!
!     to validate the 3 hrly estimated and phased real precip
!     amounts and set (or sum) into the 6 hrly values.  to
!     set the mrgd value.
!
!    \textbf{Method} \newline
!    
!     1. loop thru points in the hemisphere \newline
!       a. if 3 hrly estimated and/or real value passes validity checks \newline
!          1) put (or sum) the value into 6 hrly value \newline
!          2) put 3 hrly real value into 3 hrly mrg value \newline
!       b. else \newline
!          1) put 9999.0 into 3 hrly estimated and/or real value \newline
!          2) put 3 hrly estimated value into mrg value \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[cdfsii3]      3hr cdfsii based precip estimate (mm/3hr)
!   \item[cdfsii6]      6hr cdfsii based precip estimate (mm/6hr)
!   \item[hemi]         current hemisphere ( 1 = nh, 2= sh ) 
!   \item[est]          estimated 3hr precip array (mm/3hr)
!   \item[est6]         estimated 6hr precip array (mm/6hr)
!   \item[i]            loop index, point's i-coordinate
!   \item[ifil]	   character array holding path/name of input file(s)
!   \item[imax]         number of gridpoints in east-west direction
!   \item[j]            loop index, point's j-coordinate
!   \item[jmax]         number of gridpoints in north-south direction
!   \item[rel]          parsed 3hr real precip array (mm/3hr)
!   \item[rel6]         parsed 6hr real precip array (mm/6hr)   
!   \item[hemi]         hemisphere (1=north, 2=south)
!   \item[mrg]          merged 3hr precip array (mm/s = kg/(m2 s)) 
!   \item[pathpcp]      character array for path to precip files
!   \item[pcap]         precip cap limit (mm/3hr)
!   \item[src]          precip data source array
!   \item[srcwts]       value weight placed on each source
!  \end{description}
!EOP

  character*100                :: ifil
  integer                      :: i
  integer                      :: j
  character*10                 :: date10_03
  integer                      :: yr1,mo1,da1,hr1
  character*4                  :: fyr
  character*2                  :: fmo,fda
  logical                      :: exists
  integer                      :: ftn
  real                         :: gridDesc(6)
  real                         :: src_tmp(LIS_rc%lnc(n), LIS_rc%lnr(n))

  call AGRMET_julhr_date10( julhr, date10_03 )

!     ------------------------------------------------------------------
!     loop thru points in the hemisphere
!     ------------------------------------------------------------------

  do j = 1, LIS_rc%lnr(n)
     do i = 1, LIS_rc%lnc(n)
        
!     ------------------------------------------------------------------
!         determine if 3-hrly estimate has a valid source.
!         if not, the estimate is invalid.
!     ------------------------------------------------------------------

        if( (src(i,j) .lt. 2) .or. (src(i,j) .gt. 6) )then
           est(i,j) = LIS_rc%udef
        endif

!     ------------------------------------------------------------------
!         check estimate for gross error limit validity and set (or sum)
!         into 6-hrly estimate.  if 3-hrly estimate exceeds pcap
!         but is not invalid, set it to invalid.
!     ------------------------------------------------------------------
        
        if( est(i,j) .le. pcap .and. est(i,j).ge.0.0)then
           if( est6(i,j) .lt.-9990.0 )then
              est6(i,j) = est(i,j)
           else
              est6(i,j) = est6(i,j) + est(i,j)
           endif
        elseif( est(i,j) .gt. -9990.0 )then
           write(LIS_logunit,6100) i, j, est(i,j), pcap
           est(i,j) = LIS_rc%udef
        endif

!     ------------------------------------------------------------------
!         check cdfsii estimate for gross error limit validity and set
!         (sum) into 6-hrly estimate.  if 3-hrly estimate exceeds pcap
!         but is not invalid, set it to invalid.
!     ------------------------------------------------------------------
        
        if( cdfsii3(i,j) .le. pcap .and. cdfsii3(i,j).ge.0.0)then
           if( cdfsii6(i,j) .lt. -9990.0 )then
              cdfsii6(i,j) = cdfsii3(i,j)
           else
              cdfsii6(i,j) = cdfsii6(i,j) + cdfsii3(i,j)
           endif
        elseif( cdfsii3(i,j) .gt. -9990.0 )then
           write(LIS_logunit,6200) i, j, cdfsii3(i,j), pcap
           cdfsii3(i,j) = LIS_rc%udef
        endif

!     ------------------------------------------------------------------
!         check parsed real amt validity. set mrg value accordingly
!         and identify data source as reals. set (or sum) the
!         3 hrly real amt into 6 hrly real amt.
!     ------------------------------------------------------------------

        if( (rel(i,j) .gt. -9990.0).and.(srcwts(1) .gt. 0.0) )then
           if( rel(i,j) .le. pcap .and. rel(i,j).ge.0.0)then
              mrg(i,j) = rel(i,j)
              src(i,j) = 1
              if( rel6(i,j) .lt.-9990.0 )then
                 rel6(i,j) = rel(i,j)
              else
                 rel6(i,j) = rel6(i,j) + rel(i,j)
              endif
           else
              write(LIS_logunit,6000)  i, j, rel(i,j), pcap
              rel(i,j) = LIS_rc%udef
              if(est(i,j).gt.-9990) then 
                 mrg(i,j) = est(i,j)
              else
                 mrg(i,j) = -9999.0
              endif
           endif
        else
           if(est(i,j).gt.-9990) then 
              mrg(i,j) = est(i,j)
           else
              mrg(i,j) = -9999.0
           endif           
        endif
        
     enddo
  enddo
  return

!     ------------------------------------------------------------------
!     format statements
!     ------------------------------------------------------------------

6000 format (/, 1x, 55('-'), &
          /, 3x, 'ROUTINE VALID:  RELP VALUE EXCEEDS PCAP',&
          /, 5x, 'I/J = ',2x, i4, 2x, i4,&
          /, 5x, 'VALUE & CAP = ', f12.1, 2x, f12.1,&
          /, 1x, 55('-'))

6100 format (/, 1x, 55('-'),&
          /, 3x, 'ROUTINE VALID:  ESTP VALUE EXCEEDS PCAP',&
          /, 5x, 'I/J = ',2x, i4, 2x, i4,&
          /, 5x, 'VALUE & CAP = ', f12.1, 2x, f12.1,&
          /, 1x, 55('-'))

 6200 format (/, 1x, 55('-'),&
           /, 3x, 'ROUTINE VALID:  CDFSII VALUE EXCEEDS PCAP',&
           /, 5x, 'I/J = ',2x, i4, 2x, i4,&
           /, 5x, 'VALUE & CAP = ', f12.1, 2x, f12.1,&
           /, 1x, 55('-'))
end subroutine AGRMET_valid
  
