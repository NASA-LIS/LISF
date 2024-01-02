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
! MODULE: USAFSI_analysisMod
!
! REVISION HISTORY:
! 08 Feb 2019  Eric Kemp  First ported to LDT.
! 09 May 2019  Eric Kemp  Renamed LDTSI
! 13 Dec 2019  Eric Kemp  Renamed USAFSI
! 02 Nov 2020  Eric Kemp  Removed blacklist code at request of 557WW.
! 22 Jan 2021  Yeosang Yoon Add subroutine for new 0.1 deg snow climatology
! 13 Jan 2022  Eric Kemp Added support for GRIB1 FNMOC SST file.
! 28 Jul 2023  Eric Kemp Added support for new sfcobs file format (longer
!              station names.
!
! DESCRIPTION:
! Source code for Air Force snow depth analysis.
!-------------------------------------------------------------------------

#include "LDT_misc.h"

module USAFSI_analysisMod

   ! Defaults
   implicit none
   private

   ! Public methods
   public :: find_nearest_valid_value ! EMK
   public :: getfrac
   public :: getgeo
   public :: getobs
   public :: getsfc
   public :: getsmi
   public :: getsno
   public :: getsno_nc
   public :: getsst
   public :: getviirs
   public :: run_snow_analysis_noglacier ! EMK
   public :: run_snow_analysis_glacier ! EMK
   public :: run_seaice_analysis_ssmis ! EMK
   public :: run_seaice_analysis_gofs  ! EMK
   public :: getclimo                  ! Yeosang Yoon
 
   ! Internal constant
   real, parameter :: FILL = -1

contains
   
   ! Private subroutine
   subroutine appclm (pntclm, pntold, pntanl, pntage)

      !***********************************************************************
      !***********************************************************************
      !**
      !**  NAME: APPLY SNOW DEPTH CLIMATOLOGY
      !**
      !**  PURPOSE: APPLIES CLIMATOLOGY TO ALL REMAINING UNANALYZED POINTS.
      !**
      !**  CALLED FROM: PROCES
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  PNTAGE, PNTCLM, PNTOLD
      !**
      !**   OUTPUT: PNTAGE, PNTANL
      !**
      !**  FILES ACCESSED: NONE
      !**
      !**  SEE USAFSI_paramsMod.F90 FOR COMMON PARAMETER DESCRIPTIONS
      !**
      !**  UPDATES
      !**  =======
      !**  26 APR 95  INITIAL UNISYS VERSION......................DR KOPP/SYSM
      !**  23 OCT 97  ADDED INCLUDE OF TUNES PROC AND USE OF ITS CLMADJ
      !**             VARIABLE.................................SSGT CONRY/DNXM
      !**  22 FEB 01  INITIAL VERSION OF CODE PORTED FROM UNISYS MAINFRAME....
      !**             ........................................SSGT MILLER/DNXM
      !**  21 JUL 04  CONVERTED TO FORTRAN 90 FOR 16TH MESH...MR EYLANDER/DNXM
      !**  14 APR 05  ADDED SNOTHRESH TO REPLACE HARCODED VALUE.
      !**             CREATED UNIQUE NAMES FOR VARIABLES HOLDING ONE ELEMENT
      !**             OF ARRAY INSTEAD OF SAME NAME AS ARRAY..MR LEWISTON/DNXM
      !**  27 AUG 08  REPLACED AMIN1 FUNCTION WITH MIN....MR LEWISTON/16WS/WXE
      !**  21 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**
      !***********************************************************************
      !***********************************************************************

      ! Imports
      use LDT_usafsiMod, only: usafsi_settings
      use USAFSI_paramsMod

      ! Defaults
      implicit none

      ! Arguments
      real,      intent(in)           :: pntclm ! POINT SNOW CLIMATOLOGY IN METERS
      real,      intent(in)           :: pntold ! PREVIOUS POINT SNOW ANALYSIS IN METERS
      real,      intent(out)          :: pntanl ! CURRENT POINT SNOW ANALYSIS IN METERS
      integer,   intent(inout)        :: pntage ! POINT SNOW AGE IN DAYS

      ! Local variables
      real                            :: adjust ! ADJUSTMENT TO PREVIOUS SNOW ANALYSIS

      
      ! IF SNOW IS PRESENT, DETERMINE THE DIFFERENCE BETWEEN YESTERDAY'S
      ! ANALYSIS AND CLIMATOLOGY.  MULTIPLY BY CLMADJ (SET IN TUNING FILE)
      ! AND ADJUST THE ANALYSIS.  IF YESTERDAY'S ANALYSIS = 0, LEAVE AT 0.
      ! INCREMENT AGE FOR LOCATION WITH SNOW.  IF ADJUSTED DEPTH IS BELOW
      ! MINIMUM FOR CLIMO (SET IN TUNING FILE), SET DEPTH AND AGE TO ZERO.
      if (pntold > 0.0) then

         ! PREVIOUS DAY HAD SNOW COVER, ADJUST TOWARDS CLIMATOLOGY.
         ! INCREMENT THE AGE ONE DAY.
         adjust = (pntclm - pntold) * (usafsi_settings%clmadj)
         pntanl = pntold + adjust
         pntage = pntage + 1

         if (pntanl < snothresh) then
            ! ADJUSTED SNOW TOTAL BELOW THRESHOLD.  SET DEPTH AND AGE TO 0.
            pntanl = 0.0
            pntage = 0
         endif

         ! LIMIT AGE TO A MAXIMUM OF MAXAGE NUMBER OF DAYS.
         pntage = min(pntage,maxage)

      else

         ! THERE IS NO SNOW COVER.  SET DEPTH AND AGE TO 0.
         pntanl = 0.0
         pntage = 0

      endif

      return

   end subroutine appclm

   ! EMK Search for nearest valid point.
   ! Algorithm mimics subroutine search_extrap in GEOGRID, but code is
   ! simpler.
   function find_nearest_valid_value(nc,nr,data,missing,first_c,first_r) &
        result(answer)

      ! Defaults
      implicit none

      ! Arguments
      integer,intent(in) :: nc
      integer,intent(in) :: nr
      real, intent(in) :: data(nc,nr)
      real, intent(in) :: missing
      integer, intent(in) :: first_c
      integer, intent(in) :: first_r

      ! Result
      real :: answer

      ! Local variables
      integer :: c,r,i,j
      integer :: maxi
      integer, allocatable :: c_list(:)
      integer, allocatable :: r_list(:)
      logical, allocatable :: bitmap(:,:)
      logical :: found_valid
      real :: distance, new_distance

      ! Allocations
      found_valid = .false.
      maxi = nc*nr
      allocate(c_list(maxi))
      c_list(:) = 0
      allocate(r_list(maxi))
      r_list(:) = 0
      allocate(bitmap(nc,nr))
      bitmap(:,:) = .false.

      ! Initializations
      answer = missing
      c = first_c
      r = first_r
      i = 0
      j = 1
      c_list(1) = c
      r_list(1) = r
      bitmap(c,r) = .true.

      ! Search the grid surrounding the point in question for valid data,
      ! keeping track of what has already been searched, and keeping track of
      ! new locations to search around.
      do while (.not. found_valid .and. i .lt. maxi)
         i = i + 1
         c = c_list(i)
         r = r_list(i)

         if (data(c,r) .ne. missing) then
            found_valid = .true.
         end if
         if ( (c-1) .ge. 1 ) then
            if (j .lt. maxi .and. .not. bitmap(c-1,r)) then
               bitmap(c-1,r) = .true.
               j = j + 1
               c_list(j) = c-1
               r_list(j) = r
            end if
         end if
         if ( (c+1) .le. nc) then
            if (j .lt. maxi .and. .not. bitmap(c+1,r)) then
               bitmap(c+1,r) = .true.
               j = j + 1
               c_list(j) = c+1
               r_list(j) = r
            end if
         end if
         if ( (r-1) .ge. 1) then
            if (j .lt. maxi .and. .not. bitmap(c,r-1)) then
               bitmap(c,r-1) = .true.
               j = j + 1
               c_list(j) = c
               r_list(j) = r-1
            end if
         end if
         if ( (r+1) .le. nr) then
            if (j .lt. maxi .and. .not. bitmap(c,r+1)) then
               bitmap(c,r+1) = .true.
               j = j + 1
               c_list(j) = c
               r_list(j) = r+1
            end if
         end if         
      end do 

      if (found_valid) then
         distance = ((c-first_c)*(c-first_c)) + ((r-first_r)*(r-first_r))
         answer = data(c,r) ! First guess
         ! Search the remaining points we stored above for a closer answer.
         do while (i .lt. maxi .and. i .lt. j) 
            i = i + 1
            c = c_list(i)
            r = r_list(i)
            new_distance = &
                 ((c-first_c)*(c-first_c)) + ((r-first_r)*(r-first_r))
            if (new_distance < distance) then
               if (data(c,r) .ne. missing) then
                  distance = new_distance
                  answer = data(c,r)
               end if
            end if
         end do
      end if

      ! Clean up
      deallocate(c_list)
      deallocate(r_list)
      deallocate(bitmap)

   end function find_nearest_valid_value

   subroutine getfrac (date10, fracdir)

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME: GET FRACTIONAL SNOW
      !**
      !**  PURPOSE: READ IN FRACTIONAL SNOW DATA
      !**
      !**  CALLED FROM: SNODEP
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  DATE10, FRACDIR
      !**
      !**   OUTPUT: SNOFRAC, USEFRAC
      !**
      !**  FILES ACCESSED
      !**  ==============
      !**  FILE NAME                                R/W DESCRIPTION
      !**  ---------------------------------------- --- ------------------------------
      !**  ${FRACDIR}/snofrac_0p05deg.${DATEFR}.dat  R  FRACTIONAL SNOW
      !**
      !**  UPDATES
      !**  =======
      !**  07 NOV 12  INITIAL VERSION.............................MR LEWISTON/16WS/WXE
      !**  15 NOV 13  ADDED SEARCH BACK IF CURRENT DAY NOT FOUND..MR LEWISTON/16WS/WXE
      !**  21 Mar 19  Adapted to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
      !*******************************************************************************
      !*******************************************************************************

      ! Imports
      use LDT_coreMod, only: LDT_rc, LDT_domain
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      use LDT_usafsiMod, only: usafsi_settings
      use map_utils
      use USAFSI_arraysMod, only: USAFSI_arrays
      use USAFSI_paramsMod
      use USAFSI_utilMod

      ! Defaults
      implicit none

      ! Arguments
      character*10,  intent(in)   :: date10           ! DATE-TIME GROUP OF USAFSI CYCLE
      character*255, intent(in)   :: fracdir          ! FRACTIONAL SNOW DIRECTORY PATH
      
      ! Local constants
      character*8, parameter :: meshnp05 = '_0p05deg' ! MESH FOR 1/20 DEGREE FILE NAME

      integer,     parameter :: meshf  = 20           ! ECE GRID MESH (20 = 1/20 DEGREE)
      integer,     parameter :: igridf = 360 * meshf  ! SIZE OF GRID IN THE I-DIRECTION
      integer,     parameter :: jgridf = 180 * meshf  ! SIZE OF GRID IN THE J-DIRECTION

      ! Local variables
      character*2                 :: cyclhr           ! CYCLE HOUR

      character*10                :: datefr           ! DATE-TIME GROUP OF FRACTIONAL SNOW
      character*255               :: file_path        ! FULLY-QUALIFIED FILE NAME
      character*7                 :: iofunc           ! ACTION TO BE PERFORMED
      character*90                :: message (msglns) ! ERROR MESSAGE
      character*12                :: routine_name     ! NAME OF THIS SUBROUTINE
      integer                     :: fracnt           ! NUMBER OF FRACTIONAL POINTS
      integer                     :: i                ! SNODEP I-COORDINATE
      integer                     :: icount           ! LOOP COUNTER
      integer                     :: julhr            ! AFWA JULIAN HOUR
      integer                     :: maxdays          ! NUMBER OF DAYS TO SEARCH BACK
      integer,       allocatable  :: pntcnt  (: , :)  ! COUNT OF POINTS WITH DATA
      integer                     :: j                ! SNODEP J-COORDINATE
      logical                     :: isfile           ! FLAG INDICATING WHETHER FILE FOUND
      real,          allocatable  :: infrac_0p05deg  (: , :)  ! CDFS II FRACTIONAL SNOW DATA
      real,          allocatable  :: snocum  (: , :)  ! FRACTIONAL SNOW ACCUMULATOR
      integer :: i_0p05deg, j_0p05deg
      real :: rlat,rlon,ri,rj

      data routine_name  / 'GETFRAC     ' /

      ! ALLOCATE DATA ARRAYS.
      allocate (infrac_0p05deg (igridf, jgridf))
      allocate(pntcnt( ldt_rc%lnc(1), ldt_rc%lnr(1)))
      allocate(snocum( ldt_rc%lnc(1), ldt_rc%lnr(1)))

      ! INITIALIZE VARIABLES.
      fracnt  = 0
      icount  = 1
      iofunc  = 'READING'
      isfile  = .false.
      maxdays = 3
      message = ' '
      pntcnt  = 0
      snocum  = 0.0
      USAFSI_arrays%snofrac = -1

      ! RETRIEVE FRACTIONAL SNOW DATA.
      cyclhr = date10 (9:10)
      datefr = date10 (1:8) // '09'

      if (cyclhr .eq. '00' .or. cyclhr .eq. '06') then
         call date10_julhr (datefr, julhr, program_name, routine_name)
         julhr  = julhr  - 24
         call julhr_date10 (julhr, datefr, program_name, routine_name)
      end if

      file_search : do while ((.not. isfile) .and.(icount .le. maxdays))

         file_path = trim (fracdir) // 'snofrc' // meshnp05 // '.'       &
              // datefr // '.dat'

         inquire (file=file_path, exist=isfile)

         if (isfile) then
            write(ldt_logunit,*)'[INFO] Reading ', trim(file_path)
            write (ldt_logunit, 6000) routine_name, iofunc, file_path
            call putget_real (infrac_0p05deg, 'r', file_path, &
                 program_name,       &
                 routine_name, igridf, jgridf)

         else

            call date10_julhr (datefr, julhr, program_name, routine_name)
            julhr  = julhr  - 24
            icount = icount + 1
            call julhr_date10 (julhr, datefr, program_name, routine_name)

         end if

      end do file_search

      if (isfile) then

         ! EMK New version.  Average snow onto LDT grid.  
         do j_0p05deg = 1, jgridf
            rlat = -89.975 + (j_0p05deg-1)*0.05
            do i_0p05deg = 1, igridf
               if (.not. infrac_0p05deg(i_0p05deg, j_0p05deg) > 0) cycle
               rlon = -179.975 + (i_0p05deg-1)*0.05
               call latlon_to_ij(LDT_domain(1)%ldtproj,rlat,rlon,ri,rj)
               i = nint(ri)
               if (i .lt. 1) then
                  i = i + LDT_rc%lnc(1)
               else if (i .gt. LDT_rc%lnc(1)) then
                  i = i - LDT_rc%lnc(1)
               end if
               j = nint(rj)
               if (j .lt. 1) then
                  j = 1
               else if (j .gt. LDT_rc%lnr(1)) then
                  j = LDT_rc%lnr(1)
               end if
               pntcnt(i,j) = pntcnt(i,j) + 1
               snocum(i,j) = snocum(i,j) + infrac_0p05deg(i_0p05deg,j_0p05deg)
            end do ! i_0p05deg
         end do ! j_0p05deg

         do j = 1,LDT_rc%lnr(1)
            do i = 1, LDT_rc%lnc(1)
               if (pntcnt(i,j) > 0) then
                  fracnt = fracnt + 1
                  USAFSI_arrays%snofrac(i,j) = snocum(i,j) / real(pntcnt(i,j))
               end if
            end do ! i
         end do ! j

         write (ldt_logunit, 6200) routine_name, fracnt

      else

         usafsi_settings%usefrac = .false.
         message(1) = '[WARN]  FRACTIONAL SNOW FILE NOT FOUND'
         message(2) = '[WARN]  PATH = ' // trim(file_path)
         call error_message (program_name, routine_name, message)
         write (LDT_logunit, 6400) routine_name, file_path

      end if

      ! DEALLOCATE ARRAYS.
      deallocate (infrac_0p05deg)
      deallocate (pntcnt)
      deallocate (snocum)

      return

      ! FORMAT STATEMENTS.
6000  format (/, '[INFO] ', A7, ': ', A7, 1X, A)

6200  format (/, '[INFO]', A7, ': POINTS WITH FRACTIONAL SNOW = ', I6)

6400  format (/, '[WARN]', A7, ': FILE NOT FOUND: ', A,                    &
           /, '   WILL NOT USE FRACTIONAL SNOW DATA')

   end subroutine getfrac

   subroutine getgeo (month, static, nc, nr, elevations)

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME: GET GEOGRAPHY
      !**
      !**  PURPOSE: READ IN SNODEP STATIC DATA SETS
      !**
      !**  CALLED FROM: SNODEP
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  MONTH, STATIC
      !**
      !**   OUTPUT: ELEVAT, CLIMO, PTLAT, PTLON, SNOW_POSS
      !**
      !**  FILES ACCESSED
      !**  ==============
      !**  FILE NAME                            R/W DESCRIPTION
      !**  ------------------------------------ --- ------------------------------
      !**  ${STATIC}/snoclimo_0p25deg_jan.dat    R  SNOW CLIMATOLOGY FOR JANUARY
      !**  ${STATIC}/snoclimo_0p25deg_feb.dat    R  SNOW CLIMATOLOGY FOR FEBRUARY
      !**  ${STATIC}/snoclimo_0p25deg_mar.dat    R  SNOW CLIMATOLOGY FOR MARCH
      !**  ${STATIC}/snoclimo_0p25deg_apr.dat    R  SNOW CLIMATOLOGY FOR APRIL
      !**  ${STATIC}/snoclimo_0p25deg_may.dat    R  SNOW CLIMATOLOGY FOR MAY
      !**  ${STATIC}/snoclimo_0p25deg_jun.dat    R  SNOW CLIMATOLOGY FOR JUNE
      !**  ${STATIC}/snoclimo_0p25deg_jul.dat    R  SNOW CLIMATOLOGY FOR JULY
      !**  ${STATIC}/snoclimo_0p25deg_aug.dat    R  SNOW CLIMATOLOGY FOR AUGUST
      !**  ${STATIC}/snoclimo_0p25deg_sep.dat    R  SNOW CLIMATOLOGY FOR SEPTEMBER
      !**  ${STATIC}/snoclimo_0p25deg_oct.dat    R  SNOW CLIMATOLOGY FOR OCTOBER
      !**  ${STATIC}/snoclimo_0p25deg_nov.dat    R  SNOW CLIMATOLOGY FOR NOVEMBER
      !**  ${STATIC}/snoclimo_0p25deg_dec.dat    R  SNOW CLIMATOLOGY FOR DECEMBER
      !**  ${STATIC}/snow_mask_0p25deg.dat       R  MASK INDICATING SNOW POSSIBLE
      !**
      !**   See USAFSI_arraysMod.F90 for common array descriptions
      !**   See USAFSI_paramsMod.F90 for common parameter descriptions
      !**
      !**  UPDATES
      !**  =======
      !**  26 APR 95  INITIAL UNISYS VERSION...........................SSGT CONRY/SYSM
      !**  22 FEB 01  PORTED FROM UNISYS MAINFRAME TO UNIX............SSGT MILLER/DNXM
      !**  18 APR 03  UPGRADED TO 16TH MESH: CHANGED GEOGRAPHY RECORD LENGTH
      !**             FROM 4-BYTE TO 2-BYTE WORD AND TERRAIN TO INTEGER............AER
      !**  21 JUL O4  CONVERTED TO FORTRAN 90, MOVED MOST PASSED ARRAYS INTO
      !**             MODULES & COMPLETED 16TH MESH UPGRADE...........MR EYLANDER/DNXM
      !**  28 APR 05  ADDED PARAMETERS FOR RECORD LENGTH..............MR LEWISTON/DNXM
      !**  18 JUN 09  CONVERTED TO EQUIDISTANT CYL, ADDED AMSR-E..MR LEWISTON/2WXG/WEA
      !**  24 FEB 10  ADDED CLIMO, INPUT PATH, AND FILENAMES.
      !**             REPLACED FILE COMMANDS WITH PUTGET CALLS....MR LEWISTON/16WS/WXE
      !**  29 DEC 10  UPDATED FILENAMES...........................MR LEWISTON/16WS/WXE
      !**  08 NOV 11  REMOVED AMSR-E FILES AND PORTED TO LINUX....MR LEWISTON/16WS/WXE
      !**  04 APR 12  CHANGED CLIMO TO SEPARATE MONTHLY FILES.....MR LEWISTON/16WS/WXE
      !**  31 MAY 12  CORRECTED ORDER OF PUTGET ARGUMENTS.........MR LEWISTON/16WS/WXE
      !**  11 SEP 12  REMOVED SNOW_DENSITY........................MR LEWISTON/16WS/WXE
      !**  14 DEC 12  REMOVED OBSOLETE ERROR HANDLING FOR CLIMO...MR LEWISTON/16WS/WXE
      !**  21 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
      !**
      !*******************************************************************************
      !*******************************************************************************

      ! Imports
      use LDT_coreMod, only: LDT_domain, LDT_rc
      use LDT_logMod, only: LDT_endrun, ldt_logunit
      use map_utils
      use USAFSI_arraysMod, only: USAFSI_arrays
      use USAFSI_paramsMod
      use USAFSI_utilMod ! EMK

      ! Defaults
      implicit none

      ! Arguments
      integer,       intent(in)   :: month            ! MONTH OF YEAR (1-12)
      character*255, intent(in)   :: static           ! STATIC FILE DIRECTORY PATH
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: elevations(nc,nr)
      
      ! Local variables
      character*4                 :: cmonth  (12)     ! MONTH OF YEAR
      character*4                 :: file_ext         ! LAST PORTION OF FILE NAME
      character*255               :: file_path        ! FULLY-QUALIFIED FILE NAME
      character*90                :: message (msglns) ! ERROR MESSAGE
      character*12                :: routine_name     ! NAME OF THIS SUBROUTINE
      real, allocatable :: climo_0p25deg(:,:)
      integer*1, allocatable :: snow_poss_0p25deg(:,:)
      type(proj_info) :: snodep_0p25deg_proj
      integer :: i_0p25deg, j_0p25deg
      integer :: gindex,c,r
      real :: rlat,rlon,ri,rj

      data cmonth        / '_jan', '_feb', '_mar', '_apr', '_may', '_jun', &
           '_jul', '_aug', '_sep', '_oct', '_nov', '_dec' /

      data file_ext      / '.dat' /

      data routine_name  / 'GETGEO      ' /

      call map_set(proj_code=proj_latlon, &
           lat1=begin_lat, &
           lon1=begin_lon, &
           dx=0.25, &
           stdlon=0.25, &
           truelat1=0.25, &
           truelat2=0., &
           idim=igrid,&
           jdim=jgrid, &
           proj=snodep_0p25deg_proj)

      ! ALLOCATE ARRAYS.
      allocate(climo_0p25deg(igrid,jgrid))
      allocate(snow_poss_0p25deg(igrid,jgrid))

      ! INITIALIZE ERROR HANDLER VARIABLES.
      MESSAGE = ' '

      ! RETRIEVE THE CLIMATOLOGY FOR THE MONTH.
      ! THE CLIMO FILE CONTAINS AN ARRAY FOR EACH OF THE 12 MONTHS.
      ! EACH MONTH IS STORED CONSECUTIVELY STARTING WITH JANUARY.
      FILE_PATH = TRIM(STATIC) // 'snoclimo' // MESHNAME //             &
           CMONTH(MONTH) // FILE_EXT
      write(ldt_logunit,*)'[INFO] Reading ', trim(file_path)
      CALL PUTGET_REAL (CLIMO_0p25deg, 'r', FILE_PATH, PROGRAM_NAME,    &
           ROUTINE_NAME, IGRID, JGRID)

      USAFSI_arrays%climo(:,:) = -1
      do r = 1,nr
         do c = 1,nc
            gindex = c+(r-1)*nc
            rlat = LDT_domain(1)%lat(gindex)
            rlon = LDT_domain(1)%lon(gindex)
            call latlon_to_ij(snodep_0p25deg_proj,rlat,rlon,ri,rj)
            i_0p25deg = nint(ri)
            if (i_0p25deg .gt. igrid) then
               i_0p25deg = i_0p25deg - igrid
            else if (i_0p25deg .lt. 1) then
               i_0p25deg = i_0p25deg + igrid
            end if
            j_0p25deg = nint(rj)
            if (j_0p25deg .lt. 1) then
               j_0p25deg = 1
            else if (j_0p25deg .gt. jgrid) then
               j_0p25deg = jgrid
            end if
            if (climo_0p25deg(i_0p25deg,j_0p25deg) < -1) then
               USAFSI_arrays%climo(c,r) = -1
            else
               USAFSI_arrays%climo(c,r) = climo_0p25deg(i_0p25deg,j_0p25deg)
            end if
         end do ! c
      end do ! r

      ! EMK Copy LDT latitude data to PTLAT array
      do r = 1,nr
         do c = 1 ,nc
            gindex = c+(r-1)*LDT_rc%lnc(1)
            USAFSI_arrays%ptlat(c,r) = LDT_domain(1)%lat(gindex)
         end do ! i
      end do ! j

      ! EMK Copy LDT longitude data to PTLON array
      do r = 1,nr
         do c = 1 ,nc
            gindex = c+(r-1)*LDT_rc%lnc(1)
            USAFSI_arrays%ptlon(c,r) = LDT_domain(1)%lon(gindex)
         end do ! i
      end do ! j

      ! EMK Copy LDT terrain data to elevat array
      do r = 1,nr
         do c = 1 ,nc
            USAFSI_arrays%elevat(c,r) = elevations(c,r)
         end do ! i
      end do ! j

      ! RETRIEVE SNOW MASK DATA.
      file_path = trim(static) // 'snow_mask' // meshname // file_ext
      write(ldt_logunit,*)'[INFO] Reading ', trim(file_path)
      call putget_int1 (snow_poss_0p25deg, 'r', file_path, program_name,     &
           routine_name, igrid, jgrid)

      ! Interpolate the "snow possible" mask to the LDT grid.  
      ! For simplicity, just use the value of the 0.25 deg grid box that the
      ! LDT grid point is within.
      USAFSI_arrays%snow_poss(:,:) = 0
      do r = 1,nr
         do c = 1,nc
            gindex = c+(r-1)*LDT_rc%lnc(1)
            rlat = LDT_domain(1)%lat(gindex)
            rlon = LDT_domain(1)%lon(gindex)
            call latlon_to_ij(snodep_0p25deg_proj,rlat,rlon,ri,rj)
            i_0p25deg = nint(ri)
            if (i_0p25deg .gt. igrid) then
               i_0p25deg = i_0p25deg - igrid
            else if (i_0p25deg .lt. 1) then
               i_0p25deg = i_0p25deg + igrid
            end if
            j_0p25deg = nint(rj)
            if (j_0p25deg .lt. 1) then
               j_0p25deg = 1
            else if (j_0p25deg .gt. jgrid) then
               j_0p25deg = jgrid
            end if
            USAFSI_arrays%snow_poss(c,r) = &
                 snow_poss_0p25deg(i_0p25deg,j_0p25deg)
         end do ! i
      end do ! j

      ! DEALLOCATE ARRAYS.
      deallocate(climo_0p25deg)
      deallocate(snow_poss_0p25deg)

      return

   end subroutine getgeo

   subroutine getobs (date10, month,  sfcobs, netid,  staid, stacnt, &
        stalat, stalon, staelv, stadep, sfcobsfmt)

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME: GET OBSERVATIONS
      !**
      !**  PURPOSE: OBTAIN SNOW DEPTHS FOR ALL AVAILABLE OBSERVATIONS
      !**
      !**  CALLED FROM:  SNODEP
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  DATE10, MONTH, SFCOBS, BLKLIST
      !**
      !**   OUTPUT: STAID, STACNT, STALAT, STALON, STAELV, STADEP, OBI, OBJ
      !**
      !**  FILES ACCESSED
      !**  ==============
      !**  FILE NAME                               R/W DESCRIPTION
      !**  --------------------------------------- --- -------------------------------
      !**  ${SFCOBS}/sfcsno_nh.06hr.YYYYMMDD12.txt  R  SSMIS NH SNOW AND ICE EDR DATA
      !**  ${SFCOBS}/sfcsno_sh.06hr.YYYYMMDD12.txt  R  SSMIS SH SNOW AND ICE EDR DATA
      !**
      !**   See USAFSI_paramsMod.F90 for common parameter descriptions
      !**
      !**  UPDATES
      !**  =======
      !**  26 APR 95  INITIAL UNISYS VERSION...........................SSGT CONRY/SYSM
      !**  08 FEB 96  ADDED LOOP TO SEARCH FOR SNOW DEPTH AT 6 HOUR INTERVALS.
      !**             REMOVED CALL TO CALCSN THAT USED TO SEARCH THE PAST 24
      !**             HOURS FOR THE DEPTH......................DR KOPP,SSGT CONRY/SYSM
      !**  23 OCT 97  ADDED INCLUDE OF TUNES...........................SSGT CONRY/DNXM
      !**  23 NOV 98  MODIFIED DO LOOP TO LOOK FOR THE CORRECT NUMBER OF
      !**             ADDITIONAL WORDS................................SSGT MILLER/DNXM
      !**  09 MAR 99  ELIMINATED SNOW HOLES CAUSED BY MISINTERPRETATION OF
      !**             MISSING OBSERVATIONS AND REMOVED THE COUNTING OF
      !**             NEGATIVE OBSERVATIONS...........................SRA HERKAMP/DNXM
      !**  22 FEB 01  PORTED FROM UNISYS MAINFRAME TO UNIX............SSGT MILLER/DNXM
      !**  08 NOV 01  ADDED LOGIC TO USE THE FOLLOWING DATABASE DEPTH FLAGS -
      !**             999, 998 AND 997...................................MR GAYNO/DNXM
      !**  30 JAN 03  ADDED BLACKLIST LOGIC..............................MR GAYNO/DNXM
      !**  21 JUL 04  CONVERTED TO FORTRAN 90, MOVED MOST PASSED ARRAYS INTO
      !**             MODULES AND COMPLETED 16TH MESH ADAPTATION......MR EYLANDER/DNXM
      !**  07 MAR 05  CHANGED TO READ OBS FROM FILE INSTEAD OF CDMS...MR LEWISTON/DNXM
      !**  14 APR 05  CHANGED STADEP FROM INTEGER (TENTHS OF INCHES) TO REAL (METERS).
      !**             ADDED PARAMETERS NOCOSNO AND MINSNOW............MR LEWISTON/DNXM
      !**  03 JUN 05  DELETED JULHR INPUT VARIABLE AND CREATION OF "no_obs" FILE
      !**             (PREVIOUSLY FLAGGED GTWAPS SCRIPT TO GENERATE NONFATAL
      !**             TIVOLI MESSAGE).................................MR LEWISTON/DNXM
      !**  15 SEP 08  ADDED FILTER TO REJECT ERRONEOUSLY HIGH TEMPERATURES
      !**             IN THE ARCTIC AND ANTARCTIC REGIONS.........MR LEWISTON/2WXG/WEA
      !**  16 JUN 09  CONVERTED TO EQUIDISTANT CYLINDRICAL GRID...MR LEWISTON/2WXG/WEA
      !**  24 FEB 10  ADDED INPUT PATH AND FILENAMES..............MR LEWISTON/16WS/WXE
      !**  19 JAN 11  CHANGED LONGITUDE TO EAST POSITIVE..........MR LEWISTON/16WS/WXE
      !**  29 APR 11  MOVED BLACKLIST AND VALIDATION TO DBPULL....MR LEWISTON/16WS/WXE
      !**  08 NOV 11  CHANGED ONE IF STATEMENT TO WORK ON LINUX...MR LEWISTON/16WS/WXE
      !**  20 MAR 12  SWITCHED OBS SOURCE FROM CDMS TO JMOBS......MR LEWISTON/16WS/WXE
      !**  30 MAY 12  ADDED ITEMP TO OBS READ, REMOVED CDMS LONGITUDE SWAP,
      !**             ADDED NETID UPDATE FOR OBS WITH SNOW........MR LEWISTON/16WS/WXE
      !**  10 OCT 13  CHANGED FROM 24-HOUR TO 6-HOUR FILES........MR LEWISTON/16WS/WXE
      !**  26 JAN 17  REMOVED LEGACY JMOBS CHECKS AGAINST DEPTH.....MR PUSKAR/16WS/WXE
      !**  21 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
      !**  27 Jul 23  Added new sfcobs file format...Eric Kemp, SSAI
      !**  24 Aug 23  New global sfcsno file format...Eric Kemp, SSAI
      !**
      !*******************************************************************************
      !*******************************************************************************

      ! Imports
      use LDT_logMod, only: LDT_logunit
      use LDT_usafsiMod, only: usafsi_settings
      use map_utils ! EMK
      use USAFSI_paramsMod
      use USAFSI_utilMod ! EMK

      ! Defaults
      implicit none

      ! Arguments
      character*10,  intent(in)   :: date10                ! DATE-TIME GROUP OF CYCLE
      integer, intent(in)         :: month                 ! CURRENT MONTH (1-12)
      character*255, intent(in)   :: sfcobs                ! PATH TO DBPULL SNOW OBS DIRECTORY
      character*5,   intent(out)  :: netid       (:)       ! NETWORK ID OF AN OBSERVATION
      character*32,   intent(out)  :: staid       (:)       ! STATION ID OF AN OBSERVATION

      integer, intent(out)        :: stacnt                ! TOTAL NUMBER OF OBSERVATIONS USED
      integer, intent(out)        :: stalat      (:)       ! LATITUDE OF A STATION OBSERVATION
      integer, intent(out)        :: stalon      (:)       ! LONGITUDE OF A STATION OBSERVATION
      integer, intent(out)        :: staelv      (:)       ! ELEVATION OF A STATION OBSERVATION (METERS)
      real, intent(out)           :: stadep      (:)       ! SNOW DEPTH REPORTED AT A STATION (METERS)
      integer, intent(in) :: sfcobsfmt ! Format of sfcobs file

      ! Local variables
      character*7                 :: access_type           ! FILE ACCESS TYPE
      character*2                 :: chemi       ( 2)      ! HEMISPHERE FOR FILENAME ('nh', 'sh')
      character*2                 :: chemicap    ( 2)      ! HEMISPHERE FOR MESSAGE ('NH', 'SH')

      character*10                :: date10_hourly         ! DATE-TIME GROUP OF HOURLY DATA
      character*10                :: date10_prev           ! DATE-TIME GROUP OF LAST HOUR READ
      character*6                 :: interval              ! TIME INTERVAL FOR FILENAME
      character*4                 :: msgval                ! ERROR MESSAGE VALUE
      character*90                :: message     (msglns)  ! ERROR MESSAGE
      character*255               :: obsfile               ! NAME OF OBSERVATION TEXT FILE
      character*5                 :: obsnet                ! RETURNED OBS STATION NETWORK
      character*32                 :: obssta                ! RETURNED OBS STATION ID
      character*5,   allocatable  :: oldnet      (:)       ! ARRAY OF NETWORKS FOR OLDSTA
      character*32,   allocatable  :: oldsta      (:)       ! ARRAY OF PROCESSED STATIONS WITH SNOW DEPTHS

      character*12                :: routine_name          ! NAME OF THIS SUBROUTINE
      integer                     :: ctrgrd                ! TEMP HOLDER FOR GROUND OBS INFO
      integer                     :: ctrtmp                ! TEMP HOLDER FOR TOO WARM TEMPERATURE OBS
      integer                     :: ctrtrs                ! TEMP HOLDER FOR TEMP THRES OBS
      integer                     :: depth                 ! OBSERVATION'S SNOW DEPTH (CM)
      integer                     :: ground                ! STATE OF GROUND
      integer                     :: hemi                  ! HEMISPHERE (1 = NH, 2 = SH)
      integer                     :: istat                 ! I/O STATUS FOR OBS FILE READ
      integer                     :: itemp                 ! SURFACE TEMPERATURE IN KELVIN (SCALED BY 10)
      integer                     :: j                     ! LOOP COUNTER FOR REPEATED OBSERVATION CHECK
      integer                     :: lunsrc      (2)       ! FORTRAN UNIT NUMBER FOR INPUT FILE

      integer                     :: obcount               ! NUMBER OF OBS IN INPUT FILE
      integer                     :: obelev                ! ELEVATION OF AN OBSERVATION (METERS)
      integer                     :: oblat                 ! LATITUDE OF AN OBSERVATION (*100)
      integer                     :: oblon                 ! LONGITUDE OF AN OBSERVATION (*100)
      integer                     :: obsrtn                ! NUMBER OF OBS RETURNED FOR CURRENT HOUR
      integer                     :: obwsno                ! NUMBER OF OBSERVATIONS WITH A SNOW DEPTH
      integer                     :: printdepth            ! SNOW DEPTH FOR PRINTING TO LOG
      integer                     :: stacnt_h              ! NUMBER OF OBSERVATIONS USED PER HEMISPHERE
      integer                     :: stctp1                ! STATION ARRAY POINTER
      integer                     :: totalobs              ! TOTAL OBSERVATIONS READ FROM FILE
      logical                     :: isfile                ! FLAG INDICATING WHETHER INPUT FILE EXISTS
      logical                     :: isopen                ! FLAG INDICATING WHETHER INPUT FILE IS OPEN
      logical                     :: towarm                ! FLAG INDICATING TEMPERATURE IS TOO WARM FOR SNOW
      real                        :: printlat              ! DESCALED LATITUDE FOR PRINTING TO LOG
      real                        :: printlon              ! DESCALED LONGITUDE FOR PRINTING TO LOG
      real                        :: printtemp             ! DESCALED TEMPERATURE FOR PRINTING TO LOG

      ! DEFINE DATA VALUES
      data chemi        / 'nh', 'sh' /
      data chemicap     / 'NH', 'SH' /
      data interval     / '.06hr.' /
      data lunsrc       / 41, 42 /
      data routine_name / 'GETOBS      '/

      ! ALLOCATE PROCESSED STATION LIST TO MAX NUMBER OF OBS RETURNED.
      allocate (oldnet (usafsi_settings%maxsobs))
      allocate (oldsta (usafsi_settings%maxsobs))

      ! INITIALIZE VARIABLES.
      depth        = missing
      istat        = 0
      message      = ' '
      obcount      = 0
      obelev       = 0
      oblat        = 0
      oblon        = 0
      obsrtn       = 0
      oldnet       = ' '
      oldsta       = ' '
      stacnt       = 0
      stacnt_h     = 0
      stadep       = 0.0
      stctp1       = stacnt + 1
      towarm       = .false.

      hemi_loop : do hemi = 1, 2

         ! INITIALIZE TEMP COUNTERS TO DETERMINE HOW MANY OBS ARE NOT BEING
         ! USED BECAUSE OF UNREALISTIC/UNRELIABLE DATA INFORMATION.

         date10_prev  = ' '
         ctrgrd       = 0
         ctrtmp       = 0
         ctrtrs       = 0
         obwsno       = 0
         totalobs     = 0
         isfile       = .false.
         isopen       = .false.
         message      = ' '

         ! OPEN INPUT FILE.
         if (sfcobsfmt == 1) then
            obsfile = trim(sfcobs) // 'sfcsno_' // chemi(hemi) //      &
                 interval // date10 // '.txt'
         else if (sfcobsfmt == 2) then ! Global file
            obsfile = trim(sfcobs) // 'sfcsno_' //      &
                 '06hr_' // date10 // '.txt'
         end if
         inquire (file=obsfile, exist=isfile)
         file_check : if (isfile) then

            access_type = 'OPENING'
            open (lunsrc(hemi), file=obsfile, iostat=istat, err=5000,  &
                 form='formatted')
            isopen = .true.

            access_type = 'READING'
            write (ldt_logunit,6000) trim(routine_name), trim(obsfile)
            !read (lunsrc(hemi), 6200, iostat=istat, end=3000, err=5000) obcount
            read (lunsrc(hemi), *, iostat=istat, end=3000, err=5000) obcount

            ! LOOP THROUGH ALL OBSERVATIONS RETRIEVED FROM THE DATABASE.
            read_loop : do while (istat .eq. 0)

               if (sfcobsfmt == 1) then
                  read (lunsrc(hemi), 6400, iostat=istat, end=3000, &
                       err=5000) &
                       date10_hourly, obsnet, obssta, oblat, oblon, &
                       obelev,   &
                       itemp, depth, ground
               else if (sfcobsfmt == 2) then
                  ! New format with longer station IDs
                  read (lunsrc(hemi), 6401, iostat=istat, end=3000, &
                       err=5000) &
                       date10_hourly, obsnet, obssta, oblat, oblon, &
                       obelev,   &
                       itemp, depth, ground
               end if
               good_read : if (istat == 0) then

                  if (date10_hourly .ne. date10_prev) then
                     if (totalobs > 1) then
                        if (sfcobsfmt == 1) then
                           write(ldt_logunit,6500) &
                                trim(routine_name), chemicap(hemi),     &
                                date10_prev, obsrtn
                        else
                           write(ldt_logunit,6501) &
                                trim(routine_name),  &
                                date10_prev, obsrtn
                        end if
                        obsrtn = 0
                     end if
                  end if

                  date10_prev = date10_hourly
                  obsrtn = obsrtn + 1
                  totalobs = totalobs + 1

                  ! ENSURE OBSERVATIONS ARE NOT REPEATED.
                  duplicate_obs : do j = 1, stacnt
                     if (obsnet .eq. oldnet(j) .and. &
                          obssta .eq. oldsta(j)) then
                        cycle read_loop
                     endif
                  enddo duplicate_obs

                  ! IF TEMP NOT TOO WARM, CONTINUE PROCESSING.
                  temp_check : if (itemp <= usafsi_settings%thresh) then

                     ! IF LATITUDE IS 40 OR LESS, CHECK WHERE SNOW IS
                     ! UNLIKELY BASED ON ELEVATION, MONTH, AND LATITUDE.
                     if (abs(oblat) <= usafsi_settings%trplat(1)) then
                        call summer (obelev, hemi, oblat, month, towarm)
                     endif

                     climo_check : if (.not. towarm) then

                        ! CHECK FOR SNOW DEPTH VALUE MISSING OR NOT REPORTED.
                        valid_depth: if (depth > misval) then

                           ! INCREMENT COUNTER  FOR OBS REPORTING A DEPTH.

                           if (depth > 0) obwsno = obwsno + 1

                           ! STORE THE STATION NUMBER TO THE STATION LIST.
                           ! CONVERT REPORTED DEPTH FROM MILLIMETERS TO METERS.
                           ! INCREMENT THE STATION COUNTER AND ARRAY POINTER.

                           NETID(STCTP1)  = OBSNET
                           staid(stctp1)  = obssta
                           stadep(stctp1) = (float (depth) / 1000.0) ! convert from mm to meters

                           if (depth >= 1 .and. stadep(stctp1) < 0.001) then
                              write(ldt_logunit,6600) routine_name, depth, &
                                   stadep(stctp1)
                           end if

                           oldnet(stctp1) = obsnet
                           oldsta(stctp1) = obssta

                           stalat(stctp1) = oblat
                           stalon(stctp1) = oblon
                           staelv(stctp1) = obelev

                           stacnt         = stacnt + 1
                           stctp1         = stacnt + 1

                           
                           ! IF DEPTH IS MISSING, NOT REPORTED, OR INDETERMINATE,
                           ! AND STATE OF GROUND STATES NO SNOW, STORE STATION
                           ! NUMBER AND INCREMENT THE COUNTERS.
                        elseif ((depth <= misval)  .and.                     &
                             ((ground >= 0) .and. &
                             (ground <= 9)) ) then valid_depth

                           netid(stctp1)  = obsnet
                           staid(stctp1)  = obssta
                           stadep(stctp1) = 0.0
                           oldnet(stctp1) = obsnet
                           oldsta(stctp1) = obssta

                           stalat(stctp1) = oblat
                           stalon(stctp1) = oblon
                           staelv(stctp1) = obelev

                           stacnt         = stacnt + 1
                           stctp1         = stacnt + 1

                           ctrgrd = ctrgrd + 1

                        endif valid_depth

                     else climo_check

                        ! SNOW UNLIKELY BASED ON ELEVATION, MONTH AND LATITUDE.
                        ! SET AMOUNT TO 0, STORE THE STATION NUMBER, AND
                        ! INCREMENT THE COUNTERS.  RESET TOWARM TO FALSE.

                        netid(stctp1)  = obsnet
                        staid(stctp1)  = obssta
                        stadep(stctp1) = 0.0
                        oldnet(stctp1) = obsnet
                        oldsta(stctp1) = obssta

                        stalat(stctp1) = oblat
                        stalon(stctp1) = oblon
                        staelv(stctp1) = obelev

                        stacnt         = stacnt + 1
                        stctp1         = stacnt + 1
                        towarm         = .false.

                        ctrtmp = ctrtmp + 1

                     endif climo_check

                  else temp_check
                     
                     ! THE TEMPERATURE THRESHOLD HAS BEEN EXCEEDED, SO IF OB IS
                     ! ACCURATE, IT'S UNREALISTIC THIS STATION HAS SNOW.
                     ! IF NOT IN THE ARCTIC OR ANTARCTIC, SET AMOUNT TO 0,
                     ! STORE THE STATION NUMBER, AND INCREMENT COUNTERS.
                     if (oblat > -arctlat .and. oblat < arctlat) then

                        netid(stctp1)  = obsnet
                        staid(stctp1)  = obssta
                        stadep(stctp1) = 0.0
                        oldnet(stctp1) = obsnet
                        oldsta(stctp1) = obssta

                        stalat(stctp1) = oblat
                        stalon(stctp1) = oblon
                        staelv(stctp1) = obelev

                        stacnt         = stacnt + 1
                        stctp1         = stacnt + 1

                        ctrtrs = ctrtrs + 1

                        ! IF WE'RE IN A POLAR REGION AND TEMPERATURE EXCEEDS THE
                        ! MAXIMUM, ASSUME INVALID TEMPERATURE AND REJECT THE OB.
                     else if (itemp > usafsi_settings%arctmax) then

                        printlat = float (oblat) / 100.0
                        printlon = float (oblon) / 100.0

                        if (depth <= misval) then
                           printdepth = -99999
                        else
                           printdepth = depth
                        end if

                        if (itemp <= misval) then
                           printtemp = -9999.9
                        else
                           printtemp = itemp / 10.0
                        end if

                        write(ldt_logunit,6700) obsnet, obssta, printlat, &
                             printlon,     &
                             obelev, printtemp, ground, printdepth

                     endif

                  endif temp_check

               end if good_read

            end do read_loop

            ! DONE READING; CLOSE FILE AND PRINT TOTALS.
            ! PRINT NUMBER OF OBS WHICH REPORTED AN EXPLICIT SNOW DEPTH.
            ! THESE REPORTS HAVE A HEADER OF 43xxx, WHERE xxx IS THE DEPTH.
            ! IF NO OBSERVATIONS FOUND, SET FLAG FOR SCRIPT TO DETECT.

3000        continue

            close (lunsrc(hemi))

            if (totalobs > 0) then

               stacnt_h = stacnt - stacnt_h
               if (sfcobsfmt == 1) then
                  write (ldt_logunit,6500) trim(routine_name), &
                       chemicap(hemi), &
                       date10_prev, obsrtn
                  write (ldt_logunit,6800) trim(routine_name), &
                       chemicap(hemi), &
                       totalobs, &
                       stacnt_h, obwsno, ctrgrd, ctrtmp, ctrtrs
               else if (sfcobsfmt == 2) then
                  write (ldt_logunit,6501) trim(routine_name), &
                       date10_prev, obsrtn
                  write (ldt_logunit,6801) trim(routine_name), &
                       totalobs, &
                       stacnt_h, obwsno, ctrgrd, ctrtmp, ctrtrs
               end if
            else

               if (sfcobsfmt == 1) then
                  message(1) = &
                       '[WARN] NO SURFACE OBSERVATIONS READ FOR ' // &
                       date10 // &
                       ' ' // chemicap(hemi)
               else if (sfcobsfmt == 2) then
                  message(1) = &
                       '[WARN] NO SURFACE OBSERVATIONS READ FOR ' // &
                       date10
               end if
               call error_message (program_name, routine_name, message)

            end if

         else file_check

            if (sfcobsfmt == 1) then
               message(1) = &
                    '[WARN] NO SURFACE OBSERVATIONS FILE FOR ' // &
                    date10 // &
                    ' ' // chemicap(hemi)
            else if (sfcobsfmt == 2) then
               message(1) = &
                    '[WARN] NO SURFACE OBSERVATIONS FILE FOR ' &
                    // date10
            end if
            message(2) = '[WARN] Looked for ' // trim(obsfile)
            call error_message (program_name, routine_name, message)

         end if file_check

         ! New file format is global, so we don't need to loop again
         if (sfcobsfmt == 2) exit

      end do hemi_loop

      ! DEALLOCATE ARRAYS
      deallocate (oldsta)

      return

      ! ERROR-HANDLING SECTION.

5000  continue
      if (isopen) close (lunsrc(hemi))
      message(1) = '[ERR] ERROR ' // access_type // ' ' // trim (obsfile)
      write (msgval, '(i4)') istat
      message(2) = '[ERR] ISTAT = ' // msgval
      call abort_message (program_name, routine_name, message)
      return

      ! FORMAT STATEMENTS.
6000  format (/, '[INFO] ', A, ': READING ', A)
!6200  format (I)
6400  format (A10, 1X, A5, 1X, A10, 6(I10))
!6401  format (A10, 1X, A5, 1X, A31, 1X, 6(I10))
6401  format (A10, 1X, A5, 1X, A32, 6(I10))
6500  format (/, '[INFO] ', A6, ': SURFACE OBS READ FOR ', A2, ' DTG ',     &
           A10, ' = ', I6)
6501  format (/, '[INFO] ', A6, ': SURFACE OBS READ FOR DTG ',     &
           A10, ' = ', I6)
6600  format (1X, '**', A6, ':  DEPTH = ', I6, '   STADEP = ', I6)
6700  format (/, 1X, '[INFO] HIGH POLAR TEMP: NETW= ', A5, 1X, 'STN= ', A31,  &
           1X, 'LAT= ', F8.2, 1X, 'LON= ', F8.2,                  &
           1X, 'ELEV= ', I5, /, 6X, 'TEMP= ', F7.1,               &
           2X, 'ST OF GRND= ', I9, 2X, 'DEPTH(CM)= ', I6)
6800  format (/, 1X, 55('-'),                                           &
           /, 3X, '[INFO] SUBROUTINE:  ', A6,                               &
           /, 5X, '[INFO] TOTAL SURFACE OBS READ FOR ', A2, 9X,' = ',   I6, &
           /, 5X, '[INFO] TOTAL NON-DUPLICATE OBS PROCESSED      = ',   I6, &
           /, 5X, '[INFO] STATIONS WITH A FOUR-THREE GROUP       =   ', I4, &
           /, 5X, '[INFO] OBS NOT USED FOR STATE OF GROUND       =   ', I4, &
           /, 5X, '[INFO] OBS NOT USED FOR SEASON AND ELEVATION  =   ', I4, &
           /, 5X, '[INFO] OBS NOT USED FOR EXCEEDED TEMP THRESH  = ',   I6, &
           /, 1X, 55('-'))
6801  format (/, 1X, 55('-'),                                           &
           /, 3X, '[INFO] SUBROUTINE:  ', A6,                               &
           /, 5X, '[INFO] TOTAL SURFACE OBS READ                 = ',   I6, &
           /, 5X, '[INFO] TOTAL NON-DUPLICATE OBS PROCESSED      = ',   I6, &
           /, 5X, '[INFO] STATIONS WITH A FOUR-THREE GROUP       =   ', I4, &
           /, 5X, '[INFO] OBS NOT USED FOR STATE OF GROUND       =   ', I4, &
           /, 5X, '[INFO] OBS NOT USED FOR SEASON AND ELEVATION  =   ', I4, &
           /, 5X, '[INFO] OBS NOT USED FOR EXCEEDED TEMP THRESH  = ',   I6, &
           /, 1X, 55('-'))

   end subroutine getobs

   subroutine getsfc ( date10, stmpdir, sfctmp_found, sfctmp_lis )

      !*****************************************************************************************
      !*****************************************************************************************
      !**
      !**  NAME:  GET SURFACE TEMPERATURE DATA
      !**
      !**  PURPOSE:  GET SHELTER TEMPERATURES DEGRIBBED FROM LIS
      !**
      !**  CALLED FROM:  SNODEP
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  DATE10, STMPDIR
      !**
      !**   OUTPUT: SFCTMP_FOUND, SFCTMP_LIS
      !**
      !**  FILES ACCESSED
      !**  ==============
      !**  FILE NAME                                            R/W DESCRIPTION
      !**  ---------------------------------------------------- --- -----------------------------
      !**  ${STMPDIR}/lis_sfctmp_0p25deg.yyyymmddhh.dat         R/W LIS 2-METER SHELTER TEMP
      !**
      !**   See USAFSI_paramsMod.F90 for common parameter descriptions
      !**
      !**  UPDATES
      !**  =======
      !**  13 SEP 10  INITIAL VERSION ......................................MR LEWISTON/16WS/WXE
      !**  14 DEC 12  MOVED DEGRIB PROCESS HERE FROM SCRIPT.................MR LEWISTON/16WS/WXE
      !**  10 OCT 13  ADDED STATUS TO DEGRIB ERROR MESSAGES.................MR LEWISTON/16WS/WXE
      !**  26 SEP 17  CHANGED LIS FILENAME FROM AFWA TO 557WW.................MR PUSKAR/16WS/WXE
      !**  21 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
      !**
      !*****************************************************************************************
      !*****************************************************************************************

      ! Imports
      use LDT_coreMod, only: LDT_domain, LDT_rc
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      use map_utils
      use USAFSI_paramsMod
      use USAFSI_utilMod ! EMK

      ! Defaults
      implicit none

      ! Arguments
      character*10,  intent(in)   :: date10                ! SNODEP DATE-TIME GROUP
      character*255, intent(in)   :: stmpdir               ! SFC TEMP DIRECTORY PATH
      logical,       intent(out)  :: sfctmp_found          ! FLAG FOR SFC TEMP FILE FOUND
      real,          intent(out)  :: sfctmp_lis  ( : , : ) ! LIS SURFACE TEMPERATURE DATA

      ! Local constants
      integer, parameter          :: lis_size = igrid * jgrid_lis    ! LIS ARRAY SIZE

      ! Local variables
      character*10                :: dtglis                ! LIS DATE-TIME GROUP
      character*255               :: file_stmp             ! FULLY-QUALIFIED SFCTMP FILE NAME
      character*7                 :: iofunc                ! ACTION TO BE PERFORMED
      character*90                :: message     (msglns)  ! ERROR MESSAGE

      character*12                :: routine_name          ! NAME OF THIS ROUTINE
      integer                     :: icount                ! LOOP COUNTER
      integer                     :: julhr                 ! AFWA JULIAN HOUR
      logical                     :: isfile                ! FLAG FOR INPUT FILE FOUND
      real, allocatable :: sfctmp_lis_0p25deg(:,:)
      integer :: nc,nr
      type(proj_info) :: lis_0p25deg_proj
      real :: ri,rj,rlat,rlon
      integer :: i_0p25deg, j_0p25deg, gindex, r, c
      
      data routine_name / 'GETSFC      ' /

      allocate(sfctmp_lis_0p25deg(igrid, jgrid_lis))

      ! GET LATEST LIS SHELTER TEMPERATURES.
      dtglis       = date10
      icount       = 1
      iofunc       = 'READING'
      isfile       = .false.
      sfctmp_found = .false.
      sfctmp_lis   = -1.0

      file_search : do while ((.not. sfctmp_found) .and.(icount .le. 4))

         ! READ BINARY FILE IF FOUND.
         message   = ' '
         file_stmp = trim (stmpdir) // 'lis_sfctmp' // meshname // '.'   &
              // dtglis // '.dat'
         inquire (file=file_stmp, exist=isfile)

         if (isfile) then

            sfctmp_found = .true.
            write(ldt_logunit,*)'[INFO] Reading ', trim(file_stmp)
            write (ldt_logunit, 6000) routine_name, iofunc, file_stmp
            call putget_real (sfctmp_lis_0p25deg, 'r', file_stmp, &
                 program_name,   &
                 routine_name, igrid, jgrid_lis)

         end if

         ! IF NOT FOUND YET, BACK UP TO PREVIOUS LIS CYCLE.
         if (.not. sfctmp_found) then
            call date10_julhr (dtglis, julhr, program_name, routine_name)
            julhr  = julhr  - 6
            icount = icount + 1
            call julhr_date10 (julhr, dtglis, program_name, routine_name)
         end if

      end do file_search
      
      ! IF NOT FOUND FOR PAST 24 HOURS, SEND ERROR MESSAGE.
      if (.not. sfctmp_found) then
         message(1) = '[WARN] LIS DATA NOT FOUND FOR PAST 24 HOURS'
         call error_message (program_name, routine_name, message)
         write (ldt_logunit, 6400) routine_name
      end if

      ! Interpolate to LDT grid
      if (sfctmp_found) then
         call map_set(proj_code=proj_latlon, &
              lat1 =  -59.875, &
              lon1 = -179.875, &
              dx = 0.25, &
              stdlon = 0.25, &
              truelat1 = 0.25, &
              truelat2 = 0., &
              idim = igrid, &
              jdim = jgrid_lis, &
              proj = lis_0p25deg_proj)
         sfctmp_lis(:,:) = -1
         nr = LDT_rc%lnr(1)
         nc = LDT_rc%lnc(1)
         do r = 1, nr
            do c = 1, nc
               gindex = c + (r-1)*LDT_rc%lnc(1)
               rlat = LDT_domain(1)%lat(gindex)
               rlon = LDT_domain(1)%lon(gindex)
               call latlon_to_ij(lis_0p25deg_proj, rlat, rlon, ri, rj)
               i_0p25deg = nint(ri)
               if (i_0p25deg .gt. igrid) then
                  i_0p25deg = i_0p25deg - igrid
               else if (i_0p25deg .lt. 1) then
                  i_0p25deg = i_0p25deg + igrid
               end if
               j_0p25deg = nint(rj)
               if (j_0p25deg .lt. 1) then
                  j_0p25deg = 1
               else if (j_0p25deg .gt. jgrid_lis) then
                  j_0p25deg = jgrid_lis
               end if
               if (sfctmp_lis_0p25deg(i_0p25deg, j_0p25deg) < -1) then
                  sfctmp_lis(c,r) = -1
               else
                  sfctmp_lis(c,r) = sfctmp_lis_0p25deg(i_0p25deg, j_0p25deg)
               end if
            end do ! c
         end do ! r
                  
      end if
      
      ! Clean up
      deallocate(sfctmp_lis_0p25deg)
      
      return

      ! Format statements
6000  format (/, '[INFO] ', A6, ': ', A7, 1X, A)
!6200  format (/, '[WARN] ', A6, ': ERROR ', A7, 1X, A, /, 3X, 'STATUS = ', I6)
6400  format (/, '[WARN]', A6, ': LIS DATA NOT FOUND FOR PAST 24 HOURS')

   end subroutine getsfc

   subroutine getsmi (date10, ssmis)

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME: GET SPECIAL SENSOR MICROWAVE IMAGER/SOUNDER (SSMIS) DATA
      !**
      !**  PURPOSE: GET SSMIS SNOW DEPTH AND ICE CONCENTRATION EDR DATA
      !**
      !**  CALLED FROM: SNODEP
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  DATE10, SSMIS
      !**
      !**   OUTPUT: SSMIS_DEPTH, SSMIS_ICECON
      !**
      !**  FILES ACCESSED
      !**  ==============
      !**  FILE NAME                           R/W DESCRIPTION
      !**  ----------------------------------- --- -----------------------------------
      !**  ssmis_snoice_nh.06hr.YYYYMMDDHH.txt  R  SSMIS NH SNOW AND ICE EDR DATA
      !**  ssmis_snoice_sh.06hr.YYYYMMDDHH.txt  R  SSMIS SH SNOW AND ICE EDR DATA
      !**
      !**   See USAFSI_arraysMod.F90 for common array descriptions
      !**   See USAFSI_paramsMod.F90 for common parameter descriptions
      !**
      !**  UPDATES
      !**  =======
      !**  26 APR 95  INITIAL UNISYS VERSION...........................SSGT CONRY/SYSM
      !**  08 FEB 96  READ 5 CHANNEL TEMPERATURES AND DETERMINE SNOW COVERAGE
      !**             ........................................DR KOPP, SSGT CONRY/SYSM
      !**  05 MAY 98  UPDATED SNOW DETECTION ALGORITHM.  ADDED TWO NEW SSMI
      !**             CHANNELS (0#MIH3(8))................................DR KOPP/DNXM
      !**  10 MAR 99  ADDED CALL TO SUMMER TO CHECK FOR TOWARM AREAS..SRA HERKAMP/DNXM
      !**  22 FEB 01  PORTED FROM UNISYS MAINFRAME TO UNIX............SSGT MILLER/DNXM
      !**  11 MAR 03  MODIFIED TO USE SSM/I BRIGHTNESS TEMPERATURE SCALED BY 100.
      !**             MODIFIED SNOW DETECTION LOGIC TO USE REAL ARITHMETIC.
      !**             ADDED OUTPUT OF SNOW DEPTH ALGORITHM TO FILE.......MR GAYNO/DNXM
      !**  11 JUL 03  REPLACED CAL-VAL ICE CONCENTRATIONS FROM CDMS TERADATA
      !**             WITH NASA TEAM II ICE CONCENTRATIONS FROM FILE.....MR GAYNO/DNXM
      !**  05 MAY 04  INCORPORATED NASA TEAM II SEA ICE CALCULATIONS INTO THE
      !**             MODEL INSTEAD OF READING ICE CONCENTRATIONS FROM FILE.
      !**             MODIFIED TO BE FORTRAN 90 FREE FORM COMPLIANT.
      !**             COMPLETED 16TH MESH ADAPTATION..................MR EYLANDER/DNXM
      !**  15 MAR 05  CHANGED TO READ OBS FROM FILE INSTEAD OF CDMS...MR LEWISTON/DNXM
      !**  03 MAY 05  CHANGED TO COMPUTE SNOW OVER LAND ONLY. ADDED OPTIONAL
      !**             SAVE OF GRIDDED CHANNEL TEMPERATURES............MR LEWISTON/DNXM
      !**  03 JUN 05  DELETED JULHR INPUT VARIABLE AND CREATION OF "no_obs"
      !**             FILE (PREVIOUSLY FLAGGED GTWAPS SCRIPT TO GENERATE
      !**             NONFATAL TIVOLI MESSAGE)........................MR LEWISTON/DNXM
      !**  19 JUN 09  CONVERTED TO EQUIDISTANT CYLINDRICAL GRID...MR LEWISTON/2WXG/WEA
      !**  08 FEB 10  CONVERTED FROM SSMIS SDRS TO EDRS...........MR LEWISTON/16WS/WXE
      !**  24 FEB 10  ADDED INPUT PATH AND FILENAMES..............MR LEWISTON/16WS/WXE
      !**  26 MAY 11  COMBINED SNOW AND ICE IN SAME FILE..........MR LEWISTON/16WS/WXE
      !**  10 OCT 13  CHANGED FROM 24-HOUR TO 6-HOUR FILES........MR LEWISTON/16WS/WXE
      !**  16 JAN 18  ADDED DTG TO "NO SSMIS EDR" MESSAGES.......MR LEWISTON/16WS/WXET
      !**  21 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
      !**  28 Jan 21  Updated messages.....................Yeosang Yoon/NASA GSFC/SAIC
      !**
      !*******************************************************************************
      !*******************************************************************************

      ! Imports
      use LDT_coreMod, only: LDT_rc, LDT_domain
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      use map_utils ! EMK
      use USAFSI_arraysMod, only: USAFSI_arrays
      use USAFSI_paramsMod
      use USAFSI_utilMod ! EMK

      ! Defaults
      implicit none

      ! Arguments
      character*10,  intent(in)   :: date10                ! DATE-TIME GROUP OF CYCLE
      character*255, intent(in)   :: ssmis                 ! SSMIS FILE DIRECTORY PATH

      ! Local variables
      character*7                 :: access_type           ! FILE ACCESS TYPE
      character*2                 :: chemicap    ( 2)      ! HEMISPHERE FOR MESSAGE ('NH', 'SH')
      character*2                 :: chemifile   ( 2)      ! HEMISPHERE FOR FILENAME ('nh', 'sh')
      character*10                :: date10_hourly         ! DATE-TIME GROUP OF HOURLY DATA
      character*10                :: date10_prev           ! DATE-TIME GROUP OF LAST HOUR READ
      character*255               :: file_path             ! SSMIS SNOW OR ICE EDR TEXT FILE
      character*6                 :: interval              ! TIME INTERVAL FOR FILENAME
      character*90                :: message     (msglns)  ! ERROR MESSAGE
      character*4                 :: msgval                ! PLACEHOLDER FOR ERROR MESSAGE VALUES
      character*12                :: routine_name          ! NAME OF THIS SUBROUTINE
      integer                     :: edri16                ! EDR 16TH MESH I-COORDINATE
      integer                     :: edrj16                ! EDR 16TH MESH J-COORDINATE
      integer                     :: edrlat                ! EDR LATITUDE (100THS OF DEGREES)
      integer                     :: edrlon                ! EDR LONGITUDE (100THS OF DEGREES)
      integer                     :: hemi                  ! HEMISPHERE (1 = NH, 2 = SH)
      integer                     :: i                     ! GRID COORDINATE
      integer                     :: icerecs               ! ICE EDRS READ PER HEMI, DATE, SATELLITE
      integer                     :: iceval                ! EDR ICE CONCENTRATION VALUE
      integer                     :: irec                  ! LOOP COUNTER FOR FILE READ
      integer                     :: istat                 ! I/O STATUS
      integer                     :: j                     ! GRID COORDINATE
      integer                     :: lunsrc      ( 2)      ! FORTRAN UNIT NUMBER FOR INPUT
      integer                     :: msgline               ! ERROR MESSAGE LINE NUMBER
      integer                     :: prevsat               ! PREVIOUSLY RETRIEVED SATELLITE ID
      integer                     :: satid                 ! CURRENT SATELLITE ID
      integer                     :: snorecs               ! ICE EDRS READ PER HEMI, DATE, SATELLITE
      integer                     :: snoval                ! EDR SNOW DEPTH VALUE
      logical                     :: goodread              ! FLAG FOR SUCCESSFUL DATABASE READ
      logical                     :: isfile                ! FLAG INDICATING WHETHER INPUT FILE EXISTS
      logical                     :: isopen                ! FLAG INDICATING WHETHER INPUT FILE IS OPEN
      real                        :: rlat                  ! EDR LATITUDE (DEGREES)
      real                        :: rlon                  ! EDR LONGITUDE (DEGREES)
      real :: ri,rj
      integer :: nc,nr
      integer :: c,r
      integer, allocatable :: icecount_0p25deg(:,:)
      integer, allocatable :: icetotal_0p25deg(:,:)
      integer, allocatable :: snocount_0p25deg(:,:)
      integer, allocatable :: snototal_0p25deg(:,:)
      type(proj_info) :: snodep_0p25deg_proj
      integer :: i_0p25deg, j_0p25deg
      integer :: gindex
      integer, allocatable :: ssmis_icecon_0p25deg(:,:)
      real, allocatable :: ssmis_depth_0p25deg(:,:)

      ! DEFINE DATA VALUES
      data chemicap     / 'NH', 'SH' /
      data chemifile    / 'nh', 'sh' /
      data interval     / '.06hr.' /
      data lunsrc       /  43,  44 /
      data routine_name / 'GETSMI      '/

      ! ALLOCATE ARRAYS.
      allocate (icecount_0p25deg (igrid , jgrid))
      allocate (icetotal_0p25deg (igrid , jgrid))
      allocate (snocount_0p25deg (igrid , jgrid))
      allocate (snototal_0p25deg (igrid , jgrid))
      icecount_0p25deg(:,:) = 0
      icetotal_0p25deg(:,:) = 0
      snocount_0p25deg(:,:) = 0
      snototal_0p25deg(:,:) = 0

      call map_set(proj_code=proj_latlon, &
           lat1=begin_lat, &
           lon1=begin_lon, &
           dx=0.25, &
           stdlon=0.25, &
           truelat1=0.25, &
           truelat2=0., &
           idim=igrid,&
           jdim=jgrid, &
           proj=snodep_0p25deg_proj)

      ! INITIALIZE VARIABLES.
      date10_prev  = ' '
      goodread     = .false.
      icecount_0p25deg     = 0
      icerecs      = 0
      icetotal_0p25deg     = 0
      isfile       = .false.
      isopen       = .false.
      message      = ' '
      msgline      = 1
      prevsat      = 0
      snocount_0p25deg     = 0
      snorecs      = 0
      snototal_0p25deg     = 0
      USAFSI_arrays%ssmis_depth  = misanl
      USAFSI_arrays%ssmis_icecon = missing

      ! LOOP THROUGH HEMISPHERES.
      hemi_loop : do hemi = 1, 2

         ! OPEN INPUT FILE.
         file_path = trim(ssmis) // 'ssmis_snoice_' //                   &
              chemifile(hemi) // interval // date10 // '.txt'

         inquire (file=file_path, exist=isfile)

         file_check : if (isfile) then

            access_type = 'OPENING'
            open (lunsrc(hemi), file=file_path, iostat=istat, err=5000,   &
                 form='formatted')
            isopen = .true.

            ! LOOP THROUGH ALL EDRS RETRIEVED FROM THE DATABASE.
            write (ldt_logunit,6000) trim (routine_name), trim(file_path)
            irec = 1
            access_type = 'READING'

            read_loop : do while (istat .eq. 0)

               read (lunsrc(hemi), 6200, iostat=istat, end=3000, err=5000) &
                    date10_hourly, satid, edrlat, edrlon, edri16, edrj16, &
                    iceval, snoval

               goodread = .true.

               if (satid .ne. prevsat .or. date10_hourly .ne. date10_prev) then
                  if (irec > 1) then
                     write (ldt_logunit,6400) &
                          trim(routine_name), chemicap(hemi),      &
                          date10_prev, prevsat, icerecs, snorecs
                     icerecs = 0
                     snorecs = 0
                  end if
                  date10_prev = date10_hourly
                  prevsat = satid
               end if

               rlat = float(edrlat) / 100.0
               rlon = float(edrlon) / 100.0
               call latlon_to_ij(snodep_0p25deg_proj,rlat,rlon,ri,rj)
               i = nint(ri)
               j = nint(rj)

               if (iceval .ge. 0) then
                  icerecs = icerecs + 1
                  icetotal_0p25deg (i, j) = icetotal_0p25deg(i, j) + iceval
                  icecount_0p25deg (i, j) = icecount_0p25deg(i, j) + 1
               else if (snoval .ge. 0) then
                  snorecs = snorecs + 1
                  snototal_0p25deg (i, j) = snototal_0p25deg (i, j) + snoval
                  snocount_0p25deg (i, j) = snocount_0p25deg (i, j) + 1
               end if

               irec = irec + 1

            end do read_loop

            ! DONE READING; PRINT LAST RETRIEVAL STATS AND CLOSE FILE.
3000        continue

            close (lunsrc(hemi))
            isopen = .false.

            if (irec > 1) then

               write (ldt_logunit,6400) &
                    trim(routine_name), chemicap(hemi),          &
                    date10_prev, prevsat, icerecs, snorecs

               write (ldt_logunit,6600) &
                    trim(routine_name), chemicap(hemi), irec-1

            else

               message(msgline) = 'NO PMW READ FOR ' // date10 //   &
                    ' ' // chemicap(hemi)
               msgline = msgline + 1

            end if

         else file_check

            message(msgline) = 'NO PMW FILE FOR ' // date10 //      &
                 ' ' // chemicap(hemi)
            msgline = msgline + 1

         end if file_check

      end do hemi_loop

      ! COMPUTE AVERAGE ICE CONCENTRATIONS (WHOLE PERCENT).
      ! CONVERT AVERAGE SNOW DEPTH FROM TENTHS OF MILLIMETERS TO METERS.
      if (goodread) then
         ! Estimate the ice concentration and snow depth on the 0.25 deg
         ! grid
         allocate(ssmis_icecon_0p25deg(igrid,jgrid))
         allocate(ssmis_depth_0p25deg(igrid,jgrid))

         ssmis_icecon_0p25deg = -1
         ssmis_depth_0p25deg = -1

         do j = 1, jgrid
            do i = 1, igrid
               if (icecount_0p25deg (i, j) > 0) then
                  ssmis_icecon_0p25deg (i, j) = &
                       nint (real(icetotal_0p25deg (i, j)) /       &
                       real(icecount_0p25deg (i, j)))
               end if
               if (snocount_0p25deg (i, j) > 0) then
                  ssmis_depth_0p25deg (i, j) = &
                       real(snototal_0p25deg (i, j)) /              &
                       real(snocount_0p25deg (i, j)) / 10000.0
               end if
            end do
         end do
         
         ! Interpolate the 0.25deg data to the LDT grid
         nr = LDT_rc%lnr(1)
         nc = LDT_rc%lnc(1)
         do r = 1, nr
            do c = 1, nc
               gindex = c + (r-1)*nc
               rlat = LDT_domain(1)%lat(gindex)
               rlon = LDT_domain(1)%lon(gindex)
               call latlon_to_ij(snodep_0p25deg_proj, &
                    rlat,rlon,ri,rj)
               i_0p25deg = nint(ri)
               if (i_0p25deg < 1) then
                  i_0p25deg = i + igrid
               else if (i_0p25deg > igrid) then
                  i_0p25deg = i_0p25deg - igrid
               end if
               j_0p25deg = nint(rj)
               if (j_0p25deg < 1) then
                  j_0p25deg = 1
               else if (j_0p25deg > jgrid) then
                  j_0p25deg = jgrid
               end if

               USAFSI_arrays%ssmis_icecon(c,r) = &
                    ssmis_icecon_0p25deg(i_0p25deg, j_0p25deg)
               USAFSI_arrays%ssmis_depth(c,r) = &
                    ssmis_depth_0p25deg(i_0p25deg, j_0p25deg)

            end do ! c
         end do ! r

         deallocate(ssmis_icecon_0p25deg)
         deallocate(ssmis_depth_0p25deg)

      else

         message(msgline) = '[WARN] no ice and snow data received'
         msgline = msgline + 1

      end if

      if (msgline > 1) then

         call error_message (program_name, routine_name, message)

      end if

      deallocate (icetotal_0p25deg)
      deallocate (icecount_0p25deg)
      deallocate (snototal_0p25deg)
      deallocate (snocount_0p25deg)

      return

      ! ERROR-HANDLING SECTION.
5000  continue
      if (isopen) close (lunsrc(hemi))
      message(1) = '[ERR] ERROR ' // access_type // ' PMW FILE'
      message(2) = '[ERR] ' // trim ( file_path )
      write (msgval, '(i4)') istat
      message(3) = '[ERR] ISTAT = ' // msgval
      call abort_message (program_name, routine_name, message)
      return

      ! FORMAT STATEMENTS
6000  format (/, '[INFO] ', A, ': READING ', A)
6200  format (A10, I3, I6, I7, 2(I5), 2(I6))
6400  format (/, '[INFO] ', A, ': READ FOR ', A2, 1X,  A10,            &
           ' SATELLITE F', I2, ': ICE = ', I6, '  SNOW = ', I6)
6600  format (/, 1X, 55('-'),                                           &
           /, '[INFO] ', A6, ': TOTAL READ FOR ', A2, ' = ', I7,    &
           /, 1X, 55('-'))

   end subroutine getsmi

   subroutine getsno (date10, modif, unmod, nc, nr, landice, julhr_beg, &
        just_12z)

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME: GET SNOW DATA
      !**
      !**  PURPOSE:  RETRIEVE PREVIOUS DAY'S SNOW AND ICE DATA
      !**
      !**  CALLED FROM:  SNODEP
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  DATE10, PREVDIR
      !**
      !**   OUTPUT: ICEAGE, OLDCON, OLDDEP, OLDMASK, SNOAGE
      !**
      !**  FILES ACCESSED
      !**  ==============
      !**  FILE NAME                                     I/O DESCRIPTION
      !**  --------------------------------------------- --- ------------------
      !**  ${PREVDIR}/iceage_0p25deg.${DATE10_PREV}.dat   R  ICE AGE
      !**  ${PREVDIR}/icecon_0p25deg.${DATE10_PREV}.dat   R  ICE CONCENTRATIONS
      !**  ${PREVDIR}/icemask_0p25deg.${DATE10_PREV}.dat  R  ICE MASK
      !**  ${PREVDIR}/snoage_0p25deg.${DATE10_PREV}.dat   R  SNOW AGE
      !**  ${PREVDIR}/snodep_0p25deg.${DATE10_PREV}.dat   R  SNOW DEPTH
      !**
      !**   See USAFSI_arraysMod.F90 for common array descriptions
      !**   See USAFSI_paramsMod.F90 for common parameter descriptions
      !**
      !**  UPDATES
      !**  =======
      !**  26 APR 95  INITIAL UNISYS VERSION...........................SSGT CONRY/SYSM
      !**  05 FEB 96  REMOVED CALL TO ADDMUP...................DR KOPP,SSGT CONRY/SYSM
      !**  22 FEB 01  PORTED TO UNIX FROM UNISYS MAINFRAME............SSGT MILLER/DNXM
      !**  21 JUL O4  CONVERTED TO FORTRAN 90.  COMPLETED 16TH MESH UPGRADE.
      !**             MOVED MOST PASSED ARRAYS INTO MODULES.    MODIFIED CLIMO
      !**             READ SINCE IT'S NOW REAL, NOT INTEGER...........MR EYLANDER/DNXM
      !**  14 APR 05  CHANGED AGE (REAL) TO SNOAGE (INTEGER).
      !**             CHANGED OLDDEP FROM INCHES*10 TO METERS.........MR LEWISTON/DNXM
      !**  16 JUN 09  CONVERTED TO EQUIDISTANT CYLINDRICAL GRID...MR LEWISTON/2WXG/WEA
      !**  24 FEB 10  REPLACED FILE COMMANDS WITH PUTGET CALLS....MR LEWISTON/16WS/WXE
      !**  29 DEC 10  UPDATED FILENAMES...........................MR LEWISTON/16WS/WXE
      !**  29 MAY 12  MOVED MODIFIED DATA TEST HERE FROM SCRIPT...MR LEWISTON/16WS/WXE
      !**  10 OCT 13  CHANGED FROM 24-HOUR TO 6-HOUR CYCLES.......MR LEWISTON/16WS/WXE
      !**  21 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
      !**
      !*******************************************************************************
      !*******************************************************************************

      ! Imports
      use LDT_coreMod, only: LDT_domain, LDT_rc
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      use map_utils
      use USAFSI_arraysMod, only: USAFSI_arrays
      use USAFSI_paramsMod
      use USAFSI_utilMod ! EMK

      ! Defaults
      implicit none

      ! Arguments
      character*10,  intent(in)  :: date10           ! CURRENT CYCLE DATE-TIME GROUP
      character*255, intent(in)  :: modif            ! PATH TO MODIFIED DATA DIRECTORY
      character*255, intent(in)  :: unmod            ! PATH TO UNMODIFIED DATA DIRECTORY
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: landice(nc,nr)
      integer,       intent(out) :: julhr_beg        ! AFWA JULIAN HOUR OF PREVIOUS CYCLE
      logical, intent(in) :: just_12z

      ! Local variables
      character*10               :: date10_prev      ! PREVIOUS CYCLE DATE-TIME GROUP
      character*255              :: file_path        ! INPUT FILE PATH AND NAME
      character*90               :: message (msglns) ! ERROR MESSAGE
      character*12               :: routine_name     ! NAME OF THIS SUBROUTINE
      character*255              :: prevdir          ! PATH TO PREVIOUS CYCLE'S DATA
      integer                    :: runcycle         ! CYCLE HOUR
      integer                    :: julhr            ! AFWA JULIAN HOUR
      integer                    :: limit            ! LIMIT ON NUMBER OF CYCLES TO SEARCH
      integer                    :: tries            ! NUMBER OF CYCLES SEARCHED
      logical                    :: found            ! FLAG INDICATING DATA FOUND
      logical                    :: isfile           ! FLAG INDICATING WHETHER FILE EXISTS
      integer, allocatable :: snoage_0p25deg(:,:)
      real, allocatable :: olddep_0p25deg(:,:)
      integer, allocatable :: iceage_0p25deg(:,:)
      integer, allocatable :: oldcon_0p25deg(:,:)
      integer, allocatable :: oldmask_0p25deg(:,:)
      integer, allocatable :: snoage12z_0p25deg(:,:)
      integer, allocatable :: iceage12z_0p25deg(:,:)
      type(proj_info) :: snodep_0p25deg_proj
      integer :: i_0p25deg, j_0p25deg
      integer :: gindex
      real :: rlat,rlon,ri,rj
      integer :: c,r

      data routine_name / 'GETSNO      '/

      ! julhr must always be calculated 
      call date10_julhr (date10, julhr, program_name, routine_name)

      ! FIND THE DATE/TIME GROUP OF THE PREVIOUS CYCLE.
      if (.not. just_12z) then
         found   = .false.
         isfile  = .false.
         limit   = 20
         message = ' '
         tries   = 1

         allocate(snoage_0p25deg(igrid,jgrid))
         allocate(iceage_0p25deg(igrid,jgrid))
         allocate(olddep_0p25deg(igrid,jgrid))
         allocate(oldcon_0p25deg(igrid,jgrid))
         allocate(oldmask_0p25deg(igrid,jgrid))
         
         julhr_beg = julhr
         
         cycle_loop : do while ((.not. found) .and. (tries .le. limit))
            julhr_beg = julhr_beg - 6
            call julhr_date10 (julhr_beg, date10_prev, program_name, &
                 routine_name)
            read (date10_prev(9:10), '(i2)', err=4200) runcycle
            
            ! CHECK FOR PREVIOUS SNOW AGE.
            ! IF 12Z DETERMINE IF MODIFIED DATA AVAILABLE. IF NOT, USE 
            ! UNMODIFIED.
            if (runcycle .eq. 12) then
               
               prevdir = modif
               file_path = trim(prevdir) // 'snoage' // meshname // '.' //  &
                    date10_prev // '.dat'
               
               inquire (file=file_path, exist=isfile)
               
               if (isfile) then
                  found = .true.
               else
                  write (ldt_logunit,6000) trim (routine_name), date10_prev
               end if
               
            end if
            
            if (.not. found) then
               prevdir = unmod
               file_path = trim(prevdir) // 'snoage' // meshname // '.' //   &
                    date10_prev // '.dat'
               
               inquire (file=file_path, exist=isfile)
               
               if (isfile) found = .true.
               
            end if
            
            if (found) then
               
               ! RETRIEVE PREVIOUS SNOW AGE.
               write(LDT_logunit,*) &
                    '[INFO] Reading ',trim(file_path)
               call putget_int (snoage_0p25deg, 'r', file_path, program_name, &
                    routine_name, igrid, jgrid)
               
               ! RETRIEVE PREVIOUS SNOW DEPTH.
               file_path = trim(prevdir) // 'snodep' // meshname // '.' //    &
                    date10_prev // '.dat'
               write(LDT_logunit,*) &
                    '[INFO] Reading ',trim(file_path)
               call putget_real (olddep_0p25deg, 'r', file_path, program_name,&
                    routine_name, igrid, jgrid)
               
               ! RETRIEVE PREVIOUS ICE AGE.
               file_path = trim(prevdir) // 'iceage' // meshname // '.' //    &
                    date10_prev // '.dat'
               write(LDT_logunit,*) &
                    '[INFO] Reading ',trim(file_path)
               call putget_int (iceage_0p25deg, 'r', file_path, program_name, &
                    routine_name, igrid, jgrid)
               
               ! RETRIEVE PREVIOUS ICE CONCENTRATIONS.
               file_path = trim(prevdir) // 'icecon' // meshname // '.' //    &
                    date10_prev // '.dat'
               write(LDT_logunit,*) &
                    '[INFO] Reading ',trim(file_path)
               call putget_int (oldcon_0p25deg, 'r', file_path, program_name, &
                    routine_name, igrid, jgrid)
               
               ! RETRIEVE PREVIOUS ICE MASK.
               file_path = trim(prevdir) // 'icemask' // meshname // '.' //   &
                    date10_prev // '.dat'
               write(LDT_logunit,*) &
                    '[INFO] Reading ',trim(file_path)
               call putget_int (oldmask_0p25deg, 'r', file_path, program_name,&
                    routine_name, igrid, jgrid)
               
            else
               write (ldt_logunit,6200) trim (routine_name), date10_prev
               tries = tries + 1
            end if
            
         end do cycle_loop
      end if

      ! IF 12Z CYCLE, RETRIEVE LAST 12Z SNOW AND ICE AGE.
      read (date10(9:10), '(i2)', err=4200) runcycle
      if (just_12z) then
         found = .true.
      end if
      if (found .and. runcycle .eq. 12) then

         allocate(snoage12z_0p25deg(igrid,jgrid))
         allocate(iceage12z_0p25deg(igrid,jgrid))

         found   = .false.
         isfile  = .false.
         prevdir = modif
         limit   = 5
         tries   = 1
         
         cycle12_loop : do while ((.not. found) .and. (tries .le. limit))
            
            julhr = julhr - 24
            call julhr_date10 (julhr, date10_prev, program_name, routine_name)
            
            ! RETRIEVE PREVIOUS 12Z SNOW AGE.
            ! DETERMINE IF MODIFIED DATA AVAILABLE. IF NOT, USE UNMODIFIED.
            file_path = trim(prevdir) // 'snoage' // meshname // '.' //       &
                 date10_prev // '.dat'
            
            inquire (file=file_path, exist=isfile)

            if (isfile) then
               found = .true.
            else
               write (ldt_logunit,6000) trim (routine_name), date10_prev
               prevdir = unmod
               file_path = trim(prevdir) // 'snoage' // meshname // '.' //  &
                    date10_prev // '.dat'
               inquire (file=file_path, exist=isfile)
               if (isfile) found = .true.
            end if
            
            if (found) then
               
               ! RETRIEVE PREVIOUS 12Z SNOW AND ICE AGE.
               write(LDT_logunit,*) &
                    '[INFO] Reading ',trim(file_path)
               call putget_int (snoage12z_0p25deg, 'r', file_path, &
                    program_name, &
                    routine_name, igrid, jgrid)

               file_path = trim(prevdir) // 'iceage' // meshname // '.' //    &
                    date10_prev // '.dat'
               write(LDT_logunit,*) &
                    '[INFO] Reading ',trim(file_path)
               call putget_int (iceage12z_0p25deg, 'r', file_path, &
                    program_name, &
                    routine_name, igrid, jgrid)

            else
               write (ldt_logunit,6200) trim (routine_name), date10_prev
               tries = tries + 1
            end if

         end do cycle12_loop

      end if

      ! IF NO PREVIOUS DATA FOUND, CANNOT CONTINUE.
      if (.not. found) then
         message(1) = '[ERR] PREVIOUS USAFSI DATA MISSING'
         message(2) = '[ERR] USAFSI CANNOT CONTINUE'
         call abort_message (program_name, routine_name, message)
      end if

      ! At this point, we have data on the 0.25 deg SNODEP domain.  Interpolate
      ! to the LDT grid.
      call map_set(proj_code=proj_latlon, &
           lat1=begin_lat, &
           lon1=begin_lon, &
           dx=0.25, &
           stdlon=0.25, &    ! dlon
           truelat1=0.25, &  ! dlat
           truelat2=0., &
           idim=igrid,&
           jdim=jgrid, &
           proj=snodep_0p25deg_proj)

      ! For simplicity, just use the value of the 0.25 deg grid box that the
      ! LDT grid point is within.
      if (.not. just_12z) then
         USAFSI_arrays%olddep(:,:) = -1
         USAFSI_arrays%snoage(:,:) = -1
         USAFSI_arrays%iceage(:,:) = -1
         USAFSI_arrays%oldcon(:,:) = -1
         USAFSI_arrays%oldmask(:,:) = -1
      end if
      USAFSI_arrays%snoage12z(:,:) = -1
      USAFSI_arrays%iceage12z(:,:) = -1
      do r = 1,nr
         do c = 1,nc
            gindex = c+(r-1)*LDT_rc%lnc(1)
            rlat = LDT_domain(1)%lat(gindex)
            rlon = LDT_domain(1)%lon(gindex)
            call latlon_to_ij(snodep_0p25deg_proj,rlat,rlon,ri,rj)

            i_0p25deg = nint(ri)
            if (i_0p25deg .gt. igrid) then
               i_0p25deg = i_0p25deg - igrid
            else if (i_0p25deg .lt. 1) then
               i_0p25deg = i_0p25deg + igrid
            end if
            j_0p25deg = nint(rj)
            if (j_0p25deg .lt. 1) then
               j_0p25deg = 1
            else if (j_0p25deg .gt. jgrid) then
               j_0p25deg = jgrid
            end if
            if (.not. just_12z) then
               USAFSI_arrays%olddep(c,r) = olddep_0p25deg(i_0p25deg,j_0p25deg)
               USAFSI_arrays%snoage(c,r) = snoage_0p25deg(i_0p25deg,j_0p25deg)
               USAFSI_arrays%iceage(c,r) = iceage_0p25deg(i_0p25deg,j_0p25deg)
               USAFSI_arrays%oldcon(c,r) = oldcon_0p25deg(i_0p25deg,j_0p25deg)
               USAFSI_arrays%oldmask(c,r) = &
                    oldmask_0p25deg(i_0p25deg,j_0p25deg)
            endif
            if (runcycle .eq. 12) then
               USAFSI_arrays%snoage12z(c,r) = &
                    snoage12z_0p25deg(i_0p25deg,j_0p25deg)
               USAFSI_arrays%iceage12z(c,r) = &
                    iceage12z_0p25deg(i_0p25deg,j_0p25deg)
            end if
         end do ! c
      end do ! r

      ! Clean up
      if (allocated(snoage_0p25deg)) deallocate(snoage_0p25deg)
      if (allocated(olddep_0p25deg)) deallocate(olddep_0p25deg)
      if (allocated(iceage_0p25deg)) deallocate(iceage_0p25deg)
      if (allocated(oldcon_0p25deg)) deallocate(oldcon_0p25deg)
      if (allocated(oldmask_0p25deg)) deallocate(oldmask_0p25deg)
      if (allocated(snoage12z_0p25deg)) deallocate(snoage12z_0p25deg)
      if (allocated(iceage12z_0p25deg)) deallocate(iceage12z_0p25deg)

      ! EMK...When downscaling from 0.25 deg, some snow-free gaps appear at 
      ! high-latitudes that cannot be corrected by the snow depth analysis 
      ! due to a lack of observations.  As a kludge, we check for this and 
      ! replace with climo values.
      do r = 1,nr
         do c = 1, nc
            gindex = c+(r-1)*LDT_rc%lnc(1)
            rlat = LDT_domain(1)%lat(gindex)
            if (rlat < 40.0 .and. rlat > -40.0) cycle
            if (USAFSI_arrays%olddep(c,r) .le. 0 .and. &
                 USAFSI_arrays%climo(c,r) > 0) then
               USAFSI_arrays%olddep(c,r) = USAFSI_arrays%climo(c,r)
               USAFSI_arrays%oldmask(c,r) = 1
               if (landice(c,r) > 0.5) then
                  USAFSI_arrays%snoage(c,r) = maxage
                  USAFSI_arrays%snoage12z(c,r) = maxage                
               else
                  USAFSI_arrays%snoage(c,r) = 1
                  USAFSI_arrays%snoage12z(c,r) = 1
               endif
            end if

         end do ! c
      end do ! r

      ! DONE....RETURN TO CALLING ROUTINE
      return

      ! ERROR HANDLING SECTION.
4200  continue
      message(1) = '[ERR] ERROR CONVERTING DATA FROM CHARACTER TO INTEGER'
      message(2) = '[ERR] DATE10 = ' // date10
      call abort_message (program_name, program_name, message)
      call LDT_endrun()

      ! FORMAT STATEMENTS.
6000  format (/, '[INFO]', A, ': MODIFIED DATA FROM ', A10, ' MISSING; ',  &
           'USING UNMODIFIED DATA')
6200  format (/, '[INFO]', A, ': DATA FROM ', A10, ' MISSING; ',           &
           'CHECKING FOR PREVIOUS CYCLE')
   end subroutine getsno

   ! Fetch USAFSI data from netCDF
   subroutine getsno_nc(date10, julhr_beg, ierr)

      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      use USAFSI_netcdfMod, only: USAFSI_read_netcdf, &
           USAFSI_read_netcdf_12z
      use USAFSI_utilMod, only: abort_message, date10_julhr, &
           julhr_date10
      use USAFSI_paramsMod, only: msglns, program_name

      ! Defaults
      implicit none

      ! Arguments
      character*10, intent(in) :: date10
      integer, intent(out) :: julhr_beg
      integer, intent(out) :: ierr

      ! Local variables
      logical :: found, found_12z
      integer :: limit, tries
      integer :: runcycle
      integer :: julhr
      character*12 :: routine_name
      character*10 :: date10_prev
      character*90 :: message(msglns)

      data routine_name / 'GETSNO_NC   '/

      ! Find the date/time group of the previous cycle
      found = .false.
      found_12z = .false.
      limit = 20
      tries = 1

      call date10_julhr(date10, julhr, program_name, routine_name)
      julhr_beg = julhr

      ! Grab prior analysis
      do while ((.not. found) .and. (tries .le. limit))
         julhr_beg = julhr_beg - 6
         call julhr_date10(julhr_beg, date10_prev, program_name, routine_name)
         call USAFSI_read_netcdf(date10_prev,ierr)
         if (ierr == 0) then
            found = .true.
         else
            write (ldt_logunit,6200) trim (routine_name), date10_prev
            tries = tries + 1
         end if
      end do

      ! If 12Z cycle, retrieve last 12Z snow and ice age
      read (date10(9:10), '(i2)', err=4200) runcycle
      if (found .and. runcycle .eq. 12) then
         found_12z = .false.
         limit = 5
         tries = 1
         do while ((.not. found_12z) .and. (tries .le. limit))
            julhr = julhr - 24
            call julhr_date10 (julhr, date10_prev, program_name, routine_name)
            call USAFSI_read_netcdf_12z(date10_prev,ierr)
            if (ierr == 0) then
               found_12z = .true.
            else
               write (ldt_logunit,6200) trim (routine_name), date10_prev
               tries = tries + 1
            end if
         end do
      end if

      ! If netCDF data could not be retrieved, pass an error code back
      ! to the caller.  Caller will need to try the legacy SNODEP instead.
      if (.not. found) then
         write(LDT_logunit,*) &
              '[WARN] Cannot find prior USAFSI analysis'
         ierr = 2
      else if (.not. found_12z) then
      !EMK The above else if looks wrong, but using the below, commented-out
      !else if changes the answer.  For now, use the above.
      !else if (.not. found_12z .and. runcycle .eq. 12) then
         write(LDT_logunit,*) &
              '[WARN] Cannot find prior 12Z USAFSI analysis'
         ierr = 1
      else
         ierr = 0
      end if
      return

      ! Error handling section
4200  continue
      message(1) = '[ERR] ERROR CONVERTING DATA FROM CHARACTER TO INTEGER'
      message(2) = '[ERR] DATE10 = ' // date10
      call abort_message (program_name, program_name, message)
      call LDT_endrun()

      ! Other format statements
6200  format (/, '[INFO]', A, ': DATA FROM ', A10, ' MISSING; ',           &
           'CHECKING FOR PREVIOUS CYCLE')

   end subroutine getsno_nc

   subroutine getsst (date10, stmpdir, sstdir)

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME:  GET SEA SURFACE TEMPERATURE (SST) DATA
      !**
      !**  PURPOSE:  GET NAVY COUPLED OCEAN DATA ASSIMILATION (NCODA) SST DATA
      !**
      !**  CALLED FROM:  SNODEP
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  DATE10, SSTDIR, STMPDIR
      !**
      !**   OUTPUT: SST
      !**
      !**  FILES ACCESSED
      !**  ==============
      !**  FILE NAME                R/W             DESCRIPTION
      !**  ------------------------ --- --------------------------------------
      !**  sst.dat                   R  NAVY SEA SURFACE TEMPS
      !**
      !**   See USAFSI_arraysMod.F90 for common array descriptions
      !**   See USAFSI_paramsMod.F90 for common parameter descriptions
      !**
      !**  UPDATES
      !**  =======
      !**  18 JUN 01  INITIAL VERSION........................MR GAYNO/SSGT MILLER/DNXM
      !**  05 APR 05  CONVERTED TO FORTRAN 90 FREE FORM...............MR LEWISTON/DNXM
      !**  24 FEB 10  MOVED SST ARRAY TO MODULE; ADDED SSTDIR.....MR LEWISTON/16WS/WXE
      !**  31 MAY 12  CORRECTED ORDER OF PUTGET ARGUMENTS.........MR LEWISTON/16WS/WXE
      !**  14 DEC 12  MOVED DEGRIB PROCESS HERE FROM SCRIPT.......MR LEWISTON/16WS/WXE
      !**  13 MAY 14  CHANGED SNODEP SST RETRIEVAL PATH.............MR WRIGHT/16WS/WXE
      !**  18 SEP 14  REMOVED TRIM FROM PATH_GRIB IN COPEN CALL...MR LEWISTON/16WS/WXE
      !**  03 APR 15  UPDATED INPUT GRIB FILE NAME FORMAT.........MR LEWISTON/16WS/WXE
      !**  21 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
      !**  13 Jan 21  Added FNMOC GRIB1 file...Eric Kemp, NASA GSFC/SSAI
      !**
      !*******************************************************************************
      !*******************************************************************************

      ! Imports
      use LDT_coreMod, only: LDT_rc, LDT_domain
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      use map_utils
      use USAFSI_arraysMod, only : USAFSI_arrays
      use USAFSI_paramsMod
      use USAFSI_utilMod ! EMK

      ! Defaults
      implicit none

      ! Arguments
      character*10,  intent(in)   :: date10           ! SNODEP DATE-TIME GROUP
      character*255, intent(in)   :: stmpdir          ! SFC TEMPERATURE DIRECTORY PATH
      character*255, intent(in)   :: sstdir

      ! Local constants
      integer, parameter          :: sst_size = sst_igrid * sst_jgrid  ! SST ARRAY SIZE

      ! Local variables
      character*10                :: date10_sst       ! SST DATE-TIME GROUP
      character*255               :: file_binary      ! FULLY-QUALIFIED BINARY NAME
      character*7                 :: iofunc           ! ACTION TO BE PERFORMED
      !character*90                :: message (msglns) ! ERROR MESSAGE
      character*255                :: message (msglns) ! ERROR MESSAGE
      character*12                :: routine_name     ! NAME OF THIS SUBROUTINE
      integer                     :: runcycle         ! CYCLE TIME
      integer                     :: hrdiff           ! DIFFERENCE BETWEEN HOURS
      integer                     :: julsno           ! JULIAN HOUR OF SNODEP CYCLE
      integer                     :: julsst           ! JULIAN HOUR OF SST DATA
      integer                     :: limit            ! LIMIT OF TRIES
      integer                     :: tries            ! LIMIT OF TRIES
      logical                     :: isfile           ! FLAG INDICATING WHETHER FILE FOUND
      logical                     :: found            ! FLAG INDICATING WHETHER FILE FOUND
      real, allocatable :: sst_0p25deg(:,:)
      type(proj_info) :: sst_0p25deg_proj
      integer :: i_0p25deg, j_0p25deg
      integer :: gindex,c,r
      real :: rlat,rlon,ri,rj
      integer :: nc,nr
      character*255 :: file_grib
      integer :: grstat

      data routine_name           / 'GETSST      '/

      ! FIND THE DATE/TIME GROUP OF THE PREVIOUS CYCLE.
      ! GET SEA SURFACE TEMPERATURE DATA.
      iofunc  = 'READING'
      isfile  = .false.
      message = ' '
      found   = .false.
      limit   = 3
      tries   = 1

      call date10_julhr (date10, julsst, program_name, routine_name)

      ! LOOK FOR LATEST 12Z CYCLE.  IF IT'S NOT 12Z THEN DO BACK 6 HOURS
      ! AT A TIME FOR A MAX OF 3 TRIES UNTIL 12Z IS FOUND.
      sst_loop : do while ((mod (julsst , 24) .ne. 12) .and.            &
           (tries .le. limit))
         julsst = julsst - 6
         tries = tries + 1
      end do sst_loop

      if (mod (julsst , 24) .ne. 12) go to 5200

      allocate(sst_0p25deg(sst_igrid,sst_jgrid))
      sst_0p25deg     = -1.0

      tries = 1
      limit = 7 ! EMK Check previous 7 days

      ! LOOK FOR DEGRIBBED SST BINARY.  IF NOT FOUND LOOK IN DIFFERENT
      ! DIRECTORY FOR GR1 FILE, DEGRIB, READ, AND WRITE OUT BINARY AFTER
      ! A CHECK TO SEE IF BINARY IS ALREADY THERE FROM ANOTHER INSTANCE
      ! OF SNODEP RUNNING THIS PROGRAM.
      cycle_loop : do while ((.not. found) .and. (tries .le. limit))
         tries = tries + 1
         call julhr_date10 (julsst, date10_sst, program_name, routine_name)
         read (date10_sst(9:10), '(i2)', err=5000) runcycle
         file_binary = trim (stmpdir) // 'navyssts' // meshname // '.' &
              // date10_sst // '.dat'
         inquire (file=file_binary, exist=isfile)
         if (isfile) then
            found =.true.
            write(ldt_logunit,*)'[INFO] Reading ', trim(file_binary)
            write (ldt_logunit, 6000) routine_name, iofunc, trim (file_binary)
            call putget_real ( sst_0p25deg, 'r', file_binary, program_name,  &
                 routine_name, sst_igrid, sst_jgrid )
         else
            write(ldt_logunit,*)'[WARN] Cannot find ', trim(file_binary)
            !EMK 20220113...Reinstated GRIB1 support
            file_grib = trim(sstdir) &
                 // 'US058GOCN-GR1mdl.0043_0200_00000A0LT' &
                 // date10_sst &
                 // '_0160_000000-000000sea_temp.gr1'
            inquire(file=file_grib, exist=isfile)
            if (isfile) then
               call read_grib1_sst(file_grib, sst_igrid, sst_jgrid, &
                    sst_0p25deg, grstat)
               if (grstat .eq. 0) then
                  found = .true.
                  file_binary = trim(stmpdir) &
                       // 'navyssts' &
                       // meshname &
                       // '.' &
                       // date10_sst &
                       // '.dat'
                  inquire(file=file_binary, exist=isfile)
                  if (.not. isfile) then
                     iofunc = '[INFO] WRITING'
                     write(ldt_logunit, 6000) routine_name, iofunc, &
                          trim(file_binary)
                     call putget_real(sst_0p25deg, 'w', file_binary, &
                          program_name, routine_name, sst_igrid, sst_jgrid)
                  end if
               else
                  message(1) = '[ERR] ERROR READING FILE'
                  message(2) = '[ERR] PATH = ' // file_grib
                  call error_message(program_name, routine_name, message)
                  write(ldt_logunit, 6400) routine_name, iofunc, file_grib, &
                       grstat
               end if
            else
               message(1) = '[ERR] ERROR OPENING FILE'
               message(2) = '[ERR] PATH = ' // file_grib
               call error_message(program_name, routine_name, message)
               write(ldt_logunit, 6400) routine_name, iofunc, file_grib, grstat
            end if
         end if
         julsst = julsst - 24
      end do cycle_loop

      if (found) then
         call date10_julhr (date10, julsno, program_name, routine_name)
         if (julsno .ge. julsst) then
            hrdiff = julsno - julsst

            ! THIS CHECKS BACK 48 HOURS DUE TO THE "JULSST = JULSST - 24"
            ! ABOVE WHERE THE PROGRAM IS LOOKING BACK 24 HOURS IN CASE OF
            ! MISSING FILES.  BY SEARCHING BACK 48 HOURS THE ERROR IS
            ! REALLY LOOKING BACK ONLY 24 HOURS (IF THAT MAKES SENSE...).
            if (hrdiff .gt. 48) then
               message(1) = '  SST DATA IS MORE THAN 24 HOURS OLD'
               message(2) = '  USAFSI CYCLE = ' // date10
               message(3) = '  SEA SFC TEMP = ' // date10_sst
               call error_message (program_name, routine_name, message)

            end if

         end if

      else

         message(1) = '[WARN] SEA SURFACE TEMPERATURE DATA NOT FOUND'
         call error_message (program_name, routine_name, message)
         write (ldt_logunit, 6600) routine_name

      end if

      ! Interpolate the SST to the LDT grid
      call map_set(proj_code=proj_latlon, &
           lat1=-90., &
           lon1=0., &
           dx=0.25, &
           stdlon=0.25, &
           truelat1=0.25, &
           truelat2=0., &
           idim=sst_igrid,&
           jdim=sst_jgrid, &
           proj=sst_0p25deg_proj)

      USAFSI_arrays%sst(:,:) = -1
      nr = LDT_rc%lnr(1)
      nc = LDT_rc%lnc(1)
      do r = 1,nr
         do c = 1,nc
            gindex = c+(r-1)*nc
            rlat = LDT_domain(1)%lat(gindex)
            rlon = LDT_domain(1)%lon(gindex)
            call latlon_to_ij(sst_0p25deg_proj,rlat,rlon,ri,rj)
            i_0p25deg = nint(ri)
            if (i_0p25deg .gt. igrid) then
               i_0p25deg = i_0p25deg - igrid
            else if (i_0p25deg .lt. 1) then
               i_0p25deg = i_0p25deg + igrid
            end if
            j_0p25deg = nint(rj)
            if (j_0p25deg .lt. 1) then
               j_0p25deg = 1
            else if (j_0p25deg .gt. jgrid) then
               j_0p25deg = jgrid
            end if
            if (sst_0p25deg(i_0p25deg,j_0p25deg) < -1) then
               USAFSI_arrays%sst(c,r) = -1
            else
               USAFSI_arrays%sst(c,r) = sst_0p25deg(i_0p25deg,j_0p25deg)
            end if
         end do ! c
      end do ! r
      
      ! Clean up
      deallocate(sst_0p25deg)

      return

      ! ERROR-HANDLING SECTION.
5000  continue
      message(1) = '[ERR] ERROR CONVERTING INTEGER DATE/TIMES TO CHARACTER'
      message(2) = '[ERR] DATE10 = ' // date10_sst
      call abort_message (program_name, routine_name, message)
      if (allocated(sst_0p25deg)) deallocate(sst_0p25deg)
      return

5200  continue
      message(1) = '[ERR] INVALID CYCLE/SOMETHING IS SERIOUSLY WRONG'
      message(2) = '[ERR] DATE10 = ' // date10
      call abort_message (program_name, routine_name, message)
      if (allocated(sst_0p25deg)) deallocate(sst_0p25deg)
      return

      ! Format statements
6000  format (/, '[INFO] ', A6, ': ', A7, 1X, A)
!6200  format (/, '[INFO] ', A6, ': CURRENT SEA SURFACE TEMPERATURE DTG = ', &
!           A10)
6400  format (/, '[WARN] ', A6, ': ERROR ', A7, 1X, A, /, 3X, 'STATUS = ', I6)
6600  format (/, '[WARN] ', A6, ': SEA SURFACE TEMPERATURE DATA NOT FOUND')

   end subroutine getsst


   subroutine getviirs (date10, viirsdir)

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME: GETVIIRS
      !**
      !**  PURPOSE: READ IN VIIRS SNOW COVERED AREA MAP
      !**
      !**  CALLED FROM: SNODEP
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  DATE10, VIIRSDIR, MAXPIXAGE, MINFRAC, MINBARE
      !**
      !**   OUTPUT: VIIRSMAP, USEVIIRS
      !**
      !**  FILES ACCESSED
      !**  ==============
      !**  FILE NAME                                  R/W DESCRIPTION
      !**  ----------------------------------------   --- ----------------------------
      !**  ${VIIRSDIR}/snomap_0p005deg.${DATEFR}.tiff  R  VIIRS SNOW MAP
      !**  ${VIIRSDIR}/snoage_0p005deg.${DATEFR}.tiff  R  VIIRS SNOW AGE
      !**
      !**  UPDATES
      !**  =======
      !**  15 FEB 17  ADAPTED FROM GETFRAC ......................... M PUSKAR/16WS/WXE
      !**  22 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
      !**
      !*******************************************************************************
      !*******************************************************************************

      ! Imports
      use LDT_coreMod, only: LDT_rc, LDT_domain
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      use LDT_usafsiMod, only: usafsi_settings
      use map_utils
      use USAFSI_arraysMod, only: USAFSI_arrays
      use USAFSI_paramsMod
      use USAFSI_utilMod ! EMK

      ! Defaults
      implicit none

      ! Argments
      character(10), intent(in)   :: date10           ! DATE-TIME GROUP OF SNODEP CYCLE
      character(255), intent(in)  :: viirsdir         ! FRACTIONAL SNOW DIRECTORY PATH

      ! Local variables
      character(2)                :: cyclhr           ! CYCLE HOUR

      character(10)               :: datefr           ! DATE-TIME GROUP OF SNOW COVER
      character(255)              :: snomap_path      ! FULLY-QUALIFIED SNOMAP FILE NAME
      character(255)              :: snoage_path      ! FULLY-QUALIFIED SNOAGE FILE NAME
      character(7)                :: iofunc           ! ACTION TO BE PERFORMED
      character(90)               :: message (msglns) ! ERROR MESSAGE
      character(12)               :: routine_name     ! NAME OF THIS SUBROUTINE
      integer                     :: i                ! SNODEP I-COORDINATE
      integer                     :: icount           ! LOOP COUNTER
      integer                     :: julhr            ! AFWA JULIAN HOUR
      integer                     :: j                ! SNODEP J-COORDINATE
      logical                     :: map_exists       ! FLAG INDICATING WHETHER SNOMAP FOUND
      logical                     :: age_exists       ! FLAG INDICATING WHETHER SNOAGE FOUND
      integer :: nc,nr
      integer :: igrid_viirs, jgrid_viirs
      real :: rlat, rlon, ri, rj, ri_viirs, rj_viirs
      integer :: i_viirs, j_viirs, ierr
      integer, allocatable :: snow(:,:)
      integer, allocatable :: bare(:,:)
      integer, allocatable :: pixels(:,:)
      integer, allocatable :: mapbuf_slice(:)
      integer, allocatable :: agebuf_slice(:)
      type(proj_info) :: viirs_0p005deg_proj
      real :: min_snow_count, min_bare_count

! EMK...Must have LIBGEOTIFF support enabled!
#ifndef USE_LIBGEOTIFF
      write(LDT_logunit,*)'[ERR] LDT not compiled with LIBGEOTIFF support!'
      write(LDT_logunit,*)'[ERR] Cannot read VIIRS file!'
      write(LDT_logunit,*)&
           '[ERR] Recompile with LIBGEOTIFF support and try again!'
      call LDT_endrun()
      
#else
      external :: ztif_frac_slice ! EMK 20220113

      data routine_name  / 'GETVIIRS    ' /

      ! ALLOCATE DATA ARRAYS.
      nc = LDT_rc%lnc(1)
      nr = LDT_rc%lnr(1)
      igrid_viirs = 72000
      jgrid_viirs = 36000
      allocate(mapbuf_slice(igrid_viirs))
      mapbuf_slice(:) = 255
      allocate(agebuf_slice(igrid_viirs))
      agebuf_slice(:) = 0
      ierr = 0

      ! NOTE:  VIIRS data is upside down (origin is upper left, not lower left)
      call map_set(proj_code=proj_latlon, &
           lat1 =   89.9975, &
           lon1 = -179.9975, &
           dx = 0.005, &
           stdlon = 0.005, &    ! dlon
           truelat1 = -0.005, & ! dlat, decreases as j increases
           truelat2 = 0., &
           idim = igrid_viirs, &
           jdim = jgrid_viirs, &
           proj=viirs_0p005deg_proj)
      
      ! INITIALIZE VARIABLES.
      icount = 0
      iofunc   = 'READING'
      map_exists   = .false.
      age_exists   = .false.
      message  = ' '
      USAFSI_arrays%viirsmap = 255

      ! RETRIEVE VIIRS SNOW DATA.
      cyclhr = date10 (9:10)
      datefr = date10 (1:8) // '00'

      if (cyclhr .eq. '00' .or. cyclhr .eq. '06') then
         call date10_julhr (datefr, julhr, program_name, routine_name)
         julhr  = julhr  - 24
         call julhr_date10 (julhr, datefr, program_name, routine_name)

      end if

      file_search : do while ((.not. map_exists) .and. &
           (.not. age_exists) .and. &
           (icount .le. usafsi_settings%maxpixage))

         snomap_path = trim (viirsdir) // '/snomap_0p005deg.' // datefr &
              // '.tiff'
         snoage_path = trim (viirsdir) // '/snoage_0p005deg.' // datefr &
              // '.tiff'

         inquire (file=snomap_path, exist=map_exists)
         inquire (file=snoage_path, exist=age_exists)

         if (.not. map_exists) then
            write(LDT_logunit,*) &
                 '[WARN] Cannot find ',trim(snomap_path)
         end if
         if (.not. age_exists) then
            write(LDT_logunit,*) &
                 '[WARN] Cannot find ',trim(snoage_path)
         end if

         if (map_exists .and. age_exists) then

            write (ldt_logunit, 6000) &
                 routine_name, iofunc, trim(snomap_path)
            write (ldt_logunit, 6000) &
                 routine_name, iofunc, trim(snoage_path)

            allocate(snow(nc,nr))
            allocate(bare(nc,nr))
            allocate(pixels(nc,nr))
            snow(:,:) = 0
            bare(:,:) = 0
            pixels(:,:) = 0

            ! Read the VIIRS data at native resolution one slice at a time.
            ! For each slice, geolocate onto the LDT grid and identify
            ! as snow or bare.
            do j_viirs = 1, jgrid_viirs

               ierr = 0
               call ztif_frac_slice(mapbuf_slice, &
                    agebuf_slice, &
                    igrid_viirs, &
                    jgrid_viirs, &
                    trim(adjustl(snomap_path))//char(0), &
                    trim(adjustl(snoage_path))//char(0), &
                    icount, &
                    j_viirs, &
                    ierr)
               ! Handle errors
               if (ierr .eq. 1) then
                  write(LDT_logunit,*) &
                       '[WARN] Problem reading ',trim(snomap_path)
                  exit ! Get out of loop
               else if (ierr .eq. 2) then
                  write(LDT_logunit,*) &
                       '[WARN] Problem reading ',trim(snoage_path)
                  exit ! Get out of loop
               else

                  ! No error for this slice, so process
                  do i_viirs = 1, igrid_viirs
               
                     ! Find lat/lon of VIIRS pixel, and then determine which
                     ! LDT grid box this falls in.
                     ri_viirs = real(i_viirs)
                     rj_viirs = real(j_viirs)
                     call ij_to_latlon(viirs_0p005deg_proj, &
                          ri_viirs, &
                          rj_viirs, &
                          rlat, &
                          rlon)
                     call latlon_to_ij(LDT_domain(1)%ldtproj, &
                          rlat, &
                          rlon, &
                          ri, &
                          rj)
                     i = nint(ri)
                     if (i < 1) then
                        i = i + nc
                     else if (i > nc) then
                        i = i - nc
                     end if
                     j = nint(rj)
                     if (j < 1) then
                        j = 1
                     else if (j > nr) then
                        j = nr
                     end if                     
                     pixels(i,j) = pixels(i,j) + 1

                     ! Skip if the pixel age is too old.
                     if (agebuf_slice(i_viirs) > &
                          usafsi_settings%maxpixage) cycle
               
                     ! Increment the appropriate snow/bare counter
                     if (mapbuf_slice(i_viirs) .eq. 0) then
                        bare(i,j) = bare(i,j) + 1
                     else if (mapbuf_slice(i_viirs) .eq. 1) then
                        snow(i,j) = snow(i,j) + 1
                     end if
                  end do ! i_viirs
               end if
            end do ! j_viirs
         end if

         if (.not. map_exists .or. .not. age_exists .or. ierr .ne. 0) then

            if (allocated(snow)) deallocate(snow)
            if (allocated(bare)) deallocate(bare)
            if (allocated(pixels)) deallocate(pixels)

            call date10_julhr (datefr, julhr, program_name, routine_name)

            julhr  = julhr  - 24
            icount = icount + 1
            call julhr_date10 (julhr, datefr, program_name, routine_name)

         end if

      end do file_search

      if (map_exists .and. age_exists .and. ierr .eq. 0) then         

         ! From the geolocated data, create the final VIIRS snow cover map
         do j = 1, nr
            do i = 1, nc
               min_snow_count = usafsi_settings%minfrac * pixels(i,j)
               min_bare_count = usafsi_settings%minbare * pixels(i,j)
               if (snow(i,j) > min_snow_count) then
                  usafsi_arrays%viirsmap(i,j) = 1
               else if (bare(i,j) > min_bare_count) then
                  usafsi_arrays%viirsmap(i,j) = 0
               else
                  usafsi_arrays%viirsmap(i,j) = 255
               end if
            end do ! i
         end do ! j
         deallocate(snow)
         deallocate(bare)
         deallocate(pixels)

      else

         if (.not. map_exists) then
            usafsi_settings%useviirs = .false.
            message(1) = '[WARN] VIIRS SNOW MAP FILE NOT FOUND'
            !message(2) = '[WARN] PATH = ' // trim(snomap_path)
            call error_message (program_name, routine_name, message)
            write (ldt_logunit, 6400) routine_name, snomap_path
         end if

         if (.not. age_exists) then
            usafsi_settings%useviirs = .false.
            message(1) = '[WARN] VIIRS SNOW AGE FILE NOT FOUND'
            !message(2) = '[WARN] PATH = ' // trim(snoage_path)
            call error_message (program_name, routine_name, message)
            write (ldt_logunit, 6400) routine_name, snoage_path
         end if

      end if

      ! DEALLOCATE ARRAYS.
      deallocate(agebuf_slice)
      deallocate(mapbuf_slice)

      return

      ! FORMAT STATEMENTS.
6000  format (/, '[INFO] ', A8, ': ', A7, 1X, A, A)
6400  format (/, '[WARN] ', A8, ': FILE NOT FOUND: ', A,                    &
           /, '   WILL NOT USE VIIRS SNOW DATA')

!EMK Must have LIBGEOTIFF support enabled
#endif
   end subroutine getviirs

   ! EMK New snow analysis excluding glaciers
   subroutine run_snow_analysis_noglacier(runcycle, nc, nr, landmask, &
        landice, &
        elevations, sfctmp_found, sfctmp_lis, bratseth)

      ! Imports
      use LDT_bratsethMod
      use LDT_coreMod, only: LDT_domain
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      use LDT_usafsiMod, only: usafsi_settings
      use map_utils
      use USAFSI_arraysMod, only: USAFSI_arrays
      use USAFSI_paramsMod

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: runcycle
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: landmask(nc,nr)
      real, intent(in) :: landice(nc,nr)
      real, intent(in) :: elevations(nc,nr)
      logical, intent(in) :: sfctmp_found
      real, intent(in) :: sfctmp_lis(:,:)
      type(LDT_bratseth_t), intent(inout) :: bratseth

      ! Local variables
      real :: satdep
      real, allocatable :: gbacks(:,:)
      integer, allocatable :: snomask(:,:)
      logical, allocatable :: skip_grid_points(:,:)
      logical, allocatable :: updated(:,:)
      real :: ri,rj,lat,lon
      real :: converg_thresh
      real :: rthresh
      integer :: job,num_good_obs
      integer :: c,r
      integer :: snomask_reject_count
      integer :: bad_back_count, glacier_zone_count
      real :: ob_value
      character*32 :: new_name
      integer :: gindex
      real :: rlat

      ! Initializations
      allocate(snomask(nc,nr))
      allocate(skip_grid_points(nc,nr))
      allocate(updated(nc,nr))
      snomask(:,:) = -1 ! Unknown
      skip_grid_points(:,:) = .false.
      updated(:,:) = .false.
      rthresh = float(usafsi_settings%thresh) / 10.0

      ! Build mask to screen out areas we don't want to analyze.
      do r = 1,nr
         do c = 1,nc
            if (landmask(c,r) < 0.5) then
               skip_grid_points(c,r) = .true.
               ! EMK...Make sure prior analysis has no snow here
               USAFSI_arrays%snoanl(c,r) = -1
               USAFSI_arrays%snoage(c,r) = -1               
            else if (landice(c,r) > 0.5) then
               ! Use LDT landcover to identify glaciers
               skip_grid_points(c,r) = .true.
            else if (USAFSI_arrays%snow_poss(c,r) == 0) then
               ! Since the snow_poss mask is coarse resolution (0.25 deg),
               ! there may be artificial snow-free gaps introduced when 
               ! downscaling.  So, we will allow snow analysis to run in high 
               ! latitudes for non-glacier regions. Low and mid-latitudes
               ! will still be screened out.
               gindex = c+(r-1)*nc
               rlat = LDT_domain(1)%lat(gindex)
               if (rlat > -40.0 .and. rlat < 40.0) then
                  skip_grid_points(c,r) = .true.
                  ! No snow allowed, so zero out these fields now
                  USAFSI_arrays%snoanl(c,r) = 0
                  USAFSI_arrays%snoage(c,r) = 0               
               end if
            end if
         end do ! c
      end do ! r

      ! Initialize the analysis over nonglacier land points
      do r = 1, nr
         do c = 1, nc
            if (skip_grid_points(c,r)) cycle
            USAFSI_arrays%snoanl(c,r) = USAFSI_arrays%olddep(c,r)
         end do ! c
      end do ! r

      ! Adjust the first-guess towards gridded SSMIS.
      do r = 1, nr
         do c = 1, nc

            if (skip_grid_points(c,r)) cycle

            ! Get SSMIS value
            satdep = misanl 
            if (USAFSI_arrays%ssmis_depth(c,r) >= 0) then
               satdep = USAFSI_arrays%ssmis_depth(c,r)
            end if

            ! Handle SSMIS detection problem with very shallow snow.  If below
            ! minimum depth, preserve prior analysis if larger than SSMIS.
            if (satdep >= 0) then
               if (satdep .le. usafsi_settings%minsat) then
                  if (satdep > USAFSI_arrays%olddep(c,r)) then
                     USAFSI_arrays%snoanl(c,r) = satdep
                     updated(c,r) = .true.
                     if (satdep > 0) then
                        snomask(c,r) = 1
                     end if
                  end if
               else
                  USAFSI_arrays%snoanl(c,r) = satdep
                  updated(c,r) = .true.
                  if (satdep > 0) then
                     snomask(c,r) = 1
                  end if
               end if
            end if
            
         end do ! c
      end do ! r

      ! Update the snow mask using the fractional snow
      if (usafsi_settings%usefrac) then
         do r = 1, nr
            do c = 1, nc     
               if (skip_grid_points(c,r)) cycle
               if (USAFSI_arrays%snofrac(c,r) > usafsi_settings%minfrac) then
                  snomask(c,r) = 1
               end if
            end do ! c
         end do ! r
      end if ! usefrac

      ! Adjust snow mask based on VIIRS
      if (usafsi_settings%useviirs) then
         do r = 1,nr
            do c = 1, nc
               if (skip_grid_points(c,r)) cycle
               if (USAFSI_arrays%viirsmap(c,r) == 0) then
                  snomask(c,r) = 0
               else if (USAFSI_arrays%viirsmap(c,r) == 1) then
                  snomask(c,r) = 1
               end if
            end do ! c
         end do ! r
      end if

      ! Apply temperature check for snow.  Assumes rthresh is high
      ! enough that any VIIRS or MODIS snow are false detections.
      if (sfctmp_found) then
         do r = 1,nr
            do c = 1,nc
               if (skip_grid_points(c,r)) cycle
               if (sfctmp_lis(c,r) > rthresh) then
                  snomask(c,r) = 0
               end if
            end do ! c
         end do ! r
      end if

      ! Adjust snow depth analysis towards snow mask.  Also keep track of
      ! where snow depth was changed.
      do r = 1,nr
         do c = 1,nc
            if (skip_grid_points(c,r)) cycle
            if (snomask(c,r) .eq. 0 .and. &
                 USAFSI_arrays%snoanl(c,r) > 0) then
               USAFSI_arrays%snoanl(c,r) = 0
               USAFSI_arrays%snoage(c,r) = 0
            end if
            if (USAFSI_arrays%snoanl(c,r) .ne. USAFSI_arrays%olddep(c,r)) then
               updated(c,r) = .true.
            end if
         end do ! c
      end do ! r

      ! At this point, we have our background field and snow mask.
      ! Start QC of surface observations.

      write(LDT_logunit,*) &
           '[INFO] Reject obs that are missing elevations'
      call bratseth%run_missing_elev_qc()

      write(LDT_logunit,*)'[INFO] Reject obs over water'
      call bratseth%run_water_qc(1,nc,nr,landmask)

      write(LDT_logunit,*) &
           '[INFO] Reject obs where snow is never permitted'
      call bratseth%run_nosnow_qc(nc,nr,USAFSI_arrays%snow_poss)

      write(LDT_logunit,*) &
           '[INFO] Reject obs with elevations too different from LDT'
      call bratseth%run_elev_qc(1,nc,nr,elevations,&
           usafsi_settings%elevqc_diff_threshold)
      
      write(LDT_logunit,*)'[INFO] Check for duplicate obs'
      call bratseth%run_dup_qc()
                  
      ! Interpolate to the observations.  For simplicity, just use the value 
      ! in the grid box.
      num_good_obs = bratseth%count_good_obs()
      snomask_reject_count =0
      bad_back_count = 0
      glacier_zone_count = 0
      write(LDT_logunit,*)&
           '[INFO] Assigning background values for ',num_good_obs, &
           ' observations'
      do job = 1, num_good_obs
         call bratseth%get_lat_lon(job,lat,lon)
         call latlon_to_ij(LDT_domain(1)%ldtproj,lat,lon,ri,rj)
         c = nint(ri)
         if (c > nc) then
            c = c - nc
         else if (c < 1) then
            c = c + nc
         end if
         r = min(nr,max(1,nint(rj)))
         ! Do not use obs that are in glacier zones.
         if (skip_grid_points(c,r)) then
            !write(LDT_logunit,*) &
            !     '[INFO] Bad background, job, lat, lon, c, r, back: ', &
            !     job, lat, lon, c, r, USAFSI_arrays%snoanl(c,r)
            !write(LDT_logunit,*)'EMK: REJECT job, back: ', &
            !     job, USAFSI_arrays%snoanl(c,r)
            glacier_zone_count = glacier_zone_count + 1
            call bratseth%reject_ob(job, "GLACIER ZONE REJECT")
         else if (USAFSI_arrays%snoanl(c,r) < 0) then
            bad_back_count = bad_back_count + 1
            call bratseth%reject_ob(job, "BAD BACKGROUND REJECT")
         else if (snomask(c,r) == 0) then
            ! Reject obs that have snomask of 0, unless the ob agrees
            call bratseth%get_ob_value(job,ob_value)
            if (ob_value > 0) then
               !write(LDT_logunit,*) &
               !     '[INFO] Ob with zero snomask, job, lat, lon, ', &
               !     'c, r, back: ', &
               !     job, lat, lon, c, r, USAFSI_arrays%snoanl(c,r)
               snomask_reject_count = snomask_reject_count + 1
               call bratseth%reject_ob(job,"SAT AND T SNOMASK REJECT")
            else
               call bratseth%set_back(job,USAFSI_arrays%snoanl(c,r))
            end if
         else 
            call bratseth%set_back(job,USAFSI_arrays%snoanl(c,r))
         end if
      end do ! j
      write(LDT_logunit,*) '[INFO] Rejected ',glacier_zone_count, &
           ' obs in glacier zones'
      write(LDT_logunit,*) '[INFO] Rejected ',bad_back_count, &
           ' obs due to missing background'
      write(LDT_logunit,*) '[INFO] Rejected ',snomask_reject_count, &
           ' obs based on SAT and T snow mask'
      if ((snomask_reject_count + bad_back_count) > 0) then
         call bratseth%resort_bad_obs()
      end if

      ! Reject obs that are too low compared to the background.  This
      ! check is from Brasnett (1999).
      num_good_obs = bratseth%count_good_obs()
      write(LDT_logunit,*) &
           '[INFO] Reject obs with snow depths too low compared to background'
      call bratseth%run_skewed_back_qc(usafsi_settings%skewed_backqc_threshold)

      ! Reject obs that differ "too much" from background.  
      num_good_obs = bratseth%count_good_obs()
      write(LDT_logunit,*) &
           '[INFO] Reject obs that differ too much from background'
      call bratseth%run_back_qc(usafsi_settings%back_err_var)

      write(LDT_logunit,*) &
           '[INFO] Merge obs in the same LDT grid box'
      new_name = 'SUPEROB'
      call bratseth%run_superstat_qc(1,new_name,nc,nr)

      ! Finalize QC comments
      call bratseth%mark_good_obs()

      ! Run Bratseth scheme to analyze snow depth from ground stations.
      num_good_obs = bratseth%count_good_obs()
      if (num_good_obs > 0) then

         ! We need separate arrays for input background field and output
         ! analysis
         allocate(gbacks(nc,nr))
         gbacks(:,:) = USAFSI_arrays%snoanl(:,:)

         ! Calculate the inverse data density for each observation.
         call bratseth%calc_inv_data_dens()

         ! Calculate the observation analysis
         converg_thresh = 0.01
         call bratseth%calc_ob_anas(converg_thresh)

         ! Now interpolate to the analysis grid
         call bratseth%calc_grid_anas(1, nc, nr, gbacks, elevations, &
              USAFSI_arrays%snoanl, skip_grid_points)

         ! Clean up
         deallocate(gbacks)
      end if

      ! Remove any spurious values introducted by Bratseth scheme.  Also,
      ! keep track of changes from prior analysis.
      do r = 1,nr
         do c = 1,nc
            if (skip_grid_points(c,r)) cycle
            ! EMK...Clear out snow depth inserted where snow cover is zero.
            if (snomask(c,r) == 0) then
               USAFSI_arrays%snoanl(c,r) = 0
               USAFSI_arrays%snoage(c,r) = 0
            end if
            if (USAFSI_arrays%snoanl(c,r) < 0.01) then
               USAFSI_arrays%snoanl(c,r) = 0
               ! Leave snoage alone here, since a climatological adjustment
               ! is possible further down if snomask is positive.
            end if
            if (USAFSI_arrays%snoanl(c,r) .ne. USAFSI_arrays%olddep(c,r)) then
               updated(c,r) = .true.
            end if
         end do ! c
      end do ! r

      ! Final adjustments to snow depth and age.
      do r = 1,nr
         do c = 1,nc

            if (skip_grid_points(c,r)) cycle

            if (snomask(c,r) == 1 .and. USAFSI_arrays%snoanl(c,r) == 0) then
               ! Handle case where satellite says snow but none is analyzed.
               ! Start with prior analysis, and apply climo adjustment if 12Z.
               USAFSI_arrays%snoanl(c,r) = USAFSI_arrays%olddep(c,r)
               if (runcycle == 12) then
                  USAFSI_arrays%snoage(c,r) = &
                       min(USAFSI_arrays%snoage(c,r), &
                       USAFSI_arrays%snoage12z(c,r))
                  call appclm(USAFSI_arrays%climo(c,r), &
                       USAFSI_arrays%olddep(c,r), &
                       USAFSI_arrays%snoanl(c,r), &
                       USAFSI_arrays%snoage(c,r))
               end if
               ! If no snow in prior analysis, or if climo adjustment
               ! removes snow, put in bogus value.
               if (USAFSI_arrays%snoanl(c,r) == 0) then
                  USAFSI_arrays%snoanl(c,r) = usafsi_settings%unkdep
                  USAFSI_arrays%snoage(c,r) = 1
               end if
            else if (updated(c,r)) then
               ! Handle cases where snow depth was updated.  This
               ! includes case where satellite says there is snow and
               ! the analysis agrees.
               if (USAFSI_arrays%snoanl(c,r) == 0) then
                  USAFSI_arrays%snoage(c,r) = 0
               else if ((USAFSI_arrays%snoanl(c,r) - &
                    USAFSI_arrays%olddep(c,r)) >= snothresh) then
                  USAFSI_arrays%snoage(c,r) = 1
               else if (runcycle == 12) then
                  USAFSI_arrays%snoage(c,r) = &
                       min((USAFSI_arrays%snoage(c,r)+1), maxage)
               end if
            else
               ! No updates were made from prior analysis.  Permit
               ! adjustment to climatology.
               if (runcycle .eq. 12) then
                  USAFSI_arrays%snoage(c,r) = &
                       min(USAFSI_arrays%snoage(c,r), &
                           USAFSI_arrays%snoage12z(c,r))
                  call appclm(USAFSI_arrays%climo(c,r), &
                       USAFSI_arrays%olddep(c,r), &
                       USAFSI_arrays%snoanl(c,r), &
                       USAFSI_arrays%snoage(c,r))
               end if
            end if
         end do ! c
      end do ! r

      ! Clean up
      deallocate(snomask)
      deallocate(skip_grid_points)
      deallocate(updated)

   end subroutine run_snow_analysis_noglacier

   ! Run snow analysis over glaciers
   subroutine run_snow_analysis_glacier(runcycle, nc, nr, landmask, landice)
      
      ! Imports
      use LDT_usafsiMod, only: usafsi_settings
      use USAFSI_arraysMod, only: USAFSI_arrays
      use USAFSI_paramsMod
      
      ! Defaults
      implicit none
      
      ! Arguments
      integer, intent(in) :: runcycle
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: landmask(nc,nr)
      real, intent(in) :: landice(nc,nr)
      
      ! Local variables
      logical :: glacier
      real :: rthresh
      real :: arctlatr
      integer :: c,r
      
      ! NOTE: Do not overwrite snoanl with olddep, as that will destroy
      ! analysis over non-glaciers
      rthresh = float(usafsi_settings%thresh) / 10.0
      arctlatr = float(arctlat) / 100.0

      ! HANDLE GLACIERS.  Glaciers may occur on land points with
      ! snow_poss==0.  Greenland and Antarctica are likely to have 
      ! glaciers, and a lat/lon check is included to preliminarily put
      ! them in.  VIIRS can further remove or add points.  Once the
      ! glaciers are identified, adjustments to climatology are made at
      ! 12Z, unless the glacier first appearsx in this case, a fake
      ! depth is added.
      do r = 1, nr
         do c = 1,nc

            glacier = .false.

            ! Avoid water
            if (landmask(c,r) < 0.5) then
               ! EMK...Make sure prior analysis has no snow here
               USAFSI_arrays%snoanl(c,r) = -1
               USAFSI_arrays%snoage(c,r) = -1               
               cycle
            end if

            ! Enforce LDT landcover determination of permanent land ice
            if (landice(c,r) > 0.5) then
               glacier = .true.
            else
               cycle
            end if

            ! This is now a possible glacier point.  Initialize the analysis
            ! from the previous value.
            USAFSI_arrays%snoanl(c,r) = USAFSI_arrays%olddep(c,r)            
            ! EMK...In case of landmask disagreements between 0.25 deg
            ! legacy SNODEP and LDT grid
            if (USAFSI_arrays%snoanl(c,r) < 0) then
               USAFSI_arrays%snoanl(c,r) = USAFSI_arrays%climo(c,r)
               USAFSI_arrays%snoage(c,r) = maxage
            end if

            ! Glaciers are adjusted to climatalogy over the long term (at
            ! 12Z).  However, the adjustment only works if the prior analysis
            ! has positive depth.  So we check that first before applying.
            ! If no prior snow depth exists, the glacier is "new"; insert a 
            ! fake value for the current analysis.
            if (runcycle .eq. 12) then
               call appclm(USAFSI_arrays%climo(c,r), &
                    USAFSI_arrays%olddep(c,r), &
                    USAFSI_arrays%snoanl(c,r), &
                    USAFSI_arrays%snoage(c,r))
            end if
            ! It is possible that climo adjustment removes snow.  In that case,
            ! put in a bogus value.
            if (USAFSI_arrays%snoanl(c,r) .le. 0) then
               USAFSI_arrays%snoanl(c,r) = usafsi_settings%unkdep
               USAFSI_arrays%snoage(c,r) = 1
            end if

         end do ! c
      end do ! r

   end subroutine run_snow_analysis_glacier

   ! EMK New sea ice analysis using SSMIS.  Only use of GOFS not available.
   ! This is a refactored version of subroutine SSMICE from original SNODEP
   subroutine run_seaice_analysis_ssmis(month,runcycle,nc,nr,landmask)

      ! Imports
      use LDT_usafsiMod, only: usafsi_settings
      use USAFSI_arraysMod, only: USAFSI_arrays
      use USAFSI_paramsMod

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: month
      integer, intent(in) :: runcycle
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: landmask(nc,nr)

      ! Local variables
      integer :: hemi
      integer :: c,r

      do r = 1, nr
         do c = 1, nc

            ! Only run over water
            if (landmask(c,r) > 0.5) then
               USAFSI_arrays%icecon(c,r) = -1
               USAFSI_arrays%icemask(c,r) = -1
               USAFSI_arrays%iceage(c,r) = -1
               cycle
            end if

            ! Initialize arrays
            USAFSI_arrays%icemask(c,r) = 0
            USAFSI_arrays%icecon(c,r) = 0

            ! Find hemisphere
            if (USAFSI_arrays%ptlat(c,r) .ge. 0) then
               hemi = 1
            else
               hemi = 2
            end if

            if (abs(USAFSI_arrays%ptlat(c,r)) &
                 < usafsi_settings%latchk(month,hemi)) then
               ! Force sea ice removal in low latitudes, a function of
               ! month and hemisphere
               USAFSI_arrays%icemask(c,r) = 0
               USAFSI_arrays%icecon(c,r) = 0
               USAFSI_arrays%iceage(c,r) = 0
               
            else if (abs(USAFSI_arrays%ptlat(c,r)) >= &
                 usafsi_settings%icelat(month,hemi)) then
               ! Force sea ice in high latitudes, a function of month and
               ! hemisphere, unless valid SSMIS says otherwise.  
               if (USAFSI_arrays%ssmis_icecon(c,r) >= 0 .and. &
                    USAFSI_arrays%ssmis_icecon(c,r) .le. 100) then
                  USAFSI_arrays%icecon(c,r) = USAFSI_arrays%ssmis_icecon(c,r)
                  if (USAFSI_arrays%icecon(c,r) > usafsi_settings%minice) then
                     USAFSI_arrays%icemask(c,r) = icepnt
                  else
                     USAFSI_arrays%icemask(c,r) = 0
                     USAFSI_arrays%iceage(c,r) = 0
                  end if
               else
                  ! No SSMI available.  Force sea ice and find a reasonable
                  ! value.
                  USAFSI_arrays%icemask(c,r) = icepnt
                  if (USAFSI_arrays%oldcon(c,r) > usafsi_settings%minice) then
                     USAFSI_arrays%icecon(c,r) = &
                          USAFSI_arrays%oldcon(c,r) ! Prior analysis
                  else
                     USAFSI_arrays%icecon(c,r) = icedef ! Bogus value
                  end if
               end if

            else
               ! Mid-latitude case.  
               if ( (USAFSI_arrays%ssmis_icecon(c,r) > &
                    usafsi_settings%minice) .and. &
                    (USAFSI_arrays%ssmis_icecon(c,r) .le. 100)) then
                  ! Valid SSMIS detected ice
                  USAFSI_arrays%icemask(c,r) = icepnt
                  USAFSI_arrays%icecon(c,r) = USAFSI_arrays%ssmis_icecon(c,r)
               else if (USAFSI_arrays%ssmis_icecon(c,r) >= 0 .and. &
                    USAFSI_arrays%ssmis_icecon(c,r) .le. &
                    usafsi_settings%minice) then
                  ! Valid SSMIS detected no ice
                  USAFSI_arrays%icemask(c,r) = 0
                  USAFSI_arrays%icecon(c,r) = 0
                  USAFSI_arrays%iceage(c,r) = 0
               else
                  ! No valid SSMIS data.  Attempt to use prior analysis
                  if (USAFSI_arrays%oldcon(c,r) > usafsi_settings%minice) then
                     USAFSI_arrays%icemask(c,r) = icepnt
                     USAFSI_arrays%icecon(c,r) = USAFSI_arrays%oldcon(c,r)
                  else
                     ! Prior analysis is either ice free or missing.  Assume
                     ! ice free.
                     USAFSI_arrays%icemask(c,r) = 0
                     USAFSI_arrays%icecon(c,r) = 0
                     USAFSI_arrays%iceage(c,r) = 0
                  end if
               end if
            end if

            ! Apply SST filter
            if (USAFSI_arrays%sst(c,r) > 275) then
               USAFSI_arrays%icemask(c,r) = 0
               USAFSI_arrays%icecon(c,r) = 0
               USAFSI_arrays%iceage(c,r) = 0
            end if
            
            ! Update age if this is a 12Z analysis
            if (USAFSI_arrays%icemask(c,r) == icepnt) then
               if ( runcycle == 12 .and. &
                    USAFSI_arrays%iceage(c,r) == &
                    USAFSI_arrays%iceage12z(c,r)) then
                  USAFSI_arrays%iceage(c,r) = &
                       min( (USAFSI_arrays%iceage(c,r)+1), maxage)
               end if
            end if
               
            ! EMK...Sanity check iceage.  Make sure not missing if
            ! icemask is defined.  This can happen if USAFSI is 
            ! coldstarted at a time other than 12Z.
            if (runcycle .ne. 12) then
               if (USAFSI_arrays%icemask(c,r) .ne. -1) then
                  if (USAFSI_arrays%iceage(c,r) .eq. -1) then
                     USAFSI_arrays%iceage(c,r) = 0
                  end if
               end if
            end if

         end do ! c
      end do ! r

   end subroutine run_seaice_analysis_ssmis

   ! Update sea ice based on remapped US Navy GOFS data
   subroutine run_seaice_analysis_gofs(month, runcycle, nc, nr, landmask)

      ! Imports
      use LDT_usafsiMod, only: usafsi_settings
      use USAFSI_arraysMod, only: USAFSI_arrays
      use USAFSI_paramsMod

      ! Defaults
      implicit none
      
      ! Arguments
      integer, intent(in) :: month
      integer, intent(in) :: runcycle
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: landmask(nc,nr)

      ! Local variables
      integer :: hemi
      logical :: check_sst
      integer :: c, r

      do r = 1, nr
         do c = 1, nc
            
            ! Only run over water
            if (landmask(c,r) > 0.5) then
               USAFSI_arrays%icecon(c,r) = -1
               USAFSI_arrays%icemask(c,r) = -1
               USAFSI_arrays%iceage(c,r) = -1
               cycle
            end if
      
            ! Find hemisphere
            if (USAFSI_arrays%ptlat(c,r) >= 0) then
               hemi = 1
            else
               hemi = 2
            end if

            check_sst = .true.

            ! Use the GOFS data if available.  Otherwise, try to fall back
            ! on prior analysis subject to certain constraints.
            if (USAFSI_arrays%gofs_icecon(c,r) >= 0) then
               ! We have valid GOFS data
               USAFSI_arrays%icecon(c,r) = &
                    nint(USAFSI_arrays%gofs_icecon(c,r))
               if (USAFSI_arrays%icecon(c,r) > usafsi_settings%minice) then
                  USAFSI_arrays%icemask(c,r) = icepnt
               else
                  USAFSI_arrays%icemask(c,r) = 0
                  USAFSI_arrays%iceage(c,r) = 0
               end if
               check_sst = .false.  ! Defer to GOFS

            else if (abs(USAFSI_arrays%ptlat(c,r)) < &
                 usafsi_settings%latchk(month,hemi)) then
               ! GOFS is missing at this point, and we are in low latitudes.
               ! Force no seaice.
               USAFSI_arrays%icemask(c,r) = 0
               USAFSI_arrays%icecon(c,r) = 0
               USAFSI_arrays%iceage(c,r) = 0
               
            else if (abs(USAFSI_arrays%ptlat(c,r)) >= &
                 usafsi_settings%icelat(month,hemi)) then
               ! GOFS is missing at this point, and we are in high latitudes.
               ! Force sea ice and find a reasonable value.
               USAFSI_arrays%icemask(c,r) = icepnt
               if (USAFSI_arrays%oldcon(c,r) > usafsi_settings%minice) then
                  USAFSI_arrays%icecon(c,r) = &
                       USAFSI_arrays%oldcon(c,r) ! Prior analysis
               else
                  USAFSI_arrays%icecon(c,r) = icedef ! Bogus value
               end if

            else
               ! GOFS is missing at this point, and we are in mid-latitudes.
               ! Attempt to use prior analysis.
               if (USAFSI_arrays%oldcon(c,r) > usafsi_settings%minice) then
                  USAFSI_arrays%icemask(c,r) = icepnt
                  USAFSI_arrays%icecon(c,r) = USAFSI_arrays%oldcon(c,r)
               else
                  ! Prior analysis is either ice free or missing.  Assume
                  ! ice free.
                  USAFSI_arrays%icemask(c,r) = 0
                  USAFSI_arrays%icecon(c,r) = 0
                  USAFSI_arrays%iceage(c,r) = 0
               end if
            end if
               
            ! Apply the SST filter unless we used the GOFS sea ice value.
            if (check_sst) then
               if (USAFSI_arrays%sst(c,r) > 275) then
                  USAFSI_arrays%icemask(c,r) = 0
                  USAFSI_arrays%icecon(c,r) = 0
                  USAFSI_arrays%iceage(c,r) = 0
               end if
            end if
            
            ! Update age if this is a 12Z analysis
            if (USAFSI_arrays%icemask(c,r) == icepnt) then
               if ( runcycle == 12 .and. &
                    USAFSI_arrays%iceage(c,r) == &
                    USAFSI_arrays%iceage12z(c,r)) then
                  USAFSI_arrays%iceage(c,r) = &
                       min( (USAFSI_arrays%iceage(c,r)+1), maxage)
               end if
            end if

            ! EMK...Sanity check iceage.  Make sure not missing if
            ! icemask is defined.  This can happen if USAFSI is 
            ! coldstarted at a time other than 12Z.
            if (runcycle .ne. 12) then
               if (USAFSI_arrays%icemask(c,r) .ne. -1) then
                  if (USAFSI_arrays%iceage(c,r) .eq. -1) then
                     USAFSI_arrays%iceage(c,r) = 0
                  end if
               end if
            end if

         end do ! c
      end do ! r
         
   end subroutine run_seaice_analysis_gofs

   ! Private subroutine
   subroutine summer (obelev, hemi, oblat, month, towarm)

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME: SUMMER
      !**
      !**  PURPOSE: PLACE ZERO SNOW DEPTH WHERE SNOW IS EXTREMELY UNLIKELY
      !**
      !**  CALLED FROM: GETOBS
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  HEMI, MONTH, OBELEV, OBLAT
      !**
      !**   OUTPUT:  TOWARM
      !**
      !**  SEE USAFSI_paramsMod.F90 FOR PARAMETER DESCRIPTIONS
      !**
      !**  UPDATES
      !**  =======
      !**  26 APR 95  INITIAL UNISYS VERSION...........................SSGT CONRY/SYSM
      !**  23 OCT 97  ADDED INCLUDE OF TUNES PROC......................SSGT CONRY/DNXM
      !**  08 FEB 99  UPDATED PROLOG TO REFLECT CALL FROM GETSMI......SRA HERKAMP/DNXM
      !**  23 FEB 01  PORTED FROM UNISYS MAINFRAME TO UNIX............SSGT MILLER/DNXM
      !**  21 JUL O4  CONVERTED TO FORTRAN 90, MOVED PARMATER
      !**             AND TUNBLK DEFINITIONS TO MODULES...............MR EYLANDER/DNXM
      !**  15 SEP 10  ADDED ABS() TO LATITUDE CHECK...............MR LEWISTON/16WS/WXE
      !**  22 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**  13 Dec 19  Renamed USAFSI...Eric Kemp, NASA GSFC/SSAI
      !**
      !*******************************************************************************
      !*******************************************************************************

      ! Imports
      use LDT_usafsiMod, only: usafsi_settings
      use USAFSI_paramsMod

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in)   :: obelev             ! OBSERVATION ELEAVATION (METERS)
      integer, intent(in)   :: hemi               ! HEMISPHERE (1 = NH, 2 = SH)
      integer, intent(in)   :: oblat              ! OBSERVATION LATITUDE (DEGREES * 100)
      integer, intent(in)   :: month              ! MONTH OF YEAR (1-12)
      logical, intent(out)  :: towarm             ! FLAG INDICATING SNOW IS UNLIKELY

      ! Local variables
      integer               :: sumnth (12, 2, 2)  ! SUMMER MONTHS (MONTHS, HEMI, LATS)


      ! SUMNTH IS DIMENSIONED AS FOLLOWS:
      !
      !               J  F  M  A  M  J  J  A  S  O  N  D
      !               A  E  A  P  A  U  U  U  E  C  O  E
      !               N  B  R  R  Y  N  L  G  P  T  V  C
      !
      !           NH  x  x  x  x  x  x  x  x  x  x  x  x  -> LAT'S 20 - 30
      !           SH  x  x  x  x  x  x  x  x  x  x  x  x      "
      !           NH  x  x  x  x  x  x  x  x  x  x  x  x  -> LAT'S  0 - 20
      !           SH  x  x  x  x  x  x  x  x  x  x  x  x      "

      data sumnth / 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,                     &
           1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                     &
           0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0,                     &
           1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1 /

      towarm = .false.

      ! CHECK IF LATITUDE IS LESS THAN OR EQUAL TO 30 DEGREES AND
      ! LESS THAN OR EQUAL TO 20 DEGREES.
      if (abs (oblat) <= usafsi_settings%trplat(2)) then
         if (abs (oblat) <= usafsi_settings%trplat(3)) then

            if (obelev <= usafsi_settings%elvlim(1)) then
               ! IF THE ELEVATION IS LESS THAN OR EQUAL TO 1000 METERS,
               ! SET TOWARM FLAG TO TRUE.
               towarm = .true.
            elseif ((sumnth(month,hemi,2) == 1) .and.       &
                 (obelev <= usafsi_settings%elvlim(2)))      then
               ! IF IT'S A SUMMER MONTH AND THE ELEVATION IS LESS THAN OR
               ! EQUAL TO 1500 METERS, SET TOWARM FLAG TO TRUE.
               towarm = .true.
            endif

         elseif ((sumnth(month,hemi,2) == 1) .and.          &
              (obelev <= usafsi_settings%elvlim(3)))      then
            ! THE LATITUDE IS BETWEEN 20 AND 30 DEGRESS.  IF IT'S A SUMMER
            ! MONTH AND THE ELEVATION IS LESS THAN OR EQUAL TO 1000 METERS,
            ! SET TOWARM FLAG TO TRUE.
            towarm = .true.
         endif

      elseif ((sumnth(month,hemi,1) == 1) .and.              &
           (obelev <= usafsi_settings%elvlim(4)))      then
         ! THE LATITUDE IS BETWEEN 30 AND 40 DEGRESS.  IF IT'S A SUMMER
         ! MONTH AND THE ELEVATION IS LESS THAN OR EQUAL TO 1000 METERS,
         ! SET TOWARM FLAG TO TRUE.
         towarm = .true.
      endif

      return

   end subroutine summer

  ! Yeosang Yoon: new 10-km snow climatology
   subroutine getclimo (month, static)

      ! Imports
      use LDT_logMod, only: LDT_verify, ldt_logunit
      use USAFSI_arraysMod, only: USAFSI_arrays
      use netcdf

      ! Defaults
      implicit none

      ! Arguments
      integer,       intent(in)   :: month            ! MONTH OF YEAR (1-12)
      character*255, intent(in)   :: static           ! STATIC FILE DIRECTORY PATH

      ! Local variables
      character*4                 :: cmonth  (12)     ! MONTH OF YEAR
      character*255               :: file_path        ! FULLY-QUALIFIED FILE NAME

      data cmonth        / '_jan', '_feb', '_mar', '_apr', '_may', '_jun', &
          '_jul', '_aug', '_sep', '_oct', '_nov', '_dec' /

      integer          :: ncid, varid

      ! RETRIEVE THE CLIMATOLOGY FOR THE MONTH.
      ! THE CLIMO FILE CONTAINS AN ARRAY FOR EACH OF THE 12 MONTHS.
      ! EACH MONTH IS STORED CONSECUTIVELY STARTING WITH JANUARY.
      file_path = trim(static) //'/snoclimo_10km/'// 'snoclimo_0p10deg' &
           //cmonth(month) // '.nc'

      write(ldt_logunit,*)'[INFO] Reading ', trim(file_path)
      call LDT_verify(nf90_open(path=file_path, mode=nf90_nowrite, ncid=ncid), &
            '[ERR] Error in nf90_open for '//trim(file_path))
      call LDT_verify(nf90_inq_varid(ncid=ncid, name="snoclimo", varid=varid), &
            '[ERR] Error in nf90_inq_varid for snow climatology')
      call LDT_verify(nf90_get_var(ncid=ncid, varid=varid, values=USAFSI_arrays%climo), &
            '[ERR] Error in nf90_get_var for snow climatology')
      call LDT_verify(nf90_close(ncid), &
            '[ERR] Error in nf90_close for '//trim(file_path))

   end subroutine getclimo

   ! New routine to read FNMOC SST field from GRIB1 file, using ECCODES
   subroutine read_grib1_sst(file_grib, sst_igrid, sst_jgrid, &
        sst_0p25deg, grstat)

     ! Imports
#if (defined USE_GRIBAPI)
     use grib_api
#endif
     use LDT_logMod, only: LDT_logunit

     ! Defaults
     implicit none

     ! Arguments
     character(len=*), intent(in) :: file_grib
     integer, intent(in) :: sst_igrid
     integer, intent(in) :: sst_jgrid
     real, intent(inout) :: sst_0p25deg(sst_igrid, sst_jgrid)
     integer, intent(out) :: grstat

     ! Locals
     integer :: ftn
     integer :: igrib
     integer :: ierr
     integer :: nvars
     integer :: iedition
     integer :: igriddef
     integer :: icenter
     integer :: iparameter
     integer :: ileveltype
     integer :: ilevel
     character(len=100) :: gtype
     integer :: Ni, Nj
     real, allocatable :: dum1d(:)
     integer :: i, j, k

     grstat = 1

#if (defined USE_GRIBAPI)
     call grib_open_file(ftn, trim(file_grib), 'r', ierr)
     if (ierr .ne. 0) then
        write(ldt_logunit,*) '[WARN] Failed to open - ', trim(file_grib)
        return
     end if

     write(ldt_logunit,*)'[INFO] Reading ', trim(file_grib)

     call grib_count_in_file(ftn, nvars, ierr)
     if (ierr .ne. 0) then
        write(ldt_logunit,*) &
             '[WARN] error in grib_count_in_file for ', trim(file_grib)
        call grib_close_file(ftn)
        return
     end if

     ! Loop through the fields until we find SST
     do k = 1, nvars
        call grib_new_from_file(ftn, igrib, ierr)
        if (ierr .ne. 0) then
           write(ldt_logunit,*) '[WARN] failed to read ' // trim(file_grib)
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        endif

        call grib_get(igrib, 'editionNumber', iedition, ierr)
        if ( ierr .ne. 0 ) then
           write(ldt_logunit,*) &
                '[WARN] error in grib_get: editionNumber in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        endif
        if (iedition .ne. 1) then
           write(ldt_logunit,*) &
                '[WARN] No GRIB1 record found in read_grib1_sst!'
           call grib_release(igrib, ierr)
           cycle
        endif

        call grib_get(igrib, 'centre', icenter, ierr)
        if ( ierr .ne. 0 ) then
           write(ldt_logunit,*) &
                '[WARN] error in grib_get: centre in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        endif
        if (icenter .ne. 58) then
           write(ldt_logunit,*)'[WARN] No FNMOC message in read_grib1_sst!'
           call grib_release(igrib, ierr)
           cycle
        endif

        call grib_get(igrib, 'gridDefinition', igriddef, ierr)
        if (ierr .ne. 0) then
           write(ldt_logunit,*) &
                '[WARN] error in grib_get: gridDefinition in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        endif
        if (igriddef .ne. 200) then
           write(ldt_logunit,*)'[WARN] Wrong SST grid in read_grib1_sst!'
           call grib_release(igrib, ierr)
           cycle
        endif

        call grib_get(igrib, 'gridType', gtype, ierr)
        if (ierr .ne. 0) then
           write(ldt_logunit,*) '[WARN] error in grid_get: gridtype in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        endif
        if (gtype .ne. "regular_ll") then
           write(ldt_logunit,*) &
                '[WARN] GRIB data not on regular lat-lon grid!'
           call grib_release(igrib, ierr)
           cycle
        endif

        call grib_get(igrib, 'indicatorOfParameter', iparameter, ierr)
        if (ierr .ne. 0) then
           write(ldt_logunit,*) &
                '[WARN] error in grib_get: indicatorOfParameter in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        endif
        if (iparameter .ne. 80) then
           write(ldt_logunit,*)'[WARN] Wrong GRIB parameter in read_grib1_sst!'
           call grib_release(igrib, ierr)
           cycle
        endif

        call grib_get(igrib, 'indicatorOfTypeOfLevel', ileveltype, ierr)
        if (ierr .ne. 0) then
           write(ldt_logunit,*) &
                '[WARN] error in grib_get: indicatorOfTypeOfLevel in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        endif
        if (ileveltype .ne. 160) then
           write(ldt_logunit,*) &
                '[WARN] Wrong GRIB level type in read_grib1_sst!'
           call grib_release(igrib, ierr)
           cycle
        endif

        call grib_get(igrib, 'level', ilevel, ierr)
        if (ierr .ne. 0) then
           write(ldt_logunit,*) &
                '[WARN] error in grib_get: level in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        endif
        if (ilevel .ne. 0) then
           write(ldt_logunit,*)'[WARN] Wrong GRIB level in read_grib1_sst!'
           call grib_release(igrib, ierr)
           cycle
        endif

        call grib_get(igrib, 'Ni', Ni, ierr)
        if (ierr .ne. 0) then
           write(ldt_logunit,*) '[WARN] error in grid_get:Ni in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        endif
        if (Ni .ne. sst_igrid) then
           write(ldt_logunit,*) '[WARN] Wrong GRIB Ni dimension in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           cycle
        endif

        call grib_get(igrib, 'Nj', Nj, ierr)
        if (ierr .ne. 0) then
           write(ldt_logunit,*) '[WARN] error in grid_get:Nj in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        endif
        if (Nj .ne. sst_jgrid) then
           write(ldt_logunit,*) '[WARN] Wrong GRIB Nj dimension in ' // &
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           cycle
        endif

        ! We found the SST
        allocate(dum1d(Ni*Nj))
        call grib_get(igrib, 'values', dum1d)
        if (ierr .ne. 0) then
           write(ldt_logunit,*) &
                '[WARN] error in grib_get: values in ' //&
                'read_grib1_sst'
           call grib_release(igrib, ierr)
           call grib_close_file(ftn)
           return
        end if

        ! At this stage, we have the values of the field.
        call grib_release(igrib, ierr)
        call grib_close_file(ftn)
        do j = 1, Nj
           do i = 1, Ni
              sst_0p25deg(i,j) = dum1d(i + (j-1)*Ni)
           end do
        end do
        grstat = 0
        deallocate(dum1d)
        exit ! Get out of loop

     end do ! k

     if (grstat .ne. 0) then
        write(ldt_logunit,*) &
             '[WARN] No SST read in by ' //&
             'read_grib1_sst'
        call grib_close_file(ftn)
     end if
#endif

   end subroutine read_grib1_sst
end module USAFSI_analysisMod
