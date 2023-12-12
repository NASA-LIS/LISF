!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: AGRMET_geoest
! \label{AGRMET_geoest}
!
!
! !REVISION HISTORY:
!
!    31 mar 99  initial version ........................ mr moore/dnxm
!     2 sep 99  added option for geoswch to be 0, 1, or 2.  this allows
!               the rank values computed by geo_precip to either be 
!               used (geoswch=2) or not used (geoswch=1) in the final
!               calculation of the precip estimate.....capt hidalgo/dnxm
!     8 oct 99  ported to ibm sp-2, updated prolog, incorporated
!               FORTRAN 90 features..................capt hidalgo/agrmet
!    01 jun 00  changed grnk to integer, changed geornk to integer,
!               initialized gdgeornk to false.  Included error checking
!               for grnk values when used.  Added alert message if no
!               grnk file is found when grnks are being used.  Removed
!               variable k from passed variables and declarations, could
!               not find anywhere that variable was used................
!               ......................................ssgt campbell/dnxm
!    21 feb 01  reformatted the diagnostic prints..........mr gayno/dnxm
!    10 jun 02  changed all references to rtneph to cdfs2..mr gayno/dnxm
!    3 nov 2005 Sujay Kumar, Initial Code
!   19 Dec 2007 Marv Freimund, Simplify filename creation
!   11 Mar 2010 Changed program names in messages to LIS.  Left argument
!               to LIS_alert with old program name to keep alert numbers
!               unique........................Chris Franks/16WS/WXE/SEMS
!   10 MAR 2010 Added handling of multiple resolution grids - see in 
!               code comments for where.................Michael Shaw/WXE
!   10 MAR 2010 Add filter of corrupt GEOPRECIP files...Michael Shaw/WXE
!    9 May 2013 Bug fix for correctly reading georank files
!               ................................Ryan Ruhge/16WS/WXE/SEMS
!               
! !INTERFACE:    
subroutine AGRMET_geoest( n, j3hr, land ,gest, grnk, quad9r, &
                     imax, jmax, geoswch,&
                     alert_number, mesh) !Michael Shaw -mesh handling
                                         !of 64th polar - no mask, so
                                         !use 16th mask for now
! !USES: 
  use agrmet_forcingMod, only : agrmet_struc
  use LIS_coreMod,   only     : LIS_rc, LIS_masterproc
  use LIS_fileIOMod, only     : LIS_putget
  use LIS_timeMgrMod, only    : LIS_julhr_date
  use LIS_logMod, only        : LIS_alert, LIS_logunit

  implicit none
! !ARGUMENTS: 
  integer, intent(in)               :: n 
  integer                           :: imax
  integer                           :: jmax
  integer, intent(in)               :: j3hr  
  integer, intent(inout)            :: alert_number
  integer, intent(in)               :: geoswch
  integer                           :: mesh
! Michael Shaw - mesh handling for 64th polar; no masks, so use 16th for now
  integer, intent(in)               :: land(int(imax/mesh),int(jmax/mesh),2)
  integer, intent(inout)            :: grnk(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real, intent(out)                 :: gest(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real, intent(in)                  :: quad9r
  logical :: is_geo_corrupted

!
! !DESCRIPTION:
!
!    to retrieve a precipitation estimate based on the geo\_precip
!    technique
!
!    \textbf{Method} \newline
!    
!    - read geo\_precip file for this time/hemisphere \newline
!      - if file does not exist, return to calling routine.  Geo\_precip
!        data will not be used.  Send alert message. \newline
!      - YDT: if data are corruputed, treat as missing file (1/8/08)
!    
!    - if geoswch equals 2, read in ssmi values, if file exists. \newline
!    - loop over hemisphere \newline
!      - replace geo\_precip "missing" values w/9999.0 \newline
!      - if value is valid and over land, place precip amount in output
!        array gest. \newline
!      - if geoswch equals 2 and rank value is valid, place in output
!        array grnk. \newline
!        - if rank value invalid, set to default value of 4  \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[n]     index of the nest
!   \item[alert\_number] alert number
!   \item[exists]    a logical that indicates whether or not a file exists
!   \item[gdgeornk]  logical indicating if the rank values are usable
!   \item[geoprc]    temp array used to read in geo precip values
!   \item[goornk]    temp array used to read in rank values
!   \item[geoswch]   geo-precip estimate switch \newline
!                              0 = don't, \newline
!                              1 = process w/o grnk \newline
!                              2 = process with grnk \newline
!   \item[gest]      final geo\_precip-based precip estimates (mm/3hr)   
!   \item[grnk]      final geo\_precip rank values \newline
!                              1 = better than ssmi est. \newline
!                              2 = better than present wx based est. \newline
!                              3 = better than cdfs2 est. \newline
!                              4 = default value of geo\_precip est. \newline
!                              5 = worse than climatological est. \newline
!   \item[hemi]      hemisphere (1=nh, 2=sh)
!   \item[i]         loop counter  
!   \item[ifil]      dummy character array for input file name
!   \item[imax]      number of gridpoints in east/west direction 
!   \item[j]         loop counter  
!   \item[jmax]      number of gridpoints in north/south direction
!   \item[land]      land/sea mask (0 = water, 1 = land)
!   \item[message]   message array for alert
!   \item[quad9r]    value of 9999.0 for initializing arrays
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[julhr\_date](\ref{LIS_julhr_date}) \newline
!    converts julian hour to a date format
!  \item[agrmet\_geoprec\_filename](\ref{agrmet_geoprec_filename}) \newline
!  \item[agrmet\_geornk\_filename](\ref{agrmet_georank_filename}) \newline
!  \item[LIS\_alert](\ref{LIS_alert}) \newline
!    prints a alert message
!  \item[LIS\_putget](\ref{LIS_putget}) \newline
!    retrieves data from the specified file
!  \end{description}
!EOP  
  integer                           :: i   
  integer                           :: j,jj
! Michael Shaw - first and last hemi for mask file, e.g. 
  integer                           :: hemi,first,last
  real                              :: geoprc(imax,jmax),geopr(imax,jmax)
  integer                           :: geornk(imax,jmax)
! Michael Shaw - temp grnk and gest for flipping arrays
  real                              :: grnk_tmp(2,imax,jmax),grnk_temp(1,imax,jmax)
  real                              :: gest_tmp(2,imax,jmax),gest_temp(1,imax,jmax)
  character*100                     :: ifil
  character*100                     :: message(20)
  character*30                      :: routine_name
  logical                           :: exists
  logical                           :: gdgeornk
  integer                           :: ip
  real                              :: udef
  real                              :: grnk_r(LIS_rc%lnc(n), LIS_rc%lnr(n))
  integer                           :: yr,mo,da,hr
  real                              :: LIS_rc_gridDesc(LIS_rc%nnest,50)

  data routine_name     / 'AGRMET_geoest' /

!     ------------------------------------------------------------------
!     executable code begins here ...
!     retrieve precip amounts
!     ------------------------------------------------------------------
  gest_tmp = quad9r
  gest_temp= quad9r
  grnk_temp= quad9r
  grnk_tmp= quad9r

! Michael Shaw - if agrmet_struc(n)%global_or_hemi = 0, it's the usual sh nh reads of 
! seperate files, otherwise it's a single file lat/lon file 

  if(agrmet_struc(n)%global_or_hemi == 0)then
     first = agrmet_struc(n)%shemi  
     last  = agrmet_Struc(n)%nhemi
  else
     first = 1
     last  = 1
  endif

 
  do hemi = first,last 

     gdgeornk = .false.

     call LIS_julhr_date( j3hr, yr,mo,da,hr)
     call agrmet_geoprec_filename(ifil,agrmet_struc(n)%agrmetdir,&
          agrmet_struc(n)%geodir,agrmet_struc(n)%use_timestamp,&
          hemi,yr,mo,da,hr,n,imax,jmax)
  
     inquire( file = trim(ifil), exist = exists)
     if ( .not. exists ) then 
        write(LIS_logunit,*) '[WARN] AGRMET_geoest/precip:  error opening file ',trim(ifil)
        write(LIS_logunit,*) '[WARN]  file does not exist'
        write(LIS_logunit,*) '[WARN]  geo precip estimate will not be performed'
        message = ' '
        message(1) = 'program:  LIS'
        message(2) = '  routine:  AGRMET_geoest'
        message(3) = '  error opening file '//trim(ifil)
        message(4) = '  file does not exist'
        message(5) = '  geo precip estimate will not be performed'
        alert_number = alert_number + 1
        if(LIS_masterproc) then 
           call LIS_alert( 'precip              ', &
                alert_number, message )
        endif
        return
     endif
     write(LIS_logunit,*) '[INFO] READING ', trim(ifil)
  
     call LIS_putget( geoprc, 'r', ifil, routine_name, &
          imax, jmax )

! Michael Shaw - Need to shift the 10th degree geoprecip by 180 degrees; the dataset starts at .05, not -179.95

     if(imax==3600)then
        geopr=cshift(geoprc,1800,dim=1)
        geoprc=geopr
     endif

!     ------------------------------------------------------------------ 
!     retrieve rank values if they are going to be used in the final
!     calculation of the estimate done by the calest subroutine
!     (indicated by geoswch = 2)
!     ------------------------------------------------------------------
     if (agrmet_struc(n)%geoswch .eq. 2) then
        call agrmet_georank_filename(ifil,agrmet_struc(n)%agrmetdir,&
             agrmet_struc(n)%geodir,agrmet_struc(n)%use_timestamp,&
             hemi,yr,mo,da,hr,n,imax,jmax)
        inquire( file = trim(ifil), exist = exists)

!     ------------------------------------------------------------------ 
!     if the rank files DO exist, read in the values and set the
!     boolean gdgeornk variable to true...this will be used later.
!     ------------------------------------------------------------------ 

        if ( exists ) then
           gdgeornk = .true.
           write(LIS_logunit,*) '[INFO] READING RANK FILE', trim(ifil)
        
           call LIS_putget( geornk, 'r', ifil, routine_name, &
                imax, jmax )
        else
           write(LIS_logunit,*)'[WARN] AGRMET_GEOEST/PRECIP: ERROR OPENING FILE ',trim(ifil)
           write(LIS_logunit,*)'[WARN]  FILE DOES NOT EXIST'
           write(LIS_logunit,*)'[WARN]  DEFAULT VALUES OF 4 WILL BE USED FOR GEO RANKS'
           message = ' '
           message(1)='program:  LIS'
           message(2)='  routine:  AGRMET_geoest'
           message(3)='  error opening file '//trim(ifil)
           message(4)='  file does not exist'
           message(5)='  default values of 4 will be used for geo ranks'
           alert_number=alert_number+1
           if(LIS_masterproc) then 
              call LIS_alert('precip              ',&
                   alert_number,message)          
           endif
        endif

!     ------------------------------------------------------------------ 
!     if the geo rank files do not exist continue processing but do not
!     read in the rank values.
!     ------------------------------------------------------------------

     endif

!YDT: 1/31/08, temp disable it 
! Need to enable this in case there are any more weird geoprecips or doing an "IOC rerun"

! Michael Shaw - was decided by AFWA and Goddard teams that it could be useful to
! flag for any "Anomalous" geoprecip files (any more obvious swaths and gross
! and sudden increases in extent and intensity of precip fields.
! Could be useful to check on PDF tweaks and whether they influence "detection"
! of these cases.

     if ( is_geo_corrupted(geoprc, imax, jmax, mo, hemi) ) then 
        write(LIS_logunit,*) '[WARN] AGRMET_geoest/precip:  data corrupted - ',trim(ifil)
        write(LIS_logunit,*) '[WARN] geo precip estimate will not be performed'
        message = ' '
        message(1) = 'program:  LIS'
        message(2) = '  routine:  AGRMET_geoest'
        message(3) = '  data corrupted in file '//trim(ifil)
        message(4) = '  ' 
        message(5) = '  geo precip estimate will be assigned missing'
        alert_number = alert_number + 1
        if(LIS_masterproc) then 
           call lis_alert( 'precip              ', &
                alert_number, message )
        endif
        
        geoprc = -1.0
        
     endif

!     ------------------------------------------------------------------
!     loop thru points in the hemisphere.  in the geo-precip module,
!     bad data points are indicated by values of -1.0.  change these
!     values to quad9r to be consistent with the rest of the code
!     ------------------------------------------------------------------
  
     do j = 1, jmax
        do i = 1, imax
           if( geoprc(i,j) .lt. 0.0) geoprc(i,j) = quad9r

!     ------------------------------------------------------------------
!     for land points ...
!     ------------------------------------------------------------------

           if( land(int(i/mesh),int(j/mesh),hemi) .gt. 0 )then

!     ------------------------------------------------------------------
!           check if input geo_precip amount is valid
!     ------------------------------------------------------------------

              if( geoprc(i,j) .lt. 9990.0 )then   

!     ------------------------------------------------------------------
!             if the input amount is valid, set the output value
!     ------------------------------------------------------------------

                 gest_tmp(hemi, i,j) = geoprc(i,j)

!     ------------------------------------------------------------------
!             if rank values are needed and have been subsequently
!             read in (as indicated by the boolean gdgeornk variable), 
!             assign them to the grnk array.
!     ------------------------------------------------------------------

                 if ( (geoswch .eq. 2) .and. (gdgeornk) ) then
                    if ((geornk(i,j).lt.1).or.(geornk(i,j).gt.5)) then
                       grnk_tmp(hemi,i,j) = 4
                       write(LIS_logunit,*)' '
                       write(LIS_logunit,*)'--------------------------------------------'
                       write(LIS_logunit,*)'[WARN] Bad geo precip rank value'
                       write(LIS_logunit,*)'[WARN] This value will be set to a default of 4'
                       write(LIS_logunit,*)'[WARN] geornk(',i,',',j,') = ',geornk(i,j)
                       write(LIS_logunit,*)'[WARN] from file :'
                       write(LIS_logunit,*) trim(ifil)
                       write(LIS_logunit,*)'--------------------------------------------'
                       write(LIS_logunit,*)' '
                    else
                       grnk_tmp(hemi, i,j) = geornk(i,j)
                    endif
                 else
                    grnk_tmp(hemi, i,j) = 1 ! if not using grank file, assuming that want to use geoprecip for sure 
                 endif
              endif
           endif
        enddo
     enddo
     
  enddo

  ip = 1
  udef = -9999.0

 ! Michael Shaw - IN CASE WE NEED TO FLIP THIS THING

  if(agrmet_struc(n)%global_or_hemi==1)then
    do i=1,imax
      do j=1,jmax
        jj=jmax-(j-1)
        gest_temp(1,i,j)=gest_tmp(1,i,jj)
      enddo
    enddo
    gest_tmp(1,:,:)=gest_temp(1,:,:)
  endif

! Michael Shaw - include dimensions in pass to interp module
! for flexibility 

  call interp_agrmetvar(n,ip,gest_tmp,udef,gest,imax,jmax)
  ip = 2
  udef = quad9r

! Michael Shaw - IN CASE WE NEED TO FLIP THIS THING
 
  if(agrmet_struc(n)%global_or_hemi==1)then
    do i=1,imax
      do j=1,jmax
        jj=jmax-(j-1)
        grnk_temp(1,i,j)=grnk_tmp(1,i,jj)
      enddo
    enddo
    grnk_tmp(1,:,:)=grnk_temp(1,:,:)
  endif

! Michael Shaw - include dimensions in pass to interp module
! for flexibility

  call interp_agrmetvar(n,ip,grnk_tmp,udef,grnk_r,imax,jmax)
!Turned off by SVK for enabling global forcings
!  call agrmet_fillgaps(n,ip,grnk_r)

  where (gest == udef)
     grnk_r = udef
  end where

  do j = 1, LIS_rc%lnr(n)
     do i = 1, LIS_rc%lnc(n)
        grnk(i,j) = nint(grnk_r(i,j))
     enddo
  enddo

  return

!     ------------------------------------------------------------------
!     alert handling
!     ------------------------------------------------------------------


end subroutine AGRMET_geoest

!BOP
! 
! !ROUTINE: agrmet_geoprec_filename
! \label{agrmet_geoprec_filename}
!
! !INTERFACE: 
subroutine agrmet_geoprec_filename(fname,rootdir,dir,&
     use_timestamp,hemi,yr,mo,da,hr,n,imax,jmax) ! Michael Shaw - need
                                                 ! dimensions (or something
                                                 ! ) to distinguish resolution
                                                 ! Could do this differently -
                                                 ! Legacy from previous version
                                                 ! to ensure proper values on multiple processors
  use AGRMET_forcingMod, only : agrmet_struc
  implicit none
! !USES: 
  character(*)        :: fname
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: hemi,n
  integer, intent(in) :: yr,mo,da,hr,imax,jmax
! 
! !DESCRIPTION: 
!  This routines generates the name of the GEO precip file
!  by appending the hemisphere and timestamps to the root directory. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[use\_timestamp]
!    flag to indicate whether the directories 
!    should be timestamped or not
!   \item[rootdir]
!    path to the root directory containing the data
!   \item[dir]
!    name of the subdirectory containing the data
!   \item[fname]
!    created filename
!   \item[yr]
!    4 digit year
!   \item[mo]
!    integer value of month (1-12)
!   \item[da]
!    day of the month
!   \item[hr]
!    hour of the day
!  \end{description}
!EOP

  character( 2), parameter :: fhemi(2) = (/'nh','sh'/)
  character(10)            :: ftime1, ftime2

  write (UNIT=ftime2, FMT='(i4, i2.2, i2.2, i2.2)') yr, mo, da, hr

  if(imax==512)then
  if (use_timestamp .eq. 1) then 
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/prec08_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  else
     fname = trim(rootdir) //   '/'  // trim(dir) // '/prec08_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  endif
  elseif(imax==1024)then
  if (use_timestamp .eq. 1) then 
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/prec16_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  else
     fname = trim(rootdir) //   '/'  // trim(dir) // '/prec16_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  endif
  elseif(imax==4096)then
  if (use_timestamp .eq. 1) then
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/prec64_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  else
     fname = trim(rootdir) //   '/'  // trim(dir) // '/prec64_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  endif
  elseif(imax==3600)then
  if (use_timestamp .eq. 1) then
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/prec' // &
             '0p1deg_ll' // '.03hr.' // ftime2
  else
     fname = trim(rootdir) // '/' // trim(dir) // '/prec' // &
             '0p1deg_ll' // '.03hr.' // ftime2
  endif
  endif
end subroutine agrmet_geoprec_filename

!BOP
! 
! !ROUTINE: agrmet_georank_filename
! \label{agrmet_georank_filename}
!
! !INTERFACE: 
subroutine agrmet_georank_filename(fname,rootdir,dir,use_timestamp,&
     hemi,yr,mo,da,hr,n,imax,jmax)               ! Michael Shaw - need
                                                 ! dimensions (or something
                                                 ! ) to distinguish resolution
                                                 ! Could do this differently
                                                 ! Legacy from previous version 
  use AGRMET_forcingMod, only : agrmet_struc
  implicit none
! !ARGUMENTS: 
  character(*)        :: fname
  character(*)        :: dir
  character(*)        :: rootdir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: hemi,n
  integer, intent(in) :: yr,mo,da,hr,imax,jmax
! 
! !DESCRIPTION: 
!  This routines generates the name of the GEO rank file
!  by appending the hemisphere and timestamps to the root directory. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[use\_timestamp]
!    flag to indicate whether the directories 
!    should be timestamped or not
!   \item[rootdir]
!    path to the root directory containing the data
!   \item[dir]
!    name of the subdirectory containing the data
!   \item[fname]
!    created filename
!   \item[yr]
!    4 digit year
!   \item[mo]
!    integer value of month (1-12)
!   \item[da]
!    day of the month
!   \item[hr]
!    hour of the day
!  \end{description}
!EOP

  character(2), parameter :: fhemi(2) = (/'nh','sh'/)
  character(10)           :: ftime1, ftime2

  write (UNIT=ftime2, FMT='(i4, i2.2, i2.2, i2.2)') yr, mo, da, hr

  if(imax==512)then
  if (use_timestamp .eq. 1) then 
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/rank08_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  else
     fname = trim(rootdir) //  '/'   // trim(dir) // '/rank08_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  endif
  elseif(imax==1024)then
  if (use_timestamp .eq. 1) then
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/rank16_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  else
     fname = trim(rootdir) //  '/'   // trim(dir) // '/rank16_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  endif
  elseif(imax==4096)then
  if (use_timestamp .eq. 1) then
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/rank64_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  else
     fname = trim(rootdir) //  '/'   // trim(dir) // '/rank64_' // &
             fhemi(hemi) // '.03hr.' // ftime2
  endif
  elseif(imax==3600)then
  if (use_timestamp .eq. 1) then
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/rank' // &
             '0p1deg_ll' // '.03hr.' // ftime2
  else
     fname =  trim(rootdir) //  '/'   // trim(dir) // '/rank' // &
             '0p1deg_ll' // '.03hr.' // ftime2
  endif 
  endif

end subroutine agrmet_georank_filename

!BOP
!
! !ROUTINE: is_geo_corrupted
! \label{is_geo_corrupted}
!
! !INTERFACE:
function is_geo_corrupted(geodata, imax, jmax, mo, hemi) 


  implicit none
! !ARGUMENTS:
  logical is_geo_corrupted
  integer, intent(in) :: imax, jmax 
  integer, intent(in) :: hemi
  integer, intent(in) :: mo
  real, intent(in) :: geodata(imax, jmax) 
!
! !DESCRIPTION:
! This function does quality control of
! geoprecip raw data, to see if there are corrupted patches
! as the one on 0Z19Jan2006.
! Inspect the histogram and see if there are some outliers there.
!  returns .TRUE. if data corruption is detected. Otherwise .FALSE.
!
!  The arguments are:
!  \begin{description}
!   \item[imax]      number of gridpoints in east/west direction 
!   \item[jmax]      number of gridpoints in north/south direction
!   \item[mo]
!    integer value of month (1-12)
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[geodata]
!    imax by jmax array holding the raw geoprecip data 
!  \end{description}
!EOP

  integer, parameter :: nbins=10
  real, parameter ::  rmax=100.0, rmin=0.5
  integer :: ip

  real       :: pdf0(0:nbins) ! empirical values of normal PDF
  real       :: pdf(0:nbins), tpdf, amp
  integer    :: j, k

! only the first 5 bins are enough 
  pdf0(0) = 74.0
  pdf0(1) = 22.0
  pdf0(2) = 5.0
  pdf0(3) = 2.1
  pdf0(4) = 1.4

!add seasonal variation
       Do ip = 0, 4
         amp = pdf0(ip) * 0.13
         if (hemi.eq.1) then ! NH
           pdf0(ip) = pdf0(ip) - amp * (1.0 + cos(2.0*3.1416*(mo-1)/12) )
         else ! SH
           pdf0(ip) = pdf0(ip) - amp * (1.0 - cos(2.0*3.1416*(mo-1)/12) )
         end if
       end do

        pdf = 0.0
        Do k=1, jmax
         Do j=1, imax 
           if(geodata(j, k) .GT. rmin .and. geodata(j, k) .LE. rmax ) then
            ip = nint( (geodata(j, k) - rmin) / (rmax - rmin) * nbins )
            pdf(ip) = pdf(ip) + 1
           end if
         End Do
        End Do

        is_geo_corrupted = .FALSE.
        Do ip = 0, 4
          tpdf =  pdf(ip) * 1000 /(imax * jmax)
          if ( tpdf .GT. pdf0(ip) ) is_geo_corrupted = .TRUE. 
        End Do
    
end function is_geo_corrupted
