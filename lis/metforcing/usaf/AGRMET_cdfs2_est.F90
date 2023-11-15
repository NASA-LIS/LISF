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
! !ROUTINE: AGRMET_cdfs2_est
! \label{AGRMET_cdfs2_est}
!
! !REVISION HISTORY:
!
!
!    07 mar 97  initial version......................ssgt mccormick/sysm
!    24 mar 97  added the 'source' array which gives and sets each  
!               estimate type with a specified integer source flag; 
!               determining the values with the highest estimate
!               information.............................ssgt miller/sysm
!    19 may 97  redimensioned source array to write each 3-hour sourced 
!               estimate into its own file..............ssgt miller/sysm
!    04 dec 97  changed name of routine from crest to rtnest to more
!               correctly reflect functionality. modified to generate   
!               rtneph based estimates for all points and not climo 
!               estimates. also removed determination of a source from  
!               this routine. updated prolog and brought up to  
!               standards......................capt andrus/dnxm(agromet)
!    31 mar 99  changed box looping to hemisphere looping and added new
!               file retrieval convention for rtneph data. also added
!               abort message to cover division by zero possibility ... 
!               ......................................... mr moore/dnxm
!     7 oct 99  ported to ibm sp-2, updated prolog, incorporated
!               FORTRAN 90 features.................capt hidalgo/agrmet
!    31 mar 00  changed variable "ifil" to character*120 to be
!               consistent with routine copen.  added error checks
!               after calls to copen and gribread..........mr gayno/dnxm
!    11 may 00  changed variable rtntim from integer to real.  this
!               fixed a problem degribbing the data.....................
!               .............................capt hidalgo, mr gayno/dnxm
!    21 feb 01  reduced third dimension (number of three hourly time
!               periods) of variable rest from 4 to 2 as module precip
!               runs in 6 hourly cycles now................mr gayno/dnxm
!    10 jun 02  modified to use cdfs2 data instead of rtneph............
!               ...........................................mr gayno/dnxm
!    3 nov 2005 Sujay Kumar, incorporated into LIS
!   19 Dec 2007 Marv Freimund, Simplify filename creation and correct bug
!   11 Mar 2010 Changed program names in messages to LIS.  Left argument
!               to LIS_alert with old program name to keep alert numbers
!               unique........................Chris Franks/16WS/WXE/SEMS
!   31 MAR 2010 Add handling of 16th polar "native" mesh................
!               see in code comments for where..........Michael Shaw/WXE
!
!
! !INTERFACE:    
subroutine AGRMET_cdfs2_est( n,k, cliprc, clippd,&
     clirtn, mnpr, mxpr, mnpd, mxpd, &
     cldth, mdmv, mdpe, ovpe, & 
     cdfs2int, cdfs2est, & 
     alert_number, j3hr )
! !USES: 
  use AGRMET_forcingMod, only : agrmet_struc
  use LIS_LMLCMod,       only : LIS_LMLC
  use LIS_coreMod,       only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod,    only : LIS_julhr_date
  use LIS_logMod,        only : LIS_alert, LIS_abort, LIS_logunit, LIS_endrun

  implicit none
! !ARGUMENTS: 
  integer,       intent(in)         :: n 
  integer,       intent(in)         :: j3hr  
  integer,       intent(inout)      :: alert_number
  integer,       intent(in)         :: cdfs2int  
  integer,       intent(in)         :: k
  real,          intent(out)        :: cdfs2est(LIS_rc%lnc(n),LIS_rc%lnr(n),4) 
  real,          intent(in)         :: cldth  
  real,          intent(in)         :: clippd(LIS_rc%lnc(n),LIS_rc%lnr(n))   
  real,          intent(in)         :: cliprc(LIS_rc%lnc(n),LIS_rc%lnr(n))   
  real,          intent(in)         :: clirtn(LIS_rc%lnc(n),LIS_rc%lnr(n))   
  real,          intent(in)         :: mdmv  
  real,          intent(in)         :: mdpe  
  real,          intent(in)         :: mnpd  
  real,          intent(in)         :: mnpr  
  real,          intent(in)         :: mxpd  
  real,          intent(in)         :: mxpr  
  real,          intent(in)         :: ovpe  

!
! !DESCRIPTION:
!
!    to create a precipitation estimate based on climatological precip  
!    amount, climatological precip per precip day amounts,  
!    climatological rtneph cloud amounts and actual cdfs2 total cloud  
!    amounts.   
!
!    \textbf{Method} \newline
!
!    - retrieve current cloud data and pixel times from file. \newline
!      - if current data is missing or on a bad read, search back 
!        a maximum of two hours for a complete set of data.  if 
!        this search fails, then don't create an estimate for this
!        time period. \newline
!    - for each land point in the hemisphere:  \newline
!      - if valid time for this point is within the window .... \newline
!        - if current cdfs2 total cloud data is valid and is at least  
!          equal to the minimum cdfs2 total cloud amount threshold. \newline
!          - if the climatological precip amount for this grid point
!            meets or exceeds the minimum required climo precip amount  
!            for a cdfs2 based precip estimate. \newline
!            - if current cdfs2 cloud amount meets or exceeds climo
!              rtneph cloud amount, but is less than 100 percent. \newline
!              - determine climo precip mult factor. \newline
!              - determine precip-per-precip day mult factor based upon 
!                cloud amount. \newline
!              - calculate estimate by multiplying precip-per-precip
!                day by the climo precip mult factor, the precip-per-   
!                precip day factor, and the monthly rthneph factor. \newline
!            - else estimate zero precip \newline
!          - else estimate zero precip  \newline
!        - else estimate zero precip \newline
!      - else do not estimate precip \newline
!
!  The arguments and variables are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[alert\_number] alert number  
!  \item[cdfs2est]  cdfs2 cloud amount based precip est (mm/3hr) 
!  \item[cdfs2int]  cdfs2 time interval to look for cloud amount
!  \item[cldamt]    current cdfs2 cloud amount 
!  \item[cldth]     cloud threshold to generate a cdfs2 based
!                  precipitation estimate
!  \item[cldtim]    cdfs2 cloud pixel time on the agrmet grid 
!                  (time of newest data at gridpoint)  
!  \item[clifac]    climatological multiplication factor
!  \item[clippd]    climatological precip-per-precip day amount (mm/3hr)  
!                   interpolated to the current day
!  \item[cliprc]    climatological monthly precip amount (mm/3hr) for 
!                   interpolated to the current day
!  \item[clirtn]    climatological rtneph percent cloud cover
!                   interpolated to the current day
!  \item[date10]    3-hrly 10-digit date time group (YYYYMMDDHH)
!  \item[error\_flag] logical flag for good/bad cdf2 data ingest
!  \item[hemi]      hemisphere (1=north, 2=south) 
!  \item[i]         i-coord loop counter  
!  \item[icdfs2]    i-dimension of cdfs2 grid
!  \item[ifil]      dummy character array for input file name
!  \item[imax]      number of gridpoints in east/west direction 
!  \item[istat]     i/o status
!  \item[j]         j-coord loop counter
!  \item[jcdfs2]    j-dimension of cdfs2 grid
!  \item[jmax]      number of gridpoints in north/south direction
!  \item[j3hr]      julian hour of end of precip accumulation period
!  \item[k]         3-hrly loop counter
!  \item[land]      points of major land masses
!  \item[maxth]     max actual cdfs2 cloud cover threshold for max
!              precip estimate 
!  \item[mdmv]      median cloud cover percentage to move to for the
!              cdfs2 based precipitation estimate
!  \item[mdpe]      percent to move to median cloud cover percentage
!              to move to for the cdfs2 based precipitation estimate
!  \item[message]   message array for abort
!  \item[minth]     min actual cloud cover threshold for min
!              precip estimate 
!  \item[mnpd]      minimum precip-per-precip day multiplier used to
!              generate a non-zero cdfs2 total cloud based precip
!              estimate
!  \item[mxpd]      maximum precip-per-precip day multiplier used to
!              generate a non-zero cdfs2 total cloud based precip
!              estimate
!  \item[mnpr]      minimum 3-hour climo precip value required to
!              generate a non-zero cdfs2 total cloud based precip
!              estimate
!  \item[mxpr]      max 3-hour climo precip value at which the maximum
!              precip-per-precip day multiplier is used to
!              generate a cdfs2 total cloud based precip estimate
!  \item[ovpe]      overcast percentage to move to for cdfs2 based
!              precipitation estimate
!  \item[ppddif]    difference between mxrppd and mnrppd
!  \item[ppdfac]    precip-per-precip day mult factor based upon cdfs2
!              cloud amount 
!  \item[time]      loop index for reading in cdfs2 data 
!  \item[times]     cdfs2 pixel times
!  \item[totcld]    array of cdfs2 total cloud amounts on the agrmet grid
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[julhr\_date](\ref{LIS_julhr_date}) \newline
!   converts julian hour to a date format
!  \item[agrmet\_cdfs\_totalcld\_filename](\ref{agrmet_cdfs_totalcld_filename}) \newline
!   generates the filename for CDFS2 total cloud data
!  \item[AGRMET\_julhr\_date10](\ref{AGRMET_julhr_date10}) \newline
!    converts julian hour to a 10 character date string
!  \item[lis\_alert](\ref{LIS_alert}) \newline
!    prints a alert message
!  \item[lis\_abort](\ref{LIS_abort}) \newline
!    aborts if a fatal error occurs
!  \item[agrmet\_cdfs\_pixltime\_filename](\ref{agrmet_cdfs_pixltime_filename}) \newline
!    generates the filename to read CDFS2 pixel times
!  \end{description}
!EOP
  character*10                      :: date10
  character*120                     :: ifil
  character*255                     :: message(20)
  integer,       allocatable        :: times  ( :, : , : )
  integer*1,     allocatable        :: totalc ( :, : , : )
  real                              :: cldtim(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real                              :: cldtim_tmp(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
  integer                           :: ip
  real                              :: udef
  integer                           :: i
  integer,       parameter          :: icdfs2 = 1024
  integer                           :: istat
  integer                           :: j
  integer,       parameter          :: jcdfs2 = 1024
  integer                           :: time
  logical                           :: error_flag  
  real                              :: cldamt  
  real                              :: clifac  
  real                              :: ppddif  
  real                              :: maxth  
  real                              :: minth  
  real                              :: ppdfac  
  integer                           :: hemi
  real                              :: totcld_tmp(2,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
  real                              :: totcld(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                           :: yr,mo,da,hr

!      ------------------------------------------------------------------
!      executable code begins here ...
!      retrieve satellite-derived cdfs2 cloud analysis data for 
!      this hemisphere. note:  total cloud amount is used !!!
!      read current data.  however, if the current data does not exist or 
!      there is a read error, try to use previous hours data.  go back
!      a maximum of two hours.  if all three time periods are bad,
!      then return to calling routine and do not produce an estimate.
!      ------------------------------------------------------------------

  cdfs2est(:,:,k) = LIS_rc%udef
  error_flag = .false.
  
  allocate ( totalc (2,icdfs2, jcdfs2) )
  allocate ( times  (2,icdfs2, jcdfs2) )
  
  do hemi=1,2
     READ_DATA : do time = j3hr, (j3hr-2), -1
        
        call LIS_julhr_date( time, yr,mo,da,hr)

!      ------------------------------------------------------------------
!        read in total cloud amounts.  the data are packed into
!        one byte.
!      ------------------------------------------------------------------
     
        call agrmet_cdfs_totalcld_filename(ifil,agrmet_struc(n)%agrmetdir,&
             agrmet_struc(n)%clouddir,agrmet_struc(n)%use_timestamp,hemi,&
             yr,mo,da,hr)

        write(LIS_logunit,*)'[INFO] OPENING ', trim(ifil) 
        open(9, file=trim(ifil), access='direct', &
             recl=icdfs2*jcdfs2*1, iostat=istat, status="old")

!      ------------------------------------------------------------------
!        bad open, cycle loop and try to use previous hour's data.  
!      ------------------------------------------------------------------
     
        if (istat /= 0) then
           write(LIS_logunit,*)'[ERR] ERROR OPENING ',trim(ifil),istat
           error_flag = .true.
           cycle READ_DATA
        end if

!      ------------------------------------------------------------------
!        bad read, cycle loop and try to use previous hour's data.  
!      ------------------------------------------------------------------

        write(LIS_logunit,*)'[INFO] READING ',trim(ifil)
        read(9, rec=1, iostat=istat) totalc(hemi,:,:)
        close(9)
 
        if (istat /= 0) then
           write(LIS_logunit,*)'[ERR] ERROR READING ',trim(ifil)
           error_flag = .true.
           cycle READ_DATA
        end if

!      ------------------------------------------------------------------
!        read in pixel times.  the data are packed into
!        four bytes.
!      ------------------------------------------------------------------
        call agrmet_cdfs_pixltime_filename(ifil,agrmet_struc(n)%agrmetdir,&
             agrmet_struc(n)%clouddir,agrmet_struc(n)%use_timestamp,hemi,&
             yr,mo,da,hr)
        write(LIS_logunit,*)'[INFO] OPENING ', trim(ifil)
     
        open(9, file=trim(ifil), access='direct', &
             recl=icdfs2*jcdfs2*4, iostat=istat, status="old")

!      ------------------------------------------------------------------
!        bad open, cycle loop try to use previous hour's data.  
!      ------------------------------------------------------------------

        if (istat /= 0) then
           write(LIS_logunit,*)'[ERR] Cannot open ',trim(ifil)
           error_flag = .true.
           cycle READ_DATA
        end if
     
        write(LIS_logunit,*)'[INFO] READING ', trim(ifil)
        read(9, rec=1, iostat=istat) times(hemi, :,:)
        close(9)

!      ------------------------------------------------------------------
!        bad read, cycle loop try to use previous hour's data.  
!        if read is good, all data for this hour has been read in
!        sucessfully, therefore, exit the loop.
!      ------------------------------------------------------------------
        
        if (istat /= 0) then
           write(LIS_logunit,*)'[ERR] ERROR READING ',trim(ifil)
           error_flag = .true.
           cycle READ_DATA
        else
           exit READ_DATA
        end if
     
     enddo READ_DATA

!      ------------------------------------------------------------------
!      is the last istat is bad, then we did not successfully read
!      in the complete data for one hour.  therefore, we can't do
!      an estimate.  send an alert message to maintainer and
!      continue program.
!      ------------------------------------------------------------------

     if (istat /= 0) then
        
        write(LIS_logunit,*)
        write(LIS_logunit,*) '[WARN] IN ROUTINE AGRMET_CDFS2_EST: ERRORS WITH CDFS2 DATA.'
        write(LIS_logunit,*) '[WARN] DATA IS CORRUPT OR DOES NOT EXIST.'
        write(LIS_logunit,*) '[WARN] CDFS2 PRECIP ESTIMATE WILL NOT BE PERFORMED.'
        write(LIS_logunit,*)
        
        call AGRMET_julhr_date10 ( j3hr, date10 )
        
        message = ' '
        message(1) = 'program:  LIS'
        message(2) = '  routine:  AGRMET_cdfs2_est'
        message(3) = '  errors with cdfs2 data.'
        message(4) = '  data is corrupt or does not exist'
        message(5) = '  for time ' // date10
        message(6) = '  cdfs2 precip estimate will not be performed.'
        alert_number = alert_number + 1
        if(LIS_masterproc) then 
           call lis_alert( 'precip              ', &
                alert_number, message )
        endif
        deallocate ( totalc )
        deallocate ( times )
        
        return
     
     else

!      ------------------------------------------------------------------
!        we were able to read in a complete set of data for one hour.
!        first, convert pixel times from jul minutes to jul hours.
!      ------------------------------------------------------------------

        times(hemi,:,:) = times(hemi,:,:) / 60

!      ------------------------------------------------------------------
!        Michael Shaw - cdfs2 is 16th mesh, so decimate it down to 8th mesh by 
!        sampling every other point if native grid is 8th mesh.  However,
!        keep all data for 16th native grid..
!      ------------------------------------------------------------------
!        i*2-1 = 1,3,5,...1023 (i = 1,2,3,...512)
!        i     = 1,2,3,...1024 (i = 1,2,3,...1024)
!      ------------------------------------------------------------------

        do j = 1, agrmet_struc(n)%jmax
           do i = 1, agrmet_struc(n)%imax
              if(agrmet_struc(n)%jmax.eq.1024)then
                cldtim_tmp(hemi,i,j) = real(times(hemi, i,j))
                totcld_tmp(hemi,i,j) = real(totalc(hemi, i,j))
              elseif(agrmet_struc(n)%jmax.eq.512)then
                cldtim_tmp(hemi,i,j) = real(times(hemi, (i*2-1),(j*2-1)))
                totcld_tmp(hemi,i,j) = real(totalc(hemi, (i*2-1),(j*2-1)))
              endif
           enddo
        enddo
        
!      ------------------------------------------------------------------
!        if error flag is true, then we had problems reading the
!        most current data and had to use old data.  alert the
!        maintainer of this unfortunate development.
!      ------------------------------------------------------------------
     
        if (error_flag) then
           
           call AGRMET_julhr_date10 ( j3hr, date10 )
           
           message = ' '
           message(1) = 'program:  LIS'
           message(2) = '  routine:  AGRMET_cdfs2_est'
           message(3) = '  errors opening/reading cdfs2 data'
           message(4) = '  for time ' // date10
           message(5) = '  using previous hour/s data.'
           alert_number = alert_number + 1
           call lis_alert( 'precip              ', &
                alert_number, message )
           
        end if
        
     end if

  enddo

  deallocate ( totalc )
  deallocate ( times )
     
  ip = 1
  udef = 9999.0
  agrmet_struc(n)%global_or_hemi = 0 ! polar stereo hemispheric array versus the global lat/lon arrays also possibly interpolated with interp_agrmetvar
  call interp_agrmetvar(n,ip,cldtim_tmp,udef,cldtim,agrmet_struc(n)%imax,agrmet_struc(n)%jmax) !,agrmet_struc(n)%shemi,agrmet_struc(n)%nhemi)
  call interp_agrmetvar(n,ip,totcld_tmp,udef,totcld,agrmet_struc(n)%imax,agrmet_struc(n)%jmax) !,agrmet_struc(n)%shemi,agrmet_struc(n)%nhemi)

!      ------------------------------------------------------------------
!      now calculate the precip estimate.  first, pre-calculate
!      common cdfs2 delimiters.
!      ------------------------------------------------------------------

  ppddif = mxpd - mnpd  

!      ------------------------------------------------------------------
!      loop through the points in the hemisphere
!      ------------------------------------------------------------------

  do j = 1, LIS_rc%lnr(n)
     do i = 1, LIS_rc%lnc(n)
        
!      ------------------------------------------------------------------
!          if point is a land mass point and all the climo data is
!          valid, proceed ...
!      ------------------------------------------------------------------
        
        if( (LIS_LMLC(n)%landmask(i,j) .gt. 0) .and. &
             (clippd(i,j) .gt. -9990.0) &
             .and. (clirtn(i,j) .gt. -9990.0) )then
           
!      ------------------------------------------------------------------
!            continue processing if the pixel time is within the
!            specified time interval.  note: cdfs2int is specified
!            in the precip control file.
!      ------------------------------------------------------------------

           if ( cldtim(i,j) > (j3hr - cdfs2int)  .and. &
                cldtim(i,j) < (j3hr + cdfs2int)  ) then

!      ------------------------------------------------------------------
!              data are within the time window of interest.  check
!              if cloud amount is large enough to make a non-zero 
!              estimate.
!      ------------------------------------------------------------------

              cldamt = totcld(i,j)
              if( (cldamt .ge. cldth) .and. (cldamt .le. 100.0) )then  

!      ------------------------------------------------------------------
!                cloud amount is sufficient to generate a precip
!                estimate.  proceed to calculate the value ... 
!      ------------------------------------------------------------------

                 if( (mxpr - mnpr) .ne. 0.0 )then
                    clifac = ( cliprc(i,j) - mnpr ) &
                         / ( mxpr - mnpr )
                 else
                    message = ' '
                    message(1) = 'program:  LIS'
                    message(2) = '  routine:  AGRMET_cdfs2_est'
                    message(3) = '  error - division by zero'
                    message(4) = ' because mxpr equals mnpr.'
                    message(5) = '  fix control.precip file and restart.'
                    call lis_abort( message )
                    call LIS_endrun
                 endif
                 clifac = min( max( clifac, 0.0 ), 1.0 )
                 minth  = clirtn(i,j) &
                      + mdpe * ( mdmv - clirtn(i,j) )
                 maxth  = clirtn(i,j) &
                      + ovpe * ( 100.0 - clirtn(i,j) )
                 if( (maxth - minth) .ne. 0.0 )then
                    ppdfac = ( cldamt - minth ) / ( maxth - minth )
                 else
                    message = ' '
                    message(1) = 'program:  LIS'
                    message(2) = '  routine:  AGRMET_cdfs2_est'
                    message(3) = '  error - division by zero'
                    message(4) = ' because maxth equals minth.'
                    message(5) = '  fix control.precip file and restart.'
                    call lis_abort( message )
                    call LIS_endrun
                 endif
                 ppdfac = min( max( ppdfac, 0.0 ), 1.0 ) 
                 ppdfac = mnpd + ppddif * ppdfac 
                 cdfs2est(i,j,k) = clifac * ppdfac * clippd(i,j) 
              elseif( cldamt .ge. 0.0 )then 
                 
!      ------------------------------------------------------------------
!                cloud amount is insufficient to generate a non-zero
!                precip estimate, or exceeds 100.  since the amount is
!                valid (ge 0), set estimate to zero. 
!      ------------------------------------------------------------------
                 
                 cdfs2est(i,j,k) = 0.0
              endif
           endif
        endif
     enddo
  enddo

  return 
end subroutine AGRMET_cdfs2_est
 
!BOP
! 
! !ROUTINE: agrmet_cdfs_totalcld_filename
! \label{agrmet_cdfs_totalcld_filename}
! 
! !INTERFACE: 
subroutine agrmet_cdfs_totalcld_filename (fname,rootdir,dir,&
     use_timestamp,hemi,yr,mo,da,hr)

  implicit none
! !ARGUMENTS: 
  character(*)        :: fname
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: hemi
  integer, intent(in) :: yr,mo,da,hr
! 
! !DESCRIPTION: 
!  This routines generates the name of the CDFS2 total cloud data
!  by appending the hemisphere and timestamps to the root directory. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[fname]
!    created filename
!   \item[rootdir]
!    path to the root directory containing the data
!   \item[dir]
!    name of the subdirectory containing the data
!   \item[use\_timestamp]
!    flag to indicate whether the directories 
!    should be timestamped or not
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[yr]
!    4 digit year
!   \item[mo]
!    integer value of month (1-12)
!   \item[da]
!    day of the month (1-31)
!   \item[hr]
!    hour of the day (0-23)
!  \end{description}
!EOP
  character(3), parameter :: fhemi(2) = (/'NH_','SH_'/)
  character*10 :: ftime1, ftime2

  write (UNIT=ftime2, FMT='(i4, i2.2, i2.2, i2.2)') yr, mo, da, hr

  if (use_timestamp .eq. 1) then 
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     fname = trim(rootdir) // ftime1 // trim(dir) // '/WWMCA_TOTAL_CLOUD_PCT_' // &
             fhemi(hemi) // ftime2
  else
     fname = trim(rootdir) //  '/'   // trim(dir) // '/WWMCA_TOTAL_CLOUD_PCT_' // &
             fhemi(hemi) // ftime2
  endif

end subroutine agrmet_cdfs_totalcld_filename
