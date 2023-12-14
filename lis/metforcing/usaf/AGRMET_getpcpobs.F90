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
! !ROUTINE: AGRMET_getpcpobs
!  \label{AGRMET_getpcpobs}
!
! !REVISION HISTORY:
!   17 Aug 99  Initial version (PRECIP getobs.f)......Capt Hidalgo/Mr Moore/DNXM
!   21 Feb 01  Major rewrite to support running precip in 6-hourly
!              cycles and for new precip decoder...................Mr Gayno/DNXM
!   07 Jan 03  Expanded southern hemisphere search range for rain gauge obs
!              from synoptic time plus/minus 1 hour to -2 to +1 hours because
!              Australia reports at local time, not Zulu time......Mr Gayno/DNXM
!   29 Oct 05  Integrated into LIS and renamed.............Sujay Kumar/NASA/GSFC
!   16 Mar 07  Changed to get path names from agrmet_struc.
!              Added use_timestamp option......................Ted Lewiston/DNXM
!   19 Dec 07  Simplified filename creation ........... Marv Freimund(LMMS)/A8TM
!   27 Jan 09  Set 6 & 12 hour precip grids, which are now a single lat/lon grid,
!              to missing beofre the loop for hemispheres ..Chris Franks/2WG/WEA
!   30 Jul 10  Replaced parameter isize for obs array sizes with max_pcpobs 
!              from agrmet_struc which is configurable.  Made obs arrays 
!              allocatable and use max_pcpobs as the size..Chris Franks/16WS/WXE/SEMS
!    9 Sep 10  Modified to read JMOBS derived files........Chris Franks/16WS/WXE/SEMS
!    7 Feb 11  Enable use of either JMOBS or CDMS obs......Chris Franks/16WS/WXE/SEMS
!   11 May 11  Store obs from 3,9,15,& 21Z for India and Sri Lanka in a 
!              new array and pass to processobs............Chris Franks/16WS/WXE/SEMS
!   29 Aug 23  Call LIS_alert if a preobs file is missing..............Eric Kemp/NASA
!
! !INTERFACE:
subroutine AGRMET_getpcpobs(n, j6hr, month, prcpwe, &
     use_twelve, p6, p12, alert_number, precip6, precip12,pcp_src)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only    : LIS_tick, LIS_julhr_date
  use LIS_logMod, only        : LIS_logunit, LIS_alert
  use AGRMET_forcingMod, only : agrmet_struc
  use USAF_bratsethMod, only: USAF_ObsData, USAF_setbratsethprecipstats

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n 
  real,  intent(out)     :: prcpwe(LIS_rc%lnc(n), LIS_rc%lnr(n),4)
  real,  intent(inout)   :: p6(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real,  intent(inout)   :: p12(LIS_rc%lnc(n), LIS_rc%lnr(n))
  logical, intent(in)    :: use_twelve
  integer, intent(in)    :: j6hr
  integer, intent(in)    :: month
  integer,  intent(inout):: alert_number
  type(USAF_ObsData), intent(inout) :: precip6
  type(USAF_ObsData), intent(inout) :: precip12
  character(len=6), intent(in) :: pcp_src(4)

! !DESCRIPTION:
!  
!  This routine retrieves precip observations from the JMOBS database.
!
!    \textbf{Method}
!   
!    - Initialize necessary arrays. \newline
!      - Note: in the northern hemisphere, the database is called
!        once for current Julian hour.  In the southern hemisphere
!        it is called multiple times - at the current Julian hour;
!        one hour before and after the current Julian hour; and two hours
!        before the current Julian hour.  This is done because
!        Australia does not report at the synoptic times but rather
!        at local time (and the local reporting time changes with
!        season.  Confusing, huh?).  \newline
!    - Call getpwe to calculate rainfall based on current or past
!      weather reports. \newline
!    - At 03, 09, 15 and 21 utc, call storeobs\_offhour to perform some
!      preprocessing of the rain gauge data from India and Sri Lanka
!    - At 00, 06, 12 and 18 utc, call storeobs to perform some
!      preprocessing of the rain gauge data. then call processobs
!      to put the rain gauge data on the grid. \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[alert\_number]  number of alert messages
!   \item[network]       JMOBS network
!   \item[plat\_id]      JMOBS platform id
!   \item[bsn]           WMO block station number array
!   \item[duration]      Valid period of precipitation
!   \item[endjul]        Ending Julian hour of the one-hour loop
!   \item[hemi]          Hemisphere flag (1 = N, 2 = S)
!   \item[ierr1,ierr2,ierr3]          I/O error status
!   \item[ilat]          Array of observation latitudes
!   \item[ilon]          Array of observation longitudes
!   \item[imax]          Grid dimension - east/west direction
!   \item[jmax]          Grid dimension - north/south direction
!   \item[j1hr]          One-hourly loop index
!   \item[j3hr]          Three-hourly loop index
!   \item[j6hr]          Beginning of this 6 hourly processing period
!                       in Julian hours
!   \item[k]             Loop counter
!   \item[mscprc]        Miscellaneous precip amounts
!   \item[month]         Current month
!   \item[nsize]         Number of obs returned from database
!   \item[nsize3]        Number of obs stored from India and Sri Lanka for 3,9,15,\&21Z
!   \item[obs]           Array of processed observations \newline
!       amt6        Six hourly precip amount \newline
!       amt12       Twelve hourly precip amount \newline
!       amt24       Twenty-four hourly precip amount \newline
!       lat         Latitude \newline
!       lon         Longitude \newline
!       wmonum      WMO block station number \newline
!   \item[obs3]          Array of processed observations for 3,9,15,\&21Z
!                        for India and Sri Lanka - same structure as obs \newline
!   \item[oldd\_pwe]      Stores distances between a grid point
!                       and its nearest neighbor observation
!   \item[pastwx]        Past weather at station
!   \item[pathpcp]       Directory path of precip data directory
!   \item[prcpwe]        Present/past weather estimate on the grid
!   \item[preswx]        Present weather at station
!   \item[wmoblk]        Present weather at station
!   \item[p6]            6-hourly precip amts on the grid
!   \item[p12]           12-hourly precip amts on the grid
!   \item[quad9r]        Missing value (9999)
!   \item[sixprc]        Array of six hourly precip amounts
!   \item[startjul]      Starting julian hour of one hour loop
!   \item[stncnt]        Number of observations for this time period
!   \item[twfprc]        Array of 24 hourly precip amounts
!   \item[use\_twelve]    Logical flag to use the 12-hourly precip
!                   amts or the 6 hourly amts
!   \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[julhr\_date](\ref{LIS_julhr_date}) \newline
!    converts julian hour to a date format
!  \item[agrmetpcpobsfilename](\ref{agrmetpcpobsfilename}) \newline
!    generates the name of the precip observations file
!  \item[AGRMET\_getpwe](\ref{AGRMET_getpwe}) \newline
!    get present weather estimate
!  \item[AGRMET\_storeobs](\ref{AGRMET_storeobs}) \newline
!   store observations 
!  \item[setpathpcpobs](\ref{setpathpcpobs}) \newline
!   set the path to write the observations analysis
!  \item[AGRMET\_processobs](\ref{AGRMET_processobs}) \newline
!    process observations
!  \end{description}
!
!EOP
  real                      :: oldd_pwe(LIS_rc%lnc(n), LIS_rc%lnr(n))
  integer                   :: cdms_count
  integer                   :: hemi
  integer                   :: i
  character*100             :: filename
  real, parameter           :: quad9r = -9999.0   
  integer                   :: j1hr
  integer                   :: j3hr
  integer                   :: startjul
  integer                   :: endjul
  integer                      :: stncnt
  integer, allocatable         :: twfprc(:)
  integer, allocatable         :: duration(:)
  integer, allocatable         :: sixprc(:)
  integer, allocatable         :: mscprc(:)
  integer, allocatable         :: ilat(:)
  integer, allocatable         :: ilon(:)
  integer, allocatable         :: bsn(:)
  character*10, allocatable    :: network(:)
  character*10, allocatable    :: plat_id(:)
  integer, allocatable         :: pastwx(:)
  integer, allocatable         :: preswx(:)
  integer, allocatable         :: wmoblk(:)
  integer                      :: nsize
  integer                      :: nsize3
  integer                      :: yr,mo,da,hr
  integer                      :: ierr1, ierr2, ierr3
  integer                      :: k
  logical                      :: cdms_flag
  character(255) :: message(20)

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

  type(rain_obs), allocatable   :: obs(:)
  type(rain_obs), allocatable   :: obs3(:)

!-----------------------------------------------------------------------
! Allocate observation arrays
!-----------------------------------------------------------------------

  allocate ( twfprc (agrmet_struc(n)%max_pcpobs))
  allocate ( duration (agrmet_struc(n)%max_pcpobs))
  allocate ( sixprc (agrmet_struc(n)%max_pcpobs))
  allocate ( mscprc (agrmet_struc(n)%max_pcpobs))
  allocate ( ilat (agrmet_struc(n)%max_pcpobs))
  allocate ( ilon (agrmet_struc(n)%max_pcpobs))
  allocate ( bsn (agrmet_struc(n)%max_pcpobs))
  allocate ( network (agrmet_struc(n)%max_pcpobs))
  allocate ( plat_id (agrmet_struc(n)%max_pcpobs))
  allocate ( pastwx (agrmet_struc(n)%max_pcpobs))
  allocate ( preswx (agrmet_struc(n)%max_pcpobs))
  allocate ( wmoblk (agrmet_struc(n)%max_pcpobs))
  allocate ( obs (agrmet_struc(n)%max_pcpobs))
  allocate ( obs3 (agrmet_struc(n)%max_pcpobs))

!-----------------------------------------------------------------------
! Set the 6 & 12 hour grids missing prior to the hemisphere loop
!-----------------------------------------------------------------------

  p6 = quad9r
  p12 = quad9r
  do hemi =1,2

     if (use_twelve) then
        k = 3
     else 
        k = 1
     end if

     nsize3        = 0

     THREE_HOURLY : do j3hr = j6hr+3, j6hr+6, 3     

        ! Set Bratseth error statistics based on source of background 
        ! field.     
        call USAF_setBratsethPrecipStats(pcp_src(k),n)

        oldd_pwe      = 999.0
        obs%wmonum    = -9999
        obs%lat       = -9999.9
        obs%lon       = -9999.9
        obs%amt24     = -99999999
        obs%amt12     = -99999999
        obs%amt6      = -99999999
        obs%amtmsc    = -99999999
        stncnt        = 0
  
!-----------------------------------------------------------------------
!       Retrieve observations from JMOBS database.
!-----------------------------------------------------------------------

           call LIS_julhr_date( j3hr, yr, mo, da, hr)

           call agrmetpcpobsfilename(filename,                            &
                agrmet_struc(n)%agrmetdir,           &
                agrmet_struc(n)%cdmsdir,             &
                agrmet_struc(n)%use_timestamp,       &
                hemi, yr, mo, da, hr)
 
           write(LIS_logunit,*)'Reading PCP OBS: ',trim(filename)
           open(22,file=trim(filename),status='old', iostat=ierr1)
           if(ierr1.eq.0) then 
              nsize = 0 ! EMK
              read(22,*, iostat=ierr2) nsize
              if (ierr2 .ne. 0) nsize = 0 ! EMK

!-----------------------------------------------------------------------
!             If the number of obs is > the array sizes, write a warning
!             to the log and set back the number to avoid segfault.
!-----------------------------------------------------------------------

              if ( nsize .GT. agrmet_struc(n)%max_pcpobs ) then
     
                 write(LIS_logunit,*)' '
                 write(LIS_logunit,*)"******************************************************"
                 write(LIS_logunit,*)"* NUMBER OF RAIN GAUGE OBSERVATIONS EXCEEDS ARRAY SIZE."
                 write(LIS_logunit,*)"* NUMBER OF RAIN GAUGE OBS IS ", nsize
                 write(LIS_logunit,*)"* ARRAY SIZE IS ", agrmet_struc(n)%max_pcpobs
                 write(LIS_logunit,*)"* OBSERVATIONS BEYOND ARRAY SIZE WILL BE IGNORED."
                 write(LIS_logunit,*)"******************************************************"
                 write(LIS_logunit,*)' '

                 !EMK 20230829...Create alert file.
                 message(:) = ''
                 message(1) = '[WARN] Program:  LIS'
                 message(2) = '  Routine: AGRMET_getpcpobs'
                 message(3) = '  Too many rain gage reports in '// &
                      trim(filename)
!                 message(4) = '  Number of rain gage reports is '// nsize
                 write(message(5),'(A, I6)') &
                      '  Number of rain gage reports is ', nsize
                 !message(5) = '  Array size is '// &
                 !     agrmet_struc(n)%max_pcpobs
                 write(message(5),'(A, I6)') '  Array size is ', &
                      agrmet_struc(n)%max_pcpobs
                 message(6) = '  Observations beyond array size will be ignored'
                 message(7) = '  Increase number of AGRMET maximum precip obs in lis.config file!'
                 if (LIS_masterproc) then
                    alert_number = alert_number + 1
                    call LIS_alert('LIS.AGRMET_getpcpobs', &
                         alert_number, message)
                 end if

                 nsize = agrmet_struc(n)%max_pcpobs

     
              end if

              cdms_count = 0
              do i=1,nsize
                 read(22,*, iostat=ierr3) twfprc(i),&
                      duration(i), sixprc(i), mscprc(i),&
                      ilat(i), ilon(i), network(i), plat_id(i), pastwx(i), preswx(i), wmoblk(i)

!-----------------------------------------------------------------------
!                If this is a WMO Synoptic observation set the numeric
!                bsn field, else set it to 0.
!-----------------------------------------------------------------------

                 if (network(i) .eq. "WMO") then
                    read (plat_id(i), FMT='(I10)') bsn(i)
                 else if (network(i) .eq. "CDMS") then
                    read (plat_id(i), FMT='(I10)') bsn(i)
                    cdms_count = cdms_count + 1
                 else
                    bsn(i) = 0
                 end if
              
              enddo
              close(22)

              if (cdms_count .eq. nsize) then
                 cdms_flag = .TRUE.
              else
                 cdms_flag = .FALSE.
              end if

              if(ierr2.eq.0.and.ierr3.eq.0) then 
                 if(nsize > 0 ) then 
                    if(agrmet_struc(n)%pwswch.eq.1) then 
                       write(LIS_logunit,*)' '
                       write(LIS_logunit,*)'- CALCULATING PRES/PAST WX ESTIMATE'
                       write(LIS_logunit,*)' '
                       
                       call AGRMET_getpwe (n, nsize, agrmet_struc(n)%max_pcpobs, &
                            network, plat_id, ilat, ilon, &
                            month, pastwx, preswx, wmoblk, oldd_pwe, &
                            prcpwe(:,:,k), LIS_rc%lnc(n), LIS_rc%lnr(n))
                    end if
!-----------------------------------------------------------------------
!             At the synoptic times (00, 06, 12, and 18 utc) call
!             storeobs to preprocess the rain gauge data.
!-----------------------------------------------------------------------
                 
                    if ( mod(j3hr,6) == 0 ) then
                       
                       write(LIS_logunit,*)'- CALLING STOREOBS TO PROCESS RAIN GAUGE DATA', j3hr
                       write(LIS_logunit,*)' '
                       
                       call AGRMET_storeobs(nsize, nsize3, agrmet_struc(n)%max_pcpobs, &
                            obs, obs3, ilat, ilon, &
                            mscprc, sixprc, twfprc, network, plat_id, cdms_flag, bsn, &
                            duration, j3hr, stncnt)
                    
                    else
                       
                       write(LIS_logunit,*)'- CALLING STOREOBS_OFFHOUR TO PROCESS 3HOUR RAIN GAUGE DATA', j3hr
                       write(LIS_logunit,*)' '
                       
                       call AGRMET_storeobs_offhour(nsize, agrmet_struc(n)%max_pcpobs, &
                            obs3, ilat, ilon, &
                            mscprc, sixprc, twfprc, network, plat_id, cdms_flag, bsn, &
                            duration, nsize3)
                    
                    end if
                 
                 end if
              else

!-----------------------------------------------------------------------
!           If there is a bad database read, just move on to the
!           next time.  If all the reads end up being bad, then
!           the p6, p12 and prcpwe arrays will be filled with
!           missing flags.
!-----------------------------------------------------------------------

                 write(LIS_logunit,*)' '
                 write(LIS_logunit,*)'**********************************************'
                 write(LIS_logunit,*)'*** ERROR ON DATABASE READ.  ISTAT IS '
                 write(LIS_logunit,*)'**********************************************'
                 write(LIS_logunit,*)' '

                 !EMK 20230829...Create alert file.
                 message(:) = ''
                 message(1) = '[WARN] Program:  LIS'
                 message(2) = '  Routine: AGRMET_getpcpobs'
                 message(3) = '  Problem reading '// trim(filename)
                 if (LIS_masterproc) then
                    alert_number = alert_number + 1
                    call LIS_alert('LIS.AGRMET_getpcpobs', &
                         alert_number, message)
                 end if
              end if
           else
              write(LIS_logunit,*)' '
              write(LIS_logunit,*)'**********************************************'
              write(LIS_logunit,*)'*** ERROR ON DATABASE FILE DOES NOT EXIST'
              write(LIS_logunit,*) trim(filename)
              write(LIS_logunit,*)'**********************************************'
              write(LIS_logunit,*)' '

              !EMK 20230829...Create alert file.
              message(:) = ''
              message(1) = '[WARN] Program:  LIS'
              message(2) = '  Routine: AGRMET_getpcpobs'
              message(3) = '  Missing rain gage file '// trim(filename)
              if (LIS_masterproc) then
                 alert_number = alert_number + 1
                 call LIS_alert('LIS.AGRMET_getpcpobs', &
                      alert_number, message)
              end if

           endif
        
!-----------------------------------------------------------------------
!       Call processobs (at the synoptic times) to process and
!       place rain gauge data on the model grid.
!-----------------------------------------------------------------------

        if ( (mod(j3hr,6) == 0) .and. (stncnt > 0) ) then
           
           write(LIS_logunit,*)'- CALLING PROCESSOBS TO PLACE RAIN GAUGE DATA ON GRID'
           write(LIS_logunit,*)' '

           call LIS_julhr_date( j3hr, yr, mo, da, hr)
        
           call AGRMET_processobs(n, obs, agrmet_struc(n)%max_pcpobs, stncnt, hemi, &
                j3hr, LIS_rc%lnc(n), LIS_rc%lnr(n), p6, &
                agrmet_struc(n)%analysisdir,&
                p12, use_twelve, cdms_flag, &
                quad9r, alert_number, &
                precip6, precip12) 
           
        end if

        k = k + 1
     
     end do THREE_HOURLY
  end do

!-----------------------------------------------------------------------
! Deallocate the observation arrays
!-----------------------------------------------------------------------

  deallocate ( twfprc )
  deallocate ( duration )
  deallocate ( sixprc )
  deallocate ( mscprc )
  deallocate ( ilat )
  deallocate ( ilon )
  deallocate ( bsn )
  deallocate ( network )
  deallocate ( plat_id )
  deallocate ( pastwx )
  deallocate ( preswx )
  deallocate ( wmoblk )
  deallocate ( obs )

  return 
  
!-----------------------------------------------------------------------
!     format statements.
!-----------------------------------------------------------------------

6000 format(1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,a10,1x,a10,1x,i9,1x,i9,1x,i9)
  
end subroutine AGRMET_getpcpobs

!BOP
!
! !ROUTINE: agrmetpcpobsfilename
! \label{agrmetpcpobsfilename}
!
! !INTERFACE: 
subroutine agrmetpcpobsfilename(filename, rootdir, dir,                &
                                use_timestamp, hemi, yr, mo, da, hr)

  implicit none
! !ARGUMENTS:   
  character(*)        :: filename
  character(*)        :: rootdir
  character(*)        :: dir
  integer, intent(in) :: use_timestamp
  integer, intent(in) :: hemi
  integer, intent(in) :: yr, mo, da, hr
! 
! !DESCRIPTION: 
!  This routines generates the name of the precip observations file
!  by appending the hemisphere and timestamps to the root directory. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[filename]
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
!    day of the month
!   \item[hr]
!    hour of the day
!  \end{description}
!EOP

  character(2), parameter :: fhemi(2) = (/'nh','sh'/)
  character(10)           :: ftime1, ftime2

  write (UNIT=ftime2, FMT='(i4, i2.2, i2.2, i2.2)') yr, mo, da, hr

  if (use_timestamp .eq. 1) then 
     write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
     filename = trim(rootdir) // ftime1 // trim(dir) // '/preobs_' // &
          fhemi(hemi) // '.03hr.' // ftime2
  else
     filename = trim(rootdir) //   '/'  // trim(dir) // '/preobs_' // &
          fhemi(hemi) // '.03hr.' // ftime2
  endif

end subroutine agrmetpcpobsfilename

!BOP
! 
! !ROUTINE: setpathpcpobs
! \label{setpathpcpobs}
!
! !INTERFACE: 
subroutine setpathpcpobs(fname,dir,yr,mo,da,hr)

  implicit none
! !ARGUMENTS:   
  character(*)        :: fname
  character(*)        :: dir
  integer, intent(in) :: yr,mo,da,hr
! 
! !DESCRIPTION: 
!  This routines generates the name of the path to which processed 
!  precip observations are to be written. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[fname]
!    created filepath
!   \item[dir]
!    full path to the directory containing the data
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
!
  character(10) :: ftime1

  write (UNIT=ftime1, FMT='(a1, i4, i2.2, i2.2, a1)') '/',yr,mo,da,'/'

  fname = trim(dir) // ftime1 // 'PRECIP/'

end subroutine setpathpcpobs
