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
! !ROUTINE: get_gdas
!  \label{get_gdas}
!
! !REVISION HISTORY:
!  08 Dec 2000: Urszula Jambor; Rewrote geteta.f in fortran90 to 
!               use GDAS in GLDAS
!  09 Apr 2001: Urszula Jambor; Added capability of using
!               DAAC forcing data every
!               6 hours, rather than every 3 hours.
!  30 May 2001: Urszula Jambor; Changed forcing used: T,q,u fields taken 
!               from F00 & F03 files, radiation and precip. fields taken 
!               from F06 & F03 (F03 fields are subtracted out from F06)
!  11 Dec 2003: Sujay Kumar; Specification of elevation correction options
!  16 Mar 2008: Sujay Kumar; Enabled the capability to read the 6 and 9hr
!               forecasts as backups
!  29 Apr 2010: Sujay Kumar: Fixed the problems in mixing the instantaneous
!               and time averaged GDAS forecast files. 
!
! !INTERFACE:
subroutine get_gdas(n, findex)
! !USES:

  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,     only : LIS_tick, LIS_get_nstep
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use gdas_forcingMod,    only : gdas_struc
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, GDAS forcing. 
!  The idea is to open either the 00, 03 and 06 forecast file associated with
!  the most recent GDAS assimilation (available every 6 hours).
!  Precipitation rates and radiation fluxes will be taken from the F03 and 
!  F06 files, since averages are provided. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 3 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The strategy for missing data is to first look for 9hr forecasts and if not
!  go backwards up to 10 days to get forcing at the same time of day.
!
!  See additional notes in create\_gdasfilename(\ref{create_gdasfilename}).
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the GDAS data times
!  \item[create\_gdasfilename](\ref{create_gdasfilename}) \newline
!    Puts together appropriate file name for either 3 or 6 hour intervals
!  \item[read\_gdas](\ref{read_gdas}) \newline
!      Interpolates GDAS data to LIS grid
!   \item[read\_gdas\_elev](\ref{read_gdas_elev}) \newline
!    reads the native elevation of the GDAS data to be used
!    for topographic adjustments to the forcing 
!   \item[gdas\_reset\_interp\_input](\ref{gdas_reset_interp_input} \newline
!    resets the neighbours and weights arrays upon a grid change
!  \end{description}
!EOP

  integer :: t,f
  integer :: ferror
  integer :: try
  integer, parameter :: ndays = 10  ! # of days to look back for forcing data
  integer :: order     ! 1 indicates lesser interpolation boundary time
                       ! 2 indicates greater interpolation boundary time
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  real    :: ts1, ts2
  integer :: movetime  ! 1=move time2 into time1
  integer :: nforce    ! Number of forcing variables for model init. option
  real*8  :: timenow, time1, time2
  real*8  :: dumbtime1, dumbtime2
  real    :: gmt1, gmt2
  character(len=LIS_CONST_PATH_LEN) :: name00, name03, name06
  logical :: file_exists1, file_exists2, file_exists3
  real :: gridDesci(50)
  integer :: nstep
  logical :: F06flag
  integer :: status

  gdas_struc(n)%findtime1=0
  gdas_struc(n)%findtime2=0
  
  movetime=0
  
  !-----------------------------------------------------------------
  ! Determine the correct number of forcing variables
  !-----------------------------------------------------------------
  nstep = LIS_get_nstep(LIS_rc,n)
  nforce = gdas_struc(n)%nmif
  
  if ( LIS_rc%tscount(n).eq.1 .or. LIS_rc%rstflag(n)== 1) then
     gdas_struc(n)%findtime1=1
     gdas_struc(n)%findtime2=1
     movetime=0        ! movetime is not properly set at time-step = 1
     LIS_rc%rstflag(n) = 0
  endif
  
  !-----------------------------------------------------------------
  ! Determine required GDAS data times 
  ! (previous assimilation, current & future assimilation hours)
  ! The adjustment of the hour and the direction will be done
  ! in the subroutines that generate the names
  !-----------------------------------------------------------------
  yr1 = LIS_rc%yr  !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick( timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
  
  yr1 = LIS_rc%yr  !previous assimilation/forecast hour
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = 3*(int(real(LIS_rc%hr)/3.0))
  mn1 = 0
  ss1 = 0
  ts1 = 0
  call LIS_tick( time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
  
  yr2 = LIS_rc%yr  !next assimilation/forecast hour
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = 3*(int(real(LIS_rc%hr)/3.0))
  mn2 = 0
  ss2 = 0
  ts2 = 3*60*60
  call LIS_tick( time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )
  !-----------------------------------------------------------------
  ! Use these if need to roll back time.
  !-----------------------------------------------------------------
  dumbtime1 = time1
  dumbtime2 = time2
  if ( nstep == 0 .or. nstep == 1 ) then
     gdas_struc(n)%gdastime1 = time1
     gdas_struc(n)%gdastime2 = time2
  endif
  !-----------------------------------------------------------------
  ! Check to see if current time (timenow) has crossed past gdastime2,
  ! requiring that both gdastime2 be assigned to gdastime1 and a new
  ! gdastime2 be set 3 or 6 hours ahead of the current gdastime2.
  !-----------------------------------------------------------------
  if ( timenow > gdas_struc(n)%gdastime2 ) then
     movetime  = 1
     gdas_struc(n)%findtime2 = 1
  end if

  ! Note that the logic for determining when to change the GDAS
  ! grid must be updated.
  ! Consider the change to the T254 grid at 2002-10-29T12:00:00.
  ! This grid change does not occur until LIS' time is greater than
  ! or equal to 2002-10-29T12:00:00.  When LIS' time advances to hour 09,
  ! LIS will read ahead to get hour 12 data.  For this case, LIS needs
  ! instantaneous values from the 12f00 file.  It also needs averaged
  ! values from the 06f03 and 06f06 files.  At hour 09, LIS is still on
  ! the T170 grid, and it will properly read the 06f03 and 06f06 files,
  ! but it will not properly read the 12f00 file.
  !
  ! Updating the grid early will allow LIS to properly read the 12f00
  ! file, but then it would not read the 06f03 and 06f06 files properly.
  !
  ! The current stategy is to leave the below logic as-is.  Here
  ! LIS is not correctly read the 12f00 file, and it will roll back
  ! by one day to search for corresponding data.  Once LIS' time
  ! advances to 2002-10-29T12:00:00, and it updates the GDAS, then
  ! all GDAS files will be on the same grid, and LIS will be able
  ! read them correctly.
  !
  ! Note that there is an additional problem to consider.  LIS does not
  ! take grid changes into account when rolling back.  Say that LIS has
  ! just updated to the T254 grid.  LIS is trying to read 18f00 data, but
  ! that file is missing.  LIS will roll back by one day to read the
  ! 18f00 data from day 2002-10-28, which is on the T170 grid.  Here
  ! LIS will misread the file, causing another roll back until the 10-day
  ! roll back error condition is reached.

!= T126:
  if ( time1 >= gdas_struc(n)%griduptime1 .and. &
       time1 < gdas_struc(n)%griduptime2 .and. & 
       gdas_struc(n)%gridchange1) then 

     write(LIS_logunit,*) '[INFO] GDAS -- changing grid to 1991 - 2000 (T126)'

     gdas_struc(n)%ncold = 384
     gdas_struc(n)%nrold = 190
     gdas_struc(n)%mi = gdas_struc(n)%ncold*gdas_struc(n)%nrold
     !-------------------------------------------------------------------
     ! Reinitialize the weights and neighbors
     !-------------------------------------------------------------------
     gridDesci = 0 
     gridDesci(1) = 4
     gridDesci(2) = 384
     gridDesci(3) = 190
     gridDesci(4) = 89.277
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.277
     gridDesci(8) = -0.9375
     gridDesci(9) = 0.9375
     gridDesci(10) = 95
     gridDesci(20) = 0

! KRA - Previous:
!     if(LIS_rc%met_ecor(findex).eq."lapse-rate" .or. &
!          LIS_rc%met_ecor(findex) .eq. "lapse-rate and slope-aspect") then 
! Including MicroMet downscaling:
     if(LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
          LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" .or. &
          LIS_rc%met_ecor(findex) == "micromet" ) then
        call read_gdas_elev(n,findex, 1)
     endif

     call gdas_reset_interp_input(n, findex, gridDesci)
     gdas_struc(n)%gridchange1 = .false.

!= T170:
  elseif ( time1 >= gdas_struc(n)%griduptime2 .and. &
       time1 < gdas_struc(n)%griduptime3 .and. & 
       gdas_struc(n)%gridchange2) then 

     write(LIS_logunit,*) '[INFO] GDAS -- changing grid to 2000 - 2002 (T170)'

     gdas_struc(n)%ncold = 512
     gdas_struc(n)%nrold = 256
     gdas_struc(n)%mi = gdas_struc(n)%ncold*gdas_struc(n)%nrold
     !-------------------------------------------------------------------
     ! Reinitialize the weights and neighbors
     !-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 512
     gridDesci(3) = 256
     gridDesci(4) = 89.463
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.463
     gridDesci(8) = -0.703125
     gridDesci(9) = 0.703125
     gridDesci(10) = 128
     gridDesci(20) = 0.0

! KRA - Previous:
!     if(LIS_rc%met_ecor(findex).eq."lapse-rate" .or. &
!          LIS_rc%met_ecor(findex) .eq. "lapse-rate and slope-aspect") then 
! Including MicroMet downscaling:
     if(LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
          LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" .or. &
          LIS_rc%met_ecor(findex) == "micromet" ) then
        call read_gdas_elev(n,findex, 2)
     endif

     call gdas_reset_interp_input(n, findex, gridDesci)
     gdas_struc(n)%gridchange2 = .false.

!= T254:
  elseif ( time1 >= gdas_struc(n)%griduptime3 .and. & 
       time1 < gdas_struc(n)%griduptime4 .and. &
       gdas_struc(n)%gridchange3) then 

     write(LIS_logunit,*) '[ERR] GDAS -- changing grid to 2002 - 2005 (T254)'

     gdas_struc(n)%ncold = 768
     gdas_struc(n)%nrold = 384
     gdas_struc(n)%mi = gdas_struc(n)%ncold*gdas_struc(n)%nrold
     !-------------------------------------------------------------------
     ! Reinitialize the weights and neighbors
     !-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 768
     gridDesci(3) = 384
     gridDesci(4) = 89.642
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.642
     gridDesci(8) = -0.46875
     gridDesci(9) = 0.46875
     gridDesci(10) = 192
     gridDesci(20) = 0.0

! KRA - Previous:
!     if(LIS_rc%met_ecor(findex).eq."lapse-rate" .or. &
!          LIS_rc%met_ecor(findex) .eq. "lapse-rate and slope-aspect") then 
! Including MicroMet downscaling:
     if(LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
          LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" .or. &
          LIS_rc%met_ecor(findex) == "micromet" ) then
        call read_gdas_elev(n,findex, 3)
     endif

     call gdas_reset_interp_input(n, findex, gridDesci)
     gdas_struc(n)%gridchange3 = .false.

!= T382:
  elseif ( time1 >= gdas_struc(n)%griduptime4 .and. &
       time1 < gdas_struc(n)%griduptime5 .and. & 
       gdas_struc(n)%gridchange4) then 

     write(LIS_logunit,*) '[INFO] GDAS -- changing grid to 2005 - 2010 (T382)'

     gdas_struc(n)%ncold = 1152
     gdas_struc(n)%nrold = 576
     gdas_struc(n)%mi = gdas_struc(n)%ncold*gdas_struc(n)%nrold
     !-------------------------------------------------------------------
     ! Reinitialize the weights and neighbors
     !-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 1152
     gridDesci(3) = 576
     gridDesci(4) = 89.761
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.761
     gridDesci(8) = -0.3125
     gridDesci(9) = 0.3125
     gridDesci(10) = 288
     gridDesci(20) = 0.0

! KRA - Previous:
!     if(LIS_rc%met_ecor(findex).eq."lapse-rate" .or. &
!          LIS_rc%met_ecor(findex) .eq. "lapse-rate and slope-aspect") then 
! Including MicroMet downscaling:
     if(LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
          LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" .or. &
          LIS_rc%met_ecor(findex) == "micromet" ) then
        call read_gdas_elev(n,findex, 4)
     endif

     call gdas_reset_interp_input(n, findex, gridDesci)
     gdas_struc(n)%gridchange4 = .false.

!= T574:
  elseif ( time1 >= gdas_struc(n)%griduptime5 .and. &
       time1 < gdas_struc(n)%griduptime6 .and. &
       gdas_struc(n)%gridchange5) then

     write(LIS_logunit,*) '[INFO] GDAS -- changing grid to 2010 - 2015 (T574) --'

     gdas_struc(n)%ncold = 1760
     gdas_struc(n)%nrold = 880
     gdas_struc(n)%mi = gdas_struc(n)%ncold*gdas_struc(n)%nrold
     !-------------------------------------------------------------------
     ! Reinitialize the weights and neighbors
     !-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 1760
     gridDesci(3) = 880
     gridDesci(4) = 89.844
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.844
     gridDesci(8) = -0.204545454545455
     gridDesci(9) = 0.204545454545455
     gridDesci(10) = 440
     gridDesci(20) = 0.0

! KRA - Previous:
!     if(LIS_rc%met_ecor(findex).eq."lapse-rate" .or. &
!          LIS_rc%met_ecor(findex) .eq. "lapse-rate and slope-aspect") then 
! Including MicroMet downscaling:
     if(LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
          LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" .or. &
          LIS_rc%met_ecor(findex) == "micromet" ) then
        call read_gdas_elev(n,findex, 5)
     endif

     call gdas_reset_interp_input(n, findex, gridDesci)
     gdas_struc(n)%gridchange5 = .false.

!= T1534:
  elseif ( time1 >= gdas_struc(n)%griduptime6 .and. &
       gdas_struc(n)%gridchange6 ) then

     write(LIS_logunit,*) '[INFO] GDAS -- changing grid to 2015 (T1534) --'

     gdas_struc(n)%ncold = 3072
     gdas_struc(n)%nrold = 1536
     gdas_struc(n)%mi = gdas_struc(n)%ncold*gdas_struc(n)%nrold
     !-------------------------------------------------------------------
     ! Reinitialize the weights and neighbors
     !-------------------------------------------------------------------
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = gdas_struc(n)%ncold
     gridDesci(3) = gdas_struc(n)%nrold
     gridDesci(4) = 89.91000
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.91000
     gridDesci(8) = -0.1171875
     gridDesci(9) = 0.1171875
     gridDesci(10) = 768.0
     gridDesci(20) = 0.0

! KRA - Previous:
!     if(LIS_rc%met_ecor(findex).eq."lapse-rate" .or. &
!          LIS_rc%met_ecor(findex) .eq. "lapse-rate and slope-aspect") then 
! Including MicroMet downscaling:
     if(LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
          LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" .or. &
          LIS_rc%met_ecor(findex) == "micromet" ) then
        call read_gdas_elev(n,findex, 6)
     endif

     call gdas_reset_interp_input(n, findex, gridDesci)
     gdas_struc(n)%gridchange6 = .false.

  endif

  !-----------------------------------------------------------------
  ! Establish bookend 1: The code sets up the filenames required to 
  ! read both instantaneous and time average fields.  
  !-----------------------------------------------------------------
  if ( gdas_struc(n)%findtime1 == 1 ) then  !get new time1 from the past
     write(LIS_logunit,*) '[INFO] Getting new time1 data'

     ferror = 0
     order = 1
     try = 0
     ts1 = -24*60*60

     do
        if ( ferror /= 0 ) exit
        try = try + 1
        status = 0 

        call create_gdasfilename(order, name00, name03, name06, F06flag, &
             gdas_struc(n)%gdasdir, yr1, mo1, da1, hr1)
        inquire(file=name00,exist=file_exists1) 
        inquire(file=name03,exist=file_exists2)
        if ( .not. file_exists1 ) then
           status = 1
           write(LIS_logunit,*) '[ERR] GDAS file1 (I) ',trim(name00), &
                'does not exist. (f00) (try = ', try, ')'
        endif
        if ( .not. file_exists2 ) then
           status = 1
           write(LIS_logunit,*) '[ERR] GDAS file1 (A) ',trim(name03), &
                'does not exist. (f03) (try = ', try, ')'
        endif

        if ( ( .not. file_exists1 ) .or. ( .not. file_exists2 ) ) then 
           write(LIS_logunit,*) '[INFO] looking for backup GDAS f09 files.'
           call create_gdasf9_filename(order, name00, name03, &
                gdas_struc(n)%gdasdir, &
                yr1, mo1, da1, hr1, status)

           if ( status == 1 ) then
              write(LIS_logunit,*) '[WARN] backup GDAS f09 files are not ', &
                   'available for this hour'
           else
              inquire(file=name00,exist=file_exists1) 
              inquire(file=name03,exist=file_exists2)
              if ( .not. file_exists1 ) then
                 status = 1
                 write(LIS_logunit,*) '[ERR] backup GDAS f09 file1 (I) ', &
                      trim(name00), &
                      'does not exist. (f00) (try = ',try,')'
              endif
              if ( .not. file_exists2 ) then
                 status = 1
                 write(LIS_logunit,*) '[ERR] backup GDAS f09 file1 (A) ', &
                      trim(name03), &
                      'does not exist. (f03) (try = ',try,')'
              endif

              if ( file_exists1 .and. file_exists2 ) then 
                 status = 0 
              endif
           endif
        endif

        if ( F06flag ) then
           inquire(file=name06,exist=file_exists3)
           if ( .not. file_exists3 ) then
              status = 1
              write(LIS_logunit,*) '[ERR] GDAS file1 (A) ',trim(name06), &
                   'does not exist. (f06) (try = ', try, ')'
           endif
        endif

        if ( status == 0 ) then 
           write(LIS_logunit,*) '[INFO] Reading GDAS file1 (I) ',trim(name00)
           write(LIS_logunit,*) '[INFO] Reading GDAS file1 (A) ',trim(name03)
           if ( F06flag ) then 
              write(LIS_logunit,*) '[INFO] Reading GDAS file1 (A)',trim(name06)
           endif
           call read_gdas(order, n, findex,  &
                name00, name03, name06, F06flag, ferror, try)
           if ( ferror == 1 ) then  
              ! successfully retrieved forcing data
              gdas_struc(n)%gdastime1 = time1
           else  
              ! ferror still=0, so roll back one day and start again
              write(LIS_logunit,*) '[ERR] Error reading files. ', &
                   'Rolling back one day. (try = ', try, ')'
              call LIS_tick( dumbtime1, doy1, gmt1, &
                   yr1, mo1, da1, hr1, mn1, ss1, ts1 )
           end if
        else
           ! cannot find files, so roll back one day and start again
           write(LIS_logunit,*) '[WARN] Cannot find files. ', &
                'Rolling back one day. (try = ',try, ')'
           call LIS_tick( dumbtime1, doy1, gmt1, &
                yr1, mo1, da1, hr1, mn1, ss1, ts1 )
        endif
        if ( try > ndays ) then  
           ! data gap exceeds 10 days so stop
           write(LIS_logunit,*) '[ERR] GDAS data gap exceeds 10 days on file 1'
           call LIS_endrun
        endif
     end do
  end if

  !-----------------------------------------------------------------
  ! Establish gdastime2
  !-----------------------------------------------------------------
  if ( movetime == 1 ) then  
     gdas_struc(n)%gdastime1 = gdas_struc(n)%gdastime2
     gdas_struc(n)%findtime2 = 1  
     do f = 1, LIS_rc%met_nf(findex)
        do t = 1, LIS_rc%ngrid(n)
           gdas_struc(n)%metdata1(f,t) = gdas_struc(n)%metdata2(f,t)
        end do
     end do
  end if

  if ( gdas_struc(n)%findtime2 == 1 ) then  
     write(LIS_logunit,*) '[INFO] Getting new time2 data'

     ferror = 0
     order = 2
     try = 0
     ts2 = -24*60*60

     do
        if ( ferror /= 0 ) exit
        try = try + 1
        status = 0

        call create_gdasfilename(order, name00, name03,name06, F06flag, &
             gdas_struc(n)%gdasdir, yr1, mo1, da1, hr1)
        inquire(file=name00,exist=file_exists1) 
        inquire(file=name03,exist=file_exists2)
        if ( .not. file_exists1 ) then
           status = 1
           write(LIS_logunit,*) '[WARN] GDAS file2 (I) ',trim(name00), &
                'does not exist. (f00) (try = ', try, ')'
        endif
        if ( .not. file_exists2 ) then
           status = 1
           write(LIS_logunit,*) '[ERR] GDAS file2 (A) ',trim(name03), &
                'does not exist. (f03) (try = ', try, ')'
        endif

        if ( ( .not. file_exists1 ) .or. ( .not. file_exists2 ) ) then
           write(LIS_logunit,*) 'MSG: looking for backup GDAS f09 files.'
           call create_gdasf9_filename(order, name00, name03, &
                gdas_struc(n)%gdasdir, &
                yr1, mo1, da1, hr1, status)

           if ( status == 1 ) then
              write(LIS_logunit,*) '[WARN] backup GDAS f09 files are not ', &
                   'available for this hour'
           else
              inquire(file=name00,exist=file_exists1) 
              inquire(file=name03,exist=file_exists2)
              if ( .not. file_exists1 ) then
                 status = 1
                 write(LIS_logunit,*) '[ERR] GDAS backup f09 file2 (I) ', &
                      trim(name00), &
                      'does not exist. (f00) (try = ',try,')'
              endif
              if ( .not. file_exists2 ) then
                 status = 1
                 write(LIS_logunit,*) '[ERR] GDAS backup f09 file2 (A) ', &
                      trim(name03), &
                      'does not exist. (f03) (try = ',try,')'
              endif

              if ( file_exists1 .and. file_exists2 ) then 
                 status = 0 
              endif
           endif
        endif

        if ( F06flag ) then
           inquire(file=name06,exist=file_exists3)
           if ( .not. file_exists3 ) then
              status = 1
              write(LIS_logunit,*) '[ERR] GDAS file2 (A) ',trim(name06), &
                   'does not exist. (f06) (try = ', try, ')'
           endif
        endif

        if ( status == 0 ) then 
           write(LIS_logunit,*) '[INFO] Reading GDAS file2 (I) ',trim(name00)
           write(LIS_logunit,*) '[INFO] Reading GDAS file2 (A) ',trim(name03)
           if( F06flag ) then 
              write(LIS_logunit,*) '[INFO] Reading GDAS file2 (A)',trim(name06)
           endif
           call read_gdas(order, n, findex, &
                name00, name03, name06, F06flag, ferror, try)
           if ( ferror == 1 ) then  
              ! successfully retrieved forcing data
              gdas_struc(n)%gdastime2 = time2
           else  
              ! ferror still=0, so roll back one day and start again
              write(LIS_logunit,*) '[ERR] Error reading files. ', &
                   'Rolling back one day. (try = ', try, ')'
              call LIS_tick( dumbtime2, doy2, gmt2, yr1, mo1, &
                   da1, hr1, mn1, ss1, ts2 )
           end if
        else
           ! cannot find files, so roll back one day and start again
           write(LIS_logunit,*) '[ERR] Cannot find files. ', &
                'Rolling back one day. (try = ',try, ')'
           call LIS_tick( dumbtime2, doy2, gmt2, yr1, mo1, &
                da1, hr1, mn1, ss1, ts2 )
        endif
        if ( try > ndays ) then  
           ! data gap exceeds 10 days so stop
           write(LIS_logunit,*) '[ERR] GDAS data gap exceeds 10 days on file 2'
           call LIS_endrun
        end if
     end do
  end if

  !-----------------------------------------------------------------
  ! Check for negative values in radiation 
  !-----------------------------------------------------------------
  do f = 1,LIS_rc%met_nf(findex)  
     if ( (f == 3) .or. (f == 4) ) then
        do t=1,LIS_rc%ngrid(n)
           if ( (gdas_struc(n)%metdata2(f,t) /= -9999.9) .and.  &
                (gdas_struc(n)%metdata2(f,t) < 0)) then
              gdas_struc(n)%metdata2(f,t) = (-1) * gdas_struc(n)%metdata2(f,t)
           endif
           if ( (gdas_struc(n)%metdata1(f,t) /= -9999.9) .and.  &
                (gdas_struc(n)%metdata1(f,t) < 0)) then
              gdas_struc(n)%metdata1(f,t) = (-1) * gdas_struc(n)%metdata1(f,t)
           endif
        enddo
     endif
  enddo
  !  endif
end subroutine get_gdas

