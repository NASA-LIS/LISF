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
! !ROUTINE: get_gdasT1534
!  \label{get_gdasT1534}
!
! !REVISION HISTORY:
!  20 June 2014: Sujay Kumar; initial implementation
! !INTERFACE:
subroutine get_gdasT1534(n, findex)
! !USES:

  use LIS_coreMod,        only : LIS_rc
  use LIS_timeMgrMod,     only : LIS_tick, LIS_get_nstep
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use gdasT1534_forcingMod,    only : gdasT1534_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, GDAST1534 forcing. 
!  The idea is to open either the 00, 03 and 06 forecast file associated with
!  the most recent GDAST1534 assimilation (available every 6 hours).
!  Precipitation rates and radiation fluxes will be taken from the F03 and 
!  F06 files, since averages are provided. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 3 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The strategy for missing data is to first look for 9hr forecasts and if not
!  go backwards up to 10 days to get forcing at the same time of day.
!
!  See additional notes in create\_gdasT1534filename(\ref{create_gdasT1534filename}).
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
!    determines the GDAST1534 data times
!  \item[create\_gdasT1534filename](\ref{create_gdasT1534filename}) \newline
!    Puts together appropriate file name for either 3 or 6 hour intervals
!  \item[read\_gdasT1534](\ref{read_gdasT1534}) \newline
!      Interpolates GDAST1534 data to LIS grid
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
  real*8  :: timenow, time1, time2, ttime
  real*8  :: dumbtime1, dumbtime2
  real    :: gmt1, gmt2
  character(len=LIS_CONST_PATH_LEN) :: name
  logical :: file_exists
  real :: gridDesci(50)
  integer :: nstep
  logical :: F06flag
  integer :: status

  gdasT1534_struc(n)%findtime1=0
  gdasT1534_struc(n)%findtime2=0
  
  movetime=0
  
  !-----------------------------------------------------------------
  ! Determine the correct number of forcing variables
  !-----------------------------------------------------------------
  nstep = LIS_get_nstep(LIS_rc,n)
  
  if ( LIS_rc%tscount(n).eq.1 .or. LIS_rc%rstflag(n)== 1) then
     gdasT1534_struc(n)%findtime1=1
     gdasT1534_struc(n)%findtime2=1
     movetime=0        ! movetime is not properly set at time-step = 1
     LIS_rc%rstflag(n) = 0
  endif
  
  !-----------------------------------------------------------------
  ! Determine required GDAST1534 data times 
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
  hr1 = 1*(int(real(LIS_rc%hr)/1.0))
  mn1 = 0
  ss1 = 0
  ts1 = 0
  call LIS_tick( time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
  
  yr2 = LIS_rc%yr  !next assimilation/forecast hour
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = 1*(int(real(LIS_rc%hr)/1.0))
  mn2 = 0
  ss2 = 0
  ts2 = 1*60*60
  call LIS_tick( time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )
  !-----------------------------------------------------------------
  ! Use these if need to roll back time.
  !-----------------------------------------------------------------
  dumbtime1 = time1
  dumbtime2 = time2
  if ( nstep == 0 .or. nstep == 1 ) then
     gdasT1534_struc(n)%gdasT1534time1 = time1
     gdasT1534_struc(n)%gdasT1534time2 = time2
  endif
  !-----------------------------------------------------------------
  ! Check to see if current time (timenow) has crossed past gdasT1534time2,
  ! requiring that both gdasT1534time2 be assigned to gdasT1534time1 and a new
  ! gdasT1534time2 be set 3 or 6 hours ahead of the current gdasT1534time2.
  !-----------------------------------------------------------------
  if ( timenow > gdasT1534_struc(n)%gdasT1534time2 ) then
     movetime  = 1
     gdasT1534_struc(n)%findtime2 = 1
  end if

  !-----------------------------------------------------------------
  ! Establish bookend 1: The code sets up the filenames required to 
  ! read both instantaneous and time average fields.  
  !-----------------------------------------------------------------
  if ( gdasT1534_struc(n)%findtime1 == 1 ) then  !get new time1 from the past
     write(LIS_logunit,*) 'Getting new time1 data'

     ferror = 0
     order = 1
     try = 0
     ts1 = -24*60*60

     do
        if ( ferror /= 0 ) exit
        try = try + 1
        status = 0 

        call create_gdasT1534filename(name, &
             gdasT1534_struc(n)%gdasT1534dir, yr1, mo1, da1, hr1)

        inquire(file=name,exist=file_exists)
        if ( .not. file_exists ) then
           status = 1
           write(LIS_logunit,*) 'ERR: GDAST1534 file1 ',trim(name), &
                'does not exist. (try = ', try, ')'
        endif
        if ( status == 0 ) then 
           write(LIS_logunit,*) 'Reading GDAST1534 file1 ',trim(name)

           call read_gdasT1534(order, n, findex, name, ferror, try)
           if ( ferror == 1 ) then  
              ! successfully retrieved forcing data
              gdasT1534_struc(n)%gdasT1534time1 = time1
           else  
              ! ferror still=0, so roll back one day and start again
              write(LIS_logunit,*) 'ERR: Error reading files. ', &
                   'Rolling back one day. (try = ', try, ')'
              call LIS_tick( dumbtime1, doy1, gmt1, &
                   yr1, mo1, da1, hr1, mn1, ss1, ts1 )
           end if
        else
           ! cannot find files, so roll back one day and start again
           write(LIS_logunit,*) 'MSG: Cannot find files. ', &
                'Rolling back one day. (try = ',try, ')'
           call LIS_tick( dumbtime1, doy1, gmt1, &
                yr1, mo1, da1, hr1, mn1, ss1, ts1 )
        endif
        if ( try > ndays ) then  
           ! data gap exceeds 10 days so stop
           write(LIS_logunit,*) 'ERROR: GDAST1534 data gap exceeds 10 days on file 1'
           call LIS_endrun
        endif
     end do
  end if

  !-----------------------------------------------------------------
  ! Establish gdasT1534time2
  !-----------------------------------------------------------------
  if ( movetime == 1 ) then  
     gdasT1534_struc(n)%gdasT1534time1 = gdasT1534_struc(n)%gdasT1534time2
     gdasT1534_struc(n)%findtime2 = 1  
     do f = 1, LIS_rc%met_nf(findex)
        do t = 1, LIS_rc%ngrid(n)
           gdasT1534_struc(n)%metdata1(f,t) = gdasT1534_struc(n)%metdata2(f,t)
        end do
     end do
  end if

  if ( gdasT1534_struc(n)%findtime2 == 1 ) then  
     write(LIS_logunit,*) 'Getting new time2 data'

     ferror = 0
     order = 2
     try = 0
     ts2 = -24*60*60

     do
        if ( ferror /= 0 ) exit
        try = try + 1
        status = 0

        call create_gdasT1534filename(name, &
             gdasT1534_struc(n)%gdasT1534dir, yr2, mo2, da2, hr2)

        inquire(file=name,exist=file_exists) 
        if ( .not. file_exists ) then
           status = 1
           write(LIS_logunit,*) 'ERR: GDAST1534 file2 ',trim(name), &
                'does not exist. (f00) (try = ', try, ')'
        endif
        if ( status == 0 ) then 
           write(LIS_logunit,*) 'Reading GDAST1534 file2  ',trim(name)
           call read_gdasT1534(order, n, findex, &
                name, ferror, try)
           if ( ferror == 1 ) then  
              ! successfully retrieved forcing data
              gdasT1534_struc(n)%gdasT1534time2 = time2
           else  
              ! ferror still=0, so roll back one day and start again
              write(LIS_logunit,*) 'ERR: Error reading files. ', &
                   'Rolling back one day. (try = ', try, ')'
              call LIS_tick( dumbtime2, doy2, gmt2, yr1, mo1, &
                   da1, hr1, mn1, ss1, ts2 )
           end if
        else
           ! cannot find files, so roll back one day and start again
           write(LIS_logunit,*) 'MSG: Cannot find files. ', &
                'Rolling back one day. (try = ',try, ')'
           call LIS_tick( dumbtime2, doy2, gmt2, yr1, mo1, &
                da1, hr1, mn1, ss1, ts2 )
        endif
        if ( try > ndays ) then  
           ! data gap exceeds 10 days so stop
           write(LIS_logunit,*) 'ERROR: GDAST1534 data gap exceeds 10 days on file 2'
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
           if ( (gdasT1534_struc(n)%metdata2(f,t) /= -9999.9) .and.  &
                (gdasT1534_struc(n)%metdata2(f,t) < 0)) then
              gdasT1534_struc(n)%metdata2(f,t) = (-1) * gdasT1534_struc(n)%metdata2(f,t)
           endif
           if ( (gdasT1534_struc(n)%metdata1(f,t) /= -9999.9) .and.  &
                (gdasT1534_struc(n)%metdata1(f,t) < 0)) then
              gdasT1534_struc(n)%metdata1(f,t) = (-1) * gdasT1534_struc(n)%metdata1(f,t)
           endif
        enddo
     endif
  enddo
  !  endif
end subroutine get_gdasT1534


