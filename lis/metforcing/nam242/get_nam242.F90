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
! !ROUTINE: get_nam242
!  \label{get_nam242}
!
! !REVISION HISTORY:
!     Sep 2012: NOHRSC/NOAA: Initial specification
!  12 Mar 2013: James Geiger: Committed into LIS
!  10 Apr 2013: James Geiger: Rewrote to use GRIBAPI
!
! !INTERFACE:
subroutine get_nam242(n, findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use LIS_timeMgrMod,    only : LIS_tick, LIS_get_nstep
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN
  use nam242_forcingMod, only : nam242_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, NAM forcing. 
!  The idea is to open either the 00, 03 and 06 forecast file associated with
!  the most recent NAM assimilation (available every 6 hours).
!  Precipitation rates and radiation fluxes will be taken from the F03 and 
!  F06 files, since averages are provided. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 3 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The strategy for missing data is to first look for 9hr forecasts and if not
!  go backwards up to 10 days to get forcing at the same time of day.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the NAM data times
!  \item[create\_nam242filename](\ref{create_nam242filename}) \newline
!    Puts together appropriate file name for 3 hour intervals
!  \item[create\_nam242f9\_filename](\ref{create_nam242f9_filename}) \newline
!    Puts together appropriate file name for 9 hour intervals
!  \item[read\_nam242](\ref{read_nam242}) \newline
!      Interpolates NAM data to LIS grid
!   \item[read\_nam242\_elev](\ref{read_nam242_elev}) \newline
!    reads the native elevation of the NAM data to be used
!    for topographic adjustments to the forcing 
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
  logical :: file_exists, file_exists1, file_exists2
  integer :: option
  real :: gridDesci(50)
  integer :: nstep
  logical :: F06flag
  integer :: status

  nam242_struc(n)%findtime1=0
  nam242_struc(n)%findtime2=0

  movetime=0

!-----------------------------------------------------------------
! Determine the correct number of forcing variables
!-----------------------------------------------------------------
  nstep = LIS_get_nstep(LIS_rc,n)
  nforce = nam242_struc(n)%nmif

  if ( nstep == 0 .or. LIS_rc%rstflag(n)== 1) then
    nam242_struc(n)%findtime1=1
    nam242_struc(n)%findtime2=1
    movetime=0        ! movetime is not properly set at time-step = 1
    LIS_rc%rstflag(n) = 0
  endif 

!-----------------------------------------------------------------
! Determine required NAM data times 
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
     nam242_struc(n)%namtime1 = time1
     nam242_struc(n)%namtime2 = time2
  endif
!-----------------------------------------------------------------
! Check to see if current time (timenow) has crossed past namtime2,
! requiring that both namtime2 be assigned to namtime1 and a new
! namtime2 be set 3 or 6 hours ahead of the current namtime2.
!-----------------------------------------------------------------
  if ( timenow > nam242_struc(n)%namtime2 ) then
     movetime  = 1
     nam242_struc(n)%findtime2 = 1
  end if

!-----------------------------------------------------------------
! Establish bookend 1: The code sets up the filenames required to 
! read both instantaneous and time average fields.  
!-----------------------------------------------------------------
  if ( nam242_struc(n)%findtime1 == 1 ) then  !get new time1 from the past
     write(LIS_logunit,*) 'Getting new time1 data'
     ferror = 0
     order = 1
     try = 0
     ts1 = -24*60*60
     status = 0 
     do
        if ( ferror /= 0 ) exit
        try = try+1
        call create_nam242filename( order, name00, name03, name06, F06flag, &
             nam242_struc(n)%namdir, yr1, mo1, da1, hr1 )
        write(LIS_logunit,*) 'Reading NAM file1 (I) ',trim(name00)
        write(LIS_logunit,*) 'Reading NAM file1 (A) ',trim(name03)
        inquire(file=name00,exist=file_exists1) 
        inquire(file=name03,exist=file_exists2)
!-----------------------------------------------------------------
! look for backup f9 files
!-----------------------------------------------------------------
        if((.not.file_exists1).or.(.not.file_exists2)) then 
           call create_nam242f9_filename( order, name00, name03,&
                nam242_struc(n)%namdir, yr1, mo1, da1, hr1,status )
           inquire(file=name00,exist=file_exists1) 
           inquire(file=name03,exist=file_exists2)
           if(file_exists1.and.file_exists2) then 
              status = 0 
           else
              status = 1
           endif
        endif
           
        if(F06flag) then 
           write(LIS_logunit,*) 'Reading NAM file1 (A)',trim(name06)
        endif

        if(status.eq.0) then
           call read_nam242(n, findex, order, name00, name03, name06, &
                            F06flag, ferror, try)
           if ( ferror == 1 ) then
!-----------------------------------------------------------------
! successfully retrieved forcing data
!-----------------------------------------------------------------
              nam242_struc(n)%namtime1 = time1
           else  
!-----------------------------------------------------------------
! ferror still=0, so roll back one day & start again
!-----------------------------------------------------------------
              call LIS_tick(dumbtime1, doy1, gmt1, yr1, mo1, da1, &
                            hr1, mn1, ss1, ts1)
           end if
        endif
        if ( try > ndays ) then  
!-----------------------------------------------------------------
! data gap exceeds 10 days so stop
!-----------------------------------------------------------------
           write(LIS_logunit,*) 'ERROR: NAM data gap exceeds 10 days on file 1'
           call LIS_endrun
        endif
     end do
  end if
!-----------------------------------------------------------------
! Establish namtime2
!-----------------------------------------------------------------
  if ( movetime == 1 ) then  
     nam242_struc(n)%namtime1 = nam242_struc(n)%namtime2
     nam242_struc(n)%metdata1 = nam242_struc(n)%metdata2
  end if
  
  if ( nam242_struc(n)%findtime2 == 1 ) then  
     write(LIS_logunit,*) 'Getting new time2 data'
  
     ferror = 0
     order = 2
     try = 0
     ts2 = -24*60*60
!-----------------------------------------------------------------
! determine the required forecast files. 
!-----------------------------------------------------------------
     do
        if ( ferror /= 0 ) exit
        try = try+1
        call create_nam242filename( order, name00, name03,name06, &
             F06flag, nam242_struc(n)%namdir, yr1, mo1, da1, hr1 )
        write(LIS_logunit,*) 'First Reading NAM file2 (I) ',trim(name00)
        write(LIS_logunit,*) 'First Reading NAM file2 (A) ',trim(name03)

        inquire(file=name00,exist=file_exists1) 
        inquire(file=name03,exist=file_exists2)
        if((.not.file_exists1).or.(.not.file_exists2)) then !look for backup f9 files
           call create_nam242f9_filename( order, name00, name03,&
                nam242_struc(n)%namdir, yr1, mo1, da1, hr1,status )
           inquire(file=name00,exist=file_exists1) 
           inquire(file=name03,exist=file_exists2)
           if(file_exists1.and.file_exists2) then 
              status = 0 
              write(LIS_logunit,*) 'F9 Reading NAM file2 (I) ',trim(name00)
              write(LIS_logunit,*) 'F9 Reading NAM file2 (A) ',trim(name03)
           else
              status = 1
           endif
        endif

        if(F06flag) then 
           write(LIS_logunit,*) 'Reading NAM file2 (A)',trim(name06)
        endif
        call read_nam242(n, findex, order, name00, name03, name06, &
                         F06flag, ferror, try)
        if ( ferror == 1 ) then  
!-----------------------------------------------------------------
! successfully retrieved forcing data
!-----------------------------------------------------------------
           nam242_struc(n)%namtime2 = time2
        else  
!-----------------------------------------------------------------
! ferror still=0, so roll back one day & start again
!-----------------------------------------------------------------
           call LIS_tick(dumbtime2, doy2, gmt2, yr1, mo1, da1, &
                         hr1, mn1, ss1, ts2 )
        end if
        if ( try > ndays ) then  
!-----------------------------------------------------------------
! data gap exceeds 10 days so stop
!-----------------------------------------------------------------
           write(LIS_logunit,*) 'ERROR: NAM data gap exceeds 10 days on file 2'
           call LIS_endrun
        end if
     end do
  end if
end subroutine get_nam242
