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
! !ROUTINE: get_cmap
! \label{get_cmap}
!
! !REVISION HISTORY:
! 17 Jul 2001: Jon Gottschalck; Initial code
! 10 Oct 2001: Jon Gottschalck; Modified to adjust convective precip
!               using a ratio of the model convective / total ratio
! 16 Feb 2007: Chuck Alonge; Changed file name creation to use internal 
!               files instead of external file (avoids failure in parallel runs)
! 02 Dec 2014: KR Arsenault: Added new grid change update (~2012)
!
! !INTERFACE:
subroutine get_cmap(n,findex)
! !USES:
  use LIS_coreMod,     only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,  only : LIS_tick, LIS_get_nstep
  use LIS_logMod,      only : LIS_logunit, LIS_endrun
  use cmap_forcingMod, only : cmap_struc
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 6-hrly, CMAP forcing. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 6 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!
!   upto 2000/1/24          :   T126 (384x190)  grid
!   2001/1/24  - 2002/10/29 :   T170 (512x256)  grid
!   2002/10/29 - 2005/5/31  :   T254 (768x384)  grid
!   2005/5/31  - 2012/9/30  :   T382 (1152x576) grid ~~ CMAP only
!   2012/10/01 onwards      :   T574 (1760x880) grid ~~ CMAP only
!  Original GDAS grid change:
!   2010/7/28  onwards      :   T574 (1760x880) grid ~~ Original GDAS
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
!    determines the CMAP data times
!  \item[cmapfile](\ref{cmapfile}) \newline
!    Puts together appropriate file name for 6 hour intervals
!  \item[read\_cmap](\ref{read_cmap}) \newline
!    Interpolates CMAP data to LIS grid
!   \item[cmap\_reset\_interp\_input](\ref{cmap_reset_interp_input} \newline
!    resets the neighbours and weights arrays upon a grid change
!  \end{description}
!EOP
   
!==== Local Variables=======================
  integer :: ferror_cmap   ! Error flags for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy5, yr5, mo5, da5, hr5, mn5, ss5
  real*8  :: ctime,ftime_cmap    ! Current LDAS time and end boundary times for precip data sources 
  integer :: order
  real    :: gmt1,gmt5,ts1,ts5   ! GMT times for current LDAS time and end boundary times for precip data sources
  real    :: gridDesci(50)
  character(len=LIS_CONST_PATH_LEN) :: filename ! Filename variables for precip data sources
!=== End Variable Definition =======================

!------------------------------------------------------------------------
! Determine required observed precip data times 
! (current, accumulation end time)
! Model current time
!------------------------------------------------------------------------
  yr1 = LIS_rc%yr  !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick( ctime, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )   
!------------------------------------------------------------------------ 
! CMAP product end time
!------------------------------------------------------------------------
  yr5 = LIS_rc%yr  !end accumulation time data
  mo5 = LIS_rc%mo
  da5 = LIS_rc%da
  hr5 = 6*(LIS_rc%hr/6)
  mn5 = 0
  ss5 = 0
  ts5 = 6*60*60
  call LIS_tick( ftime_cmap, doy5, gmt5, yr5, mo5, da5, hr5, mn5, ss5, ts5 )

!------------------------------------------------------------------------
! Reinitialize the weights and neighbors for each required grid change
!------------------------------------------------------------------------

! 1991-2000 T126 Grid:
  if ( ctime > cmap_struc(n)%griduptime1 .and. &
       ctime <= cmap_struc(n)%griduptime2 .and. &
       cmap_struc(n)%gridchange1 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_cmap -- changing cmap grid to 1991 -- 2000"
     write(LIS_logunit,*) "** "

     cmap_struc(n)%ncold = 720
     cmap_struc(n)%nrold = 360
     gridDesci = 0
     gridDesci(1) = 0
     gridDesci(2) = 720
     gridDesci(3) = 360
     gridDesci(4) = -89.750
     gridDesci(5) =   0.250
     gridDesci(6) = 128
     gridDesci(7) =  89.750
     gridDesci(8) = 359.750
     gridDesci(9) =   0.500
     gridDesci(10)=   0.500
     gridDesci(20)= 0
     cmap_struc(n)%mi = gridDesci(2)*gridDesci(3)

     call cmap_reset_interp_input(n, findex, gridDesci)
     cmap_struc(n)%gridchange1 = .false.

! 2000-2002 T170 Grid:
  elseif ( ctime > cmap_struc(n)%griduptime2 .and. &
           ctime <= cmap_struc(n)%griduptime3 .and. &
           cmap_struc(n)%gridchange2 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_cmap -- changing cmap grid to 2000 -- 2002"
     write(LIS_logunit,*) "** "

     cmap_struc(n)%ncold = 512
     cmap_struc(n)%nrold = 256
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 512
     gridDesci(3) = 256
     gridDesci(4) = 89.463
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.463
     gridDesci(8) = -0.703
     gridDesci(9) = 0.703
     gridDesci(10) = 128
     gridDesci(20) = 0
     cmap_struc(n)%mi = gridDesci(2)*gridDesci(3)

     call cmap_reset_interp_input(n, findex, gridDesci)
     cmap_struc(n)%gridchange2 = .false.

! 2002-2005 T254 Grid:
  elseif ( ctime > cmap_struc(n)%griduptime3 .and. &
           ctime <= cmap_struc(n)%griduptime4 .and. &
           cmap_struc(n)%gridchange3 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_cmap -- changing cmap grid to 2002 -- 2005"
     write(LIS_logunit,*) "** "

     cmap_struc(n)%ncold = 768
     cmap_struc(n)%nrold = 384
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 768
     gridDesci(3) = 384
     gridDesci(4) = 89.462
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.462
     gridDesci(8) = -0.469
     gridDesci(9) = 0.469
     gridDesci(10) = 192
     gridDesci(20) = 0
     cmap_struc(n)%mi = gridDesci(2)*gridDesci(3)

     call cmap_reset_interp_input(n, findex, gridDesci)
     cmap_struc(n)%gridchange3 = .false.

! 2005-2012 T382 Grid:
  elseif ( ctime > cmap_struc(n)%griduptime4 .and. &
           ctime <= cmap_struc(n)%griduptime5 .and. &
           cmap_struc(n)%gridchange4 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_cmap -- changing cmap grid to 2005 -- 2012 "
     write(LIS_logunit,*) "** "

     cmap_struc(n)%ncold = 1152
     cmap_struc(n)%nrold = 576
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 1152
     gridDesci(3) = 576
     gridDesci(4) = 89.761
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.761
     gridDesci(8) = -0.313
     gridDesci(9) = 0.313
     gridDesci(10) = 288
     gridDesci(20) = 0
     cmap_struc(n)%mi = gridDesci(2)*gridDesci(3)

     call cmap_reset_interp_input(n, findex, gridDesci)
     cmap_struc(n)%gridchange4 = .false.

! 2012-- T574 Grid:
  elseif ( ctime > cmap_struc(n)%griduptime5 .and. &
           cmap_struc(n)%gridchange5 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_cmap -- changing cmap grid to 2012 -- "
     write(LIS_logunit,*) "** "

     cmap_struc(n)%ncold = 1760
     cmap_struc(n)%nrold = 880
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 1760
     gridDesci(3) = 880
     gridDesci(4) = 89.844
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = -89.844
     gridDesci(8) = -0.205
     gridDesci(9) = 0.205
     gridDesci(10) = 440
     gridDesci(20) = 0.0
     cmap_struc(n)%mi = gridDesci(2)*gridDesci(3)

     call cmap_reset_interp_input(n, findex, gridDesci)
     cmap_struc(n)%gridchange5 = .false.

  endif

!------------------------------------------------------------------------
! Ensure that data is found during first time step
!------------------------------------------------------------------------
  if( LIS_get_nstep(LIS_rc,n) == 1 .or. &
       LIS_rc%rstflag(n) == 1) then 
     LIS_rc%rstflag(n) = 0
  endif

!------------------------------------------------------------------------
! Check for and get CMAP CPC Precipitation data
!------------------------------------------------------------------------
   filename=""
   ferror_cmap = 0
   order = 2

 ! LIS timestep < CMAP time interval (6hr)
   if( LIS_rc%ts < cmap_struc(n)%ts ) then
     if( ctime > cmap_struc(n)%cmaptime ) then
      ! Put together CMAP filename:
        call cmapfile( filename, cmap_struc(n)%cmapdir, yr5, mo5, da5, hr5 )
        write(LIS_logunit,*) "Getting new CMAP CPC precip data: ",trim(filename)
        call read_cmap( n, filename, findex, order, ferror_cmap, hr5 )
        cmap_struc(n)%cmaptime = ftime_cmap
     endif  !need new time2

 ! LIS timestep >= CMAP time interval (6hr)
   elseif( LIS_rc%ts >= cmap_struc(n)%ts ) then
    ! Put together CMAP filename:
      call cmapfile( filename, cmap_struc(n)%cmapdir, yr1, mo1, da1, hr1 )
      write(LIS_logunit,*) "Getting new CMAP CPC precip data: ",trim(filename)
      call read_cmap( n, filename, findex, order, ferror_cmap, hr1 )
   endif


end subroutine get_cmap

!BOP
! !ROUTINE: cmapfile
! \label{cmapfile}
!
! !INTERFACE:
subroutine cmapfile( filename, cmapdir, yr, mo, da, hr)
  implicit none
! !ARGUMENTS: 
  character(len=*)   :: filename
  character(len=*)   :: cmapdir
  integer            :: yr, mo, da, hr
! !DESCRIPTION:
!   This subroutine puts together CMAP file name for 
!   6 hour file intervals.
! 
!  The arguments are:
!  \begin{description}
!  \item[cmapdir]
!    Name of the CMAP directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[filename]
!   file name of the timestamped CMAP file
!  \end{description}
!
!EOP

  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  character(len=6)  :: fdir
  character(len=10) :: ftime
  character(len=10), parameter :: fprefix = 'cmap_gdas_'
  character(len=4),  parameter :: fext = '.grb'

!=== End Variable Definition ===============

!------------------------------------------------------------------------
! Make variables for the time used to create the file
! We don't want these variables being passed out
!------------------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 6*(hr/6)  !hour needs to be a multiple of 6 hours
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.

  filename = ''

  write(UNIT=fdir,  fmt='(i4, i2.2)')  uyr, umo
  write(UNIT=ftime, fmt='(i4, i2.2, i2.2, i2.2)') uyr, umo, uda, uhr

  filename = trim(cmapdir) // '/' // fdir // '/' // fprefix // ftime // fext

  return
end subroutine cmapfile

!BOP
! !ROUTINE: cmapfile_old
!
! !DESCRIPTION: This subroutine puts together CMAP file name for
!               6 hour file intervals.
!
! !INTERFACE:
subroutine cmapfile_old( filename, cmapdir, yr, mo, da, hr )
!EOP
  implicit none

!==== Local Variables=======================

  character(len=*)   :: filename
  character(len=*)   :: cmapdir
  character(len=100) :: temp
  integer :: yr, mo, da, hr
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  character*1 :: fbase(80), fdir(6), ftime(10), fsubs(11), fsubs2(4)

!=== End Variable Definition ===============

!=== formats for filename segments
!------------------------------------------------------------------------
! Make variables for the time used to create the file
! We don't want these variables being passed out
!------------------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 6*(hr/6)  !hour needs to be a multiple of 6 hours
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.

  filename=''  

  write(UNIT=temp, fmt='(a40)') cmapdir
  read(UNIT=temp, fmt='(80a1)') (fbase(i), i=1,80)

  write(UNIT=temp, fmt='(a1,i4,a1)') '/', uyr, '/'
  read(UNIT=temp, fmt='(6a1)') fdir
  do i = 1, 6
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do

  write(UNIT=temp, fmt='(a11)') '_drean_smth'
  read (UNIT=temp, fmt='(11a1)') (fsubs(i), i=1,11)

  write(UNIT=temp, fmt='(i4,i2,i2,i2)') uyr, umo, uda, uhr
  read(UNIT=temp, fmt='(10a1)') ftime
  do i = 1, 10
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do

  write(UNIT=temp, fmt='(a4)') '.grb'
  read (UNIT=temp, fmt='(4a1)') (fsubs2(i), i=1,4)
  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(UNIT=temp, fmt='(80a1)') (fbase(i), i=1,c), (fdir(i), i=1,6),  &
                       (ftime(i), i=1,10),(fsubs(i), i=1,11)

  read(UNIT=temp, fmt='(a80)') filename

  return
end subroutine cmapfile_old
