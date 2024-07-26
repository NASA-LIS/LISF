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
!  10 Oct 2014: KR Arsenault: Added to LDT
!
! !INTERFACE:
subroutine get_cmap(n,findex)

! !USES:
  use LDT_coreMod, only : LDT_rc
  use LDT_logMod,  only : LDT_logunit
  use LDT_timeMgrMod, only : LDT_tick, LDT_get_nstep
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use cmap_forcingMod,only : cmap_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 6-hrly, CMAP forcing. 
!  At the beginning of a simulation, the code reads
!  the most recent past data (recent 6 hour interval), and
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
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    Determines the CMAP data times
!  \item[cmapfile](\ref{cmapfile}) \newline
!    Puts together appropriate file name for 6 hour intervals
!  \item[read\_cmap](\ref{read_cmap}) \newline
!    Interpolates CMAP data to LDT grid
!  \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    Computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
   
!==== Local Variables=======================
  integer :: ferror_cmap         ! Error flags for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy5, yr5, mo5, da5, hr5, mn5, ss5
  real*8  :: ctime, ftime_cmap   ! Current LDAS time and end boundary times for precip data sources 
  integer :: order
  real    :: gmt1,gmt5,ts1,ts5   ! GMT times for current LDAS time and end boundary times for precip data sources
  real    :: gridDesci(20)
  character(len=LDT_CONST_PATH_LEN) :: filename       ! Filename variables for precip data sources

!=== End Variable Definition =======================

!------------------------------------------------------------------------
! Determine required observed precip data times 
! (current, accumulation end time)
! Model current time
!------------------------------------------------------------------------
  yr1 = LDT_rc%yr  !current time
  mo1 = LDT_rc%mo
  da1 = LDT_rc%da
  hr1 = LDT_rc%hr
  mn1 = LDT_rc%mn
  ss1 = 0
  ts1 = 0
  call LDT_tick( ctime, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )   
!------------------------------------------------------------------------ 
! CMAP product end time
!------------------------------------------------------------------------
  yr5 = LDT_rc%yr  !end accumulation time data
  mo5 = LDT_rc%mo
  da5 = LDT_rc%da
  hr5 = 6*(LDT_rc%hr/6)
  mn5 = 0
  ss5 = 0
  ts5 = 6*60*60
  call LDT_tick( ftime_cmap, doy5, gmt5, yr5, mo5, da5, hr5, mn5, ss5, ts5 )

!------------------------------------------------------------------------
! Reinitialize the weights and neighbors for each required grid change
!------------------------------------------------------------------------

! 1991-2000 T126 Grid:
  if( ctime > cmap_struc(n)%griduptime1 .and. & 
      ctime <= cmap_struc(n)%griduptime2 .and. & 
      cmap_struc(n)%gridchange1 ) then 

     write(LDT_logunit,*) "** "
     write(LDT_logunit,*) "MSG: get_cmap -- changing cmap grid to 1991 -- 2000"
     write(LDT_logunit,*) "** "

     call update_cmapgrid( 1, cmap_struc(n)%nc, &
                           cmap_struc(n)%nr, gridDesci )

     call conserv_interp_input(n, gridDesci,&
          cmap_struc(n)%n112, cmap_struc(n)%n122,&
          cmap_struc(n)%n212, cmap_struc(n)%n222,&
          cmap_struc(n)%w112, cmap_struc(n)%w122,&
          cmap_struc(n)%w212, cmap_struc(n)%w222)

     cmap_struc(n)%gridchange1 = .false.

! 2000-2002 T170 Grid:
  elseif( ctime > cmap_struc(n)%griduptime2 .and. & 
          ctime <= cmap_struc(n)%griduptime3 .and. & 
          cmap_struc(n)%gridchange2 ) then 

     write(LDT_logunit,*) "** "
     write(LDT_logunit,*) "MSG: get_cmap -- changing cmap grid to 2000 -- 2002"
     write(LDT_logunit,*) "** "

     call update_cmapgrid( 2, cmap_struc(n)%nc, cmap_struc(n)%nr, &
                           gridDesci(:) )

     call conserv_interp_input(n, gridDesci,&
          cmap_struc(n)%n112, cmap_struc(n)%n122,&
          cmap_struc(n)%n212, cmap_struc(n)%n222,&
          cmap_struc(n)%w112, cmap_struc(n)%w122,&
          cmap_struc(n)%w212, cmap_struc(n)%w222)

     cmap_struc(n)%gridchange2 = .false.

! 2002-2005 T254 Grid:
  elseif( ctime > cmap_struc(n)%griduptime3 .and. & 
          ctime <= cmap_struc(n)%griduptime4 .and. & 
          cmap_struc(n)%gridchange3 ) then 

     write(LDT_logunit,*) "** "
     write(LDT_logunit,*) "MSG: get_cmap -- changing cmap grid to 2002 -- 2005"
     write(LDT_logunit,*) "** "

     call update_cmapgrid( 3, cmap_struc(n)%nc, &
                 cmap_struc(n)%nr, gridDesci )

     call conserv_interp_input(n, gridDesci,&
          cmap_struc(n)%n112, cmap_struc(n)%n122,&
          cmap_struc(n)%n212, cmap_struc(n)%n222,&
          cmap_struc(n)%w112, cmap_struc(n)%w122,&
          cmap_struc(n)%w212, cmap_struc(n)%w222)

     cmap_struc(n)%gridchange3 = .false.

! 2005-2012 T382 Grid:
  elseif( ctime > cmap_struc(n)%griduptime4 .and. & 
          ctime <= cmap_struc(n)%griduptime5 .and. & 
          cmap_struc(n)%gridchange4 ) then 

     write(LDT_logunit,*) "** "
     write(LDT_logunit,*) "MSG: get_cmap -- changing cmap grid to 2005 -- 2012 "
     write(LDT_logunit,*) "** "

     call update_cmapgrid( 4, cmap_struc(n)%nc, cmap_struc(n)%nr, &
                           gridDesci(:) )

     call conserv_interp_input(n, gridDesci,&
          cmap_struc(n)%n112, cmap_struc(n)%n122,&
          cmap_struc(n)%n212, cmap_struc(n)%n222,&
          cmap_struc(n)%w112, cmap_struc(n)%w122,&
          cmap_struc(n)%w212, cmap_struc(n)%w222)

     cmap_struc(n)%gridchange4 = .false.

! 2012-- T574 Grid:
  elseif( ctime > cmap_struc(n)%griduptime5 .and. &
          cmap_struc(n)%gridchange5 ) then

     write(LDT_logunit,*) "** "
     write(LDT_logunit,*) "MSG: get_cmap -- changing cmap grid to 2012 -- "
     write(LDT_logunit,*) "** "

     call update_cmapgrid( 5, cmap_struc(n)%nc, cmap_struc(n)%nr, &
                           gridDesci(:) )

     call conserv_interp_input(n, gridDesci,&
          cmap_struc(n)%n112, cmap_struc(n)%n122,&
          cmap_struc(n)%n212, cmap_struc(n)%n222,&
          cmap_struc(n)%w112, cmap_struc(n)%w122,&
          cmap_struc(n)%w212, cmap_struc(n)%w222)

     cmap_struc(n)%gridchange5 = .false.

  endif

!------------------------------------------------------------------------
! Ensure that data is found during first time step
!------------------------------------------------------------------------
  if( LDT_get_nstep(LDT_rc,n) == 1 .or. &
       LDT_rc%rstflag(n) == 1 ) then 
     LDT_rc%rstflag(n) = 0
  endif

!------------------------------------------------------------------------
! Check for and get CMAP CPC Precipitation data
!------------------------------------------------------------------------
   filename = ""
   ferror_cmap = 0
   order = 2

 ! LDT timestep < CMAP time interval (6hr)
   if( LDT_rc%ts < cmap_struc(n)%ts ) then
     if( ctime > cmap_struc(n)%cmaptime ) then
      ! Put together CMAP filename:
        call cmapfile( filename, cmap_struc(n)%cmapdir, yr5, mo5, da5, hr5 )
        write(LDT_logunit,*) "Getting new CMAP CPC precip data: ",trim(filename)
        call read_cmap( n, filename, findex, order, ferror_cmap, hr5 )
        cmap_struc(n)%cmaptime = ftime_cmap
     endif

 ! LDT timestep >= CMAP time interval (6hr)
   elseif( LDT_rc%ts >= cmap_struc(n)%ts ) then
    ! Put together CMAP filename:
      call cmapfile( filename, cmap_struc(n)%cmapdir, yr1, mo1, da1, hr1 )
      write(LDT_logunit,*) "Getting new CMAP CPC precip data: ",trim(filename)
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
  character(len=*), intent(out) :: filename
  character(len=*), intent(in)  :: cmapdir
  integer, intent(in)           :: yr, mo, da, hr
!
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
!    month
!  \item[da]
!    day of month
!  \item[hr]
!    hour of day
!  \item[filename]
!    name of the timestamped CMAP file
!  \end{description}
!
!EOP
  character*4  :: fyr
  character*2  :: fmo, fda, fhr

!=== end variable definition =============================================

  write(unit=fyr, fmt='(i4.4)')  yr
  write(unit=fmo, fmt='(i2.2)')  mo
  write(unit=fda, fmt='(i2.2)')  da
  write(unit=fhr, fmt='(i2.2)')  hr

  filename = trim(cmapdir)//"/"//fyr//fmo//"/cmap_gdas_"&
           //fyr//fmo//fda//fhr//".grb"

end subroutine cmapfile


!BOP
! !ROUTINE: update_cmapgrid
! \label{update_gmapgrid}
!
! !INTERFACE:
subroutine update_cmapgrid( gridchange, nc, nr, gridDesci )

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: gridchange
  integer, intent(out) :: nc 
  integer, intent(out) :: nr
  real,    intent(out) :: gridDesci(20)

! ___________________________________________

  select case( gridchange )

   case( 1 )  ! changing cmap grid to 1991 -- 2000

     nc = 720
     nr = 360
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

   case( 2 )  ! changing cmap grid to 2000 -- 2002

     nc = 512
     nr = 256
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

   case( 3 )  ! changing cmap grid to 2002 -- 2005

     nc = 768
     nr = 384
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

   case( 4 )  ! changing cmap grid to 2005 --

     nc = 1152
     nr = 576
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

   case( 5 )  ! changing cmap grid to 2005 --

     nc = 1760
     nr = 880
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

   case default

     print *, " No gridchange by that entry exists for CMAP ... Stopping"
!     call LDT_endrun
     stop

  end select

end subroutine update_cmapgrid
