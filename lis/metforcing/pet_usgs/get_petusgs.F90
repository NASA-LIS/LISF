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
! !ROUTINE: get_petusgs
! \label{get_petusgs}
!
!
! !REVISION HISTORY:
! 17 Jul 2001: Jon Gottschalck; Initial code
! 06 Jan 2005: Yudong Tian; Modified for LISv4.2
! 12 Mar 2012: K. Arsenault;  Updated for USGS PET data
! 25 Oct 2013: K. Arsenault;  Added PET USGS to LIS7
!
! !INTERFACE:
subroutine get_petusgs(n, findex)

! !USES:
  use LIS_coreMod,    only : LIS_rc
  use LIS_timeMgrMod, only : LIS_tick, LIS_get_nstep
  use LIS_logMod,     only : LIS_logunit
  use petusgs_forcingMod, only : petusgs_struc
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates daily USGS PET forcing. 
!  At the beginning of a simulation, the code reads the 
!  most recent ``past'' data file (current 1-day), and nearest 
!  ``future'' data file (next 1-day). These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!.
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!     determines the USGS PET data times
!  \item[petusgsfile](\ref{petusgsfile}) \newline
!     Puts together appropriate file name for 6 hour intervals
!  \item[read\_petusgs](\ref{read_petusgs}) \newline
!     Interpolates USGS PET data to LIS grid
!  \end{description}
!EOP

  integer :: ferror_pet                  ! Error flags for PETdata sources
  integer :: endtime_petusgs             ! 1=get a new file
  real*8  :: ctime, ftime_petusgs        ! Current time and end boundary times for PET data sources 
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1   ! Time parameters for current LIS time
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2   ! Time parameters for file boundary end time
  real    :: ts1, ts2
  real    :: gmt1, gmt2
  integer :: kk                          ! Forecast index

  character(len=LIS_CONST_PATH_LEN) :: filename              ! Filename variables for PET data sources

!=== End Variable Definition =======================

!------------------------------------------------------------------------
! Determine required observed USGS PET data times 
! (current, accumulation end time)
! - Model current time
!------------------------------------------------------------------------
  yr1 = LIS_rc%yr    !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick( ctime, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

!------------------------------------------------------------------------
! USGS PET end time
!------------------------------------------------------------------------
  yr2 = LIS_rc%yr    !end accumulation time data
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = 0
  mn2 = 0
  ss2 = 0
  ts2 = ((23*60)+59)*60  ! 23 hrs + 59 mins --> Seconds
  call LIS_tick( ftime_petusgs, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

!------------------------------------------------------------------------
! Ensure that data is found during the first time step
!------------------------------------------------------------------------

   endtime_petusgs = 0

!   print *, "ctime:", yr1, mo1, da1, hr1, mn1, LIS_get_nstep(LIS_rc,n)
!   print *, "ftime:", yr2, mo2, da2, hr2, mn2, LIS_get_nstep(LIS_rc,n)

   if( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1 ) then
      endtime_petusgs = 1
      LIS_rc%rstflag(n) = 0
   endif

!------------------------------------------------------------------------
! Check for and get USGS PET data file
!------------------------------------------------------------------------
    if( ctime > petusgs_struc(n)%pettime )  endtime_petusgs = 1
   
    if( endtime_petusgs == 1 ) then  ! get new time2 daily data file

      do kk= petusgs_struc(n)%st_iterid, petusgs_struc(n)%en_iterid

        ! Form USGS PET filename:
        call petusgsfile( n, kk, findex, filename, petusgs_struc(n)%pettype,&
                          petusgs_struc(n)%petdir, yr1, mo1, da1 )

        write(LIS_logunit,*) "Getting new USGS PET data (new 24 hour period): "
        write(LIS_logunit,*)  trim(filename)

        ! Read in and spatially interpolate latest USGS PET file obtained:
        ferror_pet = 0
        call read_petusgs( n, kk, findex, filename, ferror_pet )

        petusgs_struc(n)%pettime = ftime_petusgs

      end do  ! end forecast loop
    endif     ! need new time2

end subroutine get_petusgs

! --------------------------

!BOP
! !ROUTINE: petusgsfile
! \label{petusgsfile}
!
!
! !INTERFACE:
subroutine petusgsfile( n, kk, findex, filename, filetype, &
                        petusgsdir, yr, mo, da )

  use LIS_coreMod
  use LIS_forecastMod

  implicit none
! !ARGUMENTS: 
  integer         , intent(in)  :: n
  integer         , intent(in)  :: kk
  integer         , intent(in)  :: findex
  character(len=*), intent(out) :: filename
  character(len=*), intent(in)  :: petusgsdir
  character(len=*), intent(in)  :: filetype
  integer         , intent(in)  :: yr, mo, da

! !DESCRIPTION:
!   This subroutine puts together daily USGS PET filenames.
! 
!  The arguments are:
!  \begin{description}
!  \item[petusgsdir]
!    Name of the USGS PET directory
!  \item[yr]
!    year 
!  \item[mo]
!    month
!  \item[da]
!    day of month
!  \item[filename]
!    name of the timestamped USGS PET file
!  \end{description}
!
!EOP

   character(4) :: fyr4
   character(2) :: fyr2
   character(2) :: fmo, fda, fhr

!=== end variable definition =============================================

   if(LIS_rc%forecastMode.eq.0) then !hindcast run

     write(unit=fyr4,fmt='(i4.4)') yr
     write(unit=fyr2,fmt='(i2.2)') yr-2000
     write(unit=fmo, fmt='(i2.2)') mo
     write(unit=fda, fmt='(i2.2)') da

     ! Do not need to specify file extension, based on fbil_module routines ...
     if( filetype == "current" ) then 
       filename = trim(petusgsdir)//"/"//fyr4//"/et"//fyr2//fmo//fda
     elseif( filetype == "climatology" ) then
       filename = trim(petusgsdir)//"/etclim"//fmo//fda
     endif

   ! Forecast mode (e.g., ESP):
   else
     call LIS_sample_forecastDate(n,kk,findex,yr,mo,da)

     write(unit=fyr4,fmt='(i4.4)') yr
     write(unit=fyr2,fmt='(i2.2)') yr-2000
     write(unit=fmo, fmt='(i2.2)') mo
     write(unit=fda, fmt='(i2.2)') da

     ! Do not need to specify file extension, based on fbil_module routines ...
     if( filetype == "current" ) then
       filename = trim(petusgsdir)//"/"//fyr4//"/et"//fyr2//fmo//fda
     elseif( filetype == "climatology" ) then
       filename = trim(petusgsdir)//"/etclim"//fmo//fda
     endif
   endif

end subroutine petusgsfile
