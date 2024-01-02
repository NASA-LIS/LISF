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
! !ROUTINE: get_cmorph
! \label{get_cmorph}
!
!
! !REVISION HISTORY:
! 17 Jul 2001: Jon Gottschalck; Initial code
! 10 Oct 2001: Jon Gottschalck; Modified to adjust convective precip
!               using a ratio of the model convective / total ratio
! 29 Dec 2003: Luis Goncalves; Added CMORPH global observed precip data sources
! 06 Jan 2005: Yudong Tian; Modified for LDTv4.2
!
! !INTERFACE:
subroutine get_cmorph(n, findex)
! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_masterproc
  use LDT_timeMgrMod, only : LDT_tick, LDT_get_nstep
  use cmorph_forcingMod, only :cmorph_struc
  use LDT_logMod,  only : LDT_logunit
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 30 minute CMORPH forcing. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 30min interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
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
!    determines the CMORPH data times
!  \item[cmorphfile](\ref{cmorphfile}) \newline
!    Puts together appropriate file filename for 6 hour intervals
!  \item[read\_cmorph](\ref{read_cmorph}) \newline
!      Interpolates CMORPH data to LDT grid
!  \end{description}
!EOP

  integer :: ferror_cmorph      ! Error flags for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy4, yr4, mo4, da4, hr4, mn4, ss4
  integer :: endtime_cmorph     ! 1=get a new file
  integer :: sectionofcmorph
  real*8  :: ctime,ftime_cmorph ! Current LDAS time and end boundary times for precip data sources 
  real*8  :: datatime, breaktime, ffilenametime     
! Times used in HUFFMAN to determine data and file filename boundaries (see below)
  integer :: order
  real    :: gmt1,gmt4,ts1,ts4
  character(len=LDT_CONST_PATH_LEN) :: filename ! Filefilename variables for precip data sources

!=== End Variable Definition =======================

 endtime_cmorph = 0

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
! CMORPH product end time
!------------------------------------------------------------------------
  yr4 = LDT_rc%yr  !end accumulation time data
  mo4 = LDT_rc%mo
  da4 = LDT_rc%da
!!!  hr4 = 0.5*(LDT_rc%t%hr/0.5)
  hr4 = LDT_rc%hr
  mn4 = 30*(LDT_rc%mn/30)
!!!  mn4 = 0
  ss4 = 0
  ts4 = 0.5*60*60
!!!  ts4 = 1*60*60
  call LDT_tick( ftime_cmorph, doy4, gmt4, yr4, mo4, da4, hr4, mn4, ss4, ts4 )
  breaktime = ftime_cmorph - ctime
  datatime  = ftime_cmorph
  ffilenametime = ftime_cmorph

!------------------------------------------------------------------------
! Ensure that data is found during first time step
!------------------------------------------------------------------------
  endtime_cmorph = 0
  if ( LDT_get_nstep(LDT_rc,n).eq. 1 .or.LDT_rc%rstflag(n) .eq. 1) then
     endtime_cmorph = 1
     LDT_rc%rstflag(n) = 0
  endif

!------------------------------------------------------------------------
! Check for and get CMORPH precipitation data
!------------------------------------------------------------------------
   if ( ctime > cmorph_struc(n)%cmorphtime ) then
      endtime_cmorph = 1
   End If
   
    if ( LDT_rc%mn >= 0 .AND. LDT_rc%mn < 30 ) sectionofcmorph = 1
    if ( LDT_rc%mn >= 30 ) sectionofcmorph = 2
!    if(LDT_masterproc)  write(*, *) "get_cmorph 0: LDT_rc%mn=", LDT_rc%mn, &
!        " sectionofcmorph=", sectionofcmorph, "  endtime_cmorph=", endtime_cmorph
   if ( endtime_cmorph == 1 ) then  !get new time2 data first 1/2 hour
      ferror_cmorph = 0
!       print *,'dir=',cmorph_struc(n)%cmorphdir
      call cmorphfile( filename, cmorph_struc(n)%cmorphdir, yr1, mo1, da1, hr1 )
     if(LDT_masterproc) then
      if (sectionofcmorph .EQ. 1) then
        write(LDT_logunit,*) 'Getting new CMORPH precip data first 30 minutes'
        write(LDT_logunit,*) trim(filename)
      else
        write(LDT_logunit,*) 'Getting new CMORPH precip data second 30 minutes'
        write(LDT_logunit,*) trim(filename)
      endif
     end if
     order = 2
     call read_cmorph( n, filename, findex, order, ferror_cmorph, sectionofcmorph )

   endif  !need new time2

end subroutine get_cmorph

!BOP
! !ROUTINE: cmorphfile
! \label{cmorphfile}
!
!
! !INTERFACE:
subroutine cmorphfile( filename, cmorphdir, yr, mo, da, hr)

  implicit none
! !ARGUMENTS: 
  character(len=*) :: filename
  character(len=*) :: cmorphdir
  integer          :: yr, mo, da, hr
! !DESCRIPTION:
!   This subroutine puts together CMORPH file filename for 
!   30min file intervals.
! 
!  The arguments are:
!  \begin{description}
!  \item[cmorphdir]
!    Name of the CMORPH directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[filename]
!   filename of the timestamped CMORPH file
!  \end{description}
!
!EOP

!date -d "2008/06/01" +%s
!1212292800

  integer, parameter :: T2008060100 = 1212292800
  character(len=120) :: temp
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  integer tout(9), fmktime, it, ih, irec

  character*100 :: fbase, ftimedir, fstem 
  character*4 :: cyr
  character*2 :: cmo, cda, chr 

!=== End Variable Definition ===============
!=== formats for file filename segments
!------------------------------------------------------------------------
! Make variables for the time used to create the file
! We don't want these variables being passed out
!------------------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 1*(hr/1)  !hour needs to be a multiple of 1 hour
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.
 
  write(cyr, '(I4.4)') uyr
  write(cmo, '(I2.2)') umo 
  write(cda, '(I2.2)') uda 
  write(chr, '(I2.2)') uhr 

  tout=0 
  tout(6)= uyr-1900
  tout(5) = umo - 1
  tout(4) = uda 
  tout(3) = uhr
  tout(2) = umn 
  tout(1) = uss 

  fstem = '/advt-8km-intrp-prim-sat-spat-2lag-2.5+5dovlp8kmIR-'
  if ( fmktime(tout) .GE. T2008060100)  &
    fstem = '/advt-8km-interp-prim-sat-spat-2lag-2.5+5dovlp8kmIR-'
!       print *,'In cmorfile: cmorphdir=',cmorphdir
 
  filename = trim(cmorphdir) // "/" // cyr // cmo // &
                trim(fstem) // cyr // cmo // cda // chr 

  return

end subroutine cmorphfile
