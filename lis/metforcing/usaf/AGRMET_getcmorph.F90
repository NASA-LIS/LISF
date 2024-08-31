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
! !ROUTINE: AGRMET_getcmorph
! \label{AGRMET_getcmorph}
!
!
! !REVISION HISTORY:
! 17 Jul 2001: Jon Gottschalck; Initial code
! 10 Oct 2001: Jon Gottschalck; Modified to adjust convective precip
!               using a ratio of the model convective / total ratio
! 29 Dec 2003: Luis Goncalves; Added CMORPH global observed precip data sources
! 06 Jan 2005: Yudong Tian; Modified for LISv4.2
! 05 May 2013; Moved to AGRMET from suppforcing...Ryan Ruhge/16WS/WXE/SEMS
!
! !INTERFACE:
subroutine AGRMET_getcmorph(n, cmorphdata, j3hr, quad9r)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only : LIS_julhr_date
  use agrmet_forcingMod, only :agrmet_struc
  use LIS_logMod, only : LIS_logunit

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  real                :: quad9r
  real                :: cmorphdata(LIS_rc%lnc(n), LIS_rc%lnr(n))
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 30 minute CMORPH forcing. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 30min interval), and
!  the nearest future data. These two datasets are used to 
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
!  \item[cmorfile\_agrmet](\ref{cmorfile_agrmet}) \newline
!    Puts together appropriate file name for 6 hour intervals
!  \item[AGRMET\_cmorph](\ref{AGRMET_cmorph}) \newline
!      Interpolates CMORPH data to LIS grid
!  \end{description}
!EOP

  integer :: yr, mo, da, hr
  integer, intent(in)               :: j3hr  
  integer :: endtime_cmor     ! 1=get a new file
  integer :: sectionofcmor
  real*8  :: ctime,ftime_cmor       ! Current LDAS time and end boundary times for precip data sources 
  real*8  :: datatime, gap, breaktime, fnametime                    ! Times used in HUFFMAN to determine data and filename boundaries (see below)
  real    :: gmt1,gmt4
  character(len=120) :: name ! Filename variables for precip data sources

!=== End Variable Definition =======================

!------------------------------------------------------------------------
! Set parameter to measure 15 minutes time offset when using CMORPH
!------------------------------------------------------------------------
  gap = 2.8539E-5

!=== Assumption is to not find any data
 endtime_cmor = 0

!------------------------------------------------------------------------
! Determine required observed precip data times 
! (current, accumulation end time)
! Model current time
!------------------------------------------------------------------------
  call LIS_julhr_date( j3hr, yr,mo,da,hr)

!------------------------------------------------------------------------
! Check for and get CMORPH precipitation data
!------------------------------------------------------------------------
  call cmorfile_agrmet( name, agrmet_struc(n)%agrmetdir, agrmet_struc(n)%cmordir, &
    agrmet_struc(n)%use_timestamp, yr, mo, da, hr )
  call AGRMET_cmorph( n, name, cmorphdata, quad9r)
return
end subroutine AGRMET_getcmorph

!BOP
! !ROUTINE: cmorfile_agrmet
! \label{cmorfile_agrmet}
!
!
! !INTERFACE:
subroutine cmorfile_agrmet( name, agrmetdir, cmordir, use_timestamp, yr, mo, da, hr)

  implicit none
! !ARGUMENTS: 
  character(len=*) :: name
  character(len=*) :: agrmetdir
  character(len=*) :: cmordir
  integer, intent(in) :: use_timestamp
  integer            :: yr, mo, da, hr
! !DESCRIPTION:
!   This subroutine puts together CMORPH file name for 
!   30min file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[agrmetdir]
!    Name of the AGRMET forcing directory
!  \item[cmordir]
!    Name of the CMORPH directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[name]
!   name of the timestamped CMORPH file
!  \end{description}
!
!EOP

  character(len=120) :: temp
  integer :: i, c, d
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
  character*1 :: fdir(99), fbase(99), fdir2(8),ftime(10), ftimedir(10)

!=== End Variable Definition ===============
!=== formats for filename segments
90 format (a120)
91 format (120a1)
92 format (99a1)
93 format (a99)
94 format (i4, i2.2, i2.2, i2.2)
95 format (10a1)
96 format (a40)
97 format (a8)
98 format (a1, i4, i2, a1)
99 format (8a1)
89 format (a7)
88 format (a8)
86 format (a1,i4,i2.2,i2.2,a1)
87 FORMAT (10a1)
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

!       print *,'In cmorfile: cmordir=',cmordir
  write(temp, 93) agrmetdir
  read(temp, 92) (fdir(i), i=1,99)

  write(temp, 93) cmordir
  read(temp, 92) (fbase(i), i=1,99)

  write(temp, 88) '/CMORPH_'
  read(temp, 92) (fdir2(i), i=1,8)
  
 
  write(temp, 94) uyr, umo, uda, uhr
  read(temp, 95) ftime

  write(temp, 86) '/',uyr, umo, uda, '/'
  read(temp, 87) ftimedir

  do i = 1, 10
   if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do

  c = 0
  d = 0
  do i = 1, 90
    if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
    if ( (fdir(i)  == ' ') .and. (d == 0) ) d = i-1
  end do

  if (use_timestamp .eq. 1) then
    write(temp, 91) (fdir(i),i=1,d),(ftimedir(i),i=1,10),(fbase(i),i=1,c),(fdir2(i),i=1,8),  &
                       (ftime(i), i=1,10)
  else
    write(temp, 91) (fdir(i),i=1,d),(fbase(i),i=1,c),(fdir2(i),i=1,8),  &
                       (ftime(i), i=1,10)
  endif

  read(temp, 90) name
  return

end subroutine cmorfile_agrmet
