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
! !ROUTINE: get_imerg
! \label{get_imerg}
!
!
! !REVISION HISTORY:
! 17 Jul 2001: Jon Gottschalck; Initial code
! 05 Mar 2015: Jonathan Case; modified for GPM IMERG
!
! !INTERFACE:
subroutine get_imerg(n, findex)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only : LIS_tick, LIS_get_nstep
  use imerg_forcingMod, only :imerg_struc
  use LIS_logMod, only : LIS_logunit, LIS_endrun
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 30 minute IMERG forcing. 
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
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the IMERG data times
!  \item[imergfile](\ref{imergfile}) \newline
!    Puts together appropriate file name for 6 hour intervals
!  \item[read\_imerg](\ref{read_imerg}) \newline
!    Interpolates IMERG data to LIS grid
!  \end{description}
!EOP

  integer :: kk
  integer :: ferror_imerg      ! Error flags for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy4, yr4, mo4, da4, hr4, mn4, ss4
  integer :: endtime_imerg     ! 1=get a new file
  integer :: sectionofimerg
  real*8  :: ctime,ftime_imerg       ! Current LIS time and end boundary times for precip data sources 
  integer :: order
  real    :: gmt1,gmt4,ts1,ts4
  character(len=LIS_CONST_PATH_LEN) :: filename ! Filename variables for precip data sources

!=== End Variable Definition =======================

  endtime_imerg = 0

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
! Ensure that data is found during first time step
!------------------------------------------------------------------------
  endtime_imerg = 0
  if ( LIS_get_nstep(LIS_rc,n).eq. 1 .or.LIS_rc%rstflag(n) .eq. 1) then
     endtime_imerg = 1
     LIS_rc%rstflag(n) = 0
  endif

!------------------------------------------------------------------------
! Check for and get IMERG precipitation data
!------------------------------------------------------------------------
  if ( ctime > imerg_struc(n)%imergtime ) then
     endtime_imerg = 1
  endif
   
  if ( LIS_rc%mn >= 0 .AND. LIS_rc%mn < 30 ) sectionofimerg = 1
  if ( LIS_rc%mn >= 30 ) sectionofimerg = 2
  if ( endtime_imerg == 1 ) then  !get new time2 data first 1/2 hour
    ferror_imerg = 0
    order = 2
    ! Loop over number of ensemble members:
    do kk= imerg_struc(n)%st_iterid, imerg_struc(n)%en_iterid
       call imergfile(n, kk, findex, imerg_struc(n)%imergdir, &
                      yr1,mo1,da1,hr1,sectionofimerg,filename )
       call read_imerg( n, kk, filename, findex, order, ferror_imerg )
    end do
    imerg_struc(n)%imergtime = ctime
  endif  !need new time2

end subroutine get_imerg

!BOP
! !ROUTINE: imergfile
! \label{imergfile}
!
!
! !INTERFACE:
subroutine imergfile(n, kk, findex, imergdir, &
                yr, mo, da, hr, section, filename)

  use LIS_coreMod
  use LIS_forecastMod
  use LIS_logMod, only : LIS_logunit, LIS_endrun
  use imerg_forcingMod, only : imerg_struc
  implicit none

! !ARGUMENTS: 
  integer          :: n
  integer          :: kk
  integer          :: findex
  character(len=*) :: imergdir
  integer          :: yr, mo, da, hr, section
  character(len=*) :: filename
!
! !DESCRIPTION:
!   This subroutine puts together IMERG file name for 30-min file intervals.
! 
!  The arguments are:
!  \begin{description}
!  \item[imergdir]
!    Name of the IMERG directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!  \item[section]
!   section to denote the half-hourly dataset (first 30-min [1] or last 30-min [2])
!   \item[filename]
!   name of the timestamped IMERG file
!  \end{description}
!
!EOP
  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, umnadd, umnday, uss !, ts1

  character*100 :: fstem, fext
  character*4   :: cyr, cmnday, imVer
  character*2   :: cmo, cda, chr, cmn, cmnadd 

!=== End Variable Definition ===============
!------------------------------------------------------------------------
! Make variables for the time used to create the file
! We don't want these variables being passed out
!------------------------------------------------------------------------

  if(LIS_rc%forecastMode.eq.0) then !hindcast run
    uyr = yr
    umo = mo
    uda = da
    uhr = 1*(hr/1)  !hour needs to be a multiple of 1 hour
    if (section .eq. 1) then
      umn = 0
    else
      umn = 30
    endif
    umnadd = umn + 29
    umnday = uhr*60 + umn
    uss = 0
    write(cyr, '(I4.4)') uyr
    write(cmo, '(I2.2)') umo 
    write(cda, '(I2.2)') uda 
    write(chr, '(I2.2)') uhr 
    write(cmn, '(I2.2)') umn 
    write(cmnadd, '(I2.2)') umnadd 
    write(cmnday, '(I4.4)') umnday

    if(imerg_struc(n)%imergprd == 'early') then
       fstem = '/3B-HHR-E.MS.MRG.3IMERG.'
       fext  = '.RT-H5'
    elseif(imerg_struc(n)%imergprd == 'late') then
       fstem = '/3B-HHR-L.MS.MRG.3IMERG.'
       fext  = '.RT-H5'
    elseif(imerg_struc(n)%imergprd == 'final') then
       fstem = '/3B-HHR.MS.MRG.3IMERG.'
       fext  = '.HDF5'
    else
       write(LIS_logunit,*) "[ERR] Invalid IMERG product option was chosen."
       write(LIS_logunit,*) "[ERR] Please choose either 'early', 'late', or 'final'."
       call LIS_endrun()
    endif
    imVer = trim(imerg_struc(n)%imergver)
    filename = trim(imergdir)//"/"//cyr//cmo//trim(fstem)// &
          cyr//cmo//cda//"-S"//chr//cmn//"00-E"//chr//cmnadd//"59."//cmnday//"."//imVer//fext

! Forecast mode (e.g., ESP):
  else
    call LIS_sample_forecastDate(n,kk,findex,yr,mo,da)

    uyr = yr
    umo = mo
    uda = da
    uhr = 1*(hr/1)  !hour needs to be a multiple of 1 hour
    if (section .eq. 1) then
      umn = 0
    else
      umn = 30
    endif
    umnadd = umn + 29
    umnday = uhr*60 + umn
    uss = 0
    write(cyr, '(I4.4)') uyr
    write(cmo, '(I2.2)') umo
    write(cda, '(I2.2)') uda
    write(chr, '(I2.2)') uhr
    write(cmn, '(I2.2)') umn
    write(cmnadd, '(I2.2)') umnadd
    write(cmnday, '(I4.4)') umnday

    if(imerg_struc(n)%imergprd == 'early') then
       fstem = '/3B-HHR-E.MS.MRG.3IMERG.'
       fext  = '.RT-H5'
    elseif(imerg_struc(n)%imergprd == 'late') then
       fstem = '/3B-HHR-L.MS.MRG.3IMERG.'
       fext  = '.RT-H5'
    elseif(imerg_struc(n)%imergprd == 'final') then
       fstem = '/3B-HHR.MS.MRG.3IMERG.'
       fext  = '.HDF5'
    else
       write(LIS_logunit,*) "[ERR] Invalid IMERG product option was chosen."
       write(LIS_logunit,*) "[ERR] Please choose either 'early', 'late', or 'final'."
       call LIS_endrun()
    endif
    imVer = trim(imerg_struc(n)%imergver)
    filename = trim(imergdir)//"/"//cyr//cmo//trim(fstem)// &
          cyr//cmo//cda//"-S"//chr//cmn//"00-E"//chr//cmnadd//"59."//cmnday//"."//imVer//fext
  endif

end subroutine imergfile
