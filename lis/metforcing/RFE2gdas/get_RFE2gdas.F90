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
! !ROUTINE: get_RFE2gdas
!  \label{get_RFE2gdas}
!
! !REVISION HISTORY:
!  30 May 2010; Soni Yatheendradas, Initial LIS version for FEWSNET
!  25 Jan 2011: Clement Alo modified for RFE2gdas 6-hourly data
!   8 Mar 2012: Updated for LIS7
!
! !INTERFACE:
subroutine get_RFE2gdas(n, findex)

! !USES:
  use LIS_coreMod,         only : LIS_rc
  use LIS_timeMgrMod,      only : LIS_calendar, LIS_get_nstep, LIS_tick
  use LIS_logMod,          only : LIS_logunit, LIS_endrun, LIS_verify
  use LIS_constantsMod,    only : LIS_CONST_PATH_LEN
  use RFE2gdas_forcingMod, only : RFE2gdas_struc

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and upscales/interpolates 6-hourly RFE2gdas.0 forcing.
!  At the beginning of a simulation, the code
!  reads the relevant daily data (nearest daily interval).
!  Since rain value is obtained from total daily accumulation,
!  only one relevant file needs to be read, not two. 
!  Current first-cut temporal interpolation implemented is
!  even distribution over sub-6hourly model run timesteps.
!  The strategy for missing data in spatial upscaling is to retain 
!  the corresponding base precip value/s (Diurnal precip pattern 
!  then follows base precip and is not an even distribution). Upscaled 
!  values are missing only if ALL the contained RFE2gdas pixels are NODATA. 
!  The strategy for missing data in spatial interpolation is to set
!  such NODATA values in the RFE2gdas input to 0 (Ronald W Lietzow,
!  USGS, Personal E-mail Communication 05/27/2010). An upper daily 1000 mm
!  cap is also implemented as per Ron's info of current practice before
!  being input to the WRSI, MI, and BERM products 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the supplemental forcing source
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the RFE2gdas data times
!  \item[RFE2gdasfile](\ref{RFE2gdasfile}) \newline
!    Puts together appropriate file name for 6 hour intervals 
!  \item[readprecip\_RFE2gdas](\ref{readprecip_RFE2gdas}) \newline
!    reads the RFE2gdas data 
!  \end{description}
!EOP

!==== Local Variables=======================
  integer :: IsEndTime_RFE2gdas           ! 1=get a new file
  integer :: doyNow, yrNow, moNow, daNow, hrNow, mnNow, ssNow, tsNow  ! Time parameters for current LDAS time
  integer :: yr1, mo1, da1,hr1            ! Time parameters for RFE2gdas time start boundary
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2                  ! Time parameters for RFE2gdas time end boundary
  real*8  :: ctime,EndTime_RFE2gdas       ! Current LDAS time and end boundary time for precip data source
  real    :: gmtNow,gmt1,gmt2             ! GMT times for current LDAS time and begin and end boundary times for precip data sources
  integer :: ferror_RFE2gdas              ! Error flags for precip data sources
  character(len=LIS_CONST_PATH_LEN):: filename                ! Filename variables for precip data sources
  integer      :: order,kk

!=== End Variable Definition =======================

  IsEndTime_RFE2gdas = 0

!------------------------------------------------------------------------
! Determine required observed precip data times 
! (current, accumulation end, and filename timestamp start time)
! Model current time
!------------------------------------------------------------------------
  yrNow = LIS_rc%yr  !current time
  moNow = LIS_rc%mo
  daNow = LIS_rc%da
  hrNow = LIS_rc%hr
  mnNow = LIS_rc%mn
  ssNow = 0
  tsNow = 0
  call LIS_tick( ctime, doyNow, gmtNow, yrNow, moNow, daNow, hrNow, &
       mnNow, ssNow, real(tsNow))
  
!------------------------------------------------------------------------
!   RFE2gdas product end time 
!------------------------------------------------------------------------
    yr2 = LIS_rc%yr  !accumulation end time parameters
    mo2 = LIS_rc%mo 
    da2 = LIS_rc%da 
    hr2 = 6*(LIS_rc%hr/6)
    mn2 = 0
    ss2 = 0
    ts2 = 6*60*60
    call LIS_tick( EndTime_RFE2gdas, doy2, gmt2, yr2, mo2, da2, hr2, &
         mn2, ss2, real(ts2 ))

!------------------------------------------------------------------------
! Ensure that data is found during first time step
!------------------------------------------------------------------------
  if( (LIS_get_nstep(LIS_rc,n).eq. 1) .OR. &
       (LIS_rc%rstflag(n) .eq. 1) ) then
     IsEndTime_RFE2gdas = 1
     LIS_rc%rstflag(n) = 0
  endif

!------------------------------------------------------------------------
! Check for and get RFE2gdas Precipitation data
!------------------------------------------------------------------------
  filename=''
  if ( ctime > RFE2gdas_struc(n)%RFE2gdasEndTime ) IsEndTime_RFE2gdas = 1
  if ( IsEndTime_RFE2gdas == 1 ) then  !get new time2 data
     do kk=RFE2gdas_struc(n)%st_iterId, RFE2gdas_struc(n)%en_iterId
        call RFE2gdasfile(n, kk, findex, filename, &
             RFE2gdas_struc(n)%RFE2gdasDir, &
             yr2, mo2, da2, hr2 )
        ferror_RFE2gdas = 0
        order = 2
        call readprecip_RFE2gdas( n, kk,filename, mo2, findex,&
             order,ferror_RFE2gdas,hr2 )
     enddo
     RFE2gdas_struc(n)%RFE2gdasEndTime = EndTime_RFE2gdas
  endif  ! need new time2

end subroutine get_RFE2gdas


!BOP
! !ROUTINE: RFE2gdasfile
! \label{RFE2gdasfile}
!
! !INTERFACE:
subroutine RFE2gdasfile( n, kk, findex, filename, RFE2gdasDir, yr, mo, da, hr ) 

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_forecastMod

  implicit none
! !ARGUMENTS: 
  integer            :: n
  integer            :: kk
  integer            :: findex
  character(len=*)   :: filename
  character(len=*)   :: RFE2gdasDir
  integer            :: yr, mo, da, hr
!
! !DESCRIPTION:
!   This subroutine puts together RFE2-GDAS/CMAP file name for 
!   6 hour file intervals
!
!  The arguments are:
!  \begin{description}
!  \item[RFE2gdasDir]
!    Name of the RFE2gdas directory
!  \item[yr]
!    year 
!  \item[mo]
!    month
!  \item[da]
!    day of month
!  \item[hr]
!    hour of day
!  \item[filename]
!    name of the timestamped RFE2gdas file
!  \end{description}
!
!EOP

!==== Local Variables=======================
!  integer :: uyr, umo, uda,uhr, umn, uss, ts1
!  character(len=100) :: temp
!  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(13), fsubs2(4)
!  integer :: i, c
!=== End Variable Definition =======================

   character*4  :: fyr
   character*2  :: fmo, fda, fhr

  !=== end variable definition =============================================

  if(LIS_rc%forecastMode.eq.0) then ! hindcast run

    write(unit=fyr, fmt='(i4.4)')  yr
    write(unit=fmo, fmt='(i2.2)')  mo
    write(unit=fda, fmt='(i2.2)')  da
    write(unit=fhr, fmt='(i2.2)')  hr

  !=== Assemble RFE2-GDAS filename:

!    ./201009/rfe_gdas.bin.2010090706

    filename = trim(RFE2gdasDir)//"/"//fyr//fmo//"/rfe_gdas.bin."//&
               fyr//fmo//fda//fhr

  else ! Forecast run

    call LIS_sample_forecastDate(n, kk, findex, yr,mo,da)

    write(unit=fyr, fmt='(i4.4)')  yr
    write(unit=fmo, fmt='(i2.2)')  mo
    write(unit=fda, fmt='(i2.2)')  da
    write(unit=fhr, fmt='(i2.2)')  hr
 
  !=== Assemble RFE2-GDAS filename:

!   ./201009/rfe_gdas.bin.2010090706
    filename = trim(RFE2gdasDir)//"/"//fyr//fmo//"/rfe_gdas.bin."//&
              fyr//fmo//fda//fhr

  endif  ! End forecast mode


end subroutine RFE2gdasfile


