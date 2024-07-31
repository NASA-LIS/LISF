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
! !ROUTINE: get_metForcGen_filename
! \label{get_metForcGen_filename}
!
! !REVISION HISTORY:
!  06 Jan 2015: K. Arsenault; Initial Implementation
!
! !INTERFACE:
  subroutine get_metForcGen_filename(n,kk,findex,&
                 yr,mo,da,hr,mn,directory,filename)

! !USES:
   use LIS_coreMod
   use LIS_logMod
   use LIS_forecastMod

   implicit none
! !ARGUMENTS: 
   integer, intent(in)        :: n               ! Nest index
   integer, intent(in)        :: kk              ! Forecast index
   integer, intent(in)        :: findex          ! Forcing index
   integer, intent(in)        :: yr,mo,da,hr,mn  ! File timestamp
   character(len=*), intent(in)  :: directory       ! File directory
   character(len=*), intent(out) :: filename  
!
! !DESCRIPTION:
!   This subroutine puts together LDT generated forcing
!    file names.
! 
!  The arguments are:
!  \begin{description}
!  \item[directory]
!    Name of the NLDAS-2 directory
!  \item[n]
!    nest index
!  \item[kk]
!    forecast member index
!  \item[findex]
!    forcing dataset index
!  \item[yr]
!    year of filename
!  \item[mo]
!    month of filename
!  \item[da]
!    day of month (of filename)
!  \item[hr]
!    hour of day (of filename)
!  \item[mn]
!    minute of day (of filename)
!   \item[filename]
!    name of the time-stamped LDT-generated forcing file
!  \end{description}
!
!EOP
   character*4  :: fyr
   character*2  :: fmo, fda, fhr, fmn
   character*2  :: fnest

  !=== end variable definition =============================================

  !=== Assemble LDT-generated filename:
  !   ./DIRECTORY/FORCING/LDT_HIST_200505011500.d01.nc

   if(LIS_rc%forecastMode.eq.0) then ! hindcast run

     write(unit=fyr, fmt='(i4.4)')  yr
     write(unit=fmo, fmt='(i2.2)')  mo
     write(unit=fda, fmt='(i2.2)')  da
     write(unit=fhr, fmt='(i2.2)')  hr
     write(unit=fmn, fmt='(i2.2)')  mn
   !  write(unit=fnest,fmt='(i2.2)') n
 
     filename = trim(directory)//"/FORCING/"//fyr//fmo//"/LDT_HIST_"//&
                fyr//fmo//fda//fhr//fmn//".d01.nc"
   !             fyr//fmo//fda//fhr//fmn//".d"//fnest//".nc"

   else ! Forecast run
     call LIS_sample_forecastDate(n, kk, findex, yr,mo,da)

     write(unit=fyr, fmt='(i4.4)')  yr
     write(unit=fmo, fmt='(i2.2)')  mo
     write(unit=fda, fmt='(i2.2)')  da
     write(unit=fhr, fmt='(i2.2)')  hr
     write(unit=fmn, fmt='(i2.2)')  mn
   !  write(unit=fnest,fmt='(i2.2)') n

     filename = trim(directory)//"/FORCING/"//fyr//fmo//"/LDT_HIST_"//&
                fyr//fmo//fda//fhr//fmn//".d01.nc"
   !             fyr//fmo//fda//fhr//fmn//".d"//fnest//".nc"
   endif

 end subroutine get_metForcGen_filename

