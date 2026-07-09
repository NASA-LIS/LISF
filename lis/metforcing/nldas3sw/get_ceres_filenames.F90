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
! !ROUTINE: ceres_nldas3swfile
! \label{ceres_nldas3swfile}
!
! !REVISION HISTORY:
! 27 Dec 2024: David Mocko, Initial Specification
!                           (derived from get_netcdf4_filenames.F90)
!
! !INTERFACE:
subroutine ceres_nldas3swfile(n,kk,findex,filename,                    &
                              nldas3swdir,yr,mo,da,doy,hr)
! !USES:
  use LIS_coreMod
  use LIS_logMod

  implicit none
! !ARGUMENTS:
  integer                       :: n
  integer                       :: kk
  integer                       :: findex
  character(len=*), intent(out) :: filename
  character(len=*), intent(in)  :: nldas3swdir
  integer, intent(in)           :: yr,mo,da,doy,hr
!
! !DESCRIPTION:
!   This subroutine assembles CERES 4-km binary filename
!     for 1 hour file intervals.
!
!  The arguments are:
!  \begin{description}
!  \item[nldas3swdir]
!    Name of the CERES SWdown directory
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[doy]
!   Julian day of year (needed for subdirectory structure)
!  \item[hr]
!   hour of day
!   \item[filename]
!   name of the timestamped CERES 4-km binary filename
!  \end{description}
!
!EOP
  character*4  :: fyr
  character*3  :: fdoy
  character*2  :: fmo, fda, fhr
  integer      :: doy2

!=== end variable definition ===========================================

  if (LIS_rc%forecastMode.eq.0) then !hindcast run
     write(unit=fyr, fmt="(i4.4)") yr
     write(unit=fdoy,fmt="(i3.3)") doy
     write(unit=fmo, fmt="(i2.2)") mo
     write(unit=fda, fmt="(i2.2)") da
     write(unit=fhr, fmt="(i2.2)") hr

!=== Assemble CERES 4-km binary filename
     filename = trim(nldas3swdir)//"/"//fyr//"/"//                     &
                    "NLDAS3_INSOL_DAILY_"//fyr//fdoy//".dat"
  endif

end subroutine ceres_nldas3swfile

