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
! !ROUTINE: netcdf4_nldas20filea
! \label{netcdf4_nldas20filea}
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  04 Sep 2001: Brian Cosgrove; Use of NASA data enabled, updated
!               reading of data directory structure to read new format
!  24 Aug 2007: Chuck Alonge; Modified for use with NLDAS-2 data
!  22 Jan 2012: K. Arsenault; Accommodate GES DISC, NCEP/EMC filename conventions
!  11 Jul 2024: David Mocko,  Modified for use with netCDF-4 format
!
! !INTERFACE:
      subroutine netcdf4_nldas20filea(n,kk,findex,filename,            &
                                      nldas20dir,yr,mo,da,doy,hr)
! !USES:
      use LIS_coreMod
      use LIS_logMod
      use LIS_forecastMod

      implicit none
! !ARGUMENTS:
      integer                       :: n
      integer                       :: kk
      integer                       :: findex
      character(len=*), intent(out) :: filename
      character(len=*), intent(in)  :: nldas20dir
      integer, intent(in)           :: yr,mo,da,doy,hr
!
! !DESCRIPTION:
!   This subroutine assembles GES DISC NLDAS-2 netCDF-4 "A" filenames
!     for 1 hour file intervals.
!
!  The arguments are:
!  \begin{description}
!  \item[nldas20dir]
!    Name of the NLDAS-2 directory
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
!   name of the timestamped GES DISC NLDAS-2 netCDF-4 file
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

!=== Assemble GES DISC NLDAS-2 netCDF-4 FORA filename:
         filename = trim(nldas20dir)//"/"//fyr//"/"//fdoy//            &
                    "/NLDAS_FORA0125_H.A"//fyr//fmo//fda//"."//fhr//   &
                    "00.020.nc"

      else ! forecast mode
         doy2 = doy
! IF Forecast Year is a leap years, for doy:
         if ((mod(LIS_rc%yr,4).eq.0.and.mod(LIS_rc%yr,100).ne.0)       &
             .or.(mod(LIS_rc%yr,400).eq.0)) then
            if (doy.gt.59) then  ! Make sure to remove extra day
               doy2 = doy - 1
            endif
         endif

! Sample yr, mo, da
         call LIS_sample_forecastDate(n,kk,findex,yr,mo,da)

! Account for member year - leap years for doy:
         if ((mod(yr,4).eq.0.and.mod(yr,100).ne.0)                     &
             .or.(mod(yr,400).eq.0)) then
            if (doy.gt.59) then
               doy2 = doy2 + 1
            endif
         endif

         write(unit=fdoy,fmt="(i3.3)") doy2
         write(unit=fyr, fmt="(i4.4)") yr
         write(unit=fmo, fmt="(i2.2)") mo
         write(unit=fda, fmt="(i2.2)") da
         write(unit=fhr, fmt="(i2.2)") hr

!=== Assemble GES DISC NLDAS-2 netCDF-4 FORA filename:
         filename = trim(nldas20dir)//"/"//fyr//"/"//fdoy//            &
                    "/NLDAS_FORA0125_H.A"//fyr//fmo//fda//"."//fhr//   &
                    "00.020.nc"
      endif

      end subroutine netcdf4_nldas20filea

!BOP
! !ROUTINE: netcdf4_nldas20fileb
! \label{netcdf4_nldas20fileb}
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  04 Sep 2001: Brian Cosgrove; Use of NASA data enabled, updated
!               reading of data directory structure to read new format
!  24 Aug 2007: Chuck Alonge; Modified for use with NLDAS-2 data
!  22 Jan 2012: K. Arsenault; Accommodate GES DISC, GES DISC filename conventions
!  11 Jul 2024: David Mocko,  Modified for use with netCDF-4 format
!
! !INTERFACE:
      subroutine netcdf4_nldas20fileb(n,kk,findex,filename,            &
                                      nldas20dir,yr,mo,da,doy,hr)
! !USES:
      use LIS_coreMod
      use LIS_logMod
      use LIS_forecastMod

      implicit none
! !ARGUMENTS:
      integer                       :: n
      integer                       :: kk
      integer                       :: findex
      character(len=*), intent(out) :: filename
      character(len=*), intent(in)  :: nldas20dir
      integer, intent(in)           :: yr,mo,da,doy,hr
!
! !DESCRIPTION:
!   This subroutine assembles GES DISC NLDAS-2 netCDF-4 "B" filenames
!     for 1 hour file intervals.
!
!  The arguments are:
!  \begin{description}
!  \item[nldas20dir]
!    Name of the NLDAS-2 directory
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
!  \item[filename]
!   name of the timestamped GES DISC NLDAS-2 netCDF-4 file
!  \end{description}
!
!EOP
      character*4  :: fyr
      character*3  :: fdoy
      character*2  :: fmo, fda, fhr
      integer      :: doy2

!=== end variable definition =============================================

      if (LIS_rc%forecastMode.eq.0) then !hindcast run
         write(unit=fyr, fmt="(i4.4)") yr
         write(unit=fdoy,fmt="(i3.3)") doy
         write(unit=fmo, fmt="(i2.2)") mo
         write(unit=fda, fmt="(i2.2)") da
         write(unit=fhr, fmt="(i2.2)") hr

!=== Assemble GES DISC NLDAS-2 netCDF-4 FORB filename:
         filename = trim(nldas20dir)//"/"//fyr//"/"//fdoy//            &
                    "/NLDAS_FORB0125_H.A"//fyr//fmo//fda//"."//fhr//   &
                    "00.020.nc"

      else !forecast mode
         doy2 = doy
! IF Forecast Year is a leap years, for doy:
         if ((mod(LIS_rc%yr,4).eq.0.and.mod(LIS_rc%yr,100).ne.0)       &
             .or.(mod(LIS_rc%yr,400).eq.0)) then
            if (doy.gt.59) then  ! Make sure to remove extra day
               doy2 = doy - 1
            endif
         endif

! Sample yr, mo, da
         call LIS_sample_forecastDate(n,kk,findex,yr,mo,da)

! Account for member year - leap years for doy:
         if ((mod(yr,4).eq.0.and.mod(yr,100).ne.0)                     &
             .or.(mod(yr,400).eq.0)) then
            if (doy.gt.59) then
               doy2 = doy + 1
            endif
         endif

         write(unit=fdoy,fmt="(i3.3)") doy2
         write(unit=fyr, fmt="(i4.4)") yr
         write(unit=fmo, fmt="(i2.2)") mo
         write(unit=fda, fmt="(i2.2)") da
         write(unit=fhr, fmt="(i2.2)") hr

!=== Assemble GES DISC NLDAS-2 netCDF-4 FORB filename:
         filename = trim(nldas20dir)//"/"//fyr//"/"//fdoy//            &
                    "/NLDAS_FORB0125_H.A"//fyr//fmo//fda//"."//fhr//   &
                    "00.020.nc"
      endif

      end subroutine netcdf4_nldas20fileb

