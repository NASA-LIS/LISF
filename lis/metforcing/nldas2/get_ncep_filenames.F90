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
! !ROUTINE: ncep_nldas2filea
! \label{ncep_nldas2filea}
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  04 Sep 2001: Brian Cosgrove; Use of NASA data enabled, updated
!               reading of data directory structure to read new format
!  24 Aug 2007: Chuck Alonge; Modified for use with NLDAS-2 data
!  22 Jan 2012: K. Arsenault; Accommodate GES DISC, NCEP filename conventions 
!
! !INTERFACE:
 subroutine ncep_nldas2filea(n,kk,findex,filename,nldas2dir,yr,mo,da,hr)

! !USES:
   use LIS_coreMod
   use LIS_logMod
   use LIS_forecastMod

   implicit none
! !ARGUMENTS: 
   integer, intent(in)        :: n 
   integer, intent(in)        :: kk
   integer, intent(in)        :: findex
   character(len=*), intent(out) :: filename
   character(len=*), intent(in)   :: nldas2dir
   integer, intent(in)        :: yr,mo,da,hr

! !DESCRIPTION:
!   This subroutine puts together NCEP NLDAS-2 ``A'' file name for 
!   1 hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[nldas2dir]
!    Name of the NLDAS-2 directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[filename]
!   name of the timestamped NCEP NLDAS-2 file
!  \end{description}
!
!EOP

   character(len=13) :: fdir
   character(len=17), parameter :: fsubs = '.nldasforce-a.grb'
   character(len=10) :: ftime

   !=== end variable definition =============================================

   !=== put together filename

   if(LIS_rc%forecastMode.eq.0) then !hindcast run

      write(UNIT=fdir,fmt='(i4,a1,i4,i2.2,i2.2)') yr, '/', yr, mo, da

      write(UNIT=ftime,fmt='(i4,i2.2,i2.2,i2.2)') yr, mo, da, hr

   else !forecast mode
     !sample yr, mo, da

      call LIS_sample_forecastDate(n, kk, findex, yr, mo, da)
      
      write(UNIT=fdir,fmt='(i4,a1,i4,i2.2,i2.2)') yr, '/', yr, mo, da
      
      write(UNIT=ftime,fmt='(i4,i2.2,i2.2,i2.2)') yr, mo, da, hr
      
   endif
   
   filename = trim(nldas2dir) // '/' // fdir // '/' // ftime // fsubs
  
end subroutine ncep_nldas2filea


!BOP
! !ROUTINE: ncep_nldas2fileb
! \label{ncep_nldas2fileb}
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  04 Sep 2001: Brian Cosgrove; Use of NASA data enabled, updated
!               reading of data directory structure to read new format
!  24 Aug 2007: Chuck Alonge; Modified for use with NLDAS-2 data
!  22 Jan 2012: K. Arsenault; Accommodate GES DISC, NCEP filename conventions 
!
! !INTERFACE:
 subroutine ncep_nldas2fileb(n,kk,findex,filename,nldas2dir,yr,mo,da,hr)

! !USES:
   use LIS_coreMod
   use LIS_logMod
   use LIS_forecastMod

   implicit none
! !ARGUMENTS: 
   integer, intent(in)        :: n
   integer, intent(in)        :: kk
   integer, intent(in)        :: findex
   character(len=*), intent(out) :: filename
   character(len=*), intent(in)   :: nldas2dir
   integer, intent(in)        :: yr,mo,da,hr
!
! !DESCRIPTION:
!   This subroutine puts together NCEP NLDAS-2 ``B'' file name for 
!   1 hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[nldas2dir]
!    Name of the NLDAS-2 directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[filename]
!   name of the timestamped NCEP NLDAS-2 file
!  \end{description}
!
!EOP

   character(len=13) :: fdir
   character(len=17), parameter :: fsubs = '.nldasforce-b.grb'
   character(len=10) :: ftime

   !=== end variable definition =============================================

   !=== put together filename

   if(LIS_rc%forecastMode.eq.0) then !hindcast run

      write(UNIT=fdir,fmt='(i4,a1,i4,i2.2,i2.2)') yr, '/', yr, mo, da
      
      write(UNIT=ftime,fmt='(i4,i2.2,i2.2,i2.2)') yr, mo, da, hr
      
   else !forecast mode
     !sample yr, mo, da

      call LIS_sample_forecastDate(n, kk, findex, yr, mo, da)

      write(UNIT=fdir,fmt='(i4,a1,i4,i2.2,i2.2)') yr, '/', yr, mo, da
     
      write(UNIT=ftime,fmt='(i4,i2.2,i2.2,i2.2)') yr, mo, da, hr

   endif

   filename = trim(nldas2dir) // '/' // fdir // '/' // ftime // fsubs
   
end subroutine ncep_nldas2fileb

