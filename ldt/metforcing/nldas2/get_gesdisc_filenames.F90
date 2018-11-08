!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT)
!
! See RELEASE_NOTES.txt for more information.
!
! The LDT source code and documentation are not in the public domain
! and may not be freely distributed.  Only qualified entities may receive 
! the source code and documentation. 
!
! Qualified entities must be covered by a Software Usage Agreement. 
! The Software Usage Agreement contains all the terms and conditions
! regarding the release of the LDT software.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See the Software Usage Agreement for the full disclaimer of warranty.
!
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: gesdisc_nldas2filea
! \label{gesdisc_nldas2filea}
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  04 Sep 2001: Brian Cosgrove; Use of NASA data enabled, updated
!               reading of data directory structure to read new format
!  24 Aug 2007: Chuck Alonge; Modified for use with NLDAS-2 data
!  22 Jan 2012: K. Arsenault; Accommodate GES DISC, NCEP/EMC filename conventions 
!
! !INTERFACE:
 subroutine gesdisc_nldas2filea(filename,nldas2dir,yr,mo,da,doy,hr)

   implicit none
! !ARGUMENTS: 
   character*80, intent(out) :: filename
   character*40, intent(in)  :: nldas2dir
   integer, intent(in)       :: yr,mo,da,doy,hr
!
! !DESCRIPTION:
!   This subroutine puts together GES DISC NLDAS-2 "A" file name for 
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
!  \item[doy]
!   Julian day of year (needed for subdirectory structure) 
!  \item[hr]
!   hour of day
!   \item[filename]
!   name of the timestamped GES DISC NLDAS-2 file
!  \end{description}
!
!EOP
   character*4  :: fyr
   character*3  :: fdoy
   character*2  :: fmo, fda, fhr

  !=== end variable definition =============================================

   write(unit=fyr, fmt='(i4.4)')  yr
   write(unit=fdoy,fmt='(i3.3)')  doy
   write(unit=fmo, fmt='(i2.2)')  mo
   write(unit=fda, fmt='(i2.2)')  da
   write(unit=fhr, fmt='(i2.2)')  hr
 
  !=== Assemble GES DISC NLDAS-2 filename:

!  NLDAS_FOR[A]0125_H.A19850601.0000.002.grb

   filename = trim(nldas2dir)//"/"//fyr//"/"//fdoy//"/NLDAS_FORA0125_H.A"//&
              fyr//fmo//fda//"."//fhr//"00.002.grb"
   return

 end subroutine gesdisc_nldas2filea


!BOP
! !ROUTINE: gesdisc_nldas2fileb
! \label{gesdisc_nldas2fileb}
!
! !REVISION HISTORY:
!  1  Oct 1999: Jared Entin; Initial code
!  15 Oct 1999: Paul Houser; Significant F90 Revision
!  04 Sep 2001: Brian Cosgrove; Use of NASA data enabled, updated
!               reading of data directory structure to read new format
!  24 Aug 2007: Chuck Alonge; Modified for use with NLDAS-2 data
!  22 Jan 2012: K. Arsenault; Accommodate GES DISC, GES DISC filename conventions 
!
! !INTERFACE:
 subroutine gesdisc_nldas2fileb(filename,nldas2dir,yr,mo,da,doy,hr)

   implicit none
! !ARGUMENTS: 
   character*80, intent(out) :: filename
   character*40, intent(in)  :: nldas2dir
   integer, intent(in)       :: yr,mo,da,doy,hr
!
! !DESCRIPTION:
!   This subroutine puts together GES DISC NLDAS-2 "B" file name for 
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
!  \item[doy]
!   Julian day of year (needed for subdirectory structure)
!  \item[hr]
!   hour of day
!  \item[filename]
!   name of the timestamped GES DISC NLDAS-2 file
!  \end{description}
!
!EOP
   character*4  :: fyr
   character*3  :: fdoy
   character*2  :: fmo, fda, fhr

  !=== end variable definition =============================================

   write(unit=fyr, fmt='(i4.4)')  yr
   write(unit=fdoy,fmt='(i3.3)')  doy
   write(unit=fmo, fmt='(i2.2)')  mo
   write(unit=fda, fmt='(i2.2)')  da
   write(unit=fhr, fmt='(i2.2)')  hr

  !=== Assemble GES DISC NLDAS-2 filename:

!  NLDAS_FOR[B]0125_H.A19850601.0000.002.grb

   filename = trim(nldas2dir)//"/"//fyr//"/"//fdoy//"/NLDAS_FORB0125_H.A"//&
              fyr//fmo//fda//"."//fhr//"00.002.grb"

   return

 end subroutine gesdisc_nldas2fileb

