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
 subroutine ncep_nldas2filea(filename,nldas2dir,yr,mo,da,hr)

   implicit none
! !ARGUMENTS: 
   character*80, intent(out) :: filename
   character*40, intent(in)  :: nldas2dir
   integer, intent(in)       :: yr,mo,da,hr

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

   integer         :: i, c
   character*1     :: fbase(40),fsubs(17)
   character*1     :: ftime(10),fdir(15)
   character*100   :: temp

  !=== end variable definition =============================================

   !=== put together filename

   write(UNIT=temp,fmt='(a1,i4,a1,i4,i2,i2,a1)')'/',yr,'/',yr,mo,da,'/'
   read(UNIT=temp,fmt='(15a1)') fdir
   do i=1,15
      if(fdir(i).eq.(' ')) fdir(i) = '0'
   enddo
   write(UNIT=temp,fmt='(i4,i2,i2,i2)') yr,mo,da,hr
   read(UNIT=temp,fmt='(10a1)')ftime
   do i=1,10
      if(ftime(i).eq.(' ')) ftime(i) = '0'
   enddo
   write(UNIT=temp,fmt='(a17)')'.nldasforce-a.grb'
   read(UNIT=temp,fmt='(17a1)') fsubs
   write(UNIT=temp,fmt='(a40)') nldas2dir                       
   read(UNIT=temp,fmt='(40a1)') fbase
   c=0
   do i=1,40
      if(fbase(i).eq.(' ').and.c.eq.0)c=i-1
   enddo
   write(UNIT=temp,fmt='(80a1)')(fbase(i),i=1,c), (fdir(i),i=1,15), & 
        (ftime(i),i=1,10),(fsubs(i),i=1,17 ) 
   read(UNIT=temp,fmt='(a80)') filename

   return
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
 subroutine ncep_nldas2fileb(filename,nldas2dir,yr,mo,da,hr)

   implicit none
! !ARGUMENTS: 
   character*80, intent(out) :: filename
   character*40, intent(in)  :: nldas2dir
   integer, intent(in)       :: yr,mo,da,hr
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

   integer           :: i, c
   character*1       :: fbase(40),fsubs(17)
   character*1       :: ftime(10),fdir(15)
   character*100     :: temp
   !=== end variable definition =============================================

   !=== put together filename

   write(UNIT=temp,fmt='(a1,i4,a1,i4,i2,i2,a1)')'/',yr,'/',yr,mo,da,'/'
   read(UNIT=temp,fmt='(15a1)') fdir
   do i=1,15
      if(fdir(i).eq.(' ')) fdir(i) = '0'
   enddo
   write(UNIT=temp,fmt='(i4,i2,i2,i2)') yr,mo,da,hr
   read(UNIT=temp,fmt='(10a1)')ftime
   do i=1,10
      if(ftime(i).eq.(' ')) ftime(i) = '0'
   enddo
   write(UNIT=temp,fmt='(a17)')'.nldasforce-b.grb'
   read(UNIT=temp,fmt='(17a1)') fsubs
   write(UNIT=temp,fmt='(a40)') nldas2dir                       
   read(UNIT=temp,fmt='(40a1)') fbase
   c=0
   do i=1,40
      if(fbase(i).eq.(' ').and.c.eq.0)c=i-1
   enddo
   write(UNIT=temp,fmt='(80a1)')(fbase(i),i=1,c), (fdir(i),i=1,15), & 
        (ftime(i),i=1,10),(fsubs(i),i=1,17 ) 
   read(UNIT=temp,fmt='(a80)') filename

   return
 end subroutine ncep_nldas2fileb

