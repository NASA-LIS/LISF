!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: get_nldas1_filename
! \label{get_nldas1_filename}
!
! !REVISION HISTORY:
!  04 Sep 2001: Brian Cosgrove; Use of NASA data enabled, updated
!               reading of data directory structure to read new format
!  15 Feb 2012: K. Arsenault; Accommodate GES DISC, NCEP filename conventions
!  14 Mar 2014: David Mocko: Updated filename options for config file
!
! !INTERFACE:
 subroutine get_nldas1_filename(n,filename,nldas1dir,yr,mo,da,doy,hr)

   use nldas1_forcingMod,only : nldas1_struc

   implicit none
! !ARGUMENTS: 
   integer, intent(in) :: n
   character(len=*),intent(out)  :: filename
   character(len=*), intent(in)  :: nldas1dir
   integer,    intent(in)        :: yr,mo,da,doy,hr

! !DESCRIPTION:
!   This subroutine puts together NLDAS-1 filename structures for 
!   1 hour file intervals.
!   Two options for naming file conventions are given: \newline
!   1: NASA GES DISC filename type
!      (from http://disc.sci.gsfc.nasa.gov/hydrology/data-holdings) \newline
!   2: NASA-Original NLDAS-1 filename type (also NCEP-style) \newline
! 
!  The arguments are:
!  \begin{description}
!  \item[nldas1dir]
!    Name of the NLDAS-1 directory
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
!   name of the timestamped NLDAS-1 file
!  \end{description}
!
!EOP
   character*4  :: fyr
   character*3  :: fdoy
   character*2  :: fmo, fda, fhr

!=== end variable definition =============================================

   write(unit=fyr, fmt='(i4.4)') yr
   write(unit=fdoy,fmt='(i3.3)') doy
   write(unit=fmo, fmt='(i2.2)') mo
   write(unit=fda, fmt='(i2.2)') da
   write(unit=fhr, fmt='(i2.2)') hr

!- Form  NLDAS-1 File Name:
   if ( yr <= 2007 ) then
   !- GES DISC NLDAS-1 Naming Convention:
      if (trim(nldas1_struc(n)%nldas1_filesrc).eq."GES-DISC") then
         filename = trim(nldas1dir)//"/"//fyr//"/"//fdoy//&
                    "/NLDAS_FOR0125_H.A"//&
                    fyr//fmo//fda//"."//fhr//"00.001.grb"
   !- Original NLDAS-1 Naming Convention:
      elseif (trim(nldas1_struc(n)%nldas1_filesrc).eq."NCEP") then
         filename = trim(nldas1dir)//"/"//fyr//"/"//fyr//fmo//fda//&
                    "/"//fyr//fmo//fda//fhr//".FORCING.GRB"
      endif
!- NCEP NLDAS-1 filenames:
!  * (This solution is temporarily provided as files at the GES DISC are updated):
   elseif( yr >= 2008 ) then
      if (trim(nldas1_struc(n)%nldas1_filesrc).eq."NCEP") then
         filename = trim(nldas1dir)//"/"//fyr//"/"//fyr//fmo//fda//&
              "/"//fyr//fmo//fda//fhr//".FORCING.GRB"
      else
         filename = trim(nldas1dir)//"/"//fyr//"/"//fyr//fmo//fda//&
              "/"//fyr//fmo//fda//fhr//".nldasforce.grb"
      end if
   endif

   return

 end subroutine get_nldas1_filename
