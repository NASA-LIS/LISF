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
! !ROUTINE: getvicforcing
!  \label{getvicforcing412}
!
! !INTERFACE:
subroutine getvicforcing(n, findex)
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN 
  use vic_forcingMod,     only : vicforcing_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
! This is the entry point for the routines that read VIC-processed forcing
! data.
!
! Note that VIC-processed forcing data are in lock-step with either VIC's model
! time-step (for energy-balance runs) or VIC's snow time-step (for
! water-balance runs).
!
! Note that VIC-processed forcing data are considered valid at
! the given time-stamp of the corresponding forcing data files.
!EOP

   integer :: yr1, mo1, da1, hr1
   integer :: ferror
   character(len=LIS_CONST_PATH_LEN) :: fname


   yr1 = LIS_rc%yr
   mo1 = LIS_rc%mo
   da1 = LIS_rc%da
   hr1 = LIS_rc%hr

   call get_vicforcing_filename(fname, vicforcing_struc(n)%vicdir, &
                                yr1, mo1, da1, hr1)
   call vic412_read_gridded_forcing_data(n, findex, fname, ferror)

end subroutine getvicforcing

!BOP
!
! !ROUTINE: get_vicforcing_filename
!  \label{get_vicforcing412_filename}
!
! !INTERFACE:
subroutine get_vicforcing_filename(filename, dir, year, month, day, hour)
! !USES:
  implicit none
! !ARGUMENTS: 
   character(len=*), intent(out) :: filename
   character(len=*), intent(in) :: dir
   integer, intent(in) :: year
   integer, intent(in) :: month
   integer, intent(in) :: day
   integer, intent(in) :: hour
!  
! !DESCRIPTION:
!  The arguments are: 
!  \begin{description}
!  \item[filename] the name of the forcing file to be read
!  \item[dir] directory containing the VIC-processed forcing data
!  \item[year] year
!  \item[month] month
!  \item[day] day
!  \item[hour] hour
!  \end{description}
!
!  This routine generates the filename of the VIC-processed forcing file
!  to be read.
!EOP

   character(len=4) :: cyear
   character(len=2) :: cmonth
   character(len=2) :: cday
   character(len=2) :: chour
   character(len=2) :: cmin
   character(len=2) :: csec


   write(cyear,'(i4.4)')  year
   write(cmonth,'(i2.2)') month
   write(cday,'(i2.2)')   day
   write(chour,'(i2.2)')  hour
   write(cmin,'(i2.2)')   0 !LIS_rc%mn
   write(csec,'(i2.2)')   0 !LIS_rc%ss

   filename = trim(dir)//"/"//cyear//"/"//                        &
              cyear//cmonth//cday//"/"//                          &
              cyear//cmonth//cday//chour//cmin//".gd4r"

end subroutine get_vicforcing_filename
