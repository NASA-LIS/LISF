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
! !ROUTINE: get_pptEnsFcst_filename
! \label{get_pptEnsFcst_filename}
!
! !REVISION HISTORY:
!  27 Sep 2016: K. Arsenault; Initial Implementation
!
! !INTERFACE:
 subroutine get_pptEnsFcst_filename( fcsttype, fcstyr, fcstmo,&
                ensnum, yr, mo, & 
                directory, filename)

! !USES:
   use LIS_timeMgrMod

   implicit none
! !ARGUMENTS: 
   character*20,  intent(in)  :: fcsttype        ! Forecast file of origin
   integer,       intent(in)  :: fcstyr          ! Forecast year
   integer,       intent(in)  :: fcstmo          ! Forecast month - Need to convert to "3-letter month"
   integer,       intent(in)  :: ensnum          ! Forecast ensemble number
   integer,       intent(in)  :: yr, mo          ! Lead-time year, month
   character(len=*), intent(in)  :: directory    ! Dataset Directory
   character(len=*), intent(out) :: filename
!
! !DESCRIPTION:
!   This subroutine puts together ensemble forecast 
!    file names.
! 
!  The arguments are:
!  \begin{description}
!   \item[doy]
!     Integer-based day of year (DOY) input
!   \item[directory]
!     Name of the forcing directory
!   \item[filename]
!     Name of the time-stamped pptEnsFcst forcing file
!  \end{description}
!
!EOP
   character*4  :: fyr, lyr
   character*2  :: fmo, lmo
   character*3  :: fmo3
   character*2  :: fensnum

  !=== end variable definition =============================================

   write(unit=fyr, fmt='(i4.4)') fcstyr
   write(unit=fmo, fmt='(i2.2)') fcstmo
   write(unit=lyr, fmt='(i4.4)') yr
   write(unit=lmo, fmt='(i2.2)') mo
   if( ensnum < 10 ) then
     write(unit=fensnum, fmt='(i1)') ensnum
   else
     write(unit=fensnum, fmt='(i2.2)') ensnum
   endif

  !=== Assemble pptEnsFcst filename:

  select case( fcsttype )
  
    ! If forecast dataset of origin is:  GEOS-5, or CFSv2 ...
    case( "GEOS5" )
  
      !- LIS function to convert 2-digit month to 3-char month:
      call LIS_mon3char( fmo, fmo3 )

      ! Former directory structure:
      filename = trim(directory)//"/"//fyr//"/"//fmo3//"01/ens"//&
          trim(fensnum)//"/PRECTOT."//lyr//lmo//".nc4"

    case( "CFSv2" )

      call LIS_mon3char( fmo, fmo3 )

      ! New directory structure (as of Nov 30, 2022):
      filename = trim(directory)//"/"//fmo3//"01/"//fyr//"/ens"//&
          trim(fensnum)//"/PRECTOT."//lyr//lmo//".nc4"

     ! ** Will need to update later to accomodate additional start
     !     dates, like Feb 05 or Oct 03. "01" default for now.

    case default
      write(*,*) "[WARN] No other forecast datasets supported at this time "

    end select

 end subroutine get_pptEnsFcst_filename

