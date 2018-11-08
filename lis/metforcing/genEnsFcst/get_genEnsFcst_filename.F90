!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: get_genEnsFcst_filename
! \label{get_genEnsFcst_filename}
!
! !REVISION HISTORY:
!  27 Sep 2016: K. Arsenault; Initial Implementation
!
! !INTERFACE:
 subroutine get_genEnsFcst_filename( fcsttype, fcstyr, fcstmo,&
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
   character*100, intent(in)  :: directory       ! Dataset Directory
   character*140, intent(out) :: filename        
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
!     Name of the time-stamped genEnsFcst forcing file
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

  !=== Assemble genEnsFcst filename:

  select case( fcsttype )
  
    case( "GEOS5" )
  ! If forecast dataset of origin is:  GEOS-5 ...
  
     !- convert 2-digit month to 3-char month:
     ! LIS function somewhere to do that or do  somewhere here??
     ! 
!      call mon3char( fmo, fmo3 )
      call LIS_mon3char( fmo, fmo3 )

     ! What to do about the ensemble number ??  Call this routine
     !  for every individual member??  Or loop over and generate 
     !  a number of files to be passed back to main routine??

      filename = trim(directory)//"/"//fyr//"/"//fmo3//"01/ens"//&
          trim(fensnum)//"/GEOS5."//lyr//lmo//".nc4"
!          trim(fensnum)//"/GEOS5.all_forc_"//lyr//lmo//".nc4"
!            /6-hourly/1982/may01/ens1/GEOS5.all_forc_198205.nc4

     ! ** Will need to update later to accomodate additional start
     !     dates, like Feb 05 or Oct 03. "01" default for now.

    case default
      write(*,*) " No other forecast datasets supported at this time "

    end select

 end subroutine get_genEnsFcst_filename


#if 0
 subroutine mon3char( mo2digchar, mo3char )

    implicit none
    character*2, intent(in)  :: mo2digchar
    character*3, intent(out) :: mo3char

    select case( mo2digchar )
     case( "01" )
       mo3char = "jan"
     case( "02" )
       mo3char = "feb"
     case( "03" )
       mo3char = "mar"
     case( "04" )
       mo3char = "apr"
     case( "05" )
       mo3char = "may"
     case( "06" )
       mo3char = "jun"
     case( "07" )
       mo3char = "jul"
     case( "08" )
       mo3char = "aug"
     case( "09" )
       mo3char = "sep"
     case( "10" )
       mo3char = "oct"
     case( "11" )
       mo3char = "nov"
     case( "12" )
       mo3char = "dec"
     case default
       write(*,*) " Don't recognize this 2-digit month ... stopping. "
       stop
    end select

  end subroutine mon3char
#endif
