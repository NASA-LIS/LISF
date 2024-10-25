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
! !ROUTINE: agrradps_filename_sw
! \label{agrradps_filename_sw}
! 
! !INTERFACE: 
 subroutine agrradps_filename_sw(nameNH,nameSH,rootdir,yr,mo,da,hr)

   implicit none
! !ARGUMENTS:   
  
   character(len=*)          :: nameNH
   character(len=*)          :: nameSH
   character(len=*)          :: rootdir
   integer, intent(IN)       :: yr
   integer, intent(IN)       :: mo
   integer, intent(IN)       :: da
   integer, intent(IN)       :: hr
! 
! !DESCRIPTION: 
!  This routines generates the name of the shortwave data file
!  by appending the hemisphere and timestamps to the root directory. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[nameNH] filename to be generated for NH
!   \item[nameSH] filename to be generated for SH
!   \item[rootdir] name of the root directory that contains the forcing
!   \item[yr] year of data
!   \item[mo] month of  data
!   \item[da] day of data
!   \item[hr] hour of data
!  \end{description}
!
!EOP
  
   character(len=6)    :: ftime1
   character(len=10)   :: ftime2
   
   write(unit=ftime1,fmt='(i4,i2.2)') yr,mo
   write(unit=ftime2,fmt='(i4,i2.2,i2.2,i2.2)') yr,mo,da,hr
   
   nameNH = trim(rootdir)//'/SWDN/'//ftime1//'/NH/swdn_'//ftime2//'n'
   nameSH = trim(rootdir)//'/SWDN/'//ftime1//'/SH/swdn_'//ftime2//'s'
 end subroutine agrradps_filename_sw
 
!BOP
! 
! !ROUTINE: agrradps_filename_cloud
! \label{agrradps_filename_cloud}
! 
! !INTERFACE: 
 subroutine agrradps_filename_cloud(nameNH,nameSH,rootdir,yr,mo,da,hr,layer)

   implicit none
! !ARGUMENTS:   
  
   character(len=*)          :: nameNH
   character(len=*)          :: nameSH
   character(len=*)          :: rootdir
   integer, intent(IN)       :: yr
   integer, intent(IN)       :: mo
   integer, intent(IN)       :: da
   integer, intent(IN)       :: hr
   character(len=1), intent(IN)  :: layer
! 
! !DESCRIPTION: 
!  This routines generates the name of the cloud amount data file
!  by appending the hemisphere and timestamps to the root directory. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[name] filename to be generated
!   \item[rootdir] name of the root directory that contains the forcing
!   \item[yr] year of data
!   \item[mo] month of  data
!   \item[da] day of data
!   \item[hr] hour of data
!  \end{description}
!
!EOP
  
   character(len=6)    :: ftime1
   character(len=10)   :: ftime2
   
   write(unit=ftime1,fmt='(i4,i2.2)') yr,mo
   write(unit=ftime2,fmt='(i4,i2.2,i2.2,i2.2)') yr,mo,da,hr
   
   nameNH = trim(rootdir)//'/CloudAGR/'//ftime1//'/NH/cldamt' &
            //layer//'_'//ftime2//'n'
   nameSH = trim(rootdir)//'/CloudAGR/'//ftime1//'/SH/cldamt' &
            //layer//'_'//ftime2//'s'
 end subroutine agrradps_filename_cloud

