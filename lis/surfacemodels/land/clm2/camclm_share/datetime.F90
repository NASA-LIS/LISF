!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! Purpose:
!
!  A generic Date and Time routine
!
! Author: CCM Core group
!
!-----------------------------------------------------------------------
!
! $Id: datetime.F90,v 1.5 2004/05/07 22:18:33 jim Exp $
!
!-----------------------------------------------------------------------
#include "LIS_misc.h"
   subroutine datetime(cdate, ctime) 
   implicit none
!-----------------------------------------------------------------------
!
!-----------------------------Arguments---------------------------------
   character , intent(out) :: cdate*8 
   character , intent(out) :: ctime*8 
!-----------------------------------------------------------------------
!
!---------------------------Local Variables------------------------------
   integer, dimension(8) :: values 
   character :: date*8, time*10, zone*5 
!-----------------------------------------------------------------------
 
   call date_and_time (date, time, zone, values) 
   cdate(1:2) = date(5:6) 
   cdate(3:3) = '/' 
   cdate(4:5) = date(7:8) 
   cdate(6:6) = '/' 
   cdate(7:8) = date(3:4) 
   ctime(1:2) = time(1:2) 
   ctime(3:3) = ':' 
   ctime(4:5) = time(3:4) 
   ctime(6:6) = ':' 
   ctime(7:8) = time(5:6) 
 
   return  
   end subroutine datetime 
 
