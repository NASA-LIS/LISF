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
! !ROUTINE: mos_coldstart
!  \label{mos_coldstart}
!
! !REVISION HISTORY:
! 3 Jun  2003: Jon Gottschalck; Initial code
! 25 Sep 2007: Sujay Kumar; Upgraded for LIS5.0
! 
! !INTERFACE:
subroutine mos_coldstart(n)
! !USES:
   use LIS_coreMod, only : LIS_rc
   use LIS_logMod, only : LIS_logunit
   use LIS_timeMgrMod, only : LIS_date2time
   use mos_lsmMod
!
! !DESCRIPTION:
!  
!  This routine initializes the mosaic state variables with some 
!  predefined values uniformly for the entire domain. These initial values
!  will be overwritten by the values read from the supplied 
!  Mosaic model restart file. 
! 
!EOP
   implicit none
   integer, intent(in) :: n 
   integer :: t,l

   if ( LIS_rc%startcode .eq. "coldstart" ) then
      write(LIS_logunit,*)'MSG: mos_coldstart -- cold-starting mosaic', &
           '...using ics from card file'
      
      write(LIS_logunit,*)'MSG: mos_coldstart -- ntiles',LIS_rc%ntiles
      
      do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         mos_struc(n)%mos(t)%ct=mos_struc(n)%mos_it
         mos_struc(n)%mos(t)%qa=0.0
         mos_struc(n)%mos(t)%ics=0.0
         mos_struc(n)%mos(t)%snow=0.0
         mos_struc(n)%mos(t)%SoT=mos_struc(n)%mos_it
         do l=1,3
            mos_struc(n)%mos(t)%SoWET(l)=mos_struc(n)%mos_ism
         enddo
      enddo
      LIS_rc%yr=LIS_rc%syr
      LIS_rc%mo=LIS_rc%smo 
      LIS_rc%da=LIS_rc%sda
      LIS_rc%hr=LIS_rc%shr
      LIS_rc%mn=LIS_rc%smn
      LIS_rc%ss=LIS_rc%sss
      
      call LIS_date2time(LIS_rc%time,LIS_rc%doy,LIS_rc%gmt,LIS_rc%yr,&
           LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn,LIS_rc%ss) 
      write(LIS_logunit,*)'MSG: mos_coldstart -- Using the specified start time ',&
           LIS_rc%time
   endif
 end subroutine mos_coldstart
