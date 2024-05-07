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
! !ROUTINE: clm2_lairead.F90:
!
! !DESCRIPTION:
!  This program reads in AVHRR LAI data for CLM
!
! !REVISION HISTORY:
!  27 Nov 2001: Jon Gottschalck; Initial code
!  20 Feb 2002: Jon Gottschalck; Modified to use for 1/4 and 2x2.5 using 1/8 degree monthly data
!  01 Oct 2002: Jon Gottschalck; Modified to add MODIS LAI data
! 
! !INTERFACE: 
subroutine clm2_lairead (n)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_vegDataMod, only : LIS_lai, LIS_sai
  use clm2_lsmMod
!EOP
  implicit none
  integer, intent(in) :: n 
  integer :: t,tid

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
     clm2_struc(n)%clm(t)%tlai = LIS_lai(n)%tlai(tid)     
     clm2_struc(n)%clm(t)%tsai = LIS_sai(n)%tsai(tid)
!     if (clm2_struc(n)%clm(t)%itypveg .eq. 12 ) then 
     if (clm2_struc(n)%clm(t)%itypveg .eq. 0 ) then !noveg
        clm2_struc(n)%clm(t)%tlai=0.0
        clm2_struc(n)%clm(t)%tsai=0.0
        clm2_struc(n)%clm(t)%htop=0.0
        clm2_struc(n)%clm(t)%hbot=0.0  
     endif
!     if (clm2_struc(n)%clm(t)%itypveg .eq. 13 ) then 
!        clm2_struc(n)%clm(t)%tlai=0.0
!        clm2_struc(n)%clm(t)%tsai=0.0
!     endif
  enddo
  
end subroutine clm2_lairead
