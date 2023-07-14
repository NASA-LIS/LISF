!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: NoahMPnew_reset
! \label{NoahMPnew_reset}
!
! !REVISION HISTORY:
! Modified by Shugong Wang for Noah-MP.4.0.1  
! modified by Cenlin He for refactored Noah-MP v5 and later

! !INTERFACE:
subroutine NoahMPnew_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use LIS_logMod,     only : LIS_logunit
  use NoahMPnew_lsmMod

!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Noah-MP LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in Noah-MP.
!  
!EOP
  implicit none
  integer                 :: t,n
  integer                 :: status


  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) "Noah-MP.New resetting"

     ! initialize forcing variables to zeros
     do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         NoahMPnew_struc(n)%noahmpnew(t)%lwdown = 0.0
         NoahMPnew_struc(n)%noahmpnew(t)%swdown = 0.0
         NoahMPnew_struc(n)%noahmpnew(t)%psurf  = 0.0
         NoahMPnew_struc(n)%noahmpnew(t)%prcp   = 0.0
         NoahMPnew_struc(n)%noahmpnew(t)%tair   = 0.0
         NoahMPnew_struc(n)%noahmpnew(t)%qair   = 0.0
         NoahMPnew_struc(n)%noahmpnew(t)%wind_e = 0.0
         NoahMPnew_struc(n)%noahmpnew(t)%wind_n = 0.0
     enddo ! end of tile (t) loop
     NoahMPnew_struc(n)%forc_count = 0
     
  enddo ! do n=1,LIS_rc%nnest
end subroutine NoahMPnew_reset
