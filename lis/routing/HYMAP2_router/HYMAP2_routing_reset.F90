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
! !ROUTINE: HYMAP2_routing_reset
! \label{HYMAP2_routing_reset}
!
! !REVISION HISTORY:
!
! !INTERFACE:
subroutine HYMAP2_routing_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use LIS_logMod,     only : LIS_logunit

!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Noah-MP.4.0.1 LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in Noah-MP.4.0.1.
!  
!EOP
  implicit none
  integer                 :: t,n
  integer                 :: status


  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) "HYMAP2 resetting"

     
     
  enddo ! do n=1,LIS_rc%nnest
end subroutine HYMAP2_routing_reset
