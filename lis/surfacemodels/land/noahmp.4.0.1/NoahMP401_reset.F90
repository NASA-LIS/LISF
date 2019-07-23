!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: NoahMP401_reset
! \label{NoahMP401_reset}
!
! !REVISION HISTORY:
! Modified by Shugong Wang for Noah-MP.4.0.1  
! !INTERFACE:
subroutine NoahMP401_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use LIS_logMod,     only : LIS_logunit
  use NoahMP401_lsmMod

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
     write(LIS_logunit,*) "Noah-MP.4.0.1 resetting"

     ! initialize forcing variables to zeros
     do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         NOAHMP401_struc(n)%noahmp401(t)%lwdown = 0.0
         NOAHMP401_struc(n)%noahmp401(t)%swdown = 0.0
         NOAHMP401_struc(n)%noahmp401(t)%psurf  = 0.0
         NOAHMP401_struc(n)%noahmp401(t)%prcp   = 0.0
         NOAHMP401_struc(n)%noahmp401(t)%tair   = 0.0
         NOAHMP401_struc(n)%noahmp401(t)%qair   = 0.0
         NOAHMP401_struc(n)%noahmp401(t)%wind_e = 0.0
         NOAHMP401_struc(n)%noahmp401(t)%wind_n = 0.0
     enddo ! end of tile (t) loop
     NOAHMP401_struc(n)%forc_count = 0
     
  enddo ! do n=1,LIS_rc%nnest
end subroutine NoahMP401_reset
