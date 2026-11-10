!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: NoahMP52_reset
! \label{NoahMP52_reset}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; for refactored Noah-MP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine NoahMP52_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use LIS_logMod,     only : LIS_logunit
  use NoahMP52_lsmMod

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

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) "[INFO] Noah-MP.5.2 resetting"

     ! initialize forcing variables to zeros
     do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         NoahMP52_struc(n)%noahmp52(t)%lwdown = 0.0
         NoahMP52_struc(n)%noahmp52(t)%swdown = 0.0
         NoahMP52_struc(n)%noahmp52(t)%psurf  = 0.0
         NoahMP52_struc(n)%noahmp52(t)%prcp   = 0.0
         NoahMP52_struc(n)%noahmp52(t)%tair   = 0.0
         NoahMP52_struc(n)%noahmp52(t)%qair   = 0.0
         NoahMP52_struc(n)%noahmp52(t)%wind_e = 0.0
         NoahMP52_struc(n)%noahmp52(t)%wind_n = 0.0
     enddo ! end of tile (t) loop
     NoahMP52_struc(n)%forc_count = 0

  enddo ! do n=1,LIS_rc%nnest
end subroutine NoahMP52_reset
