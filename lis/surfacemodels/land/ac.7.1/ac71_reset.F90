!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: Ac71_reset
! \label{Ac71_reset}
!
! !REVISION HISTORY:
!  22 Feb 2018: Soni Yatheendradas; Initial version
!  LB: TO DO
! 
! !INTERFACE:
subroutine Ac71_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use Ac71_lsmMod
  use LIS_logMod,       only : LIS_verify, LIS_logunit

!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Ac71 LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in Ac71.
!  
!EOP
  implicit none
  integer                 :: tt,n
  integer                 :: status


  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)                        &
          'Ac71 resetting'

     ! initialize forcing variables to zeros
     do tt=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         AC71_struc(n)%ac71(tt)%lwdown = 0.0
         AC71_struc(n)%ac71(tt)%swdown = 0.0
         AC71_struc(n)%ac71(tt)%psurf = 0.0
         AC71_struc(n)%ac71(tt)%prcp = 0.0
         AC71_struc(n)%ac71(tt)%tair = 0.0
         AC71_struc(n)%ac71(tt)%qair = 0.0
         AC71_struc(n)%ac71(tt)%wind_e = 0.0
         AC71_struc(n)%ac71(tt)%wind_n = 0.0
     enddo ! end of tile (tt) loop
     AC71_struc(n)%forc_count = 0
     
  enddo ! do n=1,LIS_rc%nnest
end subroutine Ac71_reset
