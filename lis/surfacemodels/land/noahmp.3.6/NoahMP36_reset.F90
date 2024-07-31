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
! !ROUTINE: NoahMP36_reset
! \label{NoahMP36_reset}
!
! !REVISION HISTORY:
!  22 Feb 2018: Soni Yatheendradas; Initial version
! 
! !INTERFACE:
subroutine NoahMP36_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use NoahMP36_lsmMod
  use LIS_logMod,       only : LIS_verify, LIS_logunit

!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for NoahMP3.6 LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in NoahMP3.6.
!  
!EOP
  implicit none
  integer                 :: tt,n
  integer                 :: status


  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)                        &
          'NoahMP3.6 resetting'

     ! initialize forcing variables to zeros
     do tt=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         NOAHMP36_struc(n)%noahmp36(tt)%lwdown = 0.0
         NOAHMP36_struc(n)%noahmp36(tt)%swdown = 0.0
         NOAHMP36_struc(n)%noahmp36(tt)%psurf = 0.0
         NOAHMP36_struc(n)%noahmp36(tt)%prcp = 0.0
         NOAHMP36_struc(n)%noahmp36(tt)%tair = 0.0
         NOAHMP36_struc(n)%noahmp36(tt)%qair = 0.0
         NOAHMP36_struc(n)%noahmp36(tt)%wind_e = 0.0
         NOAHMP36_struc(n)%noahmp36(tt)%wind_n = 0.0
     enddo ! end of tile (tt) loop
     NOAHMP36_struc(n)%forc_count = 0
     
  enddo ! do n=1,LIS_rc%nnest
end subroutine NoahMP36_reset
