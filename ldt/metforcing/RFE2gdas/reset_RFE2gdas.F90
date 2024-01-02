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
! !MODULE: reset_RFE2gdas
! \label{reset_RFE2gdas}
! 
! !REVISION HISTORY: 
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 
! !INTERFACE:
subroutine reset_RFE2gdas
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use RFE2gdas_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for RFE2gdas forcing. 
!
!EOP  
  implicit none
  integer :: n 

  do n=1,LDT_rc%nnest
     RFE2gdas_struc(n)%RFE2gdasEndTime = 0.0
  enddo
end subroutine reset_RFE2gdas
