!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: reset_nldas1
!  \label{reset_nldas1}
!
! !REVISION HISTORY:
! 02Feb2004; Sujay Kumar, Initial Code
! 04Mar2007; Kristi Arsenault, Implemented EDAS Elevation Correction
!
! !INTERFACE:    
subroutine reset_nldas1()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use nldas1_forcingMod, only : nldas1_struc


! !DESCRIPTION:
!
!  This routine reads the options specific to NLDAS1 forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n

  do n=1,LIS_rc%nnest
     nldas1_struc(n)%nldas1time1 = 3000.0
     nldas1_struc(n)%nldas1time2 = 0.0
     
  enddo

end subroutine reset_nldas1
