!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
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
  use LDT_coreMod, only : LDT_rc
  use nldas1_forcingMod, only : nldas1_struc


! !DESCRIPTION:
!
!  This routine reads the options specific to NLDAS1 forcing from 
!  the LDT configuration file. 
!  
!EOP

  implicit none
  integer :: n

  do n=1,LDT_rc%nnest
     nldas1_struc(n)%nldas1time1 = 3000.0
     nldas1_struc(n)%nldas1time2 = 0.0
  enddo

end subroutine reset_nldas1
