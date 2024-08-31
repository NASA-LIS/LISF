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
! !ROUTINE: finalize_petusgs
! \label{finalize_petusgs}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 08Mar2012: Kristi Arsenault;  Modified for USGS PET dataset
! 25Oct2013: Kristi Arsenault;  Added PET USGS to LIS7
! 
! !INTERFACE:
subroutine finalize_petusgs(findex)

! !USES:
  use petusgs_forcingMod, only : petusgs_struc

! !DESCRIPTION: 
!  Routine to cleanup PET USGS forcing related memory allocations.   
!
!EOP
  implicit none
  integer, intent(IN) :: findex  

  deallocate( petusgs_struc )

end subroutine finalize_petusgs
