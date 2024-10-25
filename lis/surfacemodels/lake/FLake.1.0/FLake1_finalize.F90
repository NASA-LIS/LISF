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
! !ROUTINE: FLake1_finalize
! \label{FLake1_finalize}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   6/4/13: Shugong Wang; initial implementation for FLake1 with LIS-7
!
! !INTERFACE:
subroutine FLake1_finalize(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
    use FLake1_Mod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in FLake1
!
!EOP
    implicit none
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
        deallocate(FLAKE1_struc(n)%flake1)
    end do ! nest loop
  
    deallocate(FLAKE1_struc)
 
end subroutine FLake1_finalize

