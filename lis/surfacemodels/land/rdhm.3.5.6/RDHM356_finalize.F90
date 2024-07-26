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
! !ROUTINE: RDHM356_finalize
! \label{RDHM356_finalize}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   11/5/13: Shugong Wang; initial implementation for RDHM356 with LIS-7
!
! !INTERFACE:
subroutine RDHM356_finalize()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use RDHM356_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in RDHM356
!
!EOP
    implicit none
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(RDHM356_struc(n)%rdhm356(t)%SWINT)
            deallocate(RDHM356_struc(n)%rdhm356(t)%SWHINT)
            deallocate(RDHM356_struc(n)%rdhm356(t)%TSINT)
        end do  ! tile loop
 
        ! free memory for rdhm356, the data at tile level
        deallocate(RDHM356_struc(n)%rdhm356)

        ! free momory for constant parameter
        deallocate(RDHM356_struc(n)%DSINTW)
        deallocate(RDHM356_struc(n)%DSINT)

    end do ! nest loop
  
    deallocate(RDHM356_struc)
 
end subroutine RDHM356_finalize

