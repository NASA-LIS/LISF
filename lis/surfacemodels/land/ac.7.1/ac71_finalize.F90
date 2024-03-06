!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: Ac71_finalize
! \label{Ac71_finalize}
!
! !REVISION HISTORY:
!  18 JAN 2024, Louise Busschaert; initial implementation for AC71
!
! !INTERFACE:
subroutine Ac71_finalize(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
    use Ac71_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in Ac71
!
!EOP
    implicit none   
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
    !TODO deallocate all vars
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(AC71_struc(n)%ac71(t)%smc)
        end do  ! tile loop
 
        ! free memory for ac71, the data at tile level
        deallocate(AC71_struc(n)%ac71)

        ! free memory for constant parameter 
        deallocate(AC71_struc(n)%Thickness)

        ! free memory for initial state variable
        deallocate(AC71_struc(n)%init_smc)
    end do ! nest loop
  
    deallocate(AC71_struc)
 
end subroutine Ac71_finalize

