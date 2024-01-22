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
! !ROUTINE: Ac71_finalize
! \label{Ac71_finalize}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!  18 JAN 2024, Louise Busschaert; initial implementation for LIS 7 and AC71
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
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(AC71_struc(n)%ac71(t)%shdfac_monthly)
            deallocate(AC71_struc(n)%ac71(t)%smceq)
            deallocate(AC71_struc(n)%ac71(t)%sstc)
            deallocate(AC71_struc(n)%ac71(t)%sh2o)
            deallocate(AC71_struc(n)%ac71(t)%smc)
            deallocate(AC71_struc(n)%ac71(t)%zss)
            deallocate(AC71_struc(n)%ac71(t)%snowice)
            deallocate(AC71_struc(n)%ac71(t)%snowliq)
        end do  ! tile loop
 
        ! free memory for ac71, the data at tile level
        deallocate(AC71_struc(n)%ac71)

        ! free momory for constant parameter 
        deallocate(AC71_struc(n)%sldpth)
        deallocate(AC71_struc(n)%Thickness)

        ! free momory for initial state variable
        deallocate(AC71_struc(n)%init_stc)
        deallocate(AC71_struc(n)%init_sh2o)
        deallocate(AC71_struc(n)%init_smc)
        !deallocate(AC71_struc(n)%init_zss)
        !deallocate(AC71_struc(n)%init_snowice)
        !deallocate(AC71_struc(n)%init_snowliq)
    end do ! nest loop
  
    deallocate(AC71_struc)
 
end subroutine Ac71_finalize

