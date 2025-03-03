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
! !ROUTINE: RUC37_finalize
! \label{RUC37_finalize}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   1/15/15: Shugong Wang; initial implementation for RUC37 with LIS-7
!
! !INTERFACE:
subroutine RUC37_finalize()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use RUC37_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in RUC37
!
!EOP
    implicit none   
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(RUC37_struc(n)%ruc37(t)%albedo_monthly)
            deallocate(RUC37_struc(n)%ruc37(t)%shdfac_monthly)
            deallocate(RUC37_struc(n)%ruc37(t)%z0brd_monthly)
            deallocate(RUC37_struc(n)%ruc37(t)%lai_monthly)
            deallocate(RUC37_struc(n)%ruc37(t)%smc)
            deallocate(RUC37_struc(n)%ruc37(t)%sho)
            deallocate(RUC37_struc(n)%ruc37(t)%stc)
        end do  ! tile loop
 
        ! free memory for ruc37, the data at tile level
        deallocate(RUC37_struc(n)%ruc37)

        ! free momory for constant parameter 
        deallocate(RUC37_struc(n)%soil_layer_thickness)

        ! free momory for initial state variable
        deallocate(RUC37_struc(n)%init_smc)
        deallocate(RUC37_struc(n)%init_sho)
        deallocate(RUC37_struc(n)%init_stc)
    end do ! nest loop
  
    deallocate(RUC37_struc)
 
end subroutine RUC37_finalize

