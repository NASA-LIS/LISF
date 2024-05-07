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
! !ROUTINE: Crocus81_finalize
! \label{Crocus81_finalize}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   10/18/19: Mahdi Navari, Shugong Wang; initial implementation for Crocus81 with LIS-7
!
! !INTERFACE:
subroutine Crocus81_finalize()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use Crocus81_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in Crocus81
!
!EOP
    implicit none   
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(CROCUS81_struc(n)%crocus81(t)%SNOWSWE)
            deallocate(CROCUS81_struc(n)%crocus81(t)%SNOWRHO)
            deallocate(CROCUS81_struc(n)%crocus81(t)%SNOWHEAT)
            deallocate(CROCUS81_struc(n)%crocus81(t)%SNOWGRAN1)
            deallocate(CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2)
            deallocate(CROCUS81_struc(n)%crocus81(t)%SNOWHIST)
            deallocate(CROCUS81_struc(n)%crocus81(t)%SNOWAGE)
            deallocate(CROCUS81_struc(n)%crocus81(t)%ALB)
            deallocate(CROCUS81_struc(n)%crocus81(t)%SNOWLIQ)
            deallocate(CROCUS81_struc(n)%crocus81(t)%SNOWTEMP)
            deallocate(CROCUS81_struc(n)%crocus81(t)%SNOWDZ)
        end do  ! tile loop
 
        ! free memory for crocus81, the data at tile level
        deallocate(CROCUS81_struc(n)%crocus81)

        ! free momory for constant parameter 
        deallocate(CROCUS81_struc(n)%IMPWET)
        deallocate(CROCUS81_struc(n)%IMPDRY)

        ! free momory for initial state variable
        deallocate(CROCUS81_struc(n)%init_SNOWSWE)
        deallocate(CROCUS81_struc(n)%init_SNOWRHO)
        deallocate(CROCUS81_struc(n)%init_SNOWHEAT)
        deallocate(CROCUS81_struc(n)%init_SNOWGRAN1)
        deallocate(CROCUS81_struc(n)%init_SNOWGRAN2)
        deallocate(CROCUS81_struc(n)%init_SNOWHIST)
        deallocate(CROCUS81_struc(n)%init_SNOWAGE)
        deallocate(CROCUS81_struc(n)%init_SNOWLIQ)
        deallocate(CROCUS81_struc(n)%init_SNOWTEMP)
        deallocate(CROCUS81_struc(n)%init_SNOWDZ)
    end do ! nest loop
  
    deallocate(CROCUS81_struc)
 
end subroutine Crocus81_finalize

