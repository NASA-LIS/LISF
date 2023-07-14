!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!BOP
!
! !ROUTINE: NoahMPnew_finalize
! \label{NoahMPnew_finalize}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   10/25/18: Shugong Wang, Zhuo Wang; initial implementation for NoahMP401 with LIS-7
!  May 2023: Cenlin He, modified for refactored NoahMP v5 and later

! !INTERFACE:
subroutine NoahMPnew_finalize()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use NoahMPnew_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in NoahMPnew
!
!EOP
    implicit none   
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(NoahMPnew_struc(n)%noahmpnew(t)%shdfac_monthly)
            deallocate(NoahMPnew_struc(n)%noahmpnew(t)%soilcomp)
            deallocate(NoahMPnew_struc(n)%noahmpnew(t)%smc)
            deallocate(NoahMPnew_struc(n)%noahmpnew(t)%sh2o)
            deallocate(NoahMPnew_struc(n)%noahmpnew(t)%tslb)
            deallocate(NoahMPnew_struc(n)%noahmpnew(t)%tsno)
            deallocate(NoahMPnew_struc(n)%noahmpnew(t)%zss)
            deallocate(NoahMPnew_struc(n)%noahmpnew(t)%snowice)
            deallocate(NoahMPnew_struc(n)%noahmpnew(t)%snowliq)
            deallocate(NoahMPnew_struc(n)%noahmpnew(t)%smoiseq)
        end do  ! tile loop
 
        ! free memory for noahmpnew, the data at tile level
        deallocate(NoahMPnew_struc(n)%noahmpnew)

        ! free momory for constant parameter 
        deallocate(NoahMPnew_struc(n)%sldpth)

        ! free momory for initial state variable
        deallocate(NoahMPnew_struc(n)%init_smc)
        deallocate(NoahMPnew_struc(n)%init_tslb)
        
    end do ! nest loop
  
    deallocate(NoahMPnew_struc)
 
end subroutine NoahMPnew_finalize

