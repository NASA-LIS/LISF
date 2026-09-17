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
! !ROUTINE: NoahMP50_finalize
! \label{NoahMP50_finalize}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later

! !INTERFACE:
subroutine NoahMP50_finalize()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use NoahMP50_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in NoahMP50
!
!EOP
    implicit none   
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%shdfac_monthly)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%soilcomp)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%smc)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%sh2o)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%tslb)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%tsno)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%zss)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%snowice)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%snowliq)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%smoiseq)
            deallocate(NoahMP50_struc(n)%noahmp50(t)%accetrani)
        end do  ! tile loop
 
        ! free memory for noahmp50, the data at tile level
        deallocate(NoahMP50_struc(n)%noahmp50)

        ! free momory for constant parameter 
        deallocate(NoahMP50_struc(n)%sldpth)

        ! free momory for initial state variable
        deallocate(NoahMP50_struc(n)%init_smc)
        deallocate(NoahMP50_struc(n)%init_tslb)
        
    end do ! nest loop
  
    deallocate(NoahMP50_struc)
 
end subroutine NoahMP50_finalize

