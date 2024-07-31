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
! !ROUTINE: NoahMP36_finalize
! \label{NoahMP36_finalize}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   9/4/14: Shugong Wang; initial implementation for NoahMP36 with LIS-7
!
! !INTERFACE:
subroutine NoahMP36_finalize()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use NoahMP36_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in NoahMP36
!
!EOP
    implicit none   
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(NOAHMP36_struc(n)%noahmp36(t)%shdfac_monthly)
            deallocate(NOAHMP36_struc(n)%noahmp36(t)%smceq)
            deallocate(NOAHMP36_struc(n)%noahmp36(t)%sstc)
            deallocate(NOAHMP36_struc(n)%noahmp36(t)%sh2o)
            deallocate(NOAHMP36_struc(n)%noahmp36(t)%smc)
            deallocate(NOAHMP36_struc(n)%noahmp36(t)%zss)
            deallocate(NOAHMP36_struc(n)%noahmp36(t)%snowice)
            deallocate(NOAHMP36_struc(n)%noahmp36(t)%snowliq)
        end do  ! tile loop
 
        ! free memory for noahmp36, the data at tile level
        deallocate(NOAHMP36_struc(n)%noahmp36)

        ! free momory for constant parameter 
        deallocate(NOAHMP36_struc(n)%sldpth)

        ! free momory for initial state variable
        deallocate(NOAHMP36_struc(n)%init_stc)
        deallocate(NOAHMP36_struc(n)%init_sh2o)
        deallocate(NOAHMP36_struc(n)%init_smc)
        !deallocate(NOAHMP36_struc(n)%init_zss)
        !deallocate(NOAHMP36_struc(n)%init_snowice)
        !deallocate(NOAHMP36_struc(n)%init_snowliq)
    end do ! nest loop
  
    deallocate(NOAHMP36_struc)
 
end subroutine NoahMP36_finalize

