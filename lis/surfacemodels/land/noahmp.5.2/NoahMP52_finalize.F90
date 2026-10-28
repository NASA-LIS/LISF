!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!BOP
!
! !ROUTINE: NoahMP52_finalize
! \label{NoahMP52_finalize}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine NoahMP52_finalize()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use NoahMP52_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in NoahMP52
!
!EOP
    implicit none

    integer :: t, n

    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%shdfac_monthly)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%soilcomp)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%smc)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%sh2o)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%tslb)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%tsno)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%zss)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%snowice)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%snowliq)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%smoiseq)
            deallocate(NoahMP52_struc(n)%noahmp52(t)%accetrani)
        end do  ! tile loop

        ! free memory for noahmp52, the data at tile level
        deallocate(NoahMP52_struc(n)%noahmp52)

        ! free momory for constant parameter
        deallocate(NoahMP52_struc(n)%sldpth)

        ! free momory for initial state variable
        deallocate(NoahMP52_struc(n)%init_smc)
        deallocate(NoahMP52_struc(n)%init_tslb)

    end do ! nest loop

    deallocate(NoahMP52_struc)

end subroutine NoahMP52_finalize

