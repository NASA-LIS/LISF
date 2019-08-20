!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0     
!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: NoahMP401_finalize
! \label{NoahMP401_finalize}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   10/25/18: Shugong Wang, Zhuo Wang; initial implementation for NoahMP401 with LIS-7
!
! !INTERFACE:
subroutine NoahMP401_finalize(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
    use NoahMP401_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in NoahMP401
!
!EOP
    implicit none   
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%shdfac_monthly)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%soilcomp)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%smc)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%sh2o)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%tslb)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%tsno)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%zss)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%snowice)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%snowliq)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%smoiseq)
            deallocate(NOAHMP401_struc(n)%noahmp401(t)%gecros_state)
        end do  ! tile loop
 
        ! free memory for noahmp401, the data at tile level
        deallocate(NOAHMP401_struc(n)%noahmp401)

        ! free momory for constant parameter 
        deallocate(NOAHMP401_struc(n)%sldpth)

        ! free momory for initial state variable
        deallocate(NOAHMP401_struc(n)%init_smc)
        deallocate(NOAHMP401_struc(n)%init_tslb)
        deallocate(NOAHMP401_struc(n)%init_gecros_state)
        
    end do ! nest loop
  
    deallocate(NOAHMP401_struc)
 
end subroutine NoahMP401_finalize

