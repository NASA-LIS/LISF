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
! !ROUTINE: AWRAL600_finalize
! \label{AWRAL600_finalize}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   12/18/18: Wendy Sharples, Shugong Wang; initial implementation for AWRAL600 with LIS-7
!
! !INTERFACE:
subroutine AWRAL600_finalize()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use AWRAL600_lsmMod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in AWRAL600
!
!EOP
    implicit none   
    
    integer :: t, n 

    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            deallocate(AWRAL600_struc(n)%awral600(t)%s0)
            deallocate(AWRAL600_struc(n)%awral600(t)%ss)
            deallocate(AWRAL600_struc(n)%awral600(t)%sd)
            deallocate(AWRAL600_struc(n)%awral600(t)%mleaf)
            deallocate(AWRAL600_struc(n)%awral600(t)%fhru)
            deallocate(AWRAL600_struc(n)%awral600(t)%hveg)
            deallocate(AWRAL600_struc(n)%awral600(t)%laimax)
            deallocate(AWRAL600_struc(n)%awral600(t)%height)
        end do  ! tile loop
 
        ! free memory for awral600, the data at tile level
        deallocate(AWRAL600_struc(n)%awral600)

        ! free momory for constant parameter
        deallocate(AWRAL600_struc(n)%hypsperc) 
        deallocate(AWRAL600_struc(n)%alb_dry)
        deallocate(AWRAL600_struc(n)%alb_wet)
        deallocate(AWRAL600_struc(n)%cgsmax)
        deallocate(AWRAL600_struc(n)%er_frac_ref)
        deallocate(AWRAL600_struc(n)%fsoilemax)
        deallocate(AWRAL600_struc(n)%lairef)
        deallocate(AWRAL600_struc(n)%rd)
        deallocate(AWRAL600_struc(n)%s_sls)
        deallocate(AWRAL600_struc(n)%sla)
        deallocate(AWRAL600_struc(n)%tgrow)
        deallocate(AWRAL600_struc(n)%tsenc)
        deallocate(AWRAL600_struc(n)%ud0)
        deallocate(AWRAL600_struc(n)%us0)
        deallocate(AWRAL600_struc(n)%vc)
        deallocate(AWRAL600_struc(n)%w0lime)
        deallocate(AWRAL600_struc(n)%w0ref_alb)
        deallocate(AWRAL600_struc(n)%wdlimu)
        deallocate(AWRAL600_struc(n)%wslimu)

        ! free momory for initial state variable
        if ( allocated(AWRAL600_struc(n)%init_s0)) then
          deallocate(AWRAL600_struc(n)%init_s0)
        end if
        if ( allocated(AWRAL600_struc(n)%init_ss)) then
          deallocate(AWRAL600_struc(n)%init_ss)
        end if
        if ( allocated(AWRAL600_struc(n)%init_s0)) then
          deallocate(AWRAL600_struc(n)%init_sd)
        end if
        if ( allocated(AWRAL600_struc(n)%init_s0)) then
          deallocate(AWRAL600_struc(n)%init_mleaf)
        end if
    end do ! nest loop
  
    deallocate(AWRAL600_struc)
 
end subroutine AWRAL600_finalize

