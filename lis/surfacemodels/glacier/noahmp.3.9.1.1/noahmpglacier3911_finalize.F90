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
! !ROUTINE: noahmpglacier3911_finalize
! \label{noahmpglacier3911_finalize}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
!
! !INTERFACE:
subroutine noahmpglacier3911_finalize(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
    use noahmpglacier3911_Mod
!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in noahmpglacier3911
!
!EOP
    implicit none   
    
    integer :: t, n 
#if 0 
    
    do n=1, LIS_rc%nnest
        ! free memory allocated for each tile
       do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
          deallocate(Noahmpgl3911_struc(n)%noahmpgl(t)%shdfac_monthly)
          deallocate(Noahmpgl3911_struc(n)%noahmpgl(t)%smceq)
          deallocate(Noahmpgl3911_struc(n)%noahmpgl(t)%sstc)
          deallocate(Noahmpgl3911_struc(n)%noahmpgl(t)%sh2o)
          deallocate(Noahmpgl3911_struc(n)%noahmpgl(t)%smc)
          deallocate(Noahmpgl3911_struc(n)%noahmpgl(t)%zss)
          deallocate(Noahmpgl3911_struc(n)%noahmpgl(t)%snowice)
          deallocate(Noahmpgl3911_struc(n)%noahmpgl(t)%snowliq)
       end do  ! tile loop
       
        ! free memory for noahmpgl, the data at tile level
       deallocate(Noahmpgl3911_struc(n)%noahmpgl)
       
       ! free momory for constant parameter 
       deallocate(Noahmpgl3911_struc(n)%sldpth)
       
        ! free momory for initial state variable
       deallocate(Noahmpgl3911_struc(n)%init_stc)
       deallocate(Noahmpgl3911_struc(n)%init_sh2o)
       deallocate(Noahmpgl3911_struc(n)%init_smc)
       !deallocate(Noahmpgl3911_struc(n)%init_zss)
       !deallocate(Noahmpgl3911_struc(n)%init_snowice)
       !deallocate(Noahmpgl3911_struc(n)%init_snowliq)
    end do ! nest loop
    
    deallocate(Noahmpgl3911_struc)
#endif
 
  end subroutine noahmpglacier3911_finalize

