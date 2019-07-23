!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: LVT_run_optUE
!  label(LVT_run_optUE)
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
subroutine LVT_run_optUE

  use LVT_coreMod,         only : LVT_rc, LVT_endofrun, LVT_ticktime
  use LVT_optUEMod,           only : LVT_readoptUEData, LVT_computeoptUEstats, &
       LVT_optuectl
  use LVT_logMod

  implicit none
  integer :: i 

  do i=1,LVT_optuectl%maxIter
     call LVT_readoptUEData(i)
     call LVT_computeoptUEStats(i)
     call LVT_flush(LVT_logunit)
  enddo
end subroutine LVT_run_optUE
  
