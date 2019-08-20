!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: LVT_run_DAstats
!   \label(LVT_run_DAstats)
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
subroutine LVT_run_DAstats

  use LVT_coreMod
  use LVT_DAMod
  use LVT_timeMgrMod
  use LVT_statsMod
  use LVT_domainMod
  use LVT_logMod

  implicit none
  integer :: i 

  do i=1,LVT_rc%pass
     call LVT_resetClock(LVT_rc)
     do while (.NOT. LVT_endofrun())
        call LVT_ticktime   
        call LVT_readDataMask
        call LVT_readLISDAData(i)
        call LVT_computeDAStats(i)
        call LVT_flush(LVT_logunit)
     enddo    
  enddo
end subroutine LVT_run_DAstats
  
