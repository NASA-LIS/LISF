!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: LVT_init_optUE
!   label(LVT_init_optUE)
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
  subroutine LVT_init_optUE

    use LVT_domainMod,       only : LVT_domainInit
    use LVT_optUEMod,           only : LVT_optUEInit
    use LVT_logMod
    
    call LVT_domainInit()
    call LVT_optUEInit()
    call LVT_flush(LVT_logunit)

  end subroutine LVT_init_optUE
  
