!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: LVT_init_DAstats
!  \label(LVT_init_DAstats)
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
  subroutine LVT_init_DAstats

    use LVT_domainMod,       only : LVT_domainInit
    use LVT_DAMod,           only : LVT_DAInit
    use LVT_statsMod,        only : LVT_statsInit
    use LVT_logMod
    
    call LVT_domainInit()
    call LVT_statsInit()
    call LVT_DAInit()
    call LVT_flush(LVT_logunit)

  end subroutine LVT_init_DAstats
  
