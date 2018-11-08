!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: LVT_init_DAobs
! \label{LVT_init_DAobs}
!
! !INTERFACE: 
  subroutine LVT_init_DAobs
! 
! !USES: 
    use LVT_domainMod,       only : LVT_domainInit
    use LVT_LISModelDataMod, only : LVT_LISModelDataInit
    use LVT_observationsMod, only : LVT_obsDataInit
    use LVT_statsMod,        only : LVT_statsInit
    use LVT_logMod
    
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  
!   This routine defines the initialization phase of the runmode to 
!   analyse the output from a LIS data assimilation integration. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    call LVT_domainInit()
    call LVT_LISModelDataInit()
    call LVT_obsDataInit()
    call LVT_statsInit()
    call LVT_flush(LVT_logunit)

  end subroutine LVT_init_DAobs
