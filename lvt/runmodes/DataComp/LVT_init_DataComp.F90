!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: LVT_init_DataComp
! \label{LVT_init_DataComp}
!
! !INTERFACE:
  subroutine LVT_init_DataComp
! 
! !USES: 
    use LVT_coreMod
    use LVT_domainMod
    use LVT_histDataMod
    use LVT_DataStreamsMod
    use LVT_statsMod
    use LVT_logMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine defines the initialization steps for performing LVT 
!  analysis using model output from LIS. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    call LVT_domainInit()
    call LVT_histDataInit()
    call LVT_DataStreamsinit()
    call LVT_checkDatastreamSetup()
    call LVT_statsInit()
    call LVT_checkTavgSpecs()
    call LVT_core_init()
    call LVT_flush(LVT_logunit)

  end subroutine LVT_init_DataComp
  
