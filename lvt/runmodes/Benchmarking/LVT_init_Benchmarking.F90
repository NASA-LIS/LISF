!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: LVT_init_Benchmarking
! \label{LVT_init_Benchmarking}
!
! !INTERFACE:
  subroutine LVT_init_Benchmarking
! 
! !USES: 
    use LVT_coreMod
    use LVT_domainMod
    use LVT_histDataMod
    use LVT_DataStreamsMod
    use LVT_statsMod
    use LVT_logMod
    use LVT_trainingMod
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
    call LVT_trainingInit()
    call LVT_checkTavgSpecs()
    call LVT_flush(LVT_logunit)

  end subroutine LVT_init_Benchmarking
  
