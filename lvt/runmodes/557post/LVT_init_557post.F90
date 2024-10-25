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
! !ROUTINE: LVT_init_557post
! \label{LVT_init_557post}
!
! !INTERFACE:
  subroutine LVT_init_557post
! 
! !USES: 
    use ESMF
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

    integer    :: rc

    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%processHYCOM,&
         label="Process HYCOM data:",rc=rc)  
    call LVT_verify(rc,'Process HYCOM data: not defined')
  
    if(LVT_rc%processHYCOM.eq.1) then 
       call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%HYCOMdir,&
            label="HYCOM data directory:",rc=rc)  
       call LVT_verify(rc,'HYCOM data directory: not defined')
    endif

    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%applyNoiseReductionFilter,&
         label="Apply noise reduction filter:",rc=rc)  
    call LVT_verify(rc,'Apply noise reduction filter: not defined')
    if(LVT_rc%applyNoiseReductionFilter.eq.1) then 
       call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%smoothingFilterType,&
            label="Smoothing filter type:",rc=rc)  
       call LVT_verify(rc,'Smoothing filter type: not defined')       
    endif

    call LVT_domainInit()
    call LVT_histDataInit()
    call LVT_DataStreamsinit()
    call LVT_checkDatastreamSetup()
    call LVT_statsInit()
    call LVT_checkTavgSpecs()
    call LVT_core_init()
    flush(LVT_logunit)
    
  end subroutine LVT_init_557post
  
