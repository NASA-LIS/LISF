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
! !ROUTINE: LVT_init_LISpost
! \label{LVT_init_LISpost}
!
! !INTERFACE:
  subroutine LVT_init_LISpost
! 
! !USES: 
    use LVT_coreMod
    use LVT_domainMod
    use LVT_LISpostMod
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
    call LVT_initialize_LISpost()
    
    flush(LVT_logunit)

  end subroutine LVT_init_LISpost
  
