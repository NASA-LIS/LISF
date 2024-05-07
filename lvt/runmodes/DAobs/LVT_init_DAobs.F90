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
    flush(LVT_logunit)

  end subroutine LVT_init_DAobs
