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
    flush(LVT_logunit)

  end subroutine LVT_init_DAstats
  
