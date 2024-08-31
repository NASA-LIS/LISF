!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
  subroutine LDT_init_MetTimeDScale
    
    use LDT_coreMod 
    use LDT_domainMod
    use LDT_paramProcMod
    use LDT_logMod
    use LDT_metforcingMod
    use LDT_metforcDscaleMod
!    use LDT_perturbMod,   only : LDT_perturb_init

    integer  :: n 

    write(LDT_logunit,*) "------------------------------------------------"
    write(LDT_logunit,*) " Start of Metforce Temporal Downscale Processing "
    write(LDT_logunit,*) "------------------------------------------------"

!    print *, "Calling LDT_timeInit"
    call LDT_timeInit

!    print *, "Calling LDT_setDomainSpecs"
    call LDT_setDomainSpecs

!    print *, "Calling LDT_paramProcConfig"
    call LDT_paramProcConfig

!    print *, "Calling LDT_readProcParamInit"
    call LDT_readProcParamInit

!    print *, "Calling LDT_domainInit"
    call LDT_domainInit

!    print *, "Calling LDT MetforcingInit:"
    call LDT_metforcingInit

!    print *, "Calling LDT Metforce Dscale Init"
    call LDT_metforcDscaleInit

!    print *, "Calling LDT_core_init"
    call LDT_coreInit

  end subroutine LDT_init_MetTimeDScale
  
