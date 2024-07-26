!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
  subroutine LDT_init_StatDscaleMetForc
    
    use LDT_coreMod 
    use LDT_domainMod
    use LDT_paramProcMod
    use LDT_logMod
    use LDT_metforcingMod
    use LDT_statDscaleMod

    integer  :: n 

    write(LDT_logunit,*) "------------------------------------------------"
    write(LDT_logunit,*) " Start of statistical downscaling of met forcing "
    write(LDT_logunit,*) "------------------------------------------------"

    call LDT_timeInit
    call LDT_setDomainSpecs
    call LDT_paramProcConfig
    call LDT_readProcParamInit
    call LDT_domainInit
    call LDT_metforcingInit
    call LDT_statDscaleInit
    call LDT_coreInit

  end subroutine LDT_init_StatDscaleMetForc
  
