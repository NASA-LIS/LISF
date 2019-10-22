!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
  subroutine LDT_init_LISHydropreproc  

    use LDT_domainMod
    use LDT_paramProcMod
    use LDT_logMod
    use LDT_metforcingMod
    integer :: flag

    write(LDT_logunit,*) "----------------------------------"
    write(LDT_logunit,*) " Start of LDT parameter processing "
    write(LDT_logunit,*) "----------------------------------"
    
    call LDT_setDomainSpecs

    call LDT_paramProcConfig

    call LDT_paramProcInit(flag)          

    call LDT_flush(LDT_logunit)

  end subroutine LDT_init_LISHydropreproc
  
