!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
  subroutine LDT_init_rstTransformProc
    
    use LDT_coreMod
    use LDT_domainMod
    use LDT_paramProcMod
    use LDT_rstTransformProcMod
    use LDT_logMod

    integer :: n 

    write(LDT_logunit,*) "-------------------------------------------------"
    write(LDT_logunit,*) " Start of restart transformation  "
    write(LDT_logunit,*) "-------------------------------------------------"

    call LDT_setDomainSpecs
    call LDT_paramProcConfig
    call LDT_paramProcInit
    call LDT_domainInit
    call LDT_rstTransformProcInit
    call LDT_flush(LDT_logunit)

  end subroutine LDT_init_rstTransformProc
  
