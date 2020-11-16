!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
  subroutine LDT_init_climoRstproc
    
    use LDT_coreMod
    use LDT_domainMod
    use LDT_paramProcMod
    use LDT_climorstProcMod
    use LDT_logMod

    integer :: n 

    write(LDT_logunit,*) "-------------------------------------------------"
    write(LDT_logunit,*) " Start of climatological restart file processing "
    write(LDT_logunit,*) "-------------------------------------------------"

    call LDT_setDomainSpecs
    call LDT_paramProcConfig
    call LDT_paramProcInit
    call LDT_timeInit    
    call LDT_domainInit
    call LDT_climoRstProcInit
    call LDT_coreinit
    call LDT_flush(LDT_logunit)

  end subroutine LDT_init_climoRstproc
  
