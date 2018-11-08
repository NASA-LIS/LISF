!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
  subroutine LDT_init_ANNproc
    
    use LDT_coreMod 
    use LDT_domainMod
    use LDT_paramProcMod
    use LDT_ANNMod
    use LDT_logMod

    integer  :: n 

    write(LDT_logunit,*) "----------------------------------"
    write(LDT_logunit,*) " Start of ANN processing "
    write(LDT_logunit,*) "----------------------------------"

    call LDT_setDomainSpecs
    call LDT_paramProcConfig
    call LDT_paramProcInit  
    call LDT_timeInit
    call LDT_domainInit
    call LDT_ANNinit()
    call LDT_flush(LDT_logunit)

  end subroutine LDT_init_ANNproc
  
