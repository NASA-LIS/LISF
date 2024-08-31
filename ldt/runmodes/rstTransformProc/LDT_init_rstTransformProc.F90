!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
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
    flush(LDT_logunit)

  end subroutine LDT_init_rstTransformProc
  
