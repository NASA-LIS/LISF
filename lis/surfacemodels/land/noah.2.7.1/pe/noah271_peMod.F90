!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module noah271_peMod
!BOP
!
! !MODULE: noah271_peMod
!
! !DESCRIPTION:
!  This module contains the definitions of the Noah model parameters
!  used in parameter estimation. The data structure is used to expose
!  the LSM parameters to be used in opt/ue. 
!
! !REVISION HISTORY:
!  12 Jan 2012; Sujay Kumar, Initial Code
! !USES:        
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noah271_setup_pedecvars
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noah271_pe_struc

!EOP
  type, public ::  noah271_pe_dec 
     integer               :: nparams
     character*40, allocatable :: param_name(:)
     integer     , allocatable :: param_select(:)
     real        , allocatable :: param_min(:)
     real        , allocatable :: param_max(:)
  end type noah271_pe_dec

  type(noah271_pe_dec), allocatable :: noah271_pe_struc(:)

  SAVE
contains

!BOP
! !ROUTINE: noah271_setup_pedecvars
!  \label{noah271_setup_pedecvars}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine noah271_setup_pedecvars(DEC_State, Feas_State)
! !USES:
    use LIS_coreMod,       only : LIS_rc, LIS_config,LIS_vecTile
    use LIS_logMod,        only : LIS_logunit, LIS_verify
    use noah271_lsmMod,    only : noah271_struc

    implicit none
! !ARGUMENTS: 
    character(len=LIS_CONST_PATH_LEN) :: decSpaceAttribsFile
    type(ESMF_State)            :: DEC_State
    type(ESMF_State)            :: Feas_State
!
! !DESCRIPTION:
!  
!  This routine determines the list of parameters to be used in parameter
!  estimation, initializes them, and updates the LIS decision space.  
! 
!EOP

    integer                     :: n 
    type(ESMF_ArraySpec)        :: arrspec1
    type(ESMF_Field)            :: varField
    type(ESMF_Field)            :: feasField
    real, allocatable               :: vardata(:)
    integer, allocatable            :: mod_flag(:)
    integer                     :: i,t
    integer                     :: status
    

    call ESMF_StateGet(Feas_State, "Feasibility Flag", feasField, rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(feasField,localDE=0,farrayPtr=mod_flag,rc=status)
    call LIS_verify(status)
    
    mod_flag = 0

    call ESMF_ConfigGetAttribute(LIS_config,decSpaceAttribsFile,&
         label="LSM Decision space attributes files:",rc=status)
    call LIS_verify(status, "LSM Decision space attributes file: not defined")

    allocate(noah271_pe_struc(LIS_rc%nnest))
    n = 1
    noah271_pe_struc(n)%nparams = 25

    allocate(noah271_pe_struc(n)%param_name(noah271_pe_struc(n)%nparams))
    allocate(noah271_pe_struc(n)%param_select(noah271_pe_struc(n)%nparams))
    allocate(noah271_pe_struc(n)%param_min(noah271_pe_struc(n)%nparams))
    allocate(noah271_pe_struc(n)%param_max(noah271_pe_struc(n)%nparams))

    ! read the attributes file. 
    call LIS_readPEDecSpaceAttributes(decSpaceAttribsFile, &
         noah271_pe_struc(n)%nparams, &
         noah271_pe_struc(n)%param_name, &
         noah271_pe_struc(n)%param_select, &
         noah271_pe_struc(n)%param_min, &
         noah271_pe_struc(n)%param_max)

    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    do i=1,noah271_pe_struc(n)%nparams
       if(noah271_pe_struc(n)%param_select(i).eq.1) then 
          varField = ESMF_FieldCreate(arrayspec=arrspec1, &
               grid=LIS_vecTile(n),&
               name=trim(noah271_pe_struc(n)%param_name(i)),&
               rc=status)
          call LIS_verify(status, &
               'problem with fieldcreate in noah271_setup_pedecvars')
          call ESMF_AttributeSet(varField,'MinRange',&
               noah271_pe_struc(n)%param_min(i),rc=status)
          call LIS_verify(status, &
               'setting minrange to decspace obj in noah271_setup_devars')
          call ESMF_AttributeSet(varField,'MaxRange',&
               noah271_pe_struc(n)%param_max(i),rc=status)
          call LIS_verify(status, &
               'setting maxrange to decspace obj in noah271_setup_devars')

          call ESMF_StateAdd(DEC_State,(/varField/),rc=status)
          call LIS_verify(status,&
               'stateadd in noah271_setup_pedecvars')

       endif
    enddo

!initialize the decision space:
         
    do i=1,noah271_pe_struc(n)%nparams

       if(noah271_pe_struc(n)%param_select(i).eq.1) then 
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."SMCMAX")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%smcmax,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"SMCMAX",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%smcmax
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."PSISAT")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%psisat,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"PSISAT",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%psisat
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."DKSAT")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%dksat,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"DKSAT",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%dksat
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."DWSAT")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%dwsat,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"DWSAT",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%dwsat
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."BEXP")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%bexp,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"BEXP",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%bexp
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."QUARTZ")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%quartz,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"QUARTZ",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%quartz
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."RSMIN")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%rsmin,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"RSMIN",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%rsmin
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."RGL")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%rgl,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"RGL",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%rgl
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."HS")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%hs,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"HS",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%hs
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."Z0")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%z0,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"Z0",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%z0
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."LAI")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%lai,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"LAI",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%lai
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."CFACTR")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%cfactr,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"CFACTR",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%cfactr
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."CMCMAX")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%cmcmax,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"CMCMAX",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%cmcmax
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."SBETA")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%sbeta,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"SBETA",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%sbeta
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."RSMAX")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%rsmax,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"RSMAX",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%rsmax
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."TOPT")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%topt,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"TOPT",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%topt
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."REFDK")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%refdk,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"REFDK",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%refdk
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."FXEXP")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%fxexp,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"FXEXP",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%fxexp
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."REFKDT")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%refkdt,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"REFKDT",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%refkdt
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."CZIL")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%czil,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"CZIL",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%czil
             enddo

          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."CSOIL")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%csoil,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"CSOIL",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%csoil
             enddo
          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."FRZK")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%frzk,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"FRZK",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%frzk
             enddo
          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."SNUP")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%snup,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"SNUP",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%snup
             enddo
          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."SMCREF")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%smcref,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"SMCREF",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%smcref
             enddo
          endif
          if((trim(noah271_pe_struc(n)%param_name(i)).eq."SMCWLT")) then 
             call initializeDecSpace(n,&
                  noah271_struc(n)%noah(:)%smcwlt,&
                  noah271_pe_struc(n)%param_min(i),&
                  noah271_pe_struc(n)%param_max(i))

             call ESMF_StateGet(DEC_State,"SMCWLT",varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                 vardata(t) = noah271_struc(n)%noah(t)%smcwlt
             enddo
          endif
       endif
    enddo

    write(LIS_logunit,*) 'Finished setting up Noah 2.7.1 decision space '
  end subroutine noah271_setup_pedecvars

!BOP
! 
! !ROUTINE: initializeDecSpace
! \label{initializeDecSpace}
!
! !INTERFACE: 
  subroutine initializeDecSpace(n,var,varmin,varmax)
! !USES:
    use LIS_coreMod, only : LIS_rc
    use LIS_numerRecipesMod, only : LIS_rand_func
!
! !DESCRIPTION: 
!    This subroutine initialize the decision space variable. If the
!    decision space initialization mode is zero, then the population 
!    is initialized with the default LSM parameters. Otherwise the
!    population is initialized with randomly sampled values based
!    on the input range. 
! 
!EOP

    integer :: n 
    real    :: var(LIS_rc%ntiles(n))
    real    :: varmin
    real    :: varmax
    integer :: t,m
    integer :: seed
    real    :: rand

    if(LIS_rc%decSpaceInitMode.eq.1) then 
       seed = -1000
       do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
          call LIS_rand_func(seed,rand)
          
!keep the default in the first element
          do m=2,LIS_rc%nensem(n)
             call LIS_rand_func(1,rand)
             var((t-1)*LIS_rc%nensem(n)+m) = varmin + (varmax-varmin)*rand
          enddo
       enddo
    endif

  end subroutine initializeDecSpace


end module noah271_peMod
