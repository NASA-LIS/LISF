!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module noah36_peMod
!BOP
!
! !MODULE: noah36_peMod
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
  use LIS_numerRecipesMod, only : LIS_rand_func
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noah36_setup_pedecvars
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noah36_pe_struc

!EOP
  type, public ::  noah36_pe_dec 
     integer               :: nparams
     character*40, allocatable :: param_name(:)
     integer     , allocatable :: param_select(:)
     real        , allocatable :: param_min(:)
     real        , allocatable :: param_max(:)
  end type noah36_pe_dec

  type(noah36_pe_dec), allocatable :: noah36_pe_struc(:)

  SAVE
contains

!BOP
! !ROUTINE: noah36_setup_pedecvars
!  \label{noah36_setup_pedecvars}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine noah36_setup_pedecvars(DEC_State, Feas_State)
! !USES:
    use LIS_coreMod,       only : LIS_rc, LIS_config,LIS_vecPatch, LIS_surface, LIS_localPET
    use LIS_logMod,        only : LIS_logunit, LIS_verify
    use noah36_lsmMod,     only : noah36_struc

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
    real, pointer               :: vardata(:)
    integer, pointer            :: mod_flag(:)
    integer                     :: i,t,m
    integer                     :: status
    integer, parameter          :: seed_base=-1000
    integer                     :: seed
    real                        :: rand
    integer                     :: NT    
    character*100               :: vname
    integer                     :: count
    integer                     :: gid

    call ESMF_StateGet(Feas_State, "Feasibility Flag", feasField, rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(feasField,localDE=0,farrayPtr=mod_flag,rc=status)
    call LIS_verify(status)
    
!    mod_flag = 0  !now initialized centrally in LIS_optUEMod.F90 as lsms,rtms,etc. can raise flag

    call ESMF_ConfigGetAttribute(LIS_config,decSpaceAttribsFile,&
         label="LSM Decision space attributes file:",rc=status)
    call LIS_verify(status, "LSM Decision space attributes file: not defined")

    allocate(noah36_pe_struc(LIS_rc%nnest))
    n = 1
    noah36_pe_struc(n)%nparams = 39

    allocate(noah36_pe_struc(n)%param_name(noah36_pe_struc(n)%nparams))
    allocate(noah36_pe_struc(n)%param_select(noah36_pe_struc(n)%nparams))
    allocate(noah36_pe_struc(n)%param_min(noah36_pe_struc(n)%nparams))
    allocate(noah36_pe_struc(n)%param_max(noah36_pe_struc(n)%nparams))

    ! read the attributes file. 
    call LIS_readPEDecSpaceAttributes(decSpaceAttribsFile, &
         noah36_pe_struc(n)%nparams, &
         noah36_pe_struc(n)%param_name, &
         noah36_pe_struc(n)%param_select, &
         noah36_pe_struc(n)%param_min, &
         noah36_pe_struc(n)%param_max)

    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    do i=1,noah36_pe_struc(n)%nparams
       if(noah36_pe_struc(n)%param_select(i).eq.1) then 
          vname=trim(noah36_pe_struc(n)%param_name(i))
          varField = ESMF_FieldCreate(arrayspec=arrspec1, &
               grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
               name=vname,&
               rc=status)
          call LIS_verify(status, &
               'problem with fieldcreate in noah36_setup_pedecvars')
          call ESMF_AttributeSet(varField,'MinRange',&
               noah36_pe_struc(n)%param_min(i),rc=status)
          call LIS_verify(status, &
               'setting minrange to decspace obj in noah36_setup_devars')
          call ESMF_AttributeSet(varField,'MaxRange',&
               noah36_pe_struc(n)%param_max(i),rc=status)
          call LIS_verify(status, &
               'setting maxrange to decspace obj in noah36_setup_devars')

          call ESMF_StateAdd(DEC_State,(/varField/),rc=status)
          call LIS_verify(status,&
               'stateadd in noah36_setup_pedecvars')

          call ESMF_StateGet(DEC_State,vname,varField,rc=status)
          call LIS_verify(status)
          
          call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
               rc=status)
          call LIS_verify(status)
          
          !Put in vardata(:) the noah value

          NT=LIS_rc%npatch(n,LIS_rc%lsm_index)
          if(vname.eq."SMCMAX")     then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%smcmax         ;enddo ;endif 
          if(vname.eq."PSISAT")     then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%psisat         ;enddo ;endif
          if(vname.eq."DKSAT")      then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%dksat          ;enddo ;endif
          if(vname.eq."DWSAT")      then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%dwsat          ;enddo ;endif
          if(vname.eq."BEXP")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%bexp           ;enddo ;endif
          if(vname.eq."QUARTZ")     then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%quartz         ;enddo ;endif
          if(vname.eq."RSMIN")      then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%rsmin          ;enddo ;endif
          if(vname.eq."RGL")        then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%rgl            ;enddo ;endif
          if(vname.eq."HS")         then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%hs             ;enddo ;endif
          if(vname.eq."Z0")         then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%z0             ;enddo ;endif
          if(vname.eq."LAI")        then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%lai            ;enddo ;endif
          if(vname.eq."CFACTR")     then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%cfactr         ;enddo ;endif
          if(vname.eq."CMCMAX")     then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%cmcmax         ;enddo ;endif
          if(vname.eq."SBETA")      then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%sbeta          ;enddo ;endif
          if(vname.eq."RSMAX")      then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%rsmax          ;enddo ;endif
          if(vname.eq."TOPT")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%topt           ;enddo ;endif
          if(vname.eq."REFDK")      then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%refdk          ;enddo ;endif
          if(vname.eq."FXEXP")      then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%fxexp          ;enddo ;endif
          if(vname.eq."REFKDT")     then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%refkdt         ;enddo ;endif
          if(vname.eq."CZIL")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%czil           ;enddo ;endif
          if(vname.eq."CSOIL")      then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%csoil          ;enddo ;endif
          if(vname.eq."FRZK")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%frzk           ;enddo ;endif
          if(vname.eq."SNUP")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%snup           ;enddo ;endif
          if(vname.eq."SMCREF")     then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%smcref         ;enddo ;endif
          if(vname.eq."SMCDRY")     then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%smcdry         ;enddo ;endif
          if(vname.eq."SMCWLT")     then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%smcwlt         ;enddo ;endif
          if(vname.eq."F1")         then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%f1             ;enddo ;endif
          if(vname.eq."SLOPE")      then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%slope          ;enddo ;endif
          if(vname.eq."EMISS")      then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%emiss          ;enddo ;endif
          if(vname.eq."SIGMA_FLX")  then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%sigma_flx      ;enddo ;endif 
          if(vname.eq."SIGMA_SM")   then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%sigma_sm       ;enddo ;endif 
          if(vname.eq."SMC1")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%smc(1)         ;enddo ;endif
          if(vname.eq."SMC2")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%smc(2)         ;enddo ;endif
          if(vname.eq."SMC3")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%smc(3)         ;enddo ;endif
          if(vname.eq."SMC4")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%smc(4)         ;enddo ;endif
          if(vname.eq."STC1")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%stc(1)         ;enddo ;endif
          if(vname.eq."STC2")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%stc(2)         ;enddo ;endif
          if(vname.eq."STC3")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%stc(3)         ;enddo ;endif
          if(vname.eq."STC4")       then ; do t=1,NT; vardata(t) = noah36_struc(n)%noah(t)%stc(4)         ;enddo ;endif

             !Test whether any defaults are out of bounds
             count=0
             do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
                do m=1,LIS_rc%nensem(n)
                   gid=(t-1)*LIS_rc%nensem(n)+m
                   if(  (m.eq.1) &
                        .and. &
                        (count.eq.0) &
                        .and. &
                        ((vardata(gid) .lt. noah36_pe_struc(n)%param_min(i)) &
                        .or. &
                        (vardata(gid) .gt. noah36_pe_struc(n)%param_max(i))) ) then
                      count=count+1   
                      print*, '*****************************************************************', '  ', &
                           'WARNING: noah default value is out of LIS-OPT/UE bounds '                , '  ', &
                           'for ', vname                                                             , '  ', &
                           'at '                                                                     , '  ', &
                           'col: ', LIS_surface(n,LIS_rc%lsm_index)%tile(gid)%col                    , '  ', &
                           'row: ', LIS_surface(n,LIS_rc%lsm_index)%tile(gid)%row                    , '  ', &
                           'vegt class: ', noah36_struc(n)%noah(gid)%vegt                            , '  ', &
                           'soiltype: ', noah36_struc(n)%noah(gid)%soiltype                          , '  ', &
                           'default value: ', vardata(gid)                                           , '  ', &
                           'parameter min: ', noah36_pe_struc(n)%param_min(i)                        , '  ', &
                           'parameter max: ', noah36_pe_struc(n)%param_max(i)                        , '  ', &   
                           '*****************************************************************'
                      
                   endif
                enddo
             enddo
          endif
       enddo          
   
       !random initialization
       if(LIS_rc%decSpaceInitMode.eq.1) then  !random initialization 
          seed=seed_base-LIS_localPet !seed must be negative number
          call LIS_rand_func(seed,rand)  !initialize random seed with negative number
          
          do i=1,noah36_pe_struc(n)%nparams
             if(noah36_pe_struc(n)%param_select(i).eq.1) then 
                vname=trim(noah36_pe_struc(n)%param_name(i))
                
                call ESMF_StateGet(DEC_State,vname,varField,rc=status)
                call LIS_verify(status)
                
                call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                     rc=status)
                call LIS_verify(status)
                
                do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
                   do m=1,LIS_rc%nensem(n)
                      if (m.eq.1) then
                         !nothing; leave ensemble 1 with default values
                      else
                         call LIS_rand_func(1,rand)
                         vardata((t-1)*LIS_rc%nensem(n)+m) = &
                              noah36_pe_struc(n)%param_min(i) &
                              + rand * ( noah36_pe_struc(n)%param_max(i) - noah36_pe_struc(n)%param_min(i) )
                      endif
                   enddo
                enddo
             endif
          enddo
       endif

       write(LIS_logunit,*) '[INFO] Finished setting up Noah 3.3 decision space '
     end subroutine noah36_setup_pedecvars

end module noah36_peMod
