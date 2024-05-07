!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module NoahMP36_peMod
!BOP
!
! !MODULE: NoahMP36_peMod
!
! !DESCRIPTION:
!  This module contains the definitions of the NoahMP.3.6 model parameters
!  used in parameter estimation. The data structure is used to expose
!  the LSM parameters to be used in opt/ue. 
!
! !REVISION HISTORY:
!  2 Feb 2018; Soni Yatheendradas, Initial Code
! !USES:        
  use ESMF
  use LIS_numerRecipesMod, only : LIS_rand_func
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: NoahMP36_setup_pedecvars
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: NoahMP36_pe_struc

!EOP
  type, public ::  NoahMP36_pe_dec 
     integer               :: nparams
     character*40, allocatable :: param_name(:)
     integer     , allocatable :: param_select(:)
     real        , allocatable :: param_min(:)
     real        , allocatable :: param_max(:)
  end type NoahMP36_pe_dec

  type(NoahMP36_pe_dec), allocatable :: NoahMP36_pe_struc(:)

  SAVE
contains

!BOP
! !ROUTINE: NoahMP36_setup_pedecvars
!  \label{NoahMP36_setup_pedecvars}
!
! !REVISION HISTORY:
! 02 Feb 2018: Soni Yatheendradas; Initial Specification
!
! !INTERFACE:
  subroutine NoahMP36_setup_pedecvars(DEC_State, Feas_State)
! !USES:
    use LIS_coreMod,       only : LIS_rc, LIS_config,LIS_vecPatch, LIS_surface, LIS_localPET
    use LIS_logMod,        only : LIS_logunit, LIS_verify
    use NoahMP36_lsmMod,     only : NoahMP36_struc

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
    
!    write(LIS_logunit,*) '[INFO] Here 1 in NoahMP36_setup_pedecvars, mod_flag: ', mod_flag

!    mod_flag = 0  !now initialized centrally in LIS_optUEMod.F90 as lsms,rtms,etc. can raise flag

    call ESMF_ConfigGetAttribute(LIS_config,decSpaceAttribsFile,&
         label="LSM Decision space attributes file:",rc=status)
    call LIS_verify(status, "LSM Decision space attributes file: not defined")

    allocate(NoahMP36_pe_struc(LIS_rc%nnest))
    n = 1
    NoahMP36_pe_struc(n)%nparams = 62

    allocate(NoahMP36_pe_struc(n)%param_name(NoahMP36_pe_struc(n)%nparams))
    allocate(NoahMP36_pe_struc(n)%param_select(NoahMP36_pe_struc(n)%nparams))
    allocate(NoahMP36_pe_struc(n)%param_min(NoahMP36_pe_struc(n)%nparams))
    allocate(NoahMP36_pe_struc(n)%param_max(NoahMP36_pe_struc(n)%nparams))

    ! read the attributes file. 
    call LIS_readPEDecSpaceAttributes(decSpaceAttribsFile, &
         NoahMP36_pe_struc(n)%nparams, &
         NoahMP36_pe_struc(n)%param_name, &
         NoahMP36_pe_struc(n)%param_select, &
         NoahMP36_pe_struc(n)%param_min, &
         NoahMP36_pe_struc(n)%param_max)
!    write(LIS_logunit,*) 'NoahMP36_pe_struc(n)%param_name: ', &
!            NoahMP36_pe_struc(n)%param_name

    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    do i=1,NoahMP36_pe_struc(n)%nparams
       if(NoahMP36_pe_struc(n)%param_select(i).eq.1) then 
          vname=trim(NoahMP36_pe_struc(n)%param_name(i))
          varField = ESMF_FieldCreate(arrayspec=arrspec1, &
               grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
               name=vname,&
               rc=status)
          call LIS_verify(status, &
               'problem with fieldcreate in NoahMP36_setup_pedecvars')
          call ESMF_AttributeSet(varField,'MinRange',&
               NoahMP36_pe_struc(n)%param_min(i),rc=status)
          call LIS_verify(status, &
               'setting minrange to decspace obj in NoahMP36_setup_devars')
          call ESMF_AttributeSet(varField,'MaxRange',&
               NoahMP36_pe_struc(n)%param_max(i),rc=status)
          call LIS_verify(status, &
               'setting maxrange to decspace obj in NoahMP36_setup_devars')

          call ESMF_StateAdd(DEC_State,(/varField/),rc=status)
          call LIS_verify(status,&
               'stateadd in NoahMP36_setup_pedecvars')

          call ESMF_StateGet(DEC_State,vname,varField,rc=status)
          call LIS_verify(status)
          
          call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
               rc=status)
          call LIS_verify(status)
          
          !Put in vardata(:) the noah value

          NT=LIS_rc%npatch(n,LIS_rc%lsm_index)
          if(vname.eq."TOPT")       then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%topt           ;enddo ;endif
          if(vname.eq."RGL")        then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%rgl            ;enddo ;endif
          if(vname.eq."RSMAX")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%rsmax          ;enddo ;endif
          if(vname.eq."RSMIN")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%rsmin          ;enddo ;endif
          if(vname.eq."HS")         then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%hs             ;enddo ;endif
          if(vname.eq."NROOT")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%nroot          ;enddo ;endif
          if(vname.eq."CSOIL")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%csoil          ;enddo ;endif
          if(vname.eq."BEXP")       then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%bexp           ;enddo ;endif
          if(vname.eq."DKSAT")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%dksat          ;enddo ;endif
          if(vname.eq."DWSAT")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%dwsat          ;enddo ;endif
          if(vname.eq."PSISAT")     then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%psisat         ;enddo ;endif
          if(vname.eq."QUARTZ")     then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%quartz         ;enddo ;endif
          if(vname.eq."SMCMAX")     then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%smcmax         ;enddo ;endif 
          if(vname.eq."SMCREF")     then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%smcref         ;enddo ;endif
          if(vname.eq."SMCWLT")     then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%smcwlt         ;enddo ;endif
          if(vname.eq."CZIL")       then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%czil           ;enddo ;endif
          if(vname.eq."FRZK")       then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%frzk           ;enddo ;endif
          if(vname.eq."REFDK")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%refdk          ;enddo ;endif
          if(vname.eq."REFKDT")     then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%refkdt         ;enddo ;endif
          if(vname.eq."SLOPE")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%slope          ;enddo ;endif
          if(vname.eq."CH2OP")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%CH2OP;enddo ;endif
          if(vname.eq."DLEAF")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%DLEAF;enddo ;endif
          if(vname.eq."Z0MVT")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%Z0MVT;enddo ;endif
          if(vname.eq."HVT")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%HVT;enddo ;endif
          if(vname.eq."HVB")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%HVB;enddo ;endif
          if(vname.eq."RC")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%RC;enddo ;endif
          if(vname.eq."RHOL1")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%RHOL1;enddo ;endif
          if(vname.eq."RHOL2")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%RHOL2;enddo ;endif
          if(vname.eq."RHOS1")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%RHOS1;enddo ;endif
          if(vname.eq."RHOS2")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%RHOS2;enddo ;endif
          if(vname.eq."TAUL1")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%TAUL1;enddo ;endif
          if(vname.eq."TAUL2")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%TAUL2;enddo ;endif
          if(vname.eq."TAUS1")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%TAUS1;enddo ;endif
          if(vname.eq."TAUS2")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%TAUS2;enddo ;endif
          if(vname.eq."XL")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%XL;enddo ;endif
          if(vname.eq."CWPVT")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%CWPVT;enddo ;endif
          if(vname.eq."C3PSN")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%C3PSN;enddo ;endif
          if(vname.eq."KC25")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%KC25;enddo ;endif
          if(vname.eq."AKC")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%AKC;enddo ;endif
          if(vname.eq."KO25")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%KO25;enddo ;endif
          if(vname.eq."AKO")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%AKO;enddo ;endif
          if(vname.eq."AVCMX")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%AVCMX;enddo ;endif
          if(vname.eq."AQE")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%AQE;enddo ;endif
          if(vname.eq."LTOVRC")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%LTOVRC;enddo ;endif
          if(vname.eq."DILEFC")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%DILEFC;enddo ;endif
          if(vname.eq."DILEFW")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%DILEFW;enddo ;endif
          if(vname.eq."RMF25")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%RMF25;enddo ;endif
          if(vname.eq."SLA")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%SLA;enddo ;endif
          if(vname.eq."FRAGR")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%FRAGR;enddo ;endif
          if(vname.eq."TMIN")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%TMIN;enddo ;endif
          if(vname.eq."VCMX25")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%VCMX25;enddo ;endif
          if(vname.eq."TDLEF")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%TDLEF;enddo ;endif
          if(vname.eq."BP")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%BP;enddo ;endif
          if(vname.eq."MP")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%MP;enddo ;endif
          if(vname.eq."QE25")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%QE25;enddo ;endif
          if(vname.eq."RMS25")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%RMS25;enddo ;endif
          if(vname.eq."RMR25")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%RMR25;enddo ;endif
          if(vname.eq."ARM")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%ARM;enddo ;endif
          if(vname.eq."FOLNMX")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%FOLNMX;enddo ;endif
          if(vname.eq."WDPOOL")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%WDPOOL;enddo ;endif
          if(vname.eq."WRRAT")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%WRRAT;enddo ;endif
          if(vname.eq."MRP")      then ; do t=1,NT; vardata(t) = NoahMP36_struc(n)%noahmp36(t)%MRP;enddo ;endif

             !Test whether any defaults are out of bounds
             count=0
             do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
                do m=1,LIS_rc%nensem(n)
                   gid=(t-1)*LIS_rc%nensem(n)+m
                   if(  (m.eq.1) &
                        .and. &
                        (count.eq.0) &
                        .and. &
                        ((vardata(gid) .lt. NoahMP36_pe_struc(n)%param_min(i)) &
                        .or. &
                        (vardata(gid) .gt. NoahMP36_pe_struc(n)%param_max(i))) ) then
                      count=count+1   
                      print*, '*****************************************************************', '  ', &
                           'WARNING: noah default value is out of LIS-OPT/UE bounds '                , '  ', &
                           'for ', vname                                                             , '  ', &
                           'at '                                                                     , '  ', &
                           'col: ', LIS_surface(n,LIS_rc%lsm_index)%tile(gid)%col                    , '  ', &
                           'row: ', LIS_surface(n,LIS_rc%lsm_index)%tile(gid)%row                    , '  ', &
                           'vegt class: ', NoahMP36_struc(n)%noahmp36(gid)%vegetype                            , '  ', &
                           'soiltype: ', NoahMP36_struc(n)%noahmp36(gid)%soiltype                          , '  ', &
                           'default value: ', vardata(gid)                                           , '  ', &
                           'parameter min: ', NoahMP36_pe_struc(n)%param_min(i)                        , '  ', &
                           'parameter max: ', NoahMP36_pe_struc(n)%param_max(i)                        , '  ', &   
                           '*****************************************************************'
                      
                   endif
                enddo
             enddo
       endif ! if(NoahMP36_pe_struc(n)%param_select(i).eq.1) then 
    enddo  ! do i=1,NoahMP36_pe_struc(n)%nparams
   
    !random initialization
    if(LIS_rc%decSpaceInitMode.eq.1) then  !random initialization 
       seed=seed_base-LIS_localPet !seed must be negative number
       call LIS_rand_func(seed,rand)  !initialize random seed with negative number
       
       do i=1,NoahMP36_pe_struc(n)%nparams
          if(NoahMP36_pe_struc(n)%param_select(i).eq.1) then 
             vname=trim(NoahMP36_pe_struc(n)%param_name(i))
!             write(LIS_logunit,*) 'Param: ', vname
             
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
                           NoahMP36_pe_struc(n)%param_min(i) &
                           + rand * ( NoahMP36_pe_struc(n)%param_max(i) - NoahMP36_pe_struc(n)%param_min(i) )
                   endif
                enddo
             enddo
             if(vname.eq."NROOT") then !integer value

                do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
                   do m=1,LIS_rc%nensem(n)
                      vardata((t-1)*LIS_rc%nensem(n)+m) = &
                           nint(vardata((t-1)*LIS_rc%nensem(n)+m))
                   enddo
                enddo
                
             endif
          endif
       enddo
    endif

    write(LIS_logunit,*) '[INFO] Finished setting up NoahMP 3.6 decision space '
  end subroutine NoahMP36_setup_pedecvars

end module NoahMP36_peMod
