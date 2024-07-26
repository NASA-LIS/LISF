!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module NoahMP401_peMod
!BOP
!
! !MODULE: NoahMP401_peMod
!
! !DESCRIPTION:
!  This module contains the definitions of the NoahMP.4.0.1 model parameters
!  used in parameter estimation. The data structure is used to expose
!  the LSM parameters to be used in opt/ue. 
!
! !REVISION HISTORY:
!  27 Apr 2020; Sujay Kumar, Initial Code
! !USES:        
  use ESMF
  use LIS_numerRecipesMod, only : LIS_rand_func
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: NoahMP401_setup_pedecvars
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: NoahMP401_pe_struc

!EOP
  type, public ::  NoahMP401_pe_dec 
     integer               :: nparams
     character*40, allocatable :: param_name(:)
     integer     , allocatable :: param_select(:)
     real        , allocatable :: param_min(:)
     real        , allocatable :: param_max(:)
  end type NoahMP401_pe_dec

  type(NoahMP401_pe_dec), allocatable :: NoahMP401_pe_struc(:)

  SAVE
contains

!BOP
! !ROUTINE: NoahMP401_setup_pedecvars
!  \label{NoahMP401_setup_pedecvars}
!
! !REVISION HISTORY:
! 02 Feb 2018: Soni Yatheendradas; Initial Specification
!
! !INTERFACE:
  subroutine NoahMP401_setup_pedecvars(DEC_State, Feas_State)
! !USES:
    use LIS_coreMod
    use LIS_logMod
    use NoahMP401_lsmMod,     only : NoahMP401_struc

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
    
    call ESMF_ConfigGetAttribute(LIS_config,decSpaceAttribsFile,&
         label="LSM Decision space attributes file:",rc=status)
    call LIS_verify(status, "LSM Decision space attributes file: not defined")

    allocate(NoahMP401_pe_struc(LIS_rc%nnest))
    n = 1
    NoahMP401_pe_struc(n)%nparams = 76

    allocate(NoahMP401_pe_struc(n)%param_name(NoahMP401_pe_struc(n)%nparams))
    allocate(NoahMP401_pe_struc(n)%param_select(NoahMP401_pe_struc(n)%nparams))
    allocate(NoahMP401_pe_struc(n)%param_min(NoahMP401_pe_struc(n)%nparams))
    allocate(NoahMP401_pe_struc(n)%param_max(NoahMP401_pe_struc(n)%nparams))

    ! read the attributes file. 
    call LIS_readPEDecSpaceAttributes(decSpaceAttribsFile, &
         NoahMP401_pe_struc(n)%nparams, &
         NoahMP401_pe_struc(n)%param_name, &
         NoahMP401_pe_struc(n)%param_select, &
         NoahMP401_pe_struc(n)%param_min, &
         NoahMP401_pe_struc(n)%param_max)

    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    do i=1,NoahMP401_pe_struc(n)%nparams
       if(NoahMP401_pe_struc(n)%param_select(i).eq.1) then 
          vname=trim(NoahMP401_pe_struc(n)%param_name(i))
          varField = ESMF_FieldCreate(arrayspec=arrspec1, &
               grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
               name=vname,&
               rc=status)
          call LIS_verify(status, &
               'problem with fieldcreate in NoahMP401_setup_pedecvars')
          call ESMF_AttributeSet(varField,'MinRange',&
               NoahMP401_pe_struc(n)%param_min(i),rc=status)
          call LIS_verify(status, &
               'setting minrange to decspace obj in NoahMP401_setup_devars')
          call ESMF_AttributeSet(varField,'MaxRange',&
               NoahMP401_pe_struc(n)%param_max(i),rc=status)
          call LIS_verify(status, &
               'setting maxrange to decspace obj in NoahMP401_setup_devars')

          call ESMF_StateAdd(DEC_State,(/varField/),rc=status)
          call LIS_verify(status,&
               'stateadd in NoahMP401_setup_pedecvars')

          call ESMF_StateGet(DEC_State,vname,varField,rc=status)
          call LIS_verify(status)
          
          call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
               rc=status)
          call LIS_verify(status)
          
          !Put in vardata(:) the noah value

          NT=LIS_rc%npatch(n,LIS_rc%lsm_index)
          if(vname.eq."TOPT")  then  
             do t=1,NT; vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%topt           
             enddo
          endif

          if(vname.eq."RGL")  then 
             do t=1,NT; vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%rgl
             enddo
          endif
          if(vname.eq."RSMAX") then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%rsmax
             enddo
          endif
          if(vname.eq."RSMIN")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%rsmin 
             enddo
          endif
          if(vname.eq."HS") then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%hs     
             enddo
          endif

          if(vname.eq."NROOT") then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%nroot  
             enddo
          endif
          if(vname.eq."CSOIL")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%csoil    
             enddo
          endif

          if(vname.eq."BEXP")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%bexp(1)
             enddo
          endif

          if(vname.eq."DKSAT")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%dksat(1)
             enddo
          endif

          if(vname.eq."DWSAT")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%dwsat(1)
             enddo
          endif

          if(vname.eq."PSISAT")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%psisat(1)
             enddo
          endif

          if(vname.eq."QUARTZ")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%quartz(1)   
             enddo
          endif

          if(vname.eq."SMCMAX")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%smcmax(1)   
             enddo
          endif

          if(vname.eq."SMCREF")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%smcref(1)
             enddo
          endif

          if(vname.eq."SMCWLT")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%smcwlt(1)
             enddo
          endif

          if(vname.eq."CZIL") then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%czil   
             enddo
          endif
                

          if(vname.eq."SLOPE") then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%slope
             enddo
          endif

          if(vname.eq."CH2OP")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%CH2OP
             enddo
          endif

          if(vname.eq."DLEAF")      then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%DLEAF
             enddo
          endif

          if(vname.eq."Z0MVT")      then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%Z0MVT
             enddo
          endif
          
          if(vname.eq."HVT")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%HVT
             enddo
          endif

          if(vname.eq."HVB")      then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%HVB
             enddo
          endif
          if(vname.eq."RC")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%RC
             enddo
          endif
          if(vname.eq."MFSNO")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%mfsno
             enddo
          endif
          if(vname.eq."ALBSAT1")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%albsat(1)
             enddo
          endif
          if(vname.eq."ALBSAT2")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%albsat(2)
             enddo
          endif
          if(vname.eq."ALBDRY1")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%albdry(1)
             enddo
          endif
          if(vname.eq."ALBDRY2")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%albdry(2)
             enddo
          endif
          if(vname.eq."ALBICE1")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%albice(1)
             enddo
          endif
          if(vname.eq."ALBICE2")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%albice(2)
             enddo
          endif
          if(vname.eq."OMEGAS1")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%omegas(1)
             enddo
          endif
          if(vname.eq."OMEGAS2")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%omegas(2)
             enddo
          endif
          if(vname.eq."BETADS")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%betads
             enddo
          endif
          if(vname.eq."BETAIS")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%betais
             enddo
          endif
          if(vname.eq."EG1")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%eg(1)
             enddo
          endif
          if(vname.eq."EG2")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%eg(2)
             enddo
          endif
          if(vname.eq."EG2")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%eg(2)
             enddo
          endif
          if(vname.eq."Z0SNO")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%z0sno
             enddo
          endif
          if(vname.eq."SSI")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%ssi
             enddo
          endif
          if(vname.eq."SWEMX")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%swemx
             enddo
          endif
          if(vname.eq."RSURF_SNOW")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%rsurf_snow
             enddo
          endif
          if(vname.eq."MNSNALB")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%mnsnalb
             enddo
          endif
          if(vname.eq."MXSNALB")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%mxsnalb
             enddo
          endif
          if(vname.eq."SNDECAYEXP")      then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%sndecayexp
             enddo
          endif

          if(vname.eq."T_ULIMIT") then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%t_ulimit
             enddo
          endif

          if(vname.eq."T_LLIMIT") then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%t_llimit
             enddo
          endif

          if(vname.eq."T_MLIMIT") then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%t_mlimit
             enddo
          endif
          
          if(vname.eq."SNOWF_SCALEF") then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%snowf_scalef
             enddo
          endif

          if(vname.eq."RHOL1")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%RHOL(1)
             enddo
          endif

          if(vname.eq."RHOL2")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%RHOL(2)
             enddo
          endif

          if(vname.eq."RHOS1")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%RHOS(1)
             enddo
          endif

          if(vname.eq."RHOS2")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%RHOS(2)
             enddo
          endif

          if(vname.eq."TAUL1")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%TAUL(1)
             enddo
          endif

          if(vname.eq."TAUL2")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%TAUL(2)
             enddo
          endif

          if(vname.eq."TAUS1")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%TAUS(1)
             enddo
          endif

          if(vname.eq."TAUS2")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%TAUS(2)
             enddo
          endif

          if(vname.eq."XL")  then
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%XL
             enddo
          endif

          if(vname.eq."CWPVT")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%CWPVT
             enddo
          endif

          if(vname.eq."C3PSN")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%C3PSN
             enddo
          endif

          if(vname.eq."KC25")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%KC25
             enddo
          endif

          if(vname.eq."AKC")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%AKC
             enddo
          endif

          if(vname.eq."KO25")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%KO25
             enddo
          endif

          if(vname.eq."AKO")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%AKO
             enddo
          endif

          if(vname.eq."AVCMX")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%AVCMX
             enddo
          endif

          if(vname.eq."AQE")  then
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%AQE
             enddo
          endif

          if(vname.eq."LTOVRC")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%LTOVRC
             enddo
          endif
                
          if(vname.eq."DILEFC")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%DILEFC
             enddo
          endif

          if(vname.eq."DILEFW")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%DILEFW
             enddo
          endif

          if(vname.eq."RMF25")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%RMF25
             enddo
          endif

          if(vname.eq."SLA")  then
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%SLA
             enddo
          endif
                
          if(vname.eq."FRAGR")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%FRAGR
             enddo
          endif
                
          if(vname.eq."TMIN")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%TMIN
             enddo
          endif

          if(vname.eq."VCMX25")  then 
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%VCMX25
             enddo
          endif

          if(vname.eq."TDLEF")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%TDLEF
             enddo
          endif

          if(vname.eq."BP")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%BP
             enddo
          endif

          if(vname.eq."MP")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%MP
             enddo
          endif

          if(vname.eq."QE25")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%QE25
             enddo
          endif

          if(vname.eq."RMS25")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%RMS25
             enddo
          endif

          if(vname.eq."RMR25")  then
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%RMR25
             enddo
          endif

          if(vname.eq."ARM")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%ARM
             enddo
          endif

          if(vname.eq."FOLNMX")  then
             do t=1,NT
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%FOLNMX
             enddo
          endif

          if(vname.eq."WDPOOL")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%WDPOOL
             enddo
          endif

          if(vname.eq."WRRAT")  then
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%WRRAT
             enddo
          endif

          if(vname.eq."MRP")  then 
             do t=1,NT 
                vardata(t) = NoahMP401_struc(n)%noahmp401(t)%param%MRP
             enddo
          endif

          !Test whether any defaults are out of bounds
          count=0
          do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
             do m=1,LIS_rc%nensem(n)
                gid=(t-1)*LIS_rc%nensem(n)+m
                if(  (m.eq.1) &
                     .and. &
                     (count.eq.0) &
                     .and. &
                     ((vardata(gid) .lt. NoahMP401_pe_struc(n)%param_min(i)) &
                     .or. &
                     (vardata(gid) .gt. NoahMP401_pe_struc(n)%param_max(i))) ) then
                   count=count+1   
                   write(LIS_logunit,*) '*****************************************************************', '  ', &
                        'WARNING: noah default value is out of LIS-OPT/UE bounds '                , '  ', &
                        'for ', vname                                                             , '  ', &
                        'at '                                                                     , '  ', &
                        'col: ', LIS_surface(n,LIS_rc%lsm_index)%tile(gid)%col                    , '  ', &
                        'row: ', LIS_surface(n,LIS_rc%lsm_index)%tile(gid)%row                    , '  ', &
                        'vegt class: ', NoahMP401_struc(n)%noahmp401(gid)%vegetype                            , '  ', &
                        'soiltype: ', NoahMP401_struc(n)%noahmp401(gid)%soiltype                          , '  ', &
                        'default value: ', vardata(gid)                                           , '  ', &
                        'parameter min: ', NoahMP401_pe_struc(n)%param_min(i)                        , '  ', &
                        'parameter max: ', NoahMP401_pe_struc(n)%param_max(i)                        , '  ', &   
                        '*****************************************************************'
                   
                endif
             enddo
          enddo
       endif ! if(NoahMP401_pe_struc(n)%param_select(i).eq.1) then 
    enddo  ! do i=1,NoahMP401_pe_struc(n)%nparams
   
    !random initialization
    if(LIS_rc%decSpaceInitMode.eq.1) then  !random initialization 
       seed=seed_base-LIS_localPet !seed must be negative number
       call LIS_rand_func(seed,rand)  !initialize random seed with negative number
       
       do i=1,NoahMP401_pe_struc(n)%nparams
          if(NoahMP401_pe_struc(n)%param_select(i).eq.1) then 
             vname=trim(NoahMP401_pe_struc(n)%param_name(i))
             
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
                           NoahMP401_pe_struc(n)%param_min(i) &
                           + rand * ( NoahMP401_pe_struc(n)%param_max(i) - NoahMP401_pe_struc(n)%param_min(i) )
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
  end subroutine NoahMP401_setup_pedecvars

end module NoahMP401_peMod
