!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LSobjFunc_Mod
!BOP
! !MODULE: LSobjFunc_Mod
! 
! !DESCRIPTION: 
!  This module provides the objects and methods to compute a Least Square (LS) metric
!  for use in parameter estimation. 
! 
! !REVISION HISTORY: 
!  
! 15 Jul 2009: Sujay Kumar; Initial implementation
! 
! !USES:
  use ESMF

  implicit none

  PRIVATE
 
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public ::  initializeLSObjFunc

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: ls_ctl

!EOP
  type lsctl
     integer :: LSobjfunc_mode
     integer :: LSobjfunc_minobs
     character*100, allocatable :: obj_name(:)
     integer :: nobjs
     real, allocatable :: w(:)  !weights
  end type lsctl
  type(lsctl) :: ls_ctl

contains

!BOP
! !ROUTINE: intializeLSObjFunc
! \label{initializeLSObjFunc}
! 
! !INTERFACE: 
  subroutine initializeLSObjFunc()
! !USES: 
    use LIS_coreMod,         only : LIS_vecTile, LIS_config
    use LIS_optUEMod,        only : LIS_ObjectiveFunc
    use LIS_logMod,          only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify
    use LIS_PE_HandlerMod,   only : LIS_PEOBS_State
    use LIS_constantsMod,    only : LIS_CONST_PATH_LEN
! 
! !DESCRIPTION:
!  This method initializes the objects to be used in the least square 
!  computations
!EOP    
    implicit none

    type(ESMF_ArraySpec) :: arrspec1
    type(ESMF_Field)     :: objfuncField
    type(ESMF_Field)     :: minField
    type(ESMF_Field)     :: maxField
    type(ESMF_Field)     :: numobsField
    type(ESMF_Field)     :: modelvField
    type(ESMF_Field)     :: nummodelvField
    real,    pointer     :: numobs(:), nummodelv(:)
    real,    pointer     :: objfunc(:),modelv(:)
    character(len=LIS_CONST_PATH_LEN) :: LSweightsFile  !link to lis.config entry for optUE
    character*100        :: obj_name_temp
    character*200          :: line
    integer              :: n, nobjs, k, i, j 
    integer              :: status
    real                 :: weight
    logical              :: in_namelist
    integer                :: ftn

    n = 1

    call ESMF_StateGet(LIS_PEOBS_State, itemCount =ls_ctl%nobjs, rc=status)
    call LIS_verify(status)

    allocate(ls_ctl%obj_name(ls_ctl%nobjs))
    allocate(ls_ctl%w(ls_ctl%nobjs))

    ls_ctl%w=1.0  !initialize weights to 1 for single objective case

    call ESMF_StateGet(LIS_PEOBS_State, itemNameList=ls_ctl%obj_name,rc=status)
    call LIS_verify(status)

    if(ls_ctl%nobjs.gt.1) then
       call ESMF_ConfigGetAttribute(LIS_config,LSweightsFile,&
            label="Least Squares objective function weights file:",rc=status)
       call LIS_verify(status, 'Least Squares objective function weights file: not defined')

       ftn = LIS_getNextUnitNumber()
       write(LIS_logunit,*) '[INFO] Reading Least Squares objective function weights ...', &
            trim(LSweightsFile)
       open(ftn,file=trim(LSweightsFile),status='old')
       read(ftn,*)  !demarcate section
       do k=1, ls_ctl%nobjs
          read(ftn,*) obj_name_temp, weight  
          in_namelist=.false.  !initialize
          do i=1, ls_ctl%nobjs
             if(obj_name_temp.eq.ls_ctl%obj_name(i)) then
                ls_ctl%w(i)=weight
                in_namelist=.true.
             endif
          enddo
          if(.not.in_namelist) then
             write(LIS_logunit,*) '[ERR] The variable ',obj_name_temp, ' is not in list ', &
                  'of objective function names that includes'
             do j=1,ls_ctl%nobjs
                write(LIS_logunit,*) ls_ctl%obj_name(j)
             enddo
             write(LIS_logunit,*) '[ERR] Program stopping....'
             call LIS_endrun
          endif
       enddo
       
    endif

    call ESMF_ConfigGetAttribute(LIS_config,ls_ctl%LSobjfunc_mode, &
         label="Least Squares objective function mode:",rc=status)
    call LIS_verify(status, &
         'Least Squares objective function mode: not defined')

    call ESMF_ConfigGetAttribute(LIS_config,ls_ctl%LSobjfunc_minobs, &
         label="Least Squares objective function minimum number of obs:",rc=status)
    call LIS_verify(status, &
         'Least Squares objective function minimum number of obs:')

    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    objfuncField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Objective Function Value",rc=status)
    call LIS_verify(status)
    
    minField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Min Criteria Value",rc=status)
    call LIS_verify(status)
    
    maxField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Max Criteria Value",rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(objfuncField, localDE=0, farrayPtr=objfunc, rc=status)
    call LIS_verify(status)
    objfunc = 0
    
    numobsField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
         name="Number of observations",rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(numobsField, localDE=0, farrayPtr=numobs, rc=status)
    call LIS_verify(status)
    numobs = 0
    
    if(ls_ctl%LSobjfunc_mode.eq.3) then
       modelvField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
            name="Model obspred",rc=status)
       call LIS_verify(status)
       
       call ESMF_FieldGet(modelvField, localDE=0, farrayPtr=modelv, rc=status)
       call LIS_verify(status)
       modelv = 0
       
       nummodelvField = ESMF_FieldCreate(arrayspec=arrspec1,grid=LIS_vecTile(n),&
            name="Count Model obspred",rc=status)
       call LIS_verify(status)
       
       call ESMF_FieldGet(nummodelvField, localDE=0, farrayPtr=nummodelv, rc=status)
       call LIS_verify(status)
       nummodelv = 0 

       call ESMF_StateAdd(LIS_ObjectiveFunc, (/modelvField/), rc=status)
       call LIS_verify(status)  

       call ESMF_StateAdd(LIS_ObjectiveFunc, (/nummodelvField/), rc=status)
       call LIS_verify(status)  
    endif

    call ESMF_StateAdd(LIS_ObjectiveFunc, (/objfuncField/), rc=status)
    call LIS_verify(status)  
    
    call ESMF_StateAdd(LIS_ObjectiveFunc, (/minField/), rc=status)
    call LIS_verify(status)  
    
    call ESMF_StateAdd(LIS_ObjectiveFunc, (/maxField/), rc=status)
    call LIS_verify(status)  
    
    call ESMF_StateAdd(LIS_ObjectiveFunc, (/numobsField/), rc=status)
    call LIS_verify(status)  
    
  end subroutine initializeLSObjFunc
  

end module LSobjFunc_Mod
