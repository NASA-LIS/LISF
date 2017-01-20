!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#define FILENAME "LIS_NUOPC_Cap.F90"
#define MODNAME "LIS_NUOPC_Cap"
#include "LIS_NUOPC_Macros.h"

module LIS_NUOPC
!BOP
!
! !MODULE: LIS_NUOPC
!
! !DESCRIPTION:
!   This modules creates a specialized the NUOPC_Model
!   for LIS.  This is also referred to as the NUOPC Cap.
!
! !REVISION HISTORY:
!  13Oct15    Dan Rosen  Initial Specification
!
  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS        => SetServices, &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_SetClock    => label_SetClock, &
    model_label_CheckImport => label_CheckImport, &
    model_label_Advance     => label_Advance, &
    model_label_Finalize    => label_Finalize
  use LIS_NUOPC_Gluecode
  use beta_NUOPC_Fill
  use beta_NUOPC_Auxiliary
  use beta_NUOPC_Log
  use beta_NUOPC_Base, only: &
   beta_NUOPC_AddNamespace

  implicit none
  
  private
  
  public SetServices

  CHARACTER(LEN=*), PARAMETER :: label_InternalState = 'InternalState'

  type type_InternalStateStruct
    integer               :: verbosity     = VERBOSITY_LV2
    character(len=100)    :: configFile    = 'lis.config'
    logical               :: lwrite_grid   = .FALSE.
    logical               :: lwrite_imp    = .FALSE.
    logical               :: lwrite_exp    = .FALSE.
    logical               :: llog_memory   = .FALSE.
    logical               :: ltestfill_imp = .FALSE.
    logical               :: ltestfill_exp = .FALSE.
    integer               :: nnests        = 0
    integer               :: nfields       = size(LIS_FieldList)
    integer               :: timeSlice     = 1
    type(ESMF_Grid),allocatable         :: grids(:)
    type(ESMF_Clock),allocatable        :: clocks(:)
    type(ESMF_TimeInterval),allocatable :: elapsedtimes(:)
    type(ESMF_State),allocatable        :: NStateImp(:)
    type(ESMF_State),allocatable        :: NStateExp(:)
    integer,allocatable                 :: modes(:)
  end type

  type type_InternalState
    type(type_InternalStateStruct), pointer :: wrap
  end type

!EOP

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "SetServices"
  
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    integer                    :: stat
    type(type_InternalState)   :: is

    rc = ESMF_SUCCESS

    ! allocate memory for this internal state and set it in the component
    allocate(is%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='Allocation of internal state memory failed.', &
      method=METHOD, file=FILENAME, rcToReturn=rc)) return ! bail out
    call ESMF_UserCompSetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeP3, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
       specRoutine=DataInitialize, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call NUOPC_CompSpecialize(gcomp, speclabel=model_label_SetClock, &
      specRoutine=SetClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_MethodRemove(gcomp, label=model_label_CheckImport, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail ou
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
       specRoutine=CheckImport, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail ou
    call NUOPC_CompSpecialize(gcomp, speclabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
      specRoutine=ModelFinalize, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "InitializeP0"

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    ! local variables
    character(ESMF_MAXSTR)     :: cname
    integer                    :: stat
    logical                    :: configIsPresent
    type(ESMF_Config)          :: config
    type(type_InternalState)   :: is
    character(len=10)          :: value

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! Determine Verbosity
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, defaultValue="default", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","debug","default"/), &
      specialValueList=(/VERBOSITY_LV0,VERBOSITY_LV2,VERBOSITY_LV3,VERBOSITY_LV2/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! Write coupled grid files
    call ESMF_AttributeGet(gcomp, name="ConfigFile", value=is%wrap%configFile, &
      defaultValue="lis.config", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! Write coupled grid files
    call ESMF_AttributeGet(gcomp, name="WriteGrid", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%lwrite_grid = (trim(value)=="true")

    ! Write coupled import data files  
    call ESMF_AttributeGet(gcomp, name="WriteImport", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%lwrite_imp = (trim(value)=="true")

    ! Write coupled export data files
    call ESMF_AttributeGet(gcomp, name="WriteExport", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%lwrite_exp = (trim(value)=="true")

    ! Log Memory
    call ESMF_AttributeGet(gcomp, name="LogMemory", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%llog_memory = (trim(value)=="true")

    ! Test fill import fields
    call ESMF_AttributeGet(gcomp, name="TestFillImport", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%ltestfill_imp = (trim(value)=="true")

    ! Test fill export fields
    call ESMF_AttributeGet(gcomp, name="TestFillExport", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%ltestfill_exp = (trim(value)=="true")

    ! Get configuration parameters from attributes or file
    call ESMF_GridCompGet(gcomp, configIsPresent=configIsPresent, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (configIsPresent) then
      call ESMF_GridCompGet(gcomp, config=config, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      call ESMF_ConfigGetAttribute(config, is%wrap%verbosity, &
        label=TRIM(cname)//"_verbosity:", default=is%wrap%verbosity, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      call ESMF_ConfigGetAttribute(config, is%wrap%lwrite_grid, &
        label=TRIM(cname)//"_write_grid:", default=is%wrap%lwrite_grid, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      call ESMF_ConfigGetAttribute(config, is%wrap%lwrite_imp, &
        label=TRIM(cname)//"_write_imp:", default=is%wrap%lwrite_imp, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      call ESMF_ConfigGetAttribute(config, is%wrap%lwrite_exp, &
        label=TRIM(cname)//"_write_exp:", default=is%wrap%lwrite_exp, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      call ESMF_ConfigGetAttribute(config, is%wrap%llog_memory, &
        label=TRIM(cname)//"_log_memory:", default=is%wrap%llog_memory, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      call ESMF_ConfigGetAttribute(config, is%wrap%ltestfill_imp, &
        label=TRIM(cname)//"_testfill_imp:", default=is%wrap%ltestfill_imp, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      call ESMF_ConfigGetAttribute(config, is%wrap%ltestfill_exp, &
        label=TRIM(cname)//"_testfill_exp:", default=is%wrap%ltestfill_exp, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif
    
    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV2) & 
      call InternalConfigLog(trim(cname),gcomp)
    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': leaving '//METHOD, ESMF_LOGMSG_INFO)
  end subroutine
  
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "InitializeP1"

  subroutine InitializeP1(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)     :: gcomp
    type(ESMF_State)        :: importState, exportState
    type(ESMF_Clock)        :: clock
    integer,intent(out)     :: rc
    
    ! LOCAL VARIABLES
    character(ESMF_MAXSTR)     :: cname
    type(type_InternalState)   :: is
    type(ESMF_VM)              :: vm
    integer                    :: localPet, petCount
    integer                    :: stat
    integer                    :: fIndex
    integer                    :: nIndex
    character(len=10)          :: nStr

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': entered '//METHOD, ESMF_LOGMSG_INFO)

    call ESMF_GridCompGet(gcomp, vm=vm, localPet=localPet, petCount=petCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! initialize lis model for this PET
    call LIS_NUOPC_Init(vm, configFile=is%wrap%configFile, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV1) then
      call LIS_Log(trim(cname)//': '//METHOD,rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif

    call LIS_FieldDictionaryAdd(rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    is%wrap%nnests = LIS_NestCntGet(rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    allocate( &
      is%wrap%grids(is%wrap%nnests), &
      is%wrap%clocks(is%wrap%nnests), &
      is%wrap%elapsedtimes(is%wrap%nnests), &
      is%wrap%NStateImp(is%wrap%nnests), &
      is%wrap%NStateExp(is%wrap%nnests), &
      is%wrap%modes(is%wrap%nnests), &
      stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of internal state nest memory failed.", &
      line=__LINE__,file=__FILE__)) &
      return  ! bail out

    is%wrap%modes=LIS_Unknown

    if (is%wrap%nnests.le.1) then
      is%wrap%NStateImp(1) = importState
      is%wrap%NStateExp(1) = exportState
    else
      ! add namespace
      call beta_NUOPC_AddNamespace(importState, &
        domain="1", &
        nestedStateName="NestedStateImp_N1", &
        nestedState=is%wrap%NStateImp(1), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      call beta_NUOPC_AddNamespace(exportState, &
        domain="1", &
        nestedStateName="NestedStateExp_N1", &
        nestedState=is%wrap%NStateExp(1), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif

    do nIndex = 2, is%wrap%nnests
      if (nIndex > 999999999) then
        call ESMF_LogSetError(ESMF_FAILURE, &
          msg="Maximum nest size is 999,999,999.", &
          line=__LINE__,file=__FILE__,rcToReturn=rc)
        return  ! bail out
      endif
      write (nStr,"(I0)") nIndex
      call beta_NUOPC_AddNamespace(importState, &
        domain=trim(nStr), &
        nestedStateName="NestedStateImp_N"//trim(nStr), &
        nestedState=is%wrap%NStateImp(nIndex), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      call beta_NUOPC_AddNamespace(exportState, &
        domain=trim(nStr), &
        nestedStateName="NestedStateExp_N"//trim(nStr), &
        nestedState=is%wrap%NStateExp(nIndex), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    enddo

    !!
    !! advertise import and export fields in each nest
    !!
    do nIndex = 1, is%wrap%nnests
     do fIndex = 1, size(LIS_FieldList)
        if (LIS_FieldList(fIndex)%adImport) then
          call NUOPC_Advertise(is%wrap%NStateImp(nIndex), &
            standardName=trim(LIS_FieldList(fIndex)%stdname), &
            name=trim(LIS_FieldList(fIndex)%stateName), &
            rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
       endif
       if (LIS_FieldList(fIndex)%adExport) then
         call NUOPC_Advertise(is%wrap%NStateExp(nIndex), &
           standardName=trim(LIS_FieldList(fIndex)%stdname), &
           name=trim(LIS_FieldList(fIndex)%stateName), &
           rc=rc)
         if (ESMF_STDERRORCHECK(rc)) return  ! bail out
       endif
      enddo
    enddo

    if (is%wrap%verbosity >= VERBOSITY_LV2) call AdvertiseLog(trim(cname))
    if (is%wrap%verbosity >= VERBOSITY_LV2) &
      call InternalConfigLog(trim(cname),gcomp)
    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': leaving '//METHOD, ESMF_LOGMSG_INFO)
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "InitializeP3"

  subroutine InitializeP3(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    character(ESMF_MAXSTR)     :: cname
    type(type_InternalState)   :: is
    integer                    :: nIndex
    type(ESMF_Field)           :: field
    integer                    :: fIndex
    character(len=9)           :: nStr
    logical                    :: imConn,exConn

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': entered '//METHOD, ESMF_LOGMSG_INFO)

    do nIndex = 1, is%wrap%nnests
      ! Nest integer to string
      if ( nIndex > 999999999) then
        call ESMF_LogSetError(ESMF_FAILURE, &
          msg="Maximum nest size is 999,999,999.", &
          line=__LINE__,file=__FILE__,rcToReturn=rc)
        return  ! bail out
      endif
      write (nStr,"(I0)") nIndex

      ! Call gluecode to create grid.
      is%wrap%grids(nIndex) = LIS_GridCreate(nIndex, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      
      ! Write grid to NetCDF file.
      if (is%wrap%lwrite_grid) then
        call beta_NUOPC_GridWrite(is%wrap%grids(nIndex), &
          trim(cname)//'_grid_nest_'//trim(nStr)//".nc", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      do fIndex = 1, size(LIS_FieldList)
        if (LIS_FieldList(fIndex)%adImport) then
          imConn = NUOPC_IsConnected(is%wrap%NStateImp(nIndex), &
            fieldName=trim(LIS_FieldList(fIndex)%stateName),rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          if (imConn) then
!            if (associated(LIS_FieldList(fIndex)%hookup(nIndex)%importField)) then
!              field = LIS_FieldList(fIndex)%hookup(nIndex)%importField
!            else
              field = ESMF_FieldCreate(name=LIS_FieldList(fIndex)%stateName, &
                grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_R4, rc=rc)
              if (ESMF_STDERRORCHECK(rc)) return  ! bail out
!            endif
            call NUOPC_Realize(is%wrap%NStateImp(nIndex), field=field, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
            LIS_FieldList(fIndex)%realizedImport = .TRUE.
          else
            call ESMF_StateRemove(is%wrap%NStateImp(nIndex), (/trim(LIS_FieldList(fIndex)%stateName)/), &
              relaxedflag=.true.,rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          endif
        endif
      enddo
      do fIndex = 1, size(LIS_FieldList)
        if (LIS_FieldList(fIndex)%adExport) then
          exConn = NUOPC_IsConnected(is%wrap%NStateExp(nIndex), &
            fieldName=trim(LIS_FieldList(fIndex)%stateName),rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          if (exConn) then
!            if (associated(LIS_FieldList(fIndex)%hookup(nIndex)%exportArray)) then
!              field = ESMF_FieldCreate(name=LIS_FieldList(fIndex)%stateName, &
!                grid=is%wrap%grids(nIndex), farray=LIS_FieldList(fIndex)%hookup(nIndex)%exportArray, &
!                indexflag=ESMF_INDEX_DELOCAL, rc=rc)
!              if (ESMF_STDERRORCHECK(rc)) return  ! bail out
!            else
              field = ESMF_FieldCreate(name=LIS_FieldList(fIndex)%stateName, &
                grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_R4, rc=rc)
              if (ESMF_STDERRORCHECK(rc)) return  ! bail out
!            endif
            call NUOPC_Realize(is%wrap%NStateExp(nIndex), field=field,rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
            LIS_FieldList(fIndex)%realizedExport = .TRUE.
          else
            call ESMF_StateRemove(is%wrap%NStateExp(nIndex),(/trim(LIS_FieldList(fIndex)%stateName)/), &
              relaxedflag=.true.,rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          endif
        endif
      enddo
  
      call NUOPC_FillState(is%wrap%NStateImp(nIndex),value=0,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call NUOPC_FillState(is%wrap%NStateExp(nIndex),value=0,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      is%wrap%modes(nIndex) = LIS_RunModeGet(LIS_FieldList,is%wrap%NStateImp(nIndex),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    enddo

    if (is%wrap%verbosity >= VERBOSITY_LV2) call RealizeLog(trim(cname))
    if (is%wrap%verbosity >= VERBOSITY_LV2) &
      call InternalConfigLog(trim(cname),gcomp)
    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': leaving '//METHOD, ESMF_LOGMSG_INFO)
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "DataInitialize"

  subroutine DataInitialize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)                 :: cname
    type(type_InternalState)               :: is
    type(ESMF_Clock)                       :: modelClock
    integer                                :: nIndex
    character(len=10)                      :: nStr
    integer                                :: iIndex
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    integer                                :: stat

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': entered '//METHOD, ESMF_LOGMSG_INFO)

    is%wrap%timeSlice = is%wrap%timeSlice + 1

    ! query the Component for its clock
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    do nIndex=1,is%wrap%nnests
      ! Nest integer to string
      if (nIndex > 999999999) then
        nStr = '999999999+'
      else
        write (nStr,"(I0)") nIndex
      endif

      ! Initialize import and export fields
      call LIS_NUOPC_DataInit(nest=nIndex, &
        importState=is%wrap%NStateImp(nIndex), &
        exportState=is%wrap%NStateExp(nIndex),rc=rc) 
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      ! Fill import fields with test data
      if (is%wrap%ltestfill_imp) then
        call LIS_TestFillImport(nest=nIndex, &
          importState=is%wrap%NStateImp(nIndex), &
          step=is%wrap%timeSlice, &
          label=trim(cname)//": ImportFieldsFill",rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      ! Fill export fields with test data
      if (is%wrap%ltestfill_exp) then
        call LIS_TestFillExport(nest=nIndex, &
          exportState=is%wrap%NStateExp(nIndex), &
          step=is%wrap%timeSlice, &
          label=trim(cname)//": ExportFieldsFill",rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      ! Write initial import fields to file
      if (is%wrap%lwrite_imp) then
        if ( nIndex > 999999999) then
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg="Maximum nest size for writing is 999,999,999.", &
            line=__LINE__,file=__FILE__,rcToReturn=rc)
         return  ! bail out
        endif
        call beta_NUOPC_Write(is%wrap%NStateImp(nIndex), &
          fileNamePrefix="field_"//trim(cname)//"_import_nest_"//trim(nStr)//"_", &
          singleFile=.false., timeslice=is%wrap%timeSlice, &
          relaxedFlag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      ! Write initial export fields to file
      if (is%wrap%lwrite_exp) then
        if ( nIndex > 999999999) then
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg="Maximum nest size for writing is 999,999,999.", &
            line=__LINE__,file=__FILE__,rcToReturn=rc)
         return  ! bail out
        endif
        call beta_NUOPC_Write(is%wrap%NStateExp(nIndex), &
          fileNamePrefix="field_"//trim(cname)//"_export_nest_"//trim(nStr)//"_", &
          singleFile=.false., timeslice=is%wrap%timeSlice, &
          relaxedFlag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      if (is%wrap%verbosity >= VERBOSITY_LV2) then
        call NUOPC_LogState(is%wrap%NStateImp(nIndex), &
          label=trim(cname)//": ImportState Init Nest="//trim(nStr), &
          fvalues=.TRUE.,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        call NUOPC_LogState(is%wrap%NStateExp(nIndex), &
          label=trim(cname)//": ExportState Init Nest="//trim(nStr), &
          fvalues=.TRUE.,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      call ESMF_StateGet(is%wrap%NStateExp(nIndex),itemCount=itemCount, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return ! bail out

      allocate( &
        itemNameList(itemCount), &
        itemTypeList(itemCount), &
        stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Allocation of state item list memory failed.", &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_StateGet(is%wrap%NStateExp(nIndex),itemNameList=itemNameList, &
        itemTypeList=itemTypeList,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      do iIndex=1, itemCount
        if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELD) then
          call ESMF_StateGet(is%wrap%NStateExp(nIndex),field=field, &
            itemName=itemNameList(iIndex),rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        endif
      enddo

      deallocate(itemNameList, itemTypeList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg="Deallocation of state item list memory failed.", &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    enddo

    ! set InitializeDataComplete Attribute to "true", indicating to the
    ! generic code that all inter-model data dependencies are satisfied
    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV3) call LIS_FieldListLog(trim(cname))
    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': leaving '//METHOD, ESMF_LOGMSG_INFO)
  end subroutine
  
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "SetClock"

  subroutine SetClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)     :: cname
    type(type_InternalState)   :: is
    integer                    :: nIndex
    real(ESMF_KIND_R8)         :: mindt
    real(ESMF_KIND_R8)         :: ndt
    type(ESMF_Clock)           :: modelClock
    type(ESMF_TimeInterval)    :: modelTimestep
    type(ESMF_TimeInterval)    :: nestTimeStep

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': entered '//METHOD, ESMF_LOGMSG_INFO)

    ! query the Component for its clock
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! Set minTimestep to the timestep of the first nest
    mindt = LIS_TimestepGet(nest=1,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    do nIndex = 1, is%wrap%nnests
      is%wrap%clocks(nIndex) = ESMF_ClockCreate(modelClock, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      ndt = LIS_TimestepGet(nest=nIndex,rc=rc)
      call ESMF_TimeIntervalSet(nestTimestep, &
        s_r8=ndt, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      call ESMF_ClockSet(is%wrap%clocks(nIndex), &
        timeStep=nestTimestep, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      if (ndt < mindt) mindt = ndt

      call ESMF_TimeIntervalSet(is%wrap%elapsedtimes(nIndex), &
        s_r8=0._ESMF_KIND_R8, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    enddo

    call ESMF_TimeIntervalSet(modelTimestep, &
      s_r8=mindt, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call NUOPC_CompSetClock(gcomp, modelClock, &
      modelTimestep, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV2) &
      call InternalClockLog(trim(cname),gcomp)
    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': leaving '//METHOD, ESMF_LOGMSG_INFO)
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "CheckImport"

subroutine CheckImport(gcomp, rc)
    type(ESMF_GridComp) :: gcomp
    integer,intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)      :: cname
    type(type_InternalState)    :: is
    integer                     :: nIndex
    character(len=10)           :: nStr
    character(len=10)           :: sStr
    type(ESMF_Clock)            :: modelClock
    type(ESMF_Time)             :: modelCurrTime
    type(ESMF_Time)             :: modelStopTime
    logical                     :: allCurrTime
    logical                     :: allStopTime

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': entered '//METHOD, ESMF_LOGMSG_INFO)

    ! query the Component for its clock
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! get the stop time out of the clock
    call ESMF_ClockGet(modelClock, currTime=modelCurrTime, stopTime=modelStopTime, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! check that Fields in the importState show correct timestamp
    do nIndex=1,is%wrap%nnests
      ! Nest integer to string
      if (nIndex > 999999999) then
        nStr = '999999999+'
      else
        write (nStr,"(I0)") nIndex
      endif
      allCurrTime = NUOPC_IsAtTime(is%wrap%NStateImp(nIndex), modelCurrTime, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      allStopTime = NUOPC_IsAtTime(is%wrap%NStateImp(nIndex), modelStopTime, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      if (.not.(allCurrTime.or.allStopTime)) then
        call ESMF_LogSetError(ESMF_FAILURE, &
          msg=METHOD//": NUOPC INCOMPATIBILITY DETECTED: Import Fields "// &
          "Nest="//trim(nStr)//" not at correct time", &
          line=__LINE__,file=__FILE__,rcToReturn=rc)
        return  ! bail out
      endif
    enddo

    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': leaving '//METHOD, ESMF_LOGMSG_INFO)
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "ModelAdvance"

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    character(ESMF_MAXSTR)      :: cname
    type(type_InternalState)    :: is
    integer                     :: nIndex
    character(len=10)           :: nStr
    character(len=10)           :: sStr
    type(ESMF_Clock)            :: modelClock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: modelCurrTime
    type(ESMF_TimeInterval)     :: modelTimeStep
    type(ESMF_TimeInterval)     :: nestTimeStep

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': entered '//METHOD, ESMF_LOGMSG_INFO)

    is%wrap%timeSlice = is%wrap%timeSlice + 1

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call ESMF_ClockGet(modelClock, currTime=modelCurrTime, &
      timeStep=modelTimeStep, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    
    if (is%wrap%timeSlice > 999999999) then
      sStr = '999999999+'
    else
      write (sStr,"(I0)") is%wrap%timeSlice
    endif

    do nIndex=1,is%wrap%nnests
      ! Nest integer to string
      if (nIndex > 999999999) then
        nStr = '999999999+'
      else
        write (nStr,"(I0)") nIndex
      endif

      ! Fill import fields with test data
      if (is%wrap%ltestfill_imp) then
        call LIS_TestFillImport(nest=nIndex, &
          importState=is%wrap%NStateImp(nIndex), &
          step=is%wrap%timeSlice, &
          label=trim(cname)//": ImportFieldsFill",rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      ! Write import fields to file
      if (is%wrap%lwrite_imp) then
        if ( nIndex > 999999999) then
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg="Maximum nest size for writing is 999,999,999.", &
            line=__LINE__,file=__FILE__,rcToReturn=rc)
         return  ! bail out
        endif
        call beta_NUOPC_Write(is%wrap%NStateImp(nIndex), &
          fileNamePrefix="field_"//trim(cname)//"_import_nest_"//trim(nStr)//"_", &
          singleFile=.false., timeslice=is%wrap%timeSlice, &
          relaxedFlag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      if (is%wrap%verbosity >= VERBOSITY_LV2) then
        call NUOPC_LogState(is%wrap%NStateImp(nIndex), &
          label=trim(cname)//": ImportState Slice="//trim(sStr)//" Nest="//trim(nStr), &
          fvalues=.TRUE.,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      is%wrap%elapsedtimes(nIndex) = is%wrap%elapsedtimes(nIndex) + modelTimeStep

      call ESMF_ClockGet(is%wrap%clocks(nIndex),timeStep=nestTimestep,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      
      do while (is%wrap%elapsedtimes(nIndex) >= nestTimestep)
        ! Gluecode NestAdvance
        call ESMF_LogWrite( &
          trim(cname)//': '//METHOD//' Advancing Slice='//trim(sStr)//' Nest='//trim(nStr), &
          ESMF_LOGMSG_INFO)
        call LIS_NUOPC_Run(nIndex,is%wrap%modes(nIndex),is%wrap%timeSlice, &
          is%wrap%NStateImp(nIndex),is%wrap%NStateExp(nIndex), &
          is%wrap%clocks(nIndex), rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        call ESMF_ClockAdvance(is%wrap%clocks(nIndex),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        is%wrap%elapsedtimes(nIndex) = &
          is%wrap%elapsedtimes(nIndex) - nestTimestep
      enddo

      ! Fill export fields with test data
      if (is%wrap%ltestfill_exp) then
        call LIS_TestFillExport(nest=nIndex, &
          exportState=is%wrap%NStateExp(nIndex), &
          step=is%wrap%timeSlice, &
          label=trim(cname)//": ExportFieldsFill",rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      ! Write export fields to file
      if (is%wrap%lwrite_exp) then
        if ( nIndex > 999999999) then
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg="Maximum nest size for writing is 999,999,999.", &
            line=__LINE__,file=__FILE__,rcToReturn=rc)
         return  ! bail out
        endif
        call beta_NUOPC_Write(is%wrap%NStateExp(nIndex), &
          fileNamePrefix="field_"//trim(cname)//"_export_nest_"//trim(nStr)//"_", &
          singleFile=.false., timeslice=is%wrap%timeSlice, &
          relaxedFlag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      if (is%wrap%verbosity >= VERBOSITY_LV2) then
        call NUOPC_LogState(is%wrap%NStateExp(nIndex), &
          label=trim(cname)//": ExportState Slice="//trim(sStr)//" Nest="//trim(nStr), &
          fvalues=.TRUE.,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif
    enddo

    if (is%wrap%verbosity >= VERBOSITY_LV3) call LIS_FieldListLog(trim(cname))
    if (is%wrap%verbosity >= VERBOSITY_LV2) &
      call InternalClockLog(trim(cname),gcomp)
    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': leaving '//METHOD, ESMF_LOGMSG_INFO)
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "ModelFinalize"

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)     :: cname
    integer                    :: stat
    type(type_InternalState)   :: is
    integer                    :: nIndex
    type(ESMF_Clock)           :: clock

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': entered '//METHOD, ESMF_LOGMSG_INFO)

    call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! finalize the LIS model
    do nIndex=1,is%wrap%nnests
      call LIS_NUOPC_Final(nIndex,clock,rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    enddo

    if (is%wrap%verbosity >= VERBOSITY_LV1) &
      call ESMF_LogWrite(trim(cname)//': leaving '//METHOD, ESMF_LOGMSG_INFO)

    deallocate(is%wrap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='Deallocation of internal state memory failed.', &
      method=METHOD,file=FILENAME,rcToReturn=rc)) return ! bail out
  end subroutine

  !-----------------------------------------------------------------------------
  ! Utilities
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "AdvertiseLog"

  subroutine AdvertiseLog(label)
    character(len=*),intent(in) :: label

    ! local variables
    integer                    :: cntImp
    integer                    :: cntExp
    integer                    :: fIndex
    character(ESMF_MAXSTR)     :: logMsg
    integer                    :: rc

    ! Count advertised import and export fields
    cntImp = 0
    cntExp = 0
    do fIndex = 1, size(LIS_FieldList)
      if (LIS_FieldList(fIndex)%adImport) cntImp = cntImp + 1
      if (LIS_FieldList(fIndex)%adExport) cntExp = cntExp + 1
    enddo

    ! Report advertised import fields
    write(logMsg,'(a,a,i0,a)') TRIM(label)//': ', &
      'List of advertised import fields(',cntImp,'):'
    call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
    write(logMsg,'(a,a5,a,a16,a,a)') TRIM(label)//': ', &
      'index',' ','name',' ','standardName'
    call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
    cntImp = 0
    do fIndex=1, size(LIS_FieldList)
      if (.NOT.LIS_FieldList(fIndex)%adImport) cycle
      cntImp = cntImp + 1
      write(logMsg,'(a,i5,a,a16,a,a)') TRIM(label)//': ', &
        cntImp,' ',TRIM(LIS_FieldList(fIndex)%stateName), &
        ' ',TRIM(LIS_FieldList(fIndex)%stdName)
      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    enddo

    ! Report advertised export fields
    write(logMsg,'(a,a,i0,a)') TRIM(label)//': ', &
      'List of advertised export fields(',cntExp,'):'
    call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
    write(logMsg,'(a,a5,a,a16,a,a)') TRIM(label)//': ', &
      'index',' ','name',' ','standardName'
    call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
    cntExp = 0
    do fIndex=1, size(LIS_FieldList)
      if (.NOT.LIS_FieldList(fIndex)%adExport) cycle
      cntExp = cntExp + 1
      write(logMsg,'(a,i5,a,a16,a,a)') TRIM(label)//': ', &
        cntExp,' ',TRIM(LIS_FieldList(fIndex)%stateName), &
        ' ',TRIM(LIS_FieldList(fIndex)%stdName)
      call ESMF_LogWrite(trim(LogMsg), ESMF_LOGMSG_INFO)
    enddo

  end subroutine

#undef METHOD
#define METHOD "RealizeLog"

  subroutine RealizeLog(label)
    character(len=*),intent(in) :: label

    ! local variables
    integer                    :: cntImp
    integer                    :: cntExp
    integer                    :: fIndex
    character(ESMF_MAXSTR)     :: logMsg
    integer                    :: rc

    ! Count advertised import and export fields
    cntImp = 0
    cntExp = 0
    do fIndex = 1, size(LIS_FieldList)
      if (LIS_FieldList(fIndex)%realizedImport) cntImp = cntImp + 1
      if (LIS_FieldList(fIndex)%realizedExport) cntExp = cntExp + 1
    enddo

    ! Report realized import fields
    write(logMsg,'(a,a,i0,a)') TRIM(label)//': ', &
      'List of realized import fields(',cntImp,'):'
    call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
    write(logMsg,'(a,a5,a,a16,a,a)') TRIM(label)//': ', &
      'index',' ','name',' ','standardName'
    call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
    cntImp = 0
    do fIndex=1, size(LIS_FieldList)
      if (.NOT.LIS_FieldList(fIndex)%realizedImport) cycle
      cntImp = cntImp + 1
      write(logMsg,'(a,i5,a,a16,a,a)') TRIM(label)//': ', &
        cntImp,' ',TRIM(LIS_FieldList(fIndex)%stateName), &
        ' ',TRIM(LIS_FieldList(fIndex)%stdName)
      call ESMF_LogWrite(trim(LogMsg), ESMF_LOGMSG_INFO)
    enddo

    ! Report realized export fields
    write(logMsg,'(a,a,i0,a)') TRIM(label)//': ', &
      'List of realized export fields(',cntExp,'):'
    call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
    write(logMsg,'(a,a5,a,a16,a,a)') TRIM(label)//': ', &
      'index',' ','name',' ','standardName'
    call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
    cntExp = 0
    do fIndex=1, size(LIS_FieldList)
      if (.NOT.LIS_FieldList(fIndex)%realizedExport) cycle
      cntExp = cntExp + 1
      write(logMsg,'(a,i5,a,a16,a,a)') TRIM(label)//': ', &
        cntExp,' ',TRIM(LIS_FieldList(fIndex)%stateName), &
        ' ',TRIM(LIS_FieldList(fIndex)%stdName)
      call ESMF_LogWrite(trim(LogMsg), ESMF_LOGMSG_INFO)
    enddo

  end subroutine

#undef METHOD
#define METHOD "InternalConfigLog"

  subroutine InternalConfigLog(label,gcomp)
    character(len=*), intent(in) :: label
    type(ESMF_GridComp)          :: gcomp

    ! local variables
    type(type_InternalState)   :: is
    integer                    :: nIndex
    character(ESMF_MAXSTR)     :: logMsg
    character(len=64)          :: nModeStr
    integer                    :: rc

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (.NOT.(rc.eq.ESMF_SUCCESS)) then
      call ESMF_LogWrite(trim(label)// &
        ' ESMF_UserCompGetInternalState failed.',ESMF_LOGMSG_ERROR)
      return  ! bail out
    endif

    write (logMsg, "(A,(A,I0))") trim(label), &
      ': Verbosity=',is%wrap%verbosity
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,A))") trim(label), &
      ': Config File=',is%wrap%configFile
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(label), &
      ': Write Grid=',is%wrap%lwrite_grid
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(label), &
      ': Write Import=',is%wrap%lwrite_imp
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(label), &
      ': Write Export=',is%wrap%lwrite_exp
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(label), &
      ': Test Fill Import=',is%wrap%ltestfill_imp
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(label), &
      ': Test Fill Export=',is%wrap%ltestfill_exp
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(label), &
      ': Log Memory=',is%wrap%llog_memory
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,I0))") trim(label), &
      ': Nest Count=',is%wrap%nnests
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,I0))") trim(label), &
      ': Time Slice=',is%wrap%timeSlice
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)

    do nIndex=1,is%wrap%nnests
      if (allocated(is%wrap%modes)) then
        select case(is%wrap%modes(nIndex))
        case (LIS_Offline)
          nModeStr ="LIS_Offline"
        case (LIS_Coupled)
          nModeStr = "LIS_Coupled"
        case (LIS_Hybrid)
          nModeStr = "LIS_Hybrid"
        case default
          nModeStr = "LIS_Unknown"
        end select
      else
        nModeStr = "(unallocated)"
      endif
      write (logMsg, "(A,(A,I0),(A,A))") trim(label), &
        ": Nest=",nIndex, &
        " Mode=",trim(nModeStr)
      call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    enddo
  end subroutine

#undef METHOD
#define METHOD "InternalClockLog"

  subroutine InternalClockLog(label,gcomp)
    character(len=*), intent(in) :: label
    type(ESMF_GridComp)          :: gcomp

    ! local variables
    type(type_InternalState)   :: is
    integer                    :: nIndex
    character(ESMF_MAXSTR)     :: logMsg
    type(ESMF_Time)            :: nestCurrTime
    type(ESMF_TimeInterval)    :: nestTimestep
    character(len=64)          :: nCurrTimeStr
    character(len=64)          :: nTimestepStr
    integer                    :: rc

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (.NOT.(rc.eq.ESMF_SUCCESS)) then
      call ESMF_LogWrite(trim(label)// &
        ' ESMF_UserCompGetInternalState failed.',ESMF_LOGMSG_ERROR)
      return  ! bail out
    endif

    write (logMsg, "(A,(A,I0))") trim(label), &
      ': Time Slice=',is%wrap%timeSlice
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)

    do nIndex=1,is%wrap%nnests
      if (allocated(is%wrap%clocks)) then
        if (ESMF_ClockIsCreated(is%wrap%clocks(nIndex))) then
          call ESMF_ClockGet(is%wrap%clocks(nIndex), &
            currTime=nestCurrTime,timeStep=nestTimestep,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          call ESMF_TimeGet(nestCurrTime, &
            timeString=nCurrTimeStr,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          call ESMF_TimeIntervalGet(nestTimestep, &
            timeString=nTimestepStr,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        else
          nCurrTimeStr = "(not_created)"
          nTimestepStr = "(not_created)"
        endif
      else
        nCurrTimeStr = "(unallocated)"
        nTimestepStr = "(unallocated)"
      endif
      write (logMsg, "(A,(A,I0),(A,A))") trim(label), &
        ": Nest=",nIndex, &
        " CurrentTime=",trim(nCurrTimeStr)
      call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
      write (logMsg, "(A,(A,I0),(A,A))") trim(label), &
        ": Nest=",nIndex, &
        " Timestep=",trim(nTimestepStr)
      call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    enddo
  end subroutine

end module

