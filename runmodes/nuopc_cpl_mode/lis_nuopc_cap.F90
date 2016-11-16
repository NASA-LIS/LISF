!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#define ESMF_STDERRORCHECK(rc) ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)
#define FILENAME "lis_nuopc_cap"
#define MODNAME "lis_nuopc"

#define VERBOSITY_MIN 0
#define VERBOSITY_MAX 1
#define VERBOSITY_DBG 1023

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
    model_label_Advance     => label_Advance, &
    model_label_Finalize    => label_Finalize
  use LIS_NUOPC_Gluecode
  use NUOPC_FillUtility
  use NUOPC_FileWriteUtility
  use NUOPC_LogUtility
  use NUOPC_MultiNestConnector, only: &
   NUOPC_AddNamespaceWithNest 

  implicit none
  
  private
  
  public SetServices

  CHARACTER(LEN=*), PARAMETER :: label_InternalState = 'InternalState'

  type type_InternalStateStruct
    integer               :: verbosity = VERBOSITY_MAX
    logical               :: gridwrite_flag = .FALSE.
    logical               :: statewrite_flag = .FALSE.
    logical               :: profile_memory = .FALSE.
    integer               :: nnests = -1
    integer               :: nfields = size(LIS_FieldList)
    integer               :: slice = -1
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
  
  subroutine SetServices(lisGridComp, rc)
    type(ESMF_GridComp)  :: lisGridComp
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
    call ESMF_UserCompSetInternalState(lisGridComp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(lisGridComp, model_routine_SS, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(lisGridComp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(lisGridComp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call NUOPC_CompSetEntryPoint(lisGridComp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(lisGridComp, specLabel=model_label_DataInitialize, &
       specRoutine=DataInitialize, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call NUOPC_CompSpecialize(lisGridComp, speclabel=model_label_SetClock, &
      specRoutine=SetClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call NUOPC_CompSpecialize(lisGridComp, speclabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call NUOPC_CompSpecialize(lisGridComp, specLabel=model_label_Finalize, &
      specRoutine=ModelFinalize, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "InitializeP0"

  subroutine InitializeP0(lisGridComp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: lisGridComp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    ! local variables
    character(ESMF_MAXSTR)     :: cname
    integer                    :: stat
    type(ESMF_Config)          :: config
    type(type_InternalState)   :: is
    character(len=10)          :: value

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(lisGridComp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(lisGridComp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! Get configuration parameters from attributes or file
    call ESMF_GridCompGet(lisGridComp, config=config, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call ESMF_AttributeGet(lisGridComp, name="Verbosity", value=value, defaultValue="max", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","debug"/), &
      specialValueList=(/VERBOSITY_MIN,VERBOSITY_MAX,VERBOSITY_DBG/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call ESMF_ConfigGetAttribute(config, is%wrap%verbosity, &
       label=TRIM(cname)//"_verbosity:", default=is%wrap%verbosity, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call ESMF_AttributeGet(lisGridComp, name="WriteGrids", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%gridwrite_flag = (trim(value)=="true")

    call ESMF_ConfigGetAttribute(config, is%wrap%gridwrite_flag, &
       label=TRIM(cname)//"_grid_write:", default=is%wrap%gridwrite_flag, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call ESMF_AttributeGet(lisGridComp, name="WriteData", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%statewrite_flag = (trim(value)=="true")

    call ESMF_ConfigGetAttribute(config, is%wrap%statewrite_flag, &
       label=TRIM(cname)//"_write_data:", default=is%wrap%statewrite_flag, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call ESMF_AttributeGet(lisGridComp, name="ProfileMemory", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    is%wrap%profile_memory = (trim(value)=="true")

    call ESMF_ConfigGetAttribute(config, is%wrap%profile_memory, &
       label=TRIM(cname)//"_profile_memory:", default=is%wrap%profile_memory, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(lisGridComp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if ( is%wrap%verbosity >= VERBOSITY_MAX ) then
     call InternalStateLog(lisGridComp,trim(cname)//': '//METHOD//' Complete',rc)
     if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif
  end subroutine
  
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "InitializeAdvertise"

  subroutine InitializeAdvertise(lisGridComp, importState, exportState, clock, rc)
    type(ESMF_GridComp)     :: lisGridComp
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
    call ESMF_GridCompGet(lisGridComp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(lisGridComp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call ESMF_GridCompGet(lisGridComp, vm=vm, localPet=localPet, petCount=petCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! initialize lis model for this PET
    call LIS_NUOPC_Init(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_MAX) then
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

    if (is%wrap%nnests == 1) then
      is%wrap%NStateImp(1) = importState
      is%wrap%NStateExp(1) = exportState
    else
      ! add namespace
      call NUOPC_AddNamespaceWithNest(importState, &
        nest="1", &
        nestedStateName="NestedStateImp_N1", &
        nestedState=is%wrap%NStateImp(1), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      call NUOPC_AddNamespaceWithNest(exportState, &
        nest="1", &
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
      call NUOPC_AddNamespaceWithNest(importState, &
        nest=trim(nStr), &
        nestedStateName="NestedStateImp_N"//trim(nStr), &
        nestedState=is%wrap%NStateImp(nIndex), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      call NUOPC_AddNamespaceWithNest(exportState, &
        nest=trim(nStr), &
        nestedStateName="NestedStateExp_N"//trim(nStr), &
        nestedState=is%wrap%NStateExp(nIndex), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    enddo

    !!
    !! advertise import and export fields
    !!
    do nIndex = 1, is%wrap%nnests
      ! Nest integer to string
      if ( nIndex > 999999999) then
        nStr = '999999999+'
      else
        write (nStr,"(I0)") nIndex
      endif
      do fIndex = 1, size(LIS_FieldList)
        if (LIS_FieldList(fIndex)%lisForc) then
          call NUOPC_Advertise(is%wrap%NStateImp(nIndex), &
            standardName=trim(LIS_FieldList(fIndex)%stdname), &
            name=trim(LIS_FieldList(fIndex)%stateName), &
            rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          if (is%wrap%verbosity >= VERBOSITY_MAX) then
            call ESMF_LogWrite( &
              trim(cname)//": "//METHOD//" Import Nest="//trim(nStr)//" Advertised Field="// &
              trim(LIS_FieldList(fIndex)%stdname),ESMF_LOGMSG_INFO)
          endif
        endif
      enddo
      do fIndex = 1, size(LIS_FieldList)
        if (LIS_FieldList(fIndex)%lisExport) then
          call NUOPC_Advertise(is%wrap%NStateExp(nIndex), &
            standardName=trim(LIS_FieldList(fIndex)%stdname), &
            name=trim(LIS_FieldList(fIndex)%stateName), &
            rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          if (is%wrap%verbosity >= VERBOSITY_MAX) then
            call ESMF_LogWrite( &
              trim(cname)//": "//METHOD//" Export Nest="//trim(nStr)//" Advertised Field="// &
              trim(LIS_FieldList(fIndex)%stdname),ESMF_LOGMSG_INFO)
          endif
        endif
      enddo
    enddo

    if (is%wrap%verbosity >= VERBOSITY_MAX) then
      call InternalStateLog(lisGridComp,trim(cname)//': '//METHOD//' Complete',rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "InitializeRealize"

  subroutine InitializeRealize(lisGridComp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: lisGridComp
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
    call ESMF_GridCompGet(lisGridComp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(lisGridComp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

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
      if (is%wrap%gridwrite_flag) then
        call NUOPC_FileWriteGrid(is%wrap%grids(nIndex), &
          'LIS_grid_nest_'//trim(nStr)//".nc", &
          nclMap=NUOPC_MAPPRESET_GLOBAL,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      do fIndex = 1, size(LIS_FieldList)
        if (LIS_FieldList(fIndex)%lisForc) then
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
            if (is%wrap%verbosity >= VERBOSITY_MAX) then
              call ESMF_LogWrite(trim(cname)//': '//METHOD//' Import Nest='//trim(nStr)// &
                ' Realized Field='//trim(LIS_FieldList(fIndex)%stateName), &
                ESMF_LOGMSG_INFO)
            endif
          else
            call ESMF_StateRemove(is%wrap%NStateImp(nIndex), (/trim(LIS_FieldList(fIndex)%stateName)/), &
              relaxedflag=.true.,rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          endif
        endif
      enddo
      do fIndex = 1, size(LIS_FieldList)
        if (LIS_FieldList(fIndex)%lisExport) then
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
            if (is%wrap%verbosity >= VERBOSITY_MAX) then
              call ESMF_LogWrite(trim(cname)//': '//METHOD//' Export Nest='//trim(nStr)// &
                ' Realized Field='//trim(LIS_FieldList(fIndex)%stateName), &
                ESMF_LOGMSG_INFO)
            endif
          else
            call ESMF_StateRemove(is%wrap%NStateExp(nIndex),(/trim(LIS_FieldList(fIndex)%stateName)/), &
              relaxedflag=.true.,rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          endif
        endif
      enddo
  
      call NUOPC_FillState(importState,value=0,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call NUOPC_FillState(exportState,value=0,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      is%wrap%modes(nIndex) = LIS_RunModeGet(LIS_FieldList,is%wrap%NStateImp(nIndex),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    enddo

    is%wrap%slice = 0

    if (is%wrap%verbosity >= VERBOSITY_DBG) then
      call LIS_FieldListLog(trim(cname)//': '//METHOD//' Complete',rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif
    if (is%wrap%verbosity >= VERBOSITY_MAX) then
      call InternalStateLog(lisGridComp,trim(cname)//': '//METHOD//' Complete',rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "DataInitialize"

  subroutine DataInitialize(lisGridComp, rc)
    type(ESMF_GridComp)  :: lisGridComp
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)     :: cname
    type(type_InternalState)   :: is
    type(ESMF_State)           :: exportState
    integer                    :: nIndex
    character(len=10)          :: nStr

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(lisGridComp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(lisGridComp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(lisGridComp, exportState=exportState, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    do nIndex=1,is%wrap%nnests
      ! Nest integer to string
      if (nIndex > 999999999) then
        nStr = '999999999+'
      else
        write (nStr,"(I0)") nIndex
      endif

      call LIS_NUOPC_DataInit(nest=nIndex,exportState=is%wrap%NStateExp(nIndex),rc=rc) 
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      if (is%wrap%statewrite_flag) then
        if ( nIndex > 999999999) then
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg="Maximum nest size for writing is 999,999,999.", &
            line=__LINE__,file=__FILE__,rcToReturn=rc)
         return  ! bail out
        endif
        call NUOPC_Write(is%wrap%NStateExp(nIndex), &
          fileNamePrefix="field_"//trim(cname)//"_export_nest_"//trim(nStr)//"_", &
          timeslice=is%wrap%slice, relaxedFlag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif
    enddo

    if (is%wrap%verbosity >= VERBOSITY_DBG) then
      call LIS_FieldListLog(trim(cname)//': '//METHOD//' Complete',rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif
    if (is%wrap%verbosity >= VERBOSITY_MAX) then
      call InternalStateLog(lisGridComp, &
        trim(cname)//': '//METHOD//' Complete',rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif
  end subroutine
  
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "SetClock"

  subroutine SetClock(lisGridComp, rc)
    type(ESMF_GridComp)  :: lisGridComp
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
    call ESMF_GridCompGet(lisGridComp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(lisGridComp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query the Component for its clock
    call NUOPC_ModelGet(lisGridComp, modelClock=modelClock, rc=rc)
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
    call NUOPC_CompSetClock(lisGridComp, modelClock, &
      modelTimestep, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (is%wrap%verbosity >= VERBOSITY_MAX) then
      call InternalStateLog(lisGridComp,trim(cname)//': '//METHOD//' Complete',rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "ModelAdvance"

  subroutine ModelAdvance(lisGridComp, rc)
    type(ESMF_GridComp)  :: lisGridComp
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
    call ESMF_GridCompGet(lisGridComp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(lisGridComp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(lisGridComp, modelClock=modelClock, importState=importState, &
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
    
    is%wrap%slice = is%wrap%slice + 1
    if (is%wrap%slice > 999999999) then
      sStr = '999999999+'
    else
      write (sStr,"(I0)") is%wrap%slice
    endif

    do nIndex=1,is%wrap%nnests
      ! Nest integer to string
      if (nIndex > 999999999) then
        nStr = '999999999+'
      else
        write (nStr,"(I0)") nIndex
      endif
      if (is%wrap%statewrite_flag) then
        if ( nIndex > 999999999) then
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg="Maximum nest size for writing is 999,999,999.", &
            line=__LINE__,file=__FILE__,rcToReturn=rc)
         return  ! bail out
        endif
        call NUOPC_Write(is%wrap%NStateImp(nIndex), &
          fileNamePrefix="field_"//trim(cname)// &
          "_import_nest_"//trim(nStr)//"_", &
          timeslice=is%wrap%slice, relaxedFlag=.true., rc=rc)
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
        call LIS_NUOPC_Run(nIndex,is%wrap%modes(nIndex),is%wrap%slice, &
          is%wrap%NStateImp(nIndex),is%wrap%NStateExp(nIndex), &
          is%wrap%clocks(nIndex), rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        call ESMF_ClockAdvance(is%wrap%clocks(nIndex),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        is%wrap%elapsedtimes(nIndex) = &
          is%wrap%elapsedtimes(nIndex) - nestTimestep
      enddo
      if (is%wrap%statewrite_flag) then
        if ( nIndex > 999999999) then
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg="Maximum nest size for writing is 999,999,999.", &
            line=__LINE__,file=__FILE__,rcToReturn=rc)
         return  ! bail out
        endif
        call NUOPC_Write(is%wrap%NStateExp(nIndex), &
          fileNamePrefix="field_"//trim(cname)//"_export_nest_"//trim(nStr)//"_", &
          timeslice=is%wrap%slice, relaxedFlag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif
    enddo

    if (is%wrap%verbosity >= VERBOSITY_DBG) then
      call LIS_FieldListLog(trim(cname)//': '//METHOD//' Complete Slice='//trim(sStr), &
        values=.TRUE.,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif
    if (is%wrap%verbosity >= VERBOSITY_MAX) then
      call InternalStateLog(lisGridComp, &
        trim(cname)//': '//METHOD//' Complete Slice='//trim(sStr),rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif
  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "ModelFinalize"

  subroutine ModelFinalize(lisGridComp, rc)
    type(ESMF_GridComp)  :: lisGridComp
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)     :: cname
    integer                    :: stat
    type(type_InternalState)   :: is
    integer                    :: nIndex
    type(ESMF_Clock)           :: clock

    rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(lisGridComp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(lisGridComp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call NUOPC_ModelGet(lisGridComp, modelClock=clock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! finalize the LIS model
    do nIndex=1,is%wrap%nnests
      call LIS_NUOPC_Final(nIndex,clock,rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    enddo

    if ( is%wrap%verbosity >= VERBOSITY_MAX ) then
      call ESMF_LogWrite(trim(cname)//": "//METHOD//" Complete", ESMF_LOGMSG_INFO)
    endif

    deallocate(is%wrap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='Deallocation of internal state memory failed.', &
      method=METHOD,file=FILENAME,rcToReturn=rc)) return ! bail out
    
  end subroutine

  !-----------------------------------------------------------------------------
  ! Utilities
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "InternalStateLog"

  subroutine InternalStateLog(lisGridComp,label,rc)
    type(ESMF_GridComp)                     :: lisGridComp
    character(len=*), intent(in), optional  :: label
    integer, intent(out),optional           :: rc

    ! local variables
    character(ESMF_MAXSTR)     :: cname
    character(len=64)          :: llabel
    type(type_InternalState)   :: is
    integer                    :: nIndex
    character(ESMF_MAXSTR)     :: logMsg
    type(ESMF_Time)            :: nestCurrTime
    type(ESMF_TimeInterval)    :: nestTimestep
    character(len=64)          :: nCurrTimeStr
    character(len=64)          :: nTimestepStr
    character(len=64)          :: nModeStr

    if(present(rc)) rc = ESMF_SUCCESS

    ! Query component for name
    call ESMF_GridCompGet(lisGridComp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if(present(label)) then
      llabel = trim(label)
    else
      llabel = trim(cname)//': '
    endif

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(lisGridComp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    write (logMsg, "(A,(A,I0))") trim(llabel), &
      ' Verbosity=',is%wrap%verbosity
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(llabel), &
      ' Grid Write=',is%wrap%gridwrite_flag
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(llabel), &
      ' State Write=',is%wrap%statewrite_flag
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,L1))") trim(llabel), &
      ' Profile Memory=',is%wrap%profile_memory
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,I0))") trim(llabel), &
      ' Nest Count=',is%wrap%nnests
    call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    write (logMsg, "(A,(A,I0))") trim(llabel), &
      ' Slice=',is%wrap%slice
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
      write (logMsg, "(A,(A,I0),(A,A))") trim(llabel), &
        " Nest=",nIndex, &
        " Mode=",trim(nModeStr)
      call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
      write (logMsg, "(A,(A,I0),(A,A))") trim(llabel), &
        " Nest=",nIndex, &
        " CurrentTime=",trim(nCurrTimeStr)
      call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
      write (logMsg, "(A,(A,I0),(A,A))") trim(llabel), &
        " Nest=",nIndex, &
        " Timestep=",trim(nTimestepStr)
      call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
    enddo
  end subroutine

end module

