!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#define FILENAME "LIS_NUOPC_Cap.F90"
#define MODNAME "LIS_NUOPC_Cap"
#include "LIS_NUOPC_Macros.h"

!> @file LIS_NUOPC_Cap.F90 LIS NUOPC Cap interfaces
module LIS_NUOPC

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
  use LIS_NUOPC_Flags
  use LIS_ESMF_Extensions

  implicit none

  private

  public SetServices

  CHARACTER(LEN=*), PARAMETER :: label_InternalState = 'InternalState'
  INTEGER, PARAMETER :: MAXNEST = 999999999

!> @brief Custom LIS NUOPC cap internal state
  type type_InternalStateStruct
    character(len=64)     :: configFile       = 'lis.config'
    logical               :: realizeAllExport = .FALSE.
    logical               :: nestToNest       = .FALSE.
    logical               :: cplEns           = .FALSE.
    type(field_init_flag) :: init_export      = FLD_INIT_FILLV
    type(field_init_flag) :: init_import      = FLD_INIT_FILLV
    type(missingval_flag) :: misg_import      = MISSINGVAL_FAIL
    character(len=40)     :: dirOutput        = "."
    integer               :: nnests           = 0
    integer               :: nfields          = size(LIS_FieldList)
    integer,allocatable                 :: ensMemberCnt(:)
    type(ESMF_Grid),allocatable         :: grids(:)
    type(ESMF_Clock),allocatable        :: clocks(:)
    type(ESMF_TimeInterval),allocatable :: stepAccum(:)
    type(ESMF_State),allocatable        :: NStateImp(:)
    type(ESMF_State),allocatable        :: NStateExp(:)
    integer,allocatable                 :: modes(:)
  end type

!> @cond IGNORE_WRAPPERS
  type type_InternalState
    type(type_InternalStateStruct), pointer :: wrap
  end type
!> @endcond

!EOP

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

!> @brief Register NUOPC compatible phases for initialize, run, and finalize
!! @param [inout] gcomp This component object
!! @param [out]   rc    Return value for subroutine
!! @details
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(type_InternalState)   :: is
    integer                    :: stat

    rc = ESMF_SUCCESS

    ! allocate memory for this internal state and set it in the component
    allocate(is%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg='Allocation of internal state memory failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    call ESMF_UserCompSetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeP3, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
       specRoutine=DataInitialize, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call NUOPC_CompSpecialize(gcomp, speclabel=model_label_SetClock, &
      specRoutine=SetClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_MethodRemove(gcomp, label=model_label_CheckImport, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
       specRoutine=CheckImport, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call NUOPC_CompSpecialize(gcomp, speclabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
      specRoutine=ModelFinalize, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

!> @brief Set the Initialize Phase Definition (IPD). Read model configuration
!! @param [inout] gcomp       This component object
!! @param [inout] importState The coupled import state
!! @param [inout] exportState The coupled export state
!! @param [inout] clock       The clock used for coupling
!! @param [out]   rc          Return value for subroutine
!! @details
!! During initialize phase 0 the runtime configuration is read in from model
!! attributes and the initialization phase definition version is set to
!! IPDv03.
  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    ! local variables
    character(32)              :: cname
    character(*), parameter    :: rname="InitializeP0"
    integer                    :: verbosity, diagnostic
    character(len=64)          :: value
    type(type_InternalState)   :: is
    integer                    :: stat

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call LIS_AttributeGet(rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! prepare diagnostics folder
    if (btest(diagnostic,16)) then
      call ESMF_UtilIOMkDir(pathName=trim(is%wrap%dirOutput), &
        relaxedFlag=.true., rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    endif

    contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!> @brief Read configuration from component's attributes
!! @param [out] rc Return value for subroutine
!! @details
    subroutine LIS_AttributeGet(rc)
      integer, intent(out)  :: rc

      ! local variables
      logical                    :: configIsPresent
      type(ESMF_Config)          :: config
      type(NUOPC_FreeFormat)     :: attrFF
      character(ESMF_MAXSTR)     :: logMsg

      ! check gcomp for config
      call ESMF_GridCompGet(gcomp, configIsPresent=configIsPresent, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      ! read and ingest free format component attributes
      if (configIsPresent) then
        call ESMF_GridCompGet(gcomp, config=config, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        attrFF = NUOPC_FreeFormatCreate(config, &
          label=trim(cname)//"_attributes::", relaxedflag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call NUOPC_CompAttributeIngest(gcomp, attrFF, addFlag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif

      ! Realize all export fields
      call ESMF_AttributeGet(gcomp, name="realize_all_export", value=value, &
        defaultValue="false", convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      value = ESMF_UtilStringLowerCase(value, rc=rc)
      is%wrap%realizeAllExport = (trim(value)=="true")

      ! Set configuration file name
      call ESMF_AttributeGet(gcomp, name="config_file", value=value, &
        defaultValue="lis.config", convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      is%wrap%configFile=value

      ! Turn on nest coupling
      call ESMF_AttributeGet(gcomp, name="nest_to_nest", value=value, &
        defaultValue="false", convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      value = ESMF_UtilStringLowerCase(value, rc=rc)
      is%wrap%nestToNest = (trim(value)=="true")

      ! Turn on ensemble coupling
      call ESMF_AttributeGet(gcomp, name="coupled_ensemble", value=value, &
        defaultValue="false", convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      value = ESMF_UtilStringLowerCase(value, rc=rc)
      is%wrap%cplEns = (trim(value)=="true")

      ! export data initialization type
      call ESMF_AttributeGet(gcomp, name="initialize_export", &
        value=value, defaultValue="FLD_INIT_MODEL", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      is%wrap%init_export = value

      ! Initialize import setting
      call ESMF_AttributeGet(gcomp, name="initialize_import", &
        value=value, defaultValue="FLD_INIT_FILLV", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      is%wrap%init_import = value

      ! Import dependency
      call ESMF_AttributeGet(gcomp, name="import_dependency", &
        value=value, defaultValue="false", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      value = ESMF_UtilStringLowerCase(value, rc=rc)
      if (trim(value)=="true") is%wrap%init_import = FLD_INIT_IMPORT

      ! Missing import value option
      call ESMF_AttributeGet(gcomp, name="missing_import", &
        value=value, defaultValue="MISSINGVAL_FAIL", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      is%wrap%misg_import = value

      ! Get component output directory
      call ESMF_AttributeGet(gcomp, name="output_directory", &
        value=value, defaultValue=trim(cname)//"_OUTPUT", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      is%wrap%dirOutput = trim(value)

      if (btest(verbosity,16)) then
        call ESMF_LogWrite(trim(cname)//": Settings",ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'Verbosity              = ',verbosity
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'Diagnostic             = ',diagnostic
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,L1))") trim(cname)//': ', &
          'Nest To Nest           = ',is%wrap%nestToNest
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,L1))") trim(cname)//': ', &
          'Coupled Ensemble       = ',is%wrap%cplEns
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        value = is%wrap%init_export
        write (logMsg, "(A,(A,A))") trim(cname)//': ', &
          'Initialize Export      = ',trim(value)
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        value = is%wrap%init_import
        write (logMsg, "(A,(A,A))") trim(cname)//': ', &
          'Initialize Import      = ',trim(value)
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        value = is%wrap%misg_import
        write (logMsg, "(A,(A,A))") trim(cname)//': ', &
          'Missing Imports        = ',trim(value)
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,L1))") trim(cname)//': ', &
          'Realze All Exports     = ',is%wrap%realizeAllExport
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,A))") trim(cname)//': ', &
          'Config File            = ',is%wrap%configFile
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,A))") trim(cname)//': ', &
          'Output Directory       = ',is%wrap%dirOutput
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
      endif

    end subroutine

  end subroutine

  !-----------------------------------------------------------------------------

!> @brief Initialize internal model. Advertise import and export fields
!! @param [inout] gcomp       This component object
!! @param [inout] importState The coupled import state
!! @param [inout] exportState The coupled export state
!! @param [inout] clock       The clock used for coupling
!! @param [out]   rc          Return value for subroutine
!! @details
!! During initialize phase 1 the model is initialized and the import and
!! export fields are advertised in nested import and export states. Import
!! fields are configured in the forcing variables list file.
  subroutine InitializeP1(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)     :: gcomp
    type(ESMF_State)        :: importState, exportState
    type(ESMF_Clock)        :: clock
    integer,intent(out)     :: rc

    ! LOCAL VARIABLES
    character(32)              :: cname
    character(*), parameter    :: rname="InitializeP1"
    integer                    :: verbosity, diagnostic
    character(len=64)          :: value
    type(type_InternalState)   :: is
    integer                    :: stat
    type(ESMF_VM)              :: vm
    integer                    :: localPet, petCount
    integer                    :: fIndex
    integer                    :: nIndex
    character(len=9)           :: nStr

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_GridCompGet(gcomp, vm=vm, localPet=localPet, petCount=petCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! initialize lis model for this PET
    call LIS_NUOPC_Init(vm, configFile=is%wrap%configFile, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (btest(verbosity,16)) then
      call LIS_Log(trim(cname)//': '//rname,rc)
      if (ESMF_STDERRORCHECK(rc)) return
    endif

    call LIS_FieldDictionaryAdd(rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    is%wrap%nnests = LIS_NestCntGet(rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! Max nest check
    if ( is%wrap%nnests > MAXNEST ) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg="Maximum nest size is 999,999,999.", &
        line=__LINE__,file=__FILE__,rcToReturn=rc)
      return
    endif

    allocate( &
      is%wrap%ensMemberCnt(is%wrap%nnests), &
      is%wrap%grids(is%wrap%nnests), &
      is%wrap%clocks(is%wrap%nnests), &
      is%wrap%stepAccum(is%wrap%nnests), &
      is%wrap%NStateImp(is%wrap%nnests), &
      is%wrap%NStateExp(is%wrap%nnests), &
      is%wrap%modes(is%wrap%nnests), &
      stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of internal state nest memory failed.", &
      line=__LINE__,file=__FILE__)) &
      return

    is%wrap%modes=LIS_Unknown

#ifdef GSM_EXTLND
    is%wrap%NStateImp(1) = importState
    is%wrap%NStateExp(1) = exportState
#else
    if (.NOT.is%wrap%nestToNest) then
      if (is%wrap%nnests.le.1) then
        is%wrap%NStateImp(1) = importState
        is%wrap%NStateExp(1) = exportState
      else
        call ESMF_LogSetError(ESMF_FAILURE, &
          msg="Nest to nest must be turned on when multiple domains exist.", &
          line=__LINE__,file=__FILE__,rcToReturn=rc)
        return
      endif
    else
      ! add namespace
      call NUOPC_AddNestedState(importState, &
        CplSet="1", &
        nestedStateName="NestedStateImp_N1", &
        nestedState=is%wrap%NStateImp(1), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call NUOPC_AddNestedState(exportState, &
        CplSet="1", &
        nestedStateName="NestedStateExp_N1", &
        nestedState=is%wrap%NStateExp(1), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    endif

    do nIndex = 2, is%wrap%nnests
      write (nStr,"(I0)") nIndex
      call NUOPC_AddNestedState(importState, &
        CplSet=trim(nStr), &
        nestedStateName="NestedStateImp_N"//trim(nStr), &
        nestedState=is%wrap%NStateImp(nIndex), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call NUOPC_AddNestedState(exportState, &
        CplSet=trim(nStr), &
        nestedStateName="NestedStateExp_N"//trim(nStr), &
        nestedState=is%wrap%NStateExp(nIndex), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    enddo
#endif

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
          if (ESMF_STDERRORCHECK(rc)) return
       endif
       if (LIS_FieldList(fIndex)%adExport) then
         call NUOPC_Advertise(is%wrap%NStateExp(nIndex), &
           standardName=trim(LIS_FieldList(fIndex)%stdname), &
           name=trim(LIS_FieldList(fIndex)%stateName), &
           rc=rc)
         if (ESMF_STDERRORCHECK(rc)) return
       endif
      enddo
    enddo

    if (btest(verbosity,16)) call LogAdvertised()

    contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!> @brief Log advertised import and export fields
!! @details
    subroutine LogAdvertised()
      ! local variables
      integer                    :: cntImp
      integer                    :: cntExp
      integer                    :: fIndex
      character(ESMF_MAXSTR)     :: logMsg

      ! Count advertised import and export fields
      cntImp = 0
      cntExp = 0
      do fIndex = 1, size(LIS_FieldList)
        if (LIS_FieldList(fIndex)%adImport) cntImp = cntImp + 1
        if (LIS_FieldList(fIndex)%adExport) cntExp = cntExp + 1
      enddo

      ! Report advertised import fields
      write(logMsg,'(a,a,i0,a)') TRIM(cname)//': ', &
        'List of advertised import fields(',cntImp,'):'
      call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
      write(logMsg,'(a,a5,a,a16,a,a)') TRIM(cname)//': ', &
        'index',' ','name',' ','standardName'
      call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
      cntImp = 0
      do fIndex=1, size(LIS_FieldList)
        if (.NOT.LIS_FieldList(fIndex)%adImport) cycle
        cntImp = cntImp + 1
        write(logMsg,'(a,i5,a,a16,a,a)') TRIM(cname)//': ', &
          cntImp,' ',TRIM(LIS_FieldList(fIndex)%stateName), &
          ' ',TRIM(LIS_FieldList(fIndex)%stdName)
        call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
      enddo

      ! Report advertised export fields
      write(logMsg,'(a,a,i0,a)') TRIM(cname)//': ', &
        'List of advertised export fields(',cntExp,'):'
      call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
      write(logMsg,'(a,a5,a,a16,a,a)') TRIM(cname)//': ', &
        'index',' ','name',' ','standardName'
      call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
      cntExp = 0
      do fIndex=1, size(LIS_FieldList)
        if (.NOT.LIS_FieldList(fIndex)%adExport) cycle
        cntExp = cntExp + 1
        write(logMsg,'(a,i5,a,a16,a,a)') TRIM(cname)//': ', &
          cntExp,' ',TRIM(LIS_FieldList(fIndex)%stateName), &
          ' ',TRIM(LIS_FieldList(fIndex)%stdName)
        call ESMF_LogWrite(trim(LogMsg), ESMF_LOGMSG_INFO)
      enddo

    end subroutine

  end subroutine

  !-----------------------------------------------------------------------------

!> @brief Realize import and export fields
!! @param [inout] gcomp       This component object
!! @param [inout] importState The coupled import state
!! @param [inout] exportState The coupled export state
!! @param [inout] clock       The clock used for coupling
!! @param [out]   rc          Return value for subroutine
!! @details
!! During initialize phase 3 import and export fields are realized in each
!! nested import and export state if they are connected through NUOPC.
!! Realized fields are created on the LIS grid. All export fields are realized
!! if realize all export fields is turned on.
  subroutine InitializeP3(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    character(32)              :: cname
    character(*), parameter    :: rname="InitializeP3"
    integer                    :: verbosity, diagnostic
    character(len=64)          :: value
    type(type_InternalState)   :: is
    integer                    :: nIndex
    type(ESMF_Field)           :: field
    integer                    :: fIndex
    character(len=9)           :: nStr
    logical                    :: realizeImp, realizeExp
    type(ESMF_Array)           :: array

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do nIndex = 1, is%wrap%nnests
      write (nStr,"(I0)") nIndex

      ! count ensemble members
      if (is%wrap%cplEns) then
        call LIS_EnsMemberCntGet(nIndex, &
          is%wrap%ensMemberCnt(nIndex), rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      else
        is%wrap%ensMemberCnt(nIndex) = 0
      endif

      ! Call gluecode to create grid.
      is%wrap%grids(nIndex) = LIS_GridCreate(nIndex, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      if (btest(verbosity,16)) then
        call LIS_ESMF_LogGrid(is%wrap%grids(nIndex), &
          trim(cname)//"_"//rname//"_D"//trim(nStr),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif

      ! Write grid to NetCDF file.
      if (btest(diagnostic,16)) then
        call LIS_ESMF_GridWrite(is%wrap%grids(nIndex), &
          trim(is%wrap%dirOutput)//"/diag_"//trim(cname)//"_"// &
          rname//'_grid_D'//trim(nStr)//".nc", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif

      do fIndex = 1, size(LIS_FieldList)
        realizeImp = .FALSE.
        realizeExp = .FALSE.
        if (LIS_FieldList(fIndex)%adExport) then
          if (is%wrap%realizeAllExport) then
            realizeExp = .TRUE.
          else
            realizeExp = NUOPC_IsConnected(is%wrap%NStateExp(nIndex), &
              fieldName=trim(LIS_FieldList(fIndex)%stateName),rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          endif
        endif
        if (LIS_FieldList(fIndex)%adImport) then
          realizeImp = NUOPC_IsConnected(is%wrap%NStateImp(nIndex), &
            fieldName=trim(LIS_FieldList(fIndex)%stateName),rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        endif

        if (realizeImp .AND. realizeExp .AND. LIS_FieldList(fIndex)%sharedMem) then
          if (is%wrap%cplEns) then
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, &
              ungriddedLBound=(/1/), ungriddedUBound=(/is%wrap%ensMemberCnt(nIndex)/), rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          else
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          endif
          call NUOPC_Realize(is%wrap%NStateExp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedExport = .TRUE.
          call ESMF_FieldGet(field=field,array=array,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          if (is%wrap%cplEns) then
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), array=array, &
              ungriddedLBound=(/1/), ungriddedUBound=(/is%wrap%ensMemberCnt(nIndex)/), rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          else
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), array=array, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          endif
          call NUOPC_Realize(is%wrap%NStateImp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedImport = .TRUE.
        elseif (realizeImp .AND. realizeExp) then
          if (is%wrap%cplEns) then
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, &
              ungriddedLBound=(/1/), ungriddedUBound=(/is%wrap%ensMemberCnt(nIndex)/), rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          else
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          endif
          call NUOPC_Realize(is%wrap%NStateExp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedExport = .TRUE.
          if (is%wrap%cplEns) then
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, &
              ungriddedLBound=(/1/), ungriddedUBound=(/is%wrap%ensMemberCnt(nIndex)/), rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          else
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, rc=rc)
            if(ESMF_STDERRORCHECK(rc)) return
          endif
          call NUOPC_Realize(is%wrap%NStateImp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedImport = .TRUE.
        elseif (realizeExp) then
          if (is%wrap%cplEns) then
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, &
              ungriddedLBound=(/1/), ungriddedUBound=(/is%wrap%ensMemberCnt(nIndex)/), rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          else
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          endif
          call NUOPC_Realize(is%wrap%NStateExp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedExport = .TRUE.
          call ESMF_StateRemove(is%wrap%NStateImp(nIndex),(/trim(LIS_FieldList(fIndex)%stateName)/), &
            relaxedflag=.true.,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        elseif (realizeImp) then
          if (is%wrap%cplEns) then
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, &
              ungriddedLBound=(/1/), ungriddedUBound=(/is%wrap%ensMemberCnt(nIndex)/), rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          else
            field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
              grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          endif
          call NUOPC_Realize(is%wrap%NStateImp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedImport = .TRUE.
          call ESMF_StateRemove(is%wrap%NStateExp(nIndex),(/trim(LIS_FieldList(fIndex)%stateName)/), &
            relaxedflag=.true.,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        else
          call ESMF_StateRemove(is%wrap%NStateExp(nIndex),(/trim(LIS_FieldList(fIndex)%stateName)/), &
            relaxedflag=.true.,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call ESMF_StateRemove(is%wrap%NStateImp(nIndex),(/trim(LIS_FieldList(fIndex)%stateName)/), &
            relaxedflag=.true.,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        endif
      enddo

      if(.not.(is%wrap%init_import .eq. FLD_INIT_MODEL)) then
        call LIS_ESMF_FillState(is%wrap%NStateImp(nIndex), MISSINGVALUE, &
          rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif

      if(.not.(is%wrap%init_export .eq. FLD_INIT_MODEL)) then
        call LIS_ESMF_FillState(is%wrap%NStateExp(nIndex), MISSINGVALUE, &
          rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif

      is%wrap%modes(nIndex) = LIS_RunModeGet(LIS_FieldList,is%wrap%NStateImp(nIndex),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

    enddo

    if (btest(verbosity,16)) call LogRealized()

    contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!> @brief Log realized import and export fields
!! @details
    subroutine LogRealized()
      ! local variables
      integer                    :: cntImp
      integer                    :: cntExp
      integer                    :: fIndex
      character(ESMF_MAXSTR)     :: logMsg

      ! Count advertised import and export fields
      cntImp = 0
      cntExp = 0
      do fIndex = 1, size(LIS_FieldList)
        if (LIS_FieldList(fIndex)%realizedImport) cntImp = cntImp + 1
        if (LIS_FieldList(fIndex)%realizedExport) cntExp = cntExp + 1
      enddo

      ! Report realized import fields
      write(logMsg,'(a,a,i0,a)') TRIM(cname)//': ', &
        'List of realized import fields(',cntImp,'):'
      call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
      write(logMsg,'(a,a5,a,a16,a,a)') TRIM(cname)//': ', &
        'index',' ','name',' ','standardName'
      call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
      cntImp = 0
      do fIndex=1, size(LIS_FieldList)
        if (.NOT.LIS_FieldList(fIndex)%realizedImport) cycle
        cntImp = cntImp + 1
        write(logMsg,'(a,i5,a,a16,a,a)') TRIM(cname)//': ', &
          cntImp,' ',TRIM(LIS_FieldList(fIndex)%stateName), &
          ' ',TRIM(LIS_FieldList(fIndex)%stdName)
        call ESMF_LogWrite(trim(LogMsg), ESMF_LOGMSG_INFO)
      enddo

      ! Report realized export fields
      write(logMsg,'(a,a,i0,a)') TRIM(cname)//': ', &
        'List of realized export fields(',cntExp,'):'
      call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
      write(logMsg,'(a,a5,a,a16,a,a)') TRIM(cname)//': ', &
        'index',' ','name',' ','standardName'
      call ESMF_LogWrite(TRIM(logMsg), ESMF_LOGMSG_INFO)
      cntExp = 0
      do fIndex=1, size(LIS_FieldList)
        if (.NOT.LIS_FieldList(fIndex)%realizedExport) cycle
        cntExp = cntExp + 1
        write(logMsg,'(a,i5,a,a16,a,a)') TRIM(cname)//': ', &
          cntExp,' ',TRIM(LIS_FieldList(fIndex)%stateName), &
          ' ',TRIM(LIS_FieldList(fIndex)%stdName)
        call ESMF_LogWrite(trim(LogMsg), ESMF_LOGMSG_INFO)
      enddo

    end subroutine

  end subroutine

  !-----------------------------------------------------------------------------

!> @brief Initialize import and export data
!! @param [inout] gcomp This component object
!! @param [out]   rc    Return value for subroutine
!! @details
!! During data initialize this cap checks the timestamp of all import fields
!! dependent on a coupled model.  Once all dependent import fields have been
!! initialized this cap is marked initalized. The export fields are updated
!! and initialized regardless of the import field dependencies.
  subroutine DataInitialize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(32)                          :: cname
    character(*), parameter                :: rname="DataInitialize"
    integer                                :: verbosity, diagnostic
    character(len=64)                      :: value
    type(type_InternalState)               :: is
    integer                                :: stat
    type(ESMF_Clock)                       :: modelClock
    type(ESMF_Time)                        :: currTime
    character(len=32)                      :: currTimeStr
    integer                                :: nIndex
    character(len=9)                       :: nStr
    integer                                :: iIndex
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: ifield, efield
    integer                                :: eSearch
    logical                                :: importCurrent
    logical                                :: importUpdated

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query the Component for its clock
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! get the current time out of the clock
    call ESMF_ClockGet(modelClock, currTime=currTime, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_TimeGet(currTime, timeString=currTimeStr, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    importUpdated = .TRUE.
    do nIndex=1,is%wrap%nnests
      write (nStr,"(I0)") nIndex

      ! Initialize export fields
      call LIS_NUOPC_DataInit(nest=nIndex, &
        exportState=is%wrap%NStateExp(nIndex),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_StateGet(is%wrap%NStateExp(nIndex),itemCount=itemCount, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      allocate( &
        itemNameList(itemCount), &
        itemTypeList(itemCount), &
        stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Allocation of state item list memory failed.", &
        line=__LINE__, &
        file=__FILE__)) &
        return
      call ESMF_StateGet(is%wrap%NStateExp(nIndex),itemNameList=itemNameList, &
        itemTypeList=itemTypeList,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      do iIndex=1, itemCount
        if (itemTypeList(iIndex) == ESMF_STATEITEM_FIELD) then
          call ESMF_StateGet(is%wrap%NStateExp(nIndex),field=efield, &
            itemName=itemNameList(iIndex),rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call NUOPC_SetAttribute(efield, name="Updated", value="true", rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        endif
      enddo
      deallocate(itemNameList, itemTypeList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg="Deallocation of state item list memory failed.", &
        line=__LINE__, &
        file=__FILE__)) &
        return

      if (is%wrap%init_import .eq. FLD_INIT_IMPORT) then
        ! Check data dependencies
        importCurrent = NUOPC_IsAtTime(is%wrap%NStateImp(nIndex), &
          time=currTime, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        if (importCurrent) then
          call ESMF_LogWrite( &
            trim(cname)//': '//rname//' Initialize-Data-Dependency SATISFIED!!! Nest='//trim(nStr), &
            ESMF_LOGMSG_INFO)
          call LIS_ImportFieldsCopy(nIndex,is%wrap%NStateImp(nIndex),is%wrap%misg_import,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        else
          call ESMF_LogWrite( &
            trim(cname)//': '//rname//' Initialize-Data-Dependency NOT YET SATISFIED!!! Nest='//trim(nStr), &
            ESMF_LOGMSG_INFO)
          importUpdated = .FALSE.
        endif
      else
        ! Reset all import fields to export or fillv
        call ESMF_StateGet(is%wrap%NStateImp(nIndex),itemCount=itemCount, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        allocate( &
          itemNameList(itemCount), &
          itemTypeList(itemCount), &
          stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Allocation of state item list memory failed.", &
          line=__LINE__, &
          file=__FILE__)) &
          return

        call ESMF_StateGet(is%wrap%NStateImp(nIndex),itemNameList=itemNameList, &
          itemTypeList=itemTypeList,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        do iIndex=1, itemCount
          if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELD) then
            call ESMF_StateGet(is%wrap%NStateImp(nIndex),field=ifield, &
              itemName=itemNameList(iIndex),rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
            call ESMF_StateGet(is%wrap%NStateExp(nIndex), &
              itemSearch=itemNameList(iIndex), itemCount=eSearch, rc=rc)
            if (eSearch .gt. 0) then
              call ESMF_StateGet(is%wrap%NStateExp(nIndex),field=efield, &
                itemName=itemNameList(iIndex),rc=rc)
              if (ESMF_STDERRORCHECK(rc)) return
              call ESMF_FieldCopy(ifield, fieldIn=efield, rc=rc)
              if (ESMF_STDERRORCHECK(rc)) return
            else
              if (is%wrap%init_import .eq. FLD_INIT_ZERO) then
                call ESMF_FieldFill(ifield, dataFillScheme="const", &
                  const1=0.0_ESMF_KIND_R8, rc=rc)
                if (ESMF_STDERRORCHECK(rc)) return
              elseif (is%wrap%init_import .eq. FLD_INIT_FILLV) then
                call ESMF_FieldFill(ifield, dataFillScheme="const", &
                  const1=MISSINGVALUE, rc=rc)
                if (ESMF_STDERRORCHECK(rc)) return
              else
                call ESMF_LogSetError(ESMF_FAILURE, &
                  msg="Invalid import initialization option.", &
                  line=__LINE__,file=__FILE__,rcToReturn=rc)
                return
              endif
            endif
          endif
        enddo

        deallocate(itemNameList, itemTypeList, stat=stat)
        if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
          msg="Deallocation of state item list memory failed.", &
          line=__LINE__, &
          file=__FILE__)) &
          return
      endif

      if(is%wrap%init_export .eq. FLD_INIT_FILLV) then
        call LIS_ESMF_FillState(is%wrap%NStateExp(nIndex), MISSINGVALUE, &
          rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      elseif(is%wrap%init_export .eq. FLD_INIT_ZERO) then
        call LIS_ESMF_FillState(is%wrap%NStateExp(nIndex), 0.0_ESMF_KIND_R8, &
          rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      elseif(is%wrap%init_export .eq. FLD_INIT_MODEL) then
        ! do nothing
      else
        call ESMF_LogSetError(ESMF_FAILURE, &
          msg="Invalid export initialization option.", &
          line=__LINE__,file=__FILE__,rcToReturn=rc)
        return
      endif

    enddo ! enddo nnests

    ! set InitializeDataComplete Attribute to "true", indicating to the
    ! generic code that all inter-model data dependencies are satisfied
    if (importUpdated) then
      call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      ! Write initialization files
      if (btest(diagnostic,16)) then
        do nIndex=1,is%wrap%nnests
          write (nStr,"(I0)") nIndex
          call NUOPC_Write(is%wrap%NStateImp(nIndex), &
            fileNamePrefix=trim(is%wrap%dirOutput)//"/diag_"//trim(cname)//"_"// &
              rname//"_imp_D"//trim(nStr)//"_"//trim(currTimeStr)//"_", &
            overwrite=.true., status=ESMF_FILESTATUS_REPLACE, timeslice=1, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call NUOPC_Write(is%wrap%NStateExp(nIndex), &
            fileNamePrefix=trim(is%wrap%dirOutput)//"/diag_"//trim(cname)//"_"// &
              rname//"_exp_D"//trim(nStr)//"_"//trim(currTimeStr)//"_", &
            overwrite=.true., status=ESMF_FILESTATUS_REPLACE, timeslice=1, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        enddo
      endif
    endif

  end subroutine

  !-----------------------------------------------------------------------------

!> @brief Set model clock during initialization
!! @param [inout] gcomp This component object
!! @param [out]   rc    Return value for subroutine
!! @details
!! During set clock the cap creates a new clock for each nest. The time step
!! for each nest is set in LIS configuration file and initialized during LIS
!! initialization. The time accumulation tracker for each timestep is reset to
!! zero.  The cap's time step is updated to the shortest time step
!! of all nests. The restart write time step is also created and the restart
!! write time accumulation tracker is reset to zero.
  subroutine SetClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(32)              :: cname
    character(*), parameter    :: rname="SetClock"
    integer                    :: verbosity, diagnostic
    character(len=64)          :: value
    type(type_InternalState)   :: is
    integer                    :: nIndex
    real(ESMF_KIND_R8)         :: mindt
    real(ESMF_KIND_R8)         :: ndt
    type(ESMF_Clock)           :: modelClock
    type(ESMF_TimeInterval)    :: modelTimestep
    type(ESMF_TimeInterval)    :: nestTimeStep

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query the Component for its clock
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! Set minTimestep to the timestep of the first nest
    mindt = LIS_TimestepGet(nest=1,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do nIndex = 1, is%wrap%nnests
      is%wrap%clocks(nIndex) = ESMF_ClockCreate(modelClock, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      ndt = LIS_TimestepGet(nest=nIndex,rc=rc)
      call ESMF_TimeIntervalSet(nestTimestep, &
        s_r8=ndt, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_ClockSet(is%wrap%clocks(nIndex), &
        timeStep=nestTimestep, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      if (ndt < mindt) mindt = ndt

      call ESMF_TimeIntervalSet(is%wrap%stepAccum(nIndex), &
        s_r8=0._ESMF_KIND_R8, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    enddo

    call ESMF_TimeIntervalSet(modelTimestep, &
      s_r8=mindt, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call NUOPC_CompSetClock(gcomp, modelClock, &
      modelTimestep, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

!> @brief Check timestamp on import data
!! @param [inout] gcomp This component object
!! @param [out]   rc    Return value for subroutine
!! @details
!! During check import the import data is checked to verify that it is at
!! the beginning or end of the timestep.
  subroutine CheckImport(gcomp, rc)
    type(ESMF_GridComp) :: gcomp
    integer,intent(out) :: rc

    ! local variables
    character(32)               :: cname
    character(*), parameter     :: rname="CheckImport"
    integer                     :: verbosity, diagnostic
    character(len=64)           :: value
    type(type_InternalState)    :: is
    integer                     :: nIndex
    character(len=9)            :: nStr
    type(ESMF_Clock)            :: modelClock
    type(ESMF_Time)             :: modelCurrTime
    type(ESMF_Time)             :: modelStartTime
    type(ESMF_Field)            :: field
    type(ESMF_StateItem_Flag)   :: itemType
    logical                     :: allCurrTime
    logical                     :: fieldCurrTime
    integer                     :: fIndex

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query the component for its clock, importState, and exportState
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! get the stop time out of the clock
    call ESMF_ClockGet(modelClock, startTime=modelStartTime, &
      currTime=modelCurrTime, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do nIndex=1,is%wrap%nnests
      write (nStr,"(I0)") nIndex
      allCurrTime = NUOPC_IsAtTime(is%wrap%NStateImp(nIndex), modelCurrTime,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      if (.NOT.allCurrTime) then
        call ESMF_LogWrite(trim(cname)//": NUOPC INCOMPATIBILITY DETECTED: "// &
          "Import Fields Nest="//trim(nStr)//" not at correct time", &
          ESMF_LOGMSG_WARNING)
      endif
    enddo

  end subroutine

  !-----------------------------------------------------------------------------

!> @brief Advances the model by a timestep
!! @param [inout] gcomp This component object
!! @param [out]   rc    Return value for subroutine
!! @details
!! During model advance each nest time accumulation tracker is increased by
!! the timestep of the cap.  If the time accumlation tracker is greater than
!! the time step of the nest then the nest is advanced.
  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(32)               :: cname
    character(*), parameter     :: rname="ModelAdvance"
    integer                     :: verbosity, diagnostic
    character(len=64)           :: value
    type(type_InternalState)    :: is
    integer                     :: nIndex
    character(len=9)            :: nStr
    type(ESMF_Clock)            :: modelClock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: currTime, advEndTime
    character(len=32)           :: currTimeStr, advEndTimeStr
    type(ESMF_TimeInterval)     :: timeStep
    type(ESMF_TimeInterval)     :: nestTimeStep

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query the component for its clock, importState, and exportState
    call NUOPC_ModelGet(gcomp, &
      modelClock=modelClock, &
      importState=importState, &
      exportState=exportState, &
      rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query the clock for the current time and time step
    call ESMF_ClockGet(modelClock, &
      currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    advEndTime = currTime + timeStep
    call ESMF_TimeGet(currTime, timeString=currTimeStr, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_TimeGet(advEndTime, timeString=advEndTimeStr, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.

    do nIndex=1,is%wrap%nnests
      write (nStr,"(I0)") nIndex
      ! Write import files
      if (btest(diagnostic,16)) then
        call NUOPC_Write(is%wrap%NStateImp(nIndex), &
          fileNamePrefix=trim(is%wrap%dirOutput)//"/diag_"//trim(cname)//"_"// &
            rname//"_imp_D"//trim(nStr)//"_"//trim(currTimeStr)//"_", &
          overwrite=.true., status=ESMF_FILESTATUS_REPLACE, timeslice=1, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif

      is%wrap%stepAccum(nIndex) = is%wrap%stepAccum(nIndex) + timeStep

      call ESMF_ClockGet(is%wrap%clocks(nIndex),timeStep=nestTimestep,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      do while (is%wrap%stepAccum(nIndex) >= nestTimestep)
        ! Gluecode NestAdvance
        if (btest(verbosity,16)) then
          call LogAdvance(nIndex)
        endif
        call LIS_NUOPC_Run(nIndex,is%wrap%modes(nIndex), &
          is%wrap%NStateImp(nIndex),is%wrap%NStateExp(nIndex), &
          is%wrap%clocks(nIndex), is%wrap%misg_import, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call ESMF_ClockAdvance(is%wrap%clocks(nIndex),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        is%wrap%stepAccum(nIndex) = &
          is%wrap%stepAccum(nIndex) - nestTimestep
      enddo

      ! Write export files
      if (btest(diagnostic,16)) then
        call NUOPC_Write(is%wrap%NStateExp(nIndex), &
          fileNamePrefix=trim(is%wrap%dirOutput)//"/diag_"//trim(cname)//"_"// &
            rname//"_exp_D"//trim(nStr)//"_"//trim(advEndTimeStr)//"_", &
          overwrite=.true., status=ESMF_FILESTATUS_REPLACE, timeslice=1, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
    enddo

    contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!> @brief Log model advance call
!! @param [inout] nIndex Index of nest
!! @details
    subroutine LogAdvance(nIndex)
      integer                    :: nIndex
      ! local variables
      character(ESMF_MAXSTR)     :: logMsg
      character(len=64)          :: nModeStr
      type(ESMF_Time)            :: nestCurrTime
      type(ESMF_TimeInterval)    :: nestTimestep
      character(len=64)          :: nCurrTimeStr
      character(len=64)          :: nTimestepStr
      integer                    :: rc

      call ESMF_LogWrite(trim(cname)//': '//rname//&
        ' Advancing Nest='//trim(nStr),ESMF_LOGMSG_INFO)

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
      write (logMsg, "(A,(A,I0,A),(A,A))") trim(cname)//': ', &
        'Nest(',nIndex,') ', &
        'Mode = ',trim(nModeStr)
      call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)

      if (allocated(is%wrap%clocks)) then
        if (ESMF_ClockIsCreated(is%wrap%clocks(nIndex))) then
          call ESMF_ClockGet(is%wrap%clocks(nIndex), &
            currTime=nestCurrTime,timeStep=nestTimestep,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call ESMF_TimeGet(nestCurrTime, &
            timeString=nCurrTimeStr,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call ESMF_TimeIntervalGet(nestTimestep, &
            timeString=nTimestepStr,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        else
          nCurrTimeStr = "(not_created)"
          nTimestepStr = "(not_created)"
        endif
      else
        nCurrTimeStr = "(unallocated)"
        nTimestepStr = "(unallocated)"
      endif
      write (logMsg, "(A,(A,I0,A),(A,A))") trim(cname)//": ", &
        "Nest(",nIndex,") ", &
        "Current Time = ",trim(nCurrTimeStr)
      call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
      write (logMsg, "(A,(A,I0,A),(A,A))") trim(cname)//": ", &
        "Nest(",nIndex,") ", &
        "Time Step    = ",trim(nTimestepStr)
      call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)

    end subroutine

  end subroutine

  !-----------------------------------------------------------------------------

!> @brief Releases memory
!! @param [inout] gcomp This component object
!! @param [out]   rc    Return value for subroutine
!! @details
!! During model finalize LIS finalize subroutines are called and memory
!! allocated during cap initialization is released. see ::modelfinalize
  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(32)              :: cname
    character(*), parameter    :: rname="ModelFinalize"
    integer                    :: verbosity, diagnostic
    character(len=64)          :: value
    type(type_InternalState)   :: is
    integer                    :: stat
    integer                    :: nIndex
    character(len=9)           :: nStr
    type(ESMF_Clock)           :: modelClock
    type(ESMF_Time)            :: currTime
    character(len=32)          :: currTimeStr

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! query the Component for its clock
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_ClockGet(modelClock, currTime=currTime, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_TimeGet(currTime, timeString=currTimeStr, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! finalize the LIS model
    do nIndex=1,is%wrap%nnests
      call LIS_NUOPC_Final(nIndex,modelClock,rc)
      if (ESMF_STDERRORCHECK(rc)) return
    enddo

    deallocate(is%wrap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='Deallocation of internal state memory failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return

  end subroutine

end module

