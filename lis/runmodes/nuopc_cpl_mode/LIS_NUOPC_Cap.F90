!>
!! @mainpage NASA's Land Information System (LIS) NUOPC Cap
!! @author Daniel Rosen (daniel.rosen@noaa.gov)
!! @author ESMF Support (esmf_support@ucar.edu)
!! @date 03/14/2017 LIS NUOPC Cap Added to GitHub
!! @date 03/15/2017 Documentation Added
!!
!! @tableofcontents
!!
!! @section Overview Overview
!!
!! The Land Information System (LIS) is a surface model wrapper developed
!! and maintained by the National Aeronautics and Space Administration (NASA).
!! LIS abstracts several land and ocean surface models with standard
!! interfaces. The LIS cap wraps LIS with NUOPC compliant interfaces. The
!! result is configurable surface model capable of coupling with other models
!! in the National Unified Operational Prediction Capability (NUOPC).  As of
!! LIS version 7.1 only Noah v.3.3 has two-way NUOPC coupling capabilities.
!!
!! This page documents the technical design of the specialized NUOPC model and
!! the LIS gluecode.  For generic NUOPC model documentation please see the
!! [NUOPC Reference Manual] (https://www.earthsystemcog.org/projects/nuopc/refmans).
!!
!!
!! @section NuopcSpecialization NUOPC Model Specialized Entry Points
!!
!! This cap specializes the cap configuration, initialization, advertised
!! fields, realized fields, data initialization, clock, run, and finalize.
!!
!! @subsection SetServices Set Services (Register Subroutines)
!!
!! Table summarizing the NUOPC specialized subroutines registered during
!! [SetServices] (@ref LIS_NUOPC::SetServices).  The "Phase" column says
!! whether the subroutine is called during the initialization, run, or
!! finalize part of the coupled system run.
!!
!! Phase  |     Cap Subroutine                                | Description
!! -------|---------------------------------------------------|-------------------------------------------------------------
!! Init   | [InitializeP0] (@ref LIS_NUOPC::InitializeP0)     | Set the Initialize Phase Definition (IPD). Configure model
!! Init   | [InitializeP1] (@ref LIS_NUOPC::InitializeP1)     | Initialize model.  Advertize import and export fields
!! Init   | [InitializeP3] (@ref LIS_NUOPC::InitializeP3)     | Realize import and export fields
!! Init   | [DataInitialize] (@ref LIS_NUOPC::DataInitialize) | Initialize import and export data
!! Init   | [SetClock] (@ref LIS_NUOPC::SetClock)             | Set model clock during initialization
!! Run    | [CheckImport] (@ref LIS_NUOPC::CheckImport)       | Check timestamp on import data.
!! Run    | [ModelAdvance] (@ref LIS_NUOPC::ModelAdvance)     | Advances the model by a timestep
!! Final  | [ModelFinalize] (@ref LIS_NUOPC::ModelFinalize)   | Releases memory
!!
!!
!! @section Initialize Initialize
!!
!! Description of the initialization phases and internal model calls.
!! - [InitializeP0] (@ref LIS_NUOPC::InitializeP0)
!! - [InitializeP1] (@ref LIS_NUOPC::InitializeP1)
!! - [InitializeP3] (@ref LIS_NUOPC::InitializeP3)
!! - [DataInitialize] (@ref LIS_NUOPC::DataInitialize)
!! - [SetClock] (@ref LIS_NUOPC::SetClock)
!!
!! @subsection InitializeP0 InitializeP0
!!
!! During initialize phase 0 the runtime configuration is read in from model
!! attributes and the initialization phase definition version is set to
!! IPDv03.
!!
!! @subsection InitializeP1 InitializeP1
!!
!! During initialize phase 1 the model is initialized and the import and
!! export fields are advertised in nested import and export states. Import
!! fields are configured in the forcing variables list file.
!!
!! @subsection InitializeP3 InitializeP3
!!
!! During initialize phase 3 import and export fields are realized in each
!! nested import and export state if they are connected through NUOPC.
!! Realized fields are created on the LIS grid. All export fields are realized
!! if realize all export fields is turned on.
!!
!! @subsection DataInitialize DataInitialize
!!
!! During data initialize this cap checks the timestamp of all import fields
!! dependent on a coupled model.  Once all dependent import fields have been
!! initialized this cap is marked initalized. The export fields are updated
!! and initialized regardless of the import field dependencies.
!!
!! @subsection SetClock SetClock
!!
!! During set clock the cap creates a new clock for each nest. The time step
!! for each nest is set in LIS configuration file and initialized during LIS
!! initialization. The time accumulation tracker for each timestep is reset to
!! zero.  The cap's time step is updated to the shortest time step
!! of all nests. The restart write time step is also created and the restart
!! write time accumulation tracker is reset to zero.
!!
!!
!! @section Run Run
!!
!! Description of the run phase(s) and internal model calls.
!! - [CheckImport] (@ref LIS_NUOPC::CheckImport)
!! - [ModelAdvance] (@ref LIS_NUOPC::ModelAdvance)
!!
!! @subsection CheckImport CheckImport
!!
!! During check import the import data is checked to verify that it is at
!! the beginning or end of the timestep.
!!
!! @subsection ModelAdvance ModelAdvance
!!
!! During model advance each nest time accumulation tracker is increased by
!! the timestep of the cap.  If the time accumlation tracker is greater than
!! the time step of the nest then the nest is advanced.
!!
!!
!! @section Finalize Finalize
!!
!! Description of the finalize phase and internal model calls.
!! - [ModelFinalize] (@ref LIS_NUOPC::ModelFinalize)
!!
!! @subsection ModelFinalize ModelFinalize
!!
!! During model finalize LIS finalize subroutines are called and memory
!! allocated during cap initialization is released.
!!
!!
!! @subsection ModelConfiguration Model Configuration
!!
!! Model attributes are used to configure the model.
!!
!! Attribute          | Default        | Description
!! -------------------|----------------|-------------------------------------------------------------------------------------
!! Verbosity          | 0              | String, converted into an integer. Bit 16: LIS cap information will be logged.
!! Diagnostic         | 0              | String, converted into an integer. Bit 16: LIS cap diagnostics will be written.
!! realize_all_export | false          | Realize all export fields including non connected fields.
!! config_file        | lis.config     | Set the LIS configuration file.
!! nest_to_nest       | false          | Turn on nest to nest coupling. Each nest will be identified with an integer.
!! import_dependency  | false          | Data initialization will loop until all import field dependencies are satisfied.
!! output_directory   | [CNAME]_OUTPUT | Configure the LIS Cap output directory.
!!
!!
!! @section ModelFields Model Fields
!!
!! The following tables list the import and export fields.
!!
!! @subsection ImportFields Import Fields
!!
!! Import fields arelisted in the import_list parameter.
!!
!! Standard Name  | Units  | Model Variable  | Description                                | Notes
!! ---------------|--------|-----------------|--------------------------------------------|--------------------------------------
!! dummy_field_1  | Pa     | forcing_1       | field description for first import field   | |
!! dummy_field_2  | kg     | forcing_2       | field description for second import field  | |
!! dummy_field_3  | W m-2  | forcing_3       | field description for third import field   | field notes
!!
!! @subsection ExportField Export Fields
!!
!! Export fields are listed in the export_list parameter.
!!
!! Standard Name  | Units   | Model Variable  | Description                               | Notes
!! ---------------|---------|-----------------|-------------------------------------------|---------------------------
!! dummy_field_1  | m       | output_1        | field description for first export field  | field notes
!! dummy_field_2  | kg      | output_2        | field description for second export field | |
!! dummy_field_3  | m s-1   | output_3        | field description for third export field  | field notes
!!
!!
!! @section MemoryManagement Memory Management
!!
!! Model configuration is stored in a custom internal state data type. A
!! pointer to the custom internal state data type is stored in the component.
!!
!! The cap allocates new memory for each field so that 2-D coordinate points
!! can be translated into the LIS tiled field points.
!!
!! @section IO Input and Output
!!
!! Cap diagnostic output is written to the ESMF PET Logs. Cap diagnostic
!! output can be increased or decreased by setting the Verbosity attribute.
!!
!! NUOPC state restart write files are written depending on the
!! RestartInterval attribute. If set to 0 then NUOPC state restart write files
!! will never be written.
!!
!! LIS diagnostics output is written to the LIS logs configured in the LIS
!! configuration file.
!!
!! LIS output files are written to the output directory configured in the LIS
!! configuration file.  LIS output includes LIS history files and LIS restart
!! files.
!!
!! @section Dependencies Dependencies
!!
!! Dependencies
!! - [HDF5 v1.8.11+] (https://support.hdfgroup.org/HDF5/)
!! - [NetCDF v4.3.0+] (http://www.unidata.ucar.edu/software/netcdf/docs/)
!! - [NetCDF FORTRAN] (http://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html)
!! - [JasPer v1.900.1] (https://www.ece.uvic.ca/~frodo/jasper)
!! - [Grib API v1.12.3+] (https://software.ecmwf.int/wiki/display/GRIB/Home)
!!
!! @subsection HDF5 HDF5
!!
!! Configure, Build, and Install
!! $ ./configure --enable-fortran --prefix=HDF5_DIR
!! $ make
!! $ make install
!!
!! @subsection NetCDF NetCDF
!!
!! Add the following environment variables.
!! - CPPFLAGS '-IHDF5_DIR/include'
!! - LDFLAGS '-LHDF5_DIR/lib'
!!
!! Configure, Build, and Install
!! $ ./configure --prefix=NETCDF_DIR --enable-netcdf-4 --disable-dap-remote-tests
!! $ make
!! $ make install
!!
!! @subsection NetCDFFortran NetCDF Fortran
!!
!! Add the following environment variables.
!! - LD_LIBRARY_PATH 'NETCDF_DIR/lib:$LD_LIBRARY_PATH'
!! - CPPFLAGS '-INETCDF_DIR/include'
!! - LDFLAGS '-LNETCDF_DIR/lib'
!!
!! Configure, Build, and Install
!! $ ./configure --prefix=NETCDF_DIR
!! $ make
!! $ make install
!!
!! @subsection JasPer JasPer
!!
!! Configure, Build, and Install
!! $ ./configure --enable-shared --prefix=JASPER_DIR
!! $ make
!! $ make install
!!
!! @subsection GribAPI Grib API
!!
!! Configure, Build, and Install
!! $ ./configure --with-jasper=JASPER_DIR --with-netcdf=NETCDF_DIR --prefix=GRIB_API_DIR
!! $ make
!! $ make install
!!
!! @section BuildingAndInstalling Building and Installing
!!
!! Environment Variables
!! - ESMFMKFILE
!! - JASPER
!! - GRIB_API
!!
!! NUOPC Makefile Targets
!! - nuopc
!! - nuopcinstall
!! - nuopcclean
!!
!! The build system in [Makefile] (@ref Makefile) wraps the LIS build system
!! and adds the nuopc, nuopcinstall, and nuopcclean targets. Before building
!! make sure to configure the internal model.
!!
!! To build and install into the current directory run:
!!    $ make nuopc
!!
!! To install into an alternative directory run:
!!    $ make nuopcinstall DESTDIR=<INSTALL_DIR> INSTDIR=<SUBDIR>
!!
!! To build with debugging information run:
!!    $ make nuopc DEBUG=on
!!
!! @section Repository
!! The LIS NUOPC cap is maintained in a GitHub repository:
!! https://github.com/NESII/lis_cap
!!
!! @section References
!!
!! - [LIS] (https://modelingguru.nasa.gov/community/atmospheric/lis)
!! - [ESPS] (https://www.earthsystemcog.org/projects/esps)
!! - [ESMF] (https://www.earthsystemcog.org/projects/esmf)
!! - [NUOPC] (https://www.earthsystemcog.org/projects/nuopc/)

#define FILENAME "LIS_NUOPC_Cap.F90"
#define MODNAME "LIS_NUOPC_Cap"
#include "LIS_NUOPC_Macros.h"

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
  use LIS_ESMF_Extensions

  implicit none

  private

  public SetServices

  CHARACTER(LEN=*), PARAMETER :: label_InternalState = 'InternalState'
  INTEGER, PARAMETER :: MAXNEST = 999999999

  type type_InternalStateStruct
    character(len=64)     :: configFile       = 'lis.config'
    logical               :: realizeAllExport = .FALSE.
    logical               :: nestToNest       = .FALSE.
    logical               :: importDependency = .FALSE.
    character(len=40)     :: dirOutput        = "."
    integer               :: nnests           = 0
    integer               :: nfields          = size(LIS_FieldList)
    type(ESMF_Grid),allocatable         :: grids(:)
    type(ESMF_Clock),allocatable        :: clocks(:)
    type(ESMF_TimeInterval),allocatable :: stepAccum(:)
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
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out
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
!    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call LIS_AttributeGet(rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! prepare diagnostics folder
    if (btest(diagnostic,16)) then
      call ESMF_UtilIOMkDir(pathName=trim(is%wrap%dirOutput), &
        relaxedFlag=.true., rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif

    contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine LIS_AttributeGet(rc)
      integer, intent(out)  :: rc

      ! local variables
      logical                    :: configIsPresent
      type(ESMF_Config)          :: config
      type(NUOPC_FreeFormat)     :: attrFF
      character(ESMF_MAXSTR)     :: logMsg

      ! check gcomp for config
      call ESMF_GridCompGet(gcomp, configIsPresent=configIsPresent, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      ! read and ingest free format component attributes
      if (configIsPresent) then
        call ESMF_GridCompGet(gcomp, config=config, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        attrFF = NUOPC_FreeFormatCreate(config, &
          label=trim(cname)//"_attributes::", relaxedflag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        call NUOPC_CompAttributeIngest(gcomp, attrFF, addFlag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      ! Realize all export fields
      call ESMF_AttributeGet(gcomp, name="realize_all_export", value=value, &
        defaultValue="false", convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      is%wrap%realizeAllExport = (trim(value)=="true")

      ! Set configuration file name
      call ESMF_AttributeGet(gcomp, name="config_file", value=value, &
        defaultValue="lis.config", convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      is%wrap%configFile=value

      ! Turn on nest coupling
      call ESMF_AttributeGet(gcomp, name="nest_to_nest", value=value, &
        defaultValue="false", convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      is%wrap%nestToNest = (trim(value)=="true")

      ! Realize all export fields
      call ESMF_AttributeGet(gcomp, name="import_dependency", &
        value=value, defaultValue="false", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      is%wrap%importDependency = (trim(value)=="true")

      ! Get component output directory
      call ESMF_AttributeGet(gcomp, name="output_directory", &
        value=value, defaultValue=trim(cname)//"_OUTPUT", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
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
          'Import Dependency      = ',is%wrap%importDependency
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
!    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call ESMF_GridCompGet(gcomp, vm=vm, localPet=localPet, petCount=petCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! initialize lis model for this PET
    call LIS_NUOPC_Init(vm, configFile=is%wrap%configFile, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (btest(verbosity,16)) then
      call LIS_Log(trim(cname)//': '//rname,rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif

    call LIS_FieldDictionaryAdd(rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    is%wrap%nnests = LIS_NestCntGet(rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! Max nest check
    if ( is%wrap%nnests > MAXNEST ) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg="Maximum nest size is 999,999,999.", &
        line=__LINE__,file=__FILE__,rcToReturn=rc)
      return  ! bail out
    endif

    allocate( &
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
      return  ! bail out

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
        return  ! bail out
      endif
    else
      ! add namespace
      call NUOPC_AddNestedState(importState, &
        CplSet="1", &
        nestedStateName="NestedStateImp_N1", &
        nestedState=is%wrap%NStateImp(1), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      call NUOPC_AddNestedState(exportState, &
        CplSet="1", &
        nestedStateName="NestedStateExp_N1", &
        nestedState=is%wrap%NStateExp(1), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif

    do nIndex = 2, is%wrap%nnests
      write (nStr,"(I0)") nIndex
      call NUOPC_AddNestedState(importState, &
        CplSet=trim(nStr), &
        nestedStateName="NestedStateImp_N"//trim(nStr), &
        nestedState=is%wrap%NStateImp(nIndex), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      call NUOPC_AddNestedState(exportState, &
        CplSet=trim(nStr), &
        nestedStateName="NestedStateExp_N"//trim(nStr), &
        nestedState=is%wrap%NStateExp(nIndex), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
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

    if (btest(verbosity,16)) call LogAdvertised()

    contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
!    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    do nIndex = 1, is%wrap%nnests
      write (nStr,"(I0)") nIndex

      ! Call gluecode to create grid.
      is%wrap%grids(nIndex) = LIS_GridCreate(nIndex, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      if (btest(verbosity,16)) then
        call LIS_ESMF_LogGrid(is%wrap%grids(nIndex), &
          trim(cname)//"_"//rname//"_D"//trim(nStr),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      ! Write grid to NetCDF file.
      if (btest(diagnostic,16)) then
        call LIS_ESMF_GridWrite(is%wrap%grids(nIndex), &
          trim(is%wrap%dirOutput)//"/diag_"//trim(cname)//"_"// &
          rname//'_grid_D'//trim(nStr)//".nc", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
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
          field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
            grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          call NUOPC_Realize(is%wrap%NStateExp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedExport = .TRUE.
          call ESMF_FieldGet(field=field,array=array,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
            grid=is%wrap%grids(nIndex), array=array, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call NUOPC_Realize(is%wrap%NStateImp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedImport = .TRUE.
        elseif (realizeImp .AND. realizeExp) then
          field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
            grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          call NUOPC_Realize(is%wrap%NStateExp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedExport = .TRUE.
          field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
            grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, rc=rc)
          if(ESMF_STDERRORCHECK(rc)) return
          call NUOPC_Realize(is%wrap%NStateImp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedImport = .TRUE.
        elseif (realizeExp) then
          field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
            grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          call NUOPC_Realize(is%wrap%NStateExp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedExport = .TRUE.
          call ESMF_StateRemove(is%wrap%NStateImp(nIndex),(/trim(LIS_FieldList(fIndex)%stateName)/), &
            relaxedflag=.true.,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        elseif (realizeImp) then
          call ESMF_StateRemove(is%wrap%NStateExp(nIndex),(/trim(LIS_FieldList(fIndex)%stateName)/), &
            relaxedflag=.true.,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          field = ESMF_FieldCreate(name=trim(LIS_FieldList(fIndex)%stateName), &
            grid=is%wrap%grids(nIndex), typekind=ESMF_TYPEKIND_FIELD, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          call NUOPC_Realize(is%wrap%NStateImp(nIndex), field=field,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          LIS_FieldList(fIndex)%realizedImport = .TRUE.
        else
          call ESMF_StateRemove(is%wrap%NStateExp(nIndex),(/trim(LIS_FieldList(fIndex)%stateName)/), &
            relaxedflag=.true.,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call ESMF_StateRemove(is%wrap%NStateImp(nIndex),(/trim(LIS_FieldList(fIndex)%stateName)/), &
            relaxedflag=.true.,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        endif
      enddo

      call LIS_ESMF_FillState(is%wrap%NStateImp(nIndex),value=MISSINGVALUE,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call LIS_ESMF_FillState(is%wrap%NStateExp(nIndex),value=MISSINGVALUE,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      is%wrap%modes(nIndex) = LIS_RunModeGet(LIS_FieldList,is%wrap%NStateImp(nIndex),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    enddo

    if (btest(verbosity,16)) call LogRealized()

    contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
    type(ESMF_Field)                       :: field
    logical                                :: importCurrent
    logical                                :: importUpdated

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
!    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query the Component for its clock
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! get the current time out of the clock
    call ESMF_ClockGet(modelClock, currTime=currTime, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_TimeGet(currTime, timeString=currTimeStr, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    importUpdated = .TRUE.
    do nIndex=1,is%wrap%nnests
      write (nStr,"(I0)") nIndex

      if (is%wrap%importDependency) then
        importCurrent = NUOPC_IsAtTime(is%wrap%NStateImp(nIndex), &
          time=currTime, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out

        if (importCurrent) then
          call ESMF_LogWrite( &
            trim(cname)//': '//rname//' Initialize-Data-Dependency SATISFIED!!! Nest='//trim(nStr), &
            ESMF_LOGMSG_INFO)
          call LIS_ImportFieldsCopy(nIndex,is%wrap%NStateImp(nIndex),rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        else
          call ESMF_LogWrite( &
            trim(cname)//': '//rname//' Initialize-Data-Dependency NOT YET SATISFIED!!! Nest='//trim(nStr), &
            ESMF_LOGMSG_INFO)
          importUpdated = .FALSE.
        endif
      endif

      ! Initialize import and export fields
      call LIS_NUOPC_DataInit(nest=nIndex, &
        importState=is%wrap%NStateImp(nIndex), &
        exportState=is%wrap%NStateExp(nIndex),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

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
    enddo ! enddo nnests

    ! set InitializeDataComplete Attribute to "true", indicating to the
    ! generic code that all inter-model data dependencies are satisfied
    if (importUpdated) then
      call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      ! Write initialization files
      if (btest(diagnostic,16)) then
        do nIndex=1,is%wrap%nnests
          write (nStr,"(I0)") nIndex
          call NUOPC_Write(is%wrap%NStateImp(nIndex), &
            fileNamePrefix=trim(is%wrap%dirOutput)//"/diag_"//trim(cname)//"_"// &
              rname//"_imp_D"//trim(nStr)//"_"//trim(currTimeStr)//"_", &
            overwrite=.true., status=ESMF_FILESTATUS_REPLACE, timeslice=1, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
          call NUOPC_Write(is%wrap%NStateExp(nIndex), &
            fileNamePrefix=trim(is%wrap%dirOutput)//"/diag_"//trim(cname)//"_"// &
              rname//"_exp_D"//trim(nStr)//"_"//trim(currTimeStr)//"_", &
            overwrite=.true., status=ESMF_FILESTATUS_REPLACE, timeslice=1, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        enddo
      endif
    endif

  end subroutine

  !-----------------------------------------------------------------------------

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
!    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

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

      call ESMF_TimeIntervalSet(is%wrap%stepAccum(nIndex), &
        s_r8=0._ESMF_KIND_R8, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    enddo

    call ESMF_TimeIntervalSet(modelTimestep, &
      s_r8=mindt, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call NUOPC_CompSetClock(gcomp, modelClock, &
      modelTimestep, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

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
!    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query the component for its clock, importState, and exportState
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! get the stop time out of the clock
    call ESMF_ClockGet(modelClock, startTime=modelStartTime, &
      currTime=modelCurrTime, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    do nIndex=1,is%wrap%nnests
      write (nStr,"(I0)") nIndex
      allCurrTime = NUOPC_IsAtTime(is%wrap%NStateImp(nIndex), modelCurrTime,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      if (.NOT.allCurrTime) then
        call ESMF_LogWrite(trim(cname)//": NUOPC INCOMPATIBILITY DETECTED: "// &
          "Import Fields Nest="//trim(nStr)//" not at correct time", &
          ESMF_LOGMSG_WARNING)
      endif
    enddo

  end subroutine

  !-----------------------------------------------------------------------------

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
!    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query the component for its clock, importState, and exportState
    call NUOPC_ModelGet(gcomp, &
      modelClock=modelClock, &
      importState=importState, &
      exportState=exportState, &
      rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query the clock for the current time and time step
    call ESMF_ClockGet(modelClock, &
      currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    advEndTime = currTime + timeStep
    call ESMF_TimeGet(currTime, timeString=currTimeStr, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_TimeGet(advEndTime, timeString=advEndTimeStr, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

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
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif

      is%wrap%stepAccum(nIndex) = is%wrap%stepAccum(nIndex) + timeStep

      call ESMF_ClockGet(is%wrap%clocks(nIndex),timeStep=nestTimestep,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      do while (is%wrap%stepAccum(nIndex) >= nestTimestep)
        ! Gluecode NestAdvance
        if (btest(verbosity,16)) then
          call LogAdvance(nIndex)
        endif
        call LIS_NUOPC_Run(nIndex,is%wrap%modes(nIndex), &
          is%wrap%NStateImp(nIndex),is%wrap%NStateExp(nIndex), &
          is%wrap%clocks(nIndex), rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        call ESMF_ClockAdvance(is%wrap%clocks(nIndex),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
        is%wrap%stepAccum(nIndex) = &
          is%wrap%stepAccum(nIndex) - nestTimestep
      enddo

      ! Write export files
      if (btest(diagnostic,16)) then
        call NUOPC_Write(is%wrap%NStateExp(nIndex), &
          fileNamePrefix=trim(is%wrap%dirOutput)//"/diag_"//trim(cname)//"_"// &
            rname//"_exp_D"//trim(nStr)//"_"//trim(advEndTimeStr)//"_", &
          overwrite=.true., status=ESMF_FILESTATUS_REPLACE, timeslice=1, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return  ! bail out
      endif
    enddo

    contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
!    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(gcomp, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query Component for its internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! query the Component for its clock
    call NUOPC_ModelGet(gcomp, modelClock=modelClock, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call ESMF_ClockGet(modelClock, currTime=currTime, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    call ESMF_TimeGet(currTime, timeString=currTimeStr, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    ! finalize the LIS model
    do nIndex=1,is%wrap%nnests
      call LIS_NUOPC_Final(nIndex,modelClock,rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    enddo

    deallocate(is%wrap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg='Deallocation of internal state memory failed.', &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out

  end subroutine

end module

