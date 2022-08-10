module lis_nuopc_test_atm

  !-----------------------------------------------------------------------------
  ! ATM Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS        => SetServices, &
    model_label_CheckImport => label_CheckImport, &
    model_label_Advance     => label_Advance

  implicit none

  private

  public SetServices

  integer, parameter :: outputPet = 0
  integer, parameter :: nx = 64
  integer, parameter :: ny = 32
  integer, parameter :: nz = 4
  real(ESMF_KIND_R4), parameter :: filv = -1.0E34

  type fld3d
    character(len=64)           :: std        = "dummy"
    type(ESMF_Field), pointer   :: fld        => null()
    real(ESMF_KIND_R4), pointer :: ptr(:,:,:) => null()
    real(ESMF_KIND_R4)          :: dft        = filv
    logical                     :: rlz        = .false.
    real(ESMF_KIND_R4)          :: lsum(1)    = filv
    real(ESMF_KIND_R4)          :: gsum(1)    = filv
  endtype fld3d

  type fld2d
    character(len=64)           :: std      = "dummy"
    type(ESMF_Field), pointer   :: fld      => null()
    real(ESMF_KIND_R4), pointer :: ptr(:,:) => null()
    real(ESMF_KIND_R4)          :: dft      = filv
    logical                     :: rlz      = .false.
    real(ESMF_KIND_R4)          :: lsum(1)  = filv
    real(ESMF_KIND_R4)          :: gsum(1)  = filv
  endtype fld2d

  ! export fields
  type(fld2d) :: exp_tair = &
    fld2d(std="air_temperature", dft=289.5)
  type(fld2d) :: exp_qair = &
    fld2d(std="specific_humidity", dft=0.01)
  type(fld2d) :: exp_swdown = &
    fld2d(std="surface_downwelling_shortwave_flux_in_air", dft=87.5)
  type(fld2d) :: exp_lwdown = &
    fld2d(std="surface_downwelling_longwave_flux_in_air", dft=337)
  type(fld2d) :: exp_winde = &
    fld2d(std="eastward_wind", dft=-0.5)
  type(fld2d) :: exp_windn = &
    fld2d(std="northward_wind", dft=4.0)
  type(fld2d) :: exp_psurf = &
    fld2d(std="surface_air_pressure", dft=97600)
  type(fld2d) :: exp_rainf = &
    fld2d(std="rainfall_flux", dft=0.000)
  type(fld2d) :: exp_crainf = &
    fld2d(std="convective_rainfall_flux", dft=0.0000)

  interface RealizeField
    module procedure RealizeField3D
    module procedure RealizeField2D
  end interface

  interface SumField
    module procedure SumField3D
    module procedure SumField2D
  end interface

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! derive generic model phases
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! set entry points for initialization phases
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! set entry points for specialized methods
    call ESMF_MethodRemove(model, label=model_label_CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompSpecialize(model, specLabel=model_label_CheckImport, &
       specRoutine=CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    ! local variables
    character(len=64)    :: value

    rc = ESMF_SUCCESS

    ! advertise export fields
    call NUOPC_Advertise(exportState, exp_tair%std, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, exp_qair%std, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, exp_swdown%std, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, exp_lwdown%std, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, exp_winde%std, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, exp_windn%std, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, exp_psurf%std, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, exp_rainf%std, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, exp_crainf%std, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine InitializeP1

  !-----------------------------------------------------------------------------

  subroutine InitializeP2(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    ! local variables
    character(len=64)       :: value
    integer                 :: diagnostic
    character(len=64)       :: coord_type
    type(ESMF_Grid)         :: atm_grid

    rc = ESMF_SUCCESS

    ! diagnostic
    call ESMF_AttributeGet(model, name="Diagnostic", &
      value=value, defaultValue="0", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"max ","high","low ","off "/), &
      specialValueList= (/ 65535, 65535, 65535,     0/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! field geom type
    call ESMF_AttributeGet(model, name="coord_type", &
      value=coord_type, defaultValue="GRD_COORD_SPHDEG_LW", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! create grid
    if (coord_type .eq. "GRD_COORD_CARTESIAN_LW") then
      atm_grid = ESMF_GridCreateNoPeriDimUfrm(name="ATM-Grid", &
        minIndex=(/1, 1/), maxIndex=(/nx, ny/), &
        minCornerCoord=(/    0._ESMF_KIND_R8,     0._ESMF_KIND_R8/), &
        maxCornerCoord=(/64000._ESMF_KIND_R8, 32000._ESMF_KIND_R8/), &
        staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
        coordSys=ESMF_COORDSYS_CART, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    elseif (coord_type .eq. "GRD_COORD_SPHDEG_LW") then
      atm_grid = ESMF_GridCreateNoPeriDimUfrm(name="ATM-Grid", &
        minIndex=(/1, 1/), maxIndex=(/nx, ny/), &
        minCornerCoord=(/-98.5_ESMF_KIND_R8, 34.5_ESMF_KIND_R8/), &
        maxCornerCoord=(/-97.5_ESMF_KIND_R8, 35.5_ESMF_KIND_R8/), &
        staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
        coordSys=ESMF_COORDSYS_SPH_DEG, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="Unsupported coordinate type", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    ! write grid to NetCDF file
    if (btest(diagnostic,16)) then
      call Grid_Diag(atm_grid, "diagnostic_ATM_InitializeP2_grid.nc", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

    ! create export fields
    call RealizeField(exp_tair, atm_grid, exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call RealizeField(exp_qair, atm_grid, exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call RealizeField(exp_swdown, atm_grid, exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call RealizeField(exp_lwdown, atm_grid, exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call RealizeField(exp_winde, atm_grid, exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call RealizeField(exp_windn, atm_grid, exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call RealizeField(exp_psurf, atm_grid, exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call RealizeField(exp_rainf, atm_grid, exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call RealizeField(exp_crainf, atm_grid, exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine InitializeP2

  !-----------------------------------------------------------------------------

  subroutine CheckImport(model, rc)
    type(ESMF_GridComp) :: model
    integer,intent(out) :: rc
    ! local variables
    character(32)               :: cname
    type(ESMF_Clock)            :: modelClock
    type(ESMF_Time)             :: modelCurrTime
    type(ESMF_State)            :: importState
    logical                     :: allCurrTime

    rc = ESMF_SUCCESS

    ! query component for name
    call ESMF_GridCompGet(model, name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! query the component for its clock and import state
    call NUOPC_ModelGet(model, modelClock=modelClock, &
      importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! get the stop time out of the clock
    call ESMF_ClockGet(modelClock, currTime=modelCurrTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    allCurrTime = NUOPC_IsAtTime(importState, modelCurrTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (.NOT.allCurrTime) then
      call ESMF_LogWrite(trim(cname)//": "// &
        "NUOPC INCOMPATIBILITY DETECTED: Import Fields not at current time", &
        ESMF_LOGMSG_WARNING)
    endif
  end subroutine CheckImport

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    ! local variables
    type(ESMF_VM)      :: vm
    integer            :: localPet
    type(ESMF_Clock)   :: modelClock
    type(ESMF_State)   :: importState
    type(ESMF_State)   :: exportState
    character(len=160) :: clockString

    rc = ESMF_SUCCESS

    ! query component for vm and local pet
    call ESMF_GridCompGet(model, vm=vm, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! query component for import and export states
    call NUOPC_ModelGet(model, modelClock=modelClock, &
      importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_ClockPrint(modelClock, options="currTime", &
      unit=clockString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! sum export data from all PETs
    call SumField(exp_tair, vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call SumField(exp_qair, vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call SumField(exp_swdown, vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call SumField(exp_lwdown, vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call SumField(exp_winde, vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call SumField(exp_windn, vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call SumField(exp_psurf, vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call SumField(exp_rainf, vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call SumField(exp_crainf, vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! print import field sums
    if (localPet .eq. outputPet) then
      print *,"ATM Export Sums"
      print *,"  clock=",trim(clockString)
      print *,"  sum(exp_tair)=",exp_tair%gsum(1)
      print *,"  sum(exp_qair)=",exp_qair%gsum(1)
      print *,"  sum(exp_swdown)=",exp_swdown%gsum(1)
      print *,"  sum(exp_lwdown)=",exp_lwdown%gsum(1)
      print *,"  sum(exp_winde)=",exp_winde%gsum(1)
      print *,"  sum(exp_windn)=",exp_windn%gsum(1)
      print *,"  sum(exp_psurf)=",exp_psurf%gsum(1)
      print *,"  sum(exp_rainf)=",exp_rainf%gsum(1)
      print *,"  sum(exp_crainf)=",exp_crainf%gsum(1)
    end if

  end subroutine ModelAdvance

  subroutine Grid_Diag(grid, fileName, overwrite, status, timeslice, iofmt, &
  relaxedflag, rc)
    type(ESMF_Grid), intent(in)                      :: grid
    character(len=*), intent(in), optional           :: fileName
    logical, intent(in), optional                    :: overwrite
    type(ESMF_FileStatus_Flag), intent(in), optional :: status
    integer, intent(in), optional                    :: timeslice
    type(ESMF_IOFmt_Flag), intent(in), optional      :: iofmt
    logical, intent(in), optional                    :: relaxedflag
    integer, intent(out)                             :: rc
    ! local variables

    logical                 :: ioCapable
    logical                 :: doItFlag
    character(len=64)       :: lfileName
    character(len=64)       :: gridName
    type(ESMF_Array)        :: array
    type(ESMF_ArrayBundle)  :: arraybundle
    logical                 :: isPresent
    integer                 :: dimCount
    integer                 :: dimIndex
    integer,allocatable     :: coordDimCount(:)
    integer                 :: coordDimMax
    integer                 :: stat
    logical                 :: hasCorners

    rc = ESMF_SUCCESS

    ioCapable = (ESMF_IO_PIO_PRESENT .and. &
      (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT))

    doItFlag = .true. ! default
    if (present(relaxedFlag)) then
      doItFlag = .not.relaxedflag .or. (relaxedflag.and.ioCapable)
    endif

    if (doItFlag) then

      if (present(fileName)) then
        lfileName = trim(fileName)
      else
        call ESMF_GridGet(grid, name=gridName, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        lfileName = trim(gridName)//".nc"
      endif

      arraybundle = ESMF_ArrayBundleCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      ! -- centers --

      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="lon_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="lat_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      ! -- corners --

      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        isPresent=hasCorners, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (hasCorners) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rcToCheck=rc)) then
          call ESMF_ArraySet(array, name="lon_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rcToCheck=rc)) then
          call ESMF_ArraySet(array, name="lat_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
      endif

      ! -- mask --

      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="mask", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      ! -- area --

      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="area", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      call ESMF_ArrayBundleWrite(arraybundle, &
        fileName=trim(lfileName),rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call ESMF_ArrayBundleDestroy(arraybundle,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif
  end subroutine

  subroutine RealizeField3D(field, grid, layers, state, rc)
    type(fld3d), intent(inout)  :: field
    type(ESMF_Grid), intent(in) :: grid
    integer, intent(in)         :: layers
    type(ESMF_State)            :: state
    integer, intent(out)        :: rc
    ! local variables

    rc = ESMF_SUCCESS

    if (associated(field%fld)) then
      call ESMF_LogSetError(ESMF_RC_MEM_ALLOCATE, &
        msg="Field already associated: "//trim(field%std), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
    allocate(field%fld)
    field%fld = ESMF_FieldCreate(name=trim(field%std), grid=grid, &
      typekind=ESMF_TYPEKIND_R4, gridToFieldMap=(/1,3/), &
      ungriddedLBound=(/1/), ungriddedUBound=(/layers/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Realize(state, field=field%fld, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldFill(field%fld, dataFillScheme="const", &
      const1=real(field%dft,ESMF_KIND_R8), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(field%fld, farrayPtr=field%ptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    field%rlz = .true.
  end subroutine

  subroutine RealizeField2D(field, grid, state, rc)
    type(fld2d), intent(inout)  :: field
    type(ESMF_Grid), intent(in) :: grid
    type(ESMF_State)            :: state
    integer, intent(out)        :: rc
    ! local variables

    rc = ESMF_SUCCESS

    if (associated(field%fld)) then
      call ESMF_LogSetError(ESMF_RC_MEM_ALLOCATE, &
        msg="Field already associated: "//trim(field%std), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
    allocate(field%fld)
    field%fld = ESMF_FieldCreate(name=trim(field%std), grid=grid, &
      typekind=ESMF_TYPEKIND_R4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Realize(state, field=field%fld, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldFill(field%fld, dataFillScheme="const", &
      const1=real(field%dft,ESMF_KIND_R8), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(field%fld, farrayPtr=field%ptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    field%rlz = .true.
  end subroutine

  subroutine SumField3D(field, vm, rc)
    type(fld3d), intent(inout) :: field
    type(ESMF_VM), intent(in)  :: vm
    integer, intent(out)       :: rc
    ! local variables

    rc = ESMF_SUCCESS

    if (field%rlz) then
      field%lsum(1)=sum(field%ptr,field%ptr.ne.filv)
      call ESMF_VMReduce(vm=vm, sendData=field%lsum, &
        recvData=field%gsum, count=1, &
        reduceflag=ESMF_REDUCE_SUM, rootPet=outputPet, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      field%lsum(1)=filv
      field%gsum(1)=filv
    endif
  end subroutine

  subroutine SumField2D(field, vm, rc)
    type(fld2d), intent(inout) :: field
    type(ESMF_VM), intent(in)  :: vm
    integer, intent(out)       :: rc
    ! local variables

    rc = ESMF_SUCCESS

    if (field%rlz) then
      field%lsum(1)=sum(field%ptr,field%ptr.ne.filv)
      call ESMF_VMReduce(vm=vm, sendData=field%lsum, &
        recvData=field%gsum, count=1, &
        reduceflag=ESMF_REDUCE_SUM, rootPet=outputPet, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      field%lsum(1)=filv
      field%gsum(1)=filv
    endif
  end subroutine

end module lis_nuopc_test_atm
