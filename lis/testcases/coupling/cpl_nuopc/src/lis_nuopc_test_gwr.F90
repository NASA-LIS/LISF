module lis_nuopc_test_gwr

  !-----------------------------------------------------------------------------
  ! GWR Component.
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
  integer, parameter :: em = 20

  real(ESMF_KIND_R4), parameter :: filv = -1.0E34
  logical :: ensemble = .false.

  type fld4d
    character(len=64)           :: std          = "dummy"
    type(ESMF_Field), pointer   :: fld          => null()
    real(ESMF_KIND_R4), pointer :: ptr(:,:,:,:) => null()
    real(ESMF_KIND_R4)          :: dft          = filv
    logical                     :: rlz          = .false.
    real(ESMF_KIND_R4)          :: lsum(1)      = filv
    real(ESMF_KIND_R4)          :: gsum(1)      = filv
  endtype fld4d

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

  ! import fields
  type(fld2d) :: imp_flux1 = &
    fld2d(std="total_water_flux_layer_1", dft=filv)
  type(fld2d) :: imp_flux2 = &
    fld2d(std="total_water_flux_layer_2", dft=filv)
  type(fld2d) :: imp_flux3 = &
    fld2d(std="total_water_flux_layer_3", dft=filv)
  type(fld2d) :: imp_flux4 = &
    fld2d(std="total_water_flux_layer_4", dft=filv)
  ! import fields ensemble
  type(fld3d) :: imp_flux1ens = &
    fld3d(std="total_water_flux_layer_1", dft=filv)
  type(fld3d) :: imp_flux2ens = &
    fld3d(std="total_water_flux_layer_2", dft=filv)
  type(fld3d) :: imp_flux3ens = &
    fld3d(std="total_water_flux_layer_3", dft=filv)
  type(fld3d) :: imp_flux4ens = &
    fld3d(std="total_water_flux_layer_4", dft=filv)
  ! export fields
  type(fld2d) :: exp_smois1 = &
    fld2d(std="soil_moisture_fraction_layer_1          ", dft=0.20)
  type(fld2d) :: exp_smois2 = &
    fld2d(std="soil_moisture_fraction_layer_2          ", dft=0.20)
  type(fld2d) :: exp_smois3 = &
    fld2d(std="soil_moisture_fraction_layer_3          ", dft=0.20)
  type(fld2d) :: exp_smois4 = &
    fld2d(std="soil_moisture_fraction_layer_4          ", dft=0.20)
  type(fld2d) :: exp_sh2o1 = &
    fld2d(std="liquid_fraction_of_soil_moisture_layer_1", dft=0.20)
  type(fld2d) :: exp_sh2o2 = &
    fld2d(std="liquid_fraction_of_soil_moisture_layer_2", dft=0.20)
  type(fld2d) :: exp_sh2o3 = &
    fld2d(std="liquid_fraction_of_soil_moisture_layer_3", dft=0.20)
  type(fld2d) :: exp_sh2o4 = &
    fld2d(std="liquid_fraction_of_soil_moisture_layer_4", dft=0.20)
  type(fld2d) :: exp_gws = &
    fld2d(std="ground_water_storage                    ", dft=4900)
  ! export fields ensemble
  type(fld3d) :: exp_smois1ens = &
    fld3d(std="soil_moisture_fraction_layer_1          ", dft=0.20)
  type(fld3d) :: exp_smois2ens = &
    fld3d(std="soil_moisture_fraction_layer_2          ", dft=0.20)
  type(fld3d) :: exp_smois3ens = &
    fld3d(std="soil_moisture_fraction_layer_3          ", dft=0.20)
  type(fld3d) :: exp_smois4ens = &
    fld3d(std="soil_moisture_fraction_layer_4          ", dft=0.20)
  type(fld3d) :: exp_sh2o1ens = &
    fld3d(std="liquid_fraction_of_soil_moisture_layer_1", dft=0.20)
  type(fld3d) :: exp_sh2o2ens = &
    fld3d(std="liquid_fraction_of_soil_moisture_layer_2", dft=0.20)
  type(fld3d) :: exp_sh2o3ens = &
    fld3d(std="liquid_fraction_of_soil_moisture_layer_3", dft=0.20)
  type(fld3d) :: exp_sh2o4ens = &
    fld3d(std="liquid_fraction_of_soil_moisture_layer_4", dft=0.20)
  type(fld3d) :: exp_gwsens = &
    fld3d(std="ground_water_storage                    ", dft=4900)

  interface RealizeField
    module procedure RealizeField4D
    module procedure RealizeField3D
    module procedure RealizeField2D
  end interface

  interface SumField
    module procedure SumField4D
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

    ! ensemble
    call ESMF_AttributeGet(model, name="ensemble", value=value, &
      defaultValue="false", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    value = ESMF_UtilStringLowerCase(value, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ensemble = (trim(value)=="true")

    if (ensemble) then
      ! advertise import fields
      call NUOPC_Advertise(importState, imp_flux1ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(importState, imp_flux2ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(importState, imp_flux3ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(importState, imp_flux4ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ! advertise export fields
      call NUOPC_Advertise(exportState, exp_smois1ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_smois2ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_smois3ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_smois4ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_sh2o1ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_sh2o2ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_sh2o3ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_sh2o4ens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_gwsens%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      ! advertise import fields
      call NUOPC_Advertise(importState, imp_flux1%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(importState, imp_flux2%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(importState, imp_flux3%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(importState, imp_flux4%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ! advertise export fields
      call NUOPC_Advertise(exportState, exp_smois1%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_smois2%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_smois3%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_smois4%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_sh2o1%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_sh2o2%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_sh2o3%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_sh2o4%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Advertise(exportState, exp_gws%std, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

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
    type(ESMF_Grid)         :: gwr_grid

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
      gwr_grid = ESMF_GridCreateNoPeriDimUfrm(name="GWR-Grid", &
        minIndex=(/1, 1/), maxIndex=(/nx, ny/), &
        minCornerCoord=(/    0._ESMF_KIND_R8,     0._ESMF_KIND_R8/), &
        maxCornerCoord=(/64000._ESMF_KIND_R8, 32000._ESMF_KIND_R8/), &
        staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
        coordSys=ESMF_COORDSYS_CART, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    elseif (coord_type .eq. "GRD_COORD_SPHDEG_LW") then
      gwr_grid = ESMF_GridCreateNoPeriDimUfrm(name="GWR-Grid", &
        minIndex=(/1, 1/), maxIndex=(/nx, ny/), &
        minCornerCoord=(/-98.426653_ESMF_KIND_R8, 34.739932_ESMF_KIND_R8/), &
        maxCornerCoord=(/-97.718663_ESMF_KIND_R8, 35.031552_ESMF_KIND_R8/), &
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
      call Grid_Diag(gwr_grid, "diagnostic_GWR_InitializeP2_grid.nc", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif

    if (ensemble) then
      ! create import fields
      call RealizeField(imp_flux1ens, gwr_grid, importState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(imp_flux2ens, gwr_grid, importState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(imp_flux3ens, gwr_grid, importState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(imp_flux4ens, gwr_grid, importState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ! create export fields
      call RealizeField(exp_smois1ens, gwr_grid, exportState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_smois2ens, gwr_grid, exportState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_smois3ens, gwr_grid, exportState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_smois4ens, gwr_grid, exportState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_sh2o1ens, gwr_grid, exportState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_sh2o2ens, gwr_grid, exportState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_sh2o3ens, gwr_grid, exportState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_sh2o4ens, gwr_grid, exportState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_gwsens, gwr_grid, exportState, members=em, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
       ! create import fields
      call RealizeField(imp_flux1, gwr_grid, importState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(imp_flux2, gwr_grid, importState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(imp_flux3, gwr_grid, importState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(imp_flux4, gwr_grid, importState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ! create export fields
      call RealizeField(exp_smois1, gwr_grid, exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_smois2, gwr_grid, exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_smois3, gwr_grid, exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_smois4, gwr_grid, exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_sh2o1, gwr_grid, exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_sh2o2, gwr_grid, exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_sh2o3, gwr_grid, exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_sh2o4, gwr_grid, exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call RealizeField(exp_gws, gwr_grid, exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
   endif

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

    if (ensemble) then
      ! sum import data from all PETs
      call SumField(imp_flux1ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(imp_flux2ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(imp_flux3ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(imp_flux4ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ! sum export data from all PETs
      call SumField(exp_smois1ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_smois2ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_smois3ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_smois4ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_sh2o1ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_sh2o2ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_sh2o3ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_sh2o4ens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_gwsens, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ! print field sums
      if (localPet .eq. outputPet) then
        print *,"GWR Import Sums (Ensembles)"
        print *,"  clock=",trim(clockString)
        print *,"  sum(imp_flux1ens)=",imp_flux1ens%gsum(1)
        print *,"  sum(imp_flux2ens)=",imp_flux2ens%gsum(1)
        print *,"  sum(imp_flux3ens)=",imp_flux3ens%gsum(1)
        print *,"  sum(imp_flux4ens)=",imp_flux4ens%gsum(1)
        print *,"  sum(exp_smois1ens)=",exp_smois1ens%gsum(1)
        print *,"  sum(exp_smois2ens)=",exp_smois2ens%gsum(1)
        print *,"  sum(exp_smois3ens)=",exp_smois3ens%gsum(1)
        print *,"  sum(exp_smois4ens)=",exp_smois4ens%gsum(1)
        print *,"  sum(exp_sh2o1ens)=",exp_sh2o1ens%gsum(1)
        print *,"  sum(exp_sh2o2ens)=",exp_sh2o2ens%gsum(1)
        print *,"  sum(exp_sh2o3ens)=",exp_sh2o3ens%gsum(1)
        print *,"  sum(exp_sh2o4ens)=",exp_sh2o4ens%gsum(1)
        print *,"  sum(exp_gwsens)=",exp_gwsens%gsum(1)
      endif
    else
      ! sum import data from all PETs
      call SumField(imp_flux1, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(imp_flux2, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(imp_flux3, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(imp_flux4, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ! sum export data from all PETs
      call SumField(exp_smois1, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_smois2, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_smois3, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_smois4, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_sh2o1, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_sh2o2, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_sh2o3, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_sh2o4, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call SumField(exp_gws, vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      ! print field sums
      if (localPet .eq. outputPet) then
        print *,"GWR Import Sums"
        print *,"  clock=",trim(clockString)
        print *,"  sum(imp_flux1)=",imp_flux1%gsum(1)
        print *,"  sum(imp_flux2)=",imp_flux2%gsum(1)
        print *,"  sum(imp_flux3)=",imp_flux3%gsum(1)
        print *,"  sum(imp_flux4)=",imp_flux4%gsum(1)
        print *,"  sum(exp_smois1)=",exp_smois1%gsum(1)
        print *,"  sum(exp_smois2)=",exp_smois2%gsum(1)
        print *,"  sum(exp_smois3)=",exp_smois3%gsum(1)
        print *,"  sum(exp_smois4)=",exp_smois4%gsum(1)
        print *,"  sum(exp_sh2o1)=",exp_sh2o1%gsum(1)
        print *,"  sum(exp_sh2o2)=",exp_sh2o2%gsum(1)
        print *,"  sum(exp_sh2o3)=",exp_sh2o3%gsum(1)
        print *,"  sum(exp_sh2o4)=",exp_sh2o4%gsum(1)
        print *,"  sum(exp_gws)=",exp_gws%gsum(1)
      endif
    endif

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

  subroutine RealizeField4D(field, grid, state, layers, members, rc)
    type(fld4d), intent(inout)  :: field
    type(ESMF_Grid), intent(in) :: grid
    type(ESMF_State)            :: state
    integer, intent(in)         :: layers
    integer, intent(in)         :: members
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
      ungriddedLBound=(/1,1/), ungriddedUBound=(/layers,members/), rc=rc)
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

  subroutine RealizeField3D(field, grid, state, layers, members, rc)
    type(fld3d), intent(inout)    :: field
    type(ESMF_Grid), intent(in)   :: grid
    type(ESMF_State)              :: state
    integer, intent(in), optional :: layers
    integer, intent(in), optional :: members
    integer, intent(out)          :: rc
    ! local variables

    rc = ESMF_SUCCESS

    if (associated(field%fld)) then
      call ESMF_LogSetError(ESMF_RC_MEM_ALLOCATE, &
        msg="Field already associated: "//trim(field%std), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
    allocate(field%fld)
    if (present(layers)) then
      field%fld = ESMF_FieldCreate(name=trim(field%std), grid=grid, &
        typekind=ESMF_TYPEKIND_R4, gridToFieldMap=(/1,3/), &
        ungriddedLBound=(/1/), ungriddedUBound=(/layers/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    elseif (present(members)) then
      field%fld = ESMF_FieldCreate(name=trim(field%std), grid=grid, &
        typekind=ESMF_TYPEKIND_R4, gridToFieldMap=(/1,2/), &
        ungriddedLBound=(/1/), ungriddedUBound=(/members/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      call ESMF_LogSetError(ESMF_RC_MEM_ALLOCATE, &
        msg="Must provide layers or members: "//trim(field%std), &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
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

  subroutine SumField4D(field, vm, rc)
    type(fld4d), intent(inout) :: field
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

end module lis_nuopc_test_gwr
