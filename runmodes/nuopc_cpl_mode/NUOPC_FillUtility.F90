#define ESMF_STDERRORCHECK(rc) ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)
#define FILENAME "NUOPC_FillUtility.F90"
#define MODNAME "NUOPC_FillUtility"

module NUOPC_FillUtility
  use ESMF
  use NUOPC

  implicit none

  private

  public :: NUOPC_FillField
  public :: NUOPC_FillArray
  public :: NUOPC_FillFieldBundle
  public :: NUOPC_FillState

  interface NUOPC_FillState
    module procedure NUOPC_FillState_I4
    module procedure NUOPC_FillState_I8
    module procedure NUOPC_FillState_R4
    module procedure NUOPC_FillState_R8
    module procedure NUOPC_FillState_SCHEME
  end interface

  interface NUOPC_FillFieldBundle
    module procedure NUOPC_FillFieldBundle_I4
    module procedure NUOPC_FillFieldBundle_I8
    module procedure NUOPC_FillFieldBundle_R4
    module procedure NUOPC_FillFieldBundle_R8
    module procedure NUOPC_FillFieldBundle_SCHEME
  end interface

  interface NUOPC_FillField
    module procedure NUOPC_FillField_I4
    module procedure NUOPC_FillField_I8
    module procedure NUOPC_FillField_R4
    module procedure NUOPC_FillField_R8
  end interface

  interface NUOPC_FillArray
    module procedure NUOPC_FillArray_I4
    module procedure NUOPC_FillArray_I8
    module procedure NUOPC_FillArray_R4
    module procedure NUOPC_FillArray_R8
  end interface

  character(len=*),parameter :: NUOPC_COPY_FWD = 'from ESMF_Array to FORTRAN array'
  character(len=*),parameter :: NUOPC_COPY_BWD = 'from FORTRAN array to ESMF_Array'

contains

  !-----------------------------------------------------------------------------
  ! Fill ESMF State
  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillState_I4(state,value,rc)
   ! ARGUMENTS
    type(ESMF_State), intent(in)                :: state
    integer(ESMF_KIND_I4),intent(in)            :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    integer                                :: iIndex
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    integer                                :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate( &
      itemNameList(itemCount), &
      itemTypeList(itemCount), &
      stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item list memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do iIndex = 1, itemCount
      if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(state,field=field, &
          itemName=itemNameList(iIndex),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call NUOPC_FillField(field,value=value,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
    enddo

    deallocate(itemNameList, itemTypeList, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of state item list memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillState_I8(state,value,rc)
   ! ARGUMENTS
    type(ESMF_State), intent(in)                :: state
    integer(ESMF_KIND_I8),intent(in)            :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    integer                                :: iIndex
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    integer                                :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate( &
      itemNameList(itemCount), &
      itemTypeList(itemCount), &
      stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item list memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do iIndex = 1, itemCount
      if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(state,field=field, &
          itemName=itemNameList(iIndex),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call NUOPC_FillField(field,value=value,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
    enddo

    deallocate(itemNameList, itemTypeList, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of state item list memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillState_R4(state,value,rc)
   ! ARGUMENTS
    type(ESMF_State), intent(in)                :: state
    real(ESMF_KIND_R4),intent(in)               :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    integer                                :: iIndex
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    integer                                :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate( &
      itemNameList(itemCount), &
      itemTypeList(itemCount), &
      stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item list memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do iIndex = 1, itemCount
      if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(state,field=field, &
          itemName=itemNameList(iIndex),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call NUOPC_FillField(field,value=value,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
    enddo

    deallocate(itemNameList, itemTypeList, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of state item list memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillState_R8(state,value,rc)
   ! ARGUMENTS
    type(ESMF_State), intent(in)                :: state
    real(ESMF_KIND_R8),intent(in)               :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    integer                                :: iIndex
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    integer                                :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate( &
      itemNameList(itemCount), &
      itemTypeList(itemCount), &
      stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item list memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do iIndex = 1, itemCount
      if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(state,field=field, &
          itemName=itemNameList(iIndex),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call NUOPC_FillField(field,value=value,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
    enddo

    deallocate(itemNameList, itemTypeList, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of state item list memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillState_SCHEME(state,dataFillScheme,step,rc)
   ! ARGUMENTS
    type(ESMF_State), intent(in)                :: state
    character(len=*), intent(in)                :: dataFillScheme
    integer, intent(in), optional               :: step
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    integer                                :: k
    integer                                :: iIndex
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    integer                                :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate( &
      itemNameList(itemCount), &
      itemTypeList(itemCount), &
      stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item list memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    k=1 ! initialize
    do iIndex = 1, itemCount
      if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(state,field=field, &
          itemName=itemNameList(iIndex),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call ESMF_FieldFill(field, dataFillScheme=dataFillScheme, &
          member=k, step=step, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        k=k+1 ! increment the member counter
      endif
    enddo

    deallocate(itemNameList, itemTypeList, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of state item list memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------
  ! Fill ESMF Field Bundle
  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillFieldBundle_I4(fieldbundle,value,rc)
   ! ARGUMENTS
    type(ESMF_FieldBundle), intent(in)          :: fieldbundle
    integer(ESMF_KIND_I4),intent(in)            :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    integer                         :: fIndex
    integer                         :: fieldCount
    type(ESMF_Field),pointer        :: fieldList(:)
    integer                         :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(fieldbundle, &
      fieldCount=fieldCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate(fieldList(fieldCount),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of field lists failed.", &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) &
      return  ! bail out

    call ESMF_FieldBundleGet(fieldbundle, &
      fieldList=fieldList, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do fIndex=1,fieldCount
      call NUOPC_FillField(fieldList(fIndex),value=value,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    enddo

    deallocate(fieldList,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of field lists failed.", &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillFieldBundle_I8(fieldbundle,value,rc)
   ! ARGUMENTS
    type(ESMF_FieldBundle), intent(in)          :: fieldbundle
    integer(ESMF_KIND_I8),intent(in)            :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    integer                         :: fIndex
    integer                         :: fieldCount
    type(ESMF_Field),pointer        :: fieldList(:)
    integer                         :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(fieldbundle, &
      fieldCount=fieldCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate(fieldList(fieldCount),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of field lists failed.", &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) &
      return  ! bail out

    call ESMF_FieldBundleGet(fieldbundle, &
      fieldList=fieldList, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do fIndex=1,fieldCount
      call NUOPC_FillField(fieldList(fIndex),value=value,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    enddo

    deallocate(fieldList,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of field lists failed.", &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillFieldBundle_R4(fieldbundle,value,rc)
   ! ARGUMENTS
    type(ESMF_FieldBundle), intent(in)          :: fieldbundle
    real(ESMF_KIND_R4),intent(in)               :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    integer                         :: fIndex
    integer                         :: fieldCount
    type(ESMF_Field),pointer        :: fieldList(:)
    integer                         :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(fieldbundle, &
      fieldCount=fieldCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate(fieldList(fieldCount),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of field lists failed.", &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) &
      return  ! bail out

    call ESMF_FieldBundleGet(fieldbundle, &
      fieldList=fieldList, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do fIndex=1,fieldCount
      call NUOPC_FillField(fieldList(fIndex),value=value,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    enddo

    deallocate(fieldList,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of field lists failed.", &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillFieldBundle_R8(fieldbundle,value,rc)
   ! ARGUMENTS
    type(ESMF_FieldBundle), intent(in)          :: fieldbundle
    real(ESMF_KIND_R8),intent(in)               :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    integer                         :: fIndex
    integer                         :: fieldCount
    type(ESMF_Field),pointer        :: fieldList(:)
    integer                         :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(fieldbundle, &
      fieldCount=fieldCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate(fieldList(fieldCount),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of field lists failed.", &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) &
      return  ! bail out

    call ESMF_FieldBundleGet(fieldbundle, &
      fieldList=fieldList, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do fIndex=1,fieldCount
      call NUOPC_FillField(fieldList(fIndex),value=value,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    enddo

    deallocate(fieldList,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of field lists failed.", &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillFieldBundle_SCHEME(fieldbundle,dataFillScheme,step,rc)
   ! ARGUMENTS
    type(ESMF_FieldBundle), intent(in)          :: fieldbundle
    character(len=*), intent(in)                :: dataFillScheme
    integer, intent(in), optional               :: step
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    integer                         :: fIndex
    integer                         :: fieldCount
    type(ESMF_Field),pointer        :: fieldList(:)
    integer                         :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(fieldbundle, &
      fieldCount=fieldCount, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    allocate(fieldList(fieldCount),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of field lists failed.", &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) &
      return  ! bail out

    call ESMF_FieldBundleGet(fieldbundle, &
      fieldList=fieldList, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    do fIndex=1,fieldCount
      call ESMF_FieldFill(fieldList(fIndex), dataFillScheme=dataFillScheme, &
        member=fIndex, step=step, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return 
    enddo

    deallocate(fieldList,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of field lists failed.", &
      line=__LINE__, file=FILENAME, rcToReturn=rc)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------
  ! Fill ESMF Field to FORTRAN Array
  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillField_I4(field,value,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I4),intent(in)            :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array    

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_FillArray(array,value=value,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  subroutine NUOPC_FillField_I8(field,value,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I8),intent(in)            :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_FillArray(array,value=value,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  subroutine NUOPC_FillField_R4(field,value,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R4),intent(in)               :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_FillArray(array,value=value,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  subroutine NUOPC_FillField_R8(field,value,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R8),intent(in)               :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_FillArray(array,value=value,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------
  ! Fill ESMF Array
  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillArray_I4(array,value,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I4),intent(in)            :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VALUES
    integer(ESMF_KIND_I4),pointer :: farray_I41D(:)
    integer(ESMF_KIND_I4),pointer :: farray_I42D(:,:)
    integer(ESMF_KIND_I4),pointer :: farray_I43D(:,:,:)
    integer(ESMF_KIND_I8),pointer :: farray_I81D(:)
    integer(ESMF_KIND_I8),pointer :: farray_I82D(:,:)
    integer(ESMF_KIND_I8),pointer :: farray_I83D(:,:,:)
    real(ESMF_KIND_R4),pointer    :: farray_R41D(:)
    real(ESMF_KIND_R4),pointer    :: farray_R42D(:,:)
    real(ESMF_KIND_R4),pointer    :: farray_R43D(:,:,:)
    real(ESMF_KIND_R8),pointer    :: farray_R81D(:)
    real(ESMF_KIND_R8),pointer    :: farray_R82D(:,:)
    real(ESMF_KIND_R8),pointer    :: farray_R83D(:,:,:)
    type(ESMF_TypeKind_Flag)      :: typekind
    integer                       :: rank
    integer                       :: localDeCount
    integer                       :: deIndex

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (rank == 1) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I41D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I41D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I81D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I81D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R41D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R41D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R81D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R81D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    elseif (rank == 2) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I42D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I42D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I82D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I82D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R42D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R42D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R82D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R82D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    elseif (rank == 3) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I43D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I43D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I83D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I83D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R43D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R43D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R83D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R83D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot fill ESMF Array because rank is not supported", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillArray_I8(array,value,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I8),intent(in)            :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VALUES
    integer(ESMF_KIND_I4),pointer :: farray_I41D(:)
    integer(ESMF_KIND_I4),pointer :: farray_I42D(:,:)
    integer(ESMF_KIND_I4),pointer :: farray_I43D(:,:,:)
    integer(ESMF_KIND_I8),pointer :: farray_I81D(:)
    integer(ESMF_KIND_I8),pointer :: farray_I82D(:,:)
    integer(ESMF_KIND_I8),pointer :: farray_I83D(:,:,:)
    real(ESMF_KIND_R4),pointer    :: farray_R41D(:)
    real(ESMF_KIND_R4),pointer    :: farray_R42D(:,:)
    real(ESMF_KIND_R4),pointer    :: farray_R43D(:,:,:)
    real(ESMF_KIND_R8),pointer    :: farray_R81D(:)
    real(ESMF_KIND_R8),pointer    :: farray_R82D(:,:)
    real(ESMF_KIND_R8),pointer    :: farray_R83D(:,:,:)
    type(ESMF_TypeKind_Flag)      :: typekind
    integer                       :: rank
    integer                       :: localDeCount
    integer                       :: deIndex

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (rank == 1) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I41D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I41D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I81D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I81D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R41D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R41D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R81D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R81D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    elseif (rank == 2) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I42D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I42D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I82D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I82D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R42D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R42D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R82D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R82D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    elseif (rank == 3) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I43D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I43D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I83D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I83D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R43D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R43D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R83D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R83D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot fill ESMF Array because rank is not supported", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillArray_R4(array,value,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R4),intent(in)               :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VALUES
    integer(ESMF_KIND_I4),pointer :: farray_I41D(:)
    integer(ESMF_KIND_I4),pointer :: farray_I42D(:,:)
    integer(ESMF_KIND_I4),pointer :: farray_I43D(:,:,:)
    integer(ESMF_KIND_I8),pointer :: farray_I81D(:)
    integer(ESMF_KIND_I8),pointer :: farray_I82D(:,:)
    integer(ESMF_KIND_I8),pointer :: farray_I83D(:,:,:)
    real(ESMF_KIND_R4),pointer    :: farray_R41D(:)
    real(ESMF_KIND_R4),pointer    :: farray_R42D(:,:)
    real(ESMF_KIND_R4),pointer    :: farray_R43D(:,:,:)
    real(ESMF_KIND_R8),pointer    :: farray_R81D(:)
    real(ESMF_KIND_R8),pointer    :: farray_R82D(:,:)
    real(ESMF_KIND_R8),pointer    :: farray_R83D(:,:,:)
    type(ESMF_TypeKind_Flag)      :: typekind
    integer                       :: rank
    integer                       :: localDeCount
    integer                       :: deIndex

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (rank == 1) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I41D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I41D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I81D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I81D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R41D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R41D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R81D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R81D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    elseif (rank == 2) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I42D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I42D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I82D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I82D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R42D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R42D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R82D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R82D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    elseif (rank == 3) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I43D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I43D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I83D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I83D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R43D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R43D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R83D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R83D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot fill ESMF Array because rank is not supported", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_FillArray_R8(array,value,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R8),intent(in)               :: value
    integer, intent(out),optional               :: rc

    ! LOCAL VALUES
    integer(ESMF_KIND_I4),pointer :: farray_I41D(:)
    integer(ESMF_KIND_I4),pointer :: farray_I42D(:,:)
    integer(ESMF_KIND_I4),pointer :: farray_I43D(:,:,:)
    integer(ESMF_KIND_I8),pointer :: farray_I81D(:)
    integer(ESMF_KIND_I8),pointer :: farray_I82D(:,:)
    integer(ESMF_KIND_I8),pointer :: farray_I83D(:,:,:)
    real(ESMF_KIND_R4),pointer    :: farray_R41D(:)
    real(ESMF_KIND_R4),pointer    :: farray_R42D(:,:)
    real(ESMF_KIND_R4),pointer    :: farray_R43D(:,:,:)
    real(ESMF_KIND_R8),pointer    :: farray_R81D(:)
    real(ESMF_KIND_R8),pointer    :: farray_R82D(:,:)
    real(ESMF_KIND_R8),pointer    :: farray_R83D(:,:,:)
    type(ESMF_TypeKind_Flag)      :: typekind
    integer                       :: rank
    integer                       :: localDeCount
    integer                       :: deIndex

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (rank == 1) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I41D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I41D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I81D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I81D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R41D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R41D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R81D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R81D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    elseif (rank == 2) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I42D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I42D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I82D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I82D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R42D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R42D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R82D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R82D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    elseif (rank == 3) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I43D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I43D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I83D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_I83D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R43D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R43D = value
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R83D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          farray_R83D = value
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot fill ESMF Array because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot fill ESMF Array because rank is not supported", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

end module
