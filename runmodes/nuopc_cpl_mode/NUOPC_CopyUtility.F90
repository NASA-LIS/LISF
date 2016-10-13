#define ESMF_STDERRORCHECK(rc) ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)
#define FILENAME "NUOPC_CopyUtility.F90"
#define MODNAME "NUOPC_CopyUtility"

module NUOPC_CopyUtility
  use ESMF
  use NUOPC

  implicit none

  private

  public :: NUOPC_CopyFieldToFarray
  public :: NUOPC_CopyFarrayToField
  public :: NUOPC_CopyArrayToFarray
  public :: NUOPC_CopyFarrayToArray

  interface NUOPC_CopyFieldToFarray
    module procedure NUOPC_CopyFieldToFarray_I41D
    module procedure NUOPC_CopyFieldToFarray_I42D
    module procedure NUOPC_CopyFieldToFarray_I43D
    module procedure NUOPC_CopyFieldToFarray_I81D
    module procedure NUOPC_CopyFieldToFarray_I82D
    module procedure NUOPC_CopyFieldToFarray_I83D
    module procedure NUOPC_CopyFieldToFarray_R41D
    module procedure NUOPC_CopyFieldToFarray_R42D
    module procedure NUOPC_CopyFieldToFarray_R43D
    module procedure NUOPC_CopyFieldToFarray_R81D
    module procedure NUOPC_CopyFieldToFarray_R82D
    module procedure NUOPC_CopyFieldToFarray_R83D
  end interface

  interface NUOPC_CopyFarrayToField
    module procedure NUOPC_CopyFarrayToField_I41D
    module procedure NUOPC_CopyFarrayToField_I42D
    module procedure NUOPC_CopyFarrayToField_I43D
    module procedure NUOPC_CopyFarrayToField_I81D
    module procedure NUOPC_CopyFarrayToField_I82D
    module procedure NUOPC_CopyFarrayToField_I83D
    module procedure NUOPC_CopyFarrayToField_R41D
    module procedure NUOPC_CopyFarrayToField_R42D
    module procedure NUOPC_CopyFarrayToField_R43D
    module procedure NUOPC_CopyFarrayToField_R81D
    module procedure NUOPC_CopyFarrayToField_R82D
    module procedure NUOPC_CopyFarrayToField_R83D
  end interface

  interface NUOPC_CopyArrayToFarray
    module procedure NUOPC_CopyArrayToFarray_I41D
    module procedure NUOPC_CopyArrayToFarray_I42D
    module procedure NUOPC_CopyArrayToFarray_I43D
    module procedure NUOPC_CopyArrayToFarray_I81D
    module procedure NUOPC_CopyArrayToFarray_I82D
    module procedure NUOPC_CopyArrayToFarray_I83D
    module procedure NUOPC_CopyArrayToFarray_R41D
    module procedure NUOPC_CopyArrayToFarray_R42D
    module procedure NUOPC_CopyArrayToFarray_R43D
    module procedure NUOPC_CopyArrayToFarray_R81D
    module procedure NUOPC_CopyArrayToFarray_R82D
    module procedure NUOPC_CopyArrayToFarray_R83D
  end interface

  interface NUOPC_CopyFarrayToArray
    module procedure NUOPC_CopyFarrayToArray_I41D
    module procedure NUOPC_CopyFarrayToArray_I42D
    module procedure NUOPC_CopyFarrayToArray_I43D
    module procedure NUOPC_CopyFarrayToArray_I81D
    module procedure NUOPC_CopyFarrayToArray_I82D
    module procedure NUOPC_CopyFarrayToArray_I83D
    module procedure NUOPC_CopyFarrayToArray_R41D
    module procedure NUOPC_CopyFarrayToArray_R42D
    module procedure NUOPC_CopyFarrayToArray_R43D
    module procedure NUOPC_CopyFarrayToArray_R81D
    module procedure NUOPC_CopyFarrayToArray_R82D
    module procedure NUOPC_CopyFarrayToArray_R83D
  end interface

  interface NUOPC_CopyArrayFarray
    module procedure NUOPC_CopyArrayFarray_I41D
    module procedure NUOPC_CopyArrayFarray_I42D
    module procedure NUOPC_CopyArrayFarray_I43D
    module procedure NUOPC_CopyArrayFarray_I81D
    module procedure NUOPC_CopyArrayFarray_I82D
    module procedure NUOPC_CopyArrayFarray_I83D
    module procedure NUOPC_CopyArrayFarray_R41D
    module procedure NUOPC_CopyArrayFarray_R42D
    module procedure NUOPC_CopyArrayFarray_R43D
    module procedure NUOPC_CopyArrayFarray_R81D
    module procedure NUOPC_CopyArrayFarray_R82D
    module procedure NUOPC_CopyArrayFarray_R83D
  end interface

  character(len=*),parameter :: NUOPC_COPY_FWD = 'from ESMF_Array to FORTRAN array'
  character(len=*),parameter :: NUOPC_COPY_BWD = 'from FORTRAN array to ESMF_Array'

contains

  !-----------------------------------------------------------------------------
  ! Copy ESMF Field to FORTRAN Array
  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_I41D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array    

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_I42D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_I43D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_I81D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_I82D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_I83D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_R41D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_R42D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_R43D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_R81D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_R82D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFieldToFarray_R83D(field,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyArrayToFarray(array,farray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------
  ! Copy FORTRAN Array to ESMF Field
  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_I41D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array    

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_I42D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_I43D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_I81D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_I82D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_I83D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_R41D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_R42D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_R43D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_R81D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_R82D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToField_R83D(farray,field,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Field), intent(in)                :: field
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_Array)  :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CopyFarrayToArray(farray,array,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------
  ! Copy ESMF Array to FORTRAN Array.
  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_I41D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_I42D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_I43D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_I81D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_I82D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_I83D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_R41D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_R42D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_R43D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_R81D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_R82D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayToFarray_R83D(array,farray,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFarray(array,farray,reverse=.FALSE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------
  ! Copy FORTRAN Array to ESMF Array
  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_I41D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_I42D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc
    
    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_I43D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_I81D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_I82D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_I83D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_R41D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_R42D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_R43D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_R81D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_R82D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyFarrayToArray_R83D(farray,array,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:,:,:)
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    if (present(rc)) rc = ESMF_SUCCESS

    call NUOPC_CopyArrayFArray(array,farray,reverse=.TRUE.,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------
  ! Copy ESMF Array to/from FORTRAN Array - Direction controlled by reverse
  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_I41D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_I4),pointer  :: srcarray(:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 1) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_I4) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_I42D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:,:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_I4),pointer  :: srcarray(:,:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 2) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_I4) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_I43D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I4),pointer,intent(inout) :: farray(:,:,:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_I4),pointer  :: srcarray(:,:,:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 3) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_I4) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_I81D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_I8),pointer  :: srcarray(:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 1) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_I8) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_I82D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:,:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_I8),pointer  :: srcarray(:,:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 2) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_I8) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_I83D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    integer(ESMF_KIND_I8),pointer,intent(inout) :: farray(:,:,:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_I8),pointer  :: srcarray(:,:,:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 3) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_I8) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_R41D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_R4),pointer  :: srcarray(:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 1) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_R4) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_R42D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:,:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_R4),pointer  :: srcarray(:,:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 2) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_R4) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_R43D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R4),pointer,intent(inout)    :: farray(:,:,:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_R4),pointer  :: srcarray(:,:,:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 3) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_R4) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_R81D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_R8),pointer  :: srcarray(:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 1) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_R8) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_R82D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:,:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_R8),pointer  :: srcarray(:,:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 2) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_R8) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_CopyArrayFarray_R83D(array,farray,reverse,localDe,rc)
   ! ARGUMENTS
    type(ESMF_Array), intent(in)                :: array
    real(ESMF_KIND_R8),pointer,intent(inout)    :: farray(:,:,:)
    logical,intent(in)                          :: reverse
    integer,intent(in),optional                 :: localDe
    integer, intent(out),optional               :: rc

    ! LOCAL VARIABLES
    type(ESMF_TypeKind_Flag)       :: typekind
    integer                        :: rank
    integer                        :: localDeCount
    integer(ESMF_KIND_R8),pointer  :: srcarray(:,:,:)
    character(len=32),pointer      :: direction

    if (present(rc)) rc = ESMF_SUCCESS

    if (reverse) then
      direction = NUOPC_COPY_BWD
    else
      direction = NUOPC_COPY_FWD
    endif

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (localDeCount > 1 .AND. (.NOT. present(localDe))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_INCOMP,   &
        msg="Cannot copy "//trim(direction)//" without localDe.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (rank /= 3) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot copy "//trim(direction)//" unless ranks match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if (typekind /= ESMF_TYPEKIND_R8) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SAMETYPE,   &
        msg="Cannot copy "//trim(direction)//" unless typekinds match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    call ESMF_ArrayGet(array,farrayPtr=srcarray,localDe=localDe,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (size(srcarray) == size(farray)) then
      if (reverse) then
        srcarray = farray
      else
        farray = srcarray
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Cannot copy "//trim(direction)//" unless array sizes match.", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

end module
