#define ESMF_STDERRORCHECK(rc) ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)
#define FILENAME "beta_NUOPC_Fill.F90"
#define MODNAME "beta_NUOPC_Fill"

module beta_NUOPC_Fill
  use ESMF
  use NUOPC

  implicit none

  private

  public :: beta_ESMF_FieldFill
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

!------------------------------------------------------------------------------


#define ESMF_FILE "beta_NUOPC_Fill.F90"
#define ESMF_METHOD "beta_ESMF_FieldFill()"
#define ESMF_ERR_PASSTHRU msg="Internal subroutine call returned Error"
!BOP
! !IROUTINE: beta_ESMF_FieldFill - Fill data into a Field
! !INTERFACE:
  subroutine beta_ESMF_FieldFill(field, keywordEnforcer, &
    dataFillScheme, member, step, amplitude, meanValue, rc)
! !ARGUMENTS:
    type(ESMF_Field), intent(inout)           :: field
type(ESMF_KeywordEnforcer), optional:: keywordEnforcer ! must use keywords below
    character(len=*), intent(in), optional    :: dataFillScheme
    integer, intent(in), optional             :: member
    integer, intent(in), optional             :: step
    real, intent(in), optional                :: amplitude
    real, intent(in), optional                :: meanValue
    integer, intent(out), optional            :: rc
! !DESCRIPTION:
!   \label{ESMF_FieldFill}
!   Fill {\tt field} with data according to {\tt dataFillScheme}. Depending
!   on the chosen fill scheme, the {\tt member} and {\tt step} arguments are
!   used to provide differing fill data patterns.
!
!   The arguments are:
!   \begin{description}
!   \item[field]
!     The {\tt ESMF\_Field} object to fill with data.
!   \item[{[dataFillScheme]}]
!     The fill scheme. The available options are "sincos", and "one".
!     Defaults to "sincos".
!   \item[{[member]}]
!     Member incrementor. Defaults to 1.
!   \item[{[step]}]
!     Step incrementor. Defaults to 1.
!   \item[{[amplitude]}]
!     Magnitude of change. Defaults to 1.
!   \item[{[meanValue]}]
!     Mean value. Defaults to 0.
!   \item[{[rc]}]
!     Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    type(ESMF_Grid)                 :: grid
    type(ESMF_TypeKind_Flag)        :: typekind
    type(ESMF_TypeKind_Flag)        :: coordTypeKind
    integer                         :: rank
    integer, allocatable            :: coordDimCount(:)
    real(ESMF_KIND_R8), pointer     :: dataPtrR8D1(:)
    real(ESMF_KIND_R8), pointer     :: dataPtrR8D2(:,:)
    real(ESMF_KIND_R8), pointer     :: dataPtrR8D3(:,:,:)
    real(ESMF_KIND_R4), pointer     :: dataPtrR4D1(:)
    real(ESMF_KIND_R4), pointer     :: dataPtrR4D2(:,:)
    real(ESMF_KIND_R4), pointer     :: dataPtrR4D3(:,:,:)
    real(ESMF_KIND_R8), pointer     :: coord1PtrR8D1(:)
    real(ESMF_KIND_R8), pointer     :: coord2PtrR8D1(:)
    real(ESMF_KIND_R8), pointer     :: coord1PtrR8D2(:,:)
    real(ESMF_KIND_R8), pointer     :: coord2PtrR8D2(:,:)
    real(ESMF_KIND_R8), pointer     :: coord1PtrR8D3(:,:,:)
    real(ESMF_KIND_R8), pointer     :: coord2PtrR8D3(:,:,:)
    real(ESMF_KIND_R8), pointer     :: coord3PtrR8D3(:,:,:)
    real(ESMF_KIND_R4), pointer     :: coord1PtrR4D1(:)
    real(ESMF_KIND_R4), pointer     :: coord2PtrR4D1(:)
    real(ESMF_KIND_R4), pointer     :: coord1PtrR4D2(:,:)
    real(ESMF_KIND_R4), pointer     :: coord2PtrR4D2(:,:)
    real(ESMF_KIND_R4), pointer     :: coord1PtrR4D3(:,:,:)
    real(ESMF_KIND_R4), pointer     :: coord2PtrR4D3(:,:,:)
    real(ESMF_KIND_R4), pointer     :: coord3PtrR4D3(:,:,:)
    integer                         :: i, j, k

    integer                         :: l_member, l_step
    real                            :: l_amplitude, l_meanValue
    character(len=16)               :: l_dataFillScheme

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, typekind=typekind, rank=rank, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
      line=__LINE__, &
      file=ESMF_FILE)) &
      return  ! bail out

    l_member = 1
    if(present(member)) l_member = member
    l_step = 1
    if(present(step)) l_step = step
    l_dataFillScheme = "sincos"
    if(present(dataFillScheme)) l_dataFillScheme = dataFillScheme
    l_amplitude = 1.0
    if(present(amplitude)) l_amplitude = amplitude
    l_meanValue = 0.0
    if(present(meanValue)) l_meanValue = meanValue

    allocate(coordDimCount(rank))

    if (trim(l_dataFillScheme)=="sincos") then
      call ESMF_FieldGet(field, grid=grid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
        line=__LINE__, &
        file=ESMF_FILE)) &
        return  ! bail out
      call ESMF_GridGet(grid,coordTypeKind=coordTypeKind, &
        coordDimCount=coordDimCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
        line=__LINE__, &
        file=ESMF_FILE)) &
        return  ! bail out
      if (rank==1) then
        ! 1D sin pattern
        ! TODO: support Meshes
        call ESMF_GridGetCoord(grid, coordDim=1, farrayPtr=coord1PtrR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
        if (typekind==ESMF_TYPEKIND_R4) then
          call ESMF_FieldGet(field, farrayPtr=dataPtrR4D1, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
            line=__LINE__, &
            file=ESMF_FILE)) &
            return  ! bail out
          do i=lbound(dataPtrR4D1,1),ubound(dataPtrR4D1,1)
            dataPtrR4D1(i) = &
              (sin(real(l_member)*3.1416*(coord1PtrR8D1(i)+real(l_step))/180.)) * &
              l_amplitude+l_meanValue
          enddo
        elseif (typekind==ESMF_TYPEKIND_R8) then
          call ESMF_FieldGet(field, farrayPtr=dataPtrR8D1, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
            line=__LINE__, &
            file=ESMF_FILE)) &
            return  ! bail out
          do i=lbound(dataPtrR8D1,1),ubound(dataPtrR8D1,1)
            dataPtrR8D1(i) = &
              (sin(real(l_member)*3.1416*(coord1PtrR8D1(i)+real(l_step))/180.)) * &
              l_amplitude+l_meanValue
          enddo
        else
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="Unsupported typekind-rank and scheme combination requested.", &
            line=__LINE__, &
            file=ESMF_FILE, &
            rcToReturn=rc)
          return ! bail out
        endif
      elseif (rank==2) then
        ! 2D sin*cos pattern
        ! TODO: support Meshes
        if (coordTypeKind==ESMF_TYPEKIND_R4) then
          if (coordDimCount(1)==1) then
            call ESMF_GridGetCoord(grid, coordDim=1, farrayPtr=coord1PtrR4D1, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
              line=__LINE__, &
              file=ESMF_FILE)) &
              return  ! bail out
          else
            ! assume the only other choice here is 2D, if not will trigger error
            call ESMF_GridGetCoord(grid, coordDim=1, farrayPtr=coord1PtrR4D2, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
              line=__LINE__, &
              file=ESMF_FILE)) &
              return  ! bail out
          endif
          if (coordDimCount(2)==1) then
            call ESMF_GridGetCoord(grid, coordDim=2, farrayPtr=coord2PtrR4D1, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
              line=__LINE__, &
              file=ESMF_FILE)) &
              return  ! bail out
          else
            ! assume the only other choice here is 2D, if not will trigger error
            call ESMF_GridGetCoord(grid, coordDim=2, farrayPtr=coord2PtrR4D2, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
              line=__LINE__, &
              file=ESMF_FILE)) &
              return  ! bail out
          endif
        elseif (coordTypeKind==ESMF_TYPEKIND_R8) then
          if (coordDimCount(1)==1) then
            call ESMF_GridGetCoord(grid, coordDim=1, farrayPtr=coord1PtrR8D1, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
              line=__LINE__, &
              file=ESMF_FILE)) &
              return  ! bail out
          else
            ! assume the only other choice here is 2D, if not will trigger error
            call ESMF_GridGetCoord(grid, coordDim=1, farrayPtr=coord1PtrR8D2, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
              line=__LINE__, &
              file=ESMF_FILE)) &
              return  ! bail out
          endif
          if (coordDimCount(2)==1) then
            call ESMF_GridGetCoord(grid, coordDim=2, farrayPtr=coord2PtrR8D1, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
              line=__LINE__, &
              file=ESMF_FILE)) &
              return  ! bail out
          else
            ! assume the only other choice here is 2D, if not will trigger error
            call ESMF_GridGetCoord(grid, coordDim=2, farrayPtr=coord2PtrR8D2, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
              line=__LINE__, &
              file=ESMF_FILE)) &
              return  ! bail out
          endif
        else
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="Unsupported coordinate typekind.", &
            line=__LINE__, &
            file=ESMF_FILE, &
            rcToReturn=rc)
          return ! bail out
        endif

        if (typekind==ESMF_TYPEKIND_R4) then
          call ESMF_FieldGet(field, farrayPtr=dataPtrR4D2, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
            line=__LINE__, &
            file=ESMF_FILE)) &
            return  ! bail out
          if (coordDimCount(1)==1 .and. coordDimCount(2)==1) then
            if (coordTypeKind==ESMF_TYPEKIND_R4) then
            do j=lbound(dataPtrR4D2,2),ubound(dataPtrR4D2,2)
            do i=lbound(dataPtrR4D2,1),ubound(dataPtrR4D2,1)
              dataPtrR4D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR4D1(i)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR4D1(j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            elseif (coordTypeKind==ESMF_TYPEKIND_R8) then
            do j=lbound(dataPtrR4D2,2),ubound(dataPtrR4D2,2)
            do i=lbound(dataPtrR4D2,1),ubound(dataPtrR4D2,1)
              dataPtrR4D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR8D1(i)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR8D1(j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            endif
          else if (coordDimCount(1)==2 .and. coordDimCount(2)==1) then
            if (coordTypeKind==ESMF_TYPEKIND_R4) then
            do j=lbound(dataPtrR4D2,2),ubound(dataPtrR4D2,2)
            do i=lbound(dataPtrR4D2,1),ubound(dataPtrR4D2,1)
              dataPtrR4D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR4D2(i,j)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR4D1(j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            elseif (coordTypeKind==ESMF_TYPEKIND_R8) then
            do j=lbound(dataPtrR4D2,2),ubound(dataPtrR4D2,2)
            do i=lbound(dataPtrR4D2,1),ubound(dataPtrR4D2,1)
              dataPtrR4D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR8D2(i,j)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR8D1(j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            endif
          else if (coordDimCount(1)==1 .and. coordDimCount(2)==2) then
            if (coordTypeKind==ESMF_TYPEKIND_R4) then
            do j=lbound(dataPtrR4D2,2),ubound(dataPtrR4D2,2)
            do i=lbound(dataPtrR4D2,1),ubound(dataPtrR4D2,1)
              dataPtrR4D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR4D1(i)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR4D2(i,j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            elseif (coordTypeKind==ESMF_TYPEKIND_R8) then
            do j=lbound(dataPtrR4D2,2),ubound(dataPtrR4D2,2)
            do i=lbound(dataPtrR4D2,1),ubound(dataPtrR4D2,1)
              dataPtrR4D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR8D1(i)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR8D2(i,j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            endif
           else
            ! only choice left is both 2d coordinate arrays
            if (coordTypeKind==ESMF_TYPEKIND_R4) then
            do j=lbound(dataPtrR4D2,2),ubound(dataPtrR4D2,2)
            do i=lbound(dataPtrR4D2,1),ubound(dataPtrR4D2,1)
              dataPtrR4D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR4D2(i,j)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR4D2(i,j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            elseif (coordTypeKind==ESMF_TYPEKIND_R8) then
            do j=lbound(dataPtrR4D2,2),ubound(dataPtrR4D2,2)
            do i=lbound(dataPtrR4D2,1),ubound(dataPtrR4D2,1)
              dataPtrR4D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR8D2(i,j)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR8D2(i,j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            endif
          endif
        elseif (typekind==ESMF_TYPEKIND_R8) then
          call ESMF_FieldGet(field, farrayPtr=dataPtrR8D2, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
            line=__LINE__, &
            file=ESMF_FILE)) &
            return  ! bail out
          if (coordDimCount(1)==1 .and. coordDimCount(2)==1) then
            if (coordTypeKind==ESMF_TYPEKIND_R4) then
            do j=lbound(dataPtrR8D2,2),ubound(dataPtrR8D2,2)
            do i=lbound(dataPtrR8D2,1),ubound(dataPtrR8D2,1)
              dataPtrR8D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR4D1(i)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR4D1(j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            elseif (coordTypeKind==ESMF_TYPEKIND_R8) then
            do j=lbound(dataPtrR8D2,2),ubound(dataPtrR8D2,2)
            do i=lbound(dataPtrR8D2,1),ubound(dataPtrR8D2,1)
              dataPtrR8D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR8D1(i)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR8D1(j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            endif
          else if (coordDimCount(1)==2 .and. coordDimCount(2)==1) then
            if (coordTypeKind==ESMF_TYPEKIND_R4) then
            do j=lbound(dataPtrR8D2,2),ubound(dataPtrR8D2,2)
            do i=lbound(dataPtrR8D2,1),ubound(dataPtrR8D2,1)
              dataPtrR8D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR4D2(i,j)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR4D1(j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            elseif (coordTypeKind==ESMF_TYPEKIND_R8) then
            do j=lbound(dataPtrR8D2,2),ubound(dataPtrR8D2,2)
            do i=lbound(dataPtrR8D2,1),ubound(dataPtrR8D2,1)
              dataPtrR8D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR8D2(i,j)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR8D1(j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            endif
          else if (coordDimCount(1)==1 .and. coordDimCount(2)==2) then
            if (coordTypeKind==ESMF_TYPEKIND_R4) then
            do j=lbound(dataPtrR8D2,2),ubound(dataPtrR8D2,2)
            do i=lbound(dataPtrR8D2,1),ubound(dataPtrR8D2,1)
              dataPtrR8D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR4D1(i)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR4D2(i,j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            elseif (coordTypeKind==ESMF_TYPEKIND_R8) then
            do j=lbound(dataPtrR8D2,2),ubound(dataPtrR8D2,2)
            do i=lbound(dataPtrR8D2,1),ubound(dataPtrR8D2,1)
              dataPtrR8D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR8D1(i)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR8D2(i,j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            endif
          else
            ! only choice left is both 2d coordinate arrays
            if (coordTypeKind==ESMF_TYPEKIND_R4) then
            do j=lbound(dataPtrR8D2,2),ubound(dataPtrR8D2,2)
            do i=lbound(dataPtrR8D2,1),ubound(dataPtrR8D2,1)
              dataPtrR8D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR4D2(i,j)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR4D2(i,j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            elseif (coordTypeKind==ESMF_TYPEKIND_R8) then
            do j=lbound(dataPtrR8D2,2),ubound(dataPtrR8D2,2)
            do i=lbound(dataPtrR8D2,1),ubound(dataPtrR8D2,1)
              dataPtrR8D2(i,j) = &
                (sin(real(l_member)*3.1416*(coord1PtrR8D2(i,j)+real(l_step))/180.) * &
                cos(real(l_member)*3.1416*(coord2PtrR8D2(i,j)+real(l_step))/180.)) * &
                l_amplitude+l_meanValue
            enddo
            enddo
            endif
          endif
        else
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="Unsupported typekind-rank and scheme combination requested.", &
            line=__LINE__, &
            file=ESMF_FILE, &
            rcToReturn=rc)
          return ! bail out
        endif
      elseif (rank==3) then
        ! 3D sin*cos*sin pattern
        ! TODO: support Meshes
        call ESMF_GridGetCoord(grid, coordDim=1, farrayPtr=coord1PtrR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
        call ESMF_GridGetCoord(grid, coordDim=2, farrayPtr=coord2PtrR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
        call ESMF_GridGetCoord(grid, coordDim=3, farrayPtr=coord3PtrR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
      if (typekind==ESMF_TYPEKIND_R4) then
        call ESMF_FieldGet(field, farrayPtr=dataPtrR4D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
        do k=lbound(dataPtrR4D3,3),ubound(dataPtrR4D3,3)
        do j=lbound(dataPtrR4D3,2),ubound(dataPtrR4D3,2)
        do i=lbound(dataPtrR4D3,1),ubound(dataPtrR4D3,1)
          dataPtrR4D3(i,j,k) = &
            (sin(real(l_member)*3.1416*(coord1PtrR8D3(i,j,k)+real(l_step))/180.) * &
            cos(real(l_member)*3.1416*(coord2PtrR8D3(i,j,k)+real(l_step))/180.) * &
            sin(real(l_member)*3.1416*(coord3PtrR8D3(i,j,k)+real(l_step))/180.)) * &
            l_amplitude+l_meanValue
        enddo
        enddo
        enddo
        elseif (typekind==ESMF_TYPEKIND_R8) then
          call ESMF_FieldGet(field, farrayPtr=dataPtrR8D3, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
            line=__LINE__, &
            file=ESMF_FILE)) &
            return  ! bail out
          do k=lbound(dataPtrR8D3,3),ubound(dataPtrR8D3,3)
          do j=lbound(dataPtrR8D3,2),ubound(dataPtrR8D3,2)
          do i=lbound(dataPtrR8D3,1),ubound(dataPtrR8D3,1)
            dataPtrR8D3(i,j,k) = &
              (sin(real(l_member)*3.1416*(coord1PtrR8D3(i,j,k)+real(l_step))/180.) * &
              cos(real(l_member)*3.1416*(coord2PtrR8D3(i,j,k)+real(l_step))/180.) * &
              sin(real(l_member)*3.1416*(coord3PtrR8D3(i,j,k)+real(l_step))/180.)) * &
              l_amplitude+l_meanValue
          enddo
          enddo
          enddo
        else
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="Unsupported typekind-rank and scheme combination requested.", &
            line=__LINE__, &
            file=ESMF_FILE, &
            rcToReturn=rc)
          return ! bail out
        endif
      else
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg="Unsupported typekind-rank and scheme combination requested.", &
          line=__LINE__, &
          file=ESMF_FILE, &
          rcToReturn=rc)
        return ! bail out
      endif
    else if (trim(dataFillScheme)=="one") then
      if (typekind==ESMF_TYPEKIND_R8 .and. rank==1) then
        ! 1D all 1.
        call ESMF_FieldGet(field, farrayPtr=dataPtrR8D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
        ! initialize the entire array
        dataPtrR8D1 = 1._ESMF_KIND_R8
      elseif (typekind==ESMF_TYPEKIND_R4 .and. rank==1) then
        ! 1D all 1.
        call ESMF_FieldGet(field, farrayPtr=dataPtrR4D1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
        ! initialize the entire array
        dataPtrR4D1 = 1._ESMF_KIND_R4
      elseif (typekind==ESMF_TYPEKIND_R8 .and. rank==2) then
        ! 2D all 1.
        call ESMF_FieldGet(field, farrayPtr=dataPtrR8D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
        ! initialize the entire array
        dataPtrR8D2 = 1._ESMF_KIND_R8
      elseif (typekind==ESMF_TYPEKIND_R4 .and. rank==2) then
        ! 2D all 1.
        call ESMF_FieldGet(field, farrayPtr=dataPtrR4D2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
        ! initialize the entire array
        dataPtrR4D2 = 1._ESMF_KIND_R4
      elseif (typekind==ESMF_TYPEKIND_R8 .and. rank==3) then
        ! 3D all 1.
        call ESMF_FieldGet(field, farrayPtr=dataPtrR8D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
        ! initialize the entire array
        dataPtrR8D3 = 1._ESMF_KIND_R8
      elseif (typekind==ESMF_TYPEKIND_R4 .and. rank==3) then
        ! 3D all 1.
        call ESMF_FieldGet(field, farrayPtr=dataPtrR4D3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, ESMF_ERR_PASSTHRU, &
          line=__LINE__, &
          file=ESMF_FILE)) &
          return  ! bail out
        ! initialize the entire array
        dataPtrR4D3 = 1._ESMF_KIND_R4
      endif
    else
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg="Unknown dataFillScheme requested.", &
        line=__LINE__, &
        file=ESMF_FILE, &
        rcToReturn=rc)
      return ! bail out
    endif

    deallocate(coordDimCount)

  end subroutine
#undef ESMF_FILE
#undef ESMF_METHOD
#undef ESMF_ERR_PASSTHRU
!------------------------------------------------------------------------------

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
