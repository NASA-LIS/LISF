#define ESMF_STDERRORCHECK(rc) ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)
#define FILENAME "NUOPC_FileReadUtility.F90"
#define MODNAME "NUOPC_FileReadUtility"

#define VERBOSITY_MIN 0
#define VERBOSITY_MAX 255
#define VERBOSITY_DBG 1023

module NUOPC_FileReadUtility
  use ESMF
  use NUOPC
  use NETCDF

  implicit none

  private
  public :: NUOPC_NetcdfReadIXJX

  interface NUOPC_NetcdfReadIXJX
    module procedure NUOPC_NetcdfReadIXJX_Field
    module procedure NUOPC_NetcdfReadIXJX_Array
    module procedure NUOPC_NetcdfReadIXJX_I4
    module procedure NUOPC_NetcdfReadIXJX_I8
    module procedure NUOPC_NetcdfReadIXJX_R4
    module procedure NUOPC_NetcdfReadIXJX_R8
  end interface  

contains

  !-----------------------------------------------------------------------------
  ! Utilities
  !-----------------------------------------------------------------------------

  subroutine NUOPC_NetcdfReadIXJX_Field(varname,filename,start,field,rc)
    ! ARGUMENTS
    character(len=*), intent(in)          :: varname
    character(len=*), intent(in)          :: filename
    integer,intent(in)                    :: start(2)
    type(ESMF_Field),intent(inout)        :: field
    integer, intent(out),optional         :: rc

    ! LOCAL VALUES
    type(ESMF_Array)                      :: array

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_NetcdfReadIXJX(varname,filename,start,array=array,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_NetcdfReadIXJX_Array(varname,filename,start,array,rc)
    ! ARGUMENTS
    character(len=*), intent(in)          :: varname
    character(len=*), intent(in)          :: filename
    integer,intent(in)                    :: start(2)
    type(ESMF_Array),intent(inout)        :: array
    integer, intent(out),optional         :: rc

    ! LOCAL VALUES
    integer(ESMF_KIND_I4),pointer :: farray_I42D(:,:)
    integer(ESMF_KIND_I8),pointer :: farray_I82D(:,:)
    real(ESMF_KIND_R4),pointer    :: farray_R42D(:,:)
    real(ESMF_KIND_R8),pointer    :: farray_R82D(:,:)
    type(ESMF_TypeKind_Flag)      :: typekind
    integer                       :: rank
    integer                       :: localDeCount
    integer                       :: deIndex

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array,typekind=typekind,rank=rank,localDeCount=localDeCount,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (rank == 2) then
      if (typekind == ESMF_TYPEKIND_I4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I42D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call NUOPC_NetcdfReadIXJX(varname,filename,start,farray=farray_I42D, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        enddo
      elseif (typekind == ESMF_TYPEKIND_I8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_I82D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call NUOPC_NetcdfReadIXJX(varname,filename,start,farray=farray_I82D, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        enddo
      elseif (typekind == ESMF_TYPEKIND_R4) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R42D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call NUOPC_NetcdfReadIXJX(varname,filename,start,farray=farray_R42D, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        enddo
      elseif (typekind == ESMF_TYPEKIND_R8) then
        do deIndex=0,localDeCount-1
          call ESMF_ArrayGet(array,farrayPtr=farray_R82D,localDe=deIndex,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call NUOPC_NetcdfReadIXJX(varname,filename,start,farray=farray_R82D, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        enddo
      else
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
          msg="Cannot read NetCDF because typekind is not supported", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
      endif
    else
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_RANK,   &
        msg="Cannot read NetCDF because rank is not supported", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif 

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_NetcdfReadIXJX_I4(varname,filename,start,farray,rc)
    ! ARGUMENTS
    character(len=*), intent(in)          :: varname
    character(len=*), intent(in)          :: filename
    integer,intent(in)                    :: start(2)
    integer(ESMF_KIND_I4),intent(inout)   :: farray(:,:)
    integer, intent(out),optional         :: rc

    ! LOCAL VARIABLES
    integer                               :: stat
    integer                               :: varid
    integer, dimension(nf90_max_var_dims) :: dimIDs
    integer                               :: dimCnt(2)
    integer                               :: ncid

    if (present(rc)) rc = ESMF_SUCCESS

    stat = nf90_open(filename,nf90_NoWrite,ncid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_OPEN,   &
        msg="Error opening NetCDF file "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inq_varid(ncid,varname,varid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    ! How big is the netCDF variable, that is, what are the lengths of
    !   its constituent dimensions?
    stat = nf90_inquire_variable(ncid, varid, dimids = dimIDs)
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error loacating variable "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inquire_dimension(ncid, dimIDs(1), len = dimCnt(1))
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable dim1 "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inquire_dimension(ncid, dimIDs(2), len = dimCnt(2))
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable dim2 "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if ((start(1)+size(farray,1)-1 > dimCnt(1)) .OR. &
    (start(2)+size(farray,2)-1 > dimCnt(2))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Error FORTRAN array size mismatch "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_get_var(ncid, varid, values=farray, start=(/start(1),start(2)/), &
      count=(/size(farray,1),size(farray,2)/))
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable values "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_close(ncid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_CLOSE,   &
        msg="Error closing NetCDF file "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_NetcdfReadIXJX_I8(varname,filename,start,farray,rc)
    ! ARGUMENTS
    character(len=*), intent(in)          :: varname
    character(len=*), intent(in)          :: filename
    integer,intent(in)                    :: start(2)
    integer(ESMF_KIND_I8),intent(inout)   :: farray(:,:)
    integer, intent(out),optional         :: rc

    ! LOCAL VARIABLES
    integer                               :: stat
    integer                               :: varid
    integer, dimension(nf90_max_var_dims) :: dimIDs
    integer                               :: dimCnt(2)
    integer                               :: ncid

    if (present(rc)) rc = ESMF_SUCCESS

    stat = nf90_open(filename,nf90_NoWrite,ncid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_OPEN,   &
        msg="Error opening NetCDF file "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inq_varid(ncid,varname,varid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    ! How big is the netCDF variable, that is, what are the lengths of
    !   its constituent dimensions?
    stat = nf90_inquire_variable(ncid, varid, dimids = dimIDs)
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error loacating variable "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inquire_dimension(ncid, dimIDs(1), len = dimCnt(1))
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable dim1 "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inquire_dimension(ncid, dimIDs(2), len = dimCnt(2))
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable dim2 "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if ((start(1)+size(farray,1)-1 > dimCnt(1)) .OR. &
    (start(2)+size(farray,2)-1 > dimCnt(2))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Error FORTRAN array size mismatch "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_get_var(ncid, varid, values=farray, start=(/start(1),start(2)/), &
      count=(/size(farray,1),size(farray,2)/))
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable values "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_close(ncid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_CLOSE,   &
        msg="Error closing NetCDF file "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_NetcdfReadIXJX_R4(varname,filename,start,farray,rc)
    ! ARGUMENTS
    character(len=*), intent(in)          :: varname
    character(len=*), intent(in)          :: filename
    integer,intent(in)                    :: start(2)
    real(ESMF_KIND_R4),intent(inout)      :: farray(:,:)
    integer, intent(out),optional         :: rc

    ! LOCAL VARIABLES
    integer                               :: stat
    integer                               :: varid
    integer, dimension(nf90_max_var_dims) :: dimIDs
    integer                               :: dimCnt(2)
    integer                               :: ncid

    if (present(rc)) rc = ESMF_SUCCESS

    stat = nf90_open(filename,nf90_NoWrite,ncid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_OPEN,   &
        msg="Error opening NetCDF file "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inq_varid(ncid,varname,varid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    ! How big is the netCDF variable, that is, what are the lengths of
    !   its constituent dimensions?
    stat = nf90_inquire_variable(ncid, varid, dimids = dimIDs)
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error loacating variable "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inquire_dimension(ncid, dimIDs(1), len = dimCnt(1))
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable dim1 "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inquire_dimension(ncid, dimIDs(2), len = dimCnt(2))
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable dim2 "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if ((start(1)+size(farray,1)-1 > dimCnt(1)) .OR. &
    (start(2)+size(farray,2)-1 > dimCnt(2))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Error FORTRAN array size mismatch "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_get_var(ncid, varid, values=farray, start=(/start(1),start(2)/), &
      count=(/size(farray,1),size(farray,2)/))
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable values "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_close(ncid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_CLOSE,   &
        msg="Error closing NetCDF file "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NUOPC_NetcdfReadIXJX_R8(varname,filename,start,farray,rc)
    ! ARGUMENTS
    character(len=*), intent(in)          :: varname
    character(len=*), intent(in)          :: filename
    integer,intent(in)                    :: start(2)
    real(ESMF_KIND_R8),intent(inout)      :: farray(:,:)
    integer, intent(out),optional         :: rc

    ! LOCAL VARIABLES
    integer                               :: stat
    integer                               :: varid
    integer, dimension(nf90_max_var_dims) :: dimIDs
    integer                               :: dimCnt(2)
    integer                               :: ncid

    if (present(rc)) rc = ESMF_SUCCESS

    stat = nf90_open(filename,nf90_NoWrite,ncid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_OPEN,   &
        msg="Error opening NetCDF file "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inq_varid(ncid,varname,varid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    ! How big is the netCDF variable, that is, what are the lengths of
    !   its constituent dimensions?
    stat = nf90_inquire_variable(ncid, varid, dimids = dimIDs)
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error loacating variable "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inquire_dimension(ncid, dimIDs(1), len = dimCnt(1))
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable dim1 "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_inquire_dimension(ncid, dimIDs(2), len = dimCnt(2))
    if(stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable dim2 "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    if ((start(1)+size(farray,1)-1 > dimCnt(1)) .OR. &
    (start(2)+size(farray,2)-1 > dimCnt(2))) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE,   &
        msg="Error FORTRAN array size mismatch "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_get_var(ncid, varid, values=farray, start=(/start(1),start(2)/), &
      count=(/size(farray,1),size(farray,2)/))
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_READ,   &
        msg="Error reading variable values "//trim(varname)// &
          " in "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    stat = nf90_close(ncid)
    if (stat /= nf90_NoErr) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_CLOSE,   &
        msg="Error closing NetCDF file "//trim(filename)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

end module
