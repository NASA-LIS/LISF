! $Id$
!
! Earth System Modeling Framework
! Copyright 2002-2016, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!
!==============================================================================
#define FILENAME "src/addon/NUOPC/src/beta_NUOPC_Auxiliary.F90"
!==============================================================================

module beta_NUOPC_Auxiliary

  use ESMF
  use NUOPC

  implicit none

  private

  public beta_NUOPC_Write                      ! method
  public beta_NUOPC_GridWrite                      ! method
  public NUOPC_MAPPRESET_GLOBAL
  public NUOPC_MAPPRESET_CONUS
  public NUOPC_MAPPRESET_IRENE
  public NUOPC_MAPPRESET_FRONTRANGE
!==============================================================================
!
! INTERFACE BLOCKS
!
!==============================================================================

  interface beta_NUOPC_Write
    module procedure beta_NUOPC_FieldBundleWrite
    module procedure beta_NUOPC_StateWrite
  end interface

  interface beta_NUOPC_GridWrite
    module procedure beta_NUOPC_GridWrite_coords
    module procedure beta_NUOPC_GridWrite_preset
    module procedure beta_NUOPC_GridWrite_default
  end interface

  interface beta_NUOPC_NclScriptWrite
    module procedure beta_NUOPC_NclScriptWrite_coords
    module procedure beta_NUOPC_NclScriptWrite_preset
    module procedure beta_NUOPC_NclScriptWrite_default
  end interface

  interface beta_NUOPC_FerretScriptWrite
    module procedure beta_NUOPC_FerretScriptWrite_coords
    module procedure beta_NUOPC_FerretScriptWrite_preset
    module procedure beta_NUOPC_FerretScriptWrite_default
  end interface

!==============================================================================
!
! DERIVED TYPES
!
!==============================================================================

  type MapDesc
    character(len=16) :: name
    logical           :: global
    real              :: minLat
    real              :: maxLat
    real              :: minLon
    real              :: maxLon
  endtype MapDesc

!==============================================================================
!
! PRESET VALUES
!
!==============================================================================

  type(MapDesc),parameter :: &
    NUOPC_MAPPRESET_GLOBAL = MapDesc("GLOBAL",.TRUE.,-90.0,90.0,-180.0,-180.0), &
    NUOPC_MAPPRESET_CONUS = MapDesc("CONUS",.FALSE.,18.0,49.0,235.0,298.0), &
    NUOPC_MAPPRESET_IRENE = MapDesc("IRENE",.FALSE.,10.0,60.0,240.0,320.0), &
    NUOPC_MAPPRESET_FRONTRANGE = MapDesc("FRONTRANGE",.FALSE.,38.5,41.0,-107.0,-103.5)

contains

  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: NUOPC_Write - Write Field data to file
! !INTERFACE:
  ! call using generic interface: NUOPC_Write
  subroutine beta_NUOPC_FieldBundleWrite(fieldbundle, fileName, singleFile, &
    overwrite, status, timeslice, iofmt, relaxedflag, rc)
! !ARGUMENTS:
    type(ESMF_FieldBundle),     intent(in)            :: fieldbundle
    character(*),               intent(in)            :: fileName
    logical,                    intent(in),  optional :: singleFile
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    type(ESMF_IOFmt_Flag),      intent(in),  optional :: iofmt
    logical,                    intent(in),  optional :: relaxedflag
    integer,                    intent(out), optional :: rc
! !DESCRIPTION:
!   Write the data in {\tt field} to {\tt file} under the field's "StandardName"
!   attribute if supported by the {\tt iofmt}.
!
!   The arguments are:
!   \begin{description}
!   \item[fieldBundle]
!     The {\tt ESMF\_FieldBundle} object whose data is to be written.
!   \item[fileName]
!     The name of the file to write to.
!   \item[{[singleFile]}]
!     A logical flag, the default is .true., i.e., all fields in the bundle 
!     are written in one single file. If .false., each field will be written 
!     in separate files; these files are numbered with the name based on the 
!     argument "file". That is, a set of files are named: [file_name]001, 
!     [file_name]002, [file_name]003,...
!   \item[{[overwrite]}]
!      A logical flag, the default is .false., i.e., existing Field data may
!      {\em not} be overwritten. If .true., the
!      data corresponding to each field's name will be
!      be overwritten. If the {\tt timeslice} option is given, only data for
!      the given timeslice may be overwritten.
!      Note that it is always an error to attempt to overwrite a NetCDF
!      variable with data which has a different shape.
!   \item[{[status]}]
!      The file status. Valid options are {\tt ESMF\_FILESTATUS\_NEW},
!      {\tt ESMF\_FILESTATUS\_OLD}, {\tt ESMF\_FILESTATUS\_REPLACE}, and
!      {\tt ESMF\_FILESTATUS\_UNKNOWN} (default).
!   \item[{[timeslice]}]
!     Time slice counter. Must be positive. The behavior of this
!     option may depend on the setting of the {\tt overwrite} flag:
!     \begin{description}
!     \item[{\tt overwrite = .false.}:]\ If the timeslice value is
!     less than the maximum time already in the file, the write will fail.
!     \item[{\tt overwrite = .true.}:]\ Any positive timeslice value is valid.
!     \end{description}
!     By default, i.e. by omitting the {\tt timeslice} argument, no
!     provisions for time slicing are made in the output file,
!     however, if the file already contains a time axis for the variable,
!     a timeslice one greater than the maximum will be written.
!   \item[{[iofmt]}]
!    The IO format.  Valid options are  {\tt ESMF\_IOFMT\_BIN} and
!    {\tt ESMF\_IOFMT\_NETCDF}. If not present, file names with a {\tt .bin}
!    extension will use {\tt ESMF\_IOFMT\_BIN}, and file names with a {\tt .nc}
!    extension will use {\tt ESMF\_IOFMT\_NETCDF}.  Other files default to
!    {\tt ESMF\_IOFMT\_NETCDF}.
!   \item[{[relaxedflag]}]
!     If {\tt .true.}, then no error is returned even if the call cannot write
!     the file due to library limitations. Default is {\tt .false.}.
!   \item[{[rc]}]
!     Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    character(ESMF_MAXSTR)  :: standardName
    logical                 :: ioCapable
    logical                 :: doItFlag

    if (present(rc)) rc = ESMF_SUCCESS

    ioCapable = (ESMF_IO_PIO_PRESENT .and. &
      (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT))

    doItFlag = .true. ! default
    if (present(relaxedFlag)) then
      doItFlag = .not.relaxedflag .or. (relaxedflag.and.ioCapable)
    endif

    if (doItFlag) then

      call ESMF_FieldBundleWrite(fieldbundle, fileName=fileName, &
        singleFile=singleFile, overwrite=overwrite, status=status, &
        timeslice=timeslice, iofmt=iofmt, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=FILENAME)) &
        return  ! bail out

    endif

  end subroutine
  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: NUOPC_Write - Write the Fields within a State to NetCDF files
! !INTERFACE:
  ! call using generic interface: NUOPC_Write
  subroutine beta_NUOPC_StateWrite(state, fieldNameList, fileNamePrefix, singleFile, &
    overwrite, status, timeslice, iofmt, relaxedflag, rc)
! !ARGUMENTS:
    type(ESMF_State),           intent(in)            :: state
    character(len=*),           intent(in),  optional :: fieldNameList(:)
    character(len=*),           intent(in),  optional :: fileNamePrefix
    logical,                    intent(in),  optional :: singleFile
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    type(ESMF_IOFmt_Flag),      intent(in),  optional :: iofmt
    logical,                    intent(in),  optional :: relaxedflag
    integer,                    intent(out), optional :: rc
! !DESCRIPTION:
!   Write the data of the fields within a {\tt state} to NetCDF files. Each
!   field is written to an individual file using the "StandardName" attribute
!   as NetCDF attribute.
!
!   The arguments are:
!   \begin{description}
!   \item[state]
!     The {\tt ESMF\_State} object containing the fields.
!   \item[{[fieldNameList]}]
!     List of names of the fields to be written. By default write all the fields
!     in {\tt state}.
!   \item[{[fileNamePrefix]}]
!     File name prefix, common to all the files written.
!   \item[{[singleFile]}]
!     A logical flag, the default is .false., i.e., all fields in the state
!     are written in one single file. If .false., each field will be written
!     in separate files; these files are numbered with the name based on the
!     argument "file". That is, a set of files are named: [file_name]001,
!     [file_name]002, [file_name]003,...
!   \item[{[overwrite]}]
!      A logical flag, the default is .false., i.e., existing Field data may
!      {\em not} be overwritten. If .true., the
!      data corresponding to each field's name will be
!      be overwritten. If the {\tt timeslice} option is given, only data for
!      the given timeslice may be overwritten.
!      Note that it is always an error to attempt to overwrite a NetCDF
!      variable with data which has a different shape.
!   \item[{[status]}]
!      The file status. Valid options are {\tt ESMF\_FILESTATUS\_NEW},
!      {\tt ESMF\_FILESTATUS\_OLD}, {\tt ESMF\_FILESTATUS\_REPLACE}, and
!      {\tt ESMF\_FILESTATUS\_UNKNOWN} (default).
!   \item[{[timeslice]}]
!     Time slice counter. Must be positive. The behavior of this
!     option may depend on the setting of the {\tt overwrite} flag:
!     \begin{description}
!     \item[{\tt overwrite = .false.}:]\ If the timeslice value is
!     less than the maximum time already in the file, the write will fail.
!     \item[{\tt overwrite = .true.}:]\ Any positive timeslice value is valid.
!     \end{description}
!     By default, i.e. by omitting the {\tt timeslice} argument, no
!     provisions for time slicing are made in the output file,
!     however, if the file already contains a time axis for the variable,
!     a timeslice one greater than the maximum will be written.
!   \item[{[relaxedflag]}]
!     If {\tt .true.}, then no error is returned even if the call cannot write
!     the file due to library limitations. Default is {\tt .false.}.
!   \item[{[rc]}]
!     Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    integer                         :: i, itemCount
    logical                         :: lsinglefile
    logical                         :: writeflag
    type(ESMF_Field)                :: field
    type(ESMF_Array)                :: array
    type(ESMF_ArrayBundle)          :: arraybundle
    type(ESMF_StateItem_Flag)       :: itemType
    character(len=80)               :: fileName
    character(len=80)               :: stateName
    character(len=80), allocatable  :: fieldNameList_loc(:)

#ifdef DEBUG
    call ESMF_LogWrite("entered",ESMF_LOGMSG_INFO,line=__LINE__,file=__FILE__)
#endif

    if (present(rc)) rc = ESMF_SUCCESS

    if (present(singlefile)) then
      lsinglefile = singlefile
    else
      lsinglefile = .FALSE.
    endif

    if (present(fieldNameList)) then
      allocate(fieldNameList_loc(size(fieldNameList)))
      do i=1, size(fieldNameList)
        fieldNameList_loc(i) = trim(fieldNameList(i))
      enddo
    else
      call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(fieldNameList_loc(itemCount))
      call ESMF_StateGet(state, itemNameList=fieldNameList_loc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    if (lsinglefile) then
      arraybundle = ESMF_ArrayBundleCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      call ESMF_StateGet(state, name=stateName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=FILENAME)) &
        return  ! bail out

      ! Reset writeflag so that empty state is not written to file.
      writeflag = .FALSE.
      do i=1, size(fieldNameList_loc)
        call ESMF_StateGet(state, itemName=fieldNameList_loc(i), &
          itemType=itemType, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=FILENAME)) &
          return  ! bail out
        if (itemType == ESMF_STATEITEM_FIELD) then
          ! field is available in the state
          writeflag = .TRUE.
          call ESMF_StateGet(state, itemName=fieldNameList_loc(i), field=field, &
            rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=FILENAME)) return  ! bail out
          call ESMF_FieldGet(field, array=array, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=FILENAME)) return  ! bail out
          call ESMF_ArraySet(array, name=trim(fieldNameList_loc(i)), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=FILENAME)) return  ! bail out
          call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=FILENAME)) return  ! bail out
        endif
      enddo

      if (writeflag) then
        ! -> output to file
        if (present(fileNamePrefix)) then
          write (fileName,"(A)") trim(fileNamePrefix)//".nc"
        else
          write (fileName,"(A)") trim(stateName)//".nc"
        endif

        call ESMF_ArrayBundleWrite(arraybundle, fileName=trim(fileName), &
          singleFile=lsingleFile, overwrite=overwrite, status=status, &
          timeslice=timeslice, iofmt=iofmt, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg="Failed writing file: "// &
          trim(fileName), &
          line=__LINE__, &
          file=FILENAME)) &
          return  ! bail out
        call ESMF_ArrayBundleDestroy(arraybundle,rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      endif
    else
      do i=1, size(fieldNameList_loc)
        call ESMF_StateGet(state, itemName=fieldNameList_loc(i), &
          itemType=itemType, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=FILENAME)) &
          return  ! bail out
        if (itemType == ESMF_STATEITEM_FIELD) then
          ! field is available in the state
          call ESMF_StateGet(state, itemName=fieldNameList_loc(i), field=field, &
            rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
          ! -> output to file
          if (present(fileNamePrefix)) then
            write (fileName,"(A)") fileNamePrefix//trim(fieldNameList_loc(i))//".nc"
          else
            write (fileName,"(A)") trim(fieldNameList_loc(i))//".nc"
          endif
          call NUOPC_Write(field, fileName=trim(fileName), &
            overwrite=overwrite, status=status, timeslice=timeslice, &
            relaxedflag=relaxedflag, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg="Failed writing file: "// &
            trim(fileName), &
            line=__LINE__, &
            file=FILENAME)) &
            return  ! bail out
        endif
      enddo
    endif

    deallocate(fieldNameList_loc)

#ifdef DEBUG
    call ESMF_LogWrite("leaving",ESMF_LOGMSG_INFO,line=__LINE__,file=__FILE__)
#endif

  end subroutine
  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: NUOPC_Write - Write Grid data to file
! !INTERFACE:
  ! call using generic interface: NUOPC_Write
  subroutine beta_NUOPC_GridWrite_coords(grid, nclScript, mapName, minCoords, &
    maxCoords, fileName, overwrite, status, timeslice, iofmt, relaxedflag, rc)
! ! ARGUMENTS
    type(ESMF_Grid),            intent(in)            :: grid
    logical,                    intent(in)            :: nclScript
    character(len=*),           intent(in)            :: mapName
    real,                       intent(in)            :: minCoords(2)
    real,                       intent(in)            :: maxCoords(2)
    character(len=*),           intent(in),  optional :: fileName
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    type(ESMF_IOFmt_Flag),      intent(in),  optional :: iofmt
    logical,                    intent(in),  optional :: relaxedflag
    integer,                    intent(out), optional :: rc
! !DESCRIPTION:
!   Write the data in {\tt grid} to {\tt file} if supported by the
!   {\tt iofmt}.
!
!   The arguments are:
!   \begin{description}
!   \item[grid]
!     The {\tt ESMF\_Grid} object whose data is to be written.
!   \item[{[nclScript]}]
!     If {\tt .true.} then NUOPC will write an NCL script that can be used to
!     generate grid graphics. Default is {\tt .false.}.
!   \item[{[mapName]}]
!     Map name to be written to NCL script.
!   \item[{[minCoords]}]
!     Minimum map coordinates to be written to NCL script.
!   \item[{[maxCoords]}]
!     Maximum map coordinates to be written to NCL script.
!   \item[fileName]
!     The name of the file to write to. If not present then the file will
!     be written to the grid's name.
!   \item[{[overwrite]}]
!      A logical flag, the default is .false., i.e., existing Field data may
!      {\em not} be overwritten. If .true., the
!      data corresponding to each field's name will be
!      be overwritten. If the {\tt timeslice} option is given, only data for
!      the given timeslice may be overwritten.
!      Note that it is always an error to attempt to overwrite a NetCDF
!      variable with data which has a different shape.
!   \item[{[status]}]
!      The file status. Valid options are {\tt ESMF\_FILESTATUS\_NEW},
!      {\tt ESMF\_FILESTATUS\_OLD}, {\tt ESMF\_FILESTATUS\_REPLACE}, and
!      {\tt ESMF\_FILESTATUS\_UNKNOWN} (default).
!   \item[{[timeslice]}]
!     Time slice counter. Must be positive. The behavior of this
!     option may depend on the setting of the {\tt overwrite} flag:
!     \begin{description}
!     \item[{\tt overwrite = .false.}:]\ If the timeslice value is
!     less than the maximum time already in the file, the write will fail.
!     \item[{\tt overwrite = .true.}:]\ Any positive timeslice value is valid.
!     \end{description}
!     By default, i.e. by omitting the {\tt timeslice} argument, no
!     provisions for time slicing are made in the output file,
!     however, if the file already contains a time axis for the variable,
!     a timeslice one greater than the maximum will be written.
!   \item[{[iofmt]}]
!    The IO format.  Valid options are  {\tt ESMF\_IOFMT\_BIN} and
!    {\tt ESMF\_IOFMT\_NETCDF}. If not present, file names with a {\tt .bin}
!    extension will use {\tt ESMF\_IOFMT\_BIN}, and file names with a {\tt .nc}
!    extension will use {\tt ESMF\_IOFMT\_NETCDF}.  Other files default to
!    {\tt ESMF\_IOFMT\_NETCDF}.
!   \item[{[relaxedflag]}]
!     If {\tt .true.}, then no error is returned even if the call cannot write
!     the file due to library limitations. Default is {\tt .false.}.
!   \item[{[rc]}]
!     Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    type(MapDesc)  :: map

    if (present(rc)) rc = ESMF_SUCCESS

    map = MapDesc(trim(mapName),.FALSE., &
      minCoords(2),maxCoords(2),minCoords(1),maxCoords(1))

    call beta_NUOPC_GridWrite(grid, fileName=fileName, overwrite=overwrite, &
      status=status, timeslice=timeslice, iofmt=iofmt, &
      relaxedflag=relaxedflag, nclScript=nclScript, map=map, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: NUOPC_Write - Write Grid data to file
! !INTERFACE:
  ! call using generic interface: NUOPC_Write
  subroutine beta_NUOPC_GridWrite_preset(grid, nclScript, mapPreset, &
    fileName, overwrite, status, timeslice, iofmt, relaxedflag, rc)
! ! ARGUMENTS
    type(ESMF_Grid),            intent(in)            :: grid
    logical,                    intent(in)            :: nclScript
    character(len=*),           intent(in)            :: mapPreset
    character(len=*),           intent(in),  optional :: fileName
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    type(ESMF_IOFmt_Flag),      intent(in),  optional :: iofmt
    logical,                    intent(in),  optional :: relaxedflag
    integer,                    intent(out), optional :: rc
! !DESCRIPTION:
!   Write the data in {\tt grid} to {\tt file} if supported by the
!   {\tt iofmt}.
!
!   The arguments are:
!   \begin{description}
!   \item[grid]
!     The {\tt ESMF\_Grid} object whose data is to be written.
!   \item[{[nclScript]}]
!     If {\tt .true.} then NUOPC will write an NCL script that can be used to
!     generate grid graphics. Default is {\tt .false.}.
!   \item[{[mapPreset]}]
!     Map preset to use when writting NCL script.
!   \item[fileName]
!     The name of the file to write to. If not present then the file will
!     be written to the grid's name.
!   \item[{[overwrite]}]
!      A logical flag, the default is .false., i.e., existing Field data may
!      {\em not} be overwritten. If .true., the
!      data corresponding to each field's name will be
!      be overwritten. If the {\tt timeslice} option is given, only data for
!      the given timeslice may be overwritten.
!      Note that it is always an error to attempt to overwrite a NetCDF
!      variable with data which has a different shape.
!   \item[{[status]}]
!      The file status. Valid options are {\tt ESMF\_FILESTATUS\_NEW},
!      {\tt ESMF\_FILESTATUS\_OLD}, {\tt ESMF\_FILESTATUS\_REPLACE}, and
!      {\tt ESMF\_FILESTATUS\_UNKNOWN} (default).
!   \item[{[timeslice]}]
!     Time slice counter. Must be positive. The behavior of this
!     option may depend on the setting of the {\tt overwrite} flag:
!     \begin{description}
!     \item[{\tt overwrite = .false.}:]\ If the timeslice value is
!     less than the maximum time already in the file, the write will fail.
!     \item[{\tt overwrite = .true.}:]\ Any positive timeslice value is valid.
!     \end{description}
!     By default, i.e. by omitting the {\tt timeslice} argument, no
!     provisions for time slicing are made in the output file,
!     however, if the file already contains a time axis for the variable,
!     a timeslice one greater than the maximum will be written.
!   \item[{[iofmt]}]
!    The IO format.  Valid options are  {\tt ESMF\_IOFMT\_BIN} and
!    {\tt ESMF\_IOFMT\_NETCDF}. If not present, file names with a {\tt .bin}
!    extension will use {\tt ESMF\_IOFMT\_BIN}, and file names with a {\tt .nc}
!    extension will use {\tt ESMF\_IOFMT\_NETCDF}.  Other files default to
!    {\tt ESMF\_IOFMT\_NETCDF}.
!   \item[{[relaxedflag]}]
!     If {\tt .true.}, then no error is returned even if the call cannot write
!     the file due to library limitations. Default is {\tt .false.}.
!   \item[{[rc]}]
!     Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! no local variables

    if (present(rc)) rc = ESMF_SUCCESS

    select case (trim(mapPreset))
      case ('global','GLOBAL','Global')
        call beta_NUOPC_GridWrite(grid, fileName=fileName, overwrite=overwrite, &
          timeslice=timeslice, iofmt=iofmt, relaxedflag=relaxedflag, &
          nclScript=nclScript, map=NUOPC_MAPPRESET_GLOBAL, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out  ! bail out
      case ('conus','CONUS','Conus')
        call beta_NUOPC_GridWrite(grid, fileName=fileName, overwrite=overwrite, &
          timeslice=timeslice, iofmt=iofmt, relaxedflag=relaxedflag, &
          nclScript=nclScript, map=NUOPC_MAPPRESET_CONUS, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out  ! bail out
      case ('irene','IRENE','Irene')
        call beta_NUOPC_GridWrite(grid, fileName=fileName, overwrite=overwrite, &
          timeslice=timeslice, iofmt=iofmt, relaxedflag=relaxedflag, &
          nclScript=nclScript, map=NUOPC_MAPPRESET_IRENE, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out  ! bail out
      case ('frontrange','FRONTRANGE','FrontRange')
        call beta_NUOPC_GridWrite(grid, fileName=fileName, overwrite=overwrite, &
          timeslice=timeslice, iofmt=iofmt, relaxedflag=relaxedflag, &
          nclScript=nclScript, map=NUOPC_MAPPRESET_FRONTRANGE, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out  ! bail out
      case default
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_VALUE,   &
          msg="Unknown map preset value "//trim(mapPreset)//".", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
    endselect

  end subroutine  

  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: NUOPC_Write - Write Grid data to file
! !INTERFACE:
  ! call using generic interface: NUOPC_Write
  subroutine beta_NUOPC_GridWrite_default(grid, fileName, overwrite, status, &
    timeslice, iofmt, relaxedflag, nclScript, map, rc)
! ! ARGUMENTS
    type(ESMF_Grid),            intent(in)            :: grid
    character(len=*),           intent(in),  optional :: fileName
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    type(ESMF_IOFmt_Flag),      intent(in),  optional :: iofmt
    logical,                    intent(in),  optional :: relaxedflag 
    logical,                    intent(in),  optional :: nclScript
    type(MapDesc),              intent(in),  optional :: map
    integer,                    intent(out), optional :: rc
! !DESCRIPTION:
!   Write the data in {\tt grid} to {\tt file} if supported by the
!   {\tt iofmt}.
!
!   The arguments are:
!   \begin{description}
!   \item[field]
!     The {\tt ESMF\_Field} object whose data is to be written.
!   \item[fileName]
!     The name of the file to write to. If not present then the file will
!     be written to the grid's name.
!   \item[{[overwrite]}]
!      A logical flag, the default is .false., i.e., existing Field data may
!      {\em not} be overwritten. If .true., the
!      data corresponding to each field's name will be
!      be overwritten. If the {\tt timeslice} option is given, only data for
!      the given timeslice may be overwritten.
!      Note that it is always an error to attempt to overwrite a NetCDF
!      variable with data which has a different shape.
!   \item[{[status]}]
!      The file status. Valid options are {\tt ESMF\_FILESTATUS\_NEW},
!      {\tt ESMF\_FILESTATUS\_OLD}, {\tt ESMF\_FILESTATUS\_REPLACE}, and
!      {\tt ESMF\_FILESTATUS\_UNKNOWN} (default).
!   \item[{[timeslice]}]
!     Time slice counter. Must be positive. The behavior of this
!     option may depend on the setting of the {\tt overwrite} flag:
!     \begin{description}
!     \item[{\tt overwrite = .false.}:]\ If the timeslice value is
!     less than the maximum time already in the file, the write will fail.
!     \item[{\tt overwrite = .true.}:]\ Any positive timeslice value is valid.
!     \end{description}
!     By default, i.e. by omitting the {\tt timeslice} argument, no
!     provisions for time slicing are made in the output file,
!     however, if the file already contains a time axis for the variable,
!     a timeslice one greater than the maximum will be written.
!   \item[{[iofmt]}]
!    The IO format.  Valid options are  {\tt ESMF\_IOFMT\_BIN} and
!    {\tt ESMF\_IOFMT\_NETCDF}. If not present, file names with a {\tt .bin}
!    extension will use {\tt ESMF\_IOFMT\_BIN}, and file names with a {\tt .nc}
!    extension will use {\tt ESMF\_IOFMT\_NETCDF}.  Other files default to
!    {\tt ESMF\_IOFMT\_NETCDF}.
!   \item[{[relaxedflag]}]
!     If {\tt .true.}, then no error is returned even if the call cannot write
!     the file due to library limitations. Default is {\tt .false.}.
!   \item[{[nclScript]}]
!     If {\tt .true.} then NUOPC will write an NCL script that can be used to 
!     generate grid graphics. Default is {\tt .false.}.
!   \item[{[map]}]
!     Derived type including the map name and boundary coordinates.
!   \item[{[rc]}]
!     Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
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
    logical                 :: lnclScript
    logical                 :: hasCorners

    if (present(rc)) rc = ESMF_SUCCESS

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
          line=__LINE__, file=FILENAME)) return  ! bail out  ! bail out
        lfileName = trim(gridName)//".nc"
      endif

      arraybundle = ESMF_ArrayBundleCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out

      ! -- centers --

      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (isPresent) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_ArraySet(array, name="lon_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_ArraySet(array, name="lat_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      endif

      ! -- corners --

      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        isPresent=hasCorners, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (hasCorners) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) then
          call ESMF_ArraySet(array, name="lon_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=FILENAME)) return  ! bail out
          call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=FILENAME)) return  ! bail out
        endif
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) then
          call ESMF_ArraySet(array, name="lat_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=FILENAME)) return  ! bail out
          call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=FILENAME)) return  ! bail out
        endif
      endif

      ! -- mask --

      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_ArraySet(array, name="mask", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      endif

      ! -- area --

      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_ArraySet(array, name="area", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      endif

      call ESMF_ArrayBundleWrite(arraybundle, &
        fileName=trim(lfileName),rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out

      call ESMF_ArrayBundleDestroy(arraybundle,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out

      if (present(nclScript)) then
        lnclScript = nclScript
      else
        lnclScript = .FALSE.
      endif

      if (lnclScript) then
        call ESMF_GridGet(grid,dimCount=dimCount,rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out  ! bail out

        ! allocate coordDim info accord. to dimCount and tileCount
        allocate(coordDimCount(dimCount), &
          stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Allocation of coordinate dimensions memory failed.", &
          line=__LINE__, file=FILENAME)) &
          return  ! bail out

        ! get coordDim info
        call ESMF_GridGet(grid, coordDimCount=coordDimCount, &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out  ! bail out

        coordDimMax = 0
        do dimIndex=1,dimCount
          coordDimMax = MAX(coordDimMax,coordDimCount(dimIndex))
        enddo

        ! deallocate coordDim info accord. to dimCount and tileCount
        deallocate(coordDimCount, stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Deallocation of coordinate dimensions memory failed.", &
          line=__LINE__, file=FILENAME)) &
          return  ! bail out

        if (coordDimMax == 1) then
          call beta_NUOPC_NclScriptWrite(gridFile=lfileName, map=map, &
            uniformRect=.TRUE., writeCorners=hasCorners, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=FILENAME)) return  ! bail out
        else
          call beta_NUOPC_NclScriptWrite(gridFile=lfileName, map=map, &
            writeCorners=hasCorners, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=FILENAME)) return  ! bail out
        endif
      endif
    endif
  end subroutine

  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: beta_NUOPC_NclScriptWrite_coords - Write NCL Script to file
! !INTERFACE:
  ! call using generic interface: beta_NUOPC_NclScriptWrite
  subroutine beta_NUOPC_NclScriptWrite_coords(gridFile,mapName,minCoords,maxCoords, &
  title,nclFile,uniformRect,writeCorners,rc)
! ! ARGUMENTS
    character(len=*),intent(in)          :: gridFile
    character(len=*),intent(in)          :: mapName         
    real,intent(in)                      :: minCoords(2)
    real,intent(in)                      :: maxCoords(2)
    character(len=*),intent(in),optional :: title
    character(len=*),intent(in),optional :: nclFile
    logical,intent(in),optional          :: uniformRect
    logical,intent(in),optional          :: writeCorners
    integer,intent(out),optional         :: rc
! !DESCRIPTION:
!   Write NCL script that can be used to generate NCL grid visual.
!
!   The arguments are:
!   \begin{description}
!   \item[gridFile]
!     NetCDF grid file name used plot data on coordinates.
!   \item[mapName]
!     Map name.
!   \item[minCoords]
!     Minimum coordinate limits for plot.
!   \item[maxCoords]
!     Maximum coordinate limits for plot.
!   \item[title]
!     Grid plot title.
!   \item[nclFile]
!     NCL script file name
!   \item[uniformRect]
!     Repeat coordinates for uniform rectangular grids.
!   \item[writeCorners]
!     Plot corners.  Default plot center coordinates.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    type(MapDesc)  :: map

    if (present(rc)) rc = ESMF_SUCCESS

    map = MapDesc(trim(mapName),.FALSE., &
      minCoords(2),maxCoords(2),minCoords(1),maxCoords(1))

    call beta_NUOPC_NclScriptWrite(gridFile, map=map, title=title, &
      nclFile=nclFile, uniformRect=uniformRect, writeCorners=writeCorners, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: beta_NUOPC_NclScriptWrite_preset - Write NCL Script to file
! !INTERFACE:
  ! call using generic interface: beta_NUOPC_NclScriptWrite
  subroutine beta_NUOPC_NclScriptWrite_preset(gridFile,mapPreset, &
  title,nclFile,uniformRect,writeCorners,rc)
! ! ARGUMENTS
    character(len=*),intent(in)          :: gridFile
    character(len=*),intent(in)          :: mapPreset
    character(len=*),intent(in),optional :: title
    character(len=*),intent(in),optional :: nclFile
    logical,intent(in),optional          :: uniformRect
    logical,intent(in),optional          :: writeCorners
    integer,intent(out),optional         :: rc
! !DESCRIPTION:
!   Write NCL script that can be used to generate NCL grid visual.
!
!   The arguments are:
!   \begin{description}
!   \item[gridFile]
!     NetCDF grid file name used plot data on coordinates.
!   \item[mapPreset]
!     Preset coordinate limits.
!   \item[title]
!     Grid plot title.
!   \item[nclFile]
!     NCL script file name
!   \item[uniformRect]
!     Repeat coordinates for uniform rectangular grids.
!   \item[writeCorners]
!     Plot corners.  Default plot center coordinates.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables

    if (present(rc)) rc = ESMF_SUCCESS

    select case (trim(mapPreset))
      case ('global','GLOBAL','Global')
        call beta_NUOPC_NclScriptWrite(gridFile, &
          map=NUOPC_MAPPRESET_GLOBAL, title=title, nclFile=nclFile, &
          uniformRect=uniformRect, writeCorners=writeCorners, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      case ('conus','CONUS','Conus')
        call beta_NUOPC_NclScriptWrite(gridFile, &
          map=NUOPC_MAPPRESET_CONUS, title=title, nclFile=nclFile, &
          uniformRect=uniformRect, writeCorners=writeCorners, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      case ('irene','IRENE','Irene')
        call beta_NUOPC_NclScriptWrite(gridFile, &
          map=NUOPC_MAPPRESET_IRENE, title=title, nclFile=nclFile, &
          uniformRect=uniformRect, writeCorners=writeCorners, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      case ('frontrange','FRONTRANGE','FrontRange')
        call beta_NUOPC_NclScriptWrite(gridFile, &
          map=NUOPC_MAPPRESET_FRONTRANGE, title=title, nclFile=nclFile, &
          uniformRect=uniformRect, writeCorners=writeCorners, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      case default
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_VALUE,   &
          msg="Unknown map preset value "//trim(mapPreset)//".", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
    endselect

  end subroutine

  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: beta_NUOPC_NclScriptWrite_default - Write NCL Script to file
! !INTERFACE:
  ! call using generic interface: NUOPC_Write
  subroutine beta_NUOPC_NclScriptWrite_default(gridFile, map, title, &
    nclFile, uniformRect, writeCorners, rc)
! ! ARGUMENTS
    character(len=*),intent(in)          :: gridFile
    type(MapDesc),intent(in)             :: map
    character(len=*),intent(in),optional :: title
    character(len=*),intent(in),optional :: nclFile
    logical,intent(in),optional          :: uniformRect
    logical,intent(in),optional          :: writeCorners
    integer,intent(out),optional         :: rc
! !DESCRIPTION:
!   Write NCL script that can be used to generate NCL grid visual.
!
!   The arguments are:
!   \begin{description}
!   \item[gridFile]
!     NetCDF grid file name used plot data on coordinates.
!   \item[map]
!     Coordinate limits.
!   \item[title]
!     Grid plot title.
!   \item[nclFile]
!     NCL script file name
!   \item[uniformRect]
!     Repeat coordinates for uniform rectangular grids.
!   \item[writeCorners]
!     Plot corners.  Default plot center coordinates.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    type(ESMF_VM)                   :: vm
    integer                         :: lpe
    character(len=64)               :: ltitle
    character(len=64)               :: lnclFile
    character(len=64)               :: fileBase
    logical                         :: luniformRect
    logical                         :: lcorners
    integer                         :: markExt
    integer                         :: fUnit
    integer                         :: stat
    character(len=10)               :: varlat
    character(len=10)               :: varlon

    if (present(rc)) rc = ESMF_SUCCESS

    ! Get current VM and pet number
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    call ESMF_VMGet(vm, localPet=lpe, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    if (lpe /= 0) return

    markExt = index(gridFile,".",back=.TRUE.)
    if (markExt > 1) then
      fileBase = gridFile(1:markExt-1)
    elseif (markExt == 1) then
      fileBase = "grid"
    elseif (len(trim(gridFile)) > 0) then
      fileBase = trim(gridFile)
    else
      fileBase = "grid"
    endif

    if (present(nclFile)) then
      lnclFile = trim(nclFile)
    else
      lnclFile = trim(fileBase)//".ncl"
    endif

    if (present(title)) then
      ltitle = trim(title)
    else
      ltitle = trim(fileBase)//" "//trim(map%name)
    endif

    if (present(uniformRect)) then
      luniformRect = uniformRect
    else
      luniformRect = .FALSE.
    endif

    if (present(writeCorners)) then
      lcorners = writeCorners
    else
      lcorners = .FALSE.
    endif

    if (lcorners) then
      varlat = 'lat_corner'
      varlon = 'lon_corner'
    else
      varlat = 'lat_center'
      varlon = 'lon_center'
    endif

    call ESMF_UtilIOUnitGet(fUnit, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out
    open (fUnit,file=trim(lnclFile),action="write", &
      status="new",iostat=stat)
    if (stat /= 0) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_OPEN,   &
        msg="Cound not open "//trim(lnclFile)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    write (fUnit,"(A)") '; NUOPC_FileWriteMapNCL used to generate this file'
    write (fUnit,"(A)") '; execute ncl <this file> to generate grid png file'
    write (fUnit,"(A)") 'begin'
    write (fUnit,"(A)") '  in = addfile("'//trim(gridFile)//'","r")'
    write (fUnit,"(A)") '  wks = gsn_open_wks("png","'//trim(fileBase)//'")'
    write (fUnit,"(A)") '  res              = True'
    write (fUnit,"(A)") '  res@gsnMaximize  = True'
    write (fUnit,"(A)") '  res@gsnDraw      = False'
    write (fUnit,"(A)") '  res@gsnFrame     = False'
    write (fUnit,"(A)") '  res@tiMainString = "'//trim(ltitle)//'"'
    write (fUnit,"(A)") '  res@pmTickMarkDisplayMode = "Always"'
    write (fUnit,"(A)") '; '//trim(map%name)//' Map Grid'
    if (.NOT. map%global) then
      write (fUnit,"(A,F0.3)") &
        '  res@mpMinLatF    = ',map%minLat
      write (fUnit,"(A,F0.3)") &
        '  res@mpMaxLatF    = ',map%maxLat
      write (fUnit,"(A,F0.3)") &
        '  res@mpMinLonF    = ',map%minLon
      write (fUnit,"(A,F0.3)") &
        '  res@mpMaxLonF    = ',map%maxLon
    endif
    write (fUnit,"(A)") '  map = gsn_csm_map_ce(wks,res)'
    write (fUnit,"(A)") '  hgt = in->lon_center(:,:)'
    write (fUnit,"(A)") '  dimlon = getfilevardimsizes(in,"'//trim(varlon)//'")'
    write (fUnit,"(A)") '  dimlat = getfilevardimsizes(in,"'//trim(varlat)//'")'
    if (luniformRect) then
      write (fUnit,"(A)") '  hgt@lat2d = conform_dims((/dimlon(0),dimlat(1)/),'// &
        'in->'//trim(varlat)//'(:,0),1)'
      write (fUnit,"(A)") '  hgt@lon2d = conform_dims((/dimlon(0),dimlat(1)/),'// &
        'in->'//trim(varlon)//'(0,:),0)'
    else
      write (fUnit,"(A)") '  hgt@lat2d = in->'//trim(varlat)//'(:,:)'
      write (fUnit,"(A)") '  hgt@lon2d = in->'//trim(varlon)//'(:,:)'
    endif
    write (fUnit,"(A)") '  pres                   = True'
    if (lcorners) then
      write (fUnit,"(A)") '  pres@gsnCoordsAsLines  = True'
    else
      write (fUnit,"(A)") '  pres@gsnCoordsAsLines  = False'
      write (fUnit,"(A)") '  if (dimlon(0)*dimlat(1) .gt. 99) then'
      write (fUnit,"(A)") '    pres@gsMarkerIndex = 1'
      write (fUnit,"(A)") '  end if'
    endif
    write (fUnit,"(A)") '  gsn_coordinates(wks,map,hgt,pres)'
    write (fUnit,"(A)") 'end'

    close (fUnit,iostat=stat)
    if (stat /= 0) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_CLOSE,   &
        msg="Cound not close "//trim(lnclFile)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine    

  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: beta_NUOPC_FerretScriptWrite_coords - Write Ferret Script to file
! !INTERFACE:
  ! call using generic interface: beta_NUOPC_FerretScriptWrite
  subroutine beta_NUOPC_FerretScriptWrite_coords(varName, dataFile, &
    gridFile, slices, mapName, minCoords, maxCoords, scale, jnlFile, &
    uniformRect,rc)
! ! ARGUMENTS
    character(len=*),intent(in)          :: varName
    character(len=*),intent(in)          :: dataFile
    character(len=*),intent(in)          :: gridFile
    integer,intent(in)                   :: slices(:)
    character(len=*),intent(in)          :: mapName
    real,intent(in)                      :: minCoords(2)
    real,intent(in)                      :: maxCoords(2)
    real,intent(in),optional             :: scale(3)
    character(len=*),intent(in),optional :: jnlFile
    logical,intent(in),optional          :: uniformRect
    integer,intent(out),optional         :: rc
! !DESCRIPTION:
!   Write Ferret script that can be used to generate Ferret field visual.
!
!   The arguments are:
!   \begin{description}
!   \item[varName]
!     NetCDF variable name in dataFile.
!   \item[dataFile]
!     NetCDF field data file name.
!   \item[gridFile]
!     NetCDF grid file name used plot data on coordinates.
!   \item[slices]
!     Array of slice numbers to for which to create plots.
!   \item[mapName]
!     Map name used to create map description.
!   \item[minCoords]
!     Minimum coordinates used to plot data.
!   \item[maxCoords]
!     Maximum coordinates used to plot data.
!   \item[scale]
!     Field data scale (/ minimum value, maximum value, step size /)
!   \item[jnlFile]
!     Output script filename.
!   \item[uniformRect]
!     Repeat coordinates for uniform rectangular grids.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    type(MapDesc)  :: map

    if (present(rc)) rc = ESMF_SUCCESS

    map = MapDesc(trim(mapName),.FALSE., &
      minCoords(2),maxCoords(2),minCoords(1),maxCoords(1))

    call beta_NUOPC_FerretScriptWrite(varName,dataFile,gridFile,slices, &
      map=map,scale=scale,jnlFile=jnlFile, uniformRect=uniformRect,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: beta_NUOPC_FerretScriptWrite_preset - Write Ferret Script to file
! !INTERFACE:
  ! call using generic interface: beta_NUOPC_FerretScriptWrite
  subroutine beta_NUOPC_FerretScriptWrite_preset(varName, dataFile, &
    gridFile, slices, mapPreset, scale, jnlFile, &
    uniformRect,rc)
! ! ARGUMENTS
    character(len=*),intent(in)          :: varName
    character(len=*),intent(in)          :: dataFile
    character(len=*),intent(in)          :: gridFile
    integer,intent(in)                   :: slices(:)
    character(len=*),intent(in)          :: mapPreset
    real,intent(in),optional             :: scale(3)
    character(len=*),intent(in),optional :: jnlFile
    logical,intent(in),optional          :: uniformRect
    integer,intent(out),optional         :: rc
! !DESCRIPTION:
!   Write Ferret script that can be used to generate Ferret field visual.
!
!   The arguments are:
!   \begin{description}
!   \item[varName]
!     NetCDF variable name in dataFile.
!   \item[dataFile]
!     NetCDF field data file name.
!   \item[gridFile]
!     NetCDF grid file name used plot data on coordinates.
!   \item[slices]
!     Array of slice numbers to for which to create plots.
!   \item[mapPreset]
!     Map preset to use to define coordinate limits.
!   \item[scale]
!     Field data scale (/ minimum value, maximum value, step size /)
!   \item[jnlFile]
!     Output script filename.
!   \item[uniformRect]
!     Repeat coordinates for uniform rectangular grids.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables

    if (present(rc)) rc = ESMF_SUCCESS

    select case (trim(mapPreset))
      case ('global','GLOBAL','Global')
        call beta_NUOPC_FerretScriptWrite(varName, dataFile, gridFile, slices, &
          map=NUOPC_MAPPRESET_GLOBAL, scale=scale, jnlFile=jnlFile, &
          uniformRect=uniformRect, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      case ('conus','CONUS','Conus')
        call beta_NUOPC_FerretScriptWrite(varName, dataFile, gridFile, slices, &
          map=NUOPC_MAPPRESET_CONUS, scale=scale, jnlFile=jnlFile, &
          uniformRect=uniformRect, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      case ('irene','IRENE','Irene')
        call beta_NUOPC_FerretScriptWrite(varName, dataFile, gridFile, slices, &
          map=NUOPC_MAPPRESET_IRENE, scale=scale, jnlFile=jnlFile, &
          uniformRect=uniformRect, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      case ('frontrange','FRONTRANGE','FrontRange')
        call beta_NUOPC_FerretScriptWrite(varName, dataFile, gridFile, slices, &
          map=NUOPC_MAPPRESET_FRONTRANGE, scale=scale, jnlFile=jnlFile, &
          uniformRect=uniformRect, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      case default
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_VALUE,   &
          msg="Unknown map preset value "//trim(mapPreset)//".", &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return
    endselect

  end subroutine

  !-----------------------------------------------------------------------------
!BOP
! !IROUTINE: beta_NUOPC_FerretScriptWrite_default - Write Ferret Script to file
! !INTERFACE:
  ! call using generic interface: beta_NUOPC_FerretScriptWrite
  subroutine beta_NUOPC_FerretScriptWrite_default(varName, dataFile, gridFile, &
    slices, map, scale, jnlFile, uniformRect, rc)
! ! ARGUMENTS
    character(len=*),intent(in)          :: varName
    character(len=*),intent(in)          :: dataFile
    character(len=*),intent(in)          :: gridFile
    integer,intent(in)                   :: slices(:)
    type(MapDesc),intent(in)             :: map
    real,intent(in),optional             :: scale(3)
    character(len=*),intent(in),optional :: jnlFile
    logical,intent(in),optional          :: uniformRect
    integer,intent(out),optional         :: rc
! !DESCRIPTION:
!   Write Ferret script that can be used to generate Ferret field visual.
!
!   The arguments are:
!   \begin{description}
!   \item[varName]
!     NetCDF variable name in dataFile.
!   \item[dataFile]
!     NetCDF field data file name.
!   \item[gridFile]
!     NetCDF grid file name used plot data on coordinates.
!   \item[slices]
!     Array of slice numbers to for which to create plots.
!   \item[map]
!     Coordinate limits.
!   \item[scale]
!     Field data scale (/ minimum value, maximum value, step size /)
!   \item[jnlFile]
!     Output script filename.
!   \item[uniformRect]
!     Repeat coordinates for uniform rectangular grids.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    type(ESMF_VM)                   :: vm
    integer                         :: lpe
    character(len=64)               :: ljnlFile
    logical                         :: luniformRect
    character(len=64)               :: fileBase
    character(len=64)               :: limits
    character(len=64)               :: levels
    character(len=64)               :: latlon
    integer                         :: markExt
    integer                         :: sIndex
    integer                         :: fUnit
    integer                         :: stat

    if (present(rc)) rc = ESMF_SUCCESS

    ! Get current VM and pet number
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    call ESMF_VMGet(vm, localPet=lpe, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    if (lpe /= 0) return

    markExt = index(dataFile,".",back=.TRUE.)
    if (markExt > 1) then
      fileBase = dataFile(1:markExt-1)
    elseif (markExt == 1) then
      fileBase = "data"
    elseif (len(trim(dataFile)) > 0) then
      fileBase = trim(dataFile)
    else
      fileBase = "data"
    endif

    if (present(jnlFile)) then
      ljnlFile = trim(jnlFile)
    else
      ljnlFile = trim(fileBase)//".jnl"
    endif

    if (present(uniformRect)) then
      luniformRect = uniformRect
    else
      luniformRect = .FALSE.
    endif

    if (map%global) then
      limits = ''
    else
      write (limits,"(2(A,F0.3,A,F0.3))") &
        '/vlimits=',map%minLat,':',map%maxLat, &
        '/hlimits=',map%minLon,':',map%maxLon
    endif

    if (luniformRect) then
      latlon = ',lon_center[d=2,j=1:1],lat_center[d=2,i=1:1]'
    else
      latlon = ',lon_center[d=2],lat_center[d=2]'
    endif

    if (present(scale)) then
      write (levels,"(A,F0.3,A,F0.3,A,F0.3,A)") &
        '/levels=(',scale(1),',',scale(2),',',scale(3),')'
    else
      levels = ''
    endif

    call ESMF_UtilIOUnitGet(fUnit, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out
    open (fUnit,file=trim(ljnlFile),action="write", &
      status="new",iostat=stat)
    if (stat /= 0) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_OPEN,   &
        msg="Cound not open "//trim(ljnlFile)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

    write (fUnit,"(A)") ' ! NUOPC_FileWriteJNL used to generate this file'
    write (fUnit,"(A)") ' ! execute ferret -script <this file> to '
    write (fUnit,"(A)") ' ! generate data gif file for each slice listed'
    write (fUnit,"(A)") ''
    write (fUnit,"(A,A)") 'use ',trim(dataFile)
    write (fUnit,"(A,A)") 'use ',trim(gridFile)
    write (fUnit,"(A,A)") ''

    do sIndex=1,size(slices)
      write (fUnit,"((A,I0),A,A,(A,A,A),A)") &
        'shade/k=',slices(sIndex), &
        trim(limits), &
        trim(levels), &
        ' ',trim(varName),'[d=1]', &
        trim(latlon)
      write (fUnit,"(A,(A,A,I0,A))") 'FRAME/FILE=', &
        trim(fileBase),'_',slices(sIndex),'.gif'
      write (fUnit,"(A)") ''
    enddo

    write (fUnit,"(A)") 'exit'

    close (fUnit,iostat=stat)
    if (stat /= 0) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_CLOSE,   &
        msg="Cound not close "//trim(ljnlFile)//".", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return
    endif

  end subroutine

end module
