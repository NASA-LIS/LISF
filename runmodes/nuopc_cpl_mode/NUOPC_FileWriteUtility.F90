#define ESMF_STDERRORCHECK(rc) ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)
#define FILENAME "NUOPC_FileWriteUtility.F90"
#define MODNAME "NUOPC_FileWriteUtility"

module NUOPC_FileWriteUtility
  use ESMF
  use NUOPC

  implicit none

  private

  public :: NUOPC_FileWriteGrid
  public :: NUOPC_FileWriteMapNCL
  public :: NUOPC_FileWriteJNL
  public :: NUOPC_MAPPRESET_GLOBAL
  public :: NUOPC_MAPPRESET_CONUS
  public :: NUOPC_MAPPRESET_IRENE
  public :: NUOPC_MAPPRESET_FRONTRANGE

  interface NUOPC_FileWriteGrid
    module procedure NUOPC_FileWriteGrid_coords
    module procedure NUOPC_FileWriteGrid_map
  end interface

  interface NUOPC_FileWriteMapNCL
    module procedure NUOPC_FileWriteMapNCL_coords
    module procedure NUOPC_FileWriteMapNCL_map
  end interface

  interface NUOPC_FileWriteJNL
    module procedure NUOPC_FileWriteJNL_coords
    module procedure NUOPC_FileWriteJNL_map
  end interface

  type MapDesc
    character(len=16) :: name
    logical           :: global
    real              :: minLat
    real              :: maxLat
    real              :: minLon
    real              :: maxLon
  endtype MapDesc

  type(MapDesc),parameter :: &
    NUOPC_MAPPRESET_GLOBAL = MapDesc("GLOBAL",.TRUE.,-90.0,90.0,-180.0,-180.0), &
    NUOPC_MAPPRESET_CONUS = MapDesc("CONUS",.FALSE.,18.0,49.0,235.0,298.0), &
    NUOPC_MAPPRESET_IRENE = MapDesc("IRENE",.FALSE.,10.0,60.0,240.0,320.0), &
    NUOPC_MAPPRESET_FRONTRANGE = MapDesc("FRONTRANGE",.FALSE.,38.5,41.0,-107.0,-103.5)

contains

  !-----------------------------------------------------------------------------
  ! Utilities
  !-----------------------------------------------------------------------------

  subroutine NUOPC_FileWriteGrid_coords(grid,filename,nclMapName, &
  nclMinCoords,nclMaxCoords,rc)
    ! ARGUMENTS
    type(ESMF_Grid),intent(in)           :: grid
    character(len=*),intent(in),optional :: fileName
    character(len=*),intent(in)          :: nclMapName
    real,intent(in)                      :: nclMinCoords(2)
    real,intent(in)                      :: nclMaxCoords(2)
    integer,intent(out),optional         :: rc

    ! LOCAL VARIABLES
    type(MapDesc)  :: map

    if (present(rc)) rc = ESMF_SUCCESS

    map = MapDesc(trim(nclMapName),.FALSE., &
      nclMinCoords(2),nclMaxCoords(2),nclMinCoords(1),nclMaxCoords(1))

    call NUOPC_FileWriteGrid(grid,fileName=fileName,nclMap=map,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

  end subroutine

  subroutine NUOPC_FileWriteGrid_map(grid,fileName,nclMap,rc)
    ! ARGUMENTS
    type(ESMF_Grid),intent(in)           :: grid
    character(len=*),intent(in),optional :: fileName
    type(MapDesc),intent(in),optional    :: nclMap 
    integer,intent(out),optional         :: rc

    ! LOCAL VARIABLES
    character(len=64)               :: lfileName
    character(len=64)               :: gridName
    type(ESMF_Array)                :: array
    type(ESMF_ArrayBundle)          :: arrayBundle
    logical                         :: isPresent
    integer                         :: dimCount
    integer                         :: dimIndex
    integer,allocatable             :: coordDimCount(:)
    integer                         :: coordDimMax
    integer                         :: stat
    logical                         :: corners

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_GridGet(grid, name=gridName, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (present(fileName)) then
      lfileName = trim(fileName)
    else
      lfileName = trim(gridName)//".nc"
    endif

    arrayBundle = ESMF_ArrayBundleCreate(rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! -- centers --

    call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      isPresent=isPresent, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    if (isPresent) then
      call ESMF_GridGetCoord(grid, coordDim=1, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_ArraySet(array, name="lon_center", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_ArrayBundleAdd(arrayBundle,(/array/),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_GridGetCoord(grid, coordDim=2, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_ArraySet(array, name="lat_center", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_ArrayBundleAdd(arrayBundle,(/array/),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    endif

    ! -- corners --

    call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
      isPresent=isPresent, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    if (isPresent) then
      corners = .TRUE.
      call ESMF_GridGetCoord(grid, coordDim=1, &
        staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
      if (.not. ESMF_STDERRORCHECK(rc)) then
        call ESMF_ArraySet(array, name="lon_corner", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call ESMF_ArrayBundleAdd(arrayBundle,(/array/),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
      call ESMF_GridGetCoord(grid, coordDim=2, &
        staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
      if (.not. ESMF_STDERRORCHECK(rc)) then
        call ESMF_ArraySet(array, name="lat_corner", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call ESMF_ArrayBundleAdd(arrayBundle,(/array/),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
    else
      corners = .FALSE.
    endif

    ! -- mask --

    call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    if (isPresent) then
      call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_ArraySet(array, name="mask", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_ArrayBundleAdd(arrayBundle,(/array/),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    endif

    ! -- area --

    call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    if (isPresent) then
      call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_ArraySet(array, name="area", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_ArrayBundleAdd(arrayBundle,(/array/),rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
   endif

    call ESMF_ArrayBundleWrite(arrayBundle, &
     fileName=trim(lfileName),rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_ArrayBundleDestroy(arrayBundle,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (present(nclMap)) then
      call ESMF_GridGet(grid,dimCount=dimCount,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

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
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out

      coordDimMax = 0
      do dimIndex=1,dimCount
        coordDimMax = MAX(coordDimMax,coordDimCount(dimIndex))
      enddo

      if (coordDimMax == 1) then
        call NUOPC_FileWriteMapNCL(trim(lfileName),nclMap, &
          repeatCoord=.TRUE.,corners=corners,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      else
        call NUOPC_FileWriteMapNCL(trim(lfileName),nclMap, &
          corners=corners,rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
   endif

  end subroutine

  subroutine NUOPC_FileWriteMapNCL_coords(gridFile,mapName,minCoords,maxCoords, &
  title,nclFile,repeatCoord,corners,rc)
    ! ARGUMENTS
    character(len=*),intent(in)          :: gridFile
    character(len=*),intent(in)          :: mapName          
    real,intent(in)                      :: minCoords(2)
    real,intent(in)                      :: maxCoords(2)
    character(len=*),intent(in),optional :: title
    character(len=*),intent(in),optional :: nclFile
    logical,intent(in),optional          :: repeatCoord
    logical,intent(in),optional          :: corners
    integer,intent(out),optional         :: rc

    ! LOCAL VARIABLES
    type(MapDesc)  :: map

    if (present(rc)) rc = ESMF_SUCCESS

    map = MapDesc(trim(mapName),.FALSE., &
      minCoords(2),maxCoords(2),minCoords(1),maxCoords(1))

    call NUOPC_FileWriteMapNCL(gridFile,map=map,title=title, &
      nclFile=nclFile, repeatCoord=repeatCoord, corners=corners, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  subroutine NUOPC_FileWriteMapNCL_map(gridFile,map,title, &
  nclFile,repeatCoord,corners,rc)
    ! ARGUMENTS
    character(len=*),intent(in)          :: gridFile
    type(MapDesc),intent(in)             :: map
    character(len=*),intent(in),optional :: title
    character(len=*),intent(in),optional :: nclFile
    logical,intent(in),optional          :: repeatCoord
    logical,intent(in),optional          :: corners
    integer,intent(out),optional         :: rc

    ! LOCAL VARIABLES
    type(ESMF_VM)                   :: vm
    integer                         :: lpe
    character(len=64)               :: ltitle
    character(len=64)               :: lnclFile
    character(len=64)               :: fileBase
    logical                         :: lrepeatCoord
    logical                         :: lcorners
    integer                         :: markExt
    integer                         :: fUnit
    integer                         :: stat
    character(len=10)               :: varlat
    character(len=10)               :: varlon

    if (present(rc)) rc = ESMF_SUCCESS

    ! Get current VM and pet number
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMGet(vm, localPet=lpe, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

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

    if (present(repeatCoord)) then
      lrepeatCoord = repeatCoord
    else
      lrepeatCoord = .FALSE.
    endif

    if (present(corners)) then
      lcorners = corners
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
    if (ESMF_STDERRORCHECK(rc)) return
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
    if (lrepeatCoord) then
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

  subroutine NUOPC_FileWriteJNL_coords(varName,dataFile,gridFile,slices,mapName, &
  minCoords,maxCoords,scale,jnlFile,repeatCoord,rc)
    ! ARGUMENTS
    character(len=*),intent(in)          :: varName
    character(len=*),intent(in)          :: dataFile
    character(len=*),intent(in)          :: gridFile
    integer,intent(in)                   :: slices(:)
    character(len=*),intent(in)          :: mapName
    real,intent(in)                      :: minCoords(2)
    real,intent(in)                      :: maxCoords(2)
    real,intent(in),optional             :: scale(3)
    character(len=*),intent(in),optional :: jnlFile
    logical,intent(in),optional          :: repeatCoord
    integer,intent(out),optional         :: rc

    ! LOCAL VARIABLES
    type(MapDesc)  :: map

    if (present(rc)) rc = ESMF_SUCCESS

    map = MapDesc(trim(mapName),.FALSE., &
      minCoords(2),maxCoords(2),minCoords(1),maxCoords(1))

    call NUOPC_FileWriteJNL(varName,dataFile,gridFile,slices, &
      map=map,scale=scale,jnlFile=jnlFile, repeatCoord=repeatCoord,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end subroutine

  subroutine NUOPC_FileWriteJNL_map(varName,dataFile,gridFile,slices,map, &
  scale,jnlFile,repeatCoord,rc)
    ! ARGUMENTS
    character(len=*),intent(in)          :: varName
    character(len=*),intent(in)          :: dataFile
    character(len=*),intent(in)          :: gridFile
    integer,intent(in)                   :: slices(:)
    type(MapDesc),intent(in)             :: map
    real,intent(in),optional             :: scale(3)
    character(len=*),intent(in),optional :: jnlFile
    logical,intent(in),optional          :: repeatCoord
    integer,intent(out),optional         :: rc

    ! LOCAL VARIABLES
    type(ESMF_VM)                   :: vm
    integer                         :: lpe
    character(len=64)               :: ljnlFile
    logical                         :: lrepeatCoord
    character(len=64)               :: fileBase
    character(len=64)               :: limits
    character(len=64)               :: levels
    character(len=64)               :: latlon
    integer                         :: markExt
    integer                         :: sIndex
    integer                         :: fUnit
    integer                         :: stat
    character(ESMF_MAXSTR)          :: ferretCmd

    if (present(rc)) rc = ESMF_SUCCESS

    ! Get current VM and pet number
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMGet(vm, localPet=lpe, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

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

    if (present(repeatCoord)) then
      lrepeatCoord = repeatCoord
    else
      lrepeatCoord = .FALSE.
    endif

    if (map%global) then
      limits = ''
    else
      write (limits,"(2(A,F0.3,A,F0.3))") &
        '/vlimits=',map%minLat,':',map%maxLat, &
        '/hlimits=',map%minLon,':',map%maxLon
    endif

    if (lrepeatCoord) then
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
    if (ESMF_STDERRORCHECK(rc)) return
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
