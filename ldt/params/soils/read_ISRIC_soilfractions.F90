!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_ISRIC_soilfractions
! \label{read_ISRIC_soilfractions}
!
! !REVISION HISTORY:
!  09 Jun 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_ISRIC_soilfractions(n, num_bins, soilsfgrd, &
       sandave, clayave, siltave )
! !USES:
    use LDT_coreMod
    use LDT_coreMod
    use LDT_logMod


    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: n          ! nest index
    integer, intent(in)   :: num_bins   ! number of bins for tiling
    real,    intent(out)  :: soilsfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
    real,    intent(out)  :: sandave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
    real,    intent(out)  :: clayave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
    real,    intent(out)  :: siltave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine retrieves ISRIC soil, silt, clay fraction data 
!   and reprojects it to the LDT target grid. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of bins (or types/classes)
!  \item[soilsfgrd]
!   output field with merged soils fraction of grid
!  \item[sandave]
!   output field with average sand binned fraction 
!  \item[clayave]
!   output field with average clay binned fraction 
!  \item[siltave]
!   output field with average silt binned fraction 
!  \end{description}
!EOP
    integer :: ftn_sa, ftn_cl, ftn_si
    integer :: c,r
    logical :: file_exists
!________________________________________________________________________

    write(LDT_logunit,*) "[INFO] Reading ISRIC sand, clay and silt files: ",&
         trim(LDT_rc%safile(n)),", ",&
         trim(LDT_rc%clfile(n)),",",&
         trim(LDT_rc%sifile(n)) 

   sandave  = 0.; clayave = 0.; siltave = 0.
   soilsfgrd = 0.    
! -------------------------------------------------------------------
!  CHECK FOR AND OPEN SOIL FRACTION FILES
! -------------------------------------------------------------------
!- Sand file:
    inquire(file=trim(LDT_rc%safile(n)), exist=file_exists)
    if(.not.file_exists) then 
       write(LDT_logunit,*) "Sand map ",trim(LDT_rc%safile(n))," not found."
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
    endif
    
    call processISRICsoilf(n, LDT_rc%safile(n), num_bins, sandave)
    call processISRICsoilf(n, LDT_rc%clfile(n), num_bins, clayave)
    call processISRICsoilf(n, LDT_rc%sifile(n), num_bins, siltave)

    ! Make sure and negative values are set to universal undefined value:
    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)
          if( sandave(c,r,1) < 0. ) sandave(c,r,1) = LDT_rc%udef
          if( clayave(c,r,1) < 0. ) clayave(c,r,1) = LDT_rc%udef
          if( siltave(c,r,1) < 0. ) siltave(c,r,1) = LDT_rc%udef
       enddo
    enddo
    soilsfgrd(:,:,1) = 1.0

    write(LDT_logunit,*) "[INFO] Done reading ISRIC soil fraction files."
    
  end subroutine read_ISRIC_soilfractions

!BOP
! !ROUTINE: processISRICsoilf
! 
! !DESCRIPTION:
!  This subroutine reads a given ISRIC soil fraction
!  file (sand, clay, silt) and processes the data to 
!  the LDT target grid. 
!
! !INTERFACE: 
  subroutine processISRICsoilf(n, filename, num_bins, output_field)
! !USES:
    use LDT_coreMod
    use LDT_gridmappingMod
    use LDT_paramTileInputMod
    use LDT_fileIOMod
    use LDT_logMod

#if (defined USE_GDAL)
    use, INTRINSIC:: iso_c_binding
    use fortranc
    use gdal
#endif

    implicit none
! !ARGUMENTS: 
    integer,            intent(in)   :: n 
    character(len=*),   intent(in)   :: filename
    integer,            intent(in)   :: num_bins
    real                             :: output_field(LDT_rc%lnc(n),&
         LDT_rc%lnr(n),num_bins)
!EOP  

    integer, parameter      :: numtiles = 1

    integer                 :: subpnc, subpnr, glpnc, glpnr
! Total number of input param grid array points
    integer                 :: mi                            
! Total number of output LIS grid array points
    integer                 :: mo                            
! Input parameter grid desc array
    real                    :: subparam_gridDesc(20)         
    integer, allocatable    :: lat_line(:,:), lon_line(:,:)
! Maps each input grid point to output grid.
    integer, allocatable    :: n11(:)           
! input parameter 1d grid
    real,    allocatable    :: gi1(:), gi2(:)   
! input logical mask (to match gi)
    logical*1,allocatable   :: li1(:), li2(:)   

    logical                :: file_read_status
    real                   :: dres
    REAL,ALLOCATABLE       :: zval(:,:)
    REAL,ALLOCATABLE       :: zval2(:,:)
    integer                :: i, t, c, r, line

#if (defined USE_GDAL)
    TYPE(gdaldriverh)      :: driver
    TYPE(gdaldataseth)     :: ds
    TYPE(gdalrasterbandh)  :: band
    INTEGER(kind=c_int)    :: xsize, ysize
    REAL(kind=c_double)    :: x1, y1, x2, y2, gt(6)
    INTEGER(kind=c_int)    :: ierr
#endif
    
    integer               :: x_offset, y_offset
! Read input parameter
    real,    allocatable  :: read_inputparm(:,:)   
! input parameter 1d grid
    real,    allocatable  :: gi(:)      
! input logical mask (to match gi)
    logical*1,allocatable :: li(:)      
! output lis 1d grid
    real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n),numtiles)  
! output logical mask (to match go1)
    logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n),numtiles)  


#if (defined USE_GDAL)
    file_read_status = .false. 
    
    call GDALAllRegister()
    
! Use GDAL routines to open the tiff files:
    driver = gdalgetdriverbyname('Tif'//CHAR(0))
    ds = gdalopen( trim(filename)//CHAR(0), GA_ReadOnly)
    
    if( .not.gdalassociated(ds) ) then
       write(LDT_logunit,*) "[ERR] Opening dataset file, ",&
            trim(filename),", failed ..."
       stop
    end if
    
    ierr = gdalgetgeotransform(ds, gt)                 
    dres = abs(gt(6))
    
    xsize = gdalgetrasterxsize(ds)
    ysize = gdalgetrasterysize(ds)
    
    CALL gdalapplygeotransform(gt, 0.5_c_double, 0.5_c_double, x2, y1)
    
    band = gdalgetrasterband(ds, 1)
    if (.NOT.gdalassociated(band)) THEN
       write(LDT_logunit,*) '[ERR] Failed getting raster band from GeoTIFF dataset on file, ',&
             TRIM(filename)
       call LDT_endrun()
    endif
    
    allocate(zval(xsize, ysize))
    ierr = gdalrasterio_f(band, GF_Read, 0, 0, zval)
    if (ierr /= 0) THEN
       write(LDT_logunit,*) '[ERR] Reading data from GeoTIFF dataset on file, ',TRIM(filename)
       call LDT_endrun()
    endif
    call gdalclose(ds)
    
    dres = abs(gt(6))
    
    x1 = (xsize-1)*dres + x2
    y2 = y1 - (ysize-1)*dres  
    
    allocate(zval2(xsize,ysize))
    do r=1,ysize
       do c=1,xsize
          zval2(c,r) = zval(c,ysize-r+1)
       enddo
    enddo
    deallocate(zval)
    
    !- Map Parameter Grid Info to LIS Target Grid/Projection Info --
    subparam_gridDesc = 0.
    call LDT_RunDomainPts( n, LDT_rc%soils_proj, LDT_rc%soil_gridDesc(n,:), &
         glpnc, glpnr, subpnc, subpnr, subparam_gridDesc,    &
         lat_line, lon_line )
    

  !- Initialize parameter read-in array:
     allocate( read_inputparm(subpnc, subpnr) )
     read_inputparm = LDT_rc%udef

     x_offset = nint((subparam_griddesc(5)-LDT_rc%soil_gridDesc(n,5))/&
          LDT_rc%soil_gridDesc(n,9)) + 1
     y_offset = nint((subparam_griddesc(4)-LDT_rc%soil_gridDesc(n,4))/&
          LDT_rc%soil_gridDesc(n,10)) + 1
     
     read_inputparm = zval2(x_offset:x_offset+ subpnc, &
          y_offset:y_offset+subpnr)

     deallocate(zval2)

     mi = subpnc*subpnr
     allocate( li(mi), gi(mi) )
     li = .false.
     gi = LDT_rc%udef
     mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
     lo1 = .false.

  !- Assign 2-D array to 1-D for grid transformation routines:
     i = 0
     do r = 1, subpnr
        do c = 1, subpnc;  i = i + 1
           gi(i) = read_inputparm(c,r)
           if( gi(i) .ne. 255 ) li(i) = .true.
        enddo
     enddo

     call LDT_transform_paramgrid( n, LDT_rc%soils_gridtransform(n), &
          subparam_gridDesc, mi, numtiles, gi, li, mo, go1, lo1 )
     
!- Convert 1D to 2D grid arrays:
     if(  LDT_rc%soils_gridtransform(n) .ne. "none" ) then
        i = 0
        do r = 1, LDT_rc%lnr(n)
           do c = 1, LDT_rc%lnc(n);  i = i + 1
              do t = 1, numtiles
                 if(go1(i,t).ne.-9999.0) then 
                    output_field(c,r,t) = go1(i,t)/100
                 endif
              end do
           enddo
        enddo
     end if

!     open(100,file='test_out.bin',form='unformatted')
!     write(100) output_field
!     close(100)
!     stop

     deallocate(li)
     deallocate(gi)
     deallocate(read_inputparm)
     deallocate(lat_line)
     deallocate(lon_line)

#endif

   end subroutine processISRICsoilf
  
