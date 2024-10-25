!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readGIMMSMODIS_NDVIObs
! \label{readGIMMSMODIS_NDVIObs}
!
! !INTERFACE: 
subroutine readGIMMSMODIS_NDVIObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use GIMMSMODIS_NDVIobsMod

#if (defined USE_GDAL)
  use, INTRINSIC:: iso_c_binding
  use fortranc
  use gdal
#endif
          
  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!
! DATASET
!
!        Description     GIMMS MODIS Terra & Aqua NDVI 8-bit 9x9 degree Tiles
!        Coverage        Global 180W to 180E, 80N to 60S
!        Projection      Lat/Lon grid
!        Spatial Res.    0.00225 x 0.00225 degrees (~250-meter)
!        Temporal Res.   8-day composite period
!        File Format     TIFF
!        Download        http://gimms.gsfc.nasa.gov/download/MODIS/
!                        ftp://gimms.gsfc.nasa.gov/MODIS/
!
!  NEAR REAL-TIME AND STANDARD PROCESSING
!
!        Forward near real-time (nrt) processing of this dataset is typically
!        available within 12 hours from the last day of the 8-day composite
!        period.  Science quality, standard (std) processing is typically
!        available one or more days after the last day of the 8-day composite
!        period. The two processing modes differ in their upstream data source.
!
!        Processing          Data Source
!        ----------          -----------
!        nrt                 LANCE near real-time
!        std                 MODAPS science quality
!
!        * GIMMS nrt datasets are only available for 30 days after production.
!
!   Qualifier           Value       Description
!   ---------           -----       -----------
!   Dataset             GMOD09Q1    GIMMS MODIS Terra (MOD=Terra MYD=Aqua)
!   Start date          A2010001    Starting year 2010, day-of-year 001
!   Composite period    08d         8-day composite
!   Projection          latlon      Lat/Lon grid
!   9x9 Tile index      x00y02      Column 00, Row 02
!   Versions            5v3         MODAPS collection 5, GIMMS version 3
!   Layer name          NDVI        Normalized Vegetation Index
!   File format         tif         Tagged Image File
!
! PROJECTION
!
!   This dataset is mapped to a global 180W to 180E, 90N to 90S Lat/Lon grid
!    and divided into 9x9 degree tiles.
!
!   Parameter           Value       Units
!   ---------           -----       -----
!   Datum               WGS84       -
!   Upper left Lat      90.0        deg
!   Upper left Lon      -180.0      deg
!   Pixel size          0.00225     deg
!   Grid x size         160000      pixel
!   Grid y size         80000       pixel
!
!   Coordinates are mapped and measured from the upper left corner of the
!     grid cell.
!
!   Since the size of a global grid file may be too large for users,
!     the grid is divided into 40 columns and 20 rows (starting at the
!     upper left corner at 180W 90N) to create 9x9 degree tiles of
!     4000 x 4000 pixels.
!
!     Parameter                            Value
!     ---------                            -----
!     Tile size             9 / 0.00225 =  4000
!     No. columns (x)     160000 / 4000 =  40
!     No. rows (y)         80000 / 4000 =  20
!
! GIMMS-MODIS NDVI tile files:
!- To calculate the upper left corner {Lon,Lat} for tile {x,y}:
!   UL-Lon = -180 + (x * 9)
!   UL-Lat =   90 - (y * 9)
!
!  NOTES: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  30 June 2016: Sujay Kumar & Kristi Arsenault, Initial Specification
! 
!EOP
  integer, parameter     :: tile_nc=4000, tile_nr=4000
  integer                :: c,r, tindex
  integer                :: flag
  integer                :: ftn
  character*200          :: fname
  logical*1              :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: ndvi_out(LVT_rc%lnc*LVT_rc%lnr)
  logical*1, allocatable :: input_bitmap(:)
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  logical                :: file_exists
  logical                :: file_read_status
  integer                :: iret
  integer                :: nc, nr
  integer                :: ii

  real                   :: dres
  integer                :: tile_x, west_x, east_x, x_cnt
  integer                :: tile_y, south_y, north_y, y_cnt
  integer                :: Xdir_numtiles, Ydir_numtiles
  real                   :: cornerlat1, cornerlat2
  real                   :: cornerlon1, cornerlon2
  real, allocatable      :: read_ndvi_uint(:,:)   ! for unsigned integer*8
  real, allocatable      :: read_ndvi(:,:,:)   
  real, allocatable      :: mosaic_ndvi(:,:)
  real, allocatable      :: ndvi_inp(:)
  integer                :: col_cnt1, col_cnt2, row_cnt1, row_cnt2

#if (defined USE_GDAL)
  TYPE(gdaldriverh)      :: driver
  TYPE(gdaldataseth)     :: ds
  TYPE(gdalrasterbandh)  :: band
  INTEGER(kind=c_int)    :: xsize, ysize
  REAL(kind=c_double)    :: x1, y1, x2, y2, gt(6)
  INTEGER(kind=c_int)    :: i, j, k, ierr
#endif
  
  varfield = LVT_rc%udef
  
#if (defined USE_GDAL) 

  file_read_status = .false. 

  call GDALAllRegister()
!
! - USE THIS APPROACH FOR SELECTING THE TILED FILES TO OPEN -
! - To calculate the {x,y} tile index for a given {Lon,Lat}:
!   x = floor((180 + Lon) / 9)
!   y = floor(( 90 - Lat) / 9)
!
! TBD: Assumes that LVT is running in the lat/lon projection. 
  
  west_x  = floor((180 + LVT_rc%gridDesc(5)) / 9)
  south_y = floor(( 90 - LVT_rc%gridDesc(4)) / 9)
  east_x  = floor((180 + LVT_rc%gridDesc(8)) / 9)
  north_y = floor(( 90 - LVT_rc%gridDesc(7)) / 9)
  
  Xdir_numtiles = (east_x-west_x)+1
  Ydir_numtiles = (south_y-north_y)+1
 
  allocate( read_ndvi(tile_nc, tile_nr, Xdir_numtiles) )
  allocate( mosaic_ndvi(Xdir_numtiles*tile_nc, Ydir_numtiles*tile_nr) )
  
  read_ndvi = -9999.
  mosaic_ndvi = -9999.

! Derive full mosaicked GIMMS tiles domain corner points == 
!  used for setting the parameter array for reprojection/interp
!  to LVT eval or LVT runnin domain:
! Initialize lat/lon points:
  cornerlat1 =  9999. ! [LL lat of mosaicked tile domain]
  cornerlat2 = -9999. ! [UR lat of mosaicked tile domain]
  cornerlon1 =  9999. ! [LL lon of mosaicked tile domain]
  cornerlon2 = -9999. ! [UR lon of mosaicked tile domain]
!- Loop over relevant tiles and mosaic them together:
  y_cnt = 0
  do tile_y = north_y, south_y
     y_cnt = y_cnt + 1
     row_cnt1 = (tile_y-north_y)*tile_nr+1
     row_cnt2 = row_cnt1+tile_nr-1
     
     x_cnt = 0
     do tile_x = west_x, east_x
        x_cnt = x_cnt + 1 
        col_cnt1 = (tile_x-west_x)*tile_nc+1
        col_cnt2 = col_cnt1+tile_nc-1
        
        call create_GIMMSMODISndvi_filename( &
             gimmsmodisndviobs(source)%odir, &
             LVT_rc%dyr(source),             &
             LVT_rc%ddoy(source),            &
             tile_x,                         &
             tile_y,                         &
             fname)
		
        inquire(file=trim(fname),exist=file_exists) 
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading GIMMS MODIS file '&
                //trim(fname)
           ! Use GDAL routines to open the tiff files:
           driver = gdalgetdriverbyname('Tif'//CHAR(0))
           ds = gdalopen( trim(fname)//CHAR(0), GA_ReadOnly)
           
           
           if( .not.gdalassociated(ds) ) then
              write(LVT_logunit,*) "[ERR] opening dataset on file ",&
                   trim(fname)," failed ..."
              call LVT_endrun()
           end if
           
           ierr = gdalgetgeotransform(ds, gt)                 
           dres = abs(gt(6))
           
           xsize = gdalgetrasterxsize(ds)
           ysize = gdalgetrasterysize(ds)
           
           CALL gdalapplygeotransform(gt, 0.5_c_double, 0.5_c_double, x1, y1)
           
           x2 = (xsize-1)*dres + x1
           y2 = y1 - (ysize-1)*dres
           
           cornerlat1 = min(cornerlat1, y2)  ! [LL lat of mosaicked tile domain]
           cornerlat2 = max(cornerlat2, y1)  ! [UR lat of mosaicked tile domain]
           cornerlon1 = min(cornerlon1, x1)  ! [LL lon of mosaicked tile domain]
           cornerlon2 = max(cornerlon2, x2)  ! [UR lon of mosaicked tile domain]
           
           ! Retrieve NDVI band/layer (layer=1) --
           band = gdalgetrasterband(ds, 1)
           if (.NOT.gdalassociated(band)) THEN
              write(LVT_logunit,*) &
                   "[ERR] failed getting NDVI band from TIFF dataset in file, ",&
                   TRIM(fname)
              call LVT_endrun()
           endif
           
           ! Read in NDVI layer input:
           allocate( read_ndvi_uint(tile_nc, tile_nr) )
           ierr = gdalrasterio_f( band, GF_Read, 0, 0, read_ndvi_uint )
           if (ierr /= 0) then
              write(LVT_logunit,*) &
                   '[ERR] reading data from TIFF dataset on file ',TRIM(fname)
              call LVT_endrun()
           endif
           !- Convert data format (UINT8) to regular floating point::      
           do j = 1, tile_nr
              do i = 1, tile_nc
                 !- Convert the UINT8 values to reals::
                 if ( read_ndvi_uint(i,j) .le. 250 ) then
                    read_ndvi(i,j,x_cnt) = read_ndvi_uint(i,j) * 0.004
                 else
                    read_ndvi(i,j,x_cnt) = read_ndvi_uint(i,j) 
                 endif
              end do
           end do
           deallocate( read_ndvi_uint )
           
           !- Mosaic NDVI tiles together:
           mosaic_ndvi(col_cnt1:col_cnt2, row_cnt1:row_cnt2) = read_ndvi(:,:,x_cnt)

           call gdalclose(ds)
           file_read_status = .true. 
        endif
     end do
  enddo
  deallocate(read_ndvi)


  if(file_read_status) then 
     allocate( ndvi_inp(Xdir_numtiles*tile_nc*Ydir_numtiles*tile_nr) )
     allocate( input_bitmap(Xdir_numtiles*tile_nc*Ydir_numtiles*tile_nr) )
     
     ndvi_inp = LVT_rc%udef
     input_bitmap = .false. 
     ii = 0
     do r = Ydir_numtiles*tile_nr, 1, -1
        ii = ii + 1
        do c = 1, Xdir_numtiles*tile_nc
           if(mosaic_ndvi(c,r).gt.240) then !TBD - check
              mosaic_ndvi(c,r) = LVT_rc%udef
           endif
           ndvi_inp(c+(ii-1)*Xdir_numtiles*tile_nc) = mosaic_ndvi(c,r)
           if(mosaic_ndvi(c,r).ne.LVT_rc%udef) then 
              input_bitmap(c+(ii-1)*Xdir_numtiles*tile_nc) = .true.
           endif
        end do
     end do
  endif
  
  deallocate(mosaic_ndvi)
        
  if((gimmsmodisndviobs(source)%startFlag.and.file_read_status).or.&
       LVT_rc%resetFlag(source))then 

     LVT_rc%resetFlag(source) = .false. 
     gimmsmodisndviobs(source)%gridDesc = 0
     
     if(file_read_status) then 
        gimmsmodisndviobs(source)%startFlag = .false. 

        gimmsmodisndviobs(source)%nc = nint((cornerlon2 - cornerlon1)/dres)+1
        gimmsmodisndviobs(source)%nr = nint((cornerlat2 - cornerlat1)/dres)+1
     
        nc = GIMMSMODISndviobs(source)%nc
        nr = GIMMSMODISndviobs(source)%nr


        !filling the items needed by the interpolation library
        gimmsmodisndviobs(source)%gridDesc(1) = 0  
        gimmsmodisndviobs(source)%gridDesc(2) = gimmsmodisndviobs(source)%nc
        gimmsmodisndviobs(source)%gridDesc(3) = gimmsmodisndviobs(source)%nr
        gimmsmodisndviobs(source)%gridDesc(4) = cornerlat1
        gimmsmodisndviobs(source)%gridDesc(5) = cornerlon1
        gimmsmodisndviobs(source)%gridDesc(7) = cornerlat2
        gimmsmodisndviobs(source)%gridDesc(8) = cornerlon2
        gimmsmodisndviobs(source)%gridDesc(6) = 128
        gimmsmodisndviobs(source)%gridDesc(9) = dres
        gimmsmodisndviobs(source)%gridDesc(10) = dres
        gimmsmodisndviobs(source)%gridDesc(20) = 64
     
        gimmsmodisndviobs(source)%datares  = dres
     
        if(LVT_isAtAfinerResolution(gimmsmodisndviobs(source)%datares)) then
           allocate(gimmsmodisndviobs(source)%rlat(LVT_rc%lnc*LVT_rc%lnr))
           allocate(gimmsmodisndviobs(source)%rlon(LVT_rc%lnc*LVT_rc%lnr))
           allocate(gimmsmodisndviobs(source)%n11(LVT_rc%lnc*LVT_rc%lnr))
           
           call neighbor_interp_input(gimmsmodisndviobs(source)%gridDesc, &
                LVT_rc%gridDesc, &
                LVT_rc%lnc*LVT_rc%lnr, &
                gimmsmodisndviobs(source)%rlat, &
                gimmsmodisndviobs(source)%rlon, &
                gimmsmodisndviobs(source)%n11)
        else
           allocate(gimmsmodisndviobs(source)%n11(&
             gimmsmodisndviobs(source)%nc*gimmsmodisndviobs(source)%nr))
           call upscaleByAveraging_input(gimmsmodisndviobs(source)%gridDesc,&
                LVT_rc%gridDesc,gimmsmodisndviobs(source)%nc*&
                gimmsmodisndviobs(source)%nr,&
                LVT_rc%lnc*LVT_rc%lnr,gimmsmodisndviobs(source)%n11)
        endif
     endif
  endif
  if(file_read_status) then 
     if(LVT_isAtAfinerResolution(gimmsmodisndviobs(source)%datares)) then
        call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
             ndvi_inp, output_bitmap, ndvi_out, &
             GIMMSMODISndviobs(source)%nc*&
             GIMMSMODISndviobs(source)%nr, &
             LVT_rc%lnc*LVT_rc%lnr, &
             gimmsmodisndviobs(source)%rlat, & 
             gimmsmodisndviobs(source)%rlon, &
             gimmsmodisndviobs(source)%n11, &
             LVT_rc%udef, iret)
        
     else
        call upscaleByAveraging(&
             GIMMSMODISndviobs(source)%nc*&
             GIMMSMODISndviobs(source)%nr, &
             LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
             gimmsmodisndviobs(source)%n11, input_bitmap, &
             ndvi_inp, output_bitmap, ndvi_out)
        
     endif
     

     do r=1, LVT_rc%lnr
        do c=1, LVT_rc%lnc
           if(ndvi_out(c+(r-1)*LVT_rc%lnc).gt.0) then 
              varfield(c,r) = ndvi_out(c+(r-1)*LVT_rc%lnc)
           else
              varfield(c,r) = LVT_rc%udef
           endif
        enddo
     enddo
     deallocate(ndvi_inp)
     deallocate(input_bitmap)

  else
     varfield  = LVT_rc%udef
  endif
  
#endif

!if(LVT_MOC_NDVI(source).ge.1) then
  call LVT_logSingleDataStreamVar(LVT_MOC_NDVI,source,varfield,&
       vlevel=1,units="-")
!endif

end subroutine readGIMMSMODIS_NDVIObs

!BOP
! 
! !ROUTINE: create_GIMMSMODISndvi_filename
! \label{create_GIMMSMODISndvi_filename}
!
! !INTERFACE: 
subroutine create_GIMMSMODISndvi_filename(odir,yr,doy,x,y,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for GIMMS MODIS NDVI data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the GIMMSMODISNDVI_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: doy
  integer                      :: x
  integer                      :: y
  character(len=*)             :: filename
!
!EOP
  
  character*4             :: fyr
  character*3             :: fdoy
  character*2             :: fx
  character*2             :: fy

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy
  write(unit=fx, fmt='(i2.2)') x
  write(unit=fy, fmt='(i2.2)') y

  filename = trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//'/'//&
       'GMOD09Q1.A'//trim(fyr)//trim(fdoy)//'.08d.latlon.x'//&
       trim(fx)//'y'//trim(fy)//'.6v1.NDVI.tif'

end subroutine create_GIMMSMODISndvi_filename


