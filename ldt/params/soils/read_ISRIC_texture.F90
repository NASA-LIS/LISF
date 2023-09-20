!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_ISRIC_texture
! \label{read_ISRIC_texture}
!
! !REVISION HISTORY:
!  31 Jul 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_ISRIC_texture( n, num_bins, fgrd, texture_layers )

! !USES:
  use LDT_coreMod
  use LDT_logMod
  use LDT_fileIOMod
  use LDT_gridmappingMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc

#if (defined USE_GDAL)
  use, INTRINSIC:: iso_c_binding
  use fortranc
  use gdal
#endif


  implicit none
! !ARGUMENTS: 
  integer,intent(in)    :: n
  integer,intent(in)    :: num_bins   ! Number of soil types
  real,   intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,   intent(inout) :: texture_layers(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!
! !DESCRIPTION:
!  This subroutine retrieves ISRIC soil texture data and reprojects
!  it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[fgrd]
!   output field with the retrieved soil texture
!  \end{description}
!EOP

  integer :: water_class
  integer :: ftn
  real    :: temp(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real    :: gridcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  logical :: file_exists
  
  integer                :: subpnc, subpnr, glpnc, glpnr
  logical                :: file_read_status
  real                   :: dres
  REAL,ALLOCATABLE       :: zval(:,:)
  REAL,ALLOCATABLE       :: zval2(:,:)
  integer                :: i, t, c, r, line
  integer                :: mi                            
  integer                :: mo                            
  real                   :: param_gridDesc(20)    ! Input parameter grid desc 
  real                   :: subparam_gridDesc(20) ! Subsetted parameter grid desc          
  integer, allocatable   :: lat_line(:,:), lon_line(:,:)

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
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))           ! Output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))           ! Output logical mask (go1)
  real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins) ! Output lis 1d grid
  logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins) ! Output logical mask (go2)

! ___________________________________________________________________
  
   water_class = 14    ! Water class for ISRIC soil texture class

   temp = LDT_rc%udef
   gridcnt = 0.
   fgrd = 0.
   texture_layers = 0.
   param_gridDesc = 0.

   inquire(file=trim(LDT_rc%txtfile(n)), exist=file_exists)
   if(.not.file_exists) then 
     write(LDT_logunit,*) "[ERR] ISRIC texture map ",trim(LDT_rc%txtfile(n))," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
   endif
   write(LDT_logunit,*) "[INFO] Reading ISRIC texture file: ",&
                                trim(LDT_rc%txtfile(n))

#if (defined USE_GDAL)
   file_read_status = .false. 
  
   call GDALAllRegister()
  
   ! Use GDAL routines to open the tiff files:
   driver = gdalgetdriverbyname('Tif'//CHAR(0))
   ds = gdalopen( trim(LDT_rc%txtfile(n))//CHAR(0), GA_ReadOnly)
  
   if( .not.gdalassociated(ds) ) then
      write(LDT_logunit,*) "[ERR] Opening dataset on file ",&
           trim(LDT_rc%txtfile(n))," failed ..."
      call LDT_endrun
   end if
  
   ! Determine Tiff-file grid information:
   ierr = gdalgetgeotransform(ds, gt)                 
   dres = abs(gt(6))
  
   xsize = gdalgetrasterxsize(ds)
   ysize = gdalgetrasterysize(ds)

   CALL gdalapplygeotransform(gt, 0.5_c_double, 0.5_c_double, x2, y1)

   x1 = (xsize-1)*dres + x2
   y2 = y1 - (ysize-1)*dres

 ! Set parameter grid fgrd inputs:
   param_gridDesc(1)  = 0.   ! Latlon "projection" 
   param_gridDesc(2)  = xsize
   param_gridDesc(3)  = ysize
   param_gridDesc(4)  = y2 !+ (dres/2) ! LL lat, for gridcell midpoint
   param_gridDesc(5)  = x2 !+ (dres/2) ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = y1 !- (dres/2)  ! UR lat
   param_gridDesc(8)  = x1 !- (dres/2)  ! UR lon
   param_gridDesc(9)  = dres     ! dy: 0.0020833
   param_gridDesc(10) = dres     ! dx: 0.0020833
   param_gridDesc(20) = 64

   band = gdalgetrasterband(ds, 1)
   if (.NOT.gdalassociated(band)) THEN
      write(LDT_logunit,*) '[ERR] Failed getting raster band from GeoTIFF dataset on file, ',&
           TRIM(LDT_rc%txtfile(n))
      call LDT_endrun()
   endif
    
   allocate(zval(xsize, ysize))
   ierr = gdalrasterio_f(band, GF_Read, 0, 0, zval)
   if (ierr /= 0) THEN
      write(LDT_logunit,*) '[ERR] Reading data from GeoTIFF dataset on file ',&
           TRIM(LDT_rc%txtfile(n))
      call LDT_endrun()
   endif
   call gdalclose(ds)
  
   allocate(zval2(xsize,ysize))
   do r=1,ysize
     do c=1,xsize
        zval2(c,r) = zval(c,ysize-r+1)
     enddo
   enddo
   deallocate(zval)

   ! Map Parameter Grid Info to LIS Target Grid/Projection Info --
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%soiltext_proj, &
        param_gridDesc(:), &
        glpnc, glpnr, subpnc, subpnr, subparam_gridDesc,    &
        lat_line, lon_line )
  
   ! Initialize parameter read-in array:
   allocate( read_inputparm(subpnc, subpnr) )
   read_inputparm = LDT_rc%udef
  
   x_offset = nint((subparam_griddesc(5)-param_gridDesc(5))/&
        param_gridDesc(9)) + 1
   y_offset = nint((subparam_griddesc(4)-param_gridDesc(4))/&
        param_gridDesc(10)) + 1
  
   read_inputparm = zval2(x_offset:x_offset+subpnc, &
                          y_offset:y_offset+subpnr)
  
   deallocate(zval2)
   deallocate(lat_line)
   deallocate(lon_line)

!   open(100,file='test.bin',form='unformatted')
!   write(100) read_inputparm
!   close(100)
!   stop

! -------------------------------------------------------------------
!    AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
   mi = subpnc*subpnr
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
   if( mi .ne. mo .and. LDT_rc%soiltext_gridtransform(n) == "none" ) then
      write(LDT_logunit,*) "[ERR] Spatial transform, 'none', is selected, but number of"
      write(LDT_logunit,*) "  input and output points do not match. Select other spatial"
      write(LDT_logunit,*) "  option (e.g., mode, neighbor, tile, etc.)."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   allocate( li(mi), gi(mi) )
   li = .false.
   gi = LDT_rc%udef
!   gi = float(water_class)
   lo1 = .false.
  
   ! Assign 2-D array to 1-D for grid transformation routines:
   i = 0
   do r = 1, subpnr
      do c = 1, subpnc;  i = i + 1
         gi(i) = read_inputparm(c,r)
         if( gi(i) .ne. 255 ) then 
           li(i) = .true.
         endif
      enddo
   enddo
   deallocate(read_inputparm)

   ! Remap to the USDA/STATSGO classification (tile-only here)
   if( LDT_rc%soiltext_gridtransform(n) == "tile" ) then
      do i = 1, mi
         if(nint(gi(i)).eq.1) then
           gi(i) = 12  !clay 
         elseif(nint(gi(i)).eq.2) then
           gi(i) = 11  !silty clay
         elseif(nint(gi(i)).eq.3) then
           gi(i) = 10  !sandy clay
         elseif(nint(gi(i)).eq.4) then
           gi(i) = 9   !clay loam
         elseif(nint(gi(i)).eq.5) then
           gi(i) = 8   !silty clay loam
         elseif(nint(gi(i)).eq.6) then
           gi(i) = 7   !sandy clay loam
         elseif(nint(gi(i)).eq.7) then
           gi(i) = 6   !loam
         elseif(nint(gi(i)).eq.8) then
           gi(i) = 4   !silty loam
         elseif(nint(gi(i)).eq.9) then
           gi(i) = 3   !sandy loam
         elseif(nint(gi(i)).eq.10) then
           gi(i) = 5   !silt
         elseif(nint(gi(i)).eq.11) then
           gi(i) = 2   !loamy sand
         elseif(nint(gi(i)).eq.12) then
           gi(i) = 1   !sand
         endif
      enddo
   endif
   
   ! Apply the spatial transform option:
   select case( LDT_rc%soiltext_gridtransform(n) )

     ! (1) Single-layer selection:
!     case( "none", "mode", "neighbor" )
     case( "mode", "neighbor" )

     ! Transform parameter from original grid to LIS output grid:
       call LDT_transform_paramgrid(n, LDT_rc%soiltext_gridtransform(n), &
                          subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )

     ! Convert 1D count to 2D grid fgrds:
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n); i = i + 1
             if( go1(i).ne.LDT_rc%udef ) then 
                temp(c,r) = go1(i)
             endif
             ! Remap to the USDA/STATSGO classification
             if(temp(c,r).ne.LDT_rc%udef) then
               if(nint(temp(c,r)).eq.1) then
                 temp(c,r) = 12  !clay 
               elseif(nint(temp(c,r)).eq.2) then
                 temp(c,r) = 11  !silty clay
               elseif(nint(temp(c,r)).eq.3) then
                 temp(c,r) = 10  !sandy clay
               elseif(nint(temp(c,r)).eq.4) then
                 temp(c,r) = 9   !clay loam
               elseif(nint(temp(c,r)).eq.5) then
                 temp(c,r) = 8   !silty clay loam
               elseif(nint(temp(c,r)).eq.6) then
                 temp(c,r) = 7   !sandy clay loam
               elseif(nint(temp(c,r)).eq.7) then
                 temp(c,r) = 6   !loam
               elseif(nint(temp(c,r)).eq.8) then
                 temp(c,r) = 4   !silty loam
               elseif(nint(temp(c,r)).eq.9) then
                 temp(c,r) = 3   !sandy loam
               elseif(nint(temp(c,r)).eq.10) then
                 temp(c,r) = 5   !silt
               elseif(nint(temp(c,r)).eq.11) then
                 temp(c,r) = 2   !loamy sand
               elseif(nint(temp(c,r)).eq.12) then
                 temp(c,r) = 1   !sand
               else
                 print*, '[WARN] ISRIC class unmapped ',temp(c,r)
               endif
             endif
           enddo
        enddo

     ! (2) Estimate TILED soiltexture files (gridcnt):
     case( "tile" )

       ! Calculate total counts for each soil type in each coarse gridcell
       call LDT_transform_paramgrid( n, LDT_rc%soiltext_gridtransform(n), &
            subparam_gridDesc, mi, num_bins, gi, li, mo, go2, lo2 )

       ! Convert 1D to 2D grid arrays:
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n);  i = i + 1
             do t = 1, num_bins
                if( go2(i,t).ne.LDT_rc%udef ) then 
                   gridcnt(c,r,t) = go2(i,t)
                endif
             end do
          enddo
       enddo

     case default
       write(LDT_logunit,*) "[ERR] Selected grid transform option, "//&
                 LDT_rc%soiltext_gridtransform(n)
       write(LDT_logunit,*) "[ERR]  not supported for the ISRIC texture reader."
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun

   end select  ! End grid cnt aggregation method

   deallocate(li)
   deallocate(gi)

!- Bring 2-D Array to 3-D Soil tile space:
   if( LDT_rc%soiltext_gridtransform(n) == "none" .or. &  ! Non-tiled surfaces
       LDT_rc%soiltext_gridtransform(n) == "mode" .or. &
       LDT_rc%soiltext_gridtransform(n) == "neighbor" ) then

   !- Assign soil texture types of less than 0 to an actual texture value:
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
          ! Assign soil texture types of less than 0 to an actual texture value:
!            if( nint(temp(c,r)) .le. 0 ) then
!              temp(c,r) = 13   ! placeholder
!            endif
            if( (nint(temp(c,r)) .ne. water_class  ) .and. &
                (nint(temp(c,r)) .ne. LDT_rc%udef) ) then
               gridcnt(c,r,NINT(temp(c,r))) = 1.0
            endif
         enddo
      enddo
   end if

   ! Estimate fraction of grid (fgrid) represented by soil type::
   call param_index_fgrdcalc( n, LDT_rc%soiltext_proj, &
        LDT_rc%soiltext_gridtransform(n), &
        water_class, num_bins, gridcnt, fgrd )
  
! ---
  write(LDT_logunit,*) "[INFO] Done reading ISRIC soil texture file."

#endif

end subroutine read_ISRIC_texture

