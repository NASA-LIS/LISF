!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_GIA_irrigfrac
!  \label{read_GIA_irrigfrac}
!
! !REVISION HISTORY:
!  20 May 2019: H. Beaudoing;  Initial specification based on Sarith's code
!
! !INTERFACE:
 subroutine read_GIA_irrigfrac(n, fgrd) 

! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_irrigationMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc

#if (defined USE_GDAL)
    use, INTRINSIC:: iso_c_binding
    use fortranc
    use gdal
#endif
!EOP      

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n))

!
! !DESCRIPTION:
! Global Irrigated Areas (GIA), 30 arcsec (0.008333333 degree) global irrigated
! areas from 4 different sources. 
! This dataset downscaled the FAO Global Map of Irrigated Areas (GMIA) version
! 5.0 (Siebert et al., 2013) from 5 arcmin to 30 arcsec using the MERIS NDVI 
! data (ESA, 2007) at 10 arcsec global.  The ESA-CCI NDVI 7-day product 
! over 1999-2012 (ESA,2015) was used to detect additional irrigated areas based
! on the annual min/max, yearly course of NDVI, and the length of growing 
! period. Further, agricultral suitability data of Zabel et al. 2014 that 
! describes suibatility for 16 staple crops according to climate, soil, and 
! topography conditions and the potential number of crop cycles per year was 
! used to identify unknown irrigated areas.  The infomation on clopland is 
! takecn from the ESA-CCI-LC product
! (cropland rain-fed, cropland irrigated, mosaic cropland>50%; ESA, 2015) and
! from the predecessor GlobCover (ESA, 2010; post-flooding or irrigated 
! croplands, rain-fed croplands, mosaic cropland; 50–70%). 
! 
! Consider pixels haveing values 1-3 as irrigated (exclude 4 due to different
! choice of land cover data).  
!
!  The legend is:
!  0 = no irrigated area
!  1 = downscaled Siebert et al. 2013
!  2 = low agricultural suitability, high NDVI and NDVI course of vegetation
!  3 = potential multiple cropping < actual multiple cropping
!  4 = classified as cropland (according to ESA-CCI-LC and GlobCover) and 
!      low suitability
!
!  Ref: 
! Meier, J., Zabel, F., and Mauser, W.: A global approach to estimate irrigated ! areas – a comparison between different data and statistics, Hydrol. Earth 
! Syst. Sci., 22, 1119-1133, https://doi.org/10.5194/hess-22-1119-2018, 2018. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[fgrd]
!   output field MODIS irrigation type
!  \end{description}
!
!EOP      

! GIA crop-irrigation type:
!  integer, parameter :: input_cols = 43200
!  integer, parameter :: input_rows = 18000
!  real,    parameter :: IN_xres = 0.0083333333333333
!  real,    parameter :: IN_yres = 0.0083333333333333
  integer, parameter :: num_bins = 1

  integer   :: nc, nr, i, t, nrec
  integer   :: ftn
  logical   :: file_exists
  integer   :: mi                     ! Total number of input param grid array points
  integer   :: mo                     ! Total number of output LIS grid array points
  integer   :: glpnc, glpnr           ! Parameter (global) total columns and rows
  integer   :: subpnc, subpnr         ! Parameter subsetted columns and rows
  real      :: param_gridDesc(20)     ! Input parameter grid desc fgrd
  real      :: subparam_gridDesc(20)  ! Input parameter grid desc fgrd
  integer, allocatable  :: lat_line(:,:), lon_line(:,:)
  integer, allocatable  :: n11(:)     ! Array that maps the location of each input grid
                                      !   point in the output grid. 
  real,    allocatable  :: gi(:)      ! Input parameter 1d grid
  logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)

  real, allocatable :: var_read(:), var_in(:,:)
  real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins)    ! Output lis 1d grid
  logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins)    ! Output logical mask (to match go)
  real      :: irrigfrac_cnt(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real      :: irrigfrac(LDT_rc%lnc(n),LDT_rc%lnr(n))

  integer   :: noncrop 
  real      :: isum
  real      :: DXY
  integer   :: DC, DR

    logical                :: file_read_status
    real                   :: dres
    REAL,ALLOCATABLE       :: zval(:,:)
    REAL,ALLOCATABLE       :: zval2(:,:)
    integer                :: c, r, line

#if (defined USE_GDAL)
    TYPE(gdaldriverh)      :: driver
    TYPE(gdaldataseth)     :: ds
    TYPE(gdalrasterbandh)  :: band
    INTEGER(kind=c_int)    :: xsize, ysize
    REAL(kind=c_double)    :: x1, y1, x2, y2, gt(6)
    INTEGER(kind=c_int)    :: ierr
#endif

    integer               :: x_offset, y_offset
! ______________________________________________________________

   noncrop = 0

!- Check if file is present:
   inquire(file=trim(LDT_irrig_struc(n)%irrigfracfile), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Irrigation fraction map ",trim(LDT_irrig_struc(n)%irrigfracfile)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(unit=LDT_logunit,fmt=*) "[INFO] Reading GIA irrigated area map file"

!- Open file to be processed::
#if (defined USE_GDAL)
    file_read_status = .false.

    call GDALAllRegister()

! Use GDAL routines to open the tiff files:
    driver = gdalgetdriverbyname('Tif'//CHAR(0))
    ds = gdalopen( trim(LDT_irrig_struc(n)%irrigfracfile)//CHAR(0), GA_ReadOnly)

    if( .not.gdalassociated(ds) ) then
       write(LDT_logunit,*) "[ERR] Opening dataset file, ",&       
            trim(LDT_irrig_struc(n)%irrigfracfile),", failed ..."
       stop
    end if

    ierr = gdalgetgeotransform(ds, gt)
    dres = abs(gt(6))

    xsize = gdalgetrasterxsize(ds)
    ysize = gdalgetrasterysize(ds)

    CALL gdalapplygeotransform(gt, 0.5_c_double, 0.5_c_double, x2, y1)

    band = gdalgetrasterband(ds, 1)
    if (.NOT.gdalassociated(band)) THEN
       write(*,*) '[ERR] Failed getting raster band from TIF dataset on file,',&
             trim(LDT_irrig_struc(n)%irrigfracfile)
       call LDT_endrun()
       stop
    endif
    
    allocate(zval(xsize, ysize))
    ierr = gdalrasterio_f(band, GF_Read, 0, 0, zval)
    if (ierr /= 0) THEN
       write(*,*) '[ERR] Reading data from TIF dataset on file, ', &
                   trim(LDT_irrig_struc(n)%irrigfracfile)
       call LDT_endrun()
       stop
    endif
    call gdalclose(ds)

    x1 = (xsize-1)*dres + x2
    y2 = y1 - (ysize-1)*dres

!- Reverse-Y 
    allocate(zval2(xsize,ysize))
    do r=1,ysize
       do c=1,xsize
          zval2(c,r) = zval(c,ysize-r+1)
       enddo
    enddo

! _____________________________
!- Set parameter grid array inputs:
   LDT_irrig_struc(n)%irrig_proj  = "latlon"
   param_gridDesc(1)  = 0.          ! Latlon
   param_gridDesc(2)  = xsize
   param_gridDesc(3)  = ysize
   param_gridDesc(4)  = y2 ! LL lat (-59.9958353680719)
   param_gridDesc(5)  = x2 ! LL lon (-179.995833333100)
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = y1  ! UR lat
   param_gridDesc(8)  = x1  ! UR lon
   param_gridDesc(9)  = dres     ! dy: 0.008333333
   param_gridDesc(10) = dres     ! dx: 0.008333333
   param_gridDesc(20) = 64

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_irrig_struc(n)%irrig_proj, param_gridDesc(:), &
           glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )

!   print*,'glpnc,glpnr,subpnc,subpnr:',glpnc, glpnr, subpnc, subpnr
! _________
   allocate( var_in(subpnc,subpnr) )
   var_in = LDT_rc%udef

   x_offset = nint((subparam_gridDesc(5)-param_gridDesc(5))/&
          param_gridDesc(9)) + 1
   y_offset = nint((subparam_gridDesc(4)-param_gridDesc(4))/&
          param_gridDesc(10)) + 1
!   print*,'x_offset,y_offset:',x_offset,y_offset

   var_in = zval2(x_offset:x_offset+subpnc-1, y_offset:y_offset+subpnr-1)
          
   deallocate( zval )
   deallocate( zval2 )

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

   mi = subpnc*subpnr
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)

#if 0
   if( mi .ne. mo .and. LDT_rc%irrigfrac_gridtransform(n) == "none" ) then
      write(LDT_logunit,*) "[ERR] Spatial transform, 'none', is selected, but number of"
      write(LDT_logunit,*) "  input and output points do not match. Select other spatial"
      write(LDT_logunit,*) "  option (e.g., average, neighbor, etc.)."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
#endif

   allocate( li(mi), gi(mi), n11(mi) )
   gi  = noncrop
   li  = .false.
   lo2 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do nr = 1, subpnr
      do nc = 1, subpnc
         i = i + 1
!  reset pixel value to be 0=no, 1=irrigated for aggregation
         if( var_in(nc,nr) > 0. .and. var_in(nc,nr) < 4 ) then  
            gi(i) = 1.0
            li(i) = .true.
         endif 
      enddo
   enddo
   deallocate( var_in )

!- Create mapping between parameter domain and LIS grid domain:
   call upscaleByAveraging_input( subparam_gridDesc, &
                           LDT_rc%gridDesc(n,:), mi, mo, n11 )

!- Calculate total counts for inland water type in each coarse gridcell:
   call upscaleByCnt( mi, mo, 1, noncrop, n11, li, gi, lo2, go2 )

   irrigfrac_cnt = 0.
   i = 0
   do nr = 1, LDT_rc%lnr(n)
      do nc = 1, LDT_rc%lnc(n)
         i = i + 1
         irrigfrac_cnt(nc,nr) = go2(i,1)
      enddo
   enddo

!- Compute irrigration fraction :
   DXY = 360. / real(xsize)
   DC = nint (LDT_rc%gridDesc(n,10) / DXY )
   DR = nint (LDT_rc%gridDesc(n,9) / DXY )
   irrigfrac = 0.
   do nr = 1, LDT_rc%lnr(n)
      do nc = 1, LDT_rc%lnc(n)

       ! Calculate gridpoint GIA pixel total:
         isum = real ( DC * DR )
!         write(500,*) nc,nr, isum, irrigfrac_cnt(nc,nr,2),irrigfrac_cnt(nc,nr,3)

       ! Estimate 2-D fraction for just irrigation:
         if( isum > 0. ) then
            irrigfrac(nc,nr) = (irrigfrac_cnt(nc,nr) / isum) * 100.   
         else
            irrigfrac(nc,nr) = 0.
         endif
 
      enddo
   end do

   fgrd = irrigfrac

  write(LDT_logunit,fmt=*) "[INFO] Done reading GIA irrigation gridcell fractions"
#endif

end subroutine read_GIA_irrigfrac

