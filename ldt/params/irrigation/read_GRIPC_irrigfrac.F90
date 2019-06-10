!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_GRIPC_irrigfrac
!  \label{read_GRIPC_irrigfrac}
!
! !REVISION HISTORY:
!  30 Sep 2013: K. Arsenault;  Added new irrigation type map
!
! !INTERFACE:
 subroutine read_GRIPC_irrigfrac(n, fgrd) 

! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_irrigationMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
!EOP      

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n))

!
! !DESCRIPTION:
! The Global Rain-Fed, Irrigated, and Paddy Croplands (GRIPC) Dataset
!     (Salmon, 2013)
! 
!   The GRIPC dataset is based on MODIS and ancillary datasets and maps
!  three classes of land use:  rain-fed, irrigated and paddy-type crop areas.
!  The non-crop areas are identified useing a MODIS-based agricultural mask
!  derived from the Friedl et al. (2010) product.  In addition, a training
!  dataset of hundreds of worldwide insitu locations (~150 original + 210 new ones)
!  are used, which were used in other MODIS land cover product verification studies.
!
!   Additional screening of the MODIS-pixels is carried out using other ancillary
!  maps from Google Earth, the FAO-GMIA irrigation map, and other MODIS
!  vegetation index datasets. In developing this GRIPC classification map,
!  Salmon (2013) applied a univariate (supervised) decision tree classification scheme,
!  which is apparantly used by other land cover dataset development groups.
!  Several other MODIS related vegetation and temperature (like LST) data
!  were used in further identifying areas where water and vegetation sources
!  are present (based on Ozdogan et al., 2010).  Finally, Salmon (2013)
!  merged the final MODIS land-use classified map with crop inventory databases
!  that include spatial crop extents and maxiumum monthly growing areas.
!
!  Please note that the ArcGIS-converted image file, irrigfrac_salmon2013.flt,
!   is North Up - that is, the beginning point is the upper left pixel.
!
!  The legend is:
!    Undefined           = 0
!    Rain-fed croplands  = 1
!    Irrigated croplands = 2
!    Paddy croplands     = 3
!    Not cropped         = 4
!
!  Ref: 
!   Salmon, J. Meghan 2013. Using Satellite Remote Sensing and Hydrologic Modeling
!   To Improve Understanding of Crop Management and Agricultural Water Use At Regional
!   To Global Scales. PhD dissertation, Boston University. Ann Arbor: ProQuest/UMI.
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

! GRIPC crop-irrigation type:
  integer, parameter :: input_cols = 86400
  integer, parameter :: input_rows = 36000
  real,    parameter :: IN_xres = 0.0041666666662926
  real,    parameter :: IN_yres = 0.0041666666662926
  integer, parameter :: num_bins = 4

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
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))       ! Output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))       ! Output logical mask (to match go)
  real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins)    ! Output lis 1d grid
  logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins)    ! Output logical mask (to match go)
  real      :: irrigfrac_cnt(LDT_rc%lnc(n),LDT_rc%lnr(n), num_bins)
  real      :: irrigfrac(LDT_rc%lnc(n),LDT_rc%lnr(n))

  integer   :: noncrop 
  real      :: isum
! ______________________________________________________________

   noncrop = 4

!- Set parameter grid array inputs:
   LDT_irrig_struc(n)%irrig_proj  = "latlon"
   param_gridDesc(1)  = 0.          ! Latlon
   param_gridDesc(2)  = input_cols
   param_gridDesc(3)  = input_rows
   param_gridDesc(4)  = -60.0  + (IN_yres/2) ! LL lat (-59.9999999S)
   param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon (-179.9999999W)
   param_gridDesc(6)  = 128
   param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat
   param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon
   param_gridDesc(9)  = IN_yres     ! dy: 0.004166666
   param_gridDesc(10) = IN_xres     ! dx: 0.004166666
   param_gridDesc(20) = 64
! _____________________________

!- Check if file is present:
   inquire(file=trim(LDT_irrig_struc(n)%irrigfracfile), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Irrigation fraction map ",trim(LDT_irrig_struc(n)%irrigfracfile)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(unit=LDT_logunit,fmt=*) "[INFO] Reading MODIS/GRIPC crop-water source map file"

!- Open file to be processed::
   ftn = LDT_getNextUnitNumber()
   open( ftn, file = LDT_irrig_struc(n)%irrigfracfile,form='unformatted', &
          access='direct', convert='little_endian',status='old', &
          recl=(4*input_cols) )

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_irrig_struc(n)%irrig_proj, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )

! _________

   allocate( var_read(input_cols) )
   allocate( var_in(subpnc,subpnr) )
   var_in = float(noncrop)

!- Reverse-Y and READ in values:
   nrec = 0
   do nr = subpnr, 1, -1
      nrec = input_rows - lat_line(1,nr) + 1
      read(ftn,rec=nrec) var_read   ! Read in GRIPC pixels
      do nc = 1, subpnc
         var_in(nc,nr) = var_read(lon_line(nc,1))
      end do
   end do
   deallocate( var_read )
   call LDT_releaseUnitNumber(ftn)

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
   lo1 = .false.;  lo2 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do nr = 1, subpnr
      do nc = 1, subpnc
         i = i + 1
!         if( var_in(nc,nr) == 2. .or. &  ! Irrigated type
!             var_in(nc,nr) == 3. ) then  ! Rice-paddy area
         if( var_in(nc,nr) > 0. ) then  
            gi(i) = var_in(nc,nr)  
            li(i) = .true.
         endif 
      enddo
   enddo
   deallocate( var_in )

!- Create mapping between parameter domain and LIS grid domain:
   call upscaleByAveraging_input( subparam_gridDesc, &
                           LDT_rc%gridDesc(n,:), mi, mo, n11 )

!- Calculate total counts for inland water type in each coarse gridcell:
!   call upscaleByCnt( mi, mo, num_bins, noncrop, n11, li, gi, &
   call upscaleByCnt( mi, mo, num_bins, LDT_rc%udef, n11, li, gi, &
                      lo2, go2 )

   irrigfrac_cnt = 0.
   i = 0
   do nr = 1, LDT_rc%lnr(n)
      do nc = 1, LDT_rc%lnc(n)
         i = i + 1
         do t = 1, num_bins
            irrigfrac_cnt(nc,nr,t) = go2(i,t)
         end do
      enddo
   enddo

!- Convert 3-D irrigration fraction to 2-D:
   irrigfrac = 0.
   do nr = 1, LDT_rc%lnr(n)
      do nc = 1, LDT_rc%lnc(n)

       ! Calculate gridpoint column total:
         isum = sum(irrigfrac_cnt(nc,nr,1:num_bins))
!         write(500,*) nc,nr, isum, irrigfrac_cnt(nc,nr,2),irrigfrac_cnt(nc,nr,3)

       ! Estimate 2-D fraction for just irrigation:
         if( isum > 0. ) then
          ! Irrigated layer "sprinkler" only:
            if( LDT_irrig_struc(n)%irrigfrac_typeopt == "sprinkler" ) then
              irrigfrac(nc,nr) = (irrigfrac_cnt(nc,nr,2) / isum) * 100.   
          ! Rice paddy layer-only:
            elseif ( LDT_irrig_struc(n)%irrigfrac_typeopt == "paddy" ) then
              irrigfrac(nc,nr) = (irrigfrac_cnt(nc,nr,3) / isum) * 100.   
          ! Irrigated + rice paddy layers:
            elseif ( LDT_irrig_struc(n)%irrigfrac_typeopt == "sprinkler+paddy" ) then
              irrigfrac(nc,nr) = ((irrigfrac_cnt(nc,nr,2)+irrigfrac_cnt(nc,nr,3)) / isum) * 100.
            endif
         else
            irrigfrac(nc,nr) = 0.
         endif
 
      enddo
   end do

   fgrd = irrigfrac

  write(LDT_logunit,fmt=*) "[INFO] Done reading MODIS/GRIPC irrigation gridcell fractions"

end subroutine read_GRIPC_irrigfrac

