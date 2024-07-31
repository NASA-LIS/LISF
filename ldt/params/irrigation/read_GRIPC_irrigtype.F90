!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_GRIPC_irrigtype
!  \label{read_GRIPC_irrigtype}
!
! !REVISION HISTORY:
!  30 Sep 2013: K. Arsenault;  Added new irrigation type map
!
! !INTERFACE:
 subroutine read_GRIPC_irrigtype(n,fgrd,num_bins) 

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
  integer, intent(in) :: num_bins
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

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
!  Please note that the ArcGIS-converted image file, irrigtype_salmon2013.flt,
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
!  \item[num_bins]
!   number of bins or tiles
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
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))           ! Output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))           ! Output logical mask (to match go)
  real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins) ! Output lis 1d grid
  logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n), num_bins) ! Output logical mask (to match go)
  real      :: irrigtype_mode(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real      :: irrigtype_cnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

  integer   :: noncrop 
  real      :: isum

! __________________________________________________________

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
   inquire(file=trim(LDT_irrig_struc(n)%irrigtypefile), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Irrigation type map ",trim(LDT_irrig_struc(n)%irrigtypefile)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(unit=LDT_logunit,fmt=*) "[INFO] Reading MODIS/GRIPC crop-water source map file"

!- Open file to be processed::
   ftn = LDT_getNextUnitNumber()
   open ( ftn, file = LDT_irrig_struc(n)%irrigtypefile,form='unformatted',&
          access='direct', convert='little_endian', &
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

!- Reverse-Y and read in values:
   nrec = 0
   do nr = subpnr, 1, -1
      nrec = input_rows - lat_line(1,nr) + 1
      read(ftn,rec=nrec) var_read
      do nc = 1, subpnc
         var_in(nc,nr) = var_read(lon_line(nc,1))
        ! Set 0-undefined values to specified 'non-crop' index value:  
!         if( var_in(nc,nr) == 0. ) var_in(nc,nr) = float(noncrop)
        ! Set 0-undefined values to universal undefined value:  
         if( var_in(nc,nr) <= 0. ) var_in(nc,nr) = LDT_rc%udef
      end do
   end do
   deallocate( var_read )
   call LDT_releaseUnitNumber(ftn)

! -------------------------------------------------------------------

   irrigtype_mode = LDT_rc%udef
   irrigtype_cnt  = 0.
   fgrd           = 0.

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

   mi = subpnc*subpnr
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)

   if( mi .ne. mo .and. LDT_irrig_struc(n)%irrigtype_gridtransform == "none" ) then
      write(LDT_logunit,*) "[ERR] Spatial transform, 'none', is selected, but number of"
      write(LDT_logunit,*) "  input and output points do not match. Select other spatial"
      write(LDT_logunit,*) "  option (e.g., mode, neighbor, tile, etc.)."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   allocate( li(mi), gi(mi) )
   gi  = float(noncrop)
   li  = .false.
   lo1 = .false.;  lo2 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do nr = 1, subpnr
      do nc = 1, subpnc
         i = i + 1
         gi(i) = var_in(nc,nr)   ! Assign read-in 2D array to 1D
!         if( gi(i) .ne. float(noncrop) ) li(i) = .true.
         if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
      enddo
   enddo
   deallocate( var_in )

!- Apply the spatial transform option:
   select case( LDT_irrig_struc(n)%irrigtype_gridtransform )

 !- (a) Single-layer selection:
    case( "none", "mode", "neighbor" )

   !- Transform parameter from original grid to LIS output grid:
      call LDT_transform_paramgrid(n, LDT_irrig_struc(n)%irrigtype_gridtransform, &
                         subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )

   !- Convert 1D count to 2D grid fgrds:
      i = 0
      do nr = 1, LDT_rc%lnr(n)
         do nc = 1, LDT_rc%lnc(n)
            i = i + 1
            irrigtype_mode(nc,nr) = go1(i)
         enddo
      enddo

 !- (b) Estimate TILED irrigation files (irrigtype_cnt):
    case( "tile" )

   !- Calculate total counts for each irrigation type in each coarse gridcell:
      call LDT_transform_paramgrid(n, LDT_irrig_struc(n)%irrigtype_gridtransform, &
               subparam_gridDesc, mi, num_bins, gi, li, mo, go2, lo2 )

   !- Convert 1D irrigtype_cnt to 2D grid fgrds:
      i = 0
      do nr = 1, LDT_rc%lnr(n)
         do nc = 1, LDT_rc%lnc(n)
            i = i + 1
            do t = 1, num_bins
               irrigtype_cnt(nc,nr,t) = go2(i,t)
            end do
         enddo
      enddo

   end select  ! End grid cnt aggregation method
   deallocate( gi, li )

! ____________

!- Estimate Fraction of Grid (fgrid) 3-D space values:
   select case( LDT_irrig_struc(n)%irrigtype_gridtransform )

 !- Bring 2-D Array to 3-D irrig/source type tile space:
    case( "none", "mode", "neighbor" )

       do nr = 1, LDT_rc%lnr(n)
          do nc = 1, LDT_rc%lnc(n)
          !- Assign types of less than 0 and dominant types to 1:
             if( irrigtype_mode(nc,nr) .le. 0 ) then
                irrigtype_mode(nc,nr) = float(noncrop)
             endif 
             if ( (nint(irrigtype_mode(nc,nr)) .ne. noncrop)  .and. &
                  (nint(irrigtype_mode(nc,nr)) .ne. LDT_rc%udef) ) then
                irrigtype_cnt(nc,nr,NINT(irrigtype_mode(nc,nr))) = 1.0
             endif
          enddo
       enddo

     ! Estimate fraction of grid (fgrid) represented by crop-water source type:
       call param_index_fgrdcalc( n, LDT_irrig_struc(n)%irrig_proj, &
            LDT_irrig_struc(n)%irrigtype_gridtransform, &
            noncrop, num_bins, irrigtype_cnt, fgrd )

 !- Tiled case:
    case( "tile" )
      do nr = 1, LDT_rc%lnr(n)
         do nc = 1, LDT_rc%lnc(n)

          ! Calculate gridpoint column total:
            isum = sum(irrigtype_cnt(nc,nr,1:num_bins))

          ! Undefined value case:
            if( isum <= 0. ) then
              fgrd(nc,nr,:) = 0.

          ! Actual summed-value case:
            else
              do t = 1, num_bins
                 if( irrigtype_cnt(nc,nr,t) .ne. LDT_rc%udef ) then
                    fgrd(nc,nr,t) = irrigtype_cnt(nc,nr,t) / isum
                 else
                    print *, "[ERR] An irrigtype count = UDEF :",nc,nr,t,irrigtype_cnt(nc,nr,t)
                    print *, " Program stopping ... "
                    call LDT_endrun
                 endif
              enddo
            endif
         enddo
      enddo

   end select  ! End irrigtype_mode/cnt aggregation method

  write(LDT_logunit,fmt=*) "[INFO] Done reading MODIS/GRIPC crop-water source file"

end subroutine read_GRIPC_irrigtype

