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
! !ROUTINE: read_GLWD_inlandwatertype
!  \label{read_GLWD_inlandwatertype}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  16 Apr 2013: K. Arsenault; Modified for inland water bodies 
!  24 May 2014: K. Arsenault; Updated reader for additional options
!
! !INTERFACE:
subroutine read_GLWD_inlandwatertype(n, num_bins, watertype, waterfrac )

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
         LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod,     only : LDT_transform_paramgrid
  use LDT_paramTileInputMod, only: param_index_fgrdcalc

  use FLAKE_parmsMod
!EOP      
  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_bins
  real, intent(inout) :: watertype(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real, intent(inout) :: waterfrac(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
!
! !DESCRIPTION:
!  This subroutine retrieves the WWF Global Lake and Wetland Database
!  (GLWD) for each gridcell and extracts inland water body type
!  (tile) and related fraction.
!
!  Source:  http://www.worldwildlife.org/science/data/item1877.html
!
!  Description of GLWD-3 data set
!  File name: glwd_3 (folders .glwd_3. and .info., legend .glwd_3.avl.)
!  File size: 26.9 MB (8.4 MB zipped)
!  File format: Grid in ArcView/ArcInfo coverage format
!  Data format: integer values, for coding see legend below
!  Spatial resolution: 30 x 30 second
!  Projection: Geographic, degrees longitude and latitude
!  Spatial domain: Global land area (except Antarctica and glaciated Greenland)
!  Cell value Lake or Wetland Type:
!   1 Lake
!   2 Reservoir
!   3 River
!   4 Freshwater Marsh, Floodplain
!   5 Swamp Forest, Flooded Forest
!   6 Coastal Wetland (incl. Mangrove, Estuary, Delta, Lagoon)
!   7 Pan, Brackish/Saline Wetland
!   8 Bog, Fen, Mire (Peatland)
!   9 Intermittent Wetland/Lake
!  10 50-100% Wetland
!  11 25-50% Wetland
!  12 Wetland Compex (0-25% Wetland)
!
!  Ref: Lehner, B. and P. Döll (2004): Development and validation of a 
!    global database of lakes, reservoirs and wetlands. Journal of Hydrology,
!    296/1-4: 1-22.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of bins (or bands)
!  \item[watertype]
!   output field for inland water body tiles
!  \item[waterfrac]
!   output grid fractions for inland water body types
!  \end{description}
!
!EOP      
  integer :: ftn1
  logical :: file_exists
  integer :: c, r, t, i, ilon, ilat, ilon2, ilat2  ! loop indexes
  integer :: NorthPix, SouthPix, WestPix, EastPix  ! The coordinates in pixels
  integer :: ErCode

  integer, parameter :: NLonB = 43200, &  ! Number of longitude pixels of the bitmap
                        NLatB = 21600     ! Number of latitude pixels of the bitmap
  integer, parameter :: PixSize = 30      ! Pixel size of the bitmap in seconds of arc
  real, dimension(NLonB) :: LonPix        ! Pixels of data along the latitude
  real    :: IN_yres, IN_xres

  integer :: glpnc, glpnr                 ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr               ! Parameter subsetted columns and rows
  real    :: param_gridDesc(20)           ! Input parameter grid desc array
  real    :: subparam_gridDesc(20)        ! Subsetted input parameter grid desc array
  integer, allocatable  :: lat_line(:,:), lon_line(:,:)
  integer :: mi                           ! Total number of input param grid array points
  integer :: mo                           ! Total number of output LIS grid array points
  integer, allocatable  :: n11(:)         ! array that maps the location of each input grid
                                          ! point in the output grid.
  real,    allocatable  :: gi(:)          ! input parameter 1d grid
  logical*1,allocatable :: li(:)          ! input logical mask (to match gi)
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output logical mask (to match go1)
  real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n),num_bins) ! Output lis 1d grid
  logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n),num_bins) ! Output logical mask 

! Interior lake/wetland type:
  real, dimension(:,:), allocatable :: Readtype
  real :: watertype2d(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real :: watypecnt3d(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

  character(50) :: projection

! __________________________________________________________

   projection = "latlon"

   IN_yres = 1.0/120.0
   IN_xres = 1.0/120.0

!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.          ! Latlon
   param_gridDesc(2)  = NLonB       ! input_cols
   param_gridDesc(3)  = NlatB       ! input_rows
   param_gridDesc(4)  = -90.0  + (IN_yres/2) ! LL lat (-89.9960000S)
   param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon (-179.9960000W)
   param_gridDesc(6)  = 128
   param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat (89.99570000N)
   param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon (179.9960000W)
   param_gridDesc(9)  = IN_yres     ! dy: 0.0083333
   param_gridDesc(10) = IN_xres     ! dx: 0.0083333
   param_gridDesc(20) = 64

   write(LDT_logunit,*) "** MAJOR NOTE:  ABOUT THE GLWD INLAND WATERBODY TYPE MAP ..."
   write(LDT_logunit,*) "** ---------- "
   write(LDT_logunit,*) "   THIS MAP IS NOT CURRENTLY SET UP TO RUN WITH ANY"
   write(LDT_logunit,*) "   LAKE, WETLAND, OR RELATED INTERIOR WATER BODY TYPE."
   write(LDT_logunit,*) "   FUTURE EXPANSIONS OF THIS FIELD'S USE ARE COMING ..."

   inquire(file = trim(FLAKE_struc(n)%inlandwaterfile), exist=file_exists)
   if(.not. file_exists) then 
      write(LDT_logunit,*) "Inland waterbody type map (",&
              trim(FLAKE_struc(n)%inlandwaterfile),") not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(LDT_logunit,*) "[INFO] Reading in GLWD Inland Water Bodies File"

   watertype = LDT_rc%udef
   waterfrac  = 0.

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, projection, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )


    WestPix=INT((180.+MIN(subparam_gridDesc(5),179.999))*3600/PixSize)+1  
    EastPix=INT((180.+MIN(subparam_gridDesc(8),179.999))*3600/PixSize)+1    
    SouthPix=NLatB-(INT((90.+MIN(subparam_gridDesc(4),89.999))*3600/PixSize)+1)+1  
    NorthPix=NLatB-(INT((90.+MIN(subparam_gridDesc(7),89.999))*3600/PixSize)+1)+1 

   allocate( Readtype(subpnc,subpnr), stat=ErCode )
   if( ErCode.ne.0 ) STOP " -- Can't allocate the array <<Readtype>>"
   Readtype = LDT_rc%udef

!- Open and read water type files:
   ftn1 = LDT_getNextUnitNumber()
   open( ftn1, file=FLAKE_struc(n)%inlandwaterfile, &
         form='unformatted', access='direct', convert="little_endian", recl=NLonB*4 )

   ilat2 = 0
   do ilat = SouthPix, NorthPix, -1
      ilat2 = ilat2 + 1
      read(ftn1,rec=ilat) LonPix

      ilon2 = 0
      do ilon = WestPix, EastPix
         ilon2 = ilon2 + 1
         Readtype(ilon2,ilat2) = LonPix(ilon)
      end do
   end do

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

   mi = subpnc*subpnr
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
   allocate( gi(mi), li(mi), n11(mi) )
   gi = LDT_rc%udef
   li = .false.
   lo1 = .false.; lo2 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do r = 1, subpnr
      do c = 1, subpnc;  i = i + 1
         gi(i) = Readtype(c,r)
         if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
      enddo
   enddo
   deallocate( Readtype )

!- Create mapping between parameter domain and LIS grid domain:
   call upscaleByAveraging_input( subparam_gridDesc, &
                           LDT_rc%gridDesc(n,:), mi, mo, n11 )


!- Transform parameter grid to LIS run domain:
   select case ( FLAKE_struc(n)%inlandwater_gridtransform )

 !- Transforming 2-D lake depth field: 
 !- (a) Estimate NON-TILED dominant inland water body type:
    case( "none", "neighbor", "mode" )

    !- Transforming 2-D lake depth field: 
       lo1 = .false.
       call LDT_transform_paramgrid(n, FLAKE_struc(n)%inlandwater_gridtransform, &
                          subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )

    !- Convert 1D to 2D grid arrays:
       watertype2d(:,:) = LDT_rc%udef
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             watertype2d(c,r) = go1(i)
          enddo
       enddo

 !- (b) Estimate TILED inland water body types:
    case( "tile" )

    !- Calculate total counts for inland water type in each coarse gridcell:
       call upscaleByCnt( mi, mo, num_bins, LDT_rc%udef, n11, li, gi, &
                          lo2, go2 )

    !- Convert 1D inland water type to 2D grid arrays:
       watypecnt3d = 0.
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             do t = 1, num_bins
                watypecnt3d(c,r,t) = go2(i,t)
             end do
          enddo
       enddo

  ! Other spatial transforms selected:
    case default
      write(LDT_logunit,*) " This spatial transform, ",&
            trim(FLAKE_struc(n)%inlandwater_gridtransform),", for water body types is "
      write(LDT_logunit,*) "  not available at this time."
      write(LDT_logunit,*) "  Please select:  neighbor, mode or tile "
      write(LDT_logunit,*) " Program stopping ..."
      call LDT_endrun

    end select  ! End inland water type aggregation method
    deallocate( gi, li, n11 )


! ........................................................................

!- Bring 2-D inland water type to 3-D water type tile space:
   if ( FLAKE_struc(n)%inlandwater_gridtransform == "none"     .or. &
        FLAKE_struc(n)%inlandwater_gridtransform == "neighbor" .or. &
        FLAKE_struc(n)%inlandwater_gridtransform == "mode" ) then  ! -- NON-TILED SURFACES
      watypecnt3d = 0.
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            if ( watertype2d(c,r) .le. 0 ) then
               watertype2d(c,r) = LDT_rc%udef
            endif
            if ( nint(watertype2d(c,r)) .ne. LDT_rc%udef ) then
               watypecnt3d(c,r,nint(watertype2d(c,r))) = 1.0
            endif
         enddo
      end do
   endif   ! End NON-TILED inland water type option


!- Estimate fraction of grid (waterfrac) represented by inland water type::
   call param_index_fgrdcalc( n, projection,    &
                    FLAKE_struc(n)%inlandwater_gridtransform, &
                    nint(LDT_rc%udef), &
                    num_bins, watypecnt3d, waterfrac )

   watertype = watypecnt3d 

   call LDT_releaseUnitNumber(ftn1)

   write(LDT_logunit,*) "[INFO] Done reading GLWD lake/wetland type file"


end subroutine read_GLWD_inlandwatertype

