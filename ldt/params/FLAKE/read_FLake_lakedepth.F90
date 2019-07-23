!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_FLake_lakedepth
!  \label{read_FLake_lakedepth}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  16 Apr 2013: K. Arsenault; Modified for lake maps
!
! !INTERFACE:
 subroutine read_FLake_lakedepth(n, num_bins, lakedepth, &
                                 lakefrac, lakedepthQC )

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
         LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod,     only : LDT_transform_paramgrid
  use LDT_paramTileInputMod, only: param_1dbin_areacalc, &
                                   param_index_fgrdcalc
  
  use FLAKE_parmsMod

!EOP      
  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_bins
  real, intent(inout) :: lakedepth(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real, intent(inout) :: lakefrac(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real, intent(inout), optional :: lakedepthQC(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
!
! !DESCRIPTION:
!  This subroutine retrieves the lake depth for each gridcell
!   and returns the values in a latlon projection.
!
!  The Global database provides the external parameters fields for the
!  parameterization of lakes in atmospheric modeling.  It combines depth
!  information for the individual lakes from different sources with a map.
!  For mapping, the raster map of ECOCLIMAP2 dataset for ecosystems was used.
!  For some large lakes the bathymetry is included.  Additionally, the software
!  to project the lake-related information accurately onto an atmospheric
!  model grid is provided.
!
!  The global gridded datasets (in GlobalLake.tar.gz) contain the following
!  information on the geographical grid with the resolution of 30 arc sec.
!  (approx. 1 km):
!
!   1) the distributed mean lake depth values OR bathymetry data, and
!   2) the distributed flagm (QC) map to estimate reliability of the lake depth
!     information in every pixel of the grid:
!
!    = 0 - no inland water,
!    = 1 - the lake was not recognized by the automating mapping software,
!    = 2 - the lake depth value was missing in the dataset for individual lakes,
!    = 3 - the real depth value was used,
!    = 4 - a river.
!
!  * Note:  For flags 1 and 2, a default lake depth value of 10 m is assigned
! ** Note:  For flag 4, 3 m is set for lake depth.
!
!  Ref:  Kourzeneva, E., E. Martin, Y. Batrak and P.L.E. Moigne, 2012:
!        Climate data for parameterisation of lakes in Numerical Weather 
!        Prediction models, Tellus A 2012, 64: 17226. DOI: 10.3402/tellusa.v64i0.17226
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of bins (or bands)
!  \item[lakedepth]
!   output field FLake lake depth
!  \item[lakefrac]
!   output grid fractions for lake tiles 
!  \item[lakedepthQC]
!   optional output for lake depth QC output 
!  \end{description}
!
!EOP      
  integer :: ftn1, ftn2
  logical :: file_exists
  logical :: qc_file_present          ! Flag to indicate that QC file is present
  integer :: c, r, t, i, ilon, ilat, ilon2, ilat2  ! loop indexes
  integer :: NorthPix, SouthPix, WestPix, EastPix  ! The coordinates in pixels
  integer :: ErCode

  integer, parameter :: NLonB = 43200, &  ! Number of longitude pixels of the bitmap
                        NLatB = 21600     ! Number of latitude pixels of the bitmap
  integer, parameter :: PixSize = 30      ! Pixel size of the bitmap in seconds of arc
  integer(1), dimension(NLonB) :: LonPix1 ! Pixels of status data along the latitude
  integer(2), dimension(NLonB) :: LonPix2 ! Pixels of depth data along the latitude
  real    :: IN_yres, IN_xres

  integer :: glpnc, glpnr                 ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr               ! Parameter subsetted columns and rows
  real    :: param_gridDesc(20)           ! Input parameter grid desc array
  real    :: subparam_gridDesc(20)        ! Subsetted Input parameter grid desc array
  integer, allocatable  :: lat_line(:,:), lon_line(:,:)
  integer :: mi                           ! Total number of input param grid array points
  integer :: mo                           ! Total number of output LIS grid array points
  integer, allocatable  :: n11(:)         ! array that maps the location of each input grid
                                          ! point in the output grid.
  real,    allocatable  :: gi1(:), gi2(:) ! input parameter 1d grid
  logical*1,allocatable :: li1(:), li2(:)  ! input logical mask (to match gi)
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output lis 1d grid
  real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output lis 1d grid
  real      :: go3(LDT_rc%lnc(n)*LDT_rc%lnr(n),1)
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output logical mask (to match go1)
  logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output logical mask (to match go2)
!  real      :: lo3(LDT_rc%lnc(n)*LDT_rc%lnr(n),1)

  real, dimension(:,:), allocatable       :: ReadDepth
! Depth of lake(s) on the bitmap in the circumscribed rectangular
  integer(1), dimension(:,:), allocatable :: ReadStatus
! Status of lake(s) on the bitmap in the circumscribed rectangular

  real  :: lakedepth2d(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real  :: lakedepthQC2d(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real  :: lakedepthcnt(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real  :: total_cnt(LDT_rc%lnc(n)*LDT_rc%lnr(n))

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

  inquire(file = trim(FLAKE_struc(n)%lakedepthfile), exist=file_exists)
  if(.not. file_exists) then 
     write(LDT_logunit,*) "Lake depth map (",trim(FLAKE_struc(n)%lakedepthfile),") not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  inquire(file = trim(FLAKE_struc(n)%lakedepthQCfile), exist=qc_file_present)
  if(.not. qc_file_present) then
     write(LDT_logunit,*) "Lake depth QC map (",trim(FLAKE_struc(n)%lakedepthQCfile),") not present."
     write(LDT_logunit,*) " No QC applied to lake depth map ..."
  endif

  write(LDT_logunit,*) "[INFO] Reading FLake lake depth files"

  lakedepth = 10.  ! Default lake depth value for FLake model
  lakefrac  = 0.
  if( qc_file_present )  lakedepthQC = 0

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, projection, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )

!  WestPix: CALL Coor2Num(West, 1, WestPix, ErCode)
   WestPix=INT((180.+MIN(subparam_gridDesc(5),179.999))*3600/PixSize)+1  
!  EastPix: CALL Coor2Num(East, 1, EastPix, ErCode)
   EastPix=INT((180.+MIN(subparam_gridDesc(8),179.999))*3600/PixSize)+1    
!  SouthPix: CALL Coor2Num(South, 2, SouthPix, ErCode)
   SouthPix=NLatB-(INT((90.+MIN(subparam_gridDesc(4),89.999))*3600/PixSize)+1)+1  
!  NorthPix: CALL Coor2Num(North, 2, NorthPix, ErCode)
   NorthPix=NLatB-(INT((90.+MIN(subparam_gridDesc(7),89.999))*3600/PixSize)+1)+1 

   allocate( ReadDepth(subpnc,subpnr), stat=ErCode )
    IF(ErCode.NE.0) STOP " -- Can't allocate the array <<ReadDepth>>"
   allocate( ReadStatus(subpnc,subpnr), stat=ErCode )
    IF(ErCode.NE.0) STOP " -- Can't allocate the array <<ReadStatus>>"

   ReadDepth = LDT_rc%udef
   if( qc_file_present )  ReadStatus = -9

! -------------------------------------------------------------------
!    READ IN LAKE DEPTH AND QC LAKE DEPTH FILES/INFO
! -------------------------------------------------------------------

   ftn1 = LDT_getNextUnitNumber()
   if( qc_file_present )  ftn2 = LDT_getNextUnitNumber()

   open(ftn1, file=FLAKE_struc(n)%lakedepthfile, &
         form='unformatted', access='direct',&
         convert="little_endian", recl=NLonB*2)

   if( qc_file_present ) then
      open(ftn2, file=FLAKE_struc(n)%lakedepthQCfile, &
           form='unformatted', access='direct',&
           convert="little_endian", recl=NLonB)
   endif

   ilat2 = 0
   do ilat = SouthPix, NorthPix, -1
      ilat2 = ilat2 + 1
    ! Read actual lake depth file:
      read(ftn1,REC=ilat) LonPix2
    ! Read lake depth QC-file:
      if( qc_file_present ) then
        read(ftn2,REC=ilat) LonPix1
      endif

      ilon2 = 0
      do ilon = WestPix, EastPix
         ilon2 = ilon2 + 1
       ! Read-in Lake Depth 2-D Array:
         ReadDepth(ilon2,ilat2)=LonPix2(ilon)/10.
       ! Read-in QC Lake Depth 2-D Array:
         if( qc_file_present ) then
           ReadStatus(ilon2,ilat2)=LonPix1(ilon)
         endif
      end do
   end do

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

   mi = subpnc*subpnr
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
   allocate( gi1(mi), li1(mi), gi2(mi), li2(mi), n11(mi) )
   gi1 = LDT_rc%udef;  gi2 = 0.   ! LDT_rc%udef
   li1 = .false.;  li2 = .false.
   lo1 = .false.; lo2 = .false.;! lo3 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do r = 1, subpnr
      do c = 1, subpnc;  i = i + 1
         gi1(i) = ReadDepth(c,r)
         if( gi1(i) .ne. LDT_rc%udef ) li1(i) = .true.

         if( qc_file_present ) then
           gi2(i) = float(ReadStatus(c,r))
           if( gi2(i) .ne. 0. ) li2(i) = .true.
!           if( gi2(i).ne. LDT_rc%udef ) li2(i) = .true.
         endif
      enddo
   enddo
   deallocate( ReadDepth )
   if( qc_file_present) deallocate( ReadStatus )

!- Create mapping between parameter domain and LIS grid domain:
   call upscaleByAveraging_input( subparam_gridDesc, &
                           LDT_rc%gridDesc(n,:), mi, mo, n11 )


!- Transform parameter grid to LIS run domain:
   select case ( FLAKE_struc(n)%lakeparms_gridtransform )

 !- Transforming 2-D lake depth field: 
!    case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )
    case( "none", "neighbor", "average" )

   !- Transform parameter from original grid to LIS output grid:
      lo1 = .false.
      call LDT_transform_paramgrid(n, FLAKE_struc(n)%lakeparms_gridtransform, &
               subparam_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )

   !- Convert 1D to 2D grid arrays:
      lakedepth2d = LDT_rc%udef
      i = 0
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            i = i + 1
            lakedepth2d(c,r) = go1(i)
          enddo
       enddo

    !- Estimate lake fraction from lake depth map:
       go3 = 0.
       total_cnt = 0.
       do i = 1, mi
          if( li1(i) ) then
         !- Count depths:
            if( n11(i).ne.0 .and. gi1(i) > 0 ) then
               total_cnt(n11(i)) = total_cnt(n11(i)) + 1.
               go3(n11(i),1) = go3(n11(i),1) + 1.0
            endif
          endif
       enddo
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             if( total_cnt(i) .ne. 0  ) then
                lakefrac(c,r,1) = go3(i,1)/total_cnt(i)
             else
                lakefrac(c,r,1) = 0.
             endif 
          enddo
       enddo

 !- 3D Tile Case:
    case( "tile" ) 

      write(*,*)  " FLake Lake Depth 'tile' option is not available at this time ... "
      write(*,*)  " Stopping ... (please select 'average' for the time being)"
      write(*,*)
      call LDT_endrun

   !- Create mapping between parameter domain and LIS grid domain:
      call param_1dbin_areacalc( n, num_bins, mi, mo, n11, &
                                 0., gi1, lakefrac, lakedepth )

 
 !- All Other Cases:
    case default
       write(LDT_logunit,*) " This lake depth spatial transform, ",&
             trim(FLAKE_struc(n)%lakeparms_gridtransform),", is not available at this time ..."
       write(LDT_logunit,*) " Program stopping ..."
       call LDT_endrun

   end select
   deallocate( gi1, li1 )


!- Apply QC file to 2-D lake depth map:
   if( qc_file_present ) then

  !- Select primary QC value (per gridcell): 
     call upscaleByMode( mi, mo, LDT_rc%udef, n11, li2, gi2, &
                         lo2, go2 )

  !- Convert 1D to 2D grid arrays:
     lakedepthQC2d = LDT_rc%udef
     i = 0
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           i = i + 1
           if( go2(i) == LDT_rc%udef) go2(i) = 0.
           lakedepthQC2d(c,r) = go2(i)    ! Real lake depth used

        !- Using QC values can modify the final output lake depths: 
           if( lakedepthQC2d(c,r) == 1 ) lakedepth2d(c,r) = 10.0  ! Lake not orig. recognized; filled 
           if( lakedepthQC2d(c,r) == 2 ) lakedepth2d(c,r) = 10.0  ! Lake depth value missing; filled 
           if( lakedepthQC2d(c,r) == 4 ) lakedepth2d(c,r) =  3.0  ! River points; filled
        enddo
     enddo
   endif
   deallocate( gi2, li2 )


!- Bring 2-D Array to 3-D lake depth tile space:
   if( FLAKE_struc(n)%lakeparms_gridtransform == "none"     .or. &
       FLAKE_struc(n)%lakeparms_gridtransform == "neighbor" .or. &
       FLAKE_struc(n)%lakeparms_gridtransform == "average" ) then
      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
         !- Single lake layer, write to first bin:
!            lakefrac(c,r,1)  = 1.0
            lakedepth(c,r,1) = lakedepth2d(c,r)
            if( qc_file_present ) then
               lakedepthQC(c,r,1)= lakedepthQC2d(c,r)
            endif
         enddo
      enddo
   end if

  call LDT_releaseUnitNumber(ftn1)
  if( qc_file_present ) then
     call LDT_releaseUnitNumber(ftn2)
  endif

  write(LDT_logunit,*) "[INFO] Done reading FLake lake depth files."


end subroutine read_FLake_lakedepth
