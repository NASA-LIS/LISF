!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_AQUASTAT_irrigtype
!  \label{read_AQUASTAT_irrigtype}
!
! !REVISION HISTORY:
!  23 May 2019: H. Beaudoing;  Adopted the routines and dataset compiled by
!                              Sarith Mahanama.
!  01 Apr 2021: H. Beaudoing;  Added US county data of USGS to overlay.
!                              The county data and reading section are provided
!                              by Sarith.
!
! !INTERFACE:
#include "LDT_misc.h"
 subroutine read_AQUASTAT_irrigtype(n,fgrd,num_bins) 

! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use LDT_paramDataMod
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_irrigationMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

!EOP      

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_bins
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

!
! !DESCRIPTION:
!  Irrigation type information was compiled from USGS and AQUASTAT census data.
!  -US irrigated areas (1000 acre) by type (drip, flood and sprinkler) for 50 states, in 2015 (Dieter et al., 2018)
!  -Irrigated areas (1000 ha) by type (drip, flood, and sprinkler) from about 65 countries reported for year 2000 or later (AQUASTAT data base).
! 
!  The irrigation type data is used to estimate dominant of three methods to 
!  create a map.  Flood type is assumed to be the most common world-wide and
!  assigned to countries without report. 
!  
!  This routine maps out country and state/county table values onto geographical
!  2D array using the GIS political boundary data (native at 1 km) and
!  US Ceneus Boundary 2015 data, respectively.
!  The grid trensformation method applies to the handling of the political data,
!  so only none, mode, or neighbor are supported. The country and county index
!  are included in the output to be used in the distribution of irrigation 
!  type onto crop tiles in LIS irrigation. 
!
!  The irrigation type output contains values over land regardless of the grid
!  cell being irrigated or cropland land cover.
!  Note the irrigtype defined with this routine is actual irrigation method,
!  which is different from those with GRIPC dataset.
!
!  The legend is:
!    Undefined           = 0
!    Sprinkler           = 1
!    Drip                = 2
!    Flood               = 3
!  
!  Ref: 
!  Dieter, C.A., Maupin, M.A., Caldwell, R.R., Harris, M.A., Ivahnenko, T.I., 
!  Lovelace, J.K., Barber, N.L., and Linsey, K.S., 2018, Estimated use of water
!  in the United States in 2015: U.S. Geological Survey Circular 1441, 65 p.,
!  https://doi.org/10.3133/cir1441. [Supersedes USGS Open-File Report 2017–1131
!
!  AQUASTAT database - FAO, http://www.fao.org/nr/water/aquastat/data/query/index.html?lang=en 
!
!  USGS Water Use Data for each state to obtain county data accessed via 
!  National Water Information System (NWIS) at https://waterdata.usgs.gov 
!  US census boundary data was downloaded from
!  https://www2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_us_county_20m.zip
!
!  Global irrigation type data of AQUASTAT was supplemented with various
!  sources for individual country compiled by Matt Rodell.
!  USDA's county level data was obtained from www.nass.usda.gov.  
!
!  USDA conducted the Census of Agriculture in 2017, followed by Irrigation and !  Water Management Survey (IWMS) conducted in 2018. These information serve as
!  basis for the algorithm assigning irrigation type per crop.  Specifically,
!  we used Table 36: state data for Pressure (sprinkler+drip) vs Gravity
!  flood) method for selected crops, and Table 38: percentage of methods
!  within Pressure vs Gravity for selected crops above at US level.
!  Tables were downloaded from https://www.nass.usda.gov/Publications/AgCensus/2017/Online_Resources/Farm_and_Ranch_Irrigation_Survey/index.php.
!
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of bins or tiles
!  \item[fgrd]
!   output field AQUASTAT irrigation type
!  \end{description}
!
!EOP      

! AQUASTAT crop-irrigation type:
  integer, parameter :: noncrop = 0
  integer, parameter :: input_cols = 43200
  integer, parameter :: input_rows = 21600
  real, parameter    :: IN_xres = 0.008333333
  real, parameter    :: IN_yres = 0.008333333
  integer, parameter :: N_GADM = 256
  integer, parameter :: N_STATES = 50
  integer, parameter :: N_COUNTRIES = N_GADM
  integer*2, allocatable, dimension (:,:,:)      :: cnt_code
  character(len=3), dimension(N_COUNTRIES,3)     :: ABR
  character(len=48), dimension(N_STATES,48)      :: state_name
  integer*2, dimension(N_COUNTRIES)              :: loc_index
! USGS county data variables:
  integer*2, allocatable, dimension(:,:)         :: fid
  integer*2, allocatable, dimension(:,:)         :: temp
  integer                                        :: input_rows2
  integer                                        :: county
  integer                                        :: states
  integer                                        :: n_countyUS
  real, allocatable, dimension(:,:)              :: sprinklerfr
  real, allocatable, dimension(:,:)              :: dripfr
  real, allocatable, dimension(:,:)              :: floodfr
  integer, allocatable, dimension(:)             :: geoid
  integer                                        :: SS,CCC

  integer                             :: i,j,ii,jj
  integer                             :: ncid,ncstatus,varid
  integer                             :: countyid,stateid,uscountyid
  integer                             :: cindex
  integer                             :: N_METHOD ! # of avail country data
  character(len=2),dimension(N_STATES):: ST_NAME
  character(len=3),allocatable        :: CNT_ABR(:)
  real, allocatable                   :: sprink(:), drip(:), flood(:), tarea(:)
  real, dimension (N_STATES)          :: us_sprink, us_drip, us_flood, us_tarea
  integer*2, allocatable              :: lis_cnt_mask(:,:)

  integer   :: nc, nr
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

  real, allocatable :: var_in(:,:,:)
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))           ! Output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))           ! Output logical mask (to match go)

  character(len=120) :: usfile, glbfile, polfile, countyfile

! -------------------------------------------------------------------
! Irrigation type area or fraction inputs
! -------------------------------------------------------------------
!- Open report files to be processed::
!!! irrigtypefile in ldt.config point to the directory where files
!!! are found---file names are Hard Coded here!!!
!!! Latest update: change Global_IMethod.data -> Global_IMethod.txt
   glbfile = trim(LDT_irrig_struc(n)%irrigtypefile)//'Global_IMethod.txt'
   usfile = trim(LDT_irrig_struc(n)%irrigtypefile)//'US_IMethod.2015'
   polfile = trim(LDT_irrig_struc(n)%irrigtypefile)//'GADM_Country_and_USStates_codes_1km.nc4'
   countyfile = trim(LDT_irrig_struc(n)%irrigtypefile)//'cb_2015_us_county_30arcsec.nc4'   

   ftn = LDT_getNextUnitNumber()
   open (ftn, file=trim(glbfile), form = 'formatted', status = 'old')
    READ (ftn, *) N_METHOD

    allocate (CNT_ABR (1:N_METHOD))
    allocate (sprink  (1:N_METHOD))
    allocate (drip    (1:N_METHOD))
    allocate (flood   (1:N_METHOD))
    allocate (tarea   (1:N_METHOD))

    do i = 1, N_METHOD
       read (ftn, *) CNT_ABR (i),sprink(i), drip(i), flood(i)
    end do
   call LDT_releaseUnitNumber(ftn)

   ftn = LDT_getNextUnitNumber()
   open (ftn, file=trim(usfile), form = 'formatted', status = 'old')
    do i = 1, N_STATES
       read (ftn, *) ST_NAME (i),us_sprink(i), us_drip(i), us_flood(i)
    end do
   call LDT_releaseUnitNumber(ftn)

   tarea  = sprink + drip + flood
   sprink = sprink / tarea
   drip   = drip   / tarea
   flood  = flood  / tarea

   us_tarea  = us_sprink + us_drip + us_flood
   us_sprink = us_sprink / us_tarea
   us_drip   = us_drip   / us_tarea
   us_flood  = us_flood  / us_tarea
 
   ncstatus = nf90_open(trim(countyfile), NF90_NOWRITE, ncid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot open ", trim(polfile)
   ncstatus = nf90_inq_dimid(ncid,"county",countyid)
   ncstatus = nf90_inquire_dimension(ncid,countyid,len=county)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] failed inqure dimension for county"
   ncstatus = nf90_inq_dimid(ncid,"states",stateid)
   ncstatus = nf90_inquire_dimension(ncid,stateid,len=states)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] failed inqure dimension for states"
   ncstatus = nf90_inq_dimid(ncid,"n_countyUS",uscountyid)
   ncstatus = nf90_inquire_dimension(ncid,uscountyid,len=n_countyUS)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] failed inqure dimension for n_countyUS"
   print*, "county data dimensions: ",county, states, n_countyUS
   allocate(sprinklerfr(county,states))
   allocate(dripfr(county,states))
   allocate(floodfr(county,states))
   allocate(geoid(n_countyUS))
   ncstatus = nf90_inq_varid(ncid, "SPRINKLERFR", varid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] dataset SPRINKLERFR not exist in ",  trim(countyfile)
   ncstatus = nf90_get_var(ncid, varid, sprinklerfr, start=(/1,1/), count=(/county,states/))
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot access SPRINKLERFR in ",  trim(countyfile)
   ncstatus = nf90_inq_varid(ncid, "DRIPFR", varid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] dataset DRIPFR not exist in ",  trim(countyfile)
   ncstatus = nf90_get_var(ncid, varid, dripfr, start=(/1,1/), count=(/county,states/))
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot access DRIPFR in ",  trim(countyfile)
   ncstatus = nf90_inq_varid(ncid, "FLOODFR", varid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] dataset FLOODFR not exist in ",  trim(countyfile)
   ncstatus = nf90_get_var(ncid, varid, floodfr, start=(/1,1/), count=(/county,states/))
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot access FLOODFR in ",  trim(countyfile)
   ncstatus = nf90_inq_varid(ncid, "GEOID", varid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] dataset GEOID not exist in ",  trim(countyfile)
   ncstatus = nf90_get_var(ncid, varid, geoid, start=(/1/), count=(/n_countyUS/))
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot access GEOID in ",  trim(countyfile)
   ncstatus = nf90_close(ncid)

! -------------------------------------------------------------------
! Country, State, and County spatial data
! -------------------------------------------------------------------
   allocate (lis_cnt_mask (LDT_rc%lnc(n),LDT_rc%lnr(n)))

!- Set parameter grid array inputs:
   LDT_irrig_struc(n)%irrig_proj  = "latlon"
   param_gridDesc(1)  = 0.          ! Latlon
   param_gridDesc(2)  = input_cols
   param_gridDesc(3)  = input_rows
   param_gridDesc(4)  = -90.0  + (IN_yres/2) ! LL lat (-89.99583)
   param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon (-179.9958)
   param_gridDesc(6)  = 128
   param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat (89.99583)
   param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon (179.9959)
   param_gridDesc(9)  = IN_yres     ! dy: 0.008333333
   param_gridDesc(10) = IN_xres     ! dx: 0.008333333
   param_gridDesc(20) = 64

!- Check if file is present:
   inquire(file=trim(polfile), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Irrigation type political map ",trim(polfile)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(unit=LDT_logunit,fmt=*) "[INFO] Reading AQUASTAT crop-water source map file"

!- Open file to be processed::
!- Country codes are stored as 16-bit ingeter 0-999 (999 is missing)
   allocate( cnt_code(input_cols,input_rows,2) )
   ncstatus = nf90_open(trim(polfile), NF90_NOWRITE, ncid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot open ", trim(polfile)
   ncstatus = nf90_inq_varid(ncid, "COUNTRY_ABR", varid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] dataset COUNTRY_ABR not exist in ",  trim(polfile)
   ncstatus = nf90_get_var(ncid, varid, ABR)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot access COUNTRY_ABR in ",  trim(polfile)

   ncstatus = nf90_inq_varid(ncid, "COUNTRY_INDEX", varid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] dataset COUNTRY_INDEX not exist in ",  trim(polfile)
   ncstatus = nf90_get_var(ncid, varid, loc_index)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot access COUNTRY_INDEX in ",  trim(polfile)

   ncstatus = nf90_inq_varid(ncid, "UNIT_CODE", varid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] dataset UNIT_CODE not exist in ",  trim(polfile)
   ncstatus = nf90_get_var(ncid, varid, cnt_code)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot access COUNTRY_INDEX in ",  trim(polfile)

   ncstatus = nf90_inq_varid(ncid, "STATE_NAME", varid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] dataset STATE_NAME not exist in ",  trim(polfile)
   ncstatus = nf90_get_var(ncid, varid, state_name)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot access STATE_NAME in ",  trim(polfile)
   ncstatus = nf90_close(ncid)

!- County data to overlay:
   inquire(file=trim(countyfile), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Irrigation type US county map ",trim(countyfile)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(unit=LDT_logunit,fmt=*) "[INFO] Reading USGS crop-water source map file"

!- Open file to be processed::
!- County 3220 over US
   allocate( fid(input_cols,input_rows) )
   fid = LDT_rc%udef
   input_rows2 = input_rows / 2  ! 0-90N latitude
   ncstatus = nf90_open(trim(countyfile), NF90_NOWRITE, ncid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot open ", trim(polfile)
   ncstatus = nf90_inq_varid(ncid, "POLYID", varid)
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] dataset POLYID not exist in ",  trim(countyfile)
   ! read north-south input 
   do j = 1, input_rows2
   ncstatus = nf90_get_var(ncid, varid, fid(:,input_rows-j+1), start=(/1,j/), &
               count=(/input_cols,1/))
   if (ncstatus /= nf90_noerr) write(LDT_logunit,*) "[INFO] cannot access FID in ",  trim(countyfile)
   end do
   ncstatus = nf90_close(ncid)

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------

!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_irrig_struc(n)%irrig_proj, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )
!   print*,'AQUASTAT:',glpnc, glpnr, subpnc, subpnr, LDT_rc%lnc(n), LDT_rc%lnr(n)
!   print*,'subparam_gridDesc:',subparam_gridDesc
! _________
   allocate( var_in(subpnc,subpnr,3) )
   var_in = float(noncrop)

!- Reverse Y-axis and subset region
   do nr = subpnr, 1, -1
      do nc = 1, subpnc
         if ( cnt_code(lon_line(nc,nr),lat_line(nc,nr),2) .ne. 999 )then
          var_in(nc,nr,2) = float(cnt_code(lon_line(nc,nr),lat_line(nc,nr),2))
         else
          var_in(nc,nr,2) = LDT_rc%udef
         endif
         if ( cnt_code(lon_line(nc,nr),lat_line(nc,nr),1) .ne. 999 )then
          var_in(nc,nr,1) = float(cnt_code(lon_line(nc,nr),lat_line(nc,nr),1))
         else
          var_in(nc,nr,1) = LDT_rc%udef
         endif
         if ( fid(lon_line(nc,nr),lat_line(nc,nr)) .ne. -9999 )then
          var_in(nc,nr,3) = float(fid(lon_line(nc,nr),lat_line(nc,nr)))
         else
          var_in(nc,nr,3) = LDT_rc%udef
         endif
      end do
   end do

! -------------------------------------------------------------------
   fgrd           = 0.
   where (LDT_LSMparam_struc(n)%landmask%value(:,:,1) > 0.5)
       fgrd(:,:,3) = 1.   ! default is flood
   end where
! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

   mi = subpnc*subpnr
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)

   if( mi .ne. mo .and. LDT_irrig_struc(n)%irrigtype_gridtransform == "none" ) then
      write(LDT_logunit,*) "[ERR] Spatial transform, 'none', is selected, but number of"
      write(LDT_logunit,*) "  input and output points do not match. "
      write(LDT_logunit,*) "  Select  mode, or neighbor."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   allocate( li(mi), gi(mi), n11(mi) )

!- Loop over 1=country, 2=state, and 3=county 
   do j = 1, 3   

   gi  = float(noncrop)
   li  = .false.
   lo1 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do nr = 1, subpnr
      do nc = 1, subpnc
         i = i + 1
         gi(i) = var_in(nc,nr,j)   ! Assign read-in 2D array to 1D
         if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
      enddo
   enddo

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
            lis_cnt_mask(nc,nr) = go1(i)
         enddo
      enddo

    case default
      write(LDT_logunit,*) "[ERR] Options other than mode, neighbor, or none"
      write(LDT_logunit,*) "are not supported with AQUASTAT GADM data"
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
 
   end select  ! End grid cnt aggregation method

! ____________
!- Map country_code -> type: flood, drip,sprink, & noncrop
   if ( j .eq. 1 ) then
    cnt_loop : do i = 1,N_METHOD
        if(tarea (i) /= 100.) then
         do ii = 1, N_COUNTRIES
            if ( CNT_ABR(i) .eq. ABR(ii,1)) then
             cindex = ii
             exit
            endif 
         enddo
         do nr = 1, LDT_rc%lnr(n)
           do nc = 1, LDT_rc%lnc(n)
             if ((lis_cnt_mask(nc,nr) .eq. loc_index(cindex)).and. &
                 (LDT_LSMparam_struc(n)%landmask%value(nc,nr,1) > 0.5)) then
                 fgrd(nc,nr,3) = flood(i)
                 fgrd(nc,nr,2) = drip (i)
                 fgrd(nc,nr,1) = sprink(i)
             end if
           enddo
         enddo
        endif  ! tarea
    end do cnt_loop
    do nr = 1, LDT_rc%lnr(n)
       do nc = 1, LDT_rc%lnc(n)
         if (LDT_LSMparam_struc(n)%landmask%value(nc,nr,1) > 0.5) then
          LDT_irrig_struc(n)%country%value(nc,nr,1) = float(lis_cnt_mask(nc,nr))
         else
          LDT_irrig_struc(n)%country%value(nc,nr,1) = LDT_rc%udef 
         endif
       enddo
    enddo
   elseif ( j .eq. 2 ) then
   ! Overlay USA irrig methods
    st_loop : do i = 1,N_STATES
         do nr = 1, LDT_rc%lnr(n)
           do nc = 1, LDT_rc%lnc(n)
             if ((lis_cnt_mask(nc,nr) .eq. i).and. &
                 (LDT_LSMparam_struc(n)%landmask%value(nc,nr,1) > 0.5)) then
                 fgrd(nc,nr,3) = us_flood(i)
                 fgrd(nc,nr,2) = us_drip (i)
                 fgrd(nc,nr,1) = us_sprink(i)
             end if
           enddo
         enddo
    end do st_loop
   elseif ( j .eq. 3 ) then
   ! Overlay US County irrig methods
         do nr = 1, LDT_rc%lnr(n)
           do nc = 1, LDT_rc%lnc(n)
             if ((lis_cnt_mask(nc,nr) .ne. LDT_rc%udef).and. &
                 (LDT_LSMparam_struc(n)%landmask%value(nc,nr,1) > 0.5)) then
                 SS = geoid(lis_cnt_mask(nc,nr)) / 1000
                 CCC= geoid(lis_cnt_mask(nc,nr)) - SS*1000
                 fgrd(nc,nr,1) = sprinklerfr(CCC,SS)
                 fgrd(nc,nr,2) = dripfr(CCC,SS)
                 fgrd(nc,nr,3) = floodfr(CCC,SS)
                 LDT_irrig_struc(n)%county%value(nc,nr,1) = float(geoid(lis_cnt_mask(nc,nr)))
             end if
           enddo
         enddo
   endif  ! j == 1, 2, or 3

   end do   ! j 

   deallocate ( gi, li, n11 )
   deallocate ( var_in )
   deallocate (lis_cnt_mask)
   deallocate (cnt_code)
   deallocate (CNT_ABR)
   deallocate (sprink)
   deallocate (drip)
   deallocate (flood)
   deallocate (tarea)
   deallocate (fid)
   deallocate (sprinklerfr)
   deallocate (dripfr)
   deallocate (floodfr)
   deallocate (geoid)
! __________________________________________________________

   write(LDT_logunit,fmt=*) "[INFO] Done reading AQUASTAT crop-water source file"

end subroutine read_AQUASTAT_irrigtype

