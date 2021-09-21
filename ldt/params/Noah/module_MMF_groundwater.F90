!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: MMF_groundwater
! \label{MMF_groundwater}
!
! !REVISION HISTORY:
!  27 Aug 2021: Sarith Mahanama
!     with inputs from David Mocko, Zhuo Wang ; Initial Specification
!
! !INTERFACE:

#include "LDT_misc.h"
module  MMF_groundwater
  
  use LDT_logMod,         only : LDT_logunit, LDT_endrun
  use LDT_coreMod,        only : LDT_rc, LDT_domain
  use LDT_gridmappingMod, only : LDT_RunDomainPts
   
  implicit none
  private

  type, public :: MMF_BCsReader
   
     real, dimension(20)                 :: param_gridDesc, subparam_gridDesc
     integer                             :: glpnc, glpnr, subpnc, subpnr   
     integer,allocatable,dimension (:,:) :: lat_line, lon_line, local_mask
     
     
   contains

     procedure, public :: mi => mmf_init
     procedure, public :: mr => mmf_data_reader
     procedure, public :: cell_area => cell_area_girard
     
     procedure, private:: cell_area_line
     procedure, private:: cell_area_curve
     procedure, private:: cell_area_girard
     
     
  end type MMF_BCsReader
  
  integer,parameter                    :: NY_MMF = 21600
  integer,parameter                    :: NX_MMF = 43200

contains 
  
  SUBROUTINE mmf_init (MBR, nest, project)
    
    implicit none
    
    class (MMF_BCsReader), intent(inout) :: MBR
    integer, intent (in)                 :: nest
    character(*), intent (in)            :: project
    real                                 :: IN_xres, IN_yres
      
    IN_xres    = 360./REAL(NX_MMF)
    IN_yres    = 180./REAL(NY_MMF)
    MBR%param_gridDesc(1)  = 0.          ! Latlon
    MBR%param_gridDesc(2)  = real(NX_MMF)
    MBR%param_gridDesc(3)  = real(NY_MMF)
    MBR%param_gridDesc(4)  = -90.0  + (IN_yres/2) ! LL lat 
    MBR%param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon 
    MBR%param_gridDesc(6)  = 128
    MBR%param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat
    MBR%param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon
    MBR%param_gridDesc(9)  = IN_yres     
    MBR%param_gridDesc(10) = IN_xres     
    MBR%param_gridDesc(20) = 64
    
    ! ------------------------------------------------------------
    !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN BCS DATA
    ! ------------------------------------------------------------
    
    !- Map Parameter Grid Info to LIS Target Grid/Projection Info --
    
    MBR%subparam_gridDesc = 0.

    call LDT_RunDomainPts( nest, project, MBR%param_gridDesc(:),        &
         MBR%glpnc, MBR%glpnr, MBR%subpnc, MBR%subpnr, MBR%subparam_gridDesc, &
         MBR%lat_line, MBR%lon_line)
 
  END SUBROUTINE mmf_init
  
  ! ----------------------------------------------------------------
  
  SUBROUTINE mmf_data_reader (MBR, nest, varname, datadir, lisout)
    
    implicit none

    class (MMF_BCsReader), intent(inout)    :: MBR
    character(*), intent (in)               :: varname, datadir
    real, dimension (:,:,:), intent (inout) :: lisout
    integer, intent (in)                    :: nest
    
    ! Adapted from Zhuo Wang's
    ! /discover/nobackup/projects/lis_aist17/zwang9/Geogrid/geogrid2netcdf.f90
    
    integer, parameter :: iSigned = 1
    integer, parameter :: endian = 0
    integer, parameter :: wordsize = 2
    integer, parameter :: missing = -9999
    integer, parameter :: ntsteps =1
    integer, parameter :: nlat_tile = 1200
    integer, parameter :: nlon_tile = 1200
    character(len=1024):: geogridFile
    character(len=80)  :: tileName
    integer            :: num_tile_lon,num_tile_lat
    integer            :: rc,status,tx,ty,txs,txe,tys,tye    
    real,allocatable   :: tarray(:,:,:), garray(:,:,:)
    real               :: SF
    
    ! c function from Michael G. Duda, NCAR/MMM
    integer,  external :: read_geogrid

    ! set scale factor

    SF = 1.
    select case (trim(varname))

    case ('T')
       ! transmissivity
       SF = 0.01
       
    case ('R')
       ! recharge
       SF = 0.03
       
    case ('E')
       ! river bed elevation
       SF = 0.3
       
    case ('W')
       ! water table depth
       SF = -0.03
       
    case default  
       write(LDT_logunit,*) '[ERROR] Unknown MMF data field : ', trim (varname)
       call LDT_endrun
    end select

    ! Reading GEOGRID data
    ! --------------------
    
    num_tile_lon = NX_MMF / nlon_tile
    num_tile_lat = NY_MMF / nlat_tile

    allocate(tarray(nlon_tile,nlat_tile,ntsteps))
    allocate(garray(NX_MMF, NY_MMF, ntsteps))
    garray = LDT_rc%udef  ! global array (43200,21600) is constructed by assebling 36x18 # of tarrays (1200,1200) 

    TILE_COLS: do tx = 1, num_tile_lon
       
       txs = 1 + (tx-1)*nlon_tile
       txe = tx * nlon_tile

       TILE_ROWS: do ty = 1, num_tile_lat
          
          tarray = missing ! GEOGRID tile array 1200,1200
          tys = 1 + (ty-1)*nlat_tile
          tye = ty * nlat_tile

          write(tileName, fmt='(4(a1,i5.5))') '/',txs, '-', txe, '.', tys, '-', tye
          geogridFile = trim(datadir)//trim(tileName)
          
          rc = read_geogrid(trim(geogridFile),len(trim(geogridFile)),tarray, &
               nlon_tile,nlat_tile,ntsteps,isigned,endian,1.,wordsize,status)

          if (rc == 1 .or. status == 1) then
             write(LDT_logunit,*) '[ERROR] reading GEOGRID file : ',trim(geogridFile)
             call LDT_endrun
           end if
             
          where(tarray == missing)
             tarray = LDT_rc%udef
          elsewhere
             tarray = SF * tarray
          end where
          
          garray(txs:txe,tys:tye,1) = tarray (1:nlon_tile, 1:nlat_tile,1)
          
       end do TILE_ROWS
    end do TILE_COLS

    ! regrid garray and construct lisout array on the LIS grid
    
    call regrid_to_lisgrid (nest, garray, lisout)
    
    deallocate (tarray, garray)
    
  contains
      
    ! ---------------------------------------------------------------------
    
    SUBROUTINE regrid_to_lisgrid (nest, data_in, var_subset)
      
      ! adapted from params/irrigation/read_GRIPC_irrigfrac.F90 
      
      implicit none
      integer, intent (in)                  :: nest
      real, dimension(:,:,:), intent   (in) :: data_in
      real, dimension(:,:,:), intent(inout) :: var_subset
      
      ! arrays read from LL
      integer,parameter :: k10 = selected_int_kind(10)
      integer   :: mi                     ! Total number of input param grid array points
      integer   :: mo                     ! Total number of output LIS grid array points
      integer   :: nc, nr, i, j
      integer, allocatable  :: n11(:)     ! Array that maps the location of each input grid
      !   point in the output grid. 
      real,    allocatable  :: gi(:)      ! Input parameter 1d grid
      logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)
      real, allocatable, dimension (:,:)    :: var_in
      real, allocatable, dimension (:)      :: go2     ! Output lis 1d grid
      logical*1, allocatable, dimension (:) :: lo2  ! Output logical mask (to match go)
      
      mi = MBR%subpnc*MBR%subpnr
      mo = LDT_rc%lnc(nest)*LDT_rc%lnr(nest)
      
      allocate (var_in(MBR%subpnc,MBR%subpnr))                
      allocate ( li(mi), gi (mi), n11(mi))
      allocate (go2(mo), lo2(mo))
      
      !- Create mapping between parameter domain and LIS grid domain:
      call upscaleByAveraging_input( MBR%subparam_gridDesc, &
           LDT_rc%gridDesc(nest,:), mi, mo, n11)
      
      Data_Fields: do j = 1, size (data_in, 3)
         
         var_in = LDT_rc%udef
         do nr = 1, MBR%subpnr
            do nc = 1, MBR%subpnc
               var_in (nc,nr) = data_in (MBR%lon_line(nc,1),MBR%lat_line (1,nr),j)
            end do
         end do
         
         ! -------------------------------------------------------------------
         !     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
         ! -------------------------------------------------------------------
         
         gi  = LDT_rc%udef
         li  = .false.
         
         !- Assign 2-D array to 1-D for aggregation routines:
         i = 0
         do nr = 1, MBR%subpnr
            do nc = 1, MBR%subpnc
               i = i + 1            
               if( var_in(nc,nr) .NE. LDT_rc%udef) then  
                  gi(i) = var_in(nc,nr)  
                  li(i) = .true.
               endif
            end do
         enddo
         
         !- Spatial average within each coarse gridcell:
         call upscaleByAveraging ( mi, mo, LDT_rc%udef, n11, li, gi, &
              lo2(:), go2 (:))
         i = 0
         do nr = 1, LDT_rc%lnr(nest)
            do nc = 1, LDT_rc%lnc(nest)
               i = i + 1
               var_subset (nc,nr,j) = go2(i)
            enddo
         enddo
         
      end do Data_Fields
      deallocate (var_in, li, gi, n11, go2, lo2)
      
    END SUBROUTINE regrid_to_lisgrid
    
  END SUBROUTINE mmf_data_reader

  ! ----------------------------------------------------------------

  SUBROUTINE cell_area_curve (MBR, nest, area)
    
    use map_utils,        only : ij_to_latlon
    use LDT_constantsMod, ONLY : radius => LDT_CONST_REARTH, pi => LDT_CONST_PI
    
    implicit none
    
    class (MMF_BCsReader), intent(inout)    :: MBR    
    integer, intent (in)                    :: nest
    real, dimension (:,:,:), intent (inout) :: area
    integer                                 :: i,j, s
    real, parameter                         :: nstrips = 100.
    real                                    :: lat_ll, lat_ur , lat_ul, lat_lr
    real                                    :: lon_ll, lon_ur , lon_ul, lon_lr
    real                                    :: c, r, d2r, dx, dyl, dyu, w_edge, e_edge, lat1, lat2

    d2r     = PI/180.
    area    = 0.
    
    do j = 1, LDT_rc%lnr(nest)
       do i = 1, LDT_rc%lnc(nest)
           
          r = float (j)
          c = float (i)
          
          call ij_to_latlon(LDT_domain(nest)%ldtproj,c-0.5, r-0.5, lat_ll, lon_ll) ! SW corner
          call ij_to_latlon(LDT_domain(nest)%ldtproj,c-0.5, r+0.5, lat_ul, lon_ul) ! NW corner          
          call ij_to_latlon(LDT_domain(nest)%ldtproj,c+0.5, r+0.5, lat_ur, lon_ur) ! NE corner         
          call ij_to_latlon(LDT_domain(nest)%ldtproj,c+0.5, r-0.5, lat_lr, lon_lr) ! SE corner

          w_edge = (lon_ll + lon_ul) / 2.
          e_edge = (lon_lr + lon_ur) / 2.
          
          dyl= (lat_lr - lat_ll) / nstrips
          dyu= (lat_ur - lat_ul) / nstrips
          dx = (e_edge - w_edge) / nstrips

          do s = 1, nstrips
             
             lat1 = lat_ll + (s-1)*dyl + dyl/2.
             lat2 = lat_ul + (s-1)*dyu + dyu/2.
             area (i,j,1) = area (i,j,1) +  (sin(d2r* lat2) - sin(d2r*lat1))*radius*radius*dx*d2r/1000./1000. ! [km2]
             
          end do
                    
       end do
    end do
    
  END SUBROUTINE cell_area_curve

  ! ----------------------------------------------------------------

  SUBROUTINE cell_area_line (MBR, nest, area)
    
    use map_utils,        only : ij_to_latlon
    use LDT_constantsMod, ONLY : radius => LDT_CONST_REARTH, pi => LDT_CONST_PI
    
    implicit none
    
    class (MMF_BCsReader), intent(inout)    :: MBR    
    integer, intent (in)                    :: nest
    real, dimension (:,:,:), intent (inout) :: area
    integer                                 :: i,j, s
    real, parameter                         :: nstrips = 100.
    real                                    :: lat_ll, lat_ur , lat_ul, lat_lr
    real                                    :: lon_ll, lon_ur , lon_ul, lon_lr
    real                                    :: c, r, d2r, dx, dw, w_edge, e_edge, lat1, lat2

    d2r     = PI/180.
    area    = 0.
    
    do j = 1, LDT_rc%lnr(nest)
       do i = 1, LDT_rc%lnc(nest)
           
          r = float (j)
          c = float (i)
          
          call ij_to_latlon(LDT_domain(nest)%ldtproj,c-0.5, r-0.5, lat_ll, lon_ll) ! SW corner
          call ij_to_latlon(LDT_domain(nest)%ldtproj,c-0.5, r+0.5, lat_ul, lon_ul) ! NW corner          
          call ij_to_latlon(LDT_domain(nest)%ldtproj,c+0.5, r+0.5, lat_ur, lon_ur) ! NE corner         
          call ij_to_latlon(LDT_domain(nest)%ldtproj,c+0.5, r-0.5, lat_lr, lon_lr) ! SE corner

          w_edge = (lon_ll + lon_ul) / 2.
          e_edge = (lon_lr + lon_ur) / 2.
          
          dx = (e_edge - w_edge) / nstrips
          dw = 1./nstrips
          do s = 1, nstrips
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c-0.5 + dw/2. + (s-1)*dw, r-0.5, lat1, lon_ll)
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c-0.5 + dw/2. + (s-1)*dw, r+0.5, lat2, lon_ll)
             area (i,j,1) = area (i,j,1) +  (sin(d2r* lat2) - sin(d2r*lat1))*radius*radius*dx*d2r/1000./1000. ! [km2]
             
          end do
                    
       end do
    end do
    
  END SUBROUTINE cell_area_line

  ! ----------------------------------------------------------------

  SUBROUTINE cell_area_girard (MBR, nest, area)
    
    use map_utils,        only : ij_to_latlon
    use LDT_constantsMod, ONLY : radius => LDT_CONST_REARTH, pi => LDT_CONST_PI
    
    implicit none

    class (MMF_BCsReader), intent(inout)    :: MBR    
    integer, intent (in)                    :: nest
    real, dimension (:,:,:), intent (inout) :: area
    integer                                 :: i,j
    real                                    :: lat_ll, lat_ur , lat_ul, lat_lr, c, r
    real                                    :: lon_ll, lon_ur , lon_ul, lon_lr, lat, lon
    real                                    :: ab, bc, cd, da, ac  ! side lengths
    
    area    = 0.

    do j = 1, LDT_rc%lnr(nest)
       do i = 1, LDT_rc%lnc(nest)
           
          r = float (j)
          c = float (i)
             
          if(trim(LDT_rc%lis_map_proj(nest)) == 'latlon') then
             
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c , r, lat, lon) ! center
             area (i,j,1) = radius * radius * &
                  (sin(d2r(lat + 0.5*LDT_rc%gridDesc(nest,10))) - &
                  sin(d2r(lat - 0.5*LDT_rc%gridDesc(nest,10))))*  &
                  (d2r(LDT_rc%gridDesc(nest,9)))/1000./1000.    ! [km2]
             
          else

             
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c-0.5, r-0.5, lat_ll, lon_ll) ! SW corner (A)
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c-0.5, r+0.5, lat_ul, lon_ul) ! NW corner (B)         
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c+0.5, r+0.5, lat_ur, lon_ur) ! NE corner (C)        
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c+0.5, r-0.5, lat_lr, lon_lr) ! SE corner (D)
             
             ! side lengths
             ab = haversine(d2r(lat_ll), d2r(lon_ll), d2r(lat_ul), d2r(lon_ul))
             bc = haversine(d2r(lat_ul), d2r(lon_ul), d2r(lat_ur), d2r(lon_ur))
             cd = haversine(d2r(lat_ur), d2r(lon_ur), d2r(lat_lr), d2r(lon_lr))
             da = haversine(d2r(lat_ll), d2r(lon_ll), d2r(lat_lr), d2r(lon_lr))
             ac = haversine(d2r(lat_ll), d2r(lon_ll), d2r(lat_ur), d2r(lon_ur))
             
             ! area = ABC area + ACD area
             
             area (i,j,1) = radius * radius * &
                  (triangle_area (ac, ab, bc) + triangle_area(ac, cd, da))/1000./1000. ! [km2]
          endif
          
       end do
    end do

  contains

    ! *****************************************************************************
    
    real function haversine(deglat1,deglon1,deglat2,deglon2)
      ! great circle distance 
      real,intent(in) :: deglat1,deglon1,deglat2,deglon2
      real            :: a,c, dlat,dlon,lat1,lat2
        
      dlat = deglat2-deglat1
      dlon = deglon2-deglon1
      lat1 = deglat1
      lat2 = deglat2     
      a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
      if(a>=0. .and. a<=1.) then
         c = 2*atan2(sqrt(a),sqrt(1-a))
         haversine = c ! [per unit radius]
      else
         haversine = 1.e20
      endif
    end function haversine


    ! *****************************************************************************
    
    function triangle_area (c, a, b) result(ABC_area)
         
      ! degrees to radians
      real,intent(in) :: c, a, b
      real :: ACB, ABC_area

      ! The spherical law of cosines per unit radius
      ACB =  acos ((cos(c) - cos(a) * cos(b)) / sin(a) / sin (b))

      ! Area = spherical excess per unit radius
      ABC_area = ABS(2.* atan( tan(a/2.)*tan(b/2.)*sin(ACB) / &
           (1. + tan(a/2.)*tan(b/2.)*cos(ACB))))

    end function triangle_area

    ! *****************************************************************************
    
    function d2r (degree) result(rad)
         
      ! degrees to radians
      real,intent(in) :: degree
      real :: rad
   
      rad = degree*PI/180.
   
    end function d2r
      
  END SUBROUTINE cell_area_girard

  ! ----------------------------------------------------------------
  
end module MMF_groundwater

  
