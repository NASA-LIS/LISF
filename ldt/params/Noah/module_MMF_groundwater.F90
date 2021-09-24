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

  use ESMF
  use LDT_logMod,           only : LDT_logunit, LDT_endrun, LDT_verify
  use LDT_coreMod,          only : LDT_rc, LDT_domain
  use LDT_gridmappingMod,   only : LDT_RunDomainPts
  use LDT_paramMaskCheckMod,only : LDT_fillopts, LDT_contIndivParam_Fill
  use LDT_paramDataMod,     only : LDT_LSMparam_struc
  implicit none

  private

  integer,parameter :: NX_MMF=43200, NY_MMF = 21600
 
  type :: geogrid

     ! defaults    
     real   :: sf       = 1.
     real   :: undef    = -9999.
     logical:: DE       = .true.
     integer:: wordsize = 2
     integer:: tile_x   = 1200
     integer:: tile_y   = 1200
     integer:: tile_z   = 1
     integer:: tile_bdr = 0
     integer:: endian   = 0
     integer:: iSigned  = 1
     
  end type geogrid
  
  type, public :: MMF_mapping
     
     real, dimension(20)              :: param_gridDesc, subparam_gridDesc
     integer                          :: glpnc, glpnr, subpnc, subpnr   
     integer, pointer, dimension (:,:):: lat_line, lon_line
 
  end type MMF_mapping
 

  type, public, extends (MMF_mapping)  ::  MMF_BCsReader

     type(LDT_fillopts), public        :: gap_fill
     type(geogrid),      public        :: gp
     
   contains
     
     procedure, public :: mi => mmf_init
     procedure, public :: mr => mmf_data_reader
     
  end type MMF_BCsReader

  public :: cell_area
  
  interface cell_area
     module procedure cell_area_girard
  end interface cell_area
  
contains
  
  SUBROUTINE mmf_init (MBR, nest, project, DATADIR, map)
    
    implicit none
    
    class (MMF_BCsReader), intent(inout) :: MBR
    type(MMF_mapping),intent(in),optional:: map
    integer, intent (in)                 :: nest
    character(*), intent (in)            :: project
    CHARACTER(*), INTENT(IN)             :: DATADIR
    real                                 :: IN_xres, IN_yres
    real                                 :: dx, dy, known_lon
    type(geogrid)                        :: GG
    type(ESMF_Config)                    :: GCF
    integer                              :: RC
    character*3                          :: signed
    integer, save, allocatable, target, dimension (:,:):: tlat_line, tlon_line
    
    ! Reads index file for from the GEOGRID data directory for GEOGRID tile information and data conversion 
    
    GCF = ESMF_ConfigCreate(RC=RC)                                  ; call LDT_verify(rc,"mmf_init: create failed."    ) 
    CALL ESMF_ConfigLoadFile     (GCF,trim(DATADIR)//'index',rc=rc) ; call LDT_verify(rc,"mmf_init: load index failed.") 
    CALL ESMF_ConfigGetAttribute (GCF, label='dx:'           , VALUE=dx           , RC=RC ); call LDT_verify(rc,"mmf_init: dx not defined.") 
    CALL ESMF_ConfigGetAttribute (GCF, label='dy:'           , VALUE=dy           , RC=RC ); call LDT_verify(rc,"mmf_init: dy not defined.") 
    CALL ESMF_ConfigGetAttribute (GCF, label='scale_factor:' , VALUE=MBR%gp%sf       ,DEFAULT=GG%sf       , RC=RC )
    CALL ESMF_ConfigGetAttribute (GCF, label='missing_value:', VALUE=MBR%gp%undef    ,DEFAULT=GG%undef    , RC=RC )
    CALL ESMF_ConfigGetAttribute (GCF, label='wordsize:'     , VALUE=MBR%gp%wordsize ,DEFAULT=GG%wordsize , RC=RC )
    CALL ESMF_ConfigGetAttribute (GCF, label='tile_x:'       , VALUE=MBR%gp%tile_x   ,DEFAULT=GG%tile_x   , RC=RC )
    CALL ESMF_ConfigGetAttribute (GCF, label='tile_y:'       , VALUE=MBR%gp%tile_y   ,DEFAULT=GG%tile_y   , RC=RC )
    CALL ESMF_ConfigGetAttribute (GCF, label='tile_z:'       , VALUE=MBR%gp%tile_z   ,DEFAULT=GG%tile_z   , RC=RC )
    CALL ESMF_ConfigGetAttribute (GCF, label='tile_bdr:'     , VALUE=MBR%gp%tile_bdr ,DEFAULT=GG%tile_bdr , RC=RC )
    CALL ESMF_ConfigGetAttribute (GCF, label='signed:'       , VALUE=signed       ,DEFAULT= 'yes'      , RC=RC )
    CALL ESMF_ConfigGetAttribute (GCF, label='known_lon:'    , VALUE=known_lon    ,DEFAULT= -180.      , RC=RC )
    CALL ESMF_ConfigDestroy      (GCF, RC=RC)

    if(trim(signed) == 'yes') then
       MBR%gp%isigned = 1
    else
       MBR%gp%isigned = 0
    endif

    ! the western edge
    if(known_lon < -179.) then
       MBR%gp%DE = .true.
    else
       MBR%gp%DE = .false.
    endif
    
    ! Verify NX_MMF, MY_MMF match with dx, dy
    if ((NINT (360./dx) /= NX_MMF) .OR. (NINT (180./dy) /= NY_MMF)) then
       write(LDT_logunit,*) "[ERR] MMF_INIT: DX, DY are NOT consistent with NX_MMF and NY_MMF"
       write(LDT_logunit,*) dx, dy, nx_mmf, ny_mmf
       write(LDT_logunit,*) "  Programming stopping ..."
       call LDT_endrun
    endif
    
    if(.not. present (map)) then

       ! NOTE: Since source data array dimensions are the same for all MMF parameters,
       !       we run this block once with the first parameter and copy mapping arrays to other parameters.
       
       ! ------------------------------------------------------------
       !    PREPARE SUBSETTED PARAMETER GRID FOR READING IN BCS DATA
       ! ------------------------------------------------------------    
       !- Map Parameter Grid Info to LIS Target Grid/Projection Info --
       
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
           
       MBR%subparam_gridDesc = 0.

       if(allocated (tlat_line)) deallocate (tlat_line)
       if(allocated (tlon_line)) deallocate (tlon_line)
       
       call LDT_RunDomainPts( nest, project, MBR%param_gridDesc(:),        &
            MBR%glpnc, MBR%glpnr, MBR%subpnc, MBR%subpnr, MBR%subparam_gridDesc, &
            tlat_line, tlon_line)
       MBR%lat_line => tlat_line
       MBR%lon_line => tlon_line
    else
       
       MBR%param_gridDesc   = map%param_gridDesc
       MBR%glpnc            = map%glpnc     
       MBR%glpnr            = map%glpnr
       MBR%subpnc           = map%subpnc
       MBR%subpnr           = map%subpnr
       MBR%subparam_gridDesc= map%subparam_gridDesc
       MBR%lat_line => map%lat_line
       MBR%lon_line => map%lon_line
              
    endif
    
  END SUBROUTINE mmf_init
  
  ! ----------------------------------------------------------------
  
  SUBROUTINE mmf_data_reader (MBR, nest, datadir, lisout, mmf_transform)
    
    implicit none

    class (MMF_BCsReader), intent(inout)    :: MBR
    character(*), intent (in)               :: datadir, mmf_transform
    real, dimension (:,:,:), intent (inout) :: lisout
    integer, intent (in)                    :: nest    
    character(len=1024):: geogridFile
    character(len=80)  :: tileName
    integer            :: num_tile_lon,num_tile_lat
    integer            :: rc,status,tx,ty,txs,txe,tys,tye
    integer            :: nlat_tile_bdr, nlon_tile_bdr, nbins
    real,allocatable   :: tarray(:,:,:), garray(:,:,:)
      
    ! c function from Michael G. Duda, NCAR/MMM
    integer,  external :: read_geogrid
    ! Adapted from Zhuo Wang's
    ! /discover/nobackup/projects/lis_aist17/zwang9/Geogrid/geogrid2netcdf.f90


    ! Reading GEOGRID data
    ! --------------------

    num_tile_lon = NX_MMF / MBR%gp%tile_x
    num_tile_lat = NY_MMF / MBR%gp%tile_y

    nlon_tile_bdr = MBR%gp%tile_x + 2*MBR%gp%tile_bdr
    nlat_tile_bdr = MBR%gp%tile_y + 2*MBR%gp%tile_bdr
    
    allocate(tarray(nlon_tile_bdr,nlat_tile_bdr,MBR%gp%tile_z))
    allocate(garray(NX_MMF, NY_MMF, MBR%gp%tile_z))

    garray = LDT_rc%udef  ! global array (43200,21600) is constructed by assebling 36x18 # of tarrays (1200,1200) 

    TILE_COLS: do tx = 1, num_tile_lon
       
       txs = 1 + (tx-1)*MBR%gp%tile_x
       txe = tx * MBR%gp%tile_x

       TILE_ROWS: do ty = 1, num_tile_lat
          
          tarray = MBR%gp%undef 
          tys = 1 + (ty-1)*MBR%gp%tile_y
          tye = ty * MBR%gp%tile_y

          if(MBR%gp%DE) then
             ! The western edge is on the dateline
             write(tileName, fmt='(4(a1,i5.5))') '/',txs, '-', txe, '.', tys, '-', tye

          else
             ! Elevation data x-axis start from the Greenwich line
             if(txe <= 21600) then
                write(tileName, fmt='(4(a1,i5.5))') '/',txs+21600, '-', txe+21600, '.', tys, '-', tye
             else
                write(tileName, fmt='(4(a1,i5.5))') '/',txs-21600, '-', txe-21600, '.', tys, '-', tye
             endif
             
          endif
          
          geogridFile = trim(datadir)//trim(tileName)
          
          rc = read_geogrid(trim(geogridFile),len(trim(geogridFile)),tarray, &
               MBR%gp%tile_x,MBR%gp%tile_y,MBR%gp%tile_z,MBR%gp%isigned,MBR%gp%endian,1.,MBR%gp%wordsize,status)

          if (rc == 1 .or. status == 1) then
             write(LDT_logunit,*) '[ERROR] reading GEOGRID file : ',trim(geogridFile)
             call LDT_endrun
           end if
             
          where(tarray == MBR%gp%undef)
             tarray = LDT_rc%udef
          elsewhere
             tarray = MBR%gp%SF * tarray
          end where
          
          garray(txs:txe,tys:tye,1) = tarray (1 + MBR%gp%tile_bdr : MBR%gp%tile_bdr + MBR%gp%tile_x, 1 + MBR%gp%tile_bdr: MBR%gp%tile_bdr + MBR%gp%tile_y,1)
          
       end do TILE_ROWS
    end do TILE_COLS

    ! regrid garray and construct lisout array on the LIS grid
    
    call regrid_to_lisgrid (nest, garray, lisout)

    ! fill gaps

    nbins = size (lisout, 3)

    call LDT_contIndivParam_Fill( nest, LDT_rc%lnc(nest), LDT_rc%lnr(nest),  &
         mmf_transform,                                    &
         nbins,                                            &
         lisout,  LDT_rc%udef,                             &
         LDT_LSMparam_struc(nest)%landmask2%value,         &
         MBR%gap_fill%filltype, MBR%gap_fill%fillvalue,    &
         MBR%gap_fill%fillradius )
             
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

  SUBROUTINE cell_area_line (nest, area)
    
    use map_utils,        only : ij_to_latlon
    use LDT_constantsMod, ONLY : radius => LDT_CONST_REARTH, pi => LDT_CONST_PI
    
    implicit none
    
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
    
  END SUBROUTINE cell_area_line

  ! ----------------------------------------------------------------

  SUBROUTINE cell_area_curve (nest, area)
    
    use map_utils,        only : ij_to_latlon
    use LDT_constantsMod, ONLY : radius => LDT_CONST_REARTH, pi => LDT_CONST_PI
    
    implicit none
    
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
    
  END SUBROUTINE cell_area_curve

  ! ----------------------------------------------------------------

  SUBROUTINE cell_area_girard (nest, area)
    
    use map_utils,        only : ij_to_latlon
    use LDT_constantsMod, ONLY : radius => LDT_CONST_REARTH, pi => LDT_CONST_PI
    
    implicit none

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

  
