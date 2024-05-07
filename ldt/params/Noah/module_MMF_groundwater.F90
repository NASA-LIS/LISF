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
! !MODULE: MMF_groundwater
! \label{MMF_groundwater}
!
! !REVISION HISTORY:
!  27 Aug 2021: Sarith Mahanama
!     with inputs from David Mocko, Shugong Wang, Timothy Lahmers, Kristi Arsenault, and Zhuo Wang ; Initial Specification
!
! !INTERFACE:

#include "LDT_misc.h"
module  MMF_groundwater

  use ESMF
  use LDT_logMod,           only : LDT_logunit, LDT_endrun, LDT_verify
  use LDT_coreMod,          only : LDT_rc, LDT_domain, LDT_config
  use LDT_gridmappingMod,   only : LDT_RunDomainPts
  use LDT_paramMaskCheckMod,only : LDT_fillopts, LDT_contIndivParam_Fill, fill_logunit
  use LDT_paramDataMod,     only : LDT_LSMparam_struc
  implicit none

  private

  integer,parameter :: NX_MMF=43200, NY_MMF = 21600

  ! defaults parameters for GEOGRID tiles based on the FDEPTH data set
  ! ------------------------------------------------------------------
  type :: geogrid

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

  ! MMF grid to LIS grid mapping
  ! ----------------------------
  type, public :: MMF_mapping
   
     real, dimension(20)              :: param_gridDesc, subparam_gridDesc
     integer                          :: glpnc, glpnr, subpnc, subpnr   
     integer, pointer, dimension (:,:):: lat_line, lon_line
 
  end type MMF_mapping

  ! MMF data reader Fortran object
  ! ------------------------------
  type, public, extends (MMF_mapping)  ::  MMF_BCsReader

     type(LDT_fillopts), public        :: gap_fill
     type(geogrid),      public        :: gp
     real, public                      :: water_value
     
   contains
     
     procedure, public :: mi => mmf_init
     procedure, public :: mr => mmf_data_reader
     
  end type MMF_BCsReader

  ! Grid cell area
  ! --------------
  public :: cell_area
    
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
  
  SUBROUTINE mmf_data_reader (MBR, nest, datadir, lisout, mmf_transform, short_name)
    
    implicit none

    class (MMF_BCsReader), intent(inout)    :: MBR
    character(*), intent (in)               :: datadir, mmf_transform, short_name
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
    write(LDT_logunit,*) "Reading data from: ",trim(datadir)
    
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
               nlon_tile_bdr,nlat_tile_bdr,MBR%gp%tile_z,MBR%gp%isigned,MBR%gp%endian,1.,MBR%gp%wordsize,status)

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

    nbins = size (lisout, 3)
    
    ! Search and fill gaps
    ! --------------------
    write(LDT_logunit,*) "Checking/filling mask values for: ", &
                 trim(short_name)
    write(fill_logunit,*) "Checking/filling mask values for: ", &
                 trim(short_name)   

    call LDT_contIndivParam_Fill( nest, LDT_rc%lnc(nest), LDT_rc%lnr(nest),  &
         mmf_transform,                                    &
         nbins,                                            &
         lisout,  LDT_rc%udef,                             &
         LDT_LSMparam_struc(nest)%landmask2%value,         &
         MBR%gap_fill%filltype, MBR%gap_fill%fillvalue,    &
         MBR%gap_fill%fillradius, leave_good_data = .true. )

    if (trim (short_name) /= "MMF_HGTM") then
    ! fill suitable parameter values in lakes/water bodies
       where (lisout == LDT_rc%udef)
          lisout = MBR%water_value
       end where
    endif
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

  SUBROUTINE cell_area (nest, area)

    use map_utils,        only : ij_to_latlon
    use LDT_constantsMod, ONLY : radius => LDT_CONST_REARTH, pi => LDT_CONST_PI
    
    implicit none

 !   interface grid_cell_area
 !      module procedure area_latlon
 !      module procedure areaint
 !      module procedure area_wps
 !   end interface grid_cell_area
 
    integer, intent (in)                    :: nest
    real, dimension (:,:,:), intent (inout) :: area
    integer                                 :: i,j
    real                                    :: lat_ll, lat_ur , lat_ul, lat_lr, c, r
    real                                    :: lon_ll, lon_ur , lon_ul, lon_lr
   
    area    = 0.

    do j = 1, LDT_rc%lnr(nest)
       do i = 1, LDT_rc%lnc(nest)
           
          r = float (j)
          c = float (i)

          select case (LDT_domain(nest)%ldtproj%code)
          case (0)
             ! lat/lon
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c, r, lat_ll, lon_ll)
             area (i,j,1) = area_latlon (nest, lat_ll)
             
          case (3)
             ! Lambert conical follows WPS
             area (i,j,1) = area_wps ()
             
          case DEFAULT
             ! Area of a polygon areaint.m from Matlab
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c-0.5, r-0.5, lat_ll, lon_ll) ! SW corner (A)
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c-0.5, r+0.5, lat_ul, lon_ul) ! NW corner (B)         
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c+0.5, r+0.5, lat_ur, lon_ur) ! NE corner (C)        
             call ij_to_latlon(LDT_domain(nest)%ldtproj,c+0.5, r-0.5, lat_lr, lon_lr) ! SE corner (D)
             area (i,j,1) = areaint((/lat_ll, lat_ul, lat_ur, lat_lr/), (/lon_ll, lon_ul, lon_ur, lon_lr/))
             
          END select
       end do
    end do
    
  contains

    ! ----------------------------------------------------------------

    real function area_latlon (nest, lat)

      implicit none

      real, intent (in)    :: lat
      integer, intent (in) :: nest
      
      area_latlon = radius * radius * &
                 (sin(d2r(lat + 0.5*LDT_rc%gridDesc(nest,9)))  - &
                  sin(d2r(lat - 0.5*LDT_rc%gridDesc(nest,9))))*  &
                  d2r(LDT_rc%gridDesc(nest,10))/1000./1000.    ! [km2]
      
    end function area_latlon

    ! ----------------------------------------------------------------

    real function area_wps ()

      implicit none
      integer           :: rc
      real              :: DXY, MSFTX, MSFTY

      MSFTY = 1.
      MSFTX = 1.

      call ESMF_ConfigGetAttribute(LDT_config, DXY,label='Run domain resolution:', rc=rc); call LDT_verify(rc,"Run domain resolution of the Lambert grid is not defined.")
      area_wps = DXY*DXY/MSFTX/MSFTY
      
    end function area_wps
    
    ! ----------------------------------------------------------------
    
    real FUNCTION areaint (lat, lon)
      
      
      ! simplified from Matlab's areaint.m to compute area of a single polygon
      ! AREAINT Surface area of polygon on sphere 
      !   A = AREAINT(LAT,LON) calculates the spherical surface area of the
      !   polygon specified by the input vectors LAT, LON.  LAT and LON are in
      !   degrees.  The calculation uses a line integral approach.  The output,
      !   A, is the surface area fraction covered by the polygon on a unit
      !   sphere.
      
      implicit none
      real, intent(in), dimension(:)    :: lat, lon
      real, allocatable , dimension (:) :: latc, lonc, colat, az, integrands 
      real                              :: lat0,lon0,dlat,dlon,a,deltas,daz,colat2
      integer                           :: n, i
      
      n = size (lat) + 1
      allocate (latc (1:n))
      allocate (lonc (1:n))
      allocate (colat(1:n))
      allocate (az   (1:n))
      
      latc(1:n-1) = lat
      lonc(1:n-1) = lon
      latc(n)     = lat(1)
      lonc(n)     = lon(1)
      lat0 = 0.
      lon0 = 0.
      
      ! greatcircle distance, and greatcircle azimuth wrt 0.,0 (Matlab's distance.m)
      ! ----------------------------------------------------------------------------
      
      do i = 1,n
       
         latc(i) = d2r(latc(i))
         lonc(i) = d2r(lonc(i))
         dlat     = latc(i) - lat0
         dlon     = lonc(i) - lon0
         
         ! haversine
         a        = (sin(dlat/2.))**2 + cos(lat0)*cos(latc(i))*(sin(dlon/2.))**2
         if(a < 0.) a =0.
         if(a > 1.) a =1.
         
         colat(i) = 2.*atan2(sqrt(a),sqrt(1.-a))         
         az(i)    = atan2(cos(latc(i)) * sin(lonc(i)-lon0),  &
              cos(lat0) * sin(latc(i)) - sin(lat0) * cos(latc(i))* cos(lonc(i)-lon0))
         ! wrap az to the range 0-2pi
         az(i)    = az(i) - 2.*pi*floor(az(i)/2./pi)
         
      end do
      
      n = n -1
      allocate (integrands (1:n))
      
      do i = 1, n
         
         ! Calculate step sizes
         daz = az(i+1) - az(i)
         ! wrap to -pi <= daz <=pi
         daz = daz - 2.*pi*floor((daz+pi)/2./pi) 
         
         ! Determine average surface distance for each step
         deltas = (colat (i+1) - colat (i))/2.
         colat2 = colat(i) + deltas
         
         ! Integral over azimuth is 1-cos(colatitudes)
         integrands (i) = (1. - cos(colat2)) * daz
      end do
      
      areaint = abs (sum (integrands))/4./pi
      areaint = 4.*pi*radius*radius * MIN (areaint, 1. - areaint) /1000./1000.
      
      deallocate (integrands, latc, lonc, colat, az)
      
    end FUNCTION areaint

    ! ----------------------------------------------------------------
    
    function d2r (degree) result(rad)
      
      ! degrees to radians
      real,intent(in) :: degree
      real :: rad
      
      rad = degree*PI/180.
      
    end function d2r

  end SUBROUTINE cell_area
   
end module MMF_groundwater

  
