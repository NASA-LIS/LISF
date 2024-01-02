!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module tile_coord_routines
!BOP
!
!  !MODULE: tile_corrd_routines
! 
!  !DESCRIPTION:
!  This file contains types and subroutines for tile coordinates and domain.
!
!  !REVISION HISTORY: 
!  26 Jan 2005 Rolf Reichle, Initial Specification
!  14 Apr 2006 Rolf Reichle, split tile_coord.F90 into 2 files to avoid 
!                            having more than one module per file
!  10 Jul 2006 James Geiger, Implementation in LIS
! 
! !USES:
  use tile_coord_types
!EOP
  use LIS_logMod, only : LIS_logunit

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------  
  PUBLIC :: read_tile_coord
  PUBLIC :: read_tile_coord_header

contains
  
!BOP
!
! !ROUTINE: read_tile_coord_header
!  \label{read_tile_coord_header}
!
! !REVISION HISTORY:
!             Rolf Reichle, Initial Specification
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! !INTERFACE:
  subroutine read_tile_coord_header( ftn, tile_coord_file, atm_grid, N_tile )
    
    
    implicit none
! !ARGUMENTS: 

    integer,             intent(IN)            :: ftn
    character(len=*),    intent(in)            :: tile_coord_file
    type(grid_def_type), intent(out), optional :: atm_grid
    
    integer,             intent(OUT)           :: N_tile
    
!
! !DESCRIPTION:
! read tile coordinates from GEOS5 *.til file 
! optionally assemble atmospheric grid parameters
!
! inputs:
!
!  tile\_coord\_file : name with full path of *.til file
!
! outputs:
!
!  atm\_grid   : atmospheric grid parameters 
!
!  tile\_coord : coordinates of tiles (see tile\_coord\_type),
!               implemented as allocatable which is allocated in 
!               this subroutine
!               NOTE: number of tiles can be diagnosed 
!                     with size(tile\_coord)
!
! -------------------------------------------------------------
! 
!EOP
    
    ! locals
    
    integer :: i, tmpint1, tmpint2, tmpint3
    
    logical :: date_line_on_center, pole_on_center
    
    ! ---------------------------------------------------------------
    
    ! read file header 
    
    write (LIS_logunit,*) 'read_tile_coord(): reading from'
    write (LIS_logunit,*) tile_coord_file
    write (LIS_logunit,*)
    
    read (ftn,*) N_tile
    
    write (LIS_logunit,*) 'file contains coordinates for ', N_tile, ' tiles' 
    write (LIS_logunit,*)
    
    read (ftn,*)    ! some number (?)
    read (ftn,*)    ! some string describing atmospheric grid (?)
    
    read (ftn,*)    tmpint1
    read (ftn,*)    tmpint2
    
    read (ftn,*)    ! some string describing ocean grid                   (?)
    read (ftn,*)    ! # ocean grid cells in longitude direction (N_i_ocn) (?)
    read (ftn,*)    ! # ocean grid cells in latitude direction (N_j_ocn)  (?)
    
    ! ---------------------------------------------------------------
    
    if (present(atm_grid)) then       ! extract atm_grid from header
       
       atm_grid%N_lon = tmpint1
       atm_grid%N_lat = tmpint2
       
       if     (atm_grid%N_lon==144 .and. atm_grid%N_lat==91 ) then
          
          atm_grid%dlon = 2.5
          atm_grid%dlat = 2.
          
       elseif (atm_grid%N_lon==288 .and. atm_grid%N_lat==181) then
          
          atm_grid%dlon = 1.25
          atm_grid%dlat = 1.
          
       elseif (atm_grid%N_lon==360 .and. atm_grid%N_lat==180) then
          
          atm_grid%dlon = 1.
          atm_grid%dlat = 1.
          
       elseif (atm_grid%N_lon==540 .and. atm_grid%N_lat==361) then
          
          atm_grid%dlon = 2./3.
          atm_grid%dlat = .5
          
       elseif (atm_grid%N_lon==576 .and. atm_grid%N_lat==361) then
          
          atm_grid%dlon = 0.625
          atm_grid%dlat = .5

       else
          
          write (LIS_logunit,*) 'read_tile_coord: unknown atmospheric grid (N_lon=', &
               atm_grid%N_lon, ' and N_lat=', atm_grid%N_lat, ')'
          write (LIS_logunit,*) 'STOPPING.'
          stop
          
       end if
       
       ! find out whether date line is on edge or through center of grid cell
       
       if     (index(tile_coord_file, 'FV_144x91_DC_360x180_DE.til')/=0) then
          date_line_on_center = .true.
          pole_on_center      = .true.
       elseif (index(tile_coord_file, 'FV_144x91_DE_360x180_DE.til')/=0) then
          date_line_on_center = .false.
          pole_on_center      = .true.
       elseif (index(tile_coord_file, 'FV_288x181_DC_360x180_DE.til')/=0) then
          date_line_on_center = .true.
          pole_on_center      = .true.
       elseif (index(tile_coord_file, 'FV_540x361_DC_360x180_DE.til')/=0) then
          date_line_on_center = .true.
          pole_on_center      = .true.
       elseif (index(tile_coord_file, 'FV_576x361_DC_360x180_DE.til')/=0) then
          date_line_on_center = .true.
          pole_on_center      = .true.
          
          ! GSWP2 1 deg by 1 deg (original GSWP-2 resolution)
          
       elseif (index(tile_coord_file, 'FV_360x180_DE_360x180_DE.til')/=0) then
          date_line_on_center = .false.
          pole_on_center      = .false.

          ! GSWP2 1 deg by 1 deg w/ irregular-shaped sub-tiles
          
       elseif (index(tile_coord_file, 'PE_360x180_DE_288x270_DE.til')/=0) then
          date_line_on_center = .false.
          pole_on_center      = .false.
          
       else
          write (LIS_logunit,*) 'read_tile_coord: unknown atmospheric grid (file=', &
               tile_coord_file, ')'
          write (LIS_logunit,*) 'STOPPING.'
          stop
       end if
       
       
       
       if (pole_on_center) then
          atm_grid%ll_lat = -90. - atm_grid%dlat/2.
       else
          atm_grid%ll_lat = -90. 
       end if

       if (date_line_on_center) then
          atm_grid%ll_lon = -180. - atm_grid%dlon/2.
       else
          atm_grid%ll_lon = -180. 
       end if
       
       write (LIS_logunit,*) 'atmospheric grid:'
       if (date_line_on_center) then
          write (LIS_logunit,*) '  date line on center'
       else
          write (LIS_logunit,*) '  date line on edge'
       end if
       if (pole_on_center) then
          write (LIS_logunit,*) '  pole on center'
       else
          write (LIS_logunit,*) '  pole on edge'
       end if
       write (LIS_logunit,*) '  N_lon  = ',  atm_grid%N_lon  
       write (LIS_logunit,*) '  N_lat  = ',  atm_grid%N_lat 
       write (LIS_logunit,*) '  ll_lon = ',  atm_grid%ll_lon
       write (LIS_logunit,*) '  ll_lat = ',  atm_grid%ll_lat
       write (LIS_logunit,*) '  dlon   = ',  atm_grid%dlon  
       write (LIS_logunit,*) '  dlat   = ',  atm_grid%dlat  
       
    end if
    
  end subroutine read_tile_coord_header


!BOP
!
! !ROUTINE: read_tile_coord
!  \label{read_tile_coord}
!
! !REVISION HISTORY:
!             Rolf Reichle, Initial Specification
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! !INTERFACE:
  subroutine read_tile_coord( ftn, N_tile, tile_coord )
    
    
    implicit none
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: N_tile
    type(tile_coord_type)   :: tile_coord(N_tile)
    
!
! !DESCRIPTION:
! read tile coordinates from GEOS5 *.til file 
! optionally assemble atmospheric grid parameters
!
! inputs:
!
!  tile\_coord\_file : name with full path of *.til file
!
! outputs:
!
!  atm\_grid   : atmospheric grid parameters 
!
!  tile\_coord : coordinates of tiles (see tile\_coord\_type),
!               implemented as allocatable which is allocated in 
!               this subroutine
!               NOTE: number of tiles can be diagnosed 
!                     with size(tile\_coord)
!
! -------------------------------------------------------------
! 
!EOP
    
    ! locals
    
    integer :: i, tmpint1, tmpint2, tmpint3
        
    ! ---------------------------------------------------------------
       
    do i=1,N_tile
       
          ! read statement for SiB_V1 tile coordinate file
          !
          ! read(10,'(i10,i9,2f10.4,2i5,f10.6,3i8,f10.6,i8,f12.3)')  
          
          ! read statement for SiB_V2 tile coordinate file
          !
          ! read(10,'(i10,i9,2f10.4,2i5,f16.12,3i8,f16.12,i8,f13.4)') 
          
          ! avoid using exact format specification
          
       read (ftn,*)                                              &
            tile_coord(i)%typ,                                  &
            tile_coord(i)%pfaf,                                 &
            tile_coord(i)%com_lon,                              &
            tile_coord(i)%com_lat,                              &
            tile_coord(i)%i_atm,                                &
            tile_coord(i)%j_atm,                                &
            tile_coord(i)%frac_atm,                             &
            tmpint1,                                            &
            tmpint2,                                            &
            tmpint3,                                            &
            tile_coord(i)%frac_pfaf,                            &
            tile_coord(i)%tile_id,                              &
            tile_coord(i)%area
       
    end do
       
    ! -----------------------------------------------------------------

    
  end subroutine read_tile_coord
  
  ! **********************************************************************
  
!BOP
!
! !ROUTINE: get_tile_num_in_atm_ij
!  \label{get_tile_num_in_atm_ij}
!
! !REVISION HISTORY:
! 22 Jul 2005 Rolf Reichle, Initial Specification
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! !INTERFACE:
  subroutine get_tile_num_in_atm_ij( N_tile, tile_coord, atm_grid, &
       d2g_grid_off_i, d2g_grid_off_j,                             &
       N_tile_in_atm_ij, tile_num_in_atm_ij )
    
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), allocatable :: tile_coord ! input
    
    type(grid_def_type), intent(in) :: atm_grid
    
    integer, intent(in) :: d2g_grid_off_i, d2g_grid_off_j
    
    integer, dimension(atm_grid%N_lon,atm_grid%N_lat), intent(inout) :: &
         N_tile_in_atm_ij
    
    ! the allocatable is an output arguments that is allocated here
    
    integer, dimension(:,:,:), allocatable, optional :: tile_num_in_atm_ij
!
! !DESCRIPTION:
! find out how many tiles are in a given atm grid cell
!
! When called without optional arguments only counts number of tiles.
! When called with optional arguments allocates and fills allocatable.
!
!   d2g\_grid\_off\_i, d2g\_grid\_off\_j :
!         tile\_coord%atm\_i and tile\_coord%atm\_j refer to the *global*
!         atmospheric grid that corresponds to the tile definitions
!         (as obtained from the tile\_coord\_file)
!         d2g\_grid\_off\_i and d2g\_grid\_off\_j describe the offset between
!         the global "atm\_grid\_g" and a smaller "atm\_grid\_d" for the domain
!         of interest.  With these offset arguments the subroutine
!         can be used with an atm\_grid that is a subgrid of the 
!         global atm\_grid:
!         atm\_grid\_d%ll\_lon = atm_grid\_g%ll\_lon + d2g\_grid\_off\_i*dlon
!         atm\_grid\_d%ll\_lat = atm\_grid\_g%ll\_lat + d2g\_grid\_off\_j*dlat
!
!
! ----------------------------------------------------------
!EOP

    ! locals 
    
    integer, parameter :: nodata = -9999
    
    integer :: i, j, k, n
    
    ! -----------------------------------------------------------------
    !
    ! allocate and initialize allocatables if present
    
    if (present(tile_num_in_atm_ij)) then
       
       allocate(tile_num_in_atm_ij(atm_grid%N_lon,atm_grid%N_lat,    &
            maxval(N_tile_in_atm_ij)))
       
       tile_num_in_atm_ij = nodata
       
    end if
    
    ! (re-)initialize
    
    N_tile_in_atm_ij = 0
    
    do n=1,N_tile
       
       i = tile_coord(n)%i_atm - d2g_grid_off_i
       j = tile_coord(n)%j_atm - d2g_grid_off_j
       
       N_tile_in_atm_ij(i,j) = N_tile_in_atm_ij(i,j) + 1
       
       if (present(tile_num_in_atm_ij)) then
          
          k = N_tile_in_atm_ij(i,j)
          
          tile_num_in_atm_ij(i,j,k) = n 
          
       end if
       
    end do
    
    if (.not. present(tile_num_in_atm_ij)) then
       
       write (LIS_logunit,*) 'Maximum number of tiles in atm grid cell = ', &
            maxval(N_tile_in_atm_ij)
       write (LIS_logunit,*)
       
    end if

  end subroutine get_tile_num_in_atm_ij

  ! *******************************************************************
  
!BOP
!
! !ROUTINE: extract_land_tiles
!  \label{extract_land_tiles}
!
! !REVISION HISTORY:
! 28 Jan 2005 Rolf Reichle, Initial Specification
! 22 Jul 2005 Rolf Reichle
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! !INTERFACE:
  subroutine extract_land_tiles( N_tile_global, tile_coord_global, &
       N_tile_land, tile_coord_land )
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: N_tile_global
    
    type(tile_coord_type), dimension(:), allocatable :: tile_coord_global ! input
    
    integer, intent(inout) :: N_tile_land
    
    ! the allocatable is an output arguments that is allocated here
    
    type(tile_coord_type), dimension(:), allocatable, optional :: tile_coord_land
!
! !DESCRIPTION:
! extract land tiles from tile\_coord\_global
!
! When called without optional arguments only counts number of tiles.
! When called with optional arguments allocates and fills allocatable.
!
! ----------------------------------------------------------
!EOP
    
    ! locals 
    
    integer :: n
    
    ! -----------------------------------------------------------------
    !
    ! allocate and initialize allocatables if present
    
    if (present(tile_coord_land)) allocate(tile_coord_land(N_tile_land))
    
    ! (re-)initialize
    
    N_tile_land           = 0
    
    do n=1,N_tile_global
       
       ! count number of land tiles
       
       if (tile_coord_global(n)%typ == tile_typ_land) then
          
          N_tile_land = N_tile_land + 1
          
          if (present(tile_coord_land)) &          
               tile_coord_land(N_tile_land) = tile_coord_global(n)
          
       end if
       
    end do
    
    write (LIS_logunit,*) 'Number of land tiles = ', N_tile_land 
    write (LIS_logunit,*)
    
  end subroutine extract_land_tiles
  
  ! *******************************************************************

#if 0 
!BOP
!
! !ROUTINE: get_land_tile_info
!  \label{get_land_tile_info}
!
! !REVISION HISTORY:
! 28 Jan 2005 Rolf Reichle, Initial Specification
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! !INTERFACE:
  subroutine get_land_tile_info( tile_coord_file, catchment_def_file, &
       atm_grid,                                                      &
       N_tile, tile_coord )
    
    
    implicit none
! !ARGUMENTS: 
    character(len=*), intent(in) :: tile_coord_file, catchment_def_file
    
    type(grid_def_type), intent(out) :: atm_grid
    
    integer, intent(out) :: N_tile
    
    type(tile_coord_type), dimension(:), allocatable :: tile_coord ! out
!
! !DESCRIPTION:
! get land tile coordinates, atmospheric grid parameters, and 
! grid-to-tile mapping from GEOS5 *.til file
!
! note use of optional inout arguments
!
! ----------------------------------------------------------
!EOP
    
    ! locals
    
    integer :: N_tile_tmp
    
    type(tile_coord_type), dimension(:), allocatable :: tile_coord_tmp
    
    ! --------------------------------------------------------------
    
    nullify(tile_coord_tmp)
    
    ! get atmospheric grid parameters and global tiles
    
    call read_tile_coord( tile_coord_file, atm_grid=atm_grid, &
         tile_coord=tile_coord_tmp )
    
    N_tile_tmp = size(tile_coord_tmp)
    
    ! first call counts land tiles in tile_coord_file
    
    call extract_land_tiles( N_tile_tmp, tile_coord_tmp, N_tile )
    
    ! second call allocates and fills tile_coord with land tiles
    
    call extract_land_tiles( N_tile_tmp, tile_coord_tmp, N_tile, tile_coord )
    
    deallocate(tile_coord_tmp)
    
    ! add lat-lon bounding box for each tile from "catchment.def" file
    
    call read_catchment_def( catchment_def_file, N_tile, tile_coord )
    
  end subroutine get_land_tile_info

  ! *******************************************************************
#endif 
!BOP
!
! !ROUTINE: read_catchment_def
!  \label{read_catchment_def}
!
! !REVISION HISTORY:
!             Rolf Reichle, Initial Specification
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! !INTERFACE:
  subroutine read_catchment_def( catchment_def_file, N_tile, tile_coord )
    
    implicit none
! !ARGUMENTS: 
    character(len=*), intent(in) :: catchment_def_file

    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), allocatable :: tile_coord ! inout
!
! !DESCRIPTION:
! This routine reads the catchment def file.
!EOP
    
    ! locals
    
    integer :: i, tmpint1, tmpint2
    
    
    ! ---------------------------------------------------------------
    
    ! read file header 
    
    write (LIS_logunit,*) 'read_catchment_def(): reading from'
    write (LIS_logunit,*) trim(catchment_def_file)
    write (LIS_logunit,*)
    
    open (10, file=trim(catchment_def_file), form='formatted', action='read') 
    
    read (10,*) tmpint1
    
    write (LIS_logunit,*) 'file contains coordinates for ', tmpint1, ' tiles' 
    write (LIS_logunit,*)
    
    if (N_tile/=tmpint1) then
       
       write (LIS_logunit,*) 'tile_coord_file and catchment_def_file mismatch. (1)'
       write (LIS_logunit,*) 'STOPPING.'
       stop
       
    end if
    
    do i=1,N_tile
       
       ! avoid using exact format specification
       
       read (10,*) tmpint1, tmpint2,  &
            tile_coord(i)%min_lon,    &
            tile_coord(i)%max_lon,    &
            tile_coord(i)%min_lat,    &
            tile_coord(i)%max_lat
       
       if ( (tile_coord(i)%tile_id/=tmpint1) .or. &
            (tile_coord(i)%pfaf   /=tmpint2)           ) then
          
          write (LIS_logunit,*) 'tile_coord_file and catchment_def_file mismatch. (2)'
          write (LIS_logunit,*) 'STOPPING.'
          stop
          
       end if
       
    end do
    
    ! -----------------------------------------------------------------

    close(10, status='keep')
    
  end subroutine read_catchment_def
  
  ! *******************************************************************
  
!BOP
!
! !ROUTINE: get_tile_num_from_latlon
!  \label{get_tile_num_from_latlon}
!
! !REVISION HISTORY:
!             Rolf Reichle, Initial Specification
! 17 Nov 2005 Rolf Reichle, bug fix re. "check that lat/lon is inside atm_grid"
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! !INTERFACE:
  subroutine get_tile_num_from_latlon(N_catd, tile_coord,           &
       atm_grid, N_tile_in_atm_ij, tile_num_in_atm_ij, N_latlon, lat, lon,  &
       tile_num)

    implicit none
    
! !ARGUMENTS: 
    integer, intent(in) :: N_catd, N_latlon
    
    type(tile_coord_type), dimension(:), allocatable :: tile_coord ! input
    
    type(grid_def_type), intent(in) :: atm_grid    
    
    integer, dimension(atm_grid%N_lon,atm_grid%N_lat), intent(in) :: &
         N_tile_in_atm_ij
    
    integer, dimension(:,:,:), allocatable :: tile_num_in_atm_ij ! input    
    
    real, dimension(N_latlon), intent(in)  :: lat, lon
    
    integer, dimension(N_latlon), intent(out) :: tile_num
!
! !DESCRIPTION:
! This routine determines the tile number of a tile specified by its lat/lon.
!EOP
    
    ! local variables
    
    integer, parameter :: nodata_tilenum = -9999
    
    integer :: n, k, i_atm, j_atm, this_tile_num
    
    real :: tmp_dist, min_dist
    
    ! -----------------------------------------------------------
    
    tile_num = nodata_tilenum          ! initialize to negative value
    
    do n=1,N_latlon
       
       ! determine grid cell that contains lat/lon 
       
       call get_ij_atm_from_latlon( atm_grid, lat(n), lon(n), i_atm, j_atm )
       
       ! make sure lat/lon is *inside* atm_grid
       
       if ( i_atm>=1 .and. i_atm<=atm_grid%N_lon .and.  &
            j_atm>=1 .and. j_atm<=atm_grid%N_lat           )  then
          
          ! Loop through all tiles within grid cell that contain lat/lon and
          ! find minimum distance (Minkowski norm).
          ! If there are no land tiles in given atm grid cell, tile_num will
          ! not change from its initialized value.
          
          min_dist = 1e10
          
          do k=1,N_tile_in_atm_ij( i_atm, j_atm ) 
             
             this_tile_num = tile_num_in_atm_ij(i_atm,j_atm,k)
             
             tmp_dist = &
                  abs(lon(n) - tile_coord(this_tile_num)%com_lon) + &
                  abs(lat(n) - tile_coord(this_tile_num)%com_lat) 
             
             if (tmp_dist<min_dist) then
                
                min_dist    = tmp_dist
                tile_num(n) = this_tile_num
                
             end if
             
          end do
          
          ! make sure that lat/lon is inside bounding box of given tile
          ! (if domain consists of a single tile and atm_grid is coarse,
          !  it is likely that lat/lon is outside bounding box of tile
          !  even if the lat/lon is inside the atm grid cell that contains
          !  the given tile) - reichle, 18 Aug 2005

          if (tile_num(n)>0) then
             
             if ( lon(n) < tile_coord(tile_num(n))%min_lon    .or.  &
                  lon(n) > tile_coord(tile_num(n))%max_lon    .or.  &
                  lat(n) < tile_coord(tile_num(n))%min_lat    .or.  &
                  lat(n) > tile_coord(tile_num(n))%max_lat  )       &
                  tile_num(n) = nodata_tilenum
             
          end if
          
       end if

    end do
    
  end subroutine get_tile_num_from_latlon

  ! *******************************************************************

!BOP
!
! !ROUTINE: is_in_box
!  \label{is_in_box}
!
! !REVISION HISTORY:
!             Rolf Reichle, Initial Specification
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! determine whether catchment is within bouding box - reichle, 7 May 2003
!
! !INTERFACE:
  logical function is_in_box(                                    &
       this_minlon, this_minlat, this_maxlon, this_maxlat,       &
       minlon, minlat, maxlon, maxlat        )
    
! !ARGUMENTS: 
    implicit none
    
    real :: this_minlon, this_minlat, this_maxlon, this_maxlat
    real :: minlon, minlat, maxlon, maxlat
!
! !DESCRIPTION:
! This routine determines whether a given box is contained within
! another box.
!EOP
    
    if ( (this_minlon >= minlon) .and.        &
         (this_maxlon <= maxlon) .and.        &
         (this_minlat >= minlat) .and.        & 
         (this_maxlat <= maxlat)       )    then
       is_in_box = .true. 
    else
       is_in_box = .false.
    end if
    
  end function is_in_box

  ! *******************************************************************
  
!BOP
!
! !ROUTINE: get_ij_atm_from_latlon
!  \label{get_ij_atm_from_latlon}
!
! !REVISION HISTORY:
!             Rolf Reichle, Initial Specification
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! !INTERFACE:
  subroutine get_ij_atm_from_latlon( atm_grid, lat, lon, i_atm, j_atm )
    
    
    implicit none
    
! !ARGUMENTS: 
    type(grid_def_type), intent(in) :: atm_grid
    
    real, intent(in) :: lat, lon
    
    integer, intent(out) :: i_atm, j_atm
!
! !DESCRIPTION:
! This routine determines the column and row indices of the atmospheric
! grid-cell corresponding to the give lat/lon.
!
! NOTE: latlon is *outside* atm\_grid if 
!       ( i\_atm<0 .or. i\_atm>atm\_grid%N\_lon .or.
!         j\_atm<0 .or. j\_atm>atm\_grid%N\_lat      )
!EOP
    
    ! -------------------------------------------------------------
    !
    ! ll_lon and ll_lat refer to lower left corner of grid cell
    ! (as opposed to the grid point in the center of the grid cell)
    
    i_atm = ceiling( (lon - atm_grid%ll_lon)/atm_grid%dlon )
    j_atm = ceiling( (lat - atm_grid%ll_lat)/atm_grid%dlat )
    
    ! For a "date line on center" grid and (180-dlon/2) < lon < 180 
    ! we now have i_atm=(atm_grid%N_lon+1).
    ! If the grid spans all longitudes, reset i_atm=1
    
    if ( (i_atm==atm_grid%N_lon+1) .and. &
         (atm_grid%dlon*real(atm_grid%N_lon)>=360.-1e-4) )   i_atm = 1
    
    ! SAME SHOULD APPLY TO "pole on center" GRIDS - IMPLEMENT AT LATER TIME
    ! -reichle, 27 Sep 2005
    
  end subroutine get_ij_atm_from_latlon
  
  ! *******************************************************************
  
!BOP
!
! !ROUTINE: tile2grid
!  \label{tile2grid}
!
! !REVISION HISTORY:
!             Rolf Reichle, Initial Specification
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! reichle, 28 Jan 2005
! reichle, 16 May 2005 - added offset for "domain" grid
!
! !INTERFACE:
  subroutine tile2grid( N_tile, tile_coord, atm_grid, tile_data,    &
       grid_data, no_data_value, no_data_tol, echo,                 &
       d2g_grid_off_i, d2g_grid_off_j                            )
    
    implicit none
    
! !ARGUMENTS: 
    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), allocatable :: tile_coord  ! input
    
    type(grid_def_type), intent(in) :: atm_grid
    
    real, dimension(N_tile), intent(in) :: tile_data
    
    real, dimension(atm_grid%N_lon,atm_grid%N_lat), intent(out) :: grid_data
    
    real, intent(in), optional :: no_data_value, no_data_tol
    
    logical, intent(in), optional :: echo
    
    integer, intent(in), optional :: d2g_grid_off_i, d2g_grid_off_j

!
! !DESCRIPTION:
! map from tile space to atmospheric grid
!
! NOTE: tile\_coord must match tile\_data
!
! optional inputs: 
!   no\_data\_value :
!   no\_data\_tol   : tolerance when checking tile\_data 
!                     against no\_data\_value
!   echo          : echo no\_data\_value and tolerance
!   d2g\_grid\_off\_i, d2g\_grid\_off\_j :
!         tile\_coord%atm\_i and tile\_coord%atm\_j refer to the *global*
!         atmospheric grid that corresponds to the tile definitions
!         (as obtained from the tile\_coord\_file)
!         d2g\_grid\_off\_i and d2g\_grid\_off\_j describe the offset between
!         the global "atm\_grid\_g" and a smaller "atm\_grid\_d" for the domain
!         of interest.  With these offset arguments tile2grid
!         can be used together with an argument atm\_grid=atm\_grid\_d that 
!         is a subgrid of the global atm\_grid:
!         atm\_grid\_d%ll\_lon = atm\_grid\_g%ll\_lon + d2g\_grid\_off\_i*dlon
!         atm\_grid\_d%ll\_lat = atm\_grid\_g%ll\_lat + d2g\_grid\_off\_j*dlat
!EOP
    
    ! local variables
    
    integer :: n, i, j, off_i, off_j
    
    real :: w, no_data, tol
    
    real, parameter :: no_data_default     = -9999.
    real, parameter :: no_data_tol_default = 1e-4
    
    integer, parameter :: off_i_default = 0
    integer, parameter :: off_j_default = 0
    
    real, dimension(atm_grid%N_lon,atm_grid%N_lat) :: wgrid
    
    ! ------------------------------------
    
    if (size(tile_coord)/=N_tile) then
       write (LIS_logunit,*) 'tile2grid: tile_coord and tile_data do not match.'
       write (LIS_logunit,*) 'STOPPING.'
       stop
    end if
    
    if (present(no_data_value)) then
       no_data = no_data_value
    else
       no_data = no_data_default
    end if
    
    if (present(no_data_tol)) then
       tol = no_data_tol
    else
       tol = no_data_tol_default
    end if
    
    if (present(echo)) then
       if (echo) then
          write (LIS_logunit,*) 'tile2grid: using no-data-value = ', no_data , &
               ' with tolerance = ', tol
       end if
    end if
    
    if (present(d2g_grid_off_i)) then
       off_i = d2g_grid_off_i
    else
       off_i = off_i_default
    end if
       
    if (present(d2g_grid_off_j)) then
       off_j = d2g_grid_off_j
    else
       off_j = off_j_default
    end if
    
    ! initialize
    
    grid_data = 0.
    wgrid     = 0.
    
    ! loop through tile space
    
    do n=1,N_tile
       
       i = tile_coord(n)%i_atm - off_i
       j = tile_coord(n)%j_atm - off_j
       
       w = tile_coord(n)%frac_atm
       
       if (abs(tile_data(n)-no_data)>tol) then
          
          grid_data(i,j) = grid_data(i,j) + w*tile_data(n)
          
          wgrid(i,j) = wgrid(i,j) + w
          
       end if
       
    end do
    
    ! normalize and set no-data-value
    
    do i=1,atm_grid%N_lon
       do j=1,atm_grid%N_lat
          
          if (wgrid(i,j)>0.) then
             
             grid_data(i,j) = grid_data(i,j)/wgrid(i,j)
             
          else
             
             grid_data(i,j) = no_data
             
          end if
          
       end do
    end do
    
  end subroutine tile2grid
  
  ! *******************************************************************
  
!BOP
!
! !ROUTINE: grid2tile
!  \label{grid2tile}
!
! !REVISION HISTORY:
!             Rolf Reichle, Initial Specification
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! reichle, 28 Jan 2005
! reichle, 16 May 2005 - added offset for "domain" grid
!
! !INTERFACE:
  subroutine grid2tile( atm_grid, N_tile, tile_coord, grid_data, tile_data, &
       d2g_grid_off_i, d2g_grid_off_j                            )
    
    implicit none
! !ARGUMENTS: 
    type(grid_def_type), intent(in) :: atm_grid
    
    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), allocatable :: tile_coord  ! input
    
    real, dimension(atm_grid%N_lon,atm_grid%N_lat), intent(in) :: grid_data
    
    real, dimension(N_tile), intent(out) :: tile_data
    
    integer, intent(in), optional :: d2g_grid_off_i, d2g_grid_off_j
!
! !DESCRIPTION:
! map from atmospheric grid to tile space
!
! NOTE: tile\_coord must match tile\_data
!
! optional inputs
!
!   d2g\_grid\_off\_i, d2g\_grid\_off\_j :
!         tile\_coord%atm\_i and tile\_coord%atm\_j refer to the *global*
!         atmospheric grid that corresponds to the tile definitions
!         (as obtained from the tile\_coord\_file)
!         d2g\_grid\_off\_i and d2g\_grid\_off\_j describe the offset between
!         the global "atm\_grid\_g" and a smaller "atm\_grid\_d" for the domain
!         of interest.  With these offset arguments tile2grid
!         can be used together with an argument atm\_grid=atm\_grid\_d that 
!         is a subgrid of the global atm\_grid:
!         atm\_grid\_d%ll\_lon = atm\_grid\_g%ll\_lon + d2g\_grid\_off\_i*dlon
!         atm\_grid\_d%ll\_lat = atm\_grid\_g%ll\_lat + d2g\_grid\_off\_j*dlat
!EOP
    
    ! local variables
    
    integer :: n, i, j, off_i=0, off_j=0
    
    ! ------------------------------------
    
    if (size(tile_coord)/=N_tile) then
       write (LIS_logunit,*) 'tile2grid: tile_coord and tile_data do not match.'
       write (LIS_logunit,*) 'STOPPING.'
       stop
    end if
    
    if (present(d2g_grid_off_i))  off_i = d2g_grid_off_i
    if (present(d2g_grid_off_j))  off_j = d2g_grid_off_j
    
    do n=1,N_tile
       
       i = tile_coord(n)%i_atm - off_i
       j = tile_coord(n)%j_atm - off_j
       
       tile_data(n) = grid_data(i,j)
       
    end do
    
  end subroutine grid2tile
  
  ! *************************************************************************

end module tile_coord_routines



