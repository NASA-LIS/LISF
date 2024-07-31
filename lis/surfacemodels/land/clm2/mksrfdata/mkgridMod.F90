!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "preproc.h"
#include "LIS_misc.h"

module mkgridMod

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Routines to create land model grid
! 
! Method: 
!
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: mkgridMod.F90,v 1.5 2004/05/07 22:18:36 jim Exp $
!-----------------------------------------------------------------------


  use LIS_precisionMod
  use clm2_varpar          !parameters 
!  use clm2_varsur          !surface variables 
  use clm_varctl          !run control variables
  use fileutils, only : getfil
  use clm2_areaMod             !area averaging routines
  use spmdMod, only: LIS_masterproc
  implicit none

!=======================================================================
contains
!=======================================================================

#if (defined OFFLINE) 

  subroutine mkgrid_offline

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Obtain land model grid
! 
! Method: 
! If namelist variable mksrf_offline_fgrid is the empty string, then
! the corresponding file will be used to determine the land model grid
! If namelist variable mksrf_offline is not the empty string, then 
! the land model grid will be generated at run time 
!
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
#if 0 
    if (mksrf_offline_fgrid /= ' ') then
       call read_grid_offline
       offline_rdgrid = .true.
    else
       call create_grid_offline
       offline_rdgrid = .false.
    endif
#endif
  end subroutine mkgrid_offline

!=======================================================================

  subroutine read_grid_offline

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read land model grid when mode is offline.
! 
! Method: 
! If namelist variable mksrf_offline_fgrid is the empty string, then
! the corresponding file will be used to determine the land model grid
! Assume that the input data file has the land grid in the following 
! form:
!    lon                    => dimension
!    lat                    => dimension
!    lon(lsmlon)            => full grid longitudes
!    nlon(lsmlat)           => reduced grid number of lats per lon
!    rlon(lsmlon,lsmlat)    => reduced grid longitudes
!    lat(lsmlat)            => grid latitudes
!    LANDMASK(lsmlon,lsmlat)=> 2d land mask
!
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    include 'netcdf.inc'

! ------------------------ local variables ------------------------
    integer  :: i,j,k,n                 !indices
    integer  :: nlon_i                  !number of input data longitudes
    integer  :: nlat_i                  !number of input data latitudes
    integer  :: ncid                    !netCDF file id 
    integer  :: dimid                   !netCDF dimension id
    integer  :: varid                   !netCDF variable id
    integer  :: ret                     !netCDF return code
    real(r8) :: lon(lsmlon)             !input longitude array (full grid)
    real(r8) :: lat(lsmlat)             !input latitude array (full grid)
    real(r8) :: oro(lsmlon,lsmlat)      !input oro field 
    character(len=256) :: locfn         !local file name
!-----------------------------------------------------------------------

    write (6,*) 'Attempting to read land grid data .....'
#if 0 
    write (6,'(72a1)') ("-",i=1,60)
    
    call getfil (mksrf_offline_fgrid, locfn, 0)
    call wrap_open(locfn, 0, ncid)
    
    call wrap_inq_dimid  (ncid, 'lon', dimid)
    call wrap_inq_dimlen (ncid, dimid, nlon_i)
    if (nlon_i /= lsmlon) then
       write(6,*)'RDGRID_OFFLINE: parameter lsmlon= ',lsmlon, &
            'does not equal input nlon_i= ',nlon_i
       call endrun
    endif
    
    call wrap_inq_dimid  (ncid, 'lat', dimid)
    call wrap_inq_dimlen (ncid, dimid, nlat_i)
    if (nlat_i /= lsmlat) then
       write(6,*)'RDGRID_OFFLINE: parameter lsmlat= ',lsmlat, &
            'does not equal input nlat_i= ',nlat_i
       call endrun
    endif
    
! Determine grid longitudes for either full or reduced grid
! if variable 'rlon' is not on grid file then have full grid
! if variable 'rlon' is on grid file then have reduced grid

    ret = nf_inq_varid (ncid, 'rlon', dimid)
    if (ret == NF_NOERR) then
       fullgrid = .false.
    else
       fullgrid = .true.
    endif

    if (fullgrid) then
       numlon(:) = lsmlon
       call wrap_inq_varid (ncid, 'lon' , varid)
       call wrap_get_var_realx (ncid, varid, lon)
       do j = 1,lsmlat
          do i = 1,lsmlon
             longxy(i,j) = lon(i)
          end do
       end do
    else
       call wrap_inq_varid (ncid, 'nlon' , varid)
       call wrap_get_var_int (ncid, varid, numlon)
       call wrap_inq_varid (ncid, 'rlon' , varid)
       call wrap_get_var_realx (ncid, varid, longxy)
    endif

! Determine grid latitudes

    call wrap_inq_varid (ncid, 'lat' , varid)
    call wrap_get_var_realx (ncid, varid, lat)
    do j = 1,lsmlat
       do i =1,lsmlon
          latixy(i,j) = lat(j)
       end do
    end do

! Define land grid edges and grid cell areas

    call celledge (lsmlat, lsmlon, numlon, longxy, latixy, &
                   lats  , lonw  )

    call cellarea (lsmlat, lsmlon, numlon, lats, lonw, &
                   area   )

! Determine land mask

    call wrap_inq_varid (ncid, 'ORO' , varid)
    call wrap_get_var_realx (ncid, varid, oro)

    do j = 1,lsmlat
       do i = 1,numlon(j)
          if (nint(oro(i,j)) == 1) then
             landmask(i,j) = 1
             landfrac(i,j) = 1.0
          else
             landmask(i,j) = 0
             landfrac(i,j) = 0.
          endif
       end do
    end do

    write (6,'(72a1)') ("-",i=1,60)
    write (6,*) 'Successfully read land grid data'
    write (6,*)
#endif
  end subroutine read_grid_offline

!=======================================================================

  subroutine create_grid_offline

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Generate land model grid when mode is offline.
! 
! Method: 
! Surface grid edges -- Grids do not have to be global. To allow this, grids 
! must define the north, east, south, and west edges:
!
! If namelist variable mksrf_offline is not the empty string, then the land model
! grid will be generated at run time using the settings of the
! namelist variables
!    o mksrf_offline_fnavyoro : 20 min navy orography dataset
!    o mksrf_offline_edgen (edge(1)) : northern edge of grid (degrees): >  -90 and <= 90
!    o mksrf_offline_edgee (edge(2)) : eastern edge of grid (degrees) : see following notes
!    o mksrf_offline_edges (edge(3)) : southern edge of grid (degrees): >= -90 and <  90
!    o mksrf_offline_edgew (edge(4)) : western edge of grid (degrees) : see following notes
!
! For partial grids, northern and southern edges are any latitude
! between 90 (North Pole) and -90 (South Pole). Western and eastern
! edges are any longitude between -180 and 180, with longitudes
! west of Greenwich negative. That is, western edge >= -180 and < 180;
! eastern edge > western edge and <= 180.
! 
! For global grids, northern and southern edges are 90 (North Pole)
! and -90 (South Pole). The western edge of the longitude grid starts 
! at the dateline if the grid is generated (the longitudes for each grid 
! cell correspond with the edges (i.e., range from -180 to 180)). 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

! ------------------------ local variables ------------------------
    character(len=256) :: locfn                !local file name

    integer  :: i,j,k,n                        !indices
    integer  :: ii,ji,io,jo                    !indices
    integer  :: ncid                           !netCDF file id 
    integer  :: dimid                          !netCDF dimension id
    integer  :: varid                          !netCDF variable id
    integer  :: ier                            !error status   

    integer  :: nlon_i                         !input number of longitudes
    integer  :: nlat_i                         !input number of latitudes
    real(r8) :: dx                             !land model cell width
    real(r8) :: dy                             !land model cell length
    real(r8) :: edge_i(4)                      !input grid: N,E,S,W edges (degrees)
    real(r8), allocatable :: latixy_i(:,:)     !input grid: latitude (degrees)
    real(r8), allocatable :: longxy_i(:,:)     !input grid: longitude (degrees)
    integer , allocatable :: numlon_i(:)       !input grid: number longitude points by lat
    real(r8), allocatable :: lon_i(:,:)        !input grid: longitude, west edge (degrees)
    real(r8), allocatable :: lon_i_offset(:,:) !input grid: longitude, west edge (degrees)
    real(r8), allocatable :: lat_i(:)          !input grid: latitude, south edge (degrees)
    real(r8), allocatable :: area_i(:,:)       !input grid: cell area
    real(r8), allocatable :: mask_i(:,:)       !input grid: mask (0, 1)
    real(r8), allocatable :: fland_i(:,:)      !input grid: fractional land

    real(r8) :: mask_o                         !output grid: mask (0, 1)
    integer  :: novr_i2o                       !number of overlapping input cells
    integer  :: iovr_i2o(maxovr)               !lon index of overlap input cell
    integer  :: jovr_i2o(maxovr)               !lat index of overlap input cell
    real(r8) :: wovr_i2o(maxovr)               !weight    of overlap input cell
    real(r8) :: offset                         !used to shift x-grid 360 degrees
    
    real(r8) :: fld_o(lsmlon,lsmlat)           !output grid: dummy field 
    real(r8) :: fld_i                          !input grid: dummy field 
    real(r8) :: sum_fldo                       !global sum of dummy output field
    real(r8) :: sum_fldi                       !global sum of dummy input field
    real(r8) :: relerr = 0.00001               !max error: sum overlap weights ne 1
    
    real(r8) :: flandmin = 0.50                !minimum land fraction for grid cell to be called land  
! -----------------------------------------------------------------

! Set numlon to uniform grid (offline ASSUMES that never have a reduced grid)
#if 0 
    numlon(:) = lsmlon

! Determine lsm grid edges

    lsmedge(1) = mksrf_offline_edgen
    lsmedge(2) = mksrf_offline_edgee
    lsmedge(3) = mksrf_offline_edges
    lsmedge(4) = mksrf_offline_edgew

! Determine grid longitudes and latitudes in increments of dx and dy
! Global latitude grid goes from south pole to north pole
! Global longitude grid starts at Dateline with western edge on Dateline

    dx = (mksrf_offline_edgee - mksrf_offline_edgew) / lsmlon
    dy = (mksrf_offline_edgen - mksrf_offline_edges) / lsmlat
    do j = 1, lsmlat
       do i = 1, lsmlon
          longxy(i,j) = mksrf_offline_edgew + (2*i-1)*dx / 2.
          latixy(i,j) = mksrf_offline_edges + (2*j-1)*dy / 2
       end do
    end do

! Define edges and area of land model grid cells

    call celledge (lsmlat     , lsmlon     , numlon     , longxy     ,&
                   latixy     , lsmedge(1) , lsmedge(2) , lsmedge(3) ,&
                   lsmedge(4) , lats       , lonw       )

    call cellarea (lsmlat     , lsmlon     , numlon     , lats       ,&
                   lonw       , lsmedge(1) , lsmedge(2) , lsmedge(3) ,&
                   lsmedge(4) , area       )

! read navy orography data and obtain fractional land

    call getfil (mksrf_offline_fnavyoro, locfn, 0)
    call wrap_open(locfn, 0, ncid)

    call wrap_inq_dimid  (ncid, 'lon', dimid)
    call wrap_inq_dimlen (ncid, dimid, nlon_i)

    call wrap_inq_dimid  (ncid, 'lat', dimid)
    call wrap_inq_dimlen (ncid, dimid, nlat_i)

    allocate (latixy_i(nlon_i,nlat_i), stat=ier)   
    if (ier/=0) call endrun
    allocate (longxy_i(nlon_i,nlat_i), stat=ier)   
    if (ier/=0) call endrun
    allocate (numlon_i(nlat_i), stat=ier)
    if (ier/=0) call endrun
    allocate (lon_i(nlon_i+1,nlat_i), stat=ier)
    if (ier/=0) call endrun
    allocate (lon_i_offset(nlon_i+1,nlat_i), stat=ier)
    if (ier/=0) call endrun
    allocate (lat_i(nlat_i+1), stat=ier)        
    if (ier/=0) call endrun
    allocate (area_i(nlon_i,nlat_i), stat=ier)  
    if (ier/=0) call endrun
    allocate (mask_i(nlon_i,nlat_i), stat=ier)     
    if (ier/=0) call endrun
    allocate (fland_i(nlon_i,nlat_i), stat=ier)     
    if (ier/=0) call endrun

    call wrap_inq_varid (ncid, 'LATIXY', varid)
    call wrap_get_var_realx (ncid, varid, latixy_i)

    call wrap_inq_varid (ncid, 'LONGXY', varid)
    call wrap_get_var_realx (ncid, varid, longxy_i)

    call wrap_inq_varid (ncid, 'NUMLON', varid)
    call wrap_get_var_int (ncid, varid, numlon_i)

    call wrap_inq_varid (ncid, 'EDGEN', varid)
    call wrap_get_var_realx (ncid, varid, edge_i(1))

    call wrap_inq_varid (ncid, 'EDGEE', varid)
    call wrap_get_var_realx (ncid, varid, edge_i(2))

    call wrap_inq_varid (ncid, 'EDGES', varid)
    call wrap_get_var_realx (ncid, varid, edge_i(3))

    call wrap_inq_varid (ncid, 'EDGEW', varid)
    call wrap_get_var_realx (ncid, varid, edge_i(4))

    call wrap_inq_varid (ncid, 'LANDFRAC', varid)
    call wrap_get_var_realx (ncid, varid, fland_i)

    call wrap_close(ncid)

! determine maximim overlap for hight resolution and lsm grids and 
! and allocate dynamic memory for overlap arrays
! first determine input grid cell and cell areas

    numlon_i(:) = nlon_i

    call celledge (nlat_i    , nlon_i    , numlon_i  , longxy_i  , &
                   latixy_i  , edge_i(1) , edge_i(2) , edge_i(3) , &
                   edge_i(4) , lat_i     , lon_i     )

    call cellarea (nlat_i    , nlon_i    , numlon_i  , lat_i     , &
                   lon_i     , edge_i(1) , edge_i(2) , edge_i(3) , &
                   edge_i(4) , area_i    )

    do ji = 1, nlat_i
       do ii = 1, numlon_i(ji)
          mask_i(ii,ji) = 1.
       end do
    end do

! Shift x-grid to locate periodic grid intersections. This
! assumes that all lon_i(1,j) have the same value for all
! latitudes j and that the same holds for lon_o(1,j)

  if (lon_i(1,1) < lonw(1,1)) then
     offset = 360.0
  else
     offset = -360.0
  end if
  
  do ji = 1, nlat_i
     do ii = 1, numlon_i(ji) + 1
        lon_i_offset(ii,ji) = lon_i(ii,ji) + offset
     end do
  end do

! Process each cell on land model grid
! novr_i2o - number of input grid cells that overlap each land grid cell
! iovr_i2o - longitude index of overlapping input grid cell
! jovr_i2o - latitude  index of overlapping input grid cell
! wovr_i2o - fraction of land grid cell overlapped by input grid cell

!$OMP PARALLEL DO PRIVATE (io,jo,ii,ji,n,mask_o,novr_i2o,iovr_i2o,jovr_i2o,wovr_i2o,fld_i)
  do jo = 1, lsmlat
     do io = 1, numlon(jo)

! Determine areas of overlap and indices

        mask_o = 1.

        call areaini_point (io        , jo          , nlon_i  , nlat_i  , numlon_i , &
                           lon_i      , lon_i_offset, lat_i   , area_i  , mask_i   , &
                           lsmlon     , lsmlat      , numlon  , lonw    , lats     , &
                           area(io,jo), mask_o      , novr_i2o, iovr_i2o, jovr_i2o , &
                           wovr_i2o)                             

! Make area average

        landfrac(io,jo) = 0.
        do n = 1, novr_i2o   !overlap cell index
           ii = iovr_i2o(n)  !lon index (input grid) of overlap cell
           ji = jovr_i2o(n)  !lat index (input grid) of overlap cell
           landfrac(io,jo) = landfrac(io,jo) + fland_i(ii,ji) * wovr_i2o(n)
        end do
        if (landfrac(io,jo) > 100.000001) then
           write (6,*) 'MKGRID error: fland = ',landfrac(io,jo), &
                ' is greater than 100 for lon,lat = ',io,jo
           call endrun
        end if

! Global sum of output field -- must multiply by fraction of
! output grid that is land as determined by input grid

        fld_o(io,jo) = 0.
        do n = 1, novr_i2o
           ii = iovr_i2o(n)
           ji = jovr_i2o(n)
           fld_i = ((ji-1)*nlon_i + ii) 
           fld_o(io,jo) = fld_o(io,jo) + wovr_i2o(n) * fld_i 
        end do

     end do  !end of output longitude loop
  end do     !end of output latitude  loop
!$OMP END PARALLEL DO

! -----------------------------------------------------------------
! Error check1
! Compare global sum fld_o to global sum fld_i. 
! -----------------------------------------------------------------

! This check is true only if both grids span the same domain. 
! To obtain global sum of input field must multiply by 
! fraction of input grid that is land as determined by input grid

  sum_fldo = 0.
  do jo = 1,lsmlat
     do io = 1,numlon(jo)
        sum_fldo = sum_fldo + area(io,jo) * fld_o(io,jo) 
     end do
  end do

  sum_fldi = 0.
  do ji = 1, nlat_i      
     do ii = 1, numlon_i(ji)
        fld_i = ((ji-1)*nlon_i + ii) 
        sum_fldi = sum_fldi + area_i(ii,ji) * fld_i
     end do
  end do

  if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
     write (6,*) 'MKGRID error: input field not conserved'
     write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
     write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
     call endrun
  end if

! Determine land mask 

    where (landfrac(:,:) < flandmin)
       landmask(:,:) = 0     !ocean
    elsewhere
       landmask(:,:) = 1     !land
    endwhere

! deallocate dynamic memory

    deallocate (numlon_i)
    deallocate (latixy_i)
    deallocate (longxy_i)
    deallocate (lon_i)
    deallocate (lon_i_offset)
    deallocate (lat_i)
    deallocate (area_i)
    deallocate (mask_i)
    deallocate (fland_i)

    write (6,*) 'Successfully makde land grid data'
    write (6,*)
#endif
  end subroutine create_grid_offline

#endif

!=======================================================================

#if (defined COUP_CAM || defined COUP_CSM)

  subroutine mkgrid_cam(cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Generate land model grid when mode is CAM or CSM
! For CAM mode get grid AND fractional land and land mask from cam model
! For CSM mode get grid AND fractional land and land mask from coupler
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

! ------------------------ arguments -----------------------------------
    real(r8), intent(in) :: cam_longxy(:,:)      !cam longitudes 
    integer , intent(in) :: cam_numlon(:)        !number of cam longitudes per latitude
    real(r8), intent(in) :: cam_latixy(:,:)      !cam latitudes  
    real(r8), intent(in) :: cam_landfrac(:,:)    !cam land fraction
    integer , intent(in) :: cam_landmask(:,:)    !cam land mask
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
    integer i,j     ! loop indices
! -----------------------------------------------------------------

! Determine if grid has pole points - if so, make sure that north pole
! is non-land and south pole is land
#if 0 
    if (abs((cam_latixy(1,lsmlat) - 90.)) < 1.e-6) then
       pole_points = .true.
       if ( LIS_masterproc ) write(6,*)'MKGRIDMOD: model has pole_points' 
       do i = 1,cam_numlon(1)
          if (cam_landmask(i,1) /= 1) then
             write(6,*)'cam grid with pole points has non-land at south pole'
             write(6,*)'longitude index= ',i
             call endrun
          endif
       end do
       do i = 1,cam_numlon(lsmlat)
          if (cam_landmask(i,lsmlat) == 1) then
             write(6,*)'cam grid with pole points has land at north pole'
             write(6,*)'longitude index= ',i
             call endrun
          endif
       end do
    else
       pole_points = .false.
       if ( LIS_masterproc ) write(6,*)'MKGRIDMOD: model does not have pole_points' 
    endif
          
! Determine land grid, land mask and land fraction

    numlon(:) = cam_numlon(:)  

    fullgrid = .true.
    do j = 1,lsmlat
       if (cam_numlon(j) < lsmlon) fullgrid = .false.
    end do

    do j = 1,lsmlat
       do i = 1,numlon(j)
          longxy(i,j)   = cam_longxy(i,j) 
          latixy(i,j)   = cam_latixy(i,j)
          landmask(i,j) = cam_landmask(i,j)
          landfrac(i,j) = cam_landfrac(i,j)
       end do
    end do

! Define land grid edges and grid cell areas

    call celledge (lsmlat, lsmlon, numlon, longxy, latixy, &
                   lats  , lonw  )

    call cellarea (lsmlat, lsmlon, numlon, lats, lonw, &
                   area   )

! Determine land fraction and land mask
#endif
    return
  end subroutine mkgrid_cam

#endif

end module mkgridMod
