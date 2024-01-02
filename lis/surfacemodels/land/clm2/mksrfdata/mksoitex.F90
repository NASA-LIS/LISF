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

subroutine mksoitex (fsoitex, ndiag, pctgla_o, sand_o, clay_o)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! make %sand and %clay from IGBP soil data, which includes
! igbp soil 'mapunits' and their corresponding textures
! 
! Method: 
! 
! Authors: Gordon Bonan and Sam Levis
! 
!-----------------------------------------------------------------------
! $Id: mksoitex.F90,v 1.5 2004/05/07 22:18:37 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2_varpar    !parameters
  use clm2_varsur    !surface variables
  use fileutils, only : getfil
  use clm2_areaMod       !area averaging routines 
  use clm2_shr_sys_mod, only : shr_sys_flush 
  implicit none

! ------------------------ arguments ------------------------------
  character(len=*), intent(in) :: fsoitex        !soil texture dataset file name
  integer , intent(in) :: ndiag                  !unit # for diagnostic output
  real(r8), intent(in) :: pctgla_o(lsmlon,lsmlat)      !% glacier (output grid)
  real(r8), intent(out):: sand_o(lsmlon,lsmlat,nlevsoi) !% sand (output grid)
  real(r8), intent(out):: clay_o(lsmlon,lsmlat,nlevsoi) !% clay (output grid)
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  character(len=256) :: locfn            !local dataset file name

  integer :: nlon_i                      !input grid: longitude points
  integer :: nlat_i                      !input grid: latitude  points
  integer :: ncid,dimid,varid            !input netCDF id's
  integer :: ier                         !error status   

  integer :: nlay                        !number of soil layers
  integer :: mapunitmax                  !maximum value of igbp soil mapunits
  integer :: mapunittemp                 !temporary igbp soil mapunit
  integer :: ii                          !longitude index for IGBP grid
  integer :: ji                          !latitude  index for IGBP grid
  integer :: io                          !longitude index for land grid
  integer :: jo                          !latitude  index for land grid
  integer :: l,m,n                       !loop indices
  integer :: miss = 99999                !missing data indicator
  integer :: k                           !igbp soil mapunit index
  integer, parameter :: num=2            !get 1st and 2nd largest areas of overlap
  integer :: wsti(num)                   !index to 1st and 2nd largest values in wst
  integer, parameter :: nwstmax = 10000  !maximum size of overlap weights, by soil mapunit
  real(r8):: wst(0:nwstmax)              !overlap weights, by soil mapunit
  integer, parameter :: nlsm=4           !number of soil textures (sand, silt, clay)
  character(len=38) :: soil(0:nlsm)      !name of each soil texture
  character(len=38) :: typ               !soil texture based on %sand, silt, clay
  real(r8) :: gast_i(0:nlsm)             !input grid : global area, by texture type
  real(r8) :: gast_o(0:nlsm)             !output grid: global area, by texture type

  integer  :: novr_i2o                       !number of overlapping input cells
  integer  :: iovr_i2o(maxovr)               !lon index of overlap input cell
  integer  :: jovr_i2o(maxovr)               !lat index of overlap input cell
  real(r8) :: wovr_i2o(maxovr)               !weight    of overlap input cell
  real(r8) :: mask_o                         !output grid: mask (0, 1)
  real(r8) :: offset                         !used to shift x-grid 360 degrees
  real(r8) :: mapunit_o(lsmlon,lsmlat)       !not needed except as diagnostic

  real(r8) :: edge_i(4)                      !input grid: N,E,S,W edges (degrees)
  real(r8), allocatable :: sand_i(:,:)       !input grid: percent sand
  real(r8), allocatable :: clay_i(:,:)       !input grid: percent clay
  real(r8), allocatable :: landmask_i(:,:)   !input grid: land=1, no_soil_data=0
  real(r8), allocatable :: mapunit_i(:,:)    !input grid: igbp soil mapunits
  real(r8), allocatable :: latixy_i(:,:)     !input grid: latitude (degrees)
  real(r8), allocatable :: longxy_i(:,:)     !input grid: longitude (degrees)
  integer , allocatable :: numlon_i(:)       !input grid: # longitude pts by lat
  real(r8), allocatable :: lon_i(:,:)        !input grid: longitude, W edge (deg)
  real(r8), allocatable :: lon_i_offset(:,:) !input grid: offset longitude, west edge (degrees)
  real(r8), allocatable :: lat_i(:)          !input grid: latitude,  S edge (deg)
  real(r8), allocatable :: area_i(:,:)       !input grid: cell area
  real(r8), allocatable :: mask_i(:,:)   !input grid: mask (0, 1)

  real(r8) :: fld_o(lsmlon,lsmlat)           !output grid: dummy field 
  real(r8) :: fld_i                          !input grid: dummy field 
  real(r8) :: sum_fldo                       !global sum of dummy output field
  real(r8) :: sum_fldi                       !global sum of dummy input field
  real(r8) :: relerr = 0.00001               !max error: sum overlap weights ne 1
! -----------------------------------------------------------------

  write (6,*) 'Attempting to make %sand and %clay .....'
  call clm2_shr_sys_flush(6)

! -----------------------------------------------------------------
! Define the LSM surface types: 0 to nlsm
! -----------------------------------------------------------------

  soil(0) = 'no soil: ocean, glacier, lake, no data'
  soil(1) = 'clays                                 '
  soil(2) = 'sands                                 '
  soil(3) = 'loams                                 '
  soil(4) = 'silts                                 '

! -----------------------------------------------------------------
! Read input data
! -----------------------------------------------------------------

! Obtain input grid info

  call getfil (fsoitex, locfn, 0)
  call wrap_open(locfn, 0, ncid)

  call wrap_inq_dimid  (ncid, 'lon', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlon_i)

  call wrap_inq_dimid  (ncid, 'lat', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlat_i)

  call wrap_inq_dimid  (ncid, 'number_of_layers', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlay)

  call wrap_inq_dimid  (ncid, 'max_value_mapunit', dimid)
  call wrap_inq_dimlen (ncid, dimid, mapunitmax)

  !NOTE: wst must not be allocatable if it is declared private in a thread                                         
  if (nwstmax < mapunitmax) then
     write(6,*)'MKSOITEX: parameter nwstmax must be increased to ',mapunitmax
     call endrun
  endif

  allocate (sand_i(mapunitmax,nlay), stat=ier)
  if (ier/=0) call endrun
  allocate (clay_i(mapunitmax,nlay), stat=ier)
  if (ier/=0) call endrun     
  allocate (landmask_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call endrun
  allocate (mapunit_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call endrun
  allocate (mask_i(nlon_i,nlat_i), stat=ier)     
  if (ier/=0) call endrun
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

  call wrap_inq_varid (ncid, 'LATIXY', varid)
  call wrap_get_var_realx (ncid, varid, latixy_i)

  call wrap_inq_varid (ncid, 'LONGXY', varid)
  call wrap_get_var_realx (ncid, varid, longxy_i)

  call wrap_inq_varid (ncid, 'EDGEN', varid)
  call wrap_get_var_realx (ncid, varid, edge_i(1))

  call wrap_inq_varid (ncid, 'EDGEE', varid)
  call wrap_get_var_realx (ncid, varid, edge_i(2))

  call wrap_inq_varid (ncid, 'EDGES', varid)
  call wrap_get_var_realx (ncid, varid, edge_i(3))

  call wrap_inq_varid (ncid, 'EDGEW', varid)
  call wrap_get_var_realx (ncid, varid, edge_i(4))

! Obtain input data

  call wrap_inq_varid (ncid, 'LANDMASK', varid)
  call wrap_get_var_realx (ncid, varid, landmask_i)

  call wrap_inq_varid (ncid, 'MAPUNITS', varid)
  call wrap_get_var_realx (ncid, varid, mapunit_i)

  call wrap_inq_varid (ncid, 'PCT_SAND', varid)
  call wrap_get_var_realx (ncid, varid, sand_i)

  call wrap_inq_varid (ncid, 'PCT_CLAY', varid)
  call wrap_get_var_realx (ncid, varid, clay_i)

  call wrap_close(ncid)

! -----------------------------------------------------------------
! Map data from input grid to land model grid. 
! -----------------------------------------------------------------

! Determine input grid cell edges and cell areas

  numlon_i(:) = nlon_i

  call celledge (nlat_i    , nlon_i    , numlon_i  , longxy_i  ,  &
                 latixy_i  , edge_i(1) , edge_i(2) , edge_i(3) ,  &
                 edge_i(4) , lat_i     , lon_i     )

  call cellarea (nlat_i    , nlon_i    , numlon_i  , lat_i     ,  &
                 lon_i     , edge_i(1) , edge_i(2) , edge_i(3) ,  &
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

!$OMP PARALLEL DO PRIVATE (io,jo,ii,ji,n,l,k,mask_o,novr_i2o,iovr_i2o,jovr_i2o,wovr_i2o,fld_i,wst,wsti)
  do jo = 1, lsmlat
     do io = 1, numlon(jo)

! Determine areas of overlap and indices

        mask_o = 1.

        call areaini_point (io        , jo          , nlon_i  , nlat_i  , numlon_i, &
                           lon_i      , lon_i_offset, lat_i   , area_i  , mask_i  , &
                           lsmlon     , lsmlat      , numlon  , lonw    , lats    , &
                           area(io,jo), mask_o      , novr_i2o, iovr_i2o, jovr_i2o, &
                           wovr_i2o)                             

! Process each cell on land grid:
! Find dominant 5 minute x 5 minute IGBP soil mapunit.
! landmask_i=0 means no soil data is available, so assume mapunit=0.
! Map from mapunit values to corresponding %sand and %clay 
! (silt not needed by land model and therefore, not calculated).
! Sum overlap weights by igbp soil mapunit
! landmask_i=0 means no soil data and landmask_i=1 means data exist

        do k = 0, mapunitmax
           wst(k) = 0.
        end do
        do n = 1, novr_i2o      !overlap cell index
           ii = iovr_i2o(n)     !lon index (input grid) of overlap cell
           ji = jovr_i2o(n)     !lat index (input grid) of overlap cell
           k = mapunit_i(ii,ji) * landmask_i(ii,ji) !mapunit (input grid)
           wst(k) = wst(k) + wovr_i2o(n)
        end do
        
! Rank non-zero weights by soil mapunit. 
! wsti(1) is the most extensive mapunit. 
! wsti(2) is the second most extensive mapunit.

        call mkrank (mapunitmax, wst, miss, wsti, num)

! Set soil texture as follows:
! If land grid cell is ocean or 100% glacier: cell has no soil 
! Otherwise, grid cell needs soil:
!   a. Use dominant igbp soil mapunit based on area of overlap unless
!     'no data' is dominant
!   b. In this case use second most dominant mapunit so long as it has data
!   c. If this has no data or if there isn't a second most dominant
!      mapunit, use loam for soil texture

        if (landmask(io,jo) == 0) then                        !ocean
           mapunit_o(io,jo) = 0.
           do l = 1, nlay
              sand_o(io,jo,l) = 0.
              clay_o(io,jo,l) = 0.
           end do
        else if (abs(pctgla_o(io,jo)-100.) < 1.e-06) then     !glacier
           mapunit_o(io,jo) = 0.
           do l = 1, nlay
              sand_o(io,jo,l) = 0.
              clay_o(io,jo,l) = 0.
           end do
        else                                                  !need soil
           if (wsti(1) /= 0) then                             !not 'no data'
              mapunit_o(io,jo) = wsti(1)
              do l = 1, nlay
                 sand_o(io,jo,l) = sand_i(wsti(1),l)
                 clay_o(io,jo,l) = clay_i(wsti(1),l)
              end do
           else                                               !if (wsti(1) == 0) then
              if (wsti(2) == 0 .or. wsti(2) == miss) then     !no data
                 mapunit_o(io,jo) = wsti(2)
                 do l = 1, nlay
                    sand_o(io,jo,l) = 43.                     !use loam
                    clay_o(io,jo,l) = 18.
                 end do
              else                                            !if (wsti(2) /= 0 and /= miss)
                 mapunit_o(io,jo) = wsti(2)
                 do l = 1, nlay
                    sand_o(io,jo,l) = sand_i(wsti(2),l)
                    clay_o(io,jo,l) = clay_i(wsti(2),l)
                 end do
              end if       !end of wsti(2) if-block
           end if          !end of wsti(1) if-block
        end if             !end of land/ocean if-block

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
     write (6,*) 'MKSOITEX error: input field not conserved'
     write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
     write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
     call endrun
  end if

! -----------------------------------------------------------------
! Error check2
! Compare global area of each soil type on input and output grids
! -----------------------------------------------------------------

! input grid: global areas by texture class

  gast_i(:) = 0.
  do l = 1, nlay
     do ji = 1, nlat_i
        do ii = 1, nlon_i
           mapunittemp = nint(mapunit_i(ii,ji))
           if (mapunittemp==0) then
              typ = 'no soil: ocean, glacier, lake, no data'
           else if (clay_i(mapunittemp,l) >= 40.) then
              typ = 'clays'
           else if (sand_i(mapunittemp,l) >= 50.) then
              typ = 'sands'
           else if (clay_i(mapunittemp,l)+sand_i(mapunittemp,l) < 50.) then
              if (landmask_i(ii,ji) /= 0.) then
                 typ = 'silts'
              else            !if (landmask_i(ii,ji) == 0.) then !no data
                 typ = 'no soil: ocean, glacier, lake, no data'
              end if
           else
              typ = 'loams'
           end if
           do m = 0, nlsm
              if (typ == soil(m)) go to 101
           end do
           write (6,*) 'MKSOITEX error: sand = ',sand_i(mapunittemp,l), &
             ' clay = ',clay_i(mapunittemp,l), &
             ' not assigned to soil type for input grid lon,lat,layer = ',ii,ji,l
           call endrun
101        continue
           gast_i(m) = gast_i(m) + area_i(ii,ji)
        end do
     end do
  end do

! output grid: global areas by texture class

  gast_o(:) = 0.
  do l = 1, nlay
     do jo = 1, lsmlat
        do io = 1, numlon(jo)
           if (clay_o(io,jo,l)==0. .and. sand_o(io,jo,l)==0.) then
              typ = 'no soil: ocean, glacier, lake, no data'
           else if (clay_o(io,jo,l) >= 40.) then
              typ = 'clays'
           else if (sand_o(io,jo,l) >= 50.) then
              typ = 'sands'
           else if (clay_o(io,jo,l)+sand_o(io,jo,l) < 50.) then
              typ = 'silts'
           else
              typ = 'loams'
           end if
           do m = 0, nlsm
              if (typ == soil(m)) go to 102
           end do
           write (6,*) 'MKSOITEX error: sand = ',sand_o(io,jo,l), &
             ' clay = ',clay_o(io,jo,l), &
             ' not assigned to soil type for output grid lon,lat,layer = ',io,jo,l
           call endrun
102        continue
           gast_o(m) = gast_o(m) + area(io,jo)
        end do
     end do
  end do

! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',l=1,70)
  write (ndiag,*) 'Soil Texture Output'
  write (ndiag,'(1x,70a1)') ('=',l=1,70)
  write (ndiag,*)

  write (ndiag,*) 'The following table of soil texture classes is for comparison only.'
  write (ndiag,*) 'The actual data is continuous %sand, %silt and %clay not textural classes'
  write (ndiag,*)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',l=1,70)
  write (ndiag,1001)
1001 format (1x,'soil texture class',17x,' input grid area output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
  write (ndiag,'(1x,70a1)') ('.',l=1,70)
  write (ndiag,*)

  do l = 0, nlsm
     write (ndiag,1002) soil(l),gast_i(l)*1.e-6,gast_o(l)*1.e-6
1002 format (1x,a38,f16.3,f17.3)
  end do

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,a14)')area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif
     
  write (6,*) 'Successfully made %sand and %clay'
  write (6,*)
  call clm2_shr_sys_flush(6)

! Deallocate dynamic memory

  deallocate (mask_i)
  deallocate (latixy_i)
  deallocate (longxy_i)
  deallocate (numlon_i)
  deallocate (lon_i)
  deallocate (lon_i_offset)
  deallocate (lat_i)
  deallocate (area_i)
  deallocate (sand_i)
  deallocate (clay_i)
  deallocate (landmask_i)
  deallocate (mapunit_i)

  return
end subroutine mksoitex
