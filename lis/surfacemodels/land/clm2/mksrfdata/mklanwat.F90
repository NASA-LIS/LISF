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

subroutine mklanwat (flanwat, ndiag, lake_o, swmp_o)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! make %lake and %wetland from Cogley's one degree data
! 
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: mklanwat.F90,v 1.5 2004/05/07 22:18:36 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2_varpar   !parameters
!  use clm2_varsur   !surface variables
  use fileutils, only : getfil
  use clm2_areaMod      !area averaging routines 
  use clm2_shr_sys_mod, only : shr_sys_flush 
  implicit none

! ------------------------ arguments ------------------------------
  character(len=*), intent(in) :: flanwat           !input lanwat dataset file name
  integer , intent(in) :: ndiag                     !unit number for diagnostic output
  real(r8), intent(out):: lake_o(lsmlon,lsmlat)     !percent lake on output grid
  real(r8), intent(out):: swmp_o(lsmlon,lsmlat)     !percent wetland on output grid
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  character(len=256) :: locfn              !local dataset file name

  integer :: nlon_i                        !input grid : longitude points (read in)
  integer :: nlat_i                        !input grid : latitude  points (read in)
  integer :: ncid,dimid,varid              !input netCDF id's
  integer :: ier                    !error status   

  real(r8) :: wt                           !overlap weight
  real(r8) :: w_sum                        !sum of %lake and %wetland
  real(r8) :: glake_o                      !output grid: global area lakes
  real(r8) :: gswmp_o                      !output grid: global area wetlands
  real(r8) :: garea_o                      !output grid: global area
  real(r8) :: glake_i                      !input grid: global area lakes
  real(r8) :: gswmp_i                      !input grid: global area wetlands
  real(r8) :: garea_i                      !input grid: global area
  integer :: ii                            !longitude index for COGLEY grid
  integer :: ji                            !latitude  index for COGLEY grid
  integer :: io                            !longitude index for land model grid
  integer :: jo                            !latitude  index for land model grid
  integer :: k,n                           !indices

  real(r8) :: edge_i(4)                    !input grid: N,E,S,W edges (degrees)
  real(r8), allocatable :: lake_i(:,:)     !input grid: percent lake
  real(r8), allocatable :: swmp_i(:,:)     !input grid: percent wetland
  real(r8), allocatable :: landmask_i(:,:) !input grid: fraction land (not ocn) per land gridcell
  real(r8), allocatable :: latixy_i(:,:)   !input grid: latitude (degrees)
  real(r8), allocatable :: longxy_i(:,:)   !input grid: longitude (degrees)
  integer , allocatable :: numlon_i(:)     !input grid: number longitude points by lat
  real(r8), allocatable :: lon_i(:,:)      !input grid: longitude, west edge (degrees)
  real(r8), allocatable :: lon_i_offset(:,:)!input grid: offset longitude, west edge (degrees)
  real(r8), allocatable :: lat_i(:)        !input grid: latitude, south edge (degrees)
  real(r8), allocatable :: area_i(:,:)     !input grid: cell area
  real(r8), allocatable :: mask_i(:,:)     !input grid: mask (0, 1)

  real(r8) :: mask_o                       !output grid: mask (0, 1)
  integer  :: novr_i2o                     !number of overlapping input cells
  integer  :: iovr_i2o(maxovr)             !lon index of overlap input cell
  integer  :: jovr_i2o(maxovr)             !lat index of overlap input cell
  real(r8) :: wovr_i2o(maxovr)             !weight    of overlap input cell
  real(r8) :: offset                       !used to shift x-grid 360 degrees

  real(r8) :: fld_o(lsmlon,lsmlat)         !output grid: dummy field 
  real(r8) :: fld_i                        !input grid: dummy field 
  real(r8) :: sum_fldo                     !global sum of dummy output field
  real(r8) :: sum_fldi                     !global sum of dummy input field
  real(r8) :: relerr = 0.00001             !max error: sum overlap weights ne 1
! -----------------------------------------------------------------

  write (6,*) 'Attempting to make %lake and %wetland .....'
#if 0 
  call clm2_shr_sys_flush(6)

! -----------------------------------------------------------------
! Read in Cogley's input data
! -----------------------------------------------------------------

! Obtain input grid info

  call getfil (flanwat, locfn, 0)
  call wrap_open(locfn, 0, ncid)

  call wrap_inq_dimid  (ncid, 'lon', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlon_i)

  call wrap_inq_dimid  (ncid, 'lat', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlat_i)

  allocate (lake_i(nlon_i,nlat_i), stat=ier)     
  if (ier/=0) call endrun
  allocate (swmp_i(nlon_i,nlat_i), stat=ier)     
  if (ier/=0) call endrun
  allocate (landmask_i(nlon_i,nlat_i), stat=ier)
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
  allocate (mask_i(nlon_i,nlat_i), stat=ier)     
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

  call wrap_inq_varid (ncid, 'PCT_LAKE', varid)
  call wrap_get_var_realx (ncid, varid, lake_i)

  call wrap_inq_varid (ncid, 'PCT_WETLAND', varid)
  call wrap_get_var_realx (ncid, varid, swmp_i)

! Close input file

  call wrap_close(ncid)

! -----------------------------------------------------------------
! Map data from input grid to land model grid. Get:
! -----------------------------------------------------------------

! Determine input grid cell and cell areas

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

!$OMP PARALLEL DO PRIVATE (io,jo,ii,ji,n,mask_o,novr_i2o,iovr_i2o,jovr_i2o,wovr_i2o,fld_i)
  do jo = 1, lsmlat
     do io = 1, numlon(jo)

! Determine areas of overlap and indices

        mask_o = 1.

        call areaini_point (io        , jo          , nlon_i  , nlat_i  , numlon_i, &
                           lon_i      , lon_i_offset, lat_i   , area_i  , mask_i  , &
                           lsmlon     , lsmlat      , numlon  , lonw    , lats    , &
                           area(io,jo), mask_o      , novr_i2o, iovr_i2o, jovr_i2o, &
                           wovr_i2o)                             

        mask_o = 0.
        do n = 1, novr_i2o        !overlap cell index
           ii = iovr_i2o(n)       !lon index (input grid) of overlap cell
           ji = jovr_i2o(n)       !lat index (input grid) of overlap cell
           mask_o = mask_o + landmask_i(ii,ji) * wovr_i2o(n)
        end do

        call areaini_point (io        , jo          , nlon_i  , nlat_i  , numlon_i  , &
                           lon_i      , lon_i_offset, lat_i   , area_i  , landmask_i, &
                           lsmlon     , lsmlat      , numlon  , lonw    , lats      , &
                           area(io,jo), mask_o      , novr_i2o, iovr_i2o, jovr_i2o  , &
                           wovr_i2o   ) 

! Make area average

        lake_o(io,jo) = 0.
        swmp_o(io,jo) = 0.
        do n = 1, novr_i2o   !overlap cell index
           ii = iovr_i2o(n)  !lon index (input grid) of overlap cell
           ji = jovr_i2o(n)  !lat index (input grid) of overlap cell
           lake_o(io,jo) = lake_o(io,jo) + lake_i(ii,ji) * wovr_i2o(n)
           swmp_o(io,jo) = swmp_o(io,jo) + swmp_i(ii,ji) * wovr_i2o(n)
        end do

! Corrections: set oceans to zero and exclude areas less than 5% of cell

        if (landmask(io,jo) == 0) then
           lake_o(io,jo) = 0.
           swmp_o(io,jo) = 0.
        else
           if (lake_o(io,jo) < 5.) lake_o(io,jo) = 0.
           if (swmp_o(io,jo) < 5.) swmp_o(io,jo) = 0.
        end if

! Check for conservation 

        if ((lake_o(io,jo) + swmp_o(io,jo)) > 100.000001) then
           write (6,*) 'MKLANWAT error: lake = ',lake_o(io,jo), &
                ' and wetland = ',swmp_o(io,jo), &
                ' sum are greater than 100 for lon,lat = ',io,jo
           call endrun
        end if

! Global sum of output field -- must multiply by fraction of
! output grid that is land as determined by input grid

        fld_o(io,jo) = 0.
        do n = 1, novr_i2o
           ii = iovr_i2o(n)
           ji = jovr_i2o(n)
           fld_i = ((ji-1)*nlon_i + ii) * landmask_i(ii,ji)
           fld_o(io,jo) = fld_o(io,jo) + wovr_i2o(n) * fld_i * mask_o
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
        fld_i = ((ji-1)*nlon_i + ii) * landmask_i(ii,ji)
        sum_fldi = sum_fldi + area_i(ii,ji) * fld_i
     end do
  end do

  if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
     write (6,*) 'AREAINI error: input field not conserved'
     write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
     write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
     call endrun
  end if

! -----------------------------------------------------------------
! Error check2
! Compare global areas on input and output grids
! -----------------------------------------------------------------

! Input grid

  glake_i = 0.
  gswmp_i = 0.
  garea_i = 0.

  do ji = 1, nlat_i
     do ii = 1, nlon_i
        garea_i = garea_i + area_i(ii,ji)
        glake_i = glake_i + lake_i(ii,ji)*area_i(ii,ji)/100.
        gswmp_i = gswmp_i + swmp_i(ii,ji)*area_i(ii,ji)/100.
     end do
  end do

! Output grid

  glake_o = 0.
  gswmp_o = 0.
  garea_o = 0.

  do jo = 1, lsmlat
     do io = 1, numlon(jo)
        garea_o = garea_o + area(io,jo)
        glake_o = glake_o + lake_o(io,jo)*area(io,jo)/100.
        gswmp_o = gswmp_o + swmp_o(io,jo)*area(io,jo)/100.
     end do
  end do

! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Inland Water Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
  write (ndiag,2002) glake_i*1.e-06,glake_o*1.e-06
  write (ndiag,2003) gswmp_i*1.e-06,gswmp_o*1.e-06
  write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'lakes       ',f14.3,f17.3)
2003 format (1x,'wetlands    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,a14)')area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif

  write (6,*) 'Successfully made %lake and %wetland'
  write (6,*)
  call clm2_shr_sys_flush(6)

! Deallocate dynamic memory

  deallocate (latixy_i)
  deallocate (longxy_i)
  deallocate (numlon_i)
  deallocate (lon_i)
  deallocate (lon_i_offset)
  deallocate (lat_i)
  deallocate (area_i)
  deallocate (mask_i)
  deallocate (landmask_i)
  deallocate (lake_i)
  deallocate (swmp_i)
#endif
  return
end subroutine mklanwat


