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

subroutine mkpft (fpft, ndiag,  noveg,  pctlnd_o, pft, pctpft)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Make PFT data for vegetated patches (1 to maxpatch_pft)
! 
! Method: 
! This dataset consists of the %cover of the [numpft]+1 PFTs used by
! the model. The %cover pertains to the "vegetated" portion of the
! grid cell and sums to 100. The real portion of each grid cell
! covered by each PFT is the PFT cover times the fraction of the
! grid cell that is land. This is the quantity preserved when
! area-averaging from the input (1/2 degree) grid to the models grid.
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: mkpft.F90,v 1.5 2004/05/07 22:18:36 jim Exp $
!-----------------------------------------------------------------------

  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_precisionMod
  use clm2_varpar    !lsm parameters
  use clm2_varsur    !lsm surface variables
  use fileutils, only : getfil
  use clm2_areaMod      !area averaging routines 
  use clm2_shr_sys_mod, only : shr_sys_flush 
  implicit none

! ------------------------ arguments ------------------------------
  character(len=*), intent(in) :: fpft                        !input pft dataset file name
  integer , intent(in) :: ndiag                               !unit number for diagnostic output
  integer , intent(in) :: noveg                               !PFT number for no vegetation
  real(r8), intent(out):: pctlnd_o(lsmlon,lsmlat)             !output grid: % land per gridcell
  integer , intent(out):: pft(lsmlon,lsmlat,maxpatch_pft)     !output grid PFT (0 to numpft)
  real(r8), intent(out):: pctpft(lsmlon,lsmlat,maxpatch_pft)  !output grid PFT cover (% of vegetated area)
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  character(len=LIS_CONST_PATH_LEN) locfn  !local dataset file name

  integer :: nlon_i                 !input grid : longitude points (read in)
  integer :: nlat_i                 !input grid : latitude  points (read in)
  integer :: ncid,dimid,varid       !input netCDF id's
  integer :: ier                    !error status   

  integer  :: miss = 99999                    !missing data indicator
  real(r8) :: pctpft_o(lsmlon,lsmlat,0:numpft)!PFT percent on output grid
  real(r8) :: wst(0:numpft)                   !as pft_o at specific io, jo
  integer  :: wsti(maxpatch_pft)              !ranked indices of largest values in wst
  real(r8) :: wst_sum                         !sum of %pft
  real(r8) :: sumpct                          !sum of %pft over maxpatch_pft
  real(r8) :: diff                            !the difference (wst_sum - sumpct)
  real(r8) :: gpft_o(0:numpft)                !output grid: global area pfts
  real(r8) :: garea_o                         !output grid: global area
  real(r8) :: gpft_i(0:numpft)                !input grid: global area pfts
  real(r8) :: garea_i                         !input grid: global area
  integer  :: ii                              !longitude index for input grid
  integer  :: ji                              !latitude  index for input grid
  integer  :: io                              !longitude index for LSM grid
  integer  :: jo                              !latitude  index for LSM grid
  integer  :: k,n,m                           !indices

  integer numpft_i                            !number of plant types on input dataset
  real(r8) :: edge_i(4)                       !input grid: N,E,S,W edges (degrees)
  real(r8), allocatable :: pctpft_i(:,:,:)    !input grid: PFT percent 
  real(r8), allocatable :: landmask_i(:,:)    !input grid: fraction land (not ocn) per land gridcell
  real(r8), allocatable :: mask_i(:,:)        !input grid: mask (0, 1)
  real(r8), allocatable :: latixy_i(:,:)      !input grid: latitude (degrees)
  real(r8), allocatable :: longxy_i(:,:)      !input grid: longitude (degrees)
  integer , allocatable :: numlon_i(:)        !input grid: number longitude points by lat
  real(r8), allocatable :: lon_i(:,:)         !input grid: longitude, west edge (degrees)
  real(r8), allocatable :: lon_i_offset(:,:)!input grid: offset longitude, west edge (degrees)
  real(r8), allocatable :: lat_i(:)           !input grid: latitude, south edge (degrees)
  real(r8), allocatable :: area_i(:,:)        !input grid: cell area

  real(r8) :: mask_o                          !output grid: mask (0, 1)
  integer  :: novr_i2o                        !number of overlapping input cells
  integer  :: iovr_i2o(maxovr)                !lon index of overlap input cell
  integer  :: jovr_i2o(maxovr)                !lat index of overlap input cell
  real(r8) :: wovr_i2o(maxovr)                !weight    of overlap input cell
  real(r8) :: offset                          !used to shift x-grid 360 degrees

  real(r8) :: fld_o(lsmlon,lsmlat)            !output grid: dummy field 
  real(r8) :: fld_i                           !input grid: dummy field 
  real(r8) :: sum_fldo                        !global sum of dummy output field
  real(r8) :: sum_fldi                        !global sum of dummy input field
  real(r8) :: relerr = 0.00001                !max error: sum overlap weights ne 1

  character(len=35)  veg(0:numpft)            !vegetation types
  data veg( 0) /'not vegetated'                      /
  data veg( 1) /'needleleaf evergreen temperate tree'/
  data veg( 2) /'needleleaf evergreen boreal tree'   /
  data veg( 3) /'needleleaf deciduous boreal tree'   /
  data veg( 4) /'broadleaf evergreen tropical tree'  /
  data veg( 5) /'broadleaf evergreen temperate tree' /
  data veg( 6) /'broadleaf deciduous tropical tree'  /
  data veg( 7) /'broadleaf deciduous temperate tree' /
  data veg( 8) /'broadleaf deciduous boreal tree'    /
  data veg( 9) /'broadleaf evergreen shrub'          /
  data veg(10) /'broadleaf deciduous temperate shrub'/
  data veg(11) /'broadleaf deciduous boreal shrub'   /
  data veg(12) /'c3 arctic grass'                    /
  data veg(13) /'c3 non-arctic grass'                /
  data veg(14) /'c4 grass'                           / 
  data veg(15) /'corn'                               /
  data veg(16) /'wheat'                              /
! -----------------------------------------------------------------

  write (6,*) 'Attempting to make PFTs .....'
  call clm2_shr_sys_flush(6)

! -----------------------------------------------------------------
! Read input PFT file
! -----------------------------------------------------------------

! Obtain input grid info

  call getfil (fpft, locfn, 0)
  call wrap_open(locfn, 0, ncid)

  call wrap_inq_dimid  (ncid, 'lon', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlon_i)

  call wrap_inq_dimid  (ncid, 'lat', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlat_i)

  call wrap_inq_dimid  (ncid, 'pft', dimid)
  call wrap_inq_dimlen (ncid, dimid, numpft_i)

  if (numpft_i .ne. numpft+1) then
     write(6,*)'MKPFT: parameter numpft+1= ',numpft+1, &
          'does not equal input dataset numpft= ',numpft_i
     call endrun
  endif

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
  allocate (landmask_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call endrun
  allocate (mask_i(nlon_i,nlat_i), stat=ier)     
  if (ier/=0) call endrun
  allocate (area_i(nlon_i,nlat_i), stat=ier)  
  if (ier/=0) call endrun
  allocate (pctpft_i(nlon_i,nlat_i,0:numpft), stat=ier)     
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

  call wrap_inq_varid (ncid, 'PCT_PFT', varid)
  call wrap_get_var_realx (ncid, varid, pctpft_i)

  call wrap_close(ncid)

! -----------------------------------------------------------------
! Map data from input grid to land model grid. Get:
! -----------------------------------------------------------------

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

!$OMP PARALLEL DO PRIVATE (io,jo,ii,ji,n,m,mask_o,novr_i2o,iovr_i2o,jovr_i2o,wovr_i2o,fld_i, &
!$OMP &  wst, wst_sum, sumpct, wsti, diff)
  do jo = 1, lsmlat
     do io = 1, numlon(jo)

! Determine areas of overlap and indices

        mask_o = 1.

        call areaini_point (io        , jo          , nlon_i  , nlat_i  , numlon_i , &
                           lon_i      , lon_i_offset, lat_i   , area_i  , mask_i   , &
                           lsmlon     , lsmlat      , numlon  , lonw    , lats     , &
                           area(io,jo), mask_o      , novr_i2o, iovr_i2o, jovr_i2o , &
                           wovr_i2o)                             

! Make percent of grid cell that is land [pctlnd_o] and map data 
! from input grid to land model grid again, this time with the real mask_i

        pctlnd_o(io,jo) = 0.   
        do n = 1, novr_i2o  !overlap cell index
           ii = iovr_i2o(n) !lon index (input grid) of overlap cell
           ji = jovr_i2o(n) !lat index (input grid) of overlap cell
           pctlnd_o(io,jo) = pctlnd_o(io,jo) + 100.*landmask_i(ii,ji)*wovr_i2o(n)
        end do
        mask_o = pctlnd_o(io,jo)/100.

        call areaini_point (io        , jo          , nlon_i  , nlat_i  , numlon_i  , &
                           lon_i      , lon_i_offset, lat_i   , area_i  , landmask_i, &
                           lsmlon     , lsmlat      , numlon  , lonw    , lats      , &
                           area(io,jo), mask_o      , novr_i2o, iovr_i2o, jovr_i2o  , &
                           wovr_i2o   ) 


! Area-average percent cover on input grid [pctpft_i] to output grid [pctpft_o] 
! and correct [pctpft_o] according to land landmask
! Note that percent cover is in terms of total grid area. 

        wst_sum = 0.
        do m = 0, numpft
           pctpft_o(io,jo,m) = 0.
           do n = 1, novr_i2o !overlap cell index
              ii = iovr_i2o(n) !lon index (input grid) of overlap cell
              ji = jovr_i2o(n) !lat index (input grid) of overlap cell
              pctpft_o(io,jo,m) = pctpft_o(io,jo,m) + pctpft_i(ii,ji,m)*wovr_i2o(n)
           end do
           if (pctlnd_o(io,jo) > 0.) then
              if (landmask(io,jo) == 0) pctpft_o(io,jo,m) = 0.
           else
              pctpft_o(io,jo,m) = 0.
              if (landmask(io,jo) == 1) pctpft_o(io,jo,0) = 100.
           end if
           wst_sum = wst_sum + pctpft_o(io,jo,m) 
        end do

! Error check: percents should sum to 100 for land grid cells

        if (landmask(io,jo) == 1) then
           if (abs(wst_sum-100.) > 0.000001) then
              write (6,*) 'MKPFT error: pft = ', &
                   (pctpft_o(io,jo,m), m = 0, numpft), &
                   ' do not sum to 100. at column, row = ',io,jo, &
                   ' but to ', wst_sum
              call endrun
           end if
        end if

! Find the output pft and pct arrays
! Save percent cover by PFT [wst] and total percent cover [wst_sum]

        wst_sum = 0. 
        sumpct = 0
        do m = 0, numpft
           wst(m) = pctpft_o(io,jo,m)       
           wst_sum = wst_sum + pctpft_o(io,jo,m) 
        end do

! Rank [wst] in ascending order to obtain the top [maxpatch_pft] PFTs 

        if (landmask(io,jo) == 1) then
           call mkrank (numpft, wst, miss, wsti, maxpatch_pft)
        end if

! Fill in [pft] and [pctpft] with data for top [maxpatch_pft] PFTs. 
! If land model grid cell is ocean, set to no PFTs. 
! If land model grid cell is land, there are three possibilities: 
!  1. If [pctlnd_o] = 0, there is no PFT data from the input grid. 
!     Since need land data, use bare ground.
!  2. If [pctlnd_o] > 0, there is some PFT data from the input grid but:
!     a. use the chosen PFT so long as it is not a missing value
!     b. missing value means no more PFTs with cover > 0

        if (landmask(io,jo) == 1) then       !LSM grid wants land
           do m = 1, maxpatch_pft
              if (wsti(m) .ne. miss) then
                 pft(io,jo,m) = wsti(m)          
                 pctpft(io,jo,m) = wst(wsti(m))
              else
                 pft(io,jo,m) = noveg
                 pctpft(io,jo,m) = 0.
              end if
              sumpct = sumpct + pctpft(io,jo,m)
           end do
        else                                     !LSM grid wants ocean
           do m = 1, maxpatch_pft
              pft(io,jo,m) = 0
              pctpft(io,jo,m) = 0.
           end do
        end if

! Correct for the case of more than [maxpatch_pft] PFTs present 

        if (sumpct < wst_sum) then
           diff  = wst_sum - sumpct 
           sumpct = 0.
           do m = 1, maxpatch_pft
              pctpft(io,jo,m) = pctpft(io,jo,m) + diff/maxpatch_pft
              sumpct = sumpct + pctpft(io,jo,m)
           end do
        end if

! Error check: make sure have a valid PFT

        do m = 1, maxpatch_pft
           if (pft(io,jo,m)<0 .or. pft(io,jo,m)>numpft) then
              write (6,*) 'MKPFT error: invalid PFT for i,j=',io,jo,pft(io,jo,m)
              call endrun
           end if
        end do

! Error check: make sure PFTs sum to 100% cover

        if (landmask(io,jo) == 1) then    
           if (abs(sumpct - 100.) > 0.000001) then
              write(6,*) 'MKPFT error: sum(pct) over maxpatch_pft'
              write(6,*) '             is not = 100.'
              write(6,*) sumpct, io,jo
              call endrun
           end if
           if (sumpct < -0.000001) then
              write(6,*) 'MKPFT error: sum(pct) over maxpatch_pft'
              write(6,*) '             is < 0.'
              write(6,*) sumpct, io,jo
              call endrun
           end if
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
     write (6,*) 'MKGLACIER error: input field not conserved'
     write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
     write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
     call endrun
  end if

! -----------------------------------------------------------------
! Error check2
! Compare global areas on input and output grids
! -----------------------------------------------------------------

! input grid

  gpft_i(:) = 0.
  garea_i = 0.
  do ji = 1, nlat_i
     do ii = 1, nlon_i
        garea_i = garea_i + area_i(ii,ji)
        do m = 0, numpft
           gpft_i(m) = gpft_i(m) + pctpft_i(ii,ji,m)*area_i(ii,ji)
        end do
     end do
  end do

! output grid

  gpft_o(:) = 0.
  garea_o = 0.
  do jo = 1, lsmlat
     do io = 1, numlon(jo)
        garea_o = garea_o + area(io,jo)
        do m = 0, numpft
           gpft_o(m) = gpft_o(m) + pctpft_o(io,jo,m)*area(io,jo)
        end do
     end do
  end do

! comparison

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'PFTs Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,1001)
1001 format (1x,'plant type     ',20x,' input grid area',' output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)

  do m = 0, numpft
     write (ndiag,1002) veg(m), gpft_i(m)*1.e-06/100.,gpft_o(m)*1.e-06/100.
1002 format (1x,a35,f16.3,f17.3)
  end do

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,a14)')area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif
     
  write (6,*) 'Successfully made PFTs'
  write (6,*)

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
  deallocate (pctpft_i)

  return
end subroutine mkpft


