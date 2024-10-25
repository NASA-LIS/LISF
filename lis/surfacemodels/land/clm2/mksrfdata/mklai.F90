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

subroutine mklai (flai, ndiag, pft, mlai, msai, mhgtt, mhgtb)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! make LAI/SAI/height data for vegetated patches (1 to maxpatch_pft)
! 
! Method: 
! Portions of this code could be moved out of the month loop
! for improved efficiency
! 
! Author: Sam Levis
! 
!-----------------------------------------------------------------------
! $Id: mklai.F90,v 1.5 2004/05/07 22:18:36 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2_varpar    !lsm parameters
!  use clm2_varsur    !lsm surface variables  
  use fileutils, only : getfil
  use clm2_areaMod      !area averaging routines 
  use clm2_shr_sys_mod, only : shr_sys_flush 
  implicit none

! ------------------------ arguments ------------------------------
  character(len=*), intent(in) :: flai                         !input lai-sai-hgt dataset
  integer , intent(in) :: ndiag                                !unit number for diagnostic output
  integer , intent(in) :: pft(lsmlon,lsmlat,maxpatch_pft)      !PFT (0 to numpft)
  real(r8), intent(out):: mlai(lsmlon,lsmlat,maxpatch_pft,1)  !monthly lai
  real(r8), intent(out):: msai(lsmlon,lsmlat,maxpatch_pft,1)  !monthly sai
  real(r8), intent(out):: mhgtt(lsmlon,lsmlat,maxpatch_pft,1) !monthly height (top)
  real(r8), intent(out):: mhgtb(lsmlon,lsmlat,maxpatch_pft,1) !monthly height (bottom)
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  character(len=256) :: locfn       !local dataset file name

  integer :: nlon_i                 !input grid : longitude points (read in)
  integer :: nlat_i                 !input grid : latitude  points (read in)
  integer :: ncid,dimid,varid       !input netCDF id's
  integer :: beg4d(4),len4d(4)      !netCDF variable edges
  integer :: ntim                   !number of input time samples
  integer :: ier                    !error status   

  real(r8) :: glai_o(0:numpft)      !output grid: global area pfts
  real(r8) :: gsai_o(0:numpft)      !output grid: global area pfts
  real(r8) :: ghgtt_o(0:numpft)     !output grid: global area pfts
  real(r8) :: ghgtb_o(0:numpft)     !output grid: global area pfts
  real(r8) :: garea_o               !output grid: global area
  real(r8) :: glai_i(0:numpft)      !input grid: global area pfts
  real(r8) :: gsai_i(0:numpft)      !input grid: global area pfts
  real(r8) :: ghgtt_i(0:numpft)     !input grid: global area pfts
  real(r8) :: ghgtb_i(0:numpft)     !input grid: global area pfts
  real(r8) :: garea_i               !input grid: global area

  integer  :: ii                    !longitude index for input grid
  integer  :: io                    !longitude index for LSM grid
  integer  :: ji                    !latitude  index for input grid
  integer  :: jo                    !latitude  index for LSM grid
  integer  :: k,l,m,n               !indices

  integer numpft_i                            !number of plant types on input dataset
  real(r8) :: edge_i(4)                       !input grid: N,E,S,W edges (degrees)
  real(r8), allocatable :: mlai_i(:,:,:)      !monthly lai in
  real(r8), allocatable :: msai_i(:,:,:)      !monthly sai in
  real(r8), allocatable :: mhgtt_i(:,:,:)     !monthly height (top) in
  real(r8), allocatable :: mhgtb_i(:,:,:)     !monthly height (bottom) in
  real(r8), allocatable :: landmask_i(:,:)    !input grid: fraction land
  real(r8), allocatable :: latixy_i(:,:)      !input grid: latitude (degrees)
  real(r8), allocatable :: longxy_i(:,:)      !input grid: longitude (degrees)
  integer , allocatable :: numlon_i(:)        !input grid: number longitude points by lat
  real(r8), allocatable :: lon_i(:,:)         !input grid: longitude, west edge (degrees)
  real(r8), allocatable :: lon_i_offset(:,:)  !input grid: offset longitude, west edge (degrees)
  real(r8), allocatable :: lat_i(:)           !input grid: latitude, south edge (degrees)
  real(r8), allocatable :: area_i(:,:)        !input grid: cell area
  real(r8), allocatable :: mask_i(:,:)        !input grid: mask (0, 1)

  real(r8) :: mlai_o(lsmlon,lsmlat,0:numpft)  !monthly lai out
  real(r8) :: msai_o(lsmlon,lsmlat,0:numpft)  !monthly sai out
  real(r8) :: mhgtt_o(lsmlon,lsmlat,0:numpft) !monthly height (top) out
  real(r8) :: mhgtb_o(lsmlon,lsmlat,0:numpft) !monthly height (bottom) out

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
! -----------------------------------------------------------------

  write (6,*) 'Attempting to make LAIs/SAIs/heights .....'
#if 0 
  call clm2_shr_sys_flush(6)

! -----------------------------------------------------------------
! Determine input grid info
! -----------------------------------------------------------------

  call getfil (flai, locfn, 0)
  call wrap_open(locfn, 0, ncid)

  call wrap_inq_dimid  (ncid, 'lon', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlon_i)

  call wrap_inq_dimid  (ncid, 'lat', dimid)
  call wrap_inq_dimlen (ncid, dimid, nlat_i)

  call wrap_inq_dimid  (ncid, 'pft', dimid)
  call wrap_inq_dimlen (ncid, dimid, numpft_i)
  if (numpft_i .ne. numpft+1) then
     write(6,*)'MKLAI: parameter numpft+1= ',numpft+1, &
          'does not equal input dataset numpft= ',numpft_i
     call endrun
  endif

  call wrap_inq_dimid  (ncid, 'time', dimid)
  call wrap_inq_dimlen (ncid, dimid, ntim)
  if (ntim .ne. 12) then
     write(6,*)'MKLAI: must have 12 time samples on input data'
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
  allocate (area_i(nlon_i,nlat_i), stat=ier)          
  if (ier/=0) call endrun
  allocate (mask_i(nlon_i,nlat_i), stat=ier)          
  if (ier/=0) call endrun
  allocate (mlai_i(nlon_i,nlat_i,0:numpft), stat=ier) 
  if (ier/=0) call endrun
  allocate (msai_i(nlon_i,nlat_i,0:numpft), stat=ier) 
  if (ier/=0) call endrun
  allocate (mhgtt_i(nlon_i,nlat_i,0:numpft), stat=ier)
  if (ier/=0) call endrun
  allocate (mhgtb_i(nlon_i,nlat_i,0:numpft), stat=ier)
  if (ier/=0) call endrun
  allocate (landmask_i(nlon_i,nlat_i), stat=ier)      
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

! Obtain time independent input data

  call wrap_inq_varid (ncid, 'LANDMASK', varid)
  call wrap_get_var_realx (ncid, varid, landmask_i)

! -----------------------------------------------------------------
! Determine input grid cell and cell areas
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

! -----------------------------------------------------------------
! Loop over input months
! -----------------------------------------------------------------

  do m = 1, ntim

! Get input data for the month

     beg4d(1) = 1 ; len4d(1) = nlon_i
     beg4d(2) = 1 ; len4d(2) = nlat_i
     beg4d(3) = 1 ; len4d(3) = numpft_i
     beg4d(4) = m ; len4d(4) = 1

     call wrap_inq_varid (ncid, 'MONTHLY_LAI', varid)
     call wrap_get_vara_realx (ncid, varid , beg4d, len4d, mlai_i )

     call wrap_inq_varid (ncid, 'MONTHLY_SAI', varid)
     call wrap_get_vara_realx (ncid, varid , beg4d, len4d, msai_i )

     call wrap_inq_varid (ncid, 'MONTHLY_HEIGHT_TOP', varid)
     call wrap_get_vara_realx (ncid, varid, beg4d, len4d, mhgtt_i)

     call wrap_inq_varid (ncid, 'MONTHLY_HEIGHT_BOT', varid)
     call wrap_get_vara_realx (ncid, varid, beg4d, len4d, mhgtb_i)

! Process each cell on land model grid
! novr_i2o - number of input grid cells that overlap each land grid cell
! iovr_i2o - longitude index of overlapping input grid cell
! jovr_i2o - latitude  index of overlapping input grid cell
! wovr_i2o - fraction of land grid cell overlapped by input grid cell

     mlai_o(:,:,:)  = 0.
     msai_o(:,:,:)  = 0.
     mhgtt_o(:,:,:) = 0.
     mhgtb_o(:,:,:) = 0.

     mlai(:,:,:,m)  = 0.
     msai(:,:,:,m)  = 0.
     mhgtt(:,:,:,m) = 0.
     mhgtb(:,:,:,m) = 0.

!$OMP PARALLEL DO PRIVATE (io,jo,ii,ji,l,n,mask_o,novr_i2o,iovr_i2o,jovr_i2o,wovr_i2o,fld_i)
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
              
           call areaini_point (io        , jo          , nlon_i  , nlat_i  , numlon_i   , &
                              lon_i      , lon_i_offset, lat_i   , area_i  , landmask_i , &
                              lsmlon     , lsmlat      , numlon  , lonw    , lats       , &
                              area(io,jo), mask_o      , novr_i2o, iovr_i2o, jovr_i2o   , &
                              wovr_i2o   ) 

! Make area average and set oceans to zero

           if (landmask(io,jo) == 0) then
              do l = 0, numpft
                 mlai_o(io,jo,l)  = 0.
                 msai_o(io,jo,l)  = 0.
                 mhgtt_o(io,jo,l) = 0.
                 mhgtb_o(io,jo,l) = 0.
              end do
           else
              do l = 0, numpft
                 do n = 1, novr_i2o !overlap cell index
                    ii = iovr_i2o(n) !lon index (input grid) of overlap cell
                    ji = jovr_i2o(n) !lat index (input grid) of overlap cell
                    mlai_o(io,jo,l)  = mlai_o(io,jo,l)  + mlai_i(ii,ji,l)  * wovr_i2o(n)
                    msai_o(io,jo,l)  = msai_o(io,jo,l)  + msai_i(ii,ji,l)  * wovr_i2o(n)
                    mhgtt_o(io,jo,l) = mhgtt_o(io,jo,l) + mhgtt_i(ii,ji,l) * wovr_i2o(n)
                    mhgtb_o(io,jo,l) = mhgtb_o(io,jo,l) + mhgtb_i(ii,ji,l) * wovr_i2o(n)
                 end do
              end do
           endif  

! Assign lai/sai/hgtt/hgtb to the top [maxpatch_pft] PFTS
! as determined by mkpft.F

           do l = 0, numpft
              do n = 1, maxpatch_pft
                 if (l == pft(io,jo,n)) then
                    mlai(io,jo,n,m) = mlai_o(io,jo,l)
                    msai(io,jo,n,m) = msai_o(io,jo,l)
                    mhgtt(io,jo,n,m) = mhgtt_o(io,jo,l)
                    mhgtb(io,jo,n,m) = mhgtb_o(io,jo,l)
                 end if
              end do
           end do

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
! Compare global sum fld_o to global sum fld_i for month m
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
        write (6,*) 'MKLAI error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        call endrun
     end if
     
! -----------------------------------------------------------------
! Error check2
! Compare global areas on input and output grids
! -----------------------------------------------------------------

! Input grid global area

     glai_i(:)  = 0.
     gsai_i(:)  = 0.
     ghgtt_i(:) = 0.
     ghgtb_i(:) = 0.
     garea_i    = 0.

     do ji = 1, nlat_i
        do ii = 1, nlon_i
           garea_i = garea_i + area_i(ii,ji)
           do l = 0, numpft
              glai_i(l) = glai_i(l)+mlai_i(ii,ji,l)*area_i(ii,ji)
              gsai_i(l) = gsai_i(l)+msai_i(ii,ji,l)*area_i(ii,ji)
              ghgtt_i(l) = ghgtt_i(l)+mhgtt_i(ii,ji,l)*area_i(ii,ji)
              ghgtb_i(l) = ghgtb_i(l)+mhgtb_i(ii,ji,l)*area_i(ii,ji)
           end do
        end do
     end do

! Output grid global area

     glai_o(:)  = 0.
     gsai_o(:)  = 0.
     ghgtt_o(:) = 0.
     ghgtb_o(:) = 0.
     garea_o    = 0.

     do jo = 1, lsmlat
        do io = 1, numlon(jo)
           garea_o = garea_o + area(io,jo)
           do l = 0, numpft
              glai_o(l) = glai_o(l)+mlai_o(io,jo,l)*area(io,jo)
              gsai_o(l) = gsai_o(l)+msai_o(io,jo,l)*area(io,jo)
              ghgtt_o(l) = ghgtt_o(l)+mhgtt_o(io,jo,l)*area(io,jo)
              ghgtb_o(l) = ghgtb_o(l)+mhgtb_o(io,jo,l)*area(io,jo)
           end do
        end do
     end do

! Comparison

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'LAI Output for month ',m
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,1001)
1001 format (1x,'PFT input grid area output grid area',/ &
             1x,3x,'     10**6 km**2','      10**6 km**2')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)

     do l = 0, numpft
        write (ndiag,1002) l, glai_i(l)*1.e-06*1.e-02,glai_o(l)*1.e-06*1.e-02
1002    format (1x,i3,f16.3,f17.3)
     end do

     write (6,*) 'Successfully made LAIs/SAIs/heights for month ', m
     write (6,*)
     call clm2_shr_sys_flush(6)

  end do  ! end loop over time

! Close input file

  call wrap_close(ncid)

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
  deallocate (mlai_i)
  deallocate (msai_i)
  deallocate (mhgtt_i)
  deallocate (mhgtb_i)
#endif
  return
end subroutine mklai
