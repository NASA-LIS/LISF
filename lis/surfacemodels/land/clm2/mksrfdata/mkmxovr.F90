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

subroutine mkmxovr (nlon_i, nlat_i, numlon_i, lon_i, lat_i, &
                    nlon_o, nlat_o, numlon_o, lon_o, lat_o, &
                    mxovr , n_ovr  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! find maxinum numver of overlapping cells
! 
! Method: 
! For each output grid cell: find overlapping input grid cells that
! that overlap with output grid cell. Cells overlap if:
!
! southern edge of input grid < northern edge of output grid AND
! northern edge of input grid > southern edge of output grid
!
! western edge of input grid < eastern edge of output grid AND
! eastern edge of input grid > western edge of output grid
!
!           lon_o(io,jo)      lon_o(io+1,jo)
!
!              |                   |
!              --------------------- lat_o(jo+1)
!              |                   |
!              |                   |
!    xxxxxxxxxxxxxxx lat_i(ji+1)   |
!    x         |   x               |
!    x  input  |   x   output      |
!    x  cell   |   x    cell       |
!    x  ii,ji  |   x   io,jo       |
!    x         |   x               |
!    x         ----x---------------- lat_o(jo  )
!    x             x
!    xxxxxxxxxxxxxxx lat_i(ji  )
!    x             x
! lon_i(ii,ji) lon_i(ii+1,ji)
!
!
! The above diagram assumes both grids are oriented South to North. Other
! combinations of North to South and South to North grids are possible:
!
!     Input Grid    Output Grid
!     -------------------------
! (1)   S to N        S to N
! (2)   N to S        N to S
! (3)   S to N        N to S
! (4)   N to S        S to N
!
! The code has been modified to allow for North to South grids. Verification
! that these changes work are: 
!    o (1) and (4) give same results for output grid
!    o (2) and (3) give same results for output grid
!    o (2) and (4) give same results for output grid when output grid inverted
!
! WARNING: this code does not vectorize but is only called during start-up
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
!
! $Id: mkmxovr.F90,v 1.5 2004/05/07 22:18:36 jim Exp $ 
!
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2_areaMod
  implicit none

! ------------------------ arguments ---------------------------------
  integer , intent(in) :: nlon_i                 !input grid : max number of longitude points
  integer , intent(in) :: nlat_i                 !input grid : number of latitude points
  integer , intent(in) :: numlon_i(nlat_i)       !input grid : number of lon points for lat
  real(r8), intent(in) :: lon_i(nlon_i+1,nlat_i) !input grid : cell longitude, W edge (deg)
  real(r8), intent(in) :: lat_i(nlat_i+1)        !input grid : cell latitude, S edge (deg)
  integer , intent(in) :: nlon_o                 !output grid: max number of longitude points
  integer , intent(in) :: nlat_o                 !output grid: number of latitude points
  integer , intent(in) :: numlon_o(nlat_o)       !output grid: number of lon points for lat
  real(r8), intent(in) :: lon_o(nlon_o+1,nlat_o) !output grid: cell longitude, W edge (deg)
  real(r8), intent(in) :: lat_o(nlat_o+1)        !output grid: cell latitude, S edge (deg)
  integer , intent(out):: n_ovr(nlon_o,nlat_o)   !number of overlapping input cells
  integer , intent(out):: mxovr                  !maximum number of overlapping input cells
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  integer, parameter :: mxovr_ceiling = 100000 !very large value should only check for bad error
  integer :: ii          !input  grid longitude loop index
  integer :: ji          !input  grid latitude  loop index
  integer :: io          !output grid longitude loop index
  integer :: jo          !output grid latitude  loop index
  integer :: indexi1     !input  grid lat. index according to orientn
  integer :: indexi2     !input  grid lat. index according to orientn
  integer :: indexi3     !input  grid lat. index according to orientn
  integer :: indexo1     !output grid lat. index according to orientn
  integer :: indexo2     !output grid lat. index according to orientn
  integer :: indexo3     !output grid lat. index according to orientn
  real(r8) :: lonw       !west longitudes of overlap 
  real(r8) :: lone       !east longitudes of overlap 
  real(r8) :: dx         !difference in longitudes
  real(r8) :: lats       !south latitudes of overlap
  real(r8) :: latn       !north latitudes of overlap
  real(r8) :: dy         !difference in latitudes 
  real(r8) :: deg2rad    !pi/180
! --------------------------------------------------------------------

! Set number of overlapping cells to zero and initialize mxovr and deg2rad

  mxovr = 0
  deg2rad = (SHR_CONST_PI) / 180.
  n_ovr(:,:) = 0

! loop through output grid cells

  do jo = 1, nlat_o

! choose the right index according to the orientation of the data

     if (lat_o(nlat_o+1) > lat_o(1)) then
        indexo1 = jo+1        !south to north along the edges 
        indexo2 = jo          !south to north along the edges
        indexo3 = jo          !south to north at the center of cell
     else
        indexo1 = nlat_o+1-jo !north to south along the edges
        indexo2 = nlat_o+2-jo !north to south along the edges
        indexo3 = nlat_o+1-jo !north to south at the center of cell
     end if

     do io = 1, numlon_o(indexo3)

! loop through all input grid cells to find overlap with output grid

        do ji = 1, nlat_i                            

! choose the right index according to the orientation of the data

           if (lat_i(nlat_i+1) > lat_i(1)) then
              indexi1 = ji          !south to north along the edges
              indexi2 = ji+1        !south to north along the edges
              indexi3 = ji          !south to north at the center of cell
           else
              indexi1 = nlat_i+2-ji !north to south along the edges
              indexi2 = nlat_i+1-ji !north to south along the edges
              indexi3 = nlat_i+1-ji !north to south at the center of cell
           end if

! if lat and lon okay then increment number of overlapping cells 
! make sure 0 < n_ovr < mxovr_ceiling

           if (lat_i(indexi1)<lat_o(indexo1) .and. lat_i(indexi2)>lat_o(indexo2)) then

              do ii = 1, numlon_i(indexi3)

                 if (lon_i(ii,indexi3)<lon_o(io+1,indexo3) .and. &
                     lon_i(ii+1,indexi3)>lon_o(io,indexo3)) then

                    n_ovr(io,indexo3) = n_ovr(io,indexo3) + 1
                    if (n_ovr(io,indexo3) > mxovr_ceiling) then
                       write (6,100) n_ovr(io,indexo3),mxovr_ceiling,io,indexo3
                       call endrun
                    end if
                    if (n_ovr(io,indexo3) > mxovr) then
                       mxovr = n_ovr(io,indexo3)
                    endif

                 end if
              end do

           end if
        end do

     end do
  end do

100 format(' ','MKMXOVR error: n_ovr= ',i4,' exceeded mx_ceiling = ', &
         i4,' for output lon,lat = ',i4,i4)

  return
end subroutine mkmxovr
