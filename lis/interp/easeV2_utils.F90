!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

! +-======-+
!  Copyright (c) 2003-2007 United States Government as represented by
!  the Admistrator of the National Aeronautics and Space Administration. 
!  All Rights Reserved.
! 
!  THIS OPEN  SOURCE  AGREEMENT  ("AGREEMENT") DEFINES  THE  RIGHTS  OF USE,
!  REPRODUCTION,  DISTRIBUTION,  MODIFICATION AND REDISTRIBUTION OF CERTAIN
!  COMPUTER SOFTWARE ORIGINALLY RELEASED BY THE UNITED STATES GOVERNMENT AS
!  REPRESENTED BY THE GOVERNMENT AGENCY LISTED BELOW ("GOVERNMENT AGENCY"). 
!  THE UNITED STATES GOVERNMENT, AS REPRESENTED BY GOVERNMENT AGENCY, IS AN
!  INTENDED  THIRD-PARTY  BENEFICIARY  OF  ALL  SUBSEQUENT DISTRIBUTIONS OR
!  REDISTRIBUTIONS  OF THE  SUBJECT  SOFTWARE.  ANYONE WHO USES, REPRODUCES,
!  DISTRIBUTES, MODIFIES  OR REDISTRIBUTES THE SUBJECT SOFTWARE, AS DEFINED
!  HEREIN, OR ANY PART THEREOF,  IS,  BY THAT ACTION, ACCEPTING IN FULL THE
!  RESPONSIBILITIES AND OBLIGATIONS CONTAINED IN THIS AGREEMENT.
! 
!  Government Agency: National Aeronautics and Space Administration
!  Government Agency Original Software Designation: GSC-15354-1
!  Government Agency Original Software Title:  GEOS-5 GCM Modeling Software
!  User Registration Requested.  Please Visit http://opensource.gsfc.nasa.gov
!  Government Agency Point of Contact for Original Software: 
!                       Dale Hithon, SRA Assistant, (301) 286-2691
! 
! +-======-+

module easeV2_utils
 
  ! ==========================================================================
  !
  ! easeV2_conv.F90 - FORTRAN routines for converting grid coordinates
  !                   (latitude/longitude <--> row/column indices)
  !                   of the Equal Area Scalable Earth, version 2 (EASEv2) grid
  !
  !    ***** ONLY cylindrical ('M') projection implemented *****
  !
  ! Ported from Steven Chan's matlab code (smapease2inverse.m,
  ! smapease2forward.m), which has been ported from NSIDC's IDL code
  ! (wgs84_convert.pro, wgs84_inverse.pro) available from 
  ! ftp://sidads.colorado.edu/pub/tools/easegrid/geolocation_tools/
  !
  ! 04-Apr-2013 - reichle
  !
  ! ==========================================================================

  implicit none
 
  private
 
  public :: easeV2_convert
  public :: easeV2_inverse

  ! ***NEVER*** change these constants to GEOS-5 MAPL constants!!!!
 
  ! radius of the earth (m) and map eccentricity
 
  real*8, parameter :: map_equatorial_radius_m         = 6378137.0
 
  real*8, parameter :: map_eccentricity                = 0.081819190843
 
  real*8, parameter :: PI                              = 3.14159265358979323846
 
  real*8, parameter :: e2      = map_eccentricity * map_eccentricity
  real*8, parameter :: e4      = e2 * e2
  real*8, parameter :: e6      = e2 * e4
 
  real*8, parameter :: epsilon = 1.e-6
 
  real*8, parameter :: map_reference_longitude         =   0.0  ! 'M', 'N', 'S'
 
  ! constants for 'N' and 'S' (azimuthal) projections
 
  real*8, parameter :: N_map_reference_latitude        =  90.0 
  real*8, parameter :: S_map_reference_latitude        = -90.0
 
  ! constants for 'M' (cylindrical) projection
 
  real*8, parameter :: M_map_reference_latitude        =   0.0
  real*8, parameter :: M_map_second_reference_latitude =  30.0
 
  real*8, parameter :: M_sin_phi1 = sin(M_map_second_reference_latitude*PI/180.)
  real*8, parameter :: M_cos_phi1 = cos(M_map_second_reference_latitude*PI/180.)
 
  real*8, parameter :: M_kz = M_cos_phi1/sqrt(1.0-e2*M_sin_phi1*M_sin_phi1)
 
 
contains 
 
  ! *******************************************************************
 
  subroutine easeV2_convert (grid, lat, lon, col_ind, row_ind)
   
    ! convert geographic coordinates (spherical earth) to
    ! azimuthal equal area or equal area cylindrical grid coordinates
    !
    ! *** NOTE order of calling arguments:  "lat-lon-lon-lat" ***
    !
    ! useage: call easeV2_convert (grid, lat, lon, r, s)
    !
    ! input : grid - projection name '[M][xx]'
    !            where xx = approximate resolution [km]
    !               ie xx = "01", "03", "09", "36"       (SMAP)
    !         lat, lon = geo. coords. (decimal degrees)
    !
    ! output: col_ind, row_ind - column, row coordinates
    !
    ! --------------------------------------------------------------------------
       
    character*(*), intent(in)  :: grid
    real,          intent(in)  :: lat, lon
    real,          intent(out) :: col_ind, row_ind

    ! local variables
   
    integer :: cols, rows, scale
    real*8  :: dlon, phi, lam, rho, map_scale_m, r0, s0, ms, x, y, sin_phi, q
   
    ! ---------------------------------------------------------------------
   
    call easeV2_get_params( grid, map_scale_m, cols, rows, r0, s0 )

    dlon = lon
   
    if (abs(map_reference_longitude)>epsilon) then
       
       dlon = lon - map_reference_longitude
       
    end if

    if (dlon .lt. -180.0) dlon = dlon + 360.0
    if (dlon .gt.  180.0) dlon = dlon - 360.0
   
    phi =  lat*PI/180.   ! convert from degree to radians
    lam = dlon*PI/180.   ! convert from degree to radians
   
    sin_phi = sin(phi)
   
    ms      = map_eccentricity*sin_phi
   
    q = (1. - e2)*                                                     &
         (                                                             &
         (sin_phi /(1. - e2*sin_phi*sin_phi))                          &
         -                                                             &
         .5/map_eccentricity*log((1.-ms)/(1.+ms))                      &
         )
   
    ! note: "qp" only needed for 'N' and 'S' projections
   
    if      (grid(1:1).eq.'M') then
       
       x =  map_equatorial_radius_m*M_kz*lam
       
       y = (map_equatorial_radius_m*q)/(2.*M_kz)
       
    else
       
       print *,'Polar projections not implemented yet'
       stop
       
    endif
   
    row_ind = s0 - (y/map_scale_m)
    col_ind = r0 + (x/map_scale_m)
   
  end subroutine easeV2_convert
 
  ! *******************************************************************
 
  subroutine easeV2_inverse (grid, r, s, lat, lon)
   
    ! convert azimuthal equal area or equal area cylindrical
    ! grid coordinates to geographic coordinates (spherical earth)
    !
    ! *** NOTE order of calling arguments:  "lon-lat-lat-lon" ***
    !
    ! useage: call easeV1_inverse (grid, r, s, lat, lon)
    !
    ! input : grid - projection name '[M][xx]'
    !            where xx = approximate resolution [km]
    !               ie xx = "01", "03", "09", "36"       (SMAP)
    !         r, s - column, row coordinates
    !
    ! output: lat, lon = geo. coords. (decimal degrees)
    !
    ! --------------------------------------------------------------------------

    character*(*), intent(in)  :: grid
    real,          intent(in)  :: r, s
    real,          intent(out) :: lat, lon

    ! local variables
   
    integer   :: cols, rows
    real*8    :: phi, lam, map_scale_m, r0, s0, beta, x, y, qp
   
    ! ---------------------------------------------------------------------
   
    call easeV2_get_params( grid, map_scale_m, cols, rows, r0, s0 )
   
    x =  (r - r0)*map_scale_m
    y = -(s - s0)*map_scale_m
   
    qp = (1. - e2)*                                                           &
         (                                                                    &
         (1./(1.-e2))                                                         &
         -                                                                    &
         .5/map_eccentricity*log((1.-map_eccentricity)/(1.+map_eccentricity)) &
         )
   
    if      (grid(1:1).eq.'M') then
       
       beta = asin(2.*y*M_kz/(map_equatorial_radius_m*qp))
       
       lam  = x/(map_equatorial_radius_m*M_kz)
       
    else
       
       print *,'Polar projections not implemented yet'
       stop
       
    endif
   
    phi = beta                                                              &
         + ( ( e2/3.       + 31./180.*e4 + 517./ 5040.*e6 )*sin(2.*beta) )  &
         + ( (               23./360.*e4 + 251./ 3780.*e6 )*sin(4.*beta) )  &
         + ( (                             761./45360.*e6 )*sin(6.*beta) )
   
    lat = phi*180./PI                            ! convert from radians to degree
    lon = lam*180./PI + map_reference_longitude  ! convert from radians to degree
   
    if (lon .lt. -180.0) lon = lon + 360.0
    if (lon .gt.  180.0) lon = lon - 360.0
   
  end subroutine easeV2_inverse
 
  ! *******************************************************************
 
  subroutine easeV2_get_params( grid, map_scale_m, cols, rows, r0, s0 )
   
    implicit none
   
    character*(*), intent(in)  :: grid
    real*8,        intent(out) :: map_scale_m, r0, s0
    integer,       intent(out) :: cols, rows
   
   
    if (grid(1:1).eq.'M') then
       
       if      (grid .eq. 'M36') then      ! SMAP 36 km grid
         
          map_scale_m = 36032.220840584    ! nominal cell size in meters
          cols = 964
          rows = 406
          r0 = (cols-1)/2.0
          s0 = (rows-1)/2.0

       else if (grid .eq. 'M09') then      ! SMAP  9 km grid

          map_scale_m = 9008.055210146     ! nominal cell size in meters
          cols = 3856
          rows = 1624
          r0 = (cols-1)/2.0
          s0 = (rows-1)/2.0
         
       else if (grid .eq. 'M03') then      ! SMAP  3 km grid

          map_scale_m = 3002.6850700487    ! nominal cell size in meters
          cols = 11568
          rows = 4872
          r0 = (cols-1)/2.0
          s0 = (rows-1)/2.0
         
       else if (grid .eq. 'M01') then      ! SMAP  1 km grid

          map_scale_m = 1000.89502334956   ! nominal cell size in meters
          cols = 34704
          rows = 14616
          r0 = (cols-1)/2.0
          s0 = (rows-1)/2.0
       
       else
         
          print *,'easeV2_convert: unknown resolution: ',grid
          stop
       
       endif

    else if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then
       
       print *,'Polar projections not implemented yet'
       stop
       
    else
       
       print *, 'easeV2_convert: unknown projection: ', grid
       stop
       
    endif
           
  end subroutine easeV2_get_params
 
  ! *******************************************************************
 
end module easeV2_utils

