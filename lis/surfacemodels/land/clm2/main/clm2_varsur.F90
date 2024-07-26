!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

module clm2_varsur

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 2-d surface boundary data 
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: clm2_varsur.F90,v 1.6 2004/11/24 22:57:14 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod, only : r8
!  use clm2_varpar
  use clm2_varpar, only : nlevsoi, nlevlak
  implicit none

  PRIVATE
! land model grid

!  integer  numlon(lsmlat)                !longitude points for each latitude strip
!  real(r8) latixy(lsmlon,lsmlat)         !latitude of grid cell (degrees)
!  real(r8) longxy(lsmlon,lsmlat)         !longitude of grid cell (degrees)
!  real(r8) area(lsmlon,lsmlat)           !grid cell area (km**2)
!  real(r8) lats(lsmlat+1)                !grid cell latitude, southern edge (degrees)
!  real(r8) lonw(lsmlon+1,lsmlat)         !grid cell longitude, western edge (degrees)
!  real(r8) lsmedge(4)                    !North,East,South,West edges of grid (deg)
!  logical :: pole_points                 !true => grid has pole points
!  logical :: fullgrid  = .true.          !true => no grid reduction towards poles
!  logical :: offline_rdgrid              !true => read offline grid rather than creating it

! fractional land and mask

!  integer  landmask(lsmlon,lsmlat)       !land mask: 1 = land. 0 = ocean
!  real(r8) landfrac(lsmlon,lsmlat)       !fractional land

! surface boundary data 

!  integer  soic2d(lsmlon,lsmlat)         !soil color
!  real(r8) sand3d(lsmlon,lsmlat,nlevsoi) !soil texture: percent sand
!  real(r8) clay3d(lsmlon,lsmlat,nlevsoi) !soil texture: percent clay
!  real(r8) pctgla(lsmlon,lsmlat)         !percent of grid cell that is glacier
!  real(r8) pctlak(lsmlon,lsmlat)         !percent of grid cell that is lake
!  real(r8) pctwet(lsmlon,lsmlat)         !percent of grid cell that is wetland
!  real(r8) pcturb(lsmlon,lsmlat)         !percent of grid cell that is urbanized

! lake and soil levels

  real(r8), public :: zlak(1:nlevlak)            !lake z  (layers) 
  real(r8), public :: dzlak(1:nlevlak)           !lake dz (thickness)
  real(r8), public :: zsoi(1:nlevsoi)            !soil z  (layers)
  real(r8), public :: dzsoi(1:nlevsoi)           !soil dz (thickness)
  real(r8), public :: zisoi(0:nlevsoi)           !soil zi (interfaces)  


end module clm2_varsur
