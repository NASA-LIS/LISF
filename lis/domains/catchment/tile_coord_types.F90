!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module tile_coord_types  
!BOP
!
! !MODULE: tile_coord_types
! 
! !DESCRIPTION:
! type definitions for tile coordinates and domain 
!
!  !REVISION HISTORY: 
! 26 Jan 2005 Rolf Reichle, Initial Specification
! 14 Apr 2006 Rolf Reichle, split tile_coord.F90 into 2 files to avoid 
!                           having more than one module per file
! 10 Jul 2006 James Geiger, Implementation in LIS
!
! ========================================================================
!EOP
!
! type definition for tile coordinates
  
  implicit none

  PRIVATE

  ! constants that define the kinds (types) of tiles in *.til file
  ! (see field "typ" in tile_coord_type)
  
  integer, public, parameter :: tile_typ_ocean        =   0
  integer, public, parameter :: tile_typ_inlandwater  =  19 
  integer, public, parameter :: tile_typ_ice          =  20
  integer, public, parameter :: tile_typ_land         = 100 
  
  ! ------------------------------------------------------------

  type, public :: tile_coord_type
     
     integer :: tile_id    ! unique tile ID
     integer :: typ        ! (0=ocean, 100=land, 19=inland water, 20=ice)
     integer :: pfaf       ! Pfafstetter number (for land tiles, NOT unique)
     real    :: com_lon    ! center-of-mass longitude
     real    :: com_lat    ! center-of-mass latitude
     real    :: min_lon    ! minimum longitude (bounding box for tile)
     real    :: max_lon    ! maximum longitude (bounding box for tile)
     real    :: min_lat    ! minimum latitude (bounding box for tile)
     real    :: max_lat    ! maximum latitude (bounding box for tile)
     integer :: i_atm      ! AGCM i index
     integer :: j_atm      ! AGCM j index
     real    :: frac_atm   ! area fraction of AGCM cell covered by tile
     real    :: frac_pfaf  ! fraction of Pfafstetter catchment for land tiles 
     real    :: area       ! area [km^2]
     
  end type tile_coord_type
  
  ! ------------------------------------------------------------
  ! 
  ! llx/lly denote the coordinates of the lower left hand
  !  corner of the lower left grid cell (=southwestern corner of domain)
  !
  ! NOTE: convention for tile_coord is -180<=lon<=180, -90<=lat<=90
  
  type, public :: grid_def_type              
     
     integer :: N_lon   ! number of longitude nodes
     integer :: N_lat   ! number of latitude nodes
     real    :: ll_lon  ! lower left latitude of grid cell edge [deg]
     real    :: ll_lat  ! lower left longitude of grid cell edge [deg]
     real    :: dlon    ! longitude grid spacing [deg]
     real    :: dlat    ! latitude grid spacing [deg]
     
  end type grid_def_type
  
  ! *******************************************************************
  
end module tile_coord_types

! =====  EOF ==============================================================

