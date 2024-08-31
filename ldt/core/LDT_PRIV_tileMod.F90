!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_PRIV_tileMod 
!BOP
!
! !ROUTINE: LDT_PRIV_tileMod
!
! !DESCRIPTION:
!  The code in this file provides a description of the tile data structure in LIS
!
!  \subsubsection{Overview}
!  This module contains the tile data structure used in LIS. The tile space in LIS 
!  is used to simulate subgrid scale variability, where each tile represents a 
!  vegetation class within a grid cell. The tile space is created from the 
!  vegetation distribution within a grid cell. 
! 
!   The tile data structure contains the following variables:
!  \begin{description}
!   \item[col] 
!    column index of the corresponding grid cell
!   \item[row] 
!    row index of the corresponding grid cell
!   \item[index] 
!    index of the corresponding grid cell
!   \item[vegt] 
!    vegetation type of the tile
!   \item[ensem]
!    ensemble index of the tile
!   \item[tile\_id]
!    global catchment id for the tile
!   \item[fgrd]
!    fraction of grid covered by the tile 
!   \item[d2g]
!    local tile count to global tile count
!   \item[com\_lon]
!    center-of-mass longitude of the tile
!   \item[com\_lat]
!    center-of-mass latitude of the tile
!   \end{description}
!
! !REVISION HISTORY:
!  14 Nov 2008: Sujay Kumar; Optimized version of tile representation
!
!EOP  
  implicit none
  public tiledec
  type tiledec
     integer :: col        !Grid Column of Tile
     integer :: row        !Grid Row of Tile
     integer :: index      !Index of corresponding grid
     integer :: ensem      !ensemble id for the tile
     integer :: sftype      !surface type of the tile

     integer :: vegt       !Vegetation Type of Tile

     integer :: soilt      !soil texture Type of Tile
     real    :: sand
     real    :: clay
     real    :: silt

     real    :: elev
     real    :: slope
     real    :: aspect
     real    :: curv

     integer :: tile_id    !global catchment id for the tile
     integer :: d2g        !local tile count to global tile count
     real    :: fgrd       !Fraction of Grid covered by tile 
     real    :: pens       !ensemble weights
     real    :: com_lon    !center-of-mass longitude of the tile
     real    :: com_lat    !center-of-mass latitude of the tile
  end type tiledec
!EOP
end module LDT_PRIV_tileMod
