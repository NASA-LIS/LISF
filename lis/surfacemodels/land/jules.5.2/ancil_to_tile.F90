!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine ancil_to_tile(n, t, pft)
  use ancil_info
  use jules52_lsmMod
  implicit none 
  integer :: n, t, pft

  jules52_struc(n)%jules52(t)%ssi_pts                =  ssi_pts                !  Number of sea or sea-ice points
  jules52_struc(n)%jules52(t)%sea_pts                =  sea_pts                !  Number of sea points
  jules52_struc(n)%jules52(t)%sice_pts               =  sice_pts               !  Number of sea-ice points
  jules52_struc(n)%jules52(t)%ssi_index              =  ssi_index(1)              ! index of sea and sea-ice points
  jules52_struc(n)%jules52(t)%sea_index              =  sea_index(1)              ! index of sea points
  jules52_struc(n)%jules52(t)%sice_index             =  sice_index(1)             ! index of sea-ice points
  jules52_struc(n)%jules52(t)%sice_pts_ncat(:)       =  sice_pts_ncat(:)       ! number of points for each sea-ice category
  jules52_struc(n)%jules52(t)%sice_index_ncat(:)     =  sice_index_ncat(1, :)     ! index of points for each sea-ice category

  jules52_struc(n)%jules52(t)%soilt_pts              =  soilt_pts(1) 
  jules52_struc(n)%jules52(t)%soilt_index            =  soilt_index(1, 1)               !  indices of land points which include the nth surface type
  
  jules52_struc(n)%jules52(t)%l_soil_point           =  l_soil_point(1)        !  TRUE if a soil point  FALSE otherwise
  jules52_struc(n)%jules52(t)%l_lice_point           =  l_lice_point(1)        !  TRUE if a land ice point  FALSE otherwise
  jules52_struc(n)%jules52(t)%l_lice_surft(:)        =  l_lice_surft(:)        ! TRUE if a land ice (surface) tile, FALSE otherwise
  
  jules52_struc(n)%jules52(t)%fssi                   =  fssi_ij(1,1)                   !  Fraction of gridbox covered by sea or sea-ice
  jules52_struc(n)%jules52(t)%sea_frac               =  sea_frac(1)               !  Fraction of gridbox covered by sea (converted to single vector array)
  jules52_struc(n)%jules52(t)%sice_frac              =  sice_frac(1)              !  Fraction of gridbox covered by sea-ice converted to single vector array)
  jules52_struc(n)%jules52(t)%sice_frac_ncat(:)      =  sice_frac_ncat(1,:)      !  Fraction of gridbox covered by each sea-ice category (converted to single vector array)
  jules52_struc(n)%jules52(t)%frac_soilt             =  frac_soilt(1,1)        ! Fraction of gridbox for each soil tile
  
  jules52_struc(n)%jules52(t)%halo_i                 =  halo_i                 !  Size of halo in i direction
  jules52_struc(n)%jules52(t)%halo_j                 =  halo_j                 !  Size of halo in j direction
  jules52_struc(n)%jules52(t)%n_rows                 =  n_rows                 !  Number of rows in a v field
  jules52_struc(n)%jules52(t)%off_x                  =  off_x                  !  Size of small halo in i
  jules52_struc(n)%jules52(t)%off_y                  =  off_y                  !  Size of small halo in j
  jules52_struc(n)%jules52(t)%row_length             =  row_length             !  Number of points on a row
  jules52_struc(n)%jules52(t)%rows                   =  rows                   !  Number of rows in a theta field
  jules52_struc(n)%jules52(t)%co2_dim_len            =  co2_dim_len            !  Length of a CO2 field row
  jules52_struc(n)%jules52(t)%co2_dim_row            =  co2_dim_row            !  Number of CO2 field rows
  
  jules52_struc(n)%jules52(t)%land_pts               =  land_pts               !  No. of land points
  jules52_struc(n)%jules52(t)%land_pts_trif          =  land_pts_trif          !  For dimensioning land fields in TRIFFID
  jules52_struc(n)%jules52(t)%lice_pts               =  lice_pts               !  Number of land ice points
  jules52_struc(n)%jules52(t)%npft_trif              =  npft_trif              !  For dimensioning pft fields in TRIFFID =npft when TRIFFID on, otherwise =1
  jules52_struc(n)%jules52(t)%nsurft                 =  nsurft                 !  Number of surface tiles
  jules52_struc(n)%jules52(t)%soil_pts               =  soil_pts               !  Number of soil points
  jules52_struc(n)%jules52(t)%dim_cs1                =  dim_cs1                !  size of second dimension in soil carbon (cs) and related respiration variables
  jules52_struc(n)%jules52(t)%dim_cs2                =  dim_cs2                !  size used for some variables that are only used with TRIFFID. If not using TRIFFID these variables are set to be smaller to save some space
  
  jules52_struc(n)%jules52(t)%land_index             =  land_index(1)             !  index of land points
  jules52_struc(n)%jules52(t)%tile_index(pft)        =  surft_index(1, pft)          !  indices of land points which include the nth surface type
  jules52_struc(n)%jules52(t)%soil_index             =  soil_index(1)             !  index of soil points (i.e. land point number for each soil point)
  jules52_struc(n)%jules52(t)%lice_index             =  lice_index(1)          !  index of land ice points (i.e. land point number for each land ice point)
  jules52_struc(n)%jules52(t)%tile_pts(pft)          =  surft_pts(pft)            !  Number of land points which include the nth surface type
  jules52_struc(n)%jules52(t)%frac(pft)              =  frac_surft(1,pft)                !  fractional cover of each surface type
  
  jules52_struc(n)%jules52(t)%z1_tq                  =  z1_tq_ij(1,1)                  !  height of temperature data
  jules52_struc(n)%jules52(t)%z1_uv                  =  z1_uv_ij(1,1)                  !  height of wind data
  jules52_struc(n)%jules52(t)%ice_fract              =  ice_fract_ij(1,1)              !  fraction of gridbox covered by sea-ice  (decimal fraction)
  jules52_struc(n)%jules52(t)%ice_fract_ncat(:)      =  ice_fract_ncat_sicat(1,1,:)      !  fraction of gridbox covered by sea-ice on catagories
  jules52_struc(n)%jules52(t)%ti_cat(:)              =  ti_cat_sicat(1,1, :)              ! sea ice surface temperature on categories
  jules52_struc(n)%jules52(t)%pond_frac_cat(:)       =  pond_frac_cat_sicat(1,1,:)       ! Meltpond fraction on sea ice categories 
  jules52_struc(n)%jules52(t)%pond_depth_cat(:)      =  pond_depth_cat_sicat(1,1,:)      ! Meltpond depth on sea ice categories (m)  
  jules52_struc(n)%jules52(t)%sstfrz                 =  sstfrz_ij(1,1)                 ! Salinity-dependent sea surface freezing temperature (K)
  jules52_struc(n)%jules52(t)%land_mask              =  land_mask(1,1)              !     t if land, f elsewhere
end subroutine ancil_to_tile

