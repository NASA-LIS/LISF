!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine tile_to_ancil(n, t, pft)
  use ancil_info
  use jules5x_lsmMod
  use jules_surface_mod,      only: l_aggregate
  implicit none 
  integer :: n, t, pft, j

  ssi_pts               =  jules5x_struc(n)%jules5x(t)%ssi_pts                  !  Number of sea or sea-ice points
  sea_pts               =  jules5x_struc(n)%jules5x(t)%sea_pts                  !  Number of sea points
  sice_pts              =  jules5x_struc(n)%jules5x(t)%sice_pts                 !  Number of sea-ice points
  ssi_index(1)          =  jules5x_struc(n)%jules5x(t)%ssi_index                   ! index of sea and sea-ice points
  sea_index(1)          =  jules5x_struc(n)%jules5x(t)%sea_index                   ! index of sea points
  sice_index(1)         =  jules5x_struc(n)%jules5x(t)%sice_index                  ! index of sea-ice points
  sice_pts_ncat(:)      =  jules5x_struc(n)%jules5x(t)%sice_pts_ncat(:)         ! number of points for each sea-ice category
  sice_index_ncat(1, :) =  jules5x_struc(n)%jules5x(t)%sice_index_ncat(:)          ! index of points for each sea-ice category
  
  soilt_pts(1)          =  jules5x_struc(n)%jules5x(t)%soilt_pts 
  soilt_index(1, 1)     =  jules5x_struc(n)%jules5x(t)%soilt_index 
  
  l_soil_point(1)       =  jules5x_struc(n)%jules5x(t)%l_soil_point             !  TRUE if a soil point  FALSE otherwise
  l_lice_point(1)       =  jules5x_struc(n)%jules5x(t)%l_lice_point             !  TRUE if a land ice point  FALSE otherwise
  l_lice_surft(:)       =  jules5x_struc(n)%jules5x(t)%l_lice_surft(:)

  fssi_ij(1,1)             =  jules5x_struc(n)%jules5x(t)%fssi                          !  Fraction of gridbox covered by sea or sea-ice
  sea_frac(1)           =  jules5x_struc(n)%jules5x(t)%sea_frac                    !  Fraction of gridbox covered by sea (converted to single vector array)
  sice_frac(1)          =  jules5x_struc(n)%jules5x(t)%sice_frac                   !  Fraction of gridbox covered by sea-ice converted to single vector array)
  sice_frac_ncat(1,:)   =  jules5x_struc(n)%jules5x(t)%sice_frac_ncat(:)          !  Fraction of gridbox covered by each sea-ice category (converted to single vector array)
  frac_soilt(1,1)       =  jules5x_struc(n)%jules5x(t)%frac_soilt

  halo_i                =  jules5x_struc(n)%jules5x(t)%halo_i                   !  Size of halo in i direction
  halo_j                =  jules5x_struc(n)%jules5x(t)%halo_j                   !  Size of halo in j direction
  n_rows                =  jules5x_struc(n)%jules5x(t)%n_rows                   !  Number of rows in a v field
  off_x                 =  jules5x_struc(n)%jules5x(t)%off_x                    !  Size of small halo in i
  off_y                 =  jules5x_struc(n)%jules5x(t)%off_y                    !  Size of small halo in j
  row_length            =  jules5x_struc(n)%jules5x(t)%row_length               !  Number of points on a row
  rows                  =  jules5x_struc(n)%jules5x(t)%rows                     !  Number of rows in a theta field
  co2_dim_len           =  jules5x_struc(n)%jules5x(t)%co2_dim_len              !  Length of a CO2 field row
  co2_dim_row           =  jules5x_struc(n)%jules5x(t)%co2_dim_row              !  Number of CO2 field rows
  
  land_pts              =  jules5x_struc(n)%jules5x(t)%land_pts                 !  No. of land points
  land_pts_trif         =  jules5x_struc(n)%jules5x(t)%land_pts_trif            !  For dimensioning land fields in TRIFFID
  lice_pts              =  jules5x_struc(n)%jules5x(t)%lice_pts                 !  Number of land ice points
  npft_trif             =  jules5x_struc(n)%jules5x(t)%npft_trif                !  For dimensioning pft fields in TRIFFID =npft when TRIFFID on, otherwise =1
  nsurft                =  jules5x_struc(n)%jules5x(t)%nsurft                   !  Number of surface tiles
  soil_pts              =  jules5x_struc(n)%jules5x(t)%soil_pts                 !  Number of soil points
  dim_cs1               =  jules5x_struc(n)%jules5x(t)%dim_cs1                  !  size of second dimension in soil carbon (cs) and related respiration variables
  dim_cs2               =  jules5x_struc(n)%jules5x(t)%dim_cs2                  !  size used for some variables that are only used with TRIFFID. If not using TRIFFID these variables are set to be smaller to save some space
  
  land_index(1)         =  jules5x_struc(n)%jules5x(t)%land_index                  !  index of land points
  soil_index(1)         =  jules5x_struc(n)%jules5x(t)%soil_index                  !  index of soil points (i.e. land point number for each soil point)
  lice_index(1)         =  jules5x_struc(n)%jules5x(t)%lice_index               !  index of land ice points (i.e. land point number for each land ice point)
  
  ! for supporting l_aggregate 
  !surft_pts(pft)        =  jules5x_struc(n)%jules5x(t)%tile_pts(pft)              !  Number of land points which include the nth surface type
  !surft_index(1, pft)   =  jules5x_struc(n)%jules5x(t)%tile_index(pft)               !  indices of land points which include the nth surface type
  !frac_surft(1,pft)     =  jules5x_struc(n)%jules5x(t)%frac(pft)                    !  fractional cover of each surface type
  
  if(l_aggregate) then
     do j=1, jules5x_struc(n)%ntype
        if(jules5x_struc(n)%jules5x(t)%surft_frac(j)>0.0) then
          surft_pts(j)     = 1
          surft_index(1,j) = 1
          frac_surft(1,j)  = jules5x_struc(n)%jules5x(t)%surft_frac(j)
        endif
     enddo
  else
     pft = jules5x_struc(n)%jules5x(t)%pft
     surft_pts(pft)     = jules5x_struc(n)%jules5x(t)%tile_pts(pft)
     surft_index(1,pft) = jules5x_struc(n)%jules5x(t)%tile_index(pft)
     frac_surft(1, pft) = jules5x_struc(n)%jules5x(t)%frac(pft)
  endif 


  z1_tq_ij(1,1)            =  jules5x_struc(n)%jules5x(t)%z1_tq                         !  height of temperature data
  z1_uv_ij(1,1)            =  jules5x_struc(n)%jules5x(t)%z1_uv                         !  height of wind data
  ice_fract_ij(1,1)        =  jules5x_struc(n)%jules5x(t)%ice_fract                     !  fraction of gridbox covered by sea-ice  (decimal fraction)
  ice_fract_ncat_sicat(1, 1,:)   =  jules5x_struc(n)%jules5x(t)%ice_fract_ncat(:)          !  fraction of gridbox covered by sea-ice on catagories
  ti_cat_sicat(1,1, :)        =  jules5x_struc(n)%jules5x(t)%ti_cat(:)                     ! sea ice surface temperature on categories
  pond_frac_cat_sicat(1,1,:)  =  jules5x_struc(n)%jules5x(t)%pond_frac_cat(:)             ! Meltpond fraction on sea ice categories 
  pond_depth_cat_sicat(1,1,:) =  jules5x_struc(n)%jules5x(t)%pond_depth_cat(:)            ! Meltpond depth on sea ice categories (m)  
  sstfrz_ij(1,1)           =  jules5x_struc(n)%jules5x(t)%sstfrz                        ! Salinity-dependent sea surface freezing temperature (K)
  land_mask(1,1)        =  jules5x_struc(n)%jules5x(t)%land_mask                     !     t if land, f elsewhere
end subroutine tile_to_ancil 


