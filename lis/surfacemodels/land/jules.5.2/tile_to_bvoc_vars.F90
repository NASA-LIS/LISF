!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine tile_to_bvoc_vars(n, t, pft)
  use bvoc_vars
  use jules52_lsmMod
  use jules_surface_types_mod,  only: npft, nnvg, ntype
  implicit none 
  integer :: n, t, pft

  isoprene_gb(1)        = jules52_struc(n)%jules52(t)%isoprene       ! Gridbox mean isoprene emission flux (kgC/m2/s) 
  terpene_gb(1)         = jules52_struc(n)%jules52(t)%terpene        ! Gridbox mean (mono-)terpene emission flux (kgC/m2/s) 
  methanol_gb(1)        = jules52_struc(n)%jules52(t)%methanol       ! Gridbox mean methanol emission flux (kgC/m2/s) 
  acetone_gb(1)         = jules52_struc(n)%jules52(t)%acetone        ! Gridbox mean acetone emission flux (kgC/m2/s) 
  if(pft .le. npft) then
    isoprene_pft(1, pft)  = jules52_struc(n)%jules52(t)%isoprene_ft(pft) ! Isoprene emission flux on PFTs (kgC/m2/s) 
    terpene_pft(1, pft)   = jules52_struc(n)%jules52(t)%terpene_ft(pft)  ! (Mono-)Terpene emission flux on PFTs (kgC/m2/s) 
    methanol_pft(1, pft)  = jules52_struc(n)%jules52(t)%methanol_ft(pft) ! Methanol emission flux on PFTs (kgC/m2/s) 
    acetone_pft(1, pft)   = jules52_struc(n)%jules52(t)%acetone_ft(pft)  ! Acetone emission flux on PFTs (kgC/m2/s) 
  endif
end subroutine tile_to_bvoc_vars

