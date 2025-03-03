!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine bvoc_vars_to_tile(n, t, pft)
  use bvoc_vars
  use jules53_lsmMod
  use jules_surface_types_mod,  only: npft, nnvg, ntype
  implicit none 
  integer :: n, t, pft

  jules53_struc(n)%jules53(t)%isoprene        = isoprene_gb(1)       ! Gridbox mean isoprene emission flux (kgC/m2/s) 
  jules53_struc(n)%jules53(t)%terpene         = terpene_gb(1)        ! Gridbox mean (mono-)terpene emission flux (kgC/m2/s) 
  jules53_struc(n)%jules53(t)%methanol        = methanol_gb(1)       ! Gridbox mean methanol emission flux (kgC/m2/s) 
  jules53_struc(n)%jules53(t)%acetone         = acetone_gb(1)        ! Gridbox mean acetone emission flux (kgC/m2/s) 
  if(pft .le. npft) then
    jules53_struc(n)%jules53(t)%isoprene_ft(pft)  = isoprene_pft(1, pft) ! Isoprene emission flux on PFTs (kgC/m2/s) 
    jules53_struc(n)%jules53(t)%terpene_ft(pft)   = terpene_pft(1, pft)  ! (Mono-)Terpene emission flux on PFTs (kgC/m2/s) 
    jules53_struc(n)%jules53(t)%methanol_ft(pft)  = methanol_pft(1, pft) ! Methanol emission flux on PFTs (kgC/m2/s) 
    jules53_struc(n)%jules53(t)%acetone_ft(pft)   = acetone_pft(1, pft)  ! Acetone emission flux on PFTs (kgC/m2/s) 
  endif
end subroutine bvoc_vars_to_tile

