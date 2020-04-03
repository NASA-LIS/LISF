subroutine tile_to_fluxes(n, t, pft)
  use jules53_lsmMod
  use fluxes
  use jules_surface_types_mod,  only: npft, nnvg, ntype
  implicit none 
  integer, intent(in) :: n, t, pft 

  tstar_ij(1,1)             = jules53_struc(n)%jules53(t)%tstar                 
end subroutine tile_to_fluxes
