subroutine jules_internal_to_tile(n, t, pft)
  use jules_internal
  use jules5x_lsmMod
  use jules_snow_mod, only: cansnowtile
  use jules_surface_types_mod,  only: npft, nnvg, ntype
  implicit none 
  integer :: n, t, pft
  if ((any(cansnowtile(1:npft)) .eqv. .true.) .and. (pft .le. npft)) then
    jules5x_struc(n)%jules5x(t)%unload_backgrnd(pft)   = unload_backgrnd_pft(1, pft) 
  else
    jules5x_struc(n)%jules5x(t)%unload_backgrnd(:)   = unload_backgrnd_pft(1,:) 
  endif
end subroutine jules_internal_to_tile
