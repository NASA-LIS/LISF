subroutine tile_to_sf_diag(n, t)
  use sf_diags_mod
  use jules50_lsmMod
  implicit none 
  integer :: n, t

  sf_diag%u10m(1,1)  =  jules50_struc(n)%jules50(t)%u10m
  sf_diag%v10m(1,1)  =  jules50_struc(n)%jules50(t)%v10m
end subroutine tile_to_sf_diag
