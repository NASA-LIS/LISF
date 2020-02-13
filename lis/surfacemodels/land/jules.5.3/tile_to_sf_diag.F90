subroutine tile_to_sf_diag(n, t)
  use sf_diags_mod
  use jules53_lsmMod
  implicit none 
  integer :: n, t

  sf_diag%u10m(1,1)  =  jules53_struc(n)%jules53(t)%u10m
  sf_diag%v10m(1,1)  =  jules53_struc(n)%jules53(t)%v10m
end subroutine tile_to_sf_diag
