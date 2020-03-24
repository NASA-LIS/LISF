subroutine sf_diag_to_tile(n, t)
  use sf_diags_mod
  use jules50_lsmMod
  implicit none 
  integer :: n, t

  jules50_struc(n)%jules50(t)%u10m   = sf_diag%u10m(1,1)
  jules50_struc(n)%jules50(t)%v10m   = sf_diag%v10m(1,1)  
end subroutine sf_diag_to_tile
