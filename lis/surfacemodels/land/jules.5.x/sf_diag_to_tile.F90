subroutine sf_diag_to_tile(n, t)
  use sf_diags_mod
  use jules5x_lsmMod
  implicit none 
  integer :: n, t

  jules5x_struc(n)%jules5x(t)%u10m   = sf_diag%u10m(1,1)
  jules5x_struc(n)%jules5x(t)%v10m   = sf_diag%v10m(1,1)  
end subroutine sf_diag_to_tile
