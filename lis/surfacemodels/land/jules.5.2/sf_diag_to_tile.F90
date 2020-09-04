subroutine sf_diag_to_tile(n, t)
  use sf_diags_mod
  use jules52_lsmMod
  implicit none 
  integer :: n, t

  jules52_struc(n)%jules52(t)%u10m   = sf_diag%u10m(1,1)
  jules52_struc(n)%jules52(t)%v10m   = sf_diag%v10m(1,1)  
end subroutine sf_diag_to_tile
