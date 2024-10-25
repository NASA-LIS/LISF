!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine sf_diag_to_tile(n, t)
  use sf_diags_mod
  use jules5x_lsmMod
  implicit none 
  integer :: n, t

  jules5x_struc(n)%jules5x(t)%u10m   = sf_diag%u10m(1,1)
  jules5x_struc(n)%jules5x(t)%v10m   = sf_diag%v10m(1,1)  
end subroutine sf_diag_to_tile
