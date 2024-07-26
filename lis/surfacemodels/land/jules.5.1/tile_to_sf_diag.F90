!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine tile_to_sf_diag(n, t)
  use sf_diags_mod
  use jules51_lsmMod
  implicit none 
  integer :: n, t

  sf_diag%u10m(1,1)  =  jules51_struc(n)%jules51(t)%u10m
  sf_diag%v10m(1,1)  =  jules51_struc(n)%jules51(t)%v10m
end subroutine tile_to_sf_diag
