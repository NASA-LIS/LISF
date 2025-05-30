!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
MODULE module_sf_aclsm_72

  implicit none

  ! SOIL PARAMETERS
  INTEGER :: SLCATS
  INTEGER, PARAMETER :: NSLTYPE=30
  CHARACTER(LEN=256) SLTYPE
  REAL, DIMENSION (1:NSLTYPE) :: &
       OC, WP, SAT, FC, INFRATE, SD, CL, SI

END MODULE module_sf_aclsm_72
