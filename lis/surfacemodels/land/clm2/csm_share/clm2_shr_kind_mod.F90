!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

MODULE clm2_shr_kind_mod

   !----------------------------------------------------------------------------
   ! precision/kind constants add data public
   !----------------------------------------------------------------------------
   public
   integer,parameter :: CLM2_SHR_KIND_R8 = selected_real_kind( 5) ! 8 byte real
   integer,parameter :: CLM2_SHR_KIND_R4 = selected_real_kind( 5) ! 4 byte real
   integer,parameter :: CLM2_SHR_KIND_RN = kind(1.0)              ! native real
   integer,parameter :: CLM2_SHR_KIND_I8 = selected_int_kind ( 5) ! 8 byte integer
   integer,parameter :: CLM2_SHR_KIND_I4 = selected_int_kind ( 5) ! 4 byte integer
   integer,parameter :: CLM2_SHR_KIND_IN = kind(1)                ! native integer

END MODULE clm2_shr_kind_mod
