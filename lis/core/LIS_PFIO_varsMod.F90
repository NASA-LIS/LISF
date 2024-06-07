!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
module LIS_PFIO_varsMod
!
! !USES:
#if ( defined USE_PFIO )
   use MAPL
   use pFIO_UnlimitedEntityMod
!
   implicit none
!
! !DESCRIPTION:
! Defines global variables needed by MAPL/PFIO

   type pfio_t
      type(FileMetadata), allocatable :: fmd(:,:,:)
      integer,            allocatable :: hist_id(:,:,:)
      integer,            allocatable :: counter(:,:)
      logical,            allocatable :: first_time(:,:,:)
      type(FileMetadata), allocatable :: fmd_rst(:)
      integer,            allocatable :: rst_id(:)
   end type pfio_t
#else
   implicit none
   integer :: dummy_int
#endif
!EOP
!------------------------------------------------------------------------------
end module LIS_PFIO_varsMod
