!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

! NOTE:  Currently only 0.05 deg daily CHIRPSv2 data are supported
module CHIRPSv2_dataMod

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: CHIRPSv2_datainit

   ! Public types
   public :: CHIRPSv2data

   type, public :: chirpsv2datadec
      character*100 :: odir
      real          :: datares
      real, allocatable           :: rlat(:)
      real, allocatable           :: rlon(:)

      ! This is only used with upscale averaging
      integer, allocatable        :: n11(:)

      ! These are only used with budget interpolation
      integer, allocatable        :: n112(:,:)
      integer, allocatable        :: n122(:,:)
      integer, allocatable        :: n212(:,:)
      integer, allocatable        :: n222(:,:)     
      real,    allocatable        :: w112(:,:)
      real,    allocatable        :: w122(:,:)
      real,    allocatable        :: w212(:,:)
      real,    allocatable        :: w222(:,:)

      integer                     :: nc
      integer                     :: nr
   end type chirpsv2datadec
   
   type(chirpsv2datadec), allocatable :: chirpsv2data(:)

contains

   !---------------------------------------------------------------------------
   subroutine CHIRPSv2_datainit(i)

      ! Imports
      use ESMF
      use LVT_coreMod
      use LVT_logMod
      use LVT_histDataMod
      use LVT_timeMgrMod

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: i

      ! Local variables
      integer              :: status
      real                 :: gridDesci(50)
      integer              :: updoy, yr1,mo1,da1,hr1,mn1,ss1
      real                 :: upgmt
      character*10         :: time
      integer              :: ts

      if(.not.allocated(chirpsv2data)) then 
         allocate(chirpsv2data(LVT_rc%nDataStreams))
      endif

      ! Get top level CHIRPSv2 data directory
      call ESMF_ConfigGetAttribute(LVT_Config, chirpsv2data(i)%odir, &
           label='CHIRPSv2 data directory:', rc=status)
      call LVT_verify(status, 'CHIRPSv2 data directory: not defined')

      ! Allocate arrays on LVT grid
      allocate(chirpsv2data(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
      allocate(chirpsv2data(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))     

      ! Set CHIRPSv2 0.05 deg grid and map projection information.
      gridDesci(:) = 0
      gridDesci(1) = 0        ! Lat/lon
      gridDesci(2) = 7200     ! Points along latitude circle
      gridDesci(3) = 2000     ! Points along longitude circle
      gridDesci(4) =  -49.975 ! Latitude of first grid point
      gridDesci(5) = -179.975 ! Longitude of first grid point
      gridDesci(7) =   49.975 ! Latitude of last grid point
      gridDesci(8) =  179.975 ! Longitude of last grid point

      gridDesci(6) = 128       ! ???
      gridDesci(9) =  0.05     ! Longitudinal direction increment
      gridDesci(10) = 0.05     ! Latitude direction increment
      gridDesci(20) = 64       ! ???
      
      ! Set up interpolation data
      chirpsv2data(i)%datares = 0.05
      chirpsv2data(i)%nc = 7200
      chirpsv2data(i)%nr = 2000

      ! EMK...Use budget-bilinear interpolation if CHIRPSv2 is at         
      ! coarser resolution than the analysis grid; otherwise, use         
      ! upscale averaging.
      if (LVT_isAtAFinerResolution(chirpsv2data(i)%datares)) then

         ! Used only with budget interpolation
         allocate(chirpsv2data(i)%n112(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(chirpsv2data(i)%n122(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(chirpsv2data(i)%n212(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(chirpsv2data(i)%n222(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(chirpsv2data(i)%w112(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(chirpsv2data(i)%w122(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(chirpsv2data(i)%w212(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(chirpsv2data(i)%w222(LVT_rc%lnc*LVT_rc%lnr,25))

         chirpsv2data(i)%n112 = 0
         chirpsv2data(i)%n122 = 0
         chirpsv2data(i)%n212 = 0
         chirpsv2data(i)%n222 = 0
         chirpsv2data(i)%w112 = 0
         chirpsv2data(i)%w122 = 0
         chirpsv2data(i)%w212 = 0
         chirpsv2data(i)%w222 = 0

         call conserv_interp_input(gridDesci,LVT_rc%gridDesc,&
              LVT_rc%lnc*LVT_rc%lnr, &
              chirpsv2data(i)%rlat, chirpsv2data(i)%rlon,&
              chirpsv2data(i)%n112, chirpsv2data(i)%n122, &
              chirpsv2data(i)%n212, chirpsv2data(i)%n222, & 
              chirpsv2data(i)%w112, chirpsv2data(i)%w122, &
              chirpsv2data(i)%w212, chirpsv2data(i)%w222)
      else
         
         ! Used only with upscale averaging
         allocate(chirpsv2data(i)%n11(chirpsv2data(i)%nc*chirpsv2data(i)%nr))
         chirpsv2data(i)%n11 = 0

         call upscaleByAveraging_input(gridDesci,LVT_rc%gridDesc,&
              chirpsv2data(i)%nc*chirpsv2data(i)%nr, &
              LVT_rc%lnc*LVT_rc%lnr, &
              chirpsv2data(i)%n11)
      end if

   end subroutine CHIRPSv2_datainit
   
end module CHIRPSv2_dataMod
