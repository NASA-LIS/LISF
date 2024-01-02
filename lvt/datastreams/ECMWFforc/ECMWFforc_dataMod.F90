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
! 
! !MODULE: ECMWFforc_dataMod
! \label(ECMWFforc_dataMod)
!
! !INTERFACE:
module ECMWFforc_dataMod

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: ECMWFforc_datainit

   ! Public types
   public :: ECMWFforcdata

   type, public :: ecmwfforcdatadec
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
   end type ecmwfforcdatadec
   
   type(ecmwfforcdatadec), allocatable :: ecmwfforcdata(:)

contains

   !---------------------------------------------------------------------------
   subroutine ECMWFforc_datainit(i)

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

      if(.not.allocated(ecmwfforcdata)) then 
         allocate(ecmwfforcdata(LVT_rc%nDataStreams))
      endif

      ! Get top level ECMWFforc data directory
      call ESMF_ConfigGetAttribute(LVT_Config, ecmwfforcdata(i)%odir, &
           label='ECMWF forcing data directory:', rc=status)
      call LVT_verify(status, 'ECMWF forcing data directory: not defined')

      ! Allocate arrays on LVT grid
      allocate(ecmwfforcdata(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
      allocate(ecmwfforcdata(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))     

      ! Set ECMWFforc 0.05 deg grid and map projection information.
      gridDesci(:) = 0
      gridDesci(1) = 0        ! Lat/lon
      gridDesci(2) = 1440     ! Points along latitude circle
      gridDesci(3) = 601     ! Points along longitude circle
      gridDesci(4) =  90.00 ! Latitude of first grid point
      gridDesci(5) = -180.00 ! Longitude of first grid point
      gridDesci(7) =  -60.00 ! Latitude of last grid point
      gridDesci(8) =  179.75 ! Longitude of last grid point

      gridDesci(6) = 128       ! ???
      gridDesci(9) =  0.25     ! Longitudinal direction increment
      gridDesci(10) = 0.25     ! Latitude direction increment
      gridDesci(20) = 64       ! ???
      
      ! Set up interpolation data
      ecmwfforcdata(i)%datares = 0.25
      ecmwfforcdata(i)%nc = 1440
      ecmwfforcdata(i)%nr = 601

      ! EMK...Use budget-bilinear interpolation if ECMWFforc is at         
      ! coarser resolution than the analysis grid; otherwise, use         
      ! upscale averaging.
      if (LVT_isAtAFinerResolution(ecmwfforcdata(i)%datares)) then

         ! Used only with budget interpolation
         allocate(ecmwfforcdata(i)%n112(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(ecmwfforcdata(i)%n122(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(ecmwfforcdata(i)%n212(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(ecmwfforcdata(i)%n222(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(ecmwfforcdata(i)%w112(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(ecmwfforcdata(i)%w122(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(ecmwfforcdata(i)%w212(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(ecmwfforcdata(i)%w222(LVT_rc%lnc*LVT_rc%lnr,25))

         ecmwfforcdata(i)%n112 = 0
         ecmwfforcdata(i)%n122 = 0
         ecmwfforcdata(i)%n212 = 0
         ecmwfforcdata(i)%n222 = 0
         ecmwfforcdata(i)%w112 = 0
         ecmwfforcdata(i)%w122 = 0
         ecmwfforcdata(i)%w212 = 0
         ecmwfforcdata(i)%w222 = 0

         call conserv_interp_input(gridDesci,LVT_rc%gridDesc,&
              LVT_rc%lnc*LVT_rc%lnr, &
              ecmwfforcdata(i)%rlat, ecmwfforcdata(i)%rlon,&
              ecmwfforcdata(i)%n112, ecmwfforcdata(i)%n122, &
              ecmwfforcdata(i)%n212, ecmwfforcdata(i)%n222, & 
              ecmwfforcdata(i)%w112, ecmwfforcdata(i)%w122, &
              ecmwfforcdata(i)%w212, ecmwfforcdata(i)%w222)
      else
         
         ! Used only with upscale averaging
         allocate(ecmwfforcdata(i)%n11(ecmwfforcdata(i)%nc*ecmwfforcdata(i)%nr))
         ecmwfforcdata(i)%n11 = 0

         call upscaleByAveraging_input(gridDesci,LVT_rc%gridDesc,&
              ecmwfforcdata(i)%nc*ecmwfforcdata(i)%nr, &
              LVT_rc%lnc*LVT_rc%lnr, &
              ecmwfforcdata(i)%n11)
      end if

      call LVT_update_timestep(LVT_rc, 3*3600)
   end subroutine ECMWFforc_datainit
   
end module ECMWFforc_dataMod
