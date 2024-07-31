!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! NOTE:  Currently only V05B IMERG data are supported
module IMERG_dataMod

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: IMERG_datainit

   ! Public types
   public :: IMERGdata

   type, public :: imergdatadec
      character*100 :: odir
      character*10  :: imergver, imergprd
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
   end type imergdatadec
   
   type(imergdatadec), allocatable :: imergdata(:)

contains
   !---------------------------------------------------------------------------
   subroutine IMERG_datainit(i)
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

      if(.not.allocated(imergdata)) then 
         allocate(imergdata(LVT_rc%nDataStreams))
      endif
      ! Get top level IMERG data directory
      call ESMF_ConfigGetAttribute(LVT_Config, imergdata(i)%odir, &
           label='IMERG data directory:', rc=status)
      call LVT_verify(status, 'IMERG data directory: not defined')
      ! Get IMERG product
      call ESMF_ConfigFindLabel(LVT_config,"IMERG product:",rc=status)
      if(status /= 0) then
         write(LVT_logunit,*) "[WARN] --------------------------------------------------------------"
         write(LVT_logunit,*) "[WARN] IMERG product not specified. Defaulting to 'final' version."
         write(LVT_logunit,*) "[WARN] --------------------------------------------------------------"
         imergdata(i)%imergprd = 'final'
      else
         call ESMF_ConfigGetAttribute(LVT_config,imergdata(i)%imergprd,rc=status)
      endif
      ! Get IMERG version
      call ESMF_ConfigFindLabel(LVT_config,"IMERG version:",rc=status)
      if(status /= 0) then
         write(LVT_logunit,*) "[WARN] --------------------------------------------------------------"
         write(LVT_logunit,*) "[WARN] IMERG version not specified. Defaulting to 'V06B' version."
         write(LVT_logunit,*) "[WARN] --------------------------------------------------------------"
         imergdata(i)%imergver = 'V06B'
      else
         call ESMF_ConfigGetAttribute(LVT_config,imergdata(i)%imergver,rc=status)
      endif

      ! Allocate arrays on LVT grid
      allocate(imergdata(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
      allocate(imergdata(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))     
      ! Set IMERG grid and map projection information.
      gridDesci(:) = 0
      gridDesci(1) = 0        ! Lat/lon
      gridDesci(2) = 3600     ! Points along latitude circle
      gridDesci(3) = 1800     ! Points along longitude circle
      gridDesci(4) =  -89.95  ! Latitude of first grid point
      gridDesci(5) = -179.95  ! Longitude of first grid point
      gridDesci(6) = 128    
      gridDesci(7) =   89.95  ! Latitude of last grid point
      gridDesci(8) =  179.95  ! Longitude of last grid point
      gridDesci(9) =  0.1     ! Longitudinal direction increment
      gridDesci(10) = 0.1     ! Latitude direction increment
      gridDesci(20) = 64      
      ! Set up interpolation data
      imergdata(i)%datares = 0.1
      imergdata(i)%nc = 3600
      imergdata(i)%nr = 1800
      ! Use budget-bilinear interpolation if IMERG resolution (0.1 deg)         
      ! is coarser than the analysis grid. Use upscale averaging if
      ! IMERG is finer than the analysis grid.
      if (LVT_isAtAFinerResolution(imergdata(i)%datares)) then
         ! Used only with budget interpolation
         allocate(imergdata(i)%n112(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(imergdata(i)%n122(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(imergdata(i)%n212(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(imergdata(i)%n222(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(imergdata(i)%w112(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(imergdata(i)%w122(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(imergdata(i)%w212(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(imergdata(i)%w222(LVT_rc%lnc*LVT_rc%lnr,25))
         imergdata(i)%n112 = 0
         imergdata(i)%n122 = 0
         imergdata(i)%n212 = 0
         imergdata(i)%n222 = 0
         imergdata(i)%w112 = 0
         imergdata(i)%w122 = 0
         imergdata(i)%w212 = 0
         imergdata(i)%w222 = 0
         call conserv_interp_input(gridDesci,LVT_rc%gridDesc,&
              LVT_rc%lnc*LVT_rc%lnr, &
              imergdata(i)%rlat, imergdata(i)%rlon,&
              imergdata(i)%n112, imergdata(i)%n122, &
              imergdata(i)%n212, imergdata(i)%n222, & 
              imergdata(i)%w112, imergdata(i)%w122, &
              imergdata(i)%w212, imergdata(i)%w222)
      else
         ! Used only with upscale averaging
         allocate(imergdata(i)%n11(imergdata(i)%nc*imergdata(i)%nr))
         imergdata(i)%n11 = 0
         call upscaleByAveraging_input(gridDesci,LVT_rc%gridDesc,&
              imergdata(i)%nc*imergdata(i)%nr, &
              LVT_rc%lnc*LVT_rc%lnr, &
              imergdata(i)%n11)
      end if
   end subroutine IMERG_datainit
end module IMERG_dataMod
