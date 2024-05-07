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
! !MODULE: TRMM3B42V7_obsMod
! \label(TRMM3B42V7_obsMod)
!
! !INTERFACE:
module TRMM3B42V7_dataMod
   
   ! Imports
   use ESMF

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: TRMM3B42V7_datainit

   ! Public types
   public :: TRMM3B42V7data

   type, public :: trmm3b42v7datadec
      character*100               :: odir
      real, allocatable           :: rlat(:)
      real, allocatable           :: rlon(:)
      integer, allocatable        :: n11(:)
      integer, allocatable        :: n12(:)
      integer, allocatable        :: n21(:)
      integer, allocatable        :: n22(:)     
      real,    allocatable        :: w11(:)
      real,    allocatable        :: w12(:)
      real,    allocatable        :: w21(:)
      real,    allocatable        :: w22(:)
      integer                     :: nc
      integer                     :: nr
      type(ESMF_TimeInterval)     :: ts
   end type trmm3b42v7datadec

   type(trmm3b42v7datadec), allocatable :: trmm3b42v7data(:)

contains

   !---------------------------------------------------------------------------
   subroutine TRMM3B42V7_datainit(i)
      
      ! Imports
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

      if(.not.allocated(trmm3b42v7data)) then 
         allocate(trmm3b42v7data(LVT_rc%nDataStreams))
      endif

      ! Get top level TRMM3B42V7 data directory
      call ESMF_ConfigGetAttribute(LVT_Config, trmm3b42v7data(i)%odir, &
         label='TRMM3B42V7 data directory:', rc=status)
      call LVT_verify(status, 'TRMM3B42V7 data directory: not defined')

      ! Get time resolution of TRMM3B42V7 data
      call ESMF_ConfigGetAttribute(LVT_Config, time, &
           label='TRMM3B42V7 data output interval:', rc=status)
      call LVT_verify(status, 'TRMM3B42V7 data output interval: not defined')
      call LVT_parseTimeString(time,ts)
      call LVT_update_timestep(LVT_rc,ts)
      
      ! Allocate arrays on LVT grid
      allocate(trmm3b42v7data(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
      allocate(trmm3b42v7data(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))     
      allocate(trmm3b42v7data(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
      allocate(trmm3b42v7data(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
      allocate(trmm3b42v7data(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
      allocate(trmm3b42v7data(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
      allocate(trmm3b42v7data(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
      allocate(trmm3b42v7data(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
      allocate(trmm3b42v7data(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
      allocate(trmm3b42v7data(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

      ! Set TRMM3B42V7 grid and map projection information.
      ! FIXME:  Support 8km fields
      gridDesci = 0            
      gridDesci(1) = 0          ! Lat/Lon
      gridDesci(2) = 1440       ! Points along latitude circle
      gridDesci(3) =  400       ! Points along longitude meridian
!EMK...TRMM3B42V7 data is from 0 to 360 E, but the LVT interpolation code
!wants longitude to range from -180 to 180 E.  The TRMM3B42V7 data are reordered
!in readTRMM3B42V7data to accomodate LVT, and the coordinates given below reflect
!that.
!      gridDesci(4) = -59.875    ! Latitude of first grid point
!      gridDesci(5) =   0.125    ! Longitude of first grid point
!      gridDesci(7) =  59.875    ! Latitude of last grid point
!      gridDesci(8) = 359.875    ! Longitude of last grid point
      gridDesci(4) = -49.875    ! Latitude of first grid point
      gridDesci(5) = -179.875    ! Longitude of first grid point
      gridDesci(7) =   49.875    ! Latitude of last grid point
      gridDesci(8) =  179.875    ! Longitude of last grid point

      gridDesci(6) = 128        ! ???
      gridDesci(9) =  0.250     ! Longitudinal direction increment
      gridDesci(10) = 0.250     ! Latitudinal direction increment
      gridDesci(20) = 64        ! ???

      ! Set up interpolation data
      trmm3b42v7data(i)%nc = 1440
      trmm3b42v7data(i)%nr =  400
      call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         trmm3b42v7data(i)%rlat, trmm3b42v7data(i)%rlon,&
         trmm3b42v7data(i)%n11, trmm3b42v7data(i)%n12, &
         trmm3b42v7data(i)%n21, trmm3b42v7data(i)%n22, & 
         trmm3b42v7data(i)%w11, trmm3b42v7data(i)%w12, &
         trmm3b42v7data(i)%w21, trmm3b42v7data(i)%w22)
      call ESMF_TimeIntervalSet(trmm3b42v7data(i)%ts, s = 10800, &
         rc=status)
      call LVT_verify(status)

   end subroutine TRMM3B42V7_datainit

end module TRMM3B42V7_dataMod
