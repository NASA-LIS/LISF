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
! !MODULE: CMORPH_obsMod
! \label(CMORPH_obsMod)
!
! !INTERFACE:
module CMORPH_dataMod
   
   ! Imports
   use ESMF

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: CMORPH_datainit

   ! Public types
   public :: CMORPHdata

   type, public :: cmorphdatadec
      character*100               :: odir
      real*8                      :: changetime1
      real*8                      :: changetime2
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
   end type cmorphdatadec

   type(cmorphdatadec), allocatable :: cmorphdata(:)

contains

   !---------------------------------------------------------------------------
   subroutine CMORPH_datainit(i)
      
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

      if(.not.allocated(cmorphdata)) then 
         allocate(cmorphdata(LVT_rc%nDataStreams))
      endif

      ! Get top level CMORPH data directory
      call ESMF_ConfigGetAttribute(LVT_Config, cmorphdata(i)%odir, &
         label='CMORPH data directory:', rc=status)
      call LVT_verify(status, 'CMORPH data directory: not defined')

      ! Get time resolution of CMORPH data
      call ESMF_ConfigGetAttribute(LVT_Config, time, &
           label='CMORPH data output interval:', rc=status)
      call LVT_verify(status, 'CMORPH data output interval: not defined')
      call LVT_parseTimeString(time,ts)
      call LVT_update_timestep(LVT_rc,ts)
      
      ! Allocate arrays on LVT grid
      allocate(cmorphdata(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
      allocate(cmorphdata(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))     
      allocate(cmorphdata(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
      allocate(cmorphdata(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
      allocate(cmorphdata(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
      allocate(cmorphdata(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
      allocate(cmorphdata(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
      allocate(cmorphdata(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
      allocate(cmorphdata(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
      allocate(cmorphdata(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

      ! Set CMORPH grid and map projection information.
      ! FIXME:  Support 8km fields
      gridDesci = 0            
      gridDesci(1) = 0          ! Lat/Lon
      gridDesci(2) = 1440       ! Points along latitude circle
      gridDesci(3) =  480       ! Points along longitude meridian
!EMK...CMORPH data is from 0 to 360 E, but the LVT interpolation code
!wants longitude to range from -180 to 180 E.  The CMORPH data are reordered
!in readCMORPHdata to accomodate LVT, and the coordinates given below reflect
!that.
!      gridDesci(4) = -59.875    ! Latitude of first grid point
!      gridDesci(5) =   0.125    ! Longitude of first grid point
!      gridDesci(7) =  59.875    ! Latitude of last grid point
!      gridDesci(8) = 359.875    ! Longitude of last grid point
      gridDesci(4) =  -59.875    ! Latitude of first grid point
      gridDesci(5) = -179.875    ! Longitude of first grid point
      gridDesci(7) =   59.875    ! Latitude of last grid point
      gridDesci(8) =  179.875    ! Longitude of last grid point

      gridDesci(6) = 128        ! ???
      gridDesci(9) =  0.250     ! Longitudinal direction increment
      gridDesci(10) = 0.250     ! Latitudinal direction increment
      gridDesci(20) = 64        ! ???

      ! Set up interpolation data
      cmorphdata(i)%nc = 1440
      cmorphdata(i)%nr =  480
      call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         cmorphdata(i)%rlat, cmorphdata(i)%rlon,&
         cmorphdata(i)%n11, cmorphdata(i)%n12, &
         cmorphdata(i)%n21, cmorphdata(i)%n22, & 
         cmorphdata(i)%w11, cmorphdata(i)%w12, &
         cmorphdata(i)%w21, cmorphdata(i)%w22)
      call ESMF_TimeIntervalSet(cmorphdata(i)%ts, s = 10800, &
         rc=status)
      call LVT_verify(status)

   end subroutine CMORPH_datainit

end module CMORPH_dataMod
