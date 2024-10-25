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
! !MODULE: GDASforc_dataMod
! \label{GDASforc_dataMod}
!
! !INTERFACE:
module GDASforc_dataMod

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: GDASforc_datainit

   ! Public types
   public :: GDASforcdata

   type, public :: gdasforcdatadec
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

      integer                     :: ncold
      integer                     :: nrold
      integer                     :: mi
      real*8        :: griduptime1, griduptime2, griduptime3
      real*8        :: griduptime4, griduptime5, griduptime6
      logical       :: gridchange1, gridchange2, gridchange3
      logical       :: gridchange4, gridchange5, gridchange6
   end type gdasforcdatadec
   
   type(gdasforcdatadec), allocatable :: gdasforcdata(:)

contains

   !---------------------------------------------------------------------------
   subroutine GDASforc_datainit(i)

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

      if(.not.allocated(gdasforcdata)) then 
         allocate(gdasforcdata(LVT_rc%nDataStreams))
      endif

      ! Get top level GDASforc data directory
      call ESMF_ConfigGetAttribute(LVT_Config, gdasforcdata(i)%odir, &
           label='GDAS forcing data directory:', rc=status)
      call LVT_verify(status, 'GDAS forcing data directory: not defined')

      ! Allocate arrays on LVT grid
      allocate(gdasforcdata(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
      allocate(gdasforcdata(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))     

      gridDesci(:) = 0
      gridDesci(1) = 4
      gridDesci(2) = 192
      gridDesci(3) = 94
      gridDesci(4) = 88.542
      gridDesci(5) = 0
      gridDesci(6) = 128
      gridDesci(7) =  -88.542
      gridDesci(8) = -1.875
      gridDesci(9) = 1.875
      gridDesci(10) = 47
      gridDesci(20) = 0      
      ! Set up interpolation data
      gdasforcdata(i)%datares = 1.875
      gdasforcdata(i)%ncold = 192
      gdasforcdata(i)%nrold = 94

      ! EMK...Use budget-bilinear interpolation if GDASforc is at         
      ! coarser resolution than the analysis grid; otherwise, use         
      ! upscale averaging.
      if (LVT_isAtAFinerResolution(gdasforcdata(i)%datares)) then

         ! Used only with budget interpolation
         allocate(gdasforcdata(i)%n112(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(gdasforcdata(i)%n122(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(gdasforcdata(i)%n212(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(gdasforcdata(i)%n222(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(gdasforcdata(i)%w112(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(gdasforcdata(i)%w122(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(gdasforcdata(i)%w212(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(gdasforcdata(i)%w222(LVT_rc%lnc*LVT_rc%lnr,25))

         gdasforcdata(i)%n112 = 0
         gdasforcdata(i)%n122 = 0
         gdasforcdata(i)%n212 = 0
         gdasforcdata(i)%n222 = 0
         gdasforcdata(i)%w112 = 0
         gdasforcdata(i)%w122 = 0
         gdasforcdata(i)%w212 = 0
         gdasforcdata(i)%w222 = 0

         call conserv_interp_input(gridDesci,LVT_rc%gridDesc,&
              LVT_rc%lnc*LVT_rc%lnr, &
              gdasforcdata(i)%rlat, gdasforcdata(i)%rlon,&
              gdasforcdata(i)%n112, gdasforcdata(i)%n122, &
              gdasforcdata(i)%n212, gdasforcdata(i)%n222, & 
              gdasforcdata(i)%w112, gdasforcdata(i)%w122, &
              gdasforcdata(i)%w212, gdasforcdata(i)%w222)
      else
         
         ! Used only with upscale averaging
         allocate(gdasforcdata(i)%n11(gdasforcdata(i)%ncold*gdasforcdata(i)%nrold))
         gdasforcdata(i)%n11 = 0

         call upscaleByAveraging_input(gridDesci,LVT_rc%gridDesc,&
              gdasforcdata(i)%ncold*gdasforcdata(i)%nrold, &
              LVT_rc%lnc*LVT_rc%lnr, &
              gdasforcdata(i)%n11)
      end if

      call LVT_update_timestep(LVT_rc, 3*3600)

      ! This grid is good for some time in the 1990's.
      ! Look up the exact dates.
      yr1 = 1991
      mo1 = 01
      da1 = 01
      hr1 = 12
      mn1 = 0; ss1 = 0
      call LVT_date2time( gdasforcdata(i)%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
      
      yr1 = 2000
      mo1 = 01
      da1 = 24
      hr1 = 12
      mn1 = 0; ss1 = 0
      call LVT_date2time( gdasforcdata(i)%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
      
      yr1 = 2002     !grid update time ~ 0.469
      mo1 = 10
      da1 = 29
      hr1 = 12
      mn1 = 0; ss1 = 0
      call LVT_date2time(gdasforcdata(i)%griduptime3,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
      
      yr1 = 2005     !grid update time ~ 0.313
      mo1 = 05
      da1 = 31
      hr1 = 12
      mn1 = 0; ss1 = 0
      call LVT_date2time(gdasforcdata(i)%griduptime4,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
      
      yr1 = 2010     !grid update time ~ 0.205
      mo1 = 07
      da1 = 28
      hr1 = 12
      mn1 = 0; ss1 = 0
      call LVT_date2time(gdasforcdata(i)%griduptime5,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
      
      yr1 = 2015     !grid update time ~ 0.117
      mo1 = 01
      da1 = 14
      hr1 = 6
      mn1 = 0; ss1 = 0
      call LVT_date2time(gdasforcdata(i)%griduptime6,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
      
      gdasforcdata(i)%gridchange1 = .true.
      gdasforcdata(i)%gridchange2 = .true.
      gdasforcdata(i)%gridchange3 = .true.
      gdasforcdata(i)%gridchange4 = .true.
      gdasforcdata(i)%gridchange5 = .true.
      gdasforcdata(i)%gridchange6 = .true.
   end subroutine GDASforc_datainit
   
end module GDASforc_dataMod
