!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: HAR_dataMod
! \label(HAR_dataMod)
!
! !INTERFACE:
module HAR_dataMod
   implicit none
   private

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This plugin supports the processing of precipitation data from 
!  High Asia Reanalysis (https://www.klima.tu-berlin.de/)
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  20 Apr 2018  Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
   public :: HAR_datainit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
   public :: HARdata
!EOP

   type, public :: hardatadec
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
   end type hardatadec
   
   type(hardatadec), allocatable :: hardata(:)

contains

!BOP
! 
! ROUTINE: HAR_datainit
! \label{HAR_datainit}
! 
! !INTERFACE: 
   subroutine HAR_datainit(i)
! !USES: 
      use ESMF
      use LVT_coreMod
      use LVT_logMod
      use LVT_histDataMod
      use LVT_timeMgrMod
! 
! !DESCRIPTION: 
!
!   This subroutine initializes and sets up the data structures required
!   for reading the HAR 10km data, including the computation of spatial 
!   interpolation weights. 
!  
!EOP

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

      if(.not.allocated(hardata)) then 
         allocate(hardata(LVT_rc%nDataStreams))
      endif

      ! Get top level HAR data directory
      call ESMF_ConfigGetAttribute(LVT_Config, hardata(i)%odir, &
           label='HAR data directory:', rc=status)
      call LVT_verify(status, 'HAR data directory: not defined')

      ! Allocate arrays on LVT grid
      allocate(hardata(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
      allocate(hardata(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))     

      ! Set HAR 0.05 deg grid and map projection information.
      gridDesci(:) = 0
      gridDesci(1) = 3       ! Lat/lon
      gridDesci(2) = 270     ! Points along latitude circle
      gridDesci(3) = 180     ! Points along longitude circle
      gridDesci(4) = 24.35494 
      gridDesci(5) = 72.91683
      gridDesci(6) = 8       ! ???
      gridDesci(7) = 30.0
      gridDesci(8) = 10 
      gridDesci(9) = 10 
      gridDesci(10) = 30.0
      gridDesci(11) = 87.0 

      ! Set up interpolation data
      hardata(i)%datares = 0.10
      hardata(i)%nc = 270
      hardata(i)%nr = 180

      ! EMK...Use budget-bilinear interpolation if HAR is at         
      ! coarser resolution than the analysis grid; otherwise, use         
      ! upscale averaging.
      if (LVT_isAtAFinerResolution(hardata(i)%datares)) then

         ! Used only with budget interpolation
         allocate(hardata(i)%n112(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(hardata(i)%n122(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(hardata(i)%n212(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(hardata(i)%n222(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(hardata(i)%w112(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(hardata(i)%w122(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(hardata(i)%w212(LVT_rc%lnc*LVT_rc%lnr,25))
         allocate(hardata(i)%w222(LVT_rc%lnc*LVT_rc%lnr,25))

         hardata(i)%n112 = 0
         hardata(i)%n122 = 0
         hardata(i)%n212 = 0
         hardata(i)%n222 = 0
         hardata(i)%w112 = 0
         hardata(i)%w122 = 0
         hardata(i)%w212 = 0
         hardata(i)%w222 = 0

         call conserv_interp_input(gridDesci,LVT_rc%gridDesc,&
              LVT_rc%lnc*LVT_rc%lnr, &
              hardata(i)%rlat, hardata(i)%rlon,&
              hardata(i)%n112, hardata(i)%n122, &
              hardata(i)%n212, hardata(i)%n222, & 
              hardata(i)%w112, hardata(i)%w122, &
              hardata(i)%w212, hardata(i)%w222)
      else
         
         ! Used only with upscale averaging
         allocate(hardata(i)%n11(hardata(i)%nc*hardata(i)%nr))
         hardata(i)%n11 = 0

         call upscaleByAveraging_input(gridDesci,LVT_rc%gridDesc,&
              hardata(i)%nc*hardata(i)%nr, &
              LVT_rc%lnc*LVT_rc%lnr, &
              hardata(i)%n11)
      end if

   end subroutine HAR_datainit
   
end module HAR_dataMod
