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
! !MODULE: COAMPSout_dataMod
! \label(COAMPSout_dataMod)
!
! !INTERFACE:
module COAMPSout_dataMod
   
   ! Imports
   use ESMF

   ! Defaults
   implicit none
   private

   ! REVISION HISTORY: 
   ! 09 Dec 2022: Mahdi Navari; Initial Specification in LVT

   ! Public routines
   public :: COAMPSout_datainit

   ! Public types
   public ::  COAMPSoutdata 

   type, public :: COAMPSoutdatadec
      character*100               :: odir
      integer                     :: COAMPSnest_id
      real*8                      :: changetime1
      real*8                      :: changetime2
      real, allocatable           :: rlat(:)
      real, allocatable           :: rlon(:)
      integer                     :: nc
      integer                     :: nr
      type(ESMF_TimeInterval)     :: ts
   end type COAMPSoutdatadec

   type(COAMPSoutdatadec), allocatable :: COAMPSoutdata(:)

contains

   !---------------------------------------------------------------------------
   subroutine COAMPSout_datainit(i)
      
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
      !real                 :: gridDesci(50)
      !integer              :: updoy, yr1,mo1,da1,hr1,mn1,ss1
      !real                 :: upgmt
      !character*10         :: time
      !integer              :: ts

      if(.not.allocated(COAMPSoutdata)) then 
         allocate(COAMPSoutdata(LVT_rc%nDataStreams))
      endif

      ! Get top level COAMPSout data directory
      call ESMF_ConfigGetAttribute(LVT_Config, COAMPSoutdata(i)%odir, &
         label='COAMPS output forcing directory:', rc=status)
      call LVT_verify(status, 'COAMPS output forcing directory: not defined')
      call ESMF_ConfigGetAttribute(LVT_Config, COAMPSoutdata(i)%COAMPSnest_id, &
         label='COAMPS nest id:', rc=status)
      call LVT_verify(status, 'COAMPS nest id: not defined')

      call ESMF_TimeIntervalSet(COAMPSoutdata(i)%ts, s = 3600, &
         rc=status)
      call LVT_verify(status)

   end subroutine COAMPSout_datainit

end module COAMPSout_dataMod
