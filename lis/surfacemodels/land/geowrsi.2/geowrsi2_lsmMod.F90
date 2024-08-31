!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module geowrsi2_lsmMod
!BOP
!
! !MODULE: geowrsi2_lsmMod
!
! !DESCRIPTION:
!  Module for 1-D land model driver variable initialization
!  
! \begin{description}
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[outInterval]
!    output writing interval
!  \item[geowrsi2]
!   WRSI LSM specific variables
! \end{description} 
!
! !REVISION HISTORY:
! 31 Jul 2011: Brad Wind; Initial Definition
! 10 Jul 2013: KR Arsenault, JV Geiger;  Additional updates to code
! 25 Oct 2013: KR Arsenault; Added GeoWRSI2.0 model to LIS-7
!
! !USES:        
  use ESMF
  use LIS_timeMgrMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use geowrsi2_module
  use geowrsi2_physics_module, only : gCalcSOSmode
  use fbil_module, only : charN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: geowrsi2_lsm_ini
  public :: outVar
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: geowrsi2_struc
  public :: num_growing_seasons
  public :: geowrsi2_lsmRunMode
  public :: geowrsi2_CalcSOSlsmRunMode
  public :: geowrsi2_udef
!EOP

  type, public :: geowrsi2_type_dec

    real           :: ts                ! GeoWRSI Model Timestep
    integer        :: mon_dekad         ! Dekad counter per month (1-3)
    integer*4      :: year_dekad        ! Dekad counter per year (1-36)
    logical        :: end_soscalc       ! Flag for final SOS Calculation
    logical        :: eos_alarm         ! End-of-season (EOS) alarm
    
    type(charN)    :: inputParmFile     ! User Settings Input file (GeoWRSI-based)
    type(charN)    :: cropDir           ! Directory where crop parameter files reside
    character(len=LIS_CONST_PATH_LEN) :: rfile             ! Restart file path entry
    real           :: rstInterval       ! Restart file interval (e.g., dekad)

    integer        :: initTstepSeason   ! Average initial timestep of season
    integer        :: finalTstepSeason  ! Average final timestep of season
    integer        :: InitialYear       ! First growing season year of all years
    integer        :: LastCurrentYear   ! Final growing season year of all years

    integer        :: season_count      ! Total growing season counter 
    integer*2      :: lastSOScalcOfSeasonYr     ! Last SOS calc allowed in given year
    integer*4      :: lastSOScalcOfSeasonTStep  ! Last SOS calc allowed in given year
    integer*2      :: next_physics_act_year     ! Placeholder for activating physics
                                                !  for subsequent year

    type(geowrsi2dec), pointer :: wrsi(:)       ! variables per grid-cell

  end type geowrsi2_type_dec
  type(geowrsi2_type_dec), pointer :: geowrsi2_struc(:)

! Public terms:
  integer       :: num_growing_seasons  ! Total Number of growing seasons 
  character*10  :: geowrsi2_lsmRunMode  ! Type of GeoWRSI model runmode
  logical       :: geowrsi2_CalcSOSlsmRunMode  ! Same as above term
  real*4        :: geowrsi2_udef        ! Local GeoWRSI undefined value

  SAVE

contains

!BOP
! 
! !ROUTINE: geowrsi2_lsm_ini
! \label{geowrsi2_lsm_ini}
! 
! !INTERFACE:
  subroutine geowrsi2_lsm_ini()

! !USES:
   use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
   use LIS_timeMgrMod
   use LIS_coreMod,  only : LIS_rc
   use LIS_logMod,   only : LIS_logunit, LIS_verify, LIS_endrun

! !DESCRIPTION:        
!
!EOP
   implicit none
   integer                    :: n, t
   type(geowrsi2dec), pointer :: geowrsi2Pt
   character*3                :: fnest

! _____________________________________________________

   write(LIS_logunit,*) "Running WRSI LSM Option:"

   geowrsi2_udef = LIS_rc%udef
   allocate(geowrsi2_struc(LIS_rc%nnest))

   do n = 1, LIS_rc%nnest
      allocate(geowrsi2_struc(n)%inputParmFile%str)
      allocate(geowrsi2_struc(n)%cropDir%str)
   enddo

 ! Read config file option inputs:
   call geowrsi2_readcrd()

   do n = 1, LIS_rc%nnest
      geowrsi2_struc(n)%mon_dekad   = 0
      geowrsi2_struc(n)%year_dekad  = 0
      geowrsi2_struc(n)%end_soscalc = .false.
      geowrsi2_struc(n)%eos_alarm   = .false.
      geowrsi2_struc(n)%season_count= 0
      geowrsi2_struc(n)%lastSOScalcOfSeasonYr    = 0
      geowrsi2_struc(n)%lastSOScalcOfSeasonTStep = 0 
      geowrsi2_struc(n)%next_physics_act_year = 0

!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
      call LIS_update_timestep(LIS_rc, n, geowrsi2_struc(n)%ts)

!      call LIS_registerAlarm("GeoWRSI2 model alarm",&
!           geowrsi2_struc(n)%ts,&
!           geowrsi2_struc(n)%ts)

      write(fnest,'(i3.3)') n    

      call LIS_registerAlarm("GEOWRSI2 model alarm "//trim(fnest),LIS_rc%ts, &
                             864000.0, "dekad", dek_offset=0, when="end")

!      call LIS_registerAlarm("GeoWRSI2 restart alarm "//trim(fnest),&
!           geowrsi2_struc(n)%ts,&
!           geowrsi2_struc(n)%rstInterval)

   end do

!------------------------------------------------------------------------
! WRSI LSM Run Mode for Calculating SOS or WRSI:
!------------------------------------------------------------------------
   geowrsi2_CalcSOSlsmRunMode = .false.

   if( geowrsi2_lsmRunMode == "SOS" ) then
      write(LIS_logunit,*)"== WRSI Model: Running SOS Run-mode =="
!      write(*,*)"== WRSI Model: Running SOS Run-mode =="
      geowrsi2_CalcSOSlsmRunMode = .true.
      gCalcSOSmode = .true.

   elseif( geowrsi2_lsmRunMode == "WRSI" ) then
      write(LIS_logunit,*)"== WRSI Model: Running WRSI Run-mode == "
!      write(*,*)"== WRSI Model: Running WRSI Run-mode =="
      geowrsi2_CalcSOSlsmRunMode = .false.
      gCalcSOSmode = .false.
   else
      write(*,*)"== WRSI Model: INCORRECT RUN-MODE SELECTED =="
      write(*,*)"==       ... Select: SOS or WRSI           == "
      write(LIS_logunit,*)"== WRSI Model: INCORRECT RUN-MODE SELECTED == "
      write(LIS_logunit,*)"==       ... Select: SOS or WRSI           == "
      write(LIS_logunit,*)"== Stopping run-time ... "
      call LIS_endrun
   end if

!- Other WRSI LSM initializations needed:

   do n = 1, LIS_rc%nnest

      allocate(geowrsi2_struc(n)%wrsi(LIS_rc%npatch(n,LIS_rc%lsm_index)))

   !- Final SOS/SOSa fields to be written during SOS Calc Mode (KRA):
      do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
         allocate(geowrsi2_struc(n)%wrsi(t)%sos_write(num_growing_seasons))
         allocate(geowrsi2_struc(n)%wrsi(t)%sosa_write(num_growing_seasons))
         geowrsi2_struc(n)%wrsi(t)%sos_write(:)  = LIS_rc%udef
         geowrsi2_struc(n)%wrsi(t)%sosa_write(:) = LIS_rc%udef
!         geowrsi2_struc(n)%wrsi(t)%sos_write  = LIS_rc%udef
!         geowrsi2_struc(n)%wrsi(t)%sosa_write = LIS_rc%udef
      enddo

      LIS_sfmodel_struc(n)%nsm_layers = 1
      LIS_sfmodel_struc(n)%nst_layers = 1
      allocate(LIS_sfmodel_struc(n)%lyrthk(1))
      LIS_sfmodel_struc(n)%lyrthk(1) = 1.0
      LIS_sfmodel_struc(n)%ts = geowrsi2_struc(n)%ts
   enddo

  end subroutine geowrsi2_lsm_ini


  function outVar(realVar, geowrsi2)
  
     use fbil_module, only : CInt2

     integer(kind=2) :: outVar
  
     real*8,    intent(in) :: realVar
     integer*2, intent(in) :: geowrsi2
  
   ! Mask out the non-active pixels(i.e., where WRSI=253, 254 or 255)
   ! - Below code based on original GeoWRSI VB code:
     select case(geowrsi2)
       case(253)   ! "No Start"
!         outVar = -9999   ! GeoWRSI original assignment
         outVar = -9999   
       case(254)   ! "Yet to Start"
!         outVar = -9998   ! GeoWRSI original assignment
         outVar = -9999
       case(255)   ! "N/A"
!         outVar = -9997   ! GeoWRSI original assignment
         outVar = -9999
       case default
         outVar = CInt2(realVar)
     end select
  
  end function outVar

end module geowrsi2_lsmMod

