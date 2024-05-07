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
#include "vic411_atmos_forcing.h"
!BOP
!
! !ROUTINE: vic411_f2t
! \label{vic411_f2t}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 14 Aug 2013; Shugong Wang, Modified for LIS-7
! !INTERFACE:
subroutine vic411_f2t(n)

! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_surface
  use LIS_timeMgrMod,     only : LIS_isAlarmRinging
  use LIS_metforcingMod,  only : LIS_FORC_State
  use LIS_FORC_AttributesMod 
  use LIS_logMod,         only : LIS_verify
  use vic411_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
! This routine transfers the LIS provided forcing onto the VIC 4.1.1
! model tiles. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

   type(ESMF_Field) :: tempField
   real,pointer     :: tempPtr(:)
   integer          :: status, iret
   integer          :: k, t, tid, VAR, NPATCH, NTILES
   logical          :: alarmCheck

   real*8 :: vic411_fsvp, tempC
   real :: tempvar

   real                          :: wind
   real,pointer     :: snowf_tempPtr(:),    &
        pressure_tempPtr(:), &
        windePtr(:),         &
        windnPtr(:)
   real,pointer   :: rainf_temp(:),    &
        air_temp_temp(:), &
        vp_temp(:)
   character*3 :: fnest 
   write(fnest, '(i3.3)') n

   k = vic411_struc(n)%snowstep
   NTILES = LIS_rc%ntiles(n)
   NPATCH = LIS_rc%npatch(n,LIS_rc%lsm_index)

   if ( LIS_FORC_Snowf%selectOpt == 1 ) then
      call ESMF_StateGet(LIS_FORC_State(n),&
                         trim(LIS_FORC_Snowf%varname(1)),&
                         tempField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(tempField,localDE=0, farrayPtr=snowf_tempPtr,rc=status)
      call LIS_verify(status)
   endif


   VAR = ATMOS_PREC
   allocate(rainf_temp(NTILES))

   call ESMF_StateGet(LIS_FORC_State(n),&
                      trim(LIS_FORC_Rainf%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   rainf_temp = tempPtr

   do t = 1, NTILES     
      ! If there is snowf add it to precipitation.  VIC does not use
      ! separate rainf and snowf.  It determines what to do with
      ! precipitation.
      if ( LIS_FORC_Snowf%selectOpt.eq.1) then
         if ( snowf_tempPtr(t) /= LIS_rc%udef ) then
            rainf_temp(t) = rainf_temp(t) + snowf_tempPtr(t)
         endif
      endif
   enddo

   do t = 1, NPATCH     
      tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
      if ( vic411_struc(n)%debugging_convert_units == 1 ) then
      ! convert from mm s-1 to mm
      rainf_temp(tid) = rainf_temp(tid) * LIS_rc%nts(n)
      call vic411_add_atmosdata(t, k, VAR, rainf_temp(tid))
      else
      ! rainf is already in mm from VIC
      call vic411_add_atmosdata(t, k, VAR, rainf_temp(tid))
      endif
   enddo


   VAR = ATMOS_AIR_TEMP
   allocate(air_temp_temp(NTILES))

   call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   if ( vic411_struc(n)%debugging_convert_units == 1 ) then
   ! convert air_temp from K to C
   air_temp_temp = tempPtr - VIC_KELVIN
   else
   ! air_temp already in C from VIC
   air_temp_temp = tempPtr
   endif

   do t = 1, NPATCH     
      tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
      call vic411_add_atmosdata(t, k, VAR, air_temp_temp(tid))
   enddo


   VAR = ATMOS_SNOWFLAG
   if ( LIS_FORC_Snowflag%selectOpt == 1 ) then
      call ESMF_StateGet(LIS_FORC_State(n),&
                         trim(LIS_FORC_SNOWFLAG%varname(1)),&
                         tempField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
      call LIS_verify(status)

      do t = 1, NPATCH     
          tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
         call vic411_add_atmosdata(t, k, VAR, tempPtr(tid))
      enddo
   else
      do t = 1, NPATCH     
         tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
         ! Assumption: air_temp is in C and rainf is in mm.
         if ( ( air_temp_temp(tid) +                     &
                vic411_struc(n)%vic(t)%min_Tfactor ) < &
              vic411_struc(n)%MAX_SNOW_TEMP .and.      &
              rainf_temp(tid) > 0 ) then
            ! EMK...Last argument should be real
            !call vic411_add_atmosdata(t, k, VAR, 1) ! TRUE
            call vic411_add_atmosdata(t, k, VAR, 1.) ! TRUE
         else
            ! EMK...Last argument should be real
            !call vic411_add_atmosdata(t, k, VAR, 0) ! FALSE
            call vic411_add_atmosdata(t, k, VAR, 0.) ! FALSE
         endif
      enddo
   endif


   VAR = ATMOS_WIND
   if ( LIS_FORC_Wind%selectOpt == 1 ) then
      call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_WIND%varname(1)),&
                         tempField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
      call LIS_verify(status)

      do t = 1, NPATCH     
         tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
         call vic411_add_atmosdata(t, k, VAR, tempPtr(tid))
      enddo
   else
      call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_WIND_E%varname(1)),&
                         tempField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(tempField,localDE=0, farrayPtr=windePtr,rc=status)
      call LIS_verify(status)

      call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_WIND_N%varname(1)),&
                         tempField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(tempField,localDE=0, farrayPtr=windnPtr,rc=status)
      call LIS_verify(status)

      do t = 1, NPATCH     
         tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
         wind = sqrt( windePtr(tid)**2 + windnPtr(tid)**2 )
         call vic411_add_atmosdata(t, k, VAR, wind)
      enddo
   endif


   VAR = ATMOS_PRESSURE
   call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=pressure_tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, NPATCH     
      tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
      call vic411_add_atmosdata(t, k, VAR, pressure_tempPtr(tid))
   enddo


   VAR = ATMOS_SHORTWAVE
   call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, NPATCH     
      tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
      call vic411_add_atmosdata(t, k, VAR, tempPtr(tid))
   enddo


   VAR = ATMOS_LONGWAVE
   call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status)
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, NPATCH     
      tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
      call vic411_add_atmosdata(t, k, VAR, tempPtr(tid))
   enddo


   VAR = ATMOS_VP
   allocate(vp_temp(NTILES))
   if ( LIS_FORC_VAPORPRESS%selectOpt == 1 ) then
      call ESMF_StateGet(LIS_FORC_State(n),&
                         trim(LIS_FORC_VAPORPRESS%varname(1)),&
                         tempField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
      call LIS_verify(status)

      vp_temp = tempPtr

      do t = 1, NPATCH     
         tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
         call vic411_add_atmosdata(t, k, VAR, tempPtr(tid))
      enddo
   else
      call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_QAIR%varname(1)),&
                         tempField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
      call LIS_verify(status)
      do t = 1, NPATCH     
         tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
         tempvar = tempPtr(tid) * pressure_tempPtr(tid) / VIC_EPS
         call vic411_add_atmosdata(t, k, VAR, tempvar)
         vp_temp(tid) = tempvar
      enddo
   endif


   VAR = ATMOS_VPD
   if ( LIS_FORC_VAPORPRESSDEFICIT%selectOpt == 1 ) then
      call ESMF_StateGet(LIS_FORC_State(n),&
                         trim(LIS_FORC_VAPORPRESSDEFICIT%varname(1)),&
                         tempField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
      call LIS_verify(status)

      do t = 1, NPATCH     
         tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
         call vic411_add_atmosdata(t, k, VAR, tempPtr(tid))
      enddo
   else
      do t = 1, NPATCH     
         tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
         ! Assumption: air_temp is in C.
         tempC = air_temp_temp(tid)
         tempvar = vic411_fsvp(tempC) - vp_temp(tid)
         call vic411_add_atmosdata(t, k, VAR, tempvar)
      enddo
   endif


   VAR = ATMOS_DENSITY
   if ( LIS_FORC_DENSITY%selectOpt == 1 ) then
      call ESMF_StateGet(LIS_FORC_State(n),&
                         trim(LIS_FORC_DENSITY%varname(1)),&
                         tempField,rc=status)
      call LIS_verify(status)
      call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
      call LIS_verify(status)

      do t = 1, NPATCH     
         tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
         call vic411_add_atmosdata(t, k, VAR, tempPtr(tid))
      enddo
   else
      do t = 1, NPATCH     
         tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
         ! Assumption: air_temp is in C.
         tempvar = pressure_tempPtr(tid) / &
                   (VIC_Rd * (VIC_KELVIN + air_temp_temp(tid)))
         call vic411_add_atmosdata(t, k, VAR, tempvar)
      enddo
   endif


   deallocate(air_temp_temp)
   deallocate(vp_temp)
   deallocate(rainf_temp)

   alarmCheck = LIS_isAlarmRinging(LIS_rc,"VIC411 model alarm "//trim(fnest))

   if ( alarmCheck ) then 
      ! Do not turn off ringer.  It is also needed in vic411_main.
      call vic411_nr_atmosdata(NPATCH)
      call vic411_write_atmos(NPATCH)
      vic411_struc(n)%snowstep = 1
   else
      vic411_struc(n)%snowstep = vic411_struc(n)%snowstep + 1
   endif


end subroutine vic411_f2t
