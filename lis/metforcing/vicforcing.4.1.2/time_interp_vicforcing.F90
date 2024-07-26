!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "vic412_atmos_forcing.h"
!BOP
! !ROUTINE: time_interp_vicforcing
! \label{time_interp_vicforcing412}
!
! !INTERFACE:
subroutine time_interp_vicforcing(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod,  only : LIS_forc, LIS_FORC_Base_State
  use LIS_logMod,         only : LIS_verify
  use vic_forcingMod,     only : vicforcing_struc

! !ARGUMENTS: 
  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: findex

! !DESCRIPTION: 
! This routine performs the temporal interpolation of the VIC-processed
! forcing data to the LIS model time-step.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the supplemental forcing scheme
!  \end{description}
!EOP

   integer :: t, k
   integer :: tindex, vindex

   type(ESMF_Field) :: tempField
   real,pointer     :: tempPtr(:)
   integer          :: status

! As mentioned previously, VIC-processed forcing data are in lock-step with
! the model or snow time-step.  Thus, when using these forcing data with VIC,
! temporal interpolation is currently not needed.

! Note that the variable indices like ATMOS_SNOWFLAG count from 0 (like in C),
! so I add 1 to them when using them within Fortran subroutines.

   vindex = ATMOS_SNOWFLAG + 1

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_SNOWFLAG%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status, 'Error: Enable Snowflag in the forcing variables list')
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, LIS_rc%ntiles(n)     
      tindex = LIS_domain(n)%tile(t)%index
      tempPtr(t) = vicforcing_struc(n)%metdata1(vindex,tindex)     
   end do

   vindex = ATMOS_PREC + 1

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Rainf%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, LIS_rc%ntiles(n)     
      tindex = LIS_domain(n)%tile(t)%index
      tempPtr(t) = vicforcing_struc(n)%metdata1(vindex,tindex)     
   end do

   vindex = ATMOS_AIR_TEMP + 1

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Tair%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, LIS_rc%ntiles(n)     
      tindex = LIS_domain(n)%tile(t)%index
      tempPtr(t) = vicforcing_struc(n)%metdata1(vindex,tindex)     
   end do

   vindex = ATMOS_WIND + 1

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_WIND%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status, 'Error: Enable Atmos Wind Speed in the forcing variables list')
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, LIS_rc%ntiles(n)     
      tindex = LIS_domain(n)%tile(t)%index
      tempPtr(t) = vicforcing_struc(n)%metdata1(vindex,tindex)     
   end do

   vindex = ATMOS_VPD + 1

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_VAPORPRESSDEFICIT%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status, 'Error: Enable Atmos Vapor Pressure Deficit in the forcing variables list')
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, LIS_rc%ntiles(n)     
      tindex = LIS_domain(n)%tile(t)%index
      tempPtr(t) = vicforcing_struc(n)%metdata1(vindex,tindex)     
   end do

   vindex = ATMOS_VP + 1

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_VAPORPRESS%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status, 'Error: Enable Atmos Vapor Pressure in the forcing variables list')
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, LIS_rc%ntiles(n)     
      tindex = LIS_domain(n)%tile(t)%index
      tempPtr(t) = vicforcing_struc(n)%metdata1(vindex,tindex)     
   end do

   vindex = ATMOS_PRESSURE + 1

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Psurf%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, LIS_rc%ntiles(n)     
      tindex = LIS_domain(n)%tile(t)%index
      tempPtr(t) = vicforcing_struc(n)%metdata1(vindex,tindex)     
   end do

   vindex = ATMOS_DENSITY + 1

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_DENSITY%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status, 'Error: Enable Atmos Density in the forcing variables list')
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, LIS_rc%ntiles(n)     
      tindex = LIS_domain(n)%tile(t)%index
      tempPtr(t) = vicforcing_struc(n)%metdata1(vindex,tindex)     
   end do

   vindex = ATMOS_SHORTWAVE + 1

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_SWdown%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, LIS_rc%ntiles(n)     
      tindex = LIS_domain(n)%tile(t)%index
      tempPtr(t) = vicforcing_struc(n)%metdata1(vindex,tindex)     
   end do

   vindex = ATMOS_LONGWAVE + 1

   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_LWdown%varname(1)),&
                      tempField,rc=status)
   call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')
   call ESMF_FieldGet(tempField,localDE=0, farrayPtr=tempPtr,rc=status)
   call LIS_verify(status)

   do t = 1, LIS_rc%ntiles(n)     
      tindex = LIS_domain(n)%tile(t)%index
      tempPtr(t) = vicforcing_struc(n)%metdata1(vindex,tindex)     
   end do

end subroutine time_interp_vicforcing
