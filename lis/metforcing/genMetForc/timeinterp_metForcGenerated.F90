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
!BOP
! !ROUTINE: timeinterp_metForcGenerated
! \label{timeinterp_metForcGenerated}
!
! !REVISION HISTORY: 
! 21JAN2015 : K. Arsenault; implemented into LIS-7
!
! !INTERFACE:
subroutine timeinterp_metForcGenerated(n, findex)

! !USES:
  use ESMF
  use LIS_coreMod,         only : LIS_rc,LIS_domain
  use LIS_constantsMod,    only : LIS_CONST_SOLAR
  use LIS_logMod,          only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_timeMgrMod,      only : LIS_time2date
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod,   only : LIS_forc, LIS_FORC_Base_State
  use metForcGenerated_forcingMod, only : metForcGen_struc
  use metForcGen_VariablesMod, only : forcopts
  use LIS_forecastMod, only : LIS_get_iteration_index

! !ARGUMENTS: 
  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: findex

! !DESCRIPTION: 
!
!  Temporally interpolates the forcing data to the current model 
!  timestep. If needed, downward shortwave radiation is interpolated using a
!  zenith-angled based approach, depending on input forcing. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the forcing
!  \end{description}
! 
!EOP

! Time-interp weights:
  real    :: wt1, wt2
! Zterp terms:
  real    :: gmt1, gmt2
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime
  integer :: zdoy, t
  real    :: zw1,zw2,czb,cze,czm
  integer :: index1
  integer :: mfactor, m, k, kk
  integer :: status

! ESMF Fields and Pointer Arrays:
  type(ESMF_Field)   :: airtmpField, spechumField
  type(ESMF_Field)   :: swdField, lwdField
  type(ESMF_Field)   :: uwindField, vwindField
  type(ESMF_Field)   :: psurfField, prcpField, cpcpField
  real,pointer       :: airtmp(:), spechum(:)
  real,pointer       :: swd(:), lwd(:)
  real,pointer       :: uwind(:), vwind(:)
  real,pointer       :: psurf(:), prcp(:), cpcp(:)
! __________________________________________________________

! Get Meteorological Field  - ESMF State Get:
  if( forcopts%read_airtmp ) then
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),airtmpField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')
    call ESMF_FieldGet(airtmpField,localDE=0,farrayPtr=airtmp,rc=status)
    call LIS_verify(status)
  endif

  if( forcopts%read_spechum ) then
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),spechumField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')
    call ESMF_FieldGet(spechumField,localDE=0,farrayPtr=spechum,rc=status)
    call LIS_verify(status)
  endif

  if( forcopts%read_swdown ) then
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')
    call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
    call LIS_verify(status)
  endif

  if( forcopts%read_lwdown ) then
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1),lwdField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')
    call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
    call LIS_verify(status)
  endif

  if( forcopts%read_uwind ) then
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uwindField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')
    call ESMF_FieldGet(uwindField,localDE=0,farrayPtr=uwind,rc=status)
    call LIS_verify(status)
  endif

  if( forcopts%read_vwind ) then
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vwindField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')
    call ESMF_FieldGet(vwindField,localDE=0,farrayPtr=vwind,rc=status)
    call LIS_verify(status)
  endif

  if( forcopts%read_psurf ) then
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(1),psurfField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')
    call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
    call LIS_verify(status)
  endif

  if( forcopts%read_rainf ) then
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),prcpField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')
    call ESMF_FieldGet(prcpField,localDE=0,farrayPtr=prcp,rc=status)
    call LIS_verify(status)
    prcp  = 0.0
  endif

  if( forcopts%read_cpcp ) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CRainf%varname(1),cpcpField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')
     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LIS_verify(status)
     cpcp = 0.0
  endif

!- Time interpolate or apply average rate:
   wt1 = (metForcGen_struc%metforc_time2 - LIS_rc%time) / &
         (metForcGen_struc%metforc_time2 - metForcGen_struc%metforc_time1)
   wt2 = 1.0 - wt1

!- Metforcing ensemble factor:
   mfactor = LIS_rc%nensem(n)/metForcGen_struc%nIter

! --- SWDown ---
!- Perform zenith interpolation on solar radiation fields: 
  if( forcopts%read_swdown ) then
    if( metForcGen_struc%zterp_flags ) then
      btime = metForcGen_struc%metforc_time1
      call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
      btime = metForcGen_struc%metforc_time2
      call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
    endif
    do k=1,LIS_rc%ntiles(n)/mfactor
       do m=1,mfactor
          t = m + (k-1)*mfactor
          index1 = LIS_domain(n)%tile(t)%index
          kk = LIS_get_iteration_index(n, k, index1, mfactor)
          if( metForcGen_struc%zterp_flags ) then
            zdoy = LIS_rc%doy
            call zterp( 0, LIS_domain(n)%grid(index1)%lat,&
                 LIS_domain(n)%grid(index1)%lon, gmt1, gmt2, &
                 LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)
            swd(t) = metForcGen_struc%metdata1(kk,forcopts%index_swdown,index1)*zw1
          else
            swd(t) = metForcGen_struc%metdata1(kk,forcopts%index_swdown,index1)
          endif
          if( swd(t) < 0 ) then
            write(LIS_logunit,*) "[ERR] SW radiation is negative!!"
            write(LIS_logunit,*) "[ERR] sw=", swd(t), "... negative"
            write(LIS_logunit,*) "[ERR] LDT-generated SWdown=", &
                  metForcGen_struc%metdata2(kk,forcopts%index_swdown,index1)
            call LIS_endrun
          end if

          if( swd(t).gt.LIS_CONST_SOLAR ) then
             swd(t) = metForcGen_struc%metdata1(kk,forcopts%index_swdown,index1)
          endif
       enddo
    enddo
  endif

!- Time Averaged Longwave, Block Interpolation
   do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
       t = m + (k-1)*mfactor
       index1 = LIS_domain(n)%tile(t)%index
       kk = LIS_get_iteration_index(n, k, index1, mfactor)

       if( forcopts%read_airtmp ) then 
         if( forcopts%stat_airtmp == "inst" ) then
           airtmp(t) = wt1 * metForcGen_struc%metdata1(kk,forcopts%index_airtmp,index1) &
                     + wt2 * metForcGen_struc%metdata2(kk,forcopts%index_airtmp,index1)
         elseif( forcopts%stat_airtmp == "tavg" ) then
           airtmp(t) = metForcGen_struc%metdata2(kk,forcopts%index_airtmp,index1)
         endif
       endif

       if( forcopts%read_spechum ) then 
         if( forcopts%stat_spechum == "inst" ) then
          spechum(t) = wt1 * metForcGen_struc%metdata1(kk,forcopts%index_spechum,index1) &
                     + wt2 * metForcGen_struc%metdata2(kk,forcopts%index_spechum,index1)
         elseif( forcopts%stat_spechum == "tavg" ) then
           spechum(t) = metForcGen_struc%metdata2(kk,forcopts%index_spechum,index1)
         endif
       endif

       if( forcopts%read_psurf ) then 
         if( forcopts%stat_psurf == "inst" ) then
           psurf(t) = wt1 * metForcGen_struc%metdata1(kk,forcopts%index_psurf,index1) &
                    + wt2 * metForcGen_struc%metdata2(kk,forcopts%index_psurf,index1)
         elseif( forcopts%stat_psurf == "tavg" ) then
            psurf(t) = metForcGen_struc%metdata2(kk,forcopts%index_psurf,index1) 
         endif
       endif

       if( forcopts%read_lwdown ) then
         if( forcopts%stat_lwdown == "inst" ) then
           lwd(t) = wt1 * metForcGen_struc%metdata1(kk,forcopts%index_lwdown,index1) &
                  + wt2 * metForcGen_struc%metdata2(kk,forcopts%index_lwdown,index1)
         elseif( forcopts%stat_lwdown == "tavg" ) then
           lwd(t) = metForcGen_struc%metdata2(kk,forcopts%index_lwdown,index1)
         endif
       endif

       if( forcopts%read_uwind ) then 
         if( forcopts%stat_uwind == "inst" ) then
           uwind(t) = wt1 * metForcGen_struc%metdata1(kk,forcopts%index_uwind,index1) &
                    + wt2 * metForcGen_struc%metdata2(kk,forcopts%index_uwind,index1)
         elseif( forcopts%stat_uwind == "tavg" ) then
           uwind(t) = metForcGen_struc%metdata2(kk,forcopts%index_uwind,index1) 
         endif
       endif
       if( forcopts%read_vwind ) then 
          vwind(t) = 0.0   ! forcopts%index_vwind
       endif

       if( forcopts%read_rainf ) then 
         if( forcopts%stat_rainf == "inst" ) then
           prcp(t) = wt1 * metForcGen_struc%metdata1(kk,forcopts%index_rainf,index1) &
                   + wt2 * metForcGen_struc%metdata2(kk,forcopts%index_rainf,index1)
         elseif( forcopts%stat_rainf == "tavg" ) then
           prcp(t) = metForcGen_struc%metdata2(kk,forcopts%index_rainf,index1)
         elseif( forcopts%stat_rainf == "acc" ) then  ! Convert to rate
           prcp(t) = metForcGen_struc%metdata2(kk,forcopts%index_rainf,index1) &
                   / metForcGen_struc%ts
         endif
       endif
 
       if( forcopts%read_cpcp ) then 
         if( forcopts%stat_cpcp == "inst" ) then
           prcp(t) = wt1 * metForcGen_struc%metdata1(kk,forcopts%index_cpcp,index1) &
                   + wt2 * metForcGen_struc%metdata2(kk,forcopts%index_cpcp,index1)
         elseif( forcopts%stat_cpcp == "tavg" ) then
           cpcp(t) = metForcGen_struc%metdata2(kk,forcopts%index_cpcp,index1)
         elseif( forcopts%stat_rainf == "acc" ) then  ! Convert to rate
           cpcp(t) = metForcGen_struc%metdata2(kk,forcopts%index_cpcp,index1) &
                   / metForcGen_struc%ts
         endif
       endif     

     enddo 
   enddo   ! End tile loop

end subroutine timeinterp_metForcGenerated
