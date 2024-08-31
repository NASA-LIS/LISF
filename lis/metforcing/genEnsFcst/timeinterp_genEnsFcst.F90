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
! !ROUTINE: timeinterp_genEnsFcst
! \label{timeinterp_genEnsFcst}
!
! !REVISION HISTORY: 
! 27Sep2016 : K. Arsenault; implemented into LIS-7
!
! !INTERFACE:
subroutine timeinterp_genEnsFcst(n, findex)

! !USES:
  use ESMF
  use LIS_coreMod,         only : LIS_rc,LIS_domain
  use LIS_constantsMod,    only : LIS_CONST_SOLAR
  use LIS_logMod,          only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_timeMgrMod,      only : LIS_time2date
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod,   only : LIS_FORC_Base_State
  use genEnsFcst_forcingMod, only : genensfcst_struc
  use genEnsFcst_VariablesMod, only : forcopts

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
  integer :: zdoy
  real    :: zw1,zw2,czb,cze,czm

  integer :: t, m, k
  integer :: index1
  integer :: mfactor, tid
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
    call LIS_verify(status, 'Error: Enable Uwind in the forcing variables list')
    call ESMF_FieldGet(uwindField,localDE=0,farrayPtr=uwind,rc=status)
    call LIS_verify(status)
  endif

  if( forcopts%read_vwind ) then
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vwindField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Vwind in the forcing variables list')
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

#if 0
!- Time interpolate or apply average rate:
   wt1 = (genensfcst_struc%metforc_time2 - LIS_rc%time) / &
         (genensfcst_struc%metforc_time2 - genensfcst_struc%metforc_time1)
   wt2 = 1.0 - wt1
#endif 

  ! Ensemble member factor for this ensemble forcing dataset:
  mfactor = LIS_rc%nensem(n)/genensfcst_struc%max_ens_members

! --- SWDown ---
!- Perform zenith interpolation on solar radiation fields: 
  if( forcopts%read_swdown ) then
#if 0
    if( genensfcst_struc%zterp_flags ) then
      btime = genensfcst_struc%metforc_time1
      call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
      btime = genensfcst_struc%metforc_time2
      call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
    endif
#endif

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,genensfcst_struc%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index
#if 0
       if( genensfcst_struc%zterp_flags ) then
         zdoy = LIS_rc%doy
         call zterp( 0, LIS_domain(n)%grid(index1)%lat,&
              LIS_domain(n)%grid(index1)%lon, gmt1, gmt2, &
              LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)
         swd(t) = genensfcst_struc%metdata1(forcopts%index_swdown,m,index1)*zw1
       else
         swd(t) = genensfcst_struc%metdata1(forcopts%index_swdown,m,index1)
       endif
       if( swd(t) < 0 ) then
          write(LIS_logunit,*) "[ERR] SW radiation is negative!!"
          write(LIS_logunit,*) "[ERR] sw=", swd(t), "... negative"
          write(LIS_logunit,*) "[ERR] GenEnsFcst SWdown=", &
                genensfcst_struc%metdata2(forcopts%index_swdown,m,index1)
          call LIS_endrun
       end if
       if( swd(tid).gt.LIS_CONST_SOLAR ) then
          swd(tid) = genensfcst_struc%metdata1(forcopts%index_swdown,m,index1)
       endif
#endif
         swd(tid) = genensfcst_struc%metdata2(forcopts%index_swdown,m,index1)
        end do
      end do
    end do
  endif

!- Time Averaged Longwave, Block Interpolation
  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,genensfcst_struc%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index

          if( forcopts%read_airtmp ) then 
            if((genensfcst_struc%metdata1(forcopts%index_airtmp,m,index1).ne.LIS_rc%udef).and.&
               (genensfcst_struc%metdata2(forcopts%index_airtmp,m,index1).ne.LIS_rc%udef)) then
              airtmp(tid) = genensfcst_struc%metdata2(forcopts%index_airtmp,m,index1)
            endif
          endif
 
          if( forcopts%read_spechum ) then 
            if((genensfcst_struc%metdata1(forcopts%index_spechum,m,index1).ne.LIS_rc%udef).and.&
               (genensfcst_struc%metdata2(forcopts%index_spechum,m,index1).ne.LIS_rc%udef)) then
              spechum(tid) = genensfcst_struc%metdata2(forcopts%index_spechum,m,index1)
            endif
          endif

          if( forcopts%read_psurf ) then 
            if((genensfcst_struc%metdata1(forcopts%index_psurf,m,index1).ne.LIS_rc%udef).and.&
               (genensfcst_struc%metdata2(forcopts%index_psurf,m,index1).ne.LIS_rc%udef)) then
              psurf(tid) = genensfcst_struc%metdata2(forcopts%index_psurf,m,index1) 
            endif
          endif

          if( forcopts%read_lwdown ) then
            if((genensfcst_struc%metdata1(forcopts%index_lwdown,m,index1).ne.LIS_rc%udef).and.&
               (genensfcst_struc%metdata2(forcopts%index_lwdown,m,index1).ne.LIS_rc%udef)) then
              lwd(tid) = genensfcst_struc%metdata2(forcopts%index_lwdown,m,index1)
            endif
          endif

          if( forcopts%read_uwind ) then 
            if((genensfcst_struc%metdata1(forcopts%index_uwind,m,index1).ne.LIS_rc%udef).and.&
               (genensfcst_struc%metdata2(forcopts%index_uwind,m,index1).ne.LIS_rc%udef)) then
              uwind(tid) = genensfcst_struc%metdata2(forcopts%index_uwind,m,index1) 
            endif
          endif

          if( forcopts%read_vwind ) then 
            if((genensfcst_struc%metdata1(forcopts%index_vwind,m,index1).ne.LIS_rc%udef).and.&
               (genensfcst_struc%metdata2(forcopts%index_vwind,m,index1).ne.LIS_rc%udef)) then
              vwind(tid) = genensfcst_struc%metdata2(forcopts%index_vwind,m,index1) 
            endif
          endif

          if( forcopts%read_rainf ) then 
            if((genensfcst_struc%metdata1(forcopts%index_rainf,m,index1).ne.LIS_rc%udef).and.&
               (genensfcst_struc%metdata2(forcopts%index_rainf,m,index1).ne.LIS_rc%udef)) then
              prcp(tid) = genensfcst_struc%metdata2(forcopts%index_rainf,m,index1)
              if( prcp(tid) < 0.0 ) then
                prcp(tid) = 0.0
              endif
            endif
          endif
 
          if( forcopts%read_cpcp ) then 
            if((genensfcst_struc%metdata1(forcopts%index_cpcp,m,index1).ne.LIS_rc%udef).and.&
               (genensfcst_struc%metdata2(forcopts%index_cpcp,m,index1).ne.LIS_rc%udef)) then
              cpcp(tid) = genensfcst_struc%metdata2(forcopts%index_cpcp,m,index1)
              if( cpcp(tid) < 0.0 ) then
                cpcp(tid) = 0.0
              endif
            endif     
          endif

       end do
     enddo
   enddo   ! End tile loop

end subroutine timeinterp_genEnsFcst
