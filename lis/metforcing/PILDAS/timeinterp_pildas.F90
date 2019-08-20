!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: timeinterp_pildas
! \label{timeinterp_pildas}
!
! !REVISION HISTORY:
! 25 Apr 2013: Eric Kemp, initial code
! 14 Jul 2016: Mahdi Navari - Modified for PILDAS
!
! !INTERFACE:

subroutine timeinterp_pildas(n,findex)

! !USES:
  use ESMF
  use LIS_coreMod,          only : LIS_rc, LIS_domain, LIS_localPet
  use LIS_metforcingMod,   only : LIS_forc, LIS_FORC_Base_State
  use LIS_FORC_AttributesMod
  use LIS_constantsMod,     only : LIS_CONST_SOLAR
  use LIS_timeMgrMod,       only : LIS_time2date
  use LIS_logMod,           only : LIS_logunit, LIS_verify, LIS_endrun
  use pildas_forcingMod,      only : pildas_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in):: n
  integer, intent(in):: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. Downward shortwave radiation is interpolated using a
!  zenith-angled based approach. Precipitation and longwave radiation
!  are not temporally interpolated, and the previous 1 hourly value
!  is used. All other variables are linearly interpolated between 
!  the 1 hourly blocks. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!EOP
  integer :: t,zdoy
  integer :: index1
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  real*8 :: btime
  real :: wt1,wt2,czb,cze,czm,gmt1,gmt2
  real :: zw1,zw2
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
  type(ESMF_Field)   :: swnetField,pardrField,pardfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer       :: swnet(:),pardr(:), pardf(:),snowf(:)

  gmt1 = pildas_struc(n)%gmt1
  gmt2 = pildas_struc(n)%gmt2
!-----------------------------------------------------------------------
!  Interpolate Data in Time
!-----------------------------------------------------------------------

  wt1=(pildas_struc(n)%pildastime2-LIS_rc%time)/ & 
      (pildas_struc(n)%pildastime2-pildas_struc(n)%pildastime1)
  wt2=1.0-wt1

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),tmpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),q2Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1),lwdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(1),psurfField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')

  if (LIS_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_CRainf%varname(1),cpcpField,&
          rc=status)
     call LIS_verify(status, &
          'Error: Enable CRainf in the forcing variables list')
  endif

  if (LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_Snowf%varname(1),snowfField,&
          rc=status)
     call LIS_verify(status, &
          'Error: Enable Snowf in the forcing variables list')
  endif

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  if (LIS_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LIS_verify(status)
  endif

  if (LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
     call LIS_verify(status)
  endif
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     zdoy=LIS_rc%doy
     call zterp(1,LIS_domain(n)%grid(index1)%lat,&
          LIS_domain(n)%grid(index1)%lon,&
          gmt1,gmt2,LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)             

     if (pildas_struc(n)%metdata1(3,index1).ne.LIS_rc%udef.and.&
          pildas_struc(n)%metdata2(3,index1).ne.LIS_rc%udef) then 
        swd(t) = pildas_struc(n)%metdata1(3,index1)*zw1+&
             pildas_struc(n)%metdata2(3,index1)*zw2
        if (swd(t).gt.LIS_CONST_SOLAR) then
           write(unit=LIS_logunit,fmt=*)'warning, sw radiation too high!!'
           write(unit=LIS_logunit,fmt=*)'it is',swd(t)
           write(unit=LIS_logunit,fmt=*)'pildasdata1=',&
                pildas_struc(n)%metdata1(3,index1)
           write(unit=LIS_logunit,fmt=*)'pildasdata2=',&
                pildas_struc(n)%metdata2(3,index1)
           write(unit=LIS_logunit,fmt=*)'zw1=',zw1,'zw2=',zw2
           swd(t) = LIS_CONST_SOLAR
           write(unit=LIS_logunit,fmt=*)'forcing set to ',swd(t) 
           stop
        endif
     endif
     
     if ((swd(t).ne.LIS_rc%udef).and.(swd(t).lt.0)) then
        if (swd(t).gt.-0.00001) then 
           swd(t) = 0.0 
        else
           write(LIS_logunit,*) &
                'ERR: timeinterp_pildas -- Stopping because ', & 
                'forcing not udef but lt0,'
           write(LIS_logunit,*)'ERR: timeinterp_pildas -- ', & 
                t,swd(t),pildas_struc(n)%metdata2(3,index1), & 
                ' (',LIS_localPet,')'
           call LIS_endrun
        endif
     endif
     
     if (swd(t).gt.LIS_CONST_SOLAR) then
        swd(t)=pildas_struc(n)%metdata2(3,index1)
     endif
  enddo

!-----------------------------------------------------------------------
! precip variable Block Interpolation
!-----------------------------------------------------------------------  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if(LIS_FORC_Snowf%selectOpt.eq.1) then         
       snowf(t) = pildas_struc(n)%metdata2(10,index1)
       pcp(t) = pildas_struc(n)%metdata2(8,index1) - & 
            pildas_struc(n)%metdata2(10,index1)
    else
       pcp(t) = pildas_struc(n)%metdata2(8,index1)
    endif
    if (LIS_FORC_CRainf%selectOpt.eq.1) then
       cpcp(t) = pildas_struc(n)%metdata2(9,index1)
     endif
  enddo

!-----------------------------------------------------------------------
!    Got Time Averaged LW
!-----------------------------------------------------------------------
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     lwd(t) = pildas_struc(n)%metdata1(4,index1)*wt1+  &
          pildas_struc(n)%metdata2(4,index1)*wt2
  enddo

!-----------------------------------------------------------------------
!     Linearly interpolate everything else
!-----------------------------------------------------------------------
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     tmp(t) = pildas_struc(n)%metdata1(1,index1)*wt1 + &
          pildas_struc(n)%metdata2(1,index1)*wt2
     q2(t)  = pildas_struc(n)%metdata1(2,index1)*wt1 + &
          pildas_struc(n)%metdata2(2,index1)*wt2
     uwind(t) = pildas_struc(n)%metdata1( 5,index1)*wt1 + &
          pildas_struc(n)%metdata2( 5,index1)*wt2
     vwind(t) = pildas_struc(n)%metdata1(6,index1)*wt1 + &
          pildas_struc(n)%metdata2(6,index1)*wt2
     psurf(t) = pildas_struc(n)%metdata1(7,index1)*wt1 + &
          pildas_struc(n)%metdata2(7,index1)*wt2

  enddo
end subroutine timeinterp_pildas
  
