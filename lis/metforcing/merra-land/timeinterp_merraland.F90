!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: timeinterp_merraland
! \label{timeinterp_merraland}
!
! !REVISION HISTORY:
! 12 Oct 2009: Eric Kemp, initial code
! 22 Jul 2010: David Mocko, changed to hourly forcing
! 31 May 2013: David Mocko, changes to match MERRA-Land interpolation strategy
!
! !INTERFACE:

subroutine timeinterp_merraland(n,findex)

! !USES:
  use ESMF
  use LIS_coreMod,          only : LIS_rc, LIS_domain, LIS_localPet
  use LIS_metforcingMod,   only : LIS_forc, LIS_FORC_Base_State
  use LIS_FORC_AttributesMod
  use LIS_constantsMod,     only : LIS_CONST_SOLAR
  use LIS_timeMgrMod,       only : LIS_tick
  use LIS_logMod,           only : LIS_logunit, LIS_verify, LIS_endrun
  use merraland_forcingMod,      only : merraland_struc

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
  integer :: bss
  real*8  :: btime
  real    :: wt1,wt2,czb,cze,czm,gmt1,gmt2
  real    :: zw1,zw2,bts
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
  type(ESMF_Field)   :: swnetField,pardrField,pardfField,RefHField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),hlml(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer       :: swnet(:),pardr(:),pardf(:),snowf(:)

! Figuring GMT1 and 2 for the zterp routine (for SWdown):
  btime=merraland_struc(n)%merralandtime1
  byr = LIS_rc%yr
  bmo = LIS_rc%mo
  bda = LIS_rc%da
  bhr = LIS_rc%hr
  bmn = 30
  bss = 0
  if (LIS_rc%mn.lt.30) then
     bts = -(60*60)
  else
     bts = 0
  endif
  call LIS_tick(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn,bss,bts)

  btime=merraland_struc(n)%merralandtime2
  byr=LIS_rc%yr    !next hour
  bmo=LIS_rc%mo
  bda=LIS_rc%da
  bhr=LIS_rc%hr
  bmn=30
  bss=0
  if (LIS_rc%mn.lt.30) then
     bts=0
  else
     bts=60*60
  endif
  call LIS_tick(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn,bss,bts)

!-----------------------------------------------------------------------
!  Interpolate Data in Time
!-----------------------------------------------------------------------

  wt1=(merraland_struc(n)%merralandtime2-LIS_rc%time)/ & 
      (merraland_struc(n)%merralandtime2-merraland_struc(n)%merralandtime1)
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

  if (LIS_FORC_Pardr%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_Pardr%varname(1),pardrField,&
          rc=status)
     call LIS_verify(status, &
          'Error: Enable PARDR in the forcing variables list')
  endif

  if (LIS_FORC_Pardf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_Pardf%varname(1),pardfField,&
          rc=status)
     call LIS_verify(status, &
          'Error: Enable PARDF in the forcing variables list')
  endif

  if (LIS_FORC_Swnet%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_Swnet%varname(1),swnetField,&
          rc=status)
     call LIS_verify(status, &
          'Error: Enable SWnet in the forcing variables list')
  endif

  if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_Forc_Hgt%varname(1),RefHField,&
          rc=status)
     call LIS_verify(status, &
          'Error: Enable Forc_Hgt in the forcing variables list')
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

  if (LIS_FORC_Pardr%selectOpt.eq.1) then
     call ESMF_FieldGet(pardrField,localDE=0,farrayPtr=pardr,rc=status)
     call LIS_verify(status)
  endif

  if (LIS_FORC_Pardf%selectOpt.eq.1) then
     call ESMF_FieldGet(pardfField,localDE=0,farrayPtr=pardf,rc=status)
     call LIS_verify(status)
  endif

  if (LIS_FORC_Swnet%selectOpt.eq.1) then
     call ESMF_FieldGet(swnetField,localDE=0,farrayPtr=swnet,rc=status)
     call LIS_verify(status)
  endif
  
  if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
     call ESMF_FieldGet(RefHField,localDE=0,farrayPtr=hlml,rc=status)
     call LIS_verify(status)
  endif
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     zdoy=LIS_rc%doy
     call zterp(1,LIS_domain(n)%grid(index1)%lat,&
          LIS_domain(n)%grid(index1)%lon,&
          gmt1,gmt2,LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)
     if (merraland_struc(n)%metdata1(3,index1).ne.LIS_rc%udef.and.&
          merraland_struc(n)%metdata2(3,index1).ne.LIS_rc%udef) then 
        swd(t) = merraland_struc(n)%metdata1(3,index1)*zw1+&
             merraland_struc(n)%metdata2(3,index1)*zw2
        if (swd(t).gt.LIS_CONST_SOLAR) then
           write(unit=LIS_logunit,fmt=*) &
                          'warning, sw radiation too high!!'
           write(unit=LIS_logunit,fmt=*)'it is',swd(t)
           write(unit=LIS_logunit,fmt=*)'merralanddata1=',&
                merraland_struc(n)%metdata1(3,index1)
           write(unit=LIS_logunit,fmt=*)'merralanddata2=',&
                merraland_struc(n)%metdata2(3,index1)
           write(unit=LIS_logunit,fmt=*)'zw1=',zw1,'zw2=',zw2
           swd(t) = LIS_CONST_SOLAR
           write(unit=LIS_logunit,fmt=*)'forcing set to ',swd(t) 
        endif
! SWNET, PARDR, and PARDF should also be interpolated based
! on the solar zenith angle - dmm
        if(LIS_FORC_SWnet%selectOpt.eq.1) then 
           swnet(t) = merraland_struc(n)%metdata1(11,index1)*zw1+  &
             merraland_struc(n)%metdata2(11,index1)*zw2
           swnet(t) = min(swnet(t),swd(t))
        endif
        if(LIS_FORC_Pardr%selectOpt.eq.1) then 
           pardr(t) = merraland_struc(n)%metdata1(12,index1)*zw1+  &
             merraland_struc(n)%metdata2(12,index1)*zw2
        endif
        if(LIS_FORC_Pardf%selectOpt.eq.1) then 
           pardf(t) = merraland_struc(n)%metdata1(13,index1)*zw1+  &
             merraland_struc(n)%metdata2(13,index1)*zw2
        endif
     endif
     
     if ((swd(t).ne.LIS_rc%udef).and.(swd(t).lt.0)) then
        if (swd(t).gt.-0.00001) then 
           swd(t) = 0.0 
        else
           write(LIS_logunit,*) &
                'ERR: timeinterp_merraland -- Stopping because ', & 
                'forcing not udef but lt0,'
           write(LIS_logunit,*)'ERR: timeinterp_merraland -- ', & 
                t,swd(t),merraland_struc(n)%metdata2(3,index1), & 
                ' (',LIS_localPet,')'
           call LIS_endrun
        endif
     endif
     
     if (swd(t).gt.LIS_CONST_SOLAR) then
        swd(t)=merraland_struc(n)%metdata2(3,index1)
     endif
  enddo

!-----------------------------------------------------------------------
! precip variable - constant rate over the MERRA-Land hour
!-----------------------------------------------------------------------  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
! Total precip, convective precip, and snowfall are a constant rate
! based on the first/current time and not on the latter time - dmm
     pcp(t) = merraland_struc(n)%metdata1(8,index1)
     if ( pcp(t) < 0 ) then
        pcp(t) = 0
     endif
     if (LIS_FORC_CRainf%selectOpt.eq.1) then
        cpcp(t) = merraland_struc(n)%metdata1(9,index1)
     endif
     if(LIS_FORC_Snowf%selectOpt.eq.1) then 
       snowf(t) = merraland_struc(n)%metdata1(10,index1)
    endif
  enddo

!-----------------------------------------------------------------------
! LW down
!-----------------------------------------------------------------------
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
! Longwave down shouldn't be interpolated between MERRA-Land times - dmm
     lwd(t) = merraland_struc(n)%metdata1(4,index1)
  enddo

!-----------------------------------------------------------------------
! Linearly interpolate everything else (except winds)
!-----------------------------------------------------------------------
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     tmp(t) = merraland_struc(n)%metdata1(1,index1)*wt1 + &
          merraland_struc(n)%metdata2(1,index1)*wt2
     q2(t)  = merraland_struc(n)%metdata1(2,index1)*wt1 + &
          merraland_struc(n)%metdata2(2,index1)*wt2
! Wind in MERRA-Land is a time-averaged field and not interpolated - dmm
     uwind(t) = merraland_struc(n)%metdata1(5,index1)
     vwind(t) = merraland_struc(n)%metdata1(6,index1)
     psurf(t) = merraland_struc(n)%metdata1(7,index1)*wt1 + &
          merraland_struc(n)%metdata2(7,index1)*wt2
     if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then 
        hlml(t) = merraland_struc(n)%metdata1(14,index1)*wt1+  &
          merraland_struc(n)%metdata2(14,index1)*wt2
     endif
  enddo

end subroutine timeinterp_merraland
  
