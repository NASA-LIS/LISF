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
! !ROUTINE: timeinterp_merra2
! \label{timeinterp_merra2}
!
! !REVISION HISTORY:
! 12 Nov 2015: KR Arsenault, added to LDT
!
! !INTERFACE:

subroutine timeinterp_merra2(n,findex)

! !USES:
  use ESMF
  use LDT_coreMod,          only : LDT_rc, LDT_domain, LDT_localPet
  use LDT_metforcingMod,    only : LDT_forc, LDT_FORC_Base_State
  use LDT_FORC_AttributesMod
  use LDT_constantsMod,     only : LDT_CONST_SOLAR
  use LDT_timeMgrMod,       only : LDT_tick
  use LDT_logMod,           only : LDT_logunit, LDT_verify, LDT_endrun
  use merra2_forcingMod,    only : merra2_struc

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
!   \item[LDT\_time2date](\ref{LDT_time2date}) \newline
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

  btime=merra2_struc(n)%merra2time1
  byr = LDT_rc%yr
  bmo = LDT_rc%mo
  bda = LDT_rc%da
  bhr = LDT_rc%hr
  bmn = 30
  bss = 0
  if (LDT_rc%mn.lt.30) then
     bts = -(60*60)
  else
     bts = 0
  endif
  call LDT_tick(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn,bss,bts)
  btime=merra2_struc(n)%merra2time2
  byr=LDT_rc%yr    !next hour
  bmo=LDT_rc%mo
  bda=LDT_rc%da
  bhr=LDT_rc%hr
  bmn=30
  bss=0
  if (LDT_rc%mn.lt.30) then
     bts=0
  else
     bts=60*60
  endif
  call LDT_tick(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn,bss,bts)

!-----------------------------------------------------------------------
!  Interpolate Data in Time
!-----------------------------------------------------------------------

  wt1=(merra2_struc(n)%merra2time2-LDT_rc%time)/ & 
      (merra2_struc(n)%merra2time2-merra2_struc(n)%merra2time1)
  wt2=1.0-wt1

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Tair%varname(1),tmpField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Tair in the forcing variables list')

  call ESMF_AttributeSet(tmpField,"Enabled",1,rc=status)
  call LDT_verify(status,&
       'error in ESMF_AttributeSet in timeinterp_merra2')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Qair%varname(1),q2Field,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Qair in the forcing variables list')
  call ESMF_AttributeSet(q2Field,"Enabled",1,rc=status)
  call LDT_verify(status,&
       'error in ESMF_AttributeSet in timeinterp_merra2')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_SWdown%varname(1),swdField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable SWdown in the forcing variables list')
  call ESMF_AttributeSet(swdField,"Enabled",1,rc=status)
  call LDT_verify(status,&
       'error in ESMF_AttributeSet in timeinterp_merra2')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_LWdown%varname(1),lwdField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable LWdown in the forcing variables list')
  call ESMF_AttributeSet(lwdField,"Enabled",1,rc=status)
  call LDT_verify(status,&
       'error in ESMF_AttributeSet in timeinterp_merra2')


  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Wind_E%varname(1),uField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Wind_E in the forcing variables list')
  call ESMF_AttributeSet(uField,"Enabled",1,rc=status)
  call LDT_verify(status,&
       'error in ESMF_AttributeSet in timeinterp_merra2')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Wind_N%varname(1),vField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Wind_N in the forcing variables list')
  call ESMF_AttributeSet(vField,"Enabled",1,rc=status)
  call LDT_verify(status,&
       'error in ESMF_AttributeSet in timeinterp_merra2')


  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Psurf%varname(1),psurfField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Psurf in the forcing variables list')
  call ESMF_AttributeSet(psurfField,"Enabled",1,rc=status)
  call LDT_verify(status,&
       'error in ESMF_AttributeSet in timeinterp_merra2')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Rainf%varname(1),pcpField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Rainf in the forcing variables list')
  call ESMF_AttributeSet(pcpField,"Enabled",1,rc=status)
  call LDT_verify(status,&
       'error in ESMF_AttributeSet in timeinterp_merra2')


  if (LDT_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_StateGet(LDT_FORC_Base_State(n,findex),&
          LDT_FORC_CRainf%varname(1),cpcpField,&
          rc=status)
     call LDT_verify(status, &
          'Error: Enable CRainf in the forcing variables list')
     call ESMF_AttributeSet(cpcpField,"Enabled",1,rc=status)
     call LDT_verify(status,&
          'error in ESMF_AttributeSet in timeinterp_merra2')
  endif

  if (LDT_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_StateGet(LDT_FORC_Base_State(n,findex),&
          LDT_FORC_Snowf%varname(1),snowfField,&
          rc=status)
     call LDT_verify(status, &
          'Error: Enable Snowf in the forcing variables list')
     call ESMF_AttributeSet(snowfField,"Enabled",1,rc=status)
     call LDT_verify(status,&
          'error in ESMF_AttributeSet in timeinterp_merra2')
  endif

  if (LDT_FORC_Pardr%selectOpt.eq.1) then
     call ESMF_StateGet(LDT_FORC_Base_State(n,findex),&
          LDT_FORC_Pardr%varname(1),pardrField,&
          rc=status)
     call LDT_verify(status, &
          'Error: Enable PARDR in the forcing variables list')
     call ESMF_AttributeSet(pardrField,"Enabled",1,rc=status)
     call LDT_verify(status,&
          'error in ESMF_AttributeSet in timeinterp_merra2')
  endif

  if (LDT_FORC_Pardf%selectOpt.eq.1) then
     call ESMF_StateGet(LDT_FORC_Base_State(n,findex),&
          LDT_FORC_Pardf%varname(1),pardfField,&
          rc=status)
     call LDT_verify(status, &
          'Error: Enable PARDF in the forcing variables list')
     call ESMF_AttributeSet(pardfField,"Enabled",1,rc=status)
     call LDT_verify(status,&
          'error in ESMF_AttributeSet in timeinterp_merra2')
  endif

  if (LDT_FORC_Swnet%selectOpt.eq.1) then
     call ESMF_StateGet(LDT_FORC_Base_State(n,findex),&
          LDT_FORC_Swnet%varname(1),swnetField,&
          rc=status)
     call LDT_verify(status, &
          'Error: Enable SWnet in the forcing variables list')
     call ESMF_AttributeSet(swnetField,"Enabled",1,rc=status)
     call LDT_verify(status,&
          'error in ESMF_AttributeSet in timeinterp_merra2')
  endif

  if (LDT_FORC_Forc_Hgt%selectOpt.eq.1) then
     call ESMF_StateGet(LDT_FORC_Base_State(n,findex),&
          LDT_FORC_Forc_Hgt%varname(1),RefHField,&
          rc=status)
     call LDT_verify(status, &
          'Error: Enable Forc_Hgt in the forcing variables list')
     call ESMF_AttributeSet(refhField,"Enabled",1,rc=status)
     call LDT_verify(status,&
          'error in ESMF_AttributeSet in timeinterp_merra2')
  endif

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LDT_verify(status)
  
  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LDT_verify(status)

  if (LDT_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LDT_verify(status)
  endif

  if (LDT_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
     call LDT_verify(status)
  endif

  if (LDT_FORC_Pardr%selectOpt.eq.1) then
     call ESMF_FieldGet(pardrField,localDE=0,farrayPtr=pardr,rc=status)
     call LDT_verify(status)
  endif

  if (LDT_FORC_Pardf%selectOpt.eq.1) then
     call ESMF_FieldGet(pardfField,localDE=0,farrayPtr=pardf,rc=status)
     call LDT_verify(status)
  endif

  if (LDT_FORC_Swnet%selectOpt.eq.1) then
     call ESMF_FieldGet(swnetField,localDE=0,farrayPtr=swnet,rc=status)
     call LDT_verify(status)
  endif
  
  if (LDT_FORC_Forc_Hgt%selectOpt.eq.1) then
     call ESMF_FieldGet(RefHField,localDE=0,farrayPtr=hlml,rc=status)
     call LDT_verify(status)
  endif

! - Time interpolation options:
  select case( LDT_rc%met_tinterp(findex) )

   case( "linear" )   ! Time-interp option
  
    ! Read in Shortwave radiation fields and interpolate for time:
    do t=1,LDT_rc%ntiles(n)
       index1 = LDT_domain(n)%tile(t)%index
 
       ! Call zenith angle interpolation to estimate weights:
       zdoy=LDT_rc%doy
       call zterp(1,LDT_domain(n)%grid(index1)%lat,&
          LDT_domain(n)%grid(index1)%lon,&
          gmt1,gmt2,LDT_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LDT_rc)

       if( LDT_forc(n,findex)%metdata1(3,index1).ne.LDT_rc%udef .and.&
         LDT_forc(n,findex)%metdata2(3,index1).ne.LDT_rc%udef ) then 

         swd(t) = LDT_forc(n,findex)%metdata1(3,index1)*zw1+&
                 LDT_forc(n,findex)%metdata2(3,index1)*zw2

         if (swd(t).gt.LDT_CONST_SOLAR) then
           write(unit=LDT_logunit,fmt=*) &
                          'warning, sw radiation too high!!'
           write(unit=LDT_logunit,fmt=*)'it is',swd(t)
           write(unit=LDT_logunit,fmt=*)'merra2data1=',&
                LDT_forc(n,findex)%metdata1(3,index1)
           write(unit=LDT_logunit,fmt=*)'merra2data2=',&
                LDT_forc(n,findex)%metdata2(3,index1)
           write(unit=LDT_logunit,fmt=*)'zw1=',zw1,'zw2=',zw2
           swd(t) = LDT_CONST_SOLAR
           write(unit=LDT_logunit,fmt=*)'forcing set to ',swd(t) 
         endif
         ! SWNET, PARDR, and PARDF should also be interpolated based
         ! on the solar zenith angle - dmm
         if(LDT_FORC_SWnet%selectOpt.eq.1) then 
            swnet(t) = LDT_forc(n,findex)%metdata1(11,index1)*zw1+  &
                       LDT_forc(n,findex)%metdata2(11,index1)*zw2
            swnet(t) = min(swnet(t),swd(t))
         endif
         if(LDT_FORC_Pardr%selectOpt.eq.1) then 
            pardr(t) = LDT_forc(n,findex)%metdata1(12,index1)*zw1+  &
              LDT_forc(n,findex)%metdata2(12,index1)*zw2
         endif
         if(LDT_FORC_Pardf%selectOpt.eq.1) then 
            pardf(t) = LDT_forc(n,findex)%metdata1(13,index1)*zw1+  &
              LDT_forc(n,findex)%metdata2(13,index1)*zw2
         endif
      endif
      ! Check for when SWdown is "undefined":
      if( (swd(t).ne.LDT_rc%udef).and.(swd(t).lt.0) ) then
         if (swd(t).gt.-0.00001) then 
            swd(t) = 0.0 
         else
            write(LDT_logunit,*) &
                 '[ERR] timeinterp_merra2 -- Stopping because ', & 
                 'forcing not udef but lt0,'
            write(LDT_logunit,*)'[ERR] timeinterp_merra2 -- ', & 
                 t,swd(t),LDT_forc(n,findex)%metdata2(3,index1), & 
                 ' (',LDT_localPet,')'
            call LDT_endrun
         endif
      endif
      ! SWdown is > Solar constant value:
      if (swd(t).gt.LDT_CONST_SOLAR) then
         swd(t)=LDT_forc(n,findex)%metdata2(3,index1)
      endif
    enddo

!-----------------------------------------------------------------------
! Precip variable - constant rate over the MERRA2 hour
!-----------------------------------------------------------------------  
    do t=1,LDT_rc%ntiles(n)
      index1 = LDT_domain(n)%tile(t)%index
      ! Total precip, convective precip, and snowfall are a constant rate
      ! based on the first/current time and not on the latter time - dmm
      pcp(t) = LDT_forc(n,findex)%metdata1(8,index1)
      if ( pcp(t) < 0 ) then
         pcp(t) = 0
      endif
      if (LDT_FORC_CRainf%selectOpt.eq.1) then
         cpcp(t) = LDT_forc(n,findex)%metdata1(9,index1)
      endif
      if(LDT_FORC_Snowf%selectOpt.eq.1) then 
        snowf(t) = LDT_forc(n,findex)%metdata1(10,index1)
      endif
    enddo
!-----------------------------------------------------------------------
! LW down
!-----------------------------------------------------------------------
    do t=1,LDT_rc%ntiles(n)
      index1 = LDT_domain(n)%tile(t)%index
      ! Longwave down shouldn't be interpolated between MERRA2 times - dmm
      lwd(t) = LDT_forc(n,findex)%metdata1(4,index1)
    enddo
!-----------------------------------------------------------------------
! Linearly interpolate everything else (except winds)
!-----------------------------------------------------------------------
    do t=1,LDT_rc%ntiles(n)
      index1 = LDT_domain(n)%tile(t)%index
      tmp(t) = LDT_forc(n,findex)%metdata1(1,index1)*wt1 + &
               LDT_forc(n,findex)%metdata2(1,index1)*wt2
      q2(t)  = LDT_forc(n,findex)%metdata1(2,index1)*wt1 + &
               LDT_forc(n,findex)%metdata2(2,index1)*wt2
      ! Wind in MERRA2 is a time-averaged field and not interpolated - dmm
      uwind(t) = LDT_forc(n,findex)%metdata1(5,index1)
      vwind(t) = LDT_forc(n,findex)%metdata1(6,index1)
      psurf(t) = LDT_forc(n,findex)%metdata1(7,index1)*wt1 + &
                 LDT_forc(n,findex)%metdata2(7,index1)*wt2
      if(LDT_FORC_Forc_Hgt%selectOpt.eq.1) then 
         hlml(t) = LDT_forc(n,findex)%metdata1(14,index1)*wt1+  &
                   LDT_forc(n,findex)%metdata2(14,index1)*wt2
      endif
    enddo
!-----------------------------------------------------------------------

   case( "none" ) ! No time interpolation option turned on

    do t=1,LDT_rc%ntiles(n)
      index1 = LDT_domain(n)%tile(t)%index

       ! SWD:
       if( LDT_forc(n,findex)%metdata1(3,index1).ne.LDT_rc%udef ) then
         swd(t) = LDT_forc(n,findex)%metdata1(3,index1)

         if( swd(t).gt.LDT_CONST_SOLAR ) then
           swd(t) = LDT_CONST_SOLAR 
         endif
         if( swd(t) .lt. 0.0 ) then
           swd(t) = 0.0
         endif
       endif
       ! SWNET, PARDR, and PARDF should also be interpolated based
       ! on the solar zenith angle - dmm
       if(LDT_FORC_SWnet%selectOpt.eq.1) then
         swnet(t) = LDT_forc(n,findex)%metdata1(11,index1)
         swnet(t) = min(swnet(t),swd(t))
       endif
       if(LDT_FORC_Pardr%selectOpt.eq.1) then
         pardr(t) = LDT_forc(n,findex)%metdata1(12,index1)
       endif
       if(LDT_FORC_Pardf%selectOpt.eq.1) then
         pardf(t) = LDT_forc(n,findex)%metdata1(13,index1)
       endif

       ! Total precip, convective precip, and snowfall are a constant rate
       ! based on the first/current time and not on the latter time - dmm
       pcp(t) = LDT_forc(n,findex)%metdata1(8,index1)
       if ( pcp(t) < 0 ) then
          pcp(t) = 0
       endif
       if (LDT_FORC_CRainf%selectOpt.eq.1) then
          cpcp(t) = LDT_forc(n,findex)%metdata1(9,index1)
       endif
       if(LDT_FORC_Snowf%selectOpt.eq.1) then
         snowf(t) = LDT_forc(n,findex)%metdata1(10,index1)
       endif

       ! All other forcing variables:
       lwd(t) = LDT_forc(n,findex)%metdata1(4,index1)
       tmp(t) = LDT_forc(n,findex)%metdata1(1,index1)
       q2(t)  = LDT_forc(n,findex)%metdata1(2,index1)
       uwind(t) = LDT_forc(n,findex)%metdata1(5,index1)
       vwind(t) = LDT_forc(n,findex)%metdata1(6,index1)
       psurf(t) = LDT_forc(n,findex)%metdata1(7,index1)
       if(LDT_FORC_Forc_Hgt%selectOpt.eq.1) then
          hlml(t) = LDT_forc(n,findex)%metdata1(14,index1)
       endif
    end do

   case default   ! unrecognizable time_interp option
     write(LDT_logunit,*)"[ERR] Do not currently recognize the following"
     write(LDT_logunit,*)"  time interp option: ",trim(LDT_rc%met_tinterp(findex))
     write(LDT_logunit,*)" LDT stopping runtime ..."
     call LDT_endrun
  end select    ! End time interp case block

end subroutine timeinterp_merra2
  
