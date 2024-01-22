!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: timeinterp_merra2_ac
! \label{timeinterp_merra2_ac}
!
! !REVISION HISTORY:
! 01 Jun 2022: Michel Bechtold, initial code (based on merra-2 data preprocessed
! to daily data)
! 17 Jan 2024: Louise Busschaert, AC71 implementation in NASA master
!
! !INTERFACE:

subroutine timeinterp_merra2_ac(n,findex)

! !USES:
  use ESMF
  use LIS_coreMod,          only : LIS_rc, LIS_domain, LIS_localPet
  use LIS_metforcingMod,    only : LIS_forc, LIS_FORC_Base_State
  use LIS_FORC_AttributesMod
  use LIS_constantsMod,     only : LIS_CONST_SOLAR
  use LIS_timeMgrMod,       only : LIS_tick
  use LIS_logMod,           only : LIS_logunit, LIS_verify, LIS_endrun
  use merra2_ac_forcingMod,    only : merra2_ac_struc
  use LIS_forecastMod, only : LIS_get_iteration_index
  use LIS_ran2_gasdev


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
  integer :: t,zdoy,k,kk
  integer :: index1
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  integer :: bss
  real*8  :: btime
  real    :: wt1,wt2,czb,cze,czm,gmt1,gmt2
  real    :: zw1,zw2,bts
  integer            :: status
  integer            :: mfactor,m
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
  type(ESMF_Field)   :: swnetField,pardrField,pardfField,RefHField
  type(ESMF_Field)   :: PREC_ac_Field,TMIN_ac_Field,TMAX_ac_Field,ETo_ac_Field
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),hlml(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer       :: swnet(:),pardr(:),pardf(:),snowf(:)
  real,pointer       :: PREC_ac(:),TMIN_ac(:),TMAX_ac(:),ETo_ac(:)
  real               :: tmp_real(2)

  btime=merra2_ac_struc(n)%merra2time1
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
  btime=merra2_ac_struc(n)%merra2time2
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

  wt1=1.0
  wt2=0.0

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

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_PREC_ac%varname(1),PREC_ac_Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable PREC_ac in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_TMIN_ac%varname(1),TMIN_ac_Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable TMIN_ac in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_TMAX_ac%varname(1),TMAX_ac_Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable TMAX_ac in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_ETo_ac%varname(1),ETo_ac_Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable ETo_ac in the forcing variables list')

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


  if (LIS_FORC_PREC_ac%selectOpt.eq.1) then
  call ESMF_FieldGet(PREC_ac_Field,localDE=0,farrayPtr=PREC_ac,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(TMIN_ac_Field,localDE=0,farrayPtr=TMIN_ac,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(TMAX_ac_Field,localDE=0,farrayPtr=TMAX_ac,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(ETo_ac_Field,localDE=0,farrayPtr=ETo_ac,rc=status)
  call LIS_verify(status)
  endif

  mfactor = LIS_rc%nensem(n)/merra2_ac_struc(n)%nIter

  
 
!-----------------------------------------------------------------------
! precip variable - no interp
!-----------------------------------------------------------------------  
  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n, k, index1, mfactor)
     ! Total precip, convective precip, and snowfall are a constant rate
     !  based on the first/current time and not on the latter time - dmm
        pcp(t) = merra2_ac_struc(n)%metdata1(kk,8,index1)

        !rand is sampled from a normal distribution with zero mean and 
        !variance of 1. Then scale it by stdev * x + mean 
        !The last member is kept unperturbed. 

        if(merra2_ac_struc(n)%usescalef==1.and. &
          merra2_ac_struc(n)%usepcpsampling.eq.1) then 
           
           if(merra2_ac_struc(n)%refstdev_ip(index1).ne.-9999.0) then 
           
              call nr_gasdev(merra2_ac_struc(n)%rseed(:,t), tmp_real)
              if(m.ne.LIS_rc%nensem(n)) then 
                 pcp(t) = tmp_real(1) * merra2_ac_struc(n)%refstdev_ip(index1) + & 
                      pcp(t)
              endif
           endif
        endif
           
        if ( pcp(t) < 0 ) then
           pcp(t) = 0
        endif
        if (LIS_FORC_CRainf%selectOpt.eq.1) then
           cpcp(t) = merra2_ac_struc(n)%metdata1(kk,9,index1)
        endif
        if(LIS_FORC_Snowf%selectOpt.eq.1) then 
           snowf(t) = merra2_ac_struc(n)%metdata1(kk,10,index1)
        endif
     enddo
  enddo
!-----------------------------------------------------------------------
! No interpolation 
!-----------------------------------------------------------------------
  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n, k, index1, mfactor)

        tmp(t) = merra2_ac_struc(n)%metdata1(kk,1,index1)
        q2(t)  = merra2_ac_struc(n)%metdata1(kk,2,index1)
        lwd(t) = merra2_ac_struc(n)%metdata1(kk,4,index1)
        uwind(t) = merra2_ac_struc(n)%metdata1(kk,5,index1)
        vwind(t) = merra2_ac_struc(n)%metdata1(kk,6,index1)
        psurf(t) = merra2_ac_struc(n)%metdata1(kk,7,index1)
        if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then 
           hlml(t) = merra2_ac_struc(n)%metdata1(kk,14,index1)
        endif


        swd(t) = merra2_ac_struc(n)%metdata1(kk,3,index1)
        if(LIS_FORC_SWnet%selectOpt.eq.1) then 
           swnet(t) = merra2_ac_struc(n)%metdata1(kk,11,index1)
           swnet(t) = min(swnet(t),swd(t))
        endif
        if(LIS_FORC_Pardr%selectOpt.eq.1) then 
           pardr(t) = merra2_ac_struc(n)%metdata1(kk,12,index1)
        endif
        if(LIS_FORC_Pardf%selectOpt.eq.1) then 
           pardf(t) = merra2_ac_struc(n)%metdata1(kk,13,index1)
        endif
        
        if ((swd(t).ne.LIS_rc%udef).and.(swd(t).lt.0)) then
           if (swd(t).gt.-0.00001) then 
              swd(t) = 0.0 
           else
              write(LIS_logunit,*) &
                   '[ERR] timeinterp_merra2_ac -- Stopping because ', & 
                   'forcing not udef but lt0,'
              write(LIS_logunit,*)'[ERR] timeinterp_merra2_ac -- ', & 
                   t,swd(t),merra2_ac_struc(n)%metdata2(kk,3,index1), & 
                   ' (',LIS_localPet,')'
              call LIS_endrun
           endif
        endif
        
        if (swd(t).gt.LIS_CONST_SOLAR) then
           swd(t)=LIS_CONST_SOLAR
        endif

        PREC_ac(t) = merra2_ac_struc(n)%metdata1(kk,15,index1)
        TMIN_ac(t) = merra2_ac_struc(n)%metdata1(kk,16,index1)
        TMAX_ac(t) = merra2_ac_struc(n)%metdata1(kk,17,index1)
        ETo_ac(t) = merra2_ac_struc(n)%metdata1(kk,18,index1)

     enddo
  enddo

end subroutine timeinterp_merra2_ac
  
