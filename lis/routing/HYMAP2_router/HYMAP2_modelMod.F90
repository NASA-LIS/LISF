!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !REVISION HISTORY: 
! 19 Jan 2016: Augusto Getirana, Inclusion of the local inertia formulation 
!                                and adaptive time step. 
!  2 Mar 2020: Augusto Getirana, Fixed river/floodplain storage
! 27 Apr 2020: Augusto Getirana,  Added support for urban drainage
!
module HYMAP2_modelMod
  
  use LIS_logMod,     only : LIS_logunit
  implicit none

  public :: HYMAP2_dwi
  public :: HYMAP2_linear_reservoir_lis
  public :: HYMAP2_no_reservoir_lis
  public :: HYMAP2_calc_runoff
  public :: HYMAP2_calc_fldstg
  public :: HYMAP2_calc_rivout_kine
  public :: HYMAP2_calc_rivout_iner
  public :: HYMAP2_calc_stonxt
  public :: HYMAP2_dynstp
  public :: HYMAP2_date2frac
  public :: HYMAP2_calc_evap_fld  
  public :: HYMAP2_get_elevation_profile
  public :: HYMAP2_get_volume_profile   
  
contains
  ! ================================================
  ! ================================================  
  subroutine HYMAP2_dwi(runoff0,basflw0,roffdwi_ratio,&
                        basfdwi_ratio,runoff1,basflw1,roffdwi,basfdwi)
    implicit none
    real,    intent(in)    :: runoff0       !input runoff [mm.idt-1]
    real,    intent(in)    :: basflw0       !input baseflow [mm.idt-1]
    real,    intent(in)    :: basfdwi_ratio !deep water infiltration ratio from baseflow [-]
    real,    intent(in)    :: roffdwi_ratio !deep water infiltration ratio from surface runoff [-]
    real,    intent(out)   :: runoff1       !output runoff [mm.idt-1]
    real,    intent(out)   :: basflw1       !output baseflow [mm.idt-1]
    real,    intent(out)   :: basfdwi       !deep water infiltration from baseflow [mm.idt-1]
    real,    intent(out)   :: roffdwi       !deep water infiltration from surface runoff [mm.idt-1]
    
    roffdwi=runoff0*roffdwi_ratio
    runoff1=runoff0*(1.-max(0.,min(roffdwi_ratio,1.)))

    basfdwi=basflw0*basfdwi_ratio
    basflw1=basflw0*(1.-max(0.,min(basfdwi_ratio,1.)))

  end subroutine HYMAP2_dwi
  ! ================================================
  ! ================================================  
  subroutine HYMAP2_linear_reservoir_lis(dt,rmis,runoff0,basflw0, &
                              trnoff,tbsflw,cntime,roffsto,&
                              basfsto,runoff)
    ! ================================================
    ! to   This routine simulates surface and underground linear reservoirs  
    ! by   Augusto GETIRANA 
    ! on   6th Dec 2010
    ! at   LEGOS/CNES, Toulouse, France
    ! Adapted for flow routing implementation in LIS 09 Nov 2011
    ! Adapted for single grid cell computation in 19 Feb 2016
    ! ================================================
    
    implicit none
    
    ! ================================================
    real,    intent(in)    :: dt      !input time step length [s]
    real,    intent(in)    :: rmis    !missing data
    real,    intent(in)    :: runoff0 !runoff   before linear reservoir [m3/s]
    real,    intent(in)    :: basflw0 !baseflow before linear reservoir [m3/s]
    real,    intent(in)    :: trnoff  !runoff   concentration time parameter [day]
    real,    intent(in)    :: tbsflw  !baseflow concentration time parameter [day]
    real,    intent(in)    :: cntime  !concentration time [s]    
    real,    intent(inout) :: roffsto !runoff   reservoir storage [m3]
    real,    intent(inout) :: basfsto !baseflow reservoir storage [m3]
    real,    intent(out)   :: runoff  !runoff+baseflow after  linear reservoir [m3/s]

    real                   :: vroff, vbflw
    real                   :: qrun, qbas
    integer                :: i,j,jj
    ! ================================================
    vroff=0.0
    vbflw=0.0
    qrun=0.0
    qbas=0.0
    runoff=0.0
    !runoff
    !update reservoir volumes with input data (runoff)
    roffsto = roffsto + runoff0*dt
    !compute temporary volumes (m3)
    qrun=roffsto/(trnoff*cntime)  
    vroff=roffsto - qrun*dt
    !check if reservoir is empty
    if(vroff<0.)then
      qrun=roffsto/dt
      roffsto=0.
    else
      roffsto=vroff
    endif
    !baseflow
    !update reservoir volumes with input data (baseflow)
    basfsto = basfsto + basflw0*dt
    !compute temporary volumes (m3)
    qbas=basfsto/(tbsflw*86400.)
    vbflw=basfsto-qbas*dt
    !check if reservoir is empty
    if(vbflw<0.)then
      qbas=basfsto/dt
      basfsto=0.
    else
      basfsto=vbflw
    endif
    runoff=qrun+qbas
    
  end subroutine HYMAP2_linear_reservoir_lis
  ! ================================================
  ! ================================================  
  subroutine HYMAP2_no_reservoir_lis(dt,rmis,runoff0,basflw0,&
                          roffsto,basfsto,runoff)
    ! ================================================
    !      This routine sums runoff and baseflow variables
    ! by   Augusto GETIRANA 
    ! on   28th Dec 2010
    ! at   LEGOS/CNES, Toulouse, France
    ! Adapted for single grid cell computation in 19 Feb 2016
    ! ================================================

    implicit none

    ! ================================================
    real,    intent(in)  :: dt               !input time step length [s]
    real,    intent(in)  :: rmis             !missing data
    real,    intent(in)  :: runoff0 !runoff          before linear reservoir [m3/s]
    real,    intent(in)  :: basflw0 !baseflow        before linear reservoir [m3/s]
    real,    intent(out) :: roffsto !runoff   reservoir storage [m3]
    real,    intent(out) :: basfsto !baseflow reservoir storage [m3]
    real,    intent(out) :: runoff  !runoff+baseflow after  linear reservoir [m3/s]
    ! ================================================

    if(runoff0/=rmis.and.basflw0/=rmis)then
       runoff=runoff0+basflw0
    else
       runoff=0.
    endif
    roffsto=0.
    basfsto=0.
    
  end subroutine HYMAP2_no_reservoir_lis
  ! ================================================  
  ! ================================================  
  subroutine HYMAP2_calc_runoff(dt,fldfrc,runoff,rivsto,fldsto)
    ! ================================================
    ! to   calculate runoff to stolage
    ! by   Augusto Getirana after Dai YAMAZAKI
    ! on   10 Feb 2011
    ! at   LEGOS/CNES, Toulouse, France
    ! ================================================
    implicit none
    real,    intent(in)  ::  dt
    real,    intent(in)  ::  fldfrc !area fraction
    real,    intent(in)  ::  runoff !runoff of the grid [m3/s]
    real,    intent(out) ::  rivsto !river storage [m3]
    real,    intent(out) ::  fldsto !floodplain storage [m3]
    real                 ::  rrivrof,rfldrof
    
    !ag (24Mar2022)
    !Fix negative fldsto
    if(runoff<0.)then
      rfldrof=max(-fldsto,runoff*fldfrc*dt)
      rrivrof=max(-rivsto,runoff*dt-rfldrof)
    else
      rrivrof=runoff*(1.-fldfrc)*dt
      rfldrof=runoff*fldfrc*dt
    endif
    rivsto=max(0.,rivsto+rrivrof)
    fldsto=max(0.,fldsto+rfldrof)
    
  end subroutine HYMAP2_calc_runoff
  ! ================================================
  !ag (2Mar2020)
  ! ================================================
  subroutine HYMAP2_calc_fldstg(nz,grarea,rivlen,rivwth,            &
                         rivstomax,fldstomax,fldhgt,&
                         rivsto,fldsto,rivdph,flddph,flddph1,&
                         fldelv1,fldwth,fldfrc,fldare,       &
                         rivelv,rivare)
    use LIS_logMod
    ! ================================================
    ! to   calculate river and floodplain staging
    ! by   Augusto Getirana
    ! on   21 Aug 2018
    ! at   GSFC/NASA, Greenbelt, USA
    ! ================================================
    implicit none
  
    integer, intent(in)    :: nz           !number of stages in the sub-grid discretization
    real,    intent(in)    :: grarea       !area of the grid [m2]
    real,    intent(in)    :: rivlen       !channel length [m]
    real,    intent(in)    :: rivwth       !river width [m]
    real,    intent(in)    :: rivstomax
    real,    intent(in)    :: fldstomax(nz)!maximum floodplain storage [m3]
    real,    intent(in)    :: fldhgt(nz)
    real,    intent(in)    :: rivelv       !river bed elevation [m]
    
    real,    intent(in)    :: rivare       !river surface area [m2]
    real,    intent(inout) :: rivsto       !river storage [m3]
    real,    intent(inout) :: fldsto       !flood plain storage [m3]
    real,    intent(inout) :: rivdph       !river depth [m]
    real,    intent(inout) :: flddph       !floodplain depth [m]
    real,    intent(inout) :: flddph1      !mean floodplain depth [m]
    real,    intent(inout) :: fldelv1      !mean floodplain elevation [m]
    real,    intent(inout) :: fldwth       !floodplain width [m]
    real,    intent(out)   :: fldfrc       !area fraction
    real,    intent(out)   :: fldare
    integer                :: i
    real*8                 :: rstoall
    real*8                 :: h1,h2,v1,v2
    ! ================================================
    rstoall=rivsto+fldsto
    
    if(rstoall>rivstomax.and.grarea>rivare)then
      i=1
      !do while(rstoall>fldstomax(i).and.i<nz)
      do while(rstoall>fldstomax(i))
        i=i+1
        if(i>nz)then
          flddph=fldhgt(nz)+(rstoall-fldstomax(nz))/grarea
          fldfrc=1.
          goto 1
        endif
      enddo
      if(i>1)then
        h1=fldhgt(i-1);v1=fldstomax(i-1)
        h2=fldhgt(i);v2=fldstomax(i)
      elseif(i==1)then
        h1=0.;v1=rivstomax
        h2=fldhgt(1);v2=fldstomax(1)
      else
        write(LIS_logunit,*)"[ERR] Please check [HYMAP2_calc_fldstg] "
        call LIS_endrun()
      endif        
      flddph=h1+(h2-h1)*(rstoall-v1)/(v2-v1)
      fldfrc=min(max((rivare/grarea)+real(i-1)/real(nz)+(1./real(nz))*(rstoall-v1)/(v2-v1),0.),1.)
1 continue   
      rivsto=rivstomax+rivlen*rivwth*flddph
      !ag(27Apr2020) for cases where there's no river (e.g., urban areas)
      if(rivwth>0.and.rivlen>0)then
        rivdph=rivsto/rivlen/rivwth
      else
        rivdph=0
      endif
      fldsto=max(rstoall-rivsto,0.)
      fldare=grarea*fldfrc
      fldwth=fldare/rivlen
      flddph1=fldsto/fldwth/rivlen
    else
      rivsto=rstoall
      rivdph=max(0.,rstoall/rivare)
      fldsto=0.
      flddph=0.
      if(rstoall>0.)then
        fldfrc=min(max(rivare/grarea,0.),1.)
      else
        fldfrc=0.
      endif
      fldare=fldfrc*grarea
      fldwth=0.
      flddph1=0.
    endif
    fldelv1=rivelv+rivdph-flddph1
    
  end subroutine HYMAP2_calc_fldstg
  ! ================================================
  ! ================================================  
  subroutine HYMAP2_calc_rivout_kine(outlet,dt,rivelv,rivelv_down,nxtdst,&
                              rivwth,sfcelv,rivlen,manval,slpmin,&
                              rivsto,rivdph,rivout,rivvel)
    use LIS_logMod
    ! ================================================
    ! Calculate discharge, kinematic wave
    ! Augusto Getirana
    ! 10 Feb 2011
    ! LEGOS/CNES, Toulouse, France
    ! Adapted for flow routing implementation in LIS 23 Mar 2012
    ! Adapted for single grid cell computation in 19 Feb 2016
    ! ================================================
    implicit none    
    integer, intent(in)  :: outlet      !outlet flag: 0 - river; 1 - ocean
    real,    intent(in)  :: dt          !time step length [sec]
    real,    intent(in)  :: rivelv      !river bed elevation [m]
    real,    intent(in)  :: rivelv_down !downstream river bed elevation [m]
    real,    intent(in)  :: nxtdst      !distance to next grif [m]
    real,    intent(in)  :: rivwth      !river width [m]
    real,    intent(in)  :: rivsto      !floodplain storage [m3]
    real,    intent(in)  :: rivdph      !depth [m]
    real,    intent(in)  :: rivlen      !river length [m]
    real,    intent(in)  :: manval      !maning coefficient for rivers
    real,    intent(in)  :: slpmin      !minimum slope
    real,    intent(out) :: rivout      !outflow from the floodplain reservoir [m3/s]
    real,    intent(out) :: rivvel      !velocity [m/s]
    real,    intent(out) :: sfcelv      !water surface elevation
    real                 :: rslope,rarea,rvel,rhydrad
    ! ================================================
    sfcelv=rivelv+rivdph
    if(outlet==0)then
      rslope=(rivelv-rivelv_down)/nxtdst
      rslope=max(rslope,slpmin)
      rvel=(1./manval)*rslope**0.5*rivdph**(2./3)
      rivvel=rvel
      rarea=rivwth*rivdph
      rivout=rarea*rvel
      rivout=min(rivout,rivsto/dt)
    elseif(outlet==1)then
      rslope=slpmin
      rarea=rivwth*rivdph
      rvel=(1./manval)*rslope**0.5*rivdph**(2./3)
      rivvel=rvel
      rivout=rarea*rvel
      rivout=min(rivout,rivsto/dt)
    else
     write(LIS_logunit,*)"[ERR] Wrong outlet id",outlet
     call LIS_endrun()
    endif

  end subroutine HYMAP2_calc_rivout_kine
  ! ================================================
  ! ================================================      
  subroutine HYMAP2_calc_rivout_iner(outlet,dt,rivelv,rivelv_down,elevtn,nxtdst,        &
                              rivwth,rivsto,rivsto_down,rivdph,rivdph_down,rivlen,manval,    &
                              grv,rivout,rivvel,sfcelv,&
                              rivout_pre,rivdph_pre,rivdph_pre_down)
    use LIS_logMod
    ! ================================================
    ! Calculate discharge, local inertia 
    ! Augusto Getirana
    ! 16 Apr 2014
    ! Adapted for flow routing implementation in LIS on 5 Aug 2015
    ! Adapted for single grid cell computation in 19 Feb 2016
    ! ================================================
    implicit none
    integer, intent(in)  :: outlet          !outlet flag: 0 - river; 1 - ocean
    real,    intent(in)  :: dt              !time step length [sec]
    real,    intent(in)  :: rivelv          !river/floodplain bed elevation [m]
    real,    intent(in)  :: rivelv_down     !downstream river/floodplain bed elevation [m]
    real,    intent(in)  :: elevtn          !surface elevation [m]
    real,    intent(in)  :: nxtdst          !distance to next grid [m]
    real,    intent(in)  :: rivwth          !river/floodplain width [m]
    real,    intent(in)  :: rivsto          !river/floodplain storage [m3]
    real,    intent(in)  :: rivsto_down     !downstream river/floodplain storage [m3]
    real,    intent(in)  :: rivdph          !river/floodplain depth [m]
    real,    intent(in)  :: rivdph_down     !dowstream river/floodplain depth [m]
    real,    intent(in)  :: rivlen          !river/floodplain length [m]
    real,    intent(in)  :: manval          !manning coefficient for river/floodplain
    real,    intent(in)  :: grv             !gravity [m/s2]
    real,    intent(out) :: rivout          !river/floodplain outflow  [m3/s]
    real,    intent(out) :: rivvel          !river/floodplain flow velocity [m/s]
    real,    intent(out) :: sfcelv          !water surface elevation
    real,    intent(out) :: rivout_pre      !previous outflow [m3/s]
    real,    intent(out) :: rivdph_pre      !previous river depth [m]
    real,    intent(in)  :: rivdph_pre_down !previous downstream river depth [m]

    real                 :: sfcelv_down     !downstream water surface elevation
    real                 :: sfcelv_pre      !previous step water surface elevation (t-1) [m]
    real                 :: sfcelv_pre_down !previous step downstream water surface elevation (t-1) [m]
    real                 :: dslope, darea,dvel, dhydrad
    real                 :: dflw,dout_pre,dflw_pre,dflw_imp

    ! ================================================
    sfcelv=rivelv+rivdph
    sfcelv_down=rivelv_down+rivdph_down
    sfcelv_pre=rivelv+rivdph_pre
    sfcelv_pre_down=rivelv_down+rivdph_pre_down

    if(outlet==0)then
      dslope=(sfcelv-sfcelv_down)/nxtdst
      dflw=max(sfcelv,sfcelv_down)-rivelv
      darea=rivwth*dflw
      dflw_pre=max(sfcelv_pre,sfcelv_pre_down)-rivelv
      dflw_imp=sqrt(dflw*dflw_pre)
      if(dflw_imp<=0.)dflw_imp=dflw
        
      !if(dflw_imp>1.e-5.and.darea>1.e-5)then
      if(dflw_imp>0..and.darea>0.)then
        dout_pre=rivout_pre/rivwth
        rivout=rivwth*(dout_pre+grv*dt*dflw_imp*dslope)/&
             (1.+grv*dt*manval**2.*abs(dout_pre)*dflw_imp**(-7./3))
        if(rivout>=0.)then
          rivout=min(rivout,rivsto/dt)
        else
          rivout=-min(abs(rivout),rivsto_down/dt)
        endif
        rivvel=rivout/darea
      else
        rivout=0.
        rivvel=0.
      endif
      rivout_pre=rivout
      rivdph_pre=rivdph
      
    elseif(outlet==1)then
      !ag (20Feb2020)
      !dslope=max(0.,(sfcelv-elevtn)/nxtdst)
      dslope=(sfcelv-sfcelv_down)/nxtdst
      !ag(7Apr2021)
      dflw=max(sfcelv,sfcelv_down)-rivelv
      !dflw=rivdph
      darea=rivwth*dflw
      dflw_pre=max(sfcelv_pre,sfcelv_pre_down)-rivelv
      !dflw_pre=rivdph_pre
      dflw_imp=sqrt(dflw*dflw_pre)
      if(dflw_imp<=0.)dflw_imp=dflw
      !if(dflw_imp>1.e-5.and.darea>1.e-5)then
      if(dflw_imp>0..and.darea>0.)then
        dout_pre=rivout_pre/rivwth
        rivout=rivwth*(dout_pre+grv*dt*dflw_imp*dslope)/(1.+grv*dt*manval**2.*abs(dout_pre)*dflw_imp**(-7./3))
        if(rivout>=0.)then 
          rivout=min(rivout,rivsto/dt)
        else
          rivout=-min(abs(rivout),rivsto_down/dt)
        endif
        rivvel=rivout/darea
      else
        rivout=0.
        rivvel=0.
      endif
      rivout_pre=rivout
      rivdph_pre=rivdph
      
    else
     write(LIS_logunit,*)"[ERR] Wrong outlet id"
     call LIS_endrun()
    endif
    
  end subroutine HYMAP2_calc_rivout_iner
  ! ================================================
  ! ================================================      
  subroutine HYMAP2_no_floodplain_flow(fldelv1,flddph1,fldout,fldvel,sfcelv0,fldout_pre,flddph_pre)
    ! ================================================
    ! Option for floodplain acting as a simple reservoir with no flow  
    ! Augusto Getirana
    ! 29 June 2016
    ! ================================================
    implicit none
    real,    intent(in)  :: fldelv1         !floodplain bed elevation [m]
    real,    intent(in)  :: flddph1         !floodplain depth [m]
    real,    intent(out) :: fldout          !floodplain outflow  [m3/s]
    real,    intent(out) :: fldvel          !floodplain flow velocity [m/s]
    real,    intent(out) :: sfcelv0         !water surface elevation
    real,    intent(out) :: fldout_pre      !previous outflow [m3/s]
    real,    intent(out) :: flddph_pre      !previous floodplain depth [m]

    fldout=0.
    fldvel=0.
    sfcelv0=fldelv1+flddph1
    fldout_pre=0.
    flddph_pre=flddph1

  end subroutine HYMAP2_no_floodplain_flow
  ! ================================================
  ! ================================================      
  subroutine HYMAP2_calc_stonxt(dt,rivout,fldout,rivinf,fldinf,rivsto,fldsto)
    ! ================================================
    ! to   Calculate the storage in the next time step in FTCS diff. eq.
    ! by   Augusto GETIRANA after Dai YAMAZAKI
    ! on   10 Feb 2011
    ! at   LEGOS/CNES, Toulouse, France
    ! Adapted for flow routing implementation in LIS 15 Nov 2011
    ! Adapted for single grid cell computation in 25 Feb 2016
    ! ================================================
    implicit none
    real,    intent(in)  ::  dt          !time step length [sec]
    real,    intent(in)  ::  rivout      ! river outflow          [m3/s]
    real,    intent(in)  ::  fldout      ! floodplain outflow    [m3/s]
    real,    intent(in)  ::  rivinf      !river inflow         [m3/s]
    real,    intent(in)  ::  fldinf      !floodplain inflow    [m3/s]
    real,    intent(out) ::  rivsto      !river storage of current time step       [m3]
    real,    intent(out) ::  fldsto      !flood plain storage                      [m3]
    ! ================================================

    rivsto=max(0.,rivsto+(rivinf-rivout)*dt)
    fldsto=max(0.,fldsto+(fldinf-fldout)*dt)

  end subroutine HYMAP2_calc_stonxt  
  ! ================================================
  ! ================================================      
  subroutine HYMAP2_dynstp(dtin,nxtdst,rivdph,grv,zcadp,dt,nt)
    ! ================================================
    ! Augusto Getirana
    ! 31 Mar 2014
    ! Adapted for flow routing implementation in LIS on 5 Aug 2015
    ! Adapted for single grid cell computation in 22 Feb 2016
    ! ================================================
    implicit none
 
    real,    intent(in)  :: dtin   !input time step length [sec]
    real,    intent(in)  :: nxtdst !distance to next grid [m]
    real,    intent(in)  :: rivdph !depth [m]
    real,    intent(in)  :: grv    !gravity accerelation [m/s2]
    real,    intent(in)  :: zcadp  !alfa coefficient for adaptative time step as described in Bates et al., (2010) [-]
    real,    intent(out) :: dt     !adaptative time step length [sec]
    integer, intent(out) :: nt     !number of time steps
    real                 :: dt0, zdph
    ! ==========
    zdph=max(rivdph,0.01)
    dt0=zcadp*nxtdst/sqrt(grv*zdph)
    nt=int(dtin/dt0-0.001)+1
    dt=dtin/real(nt)

  end subroutine HYMAP2_dynstp
  !=============================================
  !=============================================  
  subroutine HYMAP2_date2frac(yr,mo,da,hr,mn,ss,dt1,time)  
    implicit none
    
    integer, intent(in)     :: yr,mo,da,hr,mn,ss
    real,    intent(in)  :: dt1
    real*8,    intent(out) :: time
    integer, parameter   :: nmday(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)  
    integer              :: day,nyday
    real                 :: dt2
    
    day=sum(nmday(1:mo-1))+da-1; if(mod(yr,4)==0.and.mo>2)day=day+1
    nyday=sum(nmday); if(mod(yr,4)==0)nyday=nyday+1
    
    !ag - 21Sep2017
    time=real(yr)+(real(day)+(real(hr)+(real(mn)+(real(ss)+dt1)/60.)/60.)/24.)/real(nyday)
!    time=yr*10000.+mo*100.+da + (hr+(mn+(ss+dt1)/60.)/60.)/24.

  end subroutine HYMAP2_date2frac
  ! ================================================
  ! ================================================     
  subroutine HYMAP2_calc_evap_fld(dt,grarea,fldare,rivare,fldsto,rivsto,&
                           evpdif,evpout)
  ! ================================================
  ! to   subtract open water efaporation from floodplain water storage
  ! by   Augusto GETIRANA
  ! on   07 Feb 2011 
  ! at   LEGOS/CNES, Toulouse, France
  ! Adapted for single grid cell computation in 19 Feb 2016
  ! ================================================
    implicit none

    real,    intent(in)    ::  dt
    real,    intent(in)    ::  grarea   !area of the grid [m2] 
    real,    intent(in)    ::  fldare   !floodplain surface area [m2]
    real,    intent(in)    ::  rivare   !river surface  area     [m2]
    real,    intent(in)    ::  evpdif   !max. evaporation from floodplains [kg m-2 s-1]
    real,    intent(inout) ::  fldsto   !floodplain storage   [m3]
    real,    intent(inout) ::  rivsto   !river storage        [m3]
    real,    intent(out)   ::  evpout   !effective evaporation from floodplains [kg m-2 s-1]
    real                   ::  devapriv,devapfld !evaporation from river and floodplains [m3  s-1]
    real                   ::  dfldare           !floodplain surface area [m2]
  ! ================================================
      dfldare=max(0.0,fldare-rivare)
      devapriv=min(rivsto,rivare*evpdif*dt/1000.0)!/dt
      devapfld=min(fldsto,dfldare*evpdif*dt/1000.0)!/dt
      rivsto = max(0.0, rivsto-devapriv)
      fldsto = max(0.0, fldsto-devapfld) 
      evpout=(devapriv+devapfld)*1000./grarea/dt

  end subroutine HYMAP2_calc_evap_fld 
  !=============================================
  !=============================================  
  subroutine HYMAP2_get_volume_profile(nz,elevtn,fldhgt,fldstomax,grarea,rivstomax,rivelv,rivlen,rivwth,elv,vol)    
    use LIS_logMod
    implicit none
   
    integer, intent(in)  :: nz
    real*8,     intent(in)  :: elevtn,fldhgt(nz),elv,rivelv !elevation/height are converted to integers [mm]
    real*8,    intent(in)  :: fldstomax(nz),grarea,rivstomax,rivlen,rivwth
    real*8,    intent(out) :: vol
    integer :: i
    real*8    :: h1,h2,v1,v2
    real*8    :: dph(nz)
    real*8    :: dphtmp
    real*8 :: vol1,dh1,dh2,dv,dh
    
    dph(:)=elevtn+fldhgt(:)   
    if(elv>elevtn.and.grarea>rivlen*rivwth)then
      i=1
      do while(elv>dph(i))
        i=i+1
        if(i>nz)then
          vol=fldstomax(nz)+(elv-dph(nz))*grarea
          goto 1
        endif
      enddo
      if(i>1)then
        h1=dph(i-1);v1=fldstomax(i-1)
        h2=dph(i);v2=fldstomax(i)
      elseif(i==1)then
        h1=elevtn;v1=rivstomax
        h2=dph(1);v2=fldstomax(1)
      else
         call LIS_endrun()
      endif
      vol=v1+(v2-v1)*(elv-h1)/(h2-h1)
    else
      vol=max(0.,real(elv-rivelv))*rivlen*rivwth
    endif
1   continue

  end subroutine HYMAP2_get_volume_profile 
  !=============================================
  !=============================================  
  subroutine HYMAP2_get_elevation_profile(nz,elevtn,fldhgt,fldstomax,grarea,rivstomax,rivelv,rivlen,rivwth,elv,vol)

    use LIS_logMod
    implicit none
   
    integer, intent(in)  :: nz
    real*8, intent(in)  :: elevtn,rivelv,fldhgt(nz)
    real*8,    intent(in)  :: fldstomax(nz),grarea,rivstomax,rivlen,rivwth
    real*8,    intent(in)  :: vol
    real*8, intent(out) :: elv
    integer              :: i
    real*8               :: h1,h2,v1,v2
    real*8               :: dph(nz)
    real*8               :: dphtmp
    
    dph(:)=elevtn+fldhgt(:)
    
    if(vol>rivstomax.and.grarea>rivlen*rivwth)then
      i=1
      do while(vol>fldstomax(i))
        i=i+1
        if(i>nz)then
          elv=dph(nz)+(vol-fldstomax(nz))/(grarea-rivlen*rivwth)
          goto 1
        endif
      enddo
      if(i>1)then
        h1=dph(i-1);v1=fldstomax(i-1)
        h2=dph(i);v2=fldstomax(i)
      elseif(i==1)then
        h1=elevtn;v1=rivstomax
        h2=dph(1);v2=fldstomax(1)
      else
        call LIS_endrun()
      endif
      elv=h1+(h2-h1)*(vol-v1)/(v2-v1)
   
    else
      elv=rivelv+(vol/rivlen/rivwth)
    endif
1   continue
  end subroutine HYMAP2_get_elevation_profile
  
end module HYMAP2_modelMod
