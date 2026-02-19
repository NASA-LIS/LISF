!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! !MODULE: HYMAP3_modelMod
!
! !DESCRIPTION
!
! !REVISION HISTORY:
! 19 Jan 2016: Augusto Getirana, Inclusion of the local inertia
!   formulation and adaptive time step.
! 02 Mar 2020: Augusto Getirana, Fixed river/floodplain storage
! 27 Apr 2020: Augusto Getirana,  Added support for urban drainage
!
module HYMAP3_modelMod

  use LIS_logMod,     only : LIS_logunit

  implicit none
  private

  public :: HYMAP3_dwi
  public :: HYMAP3_linear_reservoir_lis
  public :: HYMAP3_no_reservoir_lis
  public :: HYMAP3_calc_runoff
  public :: HYMAP3_calc_fldstg
  public :: HYMAP3_calc_rivout_kine
  public :: HYMAP3_calc_rivout_iner
  public :: HYMAP3_calc_stonxt
  public :: HYMAP3_dynstp
  public :: HYMAP3_date2frac
  public :: HYMAP3_calc_evap_fld
  public :: HYMAP3_get_elevation_profile
  public :: HYMAP3_get_volume_profile
  public :: HYMAP3_merge_river_floodplain
  public :: HYMAP3_no_floodplain_flow
  public :: HYMAP3_calc_fldstg_lev

contains
  ! ================================================
  ! ================================================
  subroutine HYMAP3_dwi(runoff0, basflw0, grarea, roffdwi_ratio, &
       basfdwi_ratio, roffdwi, basfdwi, runoff)

    implicit none

    real,    intent(inout) :: runoff0       !input runoff [mm.idt-1]
    real,    intent(inout) :: basflw0       !input baseflow [mm.idt-1]
    real,    intent(in)    :: grarea        !grid area [m2]
    real,    intent(in)    :: basfdwi_ratio !deep water infiltration ratio from baseflow [-]
    real,    intent(in)    :: roffdwi_ratio !deep water infiltration ratio from surface runoff [-]
    real,    intent(out)   :: basfdwi       !deep water infiltration from baseflow [mm.idt-1]
    real,    intent(out)   :: roffdwi       !deep water infiltration from surface runoff [mm.idt-1]
    real,    intent(out)   :: runoff        !remaining total runoff floding to the river network [m3.idt-1]

    roffdwi=runoff0*roffdwi_ratio
    runoff0=runoff0*(1.-max(0.,min(roffdwi_ratio,1.)))

    basfdwi=basflw0*basfdwi_ratio
    basflw0=basflw0*(1.-max(0.,min(basfdwi_ratio,1.)))

    runoff=(runoff0+basflw0)*grarea/1e3

  end subroutine HYMAP3_dwi
  ! ================================================
  ! ================================================
  subroutine HYMAP3_linear_reservoir_lis(dt, rmis, runoff0, basflw0, &
       grarea, trnoff, tbsflw, cntime, roffsto,&
       basfsto, runoff)
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
    real,    intent(in)    :: runoff0 !runoff   before linear reservoir [mm/s]
    real,    intent(in)    :: basflw0 !baseflow before linear reservoir [mm/s]
    real,    intent(in)    :: grarea  !grid area [m2]
    real,    intent(in)    :: trnoff  !runoff   concentration time parameter [day]
    real,    intent(in)    :: tbsflw  !baseflow concentration time parameter [day]
    real,    intent(in)    :: cntime  !concentration time [s]
    real,    intent(inout) :: roffsto !runoff   reservoir storage [m3]
    real,    intent(inout) :: basfsto !baseflow reservoir storage [m3]
    real,    intent(out)   :: runoff  !runoff+baseflow after  linear reservoir [m3/s]

    real                   :: vroff, vbflw
    real                   :: qrun, qbas
    ! ================================================
    vroff=0.0
    vbflw=0.0
    qrun=0.0
    qbas=0.0
    runoff=0.0
    !runoff
    !update reservoir volumes with input data (runoff)
    roffsto = roffsto + runoff0*dt*grarea/1e3
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
    basfsto = basfsto + basflw0*dt*grarea/1e3
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

  end subroutine HYMAP3_linear_reservoir_lis
  ! ================================================
  ! ================================================
  subroutine HYMAP3_no_reservoir_lis(dt, rmis, runoff0, basflw0, &
       roffsto, basfsto, runoff)
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

  end subroutine HYMAP3_no_reservoir_lis
  ! ================================================
  ! ================================================
  subroutine HYMAP3_calc_runoff(dt, fldfrc, runoff, rivsto, fldsto)
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

  end subroutine HYMAP3_calc_runoff
  ! ================================================
  !ag (2Mar2020)
  ! ================================================
  subroutine HYMAP3_calc_fldstg(nz, grarea, rivlen, rivwth,       &
                         rivstomax, fldstomax, fldhgt,            &
                         rivsto, fldsto, rivdph, flddph, flddph1, &
                         fldelv1, fldwth, fldfrc, fldare,         &
                         rivelv, rivare, levhgt)
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
    real,    intent(in)    :: levhgt       !levee height [m]
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
          write(LIS_logunit,*) '[ERR] Please check [HYMAP3_calc_fldstg] '
          call LIS_endrun()
       endif
       flddph=h1+(h2-h1)*(rstoall-v1)/(v2-v1)
       fldfrc=min(max((rivare/grarea)+real(i-1)/real(nz)+ &
            (1./real(nz))*(rstoall-v1)/(v2-v1),0.),1.)
1      continue
       !ag(23Feb2023) adapt line below to account for levees
       rivsto=rivstomax+rivlen*rivwth*max(0.,flddph-levhgt)
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

  end subroutine HYMAP3_calc_fldstg
  ! ================================================
  ! ================================================
  subroutine HYMAP3_calc_rivout_kine(outlet, dt, rivelv, rivelv_down, &
       nxtdst, &
       rivwth, sfcelv, rivlen, manval, slpmin, &
       rivsto, rivdph, rivout, rivvel)

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

    real                 :: rslope,rarea,rvel

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
       write(LIS_logunit,*)"[ERR] Wrong outlet id", outlet
       call LIS_endrun()
    endif

  end subroutine HYMAP3_calc_rivout_kine
  ! ================================================
  ! ================================================
  subroutine HYMAP3_calc_rivout_iner(time, ic, flowid, slpmax,           &
       bckslpmax,                                                        &
       refout, outlet, dt, rivelv, rivelv_down, elevtn, nxtdst,          &
       rivwth, rivsto, rivsto_down, rivdph, rivdph_down, rivlen, manval, &
       grv, rivout, rivvel, sfcelv,                                      &
       rivout_pre, rivdph_pre, rivdph_pre_down)

    use LIS_logMod

    ! ================================================
    ! Calculate discharge, local inertia
    ! Augusto Getirana
    ! 16 Apr 2014
    ! Adapted for flow routing implementation in LIS on 5 Aug 2015
    ! Adapted for single grid cell computation in 19 Feb 2016
    ! ================================================

    implicit none

    real*8,  intent(in)  :: time         !current time [year]
    integer, intent(in)  :: ic
    integer, intent(in)  :: flowid          !identifies river (1) and floodplain (2) dynamics
    real,    intent(in)  :: slpmax          !maximum slope allowed [m/m]
    real,    intent(in)  :: bckslpmax       !maximum backwater slope allowed [m/m]
    real,    intent(in)  :: refout          !reference outflow - only applied to floodplain
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
    real                 :: avgslp, darea
    real                 :: dflw,dout_pre,dflw_pre,dflw_imp

    ! ================================================
    sfcelv=rivelv+rivdph
    sfcelv_down=rivelv_down+rivdph_down
    sfcelv_pre=rivelv+rivdph_pre
    sfcelv_pre_down=rivelv_down+rivdph_pre_down

    if(outlet==0)then
       avgslp=min(max((sfcelv-sfcelv_down)/nxtdst,bckslpmax),slpmax)
       dflw=max(0.,max(sfcelv,sfcelv_down)-rivelv)
       darea=rivwth*dflw
       dflw_pre=max(0.,max(sfcelv_pre,sfcelv_pre_down)-rivelv)
       dflw_imp=sqrt(dflw*dflw_pre)
       if(dflw_imp<=0.)dflw_imp=dflw

       if(dflw_imp>1e-10.and.darea>1e-10)then
          dout_pre=rivout_pre/rivwth
          if(dout_pre<0)then
             dout_pre=max(rivout_pre,-4e5)/rivwth
          else
             dout_pre=min(rivout_pre,4e5)/rivwth
          endif
          rivout=rivwth*(dout_pre+grv*dt*dflw_imp*avgslp)/ &
               (1.+grv*dt*manval**2.*abs(dout_pre)*dflw_imp**(-7./3))

          if(flowid==2.and.((rivout<0.and.refout>0).or. &
               (rivout>0.and.refout<0)))rivout=0.
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
       avgslp=min(max((sfcelv-sfcelv_down)/nxtdst,-1.e-3),1.e-3)
       !ag(7Apr2021)
       dflw=max(0.,max(sfcelv,sfcelv_down)-rivelv)
       darea=rivwth*dflw
       dflw_pre=max(0.,max(sfcelv_pre,sfcelv_pre_down)-rivelv)
       dflw_imp=sqrt(dflw*dflw_pre)
       if(dflw_imp<=0.)dflw_imp=dflw

       if(dflw_imp>1e-10.and.darea>1e-10)then
          dout_pre=rivout_pre/rivwth
          if(dout_pre<0)then
             dout_pre=max(rivout_pre,-4e5)/rivwth
          else
             dout_pre=min(rivout_pre,4e5)/rivwth
          endif
          rivout=rivwth*(dout_pre+grv*dt*dflw_imp*avgslp)/ &
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

    else
       write(LIS_logunit,*)"[ERR] Wrong outlet id"
       call LIS_endrun()
    endif

  end subroutine HYMAP3_calc_rivout_iner
  ! ================================================
  ! ================================================
  subroutine HYMAP3_calc_rivfldout_iner(ic, time, outlet, dt, elevtn, &
       nxtdst,                                                   &
       rivlen, grv,                                              &
       rivelv, rivelv_down, rivwth, rivsto, rivsto_down, rivdph, &
       rivdph_down,                                              &
       rivdph_pre_down, rivman,                                  &
       fldelv, fldelv_down, fldwth, fldsto, fldsto_down, flddph, &
       flddph_down,                                              &
       flddph_pre_down, fldman,                                  &
       rivout, rivvel, rivsfcelv, rivout_pre, rivdph_pre,        &
       fldout, fldvel, fldsfcelv, fldout_pre, flddph_pre)

    use LIS_logMod
    ! ================================================
    ! Calculate weighted average discharge, local inertia
    ! Augusto Getirana
    ! 19 Apr 2024
    ! ================================================
    implicit none

    integer,  intent(in)   :: ic
    real*8,   intent(in)   :: time
    integer, intent(in)    :: outlet          !outlet flag: 0 - river; 1 - ocean
    real,    intent(in)    :: dt              !time step length [sec]
    real,    intent(in)    :: elevtn          !surface elevation [m]
    real,    intent(in)    :: nxtdst          !distance to next grid [m]
    real,    intent(in)    :: rivlen          !river/floodplain length [m]
    real,    intent(in)    :: grv             !gravity [m/s2]
    real,    intent(in)    :: rivelv          !river/floodplain bed elevation [m]
    real,    intent(in)    :: rivelv_down     !downstream river/floodplain bed elevation [m]
    real,    intent(in)    :: rivwth          !river/floodplain width [m]
    real,    intent(in)    :: rivsto          !river/floodplain storage [m3]
    real,    intent(in)    :: rivsto_down     !downstream river/floodplain storage [m3]
    real,    intent(in)    :: rivdph          !river/floodplain depth [m]
    real,    intent(in)    :: rivdph_down     !dowstream river/floodplain depth [m]
    real,    intent(in)    :: rivdph_pre_down !previous downstream river depth [m]
    real,    intent(in)    :: rivman          !manning coefficient for river/floodplain
    real,    intent(in)    :: fldelv          !river/floodplain bed elevation [m]
    real,    intent(in)    :: fldelv_down     !downstream river/floodplain bed elevation [m]
    real,    intent(in)    :: fldwth          !river/floodplain width [m]
    real,    intent(in)    :: fldsto          !river/floodplain storage [m3]
    real,    intent(in)    :: fldsto_down     !downstream river/floodplain storage [m3]
    real,    intent(in)    :: flddph          !river/floodplain depth [m]
    real,    intent(in)    :: flddph_down     !dowstream river/floodplain depth [m]
    real,    intent(in)    :: flddph_pre_down !previous downstream river depth [m]
    real,    intent(in)    :: fldman          !manning coefficient for river/floodplain
    real,    intent(out)   :: rivout          !river/floodplain outflow  [m3/s]
    real,    intent(out)   :: rivvel          !river/floodplain flow velocity [m/s]
    real,    intent(out)   :: rivsfcelv       !river water surface elevation
    real,    intent(inout) :: rivout_pre      !previous outflow [m3/s]
    real,    intent(out)   :: rivdph_pre      !previous river depth [m]
    real,    intent(out)   :: fldout          !river/floodplain outflow  [m3/s]
    real,    intent(out)   :: fldvel          !river/floodplain flow velocity [m/s]
    real,    intent(out)   :: fldsfcelv       !floodplain water surface elevation
    real,    intent(inout) :: fldout_pre      !previous outflow [m3/s]
    real,    intent(out)   :: flddph_pre      !previous river depth [m]

    real                   :: rivsfcelv_down  !downstream water surface elevation [m]
    real                   :: fldsfcelv_down  !downstream water surface elevation [m]
    real                   :: avgelv          !weighted average river/floodplain bed elevation [m]
    real                   :: avgelv_down     !weighted average downstream river/floodplain bed elevation [m]
    real                   :: avgwth          !weighted average river/floodplain width [m]
    real                   :: avgdph          !weighted average river/floodplain depth [m]
    real                   :: avgdph_down     !weighted average dowstream river/floodplain depth [m]
    real                   :: avgdph_pre_down !weighted average previous downstream river depth [m]
    real                   :: avgman          !weighted average manning coefficient for river/floodplain
    real                   :: totsto          !total water storage in river and floodplain [m3]
    real                   :: totsto_down     !total water storage in downstream river and floodplain [m3]
    real                   :: avgout          !weighted average river/floodplain outflow  [m3/s]
    real                   :: avgvel          !weighted average river/floodplain flow velocity [m/s]
    real                   :: avgsfcelv       !weighted average river water surface elevation
    real                   :: avgout_pre      !weighted average previous outflow [m3/s]
    real                   :: avgdph_pre      !weighted average previous river depth [m]
    real                   :: avgsfcelv_down     !weighted average downstream water surface elevation
    real                   :: avgsfcelv_pre      !weighted average previous step water surface elevation (t-1) [m]
    real                   :: avgsfcelv_pre_down !weighted average previous step downstream water surface elevation (t-1) [m]
    real                   :: rivslp,fldslp,avgslp, darea
    real                   :: dflw,dout_pre,dflw_pre,dflw_imp

    ! ================================================
    totsto             = rivsto + fldsto
    totsto_down        = rivsto_down + fldsto_down

    rivsfcelv          = rivelv + rivdph
    fldsfcelv          = fldelv + flddph

    rivsfcelv_down     = rivelv_down + rivdph_down
    fldsfcelv_down     = fldelv_down + flddph_down

    rivslp             = (rivsfcelv - rivsfcelv_down) / nxtdst
    fldslp             = (fldsfcelv - fldsfcelv_down) / nxtdst

    if(totsto>0)then
       avgelv             = rivelv
       avgwth             = (rivwth*rivsto + fldwth*fldsto) / totsto
       avgdph             = (rivdph*rivsto + flddph*fldsto) / totsto
       avgman             = (rivman*rivsto + fldman*fldsto) / totsto
       avgslp             = (rivslp*rivsto + fldslp*fldsto) / totsto
    else
       avgelv             = rivelv
       avgwth             = rivwth
       avgdph             = rivdph
       avgman             = rivman
       avgslp             = rivslp
    endif

    if(totsto_down>0.and.totsto_down<=1e20)then
       avgelv_down        = rivelv_down
       avgdph_down        = (rivdph_down*rivsto_down + &
            flddph_down*fldsto_down) / totsto_down
       avgdph_pre_down    = (rivdph_pre_down*rivsto_down + &
            flddph_pre_down*fldsto_down) / totsto_down
    else
       avgelv_down        = rivelv_down
       avgdph_down        = rivdph_down
       avgdph_pre_down    = rivdph_pre_down
    endif

    avgsfcelv          = avgelv + avgdph
    avgsfcelv_down     = avgelv_down + avgdph_down
    avgsfcelv_pre      = avgelv + avgdph_pre
    avgsfcelv_pre_down = avgelv_down + avgdph_pre_down

    if(outlet==0)then
       dflw=max(avgsfcelv,avgsfcelv_down)-avgelv
       darea=avgwth*dflw
       dflw_pre=max(0.,max(avgsfcelv_pre,avgsfcelv_pre_down)-avgelv)
       dflw_imp=sqrt(dflw*dflw_pre)
       if(dflw_imp<=0.)dflw_imp=dflw

       if(dflw_imp>0..and.darea>0.)then
          dout_pre=avgout_pre/avgwth
          avgout=avgwth*(dout_pre+grv*dt*dflw_imp*avgslp)/&
               (1.+grv*dt*avgman**2.*abs(dout_pre)*dflw_imp**(-7./3))
          if(avgout>=0.)then
             avgout=min(avgout,totsto/dt)
          else
             avgout=-min(abs(avgout),totsto_down/dt)
          endif
          avgvel=avgout/darea
       else
          avgout=0.
          avgvel=0.
       endif
       avgout_pre=avgout
       avgdph_pre=avgdph

    elseif(outlet==1)then
       !prevent numerical instability from steep slopes at the outlet
       if(avgslp<0)then
          avgslp=max(avgslp,-0.001)
       else
          avgslp=min(avgslp,0.001)
       endif
       dflw=max(avgsfcelv,avgsfcelv_down)-avgelv
       darea=avgwth*dflw
       dflw_pre=max(avgsfcelv_pre,avgsfcelv_pre_down)-avgelv
       dflw_imp=sqrt(dflw*dflw_pre)
       if(dflw_imp<=0.)dflw_imp=dflw
       if(dflw_imp>0..and.darea>0.)then
          dout_pre=avgout_pre/avgwth
          avgout=avgwth*(dout_pre+grv*dt*dflw_imp*avgslp)/ &
               (1.+grv*dt*avgman**2.*abs(dout_pre)*dflw_imp**(-7./3))
          if(avgout>=0.)then
             avgout=min(avgout,totsto/dt)
          else
             avgout=-min(abs(avgout),totsto_down/dt)
          endif
          avgvel=avgout/darea
       else
          avgout=0.
          avgvel=0.
       endif
       avgout_pre=avgout
       avgdph_pre=avgdph

    else
       write(LIS_logunit,*)"[ERR] Wrong outlet id"
       call LIS_endrun()
    endif

    !all outflow is attributed to river; floodplain has no flow
    rivout     = avgout
    rivvel     = avgvel
    rivout_pre = rivout
    rivdph_pre = rivdph

    fldout     = 0.
    fldvel     = 0.
    fldout_pre = fldout
    flddph_pre = flddph

  end subroutine HYMAP3_calc_rivfldout_iner
  ! ================================================
  ! ================================================
  subroutine HYMAP3_no_floodplain_flow(fldelv1, flddph1, fldout, fldvel, &
       sfcelv0, fldout_pre, flddph_pre)
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

  end subroutine HYMAP3_no_floodplain_flow
  ! ================================================
  ! ================================================
  subroutine HYMAP3_calc_stonxt(dt, rivout, fldout, rivinf, fldinf, &
       rivsto, fldsto)

    use LIS_coreMod
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
    real,    intent(out) ::  fldsto      !flood plain storage [m3]
    ! ================================================

    rivsto=max(0.,rivsto+(rivinf-rivout)*dt)
    fldsto=max(0.,fldsto+(fldinf-fldout)*dt)

  end subroutine HYMAP3_calc_stonxt
  ! ================================================
  ! ================================================
  subroutine HYMAP3_dynstp(dtin, nxtdst, rivdph, grv, zcadp, dt, nt)
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
    dt0=max(10.,zcadp*nxtdst/sqrt(grv*zdph))
    nt=int(dtin/dt0-0.001)+1
    dt=dtin/real(nt)

  end subroutine HYMAP3_dynstp
  !=============================================
  !=============================================
  subroutine HYMAP3_date2frac(yr, mo, da, hr, mn, ss, dt1, time)

    implicit none

    integer, intent(in)  :: yr,mo,da,hr,mn,ss
    real,    intent(in)  :: dt1
    real*8,  intent(out) :: time

    integer, parameter   :: &
         nmday(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    integer              :: day,nyday

    day=sum(nmday(1:mo-1))+da-1; if(mod(yr,4)==0.and.mo>2)day=day+1
    nyday=sum(nmday); if(mod(yr,4)==0)nyday=nyday+1

    !ag - 21Sep2017
    time=real(yr)+(real(day)+(real(hr)+ &
         (real(mn)+(real(ss)+dt1)/60.)/60.)/24.)/real(nyday)

  end subroutine HYMAP3_date2frac
  ! ================================================
  ! ================================================
  subroutine HYMAP3_calc_evap_fld(dt, grarea, fldare, rivare, fldsto, &
       rivsto, evpdif, evpout)
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
    devapriv=min(rivsto,rivare*evpdif*dt/1000.0)
    devapfld=min(fldsto,dfldare*evpdif*dt/1000.0)
    rivsto = max(0.0, rivsto-devapriv)
    fldsto = max(0.0, fldsto-devapfld)
    evpout=(devapriv+devapfld)*1000./grarea/dt

  end subroutine HYMAP3_calc_evap_fld
  !=============================================
  !=============================================
  subroutine HYMAP3_get_volume_profile(nz, elevtn, fldhgt, fldstomax, &
       grarea, rivstomax, rivelv, rivlen, rivwth, elv, vol)

    use LIS_logMod

    implicit none

    integer, intent(in)  :: nz
    real*8,  intent(in)  :: elevtn,fldhgt(nz),elv,rivelv !elevation/height are converted to integers [mm]
    real*8,    intent(in)  :: fldstomax(nz),grarea,rivstomax,rivlen,rivwth
    real*8,    intent(out) :: vol

    integer :: i
    real*8    :: h1,h2,v1,v2
    real*8    :: dph(nz)

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
          write(LIS_logunit,*) &
               '[ERR] Internal error in HYMAP3_get_volume_profile'
          call LIS_endrun()
       endif
       vol=v1+(v2-v1)*(elv-h1)/(h2-h1)
    else
       vol=max(0.,real(elv-rivelv))*rivlen*rivwth
    endif
1   continue

  end subroutine HYMAP3_get_volume_profile
  !=============================================
  !=============================================
  subroutine HYMAP3_get_elevation_profile(nz, elevtn, fldhgt, fldstomax, &
       grarea, rivstomax, rivelv, rivlen, rivwth, elv, vol)

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
          write(LIS_logunit,*) &
               '[ERR] Internal error in HYMAP3_get_elevation_profile'
          call LIS_endrun()
       endif
       elv=h1+(h2-h1)*(vol-v1)/(v2-v1)

    else
       elv=rivelv+(vol/rivlen/rivwth)
    endif
1   continue
  end subroutine HYMAP3_get_elevation_profile

  !=============================================
  subroutine HYMAP3_merge_river_floodplain(rivsto, rivman, rivdph, &
       rivwth,                                                     &
       fldsto, fldman, flddph, rivlen, rivdph_down, rivsto_down,   &
       fldsto_down,                                                &
       tmpsto, tmpman, tmpdph, tmpwth, tmpsto_down, tmpdph_down)

    use LIS_logMod

    implicit none

    real, intent(in)  :: rivsto,rivman,rivdph,rivwth,fldsto,fldman,&
         flddph,rivlen,rivdph_down,rivsto_down,fldsto_down
    real, intent(out) :: tmpsto,tmpman,tmpdph,tmpwth,tmpsto_down,&
         tmpdph_down

    tmpsto=rivsto+fldsto
    tmpdph=rivdph
    if(tmpsto>0)then
       tmpman=(rivman*rivsto+fldman*fldsto)/tmpsto
       tmpwth=tmpsto/tmpdph/rivlen
    else
       tmpman=rivman
       tmpwth=rivwth
    endif
    tmpdph_down=rivdph_down
    tmpsto_down=min(rivsto_down+fldsto_down,1e20)

  end subroutine HYMAP3_merge_river_floodplain
  !=============================================
  !ag (22May2024)
  ! ================================================
  subroutine HYMAP3_calc_fldstg_lev(nz, grarea, rivlen, rivwth, &
       rivstomax,                                               &
       fldstomax, fldhgt, rivare, levhgt, levstomax,            &
       fldonlystomax, fldstoatlev, elevtn,                      &
       rivsto, fldsto, rivdph, flddph, flddph1, fldelv1,        &
       fldwth, fldfrc, fldare)

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
    real,    intent(in)    :: rivare       !river surface area [m2]
    real,    intent(in)    :: levhgt       !levee height [m]
    real,    intent(in)    :: levstomax    !maximum levee water storage [m3]
    real,    intent(in)    :: fldonlystomax(nz)!maximum floodplain-only storage - excludes river storage [m3]
    real,    intent(in)    :: fldstoatlev  !floodplain-only storage at levee elevation - excludes river storage [m3]
    real,    intent(in)    :: elevtn       !grid elevation [m]
    real,    intent(inout) :: rivsto       !river storage [m3]
    real,    intent(inout) :: fldsto       !flood plain storage [m3]
    real,    intent(out)   :: rivdph       !river depth [m]
    real,    intent(out)   :: flddph       !floodplain depth [m]
    real,    intent(out)   :: flddph1      !mean floodplain depth [m]
    real,    intent(out)   :: fldelv1      !mean floodplain elevation [m]
    real,    intent(out)   :: fldwth       !floodplain width [m]
    real,    intent(out)   :: fldfrc       !area fraction
    real,    intent(out)   :: fldare

    integer                :: i
    real*8                 :: rstoall,fill
    real*8                 :: h1,h2,v1,v2
    ! ================================================
    rstoall=max(0.,rivsto+fldsto)
    !case 1: levee exists
    if(levhgt>0.)then

      !case 1.1: rivsto+fldsto<=rivstomax - all water fits in main channel
       if(rstoall<=rivstomax.or.grarea>=rivare)then
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
          fldelv1=elevtn !rivelv+rivdph

       !case 1.2: total water storage doesn't fit in main channel but
       ! fits in space below levee (riv+fld)
       !in this case, there's no exchange between river and floodplains,
       ! except if rivsto<rivstomax (i.e., case 1.2.1)
       elseif(rstoall>rivstomax.and. &
            rstoall<=fldstoatlev+rivstomax+levstomax)then

        !case 1.2.1: river bank isn't full; fill it with fldsto, if available
          if(rivsto<rivstomax)then
             fill=max(0.,min(fldsto,rivstomax-rivsto))
             fldsto=max(0.,fldsto-fill)
             rivsto=max(0.,rivsto+fill)
          endif
          if(rivare>0)then
             rivdph=rivsto/rivare
          else
             rivdph=0.
          endif

          !case 1.2.2: no floodplain water storage
          if(fldsto==0.)then
             flddph=0.
             fldwth=0.
             flddph1=0.
             fldelv1=elevtn
             fldfrc=min(max(rivare/grarea,0.),1.)
             fldare=fldfrc*grarea

          !case 1.2.3: there is floodplain, but level below levee
          ! elevation, i.e., levee isn't flooded
          else
             i=1
             do while(fldsto>fldonlystomax(i))
                i=i+1
                if(i>nz)then
                   flddph=fldhgt(nz)+ &
                        (fldsto-fldonlystomax(nz))/(grarea-rivare)
                   fldfrc=1.
                   goto 1
                endif
             enddo
             if(i>1)then
                h1=fldhgt(i-1);v1=fldonlystomax(i-1)
                h2=fldhgt(i);v2=fldonlystomax(i)
             elseif(i==1)then
                h1=0.;v1=0.
                h2=fldhgt(1);v2=fldonlystomax(1)
             else
                write(LIS_logunit,*) &
                     '[ERR] Please check [HYMAP3_calc_fldstg_lev] '
                call LIS_endrun()
             endif
             flddph=h1+(h2-h1)*(fldsto-v1)/(v2-v1)
             fldfrc=min(max((rivare/grarea)+ &
                  real(i-1)/real(nz)+(1./real(nz))*(fldsto-v1)/(v2-v1),0.),1.)
1            continue
             fldare=grarea*fldfrc
             fldwth=(fldare-rivare)/rivlen !fldare/rivlen
             flddph1=flddph !fldsto/fldwth/rivlen
             fldelv1=elevtn+flddph1 !rivelv+rivdph-flddph1
          endif
       !case 1.3: levee is flooded (solution is similar to case without levees)
       else
          i=1
          do while(rstoall>fldstomax(i))
             i=i+1
             if(i>nz)then
                flddph=fldhgt(nz)+(rstoall-fldstomax(nz))/grarea
                fldfrc=1.
                goto 2
             endif
          enddo
          if(i>1)then
             h1=fldhgt(i-1);v1=fldstomax(i-1)
             h2=fldhgt(i);v2=fldstomax(i)
          elseif(i==1)then
             h1=0.;v1=rivstomax
             h2=fldhgt(1);v2=fldstomax(1)
          else
             write(LIS_logunit,*) '[ERR] Please check [HYMAP3_calc_fldstg] '
             call LIS_endrun()
          endif
          flddph=h1+(h2-h1)*(rstoall-v1)/(v2-v1)
          fldfrc=min(max((rivare/grarea)+ &
               real(i-1)/real(nz)+(1./real(nz))*(rstoall-v1)/(v2-v1),0.),1.)
2         continue
          rivsto=rivstomax+rivlen*rivwth*max(0.,flddph)
          !ag(27Apr2020) for cases where there's no river (e.g., urban areas)
          if(rivwth>0.and.rivlen>0)then
             rivdph=rivsto/rivlen/rivwth
          else
             rivdph=0
          endif
          fldsto=max(rstoall-rivsto,0.)
          fldare=grarea*fldfrc
          fldwth=(fldare-rivare)/rivlen
          flddph1=flddph
          fldelv1=elevtn+flddph1
       endif

    !case 2: levees do not exist
    else
       if(rstoall>rivstomax.and.grarea>rivare)then
          i=1
          do while(rstoall>fldstomax(i))
             i=i+1
             if(i>nz)then
                flddph=fldhgt(nz)+(rstoall-fldstomax(nz))/grarea
                fldfrc=1.
                goto 3
             endif
          enddo
          if(i>1)then
             h1=fldhgt(i-1);v1=fldstomax(i-1)
             h2=fldhgt(i);v2=fldstomax(i)
          elseif(i==1)then
             h1=0.;v1=rivstomax
             h2=fldhgt(1);v2=fldstomax(1)
          else
             write(LIS_logunit,*) '[ERR] Please check [HYMAP3_calc_fldstg] '
             call LIS_endrun()
          endif
          flddph=h1+(h2-h1)*(rstoall-v1)/(v2-v1)
          fldfrc=min(max((rivare/grarea)+ &
               real(i-1)/real(nz)+(1./real(nz))*(rstoall-v1)/(v2-v1),0.),1.)
3         continue
          rivsto=rivstomax+rivlen*rivwth*max(0.,flddph)
          !ag(27Apr2020) for cases where there's no river (e.g., urban areas)
          if(rivwth>0.and.rivlen>0)then
             rivdph=rivsto/rivlen/rivwth
          else
             rivdph=0
          endif
          fldsto=max(rstoall-rivsto,0.)
          fldare=grarea*fldfrc
          fldwth=(fldare-rivare)/rivlen
          flddph1=flddph
          fldelv1=elevtn+flddph1
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
          fldelv1=elevtn
       endif
    endif

  end subroutine HYMAP3_calc_fldstg_lev
  ! ================================================

end module HYMAP3_modelMod
