!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !REVISION HISTORY: 
! 19 Jan 2016: Augusto Getirana, Inclusion of the local inertia formulation 
!                                and adaptive time step. 
! 13 Apr 2016: Augusto Getirana, Inclusion of option for hybrid runs with a 
!                                river flow map. 
!
subroutine HYMAP2_model_core(n,mis,nseqall,nz,time,dt,  &
     flowmap,linres,evapflag,resopflag,floodflag,dwiflag, &
     rivout_pre,rivdph_pre,grv,                         &
     fldout_pre,flddph_pre,fldelv1,                     &
     outlet,next,elevtn,nxtdst,grarea,                  &
     fldgrd,fldman,fldhgt,fldstomax,                    &
     rivman,rivelv,rivstomax,                           &
     rivlen,rivwth,rivhgt,rivare,slpmin,                &
     trnoff,tbsflw,cntime,                              &
     runoff0,basflw0,evpdif,rnfdwi_ratio,bsfdwi_ratio,&
     rivsto,rivdph,rivvel,rivout,evpout,                &
     fldout,fldsto,flddph,fldvel,fldfrc,                &
     fldare,sfcelv,roffsto,basfsto,rnfdwi,bsfdwi,surfws)
  ! ================================================
  !HyMAP core 
  !Augusto Getirana 
  !6th Dec 2010
  !LEGOS/CNES, Toulouse, France
  !Adapted for flow routing implementation in LIS 09 Nov 2011
  ! ================================================
  
  use HYMAP2_modelMod
  use HYMAP2_resopMod
  use LIS_logMod,     only : LIS_logunit
  
  implicit none

  real*8,  intent(in)  :: time         !current time [year]
  real,    intent(in)  :: dt           !time step length [s]
  integer, intent(in)  :: n
  integer, intent(in)  :: nseqall               !length of 1D sequnece for river and mouth  
  integer, intent(in)  :: nz           !number of stages in the sub-grid discretization

  integer, intent(in)  :: linres       !linear reservoir flag: 1 - use linear reservoirs; or 0 - do not use
  integer, intent(in)  :: evapflag     !evaporation flag: 1 - do not compute ; or 2 - compute evapotation in floodplains
  integer, intent(in)  :: resopflag    !reservoir operation flag: 1 - simulate reservoirs (it requires additional data)
  integer, intent(in)  :: floodflag    !floodplain dynamics flag: 1 - simulate floodplain dynamics; or 2 - do not simulate floodplain dynamics
  integer, intent(in)  :: dwiflag      !deep water infiltration flag
  real,    intent(in)  :: mis         !real undefined value
  real,    intent(in)  :: grv
  real,    intent(in)  :: slpmin                !minimum river slope [-]

  integer, intent(in)  :: outlet(nseqall)       !outlet flag: 0 - river; 1 - ocean
  integer, intent(in)  :: next(nseqall)         !point downstream horizontal
  real,    intent(in)  :: elevtn(nseqall)       !river bed elevation [m]
  real,    intent(in)  :: nxtdst(nseqall)       !distance to the next grid [m]
  real,    intent(in)  :: grarea(nseqall)       !area of the grid [m^2] 
  real,    intent(in)  :: fldgrd(nseqall,nz)    !floodplain gradient [-]
  real,    intent(in)  :: fldman(nseqall)       !maning coefficient for floodplains [-]
  real,    intent(in)  :: fldhgt(nseqall,nz)    !floodplain height [m]
  real,    intent(in)  :: fldstomax(nseqall,nz) !maximum floodplain storage [m3]
  real,    intent(in)  :: rivman(nseqall)       !maning coefficient for rivers [-]
  real,    intent(in)  :: rivelv(nseqall)       !elevation of river bed [m]
  real,    intent(in)  :: rivstomax(nseqall)    !maximum river storage [m3]
  real,    intent(in)  :: rivlen(nseqall)       !river length [m] 
  real,    intent(in)  :: rivwth(nseqall)       !river wiidth [m]
  real,    intent(in)  :: rivhgt(nseqall)       !river heihgt [m]
  real,    intent(in)  :: rivare(nseqall)       !river surface area [m2]
  real,    intent(in)  :: trnoff(nseqall)       !runoff   concentration time parameter [day]
  real,    intent(in)  :: tbsflw(nseqall)       !baseflow concentration time parameter [day]
  real,    intent(in)  :: cntime(nseqall)       !concentration time [s]
  real,    intent(in)  :: runoff0(nseqall)      !input runoff [mm.idt-1]
  real,    intent(in)  :: basflw0(nseqall)      !input baseflow [mm.idt-1]
  real,    intent(in)  :: flowmap(nseqall)      !river flow type map: 1 - kinematic; 2 - diffusive; 3 - inertial
  real,    intent(in)  :: evpdif(nseqall)
  real,    intent(in)  :: bsfdwi_ratio(nseqall) !deep water infiltration ratio from baseflow [-]
  real,    intent(in)  :: rnfdwi_ratio(nseqall) !deep water infiltration ratio from surface runoff [-]

  real,    intent(out) :: roffsto(nseqall)      !runoff   reservoir storage [m3]
  real,    intent(out) :: basfsto(nseqall)      !baseflow reservoir storage [m3]
  real,    intent(out) :: rivsto(nseqall)       !floodplain storage [m3]
  real,    intent(out) :: fldsto(nseqall)       !floodplain storage   [m3]
  real,    intent(out) :: surfws(nseqall)       !surface water storage   [mm]

  real,    intent(out) :: rivdph(nseqall)       !river depth [m]
  real,    intent(out) :: rivvel(nseqall)       !river flow velocity [m/s]
  real,    intent(out) :: rivout(nseqall)       !average river outflow  [m3/s]
  real,    intent(out) :: evpout(nseqall)       !effective evaporation from floodplains [m3]
  real,    intent(out) :: fldout(nseqall)       !floodplain outflow  [m3/s]
  real,    intent(out) :: flddph(nseqall)       !floodplain depth [m]
  real,    intent(out) :: fldvel(nseqall)       !floodplain flow velocity [m/s]
  real,    intent(out) :: fldfrc(nseqall)       !area fraction [-]
  real,    intent(out) :: fldare(nseqall)       !flooded area [m2]
  real,    intent(out) :: sfcelv(nseqall)       !water surface elevation [m]  
  
  real,    intent(out) :: rivout_pre(nseqall)   !previous river outflow [m3/s]
  real,    intent(out) :: rivdph_pre(nseqall)   !previous river depth [m]
  real,    intent(out) :: fldout_pre(nseqall)   !previous flood outflow [m3/s]
  real,    intent(out) :: flddph_pre(nseqall)   !previous flood depth [m]
  real,    intent(out) :: fldelv1(nseqall)      !floodplain elevation [m]
  real,    intent(out) :: bsfdwi(nseqall)      !deep water infiltration from baseflow [mm.idt-1]
  real,    intent(out) :: rnfdwi(nseqall)      !deep water infiltration from surface runoff [mm.idt-1]

  real                 :: flddph1(nseqall)
  real                 :: fldwth(nseqall)
  real                 :: runoff(nseqall)     !runoff+baseflow after linear reservoirs [m3/s]
  real                 :: runoff1(nseqall)    !input runoff [mm.dt-1]
  real                 :: basflw1(nseqall)    !input baseflow [mm.dt-1]
  
  real                 :: sfcelv0(nseqall)    !averaged floodplain surface elevation [m]
  integer              :: ic,ic_down
  
  real                 :: rivelv_down
  real                 :: rivsto_down
  real                 :: rivdph_down
  real                 :: rivdph_pre_down
  real                 :: fldelv1_down
  real                 :: fldsto_down
  real                 :: flddph1_down
  real                 :: flddph_pre_down
! ================================================
  !23 Nov 2016
  !Deep water infiltration
  do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    if(dwiflag==1)then
      call HYMAP2_dwi(runoff0(ic),basflw0(ic),rnfdwi_ratio(ic),&
                      bsfdwi_ratio(ic),runoff1(ic),basflw1(ic),rnfdwi(ic),bsfdwi(ic))
    else
      runoff1(ic)=runoff0(ic)
      basflw1(ic)=basflw0(ic)
    endif
  enddo
! ================================================
  !convert units [kg m-2 sec-1] to [m3 sec-1]
  where(runoff0/=mis.and.grarea/=mis)
     runoff1=runoff1*grarea/1e3
     basflw1=basflw1*grarea/1e3
  else where
     runoff1=0.
     basflw1=0.
  end where
! ================================================
  !Store runoff and baseflow in linear reservoirs
  do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    !ic_down=next(ic)
    if(linres==1)then
       !print*,ic,basfsto(ic),basflw1(ic),tbsflw(ic)
       call HYMAP2_linear_reservoir_lis(dt,mis,runoff1(ic),basflw1(ic),trnoff(ic),tbsflw(ic),&
                                 cntime(ic),roffsto(ic),basfsto(ic),runoff(ic))
    elseif(linres==0)then
       call HYMAP2_no_reservoir_lis(dt,mis,runoff1(ic),basflw1(ic),&
                             roffsto(ic),basfsto(ic),runoff(ic))      
    else
       write(LIS_logunit,*)"HYMAP routing model linear reservoir flag: unknown value"
       stop
    endif
  enddo
! ================================================
  do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    !ic_down=next(ic)
    if(evapflag.ne.0)then
      !calculate evaporation from floodplains
      call HYMAP2_calc_evap_fld(dt,grarea(ic),fldare(ic),rivare(ic),fldsto(ic),rivsto(ic),evpdif(ic),evpout(ic))
    else
      !does not calculate evaporation from floodplains
      evpout(ic)=0.
    endif
  enddo
! ================================================
  do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    !ic_down=next(ic)
    !Calculate flood and river storage
    call HYMAP2_calc_fldstg(nz,grarea(ic),rivlen(ic),rivwth(ic),       &
         rivstomax(ic),fldstomax(ic,:),fldgrd(ic,:),            &
         rivsto(ic),fldsto(ic),rivdph(ic),flddph(ic),flddph1(ic),   &
         fldelv1(ic),fldwth(ic),fldfrc(ic),fldare(ic),          &
         rivelv(ic),nxtdst(ic),rivare(ic))
  enddo
! ================================================
  do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
        
    if(outlet(ic)==0)then
      ic_down=next(ic)
      rivelv_down=rivelv(ic_down)
      rivsto_down=rivsto(ic_down)
      rivdph_down=rivdph(ic_down)
      rivdph_pre_down=rivdph_pre(ic_down)
      fldelv1_down=fldelv1(ic_down)
      fldsto_down=fldsto(ic_down)
      flddph1_down=flddph1(ic_down)
      flddph_pre_down=flddph_pre(ic_down)
    elseif(outlet(ic)==1)then
      !temporarily defining outlet's downstream pixel dynamics as:
      rivelv_down=rivelv(ic)
      rivsto_down=rivsto(ic)
      rivdph_down=-rivelv(ic)
      rivdph_pre_down=-rivelv(ic)
      fldelv1_down=fldelv1(ic)
      fldsto_down=fldsto(ic)
      flddph1_down=-fldelv1(ic)
      flddph_pre_down=-fldelv1(ic)
    else
      write(LIS_logunit,*)"Wrong outlet id"
      stop
    endif
    
    if(flowmap(ic)==1)then
      !Calculate river flow based on the kinematic wave equation
      call HYMAP2_calc_rivout_kine(outlet(ic),dt,rivelv(ic),rivelv_down,nxtdst(ic),&
                            rivwth(ic),sfcelv(ic),rivlen(ic),rivman(ic),slpmin, &
                            rivsto(ic),rivdph(ic),rivout(ic),rivvel(ic))
      !Calculate floodplain
      if(floodflag==1)then
        !Calculate floodplain flow based on the kinematic wave equation
        call HYMAP2_calc_rivout_kine(outlet(ic),dt,fldelv1(ic),fldelv1_down,nxtdst(ic),&
                              fldwth(ic),sfcelv(ic),rivlen(ic),fldman(ic),slpmin,&
                              fldsto(ic),flddph1(ic),fldout(ic),fldvel(ic))
      elseif(floodflag==0)then
        !No floodplain dynamics
        call HYMAP2_no_floodplain_flow(fldelv1(ic),flddph1(ic),fldout(ic),fldvel(ic),sfcelv0(ic),fldout_pre(ic),flddph_pre(ic))
      else
        write(LIS_logunit,*)"HYMAP floodplain dynamics: unknown value"
      endif
    elseif(flowmap(ic)==2.or.flowmap(ic)==3.or.flowmap(ic)==4.or.flowmap(ic)==5)then
      !Calculate river flow based on the local inertia wave equation
      call HYMAP2_calc_rivout_iner(outlet(ic),dt,rivelv(ic),rivelv_down,elevtn(ic),nxtdst(ic),    &
                            rivwth(ic),rivsto(ic),rivsto_down,rivdph(ic),rivdph_down,rivlen(ic),rivman(ic),&
                            grv,rivout(ic),rivvel(ic),sfcelv(ic), &
                            rivout_pre(ic),rivdph_pre(ic),rivdph_pre_down)
      !Calculate floodplain
      if(floodflag==1)then
        !Calculate floodplain flow based on the local inertia wave equation
        call HYMAP2_calc_rivout_iner(outlet(ic),dt,fldelv1(ic),fldelv1_down,elevtn(ic),nxtdst(ic),    &
                              fldwth(ic),fldsto(ic),fldsto_down,flddph1(ic),flddph1_down,rivlen(ic),fldman(ic),&
                              grv,fldout(ic),fldvel(ic),sfcelv0(ic),  &
                              fldout_pre(ic),flddph_pre(ic),flddph_pre_down)
      elseif(floodflag==0)then
        !No floodplain dynamics
        call HYMAP2_no_floodplain_flow(fldelv1(ic),flddph1(ic),fldout(ic),fldvel(ic),sfcelv0(ic),fldout_pre(ic),flddph_pre(ic))
      else
        write(LIS_logunit,*)"HYMAP floodplain dynamics: unknown value"
      endif
    else
       write(LIS_logunit,*)"HYMAP routing method: unknown value",ic,flowmap(ic)
       stop
    endif
  enddo
! ================================================
    !reservoir simulation
    if(resopflag==1)call HYMAP2_resop_run(n,mis,nseqall,nz,time,dt,next,rivsto,fldsto,rivout,fldout,rivvel,fldvel,&
                     elevtn,fldhgt,fldstomax,grarea,rivstomax,rivelv,rivlen,rivwth)
! ================================================
  do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    if(outlet(ic)==0)then
      ic_down=next(ic)
      rivsto_down=rivsto(ic_down)
      fldsto_down=fldsto(ic_down)
    elseif(outlet(ic)==1)then
      rivsto_down=1e20 !rivsto(ic)
      fldsto_down=1e20 !fldsto(ic)
    else
      write(LIS_logunit,*)"Wrong outlet id"
      stop
    endif
    !Update water storage in river and floodplain reservoirs
    call HYMAP2_calc_stonxt(outlet(ic),dt,rivout(ic),fldout(ic),rivsto(ic),fldsto(ic),rivsto_down,fldsto_down)
    if(outlet(ic)==0)then
      rivsto(ic_down)=rivsto_down
      fldsto(ic_down)=fldsto_down
    endif
  enddo
! ================================================
! Calculate runoff in the river network for the current time step
  do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    !ic_down=next(ic)
    call HYMAP2_calc_runoff(dt,fldfrc(ic),runoff(ic),rivsto(ic),fldsto(ic))
  enddo
! ================================================
! Calculate surface water storage [mm]
  do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    surfws(ic)=1e3*(rivsto(ic)+fldsto(ic))/grarea(ic)
  enddo  
  
end subroutine HYMAP2_model_core
