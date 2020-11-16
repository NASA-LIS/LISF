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
! 13 May 2019: Augusto Getirana, Fixed reservoir operation module and replaced 
!                                routine calc_fldstg with calc_fldstg 
!                                
#include "LIS_misc.h"
subroutine HYMAP2_model_core(n,it,mis,nseqall,nz,time,dt,  &
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
  use LIS_logMod
  use LIS_coreMod
  use HYMAP2_routingMod

  use LIS_coreMod
  use LIS_mpiMod

  implicit none

  real*8,  intent(in)  :: time         !current time [year]
  real,    intent(in)  :: dt           !time step length [s]
  integer, intent(in)  :: n
  integer, intent(in)  :: it           !sub time step
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
  real                 :: inflow              !inflow from upstream grid cells [m3/s]
  
  real                 :: sfcelv0(nseqall)    !averaged floodplain surface elevation [m]
  integer              :: ic,ic_down,icg
  
  real                 :: rivelv_down
  real                 :: rivsto_down
  real                 :: rivdph_down
  real                 :: rivdph_pre_down
  real                 :: fldelv1_down
  real                 :: fldsto_down
  real                 :: flddph1_down
  real                 :: flddph_pre_down

  real                 :: uprivout
  real                 :: upfldout
  real                 :: uprivout_down
  real                 :: upfldout_down
  real                 :: rivout_down
  real                 :: fldout_down
  
  real, allocatable    :: rivelv_glb(:)
  real, allocatable    :: rivsto_glb(:)
  real, allocatable    :: rivdph_glb(:)
  real, allocatable    :: rivdph_pre_glb(:)
  real, allocatable    :: fldelv1_glb(:)
  real, allocatable    :: fldsto_glb(:)
  real, allocatable    :: flddph1_glb(:)
  real, allocatable    :: flddph_pre_glb(:)
  real, allocatable    :: rivout_glb(:)
  real, allocatable    :: fldout_glb(:)
  real, allocatable    :: next_glb(:)

  real, allocatable    :: uprivout_glb(:)
  real, allocatable    :: upfldout_glb(:)


  integer              :: ic1,ic2
  integer              :: iloc(1)
  integer              :: status
  real                 :: tmp_value

  integer              :: ix,iy,ix1,iy1
! ================================================
  !If 2-way coupling, update river/floodplain water storage with LSM outputs at the first sub time step (i.e. it==1)
  if (HYMAP2_routing_struc(n)%enable2waycpl==1.and.it==1) then
    do ic=1,nseqall 
      if(outlet(ic)==mis)cycle
      if(runoff0(ic)/=mis.and.basflw0(ic)/=mis.and.grarea(ic)/=mis)then
        !convert unit from [kg m-2] to [m3]
        rivsto(ic)=(runoff0(ic)+basflw0(ic))*grarea(ic)/1e3
        fldsto(ic)=0.
      endif
    enddo

  !If 1-way coupled, allow to compute deep water infiltration, evaporation from open waters and delay
  elseif(HYMAP2_routing_struc(n)%enable2waycpl==0) then
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

 !convert units [kg m-2 sec-1] to [m3 sec-1]
 where(runoff0/=mis.and.grarea/=mis)
    runoff1=runoff1*grarea/1e3
    basflw1=basflw1*grarea/1e3
 else where
    runoff1=0.
    basflw1=0.
 end where

 !Store runoff and baseflow in linear reservoirs
 do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    if(linres==1)then
       call HYMAP2_linear_reservoir_lis(dt,mis,runoff1(ic),basflw1(ic),&
            trnoff(ic),tbsflw(ic),cntime(ic),roffsto(ic),&
            basfsto(ic),runoff(ic))
    elseif(linres==0)then
       call HYMAP2_no_reservoir_lis(dt,mis,runoff1(ic),basflw1(ic),&
            roffsto(ic),basfsto(ic),runoff(ic))      
    else
       write(LIS_logunit,*)"HYMAP routing model linear reservoir flag: unknown value"
       call LIS_endrun()
    endif
 enddo

 do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    if(evapflag.ne.0)then
       !calculate evaporation from floodplains
       call HYMAP2_calc_evap_fld(dt,grarea(ic),fldare(ic),rivare(ic),&
            fldsto(ic),rivsto(ic),evpdif(ic),evpout(ic))
    else
       !does not calculate evaporation from floodplains
       evpout(ic)=0.
    endif
 enddo
  elseif(HYMAP2_routing_struc(n)%enable2waycpl/=1) then
       write(LIS_logunit,*)"Please check 2-way coupling option"
       call LIS_endrun()
  endif
 ! ================================================
 do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    call HYMAP2_calc_fldstg_1(nz,grarea(ic),rivlen(ic),rivwth(ic),       &
         rivstomax(ic),fldstomax(ic,:),fldhgt(ic,:),            &
         rivsto(ic),fldsto(ic),rivdph(ic),flddph(ic),flddph1(ic),   &
         fldelv1(ic),fldwth(ic),fldfrc(ic),fldare(ic),          &
         rivelv(ic),rivare(ic)) 
 enddo

 allocate(rivelv_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(rivsto_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(rivdph_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(rivdph_pre_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(fldelv1_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(fldsto_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(flddph1_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(flddph_pre_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(next_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(rivout_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(fldout_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(uprivout_glb(LIS_rc%glbnroutinggrid(n)))
 allocate(upfldout_glb(LIS_rc%glbnroutinggrid(n)))

 call HYMAP2_gather_tiles(n,rivelv,rivelv_glb)
 call HYMAP2_gather_tiles(n,rivsto,rivsto_glb)
 call HYMAP2_gather_tiles(n,rivdph,rivdph_glb)
 call HYMAP2_gather_tiles(n,rivdph_pre,rivdph_pre_glb)
 call HYMAP2_gather_tiles(n,fldelv1,fldelv1_glb)
 call HYMAP2_gather_tiles(n,fldsto,fldsto_glb)
 call HYMAP2_gather_tiles(n,flddph1,flddph1_glb)
 call HYMAP2_gather_tiles(n,flddph_pre,flddph_pre_glb)
 call HYMAP2_gather_tiles(n,real(next),next_glb)
  
 rivout_glb=0.
 fldout_glb=0.
 uprivout_glb=0.
 upfldout_glb=0.
  
 ! ================================================
 do ic=1,nseqall 
    if(outlet(ic)==mis)cycle
    
    if(outlet(ic)==0)then
       ic_down=next(ic)
       rivelv_down=rivelv_glb(ic_down)
       rivsto_down=rivsto_glb(ic_down)
       rivdph_down=rivdph_glb(ic_down)
       rivdph_pre_down=rivdph_pre_glb(ic_down)
       fldelv1_down=fldelv1_glb(ic_down)
       fldsto_down=fldsto_glb(ic_down)
       flddph1_down=flddph1_glb(ic_down)
       flddph_pre_down=flddph_pre_glb(ic_down)

       rivout_down=rivout_glb(ic_down)
       fldout_down=fldout_glb(ic_down)
       uprivout=uprivout_glb(ic)
       upfldout=upfldout_glb(ic)
       uprivout_down=uprivout_glb(ic_down)
       upfldout_down=upfldout_glb(ic_down)
 
    elseif(outlet(ic)==1)then

       rivelv_down=rivelv(ic)
       rivsto_down=rivsto(ic)
       rivdph_down=rivhgt(ic)
       rivdph_pre_down=rivhgt(ic)
       fldelv1_down=fldelv1(ic)
       fldsto_down=fldsto(ic)
       flddph1_down=fldelv1(ic)
       flddph_pre_down=fldelv1(ic)

       rivout_down=0.
       fldout_down=0.
       uprivout=uprivout_glb(ic)
       upfldout=upfldout_glb(ic)
       uprivout_down=0.
       upfldout_down=0.

    else
       write(LIS_logunit,*)"[ERR] Wrong outlet id"
       call LIS_endrun()
    endif
    if(flowmap(ic)==1)then
       !Calculate river flow based on the kinematic wave equation
       call HYMAP2_calc_rivout_kine(outlet(ic),dt,rivelv(ic),rivelv_down,nxtdst(ic),&
            rivwth(ic),sfcelv(ic),rivlen(ic),rivman(ic),slpmin, &
            rivsto(ic),rivdph(ic),rivout(ic),rivvel(ic))
       !Calculate floodplain
       if(floodflag==1)then
          !Calculate floodplain flow based on the kinematic wave equation
          call HYMAP2_calc_rivout_kine(outlet(ic),dt,fldelv1(ic),fldelv1_down,&
               nxtdst(ic),fldwth(ic),sfcelv(ic),rivlen(ic),fldman(ic),slpmin,&
               fldsto(ic),flddph1(ic),fldout(ic),fldvel(ic))
       elseif(floodflag==0)then
          !No floodplain dynamics
          call HYMAP2_no_floodplain_flow(fldelv1(ic),flddph1(ic),fldout(ic),&
               fldvel(ic),sfcelv0(ic),fldout_pre(ic),flddph_pre(ic))
       else
          write(LIS_logunit,*)"HYMAP floodplain dynamics: unknown value"
       endif
    elseif(flowmap(ic)==2.or.flowmap(ic)==3.or.flowmap(ic)==4.or.flowmap(ic)==5)then
       !Calculate river flow based on the local inertia wave equation
       call HYMAP2_calc_rivout_iner(outlet(ic),dt,rivelv(ic),rivelv_down,&
            elevtn(ic),nxtdst(ic), rivwth(ic),rivsto(ic),rivsto_down,&
            rivdph(ic),rivdph_down,rivlen(ic),rivman(ic),&
            grv,rivout(ic),rivvel(ic),sfcelv(ic), &
            rivout_pre(ic),rivdph_pre(ic),rivdph_pre_down)!,&
            !uprivout,uprivout_down,rivout_down)

        !Calculate floodplain
       if(floodflag==1)then
          !Calculate floodplain flow based on the local inertia wave equation
          call HYMAP2_calc_rivout_iner(outlet(ic),dt,fldelv1(ic),fldelv1_down,&
               elevtn(ic),nxtdst(ic),fldwth(ic),fldsto(ic),fldsto_down,&
               flddph1(ic),flddph1_down,rivlen(ic),fldman(ic),&
               grv,fldout(ic),fldvel(ic),sfcelv0(ic),  &
               fldout_pre(ic),flddph_pre(ic),flddph_pre_down)!,&
               !upfldout,upfldout_down,fldout_down)
       elseif(floodflag==0)then
          !No floodplain dynamics
          call HYMAP2_no_floodplain_flow(fldelv1(ic),flddph1(ic),&
               fldout(ic),fldvel(ic),sfcelv0(ic),fldout_pre(ic),flddph_pre(ic))
       else
          write(LIS_logunit,*)"HYMAP floodplain dynamics: unknown value"
       endif
    else
       write(LIS_logunit,*)"HYMAP routing method: unknown value",ic,flowmap(ic)
       call LIS_endrun()
    endif
    
    !set the updated variable back into the global one
    call HYMAP2_map_l2g_index(n, ic,icg)
    rivdph_pre_glb(icg) = rivdph_pre(ic)
    rivelv_glb(icg) = rivelv(ic)
    rivsto_glb(icg) = rivsto(ic)
    rivdph_glb(icg) = rivdph(ic)
    rivdph_pre_glb(icg) = rivdph_pre(ic)
    fldelv1_glb(icg) = fldelv1(ic)
    fldsto_glb(icg) = fldsto(ic)
    flddph1_glb(icg) = flddph1(ic)
    flddph_pre_glb(icg) = flddph_pre(ic)
    rivout_glb(icg) = rivout(ic)
    fldout_glb(icg) = fldout(ic)

    if(outlet(ic)==0)then
      uprivout_glb(ic_down)=uprivout_glb(ic_down)+uprivout_down
      upfldout_glb(ic_down)=upfldout_glb(ic_down)+upfldout_down
    endif
    
  enddo
! ================================================
  !reservoir simulation
  if(resopflag==1)then
    do ic=1,nseqall 
      call HYMAP2_map_l2g_index(n, ic,icg)
      if(minval(abs(HYMAP2_routing_struc(n)%resopaltloc-icg))==0)then
        iloc=minloc(abs(HYMAP2_routing_struc(n)%resopaltloc-icg))

        inflow=sum(rivout_glb+fldout_glb,next_glb==icg)
        call HYMAP2_resop_main_1(mis,nseqall,nz,time,dt,next(ic),rivsto(ic),fldsto(ic),runoff(ic),&
                        rivout(ic),fldout(ic),rivvel(ic),fldvel(ic),&
                        elevtn(ic),fldhgt(ic,:),fldstomax(ic,:),grarea(ic),rivstomax(ic),rivelv(ic),rivlen(ic),rivwth(ic),inflow,&
                        HYMAP2_routing_struc(n)%ntresop,iloc(1),HYMAP2_routing_struc(n)%tresopalt(iloc(1),:),HYMAP2_routing_struc(n)%resopalt(iloc(1),:),&
                        HYMAP2_routing_struc(n)%resopoutmin(iloc(1)))
      endif
    enddo
  endif
  ! ================================================

  call HYMAP2_gather_tiles(n,rivout,rivout_glb)
  call HYMAP2_gather_tiles(n,fldout,fldout_glb)
  call HYMAP2_gather_tiles(n,rivsto,rivsto_glb)
  call HYMAP2_gather_tiles(n,fldsto,fldsto_glb)

  if(LIS_masterproc) then 
     do ic=1,LIS_rc%glbnroutinggrid(n) 
        if(HYMAP2_routing_struc(n)%outlet_glb(ic)==mis)cycle
        
        if(HYMAP2_routing_struc(n)%outlet_glb(ic)==0)then
           ic_down=HYMAP2_routing_struc(n)%next_glb(ic)
           rivsto_down=rivsto_glb(ic_down)
           fldsto_down=fldsto_glb(ic_down)          
        elseif(HYMAP2_routing_struc(n)%outlet_glb(ic)==1)then
           rivsto_down=1e20 !rivsto(ic)
           fldsto_down=1e20 !fldsto(ic)
        else
           write(LIS_logunit,*)"Wrong outlet id"
           call LIS_endrun()
        endif

        !Update water storage in river and floodplain reservoirs
        call HYMAP2_calc_stonxt(HYMAP2_routing_struc(n)%outlet_glb(ic),&
             dt,rivout_glb(ic),fldout_glb(ic),&
             rivsto_glb(ic),fldsto_glb(ic),rivsto_down,fldsto_down)
        if(HYMAP2_routing_struc(n)%outlet_glb(ic)==0)then
           rivsto_glb(ic_down)=rivsto_down
           fldsto_glb(ic_down)=fldsto_down
        endif
     enddo 
  endif

#if (defined SPMD)
  call MPI_BCAST(rivout_glb, &
       LIS_rc%glbnroutinggrid(n), &
       MPI_REAL,0, &
       LIS_mpi_comm, status)
  
  call MPI_BCAST(fldout_glb, &
       LIS_rc%glbnroutinggrid(n), &
       MPI_REAL,0, &
       LIS_mpi_comm, status)

  call MPI_BCAST(rivsto_glb, &
       LIS_rc%glbnroutinggrid(n), &
       MPI_REAL,0, &
       LIS_mpi_comm, status)

  call MPI_BCAST(fldsto_glb, &
       LIS_rc%glbnroutinggrid(n), &
       MPI_REAL,0, &
       LIS_mpi_comm, status)

#endif

  call HYMAP2_map_g2l(n, rivout_glb,rivout)
  call HYMAP2_map_g2l(n, fldout_glb,fldout)
  call HYMAP2_map_g2l(n, rivsto_glb,rivsto)
  call HYMAP2_map_g2l(n, fldsto_glb,fldsto)

! ================================================
  !If 1-way coupled, add LSM outputs (surface runoff and baseflow) at the end of the time step
  if (HYMAP2_routing_struc(n)%enable2waycpl==0) then
    !Calculate runoff in the river network for the current time step
    do ic=1,nseqall 
       if(outlet(ic)==mis)cycle
       call HYMAP2_calc_runoff(dt,fldfrc(ic),runoff(ic),rivsto(ic),fldsto(ic))
     enddo
  endif
  ! ================================================
  ! Calculate surface water storage [mm]
  do ic=1,nseqall 
     if(outlet(ic)==mis)cycle
     surfws(ic)=1e3*(rivsto(ic)+fldsto(ic))/grarea(ic)
  enddo

  deallocate(rivelv_glb)
  deallocate(rivdph_glb)
  deallocate(rivdph_pre_glb)
  deallocate(fldelv1_glb)
  deallocate(flddph1_glb)
  deallocate(flddph_pre_glb)
  deallocate(next_glb)
  deallocate(rivout_glb)
  deallocate(rivsto_glb)
  deallocate(fldsto_glb)
  deallocate(fldout_glb)

end subroutine HYMAP2_model_core
! ================================================
! ================================================
!BOP
!
! !ROUTINE: HYMAP2_gather_tiles
! \label{HYMAP2_gather_tiles}
! 
! !INTERFACE:
subroutine HYMAP2_gather_tiles(n,var,var_glb)
! !USES:
  use LIS_coreMod
  use LIS_routingMod
  use LIS_mpiMod
  use HYMAP2_routingMod
!
! !DESCRIPTION: 
!  This subroutine gathers an individual variable
!  across different processors into a global array
!EOP

  implicit none

  integer        :: n 
  real           :: var(LIS_rc%nroutinggrid(n))
  real           :: var_glb(LIS_rc%glbnroutinggrid(n))

  real           :: tmpvar(LIS_rc%glbnroutinggrid(n))
  integer        :: i,l,ix,iy,ix1,iy1
  integer        :: status

#if (defined SPMD)
  call MPI_ALLGATHERV(var,&
       LIS_rc%nroutinggrid(n),&
       MPI_REAL,tmpvar,&
       LIS_routing_gdeltas(n,:),&
       LIS_routing_goffsets(n,:),&
       MPI_REAL,LIS_mpi_comm,status)
#endif
  !rearrange them to be in correct order.
  do l=1,LIS_npes
     do i=1,LIS_routing_gdeltas(n,l-1)
        ix = HYMAP2_routing_struc(n)%seqx_glb(i+&
             LIS_routing_goffsets(n,l-1))
        iy = HYMAP2_routing_struc(n)%seqy_glb(i+&
             LIS_routing_goffsets(n,l-1))
        ix1 = ix + LIS_ews_halo_ind(n,l) - 1
        iy1 = iy + LIS_nss_halo_ind(n,l)-1
        var_glb(HYMAP2_routing_struc(n)%sindex(ix1,iy1)) = &
             tmpvar(i+LIS_routing_goffsets(n,l-1))
     enddo
  enddo

end subroutine HYMAP2_gather_tiles

!BOP
! !ROUTINE: HYMAP2_map_g2l
! \label{HYMAP2_map_g2l}
! 
! !INTERFACE:
subroutine HYMAP2_map_g2l(n, var_glb,var_local)
! !USES:
  use LIS_coreMod
  use HYMAP2_routingMod
! 
! !DESCRIPTION:
! This subroutine maps a global array in the HYMAP2
! tile space to the local processor space. 
! 
!EOP
  implicit none

  integer             :: n 
  real                :: var_glb(LIS_rc%glbnroutinggrid(n))
  real                :: var_local(LIS_rc%nroutinggrid(n))

  integer             :: i, ix,iy,ix1,iy1,jx,jy

  do i=1,LIS_rc%nroutinggrid(n)
     ix = HYMAP2_routing_struc(n)%seqx(i)
     iy = HYMAP2_routing_struc(n)%seqy(i)
     ix1 = ix + LIS_ews_halo_ind(n,LIS_localPet+1) -1
     iy1 = iy + LIS_nss_halo_ind(n,LIS_localPet+1) -1
     var_local(i)  = var_glb(HYMAP2_routing_struc(n)%sindex(ix1,iy1))
  enddo

end subroutine HYMAP2_map_g2l

!BOP
! !ROUTINE: HYMAP2_map_l2g_index
! \label{HYMAP2_map_l2g_index}
! 
! !INTERFACE: 
subroutine HYMAP2_map_l2g_index(n, local_index,glb_index)
! !USES:
  use LIS_coreMod
  use HYMAP2_routingMod
!
! !DESCRIPTION: 
!  This subroutine converts the local tile index into the 
!  the global index. 
! 
!EOP
  implicit none

  integer             :: n 
  integer             :: local_index
  integer             :: glb_index

  integer             :: i, ix,iy,ix1,iy1,jx,jy

  ix = HYMAP2_routing_struc(n)%seqx(local_index)
  iy = HYMAP2_routing_struc(n)%seqy(local_index)
  ix1 = ix + LIS_ews_halo_ind(n,LIS_localPet+1) -1
  iy1 = iy + LIS_nss_halo_ind(n,LIS_localPet+1) -1
  glb_index = HYMAP2_routing_struc(n)%sindex(ix1,iy1)

end subroutine HYMAP2_map_l2g_index
