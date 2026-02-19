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
! !SUBROUTINE: HYMAP3_model_core
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
! 19 Jan 2016: Augusto Getirana, Inclusion of the local inertia formulation
!                                and adaptive time step.
! 13 Apr 2016: Augusto Getirana, Inclusion of option for hybrid runs with a
!                                river flow map.
! 13 May 2019: Augusto Getirana, Fixed reservoir operation module
!  7 Sep 2019: Augusto Getirana,  Added support for 2-way coupling
! 10 Mar 2020: Augusto Getirana, Fixed river/floodplain inflow
! 27 Apr 2020: Augusto Getirana,  Added support for urban drainage
!  8 Apr 2022: Augusto Getirana,  Added support for bifurcation
!
#include "LIS_misc.h"
subroutine HYMAP3_model_core(n, it, mis, nseqall, nz, time, dt,  &
     flowmap, linres, evapflag, resopflag, floodflag, dwiflag,   &
     rivout_pre, rivdph_pre, grv,                                &
     fldout_pre, flddph_pre, fldelv1,                            &
     outlet, next, elevtn, nxtdst, grarea,                       &
     fldgrd, fldman, fldhgt, fldstomax,                          &
     rivman, rivelv, rivstomax,                                  &
     rivlen, rivwth, rivhgt, rivare, slpmin,                     &
     trnoff, tbsflw, cntime,                                     &
     runoff0, basflw0, evpdif, rnfdwi_ratio, bsfdwi_ratio,       &
     rivsto, rivdph, rivvel, rivout,evpout,                      &
     fldout, fldsto, flddph, fldvel, fldfrc,                     &
     fldare, sfcelv, roffsto, basfsto, rnfdwi, bsfdwi, surfws,   &
     !ag (27Apr2020)
     !urban drainage variables/parameters
     flowtype, drvel, drtotwth, drrad, drman, drslp, drtotlgh,   &
     drnoutlet, drstomax, drout, drsto,                          &
     !ag(25feb2023)
     levhgt, levstomax, fldonlystomax, fldstoatlev)
  ! ================================================
  !HYMAP3 core
  !Augusto Getirana
  !6th Dec 2010
  !LEGOS/CNES, Toulouse, France
  !Adapted for flow routing implementation in LIS 09 Nov 2011
  ! ================================================
  !ag(1May2021)
  !ag(8Apr2022)
  use HYMAP3_bifMod
  use HYMAP3_managMod
  use HYMAP3_modelMod
  use HYMAP3_resopMod
  !ag(4Apr2025)
  use HYMAP3_resopMod_Yassin
  use HYMAP3_routingMod
  !ag(27Apr2020)
  use HYMAP3_urbanMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_mpiMod

  implicit none

  real*8,  intent(in)  :: time         !current time [year]
  real,    intent(in)  :: dt           !time step length [s]
  integer, intent(in)  :: n
  integer, intent(in)  :: it           !sub time step
  integer, intent(in)  :: nseqall      !length of 1D sequnece for river and mouth
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
  real,    intent(in)  :: rivwth(nseqall)       !river width [m]
  real,    intent(in)  :: rivhgt(nseqall)       !river heihgt [m]
  real,    intent(in)  :: rivare(nseqall)       !river surface area [m2]
  real,    intent(in)  :: trnoff(nseqall)       !runoff   concentration time parameter [day]
  real,    intent(in)  :: tbsflw(nseqall)       !baseflow concentration time parameter [day]
  real,    intent(in)  :: cntime(nseqall)       !concentration time [s]
  real,    intent(inout) :: runoff0(nseqall)      !input runoff [mm.idt-1]
  real,    intent(inout) :: basflw0(nseqall)      !input baseflow [mm.idt-1]
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
  !ag (27Apr2020)
  !urban drainage variables/parameters
  integer, intent(in)  :: flowtype     !urban drainage flag: 4 - compute urban drainage
  real,    intent(in)  :: drvel
  real,    intent(in)  :: drtotwth(nseqall)
  real,    intent(in)  :: drrad
  real,    intent(in)  :: drman
  real,    intent(in)  :: drslp
  real,    intent(in)  :: drtotlgh(nseqall)
  real,    intent(in)  :: drnoutlet(nseqall)
  real,    intent(in)  :: drstomax(nseqall)    !maximum urban drainage water storage capacity [m3]
  real,    intent(inout)    :: drsto(nseqall)   !urban drainage water storage [m3]
  real,    intent(inout)    :: drout(nseqall)   !urban drainage outflow [m3/s]
  !ag(25feb2023)
  !levee parameters
  real,    intent(in)  :: levhgt(nseqall)
  real,    intent(in)  :: levstomax(nseqall)
  real,    intent(in)  :: fldonlystomax(nseqall,nz)
  real,    intent(in)  :: fldstoatlev(nseqall)

  real                 :: flddph1(nseqall)
  real                 :: fldwth(nseqall)
  real                 :: runoff(nseqall)     !runoff+baseflow after linear reservoirs [m3/s]
  real                 :: inflow              !inflow from upstream grid cells [m3/s]
  real                 :: sfcelv0(nseqall)    !averaged floodplain surface elevation [m]
  integer              :: ic,ic_down,icg,icg_down,icg_down1
  real                 :: rivelv_down
  real                 :: rivsto_down
  real                 :: rivdph_down
  real                 :: rivdph_pre_down
  real                 :: fldelv1_down
  real                 :: fldsto_down
  real                 :: flddph1_down
  real                 :: flddph_pre_down
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

  !ag (14Jun2024)
  real, allocatable    :: runoff_glb(:)

  !ag (2Mar2020)
  real                             :: rivinf(nseqall)   !local river inflow  [m3/s]
  real                             :: fldinf(nseqall)   !local floodplain inflow  [m3/s]
  real, allocatable    :: rivinf_glb(:)   !global river inflow  [m3/s]
  real, allocatable    :: fldinf_glb(:)   !global floodplain inflow  [m3/s]
  integer              :: iloc(1)
  integer              :: status

  real                      :: drinf(nseqall)   !urban drainage inflow  [m3/s]
  real, allocatable    :: drsto_glb(:)  !global urban drainage water storage  [m3]
  real, allocatable    :: drout_glb(:)  !global urban drainage outflow  [m3/s]
  real, allocatable    :: drinf_glb(:)   !global urban drainage inflow  [m3/s]

  !ag (1sep2020)
  !direct insertion
  real                 :: disdif
  real*8               :: disobs

  !ag(1apr2021)
  real*8, allocatable  :: sealevelobs(:)
  integer              :: isea

  !ag(03Aug2022)
  real                 :: tmpsto,tmpman,tmpdph,tmpwth,tmpsto_down,tmpdph_down

  !ag(8Apr2022)
  integer              :: ibif,ielv
  real                 :: bifsto1,bifsto1_down
  real                 :: bifelv1,sfcelv1,bifout1

  !ag(10Oct2024)
  !resop parallel
  integer              :: iresop

  !ag(8Nov2024)
  !slope constrainsts for model stability
  real, parameter :: rslpmax    = 1.
  real, parameter :: bckslpmax  = -1.e-5
  real, parameter :: fslpmax    = 1.e-5

  !ag (5Apr2025)
  real, allocatable    :: sfcelv_glb(:)

#if (defined SPMD)
  external :: MPI_BCAST
#endif

  ! ================================================
  !convert units [kg m-2 sec-1] to [m3 sec-1] and aggregate surface
  !runoff and baseflow
  do ic=1,nseqall
     if(runoff0(ic)/=mis.and.basflw0(ic)/=mis.and.grarea(ic)/=mis)then
        runoff(ic)=(runoff0(ic)+basflw0(ic))*grarea(ic)/1e3
     else
        runoff(ic)=0.
     endif
  enddo

  !23 Nov 2016
  !Deep water infiltration
  if(dwiflag==1)then
     do ic=1,nseqall
        call HYMAP3_dwi(runoff0(ic),basflw0(ic),grarea(ic), &
             rnfdwi_ratio(ic), &
             bsfdwi_ratio(ic),rnfdwi(ic),bsfdwi(ic),runoff(ic))
     enddo
  endif

  !Store runoff and baseflow in linear reservoirs
  if(linres==1)then
     do ic=1,nseqall
        call HYMAP3_linear_reservoir_lis(dt,mis,runoff0(ic),basflw0(ic),&
             grarea(ic),trnoff(ic),tbsflw(ic),cntime(ic),roffsto(ic),&
             basfsto(ic),runoff(ic))
     enddo
  endif

  !calculate evaporation from floodplains
  if(evapflag.ne.0)then
     do ic=1,nseqall
        call HYMAP3_calc_evap_fld(dt,grarea(ic),fldare(ic),rivare(ic),&
             fldsto(ic),rivsto(ic),evpdif(ic),evpout(ic))
     enddo
  endif

  ! ================================================
  do ic=1,nseqall
     if(outlet(ic)==mis)cycle
     call HYMAP3_calc_fldstg_lev(nz,grarea(ic),rivlen(ic),rivwth(ic), &
          rivstomax(ic),&
          fldstomax(ic,:),fldhgt(ic,:),rivare(ic),levhgt(ic), &
          levstomax(ic),&
          fldonlystomax(ic,:),fldstoatlev(ic),elevtn(ic),&
          rivsto(ic),fldsto(ic),rivdph(ic),flddph(ic),flddph1(ic),&
          fldelv1(ic),&
          fldwth(ic),fldfrc(ic),fldare(ic))
  enddo

  allocate(rivelv_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(rivsto_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(rivdph_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(rivdph_pre_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(fldelv1_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(fldsto_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(flddph1_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(flddph_pre_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(rivout_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(fldout_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(rivinf_glb(LIS_rc%glbnroutinggrid(n)))
  allocate(fldinf_glb(LIS_rc%glbnroutinggrid(n)))

  call HYMAP3_gather_tiles(n,rivelv,rivelv_glb)
  call HYMAP3_gather_tiles(n,rivsto,rivsto_glb)
  call HYMAP3_gather_tiles(n,rivdph,rivdph_glb)
  call HYMAP3_gather_tiles(n,rivdph_pre,rivdph_pre_glb)
  call HYMAP3_gather_tiles(n,fldelv1,fldelv1_glb)
  call HYMAP3_gather_tiles(n,fldsto,fldsto_glb)
  call HYMAP3_gather_tiles(n,flddph1,flddph1_glb)
  call HYMAP3_gather_tiles(n,flddph_pre,flddph_pre_glb)

  rivinf_glb=0.
  fldinf_glb=0.

  !ag (27Apr2020)
  !urban drainage global variables
  if(flowtype==4)then
     allocate(drsto_glb(LIS_rc%glbnroutinggrid(n)))
     allocate(drout_glb(LIS_rc%glbnroutinggrid(n)))
     allocate(drinf_glb(LIS_rc%glbnroutinggrid(n)))
     call HYMAP3_gather_tiles(n,drsto,drsto_glb)
     drout_glb=0.
     drinf_glb=0.
  endif
  !================================================
  !ag(13May2021)
  !get sealevel at outlets
  if(HYMAP3_routing_struc(n)%sealevelflag==1)then
     if(allocated(sealevelobs))deallocate(sealevelobs)
     allocate(sealevelobs(HYMAP3_routing_struc(n)%nsealevel))
     do isea=1,HYMAP3_routing_struc(n)%nsealevel
        call HYMAP3_resop_inter(time,HYMAP3_routing_struc(n)%ntsealevel, &
             HYMAP3_routing_struc(n)%tsealevel(isea,:),     &
             dble(HYMAP3_routing_struc(n)%sealevel(isea,:)),&
             mis,sealevelobs(isea))
     enddo
  endif
  !================================================
  do ic=1,nseqall
     call HYMAP3_map_l2g_index(n, ic,icg)
     icg_down=HYMAP3_routing_struc(n)%next_glb(icg)
     if(outlet(ic)==0)then
        rivelv_down=rivelv_glb(icg_down)
        rivsto_down=rivsto_glb(icg_down)
        rivdph_down=rivdph_glb(icg_down)
        rivdph_pre_down=rivdph_pre_glb(icg_down)
        fldelv1_down=fldelv1_glb(icg_down)
        fldsto_down=fldsto_glb(icg_down)
        flddph1_down=flddph1_glb(icg_down)
        flddph_pre_down=flddph_pre_glb(icg_down)
     elseif(outlet(ic)>0)then
        rivelv_down=min(-1.,rivelv(ic)-slpmin*nxtdst(ic))
        fldelv1_down=rivelv_down
        rivdph_pre_down=max(0.,-rivelv_down)
        flddph_pre_down=rivdph_pre_down
        rivdph_down=max(0.,-rivelv_down)
        flddph1_down=rivdph_down
        rivsto_down=1e20
        fldsto_down=1e20
        ! ================================================
        !temporary solution to speed up the run
        !ag(30Mar2021)
        !get sealevel at outlets
        if(HYMAP3_routing_struc(n)%sealevelflag==1.and. &
             HYMAP3_routing_struc(n)%outletid(ic)/= &
             HYMAP3_routing_struc(n)%imis)then
           isea=HYMAP3_routing_struc(n)%outletid(ic)
           rivdph_down=max(0.,sealevelobs(isea)-rivelv_down)
           flddph1_down=max(0.,sealevelobs(isea)-fldelv1_down)
           rivdph_pre_down=rivdph_down
           flddph_pre_down=flddph1_down
        endif
     ! ================================================
     else
        write(LIS_logunit,*)"[ERR] Wrong outlet id",outlet(ic)
        call LIS_endrun()
     endif
     ! ================================================
     if(flowmap(ic)==1)then
        !Calculate river flow based on the kinematic wave equation
        call HYMAP3_calc_rivout_kine(outlet(ic),dt,rivelv(ic), &
             rivelv_down,nxtdst(ic),&
             rivwth(ic),sfcelv(ic),rivlen(ic),rivman(ic),slpmin, &
             rivsto(ic),rivdph(ic),rivout(ic),rivvel(ic))
        !Calculate floodplain
        if(floodflag==1)then
           !Calculate floodplain flow based on the kinematic wave equation
           call HYMAP3_calc_rivout_kine(outlet(ic),dt,fldelv1(ic), &
                fldelv1_down,&
                nxtdst(ic),fldwth(ic),sfcelv(ic),rivlen(ic),fldman(ic), &
                slpmin,&
                fldsto(ic),flddph1(ic),fldout(ic),fldvel(ic))
        elseif(floodflag==0)then
           !No floodplain dynamics
           call HYMAP3_no_floodplain_flow(fldelv1(ic),flddph1(ic), &
                fldout(ic),&
                fldvel(ic),sfcelv0(ic),fldout_pre(ic),flddph_pre(ic))
        else
           write(LIS_logunit,*) &
                "[ERR] HYMAP3 floodplain dynamics: unknown value"
           call LIS_endrun()
        endif
        !ag (27Apr2020)
     elseif(flowmap(ic)==2.or.flowmap(ic)==3)then
        !ag(03Aug2022) Add option to merge the computation of river and
        !floodplain dynamics
        if(floodflag==2)then
           call HYMAP3_merge_river_floodplain(rivsto(ic),rivman(ic), &
                rivdph(ic),rivwth(ic),&
                fldsto(ic),fldman(ic),flddph1(ic),rivlen(ic), &
                rivdph_down,rivsto_down,fldsto_down,&
                tmpsto,tmpman,tmpdph,tmpwth,tmpsto_down,tmpdph_down)
           !Calculate river flow based on the local inertia wave equation
           call HYMAP3_calc_rivout_iner(time,ic,2,rslpmax,bckslpmax,mis,&
                outlet(ic),dt,rivelv(ic),rivelv_down,&
                elevtn(ic),nxtdst(ic), tmpwth,tmpsto,tmpsto_down,&
                tmpdph,tmpdph_down,rivlen(ic),tmpman,&
                grv,rivout(ic),rivvel(ic),sfcelv(ic), &
                rivout_pre(ic),rivdph_pre(ic),rivdph_pre_down)
           !No floodplain dynamics
           call HYMAP3_no_floodplain_flow(fldelv1(ic),flddph1(ic),&
                fldout(ic),fldvel(ic),sfcelv0(ic),fldout_pre(ic),&
                flddph_pre(ic))
        else
           !Calculate river flow based on the local inertia wave equation
           call HYMAP3_calc_rivout_iner(time,ic,1,rslpmax,bckslpmax,mis,&
                outlet(ic),dt,rivelv(ic),rivelv_down,&
                elevtn(ic),nxtdst(ic), rivwth(ic),rivsto(ic),rivsto_down,&
                rivdph(ic),rivdph_down,rivlen(ic),rivman(ic),&
                grv,rivout(ic),rivvel(ic),sfcelv(ic), &
                rivout_pre(ic),rivdph_pre(ic),rivdph_pre_down)

           !Calculate floodplain
           if(floodflag==1)then
              !Calculate floodplain flow based on the local inertia wave
              !equation
              call HYMAP3_calc_rivout_iner(time,ic,2,fslpmax,-fslpmax, &
                   rivout(ic),outlet(ic),dt,fldelv1(ic),fldelv1_down,&
                   elevtn(ic),nxtdst(ic),fldwth(ic),fldsto(ic),&
                   fldsto_down,&
                   flddph1(ic),flddph1_down,rivlen(ic),fldman(ic),&
                   grv,fldout(ic),fldvel(ic),sfcelv0(ic),  &
                   fldout_pre(ic),flddph_pre(ic),flddph_pre_down)
           elseif(floodflag==0)then
              !No floodplain dynamics
              call HYMAP3_no_floodplain_flow(fldelv1(ic),flddph1(ic),&
                   fldout(ic),fldvel(ic),sfcelv0(ic),fldout_pre(ic),&
                   flddph_pre(ic))
           else
              write(LIS_logunit,*)"HYMAP3 floodplain dynamics: unknown value"
              call LIS_endrun()
           endif
        endif
     !ag (27Apr2020)
     !Urban flood modeling
     elseif(flowmap(ic)==4)then
        !compute urban drainage outflow
        call HYMAP3_calc_urb_drain_out(dt,drrad,drman,drslp,drstomax(ic),&
             drtotlgh(ic),drnoutlet(ic),drsto(ic),drout(ic))

        !Calculate river flow based on the local inertia wave equation
        call HYMAP3_calc_rivout_iner(time,ic,1,rslpmax,bckslpmax,mis, &
             outlet(ic),dt,rivelv(ic),rivelv_down,&
             elevtn(ic),nxtdst(ic), rivwth(ic),rivsto(ic),rivsto_down,&
             rivdph(ic),rivdph_down,rivlen(ic),rivman(ic),&
             grv,rivout(ic),rivvel(ic),sfcelv(ic), &
             rivout_pre(ic),rivdph_pre(ic),rivdph_pre_down)

        !Calculate floodplain flow based on the local inertia wave equation
        call HYMAP3_calc_rivout_iner(time,ic,2,fslpmax,-fslpmax, &
             rivout(ic),outlet(ic),dt,fldelv1(ic),fldelv1_down,&
             elevtn(ic),nxtdst(ic),fldwth(ic),fldsto(ic),fldsto_down,&
             flddph1(ic),flddph1_down,rivlen(ic),fldman(ic),&
             grv,fldout(ic),fldvel(ic),sfcelv0(ic),  &
             fldout_pre(ic),flddph_pre(ic),flddph_pre_down)!,&
        !ag(07Jan2021)
        !Urban flood modeling over kinematic-dominated flow (use
        !kinematic wave)
     elseif(flowmap(ic)==5)then
        !compute urban drainage outflow
        call HYMAP3_calc_urb_drain_out(dt,drrad,drman,drslp,drstomax(ic),&
             drtotlgh(ic),drnoutlet(ic),drsto(ic),drout(ic))

        !Calculate river flow based on the kinematic wave equation
        call HYMAP3_calc_rivout_kine(outlet(ic),dt,rivelv(ic), &
             rivelv_down,nxtdst(ic),&
             rivwth(ic),sfcelv(ic),rivlen(ic),rivman(ic),slpmin, &
             rivsto(ic),rivdph(ic),rivout(ic),rivvel(ic))

        !Calculate floodplain flow based on the kinematic wave equation
        call HYMAP3_calc_rivout_kine(outlet(ic),dt,fldelv1(ic), &
             fldelv1_down,&
             nxtdst(ic),fldwth(ic),sfcelv(ic),rivlen(ic),fldman(ic), &
             slpmin,&
             fldsto(ic),flddph1(ic),fldout(ic),fldvel(ic))
     else
        write(LIS_logunit,*) &
             "HYMAP3 routing method: unknown value", ic, flowmap(ic)
        call LIS_endrun()
     endif
  enddo
  ! ================================================
  call HYMAP3_gather_tiles(n,rivout,rivout_glb)
  call HYMAP3_gather_tiles(n,fldout,fldout_glb)
  call HYMAP3_gather_tiles(n,rivdph_pre,rivdph_pre_glb)

  !get global grid inflow
  do icg=1,LIS_rc%glbnroutinggrid(n)
     if(HYMAP3_routing_struc(n)%outlet_glb(icg)==0)then
        icg_down=HYMAP3_routing_struc(n)%next_glb(icg)
        rivinf_glb(icg_down)=rivinf_glb(icg_down)+rivout_glb(icg)
        fldinf_glb(icg_down)=fldinf_glb(icg_down)+fldout_glb(icg)
     endif
  enddo
  call HYMAP3_map_g2l(n, rivinf_glb,rivinf)
  call HYMAP3_map_g2l(n, fldinf_glb,fldinf)
  ! ================================================
  !ag (8Apr2022): river bifurcation
  if(HYMAP3_routing_struc(n)%bifflag==1)then
     do ibif=1,HYMAP3_routing_struc(n)%nbif
        icg=HYMAP3_routing_struc(n)%bifloc(ibif,1)
        icg_down=HYMAP3_routing_struc(n)%bifloc(ibif,2)
        HYMAP3_routing_struc(n)%bifout(ibif)=0.
        do ielv=1,HYMAP3_routing_struc(n)%nbifelv
           bifelv1=HYMAP3_routing_struc(n)%bifelv(ibif)+ &
                HYMAP3_routing_struc(n)%bifdelv(ielv)
           sfcelv1=rivelv_glb(icg)+rivdph_glb(icg)
           if(HYMAP3_routing_struc(n)%bifwth(ibif,ielv)>0.and. &
                sfcelv1>bifelv1)then
              bifsto1=max(0.,rivsto_glb(icg)+fldsto_glb(icg)- &
                   HYMAP3_routing_struc(n)%bifsto(ibif,ielv)-&
                   (rivout_glb(icg)+fldout_glb(icg)+ &
                   HYMAP3_routing_struc(n)%bifout(ibif))*dt)
              bifsto1_down=max(0.,rivsto_glb(icg_down)+ &
                   fldsto_glb(icg_down)-(max(0.,rivout_glb(icg_down))+&
                   max(0.,fldout_glb(icg_down)))*dt)
              call HYMAP3_calc_bifout_iner(dt,rivelv_glb(icg),&
                   rivelv_glb(icg_down),&
                   bifelv1,HYMAP3_routing_struc(n)%biflen(ibif),&
                   HYMAP3_routing_struc(n)%bifwth(ibif,ielv),bifsto1,&
                   bifsto1_down,rivdph_glb(icg),rivdph_glb(icg_down),&
                   HYMAP3_routing_struc(n)%bifman(ielv),grv,bifout1, &
                   HYMAP3_routing_struc(n)%bifout_pre(ibif),&
                   rivdph_pre_glb(icg),rivdph_pre_glb(icg_down))
              HYMAP3_routing_struc(n)%bifout(ibif)= &
                   HYMAP3_routing_struc(n)%bifout(ibif)+bifout1
           endif
        enddo
        rivout_glb(icg)=rivout_glb(icg)+ &
             HYMAP3_routing_struc(n)%bifout(ibif)
        rivinf_glb(icg_down)=rivinf_glb(icg_down)+ &
             HYMAP3_routing_struc(n)%bifout(ibif)
     enddo
#if (defined SPMD)
     call MPI_BCAST(rivout_glb, &
          LIS_rc%glbnroutinggrid(n), &
          MPI_REAL,0, &
          LIS_mpi_comm, status)

     call MPI_BCAST(rivinf_glb, &
          LIS_rc%glbnroutinggrid(n), &
          MPI_REAL,0, &
          LIS_mpi_comm, status)
#endif
     call HYMAP3_map_g2l(n, rivout_glb,rivout)
     call HYMAP3_map_g2l(n, rivinf_glb,rivinf)
  endif
  ! ================================================
  !ag (27Apr2020)
  !Urban flood modeling
  if(flowtype==4)then
     call HYMAP3_gather_tiles(n,drout,drout_glb)
     if(LIS_masterproc) then
        !get global grid inflow
        do icg=1,LIS_rc%glbnroutinggrid(n)
           icg_down=HYMAP3_routing_struc(n)%next_glb(icg)
           if(icg_down>0)then
              if(HYMAP3_routing_struc(n)%droutlet_glb(icg)==0)then
                 drinf_glb(icg_down)=drinf_glb(icg_down)+drout_glb(icg)
              else
                 rivinf_glb(icg_down)=rivinf_glb(icg_down)+drout_glb(icg)
              endif
           endif
        enddo
     endif
#if (defined SPMD)
     call MPI_BCAST(drinf_glb, &
          LIS_rc%glbnroutinggrid(n), &
          MPI_REAL,0, &
          LIS_mpi_comm, status)
#endif
     call HYMAP3_map_g2l(n, drinf_glb,drinf)

     !Update water storage in urban drainage network
     do ic=1,nseqall
        if(flowmap(ic)==4)then
           call HYMAP3_calc_urban_drain_stonxt(dt,drvel,drtotwth(ic), &
                drstomax(ic),drout(ic),drinf(ic),rivdph(ic),drsto(ic), &
                rivsto(ic))
        endif
     enddo

  endif
  ! ================================================
  !reservoir simulation
  if(HYMAP3_routing_struc(n)%resopflag>0)then
     allocate(runoff_glb(LIS_rc%glbnroutinggrid(n)))
     call HYMAP3_gather_tiles(n,runoff,runoff_glb)
     call HYMAP3_gather_tiles(n,rivsto,rivsto_glb)
     call HYMAP3_gather_tiles(n,fldsto,fldsto_glb)
     do iresop=1,HYMAP3_routing_struc(n)%nresop
        if(HYMAP3_routing_struc(n)%resoploc(iresop)== &
             HYMAP3_routing_struc(n)%imis.and.&
             HYMAP3_routing_struc(n)%resoploc_dwn(iresop)== &
             HYMAP3_routing_struc(n)%imis)cycle
        icg=HYMAP3_routing_struc(n)%resoploc_glb(iresop)
        icg_down=HYMAP3_routing_struc(n)%resoploc_dwn_glb(iresop)
        inflow=rivinf_glb(icg)+fldinf_glb(icg)
        rivout_glb(icg)=rivout_glb(icg)+fldout_glb(icg)
        rivinf_glb(icg_down)=rivinf_glb(icg_down)+fldinf_glb(icg_down)
        fldout_glb(icg)=0.
        fldinf_glb(icg_down)=0.
        if(HYMAP3_routing_struc(n)%resopflag==1)then
           call HYMAP3_resop_main_glb(mis,nz,time,dt,inflow,&
                rivsto_glb(icg),&
                fldsto_glb(icg),&
                runoff_glb(icg),&
                rivelv_glb(icg_down),&
                rivdph_glb(icg_down),&
                HYMAP3_routing_struc(n)%elevtn_resop(iresop),&
                HYMAP3_routing_struc(n)%fldhgt_resop(iresop,:),&
                HYMAP3_routing_struc(n)%fldstomax_resop(iresop,:),&
                HYMAP3_routing_struc(n)%grarea_resop(iresop),&
                HYMAP3_routing_struc(n)%rivstomax_resop(iresop),&
                HYMAP3_routing_struc(n)%rivelv_resop(iresop),&
                HYMAP3_routing_struc(n)%rivlen_resop(iresop),&
                HYMAP3_routing_struc(n)%rivwth_resop(iresop),&
                HYMAP3_routing_struc(n)%resoptype(iresop),&
                HYMAP3_routing_struc(n)%ntresop,&
                HYMAP3_routing_struc(n)%tresopalt(iresop,:),&
                HYMAP3_routing_struc(n)%resopalt(iresop,:),&
                HYMAP3_routing_struc(n)%resopoutmin(iresop),&
                rivout_glb(icg),&
                rivinf_glb(icg_down))
        elseif(HYMAP3_routing_struc(n)%resopflag==2)then
           allocate(sfcelv_glb(LIS_rc%glbnroutinggrid(n)))
           call HYMAP3_gather_tiles(n,sfcelv,sfcelv_glb)
           call HYMAP3_resop_main_glb_Yassin(mis,nz,time,dt,inflow,&
                rivsto_glb(icg),&
                fldsto_glb(icg),&
                runoff_glb(icg),&
                sfcelv_glb(icg),&
                HYMAP3_routing_struc(n)%maxsto_resop(iresop),&
                HYMAP3_routing_struc(n)%inidis_resop(iresop),&
                HYMAP3_routing_struc(n)%inisto_resop(iresop),&
                HYMAP3_routing_struc(n)%dwndis_resop(iresop),&
                HYMAP3_routing_struc(n)%deadis_resop(iresop),&
                HYMAP3_routing_struc(n)%minsto_mo_resop(iresop,:),&
                HYMAP3_routing_struc(n)%nupsto_mo_resop(iresop,:),&
                HYMAP3_routing_struc(n)%uppsto_mo_resop(iresop,:),&
                HYMAP3_routing_struc(n)%mindis_mo_resop(iresop,:),&
                HYMAP3_routing_struc(n)%nupdis_mo_resop(iresop,:),&
                HYMAP3_routing_struc(n)%uppdis_mo_resop(iresop,:),&
                HYMAP3_routing_struc(n)%reg1_resop(iresop),&
                HYMAP3_routing_struc(n)%reg2_resop(iresop),&
                HYMAP3_routing_struc(n)%reg3_resop(iresop),&
                rivout_glb(icg),&
                rivinf_glb(icg_down))
        endif
        if(HYMAP3_routing_struc(n)%resoploc(iresop)/= &
             HYMAP3_routing_struc(n)%imis)then
           ic=HYMAP3_routing_struc(n)%resoploc(iresop)
           rivout(ic)=rivout_glb(icg)
           fldout(ic)=fldout_glb(icg)
        endif
        if(HYMAP3_routing_struc(n)%resoploc_dwn(iresop)/= &
             HYMAP3_routing_struc(n)%imis)then
           ic_down=HYMAP3_routing_struc(n)%resoploc_dwn(iresop)
           rivinf(ic_down)=rivinf_glb(icg_down)
           fldinf(ic_down)=fldinf_glb(icg_down)
        endif
     enddo
  endif

  ! ================================================
  !ag (30Apr2021): water management rules
  if(HYMAP3_routing_struc(n)%managflag==1)then
     if(LIS_masterproc) then
        do ic=1,HYMAP3_routing_struc(n)%nmanag
           icg=HYMAP3_routing_struc(n)%managloc(ic,1)
           icg_down1=HYMAP3_routing_struc(n)%managloc(ic,2)
           icg_down=HYMAP3_routing_struc(n)%next_glb(icg)
           rivout_glb(icg)=rivout_glb(icg)+fldout_glb(icg)
           rivinf_glb(icg_down1)=rivinf_glb(icg_down1)+ &
                fldinf_glb(icg_down1)
           rivinf_glb(icg_down)=rivinf_glb(icg_down)+fldinf_glb(icg_down)
           fldout_glb(icg)=0.
           fldinf_glb(icg_down1)=0.
           fldinf_glb(icg_down)=0.
           call HYMAP3_manag_rules(time, &
                HYMAP3_routing_struc(n)%managtype(ic),&
                HYMAP3_routing_struc(n)%nmanagcoef,&
                HYMAP3_routing_struc(n)%managqmax(ic),&
                HYMAP3_routing_struc(n)%managact(ic),&
                HYMAP3_routing_struc(n)%managcoef(ic,:),&
                rivout_glb(icg),rivinf_glb(icg_down), &
                rivinf_glb(icg_down1))
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
     call MPI_BCAST(rivinf_glb, &
          LIS_rc%glbnroutinggrid(n), &
          MPI_REAL,0, &
          LIS_mpi_comm, status)
     call MPI_BCAST(fldinf_glb, &
          LIS_rc%glbnroutinggrid(n), &
          MPI_REAL,0, &
          LIS_mpi_comm, status)
#endif
     call HYMAP3_map_g2l(n, rivout_glb,rivout)
     call HYMAP3_map_g2l(n, fldout_glb,fldout)
     call HYMAP3_map_g2l(n, rivinf_glb,rivinf)
     call HYMAP3_map_g2l(n, fldinf_glb,fldinf)
  endif
  ! ================================================
  !ag (9Aug2020)
  !ag (20Apr2021): fixes for global loop and rivout/fldout updates
  !direct streamflow insertion
  if(HYMAP3_routing_struc(n)%insertflag==1)then
     if(LIS_masterproc) then
        do icg=1,LIS_rc%glbnroutinggrid(n)
           if(minval(abs(HYMAP3_routing_struc(n)%insertloc-icg))==0)then
              iloc=minloc(abs(HYMAP3_routing_struc(n)%insertloc-icg))
              icg_down=HYMAP3_routing_struc(n)%next_glb(icg)
              rivout_glb(icg)=rivout_glb(icg)+fldout_glb(icg)
              rivinf_glb(icg_down)=rivinf_glb(icg_down)+ &
                   fldinf_glb(icg_down)
              fldout_glb(icg)=0.
              fldinf_glb(icg_down)=0.
              !get linearly intepolated observed streamflow
              call HYMAP3_resop_inter(time, &
                   HYMAP3_routing_struc(n)%ntinsert,&
                   HYMAP3_routing_struc(n)%tinsert(iloc(1),:),&
                   dble(HYMAP3_routing_struc(n)%insertdis(iloc(1),:)),&
                   mis,disobs)
              disdif=rivout_glb(icg)-real(disobs)
              rivout_glb(icg)=disobs
              rivinf_glb(icg_down)=rivinf_glb(icg_down)-disdif
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
     call MPI_BCAST(rivinf_glb, &
          LIS_rc%glbnroutinggrid(n), &
          MPI_REAL,0, &
          LIS_mpi_comm, status)
     call MPI_BCAST(fldinf_glb, &
          LIS_rc%glbnroutinggrid(n), &
          MPI_REAL,0, &
          LIS_mpi_comm, status)
#endif
     call HYMAP3_map_g2l(n, rivout_glb,rivout)
     call HYMAP3_map_g2l(n, fldout_glb,fldout)
     call HYMAP3_map_g2l(n, rivinf_glb,rivinf)
     call HYMAP3_map_g2l(n, fldinf_glb,fldinf)
  endif

  ! ================================================

  !Update water storage in rivers and floodplains
  do ic=1,nseqall
     call HYMAP3_calc_stonxt(dt,rivout(ic),fldout(ic),rivinf(ic), &
          fldinf(ic),rivsto(ic),fldsto(ic))
  enddo

  ! ================================================
  !Calculate runoff in the river network for the current time step
  do ic=1,nseqall
     if(outlet(ic)==mis)cycle
     call HYMAP3_calc_runoff(dt,fldfrc(ic),runoff(ic),rivsto(ic),&
          fldsto(ic))
  enddo

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
  deallocate(rivout_glb)
  deallocate(rivsto_glb)
  deallocate(fldsto_glb)
  deallocate(fldout_glb)
  !ag(27Apr2020)
  if(allocated(drinf_glb))deallocate(drinf_glb)
  !ag(11Oct2024)
  if(allocated(runoff_glb))deallocate(runoff_glb)
  !ag(5Apr2025)
  if(allocated(sfcelv_glb))deallocate(sfcelv_glb)

end subroutine HYMAP3_model_core
