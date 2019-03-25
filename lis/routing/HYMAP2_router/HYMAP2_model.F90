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
#include "LIS_misc.h"
subroutine HYMAP2_model(n,mis,nx,ny,yr,mo,da,hr,mn,ss,&
     nseqall,nz,dt,flowmap,linres,evapflag,  &
     rivout_pre,rivdph_pre,                  &
     fldout_pre,flddph_pre,fldelv1,          &
     grv,cadp,steptype,resopflag,floodflag,  &
     outlet,next,elevtn,nxtdst,grarea,       &
     fldgrd,fldman,fldhgt,fldstomax,         &
     rivman,rivelv,rivstomax,                &
     rivlen,rivwth,rivhgt,rivare,slpmin,     &
     trnoff,tbsflw,cntime,                   &
     dwiflag,rnfdwi_ratio,bsfdwi_ratio,      &
     runoff0,basflw0,evpdif,                 &
     rivsto,rivdph,rivvel,rivout,evpout,     &
     fldout,fldsto,flddph,fldvel,fldfrc,     &
     fldare,sfcelv,roffsto,basfsto,          &
     rnfdwi,bsfdwi,surfws,                   &
     dtaout                                  )                   
                   
  use HYMAP2_modelMod
  use HYMAP2_routingMod
  use LIS_mpiMod
  use LIS_coreMod
  
  implicit none

  real,    intent(in)    :: dt                    !time step length [s]
  integer, intent(in)    :: n
  integer, intent(in)    :: nx                    !number of grids in horizontal
  integer, intent(in)    :: ny                    !number of grids in vertical
  integer, intent(in)    :: yr                    !current year
  integer, intent(in)    :: mo                    !current month
  integer, intent(in)    :: da                    !current day
  integer, intent(in)    :: hr                    !current hour
  integer, intent(in)    :: mn                    !current min
  integer, intent(in)    :: ss                    !current second
  integer, intent(in)    :: nz                    !number of stages in the sub-grid discretization

  integer, intent(in)    :: linres                !linear reservoir flag: 1 - use linear reservoirs; or 0 - do not use
  integer, intent(in)    :: evapflag              !evaporation flag: 1 - do not compute ; or 2 - compute evapotation in floodplains

  real,    intent(in)    :: mis                   !real undefined value

  integer, intent(in)    :: nseqall               !length of 1D sequnece for river and mouth    
  integer, intent(in)    :: outlet(nseqall)       !outlet flag: 0 - river; 1 - ocean
  integer, intent(in)    :: next(nseqall)         !point downstream horizontal
  real,    intent(in)    :: elevtn(nseqall)       !river bed elevation [m]
  real,    intent(in)    :: nxtdst(nseqall)       !distance to the next grid [m]
  real,    intent(in)    :: grarea(nseqall)       !area of the grid [m^2] 
  real,    intent(in)    :: fldgrd(nseqall,nz)    !floodplain gradient [-]
  real,    intent(in)    :: fldman(nseqall)       !maning coefficient for floodplains [-]
  real,    intent(in)    :: fldhgt(nseqall,nz)    !floodplain height [m]
  real,    intent(in)    :: fldstomax(nseqall,nz) !maximum floodplain storage [m3]
  real,    intent(in)    :: rivman(nseqall)       !maning coefficient for rivers [-]
  real,    intent(in)    :: rivelv(nseqall)       !elevation of river bed [m]
  real,    intent(in)    :: rivstomax(nseqall)    !maximum river storage [m3]
  real,    intent(in)    :: rivlen(nseqall)       !river length [m] 
  real,    intent(in)    :: rivwth(nseqall)       !river wiidth [m]
  real,    intent(in)    :: rivhgt(nseqall)       !river heihgt [m]
  real,    intent(in)    :: rivare(nseqall)       !river surface area [m2]
  real,    intent(in)    :: trnoff(nseqall)       !runoff   concentration time parameter [day]
  real,    intent(in)    :: tbsflw(nseqall)       !baseflow concentration time parameter [day]
  real,    intent(in)    :: cntime(nseqall)       !concentration time [s]
  real,    intent(in)    :: slpmin                !minimum river slope [-]
  real,    intent(inout) :: runoff0(nseqall)      !input runoff [mm.idt-1]
  real,    intent(inout) :: basflw0(nseqall)      !input baseflow [mm.idt-1]
  real,    intent(in)    :: flowmap(nseqall)     !river flow type map: 1 - kinematic; 2 - diffusive; 3 - inertial

  real,    intent(inout) :: roffsto(nseqall)      !runoff   reservoir storage [m3]
  real,    intent(inout) :: basfsto(nseqall)      !baseflow reservoir storage [m3]
  real,    intent(inout) :: rivsto(nseqall)       !floodplain storage [m3]
  real,    intent(inout) :: fldsto(nseqall)       !floodplain storage   [m3]
  real,    intent(inout) :: surfws(nseqall)       !surface water storage   [mm]

  real,    intent(out)   :: rivdph(nseqall)       !river depth [m]
  real,    intent(out)   :: rivvel(nseqall)       !river flow velocity [m/s]
  real,    intent(out)   :: rivout(nseqall)       !average river outflow  [m3/s]
  real,    intent(inout) :: evpout(nseqall)       !effective evaporation from floodplains [m3]
  real,    intent(out)   :: fldout(nseqall)       !floodplain outflow  [m3/s]
  real,    intent(out)   :: flddph(nseqall)       !floodplain depth [m]
  real,    intent(out)   :: fldvel(nseqall)       !floodplain flow velocity [m/s]
  real,    intent(out)   :: fldfrc(nseqall)       !area fraction [-]
  real,    intent(out)   :: fldare(nseqall)       !flooded area [m2]
  real,    intent(out)   :: sfcelv(nseqall)       !water surface elevation [m]  

  real                   :: flddph1(nseqall)
  real                   :: fldinf(nseqall)
  real                   :: fldwth(nseqall)
  real                   :: runoff(nseqall)       !runoff+baseflow after linear reservoirs [m3/s]
  real                   :: etrans(nseqall)
  real                   :: evpdif(nseqall)
  real                   :: evpwat(nseqall)
  real                   :: rivout0(nseqall)       !adaptive time step river outflow  [m3/s]
  real                   :: rivvel0(nseqall)       !adaptive time step river flow velocity  [m/s]
  real                   :: fldout0(nseqall)       !adaptive time step floodplain outflow  [m3/s]
  real                   :: fldvel0(nseqall)       !adaptive time step floodplain flow velocity  [m/s]
  real                   :: evpout0(nseqall)       !adaptive time step evaporation from open waters  [mm]

  ! Local inertia variables
  real,    intent(inout) :: rivout_pre(nseqall)   !previous river outflow [m3/s]
  real,    intent(inout) :: rivdph_pre(nseqall)   !previous river depth [m]
  real,    intent(inout) :: fldout_pre(nseqall)   !previous flood outflow [m3/s]
  real,    intent(inout) :: flddph_pre(nseqall)   !previous flood depth [m]
  real,    intent(inout) :: fldelv1(nseqall)      !floodplain elevation [m]
  real,    intent(in)    :: grv                   !gravity accerelation [m/s2]
  real,    intent(in)    :: cadp                  !alfa coefficient for adaptative time step as described in Bates et al., (2010) [-]
  !ag (19Jan2016)
  !Adaptive time step 
  integer, intent(in)    :: steptype            !time step type: 0 - constant time step; 1 - adaptive time step
  integer                :: nt(nseqall)         !local number of time steps (>=1 if adaptive time step is active)
  real                   :: dta(nseqall)        !local time step [s] (<=dt if adaptive time step is active)
  real,    intent(inout) :: dtaout(nseqall)     !minimum local time step [s] (<=dt if adaptive time step is active)
  integer                :: it
  
  !Reservoir operation
  integer, intent(in)    :: resopflag
  
  !Floodplain dynamics
  integer, intent(in)    :: floodflag    !floodplain dynamics flag: 1 - simulate floodplain dynamics; or 2 - do not simulate floodplain dynamics

  !Deep water infiltration
  integer, intent(in)    :: dwiflag               !deep water infiltration flag
  real,    intent(in)    :: bsfdwi_ratio(nseqall) !deep water infiltration ratio from baseflow [-]
  real,    intent(in)    :: rnfdwi_ratio(nseqall) !deep water infiltration ratio from surface runoff [-]
  real,    intent(inout) :: bsfdwi(nseqall)       !deep water infiltration from baseflow [mm.idt-1]
  real,    intent(inout) :: rnfdwi(nseqall)       !deep water infiltration from surface runoff [mm.idt-1]

  !local variables
  real*8                   :: time                !current time
  integer                :: ic,ic_down,i,iloc(1)

  !ag (29Jan2016)
  integer                :: nt_local
  real                   :: dta_local
  integer                :: counti, countf, count_rate
  integer                :: maxi
  real                   :: mindt
  real                   :: tmp
  real*8                 :: dt1
  integer                :: status

#if 0 
  call system_clock(counti,count_rate)
#endif

! ================================================
  !ag 3 Apr 2014 
  !define sub time step
! ================================================
!  !ag 12 Dec 2016: get minimum dta values (dtaout)
!  do ic=1,nseqall
!    call HYMAP2_dynstp(dt,nxtdst(ic),rivdph(ic),grv,cadp,dta(ic),nt(ic))
!    if(dta(ic)>dt)then
!      dta(ic)=dt
!      nt(ic)=1
!    endif
!  enddo
!  where(dtaout>dta)dtaout=dta
! ================================================

  if(steptype==1)then
    dta=dt
    nt(:)=1
  elseif(steptype==2)then
    do ic=1,nseqall
      call HYMAP2_dynstp(dt,nxtdst(ic),rivdph(ic),grv,cadp,dta(ic),nt(ic))
      if(dta(ic)>dt)then
        dta(ic)=dt
        nt(ic)=1
      endif
    enddo
  endif
  !where(dtaout>dta)
  dtaout=dta
  iloc(:)=minloc(dta)
  i=iloc(1)

  if(nseqall.gt.0) then 
     nt_local = nt(i)
     dta_local = dta(i)
  else
     nt_local = 1
     dta_local = dt
  endif

!find the maximum of i and minimum dt across all processors
#if (defined SPMD)
  call MPI_ALLREDUCE(nt_local,maxi, 1, MPI_INTEGER, MPI_MAX, &
       LIS_mpi_comm,status)

  call MPI_ALLREDUCE(dta_local,mindt, 1, MPI_REAL, MPI_MIN, &
       LIS_mpi_comm,status)
#else 
  maxi = nt_local
  mindt = dta_local
#endif
  rivout0=rivout
  rivvel0=rivvel
  fldout0=fldout
  fldvel0=fldvel
  evpout0=0.
  rivout=0.
  rivvel=0.
  fldout=0.
  fldvel=0. 
  evpout=0. 

  do it=1,maxi
    !call HYMAP2_date2frac(yr,mo,da,hr,mn,ss,real(dta(i)*(it-1)-dt),time)

    dt1=real(mindt*it-dt)

    time=dble(yr*10000+mo*100+da) +&
         (dble(hr)+(dble(mn)+(dble(ss)+dt1)/60.)/60.)/24.

    call HYMAP2_model_core(n,mis,nseqall,nz,time,&
         mindt,flowmap,linres,evapflag, &
         resopflag,floodflag,dwiflag,                               &
         rivout_pre,rivdph_pre,grv,                 &
         fldout_pre,flddph_pre,fldelv1,                &
         outlet,next,elevtn,nxtdst,grarea,     & 
         fldgrd,fldman,fldhgt,fldstomax,     &
         rivman,rivelv,rivstomax,      &
         rivlen,rivwth,rivhgt,rivare,slpmin,       &
         trnoff,tbsflw,cntime,                         &
         runoff0,basflw0,evpdif,rnfdwi_ratio,bsfdwi_ratio,      &
         rivsto,rivdph,rivvel0,rivout0,evpout0,&
         fldout0,fldsto,flddph,fldvel0,fldfrc,     &
         fldare,sfcelv,roffsto,basfsto,rnfdwi,bsfdwi,surfws)
        
    do ic=1,nseqall 
      rivout0(ic)=rivout0(ic)+fldout0(ic)
!      rivout(ic)=rivout(ic)+rivout0(ic)/real(nt(i))
!      rivvel(ic)=rivvel(ic)+rivvel0(ic)/real(nt(i))
!      fldout(ic)=fldout(ic)+fldout0(ic)/real(nt(i))
!      fldvel(ic)=fldvel(ic)+fldvel0(ic)/real(nt(i))
!      evpout(ic)=evpout(ic)+evpout0(ic)/real(nt(i))

      rivout(ic)=rivout(ic)+rivout0(ic)/real(maxi)
      rivvel(ic)=rivvel(ic)+rivvel0(ic)/real(maxi)
      fldout(ic)=fldout(ic)+fldout0(ic)/real(maxi)
      fldvel(ic)=fldvel(ic)+fldvel0(ic)/real(maxi)
      evpout(ic)=evpout(ic)+evpout0(ic)/real(maxi)
    enddo

  enddo
  !do ic=1,nseqall 
  !  rivout0(ic)=rivout0(ic)+fldout0(ic)
  !  rivout(ic)=rivout0(ic)
  !  rivvel(ic)=rivvel0(ic)
  !  fldout(ic)=fldout0(ic)
  !  fldvel(ic)=fldvel0(ic)
  !  evpout(ic)=evpout0(ic)
  !enddo

#if 0 
  call system_clock(countf)
  HYMAP2_routing_struc(n)%dt_proc=HYMAP2_routing_struc(n)%dt_proc+real(countf-counti)/real(count_rate)
#endif
  tmp=minval(dtaout)
  !write(unit=HYMAP2_logunit,fmt='(8i5,10f15.4)')n,yr,mo,da,hr,mn,ss,nt(i),dta(i),sum(HYMAP2_routing_struc(:)%dt_proc),minval(dtaout,dtaout>tmp),maxval(dtaout),real(count(dtaout==minval(dtaout)))
!  write(unit=LIS_logunit,fmt='(a,i5,f10.2,f10.4,3f10.2)')'[INFO] HYMAP2_log: ',nt(i),dta(i),sum(HYMAP2_routing_struc(:)%dt_proc),minval(dtaout,dtaout>tmp),maxval(dtaout),real(count(dtaout==minval(dtaout)))

end subroutine HYMAP2_model
