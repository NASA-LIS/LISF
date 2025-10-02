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
! !MODULE HYMAP3_resopMod
!
! !DESCRIPTION
!
! !REVISION HISTORY:
! 27 Sep 2014: Augusto Getirana, Initial implementation
! 30 Jan 2016: Augusto Getirana, Adapt to run in LIS v7.1
! 13 May 2019: Augusto Getirana, Fixed bug in HYMAP3_get_volume_profile
!
module HYMAP3_resopMod
  !module defining reservoir operation variables
  !Augusto Getirana
  !NASA Goddard Space Flight Center
  !27 Sep 2014
  use HYMAP3_modelMod

  implicit none
  private

  public :: HYMAP3_resop_main_glb
  public :: HYMAP3_resop_inter

contains
  !=============================================
  !=============================================
  subroutine HYMAP3_resop_main_glb(mis, nz, time, dt, inflow, &
       rivsto, fldsto, runoff, rivelv_down, rivdph_down,      &
       elevtn, fldhgt, fldstomax, grarea, rivstomax,          &
       rivelv, rivlen, rivwth, resoptype, ntresop,            &
       tresopalt, resopalt, resopoutmin, outflow, inflow_down)

    use LIS_logMod

    implicit none

    real,    intent(in)    :: mis                 !missing data
    integer, intent(in)    :: nz                  !number of grids in vertical
    real*8,  intent(in)    :: time                !current time [year]
    real,    intent(in)    :: dt                  !time step length [s]
    real,    intent(in)    :: inflow              !grid cell inflow [m3/s]
    real,    intent(in)    :: rivsto              !floodplain storage [m3]
    real,    intent(in)    :: fldsto              !floodplain storage   [m3]
    real,    intent(in)    :: runoff              !total runoff in grid cell  [m3/s]
    real,    intent(in)    :: elevtn              !river bed elevation [m]
    real,    intent(in)    :: fldhgt(nz)          !floodplain height
    real,    intent(in)    :: fldstomax(nz)       !maximum floodplain storage [m3]
    real,    intent(in)    :: grarea              !area of the grid [m^2]
    real,    intent(in)    :: rivstomax           !maximum river storage [m3]
    real,    intent(in)    :: rivelv              !elevation of river bed [m]
    real,    intent(in)    :: rivlen              !river length [m]
    real,    intent(in)    :: rivwth              !river width [m]
    integer, intent(in)    :: resoptype           !reservoir operation flag (1: water level, 2: beaver dam)
    integer, intent(in)    :: ntresop
    real*8,  intent(in)    :: tresopalt(ntresop)
    real,    intent(in)    :: resopalt(ntresop)
    real,    intent(in)    :: resopoutmin
    real,    intent(inout) :: outflow              !grid cell outflow (river+floodplain) [m3/s]
    real,    intent(inout) :: inflow_down              !downstream grid cell inflow [m3/s]

    real*8                 :: elvobs,elvsim
    real*8                 :: volobs,volsim,volsim0
    real*8                 :: outsim

    !ag(17Apr2024)
    real*8,  parameter     :: man = 0.1          !manning roughness for beaver dam overflow
    real,    intent(in)    :: rivelv_down     !downstream river/floodplain bed elevation [m]
    real,    intent(in)    :: rivdph_down     !dowstream river/floodplain depth [m]
    real                   :: sfcelv_down     !downstream water surface elevation
    real*8                 :: outfrc,maxdamvol
    real*8                 :: slp,hgt

    if(resoptype==1)then
       volsim=max(0.,rivsto+fldsto+(inflow+runoff)*dble(dt))

       !get observed elevation from linear intepolation
       call HYMAP3_resop_inter(time,ntresop,tresopalt,dble(resopalt),mis,&
            elvobs)

       !calculate water volume as a function of observed elevation
       call HYMAP3_get_volume_profile(nz,dble(elevtn),dble(fldhgt), &
            dble(fldstomax),dble(grarea),dble(rivstomax),&
            dble(rivelv),dble(rivlen),dble(rivwth),elvobs,volobs)

       !update model variables based on water volume availability
       if(volsim>=volobs+dble(resopoutmin*dble(dt)))then
          outsim=max((volsim-volobs)/dble(dt),resopoutmin)
       else
          outsim=max(0.,min(volsim/dble(dt),resopoutmin))
       endif

       !calculate water volume and elevation as a function of new
       ! simulated outflow
       volsim0=max(0.,volsim-outsim*dble(dt))

       call HYMAP3_get_elevation_profile(nz,dble(elevtn),dble(fldhgt), &
            dble(fldstomax),dble(grarea),dble(rivstomax),&
            dble(rivelv),dble(rivlen),dble(rivwth),elvsim,volsim0)

    !ag(17Apr2024)
    elseif(resoptype==2)then

       !for beaver dam:
       ! resopalt stands for dam height, NOT dam elevation;
       ! resopoutmin stands for dam leakage fraction, not minimum outflow
       outfrc=resopoutmin !outflow fraction (beaver dam leakage fraction)
       elvobs=dble(rivelv)+dble(resopalt(1)) !riverbed elevation + dam height
       volsim=max(0.,rivsto+fldsto+(inflow+runoff)*dble(dt)) !reservoir storage
       !reservoir elevation
       call HYMAP3_get_elevation_profile(nz,dble(elevtn),dble(fldhgt), &
            dble(fldstomax),dble(grarea),dble(rivstomax),&
            dble(rivelv),dble(rivlen),dble(rivwth),elvsim,volsim)
       if(elvsim<=elvobs)then
          outsim=min(outflow,outflow*outfrc) !dam leak
       else
          sfcelv_down=rivelv_down+rivdph_down
          slp=max(1e-5,(elvsim-dble(sfcelv_down))/dble(rivlen)) !compute slope
          hgt=elvsim-elvobs !water height above dam
          maxdamvol=dble(resopalt(1))*dble(rivwth)*dble(rivlen)
          !dam overflow (manning formula) + dam leak
          outsim=(1./man)*(slp**0.5)*dble(rivwth)*(hgt**(5./3))
          outsim=min(max(0.,(volsim-maxdamvol))/dble(dt),outsim) !overflow is limited by the water volume above dam
          outsim=outsim + outflow*outfrc
          outsim=max(0.,min(volsim/dble(dt),outsim)) !total outflow is limited by the total water volume in the river reach
       endif
    else
       write(LIS_logunit,*) '[ERR] Please check reservoir operation type'
       write(LIS_logunit,*) '[ERR] Ending run at HYMAP3_resop_main'
       call LIS_endrun()
    endif
    !update HyMAP variables
    inflow_down=inflow_down-outflow+outsim
    outflow=outsim

  end subroutine HYMAP3_resop_main_glb
  !=============================================
  !=============================================
  subroutine HYMAP3_resop_inter(time, ntvar, tvar, var, mis, varout)

    implicit none

    integer, intent(in)  :: ntvar
    real*8,  intent(in)  :: time
    real,    intent(in)  :: mis
    real*8,  intent(in)  :: tvar(ntvar)
    real*8,  intent(in)  :: var(ntvar)
    real*8,  intent(out) :: varout

    integer, parameter :: &
         indays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    integer            :: iloc(1),ipre,ipos
    real*8             :: ztpre,ztpos,ztime

    if(time<minval(tvar,var/=mis))then
       ipre=0
    else
       iloc(:)=minloc(abs(tvar-time),tvar<=time.and.var/=mis)
       ipre=iloc(1)
    endif

    if(time>maxval(tvar,var/=mis))then
       ipos=0
    else
       iloc(:)=minloc(abs(tvar-time),tvar>time.and.var/=mis)
       ipos=iloc(1)
    endif

    if(ipre==0.and.ipos>0)then
       varout=var(ipos)
    elseif(ipos==0.and.ipre>0)then
       varout=var(ipre)
    elseif(ipre>0.and.ipos>0)then
       if(time==tvar(ipre))then
          varout=var(ipre)
       elseif(time==tvar(ipos))then
          varout=var(ipos)
       else
          call convert_date2fraction_res(tvar(ipre),ztpre)
          call convert_date2fraction_res(tvar(ipos),ztpos)
          call convert_date2fraction_res(time,ztime)
          varout=real(var(ipre))+real(var(ipos)-var(ipre))* &
               (ztime-ztpre)/(ztpos-ztpre)
       endif
    endif

  end subroutine HYMAP3_resop_inter
  !=============================================
  !=============================================
  subroutine convert_date2fraction_res(zdate0, zdate)

    implicit none

    real*8,       intent(in)  :: zdate0
    real*8,       intent(out) :: zdate

    integer,      parameter   :: &
         indays(12)    = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    integer     :: ii,inday0,inday1
    integer     :: iyear,imon,iday

    iyear=int(zdate0/10000d0)
    imon=int(mod(zdate0,10000d0)/100)
    iday=int(mod(zdate0,100d0))
    inday0=365;if(mod(iyear,4)==0)inday0=inday0+1
    inday1=0
    do ii=1,imon-1
       inday1=inday1+indays(ii)
       if(ii==2.and.mod(iyear,4)==0)inday1=inday1+1
    enddo
    inday1=inday1+iday
    zdate=iyear+real(inday1)/inday0+mod(zdate0,1d0)/inday0
  end subroutine convert_date2fraction_res

  !=============================================
  !=============================================
end module HYMAP3_resopMod
