!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT----------------------- 
! !REVISION HISTORY: 
! 27 Sep 2014: Augusto Getirana, Initial implementation
! 30 Jan 2016: Augusto Getirana, Adapt to run in LIS v7.1 
! 13 May 2019: Augusto Getirana, Fixed bug in HYMAP2_get_volume_profile 
!                                
!
module HYMAP2_resopMod
  !module defining reservoir operation variables
  !Augusto Getirana
  !NASA Goddard Space Flight Center
  !27 Sep 2014

  use HYMAP2_modelMod
contains
  !=============================================
  !=============================================  
  subroutine HYMAP2_resop_main_1(mis,nseqall,nz,time,dt,next,rivsto,fldsto,runoff,&
                        rivout,fldout,rivvel,fldvel,&
                        elevtn,fldhgt,fldstomax,grarea,rivstomax,rivelv,rivlen,rivwth,inflow,&
                        ntresop,resopaltloc,tresopalt,resopalt,resopoutmin)
    implicit none
   
    real*8,  intent(in)    :: time                !current time [year]
    real,    intent(in)    :: dt                  !time step length [s]
    real,    intent(in)    :: mis                 !missing data
    integer, intent(in)    :: nseqall             !length of 1D sequnece for river and mouth  
    integer, intent(in)    :: nz                  !number of grids in vertical
    integer, intent(in)    :: next                !point downstream horizontal
    real,    intent(in)    :: rivsto              !floodplain storage [m3]
    real,    intent(in)    :: fldsto              !floodplain storage   [m3]
    real,    intent(in)    :: runoff              !total runoff in grid cell  [m3/s]
    real,    intent(inout) :: rivout              !river outflow  [m3/s]
    real,    intent(inout) :: fldout              !floodplain outflow  [m3/s]
    real,    intent(inout) :: rivvel              !river flow velocity [m/s]
    real,    intent(inout) :: fldvel              !floodplain flow velocity [m/s]
    real,    intent(in)    :: elevtn              !river bed elevation [m]
    real,    intent(in)    :: fldhgt(nz)          !floodplain height
    real,    intent(in)    :: fldstomax(nz)       !maximum floodplain storage [m3]
    real,    intent(in)    :: grarea              !area of the grid [m^2] 
    real,    intent(in)    :: rivstomax           !maximum river storage [m3]
    real,    intent(in)    :: rivelv              !elevation of river bed [m]
    real,    intent(in)    :: rivlen              !river length [m] 
    real,    intent(in)    :: rivwth              !river width [m]
    real,    intent(in)    :: inflow              !inflow from upstream grid cells [m3/s]
    !integer, intent(in)    :: nresop
    integer, intent(in)    :: ntresop
    integer, intent(in)    :: resopaltloc
    real*8,  intent(in)    :: tresopalt(ntresop)
    real,    intent(in)    :: resopalt(ntresop)
    real,    intent(in)    :: resopoutmin
 
    integer                :: itime,iresop,ic,ic_down
    integer                :: elvobs,elvsim
    real*8                 :: volobs,volsim,volsim0
    real*8                 :: outsim

    integer                :: ielevtn
    integer                :: ifldhgt(nz)
    integer                :: irivelv
    integer                :: iresopalt(ntresop)

      ielevtn=int(elevtn*1000.)
      ifldhgt=int(fldhgt*1000.)
      irivelv=int(rivelv*1000.)
      iresopalt=int(mis); where(resopalt/=mis)iresopalt=int(resopalt*1000.)

      !water storage in reservoir considering inflow 
      volsim=rivsto+fldsto+(inflow+runoff)*dble(dt)
      
!!write(12,'(10f18.3)'),time,volsim,rivsto,fldsto,(inflow+runoff)*dble(dt)

!!write(12,'(10f18.3)')fldstomax
!!write(12,'(10f18.3)')real(ielevtn)/1000,real(ifldhgt)/1000,grarea,rivstomax,real(irivelv)/1000,rivlen,rivwth,real(elvobs)/1000.

      !get observed elevation from linear intepolation
      call HYMAP2_resop_inter(time,ntresop,tresopalt,iresopalt,mis,elvobs)

#if 0       
      !calculate water volume as a function of observed elevation
      call HYMAP2_get_volume_profile(nz,ielevtn,ifldhgt,dble(fldstomax),dble(grarea),dble(rivstomax),&
                              irivelv,dble(rivlen),dble(rivwth),elvobs,volobs)
#endif
!write(12,'(i,10f18.3)')elvobs,volobs
      
      !update model variables based on water volume availability
      if(volsim>=volobs+dble(resopoutmin*dble(dt)))then
        outsim=max((volsim-volobs)/dble(dt),resopoutmin)
      else
        outsim=max(0.,min(volsim/dble(dt),resopoutmin))
      endif

      !calculate new water volume and elevation as a function of new simulated outflow
      volsim0=max(0.,volsim-outsim*dble(dt)) 

#if 0
      call HYMAP2_get_elevation_profile(nz,ielevtn,ifldhgt,fldstomax,grarea,rivstomax,&
                            irivelv,rivlen,rivwth,elvsim,volsim0)
#endif                                                       
!write(12,'(i,10f18.3)')elvsim,volsim,volsim0,outsim

      !update HyMAP variables
      if(rivout>=0)then
        rivout=outsim
        rivvel=rivout/dble(dt)
        fldout=0.
        fldvel=0.
      endif

!write(11,'(f18.6,2i,10f18.3)')time,elvsim,elvobs,dt,volsim,volobs,volobs-volsim,outsim,rivout

  end subroutine HYMAP2_resop_main_1
  !=============================================
  !=============================================  
  subroutine HYMAP2_resop_main(mis,nseqall,nz,time,dt,next,rivsto,fldsto,runoff,&
                        rivout,fldout,rivvel,fldvel,&
                        elevtn,fldhgt,fldstomax,grarea,rivstomax,rivelv,rivlen,rivwth,&
                        nresop,ntresop,resopaltloc,tresopalt,resopalt,resopoutmin)
    implicit none
   
    real*8,  intent(in)    :: time         !current time [year]
    real,    intent(in)    :: dt           !time step length [s]
    real,    intent(in)    :: mis          !missing data
    integer, intent(in)    :: nseqall               !length of 1D sequnece for river and mouth  
    integer, intent(in)    :: nz           !number of grids in vertical
    integer, intent(in)    :: next(nseqall)         !point downstream horizontal
    real,    intent(in)    :: rivsto(nseqall)     !floodplain storage [m3]
    real,    intent(in)    :: fldsto(nseqall)     !floodplain storage   [m3]
    real,    intent(in)    :: runoff(nseqall)     !total runoff in grid cell  [m3/s]
    real,    intent(inout) :: rivout(nseqall)     !river outflow  [m3/s]
    real,    intent(inout) :: fldout(nseqall)     !floodplain outflow  [m3/s]
    real,    intent(inout) :: rivvel(nseqall)     !river flow velocity [m/s]
    real,    intent(inout) :: fldvel(nseqall)     !floodplain flow velocity [m/s]
    real,    intent(in)    :: elevtn(nseqall)     !river bed elevation [m]
    real,    intent(in)    :: fldhgt(nseqall,nz)  !floodplain height
    real,    intent(in)    :: fldstomax(nseqall,nz) !maximum floodplain storage [m3]
    real,    intent(in)    :: grarea(nseqall)     !area of the grid [m^2] 
    real,    intent(in)    :: rivstomax(nseqall)  !maximum river storage [m3]
    real,    intent(in)    :: rivelv(nseqall)     !elevation of river bed [m]
    real,    intent(in)    :: rivlen(nseqall)     !river length [m] 
    real,    intent(in)    :: rivwth(nseqall)     !river width [m]
    integer, intent(in)    :: nresop
    integer, intent(in)    :: ntresop
    integer, intent(in)    :: resopaltloc(nresop)
    real*8,  intent(in)    :: tresopalt(nresop,ntresop)
    real,    intent(in)    :: resopalt(nresop,ntresop)
    real,    intent(in)    :: resopoutmin(nresop)
 
    integer                :: itime,iresop,ic,ic_down
    integer                :: elvobs,elvsim
    real*8                 :: volobs,volsim
    real*8                 :: inflow,outsim

    integer                :: ielevtn
    integer                :: ifldhgt(nz)
    integer                :: irivelv
    integer                :: iresopalt(ntresop)

    do iresop=1,nresop
      ic=resopaltloc(iresop)
      ic_down=next(ic)
      inflow=sum(rivout+fldout,next==ic)

      ielevtn=int(elevtn(ic)*1000.)
      ifldhgt=int(fldhgt(ic,:)*1000.)
      irivelv=int(rivelv(ic)*1000.)
      iresopalt=int(mis); where(resopalt(iresop,:)/=mis)iresopalt=int(resopalt(iresop,:)*1000.)

      !water storage in reservoir considering inflow 
      volsim=rivsto(ic)+fldsto(ic)+(inflow+runoff(ic))*dble(dt)
      
      !get observed elevation from linear intepolation
      call HYMAP2_resop_inter(time,ntresop,tresopalt(iresop,:),iresopalt,mis,elvobs)
#if 0       
      !calculate water volume as a function of observed elevation
      call HYMAP2_get_volume_profile(nz,ielevtn,ifldhgt,dble(fldstomax(ic,:)),dble(grarea(ic)),dble(rivstomax(ic)),&
                              irivelv,dble(rivlen(ic)),dble(rivwth(ic)),elvobs,volobs)

#endif
      !update model variables based on water volume availability
      if(volsim>=volobs+dble(resopoutmin(iresop)*dble(dt)))then
        outsim=max((volsim-volobs)/dble(dt),resopoutmin(iresop))
      else
        outsim=max(0.,min(volsim/dble(dt),resopoutmin(iresop)))
      endif

      !calculate new water volume and elevation as a function of new simulated outflow
      volsim=max(0.,volsim-outsim*dble(dt)) 
#if 0 
      call HYMAP2_get_elevation_profile(nz,ielevtn,ifldhgt,fldstomax(ic,:),grarea(ic),rivstomax(ic),&
                            irivelv,rivlen(ic),rivwth(ic),elvsim,volsim)
#endif
      !update HyMAP variables
      rivout(ic)=outsim
      rivvel(ic)=rivout(ic)/dble(dt)
      fldout(ic)=0.
      fldvel(ic)=0.

    enddo
  end subroutine HYMAP2_resop_main
  !=============================================

  !=============================================
  !=============================================  
  subroutine HYMAP2_resop_inter(time,ntvar,tvar,var,mis,varout)
    implicit none
    
    integer, intent(in)  :: ntvar
    real*8,  intent(in)  :: time
    real,    intent(in)  :: mis
    real*8,  intent(in)  :: tvar(ntvar)
    integer, intent(in)  :: var(ntvar)
    integer, intent(out) :: varout
    
    integer, parameter :: indays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)  
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
      call convert_date2fraction_res(tvar(ipre),ztpre) 
      call convert_date2fraction_res(tvar(ipos),ztpos) 
      call convert_date2fraction_res(time,ztime) 
      varout=int(real(var(ipre))+real(var(ipos)-var(ipre))*(ztime-ztpre)/(ztpos-ztpre))
    endif

   end subroutine HYMAP2_resop_inter
  !=============================================
  !=============================================  
  subroutine convert_date2fraction_res(zdate0,zdate)
    implicit none
    real*8,       intent(in)  :: zdate0
    real*8,       intent(out) :: zdate
    integer,      parameter   :: indays(12)    = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    integer     :: ii,inday0,inday1
    integer     :: iyear,imon,iday
    
    iyear=int(zdate0/10000d0)
    imon=int(mod(zdate0,10000d0)/100)
    iday=int(mod(zdate0,100d0))
    inday0=365;if(mod(iyear,4)==0)inday0=inday0+1
    inday1=0
    do ii=1,imon-1
      inday1=inday1+indays(ii)
    enddo
    inday1=inday1+iday
    zdate=iyear+real(inday1)/inday0+mod(zdate0,1d0)/inday0
  end subroutine convert_date2fraction_res
  !=============================================
  !=============================================  
end module HYMAP2_resopMod
