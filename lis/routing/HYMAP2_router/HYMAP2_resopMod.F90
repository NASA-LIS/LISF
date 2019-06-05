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
!
module HYMAP2_resopMod
  !module defining reservoir operation variables
  !Augusto Getirana
  !NASA Goddard Space Flight Center
  !27 Sep 2014

contains
  !=============================================
  !=============================================  
  subroutine HYMAP2_resop_run(n,mis,nseqall,nz,time,dt,next,rivsto,fldsto,rivout,fldout,rivvel,fldvel,&
                       elevtn,fldhgt,fldstomax,grarea,rivstomax,rivelv,rivlen,rivwth)

    use HYMAP2_routingMod, only : HYMAP2_routing_struc
    
    implicit none
    
    real*8,  intent(in)    :: time                !current time [year]
    real,    intent(in)    :: dt                  !time step length [s]
    real,    intent(in)    :: mis                 !missing data
    integer, intent(in)    :: n
    !integer, intent(in)    :: nx                  !number of grids in horizontal
    !integer, intent(in)    :: ny                  !number of grids in vertical
    integer, intent(in)    :: nseqall               !length of 1D sequnece for river and mouth  
    integer, intent(in)    :: nz                  !number of stages in the sub-grid discretization
    integer, intent(in)    :: next(nseqall)         !point downstream horizontal
    !integer, intent(in)    :: nextx(nseqall)        !point downstream horizontal
    !integer, intent(in)    :: nexty(nseqall)        !point downstream vertical
    real,    intent(in)    :: rivsto(nseqall)       !floodplain storage [m3]
    real,    intent(in)    :: fldsto(nseqall)       !floodplain storage   [m3]
    real,    intent(inout) :: rivout(nseqall)       !river outflow  [m3/s]
    real,    intent(inout) :: fldout(nseqall)       !floodplain outflow  [m3/s]
    real,    intent(inout) :: rivvel(nseqall)       !river flow velocity [m/s]
    real,    intent(inout) :: fldvel(nseqall)       !floodplain flow velocity [m/s]
    !real,    intent(inout) :: rivinf(nseqall)       !river inflow   [m3/s]
    !real,    intent(inout) :: fldinf(nseqall)
    real,    intent(in)    :: elevtn(nseqall)       !river bed elevation [m]
    real,    intent(in)    :: fldhgt(nseqall,nz)    !floodplain height
    real,    intent(in)    :: fldstomax(nseqall,nz) !maximum floodplain storage [m3]
    real,    intent(in)    :: grarea(nseqall)       !area of the grid [m^2] 
    real,    intent(in)    :: rivstomax(nseqall)    !maximum river storage [m3]
    real,    intent(in)    :: rivelv(nseqall)       !elevation of river bed [m]
    real,    intent(in)    :: rivlen(nseqall)       !river length [m] 
    real,    intent(in)    :: rivwth(nseqall)       !river width [m]

    call HYMAP2_resop_main(mis,nseqall,nz,time,dt,next,rivsto,fldsto,&
                    rivout,fldout,rivvel,fldvel,&
                    elevtn,fldhgt,fldstomax,grarea,rivstomax,rivelv,rivlen,rivwth,&
                    HYMAP2_routing_struc(n)%nresop,HYMAP2_routing_struc(n)%ntresop,&
                    HYMAP2_routing_struc(n)%resopaltloc,&
                    HYMAP2_routing_struc(n)%tresopalt,HYMAP2_routing_struc(n)%resopalt,&
                    HYMAP2_routing_struc(n)%resopoutmin)
         
  end subroutine HYMAP2_resop_run
  !=============================================
  !=============================================  
  subroutine HYMAP2_resop_main(mis,nseqall,nz,time,dt,next,rivsto,fldsto,&
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
    real                   :: elvobs,volobs,elvsim,volsim,outsim
    real                   :: inflow

    do iresop=1,nresop
      ic=resopaltloc(iresop)
      ic_down=next(ic)
      inflow=sum(rivout+fldout,next==ic)
      
      !water storage in reservoir considering inflow 
      volsim=rivsto(ic)+fldsto(ic)+inflow*dt
      
      !get observed elevation from linear intepolation
      call HYMAP2_resop_inter(time,ntresop,tresopalt(iresop,:),resopalt(iresop,:),mis,elvobs)
      
      !calculate water volume as a function of observed elevation
      call HYMAP2_get_volume_profile(nz,elevtn(ic),fldhgt(ic,:),fldstomax(ic,:),grarea(ic),rivstomax(ic),&
                              rivelv(ic),rivlen(ic),rivwth(ic),elvobs,volobs)

      !update model variables based on water volume availability
      if(volsim>=volobs+resopoutmin(iresop)*dt)then
        outsim=max((volsim-volobs)/dt,resopoutmin(iresop))
      else
        outsim=max(0.,min(volsim/dt,resopoutmin(iresop)))
      endif
      !calculate new water volume and elevation as a function of new simulated outflow
      volsim=max(0.,volsim-outsim*dt) 
      call HYMAP2_get_elevation_profile(nz,elevtn(ic),fldhgt(ic,:),fldstomax(ic,:),grarea(ic),rivstomax(ic),&
                            rivelv(ic),rivlen(ic),rivwth(ic),elvsim,volsim)
      !update HyMAP variables
      rivout(ic)=outsim
      rivvel(ic)=rivout(ic)/dt
      fldout(ic)=0.
      fldvel(ic)=0.
    enddo
  end subroutine HYMAP2_resop_main
  !=============================================
  !=============================================  
  subroutine HYMAP2_get_volume_profile(nz,elevtn,fldhgt,fldstomax,grarea,rivstomax,rivelv,rivlen,rivwth,elv,vol)
   
    use HYMAP2_routingMod, only : HYMAP2_routing_struc
    use LIS_logMod,        only : LIS_logunit, LIS_endrun
    
    implicit none
   
    integer, intent(in)  :: nz
    real,    intent(in)  :: elevtn,fldhgt(nz)
    real,    intent(in)  :: fldstomax(nz),grarea,rivstomax,rivelv,rivlen,rivwth
    real,    intent(in)  :: elv
    real,    intent(out) :: vol
    integer :: i
    real    :: h1,h2,v1,v2
    real    :: dph(nz)!,zwth(nz)
    real    :: dphtmp!,zwthtmp
    
    dph(:)=elevtn+fldhgt(:)
    
    !if(elv>flddph(ix,iy,1))then
    if(elv>elevtn)then
      i=1
      do while(elv>dph(i).and.i<nz)
        i=i+1
      enddo
      if(i>=nz)then
        vol=fldstomax(nz)+(elv-dph(nz))*grarea
      else
        if(i>1)then
          h1=dph(i-1);v1=fldstomax(i-1)
          h2=dph(i);v2=fldstomax(i)
        elseif(i==1)then
          h1=elevtn;v1=rivstomax
          h2=dph(1);v2=fldstomax(1)
        else
          write(LIS_logunit,*) '[ERR] HYMAP2_get_volume_profile: ' // &
                               'Please check Reservoir elevation'
          call LIS_endrun
        endif
        vol=v1+(v2-v1)*(elv-h1)/(h2-h1)
      endif   
    else
      vol=max(0.,elv-rivelv)*rivlen*rivwth
    endif
  end subroutine HYMAP2_get_volume_profile 
  !=============================================
  !=============================================  
  subroutine HYMAP2_get_elevation_profile(nz,elevtn,fldhgt,fldstomax,grarea,rivstomax,rivelv,rivlen,rivwth,elv,vol)
   
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
   
    integer, intent(in)  :: nz
    real,    intent(in)  :: elevtn,fldhgt(nz)
    real,    intent(in)  :: fldstomax(nz),grarea,rivstomax,rivelv,rivlen,rivwth
    real,    intent(in)  :: vol
    real,    intent(out) :: elv
    integer :: i
    real    :: h1,h2,v1,v2
    real   :: dph(nz)!,zwth(nz)
    real   :: dphtmp!,zwthtmp
    
    dph(:)=elevtn+fldhgt(:)
    
    if(vol>rivstomax)then
      i=1
      do while(vol>fldstomax(i).and.i<nz)
        i=i+1
      enddo
      if(i>=nz)then
        elv=dph(nz)+(vol-fldstomax(nz))/grarea
      else
        if(i>1)then
          h1=dph(i-1);v1=fldstomax(i-1)
          h2=dph(i);v2=fldstomax(i)
        elseif(i==1)then
          h1=elevtn;v1=rivstomax
          h2=dph(1);v2=fldstomax(1)
        else
          write(LIS_logunit,*) '[ERR] HYMAP2_get_elevation_profile: ' // &
                               'Please check Reservoir elevation'
          call LIS_endrun
        endif
        elv=h1+(h2-h1)*(vol-v1)/(v2-v1)
      endif   
    else
      elv=rivelv+vol/rivlen/rivwth
    endif
  end subroutine HYMAP2_get_elevation_profile
  !=============================================
  !=============================================  
  subroutine HYMAP2_resop_inter(time,ntvar,tvar,var,mis,varout)
    implicit none
    
    integer, intent(in)  :: ntvar
    real*8,  intent(in)  :: time
    real,    intent(in)  :: mis
    real*8,  intent(in)  :: tvar(ntvar)
    real,    intent(in)  :: var(ntvar)
    real,    intent(out) :: varout
    
    integer, parameter :: indays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)  
    integer :: iloc(1),ipre,ipos
   
    iloc(:)=minloc(abs(tvar-time),tvar<=time.and.var/=mis)
    ipre=iloc(1)    
    iloc(:)=minloc(abs(tvar-time),tvar>time.and.var/=mis)
    ipos=iloc(1)

    if(ipre==0.and.ipos>0)then
      varout=var(ipos)    
    elseif(ipos==0.and.ipre>0)then
      varout=var(ipre)
    elseif(ipre>0.and.ipos>0)then
      varout=var(ipre)+(var(ipos)-var(ipre))*(time-tvar(ipre))/(tvar(ipos)-tvar(ipre))
    endif
   end subroutine HYMAP2_resop_inter
 
end module HYMAP2_resopMod
