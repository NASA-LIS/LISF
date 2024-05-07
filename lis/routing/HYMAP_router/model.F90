!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine model(rmis,inx,iny,inz,idt,iflowtype,ilinres,ievapflag,&
     i1seqx,i1seqy,inseqriv,inseqall, &
     i2nextx,i2nexty,r2elevtn,r2nxtdst,r2grarea,    &
     r2fldgrd,r2fldman,r2fldstomax,                 &
     r2rivman,r2rivelv,r2rivstomax,                 &
     r2rivlen,r2rivwth,r2rivhgt,r2rivare,rslpmin,   &
     r2trnoff,r2tbsflw,r2cntime,r2mask,&
     r2runoff0,r2basflw0,&
     !==== TEMPORARILY COMMENTED ====
     !                 !input evap				 
     !                 r2tsfc,r2psur,r2pres,r2tair,r2relh,          &
     !                 r2wind,r2zref,                               &
     
     !outputs
     r2rivsto,r2rivdph,r2rivvel,r2rivout,r2evpout,&
     r2fldout,r2fldsto,r2flddph,r2fldvel,r2fldfrc,&
     r2fldare,r2sfcelv,r2roffsto,r2basfsto)
  ! ================================================
  ! to   HyMAP core 
  ! by   Augusto Getirana 
  ! on   6th Dec 2010
  ! at   LEGOS/CNES, Toulouse, France
  ! Adapted for flow routing implementation in LIS 09 Nov 2011
  ! ================================================
  
  use HYMAP_modelMod
  use LIS_logMod,     only : LIS_logunit
  
  implicit none

! === Resolution ========================================
  integer, intent(in)                           :: idt           !! time step length [s]
  integer, intent(in)                           :: inx           !! number of grids in horizontal
  integer, intent(in)                           :: iny           !! number of grids in vertical
  integer, intent(in)                           :: inz             !! number of stages in the sub-grid discretization
! === Flags =============================================
  integer, intent(in)                           :: iflowtype     !! routing method: 1 - kinematic; 2 - diffusive
  integer, intent(in)                           :: ilinres       !! linear reservoir flag: 1 - use linear reservoirs; or 2 - do not use
  integer, intent(in)                           :: ievapflag     !! evaporation flag: 1 - do not compute ; or 2 - compute evapotation in floodplains
! === Undefined Values ==================================
  real,    intent(in)                           :: rmis         !! real undefined value
! === River sequence ====================================
  integer, intent(in),    dimension(inx*iny)    :: i1seqx      !! 1D sequence horizontal
  integer, intent(in),    dimension(inx*iny)    :: i1seqy      !! 1D sequence vertical
  integer, intent(in)                           :: inseqriv    !! length of 1D sequnece for river
  integer, intent(in)                           :: inseqall    !! length of 1D sequnece for river and mouth  
! === Map ===============================================
  integer, intent(in),    dimension(inx,iny)    :: i2nextx      !! point downstream horizontal
  integer, intent(in),    dimension(inx,iny)    :: i2nexty      !! point downstream vertical
  real,    intent(in),    dimension(inx,iny)    :: r2elevtn     !! river bed elevation [m]
  real,    intent(in),    dimension(inx,iny)    :: r2nxtdst     !! distance to the next grid [m]
  real,    intent(in),    dimension(inx,iny)    :: r2grarea     !! area of the grid [m^2] 
  real,    intent(in),    dimension(inx,iny,inz) :: r2fldgrd    !! floodplain gradient [-]
  real,    intent(in),    dimension(inx,iny)    :: r2fldman     !! maning coefficient for floodplains [-]
  real,    intent(in),    dimension(inx,iny,inz) :: r2fldstomax !! maximum floodplain storage [m3]
  real,    intent(in),    dimension(inx,iny)    :: r2rivman     !! maning coefficient for rivers [-]
  real,    intent(in),    dimension(inx,iny)    :: r2rivelv     !! elevation of river bed [m]
  real,    intent(in),    dimension(inx,iny)    :: r2rivstomax  !! maximum river storage [m3]
  real,    intent(in),    dimension(inx,iny)    :: r2rivlen     !! river length [m] 
  real,    intent(in),    dimension(inx,iny)    :: r2rivwth     !! river wiidth [m]
  real,    intent(in),    dimension(inx,iny)    :: r2rivhgt     !! river heihgt [m]
  real,    intent(in),    dimension(inx,iny)    :: r2rivare     !! river surface area [m2]
  real,    intent(in),    dimension(inx,iny)    :: r2trnoff     !! runoff   concentration time parameter [day]
  real,    intent(in),    dimension(inx,iny)    :: r2tbsflw     !! baseflow concentration time parameter [day]
  real,    intent(in),    dimension(inx,iny)    :: r2cntime     !! concentration time [s]
  real,    intent(in)                           :: rslpmin      !! minimum river slope [-]
  real,    intent(inout), dimension(inx,iny)    :: r2runoff0    !! input runoff [mm.idt-1]
  real,    intent(inout), dimension(inx,iny)    :: r2basflw0    !! input baseflow [mm.idt-1]
  integer, intent(in),    dimension(inx,iny)    :: r2mask

 ! === InOut ============================================
  real,    intent(inout), dimension(inx,iny)    :: r2roffsto    !! runoff   reservoir storage [m3]
  real,    intent(inout), dimension(inx,iny)    :: r2basfsto    !! baseflow reservoir storage [m3]
  real,    intent(inout), dimension(inx,iny)    :: r2rivsto     !! floodplain storage [m3]
  real,    intent(inout), dimension(inx,iny)    :: r2fldsto     !! floodplain storage   [m3]
 ! === Outputs ==========================================
  real,    intent(out),   dimension(inx,iny)    :: r2rivdph     !! river depth [m]
  real,    intent(out),   dimension(inx,iny)    :: r2rivvel     !! river flow velocity [m/s]
  real,    intent(out),   dimension(inx,iny)    :: r2rivout     !! river outflow  [m3/s]
  real,    intent(inout), dimension(inx,iny)    :: r2evpout     !! effective evaporation from floodplains [m3]
  real,    intent(out),   dimension(inx,iny)    :: r2fldout     !! floodplain outflow  [m3/s]
  real,    intent(out),   dimension(inx,iny)    :: r2flddph     !! floodplain depth [m]
  real,    intent(out),   dimension(inx,iny)    :: r2fldvel     !! floodplain flow velocity [m/s]
  real,    intent(out),   dimension(inx,iny)    :: r2fldfrc     !! area fraction [-]
  real,    intent(out),   dimension(inx,iny)    :: r2fldare     !! flooded area [m2]
  real,    intent(out),   dimension(inx,iny)    :: r2sfcelv     !! water surface elevation [m]  
  ! === Local variables ============================
  real,                   dimension(inx,iny)    :: r2flddph1
  real,                   dimension(inx,iny)    :: r2fldinf
  real,                   dimension(inx,iny)    :: r2fldwth
  real,                   dimension(inx,iny)    :: r2runoff     !! runoff+baseflow after linear reservoirs [m3/s]
  real,                   dimension(inx,iny)    :: r2etrans
  real,                   dimension(inx,iny)    :: r2evpdif
  real,                   dimension(inx,iny)    :: r2evpwat
  real,                   dimension(inx,iny)    :: r2rivinf     !! river inflow   [m3/s]
  real                                          :: rglbfldare

integer :: ii,x,y
!integer, parameter :: LIS_logunit = 6

! ================================================
  !CONVERT UNITS [kg m-2 sec-1] to [m3 sec-1]
  where(r2runoff0.ne.rmis)
     r2runoff0=r2runoff0*r2grarea/1e3
     r2basflw0=r2basflw0*r2grarea/1e3
  else where
     r2runoff0=0.
     r2basflw0=0.
  endwhere
! ================================================
! Linear reservoirs
! Store runoff and baseflow in linear reservoirs
! Concentration times are parameters provided as input files 
! ================================================
  if(ilinres==1)then
     call linear_reservoir_lis(inx,iny,i1seqx,i1seqy,inseqall,idt,    &
                           rmis,r2runoff0,r2basflw0,r2trnoff,r2tbsflw,&
                           r2cntime,r2roffsto,r2basfsto,r2runoff      )
  elseif(ilinres==2)then
     call no_reservoir_lis(inx,iny,idt,rmis,r2runoff0,r2basflw0,&
                       r2roffsto,r2basfsto,r2runoff)      
  elseif(ilinres==3)then
     !runs linear reservoir at the daily time step in the model.f90  routine
  else
     write(LIS_logunit,*)"HYMAP routing model linear reservoir flag: unknown value"
     stop
  endif
! ================================================
!
! ================================================
!==== TEMPORARILY COMMENTED =====
!! read meteorological forcings
!  call get_evap(cevapflag,inx,iny,inxin,inyin,i1seqx,i1seqy,inseqall,rmis,iyear,imon,iday,csufinp,r2evpdif)
! ================================================
!
! ================================================
!==== TEMPORARILY COMMENTED =====
!! calculate evaporation from floodplains
!  call calc_evap_fld(inx,iny,i1seqx,i1seqy,inseqall,idt,idt,&
!       r2fldare,r2rivare,r2fldsto,r2rivsto,            &
!       r2evpdif,r2evpout                               )
! ================================================
!
! ================================================
! Calculate flood and river storage
  call calc_fldstg(inx,iny,i1seqx,i1seqy,inseqall,        &
       r2grarea,r2rivlen,r2rivwth,r2rivinf,               &
       r2rivstomax,r2fldstomax,r2fldgrd,                  &
       r2rivsto,r2fldsto,r2rivdph,r2flddph,r2flddph1,     &
       r2fldwth,r2fldfrc,r2fldare,rglbfldare,             &
       i2nextx,i2nexty,inseqriv,r2rivelv,r2nxtdst,r2rivare,r2mask)
       
! ================================================
!
! ================================================
  if(iflowtype==1)then
     !Calculate river flow based on the kinematic wave equation
     call calc_rivout_kine(1,inx,iny,idt,i1seqx,i1seqy,inseqriv,inseqall,&
          i2nextx,i2nexty,r2rivelv,r2elevtn,r2nxtdst,r2rivwth,    &
          r2rivsto,r2rivdph,                                      &
          r2rivinf,r2rivout,r2rivvel,r2sfcelv,                    &
          r2rivlen,r2rivman,rslpmin                               )
  elseif(iflowtype==2)then
     !Calculate river flow based on the diffusive wave equation
     call calc_rivout_diff(1,inx,iny,idt,i1seqx,i1seqy,inseqriv,inseqall,&
          i2nextx,i2nexty,r2rivelv,r2elevtn,r2nxtdst,r2rivwth,    &
          r2rivsto,r2rivdph,                                      &
          r2rivinf,r2rivout,r2rivvel,r2sfcelv,                    &
          r2rivlen,r2rivman,rslpmin                               )
  else
     write(LIS_logunit,*)"HYMAP routing method: unknown value"
     stop
  endif
! ================================================
!
! ================================================
! Calculate floodplain flow based on the kinematic wave equation
  call calc_rivout_kine(2,inx,iny,idt,i1seqx,i1seqy,inseqriv,inseqall,&
       i2nextx,i2nexty,r2rivelv,r2elevtn,r2nxtdst,r2fldwth,    &
       r2fldsto,r2flddph1,                                     &
       r2fldinf,r2fldout,r2fldvel,r2sfcelv,                    &
       r2rivlen,r2fldman,rslpmin                               )
! ================================================
!
! ================================================
! Update water storage in river reservoir
     call calc_stonxt(inx,iny,i1seqx,i1seqy,inseqall,&
          idt,r2rivinf,r2rivout,r2fldinf,r2fldout,r2rivsto,r2fldsto )
! ================================================
!
! ================================================
! Calculate runoff in the river network for the current time step
  call calc_runoff(inx,iny,i1seqx,i1seqy,inseqall,idt,&
       r2fldfrc,r2runoff,                             &
       r2rivsto,r2fldsto,r2rivinf,                    &
       i2nextx,i2nexty                                )
! ================================================

return
end subroutine model
