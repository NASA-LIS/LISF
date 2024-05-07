!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module HYMAP_modelMod
  
  implicit none
  
contains
  
  subroutine linear_reservoir_lis(nx,ny,i1seqx,i1seqy,nseqall,dt,rmis,r2runoff0,r2basflw0, &
                              r2trnoff,r2tbsflw,r2cntime,r2roffsto,&
                              r2basfsto,r2runoff)
    ! ================================================
    ! to   This routine simulates surface and underground linear reservoirs  
    ! by   Augusto GETIRANA 
    ! on   6th Dec 2010
    ! at   LEGOS/CNES, Toulouse, France
    ! Adapted for flow routing implementation in LIS 09 Nov 2011
    ! ================================================
    
    implicit none
    
    ! ================================================
    !Dummy variables
    integer, intent(in)                   :: nx,ny       !! grid size
    integer, intent(in), dimension(nx*ny) :: i1seqx      !! 1D sequence horizontal
    integer, intent(in), dimension(nx*ny) :: i1seqy      !! 1D sequence vertical
    integer, intent(in)                   :: nseqall     !! length of 1D sequnece for river and mouth
    integer, intent(in)                   :: dt          !! input time step length [s]
    real,  intent(in)                     :: rmis        !! missing data
    real,  intent(in),   dimension(nx,ny) :: r2runoff0   !! runoff   before linear reservoir [m3/s]
    real,  intent(in),   dimension(nx,ny) :: r2basflw0   !! baseflow before linear reservoir [m3/s]
    real,  intent(in),   dimension(nx,ny) :: r2trnoff    !! runoff   concentration time parameter [day]
    real,  intent(in),   dimension(nx,ny) :: r2tbsflw    !! baseflow concentration time parameter [day]
    real,  intent(in),   dimension(nx,ny) :: r2cntime    !! concentration time [s]    
    real,  intent(inout),dimension(nx,ny) :: r2roffsto   !! runoff   reservoir storage [m3]
    real,  intent(inout),dimension(nx,ny) :: r2basfsto   !! baseflow reservoir storage [m3]
    real,  intent(out),  dimension(nx,ny) :: r2runoff    !! runoff+baseflow after  linear reservoir [m3/s]
    ! ================================================
    !Local variables
    real,           dimension(nx,ny) :: vroff, vbflw
    real,           dimension(nx,ny) :: qrun, qbas
  
    integer :: i,j,jj
    integer :: ix, iy, iseq
    ! ================================================
    !   
    ! ================================================
    vroff=0.0
    vbflw=0.0
    qrun=0.0
    qbas=0.0
    !
    do iseq=1, nseqall
       ix=i1seqx(iseq)
       iy=i1seqy(iseq)
       !runoff
       !update reservoir volumes with input data (runoff)
       r2roffsto(ix,iy) = r2roffsto(ix,iy) + r2runoff0(ix,iy) * real(dt)
       !compute temporary volumes (m3)
       qrun(ix,iy)=r2roffsto(ix,iy)/(r2trnoff(ix,iy)*r2cntime(ix,iy))  
       vroff(ix,iy)=r2roffsto(ix,iy) - qrun(ix,iy) * real(dt)
       !check if reservoir is empty
       if(vroff(ix,iy)<0.)then
          qrun(ix,iy)=r2roffsto(ix,iy)/real(dt)
          r2roffsto(ix,iy)=0.
       else
          r2roffsto(ix,iy)=vroff(ix,iy)
       endif
       !baseflow
       !update reservoir volumes with input data (baseflow)
       r2basfsto(ix,iy) = r2basfsto(ix,iy) + r2basflw0(ix,iy) * real(dt)
       !compute temporary volumes (m3)
       qbas(ix,iy)=r2basfsto(ix,iy)/(r2tbsflw(ix,iy)*86400.)
       vbflw(ix,iy)=r2basfsto(ix,iy)-qbas(ix,iy) * real(dt)
       !check if reservoir is empty
       if(vbflw(ix,iy)<0.)then
          qbas(ix,iy)=r2basfsto(ix,iy)/real(dt)
          r2basfsto(ix,iy)=0.
       else
          r2basfsto(ix,iy)=vbflw(ix,iy)
       endif
       r2runoff(ix,iy)=qrun(ix,iy)+qbas(ix,iy)
    enddo
    
    return
  endsubroutine linear_reservoir_lis
  ! ================================================
  ! ================================================  
  subroutine no_reservoir_lis(nx,ny,dt,rmis,r2runoff0,r2basflw0,&
                          r2roffsto,r2basfsto,r2runoff)
    ! ================================================
    !      This routine sums runoff and baseflow variables
    ! by   Augusto GETIRANA 
    ! on   28th Dec 2010
    ! at   LEGOS/CNES, Toulouse, France
    ! ================================================

    implicit none

    ! ================================================
    integer,                intent(in)  :: nx,ny       !grid size
    integer,                intent(in)  :: dt          !input time step length [s]
    real,                   intent(in)  :: rmis        !missing data
    real, dimension(nx,ny), intent(in)  :: r2runoff0   !runoff          before linear reservoir [m3/s]
    real, dimension(nx,ny), intent(in)  :: r2basflw0   !baseflow        before linear reservoir [m3/s]
    real, dimension(nx,ny), intent(out) :: r2roffsto   !! runoff   reservoir storage [m3]
    real, dimension(nx,ny), intent(out) :: r2basfsto   !! baseflow reservoir storage [m3]
    real, dimension(nx,ny), intent(out) :: r2runoff    !! runoff+baseflow after  linear reservoir [m3/s]
    ! ================================================
    !   
    ! ================================================
    where(r2runoff0/=rmis.and.r2basflw0/=rmis)
       r2runoff=r2runoff0+r2basflw0
    elsewhere 
       r2runoff=0.d0
    end where
    r2roffsto=0.d0
    r2basfsto=0.d0



  endsubroutine no_reservoir_lis
  ! ================================================  
  ! ================================================  
  subroutine calc_runoff(nx,ny,i1seqx,i1seqy,nseqall,dt,&
                         r2fldfrc,r2runoff,             &
                         r2rivsto,r2fldsto,r2rivinf,    &
                         i2nextx,i2nexty)
    ! ================================================
    ! to   calculate runoff to storage
    ! by   Augusto Getirana after Dai YAMAZAKI
    ! on   10 Feb 2011
    ! at   LEGOS/CNES, Toulouse, France
    ! ================================================
    implicit none
    !Dummy variables
    integer, intent(in)                      ::  nx                  !! number of grids in horizontal
    integer, intent(in)                      ::  ny                  !! number of grids in vertical
    integer, intent(in),    dimension(nx*ny) ::  i1seqx       !! 1D sequence horizontal
    integer, intent(in),    dimension(nx*ny) ::  i1seqy       !! 1D sequence vertical
    integer, intent(in)                      ::  nseqall             !! length of 1D sequnece for river and mouth
    integer, intent(in)                      ::  dt
    real,    intent(in),    dimension(nx,ny) ::  r2fldfrc     !! area fraction
    real,    intent(in),    dimension(nx,ny) ::  r2runoff     !! runoff of the grid       [m3/s]
    real,    intent(inout), dimension(nx,ny) ::  r2rivsto     !! river storage        [m3]
    real,    intent(inout), dimension(nx,ny) ::  r2fldsto     !! floodplain storage   [m3]
    !the three variables below are used in a similar sub-routine adapted for the Muskingum-Cunge method
    integer, intent(in),    dimension(nx,ny) ::  i2nextx  !! point downstream horizontal
    integer, intent(in),    dimension(nx,ny) ::  i2nexty  !! point downstream vertical
    real,    intent(inout), dimension(nx,ny) :: r2rivinf        !! total inflow to the grid              [m3/s]
    !Local variables
    integer ::  ix, iy, iseq
    real    ::  rrivrof
    real    ::  rfldrof
    
    ! flow calculation
    do iseq=1, nseqall
       ix=i1seqx(iseq)
       iy=i1seqy(iseq)
       rrivrof = r2runoff(ix,iy) * (1.d0-r2fldfrc(ix,iy)) * real(dt)
       rfldrof = r2runoff(ix,iy) *       r2fldfrc(ix,iy)  * real(dt)
       r2rivsto(ix,iy) = r2rivsto(ix,iy) + rrivrof
       r2fldsto(ix,iy) = r2fldsto(ix,iy) + rfldrof
    end do
    return
  end subroutine calc_runoff
  ! ================================================
  ! ================================================
  subroutine calc_fldstg(nx,ny,i1seqx,i1seqy,nseqall,                      &
                         r2areamat,r2rivlen,r2rivwth,r2rivinf,             &
                         r2rivstomax,r2fldstomax,r2fldgrd,                 &
                         r2rivsto,r2fldsto,r2rivdph,r2flddph,r2flddph1,    &
                         r2fldwth,r2fldfrc,r2fldare,rglbfldarea,           &
                         i2nextx,i2nexty,nseqriv,r2rivelv,r2nxtdst,r2rivare,r2mask)
    ! ================================================
    ! to   calculate river and floodplain staging
    ! by   Augusto Getirana after Dai YAMAZAKI
    ! on   27 Oct 2010
    ! at   LEGOS/CNES, Toulouse, France
    ! ================================================
    implicit none
    ! index
    integer             ::  nx                  !! number of grids in horizontal
    integer             ::  ny                  !! number of grids in vertical
    integer             ::  i1seqx(nx*ny)       !! 1D sequence horizontal
    integer             ::  i1seqy(nx*ny)       !! 1D sequence vertical
    integer             ::  nseqall             !! length of 1D sequnece for river and mouth
    ! boudary

    integer    ::  r2mask(nx,ny)

    real    ::  r2areamat(nx,ny)    !! area of the grid [m2]
    real    ::  r2rivlen(nx,ny)     !! channel length [m]
    real    ::  r2rivwth(nx,ny)     !! river width [m]
    real    ::  r2rivstomax(nx,ny)
    real    ::  r2fldstomax(nx,ny,10)  !! maximum floodplain storage [m3]
    real    ::  r2fldgrd(nx,ny,10)  !! floodplain gradient [-]
    ! modify
    real    ::  r2rivsto(nx,ny)     !! river storage [m3]
    real    ::  r2fldsto(nx,ny)     !! flood plain storage [m3]
    real    ::  r2rivdph(nx,ny)     !! river depth [m]
    real    ::  r2flddph(nx,ny)     !! floodplain depth [m]
    real    ::  r2flddph1(nx,ny)    !! mean floodplain depth [m]
    real    ::  r2fldwth(nx,ny)     !! floodplain width [m]
    ! output
    real    ::  r2fldfrc(nx,ny)     !! area fraction
    real    ::  r2fldare(nx,ny)
    real    ::  rglbfldarea
    ! index 1
    integer, intent(in),    dimension(nx,ny) :: i2nextx  !! point downstream horizontal
    integer, intent(in),    dimension(nx,ny) :: i2nexty  !! point downstream vertical
    integer, intent(in)                      :: nseqriv !! length of 1D sequnece for river
    ! boudary 1
    real,  intent(in),    dimension(nx,ny)    :: r2rivelv !! river bed elevation [m]
    real,  intent(in),    dimension(nx,ny)    :: r2nxtdst !! distance to next grid [m]
    real,  intent(in),    dimension(nx,ny)    :: r2rivare !! river surface area [m2]
    ! local
    integer             ::  ix, iy, iseq, i
    real    ::  rstoall
    real    ::  rstonow
    real    ::  rwthnow
    real    ::  rstopre
    real    ::  rwthpre
    real    ::  rdphpre
    real    ::  rwthinc

    real :: rvel
    real :: r2fldstoold
    
    real,  intent(inout), dimension(nx,ny) :: r2rivinf        !! total inflow to the grid              [m3/s]
    
    ! ================================================
    ! flow calculation
    ! ================================================
    rglbfldarea=0
    do iseq=1, nseqall
       ix=i1seqx(iseq)
       iy=i1seqy(iseq)
                 
       !
       r2fldstoold=r2fldsto(ix,iy)
       rstoall = r2rivsto(ix,iy) + r2fldsto(ix,iy)
       if( rstoall > r2rivstomax(ix,iy) )then
          i=1
          rstopre = r2rivstomax(ix,iy)
          rwthpre = r2rivwth(ix,iy)
          rdphpre = 0.d0
          rwthinc = r2areamat(ix,iy) / r2rivlen(ix,iy) * 0.1d0
          do while( rstoall > r2fldstomax(ix,iy,i) )
             rstopre = r2fldstomax(ix,iy,i)
             rwthpre = rwthpre + rwthinc
             rdphpre = rdphpre + r2fldgrd(ix,iy,i) * rwthinc
             i=i+1
             if( i>10 )then
                rstonow = rstoall - rstopre
                rwthnow = 0.d0
                r2flddph(ix,iy) = rdphpre + rstonow / rwthpre / r2rivlen(ix,iy)
                goto 1000
             endif
          end do
          rstonow =  rstoall - rstopre
          rwthnow = -rwthpre + ( rwthpre**2.d0 + 2.d0 * rstonow/r2rivlen(ix,iy)/r2fldgrd(ix,iy,i) )**0.5d0
          r2flddph(ix,iy) = rdphpre + r2fldgrd(ix,iy,i) * rwthnow
1000      continue
          r2rivsto(ix,iy) = r2rivstomax(ix,iy) + r2rivlen(ix,iy) * r2rivwth(ix,iy) * r2flddph(ix,iy)
          r2rivdph(ix,iy) = r2rivsto(ix,iy) / r2rivlen(ix,iy) / r2rivwth(ix,iy)
          r2fldsto(ix,iy) = rstoall - r2rivsto(ix,iy)
          r2fldsto(ix,iy) = max( r2fldsto(ix,iy), 0.d0 )
          r2fldfrc(ix,iy) = ((-r2rivwth(ix,iy) + rwthpre + rwthnow ) / (rwthinc*10.d0)) + r2rivare(ix,iy) / r2areamat(ix,iy)
          r2fldfrc(ix,iy) = max( r2fldfrc(ix,iy),0.d0)
          r2fldfrc(ix,iy) = min( r2fldfrc(ix,iy),1.d0)
          r2fldare(ix,iy)= r2areamat(ix,iy)*r2fldfrc(ix,iy)
          rglbfldarea = rglbfldarea + r2fldare(ix,iy)
          r2fldwth(ix,iy) = r2fldare(ix,iy)/r2rivlen(ix,iy)
          r2flddph1(ix,iy)= r2fldsto(ix,iy)/r2fldwth(ix,iy)/r2rivlen(ix,iy)
       else
          r2rivsto(ix,iy) = rstoall
          r2rivdph(ix,iy) = max(0.d0,rstoall / r2rivare(ix,iy))
          r2fldsto(ix,iy) = 0.d0
          r2flddph(ix,iy) = 0.d0
          if(rstoall>0d0)then
             r2fldfrc(ix,iy) = r2rivare(ix,iy) / r2areamat(ix,iy)
          else
             r2fldfrc(ix,iy)=0.0
          endif
          r2fldare(ix,iy) = r2fldfrc(ix,iy)*r2areamat(ix,iy)
          r2fldwth(ix,iy) = 0.d0
          r2flddph1(ix,iy)= 0.d0
       endif
    end do

    return
  end subroutine calc_fldstg
  ! ================================================
  ! ================================================  
  subroutine calc_rivout_kine(flag,nx,ny,dt,i1seqx,i1seqy,nseqriv,nseqall,       &
                              i2nextx,i2nexty,r2rivelv,r2elevtn,r2nxtdst,r2rivwth,&
                              r2rivsto,r2rivdph,                                  &
                              r2rivinf,r2rivout,r2rivvel,r2sfcelv,                   &
                              r2rivlen,r2manval,rslpmin)
    ! ================================================
    ! to   Calculate discharge, kinematic wave
    ! by   Augusto Getirana after Dai YAMAZAKI
    ! on   10 Feb 2011
    ! at   LEGOS/CNES, Toulouse, France
    ! Adapted for flow routing implementation in LIS 23 Mar 2012
    ! ================================================
    implicit none
    ! index
    integer             ::  flag             !! flag for rivers (1) of floodplains (2)
    integer             ::  nx               !! number of grids in horizontal
    integer             ::  ny               !! number of grids in vertical
    integer             ::  i1seqx(nx*ny)    !! 1D sequence horizonal
    integer             ::  i1seqy(nx*ny)    !! 1D sequence vertical
    integer             ::  nseqriv          !! length of 1D sequnece for river
    integer             ::  nseqall          !! length of 1D sequnece for river and mouth
    ! boundary
    integer             ::  i2nextx(nx,ny)   !! point downstream horizontal
    integer             ::  i2nexty(nx,ny)   !! point downstream vertical
    integer             ::  dt               !! time step length [sec]
    real    ::  r2rivelv(nx,ny)     !! river bed elevation [m]
    real    ::  r2elevtn(nx,ny)  !! surface elevation [m]
    real    ::  r2nxtdst(nx,ny)  !! distance to next grif [m]
    real    ::  r2rivwth(nx,ny)     !! river width [m]
    ! input
    real    ::  r2rivsto(nx,ny)     !! floodplain storage [m3]
    real    ::  r2rivdph(nx,ny)     !! depth [m]
    ! output
    real    ::  r2rivinf(nx,ny)     !! total inflow to the grid              [m3/s]
    real    ::  r2rivout(nx,ny)     !! outflow from the floodplain reservoir [m3/s]
    real    ::  r2rivvel(nx,ny)     !! velocity [m/s]
    real    ::  r2sfcelv(nx,ny)  !! water surface elevation
    ! boundary 1
    real    ::  r2rivlen(nx,ny)     !! river length [m]
    !! parameter
    real    ::  r2manval(nx,ny)   !! maning coefficient for rivers
    real    ::  rslpmin           !! minimum slope
    ! local
    real    ::  rslope, rarea,rvel, rhydrad
    integer ::  ix,iy,jx,jy,iseq
    ! ==========
    do iseq=1, nseqall
       ix=i1seqx(iseq)
       iy=i1seqy(iseq)
       r2sfcelv(ix,iy) = r2rivelv(ix,iy) + r2rivdph(ix,iy)
       r2rivinf(ix,iy) = 0.d0
    end do
    !
    do iseq=1, nseqriv
       ix=i1seqx(iseq)
       iy=i1seqy(iseq)
       jx=i2nextx(ix,iy)
       jy=i2nexty(ix,iy)   
       rslope = ( r2rivelv(ix,iy)-r2rivelv(jx,jy) ) / r2nxtdst(ix,iy)
       rslope = max( rslope, rslpmin )
       rvel = (1.d0 / r2manval(ix,iy)) * rslope**0.5d0 * r2rivdph(ix,iy)**0.67d0
       
       r2rivvel(ix,iy) = rvel
       rarea  = r2rivwth(ix,iy) * r2rivdph(ix,iy)
       r2rivout(ix,iy) = rarea * rvel
       r2rivout(ix,iy) = min( r2rivout(ix,iy), r2rivsto(ix,iy)/real(dt) )
       r2rivinf(jx,jy) = r2rivinf(jx,jy) + r2rivout(ix,iy)
    end do
    do iseq=nseqriv+1, nseqall
       ix=i1seqx(iseq)
       iy=i1seqy(iseq)
       rslope = rslpmin
       rarea  = r2rivwth(ix,iy) * r2rivdph(ix,iy)
       rvel = (1.d0 / r2manval(ix,iy)) * rslope**0.5d0 * r2rivdph(ix,iy)**0.67d0
       
       r2rivvel(ix,iy) = rvel
       r2rivout(ix,iy) = rarea * rvel
       r2rivout(ix,iy) = min( r2rivout(ix,iy), r2rivsto(ix,iy)/real(dt) )
    end do
    
    return
  end subroutine calc_rivout_kine
  ! ================================================
  ! ================================================   
      subroutine calc_rivout_diff(flag,nx,ny,dt,i1seqx,i1seqy,nseqriv,nseqall,             &
                             i2nextx,i2nexty,r2rivelv,r2elevtn,r2nxtdst,r2rivwth,&
                             r2rivsto,r2rivdph,                                  &
                             r2rivinf,r2rivout,r2rivvel,r2sfcelv,                &
                             r2rivlen,dmanval,dslpmin)
! ================================================
! to   Calculate discharge, diffusive wave
! by   Dai YAMAZAKI
! on   1st dec 2008
! at   IIS,UT
! ================================================
      !use mod_main, only : dmanval,dslpmin
      implicit none
! index
    integer               ::  flag             !! flag for rivers (1) of floodplains (2)
      integer             ::  nx                  !! number of grids in horizontal
      integer             ::  ny                  !! number of grids in vertical
      integer             ::  i1seqx(nx*ny)       !! 1D sequence horizontal
      integer             ::  i1seqy(nx*ny)       !! 1D sequence vertical
      integer             ::  nseqriv             !! length of 1D sequnece for river
      integer             ::  nseqall             !! length of 1D sequnece for river and mouth
! boundary
      integer             ::  i2nextx(nx,ny)      !! point downstream horizontal
      integer             ::  i2nexty(nx,ny)      !! point downstream vertical
      integer             ::  dt                  !! time step length [sec]
      real                ::  r2rivelv(nx,ny)     !! river bed elevation [m]
      real                ::  r2elevtn(nx,ny)     !! surface elevation [m]
      real                ::  r2nxtdst(nx,ny)     !! distance to next grif [m]
      real                ::  r2rivwth(nx,ny)     !! river width [m]
! input
      real                ::  r2rivsto(nx,ny)     !! floodplain storage [m3]
      real                ::  r2rivdph(nx,ny)     !! depth [m]
! output
      real                ::  r2rivinf(nx,ny)     !! total inflow to the grid              [m3/s]
      real                ::  r2rivout(nx,ny)     !! outflow from the floodplain reservoir [m3/s]
      real                ::  r2rivvel(nx,ny)     !! velocity [m/s]
      real                ::  r2sfcelv(nx,ny)     !! water surface elevation

      real                ::  dmanval(nx,ny)      !! maning coefficient for rivers
      real                ::  dslpmin             !! minimum slope
!acvg - 29/12/2010
! boundary 1
      real                ::  r2rivlen(nx,ny) !! river length [m]
! parameter
      real                ::  ddstmth             !! distance at mouth
      parameter              (ddstmth=10000)
! local
      real                ::  r2stoout(nx,ny)
      real                ::  r2rate(nx,ny)
      real                ::  r2velpre(nx,ny)
      real                ::  dslope, darea,dvel
      integer             ::  ix,iy,jx,jy,iseq
! tide
      real                ::  tide
      real                ::  pi
      parameter              (pi=3.14159265358979)
      real                ::  a
      real                ::  b1, b2, b3, b4
      real                ::  c1, c2, c3, c4
      parameter              (a =6.07954307)
      parameter              (b1=3.826602102)
      parameter              (b2=-0.182683388)
      parameter              (b3=1.071233145)
      parameter              (b4=0.185963106)
      parameter              (c1=12.42060122)
      parameter              (c2=23.93446966)
      parameter              (c3=012)
      parameter              (c4=25.81934166)
      real                ::  t1, t2, t3, t4
      data                    t1, t2, t3, t4 / -3.207496885, -7.027942709, -0.855564011, -4.304026087 /
      save                    t1, t2, t3, t4
! ================================================
!! for ocean tide
!!      t1 = t1+dt/60/60
!!      if( t1>c1 ) t1=t1-c1
!!      t2 = t2+dt/60/60
!!      if( t2>c2 ) t2=t2-c2
!!      t3 = t3+dt/60/60
!!      if( t3>c3 ) t3=t3-c3
!!      t4 = t4+dt/60/60
!!      if( t4>c4 ) t4=t4-c4
!!      tide = a + b1*cos(t1/c1*2*pi) + b2*cos(t2/c2*2*pi) + b3*cos(t3/c3*2*pi) + b4*cos(t4/c4*2*pi)
!!      tide = tide + 2 
!!      tide = a
! ==========
      do iseq=1, nseqall
        ix=i1seqx(iseq)
        iy=i1seqy(iseq)
        r2sfcelv(ix,iy) = r2rivelv(ix,iy) + r2rivdph(ix,iy)          !! Diffusive
        r2rivinf(ix,iy) = 0.d0
        r2stoout(ix,iy) = 0.d0
      end do
!
      do iseq=1, nseqriv
        ix=i1seqx(iseq)
        iy=i1seqy(iseq)
        jx=i2nextx(ix,iy)
        jy=i2nexty(ix,iy)
        if( r2sfcelv(ix,iy)-r2sfcelv(jx,jy) >= 0 )then
          dslope = ( r2sfcelv(ix,iy)-r2sfcelv(jx,jy) ) / r2nxtdst(ix,iy)
          dvel =   1.d0 / dmanval(jx,jy) * dslope**0.5d0 * r2rivdph(ix,iy)**0.67d0
          darea  = r2rivwth(ix,iy) * r2rivdph(ix,iy)
        else
          dslope = ( r2sfcelv(jx,jy)-r2sfcelv(ix,iy) ) / r2nxtdst(ix,iy)        
          dvel = - 1.d0 / dmanval(jx,jy) * dslope**0.5d0 * (r2sfcelv(jx,jy)-r2rivelv(ix,iy))**0.67d0
          darea  = r2rivwth(ix,iy) * (r2sfcelv(jx,jy)-r2rivelv(ix,iy))
        endif
        !ag - 05/09/2012
        !r2rivvel(ix,iy) = r2rivvel(ix,iy) + dvel
        r2rivvel(ix,iy) = dvel
        r2rivout(ix,iy) = (r2rivout(ix,iy) + darea * dvel) * 0.5
        if( r2rivout(ix,iy) >= 0 )then
          r2stoout(ix,iy) = r2stoout(ix,iy) + r2rivout(ix,iy)*dt
        else
          r2stoout(jx,jy) = r2stoout(jx,jy) - r2rivout(ix,iy)*dt
        endif
      end do

      do iseq=nseqriv+1, nseqall
        ix=i1seqx(iseq)
        iy=i1seqy(iseq)
        dslope = ( r2sfcelv(ix,iy)-r2elevtn(ix,iy) ) / ddstmth
        dslope = max( dslope, 0.d0 )
        darea  = r2rivwth(ix,iy) * r2rivdph(ix,iy)
        dvel = 1.d0 / dmanval(jx,jy) * dslope**0.5d0 * r2rivdph(ix,iy)**0.67d0
        !ag - 05/09/2012
        !r2rivvel(ix,iy) = r2rivvel(ix,iy) + dvel
        r2rivvel(ix,iy) = dvel
        r2rivout(ix,iy) = (r2rivout(ix,iy) + darea * dvel) * 0.5
        r2stoout(ix,iy) = r2stoout(ix,iy) + r2rivout(ix,iy)*dt
      end do

      do iseq=1, nseqall
        ix=i1seqx(iseq)
        iy=i1seqy(iseq)
        if( r2stoout(ix,iy)> 0.D0 )then
          r2rate(ix,iy) = min( r2rivsto(ix,iy)/r2stoout(ix,iy), 1.d0 )
        else
          r2rate(ix,iy) = 1.d0
        endif
      end do
!
      do iseq=1, nseqriv
        ix=i1seqx(iseq)
        iy=i1seqy(iseq)
        jx=i2nextx(ix,iy)
        jy=i2nexty(ix,iy)
        if( r2rivout(ix,iy) >= 0 )then
          r2rivout(ix,iy) = r2rivout(ix,iy)*r2rate(ix,iy)
        else
          r2rivout(ix,iy) = r2rivout(ix,iy)*r2rate(jx,jy)
        endif
        r2rivinf(jx,jy) = r2rivinf(jx,jy) + r2rivout(ix,iy)
      end do

      do iseq=nseqriv+1, nseqall
        ix=i1seqx(iseq)
        iy=i1seqy(iseq)
        r2rivout(ix,iy) = r2rivout(ix,iy)*r2rate(ix,iy)
      end do
!
      return
      end subroutine calc_rivout_diff
  ! ================================================
  ! ================================================      
  subroutine calc_stonxt(nx,ny,i1seqx,i1seqy,nseqall,           &
                         dt,r2rivinf,r2rivout,r2fldinf,r2fldout,&
                         r2rivsto,r2fldsto)
    ! ================================================
    ! to   Calculate the storage in the next time step in FTCS diff. eq.
    ! by   Augusto GETIRANA after Dai YAMAZAKI
    ! on   10 Feb 2011
    ! at   LEGOS/CNES, Toulouse, France
    ! Adapted for flow routing implementation in LIS 15 Nov 2011
    ! ================================================
    implicit none
    !Dummy vcariables
    integer, intent(in)                      ::  nx           !! number of grids in horizontal
    integer, intent(in)                      ::  ny           !! number of grids in vertical
    integer, intent(in),    dimension(nx*ny) ::  i1seqx       !! 1D sequence horizontal
    integer, intent(in),    dimension(nx*ny) ::  i1seqy       !! 1D sequence vertical
    integer, intent(in)                      ::  nseqall      !! length of 1D sequnece for river and mouth
    integer, intent(in)                      ::  dt           !! time step length [sec]
    real,    intent(in),    dimension(nx,ny) ::  r2rivinf     !! total inflow to river channel            [m3/s]
    real,    intent(inout), dimension(nx,ny) ::  r2rivout     !! outflow from the river reservoir         [m3/s]
    real,    intent(in),    dimension(nx,ny) ::  r2fldinf     !! total inflow to floodplain               [m3/s]
    real,    intent(in),    dimension(nx,ny) ::  r2fldout     !! outflow from the floodplain reservoir    [m3/s]
    real,    intent(inout), dimension(nx,ny) ::  r2rivsto     !! river storage of current time step       [m3]
    real,    intent(inout), dimension(nx,ny) ::  r2fldsto     !! flood plain storage                      [m3]
    
    ! Local variables
    integer             ::  ix, iy, jx, jy, iseq
    ! ================================================
    ! flow calculation
    ! ================================================
    do iseq=1, nseqall
       ix=i1seqx(iseq)
       iy=i1seqy(iseq)
       r2rivsto(ix,iy) = r2rivsto(ix,iy) + r2rivinf(ix,iy)*real(dt) - r2rivout(ix,iy)*real(dt)
       r2rivsto(ix,iy) = max( r2rivsto(ix,iy), 0.d0 )
       r2fldsto(ix,iy) = r2fldsto(ix,iy) + r2fldinf(ix,iy)*real(dt) - r2fldout(ix,iy)*real(dt)
       r2fldsto(ix,iy) = max( r2fldsto(ix,iy), 0.d0 )
       r2rivout(ix,iy) = r2rivout(ix,iy) + r2fldout(ix,iy)
    end do
   
    return
  end subroutine calc_stonxt
  
end module HYMAP_modelMod
