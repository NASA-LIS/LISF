!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module HYMAP3_bifMod
!BOP
!
! !MODULE: HYMAP3_routingMod
!
! !DESCRIPTION:
!  This module contains routines for HyMAP's river bifurcation module
!
! !REVISION HISTORY:
! 1 May 2022: Augusto Getirana, Initial implementation
!
contains

  !=============================================
  subroutine HYMAP3_get_bifurcation_pathways(biffile, nx, ny, sindex, &
       nbif, nbifdelv, bifloc, bifdelv, bifman, bifelv, biflen, bifwth)
    use LIS_logMod
    implicit none
    character(*), intent(in)  :: biffile
    integer,      intent(in)  :: nx,ny
    integer,      intent(in)  :: sindex(nx,ny)
    integer,      intent(in)  :: nbif,nbifdelv
    integer,      intent(out) :: bifloc(nbif,2)
    real,         intent(out) :: bifdelv(nbifdelv),bifman(nbifdelv), &
         bifelv(nbif),biflen(nbif),bifwth(nbif,nbifdelv)
    integer                   :: ii,ix,iy,jx,jy
    logical                   :: file_exists
    integer :: ftn

    inquire(file=trim(biffile),exist=file_exists)
    if(file_exists) then
       ftn = LIS_getNextUnitNumber()
       write(LIS_logunit,*) &
            '[INFO] Get bifurcation topology '//trim(biffile)
       open(ftn,file=trim(biffile), status='old')
       read(ftn,*)ii
       read(ftn,*)bifdelv(:)
       read(ftn,*)bifman(:)
       do ii=1,nbif
          read(ftn,*,end=10)ix,iy,jx,jy,bifelv(ii),biflen(ii),bifwth(ii,:)
          bifloc(ii,1)=sindex(ix,iy)
          bifloc(ii,2)=sindex(jx,jy)
       enddo
       close(ftn)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) &
            '[ERR] Failed to open file '//trim(biffile)// &
            ' in HYMAP3_bifMod'
       call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) '[ERR] Failed in getting bifement rules'
    write(LIS_logunit,*) '[ERR] Check header file '//trim(biffile)
    call LIS_endrun()

  end subroutine HYMAP3_get_bifurcation_pathways
  ! ================================================
  subroutine HYMAP3_calc_bifout_iner(dt, rivelv, rivelv_down, bifelv, &
       biflen, bifwth, bifsto, bifsto_down, rivdph, rivdph_down, manval, &
       grv, bifout, bifout_pre, rivdph_pre, rivdph_pre_down)
    use LIS_logMod
    ! ================================================
    ! Calculate discharge, local inertia
    ! Augusto Getirana
    ! 16 Apr 2014
    ! Adapted for flow routing implementation in LIS on 5 Aug 2015
    ! Adapted for single grid cell computation in 19 Feb 2016
    ! Adapted for bifurcation simulation in 8 Apr 2022
    ! ================================================
    implicit none
    real,    intent(in)    :: dt              !time step length [sec]
    real,    intent(in)    :: rivelv          !river/floodplain bed elevation [m]
    real,    intent(in)    :: rivelv_down     !downstream river/floodplain bed elevation [m]
    real,    intent(in)    :: biflen          !distance to next grid [m]
    real,    intent(in)    :: bifwth          !river/floodplain width [m]
    real,    intent(in)    :: bifsto          !river/floodplain storage [m3]
    real,    intent(in)    :: bifsto_down     !downstream river/floodplain storage [m3]
    real,    intent(in)    :: rivdph          !river/floodplain depth [m]
    real,    intent(in)    :: rivdph_down     !dowstream river/floodplain depth [m]
    real,    intent(in)    :: manval          !manning coefficient for river/floodplain
    real,    intent(in)    :: grv             !gravity [m/s2]
    real,    intent(out)   :: bifout          !river/floodplain outflow  [m3/s]
    real,    intent(inout) :: bifout_pre      !previous outflow [m3/s]
    real,    intent(in)    :: rivdph_pre      !previous river depth [m]
    real,    intent(in)    :: rivdph_pre_down !previous downstream river depth [m]
    real,    intent(in)    ::  bifelv          !bifurcation bed elevation [m]


    real                   :: sfcelv          !upstream water surface elevation
    real                   :: sfcelv_down     !downstream water surface elevation
    real                   :: sfcelv_pre      !previous step water surface elevation (t-1) [m]
    real                   :: sfcelv_pre_down !previous step downstream water surface elevation (t-1) [m]
    real                   :: dslope, darea
    real                   :: dflw,dout_pre,dflw_pre,dflw_imp

    ! ================================================

    if(bifsto==0)then
       bifout=0.
       goto 1
    endif

    sfcelv=rivelv+rivdph
    sfcelv_down=rivelv_down+rivdph_down
    sfcelv_pre=rivelv+rivdph_pre
    sfcelv_pre_down=rivelv_down+rivdph_pre_down

    dslope=(sfcelv-sfcelv_down)/biflen
    dflw=max(sfcelv,sfcelv_down)-bifelv

    !this imposes unidirectional flow; no negative flow
    if(dflw<=0.or.dslope<=0)then
       bifout=0.
       goto 1
    else!(dslope>0)then
       dflw=min(dflw,rivdph)
    endif

    darea=bifwth*dflw
    dflw_pre=max(0.,max(sfcelv_pre,sfcelv_pre_down)-rivelv)
    dflw_imp=sqrt(dflw*dflw_pre)
    if(dflw_imp==0.)dflw_imp=dflw

    if(dflw_imp>0..and.darea>0.)then
       dout_pre=bifout_pre/bifwth
       bifout=bifwth*(dout_pre+grv*dt*dflw_imp*dslope)/&
            (1.+grv*dt*manval**2.*abs(dout_pre)*dflw_imp**(-7./3))
       if(bifout>=0.)then
          bifout=min(bifout,bifsto/dt)
       else
          bifout=-min(abs(bifout),bifsto_down/dt)
       endif
    else
       bifout=0.
    endif
1   continue
    bifout_pre=bifout
  end subroutine HYMAP3_calc_bifout_iner
  !=============================================
end module HYMAP3_bifMod
