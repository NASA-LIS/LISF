!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module HYMAP3_resopMod_Yassin
  !module defining Yassin's reservoir operation variables
  !Augusto Getirana
  !NASA Goddard Space Flight Center
  !4 Apr 2025
  use HYMAP3_modelMod

contains
  !=============================================
  !=============================================  
  subroutine HYMAP3_resop_main_glb_Yassin(mis,nz,time,dt,inflow,&
               rivsto,fldsto,runoff,sfcelv1,&
               maxsto,inidis,inisto,dwndis,deadis,&
               minsto_mo,nupsto_mo,uppsto_mo,&
               mindis_mo,nupdis_mo,uppdis_mo,&
               reg1,reg2,reg3,&
               outflow,inflow_down)

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
    real,    intent(in)    :: sfcelv1              !river bed elevation [m]
    real*8,  intent(in)    :: maxsto              !Max reservoir capacity (m3)
    real*8,  intent(in)    :: inidis              !Initial discharge (m3 s-1)
    real*8,  intent(in)    :: inisto              !Initial storage (m3)
    real*8,  intent(in)    :: dwndis              !Downstream channel capacity
    real*8,  intent(in)    :: deadis              !Dead storage fraction of maximum storage ( 0.1 means 10% of max storage)
    real*8,  intent(in)    :: minsto_mo(12)       !Monthly min storage 12 values one per month (m3)
    real*8,  intent(in)    :: nupsto_mo(12)       !Monthly normal upper storage (m3)
    real*8,  intent(in)    :: uppsto_mo(12)       !Monthly upper storage (m3)
    real*8,  intent(in)    :: mindis_mo(12)       !Monthly min release (m3 s-1)
    real*8,  intent(in)    :: nupdis_mo(12)       !Monthly normal upper release (m3 s-1)
    real*8,  intent(in)    :: uppdis_mo(12)       !Monthly upper release (m3 s-1)
    real*8,  intent(in)    :: reg1                !Polynomial equation coef. a (y=ax^2 + bx +c)
    real*8,  intent(in)    :: reg2                !Polynomial equation coef. b (y=ax^2 + bx +c)
    real*8,  intent(in)    :: reg3                !Polynomial equation coef. c (y=ax^2 + bx +c)
    
    real,    intent(inout) :: outflow              !grid cell outflow (river+floodplain) [m3/s]
    real,    intent(inout) :: inflow_down              !downstream grid cell inflow [m3/s]

    real*8                 :: volsim
    real                   :: outsim

    real*8                 :: volres
    real*8                 :: Rx
    real*8                 :: sfcelv

    integer                :: imon

sfcelv=dble(sfcelv1*1d4)/1d4
Rx=0.

#if 0
print*,'inflow ',inflow
print*,'rivsto ',rivsto
print*,'fldsto ',fldsto
print*,'runoff ',runoff
print*,'sfcelv ',sfcelv
print*,'maxsto ',maxsto
print*,'inidis ',inidis
print*,'inisto ',inisto
print*,'dwndis ',dwndis
print*,'deadis ',deadis
print*,'minsto_mo ',minsto_mo
print*,'nupsto_mo ',nupsto_mo
print*,'uppsto_mo ',uppsto_mo
print*,'mindis_mo ',mindis_mo
print*,'nupdis_mo ',nupdis_mo
print*,'uppdis_mo ',uppdis_mo
print*,'reg1 ',reg1
print*,'reg2 ',reg2
print*,'reg3 ',reg3

stop
#endif

    imon=mod(int(time),1000)/100
    volsim=max(0.,rivsto+fldsto+(inflow+runoff)*dble(dt))
    volres = max(0.,(reg1 * sfcelv**2 + reg2 * sfcelv + reg3)) !/ maxsto
#if 0
print*,reg1,reg2,reg3,sfcelv
print*,(dble(reg1) * dble(sfcelv)**2 + dble(reg2) * dble(sfcelv) + dble(reg3))
print*,reg1 * sfcelv**2,reg2 * sfcelv,reg3
print*,sfcelv**2
#endif
!print*,volres
    call SIMRES15(dble(inflow),volres,deadis,maxsto,dwndis,(minsto_mo(imon)/maxsto),&
                 (nupsto_mo(imon)/maxsto),(uppsto_mo(imon)/maxsto),mindis_mo(imon),&
                 nupdis_mo(imon),uppdis_mo(imon),Rx,dble(dt),outsim)

    if(outsim*dt>volsim)outsim=volsim/dt
write(LIS_logunit,*) int(time),volsim,volres,sfcelv,outsim
!stop
    !update HyMAP variables
    inflow_down=inflow_down-outflow+outsim
    outflow=outsim

  end subroutine HYMAP3_resop_main_glb_Yassin
  !=============================================
  subroutine SIMRES15(Inflow, stoflo, Ld, smax, qds, Lc, Ln, Lf, qmin, qnorm, qnd, Rx, dt, flowSIM)
    implicit none   
    !> Input variables.
    real*8, intent(in) :: stoflo
    real*8, intent(in) :: Ld                  ! = resrv%deadst
    real*8, intent(in) :: Inflow              ! = resrv%flowINF(1:2)          ! observed u/s inflow
    real*8, intent(in) :: smax                ! = resrv%SMAX                  ! maximum reservoir storage (m3)
    real*8, intent(in) :: qds                 ! = resrv%qmaxmax               ! downstream channel capacity (m3/s)
    real*8, intent(in) :: Lc                  ! = resrv%dsto(mId)/resrv%SMAX  ! monthly lower storage threshold
    real*8, intent(in) :: Ln                  ! = resrv%nsto(mId)/resrv%SMAX  ! monthly upper storage threshold
    real*8, intent(in) :: Lf                  ! = resrv%ndsto(mId)/resrv%SMAX ! inflow starts affecting release
    real*8, intent(in) :: qmin                ! = resrv%Qmin(mId)
    real*8, intent(in) :: qnorm               ! = resrv%Qnor(mId)
    real*8, intent(in) :: qnd                 ! = resrv%Qnd(mId)
    real*8, intent(in) :: Rx                  ! = 2*resrv%RXN*rnd + resrv%RXN ! Rx not inclued in the Matlab code, but added based on MESH code
    real*8, intent(in) :: dt                  ! = dt                          ! time-step length/duration

    !> Local variables.
    real*8  Fu,rnd
    integer it

    !> Output variables.
    real, intent(out) :: flowSIM  ! = resrv%flowSIM(1:2)
    !real, intent(out) :: stoSIM   ! = resrv%stoSIM(1:2)
    real :: stoSIM

    real,    parameter :: icor  = 1.
    integer, parameter :: cse   = 1
    integer, parameter :: niter = 1

        call random_number(rnd)

        !> Determine release and storage at current time-step.
!        stoflo = stoSIM(t - 1) + (dt/2.0)*(flowI(t)*icor + flowI(t - 1)*icor - flowSIM(t - 1))
        ! print*,'stoflo is ',stoflo
!        Inflow = flowI(t)
         Fu = stoflo/smax
        ! print*,'Fu is ',Fu
!print*,Fu,stoflo,smax,Rx
!        flowSIM(t) = S2Q1JAN(Fu, Ld, Lc, Ln, Lf, qmin, qnorm, qnd, dt, smax, Inflow, qds, cse, Rx)
        flowSIM = S2Q1JAN(Fu, Ld, Lc, Ln, Lf, qmin, qnorm, qnd, dt, smax, Inflow, qds, cse, Rx)
!print*,flowSIM
        do it = 1, niter
            stoSIM  = stoflo - (dt/2.0)*(flowSIM)
            Fu = stoSIM/smax
            flowSIM = S2Q1JAN(Fu, Ld, Lc, Ln, Lf, qmin, qnorm, qnd, dt, smax, Inflow, qds, cse, Rx)
!print*,flowSIM
        end do
        if (stoSIM > smax) then
            flowSIM = flowSIM + ((stoSIM - smax)/dt)
            stoSIM = smax
        end if

        ! print*,'flowSIM(t-1) is ',flowSIM(t-1)
        ! print*,'flowSIM(t) is ',flowSIM(t)

    end subroutine SIMRES15
  !=============================================
  real function S2Q1JAN(Fu, Ld, Lc, Ln, Lf, qmin, qnorm, qnd, dt, smax, Inflow, qds, cse, Rx) result(qout)
    implicit none
    !> Input variables.
    real*8, intent(in) :: Fu, Ld, Lc, Ln, Lf, qmin, qnorm, qnd, dt, smax, Inflow, qds, Rx
    integer, intent(in) :: cse

        ! print*,'Fu is ',Fu
        ! print*,'Ld is ',Ld
        ! print*,'Lc is ',Lc
        ! print*,'Ln is ',Ln
        ! print*,'Lf is ',Lf
        ! print*,'qmin is ',qmin
        ! print*,'qnorm is ',qnorm
        ! print*,'qnd is ',qnd
        ! print*,'dt is ',dt
        ! print*,'smax is ',smax
        ! print*,'Inflow is ',Inflow
        ! print*,'qds is ',qds
        ! print*,'Rx is ',Rx

        !> Initialize the output variables.
    qout = 0.0
!print*,Fu,Ld,Lc,Ln,Lf
    !> Determine release at current time-step.
    if (Fu <= Ld) then
      qout = 0.0
    else if (Fu > Ld .and. Fu <= Lc) then
      qout = min(qmin,((Fu - Ld)*smax/(dt))) + Rx
    else if (Fu > Lc .and. Fu <= Ln) then
      qout = qmin + (qnorm - qmin)*((Fu - Lc)/(Ln - Lc)) + Rx
    else if (Fu > Ln .and. Fu <= Lf) then
      if (cse == 1) then
        qout = qnorm + ((Fu - Ln)/(Lf - Ln))*(qnd - qnorm) + Rx
      else
        qout = qnorm + ((Fu - Ln)/(Lf - Ln))*max((Inflow - qnorm), (qnd - qnorm)) + Rx
      end if
    else
      qout = min(max(((Fu - Lf)*smax/dt), qnd) + Rx, qds)
    end if

    !> Check for negative value.
    qout = max(qout, 0.0)
    ! print*,'qout is ',qout

  end function
  !=============================================

end module HYMAP3_resopMod_Yassin
