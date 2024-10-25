!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: clsmf25_turb
!
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
! 23 Nov 2012: David Mocko, Added Catchment Fortuna-2.5
!
! !INTERFACE:
subroutine clsmf25_turb(N_cat, met_force, cat_progn, &
     cat_diagn, zol, dze, lai,  &
     RA1, ETURB1, DEDQA1, DEDTC1, HSTURB1, DHSDQA1, DHSDTC1,                &
     RA2, ETURB2, DEDQA2, DEDTC2, HSTURB2, DHSDQA2, DHSDTC2,                &
     RA4, ETURB4, DEDQA4, DEDTC4, HSTURB4, DHSDQA4, DHSDTC4,                &
     RAS, ETURBS, DEDQAS, DEDTCS, HSTURBS, DHSDQAS, DHSDTCS,                & 
     sfc_turb_scheme )
    

    ! qliu+reichle: 14 Aug 2008 - land interface to subroutine louissurface() 
    ! qliu+reichle: 30 Oct 2008 - force derivatives to be non-negative in case
    !                             extra derivatives are included
    ! reichle     : 29 Nov 2010 - added Helfand Monin-Obukhov scheme
    ! reichle     :  2 Feb 2011 - added T2m, Q2m diagnostics
!EOP   
    use LIS_coreMod
    use LIS_logMod
    use clsmf25_surfacelayer
    use clsmf25_types
    use clsmf25_drv_types
    use clsmf25_esat_qsat
    use clsmf25_MAPL_constants, ONLY: &
       GRAV => MAPL_GRAV, CPAIR => MAPL_CP, VIREPS => MAPL_VIREPS, &
       RGAS => MAPL_RGAS


    implicit none
    
    integer,  intent(in) :: N_cat
    integer,  intent(in) :: sfc_turb_scheme
    
    type(met_force_type), dimension(N_cat), intent(in) :: met_force
    type(cat_progn_type), dimension(N_cat), intent(in) :: cat_progn
    type(cat_diagn_type), dimension(N_cat), intent(in) :: cat_diagn
    
    real, dimension(N_cat), intent(in) :: zol, dze, lai
    
    real, dimension(N_cat), intent(out) :: RA1, RA2, RA4, RAS
    real, dimension(N_cat), intent(out) :: ETURB1, DEDQA1, DEDTC1
    real, dimension(N_cat), intent(out) :: ETURB2, DEDQA2, DEDTC2
    real, dimension(N_cat), intent(out) :: ETURB4, DEDQA4, DEDTC4
    real, dimension(N_cat), intent(out) :: ETURBS, DEDQAS, DEDTCS
    real, dimension(N_cat), intent(out) :: HSTURB1,DHSDQA1,DHSDTC1
    real, dimension(N_cat), intent(out) :: HSTURB2,DHSDQA2,DHSDTC2
    real, dimension(N_cat), intent(out) :: HSTURB4,DHSDQA4,DHSDTC4
    real, dimension(N_cat), intent(out) :: HSTURBS,DHSDQAS,DHSDTCS
    
    
    ! local variables

    integer :: k, n
    
    real    :: rhoair

    real, dimension(N_cat,4) :: tc, qa, z0, ch, cq
    
    real, dimension(      4) :: dt, dq, dcqdq, dcqdt, dchdt, dchdq
    
    ! Louis

    real,    dimension(:),   allocatable :: cn, ri, zt, zq, uuu, ucn, re
    real,    dimension(:,:), allocatable :: cm, ww, dch, dcq 
    
    ! Helfand Monin-Obukhov

    integer, dimension(:),   allocatable :: tmpintvec
    
    real,    dimension(:),   allocatable :: tmprealvec, psmb, psl
    real,    dimension(:),   allocatable :: RHOH,VKM,USTAR,XX,YY,CU,CT,RIB,ZETA,WS

    real,    dimension(:,:), allocatable :: t2m,q2m
    real,    dimension(:),   allocatable :: u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m
    
    ! definition of land in sfclayer module  
  
    integer, parameter :: SFCTYPE     = 3     

    ! parameters for the Helfand Monin-Obuhkhov surface layer routine
  
    integer, parameter :: HELF_NITER  = 6  ! number of internal iterations
    integer, parameter :: HELF_Z0     = 1  ! CHOOSEZ0 in GEOS_CatchGridComp.F90 (ocean only)
  
    ! additional option for *Louis* surface layer turbulence scheme
    !
    ! turn on/off extra derivatives in turbulence calcuations
    ! (derivatives of drag coeffs w.r.t. T and Q and cross-derivatives
    !  of fluxes, eg. d_SensibleHeat/d_Q)

    logical, parameter :: incl_Louis_extra_derivs = .true.

    ! turn on/off simplistic estimate of ww in turb calculations

    logical, parameter :: incl_Louis_ww           = .false.

    real,    parameter :: HPBL        = 1000.                 

    ! prepare inputs that are common to all sfc_turb_scheme choices
    
    do k=1,N_cat
       
       tc(k,1) = cat_progn(k)%tc1
       tc(k,2) = cat_progn(k)%tc2
       tc(k,3) = cat_progn(k)%tc4
       tc(k,4) = cat_diagn(k)%tpsn(1)

       qa(k,1) = cat_progn(k)%qa1
       qa(k,2) = cat_progn(k)%qa2
       qa(k,3) = cat_progn(k)%qa4

       if(cat_diagn(k)%tpsn(1).lt.200.or.cat_diagn(k)%tpsn(1).gt.500) then 
          write(LIS_logunit,*) 'tpsn ',k,cat_diagn(k)%tpsn(1)
       endif
       qa(k,4) = qsat(cat_diagn(k)%tpsn(1),met_force(k)%Psurf/100.)
       
       z0(k,1:4) = zol(k)
       
    end do

    ! allocate and initialize variables as needed
    
    select case (sfc_turb_scheme)
       
    case (0)     ! Louis
       
       allocate(cn( N_cat))
       allocate(ri( N_cat))
       allocate(zt( N_cat))
       allocate(zq( N_cat))
       allocate(uuu(N_cat))
       allocate(ucn(N_cat))
       allocate(re( N_cat))
       
       allocate(cm( N_cat,4))
       allocate(ww( N_cat,4))
       allocate(dch(N_cat,4))
       allocate(dcq(N_cat,4))

       ! note that cm and z0 are intent(inout) in louissurface()

       cm = 0.
       ww = 0.
       
    case (1)     ! Helfand Monin-Obukhov
       
       allocate(tmpintvec( N_cat))       

       allocate(tmprealvec(N_cat))
       allocate(psmb(      N_cat)) 
       allocate(psl(       N_cat)) 
       allocate(RHOH(      N_cat))
       allocate(VKM(       N_cat))
       allocate(USTAR(     N_cat))
       allocate(XX(        N_cat))
       allocate(YY(        N_cat))
       allocate(CU(        N_cat))
       allocate(CT(        N_cat))
       allocate(RIB(       N_cat))
       allocate(ZETA(      N_cat))
       allocate(WS(        N_cat))
                           
       allocate(t2m(       N_cat,4))
       allocate(q2m(       N_cat,4))

       allocate(u2m(       N_cat))
       allocate(v2m(       N_cat))
       allocate(t10m(      N_cat))
       allocate(q10m(      N_cat))
       allocate(u10m(      N_cat))
       allocate(v10m(      N_cat))
       allocate(u50m(      N_cat))
       allocate(v50m(      N_cat))
              
    case default
       
       write(LIS_logunit,*) '[ERR] clsmf25_turb(): '
       write(LIS_logunit,*) 'Unknown selection for sfc_turb_scheme.'
       call LIS_endrun()
       
    end select
    
    
    ! ----------------------------------------------------------------------------
    !
    ! call surface layer subroutine
    
    do n=1,4
       
       select case (sfc_turb_scheme)
          
          ! ----------------------------------------------------------------------
          
       case (0)     ! Louis
          
          ! note that cm and z0 are intent(inout) in louissurface()
          
          call louissurface(SFCTYPE, n, met_force%Wind, ww, met_force%Psurf,    &
               met_force%Tair, tc, met_force%Qair, qa, met_force%Rainf_C,       &
               lai, z0, dze,                                                    &
               cm, cn, ri, zt, zq, ch, cq, uuu, ucn, re, dch, dcq )
          
          if (incl_Louis_ww) then
             
             ! simplistic estimate of ww for off-line integrations
             ! (for a better job and better consistency with coupled integrations,
             !  ww would probably have to be included in cat_progn)
             ! reichle+qliu,  9 Oct 2008
             
             ww(:,n) = ch(:,n)*(tc(:,n)-met_force(:)%Tair-                      & 
                  (GRAV/CPAIR)*dze(:))/met_force(:)%Tair + VIREPS*cq(:,n)*      &
                  (qa(:,n)-met_force(:)%Qair)
             ww(:,n) = max(ww(:,n),0.0)
             ww(:,n) = (HPBL*GRAV*ww(:,n))**(2./3.)  
             
             call louissurface(SFCTYPE, n, met_force%Wind, ww, met_force%Psurf, &
                  met_force%Tair, tc, met_force%Qair, qa, met_force%Rainf_C,    &
                  lai, z0, dze,                                                 &
                  cm, cn, ri, zt, zq, ch, cq, uuu, ucn, re, dch, dcq )
             
          end if

          ! ----------------------------------------------------------------------

       case (1)     ! Helfand Monin-Obukhov
          
          psmb = met_force%Psurf * 0.01            ! convert to millibar
          
          ! Approximate pressure at top of surface layer: 
          ! hydrostatic, eqn of state using avg temp and press
          
          tmprealvec = (dze*GRAV) / (RGAS*(met_force%Tair + tc(:,n)))
          
          psl = psmb * (1. - tmprealvec) / (1. + tmprealvec)
          
          ! surface type
          
          tmpintvec(1:N_cat) = SFCTYPE
          
          ! for wind use u-component=wind-speed, v-component=0.
          
          tmprealvec = 0. 

          call helfsurface( met_force%Wind, tmprealvec, met_force%Tair, tc(:,n), &
               met_force%Qair, qa(:,n), psl, psmb, z0(:,n), lai,                 &
               tmpintvec, dze, HELF_NITER, N_cat,                                &
               RHOH,ch(:,n),VKM,USTAR,XX,YY,CU,CT,RIB,ZETA,WS,                   &
               t2m(:,n),q2m(:,n),u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m,HELF_Z0)
          
          cq(:,n)  = ch(:,n) 

          ! ----------------------------------------------------------------------

       end select
       
    end do
    
    ! -------------------------------------------------------------------------------
    !
    ! finalize output arguments
    
    do k=1,N_cat
       
       ! rhoair from GEOS_CatchGridComp.F90
       
       rhoair     = RGAS*met_force(k)%Tair*(1.+VIREPS*met_force(k)%Qair)
       rhoair     = met_force(k)%Psurf/rhoair
       
       RA1(k)     = rhoair/ch(k,1)
       RA2(k)     = rhoair/ch(k,2)
       RA4(k)     = rhoair/ch(k,3)
       RAS(k)     = rhoair/ch(k,4)
       
       dq         = qa(k,:) - met_force(k)%Qair
       
       ETURB1(k)  = cq(k,1)*dq(1)
       ETURB2(k)  = cq(k,2)*dq(2)
       ETURB4(k)  = cq(k,3)*dq(3)
       ETURBS(k)  = cq(k,4)*dq(4)
       
       dt         = tc(k,:) - met_force(k)%Tair
       
       HSTURB1(k) = CPAIR*ch(k,1)*dt(1)
       HSTURB2(k) = CPAIR*ch(k,2)*dt(2)
       HSTURB4(k) = CPAIR*ch(k,3)*dt(3)
       HSTURBS(k) = CPAIR*ch(k,4)*dt(4)
       
       ! derivatives
       
       DEDQA1(k)  = cq(k,1)
       DEDQA2(k)  = cq(k,2)
       DEDQA4(k)  = cq(k,3)
       DEDQAS(k)  = cq(k,4)
       
       DHSDQA1(k) = 0.
       DHSDQA2(k) = 0.
       DHSDQA4(k) = 0.
       DHSDQAS(k) = 0.
       
       DHSDTC1(k) = CPAIR*ch(k,1)
       DHSDTC2(k) = CPAIR*ch(k,2)
       DHSDTC4(k) = CPAIR*ch(k,3)
       DHSDTCS(k) = CPAIR*ch(k,4)
       
       DEDTC1(k)  = 0.
       DEDTC2(k)  = 0.
       DEDTC4(k)  = 0.
       DEDTCS(k)  = 0.
       

       if ((sfc_turb_scheme==0) .and. incl_Louis_extra_derivs) then

          ! overwrite derivatives with cross-terms
          
          ! include derivatives of drag coefficients and cross-derivatives
          ! (evap w.r.t. temperature, sensible w.r.t. humidity)
          
          ! "dcq" from louissurface() is d(cq)/d(TVA) where TVA = T*(1+VIREPS*Q)
          !
          ! moreover, cq=cq(Ri) where Ri is proportional to deltaTVA=TVA-TVS
          ! (where TVS is the virtual surface temperature) -> this produces a 
          ! minus sign
          !
          ! taken together, we need d(cq)/dq = d(cq)/dTV * dTV/dQ = -dcq * VIREPS*T
          
          dcqdq      = -dcq(k,:)*VIREPS*tc(k,:)
          
          DEDQA1(k)  = max( 0., cq(k,1) + dcqdq(1)*dq(1) )
          DEDQA2(k)  = max( 0., cq(k,2) + dcqdq(2)*dq(2) )
          DEDQA4(k)  = max( 0., cq(k,3) + dcqdq(3)*dq(3) )
          DEDQAS(k)  = max( 0., cq(k,4) + dcqdq(4)*dq(4) )

          ! similarly, d(cq)/dt = d(cq)/dTV * dTV/dT = -dcq * (1+VIREPS*Q)
          
          dcqdt      = -dcq(k,:)*(1.+VIREPS*qa(k,:))
          
          DEDTC1(k)  = max( 0., dcqdt(1)*dq(1) )
          DEDTC2(k)  = max( 0., dcqdt(2)*dq(2) ) 
          DEDTC4(k)  = max( 0., dcqdt(3)*dq(3) )
          DEDTCS(k)  = max( 0., dcqdt(4)*dq(4) )
          
          ! "dch" from louissurface() is d(ch)/d(TVA) where TVA = T*(1+VIREPS*Q)
          !
          ! moreover, ch=ch(Ri) where Ri is proportional to deltaTVA=TVA-TVS
          ! (where TVS is the virtual surface temperature) -> this produces a 
          ! minus sign
          !
          ! taken together, we need d(ch)/dT = d(ch)/dTV * dTV/dT = -dch * (1+VIREPS*Q)
          
          dchdt      = -dch(k,:)*(1.+VIREPS*qa(k,:))
          
          DHSDTC1(k) = max( 0., CPAIR*(ch(k,1) + dchdt(1)*dt(1)) )
          DHSDTC2(k) = max( 0., CPAIR*(ch(k,2) + dchdt(2)*dt(2)) )
          DHSDTC4(k) = max( 0., CPAIR*(ch(k,3) + dchdt(3)*dt(3)) )
          DHSDTCS(k) = max( 0., CPAIR*(ch(k,4) + dchdt(4)*dt(4)) )

          ! similarly, d(ch)/dQ = d(ch)/dTV * dTV/dQ = -dch * VIREPS*T
          
          dchdq      = -dch(k,:)*VIREPS*tc(k,:)

          DHSDQA1(k) = max( 0., CPAIR*dchdq(1)*dt(1) )  
          DHSDQA2(k) = max( 0., CPAIR*dchdq(2)*dt(2) )
          DHSDQA4(k) = max( 0., CPAIR*dchdq(3)*dt(3) )
          DHSDQAS(k) = max( 0., CPAIR*dchdq(4)*dt(4) )
          
       end if
       
    end do
    
    ! -----------------------------------------------------------------
    !
    ! deallocate 
    
    select case (sfc_turb_scheme)
       
    case (0)     ! Louis
       
       deallocate(cn)
       deallocate(ri)
       deallocate(zt)
       deallocate(zq)
       deallocate(uuu)
       deallocate(ucn)
       deallocate(re)
       
       deallocate(cm)
       deallocate(ww)
       deallocate(dch)
       deallocate(dcq)
       
    case (1)     ! Helfand Monin-Obukhov

       deallocate(tmpintvec)
              
       deallocate(tmprealvec)
       deallocate(psmb)   
       deallocate(psl)    
       deallocate(RHOH)  
       deallocate(VKM)   
       deallocate(USTAR) 
       deallocate(XX)    
       deallocate(YY)    
       deallocate(CU)    
       deallocate(CT)    
       deallocate(RIB)   
       deallocate(ZETA)  
       deallocate(WS)    

       deallocate(t2m)   
       deallocate(q2m)   

       deallocate(u2m)   
       deallocate(v2m)   
       deallocate(t10m)  
       deallocate(q10m)  
       deallocate(u10m)  
       deallocate(v10m)  
       deallocate(u50m)  
       deallocate(v50m)  
    
    end select
    
  return
  end subroutine clsmf25_turb
