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
!
! !ROUTINE: hyssib_main
! \label{hyssib_main}
!
! !REVISION HISTORY:
!  21 Apr 2004: David Mocko, Initial Code
!  30 Apr 2004: David Mocko, Revisions to fix MPI version
!  25 Nov 2008: Chuck Alonge, Updates for LIS 5.0
!  27 Oct 2010: David Mocko, changes for HY-SSiB in LIS6.1
!
! !INTERFACE:
subroutine hyssib_main(n)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_surface,LIS_domain
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_albedoMod, only : LIS_alb
  use LIS_logMod,    only : LIS_logunit, LIS_endrun
  use LIS_histDataMod
  use hyssib_lsmMod  ! HY-SSiB tile variables
  use hyssibalb_module

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
! !DESCRIPTION:
! This code is the offline HY-SSiB Model (Hydrology with Simple SiB)
! Authors: Y.C. Sud and David M. Mocko, NASA Goddard Space Flight Center
! Code Maintenance by: David Mocko <mocko@climate.gsfc.nasa.gov>
!
! Referred Papers:
!
! SMP1:  Sud, Y.C., and D.M. Mocko, 1999: New snow-physics to complement
!        SSiB.  Part I: Design and evaluation with ISLSCP Initiative I
!        datasets.  J. Meteor. Soc. Japan, Vol. 77, No. 1B, 335-348.
!
! SMP2:  Mocko, D.M., Walker, G.K., and Y.C. Sud, 1999: New snow-physics
!        to complement SSiB.  Part II: Effects on soil moisture
!        initialization and simulated surface fluxes, precipitation and
!        hydrology of GEOS II GCM.
!        J. Meteor. Soc. Japan, Vol. 77, No. 1B, 349-366.
!
! SMP3:  Mocko, D.M., and Y.C. Sud, 2001: Refinements to SSiB with an
!        emphasis on snow-physics: Evaluation and validation using GSWP
!        and Valdai data.  Earth Interactions, Vol. 5, (5-001), 31 pp.
!
! Other papers in comments:
!
! Sellers et al (two sets of authors), 1996: J. Climate, Vol. 6, No. 4,
!   A revised land surface parameterization (SiB2) for atmospheric GCMs
!     Part I  - Model formulation.  p. 676-705.
!     Part II - The generation of global fields of terrestrial
!               biophysical parameters from satellite data.  p. 706-737.
!
! GEWEX/GSWP, 1995: International GEWEX Project Office.
!                   Version 1.0, 1 Dec 1995.  46 pp.
!
! Douville et al, 1995: A new snow parameterization for the Meteo-France
!                       climate model.  Parts I and II.
!                       Clim. Dyn., Vol. 12, 21-52.
! Mahrt and Sun, 1995: The subgrid velocity scale in the bulk
!                      aerodynamic relationship for spatially
!                      averaged scalar fluxes.
!                      Mon. Wea. Rev., Vol. 123, No. 10, 3032-3041.
!
! Xue et al, 1991: A simplified biosphere model for global climate
!                  studies.  J. Climate, Vol. 4, No. 3, 345-364.
!
! Sellers et al, 1986: A simple biosphere model (SiB) for use within
!                      general circulation models.
!                      J. Atmos Sci., Vol. 43, No. 6, 505-531.
!
! Milly and Eagleson, 1982: Parameterization of moisture and heat fluxes
!                           across the land surface for use in
!                           atmospheric general circulation models.
!                           Dept. of Engineering, Massachusetts Inst.
!                           of Technology, Rep. 279, 159 pp.
!
! Deardorff, 1972: Parameterization of planetary boundary-layer for use
!                  in general circulation models.
!                  Mon. Wea. Rev., Vol. 100, No. 2, 93-106.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  integer, parameter :: istrip=1  ! Number of tiles in 1-D version

  integer            :: index1, i
  integer            :: flag
  ! Declare constants
  ! -----------------
  real grav, pie, cpair, gasr, hlat, snomel
  real clai, cw, cs, dkfac, epsfac, stefan
  real thres, sskin, tf, tf_snow, tf_drain, delsig

  ! Declare description of point and time
  ! -------------------------------------
  real lonco(istrip), latco(istrip)
  integer ityp(istrip)
  integer nymd, nhms, month

  ! Declare vegetation parameters
  ! -----------------------------
  ! Monthly-varying:
  real z0(istrip), z1(istrip), z2(istrip), dd(istrip)
  real vcover(istrip,2), zlt(istrip,2), green(istrip)
  real rbc(istrip), rdc(istrip)

  ! Non-varying:
  real chil(istrip), topt(istrip), tll(istrip), tu(istrip)
  real defac(istrip), ph1(istrip), ph2(istrip), rootd(istrip,2)
  real rstpar(istrip,3)

  ! Declare soil parameters
  ! -----------------------
  real phsat(istrip), poros(istrip), bee(istrip), satco(istrip)
  real slope(istrip), zdepth(istrip,3), wiltp(istrip)
  real psilow(istrip), topostd(istrip), tdd(istrip)

  ! Declare radiation parameters
  ! ----------------------------
  !Removed - now included in hyssibalb_module to reduce model i/o
  !real*4 cedfu(13,12,9), cedir(13,12,9,3), cedfu1(2,13,12,6,3)
  !real*4 cedir1(2,13,12,6,3,3), cedfu2(2,13,12,6,3), xmiu(12,3)
  !real*4 cedir2(2,13,12,6,3,3), cledfu(13,12,9), cledir(13,12,9,3)
  !real*4 cether(13,12,2), xmiw(12,3)
  real visalb(istrip), niralb(istrip), thermk(istrip)

  ! Declare forcing variables
  ! -------------------------
  real tm(istrip), ps(istrip), sh(istrip), spdm(istrip)
  real ppl(istrip), ppc(istrip), swdown(istrip), lwdown(istrip)

  ! Declare other input variables
  ! -----------------------------
  real zs(istrip), bps(istrip), psb(istrip)
  real ros(istrip), thm(istrip), zb(istrip)
  real zfac, refzh

  ! Declare diagnostic variables
  ! ----------------------------
  real sgfg(istrip), sdens(istrip), satcap(istrip,2)
  real radt(istrip,2), etmass(istrip), hflux(istrip), roff(istrip)

  ! Declare ALMA variables
  ! ----------------------
  real swnet(istrip), lwnet(istrip), qle(istrip), qh(istrip)
  real qg(istrip), qf(istrip), qv(istrip), qtau(istrip), qa(istrip)
  real delsurfheat(istrip), delcoldcont(istrip), snowf(istrip)
  real rainf(istrip), evap(istrip), qs(istrip), qrec(istrip)
  real qsb(istrip), qsm(istrip), qfz(istrip), qst(istrip)
  real delsoilmoist(istrip), delswe(istrip), delintercept(istrip)
  real snowt(istrip), vegtc(istrip), baresoilt(istrip)
  real avgsurft(istrip), radteff(istrip), albedo(istrip)
  real swe(istrip), sweveg(istrip), soilmoist1(istrip)
  real soilmoist2(istrip), soilmoist3(istrip), soiltemp(istrip)
  real soilwet(istrip), potevap(istrip), ecanop(istrip)
  real tveg(istrip), esoil(istrip), rootmoist(istrip)
  real canopint(istrip), subsnow(istrip), subsurf(istrip)
  real acond(istrip), ccond(istrip), snowfrac(istrip)
  real snowdepth(istrip), sliqfrac(istrip)
  real ecanopt(istrip), ecanopi(istrip), eground(istrip)
  real hfluxc(istrip), hfluxg(istrip), watcan(istrip)
  real watgrd(istrip), snocan(istrip), snogrd(istrip)
  real wet1(istrip), wet2(istrip), wet3(istrip)

  ! Declare state variables
  ! -----------------------
  real tg(istrip), tc(istrip), tsn(istrip), td(istrip)
  real capac1(istrip), capac2(istrip), sn1(istrip), sn2(istrip)
  real gw1(istrip), gw2(istrip), gw3(istrip)

  ! Declare radiation variables
  ! ---------------------------
  real radc3(istrip,2), radn(istrip,3,2), extk(istrip,2,3,2)
  real cosz(istrip), salb(istrip,2,2)

  ! Declarations to get LIS to compile
  integer t
  real dtt, uwind, vwind

  ! Declarations to for diagnostic calculations
  real depth, factor, temp1, temp2, soilwm, soilww

  logical prin
  logical             :: alarmCheck
  integer             :: iret
  character*3   :: fnest

  write(fnest,'(i3.3)') n


!      if (first) then
!         print *,'setting HY-SSiB constants for the first time'
! Set constants

  grav = 9.81
  pie = 3.1415926
  cpair = 1004.16
  gasr = 287.04
  hlat = 2499000.
  snomel = 3.336242E+08
  clai = 840.0
  cw = 4200000.
  cs = 2106000.
  dkfac = 25.0
  epsfac = 0.622
  stefan = 5.67000E-08
  thres = 5.0E-03
  sskin = 4.0E-03
  tf = 273.15
  tf_snow = hyssib_struc(n)%rstemp
  tf_drain = tf - 1.0
  delsig = 0.012
      
      !Removed - now read in once in hyssibalb_module
      !open(unit=11,file=hyssib_struc(n)%afile,status='old',           &
      !             form='unformatted')
      !read(11) cedfu, cedir, cedfu1, cedir1, cedfu2, cedir2,           &
      !         cledfu, cledir, xmiu, cether, xmiw
      !close(11)
!         first = .false.
!      endif
  alarmCheck = LIS_isAlarmRinging(LIS_rc,"Hyssib model alarm "//trim(fnest))
  if(alarmCheck) then 
     
!=== Convert LIS Timestep varname to HY-SSiB timestep varname (DT) (sec)
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!        dtt = LIS_rc%ts
        dtt = LIS_rc%nts(n) ! EMK
        
!==== At this point start to assign HY-SSiB parameters to the respective variables
!*************************************
!========static parameters
         chil(1)     = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(1)
         topt(1)     = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(2)
         tll(1)      = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(3)
         tu(1)       = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(4)
         defac(1)    = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(5)
         ph1(1)      = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(6)
         ph2(1)      = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(7)
         rootd(1,1)  = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(8)
         rootd(1,2)  = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(9)
         rstpar(1,1) = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(10)
         rstpar(1,2) = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(11)
         rstpar(1,3) = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(12)
         phsat(1)    = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(13)
         poros(1)    = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(14)
         bee(1)      = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(15)
         satco(1)    = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(16)
         slope(1)    = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(17)
         zdepth(1,1) = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(18)
         zdepth(1,2) = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(19)
         zdepth(1,3) = HYSSIB_STRUC(N)%HYSSIB(T)%VEGP(20)
         
         index1       = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
         latco(1)    = LIS_DOMAIN(n)%grid(index1)%lat
         lonco(1)    = LIS_DOMAIN(n)%grid(index1)%lon
         wiltp(1)    = (-exp(ph2(1)) / phsat(1)) ** (-1 / bee(1))
         psilow(1)   = phsat(1) * (wiltp(1) ** (-bee(1)))
         topostd(1)  = HYSSIB_STRUC(N)%HYSSIB(T)%tempstd
         tdd(1)      = HYSSIB_STRUC(N)%HYSSIB(T)%tempbot
         visalb(1)   = LIS_alb(n)%albsf(LIS_domain(n)%gindex(&
              LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col, &
              LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row))
         niralb(1)   = LIS_alb(n)%albsf(LIS_domain(n)%gindex(&
              LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col, &
              LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row))
         zs(1)       = 0.0

!========monthly parameters
         z0(1)       = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(1)
         z1(1)       = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(2)
         z2(1)       = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(3)
         dd(1)       = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(4)
         vcover(1,1) = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(5)
         vcover(1,2) = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(6)
         zlt(1,1)    = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(7)
         zlt(1,2)    = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(8)
         green(1)    = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(9)
         rbc(1)      = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(10)
         rdc(1)      = HYSSIB_STRUC(N)%HYSSIB(T)%VEGIP(11)

         if (vcover(1,1).lt.0.0001) vcover(1,1) = 0.0001
         if (vcover(1,2).lt.0.0001) vcover(1,1) = 0.0001

!========state variables
         tc(1)       = HYSSIB_STRUC(N)%HYSSIB(T)%TC
         tg(1)       = HYSSIB_STRUC(N)%HYSSIB(T)%TG
         tsn(1)      = HYSSIB_STRUC(N)%HYSSIB(t)%TSN
         td(1)       = HYSSIB_STRUC(N)%HYSSIB(t)%TD
         gw1(1)      = HYSSIB_STRUC(N)%HYSSIB(t)%WWW(1)
         gw2(1)      = HYSSIB_STRUC(N)%HYSSIB(t)%WWW(2)
         gw3(1)      = HYSSIB_STRUC(N)%HYSSIB(t)%WWW(3)
         capac1(1)   = HYSSIB_STRUC(N)%HYSSIB(t)%CAPAC(1)
         capac2(1)   = HYSSIB_STRUC(N)%HYSSIB(t)%CAPAC(2)
         sn1(1)      = HYSSIB_STRUC(N)%HYSSIB(t)%SNOW(1)
         sn2(1)      = HYSSIB_STRUC(N)%HYSSIB(t)%SNOW(2)
         sgfg(1)     = HYSSIB_STRUC(N)%HYSSIB(t)%SGFG
         sdens(1)    = HYSSIB_STRUC(N)%HYSSIB(t)%SDENS

!========more parameters
         ityp(1)     = HYSSIB_STRUC(N)%HYSSIB(T)%VEGT
         month       = LIS_rc%mo
         nymd        = int((mod(LIS_rc%yr,100) * 10000) +              &
                           (LIS_rc%mo * 100) + LIS_rc%da)
         nhms        = int((LIS_rc%hr * 10000) + (LIS_rc%mn * 100) +   &
                            LIS_rc%ss)
         prin = .false.
! Amazon
!         if ((latco(1).eq.0.5).and.(lonco(1).eq.-59.5)) prin = .true.
! Valdai
!         if ((latco(1).eq.57.5).and.(lonco(1).eq.33.5)) prin = .true.
! New Zealand
!         if ((latco(1).eq.-37.5).and.(lonco(1).eq.175.5)) prin = .true.
! Greenland
!         if ((latco(1).eq.61.5).and.(lonco(1).eq.-44.5)) prin = .true.
! First point
!         if ((latco(1).eq.-54.5).and.(lonco(1).eq.-70.5)) prin = .true.
! First point for RHONE
!         if ((latco(1).eq.43.5).and.(lonco(1).eq.4.5)) goto 69
! Russia
!         if ((latco(1).eq.73.5).and.(lonco(1).eq.95.5)) prin = .true.
! Lat 40 and Lon -100
!         if ((latco(1).eq.40.125).and.(lonco(1).eq.-99.875)) prin = .true.

!======end of assigning parameters

!=== THE FOLLOWING BREAKS DOWN THE FORCING VARIABLES

         tm(1)     = hyssib_struc(n)%HYSSIB(t)%tair/hyssib_struc(n)%forc_count
         sh(1)     = hyssib_struc(n)%HYSSIB(t)%qair/hyssib_struc(n)%forc_count
         swdown(1) = amax1(&
              hyssib_struc(n)%HYSSIB(t)%swdown/hyssib_struc(n)%forc_count,0.0)
         lwdown(1) = hyssib_struc(n)%HYSSIB(t)%lwdown/hyssib_struc(n)%forc_count
         uwind     = (hyssib_struc(n)%HYSSIB(t)%uwind) *               &
                     (hyssib_struc(n)%HYSSIB(t)%uwind)
         vwind     = (hyssib_struc(n)%HYSSIB(t)%vwind) *               &
                     (hyssib_struc(n)%HYSSIB(t)%vwind)
         spdm(1)   = amax1(sqrt(uwind + vwind),0.25)/hyssib_struc(n)%forc_count
         ps(1)     = (hyssib_struc(n)%HYSSIB(t)%psurf / 100.0) &
              /hyssib_struc(n)%forc_count
         ppl(1)    = hyssib_struc(n)%HYSSIB(t)%rainf_in/hyssib_struc(n)%forc_count -              &
                     hyssib_struc(n)%HYSSIB(t)%rainf_cp/hyssib_struc(n)%forc_count
         ppc(1)    = hyssib_struc(n)%HYSSIB(t)%rainf_cp/hyssib_struc(n)%forc_count

         if ((hyssib_struc(n)%zh.lt.0.0).or.                           &
             (hyssib_struc(n)%zm.lt.0.0)) then
            hyssib_struc(n)%zh = hyssib_struc(n)%hyssib(t)%fheight
            hyssib_struc(n)%zm = hyssib_struc(n)%hyssib(t)%fheight
            if (hyssib_struc(n)%hyssib(t)%fheight.le.0.0) then
               write(LIS_logunit,*) 'Forcing height less than or'
               write(LIS_logunit,*) 'equal to zero!  Stopping run'
               call LIS_endrun()
            endif
         else
            if (hyssib_struc(n)%zh.lt.hyssib_struc(n)%zm) then
               zfac = LOG((hyssib_struc(n)%zh+z0(1))/z0(1)) /          &
                      LOG((hyssib_struc(n)%zm+z0(1))/z0(1))
               spdm(1) = amax1(zfac*spdm(1),0.25)
            endif
         endif
         refzh = hyssib_struc(n)%zh

         if (prin) then
            print *,' '
            print *,'Latitude ',latco(1),'; Longitude ',lonco(1)
            print *,' '
            print *,'Constants: '
            print *,'vegtype: ',ityp(1)
            print *,' '
            print *,'Forcing data: '
            print *,'tm: ',tm(1)
            print *,'ps: ',ps(1)
            print *,'sh: ',sh(1)
            print *,'spdm: ',spdm(1)
            print *,'ppl: ',ppl(1)
            print *,'ppc: ',ppc(1)
            print *,'swdown: ',swdown(1)
            print *,'lwdown: ',lwdown(1)
            print *,' '
         endif

         if (prin) then
            print *,'Visible albedo: ',visalb(1)
            print *,'Near-IR albedo: ',niralb(1)
            print *,' '
         endif

!===Enter here any final adjustment of variables before calling the LSM
         bps(1) = (ps(1) / 1000.) ** (gasr / cpair)
         thm(1) = tm(1) / bps(1)
         ros(1) = (ps(1) * 100.0) / (gasr * tm(1))
         ppl(1) = ppl(1) * dtt
         ppc(1) = ppc(1) * dtt
         psb(1) = ps(1) * delsig
         zb(1) = zs(1) + (100.0 * psb(1) / (ros(1) * grav))

!===End of adjustments
!
!=== ===============================================================
!=== ===============================================================
!=== ============BEGIN OF SIMULATION===============================

! Call radiation
         call hyssib_radsplit(istrip,tf,nymd,nhms,dtt,month,lonco,     &
              latco,ityp,hyssibalb(n)%xmiu,hyssibalb(n)%cledir,        &
              hyssibalb(n)%cledfu,hyssibalb(n)%cedir,                  &
              hyssibalb(n)%cedfu,hyssibalb(n)%cedir2,                  &
              hyssibalb(n)%cedfu2,hyssibalb(n)%cedir1,                 & 
              hyssibalb(n)%cedfu1,hyssibalb(n)%cether,                 &
              z1,z2,zlt,vcover,visalb,niralb,tc,tg,tsn,sgfg,           &
              sn1,sn2,satcap,sdens,swdown,lwdown,snowfrac,             &
              thermk,extk,cosz,salb,radn,radc3,prin)

! Call the LSM
         call hyssib_routines(istrip,grav,refzh,pie,dtt,cpair,gasr,    &
              hlat,snomel,clai,cw,cs,dkfac,epsfac,stefan,thres,sskin,  &
              tf,tf_snow,tf_drain,delsig,lonco,latco,ityp,z2,dd,z0,    &
              phsat,poros,bee,satco,slope,zdepth,vcover,sgfg,sdens,    &
              thermk,extk,cosz,zlt,green,chil,rbc,rdc,topt,tll,tu,     &
              defac,ph1,ph2,rootd,rstpar,satcap,wiltp,psilow,topostd,  &
              tdd,bps,psb,ros,thm,zb,tm,ps,sh,spdm,tg,capac1,sn1,      &
              capac2,sn2,gw1,gw2,gw3,tc,td,tsn,ppl,ppc,radn,radc3,     &
              radt,etmass,hflux,roff,swnet,lwnet,qle,qh,qg,qf,qv,qtau, &
              qa,delsurfheat,delcoldcont,snowf,rainf,evap,qs,qrec,qsb, &
              qsm,qfz,qst,delsoilmoist,delswe,delintercept,snowt,      &
              vegtc,baresoilt,avgsurft,radteff,swe,sweveg,soilmoist1,  &
              soilmoist2,soilmoist3,soiltemp,soilwet,potevap,ecanop,   &
              tveg,esoil,rootmoist,canopint,subsnow,subsurf,acond,     &
              ccond,snowdepth,sliqfrac,ecanopt,ecanopi,eground,hfluxc, &
              hfluxg,watcan,watgrd,snocan,snogrd,wet1,wet2,wet3,prin)

!=== ===============END OF SIMULATION===============================
!=== ===============================================================
!=== ===============================================================

!========Collect state variables
         HYSSIB_STRUC(N)%HYSSIB(T)%TC       = tc(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%TG       = tg(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%TSN      = tsn(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%TD       = td(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%WWW(1)   = gw1(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%WWW(2)   = gw2(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%WWW(3)   = gw3(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%CAPAC(1) = capac1(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%CAPAC(2) = capac2(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%SNOW(1)  = sn1(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%SNOW(2)  = sn2(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%SGFG     = sgfg(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%SDENS    = sdens(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%GREEN    = green(1)

         !Calculate MSTAVRZ (in excess of wilting point)
         soilwm=0.0
         soilww=0.0
         soilwm=(poros(1)-(poros(1)*wiltp(1)))*zdepth(1,1)
         soilwm=soilwm+((poros(1) - (poros(1)*wiltp(1)))*zdepth(1,2))
         soilww=((hyssib_struc(n)%hyssib(t)%www(1)*poros(1)) -         &
                 (poros(1)*wiltp(1)))*poros(1)
         soilww=soilww+((hyssib_struc(n)%hyssib(t)%www(2)*poros(1)) -  &
                        (poros(1)*wiltp(1)))*zdepth(1,2)
         hyssib_struc(n)%hyssib(t)%soilwetrz = soilww/soilwm

         !Calculate soil moisture (kg/m2) in top 1meter layer of soil
         !If soil is less than 1 meter deep, scale the soil moisture  
         !up so that it has an effective depth of 1 meter of soil
         DEPTH=0.0
         FLAG=0
         DEPTH=DEPTH+ZDEPTH(1,1)
         IF (DEPTH.GE.1) THEN
            FACTOR=(1.0-((DEPTH-1.0)/DEPTH))
            IF (DEPTH.EQ.1.0) THEN
               hyssib_struc(n)%hyssib(t)%soilmoist1m=hyssib_struc(n)%hyssib(t)%soilmoist(1)
               FLAG=1
            ELSE
               hyssib_struc(n)%hyssib(t)%soilmoist1m=hyssib_struc(n)%hyssib(t)%soilmoist(1)*FACTOR
               FLAG=1
            ENDIF
         ENDIF
         IF (FLAG.EQ.0) THEN
            DEPTH=DEPTH+ZDEPTH(1,2)
            IF (DEPTH.GE.1.0) THEN
               IF (DEPTH.EQ.1.0) THEN
                  hyssib_struc(n)%hyssib(t)%soilmoist1m =hyssib_struc(n)%hyssib(t)%soilmoist(1) + &
                  hyssib_struc(n)%hyssib(t)%soilmoist(2)
                  FLAG=1
               ELSE
                  TEMP1=DEPTH-1.0
                  FACTOR=(1.0-(TEMP1/ZDEPTH(1,2)))
                  TEMP2=hyssib_struc(n)%hyssib(t)%soilmoist(2)*FACTOR
                  hyssib_struc(n)%hyssib(t)%soilmoist1m=hyssib_struc(n)%hyssib(t)%soilmoist(1)+TEMP2
                  FLAG=1
               ENDIF
            ENDIF
         ENDIF
         IF (FLAG.EQ.0) THEN
            DEPTH=DEPTH+ZDEPTH(1,3)
            IF (DEPTH.GE.1.0) THEN
               IF (DEPTH.EQ.1.0) THEN
                  hyssib_struc(n)%hyssib(t)%soilmoist1m = hyssib_struc(n)%hyssib(t)%soilmoist(1) + &
                  hyssib_struc(n)%hyssib(t)%soilmoist(2) + hyssib_struc(n)%hyssib(t)%soilmoist(3)
                  FLAG=1
               ELSE
                  TEMP1=DEPTH-1.0
                  FACTOR=(1.0-(TEMP1/ZDEPTH(1,3)))
                  TEMP2=hyssib_struc(n)%hyssib(t)%soilmoist(3)*FACTOR
                  hyssib_struc(n)%hyssib(t)%soilmoist1m = hyssib_struc(n)%hyssib(t)%soilmoist(1) + &
                  hyssib_struc(n)%hyssib(t)%soilmoist(2) + TEMP2
                  FLAG=1
               ENDIF
            ELSE
               hyssib_struc(n)%hyssib(t)%soilmoist1m=(hyssib_struc(n)%hyssib(t)%soilmoist(1)+ &
               hyssib_struc(n)%hyssib(t)%soilmoist(2) + hyssib_struc(n)%hyssib(t)%soilmoist(3))/DEPTH
               FLAG=1
            ENDIF
          ENDIF

!========Collect the ALMA output variables
         HYSSIB_STRUC(N)%HYSSIB(T)%swnet        =  swnet(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%lwnet        =  lwnet(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qle          =  qle(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qh           =  qh(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qg           =  qg(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qf           =  qf(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qv           =  qv(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qtau         =  qtau(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qa           =  qa(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%delsurfheat  =  delsurfheat(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%delcoldcont  =  delcoldcont(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%snowf        =  snowf(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%rainf        =  rainf(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%evap         =  evap(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qs           =  qs(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qrec         =  qrec(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qsb          =  qsb(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qsm          =  qsm(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qfz          =  qfz(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%qst          =  qst(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%delsoilmoist =  delsoilmoist(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%delswe       =  delswe(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%delintercept =  delintercept(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%swdown_out   =  swdown(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%lwdown_out   =  lwdown(1)

         if (swdown(1).gt.0.0) then
            albedo(1) = ((radn(1,1,1) * salb(1,1,1)) +                 &
                         (radn(1,1,2) * salb(1,1,2)) +                 &
                         (radn(1,2,1) * salb(1,2,1)) +                 &
                         (radn(1,2,2) * salb(1,2,2))) / swdown(1)
         else
            albedo(1) = LIS_RC%UDEF
         endif
         if (sgfg(1).gt.0.5) then
            HYSSIB_STRUC(N)%HYSSIB(T)%snowt     = snowt(1)
         else
            HYSSIB_STRUC(N)%HYSSIB(T)%snowt     = LIS_RC%UDEF
         endif
         HYSSIB_STRUC(N)%HYSSIB(T)%vegtc        = vegtc(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%baresoilt    = baresoilt(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%avgsurft     = avgsurft(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%radteff      = radteff(1)
         if (swdown(1).gt.0.0) then
            HYSSIB_STRUC(N)%HYSSIB(T)%albedo    = albedo(1)
         else
            HYSSIB_STRUC(N)%HYSSIB(T)%albedo    = LIS_RC%UDEF
         endif
         HYSSIB_STRUC(N)%HYSSIB(T)%swe          = swe(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%sweveg       = sweveg(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%soilmoist(1) = soilmoist1(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%soilmoist(2) = soilmoist2(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%soilmoist(3) = soilmoist3(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%soiltemp     = soiltemp(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%soilwet      = soilwet(1)
         if ((sn2(1)+capac2(1)).gt.0.0) then
            HYSSIB_STRUC(N)%HYSSIB(T)%sliqfrac  = sliqfrac(1)
         else
            HYSSIB_STRUC(N)%HYSSIB(T)%sliqfrac  = LIS_RC%UDEF
         endif

         HYSSIB_STRUC(N)%HYSSIB(T)%potevap      = potevap(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%ecanop       = ecanop(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%tveg         = tveg(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%esoil        = esoil(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%rootmoist    = rootmoist(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%canopint     = canopint(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%subsnow      = subsnow(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%subsurf      = subsurf(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%acond        = acond(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%ccond        = ccond(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%snowfrac     = snowfrac(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%snowdepth    = snowdepth(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%ect          = ecanopt(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%eci          = ecanopi(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%egs          = eground(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%hc           = hfluxc(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%hg           = hfluxg(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%radnvisdir   = radn(1,1,1)
         HYSSIB_STRUC(N)%HYSSIB(T)%radnvisdif   = radn(1,1,2)
         HYSSIB_STRUC(N)%HYSSIB(T)%radnnirdir   = radn(1,2,1)
         HYSSIB_STRUC(N)%HYSSIB(T)%radnnirdif   = radn(1,2,2)
         HYSSIB_STRUC(N)%HYSSIB(T)%salbvisdir   = salb(1,1,1)
         HYSSIB_STRUC(N)%HYSSIB(T)%salbvisdif   = salb(1,1,2)
         HYSSIB_STRUC(N)%HYSSIB(T)%salbnirdir   = salb(1,2,1)
         HYSSIB_STRUC(N)%HYSSIB(T)%salbnirdif   = salb(1,2,2)
         HYSSIB_STRUC(N)%HYSSIB(T)%radtc        = radt(1,1)
         HYSSIB_STRUC(N)%HYSSIB(T)%radtg        = radt(1,2)
         HYSSIB_STRUC(N)%HYSSIB(T)%watcan       = watcan(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%watgrd       = watgrd(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%snocan       = snocan(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%snogrd       = snogrd(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%wet1         = wet1(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%wet2         = wet2(1)
         HYSSIB_STRUC(N)%HYSSIB(T)%wet3         = wet3(1)

         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWNET,value=           &
              hyssib_struc(n)%hyssib(t)%swnet,vlevel=1,unit="W m-2",&
              direction="DN",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LWNET,value=           &
              hyssib_struc(n)%hyssib(t)%lwnet,vlevel=1,unit="W m-2",&
              direction="DN",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QLE,value=             &
              hyssib_struc(n)%hyssib(t)%qle,vlevel=1,unit="W m-2",&
              direction="UP",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QH,value=              &
              hyssib_struc(n)%hyssib(t)%qh,vlevel=1,unit="W m-2",&
              direction="UP",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QG,value=              &
              hyssib_struc(n)%hyssib(t)%qg,vlevel=1,unit="W m-2",&
              direction="DN",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QV,value=              &
              hyssib_struc(n)%hyssib(t)%qv,vlevel=1,unit="W m-2",&
              direction="S2V",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QA,value=              &
              hyssib_struc(n)%hyssib(t)%qa,vlevel=1,unit="W m-2",&
              direction="DN",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QF,value=              &
              hyssib_struc(n)%hyssib(t)%qf,vlevel=1,unit="W m-2",&
              direction="S2L",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QTAU,value=            &
              hyssib_struc(n)%hyssib(t)%qtau,vlevel=1,unit="N m-2",&
              direction="DN",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSURFHEAT,value=     &
              hyssib_struc(n)%hyssib(t)%delsurfheat,vlevel=1,unit=     &
              "J m-2",direction="INC",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELCOLDCONT,value=     &
              hyssib_struc(n)%hyssib(t)%delcoldcont,vlevel=1,unit=     &
              "J m-2",direction="INC",surface_type=LIS_rc%lsm_index)
         
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=           &
              hyssib_struc(n)%hyssib(t)%snowf,vlevel=1,unit="kg m-2 s-1",&
              direction="DN",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=           &
              hyssib_struc(n)%hyssib(t)%rainf,vlevel=1,unit="kg m-2 s-1",&
              direction="DN",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_EVAP,value=            &
              hyssib_struc(n)%hyssib(t)%evap,vlevel=1,unit="kg m-2 s-1",&
              direction="UP",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QS,value=              &
              hyssib_struc(n)%hyssib(t)%qs,vlevel=1,unit="kg m-2 s-1",&
              direction="OUT",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSB,value=             &
              hyssib_struc(n)%hyssib(t)%qsb,vlevel=1,unit="kg m-2 s-1",&
              direction="OUT",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QREC,value=            &
              hyssib_struc(n)%hyssib(t)%qrec,vlevel=1,unit="kg m-2 s-1",&
              direction="IN",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QFZ,value=             &
              hyssib_struc(n)%hyssib(t)%qfz,vlevel=1,unit="kg m-2 s-1",&
              direction="L2S",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSM,value=             &
              hyssib_struc(n)%hyssib(t)%qsm,vlevel=1,unit="kg m-2 s-1",&
              direction="S2L",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QST,value=             &
              hyssib_struc(n)%hyssib(t)%qst,vlevel=1,unit="kg m-2 s-1",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSOILMOIST,value=    &
              hyssib_struc(n)%hyssib(t)%delsoilmoist,vlevel=1,unit=    &
              "kg m-2",direction="INC",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSWE,value=          &
              hyssib_struc(n)%hyssib(t)%delswe,vlevel=1,unit="kg m-2",&
              direction="INC",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELINTERCEPT,value=    &
              hyssib_struc(n)%hyssib(t)%delintercept,vlevel=1,unit=    &
              "kg m-2",direction="INC",surface_type=LIS_rc%lsm_index)
         
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWT,value=           &
              hyssib_struc(n)%hyssib(t)%snowt,vlevel=1,unit="K",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_VEGT,value=            &
              hyssib_struc(n)%hyssib(t)%vegtc,vlevel=1,unit="K",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_BARESOILT,value=       &
              hyssib_struc(n)%hyssib(t)%baresoilt,vlevel=1,unit="K",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_AVGSURFT,value=        &
              hyssib_struc(n)%hyssib(t)%avgsurft,vlevel=1,unit="K",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RADT,value=            &
              hyssib_struc(n)%hyssib(t)%radteff,vlevel=1,unit="K",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ALBEDO,value=          &
              hyssib_struc(n)%hyssib(t)%albedo,vlevel=1,unit="-",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWE,value=             &
              hyssib_struc(n)%hyssib(t)%swe,vlevel=1,unit="kg m-2",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWEVEG,value=          &
              hyssib_struc(n)%hyssib(t)%sweveg,vlevel=1,unit="kg m-2",&
              direction="-",surface_type=LIS_rc%lsm_index)
         
         do i=1,3
            call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILMOIST,vlevel=i, &
                value=hyssib_struc(n)%hyssib(t)%soilmoist(i),unit=     &
                "kg m-2",direction="-",surface_type=LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILMOIST,vlevel=i, &
                value=hyssib_struc(n)%hyssib(t)%soilmoist(i),unit=     &
                "m^3 m-3",direction="-",surface_type=LIS_rc%lsm_index)
         enddo
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILTEMP,vlevel=1,     &
             value=hyssib_struc(n)%hyssib(t)%soiltemp,unit="K",&
             direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILWET,value=         &
              hyssib_struc(n)%hyssib(t)%soilwet,vlevel=1,unit="-",&
              direction="-",surface_type=LIS_rc%lsm_index)
         
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POTEVAP,value=         &
              hyssib_struc(n)%hyssib(t)%potevap,vlevel=1,unit="kg m-2 s-1",&
              direction="UP",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ECANOP,value=          &
              hyssib_struc(n)%hyssib(t)%ecanop,vlevel=1,unit="kg m-2 s-1",&
              direction="UP",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TVEG,value=            &
              hyssib_struc(n)%hyssib(t)%tveg,vlevel=1,unit="kg m-2 s-1",&
              direction="UP",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ESOIL,value=           &
              hyssib_struc(n)%hyssib(t)%esoil,vlevel=1,unit="kg m-2 s-1",&
              direction="UP",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ROOTMOIST,value=       &
              hyssib_struc(n)%hyssib(t)%rootmoist,vlevel=1,unit="kg m-2",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CANOPINT,value=        &
              hyssib_struc(n)%hyssib(t)%canopint,vlevel=1,unit="kg m-2",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SUBSNOW,value=         &
              hyssib_struc(n)%hyssib(t)%subsnow,vlevel=1,unit="kg m-2 s-1",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SUBSURF,value=         &
              hyssib_struc(n)%hyssib(t)%subsurf,vlevel=1,unit="kg m-2 s-1",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ACOND,value=           &
              hyssib_struc(n)%hyssib(t)%acond,vlevel=1,unit="m s-1",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CCOND,value=           &
              hyssib_struc(n)%hyssib(t)%ccond,vlevel=1,unit="m s-1",&
              direction="-",surface_type=LIS_rc%lsm_index)
         
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWCOVER,value=       &
              hyssib_struc(n)%hyssib(t)%snowfrac,vlevel=1,unit="-",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWDEPTH,value=       &
              hyssib_struc(n)%hyssib(t)%snowdepth,vlevel=1,unit="m",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SLIQFRAC,value=        &
              hyssib_struc(n)%hyssib(t)%sliqfrac,vlevel=1,unit="-",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_GREENNESS,value=       &
              hyssib_struc(n)%hyssib(t)%green,vlevel=1,unit="-",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_GREENNESS,value=       &
              hyssib_struc(n)%hyssib(t)%green*100.0,vlevel=1,unit="%",&
              direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ROUGHNESS,value=       &
              z0(1),vlevel=1,unit="m",direction="-",surface_type=LIS_rc%lsm_index)
         call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LAI,value=             &
              zlt(1,1),vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)

         !      if (prin) then
         !         continue
         !      else
         !         HYSSIB_STRUC(N)%HYSSIB(T)%swnet = LIS_rc%udef
         !      endif

!reset 
        hyssib_struc(n)%hyssib(t)%tair = 0 
        hyssib_struc(n)%hyssib(t)%qair = 0 
        hyssib_struc(n)%hyssib(t)%swdown = 0 
        hyssib_struc(n)%hyssib(t)%lwdown = 0
        hyssib_struc(n)%hyssib(t)%uwind = 0
        hyssib_struc(n)%hyssib(t)%vwind = 0 
        hyssib_struc(n)%hyssib(t)%psurf = 0 
        hyssib_struc(n)%hyssib(t)%rainf_in = 0 
        hyssib_struc(n)%hyssib(t)%rainf_cp = 0 
      enddo
!>>> end of hyssib_main <<<<<<<<
      hyssib_struc(n)%forc_count = 0 
   endif
   
 end subroutine hyssib_main

!*** HY-SSiB SUBROUTINES ***********************************************

      subroutine hyssib_radsplit(istrip,tf,nymd,nhms,dtc3x,month,lonco,&
           latco,ityp,xmiu,cledir,cledfu,cedir,cedfu,cedir2,cedfu2,    &
           cedir1,cedfu1,cether,z1,z2,zlt,vcover,visalb,niralb,tc,     &
           tg,tsn,sgfg,sn1,sn2,satcap,sdens,swdown,lwdown,areas,       &
           thermk,extk,cosz,salb,radn,radc3,prin)

      implicit none
!***********************************************************************
!***  split swdown radiation into 4 parts; calculate/use albedos
!-----------------------------------------------------------------------
!     *** indices ***
!-----------------------------------------------------------------------
! cg =1...canopy
! cg =2...ground cover
! vn =1...visible      (0.0-0.7 micron)
! vn =2...near-infrared(0.7-3.0 micron)
! bd =1...beam
! bd =2...diffuse
! ld =1...live leaves
! ld =2...dead leaves
! vnt=1...visible      (0.0-0.7 micron)
! vnt=2...near-infrared(0.7-3.0 micron)
! vnt=3...thermal
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! tf           (K)       freezing temperature
! nymd,nhms    (-)       date/time
! dtc3x        (sec)     timestep
! month        (1-12)    current month
! lonco        (deg)     longitude (-180 west to 180 east)
! latco        (deg)     latitude (-90 south to 90 north)
! ityp         (1-13)    vegetation index
! xmiu-cether  (-)       coefficients of radiation absorption & albedo
! z1           (m)       canopy base height
! z2           (m)       canopy top height
! zlt(cg)      (-)       leaf area index
! vcover(cg)   (0-1)     fraction of vegetation cover [fpar/greenness]
! visalb       (0-1)     visible albedo
! niralb       (0-1)     near-ir albedo
! tc           (K)       canopy temperature
! tg           (K)       ground temperature
! tsn          (K)       snow on ground temperature
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
! sn1,sn2      (mm)      liquid equiv. snow on canopy/ground
! sdens        (kg m-3)  density of bulk snow layer
! swdown       (W m-2)   forcing sw downward radiation at the surface
! lwdown       (W m-2)   forcing lw downward radiation at the surface
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! satcap(cg)   (m)       interception capacity
! areas        (0-1)     snow fraction on ground
! thermk       (0-1)     canopy emissivity
! extk(cg,vnt,bd)  (-)   extinction coefficient
! cosz         (radians) cosine of solar zenith angle
! salb(vn,bd)  (0-1)     four components of surface albedo
! radn(vnt,bd) (W m-2)   downward sw/lw radiation at the surface
! radc3(cg)    (W m-2)   absorbed sw    radiation at the surface
!-----------------------------------------------------------------------
      integer istrip, nymd, nhms, month, ityp(istrip)
      integer iveg, iwave, irad, i
      real lonco(istrip), latco(istrip), dtc3x
      real tf, xvisb, xvisd, xnirb, xnird, rfact
      real*4 cedfu(13,12,9), cedir(13,12,9,3), cedfu1(2,13,12,6,3)
      real*4 cedir1(2,13,12,6,3,3), cedfu2(2,13,12,6,3), xmiu(12,3)
      real*4 cedir2(2,13,12,6,3,3), cledfu(13,12,9), cledir(13,12,9,3)
      real*4 cether(13,12,2)
      real z1(istrip), z2(istrip), zlt(istrip,2)
      real vcover(istrip,2), visalb(istrip), niralb(istrip)
      real tc(istrip), tg(istrip), tsn(istrip)
      real sgfg(istrip), sn1(istrip), sn2(istrip)
      real satcap(istrip,2), sdens(istrip)
      real swdown(istrip), lwdown(istrip)
      real areas(istrip), thermk(istrip)
      real extk(istrip,2,3,2), cosz(istrip)
      real salb(istrip,2,2), radn(istrip,3,2), radc3(istrip,2)
      real radfac(istrip,2,2,2), snoww(istrip,2)
      logical prin

      call swtop(istrip,nymd,nhms,dtc3x,lonco,latco,cosz)

      do i = 1,istrip
         call radspec(swdown(i),cosz(i),xvisb,xvisd,xnirb,xnird)
         radn(i,1,1) = xvisb
         radn(i,1,2) = xvisd
         radn(i,2,1) = xnirb
         radn(i,2,2) = xnird
         radn(i,3,1) = 0.0
         radn(i,3,2) = lwdown(i)
         radc3(i,1) = 0.0
         radc3(i,2) = 0.0
         snoww(i,1) = sn1(i) * 0.001
         snoww(i,2) = sn2(i) * 0.001
         if (prin) then
            print *,'Radsplit: '
            print *,'radn_11: ',radn(i,1,1)
            print *,'radn_12: ',radn(i,1,2)
            print *,'radn_21: ',radn(i,2,1)
            print *,'radn_22: ',radn(i,2,2)
            print *,'cosz: ',cosz(i)
            print *,' '
         endif
      enddo

      call intcap(istrip,zlt,z1,z2,snoww,satcap)

      call radalb_hyssib(istrip,month,tf,ityp,snoww,areas,satcap,      &
           tg,cosz,xmiu,cledir,cledfu,cedir,cedfu,cedir2,cedfu2,       &
           cedir1,cedfu1,cether,latco,tc,tsn,sgfg,sdens,zlt,vcover,    &
           thermk,extk,visalb,niralb,radfac,salb,prin)

      do iveg = 1,2
         do iwave = 1,2
            do irad = 1,2
               do i = 1,istrip
                  rfact = radfac(i,iveg,iwave,irad) * (1.0 -           &
                       salb(i,iwave,irad))
                  radc3(i,iveg) = rfact * radn(i,iwave,irad) +         &
                       radc3(i,iveg)
               enddo
            enddo
         enddo
      enddo

      if (prin) then
         print *,'radc3: ',radc3(1,1),radc3(1,2)
         print *,'salb: ',salb(1,1,1),salb(1,1,2),salb(1,2,1),salb(1,2,2)
         print *,'radn: ',radn(1,1,1),radn(1,1,2),radn(1,2,1),radn(1,2,2)
         print *,'radfac_1_X_X: ',radfac(1,1,1,1),radfac(1,1,1,2), &
                                  radfac(1,1,2,1),radfac(1,1,2,2)
         print *,'radfac_2_X_X: ',radfac(1,2,1,1),radfac(1,2,1,2), &
                                  radfac(1,2,2,1),radfac(1,2,2,2)
      endif

      return
      end

!***********************************************************************
      subroutine radalb_hyssib(istrip,mont,tf,ityp,snoww,scov2,satcap, &
           tg,cosz,xmiu,cledir,cledfu,cedir,cedfu,cedir2,cedfu2,       &
           cedir1,cedfu1,cether,latco,tc,tsn,sgfg,sdens,zlt,vcover,    &
           thermk,extk,visalb,niralb,radfac,salb,prin)

      implicit none
!***********************************************************************
!***  surface albedos via two stream approximation (direct and diffuse)
!     reference  :canopy reflectance,photosynthesis and transpiration
!                 p.j.sellers  (85) int.j.remote sensing,vol6,no.8,1135
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! month        (1-12)    current month
! tf           (K)       freezing temperature
! ityp         (1-13)    vegetation index
! snoww(cg)    (m)       liquid equiv. snow on canopy/ground
! scov2        (0-1)     snow fraction on ground
! satcap(cg)   (m)       interception capacity
! tg           (K)       ground temperature
! cosz         (radians) cosine of solar zenith angle
! xmiu-cether  (-)       coefficients of radiation absorption & albedo
! latco        (deg)     latitude (-90 south to 90 north)
! tc           (K)       canopy temperature
! tsn          (K)       snow on ground temperature
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
! sdens        (kg m-3)  density of bulk snow layer
! zlt(cg)      (-)       leaf area index
! vcover(cg)   (0-1)     fraction of vegetation cover [fpar/greenness]
! visalb       (0-1)     visible albedo
! niralb       (0-1)     near-ir albedo
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! thermk       (0-1)     canopy emissivity
! extk(cg,vnt,bd)  (-)   extinction coefficient
! radfac(cg,vn,bd) (0-1) fraction of downward solar radiation at surface
! salb(vn,bd)  (0-1)     four components of surface albedo
!-----------------------------------------------------------------------
      integer istrip, mont, mon, nk, ik, i, jj, i1, ml, k1, k2, ntyp
      integer ityp(istrip)
      real alphafac, betafac, pluscons, exppp
      real tf, xf, xf2, sc1, sc2, xm1, xm2, xtm1, xtm2
      real snoww(istrip,2), scov2(istrip), satcap(istrip,2)
      real tg(istrip), cosz(istrip)
      real*4 cedfu(13,12,9), cedir(13,12,9,3), cedfu1(2,13,12,6,3)
      real*4 cedir1(2,13,12,6,3,3), cedfu2(2,13,12,6,3), xmiu(12,3)
      real*4 cedir2(2,13,12,6,3,3), cledfu(13,12,9), cledir(13,12,9,3)
      real*4 cether(13,12,2)
      real latco(istrip), tc(istrip), tsn(istrip), sgfg(istrip)
      real sdens(istrip), zlt(istrip,2), vcover(istrip,2)
      real visalb(istrip), niralb(istrip)
      real thermk(istrip), extk(istrip,2,3,2)
      real radfac(istrip,2,2,2), salb(istrip,2,2)
      real deltg(istrip), zkat(istrip), areas(istrip)
      real f(istrip), fmelt(istrip), scov(istrip)
      real temp(18), xmi1(12,3), temp_bare(18)
      logical prin

      nk = 3
      alphafac = 3.5e4
      betafac = 2.0
! SMP3, Eq. 4
      pluscons = -(1.0-exp(2.0 * betafac)) / (1.0+exp(2.0 * betafac))

! Calculate fractional snow on canopy/ground
      do i = 1,istrip
         fmelt(i) = 1.0
         deltg(i) = tf - max((tsn(i)*sgfg(i)),tg(i))
         if ((abs(deltg(i)).lt.1.0).and.(deltg(i).gt.0.0)) then
            fmelt(i) = 0.6
         endif
         if (tc(i).le.tf) then
            scov(i) = max(min((snoww(i,1) / satcap(i,1)),0.5),0.0)
         else
            scov(i) = 0.0
         endif
         if (snoww(i,2).lt.1e-5) then
            areas(i) = 0.0
         else
! SMP3, Eq. 3
            exppp = exp(-2.0 * ((alphafac * snoww(i,2) /               &
                    max(sdens(i),100.0)) - betafac))
! SMP3, Eq. 5
            areas(i) = min(1.0, ((((1.0 - exppp) / (1.0 + exppp)) +    &
                       pluscons) / (1.0 + pluscons)))
         endif
         scov2(i) = max(min(areas(i),1.0),0.0)
      enddo

      if (prin) then
         print *,'scov: ',scov(1)
         print *,'snow1: ',snoww(1,1)
         print *,'satcap: ',satcap(1,1)
         print *,'tc: ',tc(1)
         print *,'tf: ',tf
         print *,'scov2: ',scov2(1)
         print *,'areas: ',areas(1)
         print *,'snow2: ',snoww(1,2)
         print *,'sdens: ',sdens(1)
         print *,'fmelt: ',fmelt(1)
         print *,'deltg: ',deltg(1)
         print *,'latco: ',latco(1)
         print *,'cosz: ',cosz(1)
         print *,'ityp: ',ityp(1)
      endif

! Terms which multiply incoming short wave fluxes
! to give absorption of radiation by canopy and ground
      do i = 1,istrip
         if (latco(i).ge.0.0) then
            mon = mont
         else
            mon = mont - 6
            if (mon.le.0) mon = mon + 12
         endif
         f(i) = max(cosz(i), 0.01746)
         xf = f(i)
         xf2 = xf * xf
         ntyp = ityp(i)
         if (ntyp.eq.13) ntyp = 11

! Set xmi1
         do jj = 1,nk
            xmi1(mon,jj) = xmiu(mon,jj)
         enddo

         if ((scov(i).lt.0.025).and.(scov2(i).lt.0.025)) then
            do i1 = 1,9
               temp(i1) = cledir(ntyp,mon,i1,1)       +                &
                          cledir(ntyp,mon,i1,2) * xf  +                &
                          cledir(ntyp,mon,i1,3) * xf2
               temp(i1+9) = cledfu(ntyp,mon,i1)
! Save values to scale albedo
               temp_bare(i1) = temp(i1)
               temp_bare(i1+9) = temp(i1+9)
               if ((prin).and.(i1.eq.1)) then
                  print *,'mon: ',mon
                  print *,'xf: ',xf
                  print *,'xf2: ',xf2
                  print *,'cledir1: ',cledir(ntyp,mon,i1,1)
                  print *,'cledir2: ',cledir(ntyp,mon,i1,2)
                  print *,'cledir3: ',cledir(ntyp,mon,i1,3)
                  print *,'temp',i1,': ',temp(i1)
                  print *,'temp',i1+9,': ',temp(i1+9)
                  print *,'tempbare',i1,': ',temp_bare(i1)
                  print *,'tempbare',i1+9,': ',temp_bare(i1+9)
               endif
            enddo
            goto 500
         else
! Save values to scale albedo if there is snow
            do i1 = 1,9
               temp_bare(i1) = cledir(ntyp,mon,i1,1)       +           &
                               cledir(ntyp,mon,i1,2) * xf  +           &
                               cledir(ntyp,mon,i1,3) * xf2
               temp_bare(i1+9) = cledfu(ntyp,mon,i1)
            enddo
         endif

         if (fmelt(i).eq.1.0) then
            ml = 1
         else
            ml = 2
         endif

         do i1 = 1,9
            temp(i1) = cedir(ntyp,mon,i1,1)       +                    &
                       cedir(ntyp,mon,i1,2) * xf  +                    &
                       cedir(ntyp,mon,i1,3) * xf2
            temp(i1+9) = cedfu(ntyp,mon,i1)
         enddo

! For desert
         if (ntyp.eq.11) then
            sc2 = scov2(i) * scov2(i)
            sc1 = scov2(i)
            do i1 = 1,6
               temp(i1) = cedir2(ml,ntyp,mon,i1,nk,1)       +          &
                          cedir2(ml,ntyp,mon,i1,nk,2) * sc1 +          &
                          cedir2(ml,ntyp,mon,i1,nk,3) * sc2 + temp(i1)
               temp(i1+9) = cedfu2(ml,ntyp,mon,i1,1)       +           &
                            cedfu2(ml,ntyp,mon,i1,2) * sc1 +           &
                            cedfu2(ml,ntyp,mon,i1,3) * sc2 + temp(i1+9)
            enddo
            goto 500
         endif

         k2 = 1
         k1 = 2
         do ik = nk,1,-1
            if (xf.lt.xmi1(mon,ik)) then
               k1 = ik + 1
               k2 = ik
               goto 716
            endif
         enddo

 716     xm2 = xmi1(mon,k2)
         if (k1.le.nk) xm1 = xmi1(mon,k1)

! Modify albedo if snow on canopy
         if (scov(i).gt.0.025) then
            sc2 = scov(i) * scov(i)
            sc1 = scov(i)
            do i1 = 1,6
               if ((k2.ge.nk).or.(k2.le.1)) then
                  temp(i1) = cedir1(ml,ntyp,mon,i1,k2,1)       +       &
                             cedir1(ml,ntyp,mon,i1,k2,2) * sc1 +       &
                             cedir1(ml,ntyp,mon,i1,k2,3) * sc2 +temp(i1)
               else
                  xtm1 = cedir1(ml,ntyp,mon,i1,k1,1)       +           &
                         cedir1(ml,ntyp,mon,i1,k1,2) * sc1 +           &
                         cedir1(ml,ntyp,mon,i1,k1,3) * sc2
                  xtm2 = cedir1(ml,ntyp,mon,i1,k2,1)       +           &
                         cedir1(ml,ntyp,mon,i1,k2,2) * sc1 +           &
                         cedir1(ml,ntyp,mon,i1,k2,3) * sc2
                  temp(i1) = (xtm1 * (xm2 - xf) + xtm2 * (xf - xm1)) / &
                             (xm2 - xm1) + temp(i1)
               endif
               temp(i1+9) = cedfu1(ml,ntyp,mon,i1,1)       +           &
                            cedfu1(ml,ntyp,mon,i1,2) * sc1 +           &
                            cedfu1(ml,ntyp,mon,i1,3) * sc2 + temp(i1+9)
            enddo
         endif

! Modify albedo if snow on ground
         if (scov2(i).gt.0.025) then
            sc2 = scov2(i) * scov2(i)
            sc1 = scov2(i)
            do i1 = 1,6
               if ((k2.ge.nk).or.(k2.le.1)) then
                  temp(i1) = cedir2(ml,ntyp,mon,i1,k2,1)       +       &
                             cedir2(ml,ntyp,mon,i1,k2,2) * sc1 +       &
                             cedir2(ml,ntyp,mon,i1,k2,3) * sc2 +temp(i1)
               else
                  xtm1 = cedir2(ml,ntyp,mon,i1,k1,1)       +           &
                         cedir2(ml,ntyp,mon,i1,k1,2) * sc1 +           &
                         cedir2(ml,ntyp,mon,i1,k1,3) * sc2
                  xtm2 = cedir2(ml,ntyp,mon,i1,k2,1)       +           &
                         cedir2(ml,ntyp,mon,i1,k2,2) * sc1 +           &
                         cedir2(ml,ntyp,mon,i1,k2,3) * sc2
                  temp(i1) = (xtm1 * (xm2 - xf) + xtm2 * (xf - xm1)) / &
                             (xm2 - xm1) + temp(i1)
               endif
               temp(i1+9) = cedfu2(ml,ntyp,mon,i1,1)       +           &
                            cedfu2(ml,ntyp,mon,i1,2) * sc1 +           &
                            cedfu2(ml,ntyp,mon,i1,3) * sc2 + temp(i1+9)
            enddo
         endif

 500     continue

! Calculate output parameters
         radfac(i,1,1,2) = temp(10) / (temp(10) + temp(12))
         radfac(i,1,2,2) = temp(11) / (temp(11) + temp(13))
         radfac(i,2,1,2) = temp(12) / (temp(10) + temp(12))
         radfac(i,2,2,2) = temp(13) / (temp(11) + temp(13))
         salb(i,1,2) = temp(14) / temp_bare(14)
         salb(i,2,2) = temp(15) / temp_bare(15)
         extk(i,1,1,2) = temp(17)
         extk(i,2,1,2) = temp(18)
         radfac(i,1,1,1) = temp(1) / (temp(1) + temp(3))
         radfac(i,1,2,1) = temp(2) / (temp(2) + temp(4))
         radfac(i,2,1,1) = temp(3) / (temp(1) + temp(3))
         radfac(i,2,2,1) = temp(4) / (temp(2) + temp(4))
         salb(i,1,1) = temp(5) / temp_bare(5)
         salb(i,2,1) = temp(6) / temp_bare(6)
         extk(i,1,1,1) = temp(8) / xf
         extk(i,2,1,1) = temp(9) / xf
         extk(i,1,3,1) = cether(ntyp,mon,1)
         extk(i,1,3,2) = cether(ntyp,mon,2)
         extk(i,2,3,1) = cether(ntyp,mon,1)
         extk(i,2,3,2) = cether(ntyp,mon,2)

         salb(i,1,1) = min(salb(i,1,1)*visalb(i), 0.90)
         salb(i,1,2) = min(salb(i,1,2)*visalb(i), 0.90)
         salb(i,2,1) = min(salb(i,2,1)*niralb(i), 0.90)
         salb(i,2,2) = min(salb(i,2,2)*niralb(i), 0.90)

         zkat(i) = extk(i,1,3,2) * zlt(i,1) / vcover(i,1)
         zkat(i) = min(50.0, zkat(i))
         zkat(i) = max(10.0e-5, zkat(i))
         thermk(i) = exp(-zkat(i))
      enddo

      if (prin) then
         print *,'salb11: ',salb(1,1,1)
         print *,'salb12: ',salb(1,1,2)
         print *,'salb21: ',salb(1,2,1)
         print *,'salb22: ',salb(1,2,2)
         print *,'visalb: ',visalb(1)
         print *,'niralb: ',niralb(1)
         print *,'temp14: ',temp(14)
         print *,'temp15: ',temp(15)
         print *,'temp5: ',temp(6)
         print *,'temp6: ',temp(6)
      endif

      return
      end

!***********************************************************************
      subroutine intcap(istrip,zlt,z1,z2,snoww,satcap)

      implicit none
!***********************************************************************
!***  compute saturation capacity on canopy/ground
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! zlt(cg)      (-)       leaf area index
! z1           (m)       canopy base height
! z2           (m)       canopy top height
! snoww(cg)    (m)       liquid equiv. snow on canopy/ground
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! satcap(cg)   (m)       interception capacity
!-----------------------------------------------------------------------
      integer istrip, i
      real zlt(istrip,2), z1(istrip), z2(istrip), snoww(istrip,2)
      real satcap(istrip,2)
      real depcov(istrip)

      do i = 1,istrip
         satcap(i,1) = zlt(i,1) * 1.0e-4
         satcap(i,2) = 1.0e-4 * 1.0e-4
         depcov(i) = max(0.0, ((snoww(i,2) * 5.0) - z1(i)))
         depcov(i) = min(depcov(i), (z2(i) - z1(i))*0.95)
         satcap(i,1) = satcap(i,1) * (1.0 - depcov(i) / (z2(i) - z1(i)))
      enddo

      return
      end

!***********************************************************************
      subroutine swtop(istrip,nymd,nhms,dtc3x,lonco,latco,cosz)

      implicit none
!***********************************************************************
!***  compute cosine solar zenith angle averaged over timestep
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! nymd,nhms    (-)       date/time
! lonco        (deg)     longitude (-180 west to 180 east)
! latco        (deg)     latitude (-90 south to 90 north)
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! cosz         (radians) cosine of solar zenith angle
!-----------------------------------------------------------------------
      integer istrip, nymd, nhms, numstep
      integer n, nnn, nhms1, nhms2, imadd
      real lonco(istrip), latco(istrip)
      real cosz(istrip), dtc3x
      real cosz1(istrip), cosz2(istrip), icnt(istrip), sum1(istrip)

      do n = 1,istrip
         icnt(n) = 0
         sum1(n) = 0.
      enddo

      numstep = nint(dtc3x / 120.0)
      do nnn = 1,numstep
         imadd = nnn - 1
         nhms1 = nhms + (200 * imadd)
         nhms2 = nhms + (200 * imadd) + 100
         call astro_hyssib(istrip,nymd,nhms1,lonco,latco,cosz1)
         call astro_hyssib(istrip,nymd,nhms2,lonco,latco,cosz2)

         do n = 1,istrip
            sum1(n) = sum1(n) + max(cosz1(n),1e-10) + max(cosz2(n),1e-10)
            icnt(n) = icnt(n) + 2
         enddo
      enddo

      do n = 1,istrip
         if (icnt(n).gt.0) then
            cosz(n) = sum1(n) / icnt(n)
         else
            cosz(n) = -999.
         endif
      enddo

      return
      end

!***********************************************************************
      subroutine astro_hyssib(istrip,nymd,nhms,lonco,latco,cosz)
!***********************************************************************
!***  compute cosine solar zenith angle at individual unique time
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! nymd,nhms    (-)       date/time
! lonco        (deg)     longitude (-180 west to 180 east)
! latco        (deg)     latitude (-90 south to 90 north)
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! cosz         (radians) cosine of solar zenith angle
!-----------------------------------------------------------------------
      integer istrip, nymd, nhms, km, k, kp, iday, idayp1
      integer year, month, day, sec, n, dayscy
      real pi, zero, one, two, six, dg2rd, eqnx, ob
      real daylen, fac, thm, thp, thnow, zs, zc, sj, cj
      real cosz(istrip), lonco(istrip), latco(istrip)
      parameter(pi = 3.1415926535898)
      parameter(eqnx = 80.9028, dg2rd = pi/180., daylen = 86400.)
      parameter(dayscy = 365*4+1, ob = 23.45*dg2rd)
      parameter(zero = 0., one = 1., two = 2., six = 6.)
      real th(dayscy), t0, t1, t2, t3, t4, hc, mndy(12,4)
!      data mndy /0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305,     &
!                 335, 366, 397, 34*0/
      data mndy /0,  31,  60,  91, 121, 152, 182, 213, 244, 274, 305, 335, &
          366, 397, 425, 456, 486, 517, 547, 578, 609, 639, 670, 700, &
          731, 762, 790, 821, 851, 882, 912, 943, 974,1004,1035,1065, &
          1096,1127,1155,1186,1216,1247,1277,1308,1339,1369,1400,1430/

      logical first
      data first /.true./
      real funastro
      integer nsecfcalc

! Set current time
      year = nymd / 10000
      month = mod(nymd, 10000) / 100
      day = mod(nymd, 100)
      sec = nsecfcalc(nhms)

! Compute day-angles for 4-year cycle
      if (first) then
!         print *,'doing HY-SSiB astro for the first time'
         !do i = 15,48
         !   mndy(i,1) = mndy(i-12,1) + 365
         !enddo

         km = int(eqnx) + 1
         fac = km - eqnx
         t0 = zero
         t1 = funastro(t0) * fac
         t2 = funastro(zero + (t1 / two)) * fac
         t3 = funastro(zero + (t2 / two)) * fac
         t4 = funastro(zero + t3) * fac
         th(km) = (t1 + two * (t2 + t3) + t4) / six

         do k = 2,dayscy
            t1 = funastro(th(km))
            t2 = funastro(th(km) + (t1 / two))
            t3 = funastro(th(km) + (t2 / two))
            t4 = funastro(th(km) + t3)
            kp = mod(km,dayscy) + 1
            th(kp) = th(km) + (t1 + two * (t2 + t3) + t4) / six
            km = kp
         enddo

         first = .false.
      endif

! Compute earth-sun distance to current second
      iday = day + mndy(month,mod(year, 4)+1)
      idayp1 = mod(iday, dayscy) + 1
      thm = mod(th(iday), two*pi)
      thp = mod(th(idayp1), two*pi)

      if (thp.lt.thm) thp = thp + two * pi
      fac = float(sec) / daylen
      thnow = thm * (one - fac) + thp * fac

      zs = sin(thnow) * sin(ob)
      zc = sqrt(one - zs * zs)

! Compute cosine of the zenith angle
      fac = fac * two * pi + pi

      do n = 1,istrip
         sj = sin(latco(n) * dg2rd)
         cj = sqrt(one - sj * sj)
         hc = cos(fac + lonco(n) * dg2rd)
         cosz(n) = sj * zs + cj * zc * hc
         if (cosz(n).lt.zero) cosz(n) = zero
      enddo

      return
      end

!***********************************************************************
      real function funastro(y)
!***********************************************************************
      real pi, one, two, dg2rd, yrlen, ecc, per
      real tempdmm, y

      parameter(pi = 3.1415926535898)
      parameter(dg2rd = pi/180.)
      parameter(yrlen = 365.25)
      parameter(ecc = .0167, per = 102.0*dg2rd)
      parameter(one = 1., two = 2.)

      tempdmm = (two * pi/ ((one - ecc ** 2) ** 1.5)) * (one / yrlen)  &
                * (one - ecc * cos(y - per)) ** 2

      funastro = tempdmm

      return
      end

!***********************************************************************
      integer function nsecfcalc(n)
!***********************************************************************
      integer tempdmm, n

      tempdmm = (n / 10000 * 3600) + (mod(n,10000) / 100 * 60) +       &
                (mod(n,100))

      nsecfcalc = tempdmm

      return
      end

!***********************************************************************
      subroutine radspec(sfcnsw,sunang,xvisb,xvisd,xnirb,xnird)

      implicit none
!***********************************************************************
!***  split swdown radiation into 4 parts based on expected solar
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! sfcnsw       (W m-2)   net short wave at ground, downward > 0
! sunang       (radians) cosine of zenith angle
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! xvisb        (W m-2)   visible downward short wave beam at ground
! xvisd        (W m-2)   visible downward short wave diffuse at ground
! xnirb        (W m-2)   near ir downward short wave beam at ground
! xnird        (W m-2)   near ir downward short wave diffuse at ground
!-----------------------------------------------------------------------
      real sfcnsw, sunang, xvisb, xvisd, xnirb, xnird
      real fvisb, fvisd, fnirb, fnird, cloud, difrat, swdown, vnrat

      fvisb = 0.0
      fvisd = 0.0
      fnirb = 0.0      
      fnird = 0.0
 
      swdown = max(0.0, sfcnsw)
      if (swdown.gt.0.) then
         sunang = max(sunang, 0.01764)
         cloud = (1160. * sunang - swdown) / (963. * sunang)
         cloud = max(cloud, 0.)
         cloud = min(cloud, 1.)

         difrat = 0.0604 / (sunang - 0.0223) + 0.0683
         difrat = max(difrat, 0.)
         difrat = min(difrat, 1.)

         difrat = difrat + (1. - difrat) * cloud
         vnrat = (580. - cloud * 464.) / ((580. - cloud * 499.) +      &
                 (580. - cloud * 464.))

         fvisb = (1. - difrat) * vnrat
         fvisd = difrat * vnrat
         fnirb = (1. - difrat) * (1. - vnrat)
         fnird = difrat * (1. - vnrat)
      endif

      xvisb = swdown * fvisb
      xvisd = swdown * fvisd
      xnirb = swdown * fnirb
      xnird = swdown * fnird

      return
      end

!***********************************************************************
      subroutine hyssib_routines(istrip,grav,refzh,pie,dtc3x,cpair,gasr,&
           hlat,snomel,clai,cw,cs,dkfac,epsfac,stefan,thres,sskin,     &
           tf,tf_snow,tf_drain,delsig,lonco,latco,ityp,z2,dd_in,z0_in, &
           phsat,poros,bee,satco,slope,zdepth,vcover,sgfg,sdens,       &
           thermk,extk,cosz,zlt,green,chil,rbc_in,rdc_in,topt,tll,tu,  &
           defac,ph1,ph2,rootd,rstpar,satcap,wiltp,psilow,topostd,     &
           tdd,bps,psb,ros,thm,zb,tm,ps,sh,spdm,tg,capac1,sn1,         &
           capac2,sn2,gw1,gw2,gw3,tc,td,tsn,ppl,ppc,radn,radc3,        &
           radt,etmass,hflux,roff,swnet,lwnet,qle,qh,qg,qf,qv,qtau,    &
           qa,delsurfheat,delcoldcont,snowf,rainf,evap,qs,qrec,qsb,    &
           qsm,qfz,qst,delsoilmoist,delswe,delintercept,snowt,         &
           vegtc,baresoilt,avgsurft,radteff,swe,sweveg,soilmoist1,     &
           soilmoist2,soilmoist3,soiltemp,soilwet,potevap,ecanop,      &
           tveg,esoil,rootmoist,canopint,subsnow,subsurf,acond,        &
           ccond,snowdepth,sliqfrac,ecanopt,ecanopi,eground,hfluxc,    &
           hfluxg,watcan,watgrd,snocan,snogrd,wet1,wet2,wet3,prin)

      implicit none
!***********************************************************************
!***  This code is the offline HY-SSiB Model (Hydrology with Simple SiB)
!
! Additional notes:
!
! istrip is the number of points in the 1-D vector of calculated points.
!
! Vegetation type index used: Xue et al, 1991, Table 2, Type 13 = Ice
!
! "Blended" in the documentation means that this variable is equal to
! either the diurnal snow layer temperature or the upper ground layer
! temperature, depending on if the snow model is "active" or not.
!
! fpar is the fraction of photosynthetically active radiation.
!
! 4 variables (dd, z0, rdc, rbc) use "_in" to prevent the temporary
! values calculated in routine airmod from overwriting the values.
!
! Forcing variables tm and sh are modified during the timestep.
!
! Units of (X/dt) indicate "X per timestep".
!-----------------------------------------------------------------------
!     *** indices ***
!-----------------------------------------------------------------------
! cg =1...canopy
! cg =2...ground cover
! vn =1...visible      (0.0-0.7 micron)
! vn =2...near-infrared(0.7-3.0 micron)
! bd =1...beam
! bd =2...diffuse
! ld =1...live leaves
! ld =2...dead leaves
! vnt=1...visible      (0.0-0.7 micron)
! vnt=2...near-infrared(0.7-3.0 micron)
! vnt=3...thermal
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------
! grav         (m sec-2) acceleration of gravity
! pie          (-)       3.14159....
! dtc3x        (sec)     timestep
! cpair     (J K-1 kg-1) heat capacity of dry air at constant pressure
! gasr      (J K-1 kg-1) gas constant of dry air
! hlat         (J kg-1)  latent heat of vaporization
! snomel       (J m-3)   latent heat of fusion X water density
! clai      (J m-2 K-1)  heat capacity of foliage
! cw        (J m-3 K-1)  heat capacity of liquid water X water density
! cs        (J m-3 K-1)  heat capacity of frozen water X water density
! dkfac        (m-1)     solar extinction coefficient through snow
! epsfac       (-)       ratio of dry air to water vapor gas constant
! stefan    (W m-2 K-4)  Stefan-Boltzmann constant
! thres        (m)       grnd snow depth threshold activating snow model
! sskin        (m)       depth of diurnal snow layer
! tf           (K)       freezing temperature
! tf_snow      (K)       reference-level temp. below which precip = snow
! tf_drain     (K)       temperature below which soil drainage stops
! delsig       (-)       thickness (in sigma) of lowest GEOS GCM layer
!-----------------------------------------------------------------------
! Description of point
!-----------------------------------------------------------------------
! lonco        (deg)     longitude (-180 west to 180 east)
! latco        (deg)     latitude (-90 south to 90 north)
! ityp         (1-13)    vegetation index
!-----------------------------------------------------------------------
! Vegetation parameters
!-----------------------------------------------------------------------
! z2           (m)       canopy top height
! dd           (m)       displacement height
! z0           (m)       roughness height
! vcover(cg)   (0-1)     fraction of vegetation cover [fpar/greenness]
! zlt(cg)      (-)       leaf area index
! green        (0-1)     green leaf fraction
! chil         (-)       leaf angle distribution factor
! rbc   (sqr((s m-1)))   bulk canopy boundary layer resistance coeff.
! rdc          (-)       ground to canopy air space resistance coeff.
! topt         (K)       optimum temperature for rst calculation
! tll          (K)       low temperature for rst calculation
! tu           (K)       top temperature for rst calculation
! defac        (-)       dew factor for rst calculation
! ph1          (-)       stome slope factor
! ph2          (-)       point at which stomates close
! rootd(cg)    (m)       rooting depth for canopy [1] or ground [2]
! rstpar(3)    (-)       par influence on stomatal resist. coefficients
! satcap(cg)   (m)       interception capacity
!-----------------------------------------------------------------------
! Soil parameters
!-----------------------------------------------------------------------
! phsat        (m)       soil moisture potential at saturation
! poros        (0-1)     soil porosity
! bee          (-)       Clapp-Hornberger empirical constant
! satco        (m/sec)   saturation hydraulic conductivity
! slope        (0-100)   average topographic slope in percent
! zdepth(3)    (m)       exact independent depth of 3 soil layers
! wiltp        (0-1)     wilting point wetness of soil, SMP3, Eq. 9
! psilow       (m)       lowest matric potential from evap., SMP3,Eq. 15
! topostd      (m)       topography standard deviat. 1x1deg from GTOPO30
! tdd          (K)       climatological average deep soil temperature
!-----------------------------------------------------------------------
! Radiation parameters
!-----------------------------------------------------------------------
! thermk       (0-1)     canopy emissivity
! extk(cg,vnt,bd)  (-)   extinction coefficient
! cosz         (radians) cosine of solar zenith angle
!-----------------------------------------------------------------------
! Forcing variables
!-----------------------------------------------------------------------
! tm           (K)       forcing reference-level temperature
! ps           (mb)      forcing surface pressure
! sh           (kg kg-1) forcing reference-level watervapor mixing ratio
! spdm         (m/sec)   forcing reference-level wind speed
! ppl          (mm/dt)   forcing large-scale precipitation
! ppc          (mm/dt)   forcing convective precipitation
! radn(vnt,bd) (W m-2)   downward sw/lw radiation at the surface
!-----------------------------------------------------------------------
! Other input variables
!-----------------------------------------------------------------------
! bps          (-)       (ps/p0)**(gasr/cpair); p0 = 1000.0 mb
! psb          (mb)      depth of pbl
! ros          (kg m-3)  surface air density
! thm          (K)       reference-level potential temperature
! zb           (m)       elevation of the pbl top
! radc3(cg)    (W m-2)   absorbed sw/lw radiation at the surface
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! Diagnostic variables
!-----------------------------------------------------------------------
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
! sdens        (kg m-3)  density of bulk snow layer
! radt(cg)     (W m-2)   net radiation
! etmass       (mm/sec)  evaporation to the atmosphere (mass)
! hflux        (W m-2)   sensible heat flux to the atmosphere
! roff         (m/dt)    total runoff
! swnet        (W m-2)   net shortwave at ground (incoming - outgoing)
! lwnet        (W m-2)   net longwave at ground (incoming - outgoing)
! qle          (W m-2)   evaporation to the atmosphere (energy)
! qh           (W m-2)   sensible heat flux to the atmosphere
! qg           (W m-2)   heat flux into the ground
! qf           (W m-2)   energy of fusion (melt or freeze)
! qv           (W m-2)   energy of sublimation (sublim. or deposition)
! qtau         (N m-2)   momentum flux lost to atmosphere
! qa           (W m-2)   momentum flux lost to atmosphere
! delsurfheat  (J m-2)   change in surface heat storage
! delcoldcont  (J m-2)   change in snow cold content
! snowf        (J m-2)   snowfall rate
! rainf        (J m-2)   rainfall rate
! qs      (kg m-2 sec-1) surface runoff
! qrec    (kg m-2 sec-1) recharge (from water table)
! qsb     (kg m-2 sec-1) sub-surface runoff
! qsm     (kg m-2 sec-1) snowmelt
! qfz     (kg m-2 sec-1) refreeze
! qst     (kg m-2 sec-1) liquid water flowing out of snowpack
! delsoilmoist (kg m-2)  change in soil moisture
! delswe       (kg m-2)  change in snow water equivalent
! delintercept (kg m-2)  change in interception storage
! snowt        (K)       snow on ground temperature
! vegtc        (K)       canopy temperature
! baresoilt    (K)       bare soil temperature
! avgsurft     (K)       average surface temperature
! radteff      (K)       effective radiative surface temperature
! swe          (kg m-2)  snow water equivalent
! sweveg       (kg m-2)  SWE intercepted by vegetation
! soilmoist1   (kg m-2)  Soil water content of layer 1
! soilmoist2   (kg m-2)  Soil water content of layer 2
! soilmoist3   (kg m-2)  Soil water content of layer 3
! soiltemp     (K)       Soil temperature
! soilwet      (0-1)     Total soil wetness
! potevap (kg m-2 sec-1) potential evapotranspiration
! ecanop  (kg m-2 sec-1) interception loss evaporation
! tveg    (kg m-2 sec-1) vegetation transpiration
! esoil   (kg m-2 sec-1) bare soil evaporation
! rootmoist    (kg m-2)  root zone soil moisture
! canopint     (kg m-2)  canopy interception
! subsnow (kg m-2 sec-1) snow sublimation from snowpack
! subsurf (kg m-2 sec-1) snow sublimation from vegetation
! acond        (m/sec)   aerodynamic conductance
! ccond        (m/sec)   stomatal conductance
! snowdepth    (m)       depth of snow layer
! sliqfrac     (0-1)     snow liquid fraction on ground
!-----------------------------------------------------------------------
! Prognostic variables
!-----------------------------------------------------------------------
! tg           (K)       ground temperature
! tc           (K)       canopy temperature
! tsn          (K)       snow on ground temperature
! td           (K)       deep soil temperature
! capac(cg)    (m)       liquid equiv. water&snow on canopy/ground
! snow(cg)     (m)       liquid equiv. snow on canopy/ground
! w(3)         (0-1)     wetness of surface/root/recharge zone
!-----------------------------------------------------------------------
      logical prin
      integer istrip, i, iveg, il, ipoint
      real grav, refzh, pie, dtc3x, cpair, gasr, hlat, snomel, clai, cw, cs
      real dkfac, epsfac, stefan, thres, sskin, tf, tf_snow, tf_drain
      real delsig, baseflow, timcon, hlat3i, summ, theta, difsl, water
      integer ityp(istrip)
      real lonco(istrip), latco(istrip)
      real z2(istrip), dd_in(istrip), z0_in(istrip)
      real vcover(istrip,2), zlt(istrip,2), green(istrip)
      real chil(istrip), rbc_in(istrip), rdc_in(istrip)
      real dd(istrip), z0(istrip), rbc(istrip), rdc(istrip)
      real topt(istrip), tll(istrip), tu(istrip)
      real defac(istrip), ph1(istrip), ph2(istrip)
      real rootd(istrip,2), rstpar(istrip,3), satcap(istrip,2)
      real phsat(istrip), poros(istrip), bee(istrip), satco(istrip)
      real slope(istrip), zdepth(istrip,3), wiltp(istrip)
      real psilow(istrip), topostd(istrip), tdd(istrip)
      real thermk(istrip), extk(istrip,2,3,2), cosz(istrip)
      real tm(istrip), ps(istrip), sh(istrip), spdm(istrip)
      real ppl(istrip), ppc(istrip), radn(istrip,3,2)
      real bps(istrip), psb(istrip), ros(istrip)
      real thm(istrip), zb(istrip), radc3(istrip,2), radt(istrip,2)
      real sgfg(istrip), sdens(istrip)
      real etmass(istrip), hflux(istrip), roff(istrip)
      real swnet(istrip), lwnet(istrip), qle(istrip), qh(istrip)
      real qg(istrip), qf(istrip), qv(istrip), qtau(istrip), qa(istrip)
      real delsurfheat(istrip), delcoldcont(istrip), snowf(istrip)
      real rainf(istrip), evap(istrip), qs(istrip), qrec(istrip)
      real qsb(istrip), qsm(istrip), qfz(istrip), qst(istrip)
      real delsoilmoist(istrip), delswe(istrip), delintercept(istrip)
      real snowt(istrip), vegtc(istrip), baresoilt(istrip)
      real avgsurft(istrip), radteff(istrip), swe(istrip)
      real sweveg(istrip), soilmoist1(istrip), soilmoist2(istrip)
      real soilmoist3(istrip), soiltemp(istrip), soilwet(istrip)
      real potevap(istrip), ecanop(istrip), tveg(istrip), esoil(istrip)
      real rootmoist(istrip), canopint(istrip), subsnow(istrip)
      real subsurf(istrip), acond(istrip), ccond(istrip)
      real snowdepth(istrip), sliqfrac(istrip)
      real ecanopt(istrip), ecanopi(istrip), eground(istrip)
      real hfluxc(istrip), hfluxg(istrip), watcan(istrip)
      real watgrd(istrip), snocan(istrip), snogrd(istrip)
      real wet1(istrip), wet2(istrip), wet3(istrip)

      real tg(istrip), tc(istrip), tsn(istrip), td(istrip)
      real capac1(istrip), sn1(istrip), capac2(istrip), sn2(istrip)
      real gw1(istrip), gw2(istrip), gw3(istrip)
      real capac(istrip,2), snow(istrip,2), w(istrip,3)
      real tgsb(istrip), swg(istrip), sws(istrip), fluxd(istrip)
      real fluxs(istrip), fluxmax(istrip), fluxmelt(istrip)
      real qqq(istrip,3), wattabl(istrip), q3g(istrip), roffg(istrip)
      real startcapac(istrip), startsnow(istrip), startswe(istrip)
      real startsm(istrip), starttotal(istrip), startint(istrip)
      real currtotal(istrip), currsnow(istrip), refreez(istrip)
      real cg(istrip), cc(istrip), cgsbi(istrip), zmelt(istrip)
      real zfrez(istrip), ymelt(istrip,2), yfrez(istrip,2)
      real xmelt(istrip,2), xfrez(istrip,2), smelt(istrip)
      real ewt(istrip,3), wait(istrip,3), sevap(istrip,3), ef(istrip,3)
      real anetlw(istrip,2), anetsw(istrip,2), fluxef(istrip)
      real dtc(istrip), dtg(istrip), dtd(istrip)
      real dtgsb(istrip), dtsn(istrip), dgdt(istrip)
      real chisl(istrip), zlwup(istrip), tgeff(istrip), epot(istrip)
      real sfall(istrip), absoil(istrip), egs(istrip), txsc(istrip)
      real ect(istrip), egt(istrip), totdep(istrip), dep(istrip)
      real cu(istrip), ustar(istrip), drag(istrip), aaa(istrip)
      real cbal(istrip), gbal(istrip), sbal(istrip), bbal(istrip)
      real eci(istrip), hc(istrip), chf(istrip), ghf(istrip)
      real hg(istrip), egi(istrip), shf(istrip), bhf(istrip)
      real phsoil(istrip,3), eintmass(istrip), etranmass(istrip)
      real ebaremass(istrip), etrancmass(istrip), eintcmass(istrip)
      real eintgmass(istrip), ra(istrip)
      real rst(istrip,2), par(istrip), pd(istrip)
      real gmt(istrip,3), gmq(istrip,3), gmu(istrip,4)
      real precheat(istrip), srtes_qst(istrip)

      ipoint = 1

! Modify magnitudes for SSiB
      do i = 1,istrip
         w(i,1) = gw1(i)
         w(i,2) = gw2(i)
         w(i,3) = gw3(i)
         snow(i,1) = sn1(i) * 0.001
         snow(i,2) = sn2(i) * 0.001
         capac(i,1) = capac1(i) * 0.001 + snow(i,1)
         capac(i,2) = capac2(i) * 0.001 + snow(i,2)
         gmt(i,1) = 0.0
         gmt(i,2) = 1.56
         gmt(i,3) = 0.0
         gmq(i,1) = 0.0
         gmq(i,2) = 1.56
         gmq(i,3) = 0.0
         gmu(i,1) = 0.0
         gmu(i,2) = 1.56
         gmu(i,3) = 0.0
         gmu(i,4) = 0.0
         dd(i) = dd_in(i)
         z0(i) = z0_in(i)
         rbc(i) = rbc_in(i)
         rdc(i) = rdc_in(i)
      enddo

! Set constants
      baseflow = 1.0 / 365.5 * (dtc3x / 86400.)

      do i = 1,istrip
! Initialize evaporation, sensible heat and runoff at timestep start
         etmass(i) = 0.0
         hflux(i) = 0.0
         roff(i) = 0.0
! Calculates soil water budget prior to calling pbl
         startsm(i) = ((w(i,1) * zdepth(i,1)) + (w(i,2) * zdepth(i,2)) &
                     + (w(i,3) * zdepth(i,3))) * poros(i)
         startcapac(i) = capac(i,1) + capac(i,2)
         startsnow(i) = snow(i,1) + snow(i,2)
         startint(i) = capac(i,1)
         startswe(i) = capac(i,2)
         starttotal(i) = startsm(i) + startcapac(i)
         if (prin.and.(i.eq.ipoint)) then
            print *,' '
            print *,'Initial conditions: '
            print *,'starttotal: ',starttotal(i)
            print *,'startsnow: ',startsnow(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'w_1: ',w(i,1)
            print *,'w_2: ',w(i,2)
            print *,'w_3: ',w(i,3)
            print *,'waterinlayer_1: ',(w(i,1) * poros(i) * zdepth(i,1))
            print *,'waterinlayer_2: ',(w(i,2) * poros(i) * zdepth(i,2))
            print *,'waterinlayer_3: ',(w(i,3) * poros(i) * zdepth(i,3))
         endif
! Update snow density with time based on Douville et al, 1995, Eq. 13
! SMP3, Eq. 1 (0.24 in paper should be 0.14, or 1/7)
         if ((sgfg(i).gt.0.5).and.(sdens(i).lt.300.0)) then
            sdens(i) = ((((sdens(i) - 300.0) * exp(-1./7. * dtc3x      &
                       / 86400.0)) + 300.0) * sgfg(i)) + (100.0 *      &
                       (1.0 - sgfg(i)))
         endif
      enddo

      do i = 1,istrip
! Calculating conductivity of soil, based on soil wetness and porosity
         theta = w(i,1) * poros(i)
         chisl(i) = (9.8e-4 + 1.2e-3 * theta) / (1.1 - 0.4 * theta)    &
                    * 418.6
         difsl = 5.0e-7
         cg(i) = chisl(i) * sqrt(86400.0 / (pie * difsl)) * 0.5
         if (snow(i,2).ge.thres) then
! If there was no snow the previous timestep, blending the deep and
! ground temperatures to one ground temperature, and assigning to the
! ground temperature the bulk snow temperature
            if (sgfg(i).lt.0.5) then
! Deep temperature is now a blend of ground and deep temps
               td(i) = ((cg(i) * tg(i)) + (2.0 * sqrt(pie * 365.0) *   &
                       cg(i) * td(i))) / (cg(i) + (2.0 * sqrt(pie *    &
                       365.0) * cg(i)))
! Tg now represents bulk snow temperature, and takes on ground temp
               tg(i) = tsn(i)
               sdens(i) = 100.0
            endif
            cg(i) = cg(i) + (2.0 * sqrt(365.0 * pie) * cg(i))
            sgfg(i) = 1.0
         else
! If there was snow the previous timestep, calculating the new ground
! /snow blended temp to the temp of the old separate ground and snow
! temperatures, with their capacities
            if (sgfg(i).gt.0.5) then
               water = capac(i,2) - snow(i,2)
! Ground temperatue is a blend of deep ground and bulk snow temps
               tg(i) = ((td(i) * cg(i)) + (tg(i) * ((cs * snow(i,2)) + &
                       (cw * water)))) /  (cg(i) + ((cs * snow(i,2)) + &
                       (cw * water)))
! Tsn now is the same as the ground temperature
               tsn(i) = tg(i)
            endif
            sgfg(i) = 0.0
            sdens(i) = -999.
         endif
! Setting blended temperature based on snow flag
         tgsb(i) = (tsn(i) * sgfg(i)) + (tg(i) * (1.0 - sgfg(i)))
         zmelt(i) = 0.0
         zfrez(i) = 0.0
         sevap(i,1) = 0.0
         sevap(i,2) = 0.0
         sevap(i,3) = 0.0
         do iveg = 1,2
            anetsw(i,iveg) = radc3(i,iveg)
         enddo
      enddo

      if (prin) then
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'w_1: ',w(ipoint,1)
         print *,'w_2: ',w(ipoint,2)
         print *,'w_3: ',w(ipoint,3)
         print *,'phsat: ',phsat(ipoint)
         print *,'bee: ',bee(ipoint)
         print *,'CALLING ROOT'
      endif

      call root(istrip,w,phsat,bee,phsoil)

      if (prin) then
         print *,'phsoil_1: ',phsoil(ipoint,1)
         print *,'phsoil_2: ',phsoil(ipoint,2)
         print *,'phsoil_3: ',phsoil(ipoint,3)

         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'dkfac: ',dkfac
         print *,'stefan: ',stefan
         print *,'sskin: ',sskin
         print *,'vcover_1: ',vcover(ipoint,1)
         print *,'vcover_2: ',vcover(ipoint,2)
         print *,'thermk: ',thermk(ipoint)
         print *,'snow_2: ',snow(ipoint,2)
         print *,'tc: ',tc(ipoint)
         print *,'tgsb: ',tgsb(ipoint)
         print *,'sgfg: ',sgfg(ipoint)
         print *,'radn_11: ',radn(ipoint,1,1)
         print *,'radn_12: ',radn(ipoint,1,2)
         print *,'radn_21: ',radn(ipoint,2,1)
         print *,'radn_22: ',radn(ipoint,2,2)
         print *,'radc3_1: ',radc3(ipoint,1)
         print *,'radc3_2: ',radc3(ipoint,2)
         print *,'CALLING RADUSE'
      endif

      call raduse_hyssib(istrip,prin,ipoint,dkfac,stefan,sskin,        &
           vcover,thermk,snow,tc,tgsb,sgfg,radn,radc3,radt,par,pd,sws, &
           swg)

      if (prin) then
         print *,'radt_1: ',radt(ipoint,1)
         print *,'radt_2: ',radt(ipoint,2)
         print *,'par: ',par(ipoint)
         print *,'pd: ',pd(ipoint)
         print *,'anetsw_1: ',anetsw(ipoint,1)
         print *,'anetsw_2: ',anetsw(ipoint,2)
         print *,'sws: ',sws(ipoint)
         print *,'swg: ',swg(ipoint)

         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'ityp: ',ityp(ipoint)
         print *,'zlt_1: ',zlt(ipoint,1)
         print *,'zlt_2: ',zlt(ipoint,2)
         print *,'green: ',green(ipoint)
         print *,'vcover_1: ',vcover(ipoint,1)
         print *,'vcover_2: ',vcover(ipoint,2)
         print *,'chil: ',chil(ipoint)
         print *,'rstpar_1: ',rstpar(ipoint,1)
         print *,'rstpar_2: ',rstpar(ipoint,2)
         print *,'rstpar_3: ',rstpar(ipoint,3)
         print *,'cosz: ',cosz(ipoint)
         print *,'par: ',par(ipoint)
         print *,'pd: ',pd(ipoint)
         print *,'CALLING STOMAT'
      endif

      call stomat(istrip,pie,ityp,zlt,green,vcover,chil,rstpar,        &
           extk,cosz,par,pd,rst)

      if (prin) then
         print *,'rst_1: ',rst(ipoint,1)
         print *,'rst_2: ',rst(ipoint,2)

         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'dt: ',dtc3x
         print *,'snomel: ',snomel
         print *,'clai: ',clai
         print *,'cw: ',cw
         print *,'cs: ',cs
         print *,'tf: ',tf
         print *,'tf_snow: ',tf_snow
         print *,'tf_drain: ',tf_drain
         print *,'poros: ',poros(ipoint)
         print *,'satco: ',satco(ipoint)
         print *,'zdepth_1: ',zdepth(ipoint,1)
         print *,'zdepth_2: ',zdepth(ipoint,2)
         print *,'zdepth_3: ',zdepth(ipoint,3)
         print *,'topostd: ',topostd(ipoint)
         print *,'capac_1: ',capac(ipoint,1)
         print *,'capac_2: ',capac(ipoint,2)
         print *,'snow_1: ',snow(ipoint,1)
         print *,'snow_2: ',snow(ipoint,2)
         print *,'w_1: ',w(ipoint,1)
         print *,'w_2: ',w(ipoint,2)
         print *,'w_3: ',w(ipoint,3)
         print *,'satcap_1: ',satcap(ipoint,1)
         print *,'satcap_2: ',satcap(ipoint,2)
         print *,'tm: ',tm(ipoint)
         print *,'ppl: ',ppl(ipoint)
         print *,'ppc: ',ppc(ipoint)
         print *,'tg: ',tg(ipoint)
         print *,'tc: ',tc(ipoint)
         print *,'tsn: ',tsn(ipoint)
         print *,'tgsb: ',tgsb(ipoint)
         print *,'sdens: ',sdens(ipoint)
         print *,'cg: ',cg(ipoint)
         print *,'CALLING INTERC'
      endif

      call interc_hyssib(istrip,dtc3x,snomel,clai,cw,cs,sskin,         &
           tf,tf_snow,poros,satco,zdepth,zlt,vcover,topostd,tm,ppl,    &
           ppc,tg,w,capac,snow,satcap,tc,tsn,sgfg,sdens,tgsb,cc,cg,    &
           cgsbi,extk,sfall,ymelt,yfrez,roff,precheat)

      if (prin) then
         print *,'tg: ',tg(ipoint)
         print *,'tc: ',tc(ipoint)
         print *,'tsn: ',tsn(ipoint)
         print *,'tgsb: ',tgsb(ipoint)
         print *,'capac_1: ',capac(ipoint,1)
         print *,'capac_2: ',capac(ipoint,2)
         print *,'snow_1: ',snow(ipoint,1)
         print *,'snow_2: ',snow(ipoint,2)
         print *,'sfall: ',sfall(ipoint)
         print *,'cc: ',cc(ipoint)
         print *,'cgsb: ',(1.0/cgsbi(ipoint))
         print *,'ymelt_1: ',ymelt(ipoint,1)
         print *,'ymelt_2: ',ymelt(ipoint,2)
         print *,'yfrez_1: ',yfrez(ipoint,1)
         print *,'yfrez_2: ',yfrez(ipoint,2)
         print *,'roff: ',roff(ipoint)
      endif

      do i = 1,istrip
         currtotal(i) = (((w(i,1) * zdepth(i,1)) + (w(i,2) *           &
                        zdepth(i,2)) + (w(i,3) * zdepth(i,3))) *       &
                        poros(i)) + capac(i,1) + capac(i,2) + roff(i)  &
                        - ((ppl(i) + ppc(i)) / 1000.0)
         if ((abs(currtotal(i)-starttotal(i)).gt.0.00001).or.         &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'INTERC CAPAC DIFF: ',(currtotal(i)-starttotal(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currtotal: ',currtotal(i)
            print *,'starttotal: ',starttotal(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'w_1: ',w(i,1)
            print *,'w_2: ',w(i,2)
            print *,'w_3: ',w(i,3)
            print *,'waterinlayer_1: ',(w(i,1) * poros(i) * zdepth(i,1))
            print *,'waterinlayer_2: ',(w(i,2) * poros(i) * zdepth(i,2))
            print *,'waterinlayer_3: ',(w(i,3) * poros(i) * zdepth(i,3))
            print *,'ppl: ',(ppl(i)/1000.)
            print *,'ppc: ',(ppc(i)/1000.)
            print *,'roff: ',roff(i)
            print *,'tm: ',tm(i)
            print *,'tc: ',tc(i)
            print *,'tg: ',tg(i)
            print *,'tsn: ',tsn(i)
            print *,'tgsb: ',tgsb(i)
         endif
         if (abs(currtotal(i)-starttotal(i)).gt.0.00001) then
            print *,'interc capac error!'
!            stop
         endif
         currsnow(i) = snow(i,1) + snow(i,2) - sfall(i) + ymelt(i,1) + &
                       ymelt(i,2) + yfrez(i,1) + yfrez(i,2)
         if ((abs(currsnow(i)-startsnow(i)).gt.0.00001).or.           &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'INTERC SNOW DIFF: ',(currsnow(i)-startsnow(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currsnow: ',currsnow(i)
            print *,'startsnow: ',startsnow(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'sfall: ',sfall(i)
            print *,'ymelt_1: ',ymelt(i,1)
            print *,'ymelt_2: ',ymelt(i,2)
            print *,'yfrez_1: ',yfrez(i,1)
            print *,'yfrez_2: ',yfrez(i,2)
         endif
         if (abs(currsnow(i)-startsnow(i)).gt.0.00001) then
            print *,'interc snow error!'
!            stop
         endif
      enddo

      if (prin) then
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'CALLING SFLXES'
      endif

      call sflxes(istrip,prin,ipoint,grav,refzh,pie,dtc3x,cpair,gasr,  &
           hlat,snomel,cs,epsfac,stefan,sskin,tf,delsig,ityp,lonco,    &
           latco,z2,dd,z0,phsat,poros,phsoil,psilow,bee,zdepth,zlt,    &
           vcover,rbc,rdc,topt,tll,tu,defac,ph1,ph2,rootd,chisl,bps,   &
           psb,ros,thm,zb,topostd,sws,swg,tm,ps,sh,spdm,ppl,ppc,       &
           starttotal,startsnow,cu,ustar,tg,w,capac,satcap,snow,tc,td, &
           tsn,tgsb,cc,cgsbi,sgfg,sdens,fluxd,fluxs,fluxmax,fluxmelt,  &
           sfall,ymelt,yfrez,thermk,radt,rst,ra,ewt,dtc,dtg,dtsn,      &
           dtgsb,dgdt,etmass,hflux,roff,hc,hg,eci,egi,ect,egt,egs,     &
           epot,gmt,gmq,sevap,eintmass,etranmass,ebaremass,etrancmass, &
           eintcmass,eintgmass,bhf,chf,ghf,shf)

      do i = 1,istrip
         currtotal(i) = (((w(i,1) * zdepth(i,1)) + (w(i,2) *           &
                        zdepth(i,2)) + (w(i,3) * zdepth(i,3))) *       &
                        poros(i)) + capac(i,1) + capac(i,2) + roff(i)  &
                        + (eintmass(i) / 1000.0) - ((ppl(i) + ppc(i))  &
                        / 1000.0)
         if ((abs(currtotal(i)-starttotal(i)).gt.0.00001).or.         &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'SFLXES CAPAC DIFF: ',(currtotal(i)-starttotal(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currtotal: ',currtotal(i)
            print *,'starttotal: ',starttotal(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'w_1: ',w(i,1)
            print *,'w_2: ',w(i,2)
            print *,'w_3: ',w(i,3)
            print *,'waterinlayer_1: ',(w(i,1) * poros(i) * zdepth(i,1))
            print *,'waterinlayer_2: ',(w(i,2) * poros(i) * zdepth(i,2))
            print *,'waterinlayer_3: ',(w(i,3) * poros(i) * zdepth(i,3))
            print *,'ppl: ',(ppl(i)/1000.)
            print *,'ppc: ',(ppc(i)/1000.)
            print *,'roff: ',roff(i)
            print *,'eintmass: ',(eintmass(i)/1000.)
            print *,'tm: ',tm(i)
            print *,'tc: ',tc(i)
            print *,'tg: ',tg(i)
            print *,'tsn: ',tsn(i)
            print *,'tgsb: ',tgsb(i)
         endif
         if (abs(currtotal(i)-starttotal(i)).gt.0.00001) then
            print *,'sflxes capac error!'
!            stop
         endif
         currsnow(i) = snow(i,1) + snow(i,2) - sfall(i) + ymelt(i,1) + &
                       ymelt(i,2) + yfrez(i,1) + yfrez(i,2) +          &
                       sevap(i,1) + sevap(i,2) + sevap(i,3)
         if ((abs(currsnow(i)-startsnow(i)).gt.0.00001).or.           &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'SFLXES SNOW DIFF: ',(currsnow(i)-startsnow(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currsnow: ',currsnow(i)
            print *,'startsnow: ',startsnow(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'sfall: ',sfall(i)
            print *,'ymelt_1: ',ymelt(i,1)
            print *,'ymelt_2: ',ymelt(i,2)
            print *,'yfrez_1: ',yfrez(i,1)
            print *,'yfrez_2: ',yfrez(i,2)
            print *,'sevap_1: ',sevap(i,1)
            print *,'sevap_2: ',sevap(i,2)
            print *,'sevap_3: ',sevap(i,3)
         endif
         if (abs(currsnow(i)-startsnow(i)).gt.0.00001) then
            print *,'sflxes snow error!'
!            stop
         endif
      enddo

      if (prin) then
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'RETURNING FROM SFLXES'
         print *,'etmass: ',(etmass(ipoint)/dtc3x)
         print *,'hflux: ',(hflux(ipoint)/dtc3x)
      endif

! Continue to update sib variables
      do i = 1,istrip
         anetlw(i,1) = radt(i,1) - anetsw(i,1)
         anetlw(i,2) = radt(i,2) - anetsw(i,2) + sws(i) + swg(i)
         tc(i) = tc(i) + dtc(i)
         tg(i) = tg(i) + dtg(i)
         tgsb(i) = tgsb(i) + dtgsb(i)
         tsn(i) = tsn(i) + dtgsb(i)
         zlwup(i) = -(anetlw(i,1) + anetlw(i,2) - radn(i,3,2))
         tgeff(i) = sqrt(sqrt((zlwup(i)/stefan)))
      enddo

! Dumping of small capac values onto soil surface store
      do iveg = 1,2
         do i = 1,istrip
            if (capac(i,iveg).le.1.e-6) then
               w(i,1) = w(i,1) + capac(i,iveg) /(poros(i) * zdepth(i,1))
               capac(i,iveg) = 0.0
! Making sure to zero out snow when capac = 0.  9/5/97 dmm
               ymelt(i,iveg) = ymelt(i,iveg) + snow(i,iveg)
               snow(i,iveg) = 0.0
            endif
         enddo
      enddo

      if (prin) then
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'CALLING SNOWMT'
      endif

      call snowmt(istrip,dtc3x,snomel,cw,cs,sskin,tf,poros,            &
           zdepth,topostd,tg,w,capac,snow,tc,td,tsn,tgsb,cc,cg,cgsbi,  &
           sgfg,sdens,fluxd,fluxmax,fluxmelt,fluxef,dtsn,swg,ghf,shf,  &
           smelt,refreez,xmelt,xfrez,ymelt,yfrez,zmelt,zfrez,roff,     &
           srtes_qst)

      if (prin) then
         print *,'smelt: ',smelt(ipoint)
         print *,'refreez: ',refreez(ipoint)
         print *,'xmelt: ',xmelt(ipoint,2)
         print *,'xfrez: ',xfrez(ipoint,2)
         print *,'ymelt: ',ymelt(ipoint,2)
         print *,'yfrez: ',yfrez(ipoint,2)
         print *,'zmelt: ',zmelt(ipoint)
         print *,'zfrez: ',zfrez(ipoint)
         print *,'roff: ',roff(ipoint)
      endif

      do i = 1,istrip
         currtotal(i) = (((w(i,1) * zdepth(i,1)) + (w(i,2) *           &
                        zdepth(i,2)) + (w(i,3) * zdepth(i,3))) *       &
                        poros(i)) + capac(i,1) + capac(i,2) + roff(i)  &
                        + (eintmass(i) / 1000.0) - ((ppl(i) + ppc(i))  &
                        / 1000.0)
         if ((abs(starttotal(i)-currtotal(i)).gt.0.00001).or.         &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'SNOWMT CAPAC DIFF: ',(currtotal(i)-starttotal(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currtotal: ',currtotal(i)
            print *,'starttotal: ',starttotal(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'w_1: ',w(i,1)
            print *,'w_2: ',w(i,2)
            print *,'w_3: ',w(i,3)
            print *,'waterinlayer_1: ',(w(i,1) * poros(i) * zdepth(i,1))
            print *,'waterinlayer_2: ',(w(i,2) * poros(i) * zdepth(i,2))
            print *,'waterinlayer_3: ',(w(i,3) * poros(i) * zdepth(i,3))
            print *,'ppl: ',(ppl(i)/1000.)
            print *,'ppc: ',(ppc(i)/1000.)
            print *,'roff: ',roff(i)
            print *,'eintmass: ',(eintmass(i)/1000.)
            print *,'tm: ',tm(i)
            print *,'tc: ',tc(i)
            print *,'tg: ',tg(i)
            print *,'tsn: ',tsn(i)
            print *,'tgsb: ',tgsb(i)
         endif
         if (abs(currtotal(i)-starttotal(i)).gt.0.00001) then
            print *,'snowmt capac error!'
!            stop
         endif
         currsnow(i) = snow(i,1) + snow(i,2) - sfall(i) + zmelt(i) +   &
                       zfrez(i) + sevap(i,1) + sevap(i,2) + sevap(i,3)
         if ((abs(startsnow(i)-currsnow(i)).gt.0.00001).or.           &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'SNOWMT SNOW DIFF: ',(currsnow(i)-startsnow(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currsnow: ',currsnow(i)
            print *,'startsnow: ',startsnow(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'sfall: ',sfall(i)
            print *,'zmelt: ',zmelt(i)
            print *,'zfrez: ',zfrez(i)
            print *,'sevap_1: ',sevap(i,1)
            print *,'sevap_2: ',sevap(i,2)
            print *,'sevap_3: ',sevap(i,3)
         endif
         if (abs(currsnow(i)-startsnow(i)).gt.0.00001) then
            print *,'snowmt snow error!'
!            stop
         endif
      enddo

! Update deep soil temperature using effective soil heat flux
      timcon = dtc3x / (2.0 * sqrt(pie * 365.0))
      do i = 1,istrip
         if (sgfg(i).gt.0.5) then
            dtd(i) = fluxef(i) / cg(i) * dtc3x
            td(i) = td(i) + dtd(i)
            dgdt(i) = dgdt(i) + (cg(i) * dtd(i) / dtc3x)
         else
            td(i) = td(i) + (fluxef(i) / (1.0 / cgsbi(i)) * timcon)
         endif
! Coefficient is equal to: 2pie over (86400 * 365).  SMP3, Eq. 6
         td(i) = td(i) + (1.992385e-7 * (tdd(i) - td(i))) * dtc3x
      enddo

! Bare soil evaporation loss
      hlat3i = 1. / (hlat * 1000.0)
      do i = 1,istrip
         qqq(i,1) = 0.0
         qqq(i,2) = 0.0
         qqq(i,3) = 0.0
! If there is snow above the threshold, bare soil evaporation loss
! occurs from the snow cover and not the top soil layer.  12/28/97 dmm
         if ((sgfg(i).gt.0.5).and.(capac(i,2).ge.(egs(i)*hlat3i))) then
            capac(i,2) = capac(i,2) - (egs(i) * hlat3i)
            if (egs(i).ge.0.0) then
               txsc(i) = max((snow(i,2)-capac(i,2)),0.0)
            else
               if (tgsb(i).le.tf) then
                  txsc(i) = egs(i) * hlat3i
               else
                  txsc(i) = 0.0
               endif
            endif
            snow(i,2) = snow(i,2) - txsc(i)
            sevap(i,2) = txsc(i)
         else
            summ = ewt(i,1) + ewt(i,2) + ewt(i,3)
            if ((summ.le.1.1e-15).or.(td(i).lt.tf_drain)) then
               w(i,1) = w(i,1) - (egs(i) * hlat3i / (poros(i)          &
                        * zdepth(i,1)))
            else
               do il = 1,2
                  w(i,il) = w(i,il) - (egs(i) * hlat3i * ewt(i,il) /   &
                            (summ * poros(i) * zdepth(i,il)))
               enddo
               qqq(i,3) = egs(i) * hlat3i * ewt(i,3) / (summ * dtc3x)
            endif
         endif
      enddo

      do i = 1,istrip
         currtotal(i) = (((w(i,1) * zdepth(i,1)) + (w(i,2) *           &
                        zdepth(i,2)) + (w(i,3) * zdepth(i,3))) *       &
                        poros(i)) + capac(i,1) + capac(i,2) + roff(i)  &
                        + ((eintmass(i) + ebaremass(i) - (qqq(i,3)     &
                        * 1000. * dtc3x)) / 1000.0) - ((ppl(i) +       &
                        ppc(i)) / 1000.0)
         if ((abs(starttotal(i)-currtotal(i)).gt.0.00001).or.         &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'BARESOIL CAPAC DIFF: ',(currtotal(i)-starttotal(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currtotal: ',currtotal(i)
            print *,'starttotal: ',starttotal(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'w_1: ',w(i,1)
            print *,'w_2: ',w(i,2)
            print *,'w_3: ',w(i,3)
            print *,'waterinlayer_1: ',(w(i,1) * poros(i) * zdepth(i,1))
            print *,'waterinlayer_2: ',(w(i,2) * poros(i) * zdepth(i,2))
            print *,'waterinlayer_3: ',(w(i,3) * poros(i) * zdepth(i,3))
            print *,'ppl: ',(ppl(i)/1000.)
            print *,'ppc: ',(ppc(i)/1000.)
            print *,'roff: ',roff(i)
            print *,'eintmass: ',(eintmass(i)/1000.)
            print *,'ebaremass: ',(ebaremass(i)/1000.)
            print *,'tm: ',tm(i)
            print *,'tc: ',tc(i)
            print *,'tg: ',tg(i)
            print *,'tsn: ',tsn(i)
            print *,'tgsb: ',tgsb(i)
         endif
         if (abs(currtotal(i)-starttotal(i)).gt.0.00001) then
            print *,'baresoil capac error!'
!            stop
         endif
         currsnow(i) = snow(i,1) + snow(i,2) - sfall(i) + zmelt(i) +   &
                       zfrez(i) + sevap(i,1) + sevap(i,2) + sevap(i,3)
         if ((abs(startsnow(i)-currsnow(i)).gt.0.00001).or.           &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'BARESOIL SNOW DIFF: ',(currsnow(i)-startsnow(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currsnow: ',currsnow(i)
            print *,'startsnow: ',startsnow(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'sfall: ',sfall(i)
            print *,'zmelt: ',zmelt(i)
            print *,'zfrez: ',zfrez(i)
            print *,'sevap_1: ',sevap(i,1)
            print *,'sevap_2: ',sevap(i,2)
            print *,'sevap_3: ',sevap(i,3)
         endif
         if (abs(currsnow(i)-startsnow(i)).gt.0.00001) then
            print *,'baresoil snow error!'
!            stop
         endif
      enddo

! Extraction of transpiration loss from root zone
      do iveg = 1,2
         if (iveg.eq.1) then
            do i = 1,istrip
               absoil(i) = ect(i) * hlat3i
            enddo
         else
            do i = 1,istrip
               absoil(i) = egt(i) * hlat3i
            enddo
         endif
         do i = 1,istrip
            ef(i,2) = 0.0
            ef(i,3) = 0.0
            totdep(i) = zdepth(i,1)
            ef(i,1) = zdepth(i,1)
         enddo
         do il = 2,3
            do i = 1,istrip
               totdep(i) = totdep(i) + zdepth(i,il)
               dep(i) = max(0.0,rootd(i,iveg)- totdep(i) + zdepth(i,il))
               dep(i) = min(dep(i), zdepth(i,il))
               ef(i,il) = dep(i)
            enddo
         enddo
! SMP3, Eq. 16
         do il = 1,3
            do i = 1,istrip
               wait(i,il) = max(((phsoil(i,il)-psilow(i))*ef(i,il)),0.0)
            enddo
         enddo
         do i = 1,istrip
            summ = wait(i,1) + wait(i,2) + wait(i,3)
            if ((summ.le.1.1e-15).or.(td(i).lt.tf_drain)) then
               w(i,1) = w(i,1) - (absoil(i) / (poros(i) * zdepth(i,1)))
            else
               do il = 1,2
                  w(i,il) = w(i,il) - ((wait(i,il) * absoil(i)) /      &
                            (summ * poros(i) * zdepth(i,il)))
               enddo
               qqq(i,3) = qqq(i,3) + (wait(i,3) * absoil(i) / (summ *  &
                          dtc3x))
            endif
         enddo
      enddo

      do i = 1,istrip
         currtotal(i) = (((w(i,1) * zdepth(i,1)) + (w(i,2) *           &
                        zdepth(i,2)) + (w(i,3) * zdepth(i,3))) *       &
                        poros(i)) + capac(i,1) + capac(i,2) + roff(i)  &
                        + ((etmass(i) - (qqq(i,3) * 1000. * dtc3x))    &
                        / 1000.0) - ((ppl(i) + ppc(i)) / 1000.0)
         if ((abs(starttotal(i)-currtotal(i)).gt.0.00001).or.         &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'TRANS CAPAC DIFF: ',(currtotal(i)-starttotal(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currtotal: ',currtotal(i)
            print *,'starttotal: ',starttotal(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'w_1: ',w(i,1)
            print *,'w_2: ',w(i,2)
            print *,'w_3: ',w(i,3)
            print *,'waterinlayer_1: ',(w(i,1) * poros(i) * zdepth(i,1))
            print *,'waterinlayer_2: ',(w(i,2) * poros(i) * zdepth(i,2))
            print *,'waterinlayer_3: ',(w(i,3) * poros(i) * zdepth(i,3))
            print *,'ppl: ',(ppl(i)/1000.)
            print *,'ppc: ',(ppc(i)/1000.)
            print *,'roff: ',roff(i)
            print *,'etmass: ',etmass(i)
            print *,'qqq: ',(qqq(i,3)*1000.*dtc3x)
            print *,'tm: ',tm(i)
            print *,'tc: ',tc(i)
            print *,'tg: ',tg(i)
            print *,'tsn: ',tsn(i)
            print *,'tgsb: ',tgsb(i)
         endif
         if (abs(currtotal(i)-starttotal(i)).gt.0.00001) then
            print *,'trans capac error!'
!            stop
         endif
         currsnow(i) = snow(i,1) + snow(i,2) - sfall(i) + zmelt(i) +   &
                       zfrez(i) + sevap(i,1) + sevap(i,2) + sevap(i,3)
         if ((abs(startsnow(i)-currsnow(i)).gt.0.00001).or.           &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'TRANS SNOW DIFF: ',(currsnow(i)-startsnow(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currsnow: ',currsnow(i)
            print *,'startsnow: ',startsnow(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'sfall: ',sfall(i)
            print *,'zmelt: ',zmelt(i)
            print *,'zfrez: ',zfrez(i)
            print *,'sevap_1: ',sevap(i,1)
            print *,'sevap_2: ',sevap(i,2)
            print *,'sevap_3: ',sevap(i,3)
         endif
         if (abs(currsnow(i)-startsnow(i)).gt.0.00001) then
            print *,'trans snow error!'
!            stop
         endif
      enddo

! Interflow, infiltration excess and loss to groundwater
! All losses are assigned to variable 'roff'
      do il = 1,2
         do i = 1,istrip
            if (w(i,il).le.0.0) then
               w(i,il+1) = w(i,il+1) + w(i,il) * zdepth(i,il) /        &
                           zdepth(i,il+1)
               w(i,il) = 0.0
            endif
         enddo
      enddo

      call hyssib_runoff(istrip,dtc3x,tf_drain,baseflow,ityp,phsat,    &
           poros,bee,satco,slope,wiltp,zdepth,w,td,q3g,qqq,wattabl,    &
           roffg,roff,prin)

      do i = 1,istrip
         currtotal(i) = (((w(i,1) * zdepth(i,1)) + (w(i,2) *           &
                        zdepth(i,2)) + (w(i,3) * zdepth(i,3))) *       &
                        poros(i)) + capac(i,1) + capac(i,2) + roff(i)  &
                        + (etmass(i) / 1000.0) - ((ppl(i) + ppc(i))    &
                        / 1000.0) - wattabl(i)
         if ((abs(starttotal(i)-currtotal(i)).gt.0.00001).or.         &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'RUNOFF CAPAC DIFF: ',(currtotal(i)-starttotal(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currtotal: ',currtotal(i)
            print *,'starttotal: ',starttotal(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'w_1: ',w(i,1)
            print *,'w_2: ',w(i,2)
            print *,'w_3: ',w(i,3)
            print *,'waterinlayer_1: ',(w(i,1) * poros(i) * zdepth(i,1))
            print *,'waterinlayer_2: ',(w(i,2) * poros(i) * zdepth(i,2))
            print *,'waterinlayer_3: ',(w(i,3) * poros(i) * zdepth(i,3))
            print *,'ppl: ',(ppl(i)/1000.)
            print *,'ppc: ',(ppc(i)/1000.)
            print *,'roff: ',roff(i)
            print *,'q3g: ',(q3g(i)*dtc3x)
            print *,'roffg: ',roffg(i)
            print *,'etmass: ',etmass(i)
            print *,'wattabl: ',wattabl(i)
            print *,'tm: ',tm(i)
            print *,'tc: ',tc(i)
            print *,'tg: ',tg(i)
            print *,'tsn: ',tsn(i)
            print *,'tgsb: ',tgsb(i)
         endif
         if (abs(currtotal(i)-starttotal(i)).gt.0.00001) then
            print *,'runoff capac error!'
!            stop
         endif
         currsnow(i) = snow(i,1) + snow(i,2) - sfall(i) + zmelt(i) +   &
                       zfrez(i) + sevap(i,1) + sevap(i,2) + sevap(i,3)
         if ((abs(startsnow(i)-currsnow(i)).gt.0.00001).or.           &
            (prin.and.(i.eq.ipoint))) then
            print *,' '
            print *,'RUNOFF SNOW DIFF: ',(currsnow(i)-startsnow(i))
            print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
            print *,'currsnow: ',currsnow(i)
            print *,'startsnow: ',startsnow(i)
            print *,'snow_1: ',snow(i,1)
            print *,'snow_2: ',snow(i,2)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'sfall: ',sfall(i)
            print *,'zmelt: ',zmelt(i)
            print *,'zfrez: ',zfrez(i)
            print *,'sevap_1: ',sevap(i,1)
            print *,'sevap_2: ',sevap(i,2)
            print *,'sevap_3: ',sevap(i,3)
         endif
         if (abs(currsnow(i)-startsnow(i)).gt.0.00001) then
            print *,'runoff snow error!'
!            stop
         endif
      enddo

      do i = 1,istrip
         if (w(i,1).gt.1.0) then
            w(i,2) = w(i,2) + (w(i,1) - 1.0) * zdepth(i,1) / zdepth(i,2)
            w(i,1) = 1.0
         endif
         if (w(i,2).gt.1.0) then
            w(i,3) = w(i,3) + (w(i,2) - 1.0) * zdepth(i,2) / zdepth(i,3)
            w(i,2) = 1.0
         endif
         if (w(i,3).gt.1.0) then
            roff(i) = roff(i) + (w(i,3) - 1.0) * poros(i) * zdepth(i,3)
            w(i,3) = 1.0
         endif
      enddo

! Increment prognostic variables
! Adjust theta and sh to be consistent with dew formation

! Solve implicit system for winds
      do i = 1,istrip
         drag(i) = ros(i) * cu(i) * ustar(i)
         aaa(i) = drag(i) * 0.01 * grav / psb(i)
         gmu(i,2) = gmu(i,2) + dtc3x * aaa(i)
         gmu(i,3) = (gmu(i,3) - aaa(i) * spdm(i)) / gmu(i,2)
         gmu(i,4) = (gmu(i,4) - aaa(i) * 0.0    ) / gmu(i,2)
      enddo

! Calculate soil water budget after calling pbl and compare with
! previous budget
      do i = 1,istrip
! Water balance minus precip plus evap plus runoff.  10/31/97 dmm
         currtotal(i) = (((w(i,1) * zdepth(i,1)) + (w(i,2) *           &
                        zdepth(i,2)) + (w(i,3) * zdepth(i,3))) *       &
                        poros(i)) + capac(i,1) + capac(i,2) + roff(i)  &
                        + (etmass(i) / 1000.0) - ((ppl(i) + ppc(i))    &
                        / 1000.0) - wattabl(i)
         if ((abs(starttotal(i)-currtotal(i)).gt.0.00001).or.         &
            (prin.and.(i.eq.ipoint))) then
! Printing out more variables if water imbalance.  10/23/97 dmm
            print *,' WBdiff:  ',(starttotal(i)-currtotal(i))
            print *,'starttotal: ',starttotal(i)
            print *,'currtotal: ',currtotal(i)
            print *,'latco: ',latco(i)
            print *,'lonco: ',lonco(i)
            print *,'ityp: ',ityp(i)
            print *,'w_1: ',w(i,1)
            print *,'w_2: ',w(i,2)
            print *,'w_3: ',w(i,3)
            print *,'capac_1: ',capac(i,1)
            print *,'capac_2: ',capac(i,2)
            print *,'ppl: ',ppl(i)
            print *,'ppc: ',ppc(i)
            print *,'etmass: ',etmass(i)
            print *,'roff: ',roff(i)
            print *,'tc: ',tc(i)
            print *,'tg: ',tg(i)
            print *,'td: ',td(i)
            print *,'tm: ',tm(i)
            print *,'tdd: ',tdd(i)
         endif
      enddo

! Calculate and compare energy budgets
      do i = 1,istrip
         cbal(i) = radt(i,1) - chf(i) - (ect(i) + hc(i) + eci(i)) /dtc3x
! Calculating snow and ground energy balance different for snow or
! non-snow covered conditions.  10/23/97 dmm
         if (sgfg(i).gt.0.5) then
! Ground energy balance is shortwave thru snow to ground, minus ground
! heat flux, minus heat flux to snow
            gbal(i) = swg(i) - ghf(i) - dgdt(i)
! Snow energy balance is net radiation, minus latent and sensible heat
! fluxes, minus/plus melt/freeze energy, minus change in heat capacity
! of snow, plus precipitation heat, plus flux from ground
            sbal(i) = radt(i,2) - shf(i) - ((egt(i) + egi(i) + hg(i) + &
                      egs(i)) / dtc3x) - fluxs(i)
            bbal(i) = sws(i) + fluxs(i) + fluxd(i) - bhf(i)
         else
            gbal(i) = radt(i,2) - ghf(i) - ((egt(i) + egi(i) + hg(i) + &
                      egs(i)) / dtc3x) - dgdt(i)
            sbal(i) = 0.0
            bbal(i) = 0.0
         endif
      enddo

      do i = 1,istrip
         if (((abs(cbal(i)).gt.0.5).or.(abs(gbal(i)).gt.0.5).or.       &
            ((abs(sbal(i)).gt.5.0).and.(snow(i,2).ge.thres))).or.      &
            (prin.and.(i.eq.ipoint))) then
            print *,'sgfg: ',sgfg(i)
            print *,'radt_1: ',radt(i,1)
            print *,'radt_2: ',radt(i,2)
            print *,'chf: ',chf(i)
            print *,'ghf: ',ghf(i)
            print *,'shf: ',shf(i)
            print *,'hflux: ',hflux(i)
            print *,'ect: ',ect(i)
            print *,'eci: ',eci(i)
            print *,'egt: ',egt(i)
            print *,'egi: ',egi(i)
            print *,'egs: ',egs(i)
            print *,'dgdt: ',dgdt(i)
            if (abs(sbal(i)).ge.1e-9) then
               print *,'sbal:  ',sbal(i)
            endif
            if (abs(gbal(i)).ge.1e-10) then
               print *,'gbal:  ',gbal(i)
            endif
            if (abs(bbal(i)).ge.1e-10) then
               print *,'bbal:  ',bbal(i)
            endif
            if (abs(cbal(i)).ge.1e-10) then
               print *,'cbal:  ',cbal(i)
               print *,'End of time step'
               print *,' '
            endif
!            stop
         endif
      enddo

! Return values to driver
      do i = 1,istrip
         gw1(i) = w(i,1)
         gw2(i) = w(i,2)
         gw3(i) = w(i,3)
         snow(i,1) = min(snow(i,1),capac(i,1))
         snow(i,2) = min(snow(i,2),capac(i,2))
         sn1(i) = snow(i,1) * 1000.0
         sn2(i) = snow(i,2) * 1000.0
         capac1(i) = max((capac(i,1) * 1000.0 - sn1(i)),0.0)
         capac2(i) = max((capac(i,2) * 1000.0 - sn2(i)),0.0)
         etmass(i) = etmass(i) / dtc3x
         hflux(i) = hflux(i) / dtc3x
! ALMA General Energy Balance Components
         swnet(i) = anetsw(i,1) + anetsw(i,2)
         lwnet(i) = anetlw(i,1) + anetlw(i,2)
         qle(i) = (ect(i) + eci(i) + egt(i) + egi(i) + egs(i)) / dtc3x
         qg(i) = ghf(i) + dgdt(i)
         qf(i) = (zmelt(i) + zfrez(i)) * snomel / dtc3x
         qv(i) = (sevap(i,1) + sevap(i,2) + sevap(i,3)) * (snomel +    &
                 (hlat * 1000.)) / dtc3x
         qh(i) = hflux(i)
         qtau(i) = drag(i) * spdm(i)
         qa(i) = max(precheat(i),0.0)
         delsurfheat(i) = (chf(i) - fluxd(i)) * dtc3x
         delcoldcont(i) = (bhf(i) + shf(i) + precheat(i)               &
                                           - qf(i) - qv(i)) * dtc3x
! ALMA General Water Balance Components
         snowf(i) = sfall(i) * 1000. / dtc3x
         rainf(i) = max((ppl(i) + ppc(i) - (sfall(i) * 1000.)),0.0)    &
                    / dtc3x
         evap(i) = etmass(i)
         qs(i) = max((roff(i) - (q3g(i) * dtc3x) - roffg(i)),0.0)      &
                 * 1000. / dtc3x
         qrec(i) = wattabl(i) * 1000. / dtc3x
         qsb(i) = ((q3g(i) * dtc3x) + roffg(i) - wattabl(i)) * 1000.   &
                  / dtc3x
         qsm(i) = (smelt(i) + xmelt(i,2) + ymelt(i,2)) * 1000. / dtc3x
         qfz(i) = -(refreez(i) + xfrez(i,2) + yfrez(i,2)) * 1000.      &
                  / dtc3x
         qst(i) = srtes_qst(i) * 1000. / dtc3x
         delsoilmoist(i) = ((((w(i,1) * zdepth(i,1)) + (w(i,2) *       &
                           zdepth(i,2)) + (w(i,3) * zdepth(i,3))) *    &
                           poros(i)) - startsm(i)) * 1000.
         delswe(i) = (capac(i,2) - startswe(i)) * 1000.
         delintercept(i) = (capac(i,1) - startint(i)) * 1000.
! ALMA Surface State Variables
         snowt(i) = tsn(i)
         vegtc(i) = tc(i)
         baresoilt(i) = tgsb(i)
         avgsurft(i) = (tgsb(i) * (1.0 - vcover(i,1))) +               &
                       (tc(i) *          vcover(i,1))
         radteff(i) = tgeff(i)
         swe(i) = capac(i,2) * 1000.
         sweveg(i) = snow(i,1) * 1000.
! ALMA SubSurface State Variables
         soilmoist1(i) = w(i,1) * zdepth(i,1) * poros(i) * 1000.
         soilmoist2(i) = w(i,2) * zdepth(i,2) * poros(i) * 1000.
         soilmoist3(i) = w(i,3) * zdepth(i,3) * poros(i) * 1000.
         soiltemp(i) = td(i)
         soilwet(i) = (zdepth(i,1) + zdepth(i,2) + zdepth(i,3))*wiltp(i)
         soilwet(i) = (((soilmoist1(i)+ soilmoist2(i) + soilmoist3(i)) &
                      / poros(i) /1000.) - soilwet(i)) / ((zdepth(i,1) &
                      + zdepth(i,2) + zdepth(i,3)) - soilwet(i))
! ALMA Evaporation Components
         potevap(i) = epot(i) / hlat / dtc3x
         ecanop(i) = (eintmass(i) - (sevap(i,1) * 1000.)) / dtc3x
         tveg(i) = (etranmass(i) - ((sevap(i,2) + sevap(i,3)) *1000.)) &
                   / dtc3x
         esoil(i) = ebaremass(i) / dtc3x
         rootmoist(i) = soilmoist1(i) + soilmoist2(i)
         canopint(i) = capac(i,1) * 1000.
         subsnow(i) = (sevap(i,2) + sevap(i,3)) * 1000. / dtc3x
         subsurf(i) = sevap(i,1) * 1000. / dtc3x
         acond(i) = 1.0 / ra(i)
         ccond(i) = 1.0 / rst(i,1)
! ALMA Cold Season Processes
         if ((snow(i,2).gt.sskin).and.(sdens(i).gt.0.0)) then
            snowdepth(i) = ((sskin / 119.6) + ((snow(i,2) - sskin)     &
                           / sdens(i))) * 1000.
         else
            snowdepth(i) = snow(i,2) / 119.6 * 1000.
         endif
         if (capac(i,2).gt.0.0) then
            sliqfrac(i) = (capac(i,2) - snow(i,2)) / capac(i,2)
         else
            sliqfrac(i) = 0.0
         endif
! Additional diagnostics
         ecanopt(i) = etrancmass(i) / dtc3x
         ecanopi(i) = eintcmass(i) / dtc3x
         eground(i) = eintgmass(i) / dtc3x
         hfluxc(i) = hc(i) / dtc3x
         hfluxg(i) = hg(i) / dtc3x
         watcan(i) = (capac(i,1) - snow(i,1)) * 1000.
         watgrd(i) = (capac(i,2) - snow(i,2)) * 1000.
         snocan(i) = snow(i,1) * 1000.
         snogrd(i) = snow(i,2) * 1000.
         wet1(i) = w(i,1)
         wet2(i) = w(i,2)
         wet3(i) = w(i,3)
      enddo

      return
      end

! ----------------------------------------------------------------------
! - END MAIN PROGRAM ---------------------------------------------------
! ----------------------------------------------------------------------

!***********************************************************************
      subroutine root(istrip,w,phsat,bee,phsoil)

      implicit none
!***********************************************************************
!***  soil moisture potentials in each layer.  SMP3, Eq. 14
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! w(3)         (0-1)     wetness of surface/root/recharge zone
! phsat        (m)       soil moisture potential at saturation
! bee          (-)       Clapp-Hornberger empirical constant
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! phsoil(3)    (m)       soil moisture potential of the i-th soil layer
!-----------------------------------------------------------------------
      integer istrip, isoil, i
      real w(istrip,3), phsat(istrip), bee(istrip)
      real phsoil(istrip,3), www

      do isoil = 1,3
         do i = 1,istrip
            www = max(w(i,isoil), 0.0001)
            phsoil(i,isoil) = phsat(i) * exp(-bee(i) * log(www))
         enddo
      enddo

      return
      end

!***********************************************************************
      subroutine raduse_hyssib(istrip,prin,ipoint,dkfac,stefan,sskin,  &
           vcover,thermk,snow,tc,tgsb,sgfg,radn,radc3,radt,par,pd,sws, &
           swg)

      implicit none
!***********************************************************************
!***  absorption of radiation by surface
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! dkfac        (m-1)     solar extinction coefficient through snow
! stefan    (W m-2 K-4)  Stefan-Boltzmann constant
! sskin        (m)       depth of diurnal snow layer
! vcover(cg)   (0-1)     fraction of vegetation cover [fpar/greenness]
! thermk       (0-1)     canopy emissivity
! snow(cg)     (m)       liquid equiv. snow on canopy/ground
! tc           (K)       canopy temperature
! tgsb         (K)       blended ground/snow temperature
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
! radn(vnt,bd) (W m-2)   downward sw/lw radiation at the surface
! radc3(cg)    (W m-2)   absorbed sw/lw radiation at the surface
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! radt(cg)     (W m-2)   net radiation
! par          (W m-2)   photosynthetically active radiation on canopy
! pd           (0-1)     ratio of par beam to total par
! sws          (W m-2)   sw radiation absorbed in bulk snow, SMP3, Eq.A5
! swg          (W m-2)   sw radiation absorbed in ground, SMP3, Eq. A6
!-----------------------------------------------------------------------
!     in subroutine parameters
!-----------------------------------------------------------------------
! trans(2)     (0-1)     SW trans. of diurnal[1]; bulk[2] snow,SMP1,Eq.1
! swc          (W m-2)   sw radiation absorbed in canopy
! swk          (W m-2)   sw radiation absorbed in diurnalsnow,SMP3,Eq.A4
!-----------------------------------------------------------------------
      logical prin
      integer istrip, ipoint, i
      real dkfac, stefan, sskin
      real vcover(istrip,2), thermk(istrip)
      real snow(istrip,2), tc(istrip), tgsb(istrip), sgfg(istrip)
      real radn(istrip,3,2), radc3(istrip,2), radt(istrip,2)
      real par(istrip), pd(istrip)
      real sws(istrip), swg(istrip)
      real swk(istrip), swc(istrip)
      real tc4(istrip), tg4(istrip)
      real fac1(istrip), fac2(istrip), closs(istrip), gloss(istrip)
      real trans(istrip,2)

      do i = 1,istrip
         tc4(i) = tc(i) * tc(i) * tc(i) * tc(i)
         tg4(i) = tgsb(i) * tgsb(i) * tgsb(i) * tgsb(i)
         fac1(i) = vcover(i,1) * (1.0 - thermk(i))
         fac2(i) = 1.0
         closs(i) = 2.0 * fac1(i) * stefan * tc4(i)
         closs(i) = closs(i) - fac2(i) * fac1(i) * stefan * tg4(i)
         gloss(i) = fac2(i) * stefan * tg4(i)
         gloss(i) = gloss(i) - fac1(i) * fac2(i) * stefan * tc4(i)
         if (prin.and.(i.eq.ipoint)) then
            print *,'closs: ',closs(i)
            print *,'gloss: ',gloss(i)
         endif
      enddo

      do i = 1,istrip
         if (sgfg(i).gt.0.5) then
            trans(i,1) = exp(-dkfac * sskin)
            trans(i,2) = exp(-dkfac * (snow(i,2) - sskin))
            swk(i) = radc3(i,2) * (1.0 - trans(i,1))
            sws(i) = (radc3(i,2) - swk(i)) * (1.0 - trans(i,2))
            swg(i) = radc3(i,2) - (swk(i) + sws(i))
            swc(i) = radc3(i,1)
            radt(i,1) = radc3(i,1) + radn(i,3,2) * vcover(i,1) *       &
                        (1.0 - thermk(i)) - closs(i)
            radt(i,2) = swk(i) + radn(i,3,2) * (1.0 - (vcover(i,1) *   &
                        (1.0 - thermk(i)))) - gloss(i)
         else
            trans(i,1) = -999.
            trans(i,2) = -999.
            swk(i) = 0.0
            swg(i) = 0.0
            sws(i) = 0.0
            swc(i) = 0.0
            radt(i,1) = radc3(i,1) + radn(i,3,2) *                     &
                         vcover(i,1) * (1.0 - thermk(i))   - closs(i)
            radt(i,2) = radc3(i,2) + radn(i,3,2) * (1.0 -              &
                        (vcover(i,1) * (1.0 - thermk(i)))) - gloss(i)
         endif
         if (prin.and.(i.eq.ipoint)) then
            print *,'swk: ',swk(i)
            print *,'swc: ',swc(i)
         endif
         par(i) = radn(i,1,1) + radn(i,1,2) + 0.001
         pd(i) = (radn(i,1,1) + 0.001) / par(i)
         if (par(i).le.0.00001) par(i) = 0.00001
      enddo

      return
      end

!***********************************************************************
      subroutine stomat(istrip,pie,ityp,zlt,green,vcover,chil,rstpar,  &
           extk,cosz,par,pd,rst)

      implicit none
!***********************************************************************
!***  stomatal resistance
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! ityp         (1-13)    vegetation index
! zlt(cg)      (-)       leaf area index
! green        (0-1)     green leaf fraction
! vcover(cg)   (0-1)     fraction of vegetation cover [fpar/greenness]
! chil         (-)       leaf angle distribution factor
! rstpar(3)    (-)       par influence on stomatal resist. coefficients
! extk(cg,vnt,bd)  (-)   extinction coefficient
! cosz         (radians) cosine of solar zenith angle
! par          (W m-2)   photosynthetically active radiation on canopy
! pd           (0-1)     ratio of par beam to total par
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! rst(cg)      (sec m-1) stomatal resistance
!-----------------------------------------------------------------------
      integer istrip, iveg, irad, i, ntyp
      real pie, fcon, xabc, xabd, ftemp
      integer ityp(istrip)
      real zlt(istrip,2), green(istrip), vcover(istrip,2)
      real chil(istrip), rstpar(istrip,3), extk(istrip,2,3,2)
      real cosz(istrip), par(istrip), pd(istrip), rst(istrip,2)
      real f(istrip), gamma(istrip), at(istrip), power1(istrip)
      real power2(istrip), aa(istrip), bb(istrip), zat(istrip)
      real zk(istrip), ekat(istrip), rho4(istrip), avflux(istrip)

! Define arrays to modify rst (see modification notes below).  rmin
! values from GEWEX/GSWP, 1995, Table 4; b2 values from Offline SSiB
! sheet.  Bogus values for rmin(13) and b2(11) are set to the values of
! rmin(12) and b2(10) respectively.  This shouldn't affect anything as
! rst for ntyp=11 or 13 is just set to 1.0e5.  5/8/96 dmm
      real rmin(13), b2(13)
      data rmin /80.,  80., 100., 120., 120., 110., 110., 110.,  80.,  &
                 80., 110.,  80.,  80./
      data b2 /162.147, 208.546, 226.241, 242.829, 242.829, 282.207,   &
               117.803, 1048.95, 1048.95, 208.546, 208.546, 208.546,   &
                41.47/

! Bounding of product of extinction coefficient and local LAI
      do i = 1,istrip
         f(i) = max(cosz(i), 0.01746)
      enddo

      do iveg = 1,2
         do irad = 1,2
            do i = 1,istrip
               ntyp = ityp(i)
               extk(i,iveg,1,irad) = min(extk(i,iveg,1,irad), (150.0/  &
                                     zlt(i,iveg)*vcover(i,iveg)))
            enddo
         enddo
      enddo

      fcon = 0.25 * pie + (1.0 / 3.0)

! "iveg" must be equal to one only!!!
      iveg = 1
      do i = 1,istrip
         ntyp = ityp(i)
         if (ntyp.eq.13) then
            rst(i,iveg) = 1.0e10
            goto 250
         endif

         at(i) = zlt(i,iveg) / vcover(i,iveg)

         if (par(i).le.0.00101) then
            xabc = rstpar(i,1) / rstpar(i,2) + rstpar(i,3)
            xabd = .5 / xabc * at(i)
            rst(i,iveg) = 1. / xabd
         else
            gamma(i) = (rstpar(i,1) + rstpar(i,2) * rstpar(i,3)) /     &
                        rstpar(i,3)
! Single extinction coefficient using weighted
! values of direct and diffus contributions to par
            power1(i) = at(i) * extk(i,iveg,1,1)
            power2(i) = at(i) * extk(i,iveg,1,2)
            aa(i) = 0.5 - (0.633 + 0.33 * chil(i)) * chil(i)
            bb(i) = 0.877 - 1.754 * aa(i)
            zat(i) = log((exp(-power2(i)) + 1.0) * 0.5) * pd(i) /      &
                     extk(i,iveg,1,1)
            zat(i) = zat(i) + log((exp(-power1(i)) + 1.0) * 0.5) *     &
                     (1.0 - pd(i)) / extk(i,iveg,1,2)
            zk(i) = 1.0 / zat(i) * log(pd(i) * exp(power1(i) * zat(i)  &
                    / at(i)) + (1.0 - pd(i)) * exp(power2(i) * zat(i)  &
                    / at(i)))
! Canopy and ground cover bulk resistances using ross-goudriaan leaf
! function, total par flux (avflux) and mean extinction coefficient (zk)
            ftemp = min(zk(i)*at(i), 20.0)
            ekat(i) = exp(ftemp)
            avflux(i) = par(i) * (pd(i) * (aa(i) / f(i) + bb(i)) +     &
                       (1.0 - pd(i)) * (bb(i) * fcon + aa(i) * 1.5))
            rho4(i) = gamma(i) / avflux(i)
            rst(i,iveg) = rstpar(i,2) / gamma(i) * log((rho4(i) *      &
                          ekat(i) + 1.0) / (rho4(i) + 1.0))
            rst(i,iveg) = rst(i,iveg) - log((rho4(i) + 1.0 / ekat(i))  &
                          / (rho4(i) + 1.0))
            rst(i,iveg) = rst(i,iveg) / (zk(i) * rstpar(i,3))
            rst(i,iveg) = 1.0 / (rst(i,iveg) * green(i))
         endif
! Scale the value of rst by the ration of r_min to b2 - from Offline
! SSiB Vegetation parameter sheet by P Dirmeyer.  5/8/96 dmm
         rst(i,iveg) = rst(i,iveg) * rmin(ntyp) / b2(ntyp)
 250     continue
      enddo

      do i = 1,istrip
         rst(i,2) = 1.0e10
      enddo

      return
      end

!***********************************************************************
      subroutine interc_hyssib(istrip,dtc3x,snomel,clai,cw,cs,sskin,   &
           tf,tf_snow,poros,satco,zdepth,zlt,vcover,topostd,tm,ppl,    &
           ppc,tg,w,capac,snow,satcap,tc,tsn,sgfg,sdens,tgsb,cc,cg,    &
           cgsbi,extk,sfall,ymelt,yfrez,roff,precheat)

      implicit none
!***********************************************************************
!***  calculation of (1) interception and drainage of rainfall and snow
!                    (2) specific heat terms fixed for time step
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! dtc3x        (sec)     timestep
! snomel       (J m-3)   latent heat of fusion X water density
! clai      (J m-2 K-1)  heat capacity of foliage
! cw        (J m-3 K-1)  heat capacity of liquid water X water density
! cs        (J m-3 K-1)  heat capacity of frozen water X water density
! sskin        (m)       depth of diurnal snow layer
! tf           (K)       freezing temperature
! tf_snow      (K)       reference-level temp. below which precip = snow
! poros        (0-1)     soil porosity
! satco        (m/sec)   saturation hydraulic conductivity
! zdepth(3)    (m)       exact independent depth of 3 soil layers
! zlt(cg)      (-)       leaf area index
! vcover(cg)   (0-1)     fraction of vegetation cover [fpar/greenness]
! topostd      (m)       topography standard deviat. 1x1deg from GTOPO30
! tm           (K)       forcing reference-level temperature
! ppl          (mm/dt)   forcing large-scale precipitation
! ppc          (mm/dt)   forcing convective precipitation
! satcap(cg)   (m)       interception capacity
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
! sdens        (kg m-3)  density of bulk snow layer
! cg        (J m-2 K-1)  heat capacity of the ground
! extk(cg,vnt,bd)  (-)   extinction coefficient
! precheat     (W m-2)   heat transferred to snow cover by rain
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! tg           (K)       ground temperature
! w(3)         (0-1)     wetness of surface/root/recharge zone
! capac(cg)    (m)       liquid equiv. water&snow on canopy/ground
! snow(cg)     (m)       liquid equiv. snow on canopy/ground
! tc           (K)       canopy temperature
! tsn          (K)       snow on ground temperature
! tgsb         (K)       blended ground/snow temperature
! cc        (J m-2 K-1)  heat capacity of the canopy
! cgsbi    (J-1 m^2 K^1) inverse of heat capacity of blended ground/snow
! sfall        (m/dt)    snowfall
! ymelt(cg)    (m/dt)    melted snow on canopy/ground from precip
! yfrez(cg)    (m/dt)    frozen water on canopy/ground from precip
! roff         (m/dt)    total runoff
!-----------------------------------------------------------------------
!     in subroutine parameters
!-----------------------------------------------------------------------
! roffc        (m/dt)    canopy water runoff from topographic variance
! roffo        (m/dt)    instantaneous runoff from convective precip
! roffp        (m/dt)    liquid precip runoff from topographic variance
!-----------------------------------------------------------------------
      integer istrip, iveg, i
      real dtc3x, snomel, clai, cw, cs, sskin, tf, tf_snow
      real crtes, prtes, excess
      real poros(istrip), satco(istrip), zdepth(istrip,3)
      real zlt(istrip,2), vcover(istrip,2), topostd(istrip)
      real tm(istrip), ppl(istrip), ppc(istrip), tg(istrip)
      real w(istrip,3), capac(istrip,2), snow(istrip,2)
      real satcap(istrip,2), tc(istrip), tsn(istrip)
      real sgfg(istrip), sdens(istrip), tgsb(istrip)
      real cc(istrip), cg(istrip), cgsbi(istrip), extk(istrip,2,3,2)
      real sfall(istrip), ymelt(istrip,2), yfrez(istrip,2), roff(istrip)
      real thru(istrip), roffc(istrip), roffo(istrip), roffp(istrip)
      real freeze(istrip), precheat(istrip), spwet(istrip)
      real ap(istrip), cp(istrip), totalp(istrip), fpi(istrip)
      real p0(istrip), tsp(istrip), zload(istrip), cct(istrip)
      real xsc(istrip), tti(istrip), arg(istrip), tex(istrip)
      real tsd(istrip), pinf(istrip), tsf(istrip)
      real txsc(istrip), equdep(istrip), xs(istrip)
      real snowp(istrip,2), waterp(istrip,2), water(istrip,2)
      real pcoefs(2,2), bp
      data pcoefs(1,1) /5.0/, pcoefs(1,2) /6.737946e-3/
      data pcoefs(2,1) /0.0001/, pcoefs(2,2) /0.9999/, bp /5.0/

      do i = 1,istrip
         ap(i) = pcoefs(2,1)
         cp(i) = pcoefs(2,2)
         totalp(i) = ppc(i) + ppl(i)
         if (totalp(i).ge.1.0e-8) then
            ap(i) = (ppc(i) * pcoefs(1,1) + ppl(i) * pcoefs(2,1)) /    &
                    totalp(i)
            cp(i) = (ppc(i) * pcoefs(1,2) + ppl(i) * pcoefs(2,2)) /    &
                    totalp(i)
         endif
         p0(i) = totalp(i) * 0.001
         thru(i) = 0.0
         fpi(i) = 0.0
         roffc(i) = 0.0
         roffo(i) = 0.0
         roffp(i) = 0.0
         sfall(i) = 0.0
         precheat(i) = 0.0
         ymelt(i,1) = 0.0
         ymelt(i,2) = 0.0
         yfrez(i,1) = 0.0
         yfrez(i,2) = 0.0
      enddo

      do i = 1,istrip
         xsc(i) = max(0.0, (capac(i,1)-satcap(i,1)))
         capac(i,1) = capac(i,1) - xsc(i)
         txsc(i) = max((snow(i,1)-capac(i,1)),0.0)
         snow(i,1) = snow(i,1) - txsc(i)
         capac(i,2) = capac(i,2) + xsc(i)
         snow(i,2) = snow(i,2) + txsc(i)
         water(i,2) = capac(i,2) - snow(i,2)
         xsc(i) = max(0.0, (water(i,2)-satcap(i,2)))
! crtes = Canopywater Ready To Enter Soil.  This is the excess water on
! the canopy above the capacity.  Calculate additional overland flow
! before liquid gets into the soil.  Using linefit of basins, limiting
! to half of liquid to runoff.  The remaining liquid enters the soil.
! SMP3, Eq. 8
         crtes = xsc(i)
         roffc(i) = min((2.07 * topostd(i) / 1843.47), 0.5) * crtes
         roff(i) = roff(i) + roffc(i)
         w(i,1) = w(i,1) + (crtes - roffc(i)) / (poros(i) * zdepth(i,1))
         capac(i,2) = capac(i,2) - xsc(i)
         water(i,2) = capac(i,2) - snow(i,2)
      enddo

      do iveg = 1,2
         if (iveg.eq.1) then
            do i = 1,istrip
               cct(i) = zlt(i,iveg) * clai
               tsp(i) = tc(i)
            enddo
         else
            do i = 1,istrip
               cct(i) = cg(i) * (1.0 - sgfg(i))
               tsp(i) = (tg(i) * sgfg(i)) + (tgsb(i) * (1.0 - sgfg(i)))
            enddo
         endif

         do i = 1,istrip
            waterp(i,iveg) = capac(i,iveg) - snow(i,iveg)
            snowp(i,iveg) = snow(i,iveg)
            zload(i) = capac(i,iveg)
            fpi(i) = (1.0 - exp(-extk(i,iveg,3,1) * zlt(i,iveg) /      &
                     vcover(i,iveg))) * vcover(i,iveg)
            tti(i) = p0(i) * (1.0 - fpi(i))
         enddo

! Proportional saturated area (xs) and leaf drainage (tex)
         do i = 1,istrip
            xs(i) = 1.0
            if (p0(i).ge.1.0e-9) then
               arg(i) = (satcap(i,iveg) - zload(i)) / (p0(i) * fpi(i)  &
                        * ap(i)) - cp(i) / ap(i)
               if (arg(i).ge.1.0e-9) then
                  xs(i) = -1.0 / bp * log(arg(i))
                  xs(i) = min(xs(i), 1.0)
                  xs(i) = max(xs(i), 0.0)
               endif
            endif
         enddo

! Total throughfall (thru) and store augmentation.  SMP3, Eq. 7
         do i = 1,istrip
            tex(i) = p0(i) * fpi(i) * (ap(i) / bp * (1.0 - exp(-bp     &
                     * xs(i))) + cp(i) * xs(i)) - (satcap(i,iveg) -    &
                     zload(i)) * xs(i)
            tex(i) = max(tex(i), 0.0)
            thru(i) = tti(i) + tex(i)
            if ((iveg.eq.2).and.(tgsb(i).le.tf)) then
               thru(i) = 0.0
            endif
            thru(i) = min(thru(i), p0(i))
            pinf(i) = p0(i) - thru(i)
            capac(i,iveg) = capac(i,iveg) + pinf(i)
            if (tm(i).ge.tf_snow) then
               sfall(i) = 0.0
               if (snow(i,iveg).gt.0.0) precheat(i) = precheat(i) +    &
                  ((capac(i,iveg) - snowp(i,iveg) - waterp(i,iveg)) *  &
                  (cw * (tm(i) - tf))) / dtc3x
            else
               snow(i,iveg) = snow(i,iveg) + pinf(i)
               sfall(i) = sfall(i) + pinf(i)
            endif
         enddo

         if (iveg.eq.2) then
            do i = 1,istrip
               if (tm(i).lt.tf_snow) then
                  capac(i,iveg) = waterp(i,iveg) + snowp(i,iveg) + p0(i)
                  snow(i,iveg) = snowp(i,iveg) + p0(i)
                  sfall(i) = sfall(i) + thru(i)
                  thru(i) = 0.0
               else
                  capac(i,iveg) = waterp(i,iveg) + snowp(i,iveg) +     &
                                  (p0(i) - thru(i))
                  equdep(i) = satco(i) * dtc3x
                  xs(i) = 1.0
                  if (thru(i).ge.1.0e-9) then
                     arg(i) = equdep(i) / (thru(i)*ap(i)) - cp(i)/ap(i)
                     if (arg(i).ge.1.0e-9) then
                        xs(i) = -1.0 / bp * log(arg(i))
                        xs(i) = min(xs(i), 1.0)
                        xs(i) = max(xs(i), 0.0)
                     endif
                  endif
! Instantaneous overland flow contribution.  SMP3, Eq. 7
                  roffo(i) = thru(i) * (ap(i) / bp * (1.0 - exp(-bp *  &
                             xs(i))) + cp(i) * xs(i)) - equdep(i) *    &
                             xs(i)
                  roffo(i) = max(roffo(i), 0.0)
                  roff(i) = roff(i) + roffo(i)
! prtes = Precip Ready To Enter Soil.  This is the sum of the liquid
! precipitation not instantaneously running off due to convective rain.
! Calculate additional overland flow before liquid gets into the soil.
! Using linefit of basins, limiting to half of liquid to runoff.  The
! remaining liquid enters the soil.  SMP3, Eq. 8
                  prtes = thru(i) - roffo(i)
                  roffp(i) = min((2.07 * topostd(i) / 1843.47), 0.5) * &
                             prtes
                  roff(i) = roff(i) + roffp(i)
                  w(i,1) = w(i,1) + (prtes - roffp(i)) / (poros(i) *   &
                           zdepth(i,1))
! Modified surface runoff formulation from glen liston
                  excess = max(0., (w(i,1)-1.0))
                  w(i,1) = w(i,1) - excess
                  w(i,2) = w(i,2) + excess * zdepth(i,1) / zdepth(i,2)
                  excess = max(0., (w(i,2)-1.0))
                  w(i,2) = w(i,2) - excess
                  w(i,1) = w(i,1) + excess * zdepth(i,2) / zdepth(i,1)
               endif
            enddo
         endif

! Temperature change due to addition of precipitation
         do i = 1,istrip
            water(i,iveg) = capac(i,iveg) - snow(i,iveg)
            tsf(i) = (water(i,iveg) * cw) + (snow(i,iveg) * cs)
            if (tsf(i).gt.1.e-10) then
               tsd(i) = (((waterp(i,iveg) * cw) + (snowp(i,iveg) *     &
                        cs)) * (tsp(i) - tm(i)) / tsf(i)) + tm(i)
            else
               tsd(i) = tm(i)
            endif
         enddo

! Post-phase-change temp is 0C, which is modified in case more
! melts/freezes than is available to
         do i = 1,istrip
            freeze(i) = -((((water(i,iveg) * cw) + (snow(i,iveg) *     &
                        cs)) * (tsd(i) - tf)) + (cct(i) * (tsp(i) -    &
                        tf))) / snomel
            tsd(i) = tf
            snowp(i,iveg) = snow(i,iveg)
            snow(i,iveg) = snow(i,iveg) + freeze(i)
            water(i,iveg) = capac(i,iveg) - snow(i,iveg)
            if (water(i,iveg).lt.0.0) then
               tsd(i) = tsd(i) + water(i,iveg) * snomel / (cs *        &
                        capac(i,iveg) + cct(i))
               snow(i,iveg) = capac(i,iveg)
            endif
            if (snow(i,iveg).lt.0.0) then
               tsd(i) = tsd(i) - snow(i,iveg) * snomel / (cw *         &
                        capac(i,iveg) + cct(i))
               snow(i,iveg) = 0.0
            endif
            water(i,iveg) = capac(i,iveg) - snow(i,iveg)
            freeze(i) = snowp(i,iveg) - snow(i,iveg)
         enddo

         do i = 1,istrip
            if (iveg.eq.1) then
               tc(i) = tsd(i)
               ymelt(i,1) = max(0.0, freeze(i))
               yfrez(i,1) = -max(0.0, -freeze(i))
            else
               ymelt(i,2) = max(0.0, freeze(i))
               yfrez(i,2) = (-max(0.0, -freeze(i)))
! Modification of snow density with refreeze
               if ((sgfg(i).gt.0.5).and.(yfrez(i,2).ne.0.0)) then
                  sdens(i) = ((sdens(i) * (snow(i,2) + yfrez(i,2)))    &
                             + (-917.0 * yfrez(i,2))) / snow(i,2)
               endif
               tg(i) = tsd(i)
               tsn(i) = (tsn(i) * sgfg(i)) + (tsd(i) * (1.0 - sgfg(i)))
               tgsb(i) = (tsn(i) * sgfg(i)) + (tg(i) * (1.0 - sgfg(i)))
            endif
         enddo

         do i = 1,istrip
            p0(i) = thru(i)
         enddo
      enddo

! Calculation of canopy and ground heat capacities
      do i = 1,istrip
         water(i,1) = capac(i,1) - snow(i,1)
         cc(i) = (zlt(i,1) * clai) + (water(i,1) * cw) +(snow(i,1) * cs)
         water(i,2) = capac(i,2) - snow(i,2)
         spwet(i) = (water(i,2) * cw) + (snow(i,2) * cs)
         if (sgfg(i).gt.0.5) then
            cgsbi(i) = 1.0 / (sskin * cs)
         else
            cgsbi(i) = 1.0 / (spwet(i) + cg(i))
         endif
      enddo

      return
      end

!***********************************************************************
      subroutine sflxes(istrip,prin,ipoint,grav,refzh,pie,dtc3x,cpair,gasr,  &
           hlat,snomel,cs,epsfac,stefan,sskin,tf,delsig,ityp,lonco,    &
           latco,z2,dd,z0,phsat,poros,phsoil,psilow,bee,zdepth,zlt,    &
           vcover,rbc,rdc,topt,tll,tu,defac,ph1,ph2,rootd,chisl,bps,   &
           psb,ros,thm,zb,topostd,sws,swg,tm,ps,sh,spdm,ppl,ppc,       &
           starttotal,startsnow,cu,ustar,tg,w,capac,satcap,snow,tc,td, &
           tsn,tgsb,cc,cgsbi,sgfg,sdens,fluxd,fluxs,fluxmax,fluxmelt,  &
           sfall,ymelt,yfrez,thermk,radt,rst,ra,ewt,dtc,dtg,dtsn,      &
           dtgsb,dgdt,etmass,hflux,roff,hc,hg,eci,egi,ect,egt,egs,     &
           epot,gmt,gmq,sevap,eintmass,etranmass,ebaremass,etrancmass, &
           eintcmass,eintgmass,bhf,chf,ghf,shf)

      use LIS_logMod!,    only : LIS_logunit
      implicit none
!***********************************************************************
!***  surface flux parameterization
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! grav         (m sec-2) acceleration of gravity
! pie          (-)       3.14159....
! dtc3x        (sec)     timestep
! cpair     (J K-1 kg-1) heat capacity of dry air at constant pressure
! gasr      (J K-1 kg-1) gas constant of dry air
! hlat         (J kg-1)  latent heat of vaporization
! snomel       (J m-3)   latent heat of fusion X water density
! cs        (J m-3 K-1)  heat capacity of frozen water X water density
! epsfac       (-)       ratio of dry air to water vapor gas constant
! stefan    (W m-2 K-4)  Stefan-Boltzmann constant
! sskin        (m)       depth of diurnal snow layer
! tf           (K)       freezing temperature
! delsig       (-)       thickness (in sigma) of lowest GEOS GCM layer
! ityp         (1-13)    vegetation index
! lonco        (deg)     longitude (-180 west to 180 east)
! latco        (deg)     latitude (-90 south to 90 north)
! z2           (m)       canopy top height
! dd           (m)       displacement height
! z0           (m)       roughness height
! phsat        (m)       soil moisture potential at saturation
! poros        (0-1)     soil porosity
! phsoil(3)    (m)       soil moisture potential of the i-th soil layer
! psilow       (m)       lowest matric potential from evap., SMP3,Eq. 15
! bee          (-)       Clapp-Hornberger empirical constant
! zdepth(3)    (m)       exact independent depth of 3 soil layers
! zlt(cg)      (-)       leaf area index
! rbc   (sqr((s m-1)))   bulk canopy boundary layer resistance coeff.
! rdc          (-)       ground to canopy air space resistance coeff.
! topt         (K)       optimum temperature for rst calculation
! tll          (K)       low temperature for rst calculation
! tu           (K)       top temperature for rst calculation
! defac        (-)       dew factor for rst calculation
! ph1          (-)       stome slope factor
! ph2          (-)       point at which stomates close
! rootd(cg)    (m)       rooting depth for canopy [1] or ground [2]
! chisl     (W m-1 K-1)  soil conductivity
! bps          (-)       (ps/p0)**(gasr/cpair); p0 = 1000.0 mb
! psb          (mb)      depth of pbl
! ros          (kg m-3)  surface air density
! thm          (K)       reference-level potential temperature
! zb           (m)       elevation of the pbl top
! topostd      (m)       topography standard deviat. 1x1deg from GTOPO30
! sws          (W m-2)   sw radiation absorbed in bulk snow, SMP3, Eq.A5
! swg          (W m-2)   sw radiation absorbed in ground, SMP3, Eq. A6
! ps           (mb)      forcing surface pressure
! spdm         (m/sec)   forcing reference-level wind speed
! ppl          (mm/dt)   forcing large-scale precipitation
! ppc          (mm/dt)   forcing convective precipitation
! starttotal   (m)       total water at start of timestep
! startsnow    (m)       total snow at start of timestep
! cu           (-)       friction transfer coefficient
! ustar        (m/sec)   friction velocity
! w(3)         (0-1)     wetness of surface/root/recharge zone
! capac(cg)    (m)       liquid equiv. water&snow on canopy/ground
! satcap(cg)   (m)       interception capacity
! snow(cg)     (m)       liquid equiv. snow on canopy/ground
! td           (K)       deep soil temperature
! cc        (J m-2 K-1)  heat capacity of the canopy
! cgsbi    (J-1 m^2 K^1) inverse of heat capacity of blended ground/snow
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
! sdens        (kg m-3)  density of bulk snow layer
! sfall        (m/dt)    snowfall
! ymelt(cg)    (m/dt)    melted snow on canopy/ground from precip
! yfrez(cg)    (m/dt)    frozen water on canopy/ground from precip
! thermk       (0-1)     canopy emissivity
! radt(cg)     (W m-2)   net radiation
! dgdt         (W m-2)   heat flux from tg layer to td layer
! hc           (J m-2)   canopy sensible heat flux
! hg           (J m-2)   ground sensible heat flux
! eci          (J m-2)   canopy interception loss
! egi          (J m-2)   ground interception loss
! ect          (J m-2)   canopy transpiration
! egt          (J m-2)   ground transpiration
! egs          (J m-2)   ground evaporation from the soil
! epot         (J m-2)   potential evapotranspiration
! sevap(3)     (m/dt)    snow evapor. canopy(1)/bare soil(2)/grnd cvr(3)
! eintmass     (mm)      interception loss
! etranmass    (mm)      vegetation transpiration
! ebaremass    (mm)      evaporation from the soil
! bhf          (W m-2)   heat flux into the top snow layer from the air
! chf          (W m-2)   heat flux into the canopy from the air
! ghf          (W m-2)   heat flux into the ground from the air
! shf          (W m-2)   change in heat storage of bulk snow layer
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! vcover(cg)   (0-1)     fraction of vegetation cover [fpar/greenness]
! tm           (K)       forcing reference-level temperature
! sh           (kg kg-1) forcing reference-level watervapor mixing ratio
! tg           (K)       ground temperature
! tc           (K)       canopy temperature
! tsn          (K)       snow on ground temperature
! tgsb         (K)       blended ground/snow temperature
! fluxd        (W m-2)   heat flux from ground to snow, SMP1, Eq. 2c
! fluxs        (W m-2)   heat flux from diurnal snow to bulk snow
! fluxmax      (W m-2)   max. possible snow to ground flux, SMP1, Eq. 2e
! fluxmelt     (W m-2)   bottom melt flux, SMP1, Eq. 4b
! rst(cg)      (sec m-1) stomatal resistance
! ra           (sec m-1) aerodynamic resistance
! ewt          (0-1)     bare soil evap. weight of 3 soil layers
! dtc          (K)       change in tc over timestep
! dtg          (K)       change in tg over timestep
! dtsn         (K)       change in tsn over timestep
! dtgsb        (K)       change in tgsb over timestep
! etmass       (mm)      evaporation to the atmosphere (mass)
! hflux        (J m-2)   sensible heat flux to the atmosphere
! roff         (m/dt)    total runoff
! gmt          (K/dt)    temperature change from/for GEOS GCM lowest lvl
! gmq    ((kg kg-1)/dt)  humidity change from/for GEOS GCM lowest level
!-----------------------------------------------------------------------
      logical jstneu, prin
      integer istrip, i, icmax, icount, ipoint, ncount
      real grav, refzh, pie, dtc3x, cpair, gasr, hlat
      real snomel, cs, epsfac, stefan, sskin, tf
      real delsig, alph, summ, small, dtmdt, dshdt
      integer ityp(istrip), icheck(istrip)
      real lonco(istrip), latco(istrip)
      real z2(istrip), dd(istrip), z0(istrip)
      real phsat(istrip), poros(istrip)
      real phsoil(istrip,3), psilow(istrip)
      real bee(istrip), zdepth(istrip,3), zlt(istrip,2)
      real rbc(istrip), rdc(istrip)
      real topt(istrip), tll(istrip), tu(istrip)
      real defac(istrip), ph1(istrip), ph2(istrip)
      real rootd(istrip,2), chisl(istrip)
      real bps(istrip), psb(istrip), ros(istrip)
      real thm(istrip), zb(istrip)
      real topostd(istrip), fluxd(istrip), fluxs(istrip)
      real sws(istrip), swg(istrip), ps(istrip), spdm(istrip)
      real ppl(istrip), ppc(istrip)
      real starttotal(istrip), startsnow(istrip)
      real cu(istrip), ustar(istrip), w(istrip,3), capac(istrip,2)
      real satcap(istrip,2), snow(istrip,2), td(istrip)
      real cc(istrip), cgsbi(istrip), sgfg(istrip), sdens(istrip)
      real sfall(istrip), ymelt(istrip,2), yfrez(istrip,2)
      real thermk(istrip), radt(istrip,2), dgdt(istrip)
      real hc(istrip), hg(istrip), eci(istrip), egi(istrip)
      real ect(istrip), egt(istrip), egs(istrip), sevap(istrip,3)
      real eintmass(istrip), etranmass(istrip), ebaremass(istrip)
      real etrancmass(istrip), eintcmass(istrip), eintgmass(istrip)
      real bhf(istrip), chf(istrip), ghf(istrip), shf(istrip)
      real vcover(istrip,2), tm(istrip), sh(istrip), tg(istrip)
      real tc(istrip), tsn(istrip), tgsb(istrip), rst(istrip,2)
      real ewt(istrip,3), dtc(istrip), dtg(istrip), dtsn(istrip)
      real dtgsb(istrip), etmass(istrip), hflux(istrip), roff(istrip)
      real gmt(istrip,3), gmq(istrip,3)
      real ec(istrip), eg(istrip), epot(istrip), css(istrip)
      real stm(istrip,2), wc(istrip), wg(istrip), rsoil(istrip)
      real fc(istrip), fg(istrip), hrr(istrip), hr(istrip)
      real etc(istrip), etg(istrip)
      real ta(istrip), ea(istrip), tctm(istrip), tgtm(istrip)
      real ra(istrip), rb(istrip), rd(istrip), u2(istrip)
      real btc(istrip), btg(istrip), dtm(istrip), dsh(istrip)
      real fluxmax(istrip), fluxmelt(istrip), fp1(istrip), ft1(istrip)
      real deadtc(istrip), deadtg(istrip), deadsh(istrip)
      real psit(istrip), fac(istrip), y1(istrip), y2(istrip)
      real ecf(istrip), egf(istrip), dewc(istrip), dewg(istrip)
      real tcsav(istrip), tgsav(istrip), tmsav(istrip), shsav(istrip)
      real tsav(istrip), esav(istrip), idewco(istrip), rdsav(istrip,2)
      real tsnsav(istrip), wt(istrip), dzm(istrip)
      real speedm(istrip), cuni(istrip), ctni(istrip)
      data icmax /10/, small /1.0e-20/

      do i = 1,istrip
         tcsav(i) = tc(i)
         tgsav(i) = tg(i)
         tsnsav(i) = tsn(i)
         tmsav(i) = tm(i)
         shsav(i) = sh(i)
         rdsav(i,1) = radt(i,1)
         rdsav(i,2) = radt(i,2)
         stm(i,1) = rst(i,1)
         stm(i,2) = rst(i,2)
      enddo

      if (prin) then
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'sgfg: ',sgfg(ipoint)
         print *,'z2: ',z2(ipoint)
         print *,'capac_2: ',capac(ipoint,2)
         print *,'dd: ',dd(ipoint)
         print *,'z0: ',z0(ipoint)
         print *,'rbc: ',rbc(ipoint)
         print *,'rdc: ',rdc(ipoint)
         print *,'CALLING AIRMOD'
      endif

      call airmod(istrip,z2,dd,z0,rbc,rdc,snow,sgfg,sdens)

      if (prin) then
         print *,'dd: ',dd(ipoint)
         print *,'z0: ',z0(ipoint)
         print *,'rbc: ',rbc(ipoint)
         print *,'rdc: ',rdc(ipoint)
      endif

      do i = 1,istrip
         wc(i) = min(1.0, capac(i,1)/satcap(i,1))
         wg(i) = min(1.0, capac(i,2)/satcap(i,2))

! Weights for calculating rsoil based on soil moisture.
! SMP3, Eq. 11, 12, and 13
         alph = phsat(i) * grav / (461.5 * tgsb(i))
         if (alph.eq.0) then 
            write(LIS_logunit,*) "ALPH BAD: ", alph
         endif
         if (w(i,1).gt.0.0001) then
            ewt(i,1) = exp(alph * w(i,1) ** (-bee(i)))
         else
            ewt(i,1) = 0.0
         endif
         if (w(i,2).gt.0.0001) then
            ewt(i,2) = ((2.0 * zdepth(i,1)) / ((2.0 * zdepth(i,1)) +   &
                       zdepth(i,2))) * exp(alph * w(i,2) ** (-bee(i)))
            ewt(i,2) = (1.0 - ewt(i,1)) * ewt(i,2)
         else
            ewt(i,2) = 0.0
         endif
         if (w(i,3).gt.0.0001) then
            ewt(i,3) = ((2.0 * zdepth(i,1)) / ((2.0 * zdepth(i,1)) +   &
                       (2.0 * zdepth(i,2)) + zdepth(i,3))) * exp(alph  &
                       * w(i,3) ** (-bee(i)))
            ewt(i,3) = (1.0 - ewt(i,1) - ewt(i,2)) * ewt(i,3)
         else
            ewt(i,3) = 0.0
         endif
         summ = ewt(i,1) + ewt(i,2) + ewt(i,3)
         if ((summ.gt.0.0).and.(alph.ne.0.0).and.(bee(i).ne.0.0)) then
            wt(i) = ((1.0 / alph) * log(summ)) ** (-1.0 / bee(i))
         else
            wt(i) = small
         endif
         fac(i) = (0.99999 * sgfg(i)) + (min(wt(i),0.99) * (1.0 -      &
                  sgfg(i)))
         if (fac(i).lt.1e-3) then
            rsoil(i) = 101840.0
            hrr(i) = 0.0
         else
            rsoil(i) = 101840.0 * (1.0 - exp(0.0027 * log(fac(i))))
            psit(i) = phsat(i) * exp(-bee(i) * log(fac(i)))
            hrr(i) = exp(psit(i) * grav / (461.5 * tgsb(i)))
         endif
         hr(i) = hrr(i)

         if (tgsb(i).le.tf) then
            vcover(i,2) = 1.0
            wg(i) = min(1.0, capac(i,2)/0.004)
            rst(i,2) = rsoil(i)
            stm(i,2) = rsoil(i)
         endif
         fc(i) = 1.0
         fg(i) = 1.0
      enddo

      if (prin) then
         print *,' '
         print *,'wc: ',wc(ipoint)
         print *,'wg: ',wg(ipoint)
         print *,'ewt_1: ',ewt(ipoint,1)
         print *,'ewt_2: ',ewt(ipoint,2)
         print *,'ewt_3: ',ewt(ipoint,3)
         print *,'rsoil: ',rsoil(ipoint)
         print *,'hr: ',hr(ipoint)
      endif

! Start the iteration of time integration to avoid oscillation
      ncount = 0
 7000 continue
      ncount = ncount + 1

! Calculate saturation pressure over water or over ice
      do i = 1,istrip
         icheck(i) = 1
         etc(i) = exp(21.65605 - 5418.0 / tc(i))
         if (tc(i).lt.tf) etc(i) = etc(i) * ((tc(i) / tf) ** 2.66)
         etg(i) = exp(21.65605 - 5418.0 / tgsb(i))
         if (tgsb(i).lt.tf) etg(i) = etg(i) * ((tgsb(i) / tf) ** 2.66)
      enddo

      if (prin) then
         print *,'icheck: ',icheck(ipoint)
         print *,'etc: ',etc(ipoint)
         print *,'etg: ',etg(ipoint)
      endif

! First guesses for ta and ea
      if (ncount.eq.1) then
         do i = 1,istrip
            ta(i) = tc(i)
            ea(i) = sh(i) * ps(i) / (epsfac + sh(i))
         enddo
      endif

      if (prin) then
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'ta: ',ta(ipoint)
         print *,'ea: ',ea(ipoint)
         print *,'z2: ',z2(ipoint)
         print *,'dd: ',dd(ipoint)
         print *,'z0: ',z0(ipoint)
         print *,'thm: ',thm(ipoint)
         print *,'zb: ',zb(ipoint)
         print *,'topostd: ',topostd(ipoint)
         print *,'tm: ',tm(ipoint)
         print *,'spdm: ',spdm(ipoint)
         print *,'icheck: ',icheck(ipoint)
         print *,'CALLING VNTLAX'
      endif

! The first call to vntlax gets the neutral values of ustar and ventmf

      jstneu = .true.
      call vntlax(istrip,jstneu,grav,gasr,delsig,z2,dd,z0,thm,         &
           zb,topostd,bps,tm,spdm,ta,icheck,refzh,dzm,speedm,cuni,ctni,&
           u2,cu,ustar,ra)

      if (prin) then
         print *,'dzm: ',dzm(ipoint)
         print *,'speedm: ',speedm(ipoint)
         print *,'cuni: ',cuni(ipoint)
         print *,'ctni: ',ctni(ipoint)
         print *,'u2: ',u2(ipoint)

         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'ta: ',ta(ipoint)
         print *,'CALLING VNTLAX'
      endif

      jstneu = .false.
      call vntlax(istrip,jstneu,grav,gasr,delsig,z2,dd,z0,thm,         &
           zb,topostd,bps,tm,spdm,ta,icheck,refzh,dzm,speedm,cuni,ctni,&
           u2,cu,ustar,ra)

      if (prin) then
         print *,'cu: ',cu(ipoint)
         print *,'ustar: ',ustar(ipoint)
         print *,'ra: ',ra(ipoint)
      endif

      do i = 1,istrip
         tctm(i) = tc(i) - tm(i)
         tgtm(i) = tgsb(i) - tm(i)
      enddo

      if (prin) then
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'tctm: ',tctm(ipoint)
         print *,'tgtm: ',tgtm(ipoint)
         print *,'tgsb: ',tgsb(ipoint)
         print *,'z2: ',z2(ipoint)
         print *,'u2: ',u2(ipoint)
         print *,'zlt_1: ',zlt(ipoint,1)
         print *,'rbc: ',rbc(ipoint)
         print *,'rdc: ',rdc(ipoint)
         print *,'CALLING RB_RD'
      endif

      call rb_rd(istrip,tctm,tgtm,tgsb,z2,u2,zlt,rbc,rdc,rb,rd)

      if (prin) then
         print *,'rb: ',rb(ipoint)
         print *,'rd: ',rd(ipoint)
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
      endif

! Iterate for air temperature and ventilation mass flux.  This version
! assumes dew-free conditions to estimate ea for buoyancy term in
! vntmf or ra
      icount = 0
 2000 icount = icount + 1

      do i = 1,istrip
         if (icheck(i).eq.1) then
            tsav(i) = ta(i)
            esav(i) = ea(i)
         endif
      enddo

      if (prin) then
         print *,'ta: ',ta(ipoint)
         print *,'CALLING VNTLAX'
      endif

      call vntlax(istrip,jstneu,grav,gasr,delsig,z2,dd,z0,thm,         &
           zb,topostd,bps,tm,spdm,ta,icheck,refzh,dzm,speedm,cuni,ctni,&
           u2,cu,ustar,ra)

      if (prin) then
         print *,'cu: ',cu(ipoint)
         print *,'ustar: ',ustar(ipoint)
         print *,'ra: ',ra(ipoint)

         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'icheck: ',icheck(ipoint)
         print *,'vcover_2: ',vcover(ipoint,2)
         print *,'rsoil: ',rsoil(ipoint)
         print *,'rst_1: ',rst(ipoint,1)
         print *,'rst_2: ',rst(ipoint,2)
         print *,'hr: ',hr(ipoint)
         print *,'fc: ',fc(ipoint)
         print *,'fg: ',fg(ipoint)
         print *,'wc: ',wc(ipoint)
         print *,'wg: ',wg(ipoint)
         print *,'etc: ',etc(ipoint)
         print *,'etg: ',etg(ipoint)
         print *,'ps: ',ps(ipoint)
         print *,'sh: ',sh(ipoint)
         print *,'sgfg: ',sgfg(ipoint)
         print *,'CALLING CUT'
      endif

      call cut(istrip,epsfac,icheck,vcover,rsoil,rst,ra,rb,rd,         &
           hr,fc,fg,wc,wg,etc,etg,ps,sh,sgfg,ea)

      if (prin) then
         print *,'ea: ',ea(ipoint)

         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'icount: ',icount
         print *,'icheck: ',icheck(ipoint)
         print *,'ityp: ',ityp(ipoint)
         print *,'phsoil_1: ',phsoil(ipoint,1)
         print *,'phsoil_2: ',phsoil(ipoint,2)
         print *,'phsoil_3: ',phsoil(ipoint,3)
         print *,'psilow: ',psilow(ipoint)
         print *,'topt: ',topt(ipoint)
         print *,'tll: ',tll(ipoint)
         print *,'tu: ',tu(ipoint)
         print *,'defac: ',defac(ipoint)
         print *,'ph1: ',ph1(ipoint)
         print *,'ph2: ',ph2(ipoint)
         print *,'rootd_1: ',rootd(ipoint,1)
         print *,'rootd_2: ',rootd(ipoint,2)
         print *,'w_2: ',w(ipoint,2)
         print *,'tc: ',tc(ipoint)
         print *,'ta: ',ta(ipoint)
         print *,'ea: ',ea(ipoint)
         print *,'stm_1: ',stm(ipoint,1)
         print *,'rst_1: ',rst(ipoint,1)
         print *,'rst_2: ',rst(ipoint,2)
         print *,'CALLING STRES2'
      endif

      call stres2(istrip,icount,icheck,ityp,phsoil,psilow,             &
           zdepth,vcover,topt,tll,tu,defac,ph1,ph2,rootd,w,tc,ta,ea,   &
           fp1,ft1,stm,rst)

      if (prin) then
         print *,'rst_1: ',rst(ipoint,1)
         print *,'rst_2: ',rst(ipoint,2)

         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'icheck: ',icheck(ipoint)
         print *,'vcover_2: ',vcover(ipoint,2)
         print *,'rsoil: ',rsoil(ipoint)
         print *,'rst_1: ',rst(ipoint,1)
         print *,'rst_2: ',rst(ipoint,2)
         print *,'hr: ',hr(ipoint)
         print *,'fc: ',fc(ipoint)
         print *,'fg: ',fg(ipoint)
         print *,'wc: ',wc(ipoint)
         print *,'wg: ',wg(ipoint)
         print *,'etc: ',etc(ipoint)
         print *,'etg: ',etg(ipoint)
         print *,'ps: ',ps(ipoint)
         print *,'sh: ',sh(ipoint)
         print *,'sgfg: ',sgfg(ipoint)
         print *,'CALLING CUT'
      endif

      call cut(istrip,epsfac,icheck,vcover,rsoil,rst,ra,rb,rd,         &
           hr,fc,fg,wc,wg,etc,etg,ps,sh,sgfg,ea)

      if (prin) then
         print *,'ea: ',ea(ipoint)

         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
      endif

      do i = 1,istrip
         if (icheck(i).eq.1) then
            ta(i) = (tgsb(i) / rd(i) + tc(i) / rb(i) + thm(i) / ra(i)  &
                    * bps(i)) / (1.0 / rd(i) + 1.0 / rb(i) + 1.0 /ra(i))
         endif
      enddo

      if (prin) then
         print *,'ta: ',ta(ipoint)
      endif

      do i = 1,istrip
         if (icheck(i).eq.1) then
            y1(i) = abs(ta(i)-tsav(i))
            y2(i) = abs(ea(i)-esav(i))
            if (((y1(i).le.1.0e-2).and.(y2(i).le.5.0e-3)).or.          &
               (icount.gt.icmax)) then
               icheck(i) = 0
            endif
         endif
      enddo

      if (prin) then
         print *,'y1: ',y1(ipoint)
         print *,'y2: ',y2(ipoint)
         print *,'icheck: ',icheck(ipoint)
         print *,'Might goto 2000!'
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
      endif

      do i = 1,istrip
         if (icheck(i).eq.1) goto 2000
      enddo

      do i = 1,istrip
         fc(i) = 1.0
         fg(i) = 1.0
         idewco(i) = 0
         icheck(i) = 1
      enddo

! Calculate saturation pressure over water or over ice
      do i = 1,istrip
         tc(i) = tcsav(i)
         tg(i) = tgsav(i)
         tgsb(i) = (tsnsav(i) * sgfg(i)) + (tgsav(i) * (1.0 - sgfg(i)))
         tsn(i) = tsnsav(i)
         tm(i) = tmsav(i)
         sh(i) = shsav(i)
         radt(i,1) = rdsav(i,1)
         radt(i,2) = rdsav(i,2)
         etc(i) = exp(21.65605 - 5418.0 / tc(i))
         if (tc(i).lt.tf) etc(i) = etc(i) * ((tc(i) / tf) ** 2.66)
         etg(i) = exp(21.65605 - 5418.0 / tgsb(i))
         if (tgsb(i).lt.tf) etg(i) = etg(i) * ((tgsb(i) / tf) ** 2.66)
         btc(i) = exp(30.25353 - 5418.0 / tc(i)) / (tc(i) * tc(i))
         btg(i) = exp(30.25353 - 5418.0 / tgsb(i)) / (tgsb(i) * tgsb(i))
      enddo

      if (prin) then
         print *,'icheck: ',icheck(ipoint)
         print *,'etc: ',etc(ipoint)
         print *,'etg: ',etg(ipoint)
         print *,'btc: ',btc(ipoint)
         print *,'btg: ',btg(ipoint)
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
      endif

 3000 continue

      if (prin) then
         print *,'icheck: ',icheck(ipoint)
         print *,'vcover_2: ',vcover(ipoint,2)
         print *,'rsoil: ',rsoil(ipoint)
         print *,'rst_1: ',rst(ipoint,1)
         print *,'rst_2: ',rst(ipoint,2)
         print *,'hr: ',hr(ipoint)
         print *,'fc: ',fc(ipoint)
         print *,'fg: ',fg(ipoint)
         print *,'wc: ',wc(ipoint)
         print *,'wg: ',wg(ipoint)
         print *,'etc: ',etc(ipoint)
         print *,'etg: ',etg(ipoint)
         print *,'ps: ',ps(ipoint)
         print *,'sh: ',sh(ipoint)
         print *,'sgfg: ',sgfg(ipoint)
         print *,'CALLING CUT'
      endif

      call cut(istrip,epsfac,icheck,vcover,rsoil,rst,ra,rb,rd,         &
           hr,fc,fg,wc,wg,etc,etg,ps,sh,sgfg,ea)

      if (prin) then
         print *,'ea: ',ea(ipoint)

         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
      endif

      do i = 1,istrip
         if (icheck(i).eq.1) then
            ecf(i) = sign(1.0, etc(i)-ea(i))
            egf(i) = sign(1.0, etg(i)-ea(i))
            dewc(i) = fc(i) * 2.0 - 1.0
            dewg(i) = fg(i) * 2.0 - 1.0
            ecf(i) = ecf(i) * dewc(i)
            egf(i) = egf(i) * dewg(i)
         endif
      enddo

      if (prin) then
         print *,'icheck: ',icheck(ipoint)
         print *,'ecf: ',ecf(ipoint)
         print *,'egf: ',egf(ipoint)
         print *,'dewc: ',dewc(ipoint)
         print *,'dewg: ',dewg(ipoint)
      endif

      do i = 1,istrip
         if (((ecf(i).gt.0.0).and.(egf(i).gt.0.0)).or.(idewco(i).eq.3)) then
            icheck(i) = 0
         else
            idewco(i) = idewco(i) + 1
            if (idewco(i).eq.1) then
               fc(i) = 0.0
               fg(i) = 1.0
            elseif (idewco(i).eq.2) then
               fc(i) = 1.0
               fg(i) = 0.0
            elseif (idewco(i).eq.3) then
               fc(i) = 0.0
               fg(i) = 0.0
            endif
         endif
      enddo

      if (prin) then
         print *,'icheck: ',icheck(ipoint)
         print *,'idewco: ',idewco(ipoint)
         print *,'fc: ',fc(ipoint)
         print *,'fg: ',fg(ipoint)
         print *,'Might goto 3000!'
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
      endif

      do i = 1,istrip
         if (icheck(i).eq.1) goto 3000
      enddo

      if (prin) then
         print *,'tc: ',tc(ipoint)
         print *,'tg: ',tg(ipoint)
         print *,'tsn: ',tsn(ipoint)
         print *,'tgsb: ',tgsb(ipoint)
         print *,'tm: ',tm(ipoint)
         print *,'sh: ',sh(ipoint)
         print *,'CALLING TEMRES'
      endif

      call temres(istrip,prin,ipoint,ncount,grav,pie,dtc3x,            &
           cpair,hlat,cs,epsfac,stefan,sskin,tf,lonco,latco,poros,     &
           zdepth,vcover,chisl,rsoil,rst,ra,rb,rd,hr,hrr,fc,fg,btc,    &
           btg,etc,etg,bps,psb,ros,thm,wc,wg,fluxd,fluxs,fluxmax,      &
           fluxmelt,sws,swg,tm,ps,sh,ppl,ppc,starttotal,startsnow,     &
           ta,ea,tg,w,capac,satcap,snow,tc,td,tsn,tgsb,cc,cgsbi,css,   &
           sgfg,sdens,sfall,ymelt,yfrez,thermk,radt,roff,hc,hg,ec,eg,  &
           eci,egi,ect,egt,egs,epot,gmt,gmq,dtc,dtg,dtsn,dtgsb,dtm,    &
           dsh,deadtc,deadtg,deadsh)

      if (prin) then
         print *,'hc: ',hc(ipoint)
         print *,'hg: ',hg(ipoint)
         print *,'ec: ',ec(ipoint)
         print *,'eg: ',eg(ipoint)
         print *,'eci: ',eci(ipoint)
         print *,'ect: ',ect(ipoint)
         print *,'egi: ',egi(ipoint)
         print *,'egt: ',egt(ipoint)
         print *,'egs: ',egs(ipoint)
         print *,'epot: ',epot(ipoint)
         print *,'dtc: ',dtc(ipoint)
         print *,'dtg: ',dtg(ipoint)
         print *,'dtsn: ',dtsn(ipoint)
         print *,'dtgsb: ',dtgsb(ipoint)
         print *,'dtm: ',dtm(ipoint)
         print *,'dsh: ',dsh(ipoint)
      endif

      if (ncount.le.1) then
         do i = 1,istrip
            tc(i) = tc(i) + dtc(i)
            tg(i) = tg(i) + dtg(i)
            tgsb(i) = tgsb(i) + dtgsb(i)
            tsn(i) = tsn(i) + dtgsb(i)
            tm(i) = tm(i) + dtm(i)
            sh(i) = sh(i) + dsh(i)
         enddo
         goto 7000
      endif

      if (prin) then
         print *,' '
         print *,'-----------------------------------------------------'
         print *,' '
         print *,'tc: ',tc(ipoint)
         print *,'tg: ',tg(ipoint)
         print *,'tsn: ',tsn(ipoint)
         print *,'tgsb: ',tgsb(ipoint)
         print *,'tm: ',tm(ipoint)
         print *,'sh: ',sh(ipoint)
         print *,'CALLING UPDATE'
      endif

      call update(istrip,prin,ipoint,pie,dtc3x,hlat,snomel,tf,         &
           ra,rb,rd,btc,btg,etc,etg,bps,fluxd,tm,ea,tg,capac,snow,tc,  &
           td,tgsb,cc,cgsbi,css,sgfg,hc,hg,eg,eci,egi,ect,egt,egs,dtc, &
           dtg,dtgsb,dtm,dsh,dgdt,deadtc,deadtg,deadsh,sevap,bhf,chf,  &
           ghf,shf,eintmass,etranmass,ebaremass,etrancmass,eintcmass,  &
           eintgmass,etmass,hflux)

      do i = 1,istrip
         fac(i) = grav / (100.0 * psb(i) * dtc3x)
         dtmdt = (gmt(i,3) + hflux(i) * fac(i) / (cpair * bps(i))) /   &
                  gmt(i,2)
         dshdt = (gmq(i,3) + etmass(i) * fac(i)) / gmq(i,2)
         dtm(i) = dtmdt * dtc3x
         dsh(i) = dshdt * dtc3x
         gmt(i,3) = dtmdt
         gmq(i,3) = dshdt
         tm(i) = tm(i) + dtm(i)
         sh(i) = sh(i) + dsh(i)
      enddo

      return
      end

!***********************************************************************
      subroutine airmod(istrip,z2,dd,z0,rbc,rdc,snow,sgfg,sdens)

      implicit none
!***********************************************************************
!***  alteration of aerodynamic transfer properties in case of snow
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! z2           (m)       canopy top height
! snow(cg)     (m)       liquid equiv. snow on canopy/ground
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
! sdens        (kg m-3)  density of bulk snow layer
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! dd           (m)       displacement height
! z0           (m)       roughness height
! rbc   (sqr((s m-1)))   bulk canopy boundary layer resistance coeff.
! rdc          (-)       ground to canopy air space resistance coeff.
!-----------------------------------------------------------------------
      integer istrip, i
      real z2(istrip), dd(istrip), z0(istrip), rbc(istrip), rdc(istrip)
      real snow(istrip,2), sgfg(istrip), sdens(istrip)
      real sdep(istrip)

      do i = 1,istrip
         if (sgfg(i).gt.0.5) then
            sdep(i) = snow(i,2) * 1000.0 / sdens(i)
            sdep(i) = min(sdep(i), z2(i)*0.75)
            dd(i) = z2(i) - (z2(i) - dd(i)) / z2(i) * (z2(i) - sdep(i))
            z0(i) = z0(i) / (z2(i) - dd(i)) * (z2(i) - dd(i))
            rbc(i) = rbc(i) * z2(i) / (z2(i) - sdep(i))
            rdc(i) = rdc(i) * (z2(i) - sdep(i)) / z2(i)
         endif
      enddo

      return
      end

!***********************************************************************
      subroutine vntlax(istrip,jstneu,grav,gasr,delsig,z2,dd,z0,thm,   &
           zb,topostd,bps,tm,spdm,ta,icheck,refzh,dzm,speedm,cuni,ctni,&
           u2,cu,ustar,ra)

      implicit none
!***********************************************************************
!***  ventilation mass flux, based on Deardorff, 1972
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! jstneu       (t/f)     logical: true  = returns neutral values
!                                 false = returns non-neutral values
! grav         (m sec-2) acceleration of gravity
! gasr      (J K-1 kg-1) gas constant of dry air
! delsig       (-)       thickness (in sigma) of lowest GEOS GCM layer
! z2           (m)       canopy top height
! dd           (m)       displacement height
! z0           (m)       roughness height
! thm          (K)       reference-level potential temperature
! zb           (m)       elevation of the pbl top
! topostd      (m)       topography standard deviat. 1x1deg from GTOPO30
! bps          (-)       (ps/p0)**(gasr/cpair); p0 = 1000.0 mb
! tm           (K)       forcing reference-level temperature
! spdm         (m/sec)   forcing reference-level wind speed
! ta           (K)       canopy air space temperature
! icheck       (0/1)     0 = air space temp/moist converged; 1 = not yet
! refzh        (m)       reference forcing height for T,q,U,V
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! dzm          (m)       reference-level point height
! speedm       (m/sec)   wind speed with added velocity scale
! cuni         (-)       inverse of neutral friction coefficient
! ctni         (-)       inverse of neutral heat coefficient
! u2           (m/sec)   wind speed at canopy top height (z2)
! ustar        (m/sec)   friction velocity
! ra           (sec m-1) aerodynamic resistance bet. reference & canopy
!-----------------------------------------------------------------------
      logical jstneu
      integer istrip, i
      real grav, refzh, gasr, delsig, vkrmni, g2, zl, xct1, xct2
      real xctu1, xctu2, bh, zdrg, grib, grzl, grz2
      real fvv, ftt, rzl, rz2
      integer icheck(istrip)
      real z2(istrip), dd(istrip), z0(istrip), thm(istrip), zb(istrip)
      real topostd(istrip), tm(istrip), spdm(istrip), ta(istrip)
      real ctni(istrip), cuni(istrip), u2(istrip)
      real cu(istrip), ustar(istrip), ra(istrip)
      real dzm(istrip), speedm(istrip), thvgm(istrip)
      real ct(istrip), cti(istrip), cui(istrip)
      real rib(istrip), ustarn(istrip), bps(istrip)
      data vkrmni /2.5/, g2 /0.75/

! cu and ct are the friction and heat transfer coefficients
! cun and ctn are the neutral friction and heat transfer coefficients
      if (jstneu) then
         do i = 1,istrip
            zl = z2(i) + 11.785 * z0(i)
! Assuming dzm based on the prescribed reference heights
! for T,q,U,V from hyssib_struc(n)%zh,zm
            dzm(i) = max((0.5 * gasr / grav * delsig * tm(i)),refzh)
            cuni(i) = log((dzm(i) - dd(i)) / z0(i)) * vkrmni
            if (zl.lt.dzm(i)) then
               xct1 = log((dzm(i) - dd(i)) / (zl - dd(i)))
               xct2 = log((zl - dd(i)) / z0(i))
               xctu1 = xct1
               xctu2 = log((zl - dd(i)) / (z2(i) - dd(i)))
               ctni(i) = (xct1 + g2 * xct2) * vkrmni
            else
               xct2 = log((dzm(i) - dd(i)) / z0(i))
               xctu1 = 0.0
               xctu2 = log((dzm(i) - dd(i)) / (z2(i) - dd(i)))
               ctni(i) = g2 * xct2 * vkrmni
            endif

! After Mahrt and Sun, 1995 - we effectively add to the flux by adding
! the velocity scale to the wind speed, following Eq. 11 and 12.  bh is
! a constant and we assume the pbl depth.  The wind speed is only
! increased for unstable conditions.  SMP3, Eq. 17 & 18.  2/25/99 dmm
            speedm(i) = spdm(i)
            thvgm(i) = (ta(i) / bps(i)) - thm(i)
            if (thvgm(i).gt.0.0) then
               bh = 0.00025
               speedm(i) = speedm(i) + (2.0 * bh * sqrt(grav / (ta(i)  &
                           / bps(i)) * zb(i) * thvgm(i)) * ctni(i) *   &
                           cuni(i))
            endif
            speedm(i) = max(2.0, speedm(i))
            ustarn(i) = speedm(i) / cuni(i)
! Modify cuni and ctni to account for orography.  11/7/91 gkw
            zdrg = 0.006 * topostd(i) / 2000.
            cuni(i) = cuni(i) / sqrt(1. + zdrg * cuni(i) ** 2)
            ctni(i) = ctni(i) / sqrt(1. + zdrg * ctni(i) ** 2)
! Neutral values
            u2(i) = speedm(i) - ustarn(i) * vkrmni *(xctu1 + g2 * xctu2)
         enddo

         return

      endif

! Stability branch based on bulk richardson number
      do i = 1,istrip
         if (icheck(i).eq.1) then
            zl = z2(i) + 11.785 * z0(i)
            thvgm(i) = (ta(i) / bps(i)) - thm(i)
            rib(i) = -thvgm(i) * grav * (dzm(i) - dd(i)) / (thm(i)     &
                     * (speedm(i) - u2(i)) ** 2)
            rib(i) = max(-10., rib(i))
            rib(i) = min(0.165, rib(i))
            if (rib(i).lt.0.0) then
               grib = -rib(i)
               grzl = -rib(i) * (zl - dd(i)) / (dzm(i) - dd(i))
               grz2 = -rib(i) * z0(i) / (dzm(i) - dd(i))
               fvv = 0.315 * grib
               if (zl.lt.dzm(i)) then
                  ftt = 0.904 * (grib + (g2 - 1.) * grzl - g2 * grz2)
               else
                  ftt = 0.904 * g2 * (grib - grz2)
               endif
               cui(i) = cuni(i) - fvv
               cti(i) = ctni(i) - ftt
            else
               rzl = rib(i) / (dzm(i) - dd(i)) * (zl - dd(i))
               rz2 = rib(i) / (dzm(i) - dd(i)) * z0(i)
               fvv = 66.85 * rib(i)
               if (zl.lt.dzm(i)) then
                  ftt = fvv + (g2 - 1) * 66.85 * rzl - g2 * 66.85 * rz2
               else
                  ftt = g2 * (fvv - 66.85 * rz2)
               endif
               cui(i) = cuni(i) + fvv
               cti(i) = ctni(i) + ftt
            endif
            cu(i) = 1.0 / cui(i)
            ct(i) = 1.0 / cti(i)
! Surface friction velocity
            ustar(i) = speedm(i) * cu(i)
            ra(i) = cti(i) / ustar(i)
            ra(i) = max(ra(i), 0.8)
         endif
      enddo

      return
      end

!***********************************************************************
      subroutine rb_rd(istrip,tctm,tgtm,tgsb,z2,u2,zlt,rbc,rdc,rb,rd)

      implicit none
!***********************************************************************
!***  rb and rd as functions of u2 and temperatures
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! tctm         (K)       difference between tc and tm
! tgtm         (K)       difference between tg and tm
! tgsb         (K)       blended ground/snow temperature
! z2           (m)       canopy top height
! u2           (m/sec)   wind speed at canopy top height (z2)
! zlt(cg)      (-)       leaf area index
! rbc   (sqr((s m-1)))   bulk canopy boundary layer resistance coeff.
! rdc          (-)       ground to canopy air space resistance coeff.
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! rb           (sec m-1) bulk canopy boundary layer resistance
! rd           (sec m-1) ground to canopy air space resistance
!-----------------------------------------------------------------------
      integer istrip, i
      real factg
      real tctm(istrip), tgtm(istrip), tgsb(istrip)
      real z2(istrip), u2(istrip), zlt(istrip,2)
      real rbc(istrip), rdc(istrip), rb(istrip), rd(istrip)
      real temdif(istrip), fih(istrip)

!      factc = 2.3761e-3
      factg = 88.29
      do i = 1,istrip
         if (tctm(i).gt.0.0) then
            temdif(i) = tctm(i) + 0.1
         else
            temdif(i) = 0.1
         endif
         rb(i) = 1.0 / (sqrt(u2(i)) / rbc(i) + zlt(i,1) * 0.004)
         if (tgtm(i).gt.0) then
            temdif(i) = tgtm(i) + 0.1
         else
            temdif(i) = 0.1
         endif
         fih(i) = sqrt(1.0 + factg * temdif(i) * z2(i) / (tgsb(i) *    &
                  u2(i) * u2(i)))
         rd(i) = rdc(i) / (u2(i) * fih(i))
      enddo

      return
      end

!***********************************************************************
      subroutine cut(istrip,epsfac,icheck,vcover,rsoil,rst,ra,rb,rd,   &
           hr,fc,fg,wc,wg,etc,etg,ps,sh,sgfg,ea)

      implicit none
!***********************************************************************
!***  calculate canopy air space vapor pressure
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! epsfac       (-)       ratio of dry air to water vapor gas constant
! icheck       (0/1)     0 = air space temp/moist converged; 1 = not yet
! vcover(cg)   (0-1)     fraction of vegetation cover [fpar/greenness]
! rsoil        (sec m-1) bare soil surface resistance
! rst(cg)      (sec m-1) stomatal resistance
! ra           (sec m-1) aerodynamic resistance bet. reference & canopy
! rb           (sec m-1) bulk canopy boundary layer resistance
! rd           (sec m-1) ground to canopy air space resistance
! hr           (-)       relative humidity within upper surface layer
! fc           (0/1)     flag for dew on canopy
! fg           (0/1)     flag for dew on ground
! wc           (0-1)     wetness fraction for canopy store
! wg           (0-1)     wetness fraction for ground store
! etc          (mb)      saturation vapor pressure at canopy
! etg          (mb)      saturation vapor pressure at ground
! ps           (mb)      forcing surface pressure
! sh           (kg kg-1) forcing reference-level watervapor mixing ratio
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! ea           (mb)      canopy air space vapor pressure
!-----------------------------------------------------------------------
      integer istrip, i
      real epsfac
      integer icheck(istrip)
      real vcover(istrip,2), rsoil(istrip), rst(istrip,2)
      real ra(istrip), rb(istrip), rd(istrip)
      real hr(istrip), fc(istrip), fg(istrip), wc(istrip), wg(istrip)
      real etc(istrip), etg(istrip), ps(istrip), sh(istrip)
      real sgfg(istrip), ea(istrip)
      real rc(istrip), rg(istrip), em(istrip)
      real coc(istrip), rsurf(istrip), cog1(istrip), cog2(istrip)
      real d2(istrip), top(istrip), xnum(istrip), tem(istrip)

      do i = 1,istrip
         if (icheck(i).eq.1) then
            rc(i) = rst(i,1) * fc(i) + rb(i) + rb(i) * fc(i)
            coc(i) = (1.0 - wc(i)) / rc(i) + wc(i) / (2.0 * rb(i))

! Over snow, rsurf = 0 and rg = rg + 500.  1/24/98 dmm
            rg(i) = rst(i,2) * fg(i) + (sgfg(i) * 500.0)
            rsurf(i) = rsoil(i) * fg(i) * (1.0 - sgfg(i))

            tem(i) = vcover(i,2) * (1.0 - wg(i)) / (rg(i) + rd(i))
            cog2(i) = tem(i) + (1.0 - vcover(i,2)) / (rsurf(i) +       &
                      rd(i)) + vcover(i,2) / (rsurf(i) + rd(i) + 44.)
            cog1(i) = (cog2(i) - tem(i)) * hr(i) + tem(i)

            xnum(i) = wg(i) / rd(i) * vcover(i,2)
            cog1(i) = cog1(i) + xnum(i)
            cog2(i) = cog2(i) + xnum(i)

            d2(i) = (1.0 / ra(i)) + coc(i) + cog2(i)
            em(i) = sh(i) * ps(i) / (epsfac + sh(i))
            top(i) = (coc(i) * etc(i)) + (em(i) / ra(i)) + (cog1(i) *  &
                     etg(i))
            ea(i) = top(i) / d2(i)
         endif
      enddo

      return
      end

!***********************************************************************
      subroutine stres2(istrip,icount,icheck,ityp,phsoil,psilow,       &
           zdepth,vcover,topt,tll,tu,defac,ph1,ph2,rootd,w,tc,ta,ea,   &
           fp1,ft1,stm,rst)

      implicit none
!***********************************************************************
!***  calculation of adjustment to light dependent stomatal resistance
!     by temperature, humidity and stress factors (simplified)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! icheck       (0/1)     0 = air space temp/moist converged; 1 = not yet
! ityp         (1-13)    vegetation index
! phsoil(3)    (m)       soil moisture potential of the i-th soil layer
! psilow       (m)       lowest matric potential from evap., SMP3,Eq. 15
! zdepth(3)    (m)       exact independent depth of 3 soil layers
! vcover(cg)   (0-1)     fraction of vegetation cover [fpar/greenness]
! topt         (K)       optimum temperature for rst calculation
! tll          (K)       low temperature for rst calculation
! tu           (K)       top temperature for rst calculation
! defac        (-)       dew factor for rst calculation
! ph1          (-)       stome slope factor
! ph2          (-)       point at which stomates close
! rootd(cg)    (m)       rooting depth for canopy [1] or ground [2]
! w(3)         (0-1)     wetness of surface/root/recharge zone
! tc           (K)       canopy temperature
! ta           (K)       canopy air space temperature
! ea           (mb)      canopy air space vapor pressure
! stm(cg)      (sec m-1) stomatal resistance (saved)
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! rst(cg)      (sec m-1) stomatal resistance
!-----------------------------------------------------------------------
      integer istrip, icount, i, j
      real xrot, xdr, dep(3)
      integer icheck(istrip), ityp(istrip)
      real phsoil(istrip,3), psilow(istrip), zdepth(istrip,3)
      real vcover(istrip,2), topt(istrip), tll(istrip), tu(istrip)
      real defac(istrip), ph1(istrip), ph2(istrip), rootd(istrip,2)
      real w(istrip,3), tc(istrip), ta(istrip), ea(istrip)
      real fp1(istrip), ft1(istrip)
      real stm(istrip,2), rst(istrip,2)
      real tv(istrip), d2(istrip), ft(istrip), fp(istrip)
      real drop(istrip), fd(istrip), ftpd(istrip)

! Humidity, temperature and transpiration factors
      if (icount.eq.1) then
         do i = 1,istrip
            if (icheck(i).eq.1) then
               if (ityp(i).eq.13.0) goto 777
               tv(i) = tc(i)
               tv(i) = min((tu(i) - 0.1), tv(i))
               tv(i) = max((tll(i) + 0.1), tv(i))

               d2(i) = (tu(i) - topt(i)) / (topt(i) - tll(i))
               ft(i) = (tv(i) - tll(i)) / (topt(i) - tll(i)) *         &
                       exp(d2(i) * log((tu(i) - tv(i)) / (tu(i) -      &
                       topt(i))))
               ft(i) = min(ft(i), 1.)
               ft(i) = max(ft(i), 1.e-5)
               ft1(i) = ft(i)
! Simplified calculation of leaf water potential factor, fp
               xrot = rootd(i,1)
               do j = 1,3
                  dep(j) = 0.0
               enddo
               do j = 1,3
                  dep(j) = min(zdepth(i,j), xrot)
                  xrot = xrot - zdepth(i,j)
                  if (xrot.le.0.0) goto 240
               enddo
 240           continue
               xdr = max(phsoil(i,1), phsoil(i,2), psilow(i))
               if (dep(3).gt.0.0) xdr = max(xdr, phsoil(i,3))
               xdr = -xdr
               xdr = log(xdr)
! Coefficient values from Xue et al, 1991, Table 2 (labels reversed)
               fp(i) = 1. - exp(-ph1(i) * (ph2(i) - xdr))
               if ((w(i,2).gt.0.15).and.(fp(i).lt.0.05)) fp(i) = 0.05
               fp(i) = min(fp(i), 1.)
               fp(i) = max(fp(i), 1.e-10)
               fp1(i) = fp(i)
            endif
 777        continue
         enddo
      endif

      do i = 1,istrip
         if (icheck(i).eq.1) then
            drop(i) = exp(21.65605 - 5418.0 / ta(i)) - ea(i)
            fd(i) = max(1.0e-5, (1.0 - drop(i) * defac(i)))
            fd(i) = min(fd(i), 1.)
! Stomatal control applied 9/12/91 to avoid "stomatal suicide"
            fd(i) = max(fd(i), 0.30)
         endif
      enddo

      do i = 1,istrip
         if (icheck(i).eq.1) then
            if (ityp(i).eq.13.0) then
               rst(i,1) = 1.0e10
            else
               ftpd(i) = fd(i) * ft1(i) * fp1(i)
               rst(i,1) = stm(i,1) / (ftpd(i) * vcover(i,1))
               rst(i,1) = min(rst(i,1), 1.0e10)
            endif
            rst(i,2) = 1.e10
         endif
      enddo

      return
      end

!***********************************************************************
      subroutine temres(istrip,prin,ipoint,ncount,grav,pie,dtc3x,      &
           cpair,hlat,cs,epsfac,stefan,sskin,tf,lonco,latco,poros,     &
           zdepth,vcover,chisl,rsoil,rst,ra,rb,rd,hr,hrr,fc,fg,btc,    &
           btg,etc,etg,bps,psb,ros,thm,wc,wg,fluxd,fluxs,fluxmax,      &
           fluxmelt,sws,swg,tm,ps,sh,ppl,ppc,starttotal,startsnow,     &
           ta,ea,tg,w,capac,satcap,snow,tc,td,tsn,tgsb,cc,cgsbi,css,   &
           sgfg,sdens,sfall,ymelt,yfrez,thermk,radt,roff,hc,hg,ec,eg,  &
           eci,egi,ect,egt,egs,epot,gmt,gmq,dtc,dtg,dtsn,dtgsb,dtm,    &
           dsh,deadtc,deadtg,deadsh)

      implicit none
!***********************************************************************
!***  temperature tendency equations with interception loss
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! ncount       (-)       iteration number of call to this routine
! grav         (m sec-2) acceleration of gravity
! pie          (-)       3.14159....
! dtc3x        (sec)     timestep
! cpair     (J K-1 kg-1) heat capacity of dry air at constant pressure
! hlat         (J kg-1)  latent heat of vaporization
! cs        (J m-3 K-1)  heat capacity of frozen water X water density
! epsfac       (-)       ratio of dry air to water vapor gas constant
! stefan    (W m-2 K-4)  Stefan-Boltzmann constant
! sskin        (m)       depth of diurnal snow layer
! tf           (K)       freezing temperature
! lonco        (deg)     longitude (-180 west to 180 east)
! latco        (deg)     latitude (-90 south to 90 north)
! poros        (0-1)     soil porosity
! zdepth(3)    (m)       exact independent depth of 3 soil layers
! chisl     (W m-1 K-1)  soil conductivity
! rsoil        (sec m-1) bare soil surface resistance
! ra           (sec m-1) aerodynamic resistance bet. reference & canopy
! rb           (sec m-1) bulk canopy boundary layer resistance
! rd           (sec m-1) ground to canopy air space resistance
! hrr          (-)       temp rel. humid. within upper surface layer
! fc           (0/1)     flag for dew on canopy
! fg           (0/1)     flag for dew on ground
! btc          (mb K-1)  first partial of satur. canopy vapor pressure
! etc          (mb)      saturation vapor pressure at canopy
! etg          (mb)      saturation vapor pressure at ground
! bps          (-)       (ps/p0)**(gasr/cpair); p0 = 1000.0 mb
! psb          (mb)      depth of pbl
! ros          (kg m-3)  surface air density
! thm          (K)       reference-level potential temperature
! sws          (W m-2)   sw radiation absorbed in bulk snow, SMP3, Eq.A5
! swg          (W m-2)   sw radiation absorbed in ground, SMP3, Eq. A6
! tm           (K)       forcing reference-level temperature
! ps           (mb)      forcing surface pressure
! sh           (kg kg-1) forcing reference-level watervapor mixing ratio
! ppl          (mm/dt)   forcing large-scale precipitation
! ppc          (mm/dt)   forcing convective precipitation
! starttotal   (m)       total water at start of timestep
! startsnow    (m)       total snow at start of timestep
! tg           (K)       ground temperature
! w(3)         (0-1)     wetness of surface/root/recharge zone
! capac(cg)    (m)       liquid equiv. water&snow on canopy/ground
! satcap(cg)   (m)       interception capacity
! snow(cg)     (m)       liquid equiv. snow on canopy/ground
! tc           (K)       canopy temperature
! td           (K)       deep soil temperature
! tsn          (K)       snow on ground temperature
! tgsb         (K)       blended ground/snow temperature
! cc        (J m-2 K-1)  heat capacity of the canopy
! cgsbi    (J-1 m^2 K^1) inverse of heat capacity of blended ground/snow
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
! sfall        (m/dt)    snowfall
! ymelt(cg)    (m/dt)    melted snow on canopy/ground from precip
! yfrez(cg)    (m/dt)    frozen water on canopy/ground from precip
! thermk       (0-1)     canopy emissivity
! roff         (m/dt)    total runoff
! gmt          (K/dt)    temperature change from/for GEOS GCM lowest lvl
! gmq    ((kg kg-1)/dt)  humidity change from/for GEOS GCM lowest level
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! vcover(cg)   (0-1)     fraction of vegetation cover [fpar/greenness]
! rst(cg)      (sec m-1) stomatal resistance
! hr           (-)       relative humidity within upper surface layer
! btg          (mb K-1)  first partial of satur. ground vapor pressure
! wc           (0-1)     wetness fraction for canopy store
! wg           (0-1)     wetness fraction for ground store
! fluxd        (W m-2)   heat flux from ground to snow, SMP1, Eq. 2c
! fluxs        (W m-2)   heat flux from diurnal snow to bulk snow
! fluxmax      (W m-2)   max. possible snow to ground flux, SMP1, Eq. 2e
! fluxmelt     (W m-2)   bottom melt flux, SMP1, Eq. 4b
! ta           (K)       canopy air space temperature
! ea           (mb)      canopy air space vapor pressure
! css       (J m-2 K-1)  heat capacity of bulk snow X water density
! sdens        (kg m-3)  density of bulk snow layer
! radt(cg)     (W m-2)   net radiation
! hc           (J m-2)   canopy sensible heat flux
! hg           (J m-2)   ground sensible heat flux
! ec           (J m-2)   canopy evaporation (eci + ect)
! eg           (J m-2)   ground evaporation (egi + egt + egs)
! eci          (J m-2)   canopy interception loss
! egi          (J m-2)   ground interception loss
! ect          (J m-2)   canopy transpiration
! egt          (J m-2)   ground transpiration
! egs          (J m-2)   ground evaporation from the soil
! epot         (J m-2)   potential evapotranspiration
! dtc          (K)       change in tc over timestep
! dtg          (K)       change in tg over timestep
! dtsn         (K)       change in tsn over timestep
! dtgsb        (K)       change in tgsb over timestep
! dtm          (K)       change in tm over timestep
! dsh          (kg kg-1) change in sh over timestep
! deadtc       (mb K-1)  change in ea w/r/t tc
! deadtg       (mb K-1)  change in ea w/r/t tg
! deadsh  (mb/(kg kg-1)) change in ea w/r/t sh
!-----------------------------------------------------------------------
!     in subroutine parameters
!-----------------------------------------------------------------------
! amult        (0/1)     flag to shut off some fluxes if snow melts
! snowcon   (W m-1 K-1)  snow conductivity
! condk     (W m-2 K-1)  ground/snow interface conductivity, SMP3, Eq.A7
! rsurf        (sec m-1) bare soil surface ground resistance
! rc           (sec m-1) bulk stomatal resistance of canopy
! rg           (sec m-1) bulk stomatal resistance of ground vegetation
!-----------------------------------------------------------------------
      logical prin
      integer istrip, ipoint, ncount, i, j
      logical stebt(istrip)
      real grav, pie, dtc3x, cpair, hlat
      real cs, epsfac, stefan, sskin, tf
      real dtc3xi, capi, timcon, tim, fak, fah, stb4, stb8, hlat3
      real lonco(istrip), latco(istrip), poros(istrip)
      real zdepth(istrip,3), chisl(istrip), rsoil(istrip)
      real ra(istrip), rb(istrip), rd(istrip)
      real hrr(istrip), fc(istrip), fg(istrip)
      real btc(istrip), etc(istrip), etg(istrip)
      real bps(istrip), psb(istrip), ros(istrip)
      real thm(istrip), sws(istrip), swg(istrip)
      real tm(istrip), ps(istrip), sh(istrip)
      real ppl(istrip), ppc(istrip)
      real starttotal(istrip), startsnow(istrip)
      real tg(istrip), w(istrip,3), capac(istrip,2), satcap(istrip,2)
      real snow(istrip,2), tc(istrip), td(istrip), tsn(istrip)
      real tgsb(istrip), cc(istrip), cgsbi(istrip), sgfg(istrip)
      real sfall(istrip), ymelt(istrip,2), yfrez(istrip,2)
      real thermk(istrip), roff(istrip)
      real gmt(istrip,3), gmq(istrip,3)
      real vcover(istrip,2), rst(istrip,2), hr(istrip), btg(istrip)
      real wc(istrip), wg(istrip), fluxd(istrip), fluxs(istrip)
      real fluxmax(istrip), fluxmelt(istrip)
      real ta(istrip), ea(istrip), css(istrip), sdens(istrip)
      real radt(istrip,2), hc(istrip), hg(istrip), ec(istrip)
      real eg(istrip), eci(istrip), egi(istrip), ect(istrip)
      real egt(istrip), egs(istrip), epot(istrip)
      real dtc(istrip), dtg(istrip), dtsn(istrip)
      real dtgsb(istrip), dtm(istrip), dsh(istrip)
      real deadtc(istrip), deadtg(istrip), deadsh(istrip)
      real amult(istrip), snowcon(istrip), condk(istrip)
      real rsurf(istrip), rc(istrip), rg(istrip), es(istrip)
      real pblsib(istrip,4,5), ai(istrip), coc(istrip), rcp(istrip)
      real cog1(istrip), cog2(istrip), d1(istrip), d2(istrip)
      real d1i(istrip), top(istrip), ak(istrip), ah(istrip)
      real cci(istrip), ecpot(istrip), egpot(istrip)
      real ecf(istrip), egf(istrip), coct(istrip), cogt(istrip)
      real cogs1(istrip), cogs2(istrip), psyi(istrip), fac1(istrip)
      real rcdtc(istrip), rcdtg(istrip), rgdtc(istrip), rgdtg(istrip)
      real hcdtc(istrip), hcdtg(istrip), hcdtm(istrip)
      real hgdtc(istrip), hgdtg(istrip), hgdtm(istrip)
      real ecdtc(istrip), ecdtg(istrip), ecdsh(istrip)
      real egdtc(istrip), egdtg(istrip), egdsh(istrip)
      real ecidif(istrip), egidif(istrip)
      real currtotal(istrip), currsnow(istrip)

      capi = 1.0 / 4.0e-3
      timcon = 2.0 * pie / 86400.0
      dtc3xi = 1.0 / dtc3x
      fak = 0.01 * grav / cpair
      fah = 0.01 * grav / hlat

      if (prin) then
         print *,'dtc3x: ',dtc3x
         print *,'timcon: ',timcon
      endif

! Do our simple snow model if there is snow on the ground.  7/2/97 dmm
      do i = 1,istrip
         fluxd(i) = 0.0
         fluxs(i) = 0.0
         fluxmax(i) = 0.0
         fluxmelt(i) = 0.0

         if (ncount.eq.2) then
            currtotal(i) = (((w(i,1) * zdepth(i,1)) + (w(i,2) *        &
                           zdepth(i,2)) + (w(i,3) * zdepth(i,3))) *    &
                           poros(i)) + capac(i,1) + capac(i,2) +       &
                           roff(i) - ((ppl(i) + ppc(i)) / 1000.0)
            if ((abs(currtotal(i)-starttotal(i)).gt.0.00001).or.      &
               (prin.and.(i.eq.ipoint))) then
               print *,' '
               print *,'TEMRES CAPAC DIFF:',(currtotal(i)-starttotal(i))
               print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
               print *,'currtotal: ',currtotal(i)
               print *,'starttotal: ',starttotal(i)
               print *,'snow_1: ',snow(i,1)
               print *,'snow_2: ',snow(i,2)
               print *,'capac_1: ',capac(i,1)
               print *,'capac_2: ',capac(i,2)
               print *,'w_1: ',w(i,1)
               print *,'w_2: ',w(i,2)
               print *,'w_3: ',w(i,3)
               print *,'waterinlayer_1: ',(w(i,1)*poros(i)*zdepth(i,1))
               print *,'waterinlayer_2: ',(w(i,2)*poros(i)*zdepth(i,2))
               print *,'waterinlayer_3: ',(w(i,3)*poros(i)*zdepth(i,3))
               print *,'ppl: ',(ppl(i)/1000.)
               print *,'ppc: ',(ppc(i)/1000.)
               print *,'roff: ',roff(i)
               print *,'tm: ',tm(i)
               print *,'tc: ',tc(i)
               print *,'tg: ',tg(i)
               print *,'tsn: ',tsn(i)
               print *,'tgsb: ',tgsb(i)
            endif
            if (abs(currtotal(i)-starttotal(i)).gt.0.00001) then
               print *,'temres capac error!'
!               stop
            endif
            currsnow(i) = snow(i,1) + snow(i,2) - sfall(i) +ymelt(i,1) &
                          + ymelt(i,2) + yfrez(i,1) + yfrez(i,2)
            if ((abs(currsnow(i)-startsnow(i)).gt.0.00001).or.        &
               (prin.and.(i.eq.ipoint))) then
               print *,' '
               print *,'TEMRES SNOW DIFF: ',(currsnow(i)-startsnow(i))
               print *,'point: ',i,'  ',lonco(i),'  ',latco(i)
               print *,'currsnow: ',currsnow(i)
               print *,'startsnow: ',startsnow(i)
               print *,'snow_1: ',snow(i,1)
               print *,'snow_2: ',snow(i,2)
               print *,'capac_1: ',capac(i,1)
               print *,'capac_2: ',capac(i,2)
               print *,'sfall: ',sfall(i)
               print *,'ymelt_1: ',ymelt(i,1)
               print *,'ymelt_2: ',ymelt(i,2)
               print *,'yfrez_1: ',yfrez(i,1)
               print *,'yfrez_2: ',yfrez(i,2)
            endif
            if (abs(currsnow(i)-startsnow(i)).gt.0.00001) then
               print *,'temres snow error!'
!               stop
            endif
         endif

         if (sgfg(i).gt.0.5) then
! Recalculating snow density with fresh snow.  8/19/98 dmm
            if (ncount.eq.1) then
               sdens(i) = (((snow(i,2) - sfall(i)) * sdens(i))         &
                          + (sfall(i) * 100.)) / snow(i,2)
               css(i) = cs * (snow(i,2) - sskin)
               if (prin) then
                  print *,'sdens: ',sdens(ipoint)
                  print *,'css: ',css(ipoint)
                  print *,'cs: ',cs
                  print *,'snow2: ',snow(ipoint,2)
                  print *,'sskin: ',sskin
               endif
            endif
         endif
      enddo

      do i = 1,istrip
         wc(i) = min(1.0, capac(i,1)/satcap(i,1))
         wg(i) = min(1.0, capac(i,2)/satcap(i,2))
! Using blended temperature.  8/27/97 dmm
         if (tgsb(i).le.tf) then
            vcover(i,2) = 1.0
            wg(i) = min(1.0, capac(i,2)*capi)
            rst(i,2) = rsoil(i)
         endif
         ak(i) = fak / (psb(i) * bps(i))
         ah(i) = fah / psb(i)
         cci(i) = 1.0 / cc(i)
         rcp(i) = ros(i) * cpair
         if (prin) then
            print *,'rcp: ',rcp(i)
            print *,'ros: ',ros(i)
            print *,'cpair: ',cpair
         endif
         psyi(i) = rcp(i) / (cpair * ps(i) / (hlat * epsfac))
      enddo

      if (prin) then
         print *,'wc: ',wc(ipoint)
         print *,'wg: ',wg(ipoint)
         print *,'rst: ',rst(ipoint,2)
         print *,'ak: ',ak(ipoint)
         print *,'ah: ',ah(ipoint)
         print *,'cci: ',cci(ipoint)
         print *,'rcp: ',rcp(ipoint)
         print *,'psyi: ',psyi(ipoint)
      endif

      do i = 1,istrip
! First guesses for sensible heat
         d1(i) = (1.0 / ra(i)) + (1.0 / rb(i)) + (1.0 / rd(i))
         d1i(i) = rcp(i) / d1(i)
! Using blended temperature.  8/27/97 dmm
         ta(i) = ((tgsb(i) / rd(i)) + (tc(i) / rb(i)) + (thm(i) *      &
                 bps(i) / ra(i))) / d1(i)
! First Guess
         hc(i) = rcp(i) * (tc(i) - ta(i)) / rb(i) * dtc3x
! Using blended temperature.  8/26/97 dmm
         hg(i) = rcp(i) * (tgsb(i) - ta(i)) / rd(i) * dtc3x
         if ((prin).and.(ncount.eq.2)) then
            print *,'rcp: ',rcp(i)
            print *,'ta: ',ta(i)
            print *,'tgs: ',tgsb(i)
            print *,'tc: ',tc(i)
            print *,'thm: ',thm(i)
            print *,'rd: ',rd(i)
            print *,'rb: ',rb(i)
            print *,'ra: ',ra(i)
            print *,'d1: ',d1(i)
            print *,'bps: ',bps(i)
            print *,'dtc3x: ',dtc3x
         endif

! First guesses for latent heat
         hr(i) = hrr(i) * fg(i) + 1.0 - fg(i)
         rc(i) = rst(i,1) * fc(i) + rb(i) + rb(i)
         coc(i) = (1.0 - wc(i)) / rc(i) + wc(i) / (2.0 * rb(i))
! Over snow, rsurf = 0 and rg = rg + 500.  1/24/98 dmm
         rg(i) = rst(i,2) * fg(i) + (sgfg(i) * 500.0)
         rsurf(i) = rsoil(i) * fg(i) * (1.0 - sgfg(i))
! No longer changing conductance with snow cover, resistance modified
! above instead.  1/24/98 dmm
         cog1(i) = vcover(i,2) * (1.0 - wg(i)) / (rg(i) + rd(i)) +     &
                   (1.0 - vcover(i,2)) / (rsurf(i) + rd(i)) * hr(i) +  &
                   vcover(i,2) / (rsurf(i) + rd(i) + 44.0) * hr(i)
         cog2(i) = vcover(i,2) * (1.0 - wg(i)) / (rg(i) + rd(i)) +     &
                   (1.0 - vcover(i,2)) / (rsurf(i) + rd(i)) +          &
                   vcover(i,2) / (rsurf(i) + rd(i) + 44.)
         cog1(i) = cog1(i) + wg(i) / rd(i) * vcover(i,2)
         cog2(i) = cog2(i) + wg(i) / rd(i) * vcover(i,2)
         d2(i) = 1.0 / ra(i) + coc(i) + cog2(i)
         es(i) = sh(i) * ps(i) / (epsfac + sh(i))
         top(i) = coc(i) * etc(i) + cog1(i) * etg(i) + es(i) / ra(i)
! First guess
         ea(i) = top(i) / d2(i)
         ec(i) = (etc(i) - ea(i)) * coc(i) * psyi(i) * dtc3x
         eg(i) = (etg(i) * cog1(i) - ea(i) * cog2(i)) * psyi(i) * dtc3x
! Calculate if it is likely that with the previous snow temperature and
! given fluxes, that the new snow temperature will be above freezing.
! "stebt" is a boolean variable and determines amult.  1/24/98 dmm
         stebt(i) = (sgfg(i).gt.0.5).and.(tsn(i).gt.272.16).and.       &
                    (((radt(i,2) + swg(i) + sws(i)) - ((hg(i) + eg(i)) &
                    * dtc3xi)).gt.1.17)
         if (stebt(i)) then
            amult(i) = 0.0
         else
            amult(i) = 1.0
         endif
      enddo

      if (prin) then
         print *,'d1: ',d1(ipoint)
         print *,'d1i: ',d1i(ipoint)
         print *,'ta: ',ta(ipoint)
         print *,'hc: ',hc(ipoint)
         print *,'hg: ',hg(ipoint)
         print *,'rc: ',rc(ipoint)
         print *,'d2: ',d2(ipoint)
         print *,'top: ',top(ipoint)
         print *,'coc: ',coc(ipoint)
         print *,'cog1: ',cog1(ipoint)
         print *,'cog2: ',cog2(ipoint)
         print *,'etc: ',etc(ipoint)
         print *,'etgs: ',etg(ipoint)
         print *,'em: ',es(ipoint)
         print *,'ra: ',ra(ipoint)
         print *,'wg: ',wg(ipoint)
         print *,'rd: ',rd(ipoint)
         print *,'rg: ',rg(ipoint)
         print *,'rsurf: ',rsurf(ipoint)
         print *,'hr: ',hr(ipoint)
         print *,'vcover_2: ',vcover(ipoint,2)
         print *,'fg: ',fg(ipoint)
         print *,'rsoil: ',rsoil(ipoint)
         print *,'rst_2:',rst(ipoint,2)
         print *,'ea: ',ea(ipoint)
         print *,'ec: ',ec(ipoint)
         print *,'eg: ',eg(ipoint)
         print *,'stebt: ',stebt(ipoint)
         print *,'amult: ',amult(ipoint)
      endif

! Partial derivative calculations for sensible heat
      do i = 1,istrip
         hcdtc(i) = d1i(i) / rb(i) * ((1.0 / ra(i)) + (1.0 / rd(i)))
         hcdtg(i) = -d1i(i) / (rb(i) * rd(i))
! Disallow flux in case snow goes above freezing.  1/24/98 dmm
         hcdtg(i) = hcdtg(i) * amult(i)
         hcdtm(i) = -d1i(i) / (rb(i) * ra(i)) * bps(i)

         hgdtg(i) = d1i(i) / rd(i) * ((1.0 / ra(i)) + (1.0 / rb(i)))
! Disallow flux in case snow goes above freezing.  1/24/98 dmm
         hgdtg(i) = hgdtg(i) * amult(i)
         hgdtc(i) = -d1i(i) / (rd(i) * rb(i))
         hgdtm(i) = -d1i(i) / (rd(i) * ra(i)) * bps(i)
      enddo

      if (prin) then
         print *,'hcdtc: ',hcdtc(ipoint)
         print *,'hcdtg: ',hcdtg(ipoint)
         print *,'hcdtm: ',hcdtm(ipoint)
         print *,'hgdtc: ',hgdtc(ipoint)
         print *,'hgdtg: ',hgdtg(ipoint)
         print *,'hgdtc: ',hgdtm(ipoint)
      endif

! Partial derivative calculations for longwave radiation flux
      stb4 = 4.0 * stefan
      stb8 = 8.0 * stefan
      do i = 1,istrip
         fac1(i) = vcover(i,1) * (1.0 - thermk(i))
         rcdtc(i) = fac1(i) * stb8 * tc(i) * tc(i) * tc(i)
! Using blended temperature.  8/26/97 dmm
         rcdtg(i) = -fac1(i) * stb4 * tgsb(i) * tgsb(i) * tgsb(i)
! Disallow flux in case snow goes above freezing.  1/24/98 dmm
         rcdtg(i) = rcdtg(i) * amult(i)
         rgdtc(i) = -fac1(i) * stb4 * tc(i) * tc(i) * tc(i)
! Using blended temperature.  8/26/97 dmm
         rgdtg(i) = stb4 * tgsb(i) * tgsb(i) * tgsb(i)
! Disallow flux in case snow goes above freezing.  1/24/98 dmm
         rgdtg(i) = rgdtg(i) * amult(i)
      enddo

      if (prin) then
         print *,'fac1: ',fac1(ipoint)
         print *,'rcdtc: ',rcdtc(ipoint)
         print *,'rcdtg: ',rcdtg(ipoint)
         print *,'rgdtc: ',rgdtc(ipoint)
         print *,'rgdtg: ',rgdtg(ipoint)
      endif

! Partial derivative calculations for latent heat
      do i = 1,istrip
         deadtc(i) = btc(i) * coc(i) / d2(i)
! Disallow flux in case snow goes above freezing.  1/24/98 dmm
         btg(i) = btg(i) * amult(i)
         deadtg(i) = btg(i) * cog1(i) / d2(i)
! Disallow flux in case snow goes above freezing.  1/24/98 dmm
         deadtg(i) = deadtg(i) * amult(i)
         deadsh(i) = epsfac * ps(i) / ((epsfac + sh(i)) ** 2 * ra(i) * &
                     d2(i))

         ecdtc(i) = (btc(i) - deadtc(i)) * coc(i) * psyi(i)
         ecdtg(i) = -deadtg(i) * coc(i) * psyi(i)
! Disallow flux in case snow goes above freezing.  1/24/98 dmm
         ecdtg(i) = ecdtg(i) * amult(i)
         ecdsh(i) = -deadsh(i) * coc(i) * psyi(i)

         egdtg(i) = (btg(i) * cog1(i) - deadtg(i) * cog2(i)) * psyi(i)
! Disallow flux in case snow goes above freezing.  1/24/98 dmm
         egdtg(i) = egdtg(i) * amult(i)
         egdtc(i) = -deadtc(i) * cog2(i) * psyi(i)
         egdsh(i) = -deadsh(i) * cog2(i) * psyi(i)
      enddo

      if (prin) then
         print *,'deadtc: ',deadtc(ipoint)
         print *,'deadtg: ',deadtg(ipoint)
         print *,'deadsh: ',deadsh(ipoint)
         print *,'ecdtc: ',ecdtc(ipoint)
         print *,'ecdtg: ',ecdtg(ipoint)
         print *,'ecdsh: ',ecdsh(ipoint)
         print *,'egdtc: ',egdtc(ipoint)
         print *,'egdtg: ',egdtg(ipoint)
         print *,'egdsh: ',egdsh(ipoint)
      endif

! Solve for time changes of pbl and sib variables, using a
! semi-implicit scheme.  FILL the 4 x 5 ARRAY
      do i = 1,istrip
! Using blended capacity.  8/26/96 dmm.
! Removed ifblock for faster comput.; used flag variable.  1/24/98 dmm
! Here is the tg equation
         if (sgfg(i).gt.0.5) then
            tim = 1.0 + timcon * dtc3x * (1.0 + (1.0 / (cgsbi(i) *     &
                  css(i))))
         else
            tim = 1.0 + timcon * dtc3x
         endif
         pblsib(i,1,1) = tim + (dtc3x * cgsbi(i) * (hgdtg(i) +         &
                         egdtg(i) + rgdtg(i)))
         pblsib(i,1,2) = dtc3x * cgsbi(i) * (hgdtc(i) + egdtc(i) +     &
                         rgdtc(i))
         pblsib(i,1,3) = dtc3x * cgsbi(i) * hgdtm(i)
         pblsib(i,1,4) = dtc3x * cgsbi(i) * egdsh(i)
! Here is the tc equation
         pblsib(i,2,1) = dtc3x * cci(i) * (hcdtg(i) + ecdtg(i) +       &
                         rcdtg(i))
         pblsib(i,2,2) = 1.0 + dtc3x * cci(i) * (hcdtc(i) + ecdtc(i) + &
                         rcdtc(i))
         pblsib(i,2,3) = dtc3x * cci(i) * hcdtm(i)
         pblsib(i,2,4) = dtc3x * cci(i) * ecdsh(i)
! Here is the ts equation
         pblsib(i,3,1) = -dtc3x * ak(i) * (hgdtg(i) + hcdtg(i))
         pblsib(i,3,2) = -dtc3x * ak(i) * (hgdtc(i) + hcdtc(i))
         pblsib(i,3,3) = gmt(i,2) - dtc3x * ak(i) * (hgdtm(i) +hcdtm(i))
         pblsib(i,3,4) = 0.0
! Here is the sh equation
         pblsib(i,4,1) = -dtc3x * ah(i) * (egdtg(i) + ecdtg(i))
         pblsib(i,4,2) = -dtc3x * ah(i) * (egdtc(i) + ecdtc(i))
         pblsib(i,4,3) = 0.0
         pblsib(i,4,4) = gmq(i,2) - dtc3x * ah(i) * (egdsh(i) +ecdsh(i))
! Here are the time derivatives
         if (sgfg(i).gt.0.5) then
            pblsib(i,1,5) = (radt(i,2) - (hg(i) + eg(i)) * dtc3xi)     &
                            * cgsbi(i) - (timcon * (tgsb(i) - tg(i)    &
                            - (sws(i) * dtc3x / css(i))))
         else
            pblsib(i,1,5) = (radt(i,2) - (hg(i) + eg(i)) * dtc3xi)     &
                            * cgsbi(i) - (timcon * (tgsb(i) - td(i)))
         endif
         pblsib(i,2,5) = (radt(i,1) - (hc(i) + ec(i)) * dtc3xi) * cci(i)
         pblsib(i,3,5) = gmt(i,3) + ak(i) * (hg(i) + hc(i)) * dtc3xi
         pblsib(i,4,5) = gmq(i,3) + ah(i) * (eg(i) + ec(i)) * dtc3xi
      enddo

      if (prin) then
         do i = 1,4
            do j = 1,5
               print *,'pblsib: ',i,j,pblsib(1,i,j)
            enddo
         enddo
      endif

! Solve 4 x 5 matrix equation
      do i = 1,istrip
         ai(i) = 1.0 / pblsib(i,1,1)
         pblsib(i,1,2) = pblsib(i,1,2) * ai(i)
         pblsib(i,1,3) = pblsib(i,1,3) * ai(i)
         pblsib(i,1,4) = pblsib(i,1,4) * ai(i)
         pblsib(i,1,5) = pblsib(i,1,5) * ai(i)
         pblsib(i,1,1) = 1.0

         pblsib(i,2,2) = pblsib(i,2,2) - pblsib(i,2,1) * pblsib(i,1,2)
         pblsib(i,2,3) = pblsib(i,2,3) - pblsib(i,2,1) * pblsib(i,1,3)
         pblsib(i,2,4) = pblsib(i,2,4) - pblsib(i,2,1) * pblsib(i,1,4)
         pblsib(i,2,5) = pblsib(i,2,5) - pblsib(i,2,1) * pblsib(i,1,5)
         pblsib(i,3,2) = pblsib(i,3,2) - pblsib(i,3,1) * pblsib(i,1,2)
         pblsib(i,3,3) = pblsib(i,3,3) - pblsib(i,3,1) * pblsib(i,1,3)
         pblsib(i,3,4) = pblsib(i,3,4) - pblsib(i,3,1) * pblsib(i,1,4)
         pblsib(i,3,5) = pblsib(i,3,5) - pblsib(i,3,1) * pblsib(i,1,5)
         pblsib(i,4,2) = pblsib(i,4,2) - pblsib(i,4,1) * pblsib(i,1,2)
         pblsib(i,4,3) = pblsib(i,4,3) - pblsib(i,4,1) * pblsib(i,1,3)
         pblsib(i,4,4) = pblsib(i,4,4) - pblsib(i,4,1) * pblsib(i,1,4)
         pblsib(i,4,5) = pblsib(i,4,5) - pblsib(i,4,1) * pblsib(i,1,5)

         ai(i) = 1.0 / pblsib(i,2,2)
         pblsib(i,2,3) = pblsib(i,2,3) * ai(i)
         pblsib(i,2,4) = pblsib(i,2,4) * ai(i)
         pblsib(i,2,5) = pblsib(i,2,5) * ai(i)
         pblsib(i,2,2) = 1.0

         pblsib(i,3,3) = pblsib(i,3,3) - pblsib(i,3,2) * pblsib(i,2,3)
         pblsib(i,3,4) = pblsib(i,3,4) - pblsib(i,3,2) * pblsib(i,2,4)
         pblsib(i,3,5) = pblsib(i,3,5) - pblsib(i,3,2) * pblsib(i,2,5)
         pblsib(i,4,3) = pblsib(i,4,3) - pblsib(i,4,2) * pblsib(i,2,3)
         pblsib(i,4,4) = pblsib(i,4,4) - pblsib(i,4,2) * pblsib(i,2,4)
         pblsib(i,4,5) = pblsib(i,4,5) - pblsib(i,4,2) * pblsib(i,2,5)

         ai(i) = 1.0 / pblsib(i,3,3)
         pblsib(i,3,4) = pblsib(i,3,4) * ai(i)
         pblsib(i,3,5) = pblsib(i,3,5) * ai(i)
         pblsib(i,3,3) = 1.0

         pblsib(i,4,4) = pblsib(i,4,4) - pblsib(i,4,3) * pblsib(i,3,4)
         pblsib(i,4,5) = pblsib(i,4,5) - pblsib(i,4,3) * pblsib(i,3,5)

         pblsib(i,4,5) = pblsib(i,4,5) / pblsib(i,4,4)
         pblsib(i,3,5) = pblsib(i,3,5) - pblsib(i,3,4) * pblsib(i,4,5)
         pblsib(i,2,5) = pblsib(i,2,5) - pblsib(i,2,4) * pblsib(i,4,5) &
                         - pblsib(i,2,3) * pblsib(i,3,5)
         pblsib(i,1,5) = pblsib(i,1,5) - pblsib(i,1,4) * pblsib(i,4,5) &
                         - pblsib(i,1,3) * pblsib(i,3,5)               &
                         - pblsib(i,1,2) * pblsib(i,2,5)
      enddo

      do i = 1,istrip
! Calculate blended temperature change.  8/26/97 dmm
         dtgsb(i) = pblsib(i,1,5) * dtc3x
         dtsn(i) = dtgsb(i) * sgfg(i)
         dtg(i) = (dtg(i) * sgfg(i)) + (dtgsb(i) * (1.0 - sgfg(i)))
         dtc(i) = pblsib(i,2,5) * dtc3x
         dtm(i) = pblsib(i,3,5) * dtc3x
         dsh(i) = pblsib(i,4,5) * dtc3x
      enddo

      if (prin) then
         print *,'dtgsb: ',dtgsb(ipoint)
         print *,'dtsn: ',dtsn(ipoint)
         print *,'dtg: ',dtg(ipoint)
         print *,'dtc: ',dtc(ipoint)
         print *,'dtm: ',dtm(ipoint)
         print *,'dsh: ',dsh(ipoint)
      endif

      do i = 1,istrip
         if (sgfg(i).gt.0.5) then
! Update the snow conductivity as a function of snow density, based on
! Douville et al, 1995, Eq. 14.  Scaling by density of snow to ice,
! because our snow is liquid equiv. depth.  SMP3, Eq. 2.  08/24/98 dmm
            snowcon(i) = 2.22 * (sqrt(sdens(i) * 119.6) / 1000.) ** 1.88
! First weird number is 0.5 * (1 + sqrt(365)).  SMP3, Eq. A7
            condk(i) = 1.0 / ((10.0525 * zdepth(i,1) / chisl(i)) +     &
                       (1.0 / (4.0 * stefan * ((tg(i) + td(i)) / 2.)   &
                       ** 3)) + (0.5 * (snow(i,2) - sskin) /           &
                       (snowcon(i) * sdens(i) / 1000.)))
            fluxs(i) = timcon / cgsbi(i) * ((tsn(i) + (dtsn(i) * (1.0  &
                       + (sskin / (snow(i,2) - sskin))))) - tg(i) -    &
                       (sws(i) * dtc3x / css(i)))
            dtg(i) = (sws(i) + (condk(i) * (td(i) - tg(i))) +          &
                     fluxs(i)) / ((css(i) / dtc3x) + condk(i))
! Calculate fluxd backward implicit, then calculate the maximum flux.
! Difference is energy which is used to melt from below.  1/24/98 dmm
            fluxd(i) = condk(i) * (td(i) - (tg(i) + dtg(i)))
            fluxmax(i) = snowcon(i) * sdens(i) / 1000. * (tf - (tg(i)  &
                         + dtg(i))) / ((snow(i,2) - sskin) * 0.5)
            fluxmelt(i) = max(0., (fluxd(i) - fluxmax(i)))
            if (prin) then
               print *,'fluxd: ',fluxd(i)
               print *,'fluxmax: ',fluxmax(i)
               print *,'fluxmelt: ',fluxmelt(i)
               print *,'condk: ',condk(i)
               print *,'td: ',td(i)
               print *,'tg: ',tg(i)
               print *,'dtg: ',dtg(i)
               print *,'snowcon: ',snowcon(i)
               print *,'sdens: ',sdens(i)
               print *,'snow: ',snow(i,2)
            endif
         else
            snowcon(i) = -999.
            condk(i) = -999.
            fluxs(i) = -999.
         endif
      enddo

      do i = 1,istrip
! Use blended temperature change.  8/26/97 dmm
         hc(i) = hc(i) + dtc3x * (hcdtc(i) * dtc(i) + hcdtg(i) *       &
                 dtgsb(i) + hcdtm(i) * dtm(i))
         hg(i) = hg(i) + dtc3x * (hgdtc(i) * dtc(i) + hgdtg(i) *       &
                 dtgsb(i) + hgdtm(i) * dtm(i))
! Check if interception loss term has exceeded canopy storage
         ecpot(i) = (etc(i) - ea(i)) + (btc(i) - deadtc(i)) * dtc(i)   &
                    - deadtg(i) * dtgsb(i) - deadsh(i) * dsh(i)
         egpot(i) = (etg(i) - ea(i)) + (btg(i) - deadtg(i)) * dtgsb(i) &
                    - deadtc(i) * dtc(i) - deadsh(i) * dsh(i)
         epot(i) = (ecpot(i) + egpot(i)) * psyi(i) / ra(i) * dtc3x
      enddo

      hlat3 = 1.0e+03 * hlat

      do i = 1,istrip
         eci(i) = ecpot(i) * wc(i) * psyi(i) / (2.0 * rb(i)) * dtc3x
         ecidif(i) = max(0.0, eci(i)-capac(i,1)*hlat3)
         hc(i) = hc(i) + ecidif(i)
         eci(i) = min(eci(i), capac(i,1)*hlat3)
         egi(i) = egpot(i) * vcover(i,2) * wg(i) * psyi(i) / rd(i) *   &
                  dtc3x
         egidif(i) = max(0.0, egi(i)-capac(i,2)*hlat3)
         hg(i) = hg(i) + egidif(i)
         egi(i) = min(egi(i), capac(i,2)*hlat3)
! Evaporation is given in j m-2, calculated from gradients
! Over snow, rsurf = 0.  1/24/98 dmm
         rsurf(i) = rsoil(i) * fg(i) * (1.0 - sgfg(i))
         coct(i) = (1.0 - wc(i)) / rc(i)
         cogt(i) = vcover(i,2) * (1.0 - wg(i)) / (rg(i) + rd(i))
         cogs2(i) = (1.0 - vcover(i,2)) / (rd(i) + rsurf(i)) +         &
                    vcover(i,2) / (rd(i) + rsurf(i) + 44.)
         cogs1(i) = cogs2(i) * hr(i)
         ect(i) = ecpot(i) * coct(i) * psyi(i) * dtc3x
         ec(i) = eci(i) + ect(i)
         egt(i) = egpot(i) * cogt(i) * psyi(i) * dtc3x
! Use blended temperature change.  8/26/97 dmm
         egs(i) = (etg(i) + btg(i) * dtgsb(i)) * cogs1(i) - (ea(i) +   &
                  deadtg(i) * dtgsb(i) + deadtc(i) * dtc(i) +          &
                  deadsh(i) * dsh(i)) * cogs2(i)
         egs(i) = egs(i) * psyi(i) * dtc3x
         if (prin) then
            print *,'ea: ',ea(i)
            print *,'etg: ',etg(i)
            print *,'btg: ',btg(i)
            print *,'dtg: ',dtg(i)
            print *,'dtc: ',dtc(i)
            print *,'dsh: ',dsh(i)
            print *,'psy: ',psyi(i)
            print *,'cogs1: ',cogs1(i)
            print *,'cogs2: ',cogs2(i)
            print *,'deadtg: ',deadtg(i)
            print *,'deadtc: ',deadtc(i)
            print *,'deadsh: ',deadsh(i)
            print *,'egs: ',egs(i)
!            stop
         endif
         eg(i) = egt(i) + egs(i) + egi(i)
      enddo

! Test of dew condition. recalculation ensues if necessary
      do i = 1,istrip
! Use blended temperature change.  8/26/97 dmm
         radt(i,1) = radt(i,1) - rcdtc(i) * dtc(i) - rcdtg(i) * dtgsb(i)
         radt(i,2) = radt(i,2) - rgdtc(i) * dtc(i) - rgdtg(i) * dtgsb(i)
         ecf(i) = sign(1.0, ecpot(i)) * (fc(i) * 2.0 - 1.0)
         egf(i) = sign(1.0, egpot(i)) * (fg(i) * 2.0 - 1.0)
         if (ecf(i).le.0.0) then
            hc(i) = hc(i) + eci(i) + ect(i)
            eci(i) = 0.0
            ect(i) = 0.0
            ec(i) = 0.0
         endif
         if (egf(i).le.0.0) then
            hg(i) = hg(i) + egi(i) + egt(i) + egs(i)
            egi(i) = 0.0
            egt(i) = 0.0
            egs(i) = 0.0
            eg(i) = 0.0
         endif
      enddo

      return
      end

!***********************************************************************
      subroutine update(istrip,prin,ipoint,pie,dtc3x,hlat,snomel,tf,   &
           ra,rb,rd,btc,btg,etc,etg,bps,fluxd,tm,ea,tg,capac,snow,tc,  &
           td,tgsb,cc,cgsbi,css,sgfg,hc,hg,eg,eci,egi,ect,egt,egs,dtc, &
           dtg,dtgsb,dtm,dsh,dgdt,deadtc,deadtg,deadsh,sevap,bhf,chf,  &
           ghf,shf,eintmass,etranmass,ebaremass,etrancmass,eintcmass,  &
           eintgmass,etmass,hflux)

      implicit none
!***********************************************************************
!***  updating of soil moisture stores and interception capacity
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! pie          (-)       3.14159....
! dtc3x        (sec)     timestep
! hlat         (J kg-1)  latent heat of vaporization
! snomel       (J m-3)   latent heat of fusion X water density
! tf           (K)       freezing temperature
! ra           (sec m-1) aerodynamic resistance bet. reference & canopy
! rb           (sec m-1) bulk canopy boundary layer resistance
! rd           (sec m-1) ground to canopy air space resistance
! btc          (mb K-1)  first partial of satur. canopy vapor pressure
! btg          (mb K-1)  first partial of satur. ground vapor pressure
! etc          (mb)      saturation vapor pressure at canopy
! etg          (mb)      saturation vapor pressure at ground
! bps          (-)       (ps/p0)**(gasr/cpair); p0 = 1000.0 mb
! fluxd        (W m-2)   heat flux from ground to snow, SMP1, Eq. 2c
! tm           (K)       forcing reference-level temperature
! ea           (mb)      canopy air space vapor pressure
! tg           (K)       ground temperature
! tc           (K)       canopy temperature
! td           (K)       deep soil temperature
! tgsb         (K)       blended ground/snow temperature
! cc        (J m-2 K-1)  heat capacity of the canopy
! cgsbi    (J-1 m^2 K^1) inverse of heat capacity of blended ground/snow
! css       (J m-2 K-1)  heat capacity of bulk snow X water density
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
! hc           (J m-2)   canopy sensible heat flux
! hg           (J m-2)   ground sensible heat flux
! eg           (J m-2)   ground evaporation (egi + egt + egs)
! dtc          (K)       change in tc over timestep
! dtg          (K)       change in tg over timestep
! dtgsb        (K)       change in tgsb over timestep
! dtm          (K)       change in tm over timestep
! dsh          (kg kg-1) change in sh over timestep
! deadtc       (mb K-1)  change in ea w/r/t tc
! deadtg       (mb K-1)  change in ea w/r/t tg
! deadsh  (mb/(kg kg-1)) change in ea w/r/t sh
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! capac(cg)    (m)       liquid equiv. water&snow on canopy/ground
! snow(cg)     (m)       liquid equiv. snow on canopy/ground
! eci          (J m-2)   canopy interception loss
! egi          (J m-2)   ground interception loss
! ect          (J m-2)   canopy transpiration
! egt          (J m-2)   ground transpiration
! egs          (J m-2)   ground evaporation from the soil
! dgdt         (W m-2)   heat flux from tg layer to td layer
! sevap(3)     (m/dt)    snow evapor. canopy(1)/bare soil(2)/grnd cvr(3)
! bhf          (W m-2)   heat flux into the top snow layer from the air
! chf          (W m-2)   heat flux into the canopy from the air
! ghf          (W m-2)   heat flux into the ground from the air
! shf          (W m-2)   change in heat storage of bulk snow layer
! eintmass     (mm)      interception loss
! etranmass    (mm)      vegetation transpiration
! ebaremass    (mm)      evaporation from the soil
! etmass       (mm)      evaporation to the atmosphere (mass)
! hflux        (J m-2)   sensible heat flux to the atmosphere
!-----------------------------------------------------------------------
      logical prin
      integer istrip, ipoint, i
      real pie, dtc3x, hlat, snomel, tf
      real timcon, dtc3xi, hlati, hlat3i, snofac
      real ra(istrip), rb(istrip), rd(istrip)
      real btc(istrip), btg(istrip), etc(istrip), etg(istrip)
      real bps(istrip), fluxd(istrip), tm(istrip), ea(istrip)
      real tg(istrip), tc(istrip), td(istrip), tgsb(istrip)
      real cc(istrip), cgsbi(istrip), css(istrip), sgfg(istrip)
      real hc(istrip), hg(istrip), eg(istrip)
      real dtc(istrip), dtg(istrip), dtgsb(istrip)
      real dtm(istrip), dsh(istrip)
      real deadtc(istrip), deadtg(istrip), deadsh(istrip)
      real capac(istrip,2), snow(istrip,2)
      real eci(istrip), egi(istrip)
      real ect(istrip), egt(istrip), egs(istrip)
      real dgdt(istrip), sevap(istrip,3)
      real bhf(istrip), chf(istrip), ghf(istrip), shf(istrip)
      real eintmass(istrip), etranmass(istrip), ebaremass(istrip)
      real etmass(istrip), hflux(istrip)
      real tgen(istrip), tcen(istrip), tmen(istrip), taen(istrip)
      real eaen(istrip), d1(istrip), estarc(istrip), estarg(istrip)
      real facks(istrip), xsc(istrip), txsc(istrip)
      real ecmass(istrip), egmass(istrip)
      real etrancmass(istrip), etrangmass(istrip)
      real eintcmass(istrip), eintgmass(istrip)

! Adjustment of temperatures and vapor pressure, sensible heat fluxes.
! n.b. latent heat fluxes cannot be derived from estarc, estarg, ea due
! to linear result of implicit method
      do i = 1,istrip
! Use blended temperature change.  8/26/97 dmm
! Use blended temperature.  9/15/97 dmm
         tgen(i) = tgsb(i) + dtgsb(i)
         tcen(i) = tc(i) + dtc(i)
         tmen(i) = tm(i) + dtm(i)
         d1(i) = (1.0 / ra(i)) + (1.0 / rb(i)) + (1.0 / rd(i))
! Compute the fluxes consistent with the differencing scheme
         taen(i) = (tgen(i) / rd(i) + tcen(i) / rb(i) + tmen(i) *      &
                   bps(i) / ra(i)) / d1(i)
! Use blended temperature change.  8/26/97 dmm
         eaen(i) = ea(i) + deadtc(i) * dtc(i) + deadtg(i) * dtgsb(i) + &
                   deadsh(i) * dsh(i)
! Vapor pressures within the canopy and the moss
         estarc(i) = etc(i) + btc(i) * dtc(i)
! Use blended temperature change.  8/26/97 dmm
         estarg(i) = etg(i) + btg(i) * dtgsb(i)
      enddo

      do i = 1,istrip
         if (tgen(i).le.tf) then
            egs(i) = eg(i) - egi(i)
            egt(i) = 0.0
         endif
      enddo
! Heat fluxes into the canopy and the ground, in w m-2
      timcon = 2.0 * pie / 86400.0
      dtc3xi = 1.0 / dtc3x
      hlati = 1.0 / hlat
      hlat3i = 1.0 / (1.0e3 * hlat)

      do i = 1,istrip
         shf(i) = 0.0
         bhf(i) = 0.0
         chf(i) = dtc3xi * cc(i) * dtc(i)
! Re-did ghf calculation.  1/24/98 dmm
         if (sgfg(i).gt.0.5) then
            bhf(i) = (dtc3xi * css(i) * dtg(i))
            ghf(i) = fluxd(i)
            dgdt(i) = 0.0
         else
            ghf(i) = timcon * (1.0 / cgsbi(i)) * (tg(i) + dtg(i) -td(i))
            dgdt(i) = dtc3xi * (1.0 / cgsbi(i)) * dtg(i)
         endif
      enddo

! Evaporation losses are expressed in j m-2 : when divided by
! (hlat*1000.) loss is in m m-2
      snofac = 1.0 / (1.0 + snomel * hlat3i)
      do i = 1,istrip
         facks(i) = 1.0
         if (tcen(i).le.tf) facks(i) = snofac
         if ((ect(i)+eci(i)).le.0.0) then
            eci(i) = ect(i) + eci(i)
            ect(i) = 0.0
            facks(i) = 1.0 / facks(i)
         endif
      enddo

      do i = 1,istrip
         xsc(i) = eci(i) * facks(i) * hlat3i
         capac(i,1) = max((capac(i,1)-xsc(i)),0.0)
         if (eci(i).ge.0.0) then
            txsc(i) = max((snow(i,1)-capac(i,1)),0.0)
         else
            if (tcen(i).le.tf) then
               txsc(i) = xsc(i)
            else
               txsc(i) = 0.0
            endif
         endif
         snow(i,1) = snow(i,1) - txsc(i)
         sevap(i,1) = txsc(i)
! Mass terms are in kg m-2 dt-1
         ecmass(i) = (ect(i) + eci(i) * facks(i)) * hlati
         etrancmass(i) = ect(i) * hlati
         eintcmass(i) = eci(i) * facks(i) * hlati
      enddo

      do i = 1,istrip
         facks(i) = 1.0
         if (tgen(i).le.tf) facks(i) = snofac
         if ((egt(i)+egi(i)).le.0.0) then
            egi(i) = egt(i) + egi(i)
            egt(i) = 0.0
            facks(i) = 1.0 / facks(i)
         endif
      enddo

      do i = 1,istrip
         xsc(i) = egi(i) * facks(i) * hlat3i
         capac(i,2) = max((capac(i,2)-xsc(i)),0.0)
         if (egi(i).ge.0.0) then
            txsc(i) = max((snow(i,2)-capac(i,2)),0.0)
         else
            if (tgen(i).le.tf) then
               txsc(i) = xsc(i)
            else
               txsc(i) = 0.0
            endif
         endif
         snow(i,2) = snow(i,2) - txsc(i)
         sevap(i,3) = txsc(i)
         egmass(i) = (egt(i) + egs(i) + egi(i) * facks(i)) * hlati
! Total mass of water and total sensible heat lost from the veggies
         etmass(i) = ecmass(i) + egmass(i)
         etrangmass(i) = (egt(i) * hlati)
         etranmass(i) = etrancmass(i) + etrangmass(i)
         eintgmass(i) = (egi(i) * facks(i) * hlati)
         eintmass(i) = eintcmass(i) + eintgmass(i)
         ebaremass(i) = egs(i) * hlati
         hflux(i) = hc(i) + hg(i)
      enddo

      if (prin) then
         print *,'etmass: ',(etmass(ipoint)/1000.)
         print *,'eintmass: ',(eintmass(ipoint)/1000.)
         print *,'eintcmass: ',(eintcmass(ipoint)/1000.)
         print *,'eintgmass: ',(eintgmass(ipoint)/1000.)
         print *,'etranmass: ',(etranmass(ipoint)/1000.)
         print *,'etrancmass: ',(etrancmass(ipoint)/1000.)
         print *,'etrangmass: ',(etrangmass(ipoint)/1000.)
         print *,'ebaremass: ',(ebaremass(ipoint)/1000.)
         print *,'hflux: ',hflux(ipoint)
      endif

      return
      end

!***********************************************************************
      subroutine snowmt(istrip,dtc3x,snomel,cw,cs,sskin,tf,poros,      &
          zdepth,topostd,tg,w,capac,snow,tc,td,tsn,tgsb,cc,cg,cgsbi,   &
          sgfg,sdens,fluxd,fluxmax,fluxmelt,fluxef,dtsn,swg,ghf,shf,   &
          smelt,refreez,xmelt,xfrez,ymelt,yfrez,zmelt,zfrez,roff,      &
          srtes_qst)

      implicit none
!***********************************************************************
!***  calculation of snowmelt and modification of temperatures
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! dtc3x        (sec)     timestep
! snomel       (J m-3)   latent heat of fusion X water density
! cw        (J m-3 K-1)  heat capacity of liquid water X water density
! cs        (J m-3 K-1)  heat capacity of frozen water X water density
! sskin        (m)       depth of diurnal snow layer
! tf           (K)       freezing temperature
! poros        (0-1)     soil porosity
! zdepth(3)    (m)       exact independent depth of 3 soil layers
! topostd      (m)       topography standard deviat. 1x1deg from GTOPO30
! cc        (J m-2 K-1)  heat capacity of the canopy
! cg        (J m-2 K-1)  heat capacity of the ground
! cgsbi    (J-1 m^2 K^1) inverse of heat capacity of blended ground/snow
! sgfg         (0-1)     flag if snow model is "active" (0 = no;1 = yes)
! fluxd        (W m-2)   heat flux from ground to snow, SMP1, Eq. 2c
! fluxmax      (W m-2)   max. possible snow to ground flux, SMP1, Eq. 2e
! fluxmelt     (W m-2)   bottom melt flux, SMP1, Eq. 4b
! dtsn         (K)       change in tsn over timestep
! swg          (W m-2)   sw radiation absorbed in ground, SMP3, Eq. A6
! ghf          (W m-2)   heat flux into the ground from the air
! ymelt(cg)    (m/dt)    melted snow on canopy/ground from precip
! yfrez(cg)    (m/dt)    frozen water on canopy/ground from precip
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! tg           (K)       ground temperature
! w(3)         (0-1)     wetness of surface/root/recharge zone
! capac(cg)    (m)       liquid equiv. water&snow on canopy/ground
! snow(cg)     (m)       liquid equiv. snow on canopy/ground
! tc           (K)       canopy temperature
! td           (K)       deep soil temperature
! tsn          (K)       snow on ground temperature
! tgsb         (K)       blended ground/snow temperature
! sdens        (kg m-3)  density of bulk snow layer
! fluxef       (W m-2)   intermediate flux calc. for td modification
! shf          (W m-2)   change in heat storage of bulk snow layer
! zmelt(cg)    (m/dt)    total melted snow on canopy/ground
! zfrez(cg)    (m/dt)    total frozen water on canopy/ground
! roff         (m/dt)    total runoff
! srtes_qst    (m/dt)    liquid water flowing out of snowpack
!-----------------------------------------------------------------------
!     in subroutine parameters
!-----------------------------------------------------------------------
! roffs        (m/dt)    snowmelt water runoff from topographic variance
!-----------------------------------------------------------------------
      integer istrip, iveg, i
      real dtc3x, snomel, cw, cs, sskin, tf, srtes
      real poros(istrip), zdepth(istrip,3), topostd(istrip)
      real cc(istrip), cg(istrip), cgsbi(istrip)
      real sgfg(istrip), fluxd(istrip), fluxmax(istrip)
      real fluxmelt(istrip), dtsn(istrip), swg(istrip)
      real ghf(istrip), ymelt(istrip,2), yfrez(istrip,2)
      real tg(istrip), w(istrip,3), capac(istrip,2), snow(istrip,2)
      real tc(istrip), td(istrip), tsn(istrip), tgsb(istrip)
      real sdens(istrip), fluxef(istrip), shf(istrip)
      real zmelt(istrip), zfrez(istrip), roff(istrip)
      real xmelt(istrip,2), xfrez(istrip,2), roffs(istrip)
      real smelt(istrip), refreez(istrip), srtes_qst(istrip)
      real cct(istrip), tmp(istrip), water(istrip,2), change(istrip)
      real captfrz(istrip), snowp(istrip,2), dtd(istrip)

! Main veg loop
      do iveg = 1,2
         if (iveg.eq.1) then
            do i = 1,istrip
               cct(i) = cc(i)
               tmp(i) = tc(i)
! Initializing melt diagnostic for all time steps.  10/24/97 dmm
               xmelt(i,1) = 0.0
               xmelt(i,2) = 0.0
               xfrez(i,1) = 0.0
               xfrez(i,2) = 0.0
               smelt(i) = 0.0
               refreez(i) = 0.0
               srtes_qst(i) = 0.0
            enddo
         else
            do i = 1,istrip
! If we have snow, replace tg with tsn and cg with cs.  10/23/97 dmm
               cct(i) = cg(i) * (1.0 - sgfg(i))
               tmp(i) = (tg(i) * sgfg(i)) + (tgsb(i) * (1.0 - sgfg(i)))
            enddo
         endif

! Completely re-wrote the melt/freeze after 4x5 temperature calculation.
! Post-phase-change temp is 0C, which is modified in case more melts
! /freezes than is available to.  1/24/98 dmm
         do i = 1,istrip
            water(i,iveg) = capac(i,iveg) - snow(i,iveg)
            snowp(i,iveg) = snow(i,iveg)
            change(i) = 0.0
            if ((sgfg(i).gt.0.5).and.(iveg.eq.2)) then
               if (tmp(i).gt.tf) then
                  change(i) = ((snow(i,2) - sskin) * cs) * (tmp(i) - tf)
                  tmp(i) = tf
               endif
               if (tsn(i).gt.tf) then
                  change(i) = change(i) + ((sskin * cs) * (tsn(i) - tf))
                  tsn(i) = tf
               endif
               change(i) = change(i) / (((tf - tmp(i)) * cs) + snomel)
               if (change(i).gt.0.0) then
                  snow(i,2) = snow(i,2) - change(i)
                  water(i,2) = capac(i,2) - snow(i,2)
                  if (water(i,2).lt.0.0) then
                     tmp(i) = tmp(i) + water(i,2) * (((tf - tmp(i)) *  &
                              cs) + snomel) / (cs * capac(i,2))
                     snow(i,2) = capac(i,2)
                  endif
                  if (snow(i,2).lt.0.0) then
                     tmp(i) = tmp(i) + snow(i,2) * (((tf - tmp(i)) *   &
                              cs) + snomel) / (cw * capac(i,2))
                     snow(i,2) = 0.0
                  endif
               endif
            else
               change(i) = -(((water(i,iveg) * cw) + (snow(i,iveg) *   &
                           cs) + cct(i)) * (tmp(i) - tf)) / snomel
               tmp(i) = tf
               snow(i,iveg) = snow(i,iveg) + change(i)
               water(i,iveg) = capac(i,iveg) - snow(i,iveg)
               if (water(i,iveg).lt.0.0) then
                  tmp(i) = tmp(i) + water(i,iveg) * snomel / (cs *     &
                           capac(i,iveg) + cct(i))
                  snow(i,iveg) = capac(i,iveg)
               endif
               if (snow(i,iveg).lt.0.0) then
                  tmp(i) = tmp(i) - snow(i,iveg) * snomel / (cw *      &
                           capac(i,iveg) + cct(i))
                  snow(i,iveg) = 0.0
               endif
            endif
            change(i) = snowp(i,iveg) - snow(i,iveg)
            if (iveg.eq.1) then
               tc(i) = tmp(i)
               xmelt(i,1) = max(0.0, change(i))
               xfrez(i,1) = -max(0.0, -change(i))
            else
! If we have snow, replace tg with tsn and cg with cs.  8/28/97 dmm
               xmelt(i,2) = max(0.0, change(i))
               xfrez(i,2) = (-max(0.0, -change(i)))
! Added modification of snow density with refreeze.  8/19/98 dmm
               if ((sgfg(i).gt.0.5).and.(xfrez(i,2).ne.0.0)) then
                  sdens(i) = ((sdens(i) * (snow(i,2) + xfrez(i,2))) +  &
                             (-917.0 * xfrez(i,2))) / snow(i,2)
               endif
               tg(i) = tmp(i)
               tsn(i) = (tsn(i) * sgfg(i)) + (tmp(i) * (1.0 - sgfg(i)))
               tgsb(i) = (tsn(i) * sgfg(i)) + (tg(i) * (1.0 - sgfg(i)))
            endif
         enddo
      enddo

! Do all the melt for previous diagnostics right here.  10/24/97 dmm
      do i = 1,istrip
         roffs(i) = 0.0
         if (sgfg(i).gt.0.5) then
            if (fluxd(i).gt.fluxmax(i)) then
! Smelt also includes energy required to bring snow temperature to
! freezing before it melts.  11/25/97 dmm
               smelt(i) = (fluxmelt(i) * dtc3x) / (((tf - tg(i)) * cs) &
                          + snomel)
               smelt(i) = min(smelt(i), snow(i,2))
               smelt(i) = max(0.0, smelt(i))
               snow(i,2) = snow(i,2) - smelt(i)
            endif
         endif
! Diagnostic of total melt (precip plus heat plus flux)
         zmelt(i) = ymelt(i,1) + ymelt(i,2) + xmelt(i,1) + xmelt(i,2)  &
                    + smelt(i)
         if (zmelt(i).gt.0.0) then
            water(i,2) = capac(i,2) - snow(i,2)
            if (sgfg(i).gt.0.5) then
! Allow some melt to refreeze and warm ground temperature.  9/22/97 dmm
! Calculate refreeze, remainder of which enters soil.  9/22/97 dmm
               captfrz(i) = max((cg(i) * (tf - td(i)) / snomel), 0.0)
               if (captfrz(i).ge.water(i,2)) then
                  refreez(i) = water(i,2)
               else
                  refreez(i) = captfrz(i)
               endif
               if (refreez(i).gt.0.0) then
                  dtd(i) = refreez(i) * snomel / cg(i)
                  td(i) = td(i) + dtd(i)
                  snow(i,2) = snow(i,2) + refreez(i)
                  sdens(i) = ((sdens(i) * (snow(i,2) - refreez(i))) +  &
                             (917.0 * refreez(i))) / snow(i,2)
                  refreez(i) = -refreez(i)
               endif
               water(i,2) = capac(i,2) - snow(i,2)
! srtes = Snowmelt Ready To Enter Soil.  This is the sum of all of the
! snowmelts, after any refreeze of the ground.  Calculate additional
! overland flow before liquid gets into the soil.  Using linefit of
! basins, limiting to half of liquid to runoff.  The remaining liquid
! enters the soil.  SMP3, Eq. 8
               srtes = water(i,2)
               srtes_qst(i) = srtes
               roffs(i) = min((2.07 * topostd(i) / 1843.47),0.5) * srtes
               roff(i) = roff(i) + roffs(i)
               w(i,1) = w(i,1) + (srtes - roffs(i)) / (poros(i) *      &
                        zdepth(i,1))
               capac(i,2) = snow(i,2)
            endif
         endif
         zfrez(i) = xfrez(i,1) + xfrez(i,2) + yfrez(i,1) + yfrez(i,2)  &
                    + refreez(i)
      enddo

! End Main Loop
      do i = 1,istrip
! shf calculated here, not at end of update.  12/15/97 dmm
         water(i,2) = capac(i,2) - snow(i,2)
! No longer removing fluxd from fluxef calc.  12/15/97 dmm
         if (sgfg(i).gt.0.5) then
            shf(i) = dtsn(i) / (cgsbi(i) * dtc3x)
            fluxef(i) = swg(i) - fluxd(i)
         else
            fluxef(i) = ghf(i)
         endif
      enddo

      return
      end

!***********************************************************************
      subroutine hyssib_runoff(istrip,dtc3x,tf_drain,baseflow,ityp,phsat,&
           poros,bee,satco,slope,wiltp,zdepth,w,td,q3g,qqq,wattabl,    &
           roffg,roff,prin)

      implicit none
!***********************************************************************
!***  calculation of inter-layer moisture exchanges and runoffs
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
! dtc3x        (sec)     timestep
! tf_drain     (K)       temperature below which soil drainage stops
! baseflow     (sec-1)   baseflow timescale
! ityp         (1-13)    vegetation index
! phsat        (m)       soil moisture potential at saturation
! poros        (0-1)     soil porosity
! bee          (-)       Clapp-Hornberger empirical constant
! satco        (m/sec)   saturation hydraulic conductivity
! slope        (0-100)   average topographic slope in percent
! wiltp        (0-1)     wilting point wetness of soil, SMP3, Eq. 9
! zdepth(3)    (m)       exact independent depth of 3 soil layers
! td           (K)       deep soil temperature
! qqq          (m/sec)   water flux between soil layers
! qqq(3)       (m/sec)   evaporat. from layer 3 to atmosphere [positive]
!-----------------------------------------------------------------------
!     output parameters
!-----------------------------------------------------------------------
! w(3)         (0-1)     wetness of surface/root/recharge zone
! roff         (m/dt)    total runoff
! roffg        (m/dt)    baseflow runoff
! q3g          (m/sec)   gravitationally driven drainage from layer 3
! wattabl      (m)       water into layer 3 from water table [positive]
!-----------------------------------------------------------------------
!     in subroutine parameters
!-----------------------------------------------------------------------
! qqq(1)       (m/sec)   from layer 1 to layer 2 [positive]
! qqq(2)       (m/sec)   from layer 2 to layer 3 [positive]
! wgrav        (0-1)     minumum wetness for q3g drainage, SMP3, Eq. 10
!-----------------------------------------------------------------------
      integer istrip, isoil, i
      real pows, q3x, wmax, pmax, wmin, pmin, rsame, qmax, qmin, excess
      real avk, div, avkmin, avkmax, dpdwdz, denom, rdenom, fcb, oldw3
      real dtc3x, tf_drain, baseflow
      integer ityp(istrip)
      real phsat(istrip), poros(istrip), bee(istrip)
      real satco(istrip), slope(istrip), wiltp(istrip)
      real zdepth(istrip,3), td(istrip), qqq(istrip,3)
      real w(istrip,3), roff(istrip), roffg(istrip)
      real q3g(istrip), wattabl(istrip), wgrav(istrip)
      real temw(istrip,3), temwp(istrip,3), temwpp(istrip,3)
      real aaa(istrip,2), bbb(istrip,2), ccc(istrip,2), dpdw(istrip)
      logical prin

      do i = 1,istrip
         wgrav(i) = 0.0
         if (td(i).ge.tf_drain) then
            temw(i,1) = max(0.03, w(i,1))
            temwp(i,1) = temw(i,1) ** (-bee(i))
            temwpp(i,1) = min(1., temw(i,1)) ** (2. * bee(i) + 3.)
            temw(i,2) = max(0.03, w(i,2))
            temwp(i,2) = temw(i,2) ** (-bee(i))
            temwpp(i,2) = min(1., temw(i,2)) ** (2. * bee(i) + 3.)
            temw(i,3) = max(0.03, w(i,3))
            temwp(i,3) = temw(i,3) ** (-bee(i))
            temwpp(i,3) = min(1., temw(i,3)) ** (2. * bee(i) + 3.)
! SMP3, Eq. 10
            wgrav(i) = ((phsat(i) - (0.5 * zdepth(i,3))) / phsat(i))   &
                       ** (-1.0 / bee(i))
         endif
      enddo

! Calculation of gravitationally driven drainage from w(3) :
! taken as an integral of time varying conductivity.
! q3g = Sellers et al, 1986, Eq. 62
      do i = 1,istrip
         q3g(i) = 0.0
         if (td(i).ge.tf_drain) then
            pows = 2. * bee(i) + 2.
            q3x = temw(i,3) ** (-pows) + satco(i) / zdepth(i,3) /      &
                  poros(i) * slope(i) * pows * dtc3x
            q3x = q3x ** (1./pows)
            q3x = -(1./q3x - w(i,3)) * poros(i) * zdepth(i,3)/dtc3x
            q3x = max(0., q3x)
            if (prin) then
               print *,'q3g: ',q3x
               print *,'temw_3: ',temw(i,3)
               print *,'pows: ',pows
               print *,'bee: ',bee(i)
               print *,'satco: ',satco(i)
               print *,'zdepth_3: ',zdepth(i,3)
               print *,'poros: ',poros(i)
               print *,'slope: ',slope(i)
            endif
            q3g(i) = min(q3x, (max((w(i,3) - wgrav(i)),0.0) * poros(i) &
                     * zdepth(i,3) / dtc3x))

! Calculation of inter-layer exchanges of water due to gravitation and
! hydraulic gradient.  the values of w(x) + dw(x) are used to calculate
! the potential gradients between layers.  modified calculation of mean
! conductivities follows Milly and Eagleson, 1982, reduces recharge flux
! to top layer
! dpdw = estimated derivative of soil moisture potential with respect to
!        soil wetness.  assumption of gravitational drainage used to
!        estimate likely minimum wetness over the time step
! qqq  = Sellers et al, 1986, Eq. 61
! avk  = Milly and Eagleson, 1982, Eq. 4.14
            wmax = max(w(i,1), w(i,2), w(i,3), 0.05)
            wmax = min(wmax, 1.)
            pmax = wmax ** (-bee(i))
            wmin = 0.5*zdepth(i,1) + zdepth(i,2) + 0.5*zdepth(i,3)
            wmin = (pmax - wmin / phsat(i)) ** (-1. / bee(i))
            wmin = min(w(i,1), w(i,2), w(i,3), wmin)
            wmin = max(wmin, 0.049)
            pmin = wmin ** (-bee(i))
            dpdw(i) = phsat(i) * (pmax - pmin) / (wmax - wmin)
         endif
      enddo

      do isoil = 1,2
         do i = 1,istrip
            if (td(i).ge.tf_drain) then
               rsame = 0.0
               avk = temwp(i,isoil) * temwpp(i,isoil) -                &
                     temwp(i,isoil+1) * temwpp(i,isoil+1)
               div = temwp(i,isoil+1) - temwp(i,isoil)
               if (abs(div).lt.1.e-6) rsame = 1.
               avk = satco(i) * avk / ((1. + 3./bee(i)) * div + rsame)
               avkmin =satco(i)*min(temwpp(i,isoil),temwpp(i,isoil+1))
               avkmax =satco(i)*max(temwpp(i,isoil),temwpp(i,isoil+1)) &
                        * 1.01
               avk = max(avk, avkmin)
               avk = min(avk, avkmax)
               dpdwdz = dpdw(i)*2./(zdepth(i,isoil) + zdepth(i,isoil+1))
               aaa(i,isoil) = 1. + (1. / zdepth(i,isoil) + 1. /        &
                              zdepth(i,isoil+1)) * dtc3x / poros(i)    &
                              * avk * dpdwdz
               bbb(i,isoil) = -avk * dpdwdz / zdepth(i,2) * dtc3x /    &
                              poros(i)
               ccc(i,isoil) = avk* (dpdwdz * (w(i,isoil)-w(i,isoil+1)) &
                              + 1. + (isoil - 1) * dpdwdz * q3g(i) /   &
                              zdepth(i,3) * dtc3x / poros(i))
            endif
         enddo
      enddo

      do i = 1,istrip
         if (td(i).ge.tf_drain) then
            denom = (aaa(i,1) * aaa(i,2) - bbb(i,1) * bbb(i,2))
            rdenom = 0.0
            if (abs(denom).lt.1.e-6) rdenom = 1.
            rdenom = (1. - rdenom) / (denom + rdenom)
            qqq(i,1) = (aaa(i,2)*ccc(i,1) - bbb(i,1)*ccc(i,2)) * rdenom
            qqq(i,2) = (aaa(i,1)*ccc(i,2) - bbb(i,2)*ccc(i,1)) * rdenom
            w(i,3) = w(i,3) - q3g(i) * dtc3x / (poros(i) * zdepth(i,3))
            roff(i) = roff(i) + q3g(i) * dtc3x
         endif
      enddo

      isoil = 1
      do i = 1,istrip
         if (td(i).ge.tf_drain) then
            qmax = w(i,isoil) * (poros(i) * zdepth(i,isoil) / dtc3x)
            qmin = -w(i,isoil+1) * (poros(i) * zdepth(i,isoil+1)/ dtc3x)
            qqq(i,isoil) = min(qqq(i,isoil), qmax)
            qqq(i,isoil) = max(qqq(i,isoil), qmin)
            w(i,isoil) = w(i,isoil) - qqq(i,isoil) / (poros(i) *       &
                         zdepth(i,isoil) / dtc3x)
            w(i,isoil+1) = w(i,isoil+1) + qqq(i,isoil) / (poros(i) *   &
                           zdepth(i,isoil+1) / dtc3x)
         else
            qqq(i,isoil) = 0.0
         endif
      enddo
      isoil = 2
      do i = 1,istrip
         wattabl(i) = 0.0
         if (td(i).ge.tf_drain) then
            qmax = w(i,isoil) * (poros(i) * zdepth(i,isoil) / dtc3x)
            qmin = -w(i,isoil+1) * (poros(i) * zdepth(i,isoil+1)/ dtc3x)
            qqq(i,isoil) = min(qqq(i,isoil), qmax)
            qqq(i,isoil) = max(qqq(i,isoil), qmin)
            w(i,isoil) = w(i,isoil) - qqq(i,isoil) / (poros(i) *       &
                         zdepth(i,isoil) / dtc3x)
            fcb = 1.0
            if (w(i,isoil+1).lt.wiltp(i)) fcb = 1.0 + 1.0
            oldw3 = w(i,isoil+1)
            w(i,isoil+1) = w(i,isoil+1) + (qqq(i,isoil) - qqq(i,3)) /  &
                           (fcb * poros(i) * zdepth(i,isoil+1) / dtc3x)
            wattabl(i) = -(poros(i) * zdepth(i,3) * (w(i,isoil+1) -    &
                         oldw3) * (fcb - 1.0))
         else
            qqq(i,isoil) = 0.0
         endif
      enddo

      do isoil = 1,3
         do i = 1,istrip
            excess = max(0., (w(i,isoil)-1.))
            w(i,isoil) = w(i,isoil) - excess
            roff(i) = roff(i) + excess * poros(i) * zdepth(i,isoil)
         enddo
      enddo

! Compute baseflow runoffs and adjust layer 3 soil moisture
      do i = 1,istrip
         roffg(i) = 0.0
         if ((ityp(i).ne.13).and.(w(i,3).lt.wgrav(i))) then
            if ((td(i).gt.tf_drain).and.(w(i,3).gt.wiltp(i))) then
! Also multiplying baseflow by porosity.  SMP3, Eq. 9.  05/14/98 dmm
               roffg(i) = baseflow * (w(i,3) - wiltp(i)) * zdepth(i,3) &
                          * poros(i)
               w(i,3) = w(i,3) - roffg(i) / (zdepth(i,3) * poros(i))
               roff(i) = roff(i) + roffg(i)
            endif
         endif
      enddo

      return
      end

