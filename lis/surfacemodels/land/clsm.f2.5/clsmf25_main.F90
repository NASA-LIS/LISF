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
! !ROUTINE: clsmf25_main
! \label{clsmf25_main}
!
! !REVISION HISTORY:
! 13 Jun 2005 Rolf Reichle, Initial Specification
! 17 Dec 2005 Sujay Kumar, Implementation in LIS
! 23 Nov 2012: David Mocko, Added Catchment Fortuna-2.5
! 27 Mar 2013: David Mocko, Updates for ALMA outputs
!
! !INTERFACE:
subroutine clsmf25_main(nid)

  ! reichle, 13 Jun 2005
  ! qliu+reichle, 14 Aug 2008 - major overhaul, use subroutine louissurface()
  ! reichle, 29 Nov 2010 - added Helfand Monin-Obukhov surface layer turbulence scheme
  ! reichle, 20 Dec 2011 - reinstated "PARdrct" and "PARdffs" for MERRA-Land file specs
  ! reichle, 28 Dec 2011 - removed field "totalb" from "cat_diagn" structure

! !USES:
  use LIS_coreMod
  use LIS_FORC_AttributesMod 
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_albedoMod, only : LIS_alb
  use LIS_vegDataMod, only : LIS_lai,LIS_gfrac
  use LIS_logMod,    only : LIS_logunit, LIS_endrun
  use LIS_histDataMod
  use LIS_pluginIndices

  use clsmf25_lsmMod, only : clsmf25_struc
  use clsmf25_constants
  use clsmf25_MAPL_constants, ONLY: &
       stefan_boltzmann => MAPL_STFBOL, MAPL_USMIN,                    &
       lhe => MAPL_ALHL, lhf => MAPL_ALHF
  use clsmf25_types
  use clsmf25_esat_qsat
  use clsmf25_model
  use clsmf25_diagn_routines

  implicit none

  integer, intent(in) :: nid 
!
! !DESCRIPTION:
!  This is the entry point for calling the Catchment LSM physics. This routine
!  calls {\tt catchment} routine that performs the land surface computations, 
!  to solve water and energy equations.
!
!  The arguments are: 
!  \begin{description}
!  \item[nid]
!   index of the nest
!  \end{description}
!
!EOP
  integer :: N_cat    
  integer :: tid 
  integer :: sfc_turb_scheme
  integer :: i,index,zone
  real    :: dtstep             
  real    :: dec,omega
  real    :: ltime,coszth
  real    :: bowen_ratio,evap_frac,GWS_out
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index))   :: green, lai, sunang, alat
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index))   :: hsnacc, evacc, shacc, lhacc
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index))   :: SH_SNOW, AVET_SNOW, WAT_10CM
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index))   :: TOTWAT_SOIL, TOTICE_SOIL
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index))   :: startDSH, afterDSH, CSO
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index))   :: startDCC, afterDCC
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index))   :: startDSM, afterDSM
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index))   :: startDSW, afterDSW
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index))   :: startDIN, afterDIN
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index))   :: qsm_out, qfz_out, lwc
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index),3) :: fices
  real                                                   :: totalb_out

  ! --------------------------------------------------------------

  ! the following parameters must be consistent with: subroutine clsmf25_pmonth
!  real, parameter :: nodata_generic         = -9999.
!  real, parameter :: nodata_tolfrac_generic = 1.e-4
!  real :: nodata_tol_generic = abs(nodata_generic*nodata_tolfrac_generic)
!  logical :: ignore_SWNET_for_snow = .false.

  real,    parameter :: MIN_VEG_HEIGHT = 0.01           

  real,    parameter :: Z0_BY_ZVEG  = 0.13 
  real,    parameter :: D0_BY_ZVEG  = 0.66 
  real,    parameter :: HPBL        = 1000.                 

  real,    parameter :: DZE_MIN     = 10.

  real,    parameter :: TOTALB_MIN  = 0.01
  real,    parameter :: TOTALB_MAX  = 0.90

  ! parameters for computing land surface emissivity, taken from GEOS_CatchGridComp.F90 

  ! Emissivity values from Wilber et al (1999, NATA-TP-1999-209362)
  ! Fu-Liou bands have been combined to Chou bands (though these are broadband only)
  ! IGBP veg types have been mapped to Sib-Mosaic types 
  ! Details in ~suarez/Emiss on cerebus

  real,   parameter :: EMSVEG(8) = (/  &
       0.99560, 0.99000, 0.99560, 0.99320,            &
       0.99280, 0.99180, 0.94120, 0.94120 /)
  real,   parameter :: EMSSNO        =    0.99999


  ! definition of land in sfclayer module  

  integer, parameter :: SFCTYPE     = 3     

  ! parameters for the Helfand Monin-Obuhkhov surface layer routine

  integer, parameter :: HELF_NITER  = 6  ! number of internal iterations
  integer, parameter :: HELF_Z0     = 1  ! CHOOSEZ0 in GEOS_CatchGridComp.F90 (ocean only)

  ! --------------------------------------------------------------

  real :: SFRAC
  real :: BDSNOW
  integer :: gid
  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index)) :: &
       SQSCAT, Z2, DZM, RSOIL1, RSOIL2, SATCAP, RDC, U2FAC,            &
       TRAINL, PARDIR, PARDIF,                                         &
       SWNETF, SWNETS,                                                 &
       RA1, ETURB1, DEDTC1, HSTURB1, DHSDTC1,         &
       RA2, ETURB2, DEDTC2, HSTURB2, DHSDTC2,         &
       RA4, ETURB4, DEDTC4, HSTURB4, DHSDTC4,         &
       RAS, ETURBS, DEDTCS, HSTURBS, DHSDTCS,         &
       QSAT1, DQS1, ALW1, BLW1,                                        &
       QSAT2, DQS2, ALW2, BLW2,                                        &
       QSAT4, DQS4, ALW4, BLW4,                                        &
       QSATS, DQSS, ALWS, BLWS,                                        &
       DEDQA1, DHSDQA1,                                                &
       DEDQA2, DHSDQA2,                                                &
       DEDQA4, DHSDQA4,                                                &
       DEDQAS, DHSDQAS,                                                &
       WCHANGE, ECHANGE,                                  &
       TC1_0, TC2_0, TC4_0, QA1_0, QA2_0, QA4_0 , &
       ALBVR, ALBNR, ALBVF, ALBNF,                                     & 
       SNOVR, SNONR, SNOVF, SNONF

  real, dimension(LIS_rc%npatch(nid,LIS_rc%lsm_index)) :: zol, vgd, d0, dze, tmprealvec, uuu, zvg


  real, dimension(N_gt  ,  LIS_rc%npatch(nid,LIS_rc%lsm_index)) :: GHTCNT
  real, dimension(N_snow,  LIS_rc%npatch(nid,LIS_rc%lsm_index)) :: WESN, HTSN, SNDZ

  real :: albedo, albsn, asnow
  real :: qsattc, tmpt, tmpswnet, emis

  integer :: k, n

  logical :: snow_dominates_alb
  logical :: alarmCheck 
  character*3 :: fnest

  write(fnest,'(i3.3)') nid

  alarmCheck = LIS_isAlarmRinging(LIS_rc,"CLSM F2.5 model alarm "//trim(fnest))
  N_cat = LIS_rc%npatch(nid,LIS_rc%lsm_index)
  if ((alarmCheck).and.(N_cat.gt.0)) then 

       SQSCAT = 0 
       Z2 = 0 
       DZM = 0 
       RSOIL1 = 0 
       RSOIL2 = 0 
       SATCAP = 0 
       RDC = 0 
       U2FAC = 0
       TRAINL = 0 
       PARDIR = 0 
       PARDIF = 0 
       SWNETF = 0 
       SWNETS = 0 
       RA1 = 0 
       ETURB1 = 0 
       DEDTC1 = 0 
       HSTURB1 = 0 
       DHSDTC1 = 0
       RA2 = 0 
       ETURB2 = 0 
       DEDTC2 = 0 
       HSTURB2 = 0 
       DHSDTC2 = 0 
       RA4 = 0 
       ETURB4 = 0 
       DEDTC4 = 0 
       HSTURB4 = 0 
       DHSDTC4 = 0 
       RAS = 0 
       ETURBS = 0 
       DEDTCS = 0 
       HSTURBS = 0 
       DHSDTCS = 0 
       QSAT1 = 0 
       DQS1 = 0 
       ALW1 = 0 
       BLW1 = 0 
       QSAT2 = 0 
       DQS2 = 0 
       ALW2 = 0
       BLW2 = 0 
       QSAT4 = 0 
       DQS4 = 0 
       ALW4 = 0 
       BLW4 = 0 
       QSATS = 0 
       DQSS = 0 
       ALWS = 0 
       BLWS = 0 
       DEDQA1 = 0 
       DHSDQA1 = 0 
       DEDQA2 = 0 
       DHSDQA2 = 0 
       DEDQA4 = 0 
       DHSDQA4 = 0 
       DEDQAS = 0
       DHSDQAS = 0 
       WCHANGE = 0 
       ECHANGE = 0 
       TC1_0 = 0 
       TC2_0 = 0 
       TC4_0 = 0 
       QA1_0 = 0 
       QA2_0 = 0
       QA4_0 = 0  

     dtstep = LIS_rc%ts

!     gid = LIS_domain(nid)%gindex(LIS_surface(nid,LIS_rc%lsm_index)%tile(10)%col,&
!          LIS_surface(nid,LIS_rc%lsm_index)%tile(10)%row)
!     write(LIS_logunit,*) 'gid ',gid

     do n=1,N_cat
        clsmf25_struc(nid)%met_force(n)%Tair = &
             clsmf25_struc(nid)%met_force(n)%Tair/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%Qair = &
             clsmf25_struc(nid)%met_force(n)%Qair/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%Psurf = &
             clsmf25_struc(nid)%met_force(n)%Psurf/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%Rainf = &
             clsmf25_struc(nid)%met_force(n)%Rainf/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%Rainf_C = &
             clsmf25_struc(nid)%met_force(n)%Rainf_C/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%Snowf = &
             clsmf25_struc(nid)%met_force(n)%Snowf/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%LWdown = &
             clsmf25_struc(nid)%met_force(n)%LWdown/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%SWdown = &
             clsmf25_struc(nid)%met_force(n)%SWdown/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%Wind = &
             clsmf25_struc(nid)%met_force(n)%Wind/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%RefH = & 
             clsmf25_struc(nid)%met_force(n)%RefH/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%swnet = & 
             clsmf25_struc(nid)%met_force(n)%swnet/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%pardr = & 
             clsmf25_struc(nid)%met_force(n)%pardr/&
             clsmf25_struc(nid)%forc_count
        clsmf25_struc(nid)%met_force(n)%pardf = & 
             clsmf25_struc(nid)%met_force(n)%pardf/&
             clsmf25_struc(nid)%forc_count
     enddo
     ! Set all of the forcing points to some safe values wherever the
     ! forcing is undefined.  This is to allows the catchment routine
     ! to process these points without crashing.
     clsmf25_struc(nid)%good_forcing_mask = .true.
     do n = 1,N_cat
        if (clsmf25_struc(nid)%met_force(n)%tair.eq.LIS_rc%udef) then
           clsmf25_struc(nid)%good_forcing_mask(n) = .false.
           !          write(unit=666,fmt=*) 'GREP: forcing undef',n
        endif
        if (LIS_lai(nid)%tlai(n).eq.-9999) then
           clsmf25_struc(nid)%good_forcing_mask(n) = .false.
           !          write(unit=666,fmt=*) 'GREP: lai undef',n
        endif
     enddo
     do n = 1,N_cat
        if ( .not. clsmf25_struc(nid)%good_forcing_mask(n) ) then
           clsmf25_struc(nid)%met_force(n)%Tair = 293.0
           clsmf25_struc(nid)%met_force(n)%Qair = 0.02
           clsmf25_struc(nid)%met_force(n)%Psurf = 10000.0
           clsmf25_struc(nid)%met_force(n)%Rainf_C = 0.0
           clsmf25_struc(nid)%met_force(n)%Rainf = 0.0
           clsmf25_struc(nid)%met_force(n)%Snowf = 0.0
           clsmf25_struc(nid)%met_force(n)%LWdown = 500.0
           clsmf25_struc(nid)%met_force(n)%SWdown = 500.0
           clsmf25_struc(nid)%met_force(n)%Wind = 10.0
           clsmf25_struc(nid)%met_force(n)%RefH = 10.0
        endif
     enddo

     ! Calculate sun angle for current location and time
     do n = 1,N_cat
        index = LIS_domain(nid)%gindex(&
             LIS_surface(nid,LIS_rc%lsm_index)%tile(n)%col,&
             LIS_surface(nid,LIS_rc%lsm_index)%tile(n)%row)
        call LIS_localtime(LIS_rc%gmt,LIS_domain(nid)%grid(index)%lon,  &
             ltime,zone)
        call coszenith(LIS_domain(nid)%grid(index)%lon,                 &
             LIS_domain(nid)%grid(index)%lat,                 &
             ltime,zone,LIS_rc%doy,coszth,dec,omega)
        !       sunang(n) = acos(coszth)
        sunang(n) = coszth
        sunang(n) = amax1(sunang(n),0.01)
     enddo

     ! Set the greenness and lai from LIS.
     ! vgd is not needed and zol is set in clsmf25_pmonth.

     do n = 1,N_cat
        tid = LIS_surface(nid,LIS_rc%lsm_index)%tile(n)%tile_id
!        if(clsmf25_struc(nid)%usegreennessflag.eq.1) then 
!           green(n) = clsmf25_struc(nid)%cat_param(n)%grn(LIS_rc%mo)
!        else
           green(n) = LIS_gfrac(nid)%greenness(tid)
!        endif
!       if (clsmf25_struc(nid)%uselaiflag.eq.1) then
!          lai(n) = clsmf25_struc(nid)%cat_param(n)%lai(LIS_rc%mo)
!       else
          lai(n) = LIS_lai(nid)%tlai(tid)
!       endif
       zol(n) = LIS_rc%udef
       vgd(n) = LIS_rc%udef
       alat(n) = LIS_domain(nid)%grid(&
             LIS_domain(nid)%gindex(&
             LIS_surface(nid,LIS_rc%lsm_index)%tile(n)%col,&
             LIS_surface(nid,LIS_rc%lsm_index)%tile(n)%row))%lat

        ! Set all of the forcing points to some safe values wherever the
        ! LAI is undefined/zero.  This is to allows the catchment routine
        ! to process these points without crashing.
        if ( lai(n) <= 0.0 ) then
           if ( lai(n) < 0.0 ) then
              write(unit=666,fmt=*) 'GREP: resetting lai',n,lai(n)
           endif
           lai(n) = 0.001
        endif
     enddo

     ! ----------------------------------------------------------------

     ! initialize cat_diagn fields that are needed here and in sfc_turb 

     do n=1,N_cat

        clsmf25_struc(nid)%cat_diagn(n) = -9999.

     end do

     ! diagnose top layer snow temperature in K

     call get_tf_nd( N_cat, clsmf25_struc(nid)%cat_progn%htsn(1),           &
                            clsmf25_struc(nid)%cat_progn%wesn(1),           &
                            clsmf25_struc(nid)%cat_diagn%tpsn(1), fices(:,1))
     call get_tf_nd( N_cat, clsmf25_struc(nid)%cat_progn%htsn(2),           &
                            clsmf25_struc(nid)%cat_progn%wesn(2),           &
                            tmprealvec,                           fices(:,2))
     call get_tf_nd( N_cat, clsmf25_struc(nid)%cat_progn%htsn(3),           &
                            clsmf25_struc(nid)%cat_progn%wesn(3),           &
                            tmprealvec,                           fices(:,3))

     clsmf25_struc(nid)%cat_diagn%tpsn(1) = clsmf25_struc(nid)%cat_diagn%tpsn(1) + TF   ! convert to Kelvin

!     if(LIS_localPet.eq.24) then 
!        print*, clsmf25_struc(nid)%cat_progn(7311)%htsn(1), &
!             clsmf25_struc(nid)%cat_progn(7311)%wesn(1), &
!             clsmf25_struc(nid)%cat_diagn(7311)%tpsn(1)
!     endif
     ! diagnose sub-tile area fractions

     call calc_arX( N_cat, clsmf25_struc(nid)%cat_param, &
          clsmf25_struc(nid)%cat_progn, clsmf25_struc(nid)%cat_diagn%ar1, &
          clsmf25_struc(nid)%cat_diagn%ar2, tmprealvec )

     call calc_asnow( N_cat, clsmf25_struc(nid)%cat_progn, &
          clsmf25_struc(nid)%cat_diagn%asnow )

     ! ---------------------------------------------------------------

     ! pmonth updates the seasonally varying parameters

     call clsmf25_pmonth( &
          N_cat, clsmf25_struc(nid)%cat_param%vegcls, LIS_rc%doy, alat,      &
          green, lai, zol, vgd, .false.,                                 &
          SQSCAT, Z2, DZM, &
          RSOIL1, RSOIL2, SATCAP, RDC, U2FAC )

     ! ---------------------------------------------------------------

     ! precipitation

     TRAINL = clsmf25_struc(nid)%met_force%Rainf - &
          clsmf25_struc(nid)%met_force%Rainf_C

     SFRAC = 1.

     ! ---------------------------------------------------------------

     ! turbulence

     do n=1,N_cat

        ! change units of surface pressure from Pa to mbar

        !PSUR(n)   = clsmf25_struc(nid)%met_force(n)%Psurf/100.

        ! from GEOS_CatchGridComp.F90

        d0(n)  = D0_BY_ZVEG*zol(n)/Z0_BY_ZVEG        ! zero-plane displacement height

        dze(n) = max(clsmf25_struc(nid)%met_force(n)%RefH-d0(n), DZE_MIN)   

     end do

     ! compute surface exchange coefficients etc BEFORE possibly resetting
     ! profile of Qair-QAx-Qsat(surf) -- for consistency with two-stage 
     ! run-method in GEOS_CatchGridComp.F90
     ! reichle+qliu,  9 Oct 2008

     sfc_turb_scheme = clsmf25_struc(nid)%turbscheme

     select case (sfc_turb_scheme)

     case (0)     ! Louis

        !allocate(ch_tile(N_cat))
        !allocate(cq_tile(N_cat))

        call clsmf25_turb(N_cat, clsmf25_struc(nid)%met_force,&
             clsmf25_struc(nid)%cat_progn, &
             clsmf25_struc(nid)%cat_diagn, zol, dze, lai, &
             RA1, ETURB1, DEDQA1, DEDTC1, HSTURB1, DHSDQA1, DHSDTC1,         &
             RA2, ETURB2, DEDQA2, DEDTC2, HSTURB2, DHSDQA2, DHSDTC2,         &
             RA4, ETURB4, DEDQA4, DEDTC4, HSTURB4, DHSDQA4, DHSDTC4,         &
             RAS, ETURBS, DEDQAS, DEDTCS, HSTURBS, DHSDQAS, DHSDTCS,         & 
             sfc_turb_scheme )

     case (1)     ! Helfand Monin-Obukhov

        call clsmf25_turb(N_cat, clsmf25_struc(nid)%met_force, &
             clsmf25_struc(nid)%cat_progn, &
             clsmf25_struc(nid)%cat_diagn, zol, dze, lai, &
             RA1, ETURB1, DEDQA1, DEDTC1, HSTURB1, DHSDQA1, DHSDTC1,           &
             RA2, ETURB2, DEDQA2, DEDTC2, HSTURB2, DHSDQA2, DHSDTC2,           &
             RA4, ETURB4, DEDQA4, DEDTC4, HSTURB4, DHSDQA4, DHSDTC4,           &
             RAS, ETURBS, DEDQAS, DEDTCS, HSTURBS, DHSDQAS, DHSDTCS,           &
             sfc_turb_scheme )

     case default

        write(LIS_logunit,*) '[ERR] clsmf25_main():'
        write(LIS_logunit,*) 'Unknown selection for sfc_turb_scheme.'
        call LIS_endrun()

     end select



     do n=1,N_cat

        ! make sure that Qair-QAx-Qsat(surf) profile is ok
        ! (in a nutshell, need Qair<=QAx<=Qsat(surf) or Qair>=QAx>=Qsat(surf))
        ! - reichle, 5 Mar 2008 
        ! - qliu+reichle, 12 Aug 2008

        !qsattc = qsat(clsmf25_struc(nid)%cat_progn(n)%tc1,PSUR(n),ALHE)
        qsattc = qsat(clsmf25_struc(nid)%cat_progn(n)%tc1,&
             clsmf25_struc(nid)%met_force(n)%Psurf/100.)
        if (clsmf25_struc(nid)%cat_progn(n)%qa1>qsattc .and. &
             clsmf25_struc(nid)%cat_progn(n)%qa1>&
             clsmf25_struc(nid)%met_force(n)%Qair ) then
           clsmf25_struc(nid)%cat_progn(n)%qa1 = &
                max( clsmf25_struc(nid)%met_force(n)%Qair, qsattc )
        endif
        if (clsmf25_struc(nid)%cat_progn(n)%qa1<qsattc .and. &
             clsmf25_struc(nid)%cat_progn(n)%qa1<&
             clsmf25_struc(nid)%met_force(n)%Qair ) then
           clsmf25_struc(nid)%cat_progn(n)%qa1 = &
                min( clsmf25_struc(nid)%met_force(n)%Qair, qsattc )
        endif

        !qsattc = qsat(clsmf25_struc(nid)%cat_progn(n)%tc2,PSUR(n),ALHE)
        qsattc = qsat(clsmf25_struc(nid)%cat_progn(n)%tc2,&
             clsmf25_struc(nid)%met_force(n)%Psurf/100.)
        if (clsmf25_struc(nid)%cat_progn(n)%qa2>qsattc .and.&
             clsmf25_struc(nid)%cat_progn(n)%qa2>&
             clsmf25_struc(nid)%met_force(n)%Qair ) then
           clsmf25_struc(nid)%cat_progn(n)%qa2 = &
                max( clsmf25_struc(nid)%met_force(n)%Qair, qsattc )
        endif
        if (clsmf25_struc(nid)%cat_progn(n)%qa2<qsattc .and. &
             clsmf25_struc(nid)%cat_progn(n)%qa2<&
             clsmf25_struc(nid)%met_force(n)%Qair ) then
           clsmf25_struc(nid)%cat_progn(n)%qa2 = &
                min( clsmf25_struc(nid)%met_force(n)%Qair, qsattc )
        endif

        !qsattc = qsat(clsmf25_struc(nid)%cat_progn(n)%tc4,PSUR(n),ALHE)
!        if(LIS_rc%yr.eq.1983.and.LIS_rc%mo.eq.7.and.LIS_localPet.eq.87.and.&
!             n.eq.828) then 
!        if(                clsmf25_struc(nid)%cat_progn(n)%tc4.lt.220) then 
!           print*, n, LIS_localPet, &
!                clsmf25_struc(nid)%met_force(n)%Tair, &
!                clsmf25_struc(nid)%met_force(n)%Qair, &
!                clsmf25_struc(nid)%met_force(n)%Psurf, &
!                clsmf25_struc(nid)%met_force(n)%Rainf, &
!                clsmf25_struc(nid)%met_force(n)%Rainf_C, &
!                clsmf25_struc(nid)%met_force(n)%Snowf, &
!                clsmf25_struc(nid)%met_force(n)%SWdown, &
!                clsmf25_struc(nid)%met_force(n)%LWdown, &
!                clsmf25_struc(nid)%met_force(n)%Wind, &
!                clsmf25_struc(nid)%met_force(n)%RefH, &
!                clsmf25_struc(nid)%cat_progn(n)%tc1,&
!                clsmf25_struc(nid)%cat_progn(n)%tc2,&
!                clsmf25_struc(nid)%cat_progn(n)%tc4
!           stop
!        endif
!
!        if(LIS_localPet.eq.87.and.&
!             n.eq.828) then 
!           write(112,*) n,LIS_localPet,&
!                clsmf25_struc(nid)%met_force(n)%Tair, &
!                clsmf25_struc(nid)%met_force(n)%Qair, &
!                clsmf25_struc(nid)%met_force(n)%Psurf, &
!                clsmf25_struc(nid)%met_force(n)%Rainf, &
!                clsmf25_struc(nid)%met_force(n)%Rainf_C, &
!                clsmf25_struc(nid)%met_force(n)%Snowf, &
!                clsmf25_struc(nid)%met_force(n)%SWdown, &
!                clsmf25_struc(nid)%met_force(n)%LWdown, &
!                clsmf25_struc(nid)%met_force(n)%Wind, &
!                clsmf25_struc(nid)%met_force(n)%RefH, &
!                clsmf25_struc(nid)%cat_progn(n)%tc1,&
!                clsmf25_struc(nid)%cat_progn(n)%tc2,&
!                clsmf25_struc(nid)%cat_progn(n)%tc4
!        endif
        qsattc = qsat(clsmf25_struc(nid)%cat_progn(n)%tc4,&
             clsmf25_struc(nid)%met_force(n)%Psurf/100.)
        if (clsmf25_struc(nid)%cat_progn(n)%qa4>qsattc .and. &
             clsmf25_struc(nid)%cat_progn(n)%qa4>&
             clsmf25_struc(nid)%met_force(n)%Qair ) then
           clsmf25_struc(nid)%cat_progn(n)%qa4 = &
                max( clsmf25_struc(nid)%met_force(n)%Qair, qsattc )
        endif
        if (clsmf25_struc(nid)%cat_progn(n)%qa4<qsattc .and. &
             clsmf25_struc(nid)%cat_progn(n)%qa4<&
             clsmf25_struc(nid)%met_force(n)%Qair ) then
           clsmf25_struc(nid)%cat_progn(n)%qa4 = &
                min( clsmf25_struc(nid)%met_force(n)%Qair, qsattc )
        endif

     end do

     ! convert ght and snow structures into arrays needed by sibalb() and catchment()

     do k=1,N_gt

        GHTCNT(k,:) = clsmf25_struc(nid)%cat_progn%ght(k)

     end do
     
     do k=1,N_snow

        WESN(k,:) = clsmf25_struc(nid)%cat_progn%wesn(k)
        HTSN(k,:) = clsmf25_struc(nid)%cat_progn%htsn(k)
        SNDZ(k,:) = clsmf25_struc(nid)%cat_progn%sndz(k)
     end do


    if (clsmf25_struc(nid)%usemodisalbflag.eq.1) then
       call SIBALB(N_cat, clsmf25_struc(nid)%cat_param%vegcls,        &
            LAI, green, sunang,                                   & 
            clsmf25_struc(nid)%modis_alb_param(:,LIS_rc%mo)%albvr, & ! use only diffuse
            clsmf25_struc(nid)%modis_alb_param(:,LIS_rc%mo)%albvf,  & ! use only diffuse
            clsmf25_struc(nid)%modis_alb_param(:,LIS_rc%mo)%albnr,  & ! use only diffuse
            clsmf25_struc(nid)%modis_alb_param(:,LIS_rc%mo)%albnf,  & ! use only diffuse
            WESN, SNDZ,                                           & 
            ALBVR, ALBNR, ALBVF, ALBNF,                           & ! snow-free albedos
            SNOVR, SNONR, SNOVF, SNONF  )                           ! snow albedos
    endif

     do n=1,N_cat

        ! tentative albedo for snow-free and snow-covered portion
        ! (naively assumes that each SW radiation component contributes 0.25 to 
        !  total SWdown -- these tentative albedos will NOT be used as such
        !  if SWnet is provided in forcing data)
        ! reichle+qliu,  9 Oct 2008
        if (clsmf25_struc(nid)%usemodisalbflag.eq.1) then
           albedo = min(1., 0.25*(ALBVR(n) + ALBNR(n) + ALBVF(n) + ALBNF(n)))
           albsn =  min(1., 0.25*(SNOVR(n) + SNONR(n) + SNOVF(n) + SNONF(n)))
        else
           tid = LIS_surface(nid,LIS_rc%lsm_index)%tile(n)%tile_id           
           albedo = LIS_alb(nid)%albsf(tid)
           albsn = LIS_alb(nid)%mxsnalb(tid)
        endif
        clsmf25_struc(nid)%cat_output(n)%albsn = albsn
        
        !SWNETF(n) = (1.-albedo) * clsmf25_struc(nid)%met_force(n)%SWdown
        !SWNETS(n) = (1.-albsn ) * clsmf25_struc(nid)%met_force(n)%SWdown

        ! diagnose total albedo and upward shortwave

        call calc_asnow( 1, clsmf25_struc(nid)%cat_progn(n:n), tmprealvec(1) )

        asnow = tmprealvec(1)     ! will also be needed for emissivity calcs below

        ! total albedo calculation requires a threshold "SWDN_THRESHOLD" for minimum 
        ! allowable positive SWdown as implemented in subroutine repair_forcing()
        ! qliu+reichle, 21 Aug 2008

        ! improved treatment of albedo and SWnet calculations that take advantage
        ! of SWnet forcing to obviate the need for individual SW radiation components
        ! reichle+qliu,  9 Oct 2008

        if (clsmf25_struc(nid)%met_force(n)%SWdown>0. .and. sunang(n)>0.) then

           ! sun is up and SWdown>0

           ! determine whether snow or snow-free albedo dominates
           if((asnow.eq.0.0).and.(albedo.eq.1.0)) then 
              snow_dominates_alb = .false. 
           else
              snow_dominates_alb = ( (1-asnow)*(1-albedo)<=asnow*(1-albsn) )
           endif
           ! compute SWNET for snow-free and snow conditions

           !          if ( (abs(clsmf25_struc(nid)%met_force(n)%SWnet-nodata_generic)<nodata_tol_generic) &
           !               .or.                                                        &
           !               (snow_dominates_alb .and. ignore_SWNET_for_snow) )  then

           ! SWnet is not available or wanted, use simplistic albedos (see above)
           if(LIS_FORC_SWnet%selectOpt.eq.0) then 
              
              SWNETF(n) = (1.-albedo) * clsmf25_struc(nid)%met_force(n)%SWdown
              SWNETS(n) = (1.-albsn ) * clsmf25_struc(nid)%met_force(n)%SWdown
              
              tmpswnet = (1-asnow)*SWNETF(n) + asnow* SWNETS(n)
           else
              
              
              ! use knowlegdge of SWnet to back out better estimate of albedo
              ! (works well if *all* tiles in forcing grid cell are either
              ! completely snow-covered or completely snow-free)
              
              if (.not. snow_dominates_alb) then
                 
                 ! snow-free area dominates, use simplistic snow albedo
                 ! and back out better snow-free albedo
                 
                 SWNETS(n) = (1.-albsn ) * clsmf25_struc(nid)%met_force(n)%SWdown
                 SWNETF(n) = (clsmf25_struc(nid)%met_force(n)%SWnet-asnow*SWNETS(n))/(1.-asnow)
                 
              else                                       

                 ! snow-covered area dominates, use simplistic snow-free albedo
                 ! and back out better snow-covered albedo                
                 
                 SWNETF(n) = (1.-albedo ) * clsmf25_struc(nid)%met_force(n)%SWdown
                 SWNETS(n) = (clsmf25_struc(nid)%met_force(n)%SWnet-(1.-asnow)*SWNETF(n))/asnow
                 
              end if
              
              tmpswnet = clsmf25_struc(nid)%met_force(n)%SWnet
              
           end if

           clsmf25_struc(nid)%cat_diagn(n)%totalb = &
                min(max(1. - tmpswnet/clsmf25_struc(nid)%met_force(n)%SWdown, &
                TOTALB_MIN), TOTALB_MAX)

           clsmf25_struc(nid)%cat_diagn(n)%swup   =&
                clsmf25_struc(nid)%cat_diagn(n)%totalb* &
                clsmf25_struc(nid)%met_force(n)%SWdown

        else

           SWNETS(n)           = 0.
           SWNETF(n)           = 0.

           clsmf25_struc(nid)%cat_diagn(n)%totalb = LIS_rc%udef
           clsmf25_struc(nid)%cat_diagn(n)%swup   = 0.

        end if

        ! photo-synthetically active radiation
        ! if PAR components are not available:
        !
        !   assume half of SWdown is photosynthetically active
        !   assume half of PAR is direct, half diffuse

        if(LIS_FORC_Pardr%selectOpt.eq.1) then 
           pardir(n) = clsmf25_struc(nid)%met_force(n)%pardr
        else
           PARDIR(n) = 0.5*0.5*clsmf25_struc(nid)%met_force(n)%SWdown          
              
        endif
        if(LIS_FORC_Pardf%selectOpt.eq.1) then 
           pardif(n) = clsmf25_struc(nid)%met_force(n)%pardf
        else
           PARDIF(n) = PARDIR(n)
        endif


        QSAT1(n)   = qsat(clsmf25_struc(nid)%cat_progn(n)%tc1, &
             clsmf25_struc(nid)%met_force(n)%Psurf/100.)
        QSAT2(n)   = qsat(clsmf25_struc(nid)%cat_progn(n)%tc2, &
             clsmf25_struc(nid)%met_force(n)%Psurf/100.)
        QSAT4(n)   = qsat(clsmf25_struc(nid)%cat_progn(n)%tc4, &
             clsmf25_struc(nid)%met_force(n)%Psurf/100.)

        QSATS(n)   = qsat(clsmf25_struc(nid)%cat_diagn(n)%tpsn(1), &
             clsmf25_struc(nid)%met_force(n)%Psurf/100.)


        ! compute ALWx, BLWx

        ! tile-average surface emissivity 

        emis = EMSVEG(clsmf25_struc(nid)%cat_param(n)%vegcls) + &
             (EMSVEG(8) - EMSVEG(clsmf25_struc(nid)%cat_param(n)%vegcls))*exp(-lai(n))

        emis = emis*(1.-asnow) + EMSSNO*asnow

        ! ALWx and BLWx are computed separately for each sub-tile 
        ! (unlike in GEOScatch_GridComp.F90)

        tmpt    = clsmf25_struc(nid)%cat_progn(n)%tc1
        BLW1(n) = emis*stefan_boltzmann*tmpt*tmpt*tmpt
        ALW1(n) = -3.*tmpt*BLW1(n)
        BLW1(n) =  4.*BLW1(n)

        tmpt    = clsmf25_struc(nid)%cat_progn(n)%tc2
        BLW2(n) = emis*stefan_boltzmann*tmpt*tmpt*tmpt
        ALW2(n) = -3.*tmpt*BLW2(n)
        BLW2(n) =  4.*BLW2(n)

        tmpt    = clsmf25_struc(nid)%cat_progn(n)%tc4
        BLW4(n) = emis*stefan_boltzmann*tmpt*tmpt*tmpt
        ALW4(n) = -3.*tmpt*BLW4(n)
        BLW4(n) =  4.*BLW4(n)

        tmpt    = clsmf25_struc(nid)%cat_diagn(n)%tpsn(1)
        BLWS(n) = emis*stefan_boltzmann*tmpt*tmpt*tmpt
        ALWS(n) = -3.*tmpt*BLWS(n)
        BLWS(n) =  4.*BLWS(n)

     end do


     ! get wind-related variable for input to catchment()

     zvg = zol/Z0_BY_ZVEG 

     uuu = max(clsmf25_struc(nid)%met_force%Wind, MAPL_USMIN) *&
          (log((max(zvg-d0,0.)+zol)/zol) &
          / log((max(clsmf25_struc(nid)%met_force%RefH-d0,10.)+zol)/zol))

! Calculate "before" states for DelSurfHeat, DelColdCont, DelSoilMoist
    do n = 1,N_cat
       CSO(n) = CSOIL_1
       if (clsmf25_struc(nid)%cat_param(n)%vegcls.ne.1) CSO(n) = CSOIL_2
       startDSH(n) = ((clsmf25_struc(nid)%cat_progn(n)%tc1 *           &
                      clsmf25_struc(nid)%cat_diagn(n)%ar1) +           &
                     (clsmf25_struc(nid)%cat_progn(n)%tc2 *            &
                      clsmf25_struc(nid)%cat_diagn(n)%ar2) +           &
                     (clsmf25_struc(nid)%cat_progn(n)%tc4 *            &
                      (1.0 - clsmf25_struc(nid)%cat_diagn(n)%ar1       &
                           - clsmf25_struc(nid)%cat_diagn(n)%ar2))) * CSO(n)
       startDCC(n) = sum(clsmf25_struc(nid)%cat_progn(n)%htsn(1:N_snow))
       startDSM(n) = (clsmf25_struc(nid)%cat_param(n)%dzpr *           &
                      clsmf25_struc(nid)%cat_param(n)%poros) -         &
                     clsmf25_struc(nid)%cat_progn(n)%catdef +          &
                     clsmf25_struc(nid)%cat_progn(n)%rzexc +           &
                     clsmf25_struc(nid)%cat_progn(n)%srfexc
       startDSW(n) = sum(clsmf25_struc(nid)%cat_progn(n)%wesn(1:N_snow))
       startDIN(n) = clsmf25_struc(nid)%cat_progn(n)%capac
       lwc(n) =     (clsmf25_struc(nid)%cat_progn(n)%wesn(1) *         &
                    (1.0 - fices(n,1)))                  +             &
                    (clsmf25_struc(nid)%cat_progn(n)%wesn(2) *         &
                    (1.0 - fices(n,2)))                  +             &
                    (clsmf25_struc(nid)%cat_progn(n)%wesn(3) *         &
                    (1.0 - fices(n,3)))
    enddo
    
!    if(LIS_localPet.eq.165) then 
!       write(LIS_logunit,*) 'tp1 ',clsmf25_struc(nid)%cat_diagn(22527)%tpsn(1),&
!            clsmf25_struc(nid)%cat_progn(22527)%tc1,&
!            clsmf25_struc(nid)%cat_progn(22527)%tc2,&
!            clsmf25_struc(nid)%cat_progn(22527)%tc4
!    endif
     ! call to Catchment model 
     call clsmf25( N_cat, dtstep, sfrac,                                 &
          clsmf25_struc(nid)%cat_param%cat_id,&
          clsmf25_struc(nid)%cat_param%vegcls,&
          clsmf25_struc(nid)%cat_param%dzsf,&
          clsmf25_struc(nid)%met_force%Rainf_C, &
          TRAINL, clsmf25_struc(nid)%met_force%Snowf, uuu,                &
          ETURB1, DEDQA1, DEDTC1, HSTURB1, DHSDQA1, DHSDTC1,              &
          ETURB2, DEDQA2, DEDTC2, HSTURB2, DHSDQA2, DHSDTC2,              &
          ETURB4, DEDQA4, DEDTC4, HSTURB4, DHSDQA4, DHSDTC4,              &
          ETURBS, DEDQAS, DEDTCS, HSTURBS, DHSDQAS, DHSDTCS,              &
          clsmf25_struc(nid)%met_force%Tair, &
          clsmf25_struc(nid)%met_force%Qair,   &
          RA1, RA2, RA4, RAS, sunang, PARDIR, PARDIF,                     &
          SWNETF, SWNETS, clsmf25_struc(nid)%met_force%LWdown,            &
          clsmf25_struc(nid)%met_force%Psurf/100.,                        &
          lai, green, Z2,                                                 &
          SQSCAT, RSOIL1, RSOIL2, RDC,                                    &
          QSAT1, DQS1, ALW1, BLW1,                                        &
          QSAT2, DQS2, ALW2, BLW2,                                        &
          QSAT4, DQS4, ALW4, BLW4,                                        &
          QSATS, DQSS, ALWS, BLWS,                                        &
          clsmf25_struc(nid)%cat_param%bf1,     &
          clsmf25_struc(nid)%cat_param%bf2,     &
          clsmf25_struc(nid)%cat_param%bf3,     &
          clsmf25_struc(nid)%cat_param%vgwmax,  &
          clsmf25_struc(nid)%cat_param%cdcr1, &
          clsmf25_struc(nid)%cat_param%cdcr2, &
          clsmf25_struc(nid)%cat_param%psis, &
          clsmf25_struc(nid)%cat_param%bee,    &
          clsmf25_struc(nid)%cat_param%poros,                  &
          clsmf25_struc(nid)%cat_param%wpwet, &
          clsmf25_struc(nid)%cat_param%cond,  &
          clsmf25_struc(nid)%cat_param%gnu,   &
          clsmf25_struc(nid)%cat_param%ars1, &
          clsmf25_struc(nid)%cat_param%ars2,   &
          clsmf25_struc(nid)%cat_param%ars3,   &
          clsmf25_struc(nid)%cat_param%ara1,  &
          clsmf25_struc(nid)%cat_param%ara2,   &
          clsmf25_struc(nid)%cat_param%ara3, &
          clsmf25_struc(nid)%cat_param%ara4,   &
          clsmf25_struc(nid)%cat_param%arw1, &
          clsmf25_struc(nid)%cat_param%arw2,   &
          clsmf25_struc(nid)%cat_param%arw3, &
          clsmf25_struc(nid)%cat_param%arw4,   &
          clsmf25_struc(nid)%cat_param%tsa1, &
          clsmf25_struc(nid)%cat_param%tsa2,   &
          clsmf25_struc(nid)%cat_param%tsb1, &
          clsmf25_struc(nid)%cat_param%tsb2,   &
          clsmf25_struc(nid)%cat_param%atau, &
          clsmf25_struc(nid)%cat_param%btau,   &
          .false.,                                                        &
          clsmf25_struc(nid)%cat_progn%tc1, &
          clsmf25_struc(nid)%cat_progn%tc2,     &
          clsmf25_struc(nid)%cat_progn%tc4,     &
          clsmf25_struc(nid)%cat_progn%qa1, &
          clsmf25_struc(nid)%cat_progn%qa2,     &
          clsmf25_struc(nid)%cat_progn%qa4,     &
          clsmf25_struc(nid)%cat_progn%capac,   &
          clsmf25_struc(nid)%cat_progn%catdef, &
          clsmf25_struc(nid)%cat_progn%rzexc,  &
          clsmf25_struc(nid)%cat_progn%srfexc,  &
          GHTCNT, clsmf25_struc(nid)%cat_diagn%tsurf,   &
          WESN, HTSN, SNDZ,                                               &
          clsmf25_struc(nid)%cat_diagn%evap, &
          clsmf25_struc(nid)%cat_diagn%shflux, &
          clsmf25_struc(nid)%cat_diagn%runoff,                    &
          clsmf25_struc(nid)%cat_diagn%eint, &
          clsmf25_struc(nid)%cat_diagn%esoi,   &
          clsmf25_struc(nid)%cat_diagn%eveg, &
          clsmf25_struc(nid)%cat_diagn%esno,   &
          clsmf25_struc(nid)%cat_diagn%bflow, &
          clsmf25_struc(nid)%cat_diagn%runsrf,&
          clsmf25_struc(nid)%cat_diagn%snmelt,  &
          clsmf25_struc(nid)%cat_diagn%lwup, &
          clsmf25_struc(nid)%cat_diagn%swland, &
          clsmf25_struc(nid)%cat_diagn%lhflux, &
          clsmf25_struc(nid)%cat_diagn%qinfil,&
          clsmf25_struc(nid)%cat_diagn%ar1, &
          clsmf25_struc(nid)%cat_diagn%ar2,     &
          clsmf25_struc(nid)%cat_diagn%rzeq,    &
          clsmf25_struc(nid)%cat_diagn%ghflux,     &
          clsmf25_struc(nid)%cat_diagn%tpsn(1), &
          clsmf25_struc(nid)%cat_diagn%asnow,&
          clsmf25_struc(nid)%cat_diagn%tp(1), &
          clsmf25_struc(nid)%cat_diagn%tp(2), &
          clsmf25_struc(nid)%cat_diagn%tp(3), &
          clsmf25_struc(nid)%cat_diagn%tp(4), &
          clsmf25_struc(nid)%cat_diagn%tp(5), &
          clsmf25_struc(nid)%cat_diagn%tp(6), &
          clsmf25_struc(nid)%cat_diagn%sfmc, &
          clsmf25_struc(nid)%cat_diagn%rzmc,   &
          clsmf25_struc(nid)%cat_diagn%prmc,   &
          clsmf25_struc(nid)%cat_diagn%entot, &
          clsmf25_struc(nid)%cat_diagn%wtot,  &
          wchange, echange, hsnacc, evacc, shacc, &
          SH_SNOW, AVET_SNOW, WAT_10CM, TOTWAT_SOIL, TOTICE_SOIL,         &
          lhacc,TC1_0, TC2_0, TC4_0, QA1_0, QA2_0, QA4_0, fices,&
          !ag(02Jan2021)
          clsmf25_struc(nid)%cat_route%rivsto,&
          clsmf25_struc(nid)%cat_route%fldsto,&
          clsmf25_struc(nid)%cat_route%fldfrc)

!    if(LIS_localPet.eq.165) then 
!       write(LIS_logunit,*) 'tp2 ',clsmf25_struc(nid)%cat_diagn(22527)%tpsn(1),&
!            clsmf25_struc(nid)%cat_progn(22527)%tc1,&
!            clsmf25_struc(nid)%cat_progn(22527)%tc2,&
!            clsmf25_struc(nid)%cat_progn(22527)%tc4
!    endif

     ! deal with AGCM fix at the end of subroutine catchment()
     !
     ! reichle+qliu, 9 Oct 2008

     !    if (GEOS5_coupled_style) then
     !       
     !       ! keep "modified" tc's and qa's, also modify fluxes
     !       
     !       clsmf25_struc(nid)%cat_diagn%evap   = clsmf25_struc(nid)%cat_diagn%evap   - evacc
     !       clsmf25_struc(nid)%cat_diagn%lhflux = clsmf25_struc(nid)%cat_diagn%lhflux - lhacc
     !       clsmf25_struc(nid)%cat_diagn%shflux = clsmf25_struc(nid)%cat_diagn%shflux - shacc
     !              
     !    else

     ! revert to "unmodified" tc's and qa's
     ! note that almost all diagnostic fluxes (including
     ! "evap", "shflux", "lhflux", "echange", "wchange")
     ! EXCEPT "evacc", "shacc", and "lhacc" are based on these 
     ! unmodified tc's and qa's

     clsmf25_struc(nid)%cat_progn%tc1 = TC1_0
     clsmf25_struc(nid)%cat_progn%tc2 = TC2_0
     clsmf25_struc(nid)%cat_progn%tc4 = TC4_0
     clsmf25_struc(nid)%cat_progn%qa1 = QA1_0
     clsmf25_struc(nid)%cat_progn%qa2 = QA2_0
     clsmf25_struc(nid)%cat_progn%qa4 = QA4_0

     !    end if


     ! convert ght and snow arrays into structures needed by LDASsa

     do k=1,N_gt

        clsmf25_struc(nid)%cat_progn%ght(k)=ghtcnt(k,:)

     end do

     do k=1,N_snow

        clsmf25_struc(nid)%cat_progn%wesn(k) = WESN(k,:)
        clsmf25_struc(nid)%cat_progn%htsn(k) = HTSN(k,:)
        clsmf25_struc(nid)%cat_progn%sndz(k) = SNDZ(k,:)

     end do


     ! diagnose snow temperatures  in layer 2:N_snow 
     ! ( tpsn1 already provided by catchment() )

     call get_tf_nd( N_cat, clsmf25_struc(nid)%cat_progn%htsn(1), &
          clsmf25_struc(nid)%cat_progn%wesn(1), &
          tmprealvec, fices(:,1) )

     do k=2,N_snow

        call get_tf_nd( N_cat, clsmf25_struc(nid)%cat_progn%htsn(k),&
             clsmf25_struc(nid)%cat_progn%wesn(k), &
             clsmf25_struc(nid)%cat_diagn%tpsn(k), fices(:,k) )

        clsmf25_struc(nid)%cat_diagn%tpsn(k) = &
             clsmf25_struc(nid)%cat_diagn%tpsn(k) + TF   ! convert to Kelvin

     end do

! Calculate "after" states for DelSurfHeat, DelColdCont, DelSoilMoist
    do n = 1,N_cat
       afterDSH(n) = ((clsmf25_struc(nid)%cat_progn(n)%tc1 *           &
                      clsmf25_struc(nid)%cat_diagn(n)%ar1) +           &
                     (clsmf25_struc(nid)%cat_progn(n)%tc2 *            &
                      clsmf25_struc(nid)%cat_diagn(n)%ar2) +           &
                     (clsmf25_struc(nid)%cat_progn(n)%tc4 *            &
                      (1.0 - clsmf25_struc(nid)%cat_diagn(n)%ar1       &
                           - clsmf25_struc(nid)%cat_diagn(n)%ar2))) * CSO(n)
       afterDCC(n) = sum(clsmf25_struc(nid)%cat_progn(n)%htsn(1:N_snow)) + &
                     ((((((clsmf25_struc(nid)%met_force(n)%tair-TF) *  &
                           SCONST) - lhf) *                            &
                           clsmf25_struc(nid)%met_force(n)%snowf) -    &
                     (clsmf25_struc(nid)%cat_diagn(n)%snmelt*lhf)) * dtstep)
       afterDSM(n) = (clsmf25_struc(nid)%cat_param(n)%dzpr *           &
                      clsmf25_struc(nid)%cat_param(n)%poros) -         &
                     clsmf25_struc(nid)%cat_progn(n)%catdef +          &
                     clsmf25_struc(nid)%cat_progn(n)%rzexc +           &
                     clsmf25_struc(nid)%cat_progn(n)%srfexc
       afterDSW(n) = sum(clsmf25_struc(nid)%cat_progn(n)%wesn(1:N_snow))
       afterDIN(n) = clsmf25_struc(nid)%cat_progn(n)%capac
! Prevent negative temperatures to stop model crashes
! during coldstart with globally constant temperatures.
! These messages should ONLY appear during the first
! few timesteps of a coldstart.  If they also appear
! at other times (or when using a restart file),
! you've got some other problems! - dmm
!       if(LIS_localPet.eq.165.and.n.eq.22527) then 
!          write(LIS_logunit,*) 'tp ',clsmf25_struc(nid)%cat_diagn(n)%tpsn(1),&
!               clsmf25_struc(nid)%cat_progn(n)%wesn(:),&
!               clsmf25_struc(nid)%cat_progn(n)%sndz(:)               
!       endif
       if(clsmf25_struc(nid)%cat_diagn(n)%tpsn(1).lt.200) then 
          write(LIS_logunit,*) '[WARN] tp problem ',n,&
               clsmf25_struc(nid)%cat_diagn(n)%tpsn(1)
          
       endif
       if (clsmf25_struc(nid)%cat_progn(n)%tc1.le.0.0) then
          write (LIS_logunit,*) 'Resetting negative TC1 at tile ',n
          clsmf25_struc(nid)%cat_progn(n)%tc1 = 1.0
       endif
       if (clsmf25_struc(nid)%cat_progn(n)%tc2.le.0.0) then
          write (LIS_logunit,*) 'Resetting negative TC2 at tile ',n
          clsmf25_struc(nid)%cat_progn(n)%tc2 = 1.0
       endif
       if (clsmf25_struc(nid)%cat_progn(n)%tc4.le.0.0) then
          write (LIS_logunit,*) 'Resetting negative TC4 at tile ',n
          clsmf25_struc(nid)%cat_progn(n)%tc4 = 1.0
       endif
    enddo

! Gather output for LIS
     do n = 1,N_cat

! CLSM-F2.5 below follows the traditional-positive direction:
! http://www.lmd.jussieu.fr/~polcher/ALMA/convention_output_3.html

! ALMA General Energy Balance Components
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SWNET,value=   &
                        (clsmf25_struc(nid)%met_force(n)%SWdown -      &
                         clsmf25_struc(nid)%cat_diagn(n)%swup),        &
                        vlevel=1,unit="W m-2",direction="DN",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_LWNET,value=   &
                        (clsmf25_struc(nid)%met_force(n)%LWdown -      &
                         clsmf25_struc(nid)%cat_diagn(n)%lwup),        &
                        vlevel=1,unit="W m-2",direction="DN",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QLE,value=     &
                        (clsmf25_struc(nid)%cat_diagn(n)%lhflux -      &
                         clsmf25_struc(nid)%cat_diagn(n)%esno),        &
                        vlevel=1,unit="W m-2",direction="UP",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QH,value=      &
                        clsmf25_struc(nid)%cat_diagn(n)%shflux,        &
                        vlevel=1,unit="W m-2",direction="UP",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QG,value=      &
                        clsmf25_struc(nid)%cat_diagn(n)%ghflux,        &
                        vlevel=1,unit="W m-2",direction="DN",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QF,value=      &
                        (clsmf25_struc(nid)%cat_diagn(n)%snmelt*lhf),  &
                        vlevel=1,unit="W m-2",direction="S2L",          &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QV,value=      &
                        clsmf25_struc(nid)%cat_diagn(n)%esno,          &
                        vlevel=1,unit="W m-2",direction="S2V",          &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QA,value=      &
                        (((clsmf25_struc(nid)%met_force(n)%tair-TF) *  &
                          (SCONST) - lhf) *                            &
                        clsmf25_struc(nid)%met_force(n)%snowf),        &
                        vlevel=1,unit="W m-2",direction="DN",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_DELSURFHEAT,value=    &
                        (afterDSH(n)-startDSH(n)),                     &
                        vlevel=1,unit="J m-2",direction="INC",          &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_DELCOLDCONT,value=    &
                        (afterDCC(n)-startDCC(n)),                     &
                        vlevel=1,unit="J m-2",direction="INC",          &
                        surface_type=LIS_rc%lsm_index)
        if ((clsmf25_struc(nid)%cat_diagn(n)%lhflux -                  &
             clsmf25_struc(nid)%cat_diagn(n)%esno).ne.0) then
           bowen_ratio = clsmf25_struc(nid)%cat_diagn(n)%shflux /      &
                        (clsmf25_struc(nid)%cat_diagn(n)%lhflux -      &
                         clsmf25_struc(nid)%cat_diagn(n)%esno)
        else
           bowen_ratio = 0.0
        endif
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_BR,value=      &
                        bowen_ratio,                                   &
                        vlevel=1,unit="-",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        if (((clsmf25_struc(nid)%cat_diagn(n)%lhflux -                 &
              clsmf25_struc(nid)%cat_diagn(n)%esno) +                  &
              clsmf25_struc(nid)%cat_diagn(n)%shflux).ne.0) then
           evap_frac = (clsmf25_struc(nid)%cat_diagn(n)%lhflux -       &
                        clsmf25_struc(nid)%cat_diagn(n)%esno) /        &
                      ((clsmf25_struc(nid)%cat_diagn(n)%lhflux -       &
                        clsmf25_struc(nid)%cat_diagn(n)%esno) +        &
                        clsmf25_struc(nid)%cat_diagn(n)%shflux)
        else
           evap_frac = 1.0       ! Assume that EF varies between 0 and 1
        endif
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_EF,value=      &
                        evap_frac,                                     &
                        vlevel=1,unit="-",direction="-",               &
                        surface_type=LIS_rc%lsm_index)

! ALMA General Water Balance Components
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SNOWF,value=   &
                        clsmf25_struc(nid)%met_force(n)%snowf,         &
                        vlevel=1,unit="kg m-2 s-1",direction="DN",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_RAINF,value=   &
                        clsmf25_struc(nid)%met_force(n)%rainf,         &
                        vlevel=1,unit="kg m-2 s-1",direction="DN",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_EVAP,value=    &
                        clsmf25_struc(nid)%cat_diagn(n)%evap,          &
                        vlevel=1,unit="kg m-2 s-1",direction="UP",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QS,value=      &
                        clsmf25_struc(nid)%cat_diagn(n)%runsrf,        &
                        vlevel=1,unit="kg m-2 s-1",direction="OUT",        &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QS,value=      &
                        clsmf25_struc(nid)%cat_diagn(n)%runsrf*dtstep, &
                        vlevel=1,unit="kg m-2",direction="OUT",        &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QSB,value=     &
                        clsmf25_struc(nid)%cat_diagn(n)%bflow,         &
                        vlevel=1,unit="kg m-2 s-1",direction="OUT",        &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QSB,value=     &
                        clsmf25_struc(nid)%cat_diagn(n)%bflow*dtstep,  &
                        vlevel=1,unit="kg m-2",direction="OUT",        &
                        surface_type=LIS_rc%lsm_index)
        qsm_out(n) = 0.0
        qfz_out(n) = 0.0
        if (sum(clsmf25_struc(nid)%cat_progn(n)%wesn(1:N_snow)).gt.0.0) then
           lwc(n) = lwc(n) + ((clsmf25_struc(nid)%met_force(n)%Snowf   +   &
                               clsmf25_struc(nid)%met_force(n)%Rainf   -   &
                               clsmf25_struc(nid)%cat_diagn(n)%snmelt) *   &
                               dtstep) -                               &
                             ((clsmf25_struc(nid)%cat_progn(n)%wesn(1) *   &
                              (1.0 - fices(n,1)))                  +   &
                              (clsmf25_struc(nid)%cat_progn(n)%wesn(2) *   &
                              (1.0 - fices(n,2)))                  +   &
                              (clsmf25_struc(nid)%cat_progn(n)%wesn(3) *   &
                              (1.0 - fices(n,3))))
           if (lwc(n).lt.0.0) qsm_out(n) = abs(lwc(n)) / dtstep
           if (lwc(n).gt.0.0) qfz_out(n) =     lwc(n)  / dtstep
        endif
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QSM,value=     &
                        qsm_out(n),                                    &
                        vlevel=1,unit="kg m-2 s-1",direction="S2L",        &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QSM,value=     &
                        qsm_out(n)*dtstep,                             &
                        vlevel=1,unit="kg m-2",direction="S2L",        &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QFZ,value=     &
                        qfz_out(n),                                    &
                        vlevel=1,unit="kg m-2 s-1",direction="L2S",        &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_QST,value=     &
                        clsmf25_struc(nid)%cat_diagn(n)%snmelt,        &
                        vlevel=1,unit="kg m-2 s-1",direction="-",          &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_DELSOILMOIST,value=   &
                        (afterDSM(n)-startDSM(n)),                     &
                        vlevel=1,unit="kg m-2",direction="INC",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_DELSWE,value=  &
                        (afterDSW(n)-startDSW(n)),                     &
                        vlevel=1,unit="kg m-2",direction="INC",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_DELINTERCEPT,value=   &
                        (afterDIN(n)-startDIN(n)),                     &
                        vlevel=1,unit="kg m-2",direction="INC",         &
                        surface_type=LIS_rc%lsm_index)

! ALMA Surface State Variables
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SNOWT,value=   &
                        clsmf25_struc(nid)%cat_diagn(n)%tsurf,         &
                        vlevel=1,unit="K",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_VEGT,value=    &
                        clsmf25_struc(nid)%cat_diagn(n)%tsurf,         &
                        vlevel=1,unit="K",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_BARESOILT,value=&
                        clsmf25_struc(nid)%cat_diagn(n)%tsurf,         &
                        vlevel=1,unit="K",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_AVGSURFT,value=&
                        clsmf25_struc(nid)%cat_diagn(n)%tsurf,         &
                        vlevel=1,unit="K",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_RADT,value=    &
                        clsmf25_struc(nid)%cat_diagn(n)%tsurf,         &
                        vlevel=1,unit="K",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_ALBEDO,value=  &
                        clsmf25_struc(nid)%cat_diagn(n)%totalb,        &
                        vlevel=1,unit="-",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        totalb_out = LIS_rc%udef
        if ( clsmf25_struc(nid)%cat_diagn(n)%totalb .ne. LIS_rc%udef ) &
         totalb_out = clsmf25_struc(nid)%cat_diagn(n)%totalb*100.0
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_ALBEDO,value=  &
                        totalb_out,&
                        vlevel=1,unit="%",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SWE,value=     &
                        (clsmf25_struc(nid)%cat_progn(n)%wesn(1) +     &
                         clsmf25_struc(nid)%cat_progn(n)%wesn(2) +     &
                         clsmf25_struc(nid)%cat_progn(n)%wesn(3)),     &
                        vlevel=1,unit="kg m-2",direction="-",           &
                        surface_type=LIS_rc%lsm_index)

! ALMA SubSurface State Variables
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SOILMOIST,value=&
                        (clsmf25_struc(nid)%cat_diagn(n)%sfmc *        &
                         clsmf25_struc(nid)%cat_param(n)%dzsf),        &
                        vlevel=1,unit="kg m-2",direction="-",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SOILMOIST,value=&
                         (clsmf25_struc(nid)%cat_diagn(n)%rzmc *       &
                          clsmf25_struc(nid)%cat_param(n)%dzrz),       &
                        vlevel=2,unit="kg m-2",direction="-",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SOILMOIST,value=&
                        (clsmf25_struc(nid)%cat_diagn(n)%prmc *        &
                         clsmf25_struc(nid)%cat_param(n)%dzpr),        &
                        vlevel=3,unit="kg m-2",direction="-",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SOILMOIST,value=&
                        clsmf25_struc(nid)%cat_diagn(n)%sfmc,          &
                        vlevel=1,unit="m^3 m-3",direction="-",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SOILMOIST,value=&
                        clsmf25_struc(nid)%cat_diagn(n)%rzmc,          &
                        vlevel=2,unit="m^3 m-3",direction="-",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SOILMOIST,value=&
                        clsmf25_struc(nid)%cat_diagn(n)%prmc,          &
                        vlevel=3,unit="m^3 m-3",direction="-",           &
                        surface_type=LIS_rc%lsm_index)
        do i = 1,N_gt
           call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SOILTEMP,value=&
                        clsmf25_struc(nid)%cat_diagn(n)%tp(i)+TF,      &
                        vlevel=i,unit="K",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        enddo
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SOILWET,value= &
                        (clsmf25_struc(nid)%cat_param(n)%cdcr2   +     &
                         clsmf25_struc(nid)%cat_progn(n)%srfexc  +     &
                         clsmf25_struc(nid)%cat_progn(n)%rzexc   -     &
                         clsmf25_struc(nid)%cat_progn(n)%catdef) /     &
                         clsmf25_struc(nid)%cat_param(n)%cdcr2,        &
                        vlevel=1,unit="-",direction="-",               &
                        surface_type=LIS_rc%lsm_index)

! ALMA Evaporation Components
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_ECANOP,value=  &
                        (clsmf25_struc(nid)%cat_diagn(n)%eint/lhe),    &
                        vlevel=1,unit="kg m-2 s-1",direction="UP",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_ECANOP,value=  &
                        (clsmf25_struc(nid)%cat_diagn(n)%eint),    &
                        vlevel=1,unit="W m-2",direction="UP",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_TVEG,value=    &
                        (clsmf25_struc(nid)%cat_diagn(n)%eveg/lhe),    &
                        vlevel=1,unit="kg m-2 s-1",direction="UP",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_TVEG,value=    &
                        (clsmf25_struc(nid)%cat_diagn(n)%eveg),    &
                        vlevel=1,unit="W m-2",direction="UP",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_ESOIL,value=   &
                        (clsmf25_struc(nid)%cat_diagn(n)%esoi/lhe),    &
                        vlevel=1,unit="kg m-2 s-1",direction="UP",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_ESOIL,value=   &
                        (clsmf25_struc(nid)%cat_diagn(n)%esoi),    &
                        vlevel=1,unit="W m-2",direction="UP",         &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_ROOTMOIST,value=&
                         (clsmf25_struc(nid)%cat_diagn(n)%rzmc   *     &
                          clsmf25_struc(nid)%cat_param(n)%dzrz),       &
                        vlevel=1,unit="kg m-2",direction="-",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_CANOPINT,value=&
                        clsmf25_struc(nid)%cat_progn(n)%capac,         &
                        vlevel=1,unit="kg m-2",direction="-",           &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_EVAPSNOW,value=&
                        (clsmf25_struc(nid)%cat_diagn(n)%esno/ALHS),   &
                        vlevel=1,unit="kg m-2 s-1",direction="-",          &
                        surface_type=LIS_rc%lsm_index)
        clsmf25_struc(nid)%cat_output(n)%acond =                            &
                         (1.0 - clsmf25_struc(nid)%cat_diagn(n)%asnow) *    &
                         (clsmf25_struc(nid)%cat_diagn(n)%ar1 / RA1(n) +    &
                          clsmf25_struc(nid)%cat_diagn(n)%ar2 / RA2(n) +    &
                         (1.0 - clsmf25_struc(nid)%cat_diagn(n)%ar1 -       &
                          clsmf25_struc(nid)%cat_diagn(n)%ar2 ) / RA4(n)) + &
                          clsmf25_struc(nid)%cat_diagn(n)%asnow / RAS(n)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_ACOND,value=   &
                        clsmf25_struc(nid)%cat_output(n)%acond,        &
                        vlevel=1,unit="m s-1",direction="-",             &
                        surface_type=LIS_rc%lsm_index)

! ALMA Other hydrologic variables
        if (clsmf25_struc(nid)%cat_progn(n)%catdef.le.&
             clsmf25_struc(nid)%cat_param(n)%cdcr1) then
           clsmf25_struc(nid)%cat_output(n)%watertabled =                   &
                (sqrt(1.e-20 + clsmf25_struc(nid)%cat_progn(n)%catdef /    &
                clsmf25_struc(nid)%cat_param(n)%bf1) -      &
                clsmf25_struc(nid)%cat_param(n)%bf2)
        else
           clsmf25_struc(nid)%cat_output(n)%watertabled =                   &
                (clsmf25_struc(nid)%cat_param(n)%cdcr2 /     &
                (1.0 - clsmf25_struc(nid)%cat_param(n)%wpwet)) /   &
                clsmf25_struc(nid)%cat_param(n)%poros / 1000.
        endif
        clsmf25_struc(nid)%cat_output(n)%watertabled =                 &
             max(clsmf25_struc(nid)%cat_output(n)%watertabled,0.0)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_WATERTABLED,value=  &
             clsmf25_struc(nid)%cat_output(n)%watertabled, &
             vlevel=1,unit="m",direction="-",              &
             surface_type=LIS_rc%lsm_index)
        clsmf25_struc(nid)%cat_output(n)%tws =                              &
             (clsmf25_struc(nid)%cat_param(n)%cdcr2 /     &
             (1.0 - clsmf25_struc(nid)%cat_param(n)%wpwet)) -   &
             clsmf25_struc(nid)%cat_progn(n)%catdef  +   &
             clsmf25_struc(nid)%cat_progn(n)%srfexc  +   &
             clsmf25_struc(nid)%cat_progn(n)%rzexc   +   &
             sum(clsmf25_struc(nid)%cat_progn(n)%wesn(1:N_snow)) +   &
             clsmf25_struc(nid)%cat_progn(n)%capac
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_TWS,value=     &
             clsmf25_struc(nid)%cat_output(n)%tws,          &
             vlevel=1,unit="mm",direction="-",              &
             surface_type=LIS_rc%lsm_index)
        GWS_out =     clsmf25_struc(nid)%cat_output(n)%tws -           &
             (clsmf25_struc(nid)%cat_diagn(n)%rzmc *           &
             clsmf25_struc(nid)%cat_param(n)%dzrz) -          &
             sum(clsmf25_struc(nid)%cat_progn(n)%wesn(1:N_snow))- &
             clsmf25_struc(nid)%cat_progn(n)%capac
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_GWS,value=     &
             GWS_out,                                       &
             vlevel=1,unit="mm",direction="-",              &
             surface_type=LIS_rc%lsm_index)
        
        ! ALMA Cold Season Processes
        call calc_asnow( 1, clsmf25_struc(nid)%cat_progn(n:n),         &
             clsmf25_struc(nid)%cat_output(n:n)%snocovr)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SNOWCOVER,value=&
             clsmf25_struc(nid)%cat_output(n)%snocovr,      &
             vlevel=1,unit="-",direction="-",               &
             surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SALBEDO,value= &
             clsmf25_struc(nid)%cat_output(n)%albsn,        &
             vlevel=1,unit="-",direction="-",               &
             surface_type=LIS_rc%lsm_index)
        do i = 1,N_snow
           call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SNOWTPROF,value=&
                        clsmf25_struc(nid)%cat_diagn(n)%tpsn(i),       &
                        vlevel=i,unit="K",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        enddo
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SNOWDEPTH,value=&
                        ((clsmf25_struc(nid)%cat_progn(n)%sndz(1)   +  &
                          clsmf25_struc(nid)%cat_progn(n)%sndz(2)   +  &
                          clsmf25_struc(nid)%cat_progn(n)%sndz(3))  *  &
                          clsmf25_struc(nid)%cat_output(n)%snocovr),   &
                        vlevel=1,unit="m",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SNOWDEPTH,value=      &
                        (((clsmf25_struc(nid)%cat_progn(n)%sndz(1)*1000.0)  + &
                          (clsmf25_struc(nid)%cat_progn(n)%sndz(2)*1000.0)  + &
                          (clsmf25_struc(nid)%cat_progn(n)%sndz(3)*1000.0)) * &
                           clsmf25_struc(nid)%cat_output(n)%snocovr),         &
                        vlevel=1,unit="mm",direction="-",              &
                        surface_type=LIS_rc%lsm_index)

        do i = 1,N_snow
           bdsnow = LIS_rc%udef
           if ((clsmf25_struc(nid)%cat_progn(n)%sndz(i) *              &
                clsmf25_struc(nid)%cat_output(n)%snocovr).gt.0) then 
              bdsnow = (clsmf25_struc(nid)%cat_progn(n)%wesn(i)) /       &
                   (clsmf25_struc(nid)%cat_progn(n)%sndz(i)*    &
                    clsmf25_struc(nid)%cat_output(n)%snocovr)
           endif
           call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_LAYERSNOWDENSITY,&
                value=bdsnow,                                          &
                vlevel=i,unit="kg m-3",direction="-",                  &
                surface_type=LIS_rc%lsm_index)
        enddo


! Rhae Sung added snow depth, snow ice and snow liquid for each layer
        do i = 1,N_snow
           call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_LAYERSNOWDEPTH,value=&
                        (clsmf25_struc(nid)%cat_progn(n)%sndz(i) *     &
                         clsmf25_struc(nid)%cat_output(n)%snocovr),    &
                        vlevel=i,unit="m",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        enddo

        do i = 1,N_snow
           call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SNOWICE,value=&
                        (clsmf25_struc(nid)%cat_progn(n)%sndz(i) *     &
                         clsmf25_struc(nid)%cat_output(n)%snocovr *    &
                         fices(n,i)*1000.0),                           &
                        vlevel=i,unit="mm",direction="-",              &
                        surface_type=LIS_rc%lsm_index)
        enddo

        do i = 1,N_snow
           call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_SNOWLIQ,value=&
                        (clsmf25_struc(nid)%cat_progn(n)%sndz(i) *     &
                         clsmf25_struc(nid)%cat_output(n)%snocovr *    &
                         (1.0-fices(n,i))*1000.0),                     &
                        vlevel=i,unit="mm",direction="-",              &
                        surface_type=LIS_rc%lsm_index)
        enddo    

! ALMA variables to be compared with remote sensed data
       call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_LWUP,value=     &
                        clsmf25_struc(nid)%cat_diagn(n)%lwup,          &
                        vlevel=1,unit="W m-2",direction="UP",           &
                        surface_type=LIS_rc%lsm_index)

! Parameters
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_LANDCOVER,value=&
                        (clsmf25_struc(nid)%cat_param(n)%vegcls*1.0),  &
                        vlevel=1,unit="-",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_POROSITY,value=&
                        clsmf25_struc(nid)%cat_param(n)%poros,         &
                        vlevel=1,unit="-",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_GREENNESS,value=&
                        green(n),                                      &
                        vlevel=1,unit="-",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_GREENNESS,value=&
                        green(n)*100.0,                                &
                        vlevel=1,unit="%",direction="-",               &
                        surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_LAI,value=     &
                        lai(n),                                        &
                        vlevel=1,unit="-",direction="-",               &
                        surface_type=LIS_rc%lsm_index)

! ALMA Other/Miscellaneous Variables

! Forcing variables
        call LIS_diagnoseSurfaceOutputVar(nid,n,LIS_MOC_FHEIGHTFORC,value=       &
                        clsmf25_struc(nid)%met_force(n)%refh,          &
                        vlevel=1,unit="m",direction="-",               &
                        surface_type=LIS_rc%lsm_index)

! Reset forcing variables
        clsmf25_struc(nid)%met_force(n)%tair = 0
        clsmf25_struc(nid)%met_force(n)%qair = 0
        clsmf25_struc(nid)%met_force(n)%swdown = 0
        clsmf25_struc(nid)%met_force(n)%lwdown = 0
        clsmf25_struc(nid)%met_force(n)%wind = 0
        clsmf25_struc(nid)%met_force(n)%psurf = 0
        clsmf25_struc(nid)%met_force(n)%snowf = 0
        clsmf25_struc(nid)%met_force(n)%rainf = 0
        clsmf25_struc(nid)%met_force(n)%rainf_c = 0
        clsmf25_struc(nid)%met_force(n)%refh = 0
        clsmf25_struc(nid)%met_force(n)%swnet = 0
        clsmf25_struc(nid)%met_force(n)%pardr = 0
        clsmf25_struc(nid)%met_force(n)%pardf = 0

#if 0 
        if(LIS_localPet.eq.147.and.n.eq.1882) then 
           print*, 'inter',LIS_rc%yr, LIS_rc%mo, LIS_rc%da,&
                LIS_rc%hr, LIS_rc%mn, LIS_rc%ss, &
                clsmf25_struc(nid)%cat_progn(n)%catdef,&
                clsmf25_struc(nid)%cat_progn(n)%rzexc, &
                clsmf25_struc(nid)%cat_progn(n)%srfexc, &
                clsmf25_struc(nid)%cat_progn(n)%sndz, &
                clsmf25_struc(nid)%cat_progn(n)%wesn, &
                clsmf25_struc(nid)%cat_diagn(n)%lhflux, & 
                clsmf25_struc(nid)%cat_diagn(n)%shflux
           print*,''
        endif
!checks:
        if((clsmf25_struc(nid)%cat_diagn(n)%tsurf.lt.200.or.&
             clsmf25_struc(nid)%cat_diagn(n)%tsurf.gt.500).or.&
             abs(clsmf25_struc(nid)%cat_diagn(n)%lhflux).gt.2000.or.&
             abs( clsmf25_struc(nid)%cat_diagn(n)%shflux).gt.2000) then 
           print*, 'model unstable ',LIS_rc%yr, LIS_rc%mo, LIS_rc%da,&
                LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
           print*, LIS_localPet, n, &
                clsmf25_struc(nid)%cat_progn(n)%catdef,&
                clsmf25_struc(nid)%cat_progn(n)%rzexc, &
                clsmf25_struc(nid)%cat_progn(n)%srfexc, &
                clsmf25_struc(nid)%cat_progn(n)%sndz, &
                clsmf25_struc(nid)%cat_progn(n)%wesn, &
                clsmf25_struc(nid)%cat_diagn(n)%lhflux, & 
                clsmf25_struc(nid)%cat_diagn(n)%shflux
           stop
        endif
#endif
     enddo
!     write(LIS_logunit,*) 'main ',clsmf25_struc(nid)%cat_progn(327)%catdef
     clsmf25_struc(nid)%forc_count = 0
  endif

end subroutine clsmf25_main

