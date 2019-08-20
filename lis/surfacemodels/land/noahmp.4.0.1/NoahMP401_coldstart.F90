!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0     
!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: NoahMP401_coldstart
! \label{NoahMP401_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/25/18: Shugong Wang, Zhuo Wang; initial implementation for LIS 7 and NoahMP401
!
! !INTERFACE:
subroutine NoahMP401_coldstart(mtype)
! !USES:
!   use LIS_coreMod, only: LIS_rc
!   use LIS_logMod, only: LIS_logunit
!   use LIS_timeMgrMod, only: LIS_date2time

    use LIS_coreMod
    use LIS_logMod
    use LIS_histDataMod

    use LIS_timeMgrMod, only: LIS_date2time
    use NoahMP401_lsmMod
    use module_sf_noahmpdrv_401

!
! !DESCRIPTION:
!
!  This routine initializes the NoahMP401 state variables with
!  some predefined values constantly for the entire domain.
!
!EOP
 
    implicit none

    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r

    ! Added by Zhuo Wang on 11/21/2018
    integer :: ids,ide, jds,jde, kds,kde,  &
         &     ims,ime, jms,jme, kms,kme,  &
         &     its,ite, jts,jte, kts,kte

    !------------------------------------------------------------------------
    ! End 2D variables not used in WRF
    !------------------------------------------------------------------------
    CHARACTER(LEN=256) :: LLANDUSE          ! (=USGS, using USGS landuse classification)
    !------------------------------------------------------------------------

    integer ::     row, col
    integer ::     NSOIL    ! number of soil layers
    LOGICAL ::     restart,          &
         &         allowed_to_read

    real, dimension(4) ::     DZS  ! Thickness of the soil layers [m]
    real ::     dx, dy
    real, dimension( 1, 1 ) :: msftx, msfty
    real :: wtddt, dtbl
    integer :: stepwtd

    real, dimension( 1, 1 ) ::                      &
         &                         snowxy,          &  ! snow water equivalent [mm]
         &                         snowhxy,         &  ! physical snow depth [m]
         &                         canwatxy            ! total canopy water + ice [mm]
 
    integer, dimension( 1, 1 ) ::                   &
         &                           ISLTYP,        &  ! soil type
                                     IVGTYP            ! vegetation type

    real,    dimension( 1, 4, 1 ) ::            &
         &                            tslb_3d,      &  ! soil temperature [K]
         &                           smois_3d,      &  ! volumetric soil moisture [m3/m3]
         &                            sh2o_3d          ! volumetric liquid soil moisture [m3/m3]
 
    LOGICAL                   ::     FNDSOILW,      &  ! soil water present in input
         &                           FNDSNOWH          ! snow depth present in input

    integer                   :: runoff_option
    integer                   :: crop_option
  
    real, dimension(1,1)      :: XLAT        !latitude
    real, dimension(1,1)      :: TSK         !skin temperature (k)
    real, dimension(1,1)      :: tvxy        !vegetation canopy temperature
    real, dimension(1,1)      :: tgxy        !ground surface temperature
    real, dimension(1,1)      :: canicexy    !canopy-intercepted ice (mm)
    real, dimension(1,1)      :: canliqxy    !canopy-intercepted liquid water (mm)
    real, dimension(1,1)      :: tmnxy       !deep soil temperature (k)
    real, dimension(1,1)      :: XICE        !sea ice fraction
    real, dimension(1,1)      :: eahxy       !canopy air vapor pressure (pa)
    real, dimension(1,1)      :: tahxy       !canopy air temperature (k)
    real, dimension(1,1)      :: cmxy        !momentum drag coefficient
    real, dimension(1,1)      :: chxy        !sensible heat exchange coefficient
    real, dimension(1,1)      :: fwetxy      !wetted or snowed fraction of the canopy (-)
    real, dimension(1,1)      :: sneqvoxy    !snow mass at last time step(mm h2o)
    real, dimension(1,1)      :: alboldxy    !snow albedo at last time step (-)
    real, dimension(1,1)      :: qsnowxy     !snowfall on the ground [mm/s]
    real, dimension(1,1)      :: wslakexy    !lake water storage [mm]
    real, dimension(1,1)      :: zwtxy       !water table depth [m]
    real, dimension(1,1)      :: waxy        !water in the "aquifer" [mm]
    real, dimension(1,1)      :: wtxy        !groundwater storage [mm]
    real, allocatable, dimension(:,:,:) :: tsnoxy  !snow temperature [K]
    real, allocatable, dimension(:,:,:) :: zsnsoxy !snow layer depth [m]
    real, allocatable, dimension(:,:,:) :: snicexy !snow layer ice [mm]
    real, allocatable, dimension(:,:,:) :: snliqxy !snow layer liquid water [mm]
    real, dimension(1,1)      :: lfmassxy    !leaf mass [g/m2]
    real, dimension(1,1)      :: rtmassxy    !mass of fine roots [g/m2]
    real, dimension(1,1)      :: stmassxy    !stem mass [g/m2]
    real, dimension(1,1)      :: woodxy      !mass of wood (incl. woody roots) [g/m2]
    real, dimension(1,1)      :: stblcpxy    !stable carbon in deep soil [g/m2]
    real, dimension(1,1)      :: fastcpxy    !short-lived carbon, shallow soil [g/m2]
    real, dimension(1,1)      :: saixy       !stem area index
    real, dimension(1,1)      :: laixy       !leaf area index
    real, dimension(1,1)      :: grainxy     !mass of grain [g/m2] !XING
    real, dimension(1,1)      :: gddxy       !growing degree days !XING
    integer, dimension(1,  1) :: cropcatxy

    real, dimension(1,1) :: t2mvxy        !2m temperature vegetation part (k)
    real, dimension(1,1) :: t2mbxy        !2m temperature bare ground part (k)
    real, dimension(1,1) :: chstarxy      !dummy

    ! Optional
    real, dimension(1,1) :: smcwtdxy    !deep soil moisture content [m3 m-3]
    real, dimension(1,1) :: deeprechxy  !deep recharge [m]
    real, dimension(1,1) :: rechxy      !accumulated recharge [mm]
    real, dimension(1,1) :: qrfsxy      !accumulated flux from groundwater to rivers [mm]
    real, dimension(1,1) :: qspringsxy  !accumulated seeping water [mm]
    real, dimension(1,1) :: qslatxy     !accumulated lateral flow [mm]
    real, dimension(1,1) :: areaxy      !grid cell area [m2]
    real, dimension(1,1) :: FDEPTHXY    !efolding depth for transmissivity (m)
    real, dimension(1,1) :: HT          !terrain height (m)
    real, dimension(1,1) :: RIVERBEDXY  !riverbed depth (m)
    real, dimension(1,1) :: EQZWT       !equilibrium water table depth (m)
    real, dimension(1,1) :: RIVERCONDXY !river conductance
    real, dimension(1,1) :: PEXPXY      !factor for river conductance
    real, dimension(1,1) :: rechclim
    integer            :: sf_urban_physics

   real, dimension(1,1)    :: taussxy
   real, dimension(1,1)    :: acsnomxy
   real, dimension(1,1)    :: acsnowxy
   real, dimension(1,1)    :: sfcrunoffxy
   real, dimension(1,1)    :: udrrunoffxy
   real, allocatable, dimension(:,:,:)  :: smoiseqxy
 
   integer, dimension(1,1) :: pgsxy
   real, allocatable, dimension(:,:,:) :: croptype
!  integer, dimension(1,1) :: PLANTING,HARVEST,SEASON_GDD
!  integer, dimension(1,1) :: SLOPETYP

   ! added by Shugong Wang
   integer :: isnow
   real, allocatable, dimension(:) :: zsnso
   real, allocatable, dimension(:) :: snice
   real, allocatable, dimension(:) :: snliq
   real, allocatable, dimension(:) :: zsoil

!  Modified by Zhuo Wang
   real, allocatable, dimension(:) :: tsnow 
   ! end add

   !EMK...Temporary arrays.
   real :: tmp_swe(1,1),tmp_snodep(1,1)
   integer :: isnowxy(1,1)
   integer :: tmp_isnowxy(1,1)

   real, dimension(1,60,1) :: gecros_state  ! Optional gecros crop
  
!  PLANTING(1,1)   = 126     ! default planting date
!  HARVEST(1,1)    = 290     ! default harvest date
!  SEASON_GDD(1,1) = 1605    ! default total seasonal growing degree days
!  SLOPETYP   = 2
!---------------------------------------------------------------------
    ids = 1
    ide = 2
    jds = 1
    jde = 2
    kds = 1
    kde = 1
    ims = 1
    ime = 1
    jms = 1
    jme = 1
    kms = 1
    kme = 1
    its = 1
    ite = 1
    jts = 1
    jte = 1
    kts = 1
    kte = 1

    do n=1, LIS_rc%nnest
        isnow = -NOAHMP401_struc(n)%nsnow

        ! Added by Shugong
        allocate(zsnso(-NOAHMP401_struc(n)%nsnow+1:NOAHMP401_struc(n)%nsoil))
        ! Modified by Zhuo Wang
        allocate(tsnow(-NOAHMP401_struc(n)%nsnow+1:0))
        allocate(snice(-NOAHMP401_struc(n)%nsnow+1:0))
        allocate(snliq(-NOAHMP401_struc(n)%nsnow+1:0))
        allocate(zsoil(NOAHMP401_struc(n)%nsoil))
        ! Added by David
        allocate(tsnoxy(1,-NOAHMP401_struc(n)%nsnow+1:0,1))
        allocate(zsnsoxy(1,-NOAHMP401_struc(n)%nsnow+1:NOAHMP401_struc(n)%nsoil,1))
        allocate(snicexy(1,-NOAHMP401_struc(n)%nsnow+1:0,1))
        allocate(snliqxy(1,-NOAHMP401_struc(n)%nsnow+1:0,1))
        allocate(smoiseqxy(1,NOAHMP401_struc(n)%nsoil,1))
        allocate(croptype(1,5,1)) ! Check to see if 5 should be hard-coded

        zsoil(1) = -NOAHMP401_struc(n)%sldpth(1)
        do l=2, NOAHMP401_struc(n)%nsoil
          zsoil(l) = zsoil(l-1) - NOAHMP401_struc(n)%sldpth(l) 
        enddo
        ! end add 

!------------------------------------
! Initialize Model State Variables
!----------------------------------------------------------------------
        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) &
             "[INFO] NoahMP401_coldstart -- cold-starting Noah-MP.4.0.1"

           do t=1, LIS_rc%npatch(n,mtype)

             row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
             col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
             XLAT(1,1) = LIS_domain(n)%grid(LIS_domain(n)%gindex(col,row))%lat

            ! added by shugong 
            zsnso = 0.0
            tmp_swe(1,1) = NOAHMP401_struc(n)%init_sneqv
            tgxy(1,1) = NOAHMP401_struc(n)%init_tskin
            tmp_snodep(1,1) = NOAHMP401_struc(n)%init_snowh

            call snow_init_401(1, 1, 1, 1, 1, 1, 1, 1,           & !input 
                              NOAHMP401_struc(n)%nsnow,          & !input 
                              NOAHMP401_struc(n)%nsoil,          & !input 
                              zsoil,                             & !input
                              tmp_swe,                           & !input
                              tgxy,                              & !input
                              tmp_snodep,                        & !input
                              zsnso, tsnow, snice, snliq,        & !output
                              isnowxy) ! output
             isnow = isnowxy(1,1)
             ! end add

             snicexy(1,-NOAHMP401_struc(n)%nsnow+1:0,1) = snice(-NOAHMP401_struc(n)%nsnow+1:0)
             snliqxy(1,-NOAHMP401_struc(n)%nsnow+1:0,1) = snliq(-NOAHMP401_struc(n)%nsnow+1:0)
             zsnsoxy(1,-NOAHMP401_struc(n)%nsnow+1:NOAHMP401_struc(n)%nsoil,1) = zsnso(-NOAHMP401_struc(n)%nsnow+1:NOAHMP401_struc(n)%nsoil) 
             tsnoxy(1,-NOAHMP401_struc(n)%nsnow+1:0,1) = tsnow(-NOAHMP401_struc(n)%nsnow+1:0)

             do l=1,5
               croptype(1,l,1) = 0.0
             enddo

             cropcatxy(1,1) = LIS_rc%cropclass

             if (NOAHMP401_struc(n)%crop_opt.eq.2) then
                do l=1, 60      ! TODO: check loop
                   gecros_state(1,l,1) = NOAHMP401_struc(n)%init_gecros_state(l)
                enddo
             else
                gecros_state(1,l,1) = 0.0
             endif

             canwatxy(1,1) = NOAHMP401_struc(n)%init_canwat
             IVGTYP(1,1) = NOAHMP401_struc(n)%noahmp401(t)%vegetype
             ISLTYP(1,1) = NOAHMP401_struc(n)%noahmp401(t)%soiltype

             do l=1, NOAHMP401_struc(n)%nsoil
               tslb_3d(1,l,1) = NOAHMP401_struc(n)%init_tslb(l)  ! Noah-MP.4.0.1 initial soil temperatures: 
               smois_3d(1,l,1) = NOAHMP401_struc(n)%init_smc(l)  ! Noah-MP.4.0.1 initial total soil moistures:
               sh2o_3d(1,l,1) = 0.0 ! Noah-MP.4.0.1 initial liquid soil moistures set in NOAHMP_INIT
               DZS(l) = NOAHMP401_struc(n)%sldpth(l)    ! Noah-MP.4.0.1 thickness of soil layers:
             enddo

             FNDSOILW = .true.
             FNDSNOWH = .true.

             TSK(1,1) = NOAHMP401_struc(n)%init_tskin       ! Noah-MP.4.0.1 surface skin temperature:
             tvxy(1,1) = 0.0 ! Noah-MP.4.0.1 initial vegetation temperature set in NOAHMP_INIT
             tgxy(1,1) = 0.0 ! Noah-MP.4.0.1 initial ground temperature set in NOAHMP_INIT
             canicexy(1,1) = 0.0  ! Noah-MP.4.0.1 initial canopy-intercepted ice set in NOAHMP_INIT
             canliqxy(1,1) = 0.0  ! Noah-MP.4.0.1 initial canopy-intercepted liquid water set in NOAHMP_INIT
             tmnxy(1,1) = NOAHMP401_struc(n)%noahmp401(t)%tbot ! Noah-MP.4.0.1 soil temperature lower boundary:
! Fix later (should it be zero?)
!             XICE(1,1) = NOAHMP401_struc(n)%seaice   ! Noah-MP.4.0.1 sea ice fraction:
             XICE(1,1) = 0.0
             eahxy(1,1) = 0.0
             tahxy(1,1) = 0.0
             cmxy(1,1) = 0.0   
             chxy(1,1) = 0.0   
             fwetxy(1,1) = 0.0
             sneqvoxy(1,1) = 0.0 ! Noah-MP.4.0.1 initial snow mass at last time step set in NOAHMP_INIT
             alboldxy(1,1) = 0.0 ! Noah-MP.4.0.1 initial snow albedo at last time step set in NOAHMP_INIT
             qsnowxy(1,1) =  0.0 ! Noah-MP.4.0.1 initial snowfall on the ground set in NOAHMP_INIT
             wslakexy(1,1) = 0.0 
             zwtxy(1,1) = NOAHMP401_struc(n)%init_zwt         ! Noah-MP.4.0.1 initial water table depth:
             waxy(1,1) = NOAHMP401_struc(n)%init_wa           ! Noah-MP.4.0.1 initial water in the aquifer:
             wtxy(1,1) = NOAHMP401_struc(n)%init_wt           ! Noah-MP.4.0.1 initial water in aquifer and saturated soil:
             lfmassxy(1,1) = 0.0   ! Noah-MP.4.0.1 initial leaf mass set in NOAHMP_INIT
             rtmassxy(1,1) = 0.0   ! Noah-MP.4.0.1 initial mass of fine roots set in NOAHMP_INIT
             stmassxy(1,1) = 0.0   ! Noah-MP.4.0.1 initial stem mass set in NOAHMP_INIT
             woodxy(1,1)   = 0.0   ! Noah-MP.4.0.1 initial mass of wood including woody roots set in NOAHMP_INIT 
             stblcpxy(1,1) = 0.0   ! Noah-MP.4.0.1 initial stable carbon in deep soil set in NOAHMP_INIT
             fastcpxy(1,1) = 0.0   ! Noah-MP.4.0.1 initial short-lived carbon in shallow soil set in NOAHMP_INIT
             saixy(1,1)    = 0.0   ! Noah-MP.4.0.1 initial stem area index set in NOAHMP_INIT
             laixy(1,1) = NOAHMP401_struc(n)%init_lai         ! Noah-MP.4.0.1 initial leaf area index:
             grainxy(1,1) = 0.0
             gddxy(1,1) = 0.0
             cropcatxy(1,1) = 0.0
             t2mvxy(1,1) = 0.0 ! set in NOAHMP_INIT
             t2mbxy(1,1) = 0.0 ! set in NOAHMP_INIT
             chstarxy(1,1) = 0.0   ! dummy
             NSOIL = NOAHMP401_struc(n)%nsoil 
             restart = .false. 
             allowed_to_read = .true. 
             runoff_option = NOAHMP401_struc(n)%run_opt       ! Noah-MP.4.0.1 runoff and groundwater option:
             crop_option = NOAHMP401_struc(n)%crop_opt        ! Noah-MP.4.0.1 crop model option:
             sf_urban_physics = NOAHMP401_struc(n)%urban_opt  ! Noah-MP.4.0.1 urban physics option:

             ! The following variables are optional for groundwater dynamics iopt_run=5, 
             !  so some random values are set temporarily.

             do l=1, NOAHMP401_struc(n)%nsoil
               smoiseqxy(1,l,1) = 0.0  ! Noah-MP.4.0.1 initial equilibrium soil moisture content
             end do

             smcwtdxy(1,1) = 0.0
             rechxy(1,1) = 0.0
             deeprechxy(1,1) = 0.0
             areaxy(1,1) = 100.
             dx = 10.0
             dy = 10.0
             msftx(1,1) = 0.0
             msfty(1,1) = 0.0
             wtddt = 0.0
             stepwtd = 0
             dtbl = 0.0
             qrfsxy(1,1) = 0.0
             qslatxy(1,1) = 0.0
             fdepthxy(1,1) = 0.0
! Fix terrain height later for crop model (Get from lis_input file)
             HT(1,1) = -9999.9
             riverbedxy(1,1) = 0.0
             eqzwt(1,1) = 0.0
             rivercondxy(1,1) = 0.0
             pexpxy(1,1) = 0.0
             rechclim(1,1) = 0.0
             sfcrunoffxy(1,1) = 0.0
             udrrunoffxy(1,1) = 0.0
             acsnomxy(1,1) = 0.0
             acsnowxy(1,1) = 0.0
             taussxy(1,1) = 0.0
             pgsxy(1,1) = 0

            CALL NOAHMP_INIT(    LLANDUSE, tmp_swe, tmp_snodep, canwatxy,   ISLTYP,   IVGTYP, XLAT, & 
                           tslb_3d, smois_3d,     sh2o_3d,      DZS, FNDSOILW, FNDSNOWH, &
                          TSK,  isnowxy,   tvxy, tgxy, canicexy,      tmnxy,     XICE, &
                       canliqxy,    eahxy,    tahxy,     cmxy,     chxy,                     &
                         fwetxy, sneqvoxy, alboldxy,  qsnowxy, wslakexy,    zwtxy,     waxy, &
                           wtxy,   tsnoxy,  zsnsoxy,  snicexy,  snliqxy, lfmassxy, rtmassxy, &
                       stmassxy,   woodxy, stblcpxy, fastcpxy,    saixy,    laixy,           &
                        grainxy,    gddxy,                                                   &
                       croptype,  cropcatxy,                                                 &
                         t2mvxy,   t2mbxy, chstarxy,                                         &
                          NSOIL,  restart,                                                   &
                         allowed_to_read,runoff_option, crop_option,                                 &
                         sf_urban_physics,                         &  ! urban scheme
                         ids,ide, jds,jde, kds,kde,                &  ! domain
                         ims,ime, jms,jme, kms,kme,                &  ! memory
                         its,ite, jts,jte, kts,kte                 &  ! tile
                            ,smoiseqxy,smcwtdxy ,rechxy   ,deeprechxy, areaxy ,dx, dy, msftx, msfty,&
                            wtddt    ,stepwtd  ,dtbl  ,qrfsxy ,qspringsxy  ,qslatxy,                  &
                            fdepthxy ,HT       ,riverbedxy ,eqzwt ,rivercondxy ,pexpxy,              &
                            rechclim ,gecros_state                 &
                            )

                NOAHMP401_struc(n)%noahmp401(t)%sfcrunoff = sfcrunoffxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%udrrunoff = udrrunoffxy(1,1)

                do l=1, NOAHMP401_struc(n)%nsoil
                    NOAHMP401_struc(n)%noahmp401(t)%smc(l) = smois_3d(1,l,1) 
                enddo
                do l=1, NOAHMP401_struc(n)%nsoil
                    NOAHMP401_struc(n)%noahmp401(t)%sh2o(l) = sh2o_3d(1,l,1)
                enddo
                do l=1, NOAHMP401_struc(n)%nsoil
                    NOAHMP401_struc(n)%noahmp401(t)%tslb(l) = tslb_3d(1,l,1)
                enddo
                NOAHMP401_struc(n)%noahmp401(t)%sneqv = tmp_swe(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%snowh = tmp_snodep(1,1) 
                NOAHMP401_struc(n)%noahmp401(t)%canwat = canwatxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%acsnom = acsnomxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%acsnow = acsnowxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%isnow = isnowxy(1,1) 
                NOAHMP401_struc(n)%noahmp401(t)%tv = tvxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%tg = tgxy(1,1) 
                NOAHMP401_struc(n)%noahmp401(t)%canice = canicexy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%canliq = canliqxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%fwet = fwetxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%sneqvo = sneqvoxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%albold = alboldxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%qsnow = qsnowxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%wslake = wslakexy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%zwt = zwtxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%wa = waxy(1,1) 
                NOAHMP401_struc(n)%noahmp401(t)%wt = wtxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%lfmass = lfmassxy(1,1) 
                NOAHMP401_struc(n)%noahmp401(t)%rtmass = rtmassxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%stmass = stmassxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%wood = woodxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%stblcp = stblcpxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%fastcp = fastcpxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%lai = laixy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%sai = saixy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%tauss = taussxy(1,1)

                do l=1, NOAHMP401_struc(n)%nsoil
                    NOAHMP401_struc(n)%noahmp401(t)%smoiseq(l) = smoiseqxy(1,l,1) 
                enddo
                NOAHMP401_struc(n)%noahmp401(t)%smcwtd = smcwtdxy(1,1) 
                NOAHMP401_struc(n)%noahmp401(t)%deeprech = deeprechxy(1,1) 
                NOAHMP401_struc(n)%noahmp401(t)%rech = rechxy(1,1) 
                NOAHMP401_struc(n)%noahmp401(t)%grain = grainxy(1,1) 
                NOAHMP401_struc(n)%noahmp401(t)%gdd = gddxy(1,1)
                NOAHMP401_struc(n)%noahmp401(t)%pgs = pgsxy(1,1)

                do l=1, 60 ! TODO: check loop
                   NOAHMP401_struc(n)%noahmp401(t)%gecros_state(l) = gecros_state(1,l,1)
                enddo

               NOAHMP401_struc(n)%noahmp401(t)%snowice(1:NOAHMP401_struc(n)%nsnow) = snicexy(1,-NOAHMP401_struc(n)%nsnow+1:0,1)
               NOAHMP401_struc(n)%noahmp401(t)%snowliq(1:NOAHMP401_struc(n)%nsnow) = snliqxy(1,-NOAHMP401_struc(n)%nsnow+1:0,1)
               NOAHMP401_struc(n)%noahmp401(t)%zss(1:NOAHMP401_struc(n)%nsnow+NOAHMP401_struc(n)%nsoil) = zsnsoxy(1,-NOAHMP401_struc(n)%nsnow+1:NOAHMP401_struc(n)%nsoil,1) 
               NOAHMP401_struc(n)%noahmp401(t)%isnow = isnowxy(1,1) 
               NOAHMP401_struc(n)%noahmp401(t)%tsno(1:NOAHMP401_struc(n)%nsnow) = tsnoxy(1,-NOAHMP401_struc(n)%nsnow+1:0,1)
!-----------------------------------------------------------------------
            enddo   ! t=1,1

        endif       ! coldstart

        deallocate(tsnoxy)
        deallocate(zsnsoxy)
        deallocate(snicexy)
        deallocate(snliqxy)
        deallocate(smoiseqxy)
        deallocate(croptype)
    
        LIS_rc%yr = LIS_rc%syr
        LIS_rc%mo = LIS_rc%smo
        LIS_rc%da = LIS_rc%sda
        LIS_rc%hr = LIS_rc%shr
        LIS_rc%mn = LIS_rc%smn
        LIS_rc%ss = LIS_rc%sss
        
        call LIS_date2time(LIS_rc%time, LIS_rc%doy, LIS_rc%gmt, LIS_rc%yr,      &
                           LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
        write(LIS_logunit,*) "[INFO] NoahMP401_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
       enddo        ! nnest
end subroutine NoahMP401_coldstart
