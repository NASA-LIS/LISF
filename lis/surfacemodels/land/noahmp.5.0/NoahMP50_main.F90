!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!BOP
!
! !ROUTINE: NoahMP50_main
! \label{NoahMP50_main}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit
!  developed by Shugong Wang for the NASA Land Information System V7.
!  The initial specification of the subroutine is by Sujay Kumar.
!
!   10/25/18: Shugong Wang, Zhuo Wang; initial implementation for NoahMP50 with LIS-7
!   05/15/19: Yeosang Yoon; code added for snow DA to work
!   10/29/19: David Mocko; Added RELSMC to output, and an option
!                          for different units for Qs/Qsb/Albedo
!   03/09/22: David Mocko: Fixed "input LAI" for dynamic vegetation options 7/8/9
!   05/23/23: Cenlin He: modified for refactored NoahMP v5 and later

! !INTERFACE:
subroutine NoahMP50_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod,    only : LIS_isAlarmRinging
    use LIS_constantsMod,  only : LIS_CONST_RHOFW   !New
    use LIS_vegDataMod,    only : LIS_lai, LIS_sai
    use LIS_logMod,        only : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod
    use NoahMP50_lsmMod
    use NoahmpIOVarType

    implicit none
! !ARGUMENTS:
    integer, intent(in)  :: n
    integer              :: t
    integer              :: i
    real                 :: dt
    real                 :: lat, lon
    real                 :: tempval
    integer              :: row, col, tid
    logical              :: alarmCheck

!
! !DESCRIPTION:
!  This is the entry point for calling the NoahMP50 physics.
!  This routine calls the {\tt noahmp_driver_50} routine that performs
!  the land surface computations, to solve water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for NoahMP50

    ! Code added by Zhuo Wang 02/28/2019
    real                 :: AvgSurfT_out           ! average surface temperature [K]
    real                 :: TWS_out                ! terrestrial water storage [mm]
    ! Code added by David Mocko 04/25/2019
    real                 :: startsm, startswe, startint, startgw, endsm
   
    ! EMK for 557WW
    real :: tmp_q2sat, tmp_es
    character*3 :: fnest
    REAL, PARAMETER:: LVH2O = 2.501000E+6 ! Latent heat for evapo for water  

    ! --------------------------------

    ! check NoahMP50 alarm. If alarm is ring, run model.

    write(fnest,'(i3.3)') n
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "NoahMP50 model alarm "//trim(fnest))

    if (alarmCheck) Then
        do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            dt  = LIS_rc%ts
            row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
            lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

            ! retrieve forcing data from NoahMP50_struc(n)%noahmp50(t) and assign to 1-D NoahmpIO variables
            ! T_PHY: air temperature
            NoahmpIO%T_PHY(1,1,1) = NoahMP50_struc(n)%noahmp50(t)%tair   / NoahMP50_struc(n)%forc_count
            ! Yeosang Yoon, for snow DA
            NoahMP50_struc(n)%noahmp50(t)%sfctmp = NoahmpIO%T_PHY(1,1,1)

            ! P8W: air pressure
            NoahmpIO%P8W(1,1,1) = NoahMP50_struc(n)%noahmp50(t)%psurf  / NoahMP50_struc(n)%forc_count

            ! U_PHY: U wind component
            NoahmpIO%U_PHY(1,1,1) = NoahMP50_struc(n)%noahmp50(t)%wind_e / NoahMP50_struc(n)%forc_count

            ! V_PHY: V wind component
            NoahmpIO%V_PHY(1,1,1) = NoahMP50_struc(n)%noahmp50(t)%wind_n / NoahMP50_struc(n)%forc_count

            ! QV_CURR: specific humidity
            NoahmpIO%QV_CURR(1,1,1) = NoahMP50_struc(n)%noahmp50(t)%qair   / NoahMP50_struc(n)%forc_count

            ! SWDOWN: downward solar radiation
            NoahmpIO%SWDOWN(1,1) = NoahMP50_struc(n)%noahmp50(t)%swdown / NoahMP50_struc(n)%forc_count

            ! GLW: downward longwave radiation
            NoahmpIO%GLW(1,1) = NoahMP50_struc(n)%noahmp50(t)%lwdown / NoahMP50_struc(n)%forc_count

            ! prcp: total precipitation (rainfall+snowfall)
            ! Noah-MP require total precipitation as forcing input. [mm per model timestep]
            ! T. Lahmers: Correct total precip for cases when model time step > forcing timestep. 
            ! Edit suggested by D. Mocko and K. Arsenault
            !if (NoahMP50_struc(n)%ts > LIS_rc%ts) then
                NoahmpIO%RAINBL(1,1) = NoahMP50_struc(n)%ts * &
                                       (NoahMP50_struc(n)%noahmp50(t)%prcp / NoahMP50_struc(n)%forc_count)
            !else
            !   NoahmpIO%RAINBL(1,1) = dt * (NoahMP50_struc(n)%noahmp50(t)%prcp / NoahMP50_struc(n)%forc_count)
            !endif

            !ag(05Jan2021)
            ! rivsto/fldsto: River storage and flood storage
            ! NoahMP50_struc(n)%noahmp50(t)%rivsto and NoahMP50_struc(n)%noahmp50(t)%fldsto
            ! are updated in noahmp50_getsws_hymap2.F90
            NoahmpIO%rivsto(1,1) = NoahMP50_struc(n)%noahmp50(t)%rivsto
            NoahmpIO%fldsto(1,1) = NoahMP50_struc(n)%noahmp50(t)%fldsto
            NoahmpIO%fldfrc(1,1) = NoahMP50_struc(n)%noahmp50(t)%fldfrc

            ! check validity of tair
            if(NoahmpIO%T_PHY(1,1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable tair in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of psurf
            if(NoahmpIO%P8W(1,1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable psurf in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of wind_e
            if(NoahmpIO%U_PHY(1,1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable wind_e in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of wind_n
            if(NoahmpIO%V_PHY(1,1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable wind_n in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of qair
            if(NoahmpIO%QV_CURR(1,1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable qair in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of swdown
            if(NoahmpIO%SWDOWN(1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable swdown in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of lwdown
            if(NoahmpIO%GLW(1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable lwdown in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of prcp
            if(NoahmpIO%RAINBL(1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable prcp in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            !

            !ag (05Jan2021)
            ! check validity of rivsto
            if(NoahmpIO%rivsto(1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable rivsto in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of fldsto
            if(NoahmpIO%fldsto(1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable fldsto in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of fldfrc
            if(NoahmpIO%fldfrc(1,1) .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable fldfrc in Noah-MP.5.0"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! get parameters
            NoahmpIO%xlat(1,1)          = lat
            NoahmpIO%xlon(1,1)          = lon
            NoahmpIO%msftx(1,1)         = 1.0
            NoahmpIO%msfty(1,1)         = 1.0
            NoahmpIO%year               = LIS_rc%yr
            NoahmpIO%month              = LIS_rc%mo
            NoahmpIO%day                = LIS_rc%da
            NoahmpIO%hour               = LIS_rc%hr
            NoahmpIO%minute             = LIS_rc%mn
            NoahmpIO%ttile              = t
            NoahmpIO%itimestep          = LIS_rc%tscount(n)
            NoahmpIO%dtbl               = NoahMP50_struc(n)%ts
            NoahmpIO%soiltstep          = NoahMP50_struc(n)%ts_soil
            NoahmpIO%dzs(:)             = NoahMP50_struc(n)%sldpth(:)
            NoahmpIO%nsoil              = NoahMP50_struc(n)%nsoil
            NoahmpIO%nsnow              = NoahMP50_struc(n)%nsnow
            NoahmpIO%dx                 = NoahMP50_struc(n)%dx
            NoahmpIO%dy                 = NoahMP50_struc(n)%dy
            NoahmpIO%ivgtyp(1,1)        = NoahMP50_struc(n)%noahmp50(t)%vegetype
            NoahmpIO%isltyp(1,1)        = NoahMP50_struc(n)%noahmp50(t)%soiltype
            ! Multiply shdfac by 100.0 because noahmpdrv.f90, expects it in units of percentage, not fraction.
            NoahmpIO%shdfac_monthly(1,:,1) = NoahMP50_struc(n)%noahmp50(t)%shdfac_monthly(:) * 100.0
            NoahmpIO%tmn(1,1)           = NoahMP50_struc(n)%noahmp50(t)%tbot
            NoahmpIO%urban_vegtype(1,1) = LIS_rc%urbanclass
            NoahmpIO%cropcat(1,1)       = LIS_rc%cropclass
            NoahmpIO%IOPT_DVEG          = NoahMP50_struc(n)%dveg_opt
            NoahmpIO%IOPT_CRS           = NoahMP50_struc(n)%crs_opt
            NoahmpIO%IOPT_BTR           = NoahMP50_struc(n)%btr_opt
            NoahmpIO%IOPT_RUNSRF        = NoahMP50_struc(n)%runsfc_opt
            NoahmpIO%IOPT_RUNSUB        = NoahMP50_struc(n)%runsub_opt
            NoahmpIO%IOPT_SFC           = NoahMP50_struc(n)%sfc_opt
            NoahmpIO%IOPT_FRZ           = NoahMP50_struc(n)%frz_opt
            NoahmpIO%IOPT_INF           = NoahMP50_struc(n)%inf_opt
            NoahmpIO%IOPT_RAD           = NoahMP50_struc(n)%rad_opt
            NoahmpIO%IOPT_ALB           = NoahMP50_struc(n)%alb_opt
            NoahmpIO%IOPT_SNF           = NoahMP50_struc(n)%snf_opt
            NoahmpIO%IOPT_TKSNO         = NoahMP50_struc(n)%tksno_opt
            NoahmpIO%IOPT_TBOT          = NoahMP50_struc(n)%tbot_opt
            NoahmpIO%IOPT_STC           = NoahMP50_struc(n)%stc_opt
            NoahmpIO%IOPT_GLA           = NoahMP50_struc(n)%gla_opt
            NoahmpIO%IOPT_RSF           = NoahMP50_struc(n)%rsf_opt
            NoahmpIO%IOPT_SOIL          = NoahMP50_struc(n)%soil_opt
            NoahmpIO%IOPT_PEDO          = NoahMP50_struc(n)%pedo_opt
            NoahmpIO%IOPT_CROP          = NoahMP50_struc(n)%crop_opt
            NoahmpIO%IOPT_IRR           = NoahMP50_struc(n)%irr_opt
            NoahmpIO%IOPT_IRRM          = NoahMP50_struc(n)%irrm_opt
            NoahmpIO%IOPT_INFDV         = NoahMP50_struc(n)%infdv_opt
            NoahmpIO%IOPT_TDRN          = NoahMP50_struc(n)%tdrn_opt
            NoahmpIO%iz0tlnd            = 0
            NoahmpIO%sf_urban_physics   = NoahMP50_struc(n)%urban_opt
! Multiply reference height by 2.0 because module_sf_noahmpdrv
! expects this variable to be in terms of a thickness of the
! atmospheric layers, and it later divides this value by 2.0.
! Thus, the LIS user should specify the exact height of the
! reference in lis.config, and module_sf_noahmpdrv will then
! correctly use this actual value.  This code is confirmed in
! the HRLDAS driver, which also multiplies this value by 2.0.
! 11/30/2018 - dmm
            NoahmpIO%dz8w(1,1,1)        = NoahMP50_struc(n)%dz8w * 2.0

            if (NoahmpIO%IOPT_CROP > 0) then 
               NoahmpIO%planting(1,1)   = NoahMP50_struc(n)%noahmp50(t)%planting
               NoahmpIO%harvest(1,1)    = NoahMP50_struc(n)%noahmp50(t)%harvest
               NoahmpIO%season_gdd(1,1) = NoahMP50_struc(n)%noahmp50(t)%season_gdd
            else
               NoahmpIO%planting(1,1)   = 0.0
               NoahmpIO%harvest(1,1)    = 0.0
               NoahmpIO%season_gdd(1,1) = 0.0
            endif

! Zhuo Wang tested on 11/15/2018, not read from LDT-generated netcdf input file
            if (NoahmpIO%IOPT_SOIL .eq. 2) then 
               NoahmpIO%soilcL1(1,1) = NoahMP50_struc(n)%noahmp50(t)%soilcL1
               NoahmpIO%soilcL2(1,1) = NoahMP50_struc(n)%noahmp50(t)%soilcL2
               NoahmpIO%soilcL3(1,1) = NoahMP50_struc(n)%noahmp50(t)%soilcL3
               NoahmpIO%soilcL4(1,1) = NoahMP50_struc(n)%noahmp50(t)%soilcL4
            endif
            if (NoahmpIO%IOPT_SOIL .eq. 3) then 
               NoahmpIO%soilcomp(1,:,1) = NoahMP50_struc(n)%noahmp50(t)%soilcomp(:)
            endif

            ! for irrigation
            if (NoahmpIO%IOPT_IRR > 0) then
               NoahmpIO%irnumsi(1,1) = NoahMP50_struc(n)%noahmp50(t)%irnumsi
               NoahmpIO%irnummi(1,1) = NoahMP50_struc(n)%noahmp50(t)%irnummi
               NoahmpIO%irnumfi(1,1) = NoahMP50_struc(n)%noahmp50(t)%irnumfi
               NoahmpIO%irfract(1,1) = NoahMP50_struc(n)%noahmp50(t)%irfract
               NoahmpIO%sifract(1,1) = NoahMP50_struc(n)%noahmp50(t)%sifract
               NoahmpIO%mifract(1,1) = NoahMP50_struc(n)%noahmp50(t)%mifract
               NoahmpIO%fifract(1,1) = NoahMP50_struc(n)%noahmp50(t)%fifract
               NoahmpIO%irwatsi(1,1) = NoahMP50_struc(n)%noahmp50(t)%irwatsi
               NoahmpIO%irwatmi(1,1) = NoahMP50_struc(n)%noahmp50(t)%irwatmi
               NoahmpIO%irwatfi(1,1) = NoahMP50_struc(n)%noahmp50(t)%irwatfi
               NoahmpIO%ireloss(1,1) = NoahMP50_struc(n)%noahmp50(t)%ireloss
               NoahmpIO%irrsplh(1,1) = NoahMP50_struc(n)%noahmp50(t)%irrsplh
               NoahmpIO%irsivol(1,1) = NoahMP50_struc(n)%noahmp50(t)%irsivol
               NoahmpIO%irmivol(1,1) = NoahMP50_struc(n)%noahmp50(t)%irmivol
               NoahmpIO%irfivol(1,1) = NoahMP50_struc(n)%noahmp50(t)%irfivol
            endif

            ! for tile drainage
            if (NoahmpIO%IOPT_TDRN > 0) then
               NoahmpIO%TD_FRACTION(1,1)= NoahMP50_struc(n)%noahmp50(t)%tdfract
               NoahmpIO%qtdrain(1,1)    = NoahMP50_struc(n)%noahmp50(t)%qtdrain
               NoahmpIO%qtiledrain(1,1) = NoahMP50_struc(n)%noahmp50(t)%qtdrainflx
            endif

            ! for MMF groundwater
            if (NoahmpIO%IOPT_RUNSUB == 5) then
               NoahmpIO%fdepthxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%fdepth
               NoahmpIO%eqzwt(1,1)       = NoahMP50_struc(n)%noahmp50(t)%eqzwt
               NoahmpIO%riverbedxy(1,1)  = NoahMP50_struc(n)%noahmp50(t)%riverbed
               NoahmpIO%rivercondxy(1,1) = NoahMP50_struc(n)%noahmp50(t)%rivercond
               NoahmpIO%pexpxy(1,1)      = NoahMP50_struc(n)%noahmp50(t)%pexp
               NoahmpIO%areaxy(1,1)      = NoahMP50_struc(n)%noahmp50(t)%area
               NoahmpIO%qrfsxy(1,1)      = NoahMP50_struc(n)%noahmp50(t)%qrfs
               NoahmpIO%qspringxy(1,1)   = NoahMP50_struc(n)%noahmp50(t)%qspring
               NoahmpIO%qrfxy(1,1)       = NoahMP50_struc(n)%noahmp50(t)%qrf
               NoahmpIO%qspringsxy(1,1)  = NoahMP50_struc(n)%noahmp50(t)%qsprings
               NoahmpIO%qslatxy(1,1)     = NoahMP50_struc(n)%noahmp50(t)%qslat
               NoahmpIO%rechclim(1,1)    = NoahMP50_struc(n)%noahmp50(t)%rechclim
               NoahmpIO%rivermask(1,1)   = NoahMP50_struc(n)%noahmp50(t)%rivermask
               NoahmpIO%nonriverxy(1,1)  = NoahMP50_struc(n)%noahmp50(t)%nonriver 
            endif

            ! get state variables
            NoahmpIO%sfcrunoff(1,1)   = NoahMP50_struc(n)%noahmp50(t)%sfcrunoff
            NoahmpIO%udrunoff(1,1)    = NoahMP50_struc(n)%noahmp50(t)%udrrunoff
            NoahmpIO%smois(1,:,1)     = NoahMP50_struc(n)%noahmp50(t)%smc(:)
            NoahmpIO%sh2o(1,:,1)      = NoahMP50_struc(n)%noahmp50(t)%sh2o(:)
            NoahmpIO%tslb(1,:,1)      = NoahMP50_struc(n)%noahmp50(t)%tslb(:)
            NoahmpIO%snow(1,1)        = NoahMP50_struc(n)%noahmp50(t)%sneqv
            NoahmpIO%snowh(1,1)       = NoahMP50_struc(n)%noahmp50(t)%snowh
            NoahmpIO%canwat(1,1)      = NoahMP50_struc(n)%noahmp50(t)%canwat
            NoahmpIO%acsnom(1,1)      = NoahMP50_struc(n)%noahmp50(t)%acsnom
            NoahmpIO%acsnow(1,1)      = NoahMP50_struc(n)%noahmp50(t)%acsnow
            NoahmpIO%isnowxy(1,1)     = NoahMP50_struc(n)%noahmp50(t)%isnow
            NoahmpIO%tvxy(1,1)        = NoahMP50_struc(n)%noahmp50(t)%tv
            NoahmpIO%tgxy(1,1)        = NoahMP50_struc(n)%noahmp50(t)%tg
            NoahmpIO%canicexy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%canice
            NoahmpIO%canliqxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%canliq
            NoahmpIO%eahxy(1,1)       = NoahMP50_struc(n)%noahmp50(t)%eah
            NoahmpIO%tahxy(1,1)       = NoahMP50_struc(n)%noahmp50(t)%tah
            NoahmpIO%cmxy(1,1)        = NoahMP50_struc(n)%noahmp50(t)%cm
            NoahmpIO%chxy(1,1)        = NoahMP50_struc(n)%noahmp50(t)%ch
            NoahmpIO%fwetxy(1,1)      = NoahMP50_struc(n)%noahmp50(t)%fwet
            NoahmpIO%sneqvoxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%sneqvo
            NoahmpIO%alboldxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%albold
            NoahmpIO%qsnowxy(1,1)     = NoahMP50_struc(n)%noahmp50(t)%qsnow
            NoahmpIO%wslakexy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%wslake
            NoahmpIO%zwtxy(1,1)       = NoahMP50_struc(n)%noahmp50(t)%zwt
            NoahmpIO%waxy(1,1)        = NoahMP50_struc(n)%noahmp50(t)%wa
            NoahmpIO%wtxy(1,1)        = NoahMP50_struc(n)%noahmp50(t)%wt
            NoahmpIO%tsnoxy(1,:,1)    = NoahMP50_struc(n)%noahmp50(t)%tsno(:)
            NoahmpIO%zsnsoxy(1,:,1)   = NoahMP50_struc(n)%noahmp50(t)%zss(:)
            NoahmpIO%snicexy(1,:,1)   = NoahMP50_struc(n)%noahmp50(t)%snowice(:)
            NoahmpIO%snliqxy(1,:,1)   = NoahMP50_struc(n)%noahmp50(t)%snowliq(:)
            NoahmpIO%lfmassxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%lfmass
            NoahmpIO%rtmassxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%rtmass
            NoahmpIO%stmassxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%stmass
            NoahmpIO%woodxy(1,1)      = NoahMP50_struc(n)%noahmp50(t)%wood
            NoahmpIO%stblcpxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%stblcp
            NoahmpIO%fastcpxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%fastcp

            ! additional accumulated variables
            NoahmpIO%ACC_SSOILXY(1,1)    = NoahMP50_struc(n)%noahmp50(t)%accssoil
            NoahmpIO%ACC_QINSURXY(1,1)   = NoahMP50_struc(n)%noahmp50(t)%accqinsur
            NoahmpIO%ACC_QSEVAXY(1,1)    = NoahMP50_struc(n)%noahmp50(t)%accqseva
            NoahmpIO%ACC_ETRANIXY(1,:,1) = NoahMP50_struc(n)%noahmp50(t)%accetrani(:)
            NoahmpIO%ACC_DWATERXY(1,1)   = NoahMP50_struc(n)%noahmp50(t)%accdwater
            NoahmpIO%ACC_PRCPXY(1,1)     = NoahMP50_struc(n)%noahmp50(t)%accprcp
            NoahmpIO%ACC_ECANXY(1,1)     = NoahMP50_struc(n)%noahmp50(t)%accecan
            NoahmpIO%ACC_ETRANXY(1,1)    = NoahMP50_struc(n)%noahmp50(t)%accetran
            NoahmpIO%ACC_EDIRXY (1,1)    = NoahMP50_struc(n)%noahmp50(t)%accedir

! DMM - If dynamic vegetation option DVEG = 7, 8, or 9 for "input LAI",
! then send LAI/SAI from input to the Noah-MP physics.  If any
! tile has an undefined LAI/SAI value, instead use the value from the
! MPTABLE file for that vegetation class and for the month.
            if ((NoahmpIO%IOPT_DVEG .ge. 7).and.(NoahmpIO%IOPT_DVEG .le. 9)) then
               tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
! If "LAI data source:" is set to "none" for these three Noah-MP
! input LAI vegetation options, stop the run.
               if (LIS_rc%useLAImap(n).ne."none") then
                  NoahmpIO%lai(1,1) = LIS_lai(n)%tlai(tid)
               else
                  write(LIS_logunit,*)                                 &
                       '[ERR] Attempting to use input LAI, however'
                  write(LIS_logunit,*)                                 &
                       '[ERR] "LAI data source:" is set to "none".'
                  call LIS_endrun()
               endif
! If "SAI data source:" is set to "none" for these three Noah-MP-4.0.1
! input LAI vegetation options, fill in the SAI values from MPTABLE.
               if (LIS_rc%useSAImap(n).ne."none") then
                  NoahmpIO%xsaixy(1,1) = LIS_sai(n)%tsai(tid)
               endif
! If any LAI or SAI values are undefined at a tile,
! fill in the LAI or SAI values from MPTABLE.
               if (NoahmpIO%lai(1,1) .eq. LIS_rc%udef) then
                  NoahmpIO%lai(1,1) = NoahMP50_struc(n)%noahmp50(t)%lai
               endif
               if (NoahmpIO%xsaixy(1,1) .eq. LIS_rc%udef) then
                  NoahmpIO%xsaixy(1,1) = NoahMP50_struc(n)%noahmp50(t)%sai
               endif
            else
               NoahmpIO%lai(1,1)     = NoahMP50_struc(n)%noahmp50(t)%lai
               NoahmpIO%xsaixy(1,1)  = NoahMP50_struc(n)%noahmp50(t)%sai
            endif
            NoahmpIO%taussxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%tauss
            NoahmpIO%smoiseq(1,:,1)  = NoahMP50_struc(n)%noahmp50(t)%smoiseq(:)
            NoahmpIO%smcwtdxy(1,1)   = NoahMP50_struc(n)%noahmp50(t)%smcwtd
            NoahmpIO%deeprechxy(1,1) = NoahMP50_struc(n)%noahmp50(t)%deeprech
            NoahmpIO%rechxy(1,1)     = NoahMP50_struc(n)%noahmp50(t)%rech
            NoahmpIO%grainxy(1,1)    = NoahMP50_struc(n)%noahmp50(t)%grain
            NoahmpIO%gddxy(1,1)      = NoahMP50_struc(n)%noahmp50(t)%gdd
            NoahmpIO%pgsxy(1,1)      = NoahMP50_struc(n)%noahmp50(t)%pgs
            NoahmpIO%sfcheadrt(1,1)  = NoahMP50_struc(n)%noahmp50(t)%sfcheadrt

! Calculate water storages at start of timestep
            startsm = 0.0
            do i = 1,NoahmpIO%nsoil
               startsm = startsm +                                     &
                         (NoahmpIO%smois(1,i,1) * NoahmpIO%dzs(i) * LIS_CONST_RHOFW)
            enddo
            startswe = NoahmpIO%snow(1,1)
            startint = NoahmpIO%canliqxy(1,1) + NoahmpIO%canicexy(1,1)
            startgw  = NoahmpIO%waxy(1,1)

            ! call model physics
            call noahmp_driver_50(n, NoahMP50_struc(n)%noahmp50(t)%param)

            ! save state variables from NoahmpIO 1-D variables to global variables
            NoahMP50_struc(n)%noahmp50(t)%sfcrunoff       = NoahmpIO%sfcrunoff(1,1)
            NoahMP50_struc(n)%noahmp50(t)%udrrunoff       = NoahmpIO%udrunoff(1,1)
            NoahMP50_struc(n)%noahmp50(t)%smc(:)          = NoahmpIO%smois(1,:,1)
            NoahMP50_struc(n)%noahmp50(t)%sh2o(:)         = NoahmpIO%sh2o(1,:,1)
            NoahMP50_struc(n)%noahmp50(t)%tslb(:)         = NoahmpIO%tslb(1,:,1)
            NoahMP50_struc(n)%noahmp50(t)%sneqv           = NoahmpIO%snow(1,1)
            NoahMP50_struc(n)%noahmp50(t)%snowh           = NoahmpIO%snowh(1,1)
            NoahMP50_struc(n)%noahmp50(t)%canwat          = NoahmpIO%canwat(1,1)
            NoahMP50_struc(n)%noahmp50(t)%acsnom          = NoahmpIO%acsnom(1,1)
            NoahMP50_struc(n)%noahmp50(t)%acsnow          = NoahmpIO%acsnow(1,1)
            NoahMP50_struc(n)%noahmp50(t)%isnow           = NoahmpIO%isnowxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%tv              = NoahmpIO%tvxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%tg              = NoahmpIO%tgxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%canice          = NoahmpIO%canicexy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%canliq          = NoahmpIO%canliqxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%eah             = NoahmpIO%eahxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%tah             = NoahmpIO%tahxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%cm              = NoahmpIO%cmxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%ch              = NoahmpIO%chxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%fwet            = NoahmpIO%fwetxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%sneqvo          = NoahmpIO%sneqvoxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%albold          = NoahmpIO%alboldxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%qsnow           = NoahmpIO%qsnowxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%wslake          = NoahmpIO%wslakexy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%zwt             = NoahmpIO%zwtxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%wa              = NoahmpIO%waxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%wt              = NoahmpIO%wtxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%tsno(:)         = NoahmpIO%tsnoxy(1,:,1)
            NoahMP50_struc(n)%noahmp50(t)%zss(:)          = NoahmpIO%zsnsoxy(1,:,1)
            NoahMP50_struc(n)%noahmp50(t)%snowice(:)      = NoahmpIO%snicexy(1,:,1)
            NoahMP50_struc(n)%noahmp50(t)%snowliq(:)      = NoahmpIO%snliqxy(1,:,1)
            NoahMP50_struc(n)%noahmp50(t)%lfmass          = NoahmpIO%lfmassxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%rtmass          = NoahmpIO%rtmassxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%stmass          = NoahmpIO%stmassxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%wood            = NoahmpIO%woodxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%stblcp          = NoahmpIO%stblcpxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%fastcp          = NoahmpIO%fastcpxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%lai             = NoahmpIO%lai(1,1)
            NoahMP50_struc(n)%noahmp50(t)%sai             = NoahmpIO%xsaixy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%tauss           = NoahmpIO%taussxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%smoiseq(:)      = NoahmpIO%smoiseq(1,:,1)
            NoahMP50_struc(n)%noahmp50(t)%smcwtd          = NoahmpIO%smcwtdxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%deeprech        = NoahmpIO%deeprechxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%rech            = NoahmpIO%rechxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%grain           = NoahmpIO%grainxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%gdd             = NoahmpIO%gddxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%pgs             = NoahmpIO%pgsxy(1,1)

            ! save output variables from NoahmpIO variables to global variables
            NoahMP50_struc(n)%noahmp50(t)%tsk       = NoahmpIO%tsk(1,1)
            NoahMP50_struc(n)%noahmp50(t)%hfx       = NoahmpIO%hfx(1,1)
            NoahMP50_struc(n)%noahmp50(t)%qfx       = NoahmpIO%qfx(1,1)
            NoahMP50_struc(n)%noahmp50(t)%lh        = NoahmpIO%lh(1,1)
            NoahMP50_struc(n)%noahmp50(t)%grdflx    = NoahmpIO%grdflx(1,1)
            NoahMP50_struc(n)%noahmp50(t)%albedo    = NoahmpIO%albedo(1,1)
            NoahMP50_struc(n)%noahmp50(t)%snowc     = NoahmpIO%snowc(1,1)
            NoahMP50_struc(n)%noahmp50(t)%emiss     = NoahmpIO%emiss(1,1)
            NoahMP50_struc(n)%noahmp50(t)%rs        = NoahmpIO%rs(1,1)
            NoahMP50_struc(n)%noahmp50(t)%t2mv      = NoahmpIO%t2mvxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%t2mb      = NoahmpIO%t2mbxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%q2mv      = NoahmpIO%q2mvxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%q2mb      = NoahmpIO%q2mbxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%trad      = NoahmpIO%tradxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%nee       = NoahmpIO%neexy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%gpp       = NoahmpIO%gppxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%npp       = NoahmpIO%nppxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%fveg      = NoahmpIO%fvegxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%runsf     = NoahmpIO%runsfxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%runsb     = NoahmpIO%runsbxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%ecan      = NoahmpIO%ecanxy(1,1)
! Direct soil evaporation does not include sublimation of the snowpack
! on the soil (by the strict ALMA definition of ESoil). - David Mocko
            NoahMP50_struc(n)%noahmp50(t)%edir      = NoahmpIO%edirxy(1,1) - NoahmpIO%qsnsubxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%etran     = NoahmpIO%etranxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%rainf     = NoahmpIO%rainlsm(1,1)
            NoahMP50_struc(n)%noahmp50(t)%snowf     = NoahmpIO%snowlsm(1,1)
            NoahMP50_struc(n)%noahmp50(t)%fsa       = NoahmpIO%fsaxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%fira      = NoahmpIO%firaxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%apar      = NoahmpIO%aparxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%psn       = NoahmpIO%psnxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%sav       = NoahmpIO%savxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%sag       = NoahmpIO%sagxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%rssun     = NoahmpIO%rssunxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%rssha     = NoahmpIO%rsshaxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%bgap      = NoahmpIO%bgapxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%wgap      = NoahmpIO%wgapxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%tgb       = NoahmpIO%tgbxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%tgv       = NoahmpIO%tgvxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%chv       = NoahmpIO%chvxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%chb       = NoahmpIO%chbxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%shg       = NoahmpIO%shgxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%shc       = NoahmpIO%shcxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%shb       = NoahmpIO%shbxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%evg       = NoahmpIO%evgxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%evb       = NoahmpIO%evbxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%ghv       = NoahmpIO%ghvxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%ghb       = NoahmpIO%ghbxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%irg       = NoahmpIO%irgxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%irc       = NoahmpIO%ircxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%irb       = NoahmpIO%irbxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%tr        = NoahmpIO%trxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%evc       = NoahmpIO%evcxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%chleaf    = NoahmpIO%chleafxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%chuc      = NoahmpIO%chucxy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%chv2      = NoahmpIO%chv2xy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%chb2      = NoahmpIO%chb2xy(1,1)
            NoahMP50_struc(n)%noahmp50(t)%infxs1rt  = NoahmpIO%infxsrt(1,1)
            NoahMP50_struc(n)%noahmp50(t)%soldrain1rt = NoahmpIO%soldrain(1,1)
            NoahMP50_struc(n)%noahmp50(t)%sfcheadrt    = NoahmpIO%sfcheadrt(1,1)

            ! additional accumulated variables
            NoahMP50_struc(n)%noahmp50(t)%accssoil     = NoahmpIO%ACC_SSOILXY(1,1)
            NoahMP50_struc(n)%noahmp50(t)%accqinsur    = NoahmpIO%ACC_QINSURXY(1,1)
            NoahMP50_struc(n)%noahmp50(t)%accqseva     = NoahmpIO%ACC_QSEVAXY(1,1)
            NoahMP50_struc(n)%noahmp50(t)%accetrani(:) = NoahmpIO%ACC_ETRANIXY(1,:,1)
            NoahMP50_struc(n)%noahmp50(t)%accdwater    = NoahmpIO%ACC_DWATERXY(1,1)
            NoahMP50_struc(n)%noahmp50(t)%accprcp      = NoahmpIO%ACC_PRCPXY(1,1)
            NoahMP50_struc(n)%noahmp50(t)%accecan      = NoahmpIO%ACC_ECANXY(1,1)
            NoahMP50_struc(n)%noahmp50(t)%accetran     = NoahmpIO%ACC_ETRANXY(1,1)
            NoahMP50_struc(n)%noahmp50(t)%accedir      = NoahmpIO%ACC_EDIRXY (1,1)

            ! for irrigation
            if (NoahmpIO%IOPT_IRR > 0) then
               NoahMP50_struc(n)%noahmp50(t)%irnumsi = NoahmpIO%irnumsi(1,1)
               NoahMP50_struc(n)%noahmp50(t)%irnummi = NoahmpIO%irnummi(1,1)
               NoahMP50_struc(n)%noahmp50(t)%irnumfi = NoahmpIO%irnumfi(1,1)
               NoahMP50_struc(n)%noahmp50(t)%irwatsi = NoahmpIO%irwatsi(1,1)
               NoahMP50_struc(n)%noahmp50(t)%irwatmi = NoahmpIO%irwatmi(1,1)
               NoahMP50_struc(n)%noahmp50(t)%irwatfi = NoahmpIO%irwatfi(1,1)
               NoahMP50_struc(n)%noahmp50(t)%ireloss = NoahmpIO%ireloss(1,1)
               NoahMP50_struc(n)%noahmp50(t)%irrsplh = NoahmpIO%irrsplh(1,1)
               NoahMP50_struc(n)%noahmp50(t)%irsivol = NoahmpIO%irsivol(1,1)
               NoahMP50_struc(n)%noahmp50(t)%irmivol = NoahmpIO%irmivol(1,1)
               NoahMP50_struc(n)%noahmp50(t)%irfivol = NoahmpIO%irfivol(1,1)
            endif

            ! for tile drainage
            if (NoahmpIO%IOPT_TDRN > 0) then
               NoahMP50_struc(n)%noahmp50(t)%qtdrain    = NoahmpIO%qtdrain(1,1)
               NoahMP50_struc(n)%noahmp50(t)%qtdrainflx = NoahmpIO%qtiledrain(1,1)  
            endif

            ! for MMF groundwater
            if (NoahmpIO%IOPT_RUNSUB == 5) then
               NoahMP50_struc(n)%noahmp50(t)%fdepth    = NoahmpIO%fdepthxy(1,1)
               NoahMP50_struc(n)%noahmp50(t)%eqzwt     = NoahmpIO%eqzwt(1,1)
               NoahMP50_struc(n)%noahmp50(t)%riverbed  = NoahmpIO%riverbedxy(1,1)
               NoahMP50_struc(n)%noahmp50(t)%rivercond = NoahmpIO%rivercondxy(1,1)
               NoahMP50_struc(n)%noahmp50(t)%pexp      = NoahmpIO%pexpxy(1,1)
               NoahMP50_struc(n)%noahmp50(t)%area      = NoahmpIO%areaxy(1,1)
               NoahMP50_struc(n)%noahmp50(t)%qrfs      = NoahmpIO%qrfsxy(1,1)
               NoahMP50_struc(n)%noahmp50(t)%qspring   = NoahmpIO%qspringxy(1,1)
               NoahMP50_struc(n)%noahmp50(t)%qrf       = NoahmpIO%qrfxy(1,1)  
               NoahMP50_struc(n)%noahmp50(t)%qsprings  = NoahmpIO%qspringsxy(1,1)
               NoahMP50_struc(n)%noahmp50(t)%qslat     = NoahmpIO%qslatxy(1,1)
            endif

            ! EMK Update RHMin for 557WW
            if (NoahmpIO%T_PHY(1,1,1) .lt. &
                 noahmp50_struc(n)%noahmp50(t)%tair_agl_min) then
               noahmp50_struc(n)%noahmp50(t)%tair_agl_min = NoahmpIO%T_PHY(1,1,1)
               ! Use formulation based on Noah.3.6 code, which treats
               ! q2sat as saturated specific humidity
               tmp_es = 611.0*exp(2.501E6/461.0*(1./273.15 - 1./NoahmpIO%T_PHY(1,1,1)))
               tmp_q2sat = 0.622*tmp_es/(NoahmpIO%P8W(1,1,1)-(1.-0.622)*tmp_es)
               noahmp50_struc(n)%noahmp50(t)%rhmin = NoahmpIO%QV_CURR(1,1,1) / tmp_q2sat
            endif

            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RHMIN, &
             value=noahmp50_struc(n)%noahmp50(t)%rhmin, &
             vlevel=1, unit="-", direction="-",&
             surface_type=LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RHMIN, &
             value=(100.*noahmp50_struc(n)%noahmp50(t)%rhmin), &
             vlevel=1, unit="%", direction="-",&
             surface_type=LIS_rc%lsm_index)

            ![ 1] output variable: tsk (unit=K  ). ***  surface radiative temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RADT, value = NoahMP50_struc(n)%noahmp50(t)%tsk, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 2] output variable: fsh (unit=W/m2). ***  sensible heat flux to atmosphere
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QH, value = NoahMP50_struc(n)%noahmp50(t)%hfx, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 3] output variable: lh (unit=W/m2). ***  latent heat flux
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QLE, value = NoahMP50_struc(n)%noahmp50(t)%lh, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 4] output variable: grdflx (unit=W/m2). ***  ground/snow heat flux to soil
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QG, value = NoahMP50_struc(n)%noahmp50(t)%grdflx, &
                                              vlevel=1, unit="W m-2", direction="DN", surface_type = LIS_rc%lsm_index)

            ![ 5] output variable: albedo (unit=- ). ***  total grid albedo
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO, value = NoahMP50_struc(n)%noahmp50(t)%albedo, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            !if (NoahMP50_struc(n)%noahmp50(t)%albedo .ne. LIS_rc%udef) &
            !   NoahMP50_struc(n)%noahmp50(t)%albedo = NoahMP50_struc(n)%noahmp50(t)%albedo * 100.0

            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO, value = NoahMP50_struc(n)%noahmp50(t)%albedo, &
            !                                  vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 6] output variable: snowc (unit=-). ***  snow cover fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWCOVER, value = NoahMP50_struc(n)%noahmp50(t)%snowc, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWCOVER, value = (NoahMP50_struc(n)%noahmp50(t)%snowc*100.0), &
                                              vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 7] output variable: smc (unit=m3/m3). ***  volumetric soil moisture
            do i=1, NoahMP50_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = NoahMP50_struc(n)%noahmp50(t)%smc(i), &
                                                  vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 8] output variable: sh2o (unit=m3/m3). ***  equilibrium volumetric liquid soil moisture content
            do i=1, NoahMP50_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMLIQFRAC, value = NoahMP50_struc(n)%noahmp50(t)%sh2o(i), &
                                                  vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 9] output variable: tslb (unit=K). ***  soil temperature
            do i=1, NoahMP50_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILTEMP, value = NoahMP50_struc(n)%noahmp50(t)%tslb(i), &
                                                  vlevel=i, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 10] output variable: sneqv (unit=mm ). ***  snow water equivalent
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWE, value = NoahMP50_struc(n)%noahmp50(t)%sneqv, &
                                              vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 11] output variable: snowh (unit=m ). ***  physical snow depth
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDEPTH, value = NoahMP50_struc(n)%noahmp50(t)%snowh, &
                                              vlevel=1, unit="m ", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 12] output variable: canwat (unit=kg/m2). ***  total canopy water storage
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPINT, value = NoahMP50_struc(n)%noahmp50(t)%canwat, &
                                              vlevel=1, unit="kg m-2", direction="- ", surface_type = LIS_rc%lsm_index)

            ![ 13] output variable: emiss (unit=- ). ***  surface bulk emissivity
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EMISSFORC, value = NoahMP50_struc(n)%noahmp50(t)%emiss, &
                                              vlevel=1, unit="- ", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 14] output variable: rs (unit=s/m). ***  total stomatal resistance
            ! call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RS, value = NoahMP50_struc(n)%noahmp50(t)%rs, &
            !                                  vlevel=1, unit="s/m", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 15] output variable: isnow (unit=-). ***  actual number of snow layers
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOWN_NLAYER, value = -1.0*NoahMP50_struc(n)%noahmp50(t)%isnow, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 16] output variable: tv (unit=K ). ***  vegetation leaf temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_VEGT, value = NoahMP50_struc(n)%noahmp50(t)%tv, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 17] output variable: tg (unit=K). ***  averaged ground surface temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GROUNDAVGT, value = NoahMP50_struc(n)%noahmp50(t)%tg, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 18] output variable: canice (unit=mm). ***  canopy intercepted ice
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWEVEG, value = NoahMP50_struc(n)%noahmp50(t)%canice, &
                                              vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 19] output variable: canliq (unit=mm). ***  canopy intercepted liquid water
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_INTL, value = NoahMP50_struc(n)%noahmp50(t)%canliq, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 20] output variable: eah (unit=Pa  ). ***  canopy air vapor pressure
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_VP, value = NoahMP50_struc(n)%noahmp50(t)%eah, &
                                              vlevel=1, unit="Pa", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 21] output variable: tah (unit=K  ). ***  canopy air temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_TEMP, value = NoahMP50_struc(n)%noahmp50(t)%tah, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 22] output variable: cm (unit=s/m ). ***  bulk momentum drag coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CM, value = NoahMP50_struc(n)%noahmp50(t)%cm, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 23] output variable: ch (unit=s/m ). ***  bulk sensible heat exchange coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CH, value = NoahMP50_struc(n)%noahmp50(t)%ch, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 24] output variable: fwet (unit=-  ). ***  wetted or snowed fraction of canopy
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_WF, value = NoahMP50_struc(n)%noahmp50(t)%fwet, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 25] output variable: wslake (unit=mm). ***  lake water storage
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKEWATER, value = NoahMP50_struc(n)%noahmp50(t)%wslake, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 26] output variable: zwt (unit=m). ***  water table depth
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WATERTABLED, value = NoahMP50_struc(n)%noahmp50(t)%zwt, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 27] output variable: wa (unit=mm). ***  water storage in the "aquifer"
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GWS, value = NoahMP50_struc(n)%noahmp50(t)%wa, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 28] output variable: wt (unit=mm). ***  water in aquifer and saturated soil
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WT_AQUI_SATSOIL, value = NoahMP50_struc(n)%noahmp50(t)%wt, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 29] output variable: tsno (unit=K). ***  snow layer temperature
            do i=1, NoahMP50_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWTPROF, value = NoahMP50_struc(n)%noahmp50(t)%tsno(i), &
                                                  vlevel=i, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 30] output variable: snowice (unit=mm ). ***  snow layer ice
            do i=1, NoahMP50_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWICE, value = NoahMP50_struc(n)%noahmp50(t)%snowice(i), &
                                                  vlevel=i, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 31] output variable: snowliq (unit=mm ). ***  snow layer liquid water
            do i=1, NoahMP50_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWLIQ, value = NoahMP50_struc(n)%noahmp50(t)%snowliq(i), &
                                                  vlevel=i, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ! Yeosang Yoon, for snow DA
            ! output variable: z_snow (unit=m). ***  snow layer-bottom depth from snow surface
            do i=1, NoahMP50_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW_LBDFSS, value = NoahMP50_struc(n)%noahmp50(t)%zss(i), &
                                                  vlevel=i, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 32] output variable: lfmass (unit=g/m2). ***  leaf mass
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LEAFMASS, value = NoahMP50_struc(n)%noahmp50(t)%lfmass, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 33] output variable: rtmass (unit=g/m2 ). ***  mass of fine roots
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ROOTMASS, value = NoahMP50_struc(n)%noahmp50(t)%rtmass, &
                                              vlevel=1, unit="g m-2 ", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 34] output variable: stmass (unit=g/m2 ). ***  stem mass
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_STEMMASS, value = NoahMP50_struc(n)%noahmp50(t)%stmass, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 35] output variable: wood (unit=g/m2). ***  mass of wood including woody roots
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WOODMASS, value = NoahMP50_struc(n)%noahmp50(t)%wood, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 36] output variable: stblcp (unit=g/m2). ***  stable carbon in deep soil
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CARBON_DEEPSOIL, value = NoahMP50_struc(n)%noahmp50(t)%stblcp, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 37] output variable: fastcp (unit=g/m2 ). ***  short-lived carbon in shallow soil
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CARBON_SHALLOWSOIL, value = NoahMP50_struc(n)%noahmp50(t)%fastcp, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 38] output variable: lai (unit=-). ***  leave area index
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAI, value = NoahMP50_struc(n)%noahmp50(t)%lai, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 39] output variable: sai (unit=- ). ***  stem area index
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SAI, value = NoahMP50_struc(n)%noahmp50(t)%sai, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 40] output variable: tauss (unit=- ). ***  snow aging factor
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWAGE, value = NoahMP50_struc(n)%noahmp50(t)%tauss, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 41] output variable: smcwtd (unit=m3/m3). ***  soil moisture content in the layer to the water table when deep
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BETWEENWATER, value = NoahMP50_struc(n)%noahmp50(t)%smcwtd, &
                                              vlevel=1, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 42] output variable: deeprech (unit=m). ***  recharge to the water table when groundwater is deep
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QRECTOGW, value = NoahMP50_struc(n)%noahmp50(t)%deeprech, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 43] output variable: rech (unit=m). ***  recharge from the water table when groundwater is shallow
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QRECFROMGW, value = NoahMP50_struc(n)%noahmp50(t)%rech, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 44] output variable: t2mv (unit=K). ***  2-m air temperature over vegetated part
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_VEGE2MT, value = NoahMP50_struc(n)%noahmp50(t)%t2mv, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 45] output variable: t2mb (unit=K ). ***  2-m height air temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BARE2MT, value = NoahMP50_struc(n)%noahmp50(t)%t2mb, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 46] output variable: q2mv (unit=kg/kg). ***  2-m mixing ratio of vegetation part
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_VEGE2MQ2, value = NoahMP50_struc(n)%noahmp50(t)%q2mv, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 47] output variable: q2mb (unit=kg/kg). ***  2-m mixing ratio of bare ground part
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BARE2MQ2, value = NoahMP50_struc(n)%noahmp50(t)%q2mb, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 48] output variable: trad (unit=K). ***  surface radiative temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RADT, value = NoahMP50_struc(n)%noahmp50(t)%trad, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 49] output variable: nee (unit=g/m2/s ). ***  net ecosystem exchange of CO2
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_NEE, value = NoahMP50_struc(n)%noahmp50(t)%nee, &
                                              vlevel=1, unit="g m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)

            ![ 50] output variable: gpp (unit=g/m2/s  ). ***  gross primary assimilation of carbon
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GPP, value = NoahMP50_struc(n)%noahmp50(t)%gpp, &
                                              vlevel=1, unit="g m-2 s-1", direction="IN", surface_type = LIS_rc%lsm_index)

            ![ 51] output variable: npp (unit=g/m2/s). ***  net primary productivity of carbon
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_NPP, value = NoahMP50_struc(n)%noahmp50(t)%npp, &
                                              vlevel=1, unit="g m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)

            ![ 52] output variable: fveg (unit=-). ***  Noah-MP green vegetation fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GREENNESS, value = NoahMP50_struc(n)%noahmp50(t)%fveg, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 53] output variable: runsf (unit=mm/s). ***  surface runoff
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QS, value = NoahMP50_struc(n)%noahmp50(t)%runsf/NoahMP50_struc(n)%ts,&
                                              vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)

            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QS, value = NoahMP50_struc(n)%noahmp50(t)%runsf, &
                                              vlevel=1, unit="kg m-2", direction="OUT", surface_type = LIS_rc%lsm_index)

            ![ 54] output variable: runsb (unit=mm/s ). ***  baseflow (saturation excess)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSB, value = NoahMP50_struc(n)%noahmp50(t)%runsb/NoahMP50_struc(n)%ts,&
                                              vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)

            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSB, value = NoahMP50_struc(n)%noahmp50(t)%runsb, &
                                              vlevel=1, unit="kg m-2", direction="OUT", surface_type = LIS_rc%lsm_index)

            ! Combined output variable:  qtot (unit=mm | mm/s). *** total runoff
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QTOT, &
                             value = NoahMP50_struc(n)%noahmp50(t)%runsf/NoahMP50_struc(n)%ts + &
                                     NoahMP50_struc(n)%noahmp50(t)%runsb/NoahMP50_struc(n)%ts,  &
                             vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)

            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QTOT, &
                             value = NoahMP50_struc(n)%noahmp50(t)%runsf + NoahMP50_struc(n)%noahmp50(t)%runsb, &
                             vlevel=1, unit="kg m-2", direction="OUT", surface_type = LIS_rc%lsm_index)

            ![ 55] output variable: ecan (unit=mm/s ). ***  evaporation of intercepted water
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ECANOP, value = NoahMP50_struc(n)%noahmp50(t)%ecan, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 56] output variable: edir (unit=mm/s ). ***  soil surface evaporation rate
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ESOIL, value = NoahMP50_struc(n)%noahmp50(t)%edir, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 57] output variable: etran (unit=mm/s ). ***  transpiration rate
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_TVEG, value = NoahMP50_struc(n)%noahmp50(t)%etran, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 58] output variable: fsa (unit=W/m2). ***  total absorbed solar radiation
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWNET, value = NoahMP50_struc(n)%noahmp50(t)%fsa, &
                                              vlevel=1, unit="W m-2", direction="DN", surface_type = LIS_rc%lsm_index)

            ![ 59] output variable: fira (unit=W/m2 ). ***  total net longwave radiation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LWUP, value = NoahMP50_struc(n)%noahmp50(t)%fira, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 60] output variable: apar (unit=W/m2). ***  photosynthesis active energy by canopy
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_APAR, value = NoahMP50_struc(n)%noahmp50(t)%apar, &
                                              vlevel=1, unit="W m-2", direction="IN", surface_type = LIS_rc%lsm_index)

            ![ 61] output variable: psn (unit=umol/m2/s ). ***  total photosynthesis of CO2 [+]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PSCO2, value = NoahMP50_struc(n)%noahmp50(t)%psn, &
                                              vlevel=1, unit="umol m-2 s-1", direction="IN", surface_type = LIS_rc%lsm_index)

            ![ 62] output variable: sav (unit=W/m2 ). ***  solar radiation absorbed by vegetation
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SAV, value = NoahMP50_struc(n)%noahmp50(t)%sav, &
            !                                  vlevel=1, unit="W/m2 ", direction="IN", surface_type = LIS_rc%lsm_index)

            ![ 63] output variable: sag (unit=W/m2 ). ***  solar radiation absorbed by ground
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SAG, value = NoahMP50_struc(n)%noahmp50(t)%sag, &
            !                                  vlevel=1, unit="W/m2 ", direction="IN", surface_type = LIS_rc%lsm_index)

            ![ 64] output variable: rssun (unit=s/m). ***  sunlit leaf stomatal resistance
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RSSUN, value = NoahMP50_struc(n)%noahmp50(t)%rssun, &
                                              vlevel=1, unit="s m-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 65] output variable: rssha (unit=s/m). ***  shaded leaf stomatal resistance
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RSSHA, value = NoahMP50_struc(n)%noahmp50(t)%rssha, &
                                              vlevel=1, unit="s m-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 66] output variable: bgap (unit=-). ***  between gap fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BGAP, value = NoahMP50_struc(n)%noahmp50(t)%bgap, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 67] output variable: wgap (unit=- ). ***  within gap fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WGAP, value = NoahMP50_struc(n)%noahmp50(t)%wgap, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 68] output variable: tgb (unit=K). ***  bare ground temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BARESOILT, value = NoahMP50_struc(n)%noahmp50(t)%tgb, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 69] output variable: tgv (unit=K). ***  vegetated ground surface temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GROUNDVEGT, value = NoahMP50_struc(n)%noahmp50(t)%tgv, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 70] output variable: chv (unit=s/m). ***  sensible heat exchange coefficient over vegetated fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHV, value = NoahMP50_struc(n)%noahmp50(t)%chv, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 71] output variable: chb (unit=s/m). ***  sensible heat exchange coefficient over bare-ground fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHB, value = NoahMP50_struc(n)%noahmp50(t)%chb, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 72] output variable: shg (unit=W/m2     ). ***  get ground sensible heat [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SHG, value = NoahMP50_struc(n)%noahmp50(t)%shg, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 73] output variable: shc (unit=W/m2   ). ***  canopy sensible heat [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SHC, value = NoahMP50_struc(n)%noahmp50(t)%shc, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 74] output variable: shb (unit=W/m2     ). ***  bare ground sensible heat [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SHB, value = NoahMP50_struc(n)%noahmp50(t)%shb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 75] output variable: evg (unit=W/m2  ). ***  veg ground evaporation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVG, value = NoahMP50_struc(n)%noahmp50(t)%evg, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 76] output variable: evb (unit=W/m2  ). ***  bare soil evaporation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVB, value = NoahMP50_struc(n)%noahmp50(t)%evb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 77] output variable: ghv (unit=W/m2 ). ***  vegetated ground heat flux [+ to soil]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GHV, value = NoahMP50_struc(n)%noahmp50(t)%ghv, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 78] output variable: ghb (unit=W/m2 ). ***  bare ground heat flux [+ to soil]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GHB, value = NoahMP50_struc(n)%noahmp50(t)%ghb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 79] output variable: irg (unit=W/m2 ). ***  veg ground net long wave radiation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_IRV, value = NoahMP50_struc(n)%noahmp50(t)%irg, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 80] output variable: irc (unit=W/m2 ). ***  canopy net long wave radiation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_IRC, value = NoahMP50_struc(n)%noahmp50(t)%irc, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 81] output variable: irb (unit=W/m2 ). ***  bare net long wave radiation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_IRB, value = NoahMP50_struc(n)%noahmp50(t)%irb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 82] output variable: tr (unit=W/m2 ). ***  transpiration heat [to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_HTR, value = NoahMP50_struc(n)%noahmp50(t)%tr, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 83] output variable: evc (unit=W/m2 ). ***  canopy evaporation heat [to atm]
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVC, value = NoahMP50_struc(n)%noahmp50(t)%evc, &
            !                                  vlevel=1, unit="W/m2 ", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 84] output variable: chleaf (unit=m/s). ***  leaf exchange coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHLEAF, value = NoahMP50_struc(n)%noahmp50(t)%chleaf, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 85] output variable: chuc (unit=m/s). ***  under canopy exchange coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHUC, value = NoahMP50_struc(n)%noahmp50(t)%chuc, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 86] output variable: chv2 (unit=m/s). ***  veg 2m exchange coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHV2, value = NoahMP50_struc(n)%noahmp50(t)%chv2, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 87] output variable: chb2 (unit=m/s). ***  bare 2m exchange coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHB2, value = NoahMP50_struc(n)%noahmp50(t)%chb2, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 88] output variable: evap (unit=kg/m2/s). ***  total evapotranspiration
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVAP, value = NoahMP50_struc(n)%noahmp50(t)%qfx, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)
            !PET
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_POTEVAP, &
                 value=(NoahmpIO%fgev_pet(1,1)+NoahmpIO%fcev_pet(1,1)+NoahmpIO%fctr_pet(1,1)), vlevel=1,unit="W m-2",&
                  direction="UP",surface_type=LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_POTEVAP, &
                 value=(NoahmpIO%fgev_pet(1,1)+NoahmpIO%fcev_pet(1,1)+NoahmpIO%fctr_pet(1,1))/LVH2O, &
                 vlevel=1,unit="kg m-2 s-1", direction="UP",surface_type=LIS_rc%lsm_index)

            ![ 89] output variable: rainf (unit=kg/m2). ***  precipitation rate
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RAINF, value = NoahMP50_struc(n)%noahmp50(t)%rainf, &
                                              vlevel=1, unit="kg m-2 s-1", direction="DN", surface_type = LIS_rc%lsm_index)

            ![ 90] output variable: snowf (unit=kg/m2). ***  snowfall rate
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWF, value = NoahMP50_struc(n)%noahmp50(t)%snowf, &
                                              vlevel=1, unit="kg m-2 s-1", direction="DN", surface_type = LIS_rc%lsm_index)

            ![ 91] LWnet
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LWNET,vlevel=1,  &
                  value=(-1.0 * NoahMP50_struc(n)%noahmp50(t)%fira), &
                  unit="W m-2", direction="DN", surface_type=LIS_rc%lsm_index)

            ! Code added by Zhuo Wang on 02/28/2019
            ![ 92] output variable: qsnbot (unit=kg m-2 s-1). ***  melting water out of snow bottom 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSM, value = NoahmpIO%qsnbotxy(1,1), &
                  vlevel=1, unit="kg m-2 s-1", direction="S2L", surface_type = LIS_rc%lsm_index)

            ![ 93] output variable: subsnow (unit=kg m-2 s-1). ***  snow sublimation 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SUBSNOW, value = NoahmpIO%qsnsubxy(1,1), &
                  vlevel=1, unit="kg m-2 s-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 94] output variable: AvgSurfT (unit=K). *** average surface temperature 
            AvgSurfT_out = NoahMP50_struc(n)%noahmp50(t)%fveg * NoahMP50_struc(n)%noahmp50(t)%tv + &
                  (1.0-NoahMP50_struc(n)%noahmp50(t)%fveg) * NoahMP50_struc(n)%noahmp50(t)%tgb
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AVGSURFT, value = AvgSurfT_out, &
                  vlevel=1, unit="K", direction="-",surface_type = LIS_rc%lsm_index)

            ![ 95] TWS should be SWE + CanopInt + Soil moisture + WA - David Mocko
            TWS_out = NoahMP50_struc(n)%noahmp50(t)%sneqv
            if ((NoahMP50_struc(n)%noahmp50(t)%canliq.ge.0.0).and.   &
                (NoahMP50_struc(n)%noahmp50(t)%canice.ge.0.0)) then
               TWS_out = TWS_out +                                     &
                            (NoahMP50_struc(n)%noahmp50(t)%canliq  + &
                             NoahMP50_struc(n)%noahmp50(t)%canice)
            endif
            do i = 1,NoahMP50_struc(n)%nsoil
               TWS_out = TWS_out +                                     &
                      (NoahMP50_struc(n)%noahmp50(t)%smc(i)  *       &
                      NoahmpIO%dzs(i)*LIS_CONST_RHOFW)
            enddo
            if (NoahMP50_struc(n)%noahmp50(t)%wa.ge.0.0) then
               TWS_out = TWS_out + NoahMP50_struc(n)%noahmp50(t)%wa
            endif
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_TWS, value = TWS_out, &
                   vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 96] Qa - Advective energy - Heat transferred to a snow cover by rain
            !         - (unit=W m-2) - added by David Mocko
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QA, value = NoahmpIO%pahxy(1,1), &
                  vlevel=1, unit="W m-2",direction="DN",surface_type=LIS_rc%lsm_index)

! Added water balance change terms - David Mocko
            endsm = 0.0
            do i = 1,NoahmpIO%nsoil
               endsm = endsm +                                         &
                          (NoahmpIO%smois(1,i,1) * NoahmpIO%dzs(i) * LIS_CONST_RHOFW)
            enddo
            call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSOILMOIST,&
                     value=(endsm - startsm),vlevel=1,unit="kg m-2",   &
                     direction="INC",surface_type=LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSWE,      &
                     value=(NoahMP50_struc(n)%noahmp50(t)%sneqv -    &
                            startswe),                                 &
                     vlevel=1,unit="kg m-2",direction="INC",           &
                     surface_type=LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELINTERCEPT,&
                     value=((NoahMP50_struc(n)%noahmp50(t)%canliq +  &
                             NoahMP50_struc(n)%noahmp50(t)%canice) - &
                            startint),                                 &
                     vlevel=1,unit="kg m-2",direction="INC",           &
                     surface_type=LIS_rc%lsm_index)
! For now, the ALMA standard does not provide a variable for the
! change in groundwater storage.  Instead, temporarily writing it
! to the DELSURFSTOR (which is the change in surface water storage).
! This is only a temporary fix, until LIS_MOC_DELGROUNDWATER or
! a similarly-named variable is added into LIS_histDataMod.F90.
            call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSURFSTOR, &
                     value=(NoahMP50_struc(n)%noahmp50(t)%wa -       &
                            startgw),                                  &
                     vlevel=1,unit="kg m-2",direction="INC",           &
                     surface_type=LIS_rc%lsm_index)

! David Mocko (10/29/2019) - Copy RELSMC calculation from Noah-3.X
           do i = 1,NoahmpIO%nsoil
              if (NoahmpIO%relsmc(1,i,1).gt.1.0) then
                 NoahmpIO%relsmc(1,i,1) = 1.0
              endif
              if (NoahmpIO%relsmc(1,i,1).lt.0.01) then
                 NoahmpIO%relsmc(1,i,1) = 0.01
              endif

! J.Case (9/11/2014) -- Set relative soil moisture to missing (LIS_rc%udef)
! if the vegetation type is urban class.
              if (NoahmpIO%ivgtyp(1,1) .eq. NoahmpIO%urban_vegtype(1,1)) then
                 NoahmpIO%relsmc(1,i,1) = LIS_rc%udef
              endif
              call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RELSMC,vlevel=i, &
                       value=NoahmpIO%relsmc(1,i,1),unit='-',direction="-",surface_type=LIS_rc%lsm_index)
              if (NoahmpIO%relsmc(1,i,1) .eq. LIS_rc%udef) then
                 tempval = NoahmpIO%relsmc(1,i,1)
              else
                 tempval = NoahmpIO%relsmc(1,i,1)*100.0
              endif
              call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RELSMC,vlevel=i, &
                       value=tempval,unit='%',direction="-",surface_type=LIS_rc%lsm_index)
            enddo

            ! reset forcing variables to zeros
            NoahMP50_struc(n)%noahmp50(t)%tair   = 0.0
            NoahMP50_struc(n)%noahmp50(t)%psurf  = 0.0
            NoahMP50_struc(n)%noahmp50(t)%wind_e = 0.0
            NoahMP50_struc(n)%noahmp50(t)%wind_n = 0.0
            NoahMP50_struc(n)%noahmp50(t)%qair   = 0.0
            NoahMP50_struc(n)%noahmp50(t)%swdown = 0.0
            NoahMP50_struc(n)%noahmp50(t)%lwdown = 0.0
            NoahMP50_struc(n)%noahmp50(t)%prcp   = 0.0

        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        NoahMP50_struc(n)%forc_count = 0

    endif ! end of alarmCheck loop

    ! EMK...See if noahmp50_struc(n)%noahmp50(t)%tair_agl_min needs to be 
    ! reset for calculating RHMin.  
    alarmCheck = LIS_isAlarmRinging(LIS_rc, &
         "NoahMP50 RHMin alarm "//trim(fnest))
    if (alarmCheck) then
       write(LIS_logunit,*) &
            '[INFO] Resetting tair_agl_min for RHMin calculation'
       do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          noahmp50_struc(n)%noahmp50(t)%tair_agl_min = 999.
       end do
    end if

end subroutine NoahMP50_main
