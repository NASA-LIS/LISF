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
!BOP
!
! !ROUTINE: noahmpglacier3911_main
! \label{noahmpglacier3911_main}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine noahmpglacier3911_main(n)
! !USES:
  use LIS_coreMod
  use LIS_histDataMod
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_constantsMod,  only : LIS_CONST_RHOFW
  use LIS_logMod, only     : LIS_logunit, LIS_endrun
  use LIS_FORC_AttributesMod 
  use noahmpglacier3911_Mod
  use module_sf_noahmpglacier3911, only : calc_declin
  use module_sf_noahmp_glacier
  
  implicit none
! !ARGUMENTS:
  integer, intent(in)  :: n
  integer              :: t
  integer              :: row, col
  integer              :: year, month, day, hour, minute, second
  logical              :: alarmCheck

!
! !DESCRIPTION:
!  This is the entry point for calling the noahmpglacier3911 physics.
!  This routine calls the {\tt noahmp\_driver\_36 } routine that performs the
!  land surface computations, to solve for water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP
  integer        :: i,j,k
  real           :: dt
  real           :: lat, lon
  
  real                :: t_ml,p_ml, u_ml,v_ml, q_ml,swdn
  real                :: prcp,lwdn, tbot,z_ml
  real                :: zsoil(noahmpgl3911_struc(n)%nsoil)
  real                :: cosz
  character(len=12)   :: nowdate
  integer             :: yearlen
  real                :: julian
  integer             :: nsoil, nsnow
  real                :: qsnow
  real                :: sneqvo
  real                :: ficeold(-noahmpgl3911_struc(n)%nsnow+1:0)
  integer             :: snl_idx,isnow
  real                :: albold
  real                :: ch, cm
  real                :: swe,sndpth
  real                :: ponding, ponding1, ponding2
  real                :: t2mb,q2mb
  real                :: chb2,fpice
  real                :: emissi
  real                :: qsnbot
  real                :: fsa
  real                :: fsr
  real                :: fira
  real                :: fsh
  real                :: fgev
  real                :: ssoil
  real                :: trad
  real                :: esoil
  real                :: runsf
  real                :: runsb
  real                :: sag
  real                :: salb
  real                :: tg
  real                :: tauss
  real                :: smc(noahmpgl3911_struc(n)%nsoil)
  real                :: smh2o(noahmpgl3911_struc(n)%nsoil)
  real                :: qsfc1d       ! bulk surface specific humidity
  real                :: snice(-noahmpgl3911_struc(n)%nsnow+1:0)
  real                :: snliq(-noahmpgl3911_struc(n)%nsnow+1:0 )
  real :: zsnso(-noahmpgl3911_struc(n)%nsnow+1:noahmpgl3911_struc(n)%nsoil)
  real :: stc(noahmpgl3911_struc(n)%nsoil + noahmpgl3911_struc(n)%nsnow)
  real                 :: soil_temp(noahmpgl3911_struc(n)%nsoil)   
  real                 :: snow_temp(noahmpgl3911_struc(n)%nsnow)
  real                 :: sldpth(noahmpgl3911_struc(n)%nsoil)
#ifdef WRF_HYDRO
  real                 :: sfcheadrt
#endif

  nsnow = noahmpgl3911_struc(n)%nsnow
  nsoil = noahmpgl3911_struc(n)%nsoil
  
  sldpth = noahmpgl3911_struc(n)%sldpth(:)
 ! set ZSOIL 
  zsoil(1) = -noahmpgl3911_struc(n)%sldpth(1)
  do k = 2, nsoil
     zsoil(k) = zsoil(k-1) - noahmpgl3911_struc(n)%sldpth(k)
  enddo

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "noahmpglacier3911 model alarm")

  if (alarmCheck) Then

     do t = 1, LIS_rc%npatch(n, LIS_rc%glacier_index)
        dt = LIS_rc%ts
        row = LIS_surface(n, LIS_rc%glacier_index)%tile(t)%row
        col = LIS_surface(n, LIS_rc%glacier_index)%tile(t)%col
        lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
        lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
        
        ! cosz, yearlen, and julian are calculated in subroutine calc_declin_36 
        ! be careful here!!!, LIS uses GMT. the date for calc_declin_36 should be local time. Longitude is need to cnvert GMT 
          ! into local time!!!. To be implemented after bondville test 

        write(nowdate,'(I4.4,4I2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da,&
             LIS_rc%hr, LIS_rc%mn
        call calc_declin(nowdate(1:4)//"-"//nowdate(5:6)//"-"//&
             nowdate(7:8)//"_"//nowdate(9:10)//":"//nowdate(11:12)//":00", &
             lat, lon, cosz, yearlen, julian)
        
        ! retrieve forcing data from NOAHMPGL_struc(n)%noahmpgl(t) and assign to local variables
        ! tair: air temperature
        t_ml       = noahmpgl3911_struc(n)%noahmpgl(t)%tair   / Noahmpgl3911_struc(n)%forc_count
        Noahmpgl3911_struc(n)%noahmpgl(t)%sfctmp = t_ml
        ! psurf: air pressure
        p_ml      = noahmpgl3911_struc(n)%noahmpgl(t)%psurf / &
             Noahmpgl3911_struc(n)%forc_count

        ! wind_e: eastward wind speed
        u_ml     = Noahmpgl3911_struc(n)%noahmpgl(t)%wind_e / &
             Noahmpgl3911_struc(n)%forc_count

        ! wind_n: northward wind speed
        v_ml     = Noahmpgl3911_struc(n)%noahmpgl(t)%wind_n / &
             Noahmpgl3911_struc(n)%forc_count

        ! qair: near Surface Specific Humidity
        q_ml       = Noahmpgl3911_struc(n)%noahmpgl(t)%qair / &
             Noahmpgl3911_struc(n)%forc_count

        ! swdown: downward solar radiation
        swdn     = Noahmpgl3911_struc(n)%noahmpgl(t)%swdown / &
             Noahmpgl3911_struc(n)%forc_count

        ! lwdown: downward longwave radiation
        lwdn     = Noahmpgl3911_struc(n)%noahmpgl(t)%lwdown / &
             Noahmpgl3911_struc(n)%forc_count

        ! prcp: precipitation Rate
        prcp       = Noahmpgl3911_struc(n)%noahmpgl(t)%prcp / &
             Noahmpgl3911_struc(n)%forc_count

        ! check validity of tair
        if(t_ml .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable tair in Noahmpgl"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of psurf
        if(p_ml .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable psurf in Noahmpgl"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of wind_e
        if(u_ml .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable wind_e in Noahmpgl"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of wind_n
        if(v_ml .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable wind_n in Noahmpgl"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of qair
        if(q_ml .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable qair in Noahmpgl"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of swdown
        if(swdn .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable swdown in Noahmpgl"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of lwdown
        if(lwdn .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable lwdown in Noahmpgl"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif
        ! check validity of prcp
        if(prcp .eq. LIS_rc%udef) then
           write(LIS_logunit, *) "undefined value found for forcing variable prcp in Noahmpgl"
           write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
           call LIS_endrun()
        endif

        isnow = noahmpgl3911_struc(n)%noahmpgl(t)%isnow

        z_ml = noahmpgl3911_struc(n)%noahmpgl(t)%zlvl
        tbot = noahmpgl3911_struc(n)%noahmpgl(t)%tbot

        tbot = min(tbot,263.15)    ! set deep temp to at most -10C

        qsnow  = noahmpgl3911_struc(n)%noahmpgl(t)%qsnow
        sneqvo = noahmpgl3911_struc(n)%noahmpgl(t)%sneqvo

        ficeold = 0.0
        
        zsnso(-nsnow+1:nsoil) = noahmpgl3911_struc(n)%noahmpgl(t)%zss(1:nsnow+nsoil) 
        snice(-nsnow+1:0)     = noahmpgl3911_struc(n)%noahmpgl(t)%snowice(1:nsnow)
        snliq(-nsnow+1:0)     = noahmpgl3911_struc(n)%noahmpgl(t)%snowliq(1:nsnow) 

        do snl_idx=isnow+1,0
           if(snice(snl_idx)+&
                snliq(snl_idx)>0.0) then
              ficeold(snl_idx)  =&
                   snice(snl_idx) / &
                   (snice(snl_idx)+&
                   snliq(snl_idx))
           else 
              ficeold(snl_idx)  = 0.0
           endif
        enddo

        albold = noahmpgl3911_struc(n)%noahmpgl(t)%albold
        ch = noahmpgl3911_struc(n)%noahmpgl(t)%ch
        cm = noahmpgl3911_struc(n)%noahmpgl(t)%cm
        swe = noahmpgl3911_struc(n)%noahmpgl(t)%sneqv
        sndpth = noahmpgl3911_struc(n)%noahmpgl(t)%snowh
        tg = noahmpgl3911_struc(n)%noahmpgl(t)%tg
        tauss = noahmpgl3911_struc(n)%noahmpgl(t)%tauss
        smc = noahmpgl3911_struc(n)%noahmpgl(t)%smc
        smh2o = noahmpgl3911_struc(n)%noahmpgl(t)%sh2o
        qsfc1d = q_ml
        stc = noahmpgl3911_struc(n)%noahmpgl(t)%sstc
#ifdef WRF_HYDRO
        sfcheadrt = noahmpgl3911_struc(n)%noahmpgl(t)%sfcheadrt
#endif

        call noahmp_glacier(  lis_localpet, t, &
             noahmpgl3911_struc(n)%alb_opt,&
             noahmpgl3911_struc(n)%snf_opt,&
             noahmpgl3911_struc(n)%tbot_opt,&
             noahmpgl3911_struc(n)%stc_opt,&
             noahmpgl3911_struc(n)%gla_opt,&
             cosz,nsnow,nsoil, dt, & 
             t_ml, p_ml, u_ml,  v_ml, q_ml,  swdn, & 
             prcp, lwdn, tbot,  z_ml, ficeold,  zsoil, & 
             qsnow, sneqvo, albold, cm, ch, isnow, & 
             swe,     smc,   zsnso,  sndpth,   snice,   snliq, & ! IN/OUT :
             tg,     stc,   smh2o,   tauss,  qsfc1d,          & ! IN/OUT :
             fsa,  fsr, fira, fsh, fgev, ssoil, & 
             trad, esoil,runsf,runsb,sag, salb, & 
             qsnbot,ponding,ponding1,ponding2, t2mb, q2mb, & 
             emissi, fpice,  chb2 &                          
#ifdef WRF_HYDRO
        , sfcheadrt                                           &
#endif
        )

        noahmpgl3911_struc(n)%noahmpgl(t)%isnow = isnow 

        noahmpgl3911_struc(n)%noahmpgl(t)%zlvl = z_ml
        noahmpgl3911_struc(n)%noahmpgl(t)%tbot = tbot 

        noahmpgl3911_struc(n)%noahmpgl(t)%qsnow = qsnow  
        noahmpgl3911_struc(n)%noahmpgl(t)%sneqvo = sneqvo 

        
        noahmpgl3911_struc(n)%noahmpgl(t)%zss(1:nsnow+nsoil)  = zsnso(-nsnow+1:nsoil) 
        noahmpgl3911_struc(n)%noahmpgl(t)%snowice(1:nsnow) = snice(-nsnow+1:0)    
        noahmpgl3911_struc(n)%noahmpgl(t)%snowliq(1:nsnow)  = snliq(-nsnow+1:0)     

        noahmpgl3911_struc(n)%noahmpgl(t)%albold= albold 
        noahmpgl3911_struc(n)%noahmpgl(t)%ch = ch 
        noahmpgl3911_struc(n)%noahmpgl(t)%cm= cm 
        noahmpgl3911_struc(n)%noahmpgl(t)%sneqv = swe 
        noahmpgl3911_struc(n)%noahmpgl(t)%snowh = sndpth 
        noahmpgl3911_struc(n)%noahmpgl(t)%tg = tg 
        noahmpgl3911_struc(n)%noahmpgl(t)%tauss = tauss 
        noahmpgl3911_struc(n)%noahmpgl(t)%smc = smc 
        noahmpgl3911_struc(n)%noahmpgl(t)%sh2o = smh2o 
        noahmpgl3911_struc(n)%noahmpgl(t)%sstc = stc 

        noahmpgl3911_struc(n)%noahmpgl(t)%fsa = fsa
        noahmpgl3911_struc(n)%noahmpgl(t)%fsr = fsr
        noahmpgl3911_struc(n)%noahmpgl(t)%fira = fira
        noahmpgl3911_struc(n)%noahmpgl(t)%fsh = fsh
        noahmpgl3911_struc(n)%noahmpgl(t)%fgev = fgev
        noahmpgl3911_struc(n)%noahmpgl(t)%ssoil = ssoil
        noahmpgl3911_struc(n)%noahmpgl(t)%trad = trad
        noahmpgl3911_struc(n)%noahmpgl(t)%edir = esoil
        noahmpgl3911_struc(n)%noahmpgl(t)%runsrf = runsf
        noahmpgl3911_struc(n)%noahmpgl(t)%runsub = runsb
        noahmpgl3911_struc(n)%noahmpgl(t)%sag = sag
        noahmpgl3911_struc(n)%noahmpgl(t)%albedo = salb
        noahmpgl3911_struc(n)%noahmpgl(t)%ponding = ponding
        noahmpgl3911_struc(n)%noahmpgl(t)%ponding1 = ponding1
        noahmpgl3911_struc(n)%noahmpgl(t)%ponding2 = ponding2
        noahmpgl3911_struc(n)%noahmpgl(t)%fpice = fpice
        noahmpgl3911_struc(n)%noahmpgl(t)%emissi = emissi
        noahmpgl3911_struc(n)%noahmpgl(t)%chb2 = chb2
#ifdef WRF_HYDRO
        noahmpgl3911_struc(n)%noahmpgl(t)%sfcheadrt = sfcheadrt
#endif

        ![ 1] output variable: soil_temp (unit=K). ***  soil layer temperature
        soil_temp(1:noahmpgl3911_struc(n)%nsoil) = &
             noahmpgl3911_struc(n)%noahmpgl(t)%sstc(&
             noahmpgl3911_struc(n)%nsnow+1:noahmpgl3911_struc(n)%nsoil+noahmpgl3911_struc(n)%nsnow)
        do i=1, noahmpgl3911_struc(n)%nsoil
           call LIS_diagnoseSurfaceOutputVar(n, t, &
                LIS_MOC_SOILTEMP, value = soil_temp(i), &
                vlevel=i, unit="K", direction="-", &
                surface_type = LIS_rc%glacier_index)
        end do
        
        ![ 2] output variable: snow_temp (unit=K). ***  snow layer temperature
        snow_temp(1:noahmpgl3911_struc(n)%nsnow) = &
             noahmpgl3911_struc(n)%noahmpgl(t)%sstc(1:noahmpgl3911_struc(n)%nsnow) 
        do i=1, noahmpgl3911_struc(n)%nsnow
           ! Test code to reset snow temperature to undefined
           ! when there is no corresponding snow layer - Mocko
           if ((i + abs(noahmpgl3911_struc(n)%noahmpgl(t)%isnow))     &
                .le.noahmpgl3911_struc(n)%nsnow) then
              snow_temp(i) = LIS_rc%udef
           endif
           call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWTPROF, &
                value = snow_temp(i), &
                vlevel=i, unit="K", direction="-", &
                surface_type = LIS_rc%glacier_index)
        end do

        ![ 3] output variable: sh2o (unit=m^3 m-3). ***  volumetric liquid soil moisture 
        do i=1, noahmpgl3911_struc(n)%nsoil
           call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMLIQFRAC, &
                value = noahmpgl3911_struc(n)%noahmpgl(t)%sh2o(i), &
                vlevel=i, unit="m^3 m-3", direction="-", &
                surface_type = LIS_rc%glacier_index)
        end do
        
        ![ 4] output variable: smc (unit=m^3 m-3 ). ***  volumetric soil moisture, ice + liquid 
        do i=1, noahmpgl3911_struc(n)%nsoil
           call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, &
                value = noahmpgl3911_struc(n)%noahmpgl(t)%smc(i), &
                vlevel=i, unit="m^3 m-3", direction="-", &
                surface_type = LIS_rc%glacier_index)
           call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, &
                value = noahmpgl3911_struc(n)%noahmpgl(t)%smc(i)*&
                sldpth(i)*LIS_CONST_RHOFW,  &
                vlevel=i, unit="kg m-2", direction="-", &
                surface_type = LIS_rc%glacier_index)
        enddo
        ![ 11] output variable: tg (unit=K). ***  ground averaged temperature
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GROUNDAVGT, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%tg, &
             vlevel=1, unit="K", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 12] output variable: isnow (unit=-). ***  actual number of snow layers 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOWN_NLAYER, &
             value = -1.0*noahmpgl3911_struc(n)%noahmpgl(t)%isnow, &
             vlevel=1, unit="-", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 13] output variable: z_snow (unit=m). ***  snow layer-bottom depth from snow surface
        do i=1, noahmpgl3911_struc(n)%nsnow
           call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW_LBDFSS, &
                value = noahmpgl3911_struc(n)%noahmpgl(t)%zss(i), &
                vlevel=i, unit="m", direction="-", &
                surface_type = LIS_rc%glacier_index)
        end do
        ![ 14] output variable: z_soil (unit=m). ***  soil layer-bottom depth from snow surface
        do i=1, noahmpgl3911_struc(n)%nsoil
           call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOIL_LBDFSS, &
                value = noahmpgl3911_struc(n)%noahmpgl(t)%zss(&
                i+noahmpgl3911_struc(n)%nsnow), &
                vlevel=i, unit="m", direction="-", &
                surface_type = LIS_rc%glacier_index)
        end do
        
        ![ 15] output variable: snowh (unit=m ). ***  snow height 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDEPTH, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%snowh, &
             vlevel=1, unit="m", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 16] output variable: sneqv (unit=mm ). ***  snow water equivalent 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWE, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%sneqv, &
             vlevel=1, unit="kg m-2", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 17] output variable: snowice (unit=mm ). ***  snow-layer ice 
        do i=1, noahmpgl3911_struc(n)%nsnow
           call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWICE, &
                value = noahmpgl3911_struc(n)%noahmpgl(t)%snowice(i), &
                vlevel=i, unit="mm", direction="-", &
                surface_type = LIS_rc%glacier_index)
        end do
            
        ![ 18] output variable: snowliq (unit=mm ). ***  snow-layer liquid water 
        do i=1, noahmpgl3911_struc(n)%nsnow
           call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWLIQ, &
                value = noahmpgl3911_struc(n)%noahmpgl(t)%snowliq(i), &
                vlevel=i, unit="mm", direction="-", &
                surface_type = LIS_rc%glacier_index)
        end do
        
        ![ 31] output variable: cm (unit=m s-1 ). ***  momentum drag coefficient 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CM, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%cm, &
             vlevel=1, unit="m s-1", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 32] output variable: ch (unit=m s-1 ). ***  sensible heat exchange coefficient 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CH,&
             value = noahmpgl3911_struc(n)%noahmpgl(t)%ch, &
             vlevel=1, unit="m s-1", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 33] output variable: tauss (unit=- ). ***  snow aging term 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWAGE, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%tauss, &
             vlevel=1, unit="-", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 37] output variable: fsa (unit=W m-2). ***  total absorbed solar radiation 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWNET, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%fsa, &
             vlevel=1, unit="W m-2", direction="DN", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 38] output variable: fsr (unit=W m-2 ). ***  total reflected solar radiation 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_FSR, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%fsr, &
             vlevel=1, unit="W m-2", direction="UP", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 39] output variable: fira (unit=W m-2 ). ***  total net longwave radiation to atmosphere 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LWUP, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%fira, &
             vlevel=1, unit="W m-2", direction="UP", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 40] output variable: fsh (unit=W m-2). ***  total sensible heat to atmosphere 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QH, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%fsh, &
             vlevel=1, unit="W m-2", direction="UP", &
             surface_type = LIS_rc%glacier_index)
        ![ 43] output variable: fgev (unit=W m-2  ). ***  ground evaporative heat to atmosphere 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_FGEV, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%fgev, &
             vlevel=1, unit="W m-2", direction="UP", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 41] output variable: ssoil (unit=W m-2). ***  ground heat flux to soil 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QG, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%ssoil, &
             vlevel=1, unit="W m-2", direction="DN", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 48] output variable: trad (unit=K  ). ***  surface radiative temperature 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RADT, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%trad, &
             vlevel=1, unit="K", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 47] output variable: edir (unit=kg m-2 s-1 ). ***  direct evaporation rate from surface 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ESOIL, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%edir, &
             vlevel=1, unit="kg m-2 s-1", direction="UP", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 55] output variable: runsrf (unit=kg m-2 s-1). ***  surface runoff 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QS, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%runsrf, &
             vlevel=1, unit="kg m-2 s-1", direction="OUT", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 56] output variable: runsub (unit=kg m-2 s-1 ). ***  baseflow (saturation excess) 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSB, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%runsub, &
             vlevel=1, unit="kg m-2 s-1", direction="OUT", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 60] output variable: sag (unit=W m-2 ). ***  solar radiation absorbed by ground 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SAG,&
             value = noahmpgl3911_struc(n)%noahmpgl(t)%sag, &
             vlevel=1, unit="W m-2", direction="IN", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 66] output variable: albedo (unit=- ). ***  surface albedo 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%albedo, &
             vlevel=1, unit="-", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 68] output variable: ponding (unit=mm). ***  surface ponding 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PONDING, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%ponding, &
             vlevel=1, unit="mm", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 69] output variable: ponding1 (unit=mm). ***  surface ponding1 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PONDING1, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%ponding1, &
             vlevel=1, unit="mm", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 70] output variable: ponding2 (unit=mm ). ***  surface ponding2 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PONDING2, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%ponding2, &
             vlevel=1, unit="mm", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 94] output variable: fpice (unit=- ). ***  snow fraction in precipitation 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_FPICE, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%fpice, &
             vlevel=1, unit="-", direction="-", &
             surface_type = LIS_rc%glacier_index)
        
        ![ 77] output variable: emissi (unit=- ). ***  surface emissivity 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EMISSFORC, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%emissi, &
             vlevel=1, unit="-", direction="-",&
             surface_type = LIS_rc%glacier_index)
        
        
        ![ 93] output variable: chb2 (unit=m s-1). ***  sensible heat exchange coefficient over bare-ground 
        call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHB2, &
             value = noahmpgl3911_struc(n)%noahmpgl(t)%chb2, &
             vlevel=1, unit="m s-1", direction="-", &
             surface_type = LIS_rc%glacier_index)
     end do
  endif
end subroutine noahmpglacier3911_main


