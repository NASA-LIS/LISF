!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module ac71_prep_f

contains
!BOP
!
! !ROUTINE: Ac71_ETo_calc
! \label{Ac71_ETo_calc}
!
! !REVISION HISTORY:
!  18 FEB 2024, Louise Busschaert; initial implementation

! !INTERFACE:
subroutine Ac71_ETo_calc(P, Tmax, Tmin, Tdew, ws, Rs, z, lat, eto)
! !USES:
    use LIS_tbotAdjustMod,  only: LIS_tbotTimeUtil ! consider using something else
    use LIS_coreMod,        only: LIS_rc
    use LIS_logMod, only     : LIS_logunit

    !
    ! !DESCRIPTION: 
    ! 
    !  This routine computes the reference evapotranspiration (ETo)
    !  with the Penman-Monteith equation, following the guidelines
    !  of the FAO Irrigation and Drainage Paper No. 56.
    !
    !  
    !EOP
    implicit none
    real, intent(in)      :: P ! kPa
    real, intent(in)      :: Tmax, Tmin, Tdew
    real, intent(in)      :: ws ! wind speed
    real, intent(in)      :: Rs ! radiation
    real, intent(in)      :: z
    real, intent(in)      :: lat
    real, intent(inout)   :: eto ! returns eto [mm day-1]
    
    real                  :: Tmean
    real                  :: gamma, ea, es, slope
    real                  :: phi, dr, delta, omega, Ra, Rso
    real                  :: Rns, Rnl, Rn
    real                  :: julian_in, ratio

    real, parameter       :: PI = 3.14159265359

    ! Teamn (degC)
    Tmean = (Tmin + Tmax)/2.
    write(LIS_logunit, *) "Tmax=", Tmax
    write(LIS_logunit, *) "Tmin=", Tmin

    ! Psychrometric constant [kPa degC-1]
    gamma = 0.664742 * 0.001 * P

    ! Mean saturation vapour pressure (es) [kPa]
    es = (0.6108 * EXP((17.27 * Tmin) / (237.3 + Tmin)) &
            + 0.6108 * EXP((17.27 * Tmax) / (237.3 + Tmax))) / 2
    write(LIS_logunit, *) "es=", es
    ! Actual vapor pressure (ea) [kPa]
    ea = 0.6108 * EXP((17.27 * Tdew) / (237.3 + Tdew))
    write(LIS_logunit, *) "ea=", ea
    write(LIS_logunit, *) "Tdew=", Tdew

    ! Slope of saturation vapour pressure curve [kPa degC-1]
    slope = (4098 * (0.6108 * EXP((17.27 * Tmean)/(Tmean + 237.3)))) &
            / (Tmean + 237.3)**2

    ! Extraterrestrial radiation (Ra) [MJ m-2 day-1]
    call LIS_tbotTimeUtil(julian_in,LIS_rc%yr) ! First get day of year
    phi = lat * PI/180 ! get latitude in rad
    write(LIS_logunit, *) "phi=", phi
    write(LIS_logunit, *) "doy=", julian_in
    dr = 1 + 0.033 * COS(2*PI/365 * julian_in) ! Inverse relative distance Earth-Sun
    delta = 0.409 * SIN(2*PI/365 * julian_in - 1.39) ! Solar declination
    omega = ACOS(-TAN(phi)*TAN(delta))
    Ra = (24 * 60/PI) * 0.082 * dr * (omega*SIN(phi)*SIN(delta) &
         + COS(phi)*COS(delta)*SIN(omega))
    write(LIS_logunit, *) "Ra=", Ra
    if(Ra.le.0) then ! return 0 ETo for unrealisic numbers
            eto = 0
    endif

    ! Net radiation at crop surface (Rn) [MJ m-2 day-1]
    write(LIS_logunit, *) "Rs=", Rs
    Rso = (0.75 + 2E-5 * z) * Ra
    write(LIS_logunit, *) "Rso=", Rso
    Rns = (1 - 0.23) * Rs ! fixed value of 0.23 for albedo
    write(LIS_logunit, *) "Rns=", Rns
    ratio = Rs/Rso
    if (ratio.gt.1) then
        ratio=1
    endif
    write(LIS_logunit, *) "ratio=", ratio
    write(LIS_logunit, *) "f1=", ((Tmax + 273.15)**4 + (Tmin + 273.15)**4) / 2.
    write(LIS_logunit, *) "f2=", 0.34 - 0.14 * SQRT(ea)
    write(LIS_logunit, *) "f3=", (1.35 * ratio - 0.35)
    Rnl = 4.903E-9 * (((Tmax + 273.15)**4 + (Tmin + 273.15)**4) / 2.) &
                    * (0.34 - 0.14 * SQRT(ea)) &
                    * (1.35 * ratio - 0.35) ! Net longwave radiation [MJ m-2 day-1]
    write(LIS_logunit, *) "Rnl=", Rnl
    Rn = Rns - Rnl
    write(LIS_logunit, *) "Rn=", Rn

    ! later --> add if statements to correct for 10m data (ERA5 10m T and 10m wind)
    ! Add in config file reference height for wind

    ! Penman-Monteith --> reference evapotranspiration [mm day-1]
    eto = (0.408 * slope * Rn + (gamma * (900 / (Tmean + 273.15)) * ws * (es - ea))) &
        / (slope + gamma * (1 + 0.34 * ws))
    write(LIS_logunit, *) "ws=", ws
    write(LIS_logunit, *) "eto=", eto
    if(eto.le.0) then
        eto = 0 ! avoid negative values
    endif
    
end subroutine Ac71_ETo_calc

end module ac71_prep_f
