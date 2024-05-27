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
!                                   

! !INTERFACE:
subroutine Ac71_ETo_calc(P, Tmax, Tmin, Tdew, ws, Rs, z, lat, eto)
! !USES:
    use LIS_tbotAdjustMod,  only: LIS_tbotTimeUtil
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
    real, intent(in)      :: z ! elevation
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

    ! Psychrometric constant [kPa degC-1]
    gamma = 0.664742 * 0.001 * P

    ! Mean saturation vapour pressure (es) [kPa]
    es = (0.6108 * EXP((17.27 * Tmin) / (237.3 + Tmin)) &
            + 0.6108 * EXP((17.27 * Tmax) / (237.3 + Tmax))) / 2
    ! Actual vapor pressure (ea) [kPa]
    ea = 0.6108 * EXP((17.27 * Tdew) / (237.3 + Tdew))

    ! Slope of saturation vapour pressure curve [kPa degC-1]
    slope = (4098 * (0.6108 * EXP((17.27 * Tmean)/(Tmean + 237.3)))) &
            / (Tmean + 237.3)**2

    ! Extraterrestrial radiation (Ra) [MJ m-2 day-1]
    call LIS_tbotTimeUtil(julian_in,LIS_rc%yr) ! First get day of year
    phi = lat * PI/180 ! get latitude in rad
    dr = 1 + 0.033 * COS(2*PI/365 * julian_in) ! Inverse relative distance Earth-Sun
    delta = 0.409 * SIN(2*PI/365 * julian_in - 1.39) ! Solar declination
    omega = ACOS(-TAN(phi)*TAN(delta))
    Ra = (24 * 60/PI) * 0.082 * dr * (omega*SIN(phi)*SIN(delta) &
         + COS(phi)*COS(delta)*SIN(omega))
    if(Ra.le.0) then ! return 0 ETo for unrealisic numbers
            eto = 0
    endif

    ! Net radiation at crop surface (Rn) [MJ m-2 day-1]
    Rso = (0.75 + 2E-5 * z) * Ra
    Rns = (1 - 0.23) * Rs ! fixed value of 0.23 for albedo
    ratio = Rs/Rso
    if (ratio.gt.1) then
        ratio=1
    endif
    Rnl = 4.903E-9 * (((Tmax + 273.15)**4 + (Tmin + 273.15)**4) / 2.) &
                    * (0.34 - 0.14 * SQRT(ea)) &
                    * (1.35 * ratio - 0.35) ! Net longwave radiation [MJ m-2 day-1]
    Rn = Rns - Rnl

    ! Penman-Monteith --> reference evapotranspiration [mm day-1]
    eto = (0.408 * slope * Rn + (gamma * (900 / (Tmean + 273.15)) * ws * (es - ea))) &
        / (slope + gamma * (1 + 0.34 * ws))
    if(eto.le.0) then
        eto = 0 ! avoid negative values
    endif
    
end subroutine Ac71_ETo_calc


!BOP
!
! !ROUTINE: Ac71_read_Trecord
! \label{Ac71_read_Trecord}
!
! !REVISION HISTORY:
!  24 MAY 2024, Louise Busschaert; initial implementation
!                                   

! !INTERFACE:
subroutine ac71_read_Trecord(n)
! !USES:
    use ESMF
    use LIS_metForcingMod,  only: LIS_get_met_forcing, LIS_FORC_State
    use LIS_timeMgrMod,     only: LIS_timemgr_set, LIS_advance_timestep, &
                                  LIS_update_clock
    use LIS_coreMod,        only: LIS_rc
    use LIS_logMod,         only: LIS_logunit, LIS_verify
    use LIS_FORC_AttributesMod
    use Ac71_lsmMod

    !
    ! !DESCRIPTION: 
    ! 
    !  This subroutine stores the mean temperatures for the ac71 simulation
    !  period required when AquaCrop is run with a crop calibrated in growing
    !  degree days. 
    !
    !  
    !EOP
    implicit none
    integer, intent(in)      :: n

    real, allocatable     :: tmp_arr(:,:,:)
    integer               :: i, j, m, p, status, met_ts
    integer               :: yr_saved, mo_saved, da_saved, hr_saved, mn_saved
    real*8                :: ringtime_saved, time_saved
    real                  :: gmt_tmp
    character(256)        :: spatial_interp_saved

    ! Near Surface Air Temperature [K]
    type(ESMF_Field)  :: tmpField
    real, pointer     :: tmp(:)

    write(LIS_logunit,*) "[INFO] AC71: new simulation period, reading of temperature record..."

    ! Change flags for met forcing reading
    ! Spatial interpolation flag
    do m=1,LIS_rc%nmetforc
        spatial_interp_saved = LIS_rc%met_interp(m)
        call finalmetforc(trim(LIS_rc%metforc(m))//char(0),m) ! finalize and initialize again, risky but only solution to get neighbor
        LIS_rc%met_interp(m) = "neighbor"
        call initmetforc(trim(LIS_rc%metforc(m))//char(0),m)  
    enddo

    ! Save current time
    yr_saved  = LIS_rc%yr
    mo_saved = LIS_rc%mo
    da_saved = LIS_rc%da
    hr_saved = LIS_rc%hr
    mn_saved = LIS_rc%mn
    ! Make sure it is the right met time step
    call LIS_update_clock(LIS_rc%ts)

    met_ts = int(86400./LIS_rc%nts(n))

    allocate(tmp_arr(LIS_rc%gnc(n)*LIS_rc%gnr(n),366,met_ts))

    day_loop: do i=1,366
        do j=1,met_ts
            ! read met forcing
            call LIS_get_met_forcing(n)

            ! Get Tair
            call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Tair%varname(1)), tmpField, rc=status)
            call LIS_verify(status, "Ac71_prep_f: error getting Tair")

            call ESMF_FieldGet(tmpField, localDE = 0, farrayPtr = tmp, rc = status)
            call LIS_verify(status, "Ac71_prep_f: error retrieving Tair")

            ! Store temperatures
            tmp_arr(:,i,j) = tmp

            ! Change LIS time to the next meteo time step
            call LIS_advance_timestep(LIS_rc)
        enddo
        if ((LIS_rc%da.ne.AC71_struc(n)%Sim_AnnualStartDay)&
            .and.(LIS_rc%mo.ne.AC71_struc(n)%Sim_AnnualStartMonth)&
            .and.(i.ne.1)) exit day_loop
    enddo day_loop

    ! allocate Trecord accoridng to domain grid
    allocate(AC71_struc(n)%Trecord(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
    do p=1,LIS_rc%gnc(n)*LIS_rc%gnr(n)
        allocate(AC71_struc(n)%Trecord(p)%Tmax_record(366))
        allocate(AC71_struc(n)%Trecord(p)%Tmin_record(366))
        AC71_struc(n)%Trecord(p)%Tmax_record = maxval(tmp_arr(p,:,:),2)
        AC71_struc(n)%Trecord(p)%Tmin_record = minval(tmp_arr(p,:,:),2)
    enddo

    deallocate(tmp_arr)

    ! Reset all flags
    ! Spatial interpolation flag
    do m=1,LIS_rc%nmetforc
        call finalmetforc(trim(LIS_rc%metforc(m))//char(0),m)
        LIS_rc%met_interp(m) = trim(spatial_interp_saved)
        call initmetforc(trim(LIS_rc%metforc(m))//char(0),m)  
    enddo
    ! Reset actual time
    call LIS_timemgr_set(LIS_rc,yr_saved,mo_saved,da_saved,hr_saved,mn_saved,0,0,0.0)
    LIS_rc%rstflag = 1 ! Needed for met forcings
    write(LIS_logunit,*) "[INFO] AC71: new simulation period, reading of temperature record... Done!"

end subroutine ac71_read_Trecord

end module ac71_prep_f
