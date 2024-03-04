!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: Ac71_f2t
! \label{Ac71_f2t}
!  
!  18 JAN 2024, Louise Busschaert; initial implementation for LIS 7 and AC71
!  14 FEB 2024, Louise Busschaert; Compute Tmax/Min and Averages
!
! !INTERFACE:
subroutine Ac71_f2t(n)
! !USES:
    use ESMF
    use LIS_coreMod, only       : LIS_rc , LIS_surface
    use LIS_metforcingMod, only : LIS_FORC_State
    use LIS_logMod, only        : LIS_verify
    use LIS_FORC_AttributesMod   
    use Ac71_lsmMod

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n     
!
! !DESCRIPTION:
!  This routine transfers the LIS provided forcing into the Ac71
!  model tiles.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

    integer            :: t, v, status
    integer            :: tid 
    real               :: ee, val, td
 
    ! Near Surface Air Temperature [K]
    type(ESMF_Field)  :: tmpField
    real, pointer     :: tmp(:)
 
    ! Near Surface Specific Humidity [kg kg-1]
    type(ESMF_Field)  :: q2Field
    real, pointer     :: q2(:)
 
    ! Incident Shortwave Radiation [W m-2]
    type(ESMF_Field)  :: swdField
    real, pointer     :: swd(:)
 
    ! Incident Longwave Radiation [W m-2]
    type(ESMF_Field)  :: lwdField
    real, pointer     :: lwd(:)
 
    ! Eastward Wind [m s-1]
    type(ESMF_Field)  :: uField
    real, pointer     :: uwind(:)
 
    ! Northward Wind [m s-1]
    type(ESMF_Field)  :: vField
    real, pointer     :: vwind(:)
 
    ! Surface Pressure [Pa]
    type(ESMF_Field)  :: psurfField
    real, pointer     :: psurf(:)
 
    ! Rainfall Rate [kg m-2 s-1]
    type(ESMF_Field)  :: pcpField
    real, pointer     :: pcp(:)
 
    ! Snowfall Rate [kg m-2 s-1]
    type(ESMF_Field)  :: snowField
    real, pointer     :: snowf(:)

    ! Dewpoint Temperature [K]
    real, pointer     :: tdew(:)

    ! Wind Speed [m/s]
    real, pointer     :: wndspd(:)


    integer, pointer   :: layer_h(:), layer_m(:)
 
    !!! GET FORCING FIELDS FROM LIS
    ! get near surface air temperature
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Tair%varname(1)), tmpField, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting Tair")

    call ESMF_FieldGet(tmpField, localDE = 0, farrayPtr = tmp, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving Tair")
    
    ! get near surface specific humidity
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Qair%varname(1)), q2Field, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting Qair")

    call ESMF_FieldGet(q2Field, localDE = 0, farrayPtr = q2, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving Qair")
    
    ! get incident shortwave radiation
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_SWdown%varname(1)), swdField, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting Swdown")

    call ESMF_FieldGet(swdField, localDE = 0, farrayPtr = swd, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving Swdown")
    
    ! get incident longwave radiation
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_LWdown%varname(1)), lwdField, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting Lwdown")

    call ESMF_FieldGet(lwdField, localDE = 0, farrayPtr = lwd, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving Lwdown")
    
    ! get eastward wind
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_E%varname(1)), uField, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting Wind_E")

    call ESMF_FieldGet(uField, localDE = 0, farrayPtr = uwind, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving Wind_E")
    
    ! get northward wind
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_N%varname(1)), vField, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting Wind_N")

    call ESMF_FieldGet(vField, localDE = 0, farrayPtr = vwind, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving Wind_N")
    
    ! get surface pressure
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Psurf%varname(1)), psurfField, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting Psurf")

    call ESMF_FieldGet(psurfField, localDE = 0, farrayPtr = psurf, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving Psurf")
    
    ! get rainfall rate
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Rainf%varname(1)), pcpField, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting Rainf")

    call ESMF_FieldGet(pcpField, localDE = 0, farrayPtr = pcp, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving Rainf")
    
    ! get snowfall rate
    if(LIS_Forc_Snowf%selectOpt .eq. 1) then 
        call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Snowf%varname(1)), snowField, rc=status)
        call LIS_verify(status, "Ac71_f2t: error getting Snowf")

        call ESMF_FieldGet(snowField, localDE = 0, farrayPtr = snowf, rc = status)
        call LIS_verify(status, "Ac71_f2t: error retrieving Snowf")
    endif 

    !!! set the forcing counter
    AC71_struc(n)%forc_count = AC71_struc(n)%forc_count + 1
 
    !!! pass forcing data to tiles
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id

        ! TAIR
        AC71_struc(n)%ac71(t)%tair = AC71_struc(n)%ac71(t)%tair + tmp(tid)

        if (AC71_struc(n)%forc_count.eq.1) then !First iteration set max/min 
            AC71_struc(n)%ac71(t)%tmax = tmp(tid)
            AC71_struc(n)%ac71(t)%tmin = tmp(tid)
        else
            if (tmp(tid).gt.AC71_struc(n)%ac71(t)%tmax) then
            AC71_struc(n)%ac71(t)%tmax=tmp(tid) !Replace maximum temperature
            endif
            if (tmp(tid).lt.AC71_struc(n)%ac71(t)%tmin) then
            AC71_struc(n)%ac71(t)%tmin=tmp(tid) !Replace minimum temperature
            endif
        endif

        ! QAIR
        AC71_struc(n)%ac71(t)%qair = AC71_struc(n)%ac71(t)%qair + q2(tid)

        ! SWDOWN
        AC71_struc(n)%ac71(t)%swdown = AC71_struc(n)%ac71(t)%swdown + swd(tid)

        ! LWDOWN
        AC71_struc(n)%ac71(t)%lwdown = AC71_struc(n)%ac71(t)%lwdown + lwd(tid)

        ! WIND_E
        AC71_struc(n)%ac71(t)%wind_e = AC71_struc(n)%ac71(t)%wind_e + uwind(tid)

        ! WIND_N
        AC71_struc(n)%ac71(t)%wind_n = AC71_struc(n)%ac71(t)%wind_n + vwind(tid)

        ! Calculate Magnitude of Wind Speed (m/s) 
        AC71_struc(n)%ac71(t)%wndspd = AC71_struc(n)%ac71(t)%wndspd + SQRT(uwind(tid)**2 + vwind(tid)**2)

        ! PSURF
        AC71_struc(n)%ac71(t)%psurf = AC71_struc(n)%ac71(t)%psurf + psurf(tid)

        ! RAINF
        if(pcp(tid).ne.LIS_rc%udef) then
           AC71_struc(n)%ac71(t)%prcp = AC71_struc(n)%ac71(t)%prcp + pcp(tid)
        endif

        ! SNOWF
        ! If there is snowf add it to precipitation.  Ac71 does not use
        ! separate rainf and snowf.  It determines what to do with precipitation.
        if(LIS_Forc_Snowf%selectOpt .eq. 1) then 
           if(snowf(tid).ne.LIS_rc%udef) then
              AC71_struc(n)%ac71(t)%prcp = AC71_struc(n)%ac71(t)%prcp + snowf(tid)
           endif
        endif

        ! Calculate Dewpoint

        ! Following A First Course in Atmospheric Thermodynamics, assume
        ! approximation q = epsilon*e/p

        ! Calculate vapor pressure
        ee = (q2(tid)*psurf(tid))/0.622

        ! Invert Bolton 1980 formula for saturation vapor pressure to calculate Td
        ! since es(Td) = e

        val = log(ee/611.2)
        td = (243.5 * val) / (17.67 - val) ! Dewpoint in C
        td = td + 273.15
        AC71_struc(n)%ac71(t)%tdew = AC71_struc(n)%ac71(t)%tdew + td 

    enddo
 
end subroutine Ac71_f2t
