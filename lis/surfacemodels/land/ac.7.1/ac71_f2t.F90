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
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  
!  18 JAN 2024, Louise Busschaert; initial implementation for LIS 7 and AC71
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

    integer           :: t, v, status
    integer           :: tid 
 
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
 
    ! Eastward Wind [W m-2]
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

    ! PREC_ac [mm/day]
    type(ESMF_Field)  :: PREC_ac_Field
    real, pointer     :: PREC_ac(:)
 
    ! TMIN_ac [degC]
    type(ESMF_Field)  :: TMIN_ac_Field
    real, pointer     :: TMIN_ac(:)
 
    ! TMAX_ac [degC]
    type(ESMF_Field)  :: TMAX_ac_Field
    real, pointer     :: TMAX_ac(:)
 
    ! ETo_ac [mm/d]
    type(ESMF_Field)  :: ETo_ac_Field
    real, pointer     :: ETo_ac(:)

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

    ! MB: AC_71 

    ! get PREC_ac
    if(LIS_Forc_PREC_ac%selectOpt .eq. 1) then 
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_PREC_AC%varname(1)), PREC_ac_Field, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting PREC_ac")

    call ESMF_FieldGet(PREC_ac_Field, localDE = 0, farrayPtr = PREC_ac, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving PREC_ac")
    endif
    ! get TMIN_ac
    if(LIS_Forc_TMIN_ac%selectOpt .eq. 1) then 
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_TMIN_AC%varname(1)), TMIN_ac_Field, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting TMIN_ac")

    call ESMF_FieldGet(TMIN_ac_Field, localDE = 0, farrayPtr = TMIN_ac, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving TMIN_ac")
    endif
    ! get TMAX_ac
    if(LIS_Forc_TMAX_ac%selectOpt .eq. 1) then 
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_TMAX_AC%varname(1)), TMAX_ac_Field, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting TMAX_ac")

    call ESMF_FieldGet(TMAX_ac_Field, localDE = 0, farrayPtr = TMAX_ac, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving TMAX_ac")
    endif
    ! get ETo_ac
    if(LIS_Forc_ETo_ac%selectOpt .eq. 1) then 
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_ETo_AC%varname(1)), ETo_ac_Field, rc=status)
    call LIS_verify(status, "Ac71_f2t: error getting ETo_ac")

    call ESMF_FieldGet(ETo_ac_Field, localDE = 0, farrayPtr = ETo_ac, rc = status)
    call LIS_verify(status, "Ac71_f2t: error retrieving ETo_ac")
    endif

    !!! set the forcing counter
    AC71_struc(n)%forc_count = AC71_struc(n)%forc_count + 1
 
    !!! pass forcing data to tiles
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id

        ! TAIR
        AC71_struc(n)%ac71(t)%tair = AC71_struc(n)%ac71(t)%tair + tmp(tid)

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
        
        if (trim(LIS_rc%metforc(1)) == 'MERRA2_AC') then
        ! PREC_ac
        AC71_struc(n)%ac71(t)%PREC_ac = AC71_struc(n)%ac71(t)%PREC_ac + PREC_ac(tid)
        ! TMIN
        AC71_struc(n)%ac71(t)%TMIN_ac = AC71_struc(n)%ac71(t)%TMIN_ac + TMIN_ac(tid)
        ! TMAX
        AC71_struc(n)%ac71(t)%TMAX_ac = AC71_struc(n)%ac71(t)%TMAX_ac + TMAX_ac(tid)
        ! ETo
        AC71_struc(n)%ac71(t)%ETo_ac = AC71_struc(n)%ac71(t)%ETo_ac + ETo_ac(tid)
        end if 

    enddo
 
end subroutine Ac71_f2t
