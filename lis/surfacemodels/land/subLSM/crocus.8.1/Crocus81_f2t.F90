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
! !ROUTINE: Crocus81_f2t
! \label{Crocus81_f2t}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  
!  10/18/19: Mahdi Navari, Shugong Wang; Initial implementation for LIS 7 and Crocus81
!
! !INTERFACE:
subroutine Crocus81_f2t(n)
! !USES:
    use ESMF
    use LIS_coreMod, only       : LIS_rc , LIS_surface
    use LIS_metforcingMod, only : LIS_FORC_State
    use LIS_logMod, only        : LIS_verify
    use LIS_FORC_AttributesMod   
    use Crocus81_lsmMod

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n     
!
! !DESCRIPTION:
!  This routine transfers the LIS provided forcing into the Crocus81
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
 
    ! Near Surface Specific Humidity [kg/kg]
    type(ESMF_Field)  :: q2Field
    real, pointer     :: q2(:)
 
    ! Eastward Wind [m/s]
    type(ESMF_Field)  :: uField
    real, pointer     :: uwind(:)
 
    ! Northward Wind [m/s]
    type(ESMF_Field)  :: vField
    real, pointer     :: vwind(:)
 
    ! Rainfall Rate [kg/(m2s)]
    type(ESMF_Field)  :: pcpField
    real, pointer     :: pcp(:)
 
    ! Snowfall Rate [kg/(m2s)]
    type(ESMF_Field)  :: snowField
    real, pointer     :: snowf(:)
 
    ! Incident Longwave Radiation [W/m2]
    type(ESMF_Field)  :: lwdField
    real, pointer     :: lwd(:)
 
    ! Incident Shortwave Radiation [W/m2]
    type(ESMF_Field)  :: swdField
    real, pointer     :: swd(:)
 
    ! Surface Pressure [Pa]
    type(ESMF_Field)  :: psurfField
    real, pointer     :: psurf(:)

    integer, pointer   :: layer_h(:), layer_m(:)
 
    !!! GET FORCING FIELDS FROM LIS
    ! get near surface air temperature
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Tair%varname(1)), tmpField, rc=status)
    call LIS_verify(status, "Crocus81_f2t: error getting TA")

    call ESMF_FieldGet(tmpField, localDE = 0, farrayPtr = tmp, rc = status)
    call LIS_verify(status, "Crocus81_f2t: error retrieving TA")
    
    ! get near surface specific humidity
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Qair%varname(1)), q2Field, rc=status)
    call LIS_verify(status, "Crocus81_f2t: error getting QA")

    call ESMF_FieldGet(q2Field, localDE = 0, farrayPtr = q2, rc = status)
    call LIS_verify(status, "Crocus81_f2t: error retrieving QA")
    
    ! get eastward wind
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_E%varname(1)), uField, rc=status)
    call LIS_verify(status, "Crocus81_f2t: error getting Wind_E")

    call ESMF_FieldGet(uField, localDE = 0, farrayPtr = uwind, rc = status)
    call LIS_verify(status, "Crocus81_f2t: error retrieving Wind_E")
    
    ! get northward wind
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_N%varname(1)), vField, rc=status)
    call LIS_verify(status, "Crocus81_f2t: error getting Wind_N")

    call ESMF_FieldGet(vField, localDE = 0, farrayPtr = vwind, rc = status)
    call LIS_verify(status, "Crocus81_f2t: error retrieving Wind_N")
    
    ! get rainfall rate
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Rainf%varname(1)), pcpField, rc=status)
    call LIS_verify(status, "Crocus81_f2t: error getting RRSNOW")

    call ESMF_FieldGet(pcpField, localDE = 0, farrayPtr = pcp, rc = status)
    call LIS_verify(status, "Crocus81_f2t: error retrieving RRSNOW")
    
    ! get snowfall rate
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Snowf%varname(1)), snowField, rc=status)
    call LIS_verify(status, "Crocus81_f2t: error getting SRSNOW")

    call ESMF_FieldGet(snowField, localDE = 0, farrayPtr = snowf, rc = status)
    call LIS_verify(status, "Crocus81_f2t: error retrieving SRSNOW")
    
    ! get incident longwave radiation
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_LWdown%varname(1)), lwdField, rc=status)
    call LIS_verify(status, "Crocus81_f2t: error getting LW_RAD")

    call ESMF_FieldGet(lwdField, localDE = 0, farrayPtr = lwd, rc = status)
    call LIS_verify(status, "Crocus81_f2t: error retrieving LW_RAD")
    
    ! get incident shortwave radiation
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_SWdown%varname(1)), swdField, rc=status)
    call LIS_verify(status, "Crocus81_f2t: error getting SW_RAD")

    call ESMF_FieldGet(swdField, localDE = 0, farrayPtr = swd, rc = status)
    call LIS_verify(status, "Crocus81_f2t: error retrieving SW_RAD")
    
    ! get surface pressure
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Psurf%varname(1)), psurfField, rc=status)
    call LIS_verify(status, "Crocus81_f2t: error getting PPS")

    call ESMF_FieldGet(psurfField, localDE = 0, farrayPtr = psurf, rc = status)
    call LIS_verify(status, "Crocus81_f2t: error retrieving PPS")
    
 
    !!! set the forcing counter
    CROCUS81_struc(n)%forc_count = CROCUS81_struc(n)%forc_count + 1
 
    !!! pass forcing data to tiles
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id

        ! TA
        CROCUS81_struc(n)%crocus81(t)%ta = CROCUS81_struc(n)%crocus81(t)%ta + tmp(tid)

        ! QA
        CROCUS81_struc(n)%crocus81(t)%qa = CROCUS81_struc(n)%crocus81(t)%qa + q2(tid)

        ! WIND_E
        CROCUS81_struc(n)%crocus81(t)%wind_e = CROCUS81_struc(n)%crocus81(t)%wind_e + uwind(tid)

        ! WIND_N
        CROCUS81_struc(n)%crocus81(t)%wind_n = CROCUS81_struc(n)%crocus81(t)%wind_n + vwind(tid)

        ! RRSNOW
        CROCUS81_struc(n)%crocus81(t)%rrsnow = CROCUS81_struc(n)%crocus81(t)%rrsnow + pcp(tid)

        ! SRSNOW
        CROCUS81_struc(n)%crocus81(t)%srsnow = CROCUS81_struc(n)%crocus81(t)%srsnow + snowf(tid)

        ! LW_RAD
        CROCUS81_struc(n)%crocus81(t)%lw_rad = CROCUS81_struc(n)%crocus81(t)%lw_rad + lwd(tid)

        ! SW_RAD
        CROCUS81_struc(n)%crocus81(t)%sw_rad = CROCUS81_struc(n)%crocus81(t)%sw_rad + swd(tid)

        ! PPS
        CROCUS81_struc(n)%crocus81(t)%pps = CROCUS81_struc(n)%crocus81(t)%pps + psurf(tid)
    enddo

end subroutine Crocus81_f2t
