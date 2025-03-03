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
! !ROUTINE: FLake1_f2t
! \label{FLake1_f2t}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  
!  6/4/13: Shugong Wang; Initial implementation for LIS 7 and FLake1
!
! !INTERFACE:
subroutine FLake1_f2t(n)
! !USES:
    use ESMF
    use LIS_coreMod, only       : LIS_rc , LIS_surface
    use LIS_metforcingMod, only : LIS_FORC_State
    use LIS_logMod, only        : LIS_verify
    use LIS_FORC_AttributesMod   
    use FLake1_Mod

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n     
!
! !DESCRIPTION:
!  This routine transfers the LIS provided forcing into the FLake1
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
 
    ! Incident Shortwave Radiation [W/m2]
    type(ESMF_Field)  :: swdField
    real, pointer     :: swd(:)
 
    ! Incident Longwave Radiation [W/m2]
    type(ESMF_Field)  :: lwdField
    real, pointer     :: lwd(:)
 
    ! Eastward Wind [W/m2]
    type(ESMF_Field)  :: uField
    real, pointer     :: uwind(:)
 
    ! Northward Wind [m/s]
    type(ESMF_Field)  :: vField
    real, pointer     :: vwind(:)
 
    ! Surface Pressure [Pa]
    type(ESMF_Field)  :: psurfField
    real, pointer     :: psurf(:)

    !!! GET FORCING FIELDS FROM LIS
    ! get near surface air temperature
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Tair%varname(1)), tmpField, rc=status)
    call LIS_verify(status, "FLake1_f2t: error getting Tair")

    call ESMF_FieldGet(tmpField, localDE = 0, farrayPtr = tmp, rc = status)
    call LIS_verify(status, "FLake1_f2t: error retrieving Tair")
    
    ! get near surface specific humidity
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Qair%varname(1)), q2Field, rc=status)
    call LIS_verify(status, "FLake1_f2t: error getting Qair")

    call ESMF_FieldGet(q2Field, localDE = 0, farrayPtr = q2, rc = status)
    call LIS_verify(status, "FLake1_f2t: error retrieving Qair")
    
    ! get incident shortwave radiation
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_SWdown%varname(1)), swdField, rc=status)
    call LIS_verify(status, "FLake1_f2t: error getting Swdown")

    call ESMF_FieldGet(swdField, localDE = 0, farrayPtr = swd, rc = status)
    call LIS_verify(status, "FLake1_f2t: error retrieving Swdown")
    
    ! get incident longwave radiation
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_LWdown%varname(1)), lwdField, rc=status)
    call LIS_verify(status, "FLake1_f2t: error getting Lwdown")

    call ESMF_FieldGet(lwdField, localDE = 0, farrayPtr = lwd, rc = status)
    call LIS_verify(status, "FLake1_f2t: error retrieving Lwdown")
    
    ! get eastward wind
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_E%varname(1)), uField, rc=status)
    call LIS_verify(status, "FLake1_f2t: error getting Wind_E")

    call ESMF_FieldGet(uField, localDE = 0, farrayPtr = uwind, rc = status)
    call LIS_verify(status, "FLake1_f2t: error retrieving Wind_E")
    
    ! get northward wind
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_N%varname(1)), vField, rc=status)
    call LIS_verify(status, "FLake1_f2t: error getting Wind_N")

    call ESMF_FieldGet(vField, localDE = 0, farrayPtr = vwind, rc = status)
    call LIS_verify(status, "FLake1_f2t: error retrieving Wind_N")
    
    ! get surface pressure
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Psurf%varname(1)), psurfField, rc=status)
    call LIS_verify(status, "FLake1_f2t: error getting Psurf")

    call ESMF_FieldGet(psurfField, localDE = 0, farrayPtr = psurf, rc = status)
    call LIS_verify(status, "FLake1_f2t: error retrieving Psurf")
    
 
    !!! set the forcing counter
    FLAKE1_struc(n)%forc_count = FLAKE1_struc(n)%forc_count + 1
 
    !!! pass forcing data to tiles
    do t=1, LIS_rc%npatch(n, LIS_rc%lake_index)
        tid = LIS_surface(n, LIS_rc%lake_index)%tile(t)%tile_id
        ! TAIR
        FLAKE1_struc(n)%flake1(t)%tair = FLAKE1_struc(n)%flake1(t)%tair + tmp(tid)
        ! QAIR
        FLAKE1_struc(n)%flake1(t)%qair = FLAKE1_struc(n)%flake1(t)%qair + q2(tid)
        ! SWDOWN
        FLAKE1_struc(n)%flake1(t)%swdown = FLAKE1_struc(n)%flake1(t)%swdown + swd(tid)
        ! LWDOWN
        FLAKE1_struc(n)%flake1(t)%lwdown = FLAKE1_struc(n)%flake1(t)%lwdown + lwd(tid)
        ! WIND_E
        FLAKE1_struc(n)%flake1(t)%wind_e = FLAKE1_struc(n)%flake1(t)%wind_e + uwind(tid)
        ! WIND_N
        FLAKE1_struc(n)%flake1(t)%wind_n = FLAKE1_struc(n)%flake1(t)%wind_n + vwind(tid)
        ! PSURF
        FLAKE1_struc(n)%flake1(t)%psurf = FLAKE1_struc(n)%flake1(t)%psurf + psurf(tid)
    enddo
 
end subroutine FLake1_f2t
