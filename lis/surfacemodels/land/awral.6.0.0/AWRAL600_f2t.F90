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
! !ROUTINE: AWRAL600_f2t
! \label{AWRAL600_f2t}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  
!  12/18/18: Wendy Sharples, Shugong Wang; Initial implementation for LIS 7 and AWRAL600
!
! !INTERFACE:
subroutine AWRAL600_f2t(n)
! !USES:
    use ESMF
    use LIS_coreMod, only       : LIS_rc , LIS_surface
    use LIS_metforcingMod, only : LIS_FORC_State
    use LIS_logMod, only        : LIS_verify
    use LIS_FORC_AttributesMod   
    use AWRAL600_lsmMod

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n     
!
! !DESCRIPTION:
!  This routine transfers the LIS provided forcing into the AWRAL600
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
 
    ! Near Surface Air Temperature [K]           ==> air temperature (tat, C)
    type(ESMF_Field)  :: tmpField
    real, pointer     :: tmp(:)
 
    ! Near Surface Specific Humidity [kg/kg]     ==> average vapor pressure (avpt, Pa) 
    type(ESMF_Field)  :: q2Field
    real, pointer     :: q2(:)
 
    ! Incident Shortwave Radiation [W/m2]a       ==> shortwave radiation (rgt, MJ/m2)
    type(ESMF_Field)  :: swdField
    real, pointer     :: swd(:)
 
    ! Incident Direct Shortwave Radiation [W/m2] ==> clear sky radiation (radcskyt, MJ/m2)
    type(ESMF_Field)  :: swdirField
    real, pointer     :: swdir(:)
 
    ! Eastward Wind [W/m2]                       ==> wind speed (u2t, m/s)
    type(ESMF_Field)  :: uField
    real, pointer     :: uwind(:)
 
    ! Rainfall Rate [kg/m2s]                     ==> rainfalll (pt, mm)
    type(ESMF_Field)  :: pcpField
    real, pointer     :: pcp(:)

    integer, pointer   :: layer_h(:), layer_m(:)
 
    !!! GET FORCING FIELDS FROM LIS
    ! get near surface air temperature
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Tair%varname(1)), tmpField, rc=status)
    call LIS_verify(status, "AWRAL600_f2t: error getting Tair")

    call ESMF_FieldGet(tmpField, localDE = 0, farrayPtr = tmp, rc = status)
    call LIS_verify(status, "AWRAL600_f2t: error retrieving Tair")
    
    ! get near surface specific humidity
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Qair%varname(1)), q2Field, rc=status)
    call LIS_verify(status, "AWRAL600_f2t: error getting Qair")

    call ESMF_FieldGet(q2Field, localDE = 0, farrayPtr = q2, rc = status)
    call LIS_verify(status, "AWRAL600_f2t: error retrieving Qair")
    
    ! get incident shortwave radiation
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_SWdown%varname(1)), swdField, rc=status)
    call LIS_verify(status, "AWRAL600_f2t: error getting Swdown")

    call ESMF_FieldGet(swdField, localDE = 0, farrayPtr = swd, rc = status)
    call LIS_verify(status, "AWRAL600_f2t: error retrieving Swdown")
    
    ! get incident direct shortwave radiation
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_SWdirect%varname(1)), swdirField, rc=status)
    call LIS_verify(status, "AWRAL600_f2t: error getting Swdirect")

    call ESMF_FieldGet(swdirField, localDE = 0, farrayPtr = swdir, rc = status)
    call LIS_verify(status, "AWRAL600_f2t: error retrieving Swdirect")
    
    ! get eastward wind
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_E%varname(1)), uField, rc=status)
    call LIS_verify(status, "AWRAL600_f2t: error getting Wind_E")

    call ESMF_FieldGet(uField, localDE = 0, farrayPtr = uwind, rc = status)
    call LIS_verify(status, "AWRAL600_f2t: error retrieving Wind_E")
    
    ! get rainfall rate
    call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Rainf%varname(1)), pcpField, rc=status)
    call LIS_verify(status, "AWRAL600_f2t: error getting Rainf")

    call ESMF_FieldGet(pcpField, localDE = 0, farrayPtr = pcp, rc = status)
    call LIS_verify(status, "AWRAL600_f2t: error retrieving Rainf")

    !!! set the forcing counter
    AWRAL600_struc(n)%forc_count = AWRAL600_struc(n)%forc_count + 1
 
    !!! pass forcing data to tiles
    do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id

        ! TAIR
        AWRAL600_struc(n)%awral600(t)%tair = AWRAL600_struc(n)%awral600(t)%tair + tmp(tid)

        ! QAIR
        AWRAL600_struc(n)%awral600(t)%qair = AWRAL600_struc(n)%awral600(t)%qair + q2(tid)

        ! SWDOWN
        AWRAL600_struc(n)%awral600(t)%swdown = AWRAL600_struc(n)%awral600(t)%swdown + swd(tid)

        ! SWDIRECT
        AWRAL600_struc(n)%awral600(t)%swdirect = AWRAL600_struc(n)%awral600(t)%swdirect + swdir(tid)

        ! WIND_E
        AWRAL600_struc(n)%awral600(t)%wind_e = AWRAL600_struc(n)%awral600(t)%wind_e + uwind(tid)

        ! RAINF
        AWRAL600_struc(n)%awral600(t)%rainf = AWRAL600_struc(n)%awral600(t)%rainf + pcp(tid)
    enddo
 
end subroutine AWRAL600_f2t
