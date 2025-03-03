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
! !ROUTINE: AC72_f2t
! \label{AC72_f2t}
!
!  04 NOV 2024, Louise Busschaert; initial implementation for AC72
!
! !INTERFACE:
subroutine AC72_f2t(n)
  ! !USES:
  use AC72_lsmMod
  use ESMF
  use LIS_constantsMod, only: LIS_CONST_G, LIS_CONST_TKFRZ
  use LIS_coreMod, only       : LIS_rc, LIS_surface
  use LIS_FORC_AttributesMod
  use LIS_logMod, only        : LIS_verify
  use LIS_metforcingMod, only : LIS_FORC_State

  implicit none

  ! !ARGUMENTS:
  integer, intent(in) :: n
  !
  ! !DESCRIPTION:
  !  This routine transfers the LIS provided forcing into the AC72
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

  ! For lapse rate
  real :: force_tmp,force_hum,force_prs
  real :: elevdiff
  real :: esat,qsat,rh,fesat,fqsat
  real :: tcforce,pcforce,hcforce,tbar
  real, parameter :: rdry = 287.
  real, parameter :: lapse = -0.0065

  ! For corrected wind speed
  real :: wind_tmp

  ! Near Surface Air Temperature [K]
  type(ESMF_Field)  :: tmpField
  real, pointer     :: tmp(:)

  ! Near Surface Specific Humidity [kg kg-1]
  type(ESMF_Field)  :: q2Field
  real, pointer     :: q2(:)

  ! Incident Shortwave Radiation [W m-2]
  type(ESMF_Field)  :: swdField
  real, pointer     :: swd(:)

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


!!! GET FORCING FIELDS FROM LIS
  ! get near surface air temperature
  call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Tair%varname(1)), tmpField, rc=status)
  call LIS_verify(status, "AC72_f2t: error getting Tair")

  call ESMF_FieldGet(tmpField, localDE = 0, farrayPtr = tmp, rc = status)
  call LIS_verify(status, "AC72_f2t: error retrieving Tair")

  ! get near surface specific humidity
  call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Qair%varname(1)), q2Field, rc=status)
  call LIS_verify(status, "AC72_f2t: error getting Qair")

  call ESMF_FieldGet(q2Field, localDE = 0, farrayPtr = q2, rc = status)
  call LIS_verify(status, "AC72_f2t: error retrieving Qair")

  ! get incident shortwave radiation
  call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_SWdown%varname(1)), swdField, rc=status)
  call LIS_verify(status, "AC72_f2t: error getting Swdown")

  call ESMF_FieldGet(swdField, localDE = 0, farrayPtr = swd, rc = status)
  call LIS_verify(status, "AC72_f2t: error retrieving Swdown")

  ! get eastward wind
  call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_E%varname(1)), uField, rc=status)
  call LIS_verify(status, "AC72_f2t: error getting Wind_E")

  call ESMF_FieldGet(uField, localDE = 0, farrayPtr = uwind, rc = status)
  call LIS_verify(status, "AC72_f2t: error retrieving Wind_E")

  ! get northward wind
  call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Wind_N%varname(1)), vField, rc=status)
  call LIS_verify(status, "AC72_f2t: error getting Wind_N")

  call ESMF_FieldGet(vField, localDE = 0, farrayPtr = vwind, rc = status)
  call LIS_verify(status, "AC72_f2t: error retrieving Wind_N")

  ! get surface pressure
  call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Psurf%varname(1)), psurfField, rc=status)
  call LIS_verify(status, "AC72_f2t: error getting Psurf")

  call ESMF_FieldGet(psurfField, localDE = 0, farrayPtr = psurf, rc = status)
  call LIS_verify(status, "AC72_f2t: error retrieving Psurf")

  ! get rainfall rate
  call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Rainf%varname(1)), pcpField, rc=status)
  call LIS_verify(status, "AC72_f2t: error getting Rainf")

  call ESMF_FieldGet(pcpField, localDE = 0, farrayPtr = pcp, rc = status)
  call LIS_verify(status, "AC72_f2t: error retrieving Rainf")

  ! get snowfall rate
  if(LIS_Forc_Snowf%selectOpt .eq. 1) then 
     call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Snowf%varname(1)), snowField, rc=status)
     call LIS_verify(status, "AC72_f2t: error getting Snowf")

     call ESMF_FieldGet(snowField, localDE = 0, farrayPtr = snowf, rc = status)
     call LIS_verify(status, "AC72_f2t: error retrieving Snowf")
  endif

!!! set the forcing counter
  AC72_struc(n)%forc_count = AC72_struc(n)%forc_count + 1

!!! pass forcing data to tiles
  do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id

     ! lapse-rate correction for ref height > 2 m
     ! Copied from lapse-rate correction of met forcings assuming lapse
     ! rate of 6.5 degC/km
     elevdiff = AC72_struc(n)%refz_tq - 2. ! ref 2m for AC
     if (elevdiff.gt.0.1) then          
        force_tmp = tmp(tid)
        force_hum = q2(tid)
        force_prs = psurf(tid)

        tcforce=force_tmp+(lapse*elevdiff)
        tbar=(force_tmp+tcforce)/2.
        pcforce=force_prs/(exp((LIS_CONST_G*elevdiff)/(rdry*tbar)))
        if (force_hum .eq. 0) force_hum=1e-08
        ee=(force_hum*force_prs)/0.622               
        esat=611.2*(exp((17.67*(force_tmp-LIS_CONST_TKFRZ))/&
             ((force_tmp-LIS_CONST_TKFRZ)+243.5)))
        qsat=(0.622*esat)/(force_prs-(0.378*esat))
        rh=(force_hum/qsat)*100.
        fesat=611.2*(exp((17.67*(tcforce-LIS_CONST_TKFRZ))/ &
             ((tcforce-LIS_CONST_TKFRZ)+243.5)))
        fqsat=(0.622*fesat)/(pcforce-(0.378*fesat))
        hcforce=(rh*fqsat)/100.

        tmp(tid)   = tcforce
        q2(tid)    = hcforce
        psurf(tid) = pcforce
     endif

     ! TAIR
     AC72_struc(n)%ac72(t)%tair = AC72_struc(n)%ac72(t)%tair + tmp(tid)

     if (AC72_struc(n)%forc_count.eq.1) then !First iteration set max/min 
        AC72_struc(n)%ac72(t)%tmax = tmp(tid)
        AC72_struc(n)%ac72(t)%tmin = tmp(tid)
     else
        if (tmp(tid).gt.AC72_struc(n)%ac72(t)%tmax) then
           AC72_struc(n)%ac72(t)%tmax=tmp(tid) !Replace maximum temperature
        endif
        if (tmp(tid).lt.AC72_struc(n)%ac72(t)%tmin) then
           AC72_struc(n)%ac72(t)%tmin=tmp(tid) !Replace minimum temperature
        endif
     endif

     ! SWDOWN
     AC72_struc(n)%ac72(t)%swdown = AC72_struc(n)%ac72(t)%swdown + swd(tid)

     ! Calculate Magnitude of Wind Speed (m/s)
     wind_tmp = SQRT(uwind(tid)**2 + vwind(tid)**2)
     elevdiff = AC72_struc(n)%refz_uv - 2.
     if (elevdiff.gt.0.1) then ! replace with corrected value
        wind_tmp = wind_tmp * (4.87/LOG(67.8*AC72_struc(n)%refz_uv-5.42))
     endif

     AC72_struc(n)%ac72(t)%wndspd = AC72_struc(n)%ac72(t)%wndspd + wind_tmp

     ! PSURF
     AC72_struc(n)%ac72(t)%psurf = AC72_struc(n)%ac72(t)%psurf + psurf(tid)

     ! RAINF
     if(pcp(tid).ne.LIS_rc%udef) then
        AC72_struc(n)%ac72(t)%prcp = AC72_struc(n)%ac72(t)%prcp + pcp(tid)
     endif

     ! SNOWF
     ! If there is snowf add it to precipitation.  AC72 does not use
     ! separate rainf and snowf.
     if(LIS_Forc_Snowf%selectOpt .eq. 1) then 
        if(snowf(tid).ne.LIS_rc%udef) then
           AC72_struc(n)%ac72(t)%prcp = AC72_struc(n)%ac72(t)%prcp + snowf(tid)
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
     td = td + LIS_CONST_TKFRZ
     AC72_struc(n)%ac72(t)%tdew = AC72_struc(n)%ac72(t)%tdew + td 

  enddo

end subroutine AC72_f2t
