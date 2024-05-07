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
! !ROUTINE: LIS_microMetCorrection
! \label{LIS_microMetCorrection}
!
! !REVISION HISTORY:
!
!  21 Dec 2004: Sujay Kumar; Initial Specification
!  21 Jan 2021: Kristi Arsenault; Update to the code, based on SnowModel
!  22 Mar 2022: Kristi Arsenault; Completed code updates
!
! !INTERFACE:
subroutine LIS_microMetCorrection(n, modelelev, LIS_FORC_Base_State)
! !USES:
  use LIS_logMod
  use ESMF
  use LIS_mpiMod
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_FORC_AttributesMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n   ! Nest index
  real                :: modelelev(LIS_rc%ngrid(n))
  type(ESMF_State)    :: LIS_FORC_Base_State
!
! !DESCRIPTION:
!  Uses the MicroMet topographic downscaling methods to correct
!   temperature, pressure, humidity, shortwave
!   and longwave radiation, and wind field redistribution,
!   based on:
!
!  Ref:
!  Liston, G. E., & Elder, K. (2006). A Meteorological Distribution System for High-Resolution 
!   Terrestrial Modeling (MicroMet), Journal of Hydrometeorology, 7(2), 217-234. 
!   https://journals.ametsoc.org/view/journals/hydr/7/2/jhm486_1.xml
!
! The arguments are:  
! \begin{description}
!  \item[n]
!   index of the n 
!  \item [modelelev]
!   forcing elevation
! \end{description}
!EOP

   integer :: t, k, row, col, i, j, gindex, gid
   real    :: force_tmp, tcforce   ! Input, corrected temperature
   real    :: force_hum, hcforce   ! Input, corrected spec humidity 
   real    :: force_prs, pcforce   ! Input, corrected pressure
   real    :: force_lwd, lcforce   ! Input, corrected LW down rad
   real    :: force_swd, scforce   ! Input, corrected SW down rad
   real    :: force_pcp, pptforce  ! Input, corrected precip
   real    :: force_u,   uwcforce  ! Input, corrected U-wind comp.
   real    :: force_v,   vwcforce  ! Input, corrected U-wind comp.

   real    :: elevdiff, elevdiff2
   real    :: Td, tbar, e, Tdcforce

   integer :: mbefore, mafter
   real    :: xlat, xlon, weight
   real    :: T_lapse_rate, Td_lapse_rate
   real    :: precip_lapse_rate, precip_lapse_rate_m, alfa
   real    :: theta, omegas, thetad, wind, Ww
   real    :: sigma_c, zenith, dec,tau,mu,cosi, psi_dir,psi_dif
   real    :: lhour,czenith
   real    :: ea, epsilona,epsilonb
   real    :: mee, mfe, ee, fe, ratio
   real    :: femiss,emiss
   real    :: elev700, T700, Td700, Rh700, es
   real    :: deg2rad, rad2deg
   real    :: wslope_max, wslopemax_glb, curve_max, curvemax_glb
   real    :: slope, aspect, curvature
   real    :: dirdiff, xmult
   real    :: deltax, deltay, deltaxy, inc, topo

!    B1=17.502 ==> over water, following Buck (1981) reference
!    B2=22.452 ==> over ice
   real, parameter :: A=611.21, B1=17.502, B2=22.452, C=240.97
   real, parameter :: grav = 9.81, rdry = 287., Sstar=1367.0
   real, parameter :: sigma = 5.67E-8
   real, parameter :: Tf = 273.16
   integer, parameter :: bb=2016

   ! Cloud-cover parameters:
   real, parameter :: topo_ref_cloud = 3000.
   real, parameter :: press_ratio = 0.7
   real, parameter :: dx_walcek = 80.0
   real :: f_max, one_minus_RHe, f_1
   real :: delta_topo700, Td_700, Tair_700, rh_700
   real :: e_700, es_700, cloud_frac

   ! Longwave radiation parameters:
   real, parameter :: Stef_Boltz = 5.6696e-8
   real, parameter :: alfa_lw = 1.083
   ! - Constants required for Iziomon et al. (2003).
   real, parameter :: E1 = 200.0, X1 = 0.35, Y1 = 0.100, Z1 = 0.224
   real, parameter :: E2 = 1500.0, X2 = 0.43, Y2 = 0.115, Z2 = 0.320
   real, parameter :: E3 = 3000.0, X3 = 0.51, Y3 = 0.130, Z3 = 1.100
   ! - Assume the X and Y coefficients increase linearly to 3000 m,
   !   and adjust Z to create a best fit to the CLPX data.
   real :: Xs, Ys, Zs, emiss_cloud

   ! Shortwave radiation parameters:
   real, parameter :: solar_const = 1370.
   real, parameter :: days_yr = 365.25
   real, parameter :: Trop_Can = 0.41
   real, parameter :: solstice = 173.
   real :: xhour, xxhour    ! model decimal hour
   real :: sol_dec, hr_angl, cos_Z, sin_Z, cos_i
   real :: sun_azimuth, slope_az_S0
   real :: trans_direct, trans_diffuse
   real :: Qsi_trans_dir, Qsi_trans_dif, Qsi_direct, Qsi_diffuse


   ! Month parameters:
   integer, parameter :: months = 12
   integer lastday(months)
   data lastday/31,28,31,30,31,30,31,31,30,31,30,31/

   ! MicroMet monthly lapse rates (Liston and Kelly, 2006):
   real lapse_rate(months)
   real lapse_rate_nohem(months)
   real lapse_rate_sohem(months)
   data lapse_rate_nohem /0.0044,0.0059,0.0071,0.0078,0.0081,0.0082,&
                          0.0081,0.0081,0.0077,0.0068,0.0055,0.0047/
   data lapse_rate_sohem /0.0081,0.0081,0.0077,0.0068,0.0055,0.0047,&
                          0.0044,0.0059,0.0071,0.0078,0.0081,0.0082/
   ! Atmospheric vapor pressure coeffs are in units of km-1
   !  Following Liston and Kelly (2006)
   real am(months)
   real am_nohem(months)
   real am_sohem(months)
   data am_nohem /0.41,0.42,0.40,0.39,0.38,0.36,&
                  0.33,0.33,0.36,0.37,0.40,0.40/
   data am_sohem /0.33,0.33,0.36,0.37,0.40,0.40,&
                  0.41,0.42,0.40,0.39,0.38,0.36/

   ! The precipitation lapse rate units are in km-1.
   real prec_lapse_rate(months)
   real precip_lapse_rate_nohem(months)
   real precip_lapse_rate_sohem(months)
   data precip_lapse_rate_nohem /0.35,0.35,0.35,0.30,0.25,0.20,&
  &                              0.20,0.20,0.20,0.25,0.30,0.35/
   data precip_lapse_rate_sohem /0.20,0.20,0.20,0.25,0.30,0.35,&
  &                              0.35,0.35,0.35,0.30,0.25,0.20/

   ! Avoid problems of zero (low) winds (for example, turbulence
   !   theory, log wind profile, etc., says that we must have some
   !   wind.  Thus, some equations blow up when the wind speed gets
   !   very small).  This number defines the value that any wind speed
   !   below this gets set to. From snowmodel.par file ...
   real, parameter :: windspd_min = 0.1

   ! The curvature and wind_slope values range between -0.5 and +0.5.
   !   Valid slopewt and curvewt values are between 0 and 1, with
   !   values of 0.5 giving approximately equal weight to slope and
   !   curvature. Glen Liston suggests that slopewt and curvewt be set such
   !   that slopewt + curvewt = 1.0.  This will limit the total
   !   wind weight to between 0.5 and 1.5 (but this is not required).
   real, parameter :: slopewt = 0.58
   real, parameter :: curvewt = 0.42
   integer         :: loops_windwt_smoother

   ! The curvature is used as part of the wind model.  Define a length
   !   scale that the curvature calculation will be performed on.  This
   !   has units of meters, and should be approximately one-half the
   !   wavelength of the topographic features within the domain.
   real, parameter :: curve_len_scale = 300.0

   ! 2D wind fields to apply smoother:
   real, allocatable :: windspd(:,:), winddir(:,:)
   real, allocatable :: wind_slope(:,:), windwt(:,:)
   real, allocatable :: curv2d(:,:)

   ! ESMF metforcing states:
   integer            :: status, ierr
   type(ESMF_Field)   :: tmpField, q2Field,  lwdField, psurfField
   type(ESMF_Field)   :: swdField, pcpField, uField, vField
   real, pointer      :: tmp(:),q2(:),lwd(:),psurf(:)
   real, pointer      :: swd(:),pcp(:),uwind(:),vwind(:)

! .............................................................

    ! ESMF calls for primary metforcing fields 
    ! Temperature
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Tair%varname(1)),tmpField,&
         rc=status)
    call LIS_verify(status,'Error in ESMF_StateGet for tmp in LIS_microMetCorrection')
 
    call ESMF_FieldGet(tmpField,localDE=0, farrayPtr=tmp,rc=status)
    call LIS_verify(status,'Error in ESMF_FieldGet for tmp in LIS_microMetCorrection')

    ! Surface pressure
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Psurf%varname(1)),psurfField,&
         rc=status)
    call LIS_verify(status,'Error in ESMF_StateGet for psurf in LIS_microMetCorrection')

    call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
    call LIS_verify(status,'Error in ESMF_FieldGet for psurf in LIS_microMetCorrection')

    ! Spec humidity
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Qair%varname(1)),q2Field,&
         rc=status)
    call LIS_verify(status,'Error in ESMF_StateGet for q2air in LIS_microMetCorrection')

    call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
    call LIS_verify(status,'Error in ESMF_FieldGet for q2air in LIS_microMetCorrection')

    ! LW down radiation
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_LWdown%varname(1)),lwdField,&
         rc=status)
    call LIS_verify(status,'Error in ESMF_StateGet for lwdown in LIS_microMetCorrection')

    call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
    call LIS_verify(status,'Error in ESMF_FieldGet for lwdown in LIS_microMetCorrection')

    ! SW down radiation
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_SWdown%varname(1)),swdField,&
         rc=status)
    call LIS_verify(status,'Error in ESMF_StateGet for swdown in LIS_microMetCorrection')

    call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
    call LIS_verify(status,'Error in ESMF_FieldGet for swdown in LIS_microMetCorrection')

    ! Total precip (or rainfall)
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LIS_verify(status,'Error in ESMF_StateGet for pcp in LIS_microMetCorrection')

    call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
    call LIS_verify(status,'Error in ESMF_FieldGet for pcp in LIS_microMetCorrection')

    ! U-wind component:
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Wind_E%varname(1)),uField,&
         rc=status)
    call LIS_verify(status,'Error in ESMF_StateGet for uwind in LIS_microMetCorrection')

    call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
    call LIS_verify(status,'Error in ESMF_FieldGet for uwind in LIS_microMetCorrection')

    ! V-wind component:
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Wind_N%varname(1)),vField,&
         rc=status)
    call LIS_verify(status,'Error in ESMF_StateGet for vwind in LIS_microMetCorrection')

    call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
    call LIS_verify(status,'Error in ESMF_FieldGet for vwind in LIS_microMetCorrection')

! .............................................................

   ! For topographic wind distribution calculation
   deg2rad = LIS_CONST_PI / 180.0
   rad2deg = 180.0 / LIS_CONST_PI

   ! Cloud cover fraction calculations (for LWdown and SWdown):
   ! Walcek coefficients.
   f_max = 78.0 + dx_walcek/15.5
   one_minus_RHe = 0.196 + (0.76-dx_walcek/2834.0) * (1.0 - press_ratio)
   f_1 = f_max * (press_ratio - 0.1) / 0.6 / 100.0

   ! Compute the solar declination angle (radians).
   sol_dec = Trop_Can * cos(2.*LIS_CONST_PI * (real(LIS_rc%doy) - solstice)/days_yr)

   ! Find the month before and after the day in question.
   ! - For monthly lapse rate selection.
   if (LIS_rc%da.le.15) then
      mbefore = LIS_rc%mo - 1
      if (mbefore.eq.0) mbefore = 12
      mafter = LIS_rc%mo
      weight = (real(lastday(mbefore)) - 15. + real(LIS_rc%da)) / &
              &  real(lastday(mbefore))
   else
      mbefore = LIS_rc%mo
      mafter = LIS_rc%mo + 1
      if (mafter.eq.13) mafter = 1
      weight = (real(LIS_rc%da) - 15.) / real(lastday(mbefore))
   endif

! --- 
   ! Loop over model-space tiles:
   do t=1,LIS_rc%ntiles(n)

      gindex = LIS_domain(n)%tile(t)%index
      xlat = LIS_domain(n)%grid(gindex)%lat
      xlon = LIS_domain(n)%grid(gindex)%lon
      slope  = LIS_domain(n)%tile(t)%slope
      aspect = LIS_domain(n)%tile(t)%aspect
      curvature = LIS_domain(n)%tile(t)%curvature   ! Now included in LDT output

      if(tmp(t).gt.0) then
        force_tmp = tmp(t)
        force_prs = psurf(t)
        force_hum = q2(t)
        force_pcp = pcp(t)
        force_lwd = lwd(t)
        force_swd = swd(t)
        force_u   = uwind(t)
        force_v   = vwind(t)
        
        ! Read-in lapse rate correction monthly values:
        ! (from MicroMet code: get_lapse_rates routine)
        ! Air, dewpoint temperature, and precipitation.
        do k=1,months
          if( xlat.lt.0.0 ) then   ! Southern Hemisphere
            lapse_rate(k) = lapse_rate_sohem(k)              ! Temp
            am(k) = am_sohem(k)                              ! Td
            prec_lapse_rate(k) = precip_lapse_rate_sohem(k)  ! Precip
          else                     ! Northern Hemisphere
            lapse_rate(k) = lapse_rate_nohem(k)
            am(k) = am_nohem(k)
            prec_lapse_rate(k) = precip_lapse_rate_nohem(k)
          endif
        enddo

        ! Define the temperature lapse rate (Kelvin/1000m).
        T_lapse_rate = (- (weight * lapse_rate(mafter) + &
                 (1. - weight) * lapse_rate(mbefore)))

        ! Define the dew-point temperature lapse rate (deg C/m).
        Td_lapse_rate = (- ((weight * am(mafter) + &
                 (1. - weight) * am(mbefore)) * C)) / (B1 * 1000.0)  ! MicroMet

        ! Define the precipitation lapse rate (km-1).
        precip_lapse_rate = weight * prec_lapse_rate(mafter) + &
                 (1. - weight) * prec_lapse_rate(mbefore)

        ! Calculate elev diff between hi-res elev and forcing elev
        elevdiff = LIS_domain(n)%tile(t)%elev-&
                       modelelev(gindex)
        ! ---
        ! (1) T2air (tmp | tcforce) :  Perform lapse-rate correction on temperature:
        tcforce = force_tmp+(T_lapse_rate*elevdiff)

        ! ---
        ! (2) Surface Pressure (psurf | pcforce)
        ! Pressure - elevation change (from Cosgrove et. al (2003):
        tbar=(force_tmp+tcforce)/2.
        pcforce=force_prs/(exp((LIS_CONST_G*elevdiff)/(rdry*tbar)))

        ! ---
        ! (3) Relative humidity (q2 | hcforce)
        ! Perform lapse-rate correction on humidity fields:
        ! Compute water vapor pressure from spec humidity and pressure
        e = force_hum*force_prs / (0.622+0.378*force_hum)   

        Td = C*log((e/A))/(B1-log((e/A))) + Tf      ! matches MicroMet equation
        Tdcforce = Td + (Td_lapse_rate*elevdiff)    ! final corrected dewpoint temp  

        ! Convert back to spec hum from Td (MicroMet)
        e = A*exp((B1*(Tdcforce-Tf))/(C+(Tdcforce-Tf)))  
        hcforce = (0.622*e) / (pcforce - (0.378*e))

        ! ---
        ! (4) Precipitation (pcp | pptforce)
        ! Use a precipitation "lapse rate", or adjust the precipitation on
        !  the actual elevation grid. The reason the interpolated station elevations 
        !  are used as the topographic reference surface (instead of something
        !  like sea level), is because precip adjustment factor is a non-linear
        !  function of elevation difference.
        !
        ! The adjustment factor that is used comes from: Thornton, P. E.,
        !   S. W. Running, and M. A. White, 1997: Generating surfaces of
        !   daily meteorological variables over large regions of complex
        !   terrain.  J. Hydrology, 190, 214-251.

        ! Convert the precipitation "lapse rate" (km-1) to m-1.
        precip_lapse_rate_m = precip_lapse_rate / 1000.0

        ! Apply Glen Liston's original precipitation increase with
        !  elevation scheme. Also, we could implement (in the future) 
        !  Ward van Pelt's scheme (see van Pelt et al. 2016).
        ! Don't let the elevation difference be greater than some number
        !  (like 1800 meters gives a factor of 4.4).  If it is too large
        !   you get huge precipitation adjustment, a divide by zero, or
        !   even negative adjustments for high elevations).

        elevdiff2 = min(elevdiff,1800.0)   
        alfa = precip_lapse_rate_m * elevdiff2
        pptforce = force_pcp * (1.0 + alfa)/(1.0 - alfa)
        pptforce = max(0.0,pptforce)   ! Ensure no negative precip values

        ! --- 
        ! Incoming/Downward Radiation Fields:
        ! --- 
        ! Cloud cover fraction calculations (for LWdown and SWdown):
        delta_topo700 = topo_ref_cloud - LIS_domain(n)%tile(t)%elev
        Td_700 = Tdcforce + Td_lapse_rate * delta_topo700
        Tair_700 = tcforce + T_lapse_rate * delta_topo700

        ! Convert each Td to a gridded relative humidity (0-1).
        e_700 = A * exp((B1 * (Td_700 - Tf))/(C + (Td_700 - Tf)))
        es_700 = A * exp((B1 * (Tair_700 - Tf))/(C + (Tair_700 - Tf)))
        rh_700 = e_700/es_700
        rh_700 = min(1.0,rh_700)
        rh_700 = max(0.0,rh_700)

        ! Use this RH at 700 mb to define the cloud fraction (0-1).
        cloud_frac = f_1 * exp((rh_700 - 1.0)/one_minus_RHe)

        ! ---
        ! (5) Longwave radiation (lwd | lcforce) -- MicroMet:
        ! Compute Qli following Iziomon et al. (2003).
        if( LIS_domain(n)%tile(t)%elev.lt.E1 ) then
           Xs = X1
           Ys = Y1
           Zs = Z1
        elseif( LIS_domain(n)%tile(t)%elev.gt.E2) then
           Xs = X3
           Ys = Y3
           Zs = Z3
        else
           Xs = X1 + (LIS_domain(n)%tile(t)%elev - E1) * (X3 - X1)/(E3 - E1)
           Ys = Y1 + (LIS_domain(n)%tile(t)%elev - E1) * (Y3 - Y1)/(E3 - E1)
           Zs = Z1 + (LIS_domain(n)%tile(t)%elev - E1) * (Z3 - Z1)/(E3 - E1)
        endif

        emiss_cloud = alfa_lw * &
            (1.0 - Xs * exp((- Ys) * e/tcforce)) * &
            (1.0 + Zs * cloud_frac**2)
        emiss_cloud = min(1.0,emiss_cloud)

        lcforce = emiss_cloud * Stef_Boltz * tcforce**4
#if 0
        ! From Cogsrove et al. (2003) -- 
        ee = (force_hum*force_prs)/0.622
        fe = (hcforce*pcforce)/0.622
        mee = ee/100.
        mfe = fe/100.
       !----------------------------------------------------------------------
       ! correct for negative vapor pressure at very low temperatures at
       ! high latitudes
       !----------------------------------------------------------------------
        if (mee .le. 0) mee = 1e-08
        if (mfe .le. 0) mfe = 1e-08
        emiss  = 1.08*(1-exp(-mee**(force_tmp/bb)))
        femiss = 1.08*(1-exp(-mfe**(tcforce/bb)))
        ratio  = (femiss*(tcforce**4))/(emiss*(force_tmp**4))
        lcforce= force_lwd*ratio
#endif

        ! ---
        ! (6) Shortwave radiation (swd | scforce) -- MicroMet:
        !  "solar" routine from MicroMet calculates the "local" SWdown, and
        !   only set up for radiation values for <= 3 hour model timestep.

        ! In solar radiation, adjust the time to correspond to the
        !  local time at this longitude position.
        xxhour = LIS_rc%hr + xlon / 15.0
        if (xxhour.ge.24.0) xxhour = xxhour - 24.0
        if (xxhour.lt.0.0)  xxhour = xxhour + 24.0

        ! Compute the sun's hour angle (radians).
        hr_angl = (xxhour * 15.0 - 180.0) * deg2rad

        ! Compute cos_Z.  Note that the sin of the solar elevation angle,
        !   sin_alfa, is equal to the cosine of the solar zenith angle,
        !   cos_Z.
        cos_Z = sin(sol_dec) * sin(xlat * deg2rad) + &
                cos(sol_dec) * cos(xlat * deg2rad) * cos(hr_angl)
        cos_Z = max(0.0,cos_Z)

        ! Account for clouds, water vapor, pollution, etc.
        trans_direct = (0.6 + 0.2 * cos_Z) * (1.0 - cloud_frac)
        trans_diffuse = (0.3 + 0.1 * cos_Z) * cloud_frac

        ! Compute the solar radiation transmitted through the atmosphere.
        Qsi_trans_dir = solar_const * trans_direct
        Qsi_trans_dif = solar_const * trans_diffuse

        ! COMPUTE THE CORRECTIONS TO ALLOW FOR TOPOGRAPHIC SLOPE AND ASPECT.
        ! The sine of the solar zenith angle.
        sin_Z = sqrt(1.0 - cos_Z*cos_Z)

        ! Azimuth of the sun, with south having zero azimuth for the
        !   northern hemisphere.
        sun_azimuth = asin(max(-1.0,min(1.0,cos(sol_dec)*sin(hr_angl)/sin_Z)))
        if (xlat.lt.0.0) then
          sun_azimuth = - sun_azimuth
        endif

        ! Make the corrections so that the angles below the local horizon
        !   are still measured from the normal to the slope.
        if (xlat.ge.0.0) then
          if (hr_angl.lt.0.0) then
            if (hr_angl.lt.sun_azimuth) sun_azimuth = - LIS_CONST_PI - sun_azimuth
          elseif (hr_angl.gt.0.0) then
            if (hr_angl.gt.sun_azimuth) sun_azimuth = LIS_CONST_PI - sun_azimuth
          endif
        endif

        ! Build, from the variable with north having zero azimuth, a 
        !   slope_azimuth value with south having zero azimuth.  Also
        !   make north have zero azimuth if in the southern hemsisphere.
        if (xlat.ge.0.0) then
          if( aspect.ge.180.0 ) then
             slope_az_S0 = aspect - 180.0
          else
             slope_az_S0 = aspect + 180.0
          endif
        else
           slope_az_S0 = aspect
        endif

        ! Compute the angle between the normal to the slope and the angle
        !   at which the direct solar radiation impinges on the sloping
        !   terrain (radians).
        cos_i = cos(slope * deg2rad) * cos_Z + &
                sin(slope * deg2rad) * sin_Z * &
                cos(sun_azimuth - slope_az_S0 * deg2rad)

        ! Adjust the topographic correction due to local slope so that
        !   the correction is zero if the sun is below the local horizon 
        !   (i.e., the slope is in the shade) or if the sun is below the
        !   global horizon.
        if (cos_i.lt.0.0) cos_i = 0.0
        if (cos_Z.le.0.0) cos_i = 0.0

        ! Adjust the solar radiation for slope, etc.
        Qsi_direct = cos_i * Qsi_trans_dir
        Qsi_diffuse = cos_Z * Qsi_trans_dif

        ! Combine the direct and diffuse solar components.
        scforce = Qsi_direct + Qsi_diffuse

        ! ------------------------------------
        ! Convert back to LIS-model tile space:
        tmp(t)   = tcforce
        psurf(t) = pcforce
        q2(t)    = hcforce
        pcp(t)   = pptforce
        lwd(t)   = lcforce
        swd(t)   = scforce

      endif
   enddo

   ! ---
   ! (7/8) Wind fields (uwind | uwcforce; vwind | vwcforce)

   ! Convert domain dx, dy into meters:
   if( LIS_domain(n)%dx <= 1.) then  ! Temp. solution for lat/lon dx,dy
     deltax = LIS_domain(n)%dx*100000.
     deltay = LIS_domain(n)%dy*100000.
   else
     deltax = LIS_domain(n)%dx*1000.
     deltay = LIS_domain(n)%dy*1000.
   endif

   ! Smooth the wind weighting factor to eliminate any sharp speed
   !   increases resulting from the merging of the curve wt and the
   !   slope wt.  Define the number of times this is done to be a
   !   function of the curvature length scale and the grid increment.
   !   The first 0.5 here just means that half of the caclulated
   !   number appears to be about right (additional loops with this
   !   smoother does not change the results much).  If there are
   !   unwanted wave features in the snow distribution, this 0.5
   !   factor can be increased to 1.0 or more, to get rid of these
   !   waves.  Also see "loops_snowd_smoother" in snowtran_code.f90.
   ! xmult = 0.5  |  xmult = 1.5
   xmult = 1.0     ! Based on MicroMet code

   loops_windwt_smoother = nint(xmult * curve_len_scale / &
                               (0.5 * (deltax + deltay)))

   ! Set up the 2D wind fields needed:
   allocate( windspd(LIS_rc%lnc(n),LIS_rc%lnr(n)),winddir(LIS_rc%lnc(n),LIS_rc%lnr(n)) )
   allocate( wind_slope(LIS_rc%lnc(n),LIS_rc%lnr(n)),windwt(LIS_rc%lnc(n),LIS_rc%lnr(n)) )
   allocate( curv2d(LIS_rc%lnc(n),LIS_rc%lnr(n)) )

   windspd = LIS_rc%udef
   winddir = LIS_rc%udef
   wind_slope = LIS_rc%udef
   windwt = LIS_rc%udef
   curv2d = LIS_rc%udef

   ! Loop over model-space tiles and convert key wind fields to 2D:
   do t=1,LIS_rc%ntiles(n)
      col = LIS_domain(n)%tile(t)%col
      row = LIS_domain(n)%tile(t)%row

      slope  = LIS_domain(n)%tile(t)%slope
      aspect = LIS_domain(n)%tile(t)%aspect
      curvature = LIS_domain(n)%tile(t)%curvature   ! Now included in LDT output
      curv2d(col,row) = LIS_domain(n)%tile(t)%curvature

      ! Convert these u and v components to speed and directions.
      if(tmp(t).gt.0) then
        force_u = uwind(t)
        force_v = vwind(t)

        ! Some compilers do not allow both u and v to be 0.0 in the atan2 computation.
        if( abs(force_u).lt.1e-10 ) force_u = 1e-10
        winddir(col,row) = rad2deg * atan2(force_u,force_v)
        if( winddir(col,row).ge.180.0 ) then
          winddir(col,row) = winddir(col,row) - 180.0
        else
          winddir(col,row) = winddir(col,row) + 180.0
        endif
        windspd(col,row) = sqrt(force_u**2 + force_v**2)

        ! Modify the wind speed and direction according to simple wind-topography
        ! From "topo_mod_winds" subroutine in micromet_code.f90.
        !
        ! Compute the wind modification factor which is a function of
        !  topography and wind direction following Liston and Sturm (1998).

        ! Compute the slope in the direction of the wind.
        wind_slope(col,row) = deg2rad * slope * &
             cos(deg2rad * (winddir(col,row) - aspect))

      endif
   enddo

   ! Scale the wind slope such that the max abs(wind slope) has a value
   !   of abs(0.5).  Include a 1 mm slope in slope_max to prevent
   !   divisions by zero in flat terrain where the slope is zero.
   wslope_max = 0.0 + 0.001
   do j=1,LIS_rc%lnr(n)
     do i=1,LIS_rc%lnc(n)
        wslope_max = max(wslope_max,abs(wind_slope(i,j)))
     enddo
   enddo

! Search for global "max" of wslope_max value for all parallel subdomains:
#if (defined SPMD)
        call MPI_Barrier(LIS_MPI_COMM, ierr)
        call MPI_ALLREDUCE(wslope_max, wslopemax_glb, 1,&
                 MPI_REAL, MPI_MAX, LIS_mpi_comm, ierr)
        wslope_max = wslopemax_glb
#endif
! -----

   ! Apply global slope max to 2D wind_slope calculation:
   do j=1,LIS_rc%lnr(n)
     do i=1,LIS_rc%lnc(n)

        wind_slope(i,j) = wind_slope(i,j) / (2.0 * wslope_max)

        ! Calculate the wind speed and direction adjustments. The
        !  curvature and wind_slope values range between -0.5 and +0.5.
        !  Valid slopewt and curvewt values are between 0 and 1, with
        !  values of 0.5 giving approximately equal weight to slope and
        !  curvature. Glen Liston suggests that slopewt and curvewt be set such
        !  that slopewt + curvewt = 1.0.  This will limit the total
        !  wind weight to between 0.5 and 1.5 (but this is not required).

        ! Compute the wind weighting factor -- dependent on wind direction input.
        windwt(i,j) = 1.0 + slopewt * wind_slope(i,j) + curvewt * curv2d(i,j)

     enddo
   enddo

   ! Don't apply smoothing if the domain is arbitrarily small (<=100 gnc, gnr).
   if (LIS_rc%gnc(n).gt.100 .and. LIS_rc%gnr(n).gt.100) then
     do k=1,loops_windwt_smoother
        call wind_smoother(LIS_rc%lnc(n),LIS_rc%lnr(n),windwt)  !<-- windwt(nc,nr)
     enddo
   endif

   ! Complete final calculations and translation back to 1D tilespace:
   do t=1,LIS_rc%ntiles(n)
      i = LIS_domain(n)%tile(t)%col
      j = LIS_domain(n)%tile(t)%row
      aspect = LIS_domain(n)%tile(t)%aspect

      if(tmp(t).gt.0) then

        ! Generate the terrain-modified wind speed.
        windspd(i,j) = windwt(i,j) * windspd(i,j)

        ! Modify the wind direction according to Ryan (1977).  Note that it
        !   is critical that "dirdiff" handles the cases where the slope
        !   azimuth and the wind direction are on different sides of the
        !   360-0 line.
        if( aspect.gt.270.0 .and. winddir(i,j).lt.90.0 ) then
          dirdiff = aspect - winddir(i,j) - 360.0
        elseif( aspect.lt.90.0 .and. winddir(i,j).gt.270.0 ) then
          dirdiff = aspect - winddir(i,j) + 360.0
        else
          dirdiff = aspect - winddir(i,j)
        endif
        if( abs(dirdiff).le.90.0 ) then
          winddir(i,j) = winddir(i,j) - 0.5 * min(wind_slope(i,j)*rad2deg,45.0) * &
                    sin(deg2rad * (2.0 * dirdiff))
          if( winddir(i,j).gt.360.0 ) then
             winddir(i,j) = winddir(i,j) - 360.0
          elseif( winddir(i,j).lt.0.0 ) then
             winddir(i,j) = winddir(i,j) + 360.0
          endif
        endif

        ! Avoid problems of zero (low) winds (for example, turbulence theory,
        !  log wind profile, etc., says that we must have some wind. Thus,
        !   some equations blow up when the wind speed gets very small.
        if( windspd(i,j).lt.windspd_min ) then
          windspd(i,j) = windspd_min
        endif

        ! Final conversion back to corrected uwind and vwind components:
        uwcforce = (- windspd(i,j)) * sin(deg2rad*winddir(i,j))
        vwcforce = (- windspd(i,j)) * cos(deg2rad*winddir(i,j))

        uwind(t) = uwcforce
        vwind(t) = vwcforce

      endif
   enddo

   deallocate( windspd, winddir )
   deallocate( wind_slope, windwt, curv2d )

   write(LIS_logunit,*) '[INFO] MicroMet forcing topographic downscaling completed'

end subroutine LIS_microMetCorrection

 
subroutine wind_smoother (nx, ny, weight)

   implicit none

! !ARGUMENTS: 
  integer, intent(in) :: nx, ny
  real, intent(inout) :: weight(nx,ny)
!
! !DESCRIPTION:
!
! Routine performs a 9-point smoothing operation.
!
! The result at each grid point is a weighted average of the grid
!  point and the surrounding 8 points.  The center point receives
!  a weight of 1.0, the points at each side and above and below
!  receive a weight of 0.5, and corner points receive a weight of
!  0.3.  All points are multiplied by their weights and summed,
!  then divided by the total weight.
!
!  Ref:
!  Liston, G. E., & Elder, K. (2006). A Meteorological Distribution System for High-Resolution 
!   Terrestrial Modeling (MicroMet), Journal of Hydrometeorology, 7(2), 217-234. 
!   https://journals.ametsoc.org/view/journals/hydr/7/2/jhm486_1.xml

   integer :: i,j
   real    ::  weight_tmp(nx,ny)

      ! Do the interior.
      do i=2,nx-1
        do j=2,ny-1
          weight_tmp(i,j) = (weight(i,j) + 0.5 * (weight(i,j-1) + &
            weight(i,j+1) + weight(i-1,j) + weight(i+1,j)) + 0.3 *&
            (weight(i-1,j-1) + weight(i+1,j+1) + weight(i-1,j+1) +&
            weight(i+1,j-1))) / 4.2
        enddo
      enddo

      ! Do the sides.
      j = 1
      do i=2,nx-1
        weight_tmp(i,j) = (weight(i,j) + 0.5 * (weight(i,j+1) + weight(i-1,j) +&
          weight(i+1,j)) + 0.3 * (weight(i+1,j+1) + weight(i-1,j+1))) / 3.1
      enddo

      j = ny
      do i=2,nx-1
        weight_tmp(i,j) = (weight(i,j) + 0.5 * (weight(i,j-1) + weight(i-1,j) +&
          weight(i+1,j)) + 0.3 * (weight(i+1,j-1) + weight(i-1,j-1))) / 3.1
      enddo

      i = 1
      do j=2,ny-1
        weight_tmp(i,j) = (weight(i,j) + 0.5 * (weight(i,j-1) + weight(i,j+1) +&
          weight(i+1,j)) + 0.3 * (weight(i+1,j-1) + weight(i+1,j+1))) / 3.1
      enddo

      i = nx
      do j=2,ny-1
        weight_tmp(i,j) = (weight(i,j) + 0.5 * (weight(i,j-1) + weight(i,j+1) +&
          weight(i-1,j)) + 0.3 * (weight(i-1,j-1) + weight(i-1,j+1))) / 3.1
      enddo

      ! Do the corners.
      i = 1
      j = 1
      weight_tmp(i,j) = (weight(i,j) + 0.5 * (weight(i,j+1) + weight(i+1,j)) +&
        0.3 * weight(i+1,j+1)) / 2.3

      i = nx
      j = 1
      weight_tmp(i,j) = (weight(i,j) + 0.5 * (weight(i,j+1) + weight(i-1,j)) +&
        0.3 * weight(i-1,j+1)) / 2.3

      i = 1
      j = ny
      weight_tmp(i,j) = (weight(i,j) + 0.5 * (weight(i,j-1) + weight(i+1,j)) +&
        0.3 * weight(i+1,j-1)) / 2.3

      i = nx
      j = ny
      weight_tmp(i,j) = (weight(i,j) + 0.5 * (weight(i,j-1) + weight(i-1,j)) +&
        0.3 * weight(i-1,j-1)) / 2.3

      ! Return the smoothed array.
      do i=1,nx
        do j=1,ny
          weight(i,j) = weight_tmp(i,j)
        enddo
      enddo

end subroutine 
