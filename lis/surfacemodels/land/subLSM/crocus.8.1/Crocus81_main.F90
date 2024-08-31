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
!
! !ROUTINE: Crocus81_main
! \label{Crocus81_main}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   10/18/19: Mahdi Navari, Shugong Wang; initial implementation for Crocus81 with LIS-7
!   9 Dec 2020: Mahdi Navari; edited to take into account the Crocus slope correction
!   19 Jan 2021: Mahdi Navari, edited to properly initialize precipitation 
!   21 Jan 2021: Mahdi Navari, edited to properly assign values for TG, XWG, and XWGI for the stand-alone version
!
! !INTERFACE:
subroutine Crocus81_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod 
    use Crocus81_lsmMod
   !use other modules
  
    implicit none
! !ARGUMENTS:
    integer, intent(in)  :: n
    integer              :: t
    integer              :: i
    real                 :: dt
    real                 :: lat, lon
    integer              :: row, col
    integer              :: year, month, day, hour, minute, second
    logical              :: alarmCheck

!
! !DESCRIPTION:
!  This is the entry point for calling the Crocus81 physics.
!  This routine calls the {\tt crocus_driver } routine that performs the
!  land surface computations, to solve for water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for Crocus81
    integer              :: tmp_n                  ! nest id [-]
    integer              :: tmp_nsnow              ! number of snow layer []
    integer              :: tmp_nimpur             ! number of impurtites []
    integer              :: tmp_year               ! year of the currrent time step [-]
    integer              :: tmp_month              ! month of the current time step [-]
    integer              :: tmp_day                ! day of the current time step [-]
    integer              :: tmp_hour               ! hour of the current time step [-]
    integer              :: tmp_minute             ! minute of the current time step [-]
    CHARACTER(len=3)     :: tmp_SNOWRES_opt        ! SNOWRES_opt  = ISBA-SNOW3L turbulant exchange option
!   'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!   'RIL' = Limit Richarson number under very stable
! conditions (currently testing)
!   'M98'  = Martin et Lejeune 1998 : older computation for turbulent fluxes coefficents in Crocus []
    LOGICAL              :: tmp_OMEB_BOOL          ! !True = coupled to MEB. surface fluxes are IMPOSED
! as an upper boundary condition to the explicit snow schemes. 
!If = False, then energy budget and fluxes are computed herein. [-]
    LOGICAL              :: tmp_GLACIER_BOOL       ! !True = Over permanent snow and ice, initialise WGI=WSAT, Hsnow>=10m and allow 0.8<SNOALB<0.85
! False = No specific treatment [-]
    CHARACTER(len=3)     :: tmp_HIMPLICIT_WIND_opt ! wind implicitation option  'OLD' = direct , 'NEW' = Taylor serie, order 1 [-]
    REAL*8, allocatable    :: tmp_SNOWSWE(:)         ! Snow layer(s) liquid Water Equivalent (SWE:kg m-2) [kg/m2]
    REAL*8, allocatable    :: tmp_SNOWRHO(:)         ! Snow layer(s) averaged density (kg/m3) [kg/m3]
    REAL*8, allocatable    :: tmp_SNOWHEAT(:)        ! Snow layer(s) heat content (J/m2) [J/m2]
    REAL*8                 :: tmp_SNOWALB            ! snow surface albedo [-]
    REAL*8, allocatable    :: tmp_SNOWGRAN1(:)       ! Snow layers grain feature 1 [-]
    REAL*8, allocatable    :: tmp_SNOWGRAN2(:)       ! Snow layer grain feature 2 [-]
    REAL*8, allocatable    :: tmp_SNOWHIST(:)        ! Snow layer grain historical parameter (only for non dendritic snow) (-) in {0-5} [-]
    REAL*8, allocatable    :: tmp_SNOWAGE(:)         ! Age since snowfall (day) [day]
    REAL*8                 :: tmp_PTSTEP             ! time step of the integration [s]
    REAL*8                 :: tmp_PPS                ! pressure at atmospheric model surface (Pa) [Pa]
    REAL*8                 :: tmp_SRSNOW             ! snow rate (SWE) [kg/(m2 s)] [kg/(m2 s)]
    REAL*8                 :: tmp_RRSNOW             ! rain rate [kg/(m2 s)] [kg/(m2 s)]
    REAL*8                 :: tmp_TA                 ! atmospheric temperature at level za (K) [K]
    REAL*8                 :: tmp_TG !(:)              ! Surface soil temperature (effective temperature the of layer lying below snow) (K)  (for snowcro.F90 we only use the surface layer ZP_TG(:,1))  (#nsoil depends on 2-L, 3-L DIF) [K]
    REAL*8                 :: tmp_SW_RAD             ! incoming solar radiation (W/m2) [W/m2]
    REAL*8                 :: tmp_QA                 ! atmospheric specific humidity at level za [-]
    REAL*8                 :: tmp_Wind_E             ! Eastward Wind [m/s]
    REAL*8                 :: tmp_Wind_N             ! Northward Wind [m/s]
    REAL*8                 :: tmp_LW_RAD             ! atmospheric infrared radiation (W/m2) [W/m2]
    REAL                 :: tmp_UREF               ! reference height of the wind [m]
    REAL                 :: tmp_SLOPE              ! angle between the normal to the surface and the vertical  (MN:  replaced PDIRCOSZW with slope and computed the cosine in the driver) [-]
    REAL                 :: tmp_ZREF               ! Reference height of the first atmospheric level (m) [m]
    REAL                 :: tmp_Z0NAT              ! grid box average roughness length (m) (roughness length for momentum) [m]
    REAL                 :: tmp_Z0EFF              ! roughness length for momentum (modd_diagn.F90 effective roughness length for heat(!?)) [m]
    REAL                 :: tmp_Z0HNAT             ! grid box average roughness length for heat [m]
    REAL*8, allocatable  :: tmp_ALB (:)             ! soil/vegetation albedo[-] (monthly value)
    LOGICAL              :: tmp_use_monthly_albedo_map ! if usemonalb == .true., then the alb value passed to lsmcrocus will be used 
                                                       ! as the background snow-free albedo term.  ! if usemonalb == .false., then alb will be set to 0.2
    REAL*8                 :: tmp_SOILCOND           ! soil thermal conductivity (W m-1 K-1) [W /(m K)]
    REAL                 :: tmp_D_G                ! !Assumed first soil layer thickness (m)
!Used to calculate ground/snow heat flux   (D_G(:,1)) [m]
    REAL*8, allocatable    :: tmp_SNOWLIQ(:)         ! Snow layer(s) liquid water content (m) [m]
    REAL*8, allocatable    :: tmp_SNOWTEMP(:)        ! Snow layer(s) temperature (K) [K]
    REAL*8, allocatable    :: tmp_SNOWDZ(:)          ! Snow layer(s) thickness (m) [m]
    REAL*8                 :: tmp_THRUFAL            ! Rate that liquid water leaves snow pack: paritioned into soil infiltration/runoff  by ISBA [kg/(m2 s)] [kg/(m2 s)]
    REAL*8                 :: tmp_GRNDFLUX           ! Soil/snow interface heat flux (W/m2) [W/m2]
    REAL                 :: tmp_SNDRIFT            ! Blowing snow sublimation (kg/m2/s) NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592) [kg/m2/s]
    REAL                 :: tmp_RI_n               ! Richardson number (-)  NOTE: RI has not been initialized in CALL_MODEL (If not OMED initalized to undefined in the snow3L_isba.F90) [-]
    REAL*8                 :: tmp_EMISNOW            ! Snow surface emissivity (initialize to XEMISSN
!XEMISSN  comes from MODD_SNOW_PAR) [-]
    REAL                 :: tmp_CDSNOW             ! Drag coefficient for momentum over snow (-) [-]
    REAL                 :: tmp_USTARSNOW          ! Friction velocity over snow (m/s); [m/s]
    REAL                 :: tmp_CHSNOW             ! Drag coefficient for heat over snow  (-) [-]
    REAL*8                 :: tmp_SNOWHMASS          ! Heat content change due to mass changes in snowpack: for budget calculations only. [J/m2]
    REAL*8                 :: tmp_QS                 ! Surface humidity (kg/kg) [kg/kg]
    REAL                 :: tmp_PERMSNOWFRAC       ! Fraction of permanet snow/ice [-]
    real                 :: tmp_LAT                ! Latitude in decimal degree  (latitude (degrees +North)) [degrees]
    real                 :: tmp_LON                ! Longitude in decimal year    (longitude (degrees +East)) [degrees]
    CHARACTER(len=4)     :: tmp_SNOWDRIFT_opt      ! Mechanical transformation of snow grain and compaction + effect of wind on falling snow properties
!	'NONE': No snowdrift scheme
!	'DFLT': falling snow falls as purely dendritic
! 	'GA01': Gallee et al 2001
!	'VI13': Vionnet et al 2013 [-]
    LOGICAL              :: tmp_SNOWDRIFT_SUBLIM_BOOL ! Logicals for snowdrift sublimation [-]
    LOGICAL              :: tmp_SNOW_ABS_ZENITH_BOOL ! Activate parametrization of solar absorption for polar regions (If True modify solar absorption as a function of solar zenithal angle) [-]
    CHARACTER(len=3)     :: tmp_SNOWMETAMO_opt     ! Metamorphism scheme: B92 (historical version, Brun et al 92), C13, T07, F06 (see Carmagnola et al 2014) [-]
    CHARACTER(len=3)     :: tmp_SNOWRAD_opt        ! Radiative transfer scheme. HSNOWRAD=B92 Brun et al 1992.  HSNOWRAD=T17 (Tuzet et al. 2017) (Libois et al. 2013) TARTES with impurities content scheme [-]
    LOGICAL              :: tmp_ATMORAD_BOOL       ! Activate atmotartes scheme  (default=.FALSE. # This option is not stable yet, but it is supposed to compute the direct/diffuse ratio directly from atmospheric informations (AOD, Ozone column, Water column..)) [-]
    REAL, allocatable    :: tmp_IMPWET(:)          ! [init_surf_atmn.F90   --> wet deposit coefficient for each impurity type (g/m_/s)  ,   snowcro.F90 --> Dry and wet deposit coefficient from Forcing File(g/m_/s)   (# 1553 You can either feed the model with prescribed and constant deposition fluxes or introduce a wet and dry deposition field directly in the forcing file. ) [g/m_/s]
    REAL, allocatable    :: tmp_IMPDRY(:)          ! [init_surf_atmn.F90   --> wet deposit coefficient for each impurity type (g)  ,   snowcro.F90 --> Dry and wet deposit coefficient from Forcing File(g/m_/s) [g/m_/s]
    CHARACTER(len=3)     :: tmp_SNOWFALL_opt       ! New options for multiphysics version (Cluzet et al 2016). Falling snow scheme: V12 (Vionnet et al. 2012) , A76 (Anderson 1976), S02 (Lehning and al. 2002), P75 (Pahaut 1975) [-]
    CHARACTER(len=3)     :: tmp_SNOWCOND_opt       ! Thermal conductivity scheme: Y81 (Yen 1981), I02 (Boone et al. 2002) C11 (Calonne et al. 2011) [-]
    CHARACTER(len=3)     :: tmp_SNOWHOLD_opt       ! liquid water content scheme: B92 (Brun et al. 1992) O04 (Oleson et al., 2004) S02 (SNOWPACK, Lehning et al, 2002) B02 (ISBA_ES, Boone et al. 2002) [-]
    CHARACTER(len=3)     :: tmp_SNOWCOMP_opt       ! B92 snow compaction basis version and B93 for slightly different parameters   (NOTE: IN THE CODE S14, T11, ) [-]
    CHARACTER(len=3)     :: tmp_SNOWZREF_opt       ! reference height is constant or variable from the snow surface: CST (constant from snow surface, i.e. Col de Porte) or VAR (variable from snow surface = snow depth has to be removed from reference height) [-]
    REAL*8                 :: tmp_SNOWMAK_dz         ! Snowmaking thickness (m) [m]
    LOGICAL              :: tmp_SNOWCOMPACT_BOOL   ! Snowmaking and Grooming options [-]
    LOGICAL              :: tmp_SNOWMAK_BOOL       ! Snowmaking and Grooming options [-]
    LOGICAL              :: tmp_SNOWTILLER_BOOL    ! Snowmaking and Grooming options [-]
    LOGICAL              :: tmp_SELF_PROD_BOOL     ! Snowmaking and Grooming options [-]
    LOGICAL              :: tmp_SNOWMAK_PROP_BOOL  ! Snowmaking and Grooming options [-]
    LOGICAL              :: tmp_PRODSNOWMAK_BOOL   ! Snowmaking and Grooming options [-]
    REAL                 :: tmp_SLOPE_DIR          ! !Typical slope aspect in the grid  (clockwise from N) [Radians]
    REAL*8                 :: tmp_SAND               ! Soil sand fraction (-) [-]
    REAL*8                 :: tmp_SILT               ! Soil silt fraction (-) [-]
    REAL*8                 :: tmp_CLAY               ! Soil clay fraction (-) [-]
    REAL*8                 :: tmp_POROSITY           ! Soil porosity (m3 m-3) [m3/m3]
    REAL*8                 :: tmp_XWGI               ! soil volumetric frozen water content
    REAL*8                 :: tmp_XWG                ! soil volumetric liquid water content
    LOGICAL              :: tmp_Partition_total_precip_BOOL ! Boolean option to partition total precipitation into snowfall and rainfall using Jordan 1991 [-]
    REAL                 :: tmp_SD_1D              ! Total snow depth, temporally added [m]  
    REAL                 :: tmp_SWE_1D             ! Total SWE, temporally added [kg/m2] 
    REAL                 :: FPICE
    character*3        :: fnest ! MN Bug in toolkit (added to this code)
    real*8    :: tmp
    
    allocate( tmp_SNOWSWE( CROCUS81_struc(n)%nsnow ) )
    allocate( tmp_SNOWRHO( CROCUS81_struc(n)%nsnow ) )
    allocate( tmp_SNOWHEAT( CROCUS81_struc(n)%nsnow ) )
    allocate( tmp_SNOWGRAN1( CROCUS81_struc(n)%nsnow ) )
    allocate( tmp_SNOWGRAN2( CROCUS81_struc(n)%nsnow ) )
    allocate( tmp_SNOWHIST( CROCUS81_struc(n)%nsnow ) )
    allocate( tmp_SNOWAGE( CROCUS81_struc(n)%nsnow ) )
!    allocate( tmp_TG( 12 ) )
    allocate( tmp_ALB( 12 ) )
    allocate( tmp_SNOWLIQ( CROCUS81_struc(n)%nsnow ) )
    allocate( tmp_SNOWTEMP( CROCUS81_struc(n)%nsnow ) )
    allocate( tmp_SNOWDZ( CROCUS81_struc(n)%nsnow ) )
    allocate( tmp_IMPWET( CROCUS81_struc(n)%nimpur ) )
    allocate( tmp_IMPDRY( CROCUS81_struc(n)%nimpur ) )

    ! check Crocus81 alarm. If alarm is ring, run model. 
     write(fnest,'(i3.3)') n  !MN  Bug in the toolkit 
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "CROCUS81 model alarm "// trim(fnest)) !MN  Bug in the toolkit 
    if (alarmCheck) Then
        do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)

            dt = LIS_rc%ts
            row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
            lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
           
            ! retrieve forcing data from CROCUS81_struc(n)%crocus81(t) and assign to local variables
            ! PPS: pressure at atmospheric model surface (Pa)
            tmp_PPS        = CROCUS81_struc(n)%crocus81(t)%PPS    / CROCUS81_struc(n)%forc_count
 
            ! SRSNOW: snow rate (SWE) [kg/(m2 s)]
            if (CROCUS81_struc(n)%Partition_total_precip_BOOL) then
               tmp_SRSNOW     = 0.0
            else 
               tmp_SRSNOW     = CROCUS81_struc(n)%crocus81(t)%SRSNOW / CROCUS81_struc(n)%forc_count
            endif
            ! RRSNOW: rain rate [kg/(m2 s)]
            tmp_RRSNOW     = CROCUS81_struc(n)%crocus81(t)%RRSNOW / CROCUS81_struc(n)%forc_count

            ! TA: atmospheric temperature at level za (K)
            tmp_TA         = CROCUS81_struc(n)%crocus81(t)%TA     / CROCUS81_struc(n)%forc_count
 
            ! SW_RAD: incoming solar radiation (W/m2)
            tmp_SW_RAD     = CROCUS81_struc(n)%crocus81(t)%SW_RAD / CROCUS81_struc(n)%forc_count
 
            ! QA: atmospheric specific humidity at level za
            tmp_QA         = CROCUS81_struc(n)%crocus81(t)%QA     / CROCUS81_struc(n)%forc_count
 
            ! Wind_E: Eastward Wind
            tmp_Wind_E     = CROCUS81_struc(n)%crocus81(t)%Wind_E / CROCUS81_struc(n)%forc_count
 
            ! Wind_N: Northward Wind
            tmp_Wind_N     = CROCUS81_struc(n)%crocus81(t)%Wind_N / CROCUS81_struc(n)%forc_count
 
            ! LW_RAD: atmospheric infrared radiation (W/m2)
            tmp_LW_RAD     = CROCUS81_struc(n)%crocus81(t)%LW_RAD / CROCUS81_struc(n)%forc_count

            ! Slope correction 
            ! As a convention, only vertical incoming and outgoing fluxes are 
            ! provided to and from the model; the correction of these terms
            ! according to the local slope is carried out within SURFEX.
            ! Similarly, variables such as total snow depth, total snow water
            ! equivalent, and the corresponding variables for each layer
            ! are output by the model in terms of their vertical component,
            ! i.e. projected vertically. (Vionnet et al 2012)
            
            ! For simplicity, we will use the cosine of the slope to edit the fluxes. 
            ! In the SIRFEX code, they have used sky view factor (1+cos(slope))/2 
            ! to correct the component of SW and LW radiation.
            ! See coupling_isba_orographyn.F90 in the SURFEX code
            tmp_SLOPE = CROCUS81_struc(n)%crocus81(t)%SLOPE
            if (tmp_SRSNOW .NE. LIS_rc%udef) then
               tmp_SRSNOW = tmp_SRSNOW * COS(tmp_SLOPE)
            endif
            if (tmp_RRSNOW .NE. LIS_rc%udef) then
                tmp_RRSNOW = tmp_RRSNOW * COS(tmp_SLOPE)
            endif
            if (tmp_SW_RAD .NE. LIS_rc%udef) then
               tmp_SW_RAD = tmp_SW_RAD * COS(tmp_SLOPE)
            endif
            if (tmp_LW_RAD .NE. LIS_rc%udef) then
               tmp_LW_RAD = tmp_LW_RAD * COS(tmp_SLOPE)
            endif
            ! Added lis.config option to determine whether to partition total 
            ! precipitation into snowfall and rainfall 
            ! Use Jordan model to partition total precipitation into snowfall and rainfall 
            ! Jordan (1991)
            if (CROCUS81_struc(n)%Partition_total_precip_BOOL) then 
              IF(tmp_TA > 273.16+2.5)THEN
                  FPICE = 0.
              ELSE
                IF(tmp_TA <= 273.16+0.5)THEN
                  FPICE = 1.0
                ELSE IF(tmp_TA <= 273.16+2.)THEN
                  FPICE = 1.-(-54.632 + 0.2*tmp_TA)
                ELSE
                  FPICE = 0.6
                ENDIF
              ENDIF
              ! snow rate (SWE) [kg/(m2 s)]
              tmp_SRSNOW  = tmp_RRSNOW * FPICE 

              ! rain rate (SWE) [kg/(m2 s)]
              tmp_RRSNOW  = tmp_RRSNOW - tmp_SRSNOW
            endif
            !
            ! check validity of PPS
            if(tmp_PPS .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable PPS in Crocus81"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of SRSNOW
            if(tmp_SRSNOW .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable SRSNOW in Crocus81"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of RRSNOW
            if(tmp_RRSNOW .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable RRSNOW in Crocus81"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of TA
            if(tmp_TA .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable TA in Crocus81"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of SW_RAD
            if(tmp_SW_RAD .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable SW_RAD in Crocus81"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of QA
            if(tmp_QA .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable QA in Crocus81"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Wind_E
            if(tmp_Wind_E .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable Wind_E in Crocus81"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Wind_N
            if(tmp_Wind_N .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable Wind_N in Crocus81"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of LW_RAD
            if(tmp_LW_RAD .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable LW_RAD in Crocus81"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! 
            tmp_LAT  = lat
            tmp_LON = lon
            tmp_year   = LIS_rc%yr
            tmp_month  = LIS_rc%mo
            tmp_day    = LIS_rc%da
            tmp_hour   = LIS_rc%hr
            tmp_minute = LIS_rc%mn
 
            ! get parameters 
            tmp_nsnow                               = CROCUS81_struc(n)%nsnow                           
            tmp_nimpur                              = CROCUS81_struc(n)%nimpur                          
            tmp_SNOWRES_opt                         = CROCUS81_struc(n)%SNOWRES_opt                     
            tmp_OMEB_BOOL                           = CROCUS81_struc(n)%OMEB_BOOL                       
            tmp_GLACIER_BOOL                        =  CROCUS81_struc(n)%GLACIER_BOOL   
            tmp_HIMPLICIT_WIND_opt                  = CROCUS81_struc(n)%HIMPLICIT_WIND_opt              
            tmp_PTSTEP                              = CROCUS81_struc(n)%PTSTEP                          
            tmp_UREF                                = CROCUS81_struc(n)%UREF                            
            tmp_SLOPE                               = CROCUS81_struc(n)%crocus81(t)%SLOPE ! reading from LDT output
            tmp_ZREF                                = CROCUS81_struc(n)%ZREF                            
            tmp_Z0NAT                               = CROCUS81_struc(n)%Z0NAT                           
            tmp_Z0EFF                               = CROCUS81_struc(n)%Z0EFF                           
            tmp_Z0HNAT                              = CROCUS81_struc(n)%Z0HNAT                          
            tmp_ALB(:)                              = CROCUS81_struc(n)%crocus81(t)%ALB(:) ! reading from LDT output 
            tmp_use_monthly_albedo_map              = CROCUS81_struc(n)%use_monthly_albedo_map 
            tmp_SOILCOND                            = CROCUS81_struc(n)%crocus81(t)%SOILCOND 
            tmp_D_G                                 = CROCUS81_struc(n)%D_G                             
            tmp_PERMSNOWFRAC                        = CROCUS81_struc(n)%crocus81(t)%PERMSNOWFRAC  !read from LDT output                
            tmp_SNOWDRIFT_opt                       = CROCUS81_struc(n)%SNOWDRIFT_opt                   
            tmp_SNOWDRIFT_SUBLIM_BOOL               = CROCUS81_struc(n)%SNOWDRIFT_SUBLIM_BOOL           
            tmp_SNOW_ABS_ZENITH_BOOL                = CROCUS81_struc(n)%SNOW_ABS_ZENITH_BOOL            
            tmp_SNOWMETAMO_opt                      = CROCUS81_struc(n)%SNOWMETAMO_opt                  
            tmp_SNOWRAD_opt                         = CROCUS81_struc(n)%SNOWRAD_opt                     
            tmp_ATMORAD_BOOL                        = CROCUS81_struc(n)%ATMORAD_BOOL                    
            tmp_IMPWET(:)                           = CROCUS81_struc(n)%IMPWET(:)                          
            tmp_IMPDRY(:)                           = CROCUS81_struc(n)%IMPDRY(:)                          
            tmp_SNOWFALL_opt                        = CROCUS81_struc(n)%SNOWFALL_opt                    
            tmp_SNOWCOND_opt                        = CROCUS81_struc(n)%SNOWCOND_opt                    
            tmp_SNOWHOLD_opt                        = CROCUS81_struc(n)%SNOWHOLD_opt                    
            tmp_SNOWCOMP_opt                        = CROCUS81_struc(n)%SNOWCOMP_opt                    
            tmp_SNOWZREF_opt                        = CROCUS81_struc(n)%SNOWZREF_opt                    
            tmp_SNOWCOMPACT_BOOL                    = CROCUS81_struc(n)%SNOWCOMPACT_BOOL                
            tmp_SNOWMAK_BOOL                        = CROCUS81_struc(n)%SNOWMAK_BOOL                    
            tmp_SNOWTILLER_BOOL                     = CROCUS81_struc(n)%SNOWTILLER_BOOL                 
            tmp_SELF_PROD_BOOL                      = CROCUS81_struc(n)%SELF_PROD_BOOL                  
            tmp_SNOWMAK_PROP_BOOL                   = CROCUS81_struc(n)%SNOWMAK_PROP_BOOL               
            tmp_PRODSNOWMAK_BOOL                    = CROCUS81_struc(n)%PRODSNOWMAK_BOOL                
            tmp_SLOPE_DIR                           = CROCUS81_struc(n)%crocus81(t)%SLOPE_DIR ! read from LDT output
            tmp_SAND                                = CROCUS81_struc(n)%crocus81(t)%SAND                            
            tmp_SILT                                = CROCUS81_struc(n)%crocus81(t)%SILT                            
            tmp_CLAY                                = CROCUS81_struc(n)%crocus81(t)%CLAY                            
            tmp_POROSITY                            = CROCUS81_struc(n)%crocus81(t)%POROSITY            

            if(LIS_rc%lsm.ne."none") then
               tmp_TG                                  = CROCUS81_struc(n)%crocus81(t)%TG
               tmp_XWGI                                = CROCUS81_struc(n)%crocus81(t)%XWGI
               tmp_XWG                                 = CROCUS81_struc(n)%crocus81(t)%XWG
            else
            ! For the stand-alone version, we need to provide some default values for XWGI and XWG
              tmp_XWGI = 0.0 ! MN set to zero 
              tmp_XWG  = tmp_POROSITY * 0.8 ! MN assume volumetric soil water content of the snow covered
                                            ! ground is 80% of POROSITY (For Col de Porte it is between 72-85)
              tmp_TG   = 273.15             ! no energy exchange between snow and soil
            endif
            ! get state variables
            tmp_SNOWSWE(:)    = CROCUS81_struc(n)%crocus81(t)%SNOWSWE(:)   
            tmp_SNOWRHO(:)    = CROCUS81_struc(n)%crocus81(t)%SNOWRHO(:)   
            tmp_SNOWHEAT(:)   = CROCUS81_struc(n)%crocus81(t)%SNOWHEAT(:)  
            tmp_SNOWALB       = CROCUS81_struc(n)%crocus81(t)%SNOWALB   
            tmp_SNOWGRAN1(:)  = CROCUS81_struc(n)%crocus81(t)%SNOWGRAN1(:) 
            tmp_SNOWGRAN2(:)  = CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2(:) 
            tmp_SNOWHIST(:)   = CROCUS81_struc(n)%crocus81(t)%SNOWHIST(:)  
            tmp_SNOWAGE(:)    = CROCUS81_struc(n)%crocus81(t)%SNOWAGE(:)   
            tmp_SNOWLIQ(:)    = CROCUS81_struc(n)%crocus81(t)%SNOWLIQ(:)   
            tmp_SNOWTEMP(:)   = CROCUS81_struc(n)%crocus81(t)%SNOWTEMP(:)  
            tmp_SNOWDZ(:)     = CROCUS81_struc(n)%crocus81(t)%SNOWDZ(:)    
            tmp_GRNDFLUX      = CROCUS81_struc(n)%crocus81(t)%GRNDFLUX  
            tmp_SNDRIFT       = CROCUS81_struc(n)%crocus81(t)%SNDRIFT   
            tmp_RI_n          = CROCUS81_struc(n)%crocus81(t)%RI_n      
            tmp_CDSNOW        = CROCUS81_struc(n)%crocus81(t)%CDSNOW    
            tmp_USTARSNOW     = CROCUS81_struc(n)%crocus81(t)%USTARSNOW 
            tmp_CHSNOW        = CROCUS81_struc(n)%crocus81(t)%CHSNOW    
            tmp_SNOWMAK_dz    = CROCUS81_struc(n)%crocus81(t)%SNOWMAK_dz
 

! Compute GLACIER_BOOL 
! TODO check the threshold value
! use a threshold value of 0.2  
tmp_GLACIER_BOOL = .False.
if (tmp_PERMSNOWFRAC.gt.0.2)then
      tmp_GLACIER_BOOL = .True.
endif

!------------------------------------------------------------------------------------ 
! call model physics 
            call crocus_driver(tmp_n                 , & ! IN    - nest id [-]
                               tmp_nsnow             , & ! IN    - number of snow layer []
                               tmp_nimpur            , & ! IN    - number of impurtites []
                               tmp_year              , & ! IN    - year of the currrent time step [-]
                               tmp_month             , & ! IN    - month of the current time step [-]
                               tmp_day               , & ! IN    - day of the current time step [-]
                               tmp_hour              , & ! IN    - hour of the current time step [-]
                               tmp_minute            , & ! IN    - minute of the current time step [-]
                               tmp_SNOWRES_opt       , & ! IN    - SNOWRES_opt  = ISBA-SNOW3L turbulant exchange option
!   'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!   'RIL' = Limit Richarson number under very stable
! conditions (currently testing)
!   'M98'  = Martin et Lejeune 1998 : older computation for turbulent fluxes coefficents in Crocus []
                               tmp_OMEB_BOOL         , & !       - !True = coupled to MEB. surface fluxes are IMPOSED
! as an upper boundary condition to the explicit snow schemes. 
!If = False, then energy budget and fluxes are computed herein. [-]
                               tmp_GLACIER_BOOL      , & ! IN    - !True = Over permanent snow and ice, initialise WGI=WSAT, Hsnow>=10m and allow 0.8<SNOALB<0.85
! False = No specific treatment [-]
                               tmp_HIMPLICIT_WIND_opt, & ! IN    - wind implicitation option  'OLD' = direct , 'NEW' = Taylor serie, order 1 [-]
                               tmp_SNOWSWE           , & ! INOUT - Snow layer(s) liquid Water Equivalent (SWE:kg m-2) [kg/m2]
                               tmp_SNOWRHO           , & ! INOUT - Snow layer(s) averaged density (kg/m3) [kg/m3]
                               tmp_SNOWHEAT          , & ! INOUT - Snow layer(s) heat content (J/m2) [J/m2]
                               tmp_SNOWALB           , & ! INOUT - snow surface albedo [-]
                               tmp_SNOWGRAN1         , & ! INOUT - Snow layers grain feature 1 [-]
                               tmp_SNOWGRAN2         , & ! INOUT - Snow layer grain feature 2 [-]
                               tmp_SNOWHIST          , & ! INOUT - Snow layer grain historical parameter (only for non dendritic snow) (-) in {0-5} [-]
                               tmp_SNOWAGE           , & ! INOUT - Age since snowfall (day) [day]
                               tmp_PTSTEP            , & ! IN    - time step of the integration [s]
                               tmp_PPS               , & ! IN    - pressure at atmospheric model surface (Pa) [Pa]
                               tmp_SRSNOW            , & ! IN    - snow rate (SWE) [kg/(m2 s)] [kg/(m2 s)]
                               tmp_RRSNOW            , & ! IN    - rain rate [kg/(m2 s)] [kg/(m2 s)]
                               tmp_TA                , & ! IN    - atmospheric temperature at level za (K) [K]
                               tmp_TG                , & ! IN    - Surface soil temperature (effective temperature the of layer lying below snow) (K)  (for snowcro.F90 we only use the surface layer ZP_TG(:,1))  (#nsoil depends on 2-L, 3-L DIF) [K]
                               tmp_SW_RAD            , & ! IN    - incoming solar radiation (W/m2) [W/m2]
                               tmp_QA                , & ! IN    - atmospheric specific humidity at level za [-]
                               tmp_Wind_E            , & ! IN    - Eastward Wind [m/s]
                               tmp_Wind_N            , & ! IN    - Northward Wind [m/s]
                               tmp_LW_RAD            , & ! IN    - atmospheric infrared radiation (W/m2) [W/m2]
                               tmp_UREF              , & ! IN    - reference height of the wind [m]
                               tmp_SLOPE             , & ! IN    - angle between the normal to the surface and the vertical  (MN:  replaced PDIRCOSZW with slope and computed the cosine in the driver) [-]
                               tmp_ZREF              , & ! IN    - Reference height of the first atmospheric level (m) [m]
                               tmp_Z0NAT             , & ! IN    - grid box average roughness length (m) (roughness length for momentum) [m]
                               tmp_Z0EFF             , & ! IN    - roughness length for momentum (modd_diagn.F90 effective roughness length for heat(!?)) [m]
                               tmp_Z0HNAT            , & ! IN    - grid box average roughness length for heat [m]
                               tmp_ALB               , & ! IN    - soil/vegetation albedo [-]
                               !tmp_SOILCOND          , & ! IN    - soil thermal conductivity (W m-1 K-1) [W /(m K)]! MN : for now will be computed in the driver using the sand fraction
                               tmp_D_G               , & ! IN    - !Assumed first soil layer thickness (m)
!Used to calculate ground/snow heat flux   (D_G(:,1)) [m]
                               tmp_SNOWLIQ           , & ! INOUT - Snow layer(s) liquid water content (m) [m]
                               tmp_SNOWTEMP          , & ! INOUT - Snow layer(s) temperature (K) [K]
                               tmp_SNOWDZ            , & ! INOUT - Snow layer(s) thickness (m) [m]
                               tmp_THRUFAL           , & ! OUT   - Rate that liquid water leaves snow pack: paritioned into soil infiltration/runoff  by ISBA [kg/(m2 s)] [kg/(m2 s)]
                               tmp_GRNDFLUX          , & ! INOUT - Soil/snow interface heat flux (W/m2) [W/m2]
                               tmp_SNDRIFT           , & ! OUT   - Blowing snow sublimation (kg/m2/s) NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592) [kg/m2/s]
                               tmp_RI_n              , & ! OUT   - Richardson number (-)  NOTE: RI has not been initialized in CALL_MODEL (If not OMED initalized to undefined in the snow3L_isba.F90) [-]
                               tmp_EMISNOW           , & ! OUT   - Snow surface emissivity (initialize to XEMISSN
!XEMISSN  comes from MODD_SNOW_PAR) [-]
                               tmp_CDSNOW            , & ! OUT   - Drag coefficient for momentum over snow (-) [-]
                               tmp_USTARSNOW         , & ! OUT   - Friction velocity over snow (m/s); [m/s]
                               tmp_CHSNOW            , & ! OUT   - Drag coefficient for heat over snow  (-) [-]
                               tmp_SNOWHMASS         , & ! OUT   - Heat content change due to mass changes in snowpack: for budget calculations only. [J/m2]
                               tmp_QS                , & ! OUT   - Surface humidity (kg/kg) [kg/kg]
                               tmp_PERMSNOWFRAC      , & ! IN    - Fraction of permanet snow/ice [-]
                               tmp_LAT               , & ! IN    - Latitude in decimal degree  (latitude (degrees +North)) [degrees]
                               tmp_LON               , & ! IN    - Longitude in decimal year    (longitude (degrees +East)) [degrees]
                               tmp_SNOWDRIFT_opt     , & ! IN    - Mechanical transformation of snow grain and compaction + effect of wind on falling snow properties
!	'NONE': No snowdrift scheme
!	'DFLT': falling snow falls as purely dendritic
! 	'GA01': Gallee et al 2001
!	'VI13': Vionnet et al 2013 [-]
                               tmp_SNOWDRIFT_SUBLIM_BOOL, & ! IN    - Logicals for snowdrift sublimation [-]
                               tmp_SNOW_ABS_ZENITH_BOOL, & ! IN    - Activate parametrization of solar absorption for polar regions (If True modify solar absorption as a function of solar zenithal angle) [-]
                               tmp_SNOWMETAMO_opt    , & ! IN    - Metamorphism scheme: B92 (historical version, Brun et al 92), C13, T07, F06 (see Carmagnola et al 2014) [-]
                               tmp_SNOWRAD_opt       , & ! IN    - Radiative transfer scheme. HSNOWRAD=B92 Brun et al 1992.  HSNOWRAD=T17 (Tuzet et al. 2017) (Libois et al. 2013) TARTES with impurities content scheme [-]
                               tmp_ATMORAD_BOOL      , & ! IN    - Activate atmotartes scheme  (default=.FALSE. # This option is not stable yet, but it is supposed to compute the direct/diffuse ratio directly from atmospheric informations (AOD, Ozone column, Water column..)) [-]
                               tmp_IMPWET            , & ! IN    - [init_surf_atmn.F90   --> wet deposit coefficient for each impurity type (g/m_/s)  ,   snowcro.F90 --> Dry and wet deposit coefficient from Forcing File(g/m_/s)   (# 1553 You can either feed the model with prescribed and constant deposition fluxes or introduce a wet and dry deposition field directly in the forcing file. ) [g/m_/s]
                               tmp_IMPDRY            , & ! IN    - [init_surf_atmn.F90   --> wet deposit coefficient for each impurity type (g)  ,   snowcro.F90 --> Dry and wet deposit coefficient from Forcing File(g/m_/s) [g/m_/s]
                               tmp_SNOWFALL_opt      , & ! IN    - New options for multiphysics version (Cluzet et al 2016). Falling snow scheme: V12 (Vionnet et al. 2012) , A76 (Anderson 1976), S02 (Lehning and al. 2002), P75 (Pahaut 1975) [-]
                               tmp_SNOWCOND_opt      , & ! IN    - Thermal conductivity scheme: Y81 (Yen 1981), I02 (Boone et al. 2002) C11 (Calonne et al. 2011) [-]
                               tmp_SNOWHOLD_opt      , & ! IN    - liquid water content scheme: B92 (Brun et al. 1992) O04 (Oleson et al., 2004) S02 (SNOWPACK, Lehning et al, 2002) B02 (ISBA_ES, Boone et al. 2002) [-]
                               tmp_SNOWCOMP_opt      , & ! IN    - B92 snow compaction basis version and B93 for slightly different parameters   (NOTE: IN THE CODE S14, T11, ) [-]
                               tmp_SNOWZREF_opt      , & ! IN    - reference height is constant or variable from the snow surface: CST (constant from snow surface, i.e. Col de Porte) or VAR (variable from snow surface = snow depth has to be removed from reference height) [-]
                               tmp_SNOWMAK_dz        , & ! IN    - Snowmaking thickness (m) [m]
                               tmp_SNOWCOMPACT_BOOL  , & ! IN    - Snowmaking and Grooming options [-]
                               tmp_SNOWMAK_BOOL      , & ! IN    - Snowmaking and Grooming options [-]
                               tmp_SNOWTILLER_BOOL   , & ! IN    - Snowmaking and Grooming options [-]
                               tmp_SELF_PROD_BOOL    , & ! IN    - Snowmaking and Grooming options [-]
                               tmp_SNOWMAK_PROP_BOOL , & ! IN    - Snowmaking and Grooming options [-]
                               tmp_PRODSNOWMAK_BOOL  , & ! INOUT - Snowmaking and Grooming options [-]
                               tmp_SLOPE_DIR        , &  ! IN    - !Typical slope aspect in the grid  (clockwise from N) [Radians]
                               tmp_SAND              , & ! IN    - Soil sand fraction (-) [-]
                               tmp_SILT              , & ! IN    - Soil silt fraction (-) [-]
                               tmp_CLAY              , & ! IN    - Soil clay fraction (-) [-]
                               tmp_POROSITY          ,&  ! IN    - Soil porosity (m3 m-3) [m3/m3]
                               tmp_XWGI              ,&  ! IN    - soil volumetric frozen water content
                               tmp_XWG                ,&  ! IN    - soil volumetric liquid water content
                               tmp_use_monthly_albedo_map)! ,& ! IN    - Boolean option to partition total precipitation into snowfall and rainfall using Jordan 1991 
    
            ! save state variables from local variables to global variables
            CROCUS81_struc(n)%crocus81(t)%SNOWSWE(:)    = tmp_SNOWSWE(:)   
            CROCUS81_struc(n)%crocus81(t)%SNOWRHO(:)    = tmp_SNOWRHO(:)   
            CROCUS81_struc(n)%crocus81(t)%SNOWHEAT(:)   = tmp_SNOWHEAT(:)  
            CROCUS81_struc(n)%crocus81(t)%SNOWALB       = tmp_SNOWALB     
            CROCUS81_struc(n)%crocus81(t)%SNOWGRAN1(:)  = tmp_SNOWGRAN1(:) 
            CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2(:)  = tmp_SNOWGRAN2(:) 
            CROCUS81_struc(n)%crocus81(t)%SNOWHIST(:)   = tmp_SNOWHIST(:)  
            CROCUS81_struc(n)%crocus81(t)%SNOWAGE(:)    = tmp_SNOWAGE(:)   
            CROCUS81_struc(n)%crocus81(t)%SNOWLIQ(:)    = tmp_SNOWLIQ(:)   
            CROCUS81_struc(n)%crocus81(t)%SNOWTEMP(:)   = tmp_SNOWTEMP(:)  
            CROCUS81_struc(n)%crocus81(t)%SNOWDZ(:)     = tmp_SNOWDZ(:)    
            CROCUS81_struc(n)%crocus81(t)%GRNDFLUX      = tmp_GRNDFLUX     
            CROCUS81_struc(n)%crocus81(t)%SNDRIFT       = tmp_SNDRIFT      
            CROCUS81_struc(n)%crocus81(t)%RI_n          = tmp_RI_n         
            CROCUS81_struc(n)%crocus81(t)%CDSNOW        = tmp_CDSNOW       
            CROCUS81_struc(n)%crocus81(t)%USTARSNOW     = tmp_USTARSNOW    
            CROCUS81_struc(n)%crocus81(t)%CHSNOW        = tmp_CHSNOW       
            CROCUS81_struc(n)%crocus81(t)%SNOWMAK_dz    = tmp_SNOWMAK_dz   
    
            ! save output variables from local variables to global variables
            CROCUS81_struc(n)%crocus81(t)%THRUFAL      = tmp_THRUFAL     
            CROCUS81_struc(n)%crocus81(t)%EMISNOW      = tmp_EMISNOW     
            CROCUS81_struc(n)%crocus81(t)%SNOWHMASS    = tmp_SNOWHMASS   
            CROCUS81_struc(n)%crocus81(t)%QS           = tmp_QS          
            tmp = 0.0 
            tmp = sum(tmp_SNOWDZ)
            CROCUS81_struc(n)%crocus81(t)%SD_1D        = tmp
            tmp = 0.0
            tmp = sum(tmp_SNOWSWE)
            CROCUS81_struc(n)%crocus81(t)%SWE_1D       = tmp
           

! TODO   ! Slope correction --> NOTE: the correction should be apply when we write the output
            ! As a convention, only vertical incoming and outgoing fluxes are 
            ! provided to and from the Crcous model; the correction of these terms
            ! according to the local slope is carried out here.
            ! Similarly, variables such as total snow depth, total snow water
            ! equivalent, and the corresponding variables for each layer
            ! are output by the model in terms of their vertical component,
            ! i.e. projected vertically 
 
            ![ 1] output variable: SNOWSWE (unit=kg/m2). *** Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
            do i=1, CROCUS81_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWLIQPROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWSWE(i), &
                                                  vlevel=i, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 2] output variable: SNOWHEAT (unit=J/m2). *** Snow layer(s) Heat content (J/m2)
            do i=1, CROCUS81_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWHEATCONTENTPROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWHEAT(i), &
                                                  vlevel=i, unit="J m-2", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 3] output variable: SNOWRHO (unit=kg/m3). *** Snow layer(s) Density (kg/m3)
            do i=1, CROCUS81_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWRHOPROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWRHO(i), &
                                                  vlevel=i, unit="kg m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 4] output variable: SNOWALB (unit=-). *** Snow surface albedo
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWALB, value = CROCUS81_struc(n)%crocus81(t)%SNOWALB, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 5] output variable: SNOWGRAN1 (unit=-). *** Snow layer(s) grain parameter 1
            do i=1, CROCUS81_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWGRAIN1PROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWGRAN1(i), &
                                                  vlevel=i, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 6] output variable: SNOWGRAN2 (unit=-). *** Snow layer(s) grain parameter 2
            do i=1, CROCUS81_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWGRAIN2PROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2(i), &
                                                  vlevel=i, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 7] output variable: SNOWHIST (unit=-). *** Snow layer(s) Historical parameter (-) in {0-5}
            do i=1, CROCUS81_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWHISTPROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWHIST(i), &
                                                  vlevel=i, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 8] output variable: SNOWAGE (unit=day). *** Snow layer(s) Age since snowfall (day)
            do i=1, CROCUS81_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWAGEPROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWAGE(i), &
                                                  vlevel=i, unit="day", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 9] output variable: SNOWLIQ (unit=m). *** Snow layer(s) liquid water content (m)
            do i=1, CROCUS81_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWLIQCONTENTPROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWLIQ(i), &
                                                  vlevel=i, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 10] output variable: SNOWTEMP (unit=K). *** Snow layer(s) temperature (K)
            do i=1, CROCUS81_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWTEMPPROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWTEMP(i), &
                                                  vlevel=i, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 11] output variable: SNOWDZ (unit=m). *** Snow layer(s) thickness (m)
            do i=1, CROCUS81_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWHIGHTPROF, value = CROCUS81_struc(n)%crocus81(t)%SNOWDZ(i), &
                                                  vlevel=i, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 12] output variable: THRUFAL (unit=kg/(m2 s)). *** Rate that liquid water leaves snow pack: paritioned into soil infiltration/runoff  by ISBA [kg/(m2 s)]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWQS, value = CROCUS81_struc(n)%crocus81(t)%THRUFAL, &
                                              vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
            
            ![ 13] output variable: GRNDFLUX (unit=W/m2). *** Soil/snow interface heat flux (W/m2) (If MEB input else initalize to 0)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWSOILHEATFLUX, value = CROCUS81_struc(n)%crocus81(t)%GRNDFLUX, &
                                              vlevel=1, unit="W m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 14] output variable: SNDRIFT (unit=kg/m2/s). *** Blowing snow sublimation (kg/m2/s) NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BLOWINGSNOWSUBLIM, value = CROCUS81_struc(n)%crocus81(t)%SNDRIFT, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 15] output variable: RI_n (unit=-). *** Richardson number (-)  NOTE: ZP_RI has not been initialized in CALL_MODEL (If not OMED initalized to undefined in the snow3L_isba.F90)  
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWRECHARD, value = CROCUS81_struc(n)%crocus81(t)%RI_n, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 16] output variable: EMISNOW (unit=-). *** Snow surface emissivity (initialize to XEMISSN
! XEMISSN comes from MODD_SNOW_PAR)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWEMISS, value = CROCUS81_struc(n)%crocus81(t)%EMISNOW, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 17] output variable: CDSNOW (unit=-). *** Drag coefficient for momentum over snow (-)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWCM, value = CROCUS81_struc(n)%crocus81(t)%CDSNOW, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 18] output variable: USTARSNOW (unit=m/s). *** Friction velocity over snow (m/s);
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWSHEARVLOCITY, value = CROCUS81_struc(n)%crocus81(t)%USTARSNOW, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 19] output variable: CHSNOW (unit=-). *** Drag coefficient for heat over snow  (-)                          
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWHEATDRAG, value = CROCUS81_struc(n)%crocus81(t)%CHSNOW, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 20] output variable: SNOWHMASS (unit=J/m2). *** Heat content change due to mass changes in snowpack: for budget calculations only.
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDELTAHEAT, value = CROCUS81_struc(n)%crocus81(t)%SNOWHMASS, &
                                              vlevel=1, unit="J m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 21] output variable: QS (unit=kg/kg). *** surface humidity (kg/kg)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWSURFACEQ, value = CROCUS81_struc(n)%crocus81(t)%QS, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 22] output variable: PERMSNOWFRAC (unit=-). *** Fraction of permanet snow/ice
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GLACIERFRACTION, value = CROCUS81_struc(n)%crocus81(t)%PERMSNOWFRAC, &
                                              vlevel=1, unit="-", direction="", surface_type = LIS_rc%lsm_index)

            !TODO Writing out the snow profile is a bottleneck and significantly increases the simulation time. For now 
            !     option added to the OUTPUT.TBL to be able to write out the cumulative value.  
            ![ 23] output variable: SD_1D (unit=m). *** Total snow depth, temporally added 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SnowDepth, value = CROCUS81_struc(n)%crocus81(t)%SD_1D, &                
                                                vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)                
              
            ![ 24] output variable: SWE_1D (unit=kg/m2). *** Total SWE, temporally added                                         
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWE, value = CROCUS81_struc(n)%crocus81(t)%SWE_1D, &              
                                                vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)      

            ! reset forcing variables to zeros
            CROCUS81_struc(n)%crocus81(t)%PPS = 0.0
            CROCUS81_struc(n)%crocus81(t)%SRSNOW = 0.0
            CROCUS81_struc(n)%crocus81(t)%RRSNOW = 0.0
            CROCUS81_struc(n)%crocus81(t)%TA = 0.0
            CROCUS81_struc(n)%crocus81(t)%SW_RAD = 0.0
            CROCUS81_struc(n)%crocus81(t)%QA = 0.0
            CROCUS81_struc(n)%crocus81(t)%Wind_E = 0.0
            CROCUS81_struc(n)%crocus81(t)%Wind_N = 0.0
            CROCUS81_struc(n)%crocus81(t)%LW_RAD = 0.0
        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        CROCUS81_struc(n)%forc_count = 0 
    endif ! end of alarmCheck loop 
    
    deallocate( tmp_SNOWSWE )
    deallocate( tmp_SNOWRHO )
    deallocate( tmp_SNOWHEAT )
    deallocate( tmp_SNOWGRAN1 )
    deallocate( tmp_SNOWGRAN2 )
    deallocate( tmp_SNOWHIST )
    deallocate( tmp_SNOWAGE )
    deallocate( tmp_ALB )
    deallocate( tmp_SNOWLIQ )
    deallocate( tmp_SNOWTEMP )
    deallocate( tmp_SNOWDZ )
    deallocate( tmp_IMPWET )
    deallocate( tmp_IMPDRY )
end subroutine Crocus81_main
