! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION SfcFlx_lwradatm (T_a, e_a, cl_tot, cl_low)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the long-wave radiation flux from the atmosphere
!  as function of air temperature, water vapour pressure and cloud fraction. 
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "SfcFlx".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = ireals), INTENT(IN) ::   &
  T_a                               , & ! Air temperature [K]
  e_a                               , & ! Water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
  cl_tot                            , & ! Total cloud cover [0,1]
  cl_low                                ! Lowe-level cloud cover [0,1]


!  Local parameters  

!  Coefficients in the empirical formulation  
!  developed at the Main Geophysical Observatory (MGO), St. Petersburg, Russia.
REAL (KIND = ireals), PARAMETER ::   &
  c_lmMGO_1    = 43.057924_ireals  , & ! Empirical coefficient 
  c_lmMGO_2    = 540.795_ireals        ! Empirical coefficient 
!  Temperature-dependent cloud-correction coefficients in the MGO formula
INTEGER (KIND = iintegers), PARAMETER :: &
  nband_coef = 6_iintegers                 ! Number of temperature bands
REAL (KIND = ireals), PARAMETER, DIMENSION (nband_coef) ::      &
  corr_cl_tot     = (/0.70_ireals, 0.45_ireals, 0.32_ireals,    & 
                      0.23_ireals, 0.18_ireals, 0.13_ireals/) , & ! Total clouds
  corr_cl_low     = (/0.76_ireals, 0.49_ireals, 0.35_ireals,    & 
                      0.26_ireals, 0.20_ireals, 0.15_ireals/) , & ! Low-level clouds
  corr_cl_midhigh = (/0.46_ireals, 0.30_ireals, 0.21_ireals,    & 
                      0.15_ireals, 0.12_ireals, 0.09_ireals/)     ! Mid- and high-level clouds
REAL (KIND = ireals), PARAMETER ::   &
  T_low  = 253.15_ireals           , & ! Low-limit temperature in the interpolation formula [K]
  del_T  = 10.0_ireals                 ! Temperature step in the interpolation formula [K]

!  Coefficients in the empirical water-vapour correction function 
!  (see Fung et al. 1984, Zapadka and Wozniak 2000, Zapadka et al. 2001). 
REAL (KIND = ireals), PARAMETER ::     &
  c_watvap_corr_min = 0.6100_ireals  , & ! Empirical coefficient (minimum value of the correction function)
  c_watvap_corr_max = 0.7320_ireals  , & ! Empirical coefficient (maximum value of the correction function)
  c_watvap_corr_e   = 0.0050_ireals      ! Empirical coefficient [(N m^{-2})^{-1/2}]

!  Local variables of type INTEGER
INTEGER (KIND = iintegers) :: &
  i                             ! Loop index

!  Local variables of type REAL
REAL (KIND = ireals) ::    &
  c_cl_tot_corr          , & ! The MGO cloud correction coefficient, total clouds
  c_cl_low_corr          , & ! The MGO cloud correction coefficient, low-level clouds
  c_cl_midhigh_corr      , & ! The MGO cloud correction coefficient, mid- and high-level clouds
  T_corr                 , & ! Temperature used to compute the MGO cloud correction [K]
  f_wvpres_corr          , & ! Correction function with respect to water vapour
  f_cloud_corr               ! Correction function with respect to cloudiness
 
!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Water-vapour correction function
  f_wvpres_corr = c_watvap_corr_min + c_watvap_corr_e*SQRT(e_a) 
  f_wvpres_corr = MIN(f_wvpres_corr, c_watvap_corr_max)

! Cloud-correction coefficients using the MGO formulation with linear interpolation 
IF(T_a.LT.T_low) THEN 
  c_cl_tot_corr     = corr_cl_tot(1)   
  c_cl_low_corr     = corr_cl_low(1)
  c_cl_midhigh_corr = corr_cl_midhigh(1)
ELSE IF(T_a.GE.T_low+(nband_coef-1_iintegers)*del_T) THEN
  c_cl_tot_corr     = corr_cl_tot(nband_coef)   
  c_cl_low_corr     = corr_cl_low(nband_coef)
  c_cl_midhigh_corr = corr_cl_midhigh(nband_coef)
ELSE 
  T_corr = T_low
  DO i=1, nband_coef-1
    IF(T_a.GE.T_corr.AND.T_a.LT.T_corr+del_T) THEN 
      c_cl_tot_corr = (T_a-T_corr)/del_T
      c_cl_low_corr = corr_cl_low(i) + (corr_cl_low(i+1)-corr_cl_low(i))*c_cl_tot_corr
      c_cl_midhigh_corr = corr_cl_midhigh(i) + (corr_cl_midhigh(i+1)-corr_cl_midhigh(i))*c_cl_tot_corr
      c_cl_tot_corr = corr_cl_tot(i) + (corr_cl_tot(i+1)-corr_cl_tot(i))*c_cl_tot_corr
    END IF 
    T_corr = T_corr + del_T
  END DO
END IF
! Cloud correction function
IF(cl_low.LT.0._ireals) THEN  ! Total cloud cover only 
  f_cloud_corr = 1._ireals + c_cl_tot_corr*cl_tot*cl_tot
ELSE                          ! Total and low-level cloud cover
  f_cloud_corr = (1._ireals + c_cl_low_corr*cl_low*cl_low)  &
               * (1._ireals + c_cl_midhigh_corr*(cl_tot*cl_tot-cl_low*cl_low))
END IF

! Long-wave radiation flux [W m^{-2}]

!  The MGO formulation  
!_nu The MGO formulation  
!_nu SfcFlx_lwradatm = -SfcFlx_lwradatm*c_lwrad_emis  &
!_nu                 * (c_lmMGO_1*SQRT(tpsf_C_StefBoltz*T_a**4_iintegers)-c_lmMGO_2)
!_nu 

!  "Conventional" formulation  
!  (see Fung et al. 1984, Zapadka and Wozniak 2000, Zapadka et al. 2001)  
SfcFlx_lwradatm = -c_lwrad_emis*tpsf_C_StefBoltz*T_a**4_iintegers  &
                * f_wvpres_corr*f_cloud_corr

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION SfcFlx_lwradatm

