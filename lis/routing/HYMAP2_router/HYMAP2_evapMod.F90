!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !REVISION HISTORY: 
! 14 Dec 2014: Augusto Getirana, Initial implementation
!-------------------------END NOTICE -- DO NOT EDIT----------------------- 
! Augusto Getirana
! 20 April 2017
! NASA GSFC
! This routine computes ET based on the Penman-Monteith equation, as presented in Shuttleworth (1993)  
! Shuttleworth, W.J., 1993. Evaporation. In: Maidment, D.R., Handbook of Hydrology. McGraw-Hill. New York.
! It assumes:
! 1. Temperature in degrees Celsius
! 2. Pressure in kPa
! 3. Neglects G 
! 4. Neglects vegetative surface roughness 

!LIS forcings used in this module
!Tair:          1        1     K       # Near surface air temperature
!Qair:          1        1     kg/kg   # Near surface specific humidity
!SWdown:        1        1     W/m2    # Incident shortwave radiation (total)
!LWdown:        1        1     W/m2    # Incident longwave radiation
!Wind_E:        1        1     W/m2    # Eastward wind
!Wind_N:        1        1     m/s     # Northward wind
!Psurf:         1        1     Pa      # Surface pressure

!CHECK:
!UNITS OF TEMPERATURE, PRESSURE, AIR HUMIDITY, NET RADIATION AND WIND SPEED IN LIS OUTPUTS
module HYMAP2_evapMod

  public :: HYMAP2_evap_main

  ! === Physical, Universal and Model constants ====
  real, parameter :: mesp             = 1000.0     ! water specific mass [kg m-3]
  real, parameter :: cspa             = 0.001013   ! specific heat of humid air [MJ kg-1 C-1]
  real, parameter :: xesatcf_a        = 0.6112     ! kPa
  real, parameter :: xesatcf_b        =  17.67     ! -
  real, parameter :: xesatcf_c        = 243.50     ! K
  real, parameter :: zeroC            = 273.15     ! K
  
  real, parameter :: xz0m_sea         = 0.000240   ! m initial value 
  real, parameter :: xz0m_min_sea     = 0.00001    ! m

contains
  !=======================================================================
  subroutine HYMAP2_evap_main(evaptype,n,nseqall,mis,outlet,pres,tair,qair,wind,qnet,evap,ewat,edif)

    use LIS_logMod,     only : LIS_logunit   
     
    implicit none
    
    integer,      intent(in)    :: evaptype
    integer,      intent(in)    :: n
    integer,      intent(in)    :: nseqall           ! length of 1D sequnece for river and mouth  
    !real,         intent(in)    :: dt                ! time step length [s]
    real,         intent(in)    :: mis               ! missing data
    integer,      intent(in)    :: outlet(nseqall)   ! outlet flag: 0 - river; 1 - ocean
    !Meteorological forcings
    real,         intent(in)    :: pres(nseqall)     ! surface pressure (Pa)
    real,         intent(in)    :: tair(nseqall)     ! Air temperature at reference level (K)
    real,         intent(in)    :: qair(nseqall)     ! Specific humidity at reference level (kg/kg)
    real,         intent(in)    :: wind(nseqall)     ! Wind speed at reference level (m s-1) 
    real,         intent(in)    :: qnet(nseqall)     ! Net radiation [MJ s-1 m-2]
    !LSM output
    real,         intent(in)    :: evap(nseqall)     ! total evapotranspiration (kg m-2 s-1) derived from LSM
    !HYMAP2 outputs
    real,         intent(inout) :: ewat(nseqall)     ! Surface water potential evaporation (kg m-2 s-1)
    real,         intent(out)   :: edif(nseqall)     ! Differential evaporation (kg m-2 s-1) 

    integer                     :: ic

    do ic=1,nseqall 
      if(outlet(ic)==mis)cycle
      if(evaptype==1)then
        call HYMAP2_penman_monteith(pres(ic),tair(ic),wind(ic),qair(ic),qnet(ic),ewat(ic))
      else
         write(LIS_logunit,*)"HYMAP routing model evaporation flag: unknown value"
         stop
      endif
      edif(ic)=max(0.,ewat(ic)-evap(ic))
    enddo

  end subroutine HYMAP2_evap_main
  !=======================================================================
  subroutine HYMAP2_penman_monteith(pres,tair,wind,qair,qnet,ewat)
    implicit none
    real, intent(in)  :: pres,tair,wind,qair,qnet
    real, intent(out) :: ewat
    
    real              :: esat    !vapor pressure at saturation [kPa]
    real              :: del     !change rate of vapor pressure at saturation [kPa C-1]
    real              :: mspa    !specific mass of air [kg m-3]
    real              :: ed      !vapor pressure [kPa]
    real              :: d       !vapor deficit [kPa]
    real              :: z0      !surface roughness length [m]
    real              :: wind0   !wind speed at 10-meter high
    real              :: ra      !aerodynamic roughness [s m-1]
    real              :: qle     !latent heat of vaporization [MJ/kg] 
    real              :: gama    ![kPa C-1]

    esat=xesatcf_a*exp(xesatcf_b*tair/(tair + xesatcf_c))

    del=(4098.0*esat)/((237.3+tair)**2.) ! [kPa C-1] (eq. 4.2.3)

    mspa=3.486*(pres/(tair+275.)) ! kg/m3 (eq. 4.2.4)

    ed=esat*qair

    d=esat-ed 

    z0=0.1*xz0m_sea

    wind0=max(0.001,wind*((log(10.0/z0))/log(2.0/z0))) !eq. 3.46  
    ra=(6.25/(wind0))*(log(10.0/z0))**2.

    qle=(2.501-0.002361*tair)

    gama=0.0016286*pres/qle

    ewat=1000*((del*qnet+mspa*cspa*d/ra)/(del+gama))/(qle*mesp) !kg m-2 s-1

  end subroutine HYMAP2_penman_monteith
  !=======================================================================
end module HYMAP2_evapMod

