!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module noah32_module
!BOP
!
! !MODULE: noah32_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the 
!  data structure containing the Noah3.2 1-d variables.
!  The variables specified in the data structure include: 
!
!  \begin{description}
!   \item[soiltype]
!    soil type (integer index)
!   \item[slopetype]
!    slope type (integer index)
!   \item[nroot] 
!    Number of root layers, depends on the vegetation type
!   \item[rsmin] 
!    Minimum canopy resistance (s/m)
!   \item[rgl] 
!    Solar radiation term of canopy resistance function
!   \item[hs] 
!    vapor pressure deficit term of canopy resistance function
!   \item[snup] 
!    threshold snow depth (in water equivalent m) that implies
!    100\% snow cover
!   \item[z0]
!    roughness length
!   \item[lai]
!    leaf area index
!   \item[shdfac]
!    greenness vegetation fraction
!   \item[shdmax]
!    maximum greenness vegetation fraction
!   \item[shdmin]
!    minimum greenness vegetation fraction
!   \item[smcmax]
!    maximum soil moisture content (porosity)
!   \item[slope]
!    slope coefficient for baseflow
!   \item[psisat]
!    saturated matric potential
!   \item[dksat]
!    saturated hydraulic conductivity
!   \item[bexp]
!    b parameter value
!   \item[quartz]
!    quartz fraction
!   \item[smcwlt]
!    wilting point (volumetric)
!   \item[smcref]
!    reference soil moisture (onset of soil moisture stress in
!    transpiration, volumetric)
!   \item[dwsat]
!    saturated soil diffusivity
!   \item[mxsnalb]
!    maximum albedo expected over deep snow
!   \item[tempbot]
!    bottom boundary temperature
!   \item[t1]
!    skin temperature (K)
!   \item[cmc]
!    canopy water content
!   \item[snowh]
!    actual snow depth (m)
!   \item[sneqv]
!    snow water equivalent (m)
!   \item[stc]
!    soil temperature for different layers
!   \item[smc]
!    soil moisture for different layers
!   \item[sh2o]
!    liquid-only soil moisture for different layers
!   \item[ch]
!    heat/moisture exchange coefficient
!   \item[cm]
!    momentum exchange coefficient
!   \item[vegt]
!    vegetation type of tile
!   \item[tair\_agl\_min]
!     minimum 2 meter AGL temperature
!   \item[rhmin]
!     minimum 2 meter relative humidity
!   \item[relsmc]
!    Volumetric relative soil moisture $(m^3/m^3)$
!    (smc - smcwlt) / (porosity - smcwlt)
!   \end{description}
!
! !REVISION HISTORY:
!  28 Apr 2002: K. Arsenault added Noah2.5 LSM code to LDAS. 
!  14 Nov 2002: Sujay Kumar Optimized version for LIS  
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
! 
!EOP
  implicit none

  PRIVATE
  type, public :: noah32dec     
     integer :: soiltype     
     integer :: slopetype
     real, allocatable :: lyrthk(:) !used only if distributed depths are defined
     integer :: nroot
     real    :: rsmin
     real    :: rgl
     integer :: vegt          
     real :: hs
     real :: snup
     real :: z0_old
     real :: z0brd
     real :: z0brd_old
     real :: lai
     real :: shdfac
     real :: shdmin
     real :: shdmax
     real :: slope
     real :: xice
     real :: embrd
     real :: z
     real :: q2sat
     real :: eta_kinematic
     real :: tauu
     real :: tauv
     real :: tau
     real :: soilm
     real :: smcmax
     real :: smcdry
     real :: psisat
     real :: dksat
     real :: bexp
     real :: quartz
     real :: smcwlt
     real :: smcref
     real :: dwsat
     real :: mxsnalb
     real :: tempbot
     real :: alb      ! background snow free albedo

!-------------------------------------------------------------------------
! Noah3.2-State Variables
!-------------------------------------------------------------------------
     real :: t1               
     real :: cmc              
     real :: snowh            
     real :: sneqv
     real,allocatable :: stc(:)   
     real,allocatable :: smc(:)   
     real,allocatable :: sh2o(:)  
     real :: ch               
     real :: cm    
     real :: z0
     real :: emiss
     real :: albedo
     real :: snotime1
     real :: q1
!yliu17
     real :: sndens  !snow density as fraction of water density

!-------------------------------------------------------------------------
! Noah3.2-Forcing Variables
!-------------------------------------------------------------------------
     real :: tair
     real :: qair
     real :: swdown
     real :: lwdown
     real :: uwind
     real :: vwind
     real :: psurf
     real :: rainf
     real :: rainf_c
     real :: snowf
     real :: fheight
     real :: qsfc
     real :: chs2
     real :: cqs2
     real :: q2
     real :: t2
     real :: th2

!-----------------------------------------------------------------------
! Noah3.2-Output variables
!-----------------------------------------------------------------------
     REAL :: swnet
     REAL :: lwnet
     REAL :: qle
     REAL :: qh
     REAL :: qg
     REAL :: beta
     REAL :: etanrg
     REAL :: snowfall
     REAL :: pcp
     REAL :: evap
     REAL :: qs
     REAL :: qsb
     REAL :: qsm
     REAL :: swe
     REAL :: sca
     REAL :: soilmoist(4)
     REAL :: soilwet
     REAL :: ecanop
     REAL :: tveg
     REAL :: esoil
     REAL :: canopint
     REAL :: rootmoist
     REAL :: runoff1
     REAL :: runoff2
     REAL :: tair_agl_min
     REAL :: rhmin

     real :: czil
     real :: refdk
     real :: refkdt
     real :: kdt
     real :: rsmax
     real :: topt
     real :: cfactr
     real :: frzx
     real :: frzk
     real :: fxexp
     real :: sbeta
     real :: csoil
     real :: cmcmax
     real :: zbot
     real :: salp
     real :: f1
     real :: lvcoef
     real :: sigma_flx

     real,allocatable :: relsmc(:)  !used to export relsmc to WRF
  end type noah32dec

end module noah32_module
