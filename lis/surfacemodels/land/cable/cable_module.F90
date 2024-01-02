!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
      module cable_module
!BOP
!
! !MODULE: cable_module
! \label{cable_module}
!
! !REVISION HISTORY:
!  21 Jul 2004: Sujay Kumar, Initial Specification
!  23 Oct 2007: Kristi Arsenault, Added code for LISv5.0
!  25 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!  14 Dec 2011: Claire Carouge (ccc), CABLE LSM improvements
!  23 May 2013: David Mocko, latest CABLE v1.4b version for LIS6.2
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the CABLE 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}
!   \item[forcing]
!     Array of meteorological forcing
!   \item{tair}
!     2m air temperature forcing
!   \item{qair}
!     2m specific humidity forcing
!   \item{swdown}
!     downward shortwave forcing
!   \item{lwdown}
!     downward longwave forcing
!   \item{uwind}
!     u-wind component forcing
!   \item{vwind}
!     v-wind component forcing
!   \item{psurf}
!     surface pressure forcing
!   \item{rainf}
!     total rainfall forcing
!   \item{rainf\_c}
!     convective rainfall forcing
!   \item{snowf}
!     total snowfall forcing
!   \item[lai]
!    leaf area index
!   \item[vegtype]
!    vegetation type (integer index)
!   \item[soiltype]
!    soil type (integer index)
!  \end{description}
!
! !USES:
      use cable_dimensions,   only : ms,msn,ncp,ncs
!EOP
      implicit none

      type cabledec

!-------------------------------------------------------------------------
! CABLE - State Variables
!-------------------------------------------------------------------------
      real :: cansto
      real :: rtsoil
      real :: ssdnn
      real :: snowd
      real :: osnowd
      real :: snage
! Need to check if writing/reading an integer in LIS restart file is OK - dmm
      integer :: isflag
      real :: wbice(ms)
      real :: tggsn(msn)
      real :: ssdn(msn)
      real :: smass(msn)
      real :: albsoilsn(msn)
      real :: wb(ms)
      real :: tgg(ms)
      real :: cplant(ncp)
      real :: csoil(ncs)
      real :: tk_old
!ccc
      real :: pwb_min
      real :: wbtot
      real :: gammzz(ms)

!-------------------------------------------------------------------------
! CABLE - Parameter Variables
!-------------------------------------------------------------------------
      real :: lai

      integer :: vegtype
      integer :: meth
      real :: canst1
      real :: dleaf
      real :: vcmax
      real :: ejmax
      real :: hc
      real :: xfang
      real :: vbeta
      real :: rp20
      real :: rpcoef
      real :: rs20
      real :: shelrb
      real :: wai
      real :: vegcf
      real :: extkn
      real :: tminvj
      real :: tmaxvj
      real :: frac4
      real :: cplant_init(ncp)
      real :: ratecp(ncp)
      real :: csoil_init(ncs)
      real :: ratecs(ncs)
      real :: froot(ms)

      integer :: soiltype
      real :: albsoil
      real :: silt
      real :: clay
      real :: sand
      real :: swilt
      real :: sfc
      real :: ssat
      real :: bch
      real :: hyds
      real :: sucs
      real :: rhosoil
      real :: css
      real :: zse(ms)
      real :: zshh(ms+1)

!-------------------------------------------------------------------------
! CABLE - Forcing Variables
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
      real :: co2

! Added for coupling to WRF. (ccc, 17/10/11)
      real :: za
      real :: forc_ch
      real :: forc_chs2
      real :: forc_cqs2
      real :: q2sat
      real :: coszen

!-------------------------------------------------------------------------
! CABLE - Output Variables
!-------------------------------------------------------------------------
      real :: swnet
      real :: lwnet
      real :: rnet
      real :: qle
      real :: qh
      real :: qg
      real :: qs
      real :: qsb
      real :: evap
      real :: vegt
      real :: radt
      real :: albedo
      real :: sdepth
      real :: ecanop
      real :: tveg
      real :: esoil
      real :: autoresp
      real :: heteroresp
      real :: leafresp
      real :: npp
      real :: gpp
      real :: nee
      real :: potev

! Added for coupling to WRF
      real :: z0m
      real :: smelt
      real :: trad

      end type cabledec

      end module cable_module
