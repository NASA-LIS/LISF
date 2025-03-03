!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

module pft_varcon

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Module of vegetation constants
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: pft_varcon.F90,v 1.6 2004/11/24 22:57:20 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2_varpar
  implicit none

  PRIVATE
! Vegetation type constants

  !  character(len=40) pftname(0:numpft) !PFT description
  character(len=27), public ::  pftname(0:numpft) !PFT description
  
  integer, public :: ncorn                  !value for corn
  integer, public :: nwheat                 !value for wheat
  integer, public :: noveg                  !value for not vegetated 
  integer, public :: ntree                  !value for last type of tree

  real(r8), public ::  dleaf(0:numpft)       !characteristic leaf dimension (m) 
  real(r8), public ::  c3psn(0:numpft)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(r8), public ::  vcmx25(0:numpft)      !max rate of carboxylation at 25C (umol CO2/m**2/s)
  real(r8), public ::  mp(0:numpft)          !slope of conductance-to-photosynthesis relationship
  real(r8), public ::  qe25(0:numpft)        !quantum efficiency at 25C (umol CO2 / umol photon)
  real(r8), public ::  xl(0:numpft)          !leaf/stem orientation index
  real(r8), public ::  rhol(0:numpft,numrad) !leaf reflectance: 1=vis, 2=nir 
  real(r8), public ::  rhos(0:numpft,numrad) !stem reflectance: 1=vis, 2=nir 
  real(r8), public ::  taul(0:numpft,numrad) !leaf transmittance: 1=vis, 2=nir 
  real(r8), public ::  taus(0:numpft,numrad) !stem transmittance: 1=vis, 2=nir 
  real(r8), public ::  z0mr(0:numpft)        !ratio of momentum roughness length to canopy top height (-)
  real(r8), public ::  displar(0:numpft)     !ratio of displacement height to canopy top height (-)
  real(r8), public ::  roota_par(0:numpft)   !CLM rooting distribution parameter [1/m]
  real(r8), public ::  rootb_par(0:numpft)   !CLM rooting distribution parameter [1/m]

  real(r8), public ::  sla(0:numpft)              !sp. leaf area [m2 leaf g-1 carbon]
  real(r8), public ::  pftpar(0:numpft,1:npftpar) !the rest for use with DGVM
  real(r8), public ::  lm_sapl(0:numpft)
  real(r8), public ::  sm_sapl(0:numpft)
  real(r8), public ::  hm_sapl(0:numpft)
  real(r8), public ::  rm_sapl(0:numpft)
  logical,  public ::   tree(0:numpft)
  logical,  public ::   summergreen(0:numpft)
  logical,  public ::   raingreen(0:numpft)

  real(r8), public, parameter :: reinickerp = 1.6 !parameter in allometric equation
  real(r8), public, parameter :: wooddens = 2.0e5 !wood density (gC/m3)
  real(r8), public, parameter :: latosa = 8.0e3   !ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b)
  real(r8), public, parameter :: allom1 = 100.0   !parameters in allometric
  real(r8), public, parameter :: allom2 =  40.0
  real(r8), public, parameter :: allom3 =   0.5

end module pft_varcon


















