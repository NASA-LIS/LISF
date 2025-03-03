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

subroutine SnowAlbedo (clm, coszen, nband, ind, alb)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                        M  available land surface process model.
!  M --COMMON LAND MODEL--  C
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List:
!
!-----------------------------------------------------------------------
! Purpose: 
! Determine snow albedos
! 
! Method: 
! 
! Author:
! Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: SnowAlbedo.F90,v 1.6 2004/11/24 22:56:36 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varpar,   only : numrad
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm       !CLM 1-D Module

  real(r8), intent(in) :: coszen ! cosine solar zenith angle for next time step
  integer, intent(in) :: nband   ! number of solar radiation waveband classes
  integer, intent(in) :: ind     ! 0=direct beam, 1=diffuse radiation

  real(r8), intent(out):: alb(numrad) !snow albedo by waveband

!----Local Variables----------------------------------------------------

  integer  :: ib          !waveband class

! CLM variables and constants for snow albedo calculation

  real(r8) :: snal0 = 0.95 !vis albedo of new snow for sza<60
  real(r8) :: snal1 = 0.65 !nir albedo of new snow for sza<60
  real(r8) :: conn  = 0.5  !constant for visible snow alb calculation [-]
  real(r8) :: cons  = 0.2  !constant (=0.2) for nir snow albedo calculation [-]
  real(r8) :: sl    = 2.0  !factor that helps control alb zenith dependence [-]
  real(r8) :: age          !factor to reduce visible snow alb due to snow age [-]
  real(r8) :: albs         !temporary vis snow albedo
  real(r8) :: albl         !temporary nir snow albedo
  real(r8) :: cff          !snow alb correction factor for zenith angle > 60 [-]
  real(r8) :: czf          !solar zenith correction for new snow albedo [-]

!----End Variable List--------------------------------------------------

!
! Zero albedos
!

  do ib = 1, nband
     alb(ib) = 0._r4
  end do

!
! CLM Albedo for snow cover.
! Snow albedo depends on snow-age, zenith angle, and thickness of snow,
! age gives reduction of visible radiation
!

!
! Correction for snow age

!  print*,'sn..',iam, clm%snowage
  age = 1.-1./(1.+clm%snowage)
  albs = snal0*(1.-cons*age)
  albl = snal1*(1.-conn*age)

  if (ind == 0) then

!
! czf corrects albedo of new snow for solar zenith
!

    cff    = ((1.+1./sl)/(1.+max(0.001_r4,coszen)*2.*sl )- 1./sl)
    cff    = max(cff,0._r4)
    czf    = 0.4*cff*(1.-albs)
    albs = albs+czf
    czf    = 0.4*cff*(1.-albl)
    albl = albl+czf

  endif

  alb(1) = albs
  alb(2) = albl
  return
end subroutine SnowAlbedo






