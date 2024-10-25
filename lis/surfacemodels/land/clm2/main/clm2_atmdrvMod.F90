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
module clm2_atmdrvMod
!BOP
!
! !MODULE: clm2_atmdrvMod
!  This module defines the variables and routines to enable the 
!  interaction with an atmospheric component
! 
! !USES: 

!EOP
  implicit none

  PRIVATE

  public :: clm2_atmdrv

! logical variables for file manipuation

  logical :: open_data=.true.             !true => open data file (first tstep of the run or month)
  logical :: allocated_data=.false.       !true => allocate dynamic data

! atmospheric grid data

  integer  :: atmlon                      !number of atm longitudes
  integer  :: atmlat                      !number of atm latitudes
!  real(r8) :: edge_a(4)                   !N,E,S,W edges of atm grid
!  integer , allocatable :: numlon_a(:)    !number of lon points at each lat
!  real(r8), allocatable :: latixy_a(:,:)  !latitude of grid cell (degrees)
!  real(r8), allocatable :: longxy_a(:,:)  !longitude of grid cell (degrees)

! atmospheric forcing variables on atmospheric grid

!  real(r8), allocatable :: x(:,:,:)            !temp. array in which atm data is stored
!  real(r8), allocatable :: forc_txy_a (:,:)    !atm bottom level temperature (Kelvin)
!  real(r8), allocatable :: forc_uxy_a (:,:)    !atm bottom level zonal wind (m s-1)
!  real(r8), allocatable :: forc_vxy_a (:,:)    !atm bottom level meridional wind (m s-1)
!  real(r8), allocatable :: forc_qxy_a (:,:)    !atm bottom level specific humidity (kg kg-1)
!  real(r8), allocatable :: zgcmxy_a (:,:)      !atm bottom level height above surface (m)
!  real(r8), allocatable :: prcxy_a  (:,:)      !convective precipitation rate (mm H2O/s)
!  real(r8), allocatable :: prlxy_a  (:,:)      !large-scale precipitation rate (mm H2O/s)
!  real(r8), allocatable :: flwdsxy_a(:,:)      !downward longwave rad onto surface (W/m**2)
!  real(r8), allocatable :: forc_solsxy_a (:,:) !vis direct beam solar rad onto srf (W/m**2)
!  real(r8), allocatable :: forc_sollxy_a (:,:) !nir direct beam solar rad onto srf (W/m**2)
!  real(r8), allocatable :: forc_solsdxy_a(:,:) !vis diffuse solar rad onto srf (W/m**2)
!  real(r8), allocatable :: forc_solldxy_a(:,:) !nir diffuse solar rad onto srf(W/m**2)
!  real(r8), allocatable :: forc_pbotxy_a (:,:) !atm bottom level pressure (Pa)
!  real(r8), allocatable :: forc_psrfxy_a (:,:) !atm surface pressure (Pa)

! atmospheric forcing variables on land model grid

!  real(r8) :: forc_txy (lsmlon,lsmlat)         !atm bottom level temperature (Kelvin)
!  real(r8) :: forc_uxy (lsmlon,lsmlat)         !atm bottom level zonal wind (m s-1)
!  real(r8) :: forc_vxy (lsmlon,lsmlat)         !atm bottom level meridional wind (m s-1)
!  real(r8) :: forc_qxy (lsmlon,lsmlat)         !atm bottom level specific humidity (kg kg-1)
!  real(r8) :: zgcmxy (lsmlon,lsmlat)           !atm bottom level height above surface (m)
!  real(r8) :: prcxy  (lsmlon,lsmlat)           !convective precipitation rate (mm H2O/s)
!  real(r8) :: prlxy  (lsmlon,lsmlat)           !large-scale precipitation rate (mm H2O/s)
!  real(r8) :: flwdsxy(lsmlon,lsmlat)           !downward longwave rad onto surface (W/m**2)
!  real(r8) :: forc_solsxy (lsmlon,lsmlat)      !vis direct beam solar rad onto srf (W/m**2)
!  real(r8) :: forc_sollxy (lsmlon,lsmlat)      !nir direct beam solar rad onto srf (W/m**2)
!  real(r8) :: forc_solsdxy(lsmlon,lsmlat)      !vis diffuse solar rad onto srf (W/m**2)
!  real(r8) :: forc_solldxy(lsmlon,lsmlat)      !nir diffuse solar rad onto srf(W/m**2)
!  real(r8) :: forc_pbotxy (lsmlon,lsmlat)      !atm bottom level pressure (Pa)
!  real(r8) :: forc_psrfxy (lsmlon,lsmlat)      !atm surface pressure (Pa)

! atmosphere grid to land model surface grid mapping for each land grid cell:

!  integer, parameter :: mxovr =10          !maximum number of overlapping cells
!  integer :: novr_a2s(lsmlon,lsmlat)       !number    of overlapping atm cells
!  integer :: iovr_a2s(lsmlon,lsmlat,mxovr) !lon index of overlapping atm cells
!  integer :: jovr_a2s(lsmlon,lsmlat,mxovr) !lat index of overlapping atm cells
!  real(r8):: wovr_a2s(lsmlon,lsmlat,mxovr) !weight    of overlapping atm cells

! 1d temporarys

!  real(r8), allocatable, private :: forc_sols(:)   !vis direct beam solar rad onto srf (W/m**2)
!  real(r8), allocatable, private :: forc_soll(:)   !nir direct beam solar rad onto srf (W/m**2)
!  real(r8), allocatable, private :: forc_solsd(:)  !vis diffuse solar rad onto srf (W/m**2)
!  real(r8), allocatable, private :: forc_solld(:)  !nir diffuse solar rad onto srf(W/m**2)

! file netCDF id's

  integer :: ncid                !netCDF dataset id
  integer :: nvar                !number of variables in the data file
  integer :: nlon                !number of atm longitude points
  integer :: nlat                !number of atm latitude points
  integer :: ntim                !number of atm time slices per data file
  character(len=8) :: varnam(99) !variable names of atm. fields

  SAVE

!=======================================================================
contains
!=======================================================================

!BOP
! 
! !ROUTINE: clm2_atmdrv
! \label{clm2_atmdrv}
!
! !INTERFACE: 
  subroutine clm2_atmdrv(n)
! !USES: 
    use ESMF
    use clm2_lsmMod
    use clm2_varcon  , only : rair, cpair, po2, pco2, tcrit, tfrz    
    use LIS_precisionMod
    use LIS_coreMod, only : LIS_rc, LIS_surface
    use LIS_FORC_AttributesMod 
    use LIS_metforcingMod, only : LIS_FORC_State
    use LIS_logMod,         only : LIS_verify, LIS_endrun
    
    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n 
!
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the CLM2
!  model tiles. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP
    real(r8) :: solar(LIS_rc%npatch(n,LIS_rc%lsm_index))
    real(r8) :: prcp(LIS_rc%npatch(n,LIS_rc%lsm_index))
    real(r8) :: forc_vp
    integer            :: tid
    integer            :: t,status
    type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
    type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
    real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
    real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
    
    
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),&
                       tmpField, rc=status)
    call LIS_verify(status,'clm2_atmdrv: error getting Tair')
    
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),&
                       q2Field, rc=status)
    call LIS_verify(status,'clm2_atmdrv: error getting Qair')
    
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),&
                       swdField, rc=status)
    call LIS_verify(status,'clm2_atmdrv: error getting SWdown')
    
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),&
                       lwdField, rc=status)
    call LIS_verify(status,'clm2_atmdrv: error getting LWdown')
    
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),&
                       uField, rc=status)
    call LIS_verify(status,'clm2_atmdrv: error getting Wind_E')
    
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),&
                       vField, rc=status)
    call LIS_verify(status,'clm2_atmdrv: error getting Wind_N')
    
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),&
                       psurfField, rc=status)
    call LIS_verify(status, 'clm2_atmdrv: error getting PSurf')
    
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),&
                       pcpField, rc=status)
    call LIS_verify(status,'clm2_atmdrv: error getting Rainf')
    
    if(LIS_FORC_CRainf%selectOpt.eq.1) then 
       call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CRainf%varname(1)),&
                          cpcpField, rc=status)
       call LIS_verify(status,'clm2_atmdrv: error getting CRainf')
    endif
    
    if(LIS_FORC_Snowf%selectOpt.eq.1) then 
       call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),&
                          snowfField, rc=status)
       call LIS_verify(status,'clm2_atmdrv: error getting Snowf')
    endif
    
    call ESMF_FieldGet(tmpField, localDE=0, farrayPtr=tmp,rc=status)
    call LIS_verify(status, 'clm2_atmdrv: error retrieving tmp')
    
    call ESMF_FieldGet(q2Field, localDE=0, farrayPtr=q2,rc=status)
    call LIS_verify(status, 'clm2_atmdrv: error retrieving q2')
    
    call ESMF_FieldGet(swdField, localDE=0, farrayPtr=swd,rc=status)
    call LIS_verify(status, 'clm2_atmdrv: error retrieving swd')
    
    call ESMF_FieldGet(lwdField, localDE=0, farrayPtr=lwd,rc=status)
    call LIS_verify(status, 'clm2_atmdrv: error retrieving lwd')
    
    call ESMF_FieldGet(uField, localDE=0, farrayPtr=uwind,rc=status)
    call LIS_verify(status, 'clm2_atmdrv: error retrieving uwind')
    
    call ESMF_FieldGet(vField, localDE=0, farrayPtr=vwind,rc=status)
    call LIS_verify(status, 'clm2_atmdrv: error retrieving vwind')
    
    call ESMF_FieldGet(psurfField, localDE=0, farrayPtr=psurf,rc=status)
    call LIS_verify(status, 'clm2_atmdrv: error retrieving psurf')
    
    call ESMF_FieldGet(pcpField, localDE=0, farrayPtr=pcp,rc=status)
    call LIS_verify(status, 'clm2_atmdrv: error retrieving pcp')
    
    if(LIS_FORC_CRainf%selectOpt.eq.1) then 
       call ESMF_FieldGet(cpcpField, localDE=0, farrayPtr=cpcp,rc=status)
       call LIS_verify(status, 'clm2_atmdrv: error retrieving cpcp')
    endif
    
    if(LIS_FORC_Snowf%selectOpt.eq.1) then 
       call ESMF_FieldGet(snowfField, localDE=0, farrayPtr=snowf,rc=status)
       call LIS_verify(status, 'clm2_atmdrv: error retrieving snowf')
    endif

    clm2_struc(n)%forc_count =  clm2_struc(n)%forc_count + 1
    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
       tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
#if(defined COUPLED) 
       clm2_struc(n)%clm(t)%forc_t        = tmp(tid)
       clm2_struc(n)%clm(t)%forc_q        = q2(tid)
       solar(t)                           = swd(tid)
       clm2_struc(n)%clm(t)%forc_solad(1) = solar(t)*35.0/100.0
       clm2_struc(n)%clm(t)%forc_solad(2) = solar(t)*35.0/100.0
       clm2_struc(n)%clm(t)%forc_solai(1) = solar(t)*15.0/100.0
       clm2_struc(n)%clm(t)%forc_solai(2) = solar(t)*15.0/100.0
       clm2_struc(n)%clm(t)%forc_lwrad    = lwd(tid)
       clm2_struc(n)%clm(t)%forc_u        = uwind(tid)
       clm2_struc(n)%clm(t)%forc_v        = vwind(tid)
       clm2_struc(n)%clm(t)%forc_pbot     = psurf(tid)
       prcp(t)                            = pcp(tid)
       print*,'need to translate a few more variables.. '
       call LIS_endrun()
#if JIM_FIX_THESE
       clm2_struc(n)%clm(t)%forc_hgt      = 0.5*forcing(9id)
       clm2_struc(n)%clm(t)%forc_ch       = forcing(6)
       clm2_struc(n)%clm(t)%forc_q2sat    = forcing(10)
       clm2_struc(n)%clm(t)%forc_cosz     = forcing(13)
   
       clm2_struc(n)%clm(t)%forc_hgt_u    = clm2_struc(n)%clm(t)%forc_hgt+&
           clm2_struc(n)%clm(t)%displa+clm2_struc(n)%clm(t)%z0m
       clm2_struc(n)%clm(t)%forc_hgt_t    = clm2_struc(n)%clm(t)%forc_hgt+&
           clm2_struc(n)%clm(t)%displa+clm2_struc(n)%clm(t)%z0m
       clm2_struc(n)%clm(t)%forc_hgt_q    = clm2_struc(n)%clm(t)%forc_hgt+&
           clm2_struc(n)%clm(t)%displa+clm2_struc(n)%clm(t)%z0m
   
!       clm2_struc(n)%clm(k)%forc_hgt   =10.0+clm2_struc(n)%clm(k)%displa+&
!            clm2_struc(n)%clm(k)%z0m
!       clm2_struc(n)%clm(k)%forc_hgt_t  =10.0+clm2_struc(n)%clm(k)%displa+&
!            clm2_struc(n)%clm(k)%z0m
!       clm2_struc(n)%clm(k)%forc_hgt_u  =10.0+clm2_struc(n)%clm(k)%displa+&
!            clm2_struc(n)%clm(k)%z0m
!       clm2_struc(n)%clm(k)%forc_hgt_q  =10.0+clm2_struc(n)%clm(k)%displa+&
!            clm2_struc(n)%clm(k)%z0m
#endif
#else
       clm2_struc(n)%clm(t)%forc_t        = clm2_struc(n)%clm(t)%forc_t + tmp(tid)
       clm2_struc(n)%clm(t)%forc_q        = clm2_struc(n)%clm(t)%forc_q + q2(tid)
       solar(t)                          = swd(tid)
       clm2_struc(n)%clm(t)%forc_solad(1) = &
            clm2_struc(n)%clm(t)%forc_solad(1) + solar(t)*35.0/100.0
       clm2_struc(n)%clm(t)%forc_solad(2) =  &
            clm2_struc(n)%clm(t)%forc_solad(2) + solar(t)*35.0/100.0
       clm2_struc(n)%clm(t)%forc_solai(1) = &
             clm2_struc(n)%clm(t)%forc_solai(1) + solar(t)*15.0/100.0
       clm2_struc(n)%clm(t)%forc_solai(2) = &
             clm2_struc(n)%clm(t)%forc_solai(2)+ solar(t)*15.0/100.0
       clm2_struc(n)%clm(t)%forc_lwrad    = &
             clm2_struc(n)%clm(t)%forc_lwrad + lwd(tid)
       clm2_struc(n)%clm(t)%forc_u        = &
             clm2_struc(n)%clm(t)%forc_u + uwind(tid)
       clm2_struc(n)%clm(t)%forc_v        = &
            clm2_struc(n)%clm(t)%forc_v + vwind(tid)
       clm2_struc(n)%clm(t)%forc_pbot     = &
            clm2_struc(n)%clm(t)%forc_pbot + psurf(tid)
       prcp(t)                           = pcp(tid)
#endif
    
       if (prcp(t) > 0.) then
          if (clm2_struc(n)%clm(t)%forc_t > (tfrz + tcrit)) then
             clm2_struc(n)%clm(t)%itypprc   = 1
             clm2_struc(n)%clm(t)%forc_rain = &
                  clm2_struc(n)%clm(t)%forc_rain + prcp(t)
             clm2_struc(n)%clm(t)%forc_snow = &
                  clm2_struc(n)%clm(t)%forc_snow + 0.
          else
             clm2_struc(n)%clm(t)%itypprc   = 2
             clm2_struc(n)%clm(t)%forc_rain = &
                  clm2_struc(n)%clm(t)%forc_rain + 0.
             clm2_struc(n)%clm(t)%forc_snow = &
                  clm2_struc(n)%clm(t)%forc_snow + prcp(t)
          endif
       else
          clm2_struc(n)%clm(t)%itypprc   = 0
          clm2_struc(n)%clm(t)%forc_rain = &
               clm2_struc(n)%clm(t)%forc_rain + 0.
          clm2_struc(n)%clm(t)%forc_snow = & 
               clm2_struc(n)%clm(t)%forc_snow + 0
       endif
! Derive new fields (potential temperature, vapor pressure, 
! air density, CO2, and O2) and copy solar radiations

!=== Cuurently potential temperature set to 2 m temperature since 
!=== we only get surface pressure in our forcing and elevation differences 
!=== are accounted for in for.f

!===LDAS modification: Slight change to be consistent with our forcing dataset
!          clm2_struc(n)%clm(t)%forc_th  = clm2_struc(n)%clm(t)%forc_t * (clm2_struc(n)%clm(t)%forc_psrf/clm2_struc(n)%clm(t)%forc_pbot)**(rair/cpair)
!       clm2_struc(n)%clm(t)%forc_th  = tmp(t) * &
!            (clm2_struc(n)%clm(t)%forc_pbot/clm2_struc(n)%clm(t)%forc_pbot)**(rair/cpair)
       clm2_struc(n)%clm(t)%forc_th  =  tmp(tid) * &
            (psurf(tid)/psurf(tid))**(rair/cpair)
       forc_vp  = q2(tid)*psurf(tid) / &
            (0.622+0.378*q2(tid))
       clm2_struc(n)%clm(t)%forc_rho = &
            (psurf(tid)-0.378*forc_vp) / &
            (rair*tmp(tid))

    enddo
    
  end subroutine clm2_atmdrv

end module clm2_atmdrvMod

