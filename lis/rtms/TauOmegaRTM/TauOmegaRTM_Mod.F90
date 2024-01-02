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
module TauOmegaRTM_Mod
!BOP
!
! !MODULE: TauOmegaRTM_Mod
!
! !DESCRIPTION:
!    This module provides the routines to control the execution of 
!    the tau-omega model 
!
! !REVISION HISTORY:
! 28 Aug 2012: Sujay Kumar, initial specification based on 
!              the code from Wade Crow
!
! !USES:        


#if (defined RTMS)

  use ESMF
  use LIS_coreMod
  use LIS_RTMMod
  use LIS_logMod

  implicit none
 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: TauOmegaRTM_initialize
  public :: TauOmegaRTM_f2t
  public :: TauOmegaRTM_run
  public :: TauOmegaRTM_output
  public :: TauOmegaRTM_geometry 
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: tauomega_struc
!EOP
  type, public ::  tauomega_type_dec 
     real, allocatable :: h(:)
     real, allocatable :: b(:)
     real, allocatable :: w(:)
     real, allocatable :: vwc(:)
     real, allocatable :: Tb(:)
  end type tauomega_type_dec

  type(tauomega_type_dec), allocatable :: tauomega_struc(:) 

  SAVE

contains
!BOP
! 
! !ROUTINE: TauOmegaRTM_initialize
! \label{TauOmegaRTM_initialize}
! 
! !INTERFACE:
  subroutine TauOmegaRTM_initialize()
! !USES:

! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for noah-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for noah from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readTauOmegaRTMcrd](\ref{readTauOmegaRTMcrd}) \newline
!    reads the runtime options for TauOmegaRTM EMonly
!  \end{description}
!EOP
    implicit none
    
    integer :: n,t

    real    :: h_val(13), w_val(13), vwc_val(13),b_val(13)
    data h_val /1.2,1.3,1.2,1, 1.3, 0.7, 0.7, 0.7, 0.5,0.1,0.5,0.7,0.1/
    data w_val /0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05/
    data vwc_val /1,1,1,1,1,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0/
    data b_val /0.33,0.33,0.33,0.33,0.33,0.3,0.3,0.3,0.2,0.2,0.15,0.15,0/

    allocate(tauomega_struc(LIS_rc%nnest))

    do n=1,LIS_rc%nnest
       
       allocate(tauomega_struc(n)%h(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(tauomega_struc(n)%b(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(tauomega_struc(n)%w(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(tauomega_struc(n)%vwc(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       allocate(tauomega_struc(n)%Tb(LIS_rc%npatch(n,LIS_rc%lsm_index)))

       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          tauomega_struc(n)%h(t) = & 
               h_val(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt)
          tauomega_struc(n)%b(t) = & 
               b_val(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt)
          tauomega_struc(n)%vwc(t) = & 
               vwc_val(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt)
          tauomega_struc(n)%w(t) = & 
               w_val(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt)
       enddo
       
       call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Content")
       call add_sfc_fields(n,LIS_sfcState(n), "Soil Temperature")
    enddo

  end subroutine TauOmegaRTM_initialize 

  subroutine add_sfc_fields(n, sfcState,varname)

    implicit none 

    integer            :: n 
    type(ESMF_State)   :: sfcState
    character(len=*)   :: varname

    type(ESMF_Field)     :: varField
    type(ESMF_ArraySpec) :: arrspec
    integer              :: status
    real :: sum
    call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    varField = ESMF_FieldCreate(arrayspec=arrSpec, & 
         grid=LIS_vecTile(n), name=trim(varname), &
         rc=status)
    call LIS_verify(status, 'Error in field_create of '//trim(varname))
    
    call ESMF_StateAdd(sfcState, (/varField/), rc=status)
    call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

  end subroutine add_sfc_fields


  subroutine TauOmegaRTM_f2t(n)

    implicit none

    integer, intent(in)    :: n 

  end subroutine TauOmegaRTM_f2t


  subroutine TauOmegaRTM_geometry(n)
    implicit none
    integer, intent(in)    :: n

  end subroutine TauOmegaRTM_geometry 

  subroutine TauOmegaRTM_run(n)
! !USES: 
    
    implicit none

    integer, intent(in) :: n 

    real, parameter     :: freqm=1.41
    integer             :: t
    integer             :: status
    real, pointer       :: sm(:), ts(:)
    real                :: tveg(LIS_rc%ntiles(n))
    real                :: teff(LIS_rc%ntiles(n))
    complex             :: er_r,arg,rsh
    real                :: surrh, tau, exptau,a
    real                :: theta, stheta, ctheta

#if 0 
    theta = 0.6981
    ctheta = cos(theta)
    stheta = sin(theta)

!   map surface properties to SFC    
    call getsfcvar(LIS_sfcState(n), "Soil Moisture Content",&
         sm)
    call getsfcvar(LIS_sfcState(n), "Soil Temperature", &
         ts)
    tveg = ts

!---------------------------------------------
! Tile loop 
!--------------------------------------------
    do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
       sm(t) = sm(t)*0.01d0  !convert to g/cm3?
       teff(t) = ts(t)*0.5 + tveg(t)*0.5
       
       if (LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt.ne.14 &
            .and. sm(t).gt.0.0) then
!	  call hallik(sm(t),sand(t),clay(t),er_r)
          call dobson(sm(t),&
               LIS_surface(n,LIS_rc%lsm_index)%tile(t)%sand,&
               LIS_surface(n,LIS_rc%lsm_index)%tile(t)%clay,er_r,teff(t))

          arg = csqrt(er_r - stheta**2)
          rsh = cabs((ctheta - arg)/(ctheta + arg))
          rsh = rsh * rsh
          surrh = real(rsh) * exp (-tauomega_struc(n)%h(t))
          tau = tauomega_struc(n)%b(t)*tauomega_struc(n)%vwc(t)
          exptau = exp(-tau/ctheta)
          a = tveg(t)*(1.-tauomega_struc(n)%w(t))*(1.-exptau)
          tauomega_struc(n)%Tb(t)=teff(t)*(1.0-surrh)*exptau+a*(1.+surrh*exptau)
! lc=14 means an open water surface	  
       else if(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt.eq.14&
            .and.ts(t).gt.0.0) then
          call reflct(freqm,0.0,teff(t),theta,watrrh)
          tauomega_struc(n)%Tb(t) = tveg(t)*(1.-watrrh)
!	  write(*,*)watrrh,tveg(t),tauomega_struc(n)%Tb(t),tauomega_struc(n)%Tbb(t)
       else
          tauomega_struc(n)%Tb(t)=LIS_rc%udef
       endif       
    enddo
#endif
  end subroutine TauOmegaRTM_run

#if 0 
  subroutine hallik(vsm,sf,cf,eps)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!              eps     Dielectric Constant (Real Part)
!              vsm      Soil moisture content (g cm-3)
!              sf      percent of sand content in soil
!              cf      percent of clay content in soil
!  Note:  Soil temperature must be higher than 0 C
!
    real      :: vsm, sf,cf
    complex eps
    complex aa0,aa1,aa2,bb0,bb1,bb2,cc0,cc1,cc2
    
    aa0 = cmplx(2.862,0.356)
    aa1 = cmplx(-0.012,-0.003)
    aa2 = cmplx(0.001,-0.008)
    
    bb0 = cmplx(3.803,5.507)
    bb1 = cmplx(0.462,0.044)
    bb2 = cmplx(-0.341,-0.002)
    
    cc0 = cmplx(119.006,17.753)
    cc1 = cmplx(-0.500,-0.313)
    cc2 = cmplx(0.633,0.206)

    eps = (aa0+aa1*sf+aa2*cf)+(bb0+bb1*sf+bb2*cf)*vsm & 
         + (cc0+cc1*sf+cc2*cf)*vsm*vsm
    
    return
  end subroutine hallik
#endif	

#if 0 
	
  subroutine dobson(vsm,sf,cf,eps,ts)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!              eps     Dielectric Constant (Real Part)
!              vsm      Soil moisture content (g cm-3)
!              sf      percent of sand content in soil
!              cf      percent of clay content in soil
!  Note:  Soil temperature must be higher than 0 C
!
    real vsm, sf, cf,ts
    complex :: eps
    real    :: freq
    real alpha,epso,epswinf,epss,rhos,por,fv,rhob
    real betar,omtau,epswo,fac,epsfwr,tsC
! Set physical constants and bounds
    alpha   = 0.65  
    epso    = 8.854e-12 
    epswinf = 4.9   
    freq = 1.41    
	
! Compute dielectric constant of soil solid particles
    epss = 4.70
    rhos = (sqrt(epss + 0.062) - 1.01)/0.44
    por = 0.505 - 0.142*sf*0.01 - 0.037*cf*0.01
    fv = 1 - por
    rhob = fv*rhos

! Compute optimized coefficient values
    betar = 1.27480 - 0.519*sf*0.01 - 0.152*cf*0.01
    tsC = ts - 273.15
	
	
! Compute dielectric constant of free water (expressions for
! relaxation time and static dielectric constant obtained from
! Ulaby et al., Vol. 3)
    omtau  =  freq*(1.1109e-1 - 3.824e-3*tsC + 6.938e-5*tsC**2 - 5.096e-7*tsC**3)
    epswo  =  88.045 - 0.4147*tsC + 6.295e-4*tsC**2 + 1.075e-5*tsC**3
    fac    = (epswo - epswinf) / (1.0 + omtau**2)
    epsfwr =  epswinf + fac
    eps   = ( 1.0 + (epss**alpha - 1.0)*rhob/rhos + (vsm**betar) * ( (epsfwr**alpha) - 1.0) )**(1.0/alpha)
    
    return
  end subroutine dobson
	
	
	
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!REFLCT Computes the dielectric constant of water and reflectivity of a
!  smooth water surface as a function of frequency, water temperature,
!  salinity, and incidence angle (for reflectivity)
!
!  Ref: Klein and Swift (1977): An improved model for the dielectric constant
!    of sea water at microwave frequencies, IEEE Trans. Ant. Propagat., AP-25,
!    104-111.
!
!  Dielectric model assumed valid for salinity in range 4 to 35 psu, and
!    temperature in range 0 to 40 C.
! 
!  Variables:
!    freq       Frequency (GHz)
!    s          Water salinity (psu)
!    tk         Water temperature (K) - need to be in C
!    theta      Incidence angle (deg)
!    dr         Real part of dielectric constant
!    di         Imaginary part of dielectric constant
!    rh         Reflectivity (h-pol)
!

  subroutine reflct(freq, s, tk, theta, rh)
    real :: freq, s, tk, theta, rh
    real :: t
    complex eps

! Constants:
    t = tk - 273.15
    einf = 4.9
    e0 = 8.854e-12
    dt = 25.-t
    pi = 3.1415926535897
    
! Equations for Sigma:
    beta = 2.033e-2+(1.266e-4*dt)+(2.464e-6*dt**2)-  & 
         s*(1.849e-5-2.551e-7*dt+2.551e-8*dt**2)
    sig25 = s*(0.182521-1.46192e-3*s+2.09324e-5*s**2- & 
         1.28205e-7*s**3)
    sigma = sig25*exp(-dt*beta)

! Equations for Epsilon:
    et = 87.134 - 1.949e-1*t - 1.276e-2*t**2 + 2.491e-4*t**3
    ast = 1 + 1.613e-5*s*t - 3.656e-3*s + 3.21e-5*s**2 -  & 
         4.232e-7*s**3
    ets = et*ast

! Equations for Tau:
    tt = 1.768e-11 - 6.086e-13*t + 1.104e-14*t**2 - 8.111e-17*t**3
    bst = 1. + 2.282e-5*s*t - 7.638e-4*s - 7.76e-6*s**2 + & 
         1.105e-8*s**3
    tts = tt*bst

! Calculate permittivity:
    om = 2.*pi*freq*1.e9
    dr = einf + (ets-einf)/(1.+(om*tts)**2)
    di = (om*tts) * (ets-einf)/(1.+(om*tts)**2) + sigma/(om*e0)
    eps = cmplx(dr,-di)
       
! Compute reflectivities:
    sinth = sin(theta)
    costh = cos(theta)
    epss = sqrt(eps - sinth**2)
    rh = abs((costh-epss)/(costh+epss))**2
    return
  end subroutine reflct
#endif
  subroutine getsfcvar(sfcState, varname, var)
! !USES: 
    
    implicit none
    
    type(ESMF_State)      :: sfcState
    type(ESMF_Field)      :: varField
    character(len=*)      :: varname
    real, pointer         :: var(:)
    integer               :: status

    call ESMF_StateGet(sfcState, trim(varname), varField, rc=status)
    call LIS_verify(status, 'Error in StateGet: CMEM3_handlerMod '//trim(varname))
    call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
    call LIS_verify(status, 'Error in FieldGet: CMEM3_handlerMod '//trim(varname))

  end subroutine getsfcvar

!!!!BOP
!!!! !ROUTINE: TauOmegaRTM_output
!!!! \label{TauOmegaRTM_output}
!!!!
!!!! !INTERFACE: 
  subroutine TauOmegaRTM_output(n)
    integer, intent(in) :: n 
end subroutine TauOmegaRTM_output
#endif
end module TauOmegaRTM_Mod



