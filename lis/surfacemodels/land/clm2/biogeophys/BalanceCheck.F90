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

subroutine BalanceCheck (clm) 

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
! Water and energy balance check  
!
! Method:
! This subroutine accumulates the numerical truncation errors of the water
! and energy balance calculation. It is helpful to see the performance of 
! the process of integration.
!
! The error for energy balance: 
! error = abs(Net radiation - the change of internal energy - Sensible heat
!             - Latent heat) 
! The error should be less than 0.02 W m-2 in each time integration interval;
!
! The error for water balance:
! error = abs(precipitation - change of water storage - evaporation - runoff)
! The error should be less than 0.001 mm in  each time integration interval.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! 10 November 2000: Mariana Vertenstein
!
!-----------------------------------------------------------------------
! $Id: BalanceCheck.F90,v 1.7 2004/11/24 22:56:15 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
!<debug>
!  use clm2_varpar, only : nlevsoi
!</debug>
  use clm2_varcon, only : istsoil, tfrz
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm !CLM 1-D Module

!----Local Variables----------------------------------------------------
!  integer j                    ! do loop index
!  logical :: constop = .false. ! true => stop if energy balance err too great
  real(r8) :: errh2o, errseb,errlon!,errsol
!----End Variable List--------------------------------------------------

!
! Water balance 
!

  
  errh2o = clm%endwb - clm%begwb - &
           ( clm%forc_rain  + clm%forc_snow - clm%qflx_evap_tot - clm%qflx_surf &
           - clm%qflx_qrgwl - clm%qflx_drain ) * clm%dtime
!  if (clm%kpatch .eq. 56063) then
!  print*, clm%kpatch, errh2o, clm%begwb, clm%endwb, clm%forc_rain, clm%forc_snow
!  endif
  if (abs(errh2o) > .15) then
!     print*,'h',iam,clm%forc_rain,clm%forc_snow,clm%qflx_evap_tot,clm%qflx_surf
!     print*,'h2',iam,clm%qflx_qrgwl,clm%qflx_drain,clm%dtime 
     print*, clm%h2osoi_liq, clm%h2osoi_ice, clm%t_veg
     write(6,200)'water balance error',clm%nstep,clm%kpatch,errh2o
!     write(6,*)'clm model is stopping'
!     call endrun
  endif

!
! Solar radiation energy balance
!

!  errsol = clm%fsa + clm%fsr - (clm%forc_solad(1) + clm%forc_solad(2) &
!               + clm%forc_solai(1) + clm%forc_solai(2))

!  if (abs(errsol) > .10 ) then
!     write(6,100)'solar radiation balance error',clm%nstep,clm%kpatch,errsol
!     write(6,*)'clm model is stopping'
!     call endrun
!  endif

!
! Longwave radiation energy balance
!

  errlon = clm%eflx_lwrad_out - clm%eflx_lwrad_net - clm%forc_lwrad

  if (abs(errlon) > .10 ) then
     write(6,100)'longwave energy balance error',clm%nstep,clm%kpatch,errlon
!     write(6,*)'clm model is stopping'
!     call endrun
  endif

!
! Surface energy balance
!

  errseb = clm%sabv + clm%sabg  &
             + clm%forc_lwrad - clm%eflx_lwrad_out &
             - clm%eflx_sh_tot &
             - clm%eflx_lh_tot &
             - clm%eflx_soil_grnd 
                
  if (abs(errseb) > .10 ) then
     write(6,100)'surface flux energy balance error',clm%nstep,clm%kpatch,errseb
!     write(6,*)'clm model is stopping'
!     call endrun
  endif

!
! Accumulation of water and surface energy balance error
!

!  acc_errh2o = clm%acc_errh2o + errh2o
!  acc_errseb = clm%acc_errseb + errseb

100 format (1x,a14,' nstep =',i10,' point =',i6,' imbalance =',f8.2,' W m-2') 
200 format (1x,a14,' nstep =',i10,' point =',i6,' imbalance =',f8.2,' mm') 

end subroutine BalanceCheck












