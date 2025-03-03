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

subroutine SnowAge (clm)

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
! Updates snow age
!
! Method:
! Based on BATS code.
!
! Author:
! Original Code:  Robert Dickinson
! 15 September 1999: Yongjiu Dai; Integration of code into CLM
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: SnowAge.F90,v 1.6 2004/11/24 22:56:35 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
!  use spmdMod
  use clm2_varcon, only : tfrz
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm  !CLM 1-D Module

!----Local Variables----------------------------------------------------

  real(r8) age1 ! snow aging factor due to crystal growth [-]
  real(r8) age2 ! snow aging factor due to surface growth [-]
  real(r8) age3 ! snow aging factor due to accum of other particles [-]
  real(r8) arg  ! temporary variable used in snow age calculation [-]
  real(r8) arg2 ! temporary variable used in snow age calculation [-]
  real(r8) dela ! temporary variable used in snow age calculation [-]
  real(r8) dels ! temporary variable used in snow age calculation [-]
  real(r8) sge  ! temporary variable used in snow age calculation [-]

!----End Variable List--------------------------------------------------

  if (clm%h2osno <= 0.) then

     clm%snowage = 0.

  else if (clm%h2osno > 800.) then   ! Over Antarctica

     clm%snowage = 0.

  else                               ! Away from Antarctica 

     age3  = 0.3
     arg   = 5.e3*(1./tfrz-1./clm%t_grnd)
     arg2  = min(0._r4,10.*arg)
     age2  = exp(arg2)
     age1  = exp(arg)
     dela  = 1.e-6*clm%dtime*(age1+age2+age3)
     dels  = 0.1*max(0.0_r4, clm%h2osno-clm%h2osno_old)
!	print*,'sn',iam, clm%t_grnd
     sge   = (clm%snowage+dela)*(1.0-dels)
     clm%snowage   = max(0.0_r4,sge)

  endif

end subroutine SnowAge
