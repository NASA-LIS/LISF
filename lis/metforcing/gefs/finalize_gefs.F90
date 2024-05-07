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
!BOP
! !MODULE: finalize_gefs
!  \label{finalize_gefs}
!
! !REVISION HISTORY: 
! 7 Mar 2013: Sujay Kumar, initial specification
! 1 Jul 2019: K. Arsenault, expand support for GEFS forecasts
! 
! !INTERFACE:
subroutine finalize_gefs(findex)
! !USES:
  use LIS_coreMod,     only : LIS_rc
  use gefs_forcingMod, only : gefs_struc, gefs_initialized
!
! !DESCRIPTION:
!  Routine to finalize GEFS forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex
  
  if ( gefs_initialized ) then 
     do n=1,LIS_rc%nnest
       if( LIS_rc%met_interp(findex).eq."bilinear" ) then 
          deallocate(gefs_struc(n)%n111)
          deallocate(gefs_struc(n)%n121)
          deallocate(gefs_struc(n)%n211)
          deallocate(gefs_struc(n)%n221)
          deallocate(gefs_struc(n)%w111)
          deallocate(gefs_struc(n)%w121)
          deallocate(gefs_struc(n)%w211)
          deallocate(gefs_struc(n)%w221)

       elseif( LIS_rc%met_interp(findex).eq."budget-bilinear" ) then
          deallocate(gefs_struc(n)%n111)
          deallocate(gefs_struc(n)%n121)
          deallocate(gefs_struc(n)%n211)
          deallocate(gefs_struc(n)%n221)
          deallocate(gefs_struc(n)%w111)
          deallocate(gefs_struc(n)%w121)
          deallocate(gefs_struc(n)%w211)
          deallocate(gefs_struc(n)%w221)

          deallocate(gefs_struc(n)%n112)
          deallocate(gefs_struc(n)%n122)
          deallocate(gefs_struc(n)%n212)
          deallocate(gefs_struc(n)%n222)
          deallocate(gefs_struc(n)%w112)
          deallocate(gefs_struc(n)%w122)
          deallocate(gefs_struc(n)%w212)
          deallocate(gefs_struc(n)%w222)
       endif
    enddo
    deallocate(gefs_struc)
    gefs_initialized = .false.
 endif

end subroutine finalize_gefs
