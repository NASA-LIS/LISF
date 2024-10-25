!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: finalize_geos5fcst
!  \label{finalize_geos5fcst}
!
! !REVISION HISTORY: 
! 7 Mar 2013: Sujay Kumar, initial specification
! 
! !INTERFACE:
subroutine finalize_geos5fcst(findex)
! !USES:
  use LDT_coreMod,          only : LDT_rc
  use geos5fcst_forcingMod, only : geos5fcst_struc, geos5fcst_initialized
!
! !DESCRIPTION:
!  Routine to finalize GEOS5 forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex
  
  if ( geos5fcst_initialized ) then 
     do n=1,LDT_rc%nnest
       if ( LDT_rc%met_gridtransform(findex).eq."bilinear" ) then 

          deallocate(geos5fcst_struc(n)%n111)
          deallocate(geos5fcst_struc(n)%n121)
          deallocate(geos5fcst_struc(n)%n211)
          deallocate(geos5fcst_struc(n)%n221)
          deallocate(geos5fcst_struc(n)%w111)
          deallocate(geos5fcst_struc(n)%w121)
          deallocate(geos5fcst_struc(n)%w211)
          deallocate(geos5fcst_struc(n)%w221)
       endif
    enddo
    deallocate(geos5fcst_struc)
    geos5fcst_initialized = .false.
 endif

end subroutine finalize_geos5fcst
