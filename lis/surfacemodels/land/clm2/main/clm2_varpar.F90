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

module clm2_varpar

  use LIS_precisionMod
  implicit none

  PRIVATE
!----------------------------------------------------------------------- 
! 
! Purpose: 
! land surface model array dimensions
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: clm2_varpar.F90,v 1.6 2004/11/24 22:57:13 jim Exp $
!-----------------------------------------------------------------------

! Define land surface 2-d grid. This sets the model resolution according
! to cpp directives LSMLON and LSMLAT in preproc.h. 

!  integer, parameter :: lsmlon = LSMLON  !maximum number of longitude points on lsm grid
!  integer, parameter :: lsmlat = LSMLAT  !number of latitude points on lsm grid

! Define maximum number of PFT patches per grid cell and set
! patch number for urban, lake, wetland, and glacier patches

#if (defined DGVM)
  integer, public, parameter :: maxpatch_pft = 13
#else

  integer, public :: maxpatch_pft               !maximum number of PFT subgrid patches per grid cell
#endif
  integer, public :: npatch_urban !urban   patch number: 1 to maxpatch
  integer, public :: npatch_lake  !lake    patch number: 1 to maxpatch
  integer, public :: npatch_wet   !wetland patch number: 1 to maxpatch
  integer, public :: npatch_gla   !glacier patch number: 1 to maxpatch
  integer, public :: maxpatch     !maximum number of subgrid patches per grid cell

! Define history file parameters

  integer, public , parameter :: maxhist      =   3             !max number of history files
  integer, public , parameter :: maxflds      = 200             !max number of fields in list
  integer, public , parameter :: max_slevflds =  75             !max number of active single-level fields
  integer, public , parameter :: max_mlevflds =  10             !max number of active multi-level fields (either snow or soil)
  integer, public , parameter :: maxalflds = max_slevflds + max_mlevflds !max number of active fields (all levels)

! Define number of level parameters

  integer, public, parameter :: nlevsoi     =  10   !number of soil layers
  integer, public, parameter :: nlevlak     =  10   !number of lake layers
  integer, public, parameter :: nlevsno     =   5   !maximum number of snow layers

! Define miscellaneous parameters

  integer, public, parameter :: numwat      =   5   !number of water types (soil, ice, 2 lakes, wetland)
  integer, public, parameter :: numpft      =  16   !number of plant types

! next variable used with DGVM
  integer, public, parameter :: npftpar     =  32   !number of pft parameters (in LPJ)
  integer, public, parameter :: numcol      =   8   !number of soil color types
  integer, public, parameter :: numrad      =   2   !number of solar radiation bands: vis, nir
  integer, public, parameter :: ndst        =   4   !number of dust size classes
  integer, public, parameter :: dst_src_nbr =   3   !number of size distns in src soil
  integer, public, parameter :: nvoc        =   4   !number of voc categories

! Define parameters for RTM river routing model

  integer, public, parameter :: rtmlon = 720  !# of rtm longitudes
  integer, public, parameter :: rtmlat = 360  !# of rtm latitudes

end module clm2_varpar

