!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! 
!BOP
! 
! !MODULE: enks_types
! 
! !DESCRIPTION: 
! This module contains the definition of types for the GMAO Enks  
! 
! !REVISION HISTORY: 
!   reichle, 19 Jul 2005
!EOP 
module enks_types  

  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none
  save
  
  ! everything is private by default unless made public
  
  private
  
  public :: obs_type, obs_param_type
  
  public :: update_region_ind_type, update_region_type
  
  ! -----------------------------------------------------------------------
  
  ! obs_type is basic element of vector "Observations" (length N_obs), 
  ! which contains all observations of all types that are available
  ! at a given update time
  
  type :: obs_type
     
     integer          :: species   ! identifier for type of measurement
     integer          :: catnum    ! number of catchment in domain
     real             :: lon       ! longitude
     real             :: lat       ! latitude
     real             :: value     ! value of measurement
     real             :: std       ! obs error std
     logical          :: assim     ! .T. if assimilated, .F. if innov only
     integer          :: pert_type ! 0- additive, 1-multiplicative
  end type obs_type
  
  ! ----------------------------------------------------------------------
  !
  ! vector obs_param contains information about each species of observations 
  
  type :: obs_param_type
     
     integer          :: species   ! identifier for type of measurement
     character(40)    :: descr     ! description
     logical          :: assim     ! assimilate yes/no?
     logical          :: scale     ! scale yes/no?
     logical          :: getinnov  ! compute innovations? (.T. if assim==.T.)
     real             :: nodata    ! no-data-value
     character(len=LIS_CONST_PATH_LEN)  :: path      ! path to measurements file 
     character(80)    :: name      ! name identifier for measurements 
     character(len=LIS_CONST_PATH_LEN)  :: scalepath ! path to file with scaling parameters
     character(len=LIS_CONST_PATH_LEN) :: scalename ! filename for scaling parameters
     real             :: std       ! default obs error std

     real             :: std_normal_max  ! see pert_param_type
     logical          :: zeromean        ! see pert_param_type
     real             :: xcorr           ! see pert_param_type
     real             :: ycorr           ! see pert_param_type
     
  end type obs_param_type
  
  ! ----------------------------------------------------------------------
  
  type :: update_region_ind_type
     
     integer          :: cat_id       ! ID of catchment
     integer          :: region_num   ! number of region
     
  end type update_region_ind_type
  
  ! ----------------------------------------------------------------------
  
  type :: update_region_type
     
     integer          :: region_num   ! number of region
     integer          :: n_cats       ! number of catchments in region
     real             :: min_lon      ! lat/lon bounding box of region
     real             :: max_lon
     real             :: min_lat
     real             :: max_lat
          
  end type update_region_type

  
end module enks_types

! ================== EOF ===============================================
