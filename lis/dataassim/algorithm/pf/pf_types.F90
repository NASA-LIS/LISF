!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! 
!BOP
! 
! !MODULE: pf_types
! 
! !DESCRIPTION: 
! This module contains the definition of types for the GMAO Pf  
! 
! !REVISION HISTORY: 
!   reichle, 19 Jul 2005
!EOP 
module pf_types  
  
  implicit none
  save
  
  ! everything is private by default unless made public
  
  private
  
  public :: pf_obs_type, pf_obs_param_type
  
  ! -----------------------------------------------------------------------
  
  ! obs_type is basic element of vector "Observations" (length N_obs), 
  ! which contains all observations of all types that are available
  ! at a given update time
  
  type :: pf_obs_type
     
     integer          :: species   ! identifier for type of measurement
     integer          :: catnum    ! number of catchment in domain
     real             :: lon       ! longitude
     real             :: lat       ! latitude
     real             :: value     ! value of measurement
     real             :: std       ! obs error std
     logical          :: assim     ! .T. if assimilated, .F. if innov only
     integer          :: pert_type ! 0- additive, 1-multiplicative
  end type pf_obs_type
  
  ! ----------------------------------------------------------------------
  !
  ! vector obs_param contains information about each species of observations 
  
  type :: pf_obs_param_type
     
     integer          :: species   ! identifier for type of measurement
     character(40)    :: descr     ! description
     logical          :: assim     ! assimilate yes/no?
     logical          :: scale     ! scale yes/no?
     logical          :: getinnov  ! compute innovations? (.T. if assim==.T.)
     real             :: nodata    ! no-data-value
     character(200)   :: path      ! path to measurements file 
     character(80)    :: name      ! name identifier for measurements 
     character(200)   :: scalepath ! path to file with scaling parameters
     character(80)    :: scalename ! filename for scaling parameters
     real             :: std       ! default obs error std

     real             :: std_normal_max  ! see pert_param_type
     logical          :: zeromean        ! see pert_param_type
     real             :: xcorr           ! see pert_param_type
     real             :: ycorr           ! see pert_param_type
     
  end type pf_obs_param_type
    
end module pf_types

! ================== EOF ===============================================
