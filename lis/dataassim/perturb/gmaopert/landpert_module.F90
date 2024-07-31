!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
!  !MODULE: landpert_module.F90
! 
!  !DESCRIPTION: 
!  This module contains routines that can be used to perturb forcing 
!  fields such as precipitation and radiation
!  or model prognostic variables such as soil moisture or soil temperature
!
!  MUST initialize random seed and Pert\_ntrmdt by calling 
!  get\_pert() with initialize=.true. at the start of the driver program
!  (otherwise set initialize=.false.)

!   
!  !REVISION HISTORY: 
!  24Jan05 Rolf Reichle  Initial Specification
!  07Jun05 Rolf Reichle  more init options
!  07Jul05 Sujay Kumar   Initial Specification
! 
!EOP
module landpert_module
  implicit none
  
  save

  ! everything is private by default unless made public
  private

  public :: pert_param_type
  public :: allocate_pert_param
  ! --------------------------------------------------------------------
  !
  ! parameters for each kind of perturbation (precip, radiation, 
  ! soil moisture, etc)
  
  type :: pert_param_type
     
     character(40)     :: descr    ! 'precip', 'shortwave', etc
     integer           :: typ      ! add or multiply perturbation?
     
     ! max allowed normalized perturbation (relative to N(0,1))
     
     real              :: std_normal_max  
     
     ! if .true. enforce zeromean across ensemble 
     ! (implies mean=1 for multiplicative perturbations)
     ! (not applicable if only one ensemble member is done at a time)
     
     logical           :: zeromean        ! enforce zero mean across ensemble
     
     ! Mean and std are allowed to vary in space (dimension(N_x,N_y)).
     
     real, dimension(:,:), pointer :: mean      ! mean
     real, dimension(:,:), pointer :: std       ! standard deviation
     
     ! Cross-correlations between different kinds of perturbations
     ! (eg. between precip and shortwave perturbations) are allowed to vary
     ! in space (dimension(N_pert_kind,N_x,N_y)).
     
     real, dimension(:,:,:), pointer :: ccorr
     
     ! Spatial and temporal correlation scales must be constant in space.
     ! For non-zero cross-correlations they must also be the same for
     ! all kinds for perturbations (eg. if precip and radiation
     ! perturbations are cross-correlated, their xcorr, ycorr and tcorr
     ! must be the same).
     
     real             :: xcorr  ! correlation length along latitudes   [deg]
     real             :: ycorr  ! correlation length along longitudes  [deg]
     real             :: tcorr  ! temporal correlation length          [s]
     
  end type pert_param_type
contains

!BOP
! 
! !ROUTINE: allocate_pert_param
!  \label{allocate_pert_param}
! 
! !INTERFACE:   
  subroutine allocate_pert_param(N_pert, N_x, N_y, pert_param)
    
    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: N_pert, N_x, N_y    
    type(pert_param_type), dimension(:), pointer :: pert_param
! 
! !DESCRIPTION:
!   allocates memory for the perturbation data structures. 
! 
!EOP
    
    ! local variables
    
    integer :: k
    
    ! --------------------------------------------------------
    
    nullify(pert_param)
        
    allocate(pert_param(N_pert))
    
    do k=1,N_pert
       
       allocate(pert_param(k)%mean(N_x,N_y))
       allocate(pert_param(k)%std(N_x,N_y))
       allocate(pert_param(k)%ccorr(N_pert,N_x,N_y))
       
    end do
    
  end subroutine allocate_pert_param

end module landpert_module
