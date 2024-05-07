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
!  !MODULE: landpert_routines.F90
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
!  07Jul05 Sujay Kumar   Specification in LIS
!  27Feb19 Mahdi Navari  domin id was added to the initializing the seed 
!                        this will solve the checkerboard pattern and 
!                        stripe patterns in the DA results
!
!EOP 
module landpert_routines

  use ESMF
  use landpert_module
  use LIS_ran2_gasdev
  use random_fields
  use LIS_logMod
  use LIS_coreMod
  use LIS_historyMod
  use LIS_DAobservationsMod

  implicit none

  
  ! everything is private by default unless made public
  
  private

  public :: isXYcorrEnabled
  public :: get_pert
  public :: assemble_pert_param
  public :: assemble_pert_param_obs
  public :: get_sqrt_corr_matrix
  public :: get_init_Pert_rseed

  type :: grid_def_type              
     integer          :: N_x        ! number of longitudinal nodes
     integer          :: N_y        ! number of latitudinal nodes
     integer          :: N_x_fft
     integer          :: N_y_fft
     real             :: llx        ! lower left latitude in deg
     real             :: lly        ! lower left longitude in deg
     real             :: dx         ! longitude spacing in deg
     real             :: dy         ! latitude spacing in deg
  end type grid_def_type

!BOP 
! 
! !ROUTINE: assemble_pert_param
! \label{assemble_pert_param}
! 
! !INTERFACE:
  interface assemble_pert_param
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure assemble_pert_param_global
     module procedure assemble_pert_param_local
     
  end interface

!BOP 
! 
! !ROUTINE: assemble_pert_param_obs
! \label{assemble_pert_param_obs}
! 
! !INTERFACE:
  interface assemble_pert_param_obs
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure assemble_pert_param_global_obs
     module procedure assemble_pert_param_local_obs
     
  end interface

contains

! 
!BOP
! 
! !ROUTINE: isXYcorrEnabled
! \label{isXYcorrEnabled
! 
! !INTERFACE: 
  function isXYcorrEnabled(Pert_State) result(finish)
    implicit none
! !ARGUMENTS: 
    logical          :: finish
    type(ESMF_State), intent(in) :: Pert_State
! 
! !DESCRIPTION: 
!  This routine checks to see if horizontal correlations are
!  enabled in the perturbation settings
!EOP
    integer                  :: i
    character*100, pointer   :: pertobjs(:)
    type(ESMF_Field)         :: pertField
    integer                  :: objcount
    integer                  :: status
    real                     :: xcorr 
    real                     :: ycorr 

    call ESMF_StateGet(Pert_State, itemCount=objcount,rc=status)
    call LIS_verify(status)
    
    finish = .false. 

    allocate(pertobjs(objcount))

    call ESMF_StateGet(Pert_State, itemNameList=pertobjs, rc=status)
    call LIS_verify(status)

    do i=1,objcount

       call ESMF_StateGet(Pert_State, pertobjs(i),pertField,&
            rc=status)
       call LIS_verify(status)

       call ESMF_AttributeGet(pertField,"X Correlation Scale",xcorr,&
            rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField,"Y Correlation Scale",ycorr,&
            rc=status)
       call LIS_verify(status)

       if(xcorr.gt.0.or.ycorr.gt.0) then 
          finish = .true.
       endif
    enddo
    
    deallocate(pertobjs)
    
  end function isXYcorrEnabled
!BOP
! 
! !ROUTINE: assemble_pert_param_local
!  \label{assemble_pert_param_local}
!
! !INTERFACE:  
  subroutine assemble_pert_param_local( Pert_State, N_x, N_y,               & 
       N_pert, pert_param )
! 
! !DESCRIPTION:
!
! sample subroutine that demonstrates how pert\_param can be 
! assembled for forcing perturbations
!
! return N\_force\_pert, allocate and assemble structure pert\_param
!
! make sure order of fields is compatible with your driver
!
! forcing field 1 = precip
! forcing field 2 = shortwave
! forcing field 3 = longwave
! forcing field 4 = air temperature
!
! !REVISION HISTORY:     
! reichle, 11 Feb 2005
! reichle,  8 Jun 2005
!
! !USES:     

    
    implicit none
! !ARGUMENTS:     
    type(ESMF_State), intent(in) :: Pert_State
    integer, intent(in) :: N_x, N_y
    integer, intent(in) :: N_pert 
    type(pert_param_type), dimension(:), pointer :: pert_param ! out
!EOP    
    ! -----------------------------------------------------------------
    
    ! forcing perturbation parameters
    
    ! # of forcing variables that are perturbed
    ! (currently pcp, sw, lw, tair)
    
!    integer, parameter :: N_tmp  = 4
    
    ! limit on range of random numbers:
    !  specify max absolute value allowed to be drawn from a standard
    !  normal distribution
    
    real     :: std_normal_max 
    
    ! decide whether to ensure zero mean for synthetically generated errors
    ! (IMPORTANT: this will only have an effect for N_ens>2!!!)
    
    logical  :: zeromean
    
    ! temporal correlation scale
    
    real     :: tcorr   ! 86400               ! [s]
    
    ! horizontal correlation scales 
    
    real     :: xcorr 
    real     :: ycorr 
    
    ! perturbations are either
    !
    !   typ=0: additive, mean=0
    !   typ=1: multiplicative and lognormal, mean=1
    
    integer    :: typ
    real, pointer :: std(:)
    integer             :: c,r
        
    ! correlation coefficients -1 <= rho <= 1     
    !
    ! (these numbers are made up and are not tested with any data!)
    
    real, parameter :: rho_pcp_sw   = -.5
    real, parameter :: rho_pcp_lw   =  .5
    real, parameter :: rho_pcp_tair = -.3
    real, parameter :: rho_sw_lw    = -.5
    real, parameter :: rho_sw_tair  =  .5
    real, parameter :: rho_lw_tair  =  .4
    
    ! ---------------------------------------------------------------------
    !
    ! local variables
    
    integer   :: i, j, k, l
    real      :: tmpreal
!    real, dimension(N_tmp,N_tmp) :: tmpmat1, tmpmat2
    real, dimension(N_pert,N_pert) :: tmpmat1, tmpmat2
    integer                      :: objcount
    integer                      :: status
    character*100, pointer       :: pertobjs(:)
    real,          pointer       :: ccorr(:)
    type(ESMF_Field)             :: pertField
    integer                      :: tmp
    integer                      :: n 

    n = 1 ! assuming for now -- this must be passed as arguments - TBD

    ! ---------------------------------------------------------------------
    !
    ! allocate pert_param (must not be associated at this time)

    if (associated(pert_param)) then
       write (*,*) 'assemble_pert_param(): this needs work...'
       write (*,*) 'stopping'
       stop
    end if
    
!    N_pert = N_tmp
    
    call allocate_pert_param(N_pert, N_x, N_y, pert_param)

    call ESMF_StateGet(Pert_State, itemCount=objcount,rc=status)
    call LIS_verify(status,'ESMF_StateGet1: Pert_state in landpert_routines')
    
    allocate(pertobjs(objcount))
    allocate(ccorr(objcount))

    call ESMF_StateGet(Pert_State, itemNameList=pertobjs, rc=status)
    call LIS_verify(status,'ESMF_StateGet2: Pert_state in landpert_routines')

    do i=1,objcount

       call ESMF_StateGet(Pert_State, pertobjs(i),pertField,&
            rc=status)
       call LIS_verify(status,'ESMF_StateGet3: Pert_state in landpert_routines')

       call ESMF_AttributeGet(pertField,"Perturbation Type",typ,&
            rc=status)
       call LIS_verify(status,'ESMF_AttributeGet: perttype in landpert_routines')

       call ESMF_AttributeGet(pertField,"Std Normal Max",std_normal_max,&
            rc=status)
       call LIS_verify(status,'ESMF_AttributeGet: stdnormalmax in landpert_routines')
       call ESMF_AttributeGet(pertField,"Ensure Zero Mean",tmp,&
            rc=status)
       call LIS_verify(status,'ESMF_AttributeGet: zeromean in landpert_routines')
       zeromean = .false. 
       if(tmp.eq.1) zeromean = .true.
       
       allocate(std(LIS_rc%ngrid(n)))

       call ESMF_AttributeGet(pertField,"Temporal Correlation Scale",tcorr,&
            rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField,"X Correlation Scale",xcorr,&
            rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField,"Y Correlation Scale",ycorr,&
            rc=status)
       call LIS_verify(status)
       
       if(LIS_rc%ngrid(n).gt.0) then 
          call ESMF_AttributeGet(pertField,"Standard Deviation",std,&
               rc=status)
          call LIS_verify(status)
       endif


       call ESMF_AttributeGet(pertField,"Cross Correlation Strength",&
            ccorr, rc=status)
       call LIS_verify(status)

    ! ----------------------------------------
    !
    ! copy inputs into structure
    
       pert_param(i)%descr          = pertobjs(i)
       pert_param(i)%typ            = typ
       pert_param(i)%std_normal_max = std_normal_max
       pert_param(i)%zeromean       = zeromean
       pert_param(i)%tcorr          = tcorr
       pert_param(i)%xcorr          = xcorr
       pert_param(i)%ycorr          = ycorr
       
!       pert_param(i)%std            = std

       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)             
             if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                pert_param(i)%std(c,r) = &
                     std(LIS_domain(n)%gindex(c,r))
             else
                pert_param(i)%std(c,r) = -9999.0
             endif
          enddo
       enddo

!       pert_param(i)%ccorr(:,:,:)   = 0.0
!       pert_param(i)%ccorr(i,:,:)   = 1.0
       do j=1,objcount
          pert_param(i)%ccorr(j,:,:)  = ccorr(j)
       enddo
    enddo
    ! -------------------------------------------------------------
    !
    ! set mean and (if needed) modify standard deviation according to 'typ'
    ! (additive or multiplicative perturbations)    

    
    do k=1,N_pert
       
       select case (pert_param(k)%typ)
          
       case (0)         ! additive (mean=0, std as above)
          
          pert_param(k)%mean  = 0.
          
       case (1)         ! multiplicative and lognormal (mean=1)
          
          do i=1,N_x
             do j=1,N_y
                
                tmpreal = pert_param(k)%std(i,j) 
                
                tmpreal = log( 1. + tmpreal**2)
                
                pert_param(k)%mean(i,j) = - .5*tmpreal
                pert_param(k)%std(i,j)  = sqrt(tmpreal)
                
             end do
          end do
          
       case default
          
          write (*,*) 'assemble_pert_param(): encountered unknown'
          write (*,*)  k, pert_param(k)%typ
          write (*,*) 'type of error, stopping...'
          stop
          
       end select
       
    end do
    
    
    
      ! echo part of pert_param (mean, std, and ccorr for i=1, j=1 only):
    if(LIS_rc%plevel.gt.1) then 
       do i=1,N_pert
          
          write (*,*) 'pert_param(',i,')%descr=', &
               pert_param(i)%descr
          write (*,*) 'pert_param(',i,')%typ=', &
               pert_param(i)%typ
          write (*,*) 'pert_param(',i,')%zeromean=', &
               pert_param(i)%zeromean
          write (*,*) 'pert_param(',i,')%std_normal_max=', &
               pert_param(i)%std_normal_max
          write (*,*) 'pert_param(',i,')%xcorr=', &
               pert_param(i)%xcorr
          write (*,*) 'pert_param(',i,')%ycorr=', &
               pert_param(i)%ycorr
          write (*,*) 'pert_param(',i,')%tcorr=', &
               pert_param(i)%tcorr
          
          write (*,*) 'pert_param(',i,')%mean(1,1)=', &
               pert_param(i)%mean(1,1)
          write (*,*) 'pert_param(',i,')%std(1,1)=', &
               pert_param(i)%std(1,1)
          
          do j=1,N_pert
             
             write (*,*) 'pert_param(',i,')%ccorr(',j,',1,1)=', &
                  pert_param(i)%ccorr(j,1,1)
             
          end do
          
       end do
    end if
    
    
    

    ! -------------------------------------------------------------
    !
    ! compute sqrt of correlation matrix for each grid point
    
    do i=1,N_x
       do j=1,N_y

          ! extract local correlation matrix for grid point (i,j)
          
          do k=1,N_pert
             do l=1,N_pert
                
                tmpmat1(k,l) = pert_param(k)%ccorr(l,i,j)
                
             end do
          end do
          
          ! compute sqrt of local correlation matrix
          
          call get_sqrt_corr_matrix( N_pert, tmpmat1, tmpmat2 )
          
          ! overwrite cross-correlations in forcepert_param with square 
          ! root of cross-correlation matrix
          
          do k=1,N_pert
             do l=1,N_pert
                
                pert_param(k)%ccorr(l,i,j) = tmpmat2(k,l)
                
             end do
          end do
          
       end do
    end do
    
    deallocate(pertobjs)
    deallocate(ccorr)

  end subroutine assemble_pert_param_local


!BOP
! 
! !ROUTINE: assemble_pert_param_global
!  \label{assemble_pert_param_global}
!
! !INTERFACE:  
  subroutine assemble_pert_param_global( Pert_State, N_x, N_y,  & 
       N_pert, pert_param, glob_flag)
! 
! !DESCRIPTION:
!
! sample subroutine that demonstrates how pert\_param can be 
! assembled for forcing perturbations
!
! return N\_force\_pert, allocate and assemble structure pert\_param
!
! make sure order of fields is compatible with your driver
!
! forcing field 1 = precip
! forcing field 2 = shortwave
! forcing field 3 = longwave
! forcing field 4 = air temperature
!
! !REVISION HISTORY:     
! reichle, 11 Feb 2005
! reichle,  8 Jun 2005
!
! !USES:     

    
    implicit none
! !ARGUMENTS:     
    type(ESMF_State), intent(in)                 :: Pert_State
    integer, intent(in)                          :: N_x, N_y
    integer, intent(in)                          :: N_pert 
    type(pert_param_type), dimension(:), pointer :: pert_param ! out
    integer                                      :: glob_flag
!EOP    
    ! -----------------------------------------------------------------
    
    ! forcing perturbation parameters
    
    ! # of forcing variables that are perturbed
    ! (currently pcp, sw, lw, tair)
    
!    integer, parameter :: N_tmp  = 4
    
    ! limit on range of random numbers:
    !  specify max absolute value allowed to be drawn from a standard
    !  normal distribution
    
    real     :: std_normal_max 
    
    ! decide whether to ensure zero mean for synthetically generated errors
    ! (IMPORTANT: this will only have an effect for N_ens>2!!!)
    
    logical  :: zeromean
    
    ! temporal correlation scale
    
    real     :: tcorr   ! 86400               ! [s]
    
    ! horizontal correlation scales 
    
    real     :: xcorr 
    real     :: ycorr 
    
    ! perturbations are either
    !
    !   typ=0: additive, mean=0
    !   typ=1: multiplicative and lognormal, mean=1
    
    integer    :: typ
    real, pointer :: std(:)
    integer             :: c,r
        
    ! correlation coefficients -1 <= rho <= 1     
    !
    ! (these numbers are made up and are not tested with any data!)
    
    real, parameter :: rho_pcp_sw   = -.5
    real, parameter :: rho_pcp_lw   =  .5
    real, parameter :: rho_pcp_tair = -.3
    real, parameter :: rho_sw_lw    = -.5
    real, parameter :: rho_sw_tair  =  .5
    real, parameter :: rho_lw_tair  =  .4
    
    ! ---------------------------------------------------------------------
    !
    ! local variables
    
    integer   :: i, j, k, l
    real      :: tmpreal
!    real, dimension(N_tmp,N_tmp) :: tmpmat1, tmpmat2
    real, dimension(N_pert,N_pert) :: tmpmat1, tmpmat2
    integer                      :: objcount
    integer                      :: status
    character*100, pointer       :: pertobjs(:)
    real,          pointer       :: ccorr(:)
    type(ESMF_Field)             :: pertField
    integer                      :: tmp
    integer                      :: n 
    real,          allocatable   :: global_std(:,:)
    n = 1 ! assuming for now -- this must be passed as arguments - TBD

    ! ---------------------------------------------------------------------
    !
    ! allocate pert_param (must not be associated at this time)

    if (associated(pert_param)) then
       write (*,*) 'assemble_pert_param(): this needs work...'
       write (*,*) 'stopping'
       stop
    end if
    
!    N_pert = N_tmp
    
    if(LIS_masterproc) then 
       call allocate_pert_param(N_pert, N_x, N_y, pert_param)
    endif

    call ESMF_StateGet(Pert_State, itemCount=objcount,rc=status)
    call LIS_verify(status)
    
    allocate(pertobjs(objcount))
    allocate(ccorr(objcount))

    call ESMF_StateGet(Pert_State, itemNameList=pertobjs, rc=status)
    call LIS_verify(status)

    do i=1,objcount

       call ESMF_StateGet(Pert_State, pertobjs(i),pertField,&
            rc=status)
       call LIS_verify(status)

       call ESMF_AttributeGet(pertField,"Perturbation Type",typ,&
            rc=status)
       call LIS_verify(status)

       call ESMF_AttributeGet(pertField,"Std Normal Max",std_normal_max,&
            rc=status)
       call LIS_verify(status)
       call ESMF_AttributeGet(pertField,"Ensure Zero Mean",tmp,&
            rc=status)
       call LIS_verify(status)
       zeromean = .false. 
       if(tmp.eq.1) zeromean = .true.
       
       allocate(std(LIS_rc%ngrid(n)))

       call ESMF_AttributeGet(pertField,"Temporal Correlation Scale",tcorr,&
            rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField,"X Correlation Scale",xcorr,&
            rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField,"Y Correlation Scale",ycorr,&
            rc=status)
       call LIS_verify(status)
       
       if(LIS_rc%ngrid(n).gt.0) then 
          call ESMF_AttributeGet(pertField,"Standard Deviation",std,&
               rc=status)
          call LIS_verify(status)
       endif


       call ESMF_AttributeGet(pertField,"Cross Correlation Strength",&
            ccorr, rc=status)
       call LIS_verify(status)

       call LIS_gather_1dgrid_to_2dgrid(n, global_std, std)
    ! ----------------------------------------
    !
    ! copy inputs into structure
    
       if(LIS_masterproc) then 
          pert_param(i)%descr          = pertobjs(i)
          pert_param(i)%typ            = typ
          pert_param(i)%std_normal_max = std_normal_max
          pert_param(i)%zeromean       = zeromean
          pert_param(i)%tcorr          = tcorr
          pert_param(i)%xcorr          = xcorr
          pert_param(i)%ycorr          = ycorr
       
          pert_param(i)%std            = global_std

          
          do j=1,objcount
             pert_param(i)%ccorr(j,:,:)  = ccorr(j)
          enddo
       endif
       deallocate(global_std)
    enddo
    ! -------------------------------------------------------------
    !
    ! set mean and (if needed) modify standard deviation according to 'typ'
    ! (additive or multiplicative perturbations)    

    if(LIS_masterproc) then 
       do k=1,N_pert
          
          select case (pert_param(k)%typ)
             
          case (0)         ! additive (mean=0, std as above)
             
             pert_param(k)%mean  = 0.
             
          case (1)         ! multiplicative and lognormal (mean=1)
             
             do i=1,N_x
                do j=1,N_y
                   
                   tmpreal = pert_param(k)%std(i,j) 
                   
                   tmpreal = log( 1. + tmpreal**2)
                   
                   pert_param(k)%mean(i,j) = - .5*tmpreal
                   pert_param(k)%std(i,j)  = sqrt(tmpreal)
                   
                end do
             end do
             
          case default
             
             write (*,*) 'assemble_pert_param(): encountered unknown'
             write (*,*)  k, pert_param(k)%typ
             write (*,*) 'type of error, stopping...'
             stop
          
          end select
       enddo        

    ! -------------------------------------------------------------
    !
    ! compute sqrt of correlation matrix for each grid point
    
       do i=1,N_x
          do j=1,N_y

          ! extract local correlation matrix for grid point (i,j)
          
             do k=1,N_pert
                do l=1,N_pert
                   
                   tmpmat1(k,l) = pert_param(k)%ccorr(l,i,j)
                   
                end do
             end do
          
          ! compute sqrt of local correlation matrix
          
             call get_sqrt_corr_matrix( N_pert, tmpmat1, tmpmat2 )
          
          ! overwrite cross-correlations in forcepert_param with square 
          ! root of cross-correlation matrix
          
             do k=1,N_pert
                do l=1,N_pert
                
                   pert_param(k)%ccorr(l,i,j) = tmpmat2(k,l)
                   
                end do
             end do
             
          end do
       end do
    endif

    deallocate(pertobjs)
    deallocate(ccorr)

  end subroutine assemble_pert_param_global

!BOP
! 
! !ROUTINE: assemble_pert_param_local_obs
!  \label{assemble_pert_param_local_obs}
!
! !INTERFACE:  
  subroutine assemble_pert_param_local_obs( kk, Pert_State, N_x, N_y,    & 
       N_pert, pert_param )
! 
! !DESCRIPTION:
!
! sample subroutine that demonstrates how pert\_param can be 
! assembled for forcing perturbations
!
! return N\_force\_pert, allocate and assemble structure pert\_param
!
! make sure order of fields is compatible with your driver
!
! forcing field 1 = precip
! forcing field 2 = shortwave
! forcing field 3 = longwave
! forcing field 4 = air temperature
!
! !REVISION HISTORY:     
! reichle, 11 Feb 2005
! reichle,  8 Jun 2005
!
! !USES:     

    
    implicit none
! !ARGUMENTS:     
    integer, intent(in)                          :: kk 
    type(ESMF_State), intent(in)                 :: Pert_State
    integer, intent(in)                          :: N_x, N_y
    integer, intent(in)                          :: N_pert 
    type(pert_param_type), dimension(:), pointer :: pert_param ! out
!EOP    
    ! -----------------------------------------------------------------
    
    ! forcing perturbation parameters
    
    ! # of forcing variables that are perturbed
    ! (currently pcp, sw, lw, tair)
    
!    integer, parameter :: N_tmp  = 4
    
    ! limit on range of random numbers:
    !  specify max absolute value allowed to be drawn from a standard
    !  normal distribution
    
    real     :: std_normal_max 
    
    ! decide whether to ensure zero mean for synthetically generated errors
    ! (IMPORTANT: this will only have an effect for N_ens>2!!!)
    
    logical  :: zeromean
    
    ! temporal correlation scale
    
    real     :: tcorr   ! 86400               ! [s]
    
    ! horizontal correlation scales 
    
    real     :: xcorr 
    real     :: ycorr 
    
    ! perturbations are either
    !
    !   typ=0: additive, mean=0
    !   typ=1: multiplicative and lognormal, mean=1
    
    integer    :: typ
    real, pointer :: std(:)
    integer             :: c,r
        
    ! correlation coefficients -1 <= rho <= 1     
    !
    ! (these numbers are made up and are not tested with any data!)
    
    real, parameter :: rho_pcp_sw   = -.5
    real, parameter :: rho_pcp_lw   =  .5
    real, parameter :: rho_pcp_tair = -.3
    real, parameter :: rho_sw_lw    = -.5
    real, parameter :: rho_sw_tair  =  .5
    real, parameter :: rho_lw_tair  =  .4
    
    ! ---------------------------------------------------------------------
    !
    ! local variables
    
    integer   :: i, j, k, l
    real      :: tmpreal
!    real, dimension(N_tmp,N_tmp) :: tmpmat1, tmpmat2
    real, dimension(N_pert,N_pert) :: tmpmat1, tmpmat2
    integer                      :: objcount
    integer                      :: status
    character*100, pointer       :: pertobjs(:)
    real,          pointer       :: ccorr(:)
    type(ESMF_Field)             :: pertField
    integer                      :: tmp
    integer                      :: n 

    n = 1 ! assuming for now -- this must be passed as arguments - TBD
    
    ! ---------------------------------------------------------------------
    !
    ! allocate pert_param (must not be associated at this time)

    if (associated(pert_param)) then
       write (*,*) 'assemble_pert_param(): this needs work...'
       write (*,*) 'stopping'
       stop
    end if
    
!    N_pert = N_tmp
    
    call allocate_pert_param(N_pert, N_x, N_y, pert_param)

    call ESMF_StateGet(Pert_State, itemCount=objcount,rc=status)
    call LIS_verify(status,'ESMF_StateGet1: Pert_state in landpert_routines')
    
    allocate(pertobjs(objcount))
    allocate(ccorr(objcount))

    call ESMF_StateGet(Pert_State, itemNameList=pertobjs, rc=status)
    call LIS_verify(status,'ESMF_StateGet2: pert_state in landpert_routines')

    do i=1,objcount

       call ESMF_StateGet(Pert_State, pertobjs(i),pertField,&
            rc=status)
       call LIS_verify(status,'ESMF_StateGet3: pert_state in landpert_routines')

       call ESMF_AttributeGet(pertField,"Perturbation Type",typ,&
            rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeGet failed for perturbation type')

       call ESMF_AttributeGet(pertField,"Std Normal Max",std_normal_max,&
            rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeGet failed for std normal max')
       call ESMF_AttributeGet(pertField,"Ensure Zero Mean",tmp,&
            rc=status)
       call LIS_verify(status, &
            'ESMF_AttributeGet failed for ensure zero mean')
       zeromean = .false. 
       if(tmp.eq.1) zeromean = .true.
       
       allocate(std(LIS_rc%obs_ngrid(kk)))

       call ESMF_AttributeGet(pertField,"Temporal Correlation Scale",tcorr,&
            rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeGet failed for temporal correlation scale')
       
       call ESMF_AttributeGet(pertField,"X Correlation Scale",xcorr,&
            rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeGet failed for x correlation scale')
       
       call ESMF_AttributeGet(pertField,"Y Correlation Scale",ycorr,&
            rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeGet failed for y correlation scale')

       if(LIS_rc%obs_ngrid(kk).gt.0) then 
          call ESMF_AttributeGet(pertField,"Standard Deviation",std,&
               rc=status)
          call LIS_verify(status,&
               'ESMF_AttributeGet failed for standard deviation')
       endif


       call ESMF_AttributeGet(pertField,"Cross Correlation Strength",&
            ccorr, rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeGet failed for cross correlation strength')

    ! ----------------------------------------
    !
    ! copy inputs into structure
    
       pert_param(i)%descr          = pertobjs(i)
       pert_param(i)%typ            = typ
       pert_param(i)%std_normal_max = std_normal_max
       pert_param(i)%zeromean       = zeromean
       pert_param(i)%tcorr          = tcorr
       pert_param(i)%xcorr          = xcorr
       pert_param(i)%ycorr          = ycorr
       
!       pert_param(i)%std            = std
       do r=1,LIS_rc%obs_lnr(kk)
          do c=1,LIS_rc%obs_lnc(kk)             
             if(LIS_obs_domain(n,kk)%gindex(c,r).ne.-1) then
                pert_param(i)%std(c,r) = &
                     std(LIS_obs_domain(n,kk)%gindex(c,r))
             else
                pert_param(i)%std(c,r) = -9999.0
             endif
          enddo
       enddo
!       pert_param(i)%ccorr(:,:,:)   = 0.0
!       pert_param(i)%ccorr(i,:,:)   = 1.0
       do j=1,objcount
          pert_param(i)%ccorr(j,:,:)  = ccorr(j)
       enddo
    enddo
    ! -------------------------------------------------------------
    !
    ! set mean and (if needed) modify standard deviation according to 'typ'
    ! (additive or multiplicative perturbations)    

    
    do k=1,N_pert
       
       select case (pert_param(k)%typ)
          
       case (0)         ! additive (mean=0, std as above)
          
          pert_param(k)%mean  = 0.
          
       case (1)         ! multiplicative and lognormal (mean=1)
          
          do i=1,N_x
             do j=1,N_y
                
                tmpreal = pert_param(k)%std(i,j) 
                
                tmpreal = log( 1. + tmpreal**2)
                
                pert_param(k)%mean(i,j) = - .5*tmpreal
                pert_param(k)%std(i,j)  = sqrt(tmpreal)
                
             end do
          end do
          
       case default
          
          write (*,*) 'assemble_pert_param(): encountered unknown'
          write (*,*)  k, pert_param(k)%typ
          write (*,*) 'type of error, stopping...'
          stop
          
       end select
       
    end do
    
    
    
      ! echo part of pert_param (mean, std, and ccorr for i=1, j=1 only):
    if(LIS_rc%plevel.gt.1) then 
       do i=1,N_pert
          
          write (*,*) 'pert_param(',i,')%descr=', &
               pert_param(i)%descr
          write (*,*) 'pert_param(',i,')%typ=', &
               pert_param(i)%typ
          write (*,*) 'pert_param(',i,')%zeromean=', &
               pert_param(i)%zeromean
          write (*,*) 'pert_param(',i,')%std_normal_max=', &
               pert_param(i)%std_normal_max
          write (*,*) 'pert_param(',i,')%xcorr=', &
               pert_param(i)%xcorr
          write (*,*) 'pert_param(',i,')%ycorr=', &
               pert_param(i)%ycorr
          write (*,*) 'pert_param(',i,')%tcorr=', &
               pert_param(i)%tcorr
          
          write (*,*) 'pert_param(',i,')%mean(1,1)=', &
               pert_param(i)%mean(1,1)
          write (*,*) 'pert_param(',i,')%std(1,1)=', &
               pert_param(i)%std(1,1)
          
          do j=1,N_pert
             
             write (*,*) 'pert_param(',i,')%ccorr(',j,',1,1)=', &
                  pert_param(i)%ccorr(j,1,1)
             
          end do
          
       end do
    end if
    
    
    

    ! -------------------------------------------------------------
    !
    ! compute sqrt of correlation matrix for each grid point
    
    do i=1,N_x
       do j=1,N_y

          ! extract local correlation matrix for grid point (i,j)
          
          do k=1,N_pert
             do l=1,N_pert
                
                tmpmat1(k,l) = pert_param(k)%ccorr(l,i,j)
                
             end do
          end do
          
          ! compute sqrt of local correlation matrix
          
          call get_sqrt_corr_matrix( N_pert, tmpmat1, tmpmat2 )
          
          ! overwrite cross-correlations in forcepert_param with square 
          ! root of cross-correlation matrix
          
          do k=1,N_pert
             do l=1,N_pert
                
                pert_param(k)%ccorr(l,i,j) = tmpmat2(k,l)
                
             end do
          end do
          
       end do
    end do
    
    deallocate(pertobjs)
    deallocate(ccorr)

  end subroutine assemble_pert_param_local_obs


!BOP
! 
! !ROUTINE: assemble_pert_param_global_obs
!  \label{assemble_pert_param_global_obs}
!
! !INTERFACE:  
  subroutine assemble_pert_param_global_obs(kk, Pert_State, N_x, N_y,  & 
       N_pert, pert_param, glob_flag)
! 
! !DESCRIPTION:
!
! sample subroutine that demonstrates how pert\_param can be 
! assembled for forcing perturbations
!
! return N\_force\_pert, allocate and assemble structure pert\_param
!
! make sure order of fields is compatible with your driver
!
! forcing field 1 = precip
! forcing field 2 = shortwave
! forcing field 3 = longwave
! forcing field 4 = air temperature
!
! !REVISION HISTORY:     
! reichle, 11 Feb 2005
! reichle,  8 Jun 2005
!
! !USES:     

    
    implicit none
! !ARGUMENTS:     
    integer,          intent(in)                 :: kk
    type(ESMF_State), intent(in)                 :: Pert_State
    integer, intent(in)                          :: N_x, N_y
    integer, intent(in)                          :: N_pert 
    type(pert_param_type), dimension(:), pointer :: pert_param ! out
    integer                                      :: glob_flag
!EOP    
    ! -----------------------------------------------------------------
    
    ! forcing perturbation parameters
    
    ! # of forcing variables that are perturbed
    ! (currently pcp, sw, lw, tair)
    
!    integer, parameter :: N_tmp  = 4
    
    ! limit on range of random numbers:
    !  specify max absolute value allowed to be drawn from a standard
    !  normal distribution
    
    real     :: std_normal_max 
    
    ! decide whether to ensure zero mean for synthetically generated errors
    ! (IMPORTANT: this will only have an effect for N_ens>2!!!)
    
    logical  :: zeromean
    
    ! temporal correlation scale
    
    real     :: tcorr   ! 86400               ! [s]
    
    ! horizontal correlation scales 
    
    real     :: xcorr 
    real     :: ycorr 
    
    ! perturbations are either
    !
    !   typ=0: additive, mean=0
    !   typ=1: multiplicative and lognormal, mean=1
    
    integer    :: typ
    real, pointer :: std(:)
    integer             :: c,r
        
    ! correlation coefficients -1 <= rho <= 1     
    !
    ! (these numbers are made up and are not tested with any data!)
    
    real, parameter :: rho_pcp_sw   = -.5
    real, parameter :: rho_pcp_lw   =  .5
    real, parameter :: rho_pcp_tair = -.3
    real, parameter :: rho_sw_lw    = -.5
    real, parameter :: rho_sw_tair  =  .5
    real, parameter :: rho_lw_tair  =  .4
    
    ! ---------------------------------------------------------------------
    !
    ! local variables
    
    integer   :: i, j, k, l
    real      :: tmpreal
!    real, dimension(N_tmp,N_tmp) :: tmpmat1, tmpmat2
    real, dimension(N_pert,N_pert) :: tmpmat1, tmpmat2
    integer                      :: objcount
    integer                      :: status
    character*100, pointer       :: pertobjs(:)
    real,          pointer       :: ccorr(:)
    type(ESMF_Field)             :: pertField
    integer                      :: tmp
    integer                      :: n 
    real,          allocatable   :: global_std(:,:)
    n = 1 ! assuming for now -- this must be passed as arguments - TBD

    ! ---------------------------------------------------------------------
    !
    ! allocate pert_param (must not be associated at this time)

    if (associated(pert_param)) then
       write (*,*) 'assemble_pert_param(): this needs work...'
       write (*,*) 'stopping'
       stop
    end if
    
!    N_pert = N_tmp
    
    if(LIS_masterproc) then 
       call allocate_pert_param(N_pert, N_x, N_y, pert_param)
    endif

    call ESMF_StateGet(Pert_State, itemCount=objcount,rc=status)
    call LIS_verify(status)
    
    allocate(pertobjs(objcount))
    allocate(ccorr(objcount))

    call ESMF_StateGet(Pert_State, itemNameList=pertobjs, rc=status)
    call LIS_verify(status)

    do i=1,objcount

       call ESMF_StateGet(Pert_State, pertobjs(i),pertField,&
            rc=status)
       call LIS_verify(status)

       call ESMF_AttributeGet(pertField,"Perturbation Type",typ,&
            rc=status)
       call LIS_verify(status)

       call ESMF_AttributeGet(pertField,"Std Normal Max",std_normal_max,&
            rc=status)
       call LIS_verify(status)
       call ESMF_AttributeGet(pertField,"Ensure Zero Mean",tmp,&
            rc=status)
       call LIS_verify(status)
       zeromean = .false. 
       if(tmp.eq.1) zeromean = .true.
       
       allocate(std(LIS_rc%obs_ngrid(kk)))

       call ESMF_AttributeGet(pertField,"Temporal Correlation Scale",tcorr,&
            rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField,"X Correlation Scale",xcorr,&
            rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeGet(pertField,"Y Correlation Scale",ycorr,&
            rc=status)
       call LIS_verify(status)
       
       if(LIS_rc%obs_ngrid(kk).gt.0) then 
          call ESMF_AttributeGet(pertField,"Standard Deviation",std,&
               rc=status)
          call LIS_verify(status)
       endif


       call ESMF_AttributeGet(pertField,"Cross Correlation Strength",&
            ccorr, rc=status)
       call LIS_verify(status)

       call LIS_gather_1dgrid_to_2dgrid_obs(n,kk,global_std, std)
    ! ----------------------------------------
    !
    ! copy inputs into structure
    
       if(LIS_masterproc) then 
          pert_param(i)%descr          = pertobjs(i)
          pert_param(i)%typ            = typ
          pert_param(i)%std_normal_max = std_normal_max
          pert_param(i)%zeromean       = zeromean
          pert_param(i)%tcorr          = tcorr
          pert_param(i)%xcorr          = xcorr
          pert_param(i)%ycorr          = ycorr
       
          pert_param(i)%std            = global_std

          
          do j=1,objcount
             pert_param(i)%ccorr(j,:,:)  = ccorr(j)
          enddo
       endif
       deallocate(global_std)
    enddo
    ! -------------------------------------------------------------
    !
    ! set mean and (if needed) modify standard deviation according to 'typ'
    ! (additive or multiplicative perturbations)    

    if(LIS_masterproc) then 
       do k=1,N_pert
          
          select case (pert_param(k)%typ)
             
          case (0)         ! additive (mean=0, std as above)
             
             pert_param(k)%mean  = 0.
             
          case (1)         ! multiplicative and lognormal (mean=1)
             
             do i=1,N_x
                do j=1,N_y
                   
                   tmpreal = pert_param(k)%std(i,j) 
                   
                   tmpreal = log( 1. + tmpreal**2)
                   
                   pert_param(k)%mean(i,j) = - .5*tmpreal
                   pert_param(k)%std(i,j)  = sqrt(tmpreal)
                   
                end do
             end do
             
          case default
             
             write (*,*) 'assemble_pert_param(): encountered unknown'
             write (*,*)  k, pert_param(k)%typ
             write (*,*) 'type of error, stopping...'
             stop
          
          end select
       enddo        

    ! -------------------------------------------------------------
    !
    ! compute sqrt of correlation matrix for each grid point
    
       do i=1,N_x
          do j=1,N_y

          ! extract local correlation matrix for grid point (i,j)
          
             do k=1,N_pert
                do l=1,N_pert
                   
                   tmpmat1(k,l) = pert_param(k)%ccorr(l,i,j)
                   
                end do
             end do
          
          ! compute sqrt of local correlation matrix
          
             call get_sqrt_corr_matrix( N_pert, tmpmat1, tmpmat2 )
          
          ! overwrite cross-correlations in forcepert_param with square 
          ! root of cross-correlation matrix
          
             do k=1,N_pert
                do l=1,N_pert
                
                   pert_param(k)%ccorr(l,i,j) = tmpmat2(k,l)
                   
                end do
             end do
             
          end do
       end do
    endif

    deallocate(pertobjs)
    deallocate(ccorr)

  end subroutine assemble_pert_param_global_obs

!BOP
! 
! !ROUTINE: get_sqrt_corr_matrix
! \label{get_sqrt_corr_matrix}
!
! !INTERFACE:
  subroutine get_sqrt_corr_matrix( N, A, S )
! !USES:    
    use nr_jacobi

    implicit none
! !ARGUMENTS:     
    integer, intent(in)                  :: N
    real,    intent(in),  dimension(N,N) :: A
    real,    intent(out), dimension(N,N) :: S
! !DESCRIPTION:
!    
! get sqrt S of real, symmetric (correlation) matrix A
!
! NOTE: there is no check that A is indeed symmetric!
!
! A = S*transpose(S)
!
! reichle, 7 Jun 2005
!
!EOP

    ! local
    
    integer :: j

    real, dimension(N) :: D
    
    ! ------------------------------------------------------------------
    
    call jacobi(A,N,D,S)
        
    do j=1,N
       
       if (D(j)<0.) then
          
          write (*,*) 'get_sqrt_corr_matrix(): negative eigenvalue found'
          write (*,*) '                        invalid correlation matrix'
          write (*,*) 'STOPPING.'
          stop
          
       end if

       S(:,j) = S(:,j)*sqrt(D(j))
       
    end do
    
  end subroutine get_sqrt_corr_matrix

!BOP
! 
! !ROUTINE: get_init_Pert_rseed
! \label{get_init_Pert_rseed}
! 
! !INTERFACE: 
  subroutine get_init_Pert_rseed( N_ens, N_domain, gnc, gnr, &
       ens_id, domain_id, &
       init_Pert_rseed )

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: N_ens, N_domain, gnc, gnr
    integer, dimension(N_ens),    intent(in) :: ens_id
    integer, intent(in) :: domain_id
    integer, dimension(N_ens,N_domain), intent(out) :: init_Pert_rseed

! !DESCRIPTION:    
! get initial random seed "init\_Pert\_rseed" for initializing 
! Pert\_rseed within get\_pert()
!
! A different random seed is necessary for each ensemble member.
! In parallel applications, the random seed must also be different
! for each subdomain.
!
! This subroutine is meant as a sample for how the initial Pert\_rseed 
! can be set.  
! In this example, ens\_id and domain\_id are meant to be 
! nonnegative small integers.
!
!EOP    
    ! --------------------------------------------------
    !
    ! local variables
    
    integer :: m, n, i, j
    
    character(200) :: tmpstring
    
    ! -------------------------------------------------------------
    
    do m=1,N_domain
       do n=1,N_ens

          ! EMK Replaced equation for picking unique negative number for
          ! each ensemble member of each global grid location.
          init_Pert_rseed(n,m) = (domain_id) + &
               ( ens_id(n)*gnc*gnr   ) + &
               ( (m-1)*gnc*gnr*N_ens )
          init_Pert_rseed(n,m) = -1 * init_Pert_rseed(n,m)

       end do
    end do
    
    ! make sure init_Pert_rseed is negative and no two numbers are the same
    
    do m=1,N_domain
       do n=1,N_ens
          
          if (init_Pert_rseed(n,m)>=0) then
             
             tmpstring = 'get_init_Pert_rseed(): found nonnegative'
             tmpstring = trim(tmpstring) // ' component of init_Pert_rseed'
             tmpstring = trim(tmpstring) // ' - STOPPING.'
             write (*,*) trim(tmpString)
             write (*,*) n,m, init_pert_rseed(n,m)
             stop
             
          end if
          
          do i=m+1,N_domain
             do j=n+1,N_ens
                
                if (init_Pert_rseed(j,i)==init_Pert_rseed(n,m)) then
                   
                   tmpstring = 'get_init_Pert_rseed(): found identical'
                   tmpstring = trim(tmpstring) // ' components'
                   tmpstring = trim(tmpstring) // ' of init_Pert_rseed'
                   tmpstring = trim(tmpstring) // ' - STOPPING.'
                   write (*,*) tmpstring
                   stop
                   
                end if
                
             end do
          end do
          
       end do
    end do
    
  end subroutine get_init_Pert_rseed
!BOP
! 
! !ROUTINE: get_pert
!  \label{get_pert}
!
! !INTERFACE:    
  subroutine get_pert(                                &
       N_pert, N_ens, N_x, N_y,                       &
       N_x_fft, N_y_fft,                       &
       dx, dy, dtstep,                                &
       std,                                           & 
       pert_param,                                    &
       Pert_rseed,                                    &
       Pert_ntrmdt,                                   &
       Pert,                                          &
       initialize_rseed,                              &
       initialize_ntrmdt               )
    
    implicit none
! !ARGUMENTS:
    integer, intent(in) :: N_pert   ! # different perturbations
    
    integer, intent(in) :: N_ens  ! # ensemble members
    
    integer, intent(in) :: N_x    ! # grid cells in longitudinal direction
    integer, intent(in) :: N_y    ! # grid cells in latitudinal direction

    integer, intent(in) :: N_x_fft
    integer, intent(in) :: N_y_fft
    
    real, intent(in) :: dx   ! grid cell size in longitudinal direction [deg]
    real, intent(in) :: dy   ! grid cell size in latitudinal direction [deg]
    
    real, intent(in) :: dtstep        ! perturbation time step in seconds
   
    real, intent(in) :: std(N_pert, N_x, N_y)         ! temporally varying standard deviation
    ! Parameter structure for perturbations (see type definition for details).
    
    type(pert_param_type), dimension(:), pointer :: pert_param
    
    ! Pert_ntrmdt are intermediate perturbation fields
    ! that need to be remembered between calls to this subroutine.
    ! In essence, they store N_pert mutually uncorrelated
    ! perturbation fields of standard-normal distribution.
    ! Pert_rseed is the random seed for the generation of
    ! Pert_ntrmdt and is treated similarly to a prognostic variable.
    ! Each ensemble member has its own random seed.

    integer, dimension(NRANDSEED,N_x_fft, N_y_fft, N_ens), intent(inout) :: Pert_rseed
    
    real, dimension(N_pert,N_x,N_y,N_ens), intent(inout) :: Pert_ntrmdt
    
    ! Pert are N_pert cross-correlated perturbation
    ! fields that are rotated and scaled versions of Pert_ntrmdt
    ! so that Pert has the the mean values, standard deviations and
    ! cross-correlations specified in pert_param.  
    ! The distribution is lognormal for multiplicative perturbations.  
    ! Pert should be used as follows for field F (eg. large-scale
    ! precip, convective precip, lw radiation, ...)
    !
    ! Fprime = F+Pert   for additive perturbations
    ! Fprime = F*Pert   for multiplicative perturbations
    !
    ! Note that this subroutine does NOT ensure physically meaningful
    ! perturbed fields.  This is best done outside this subroutine
    ! after the perturbations have been applied.    
    
    real, dimension(N_pert,N_x,N_y,N_ens), intent(out) :: Pert

    ! If initialize_rseed==.true., set initial random seed vector.
    ! If initialize_rseed==.true., the first row of Pert_rseed must be
    ! filled with a different negative integer for each ensemble member.
    ! Note that in parallel applications initial Pert_rseed must also differ
    ! for each subdomain.  See sample subroutine get_init_Pert_rseed().
    !
    ! If initialize_ntrmdt==.true., generate initial Pert_ntrmdt (must be 
    ! allocated!!!).  
    !
    ! Note that when get_pert() is used for generating independent 
    ! perturbations to forcings and prognostic variables, only one
    ! common Pert_rseed should be used, so one of the "ntrmdt" fields
    ! must be initialized without initializing "rseed" again.
    
    ! If initialize_*==.false. or absent, Pert_rseed and 
    ! Pert_ntrmdt must be what was obtained as output from the last call 
    ! to get_pert().  There is only one exception: if there are no temporal 
    ! correlations, it is not necessary to remember Pert_ntrmdt.
    
    logical, intent(in), optional :: initialize_rseed
    logical, intent(in), optional :: initialize_ntrmdt

! !DESCRIPTION:
! get perturbations
!
! reichle, 11 Feb 2005
!
! N\_pert is the number of *perturbation* fields and is not 
! necessarily equal to the number of forcing fields.
! E.g. generate N\_pert=3 perturbations for precip, 
! shortwave radiation, and longwave radiation. These can then be
! applied to N\_force forcing fields, possibly various precip fields 
! (incl large-scale \& convective precip and snow) and to radiation 
! fields.
!EOP    
    
    ! --------------------------------------------------
    !
    ! local variables
    
    type(grid_def_type) :: loc_grid
    
    integer :: i, j, m, mm, n
    
    real, dimension(N_x,N_y)  :: tmp_grid
    
    real :: tmpreal
    
    logical :: init_rseed, init_ntrmdt
   
    
    do m=1,N_pert

       pert_param(m)%std(:,:) = std(m,:,:)
       select case (pert_param(m)%typ)

       case (0)         ! additive (mean=0, std as above)
          
          pert_param(m)%mean  = 0.

       case (1)         ! multiplicative and lognormal (mean=1)
          
          do i=1,N_x
             do j=1,N_y
                tmpreal = pert_param(m)%std(i,j) 
                
                tmpreal = log( 1. + tmpreal**2)
                
                pert_param(m)%mean(i,j) = - .5*tmpreal
                pert_param(m)%std(i,j)  = sqrt(tmpreal)

             end do
          end do
          
       case default
          
          write (*,*) 'get_pert: encountered unknown'
          write (*,*)  m, pert_param(m)%typ
          write (*,*) 'type of error, stopping...'
          stop
          
       end select
       
    end do
    

    ! --------------------------------------------------------------------
    
    ! assemble structure for grid info
    
    loc_grid%N_x = N_x
    loc_grid%N_y = N_y
    loc_grid%N_x_fft = N_x_fft
    loc_grid%N_y_fft = N_y_fft
    loc_grid%llx = -9999.             ! lower left corner not needed
    loc_grid%lly = -9999.             ! lower left corner not needed
    loc_grid%dx  = dx
    loc_grid%dy  = dy
    
    ! ------------------------------------------------------------------
    !
    ! initialize random seed if necessary
    
    init_rseed  = .false.
    init_ntrmdt = .false.
    
    if (present(initialize_rseed))   init_rseed  = initialize_rseed
    if (present(initialize_ntrmdt))  init_ntrmdt = initialize_ntrmdt
    
    if (init_rseed) then
       do i=1,N_x_fft
          do j=1,N_y_fft
             do n=1,N_ens
                
                call init_randseed(Pert_rseed(:,i,j,n))
             enddo
          enddo
       end do
       
    end if
        
    ! ------------------------------------------------------------------
    !
    ! Pert_ntrmdt are standard-normal with desired 
    ! temporal and spatial correlation structure.
    ! Cross-correlations between different fields and scaling to desired
    ! mean and variance is NOT included in Pert_ntrmdt.
    !
    ! Propagate perturbation fields:
    ! On input, Pert_ntrmdt must contain fields from last time step.
    ! (If init_ntrmdt=.true. Pert_ntrmdt is initialized to a 
    ! standard-normal field with the desired spatial correlation structure)
    
    call propagate_pert(                             &
         N_pert, N_ens, loc_grid, dtstep,            &
         Pert_rseed,                                 &
         pert_param,                                 &
         Pert_ntrmdt,                                &
         init_ntrmdt                           )
    
    ! compute diagnostic "Pert" 
    
    ! ensure that ensemble mean perturbation is zero
    ! (must have N_ens>2 for this to make sense).
    !
    ! NOTE: since the sample mean model error varies spatially,
    !       this adjustment slightly changes the *spatial* mean and 
    !       covariance of the model error fields
    !       likely, the benefits of the adjustments for small
    !       ensemble sizes outweigh the disadvantages of altering
    !       the statistical properties
    
    do m=1,N_pert
       
       if ( (pert_param(m)%zeromean) .and. (N_ens>2)) then
          
          do i=1,loc_grid%N_x
             
             call adjust_mean(loc_grid%N_y, N_ens, Pert_ntrmdt(m,i,:,:) )
             
          end do
          
       end if
       
    end do
    
    ! compute rotated fields to get desired cross-correlations between
    ! different fields, then scale to desired mean and std
    
    do m=1,N_pert
       
       do n=1,N_ens
          
          ! rotate to get desired multivariate correlations
          
          do i=1,loc_grid%N_x
             do j=1,loc_grid%N_y
                
                tmp_grid(i,j) = 0.
                
                do mm=1,N_pert
                   
                   tmp_grid(i,j) = tmp_grid(i,j) + & 
                        pert_param(m)%ccorr(mm,i,j) * Pert_ntrmdt(mm,i,j,n)
                   
                end do
                
             end do
          end do
          
          ! scale back freak outliers
          
          call truncate_std_normal( loc_grid%N_x, loc_grid%N_y, &
               pert_param(m)%std_normal_max, tmp_grid )           

          ! scale
          
          do i=1,loc_grid%N_x
             do j=1,loc_grid%N_y
                
                tmpreal = pert_param(m)%mean(i,j) + &
                     pert_param(m)%std(i,j) * tmp_grid(i,j)

                select case (pert_param(m)%typ) 
                   
                case (0)        ! additive
                   
                   Pert(m,i,j,n) = tmpreal
                   
                case (1)        ! multiplicative and lognormal

                   Pert(m,i,j,n) = exp(tmpreal)

                case default
                   
                   write (*,*) 'generate_pert(): encountered unknown'
                   write (*,*) 'typ_pert, stopping...'
                   stop
                   
                end select
             end do
          end do
          

       end do  ! end loop through ensemble members (1:N_ens)
    end do     ! end loop through different fields (1:N_pert)
  end subroutine get_pert
  

!BOP
!
! !ROUTINE: propagate_pert
! \label{propagate_pert}
! 
! !INTERFACE:   
  subroutine propagate_pert(                       &
       N_pert, N_ens, loc_grid, dtstep,            &
       Pert_rseed,                                 &
       pert_param,                                 &
       Pert_ntrmdt,                                &
       initialize              )

! !USES:

    implicit none
    
! !ARGUMENTS:     
    integer, intent(in) :: N_ens, N_pert
    
    type(grid_def_type), intent(in) :: loc_grid   ! local grid definition
    
    type(pert_param_type), dimension(N_pert), intent(in) :: pert_param
    
    real, intent(in) :: dtstep  ! time step of generation of error fields [s]
    
    integer, dimension(NRANDSEED,&
         loc_grid%N_x_fft, loc_grid%N_y_fft,N_ens), intent(inout) :: Pert_rseed
    
    real, dimension(N_pert,loc_grid%N_x,loc_grid%N_y,N_ens), &
         intent(inout) :: Pert_ntrmdt
    
    logical, intent(in) :: initialize   ! switch
 
    
! !DESCRIPTION:
! generate zero-mean, unit-variance (!!) time series 
!  of N\_pert 2d perturbation fields
!
! can also be used just to get a set of 2d random fields (set dtstep 
!  to arbitrary number and "initialize" to .true.)
!
! on input, Pert\_ntrmdt must contain the corresponding 
!  perturbations from the previous time step 
!
! accounts for temporal correlation with AR(1) approach
!  (if pert\_param\%tcorr==0 then error is white in time) 
!
! adapted from off-line EnKF, subroutine propagate\_err() 
!  from NCAT\_59124\_tskin in enkf\_catchment.f90
!
! reichle, 14 Feb 2005
! 
!EOP   
    ! ---------------------------    
    
    ! locals
    
!    integer :: i, j, m, n
    integer :: m, n
    
    real    :: cc, dd
    
    real, dimension(loc_grid%N_x,loc_grid%N_y) :: rfield, rfield2
    
    logical :: white_in_time, white_in_space, stored_field
    
    ! ---------------------------
    
    do m=1,N_pert
       
       ! get parameters for temporal correlation
       
       if ((.not. initialize) .and. (pert_param(m)%tcorr>0.0)) then
          
          white_in_time = .false.
          
          cc = exp( - dtstep / pert_param(m)%tcorr )
          
          dd = sqrt( 1 - cc**2 )
          
       else
          
          cc = 0.
          dd = 1.
          
          white_in_time = .true.
          
       end if
       
       ! find out whether there are spatial correlations
       
       if ( (pert_param(m)%xcorr>0.0) .or. (pert_param(m)%ycorr>0.0) ) then
          
          white_in_space = .false.
       ! disable spatial correlations for now
!          write(LIS_logunit,*) 'Spatial correlations in perturbations are not currently supported..'
!          write(LIS_logunit,*) 'Program stopping....'
!          call LIS_endrun()
          
       else
          
          white_in_space = .true.

       end if
       
       
       ! generate new random fields and propagate AR(1)
       !
       ! Note that rfg2d always produces a pair of random fields!
       !
       ! Use logical variable "stored_field" to figure out whether a second
       ! standard-normal random field is available for next ensemble member.
       ! (in other words, this subroutine is most efficient if N_ens=even 
       ! number, and it is least efficient if N_ens=1)
       
       stored_field = .false.              
       
       do n=1,N_ens
          
          ! generate a random field
          
          if (white_in_space) then 
             
             call generate_white_field(loc_grid%N_x, loc_grid%N_y, &
                  Pert_rseed(:,:,:,n), rfield )
             
          else       ! spatially correlated random fields
             
             ! NOTE: rfg2d_fft() relies on CXML math library (22 Feb 05)
             
             if (.not. stored_field) then
                
                call rfg2d_fft( Pert_rseed(:,:,:,n), 1.,                &
                     loc_grid%N_x, loc_grid%N_y,                    &
                     loc_grid%N_x_fft, loc_grid%N_y_fft,            &
                     loc_grid%dx, loc_grid%dy,                      &
                     pert_param(m)%xcorr, pert_param(m)%ycorr,      &
                     rfield, rfield2 )
                
                stored_field = .true.
                
             else
                
                rfield = rfield2
                
                stored_field = .false.
                
             end if
             
          end if
          
          !! -----------------------------------------------------------
          !!
          !! adjust std of fields to match exactly 1.0
          !!
          !! WARNING: before doing this should check that 
          !!          N_x*dx>>xcorr .and. N_y*dy>>ycorr .and. N_x*N_y>>1 
          !!
          !! Cannot use adjust_std if the field is small relative to 
          !! its spatial correlation scales or if it contains only a few
          !! grid cells (even for white noise in space)
          !!
          !! reichle, 25 Jan 2005
          !!
          !!call adjust_std( loc_grid%N_x, loc_grid%N_y, rfield )
          !!
          !! -----------------------------------------------------------
          
          ! propagate AR(1) 
          
          if (white_in_time) then
             
             Pert_ntrmdt(m,:,:,n) = rfield 
             
          else
             
             Pert_ntrmdt(m,:,:,n) = cc*Pert_ntrmdt(m,:,:,n) + dd*rfield 
          end if
          
       end do
       
    end do
    
  end subroutine propagate_pert
  
!BOP
! 
! !ROUTINE: truncate_std_normal
! \label{truncate_std_normal}
! 
! !INTERFACE: 
  subroutine truncate_std_normal( N_x, N_y, std_normal_max, grid_data )

    implicit none
    
! !ARGUMENTS: 
    integer, intent(in) :: N_x, N_y
    
    real, intent(in) :: std_normal_max
    
    real, dimension(N_x,N_y), intent(inout) :: grid_data
!
! !DESCRIPTION:    
! truncate a realization of standard normal variables
! (scale back freak outliers)
!EOP    
    
    ! local variables
    
    integer :: i,j
    
    ! --------------------------------------------------------
    
    do i=1,N_x
       do j=1,N_y
          
          ! want: -std_normal_max < cat_data < std_normal_max
          
          grid_data(i,j) = &
               sign( min(abs(grid_data(i,j)),std_normal_max), grid_data(i,j) )
          
       end do
    end do
    
  end subroutine truncate_std_normal
  
  ! **************************************************************
!BOP
! 
! !ROUTINE: adjust_mean
! \label{gmaopert_adjust_mean}
! 
! !INTERFACE:   
  subroutine adjust_mean( N_row, N_col, A, M )
    
    implicit none
    
! !ARGUMENTS: 
    integer, intent(in) :: N_row, N_col    
    real, intent(inout), dimension(N_row,N_col)  :: A    
    real, intent(in), optional, dimension(N_row) :: M
!
! !DESCRIPTION:    
! adjust N\_row by N\_col matrix A such that 
! mean over columns for each row is given by the
! corresponding element in vector M of length N\_row
! 
! vector of mean values M is optional input, if not present 
! zero mean is assumed
!EOP    
    ! ----------------------------
    
    ! locals
    
    integer i
    
    real, dimension(N_row) :: correction
    
    ! ------------------------------------------------------------
    
    if (present(M)) then
       correction = M - sum(A,2)/real(N_col) 
    else
       correction = - sum(A,2)/real(N_col) 
    end if
    
    do i=1,N_col
       A(:,i) = A(:,i) + correction
    end do
    
  end subroutine adjust_mean
end module landpert_routines


