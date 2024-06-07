!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module geis_forcingMod
!BOP
! !MODULE: geis_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the global forcing data used in the 
!  Earth Information System (EIS).  The variables are produced 
!  at 5km spatial resolution, and at hourly intervals. 
!
!  NOTE: no spatial interpolation of the code is supported currently
!
! !REVISION HISTORY: 
! 15 Nov 2023: Sujay Kumar; Initial Specification
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GEIS      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: geis_struc
!EOP

  type, public ::  geis_type_dec 
     real          :: ts
     integer       :: ncold, nrold   ! AWIPS 212 dimensions
     character*50  :: geis_filesrc
     character(len=LIS_CONST_PATH_LEN) :: geisdir ! NLDAS-2 Forcing Directory
     real*8        :: geistime1,geistime2

     integer                :: findtime1, findtime2
     integer                :: niter
     integer                :: c_off,r_off
     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type geis_type_dec

  type(geis_type_dec), allocatable :: geis_struc(:)
!EOP
contains
  
!BOP
!
! !ROUTINE: init_GEIS
! \label{init_GEIS}
!
! !INTERFACE:
  subroutine init_GEIS(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit,LIS_endrun
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use map_utils,      only : proj_latlon
    use LIS_forecastMod

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for NLDAS-2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_geis](\ref{readcrd_geis}) \newline
!     reads the runtime options specified for NLDAS-2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[read\_geis\_elev](\ref{read_geis_elev}) \newline
!    reads the native elevation of the NLDAS-2 data to be used
!    for topographic adjustments to the forcing 
!  \end{description}
!EOP
    
    integer :: n
    real    :: lat1, lon1
    
    allocate(geis_struc(LIS_rc%nnest))
    call readcrd_geis()

    do n=1, LIS_rc%nnest
       geis_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc, n, geis_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 7

  ! Set NLDAS-2 grid dimensions and extent information:
    geis_struc(:)%ncold = 7200
    geis_struc(:)%nrold = 3000

    do n=1,LIS_rc%nnest

       allocate(geis_struc(n)%metdata1(1,&
            LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(geis_struc(n)%metdata2(1,&
            LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       
       geis_struc(n)%metdata1 = 0
       geis_struc(n)%metdata2 = 0

       geis_struc(n)%findtime1 = 0 
       geis_struc(n)%findtime2 = 0

       geis_struc(n)%nIter = 1

       lat1 = LIS_domain(n)%lat(1)
       lon1 = LIS_domain(n)%lon(1)
       
       geis_struc(n)%r_off = nint((lat1 + 59.975)/0.05) + 1
       geis_struc(n)%c_off = nint((lon1 + 179.975)/0.05) + 1
       
    enddo
    
  end subroutine init_GEIS
end module geis_forcingMod
