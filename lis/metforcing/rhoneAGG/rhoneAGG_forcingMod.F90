!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: rhoneAGG_forcingMod
!
! !DESCRIPTION:
!  Contains routines and variables that define the native domain
!  for RHONEAGG model forcing.
!
!EOP
module rhoneAGG_forcingMod

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_RHONEAGG      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: rhoneAGG_struc

  type, public :: rhoneAGG_type_struc
     real           :: ts
     character*50   :: rhoneAGGdir
     integer        :: mi
     integer        :: ncold
     integer        :: nrold
     real*8         :: rhoneAGGtime1, rhoneAGGtime2
     integer        :: nmif
     integer        :: findtime1, findtime2

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type rhoneAGG_type_struc

  type(rhoneAGG_type_struc), allocatable :: rhoneAGG_struc(:)

contains

!BOP
!
! !ROUTINE: init_RHONEAGG
!
! !REVISION HISTORY:
! 11 Dec 2003: Sujay Kumar, Initial Specification
!
! !INTERFACE:
  subroutine init_RHONEAGG(findex)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
    integer,   intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for RHONEAGG
!  AGG data. 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_rhoneAGG](\ref{readcrd_rhoneAGG}) \newline
!     reads the runtime options specified for RhoneAGGAGG data
!  \end{description}
!EOP
    real    :: gridDesci(50)
    integer :: n 
    allocate(rhoneAGG_struc(LIS_rc%nnest))

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the Rhone-AGG experiment forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    do n=1, LIS_rc%nnest
       rhoneAGG_struc(n)%ts = 3*3600
       call LIS_update_timestep(LIS_rc, n, rhoneAGG_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 9 

    call readcrd_rhoneAGG()

    do n=1,LIS_rc%nnest

       allocate(rhoneAGG_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(rhoneAGG_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       rhoneAGG_struc(n)%metdata1 = 0
       rhoneAGG_struc(n)%metdata2 = 0

       gridDesci(1) = 0
       gridDesci(2) = rhoneAGG_struc(n)%ncold
       gridDesci(3) = rhoneAGG_struc(n)%nrold
       gridDesci(4) = 43.500
       gridDesci(5) =  3.500
       gridDesci(7) = 48.500
       gridDesci(8) =  7.500
       gridDesci(6) = 128
       gridDesci(9) = 1.000
       gridDesci(10) = 1.000
       gridDesci(20) = 255
    
       rhoneAGG_struc(n)%mi = 0!rhoneAGGdrv%ncold*rhoneAGGdrv%nrold
    enddo

  end subroutine init_RHONEAGG
end module rhoneAGG_forcingMod
  
