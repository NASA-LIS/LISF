!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module COAMPSout_forcingMod
!BOP
! !MODULE: COAMPSout_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data extracted from COAMPS output
!  files.  Here COAMPS output files are consider input data.
!
!  The implementation in LIS has the derived data type {\tt COAMPSout\_struc} that
!  includes the variables that specify the runtime options.
!  They are described below: 
!  \begin{description}
!  \item[nest\_id]
!    Value of the COAMPS nest
!  \item[COAMPSoutdir]
!    Directory containing the input data
!  \item[ts]
!    Frequency in seconds of the forcing data
!  \item[COAMPSouttime1]
!    The nearest, previous 1 hour instance of the incoming 
!    data (as a real time). 
!  \item[COAMPSouttime2]
!    The nearest, next 1 hour instance of the incoming 
!    data (as a real time).
!  \item[findtime1, findtime2]
!    boolean flags to indicate which time is to be read for 
!    temporal interpolation.
!  \end{description}

! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_COAMPSout      !defines the native resolution of 
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: COAMPSout_struc
!EOP

  type, public    :: COAMPSout_type_dec
     integer      :: nest_id
     character(len=LIS_CONST_PATH_LEN) :: COAMPSoutdir
     real         :: ts
     real*8       :: COAMPSouttime1,COAMPSouttime2
     integer      :: findtime1,findtime2

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type COAMPSout_type_dec

  type(COAMPSout_type_dec), allocatable :: COAMPSout_struc(:)
!EOP
contains
  
!BOP
!
! !ROUTINE: init_COAMPSout
! \label{init_COAMPSout}
!
! !REVISION HISTORY: 
! 
! !INTERFACE:
  subroutine init_COAMPSout(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc
    use LIS_timeMgrMod, only : LIS_update_timestep

    implicit none
! !ARGUMENTS:  
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for LIS output
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!  Note that no interpolation is performed for
!  this forcing. The data is expected to be in the same map projection
!  and resolution as that of the current LIS run.
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_COAMPSout](\ref{readcrd_COAMPSout}) \newline
!     reads the runtime options specified for COAMPS output data
!  \end{description}
!EOP
    
    integer :: n
    
    allocate(COAMPSout_struc(LIS_rc%nnest))

    call readcrd_COAMPSout()

    LIS_rc%met_nf(findex) = 11

    do n=1, LIS_rc%nnest

       allocate(COAMPSout_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(COAMPSout_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       COAMPSout_struc(n)%metdata1 = 0
       COAMPSout_struc(n)%metdata2 = 0

       COAMPSout_struc(n)%ts = 60*60 
       call LIS_update_timestep(LIS_rc, n, COAMPSout_struc(n)%ts)
    enddo

    !tmp     <--> LIS_forc%metdata1(1,:)
    !q2      <--> LIS_forc%metdata1(2,:)
    !swd     <--> LIS_forc%metdata1(3,:)
    !lwd     <--> LIS_forc%metdata1(4,:)
    !uwind   <--> LIS_forc%metdata1(5,:)
    !vwind   <--> LIS_forc%metdata1(6,:)
    !psurf   <--> LIS_forc%metdata1(7,:)
    !pcp     <--> LIS_forc%metdata1(8,:)
    ! optional forcing
    !fheight <--> LIS_forc%metdata1(10,:)
    !acond   <--> LIS_forc%metdata1(11,:)

  end subroutine init_COAMPSout
end module COAMPSout_forcingMod
