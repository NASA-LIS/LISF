!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module cmorph_forcingMod
!BOP
! !MODULE: cmorph_forcingMod
! 
! !REVISION HISTORY:
! 05Jan2006: Yudong Tian; start from Jon G.'s code for older versions of LDT. 
!
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  Climate Prediction Center (CPC)'s MORPHing technique 
!  (CMORPH). CMORPH products are produced as global fields from 
!  60 N-60S at 8km resolution at 30 minute intervals. 
!
!  The implementation in LDT has the derived data type {\tt cmorph\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[cmordir]
!    Directory containing the input data
!  \item[cmortime]
!    The nearest, hourly instance of the incoming 
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch the input resolution to T170
!  \item[mi]
!    Number of points in the input grid
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for conservative interpolation.
!  \end{description}
!
! !USES: 
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_cmorph      !defines the native resolution of 
                                   !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: cmorph_struc
!EOP

  type, public :: cmorph_type_dec
     real    :: ts
     integer :: nc, nr           ! AWIPS 212 dimensions
     character(len=LDT_CONST_PATH_LEN) :: cmorphdir  ! CMOR Forcing Directory
     real*8  :: cmorphtime
     real*8  :: griduptime1
     logical :: gridchange1
     integer :: mi
     integer, allocatable :: n112(:,:)
     integer, allocatable :: n122(:,:)
     integer, allocatable :: n212(:,:)
     integer, allocatable :: n222(:,:)
     real, allocatable ::  w112(:,:),w122(:,:)
     real, allocatable ::  w212(:,:),w222(:,:)
  end type cmorph_type_dec

  type(cmorph_type_dec), allocatable :: cmorph_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_cmorph
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 06Jan2006: Yudong Tian; modification for LDTv4.2
! 
! !INTERFACE:
  subroutine init_cmorph(findex)
! !USES: 
   use LDT_coreMod,    only : LDT_rc, LDT_domain
   use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep
   use LDT_logMod,     only : LDT_logunit

   implicit none
   integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for CMORPH
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes \ref{interp}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_cmorph](\ref{readcrd_cmorph}) \newline
!     reads the runtime options specified for CMORPH data
!   \item[LDT\_date2time](\ref{LDT_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP
   integer :: updoy,yr1,mo1,da1,hr1,mn1,ss1
   real    :: upgmt
   integer :: n
   real    :: gridDesci(20)

    allocate(cmorph_struc(LDT_rc%nnest))

    write(LDT_logunit,fmt=*)"MSG: Initializing CMORPH forcing grid ... "

    call readcrd_cmorph()

    LDT_rc%met_nf(findex) = 2
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .false.
    LDT_rc%met_proj(findex)  = "latlon"

    cmorph_struc%nc = 4948
    cmorph_struc%nr = 1649
    LDT_rc%met_nc(findex) = cmorph_struc(1)%nc
    LDT_rc%met_nr(findex) = cmorph_struc(1)%nr

 !- CMORPH Grid description:
    gridDesci(1)  = 0
    gridDesci(2)  = cmorph_struc(1)%nc
    gridDesci(3)  = cmorph_struc(1)%nr
    gridDesci(4)  = -59.963614
    gridDesci(5)  = -179.963622
    gridDesci(6)  = 128
    gridDesci(7)  = 59.963614
    gridDesci(8)  = 179.963622
    gridDesci(9)  = 0.072756669
    gridDesci(10) = 0.072771377
    gridDesci(20) = 64

    LDT_rc%met_gridDesc(findex,:)  = gridDesci(:)

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    do n=1, LDT_rc%nnest
       cmorph_struc(n)%ts = 3600 !check 
       call LDT_update_timestep(LDT_rc, n, cmorph_struc(n)%ts)
    enddo
    
    do n=1,LDT_rc%nnest       

       allocate(cmorph_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(cmorph_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(cmorph_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(cmorph_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(cmorph_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(cmorph_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(cmorph_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       allocate(cmorph_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
       
       call conserv_interp_input(n,gridDesci(:),&
            cmorph_struc(n)%n112,cmorph_struc(n)%n122,cmorph_struc(n)%n212,&
            cmorph_struc(n)%n222,cmorph_struc(n)%w112,cmorph_struc(n)%w122,&
            cmorph_struc(n)%w212,cmorph_struc(n)%w222)
       
       yr1 = 2012     !grid update time
       mo1 = 10
       da1 = 29
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LDT_date2time(cmorph_struc(n)%griduptime1,&
            updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
       cmorph_struc(n)%gridchange1 = .true.
    enddo

  end subroutine init_cmorph

end module cmorph_forcingMod
