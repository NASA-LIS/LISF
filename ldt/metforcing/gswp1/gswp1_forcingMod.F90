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
! !MODULE: gswp1_forcingMod
!
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Soil
!  Wetness Project (GSWP). GSWP forcing variables are produced
!  on a latlon 1degree grid at 3 hour intervals. 
!
!  The implemenatation in LDT has the derived data type {\tt gswp\_struc} that
!  includes the variables that specify the runtime options
!  They are described below: 
!  \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the GSWP data
!  \item[gswptime1]
!    The nearest, previous 3 hour instance of the incoming 
!    data (as a real time). 
!  \item[gswptime2]
!    The nearest, next 3 hour instance of the incoming 
!    data (as a real time).
!  \item[tair]
!    Directory containing the 2m air temperature data
!  \item[qair]
!    Directory containing the 2m specific humidity data
!  \item[psurf]
!    Directory containing the surface pressure data
!  \item[wind]
!    Directory containing the wind data
!  \item[rainf]
!    Directory containing the total precipitation data
!  \item[snowf]
!    Directory containing the total snowfall data
!  \item[swdown]
!    Directory containing the downward shortwave radiation data
!  \item[swdown]
!    Directory containing the downward longwave radiation data
!  \item[mi]
!    Number of points in the input grid
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
!
!EOP
module gswp1_forcingMod

  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_gswp1      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gswp1_struc

  type, public :: gswp1_type_dec
     real          :: ts
     integer       :: mi 
     integer       :: nmif
     integer       :: nc
     integer       :: nr
     character(len=LDT_CONST_PATH_LEN) :: gswp1dir
     real*8        :: gswp1time1
     real*8        :: gswp1time2
  !Suffixes 1 are for bilinear 
     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)
     
  !Suffixes 2 are for conservative 
     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)
     integer            :: findtime1, findtime2
     
  end type gswp1_type_dec

  type(gswp1_type_dec), allocatable :: gswp1_struc(:)

contains

!BOP
!
! !ROUTINE: init_gswp1
! \label{init_gswp1}
!
! !REVISION HISTORY:
! 11 Dec 2003: Sujay Kumar, Initial Specification
!
! !INTERFACE:
  subroutine init_gswp1(findex)
! !USES:
    use LDT_coreMod, only : LDT_rc
    use LDT_timeMgrMod, only : LDT_update_timestep
    use LDT_logMod,  only : LDT_logunit, LDT_endrun

    implicit none
! !AGRUMENTS: 
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GSWP1
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes \ref{interp}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readgswp1crd](\ref{readgswp1crd}) \newline
!     reads the runtime options specified for GSWP1 data
!  \end{description}
!EOP

    integer  :: n 
    real     :: gridDesci(20)

    allocate(gswp1_struc(LDT_rc%nnest))    

   write(LDT_logunit,fmt=*)"MSG: Initializing GSWP-1 forcing grid ... "

!    call readcrd_gswp1()

    LDT_rc%met_nf(findex) = 10
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .true.

    gswp1_struc%nc = 360
    gswp1_struc%nr = 150
    LDT_rc%met_nc(findex) = gswp1_struc(1)%nc
    LDT_rc%met_nr(findex) = gswp1_struc(1)%nr

 !- GSWP1 Grid description:
    LDT_rc%met_proj(findex)        = "latlon"
    LDT_rc%met_gridDesc(findex,1)  = 0
    LDT_rc%met_gridDesc(findex,2)  = gswp1_struc(1)%nc
    LDT_rc%met_gridDesc(findex,3)  = gswp1_struc(1)%nr
    LDT_rc%met_gridDesc(findex,4)  = -59.500
    LDT_rc%met_gridDesc(findex,5)  = -179.500
    LDT_rc%met_gridDesc(findex,6)  = 128
    LDT_rc%met_gridDesc(findex,7)  = 89.500
    LDT_rc%met_gridDesc(findex,8)  = 179.500
    LDT_rc%met_gridDesc(findex,9)  = 1.000
    LDT_rc%met_gridDesc(findex,10) = 1.000
    LDT_rc%met_gridDesc(findex,20) = 0.

    gridDesci(:) = LDT_rc%met_gridDesc(findex,:)

    gswp1_struc%mi = gswp1_struc%nc*gswp1_struc%nr

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

#if 0
    do n=1,LDT_rc%nnest

       gswp1_struc(n)%ts = 3600 
       call LDT_update_timestep(LDT_rc, n, gswp1_struc(n)%ts)

       gswp1_struc(n)%gswp1time1 = 1000.0
       gswp1_struc(n)%gswp1time2 = 0.0

!Setting up weights for Interpolation
       if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 
          allocate(gswp1_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(:),&
               gswp1_struc(n)%n111,gswp1_struc(n)%n121,&
               gswp1_struc(n)%n211,gswp1_struc(n)%n221,&
               gswp1_struc(n)%w111,gswp1_struc(n)%w121,&
               gswp1_struc(n)%w211,gswp1_struc(n)%w221)

       elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
          allocate(gswp1_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gswp1_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(:),&
               gswp1_struc(n)%n111,gswp1_struc(n)%n121,&
               gswp1_struc(n)%n211,gswp1_struc(n)%n221,&
               gswp1_struc(n)%w111,gswp1_struc(n)%w121,&
               gswp1_struc(n)%w211,gswp1_struc(n)%w221)
          allocate(gswp1_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gswp1_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gswp1_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gswp1_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gswp1_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gswp1_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gswp1_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gswp1_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          call conserv_interp_input(n, gridDesci(:), &
               gswp1_struc(n)%n112,gswp1_struc(n)%n122,&
               gswp1_struc(n)%n212,gswp1_struc(n)%n222,&
               gswp1_struc(n)%w112,gswp1_struc(n)%w122,&
               gswp1_struc(n)%w212,gswp1_struc(n)%w222)

       elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then 
          write(LDT_logunit,*) 'Neighbor interpolation is not supported'
          write(LDT_logunit,*) 'for GSWP1 forcing... Program stopping..'
          call LDT_endrun()
       endif

    enddo
#endif
    
  end subroutine init_gswp1
end module gswp1_forcingMod

