!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gldas_forcingMod
!BOP
! !MODULE: gldas_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the GLDAS forcing data produced by merging
!  GEOS, AGRMET radiation and CMAP precipitation
!  
!  The forcing data consists of the following meteorological variables. 
!  \begin{description}
!   \item[2m air temperature]
!   \item[2m specific humidity]
!   \item[incident downward shortwave radiation]
!   \item[incident downward longwave radiation]
!   \item[eastward wind]
!   \item[northward wind]
!   \item[surface pressure]
!   \item[rainfall rate]
!   \item[snowfall rate]
!   \item[convective rainfall rate]
!  \end{description}
!
! !USES: 
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_gldas      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gldas_struc
!EOP

  type, public :: gldas_type_dec
     real         :: ts
     integer      :: nc, nr   !AWIPS 212 dimensions
     integer      :: nmif
     character(len=LDT_CONST_PATH_LEN) :: gldasdir   !GLDAS Forcing Directory
     real*8       :: gldastime1,gldastime2
     integer :: mi
  !Suffixes 1 are for bilinear 
  !Suffixes 2 are for conservative 
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
     real, allocatable      ::  w112(:,:),w122(:,:)
     real, allocatable      ::  w212(:,:),w222(:,:)

     integer, allocatable   :: n113(:)

     integer, allocatable   :: smask1(:,:)
     integer, allocatable   :: colfill(:,:)
     integer, allocatable   :: rowfill(:,:)
     logical            :: fillflag1
     integer            :: findtime1, findtime2
  end type gldas_type_dec
  
  type(gldas_type_dec), allocatable :: gldas_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_gldas
!  \label{init_gldas}
!
! !REVISION HISTORY: 
!  19 Sept 2008: Sujay Kumar: Initial Implementation
! 
! !INTERFACE:
  subroutine init_gldas(findex)
! !USES: 
    use LDT_coreMod,    only : LDT_rc
    use LDT_timeMgrMod, only : LDT_update_timestep
    use LDT_logMod,     only : LDT_logunit,LDT_endrun

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GLDAS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes \ref{interp}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readgldascrd](\ref{readgldascrd}) \newline
!     reads the runtime options specified for GLDAS data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[neighbor\_interp\_input](\ref{neighbor_interp_input}) \newline
!    computes the computational weights for neighbor interpolation
!  \end{description}
!EOP

    integer :: n 
    real    :: gridDesci(20)

    allocate(gldas_struc(LDT_rc%nnest))

   write(LDT_logunit,fmt=*)"MSG: Initializing GLDAS forcing grid ... "

!    call readcrd_gldas()

    LDT_rc%met_nf(findex) = 10
    LDT_rc%met_ts(findex) = 3*3600
    LDT_rc%met_zterp(findex) = .true.

    gldas_struc%nc = 181
    gldas_struc%nr = 61
    LDT_rc%met_nc(findex) = gldas_struc(1)%nc
    LDT_rc%met_nr(findex) = gldas_struc(1)%nr

 !- AGRMET Grid description:
    LDT_rc%met_proj(findex)        = "latlon"
    LDT_rc%met_gridDesc(findex,1)  = 0
    LDT_rc%met_gridDesc(findex,2)  = gldas_struc(1)%nc
    LDT_rc%met_gridDesc(findex,3)  = gldas_struc(1)%nr
    LDT_rc%met_gridDesc(findex,4)  = -60.0
    LDT_rc%met_gridDesc(findex,5)  = -180.0
    LDT_rc%met_gridDesc(findex,6)  = 128
    LDT_rc%met_gridDesc(findex,7)  = 90.0
    LDT_rc%met_gridDesc(findex,8)  = 180.0
    LDT_rc%met_gridDesc(findex,9)  = 2.00
    LDT_rc%met_gridDesc(findex,10) = 2.50
    LDT_rc%met_gridDesc(findex,20) = 0

    gridDesci(:) = LDT_rc%met_gridDesc(findex,:)

    gldas_struc%mi = gldas_struc%nc*gldas_struc%nr

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return
    
#if 0
    do n=1,LDT_rc%nnest

       gldas_struc(n)%ts = 3*3600
       call LDT_update_timestep(LDT_rc, n, gldas_struc(n)%ts)

     ! Setting up weights for Interpolation
       if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 
          allocate(gldas_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci, &
               gldas_struc(n)%n111,gldas_struc(n)%n121,&
               gldas_struc(n)%n211,gldas_struc(n)%n221,&
               gldas_struc(n)%w111,gldas_struc(n)%w121,&
               gldas_struc(n)%w211,gldas_struc(n)%w221)

          allocate(gldas_struc(n)%smask1(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%colfill(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%rowfill(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          gldas_struc(n)%fillflag1 = .true. 

       elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
          allocate(gldas_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci,&
               gldas_struc(n)%n111,gldas_struc(n)%n121,&
               gldas_struc(n)%n211,gldas_struc(n)%n221,&
               gldas_struc(n)%w111,gldas_struc(n)%w121,&
               gldas_struc(n)%w211,gldas_struc(n)%w221)
          allocate(gldas_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gldas_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gldas_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gldas_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gldas_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gldas_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gldas_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(gldas_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          call conserv_interp_input(n, gridDesci,&
               gldas_struc(n)%n112,gldas_struc(n)%n122,&
               gldas_struc(n)%n212,gldas_struc(n)%n222,&
               gldas_struc(n)%w112,gldas_struc(n)%w122,&
               gldas_struc(n)%w212,gldas_struc(n)%w222)
          allocate(gldas_struc(n)%smask1(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%colfill(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%rowfill(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          gldas_struc(n)%fillflag1 = .true. 

       elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then 
          allocate(gldas_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input(n, gridDesci, gldas_struc(n)%n113)

          allocate(gldas_struc(n)%smask1(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%colfill(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          allocate(gldas_struc(n)%rowfill(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          gldas_struc(n)%fillflag1 = .true. 
       endif
    enddo
#endif
  end subroutine init_gldas

end module gldas_forcingMod

