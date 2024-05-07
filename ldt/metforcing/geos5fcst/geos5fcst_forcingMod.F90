!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module geos5fcst_forcingMod
!BOP
! !MODULE: geos5fcst_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of GEOS5 ensemble forecast data used as forcing
!  within LDT. The input spatial resolution of the data is at 1.25x1 deg. 
!  This plugin supports either the use of any chosen ensemble member within the
!  GEOS5 forecast output or the full ensemble. 
!
! !USES: 
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GEOS5FCST      !defines the native resolution of 
                                !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: geos5fcst_struc
!EOP

  type, public ::  geos5fcst_type_dec 

     real               :: ts
     integer            :: nc, nr
     integer            :: max_ens_members
     character(len=LDT_CONST_PATH_LEN) :: geos5fcstdir
     real*8             :: fcsttime1,fcsttime2
     
     real               :: gridDesc(50)
     integer            :: mi

     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

     integer            :: findtime1, findtime2
  end type geos5fcst_type_dec

  type(geos5fcst_type_dec), allocatable :: geos5fcst_struc(:)
  logical, public :: geos5fcst_initialized = .false.
!EOP
contains
  
!BOP
!
! !ROUTINE: init_GEOS5FCST
! \label{init_GEOS5FCST}
!
! !REVISION HISTORY: 
! 7 Mar 2013: Sujay Kumar, initial specification
! 
! !INTERFACE:
  subroutine init_GEOS5FCST(findex)
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_domain
    use LDT_timeMgrMod, only : LDT_update_timestep
    use LDT_logMod,     only : LDT_logunit,LDT_endrun
!    use LDT_spatialDownscalingMod, only : LDT_init_pcpclimo_native

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GEOS5 
!  forecast data (which is assumed to be at 1.25x1 deg). 
!  The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}). Currently only bilinear interpolation 
!  is supported for this forcing plugin. 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_geos5fcst](\ref{readcrd_geos5fcst}) \newline
!     reads the runtime options specified for GEOS5 forecast data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!  \end{description}
!EOP
    
    integer :: n

    LDT_rc%met_nf(findex) = 14
    
    if(.not.geos5fcst_initialized) then 
       geos5fcst_initialized = .true. 

       allocate(geos5fcst_struc(LDT_rc%nnest))

       geos5fcst_struc(:)%nc = 288
       geos5fcst_struc(:)%nr = 181

       call readcrd_geos5fcst()
       
       do n=1, LDT_rc%nnest
          LDT_rc%met_nensem(findex) = geos5fcst_struc(n)%max_ens_members
          geos5fcst_struc(n)%ts = 3600
          call LDT_update_timestep(LDT_rc, n, geos5fcst_struc(n)%ts)
          allocate(geos5fcst_struc(n)%metdata1(LDT_rc%met_nf(findex),&
               geos5fcst_struc(n)%max_ens_members,&
               LDT_rc%ngrid(n)))
          allocate(geos5fcst_struc(n)%metdata2(LDT_rc%met_nf(findex),&
               geos5fcst_struc(n)%max_ens_members,&
               LDT_rc%ngrid(n)))
       enddo
       do n=1,LDT_rc%nnest
     
          geos5fcst_struc(n)%gridDesc  = 0        
          geos5fcst_struc(n)%findtime1 = 0 
          geos5fcst_struc(n)%findtime2 = 0 

          geos5fcst_struc(n)%gridDesc(1) = 0
          geos5fcst_struc(n)%gridDesc(2) = geos5fcst_struc(n)%nc
          geos5fcst_struc(n)%gridDesc(3) = geos5fcst_struc(n)%nr
          geos5fcst_struc(n)%gridDesc(4) = -90.0
          geos5fcst_struc(n)%gridDesc(5) = -180.0
          geos5fcst_struc(n)%gridDesc(6) = 128
          geos5fcst_struc(n)%gridDesc(7) = 90.0
          geos5fcst_struc(n)%gridDesc(8) = 178.75
          geos5fcst_struc(n)%gridDesc(9) =  1.25
          geos5fcst_struc(n)%gridDesc(10) = 1.0
          geos5fcst_struc(n)%gridDesc(20) = 64

          geos5fcst_struc(n)%mi = geos5fcst_struc(n)%nc*geos5fcst_struc(n)%nr
          !Setting up weights for Interpolation
          if ( LDT_rc%met_gridtransform(findex).eq."bilinear" ) then 

             allocate(geos5fcst_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
             allocate(geos5fcst_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
             allocate(geos5fcst_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
             allocate(geos5fcst_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
             allocate(geos5fcst_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
             allocate(geos5fcst_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
             allocate(geos5fcst_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
             allocate(geos5fcst_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

             call bilinear_interp_input(n,geos5fcst_struc(n)%gridDesc(:),&
                  geos5fcst_struc(n)%n111,geos5fcst_struc(n)%n121,&
                  geos5fcst_struc(n)%n211,&
                  geos5fcst_struc(n)%n221,geos5fcst_struc(n)%w111,&
                  geos5fcst_struc(n)%w121,&
                  geos5fcst_struc(n)%w211,geos5fcst_struc(n)%w221)
          else
             write(LDT_logunit,*) &
                  '[ERR] Interpolation option not specified for GEOS5FCST'
             write(LDT_logunit,*) '[ERR] Program stopping...'
             call LDT_endrun()
          endif
       enddo
    end if
  end subroutine init_GEOS5FCST
end module geos5fcst_forcingMod
