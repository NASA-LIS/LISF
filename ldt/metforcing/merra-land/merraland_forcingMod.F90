!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module merraland_forcingMod
!BOP
! !MODULE: merraland_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the MERRA-Land forcing data. 
!  The data is global 1 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt merraland\_struc}
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the ECMWF data
!  \item[merralandtime1]
!    The nearest, previous 1 hour instance of the incoming 
!    data (as a real time). 
!  \item[merralandtime2]
!    The nearest, next 1 hour instance of the incoming 
!    data (as a real time).
!  \item[merralanddir]
!    Directory containing the input data
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for n. neighbor interpolation. 
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
!
! !USES: 
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_MERRALAND      !defines the native resolution of 
                                !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: merraland_struc

!EOP
  type, public ::  merraland_type_dec
     real      :: ts
     integer   :: nc, nr
     integer   :: mi
     character*80 :: merralanddir   !MERRA-Land Forcing Directory
     real*8    :: merralandtime1,merralandtime2

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
     integer, allocatable   :: n113(:)

     integer                :: findtime1, findtime2
     logical                :: startFlag, dayFlag
     real, allocatable      :: merraforc1(:,:,:), merraforc2(:,:,:)

     integer            :: nvars
     integer            :: uselml

  end type merraland_type_dec
  
  type(merraland_type_dec), allocatable :: merraland_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_MERRALAND
! \label{init_MERRALAND}
!
! !REVISION HISTORY: 
! 12 Oct 2009: Eric Kemp
! 22 Jul 2010: David Mocko, changed to hourly forcing
! 22 Jan 2015: KR Arsenault, added to LDT
! 
! !INTERFACE:
  subroutine init_MERRALAND(findex)

! !USES: 
    use LDT_coreMod
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !AGRUMENTS: 
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for MERRA-Land
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes \ref{interp}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readmerralandcrd](\ref{readmerralandcrd}) \newline
!     reads the runtime options specified for MERRA-Land data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    integer :: n
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    real    :: gridDesci(20)

    allocate(merraland_struc(LDT_rc%nnest))

    write(LDT_logunit,fmt=*)"MSG: Initializing MERRA-Land forcing grid ... "

 !- Read in config file entries:
    call readcrd_merraland()
    
    LDT_rc%met_nf(findex) = 14
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .true.

    merraland_struc%nc = 540
    merraland_struc%nr = 361
    LDT_rc%met_nc(findex) = merraland_struc(1)%nc
    LDT_rc%met_nr(findex) = merraland_struc(1)%nr

 !- GSWP1 Grid description:
    LDT_rc%met_proj(findex)        = "latlon"
    LDT_rc%met_gridDesc(findex,1)  = 0
    LDT_rc%met_gridDesc(findex,2)  = LDT_rc%met_nc(findex)
    LDT_rc%met_gridDesc(findex,3)  = LDT_rc%met_nr(findex)
    LDT_rc%met_gridDesc(findex,4)  = -90.000
    LDT_rc%met_gridDesc(findex,5)  = -180.000
    LDT_rc%met_gridDesc(findex,6)  = 128
    LDT_rc%met_gridDesc(findex,7)  = 90.000
    LDT_rc%met_gridDesc(findex,8)  = 179.333333
    LDT_rc%met_gridDesc(findex,9)  = 0.66666666667
    LDT_rc%met_gridDesc(findex,10) = 0.5
    LDT_rc%met_gridDesc(findex,20) = 0.

    gridDesci(:) = LDT_rc%met_gridDesc(findex,:)
     
    merraland_struc%mi = merraland_struc%nc*merraland_struc%nr

    do n=1,LDT_rc%nnest

       merraland_struc(n)%ts = 3600  !check
       call LDT_update_timestep(LDT_rc, n, merraland_struc(n)%ts)

     ! Setting up weights for Interpolation
       select case ( LDT_rc%met_gridtransform(findex) )
  
         case( "bilinear" )
          allocate(merraland_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(:),&
               merraland_struc(n)%n111,merraland_struc(n)%n121,&
               merraland_struc(n)%n211,merraland_struc(n)%n221,&
               merraland_struc(n)%w111,merraland_struc(n)%w121,&
               merraland_struc(n)%w211,merraland_struc(n)%w221)

         case( "budget-bilinear" )
          allocate(merraland_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merraland_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(:), &
               merraland_struc(n)%n111,merraland_struc(n)%n121,&
               merraland_struc(n)%n211,merraland_struc(n)%n221,&
               merraland_struc(n)%w111,merraland_struc(n)%w121,&
               merraland_struc(n)%w211,merraland_struc(n)%w221)

          allocate(merraland_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merraland_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merraland_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merraland_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merraland_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merraland_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merraland_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merraland_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          call conserv_interp_input(n, gridDesci(:), &
               merraland_struc(n)%n112,merraland_struc(n)%n122,&
               merraland_struc(n)%n212,merraland_struc(n)%n222,&
               merraland_struc(n)%w112,merraland_struc(n)%w122,&
               merraland_struc(n)%w212,merraland_struc(n)%w222)

         case( "neighbor" )
          allocate(merraland_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input(n, gridDesci(:),&
               merraland_struc(n)%n113)

         case default
          write(LDT_logunit,*) 'Interpolation option '// &
               trim(LDT_rc%met_gridtransform(findex))//&
               ' for MERRA-Land forcing is not supported.'
          call LDT_endrun()
       end select

       call LDT_registerAlarm("MERRA-Land forcing alarm",&
            86400.0,86400.0)
       merraland_struc(n)%startFlag = .true. 
       merraland_struc(n)%dayFlag = .true. 

       merraland_struc(n)%nvars = 14

       allocate(merraland_struc(n)%merraforc1(&
            merraland_struc(n)%nvars, 24, &
            LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(merraland_struc(n)%merraforc2(&
            merraland_struc(n)%nvars, 24, &
            LDT_rc%lnc(n)*LDT_rc%lnr(n)))

       merraland_struc(n)%merraforc1 = LDT_rc%udef
       merraland_struc(n)%merraforc2 = LDT_rc%udef

    enddo

  end subroutine init_MERRALAND

end module merraland_forcingMod
 
