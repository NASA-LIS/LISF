!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: genEnsFcst_SpatialInterpMod
! \label{genEnsFcst_SpatialInterpMod}
! 
! !DESCRIPTION: 
!  This module contains routines that provide a general 
!   interface with the interp routines, allowing users to
!   easily plug-in interp options and keep account of
!   associated arrays.  
!   
! !REVISION HISTORY: 
!  12Jan2015 -- KR Arsenault;  Initial Specification
! 
! !INTERFACE:
module genEnsFcst_SpatialInterpMod
!
! !USES:

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: genEnsFcst_interp_init      ! initialize interp routine inputs/outputs
  public :: genEnsFcst_interp_data      ! transform input data/field to LIS target grid
  public :: genEnsFcst_interp_finalize  ! deallocate/reinit any interp arrays/inputs
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: sinterp_struc
!EOP

  type spatialinterp_type_dec

   ! Subset parameters for "none" reprojection type:
     integer, allocatable  :: lat_line(:,:)
     integer, allocatable  :: lon_line(:,:)
     real                  :: subset_gridDesc(20)
     integer               :: subset_nc, subset_nr
 
   ! Bilinear pts and wts (n111, also used for averaging):
     integer, allocatable  :: n111(:), n121(:)
     integer, allocatable  :: n211(:), n221(:)
     real,    allocatable  :: w111(:), w121(:)
     real,    allocatable  :: w211(:), w221(:)

   ! Budget-bilinear pts and wts (25 radius points):
     integer, allocatable  :: n112(:,:), n122(:,:)
     integer, allocatable  :: n212(:,:), n222(:,:)
     real,    allocatable  :: w112(:,:), w122(:,:)
     real,    allocatable  :: w212(:,:), w222(:,:)

   ! Neighbor points:
     integer, allocatable  :: n113(:)

  end type spatialinterp_type_dec

  type(spatialinterp_type_dec), allocatable :: sinterp_struc(:)

contains

!BOP
! !ROUTINE: genEnsFcst_interp_init 
!  \label{genEnsFcst_interp_init}
! 
! !INTERFACE: 
!
 subroutine genEnsFcst_interp_init( n, findex, input_gridDesc )
!
! !USES: 
   use LIS_coreMod, only : LIS_rc, LIS_domain
   use LIS_logMod,  only : LIS_logunit, LIS_endrun
   use LIS_gridmappingMod, only : LIS_RunDomainPts

   implicit none
!
! !ARGUMENTS: 
   integer, intent(in) :: n                  ! Nest index
   integer, intent(in) :: findex             ! Forcing index
   real,    intent(in) :: input_gridDesc(50) ! Input forcing griddesc array
!
! !DESCRIPTION: 
!  This routine allocates and generates the points and
!   weights for the user-specified interp or grid transform
!   option to the target LIS output grid.
!EOP  
!
   integer      :: input_nc, input_nr
   integer      :: numinpts
   integer      :: numoutpts
   integer      :: radiuspts
   integer      :: iret
   integer      :: glpnc, glpnr
! _______________________________________________


   allocate( sinterp_struc(LIS_rc%nnest) )

   write(LIS_logunit,*) "[INFO] "
   write(LIS_logunit,*) "[INFO] Reprojection option selected: ",&
         trim(LIS_rc%met_interp(findex))," (interpolation)"

 ! Determine LIS output grid projection type:
   select case( nint(LIS_rc%gridDesc(n,1)) )
    case ( 0 )  ! Lat-long
      input_nc = input_gridDesc(2)
      input_nr = input_gridDesc(3)
    case default
      write(LIS_logunit,*) "[ERR] OTHER OUTPUT GRIDS ARE NOT CURRENTLY "
      write(LIS_logunit,*) "[ERR] SUPPORTED BY THE GENERATED MET FORCING READER."
      write(LIS_logunit,*) "[ERR] Program stopping ..."
      call LIS_endrun()
   end select
  
   numinpts  = input_nc * input_nr
   numoutpts = LIS_rc%lnc(n) * LIS_rc%lnr(n)
   radiuspts = 25
  
!- Determine if reprojection option:
   select case( LIS_rc%met_interp(findex) )

   ! Spatial averaging ("upscaling"):
     case( "average" )  

     ! Lat-lon grid check:
       if( nint(LIS_rc%gridDesc(n,1)) .ne. 0 ) THEN 
          write(LIS_logunit,*) "[ERR] Averaging only supported for"
          write(LIS_logunit,*) "[ERR]  lat/lon run domains. "
          write(LIS_logunit,*) "[ERR] Program stopping ..."
          call LIS_endrun()
       endif

     ! Confirm if LIS run domain resolution really not less than forcing
     ! resolution for lat/lon projection: 
       if( input_gridDesc(10) .gt. LIS_rc%gridDesc(n,10)  ) then
          write(LIS_logunit,*) "[ERR] For 'average' reprojection option, "
          write(LIS_logunit,*) "[ERR] The LIS lat/lon run domain resolution must be "
          write(LIS_logunit,*) "[ERR]  > the generated ensemble forecast forcing resolution."
          write(LIS_logunit,*) "[ERR]  Program stopping ..."
          call LIS_endrun()
       endif

       allocate( sinterp_struc(n)%n111(numinpts) )
       call upscaleByAveraging_input( input_gridDesc(:),&
                   LIS_rc%gridDesc(n,:), numinpts, numoutpts,&
                   sinterp_struc(n)%n111)

   ! Nearest neighbor search ("downscaling"):
     case( "neighbor" )

       allocate( sinterp_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)) )
       call neighbor_interp_input( n, input_gridDesc, &
                     sinterp_struc(n)%n113 )

   ! Bilinear interpolation ("downscaling"):
     case( "bilinear" )
       allocate( sinterp_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)) )
       allocate( sinterp_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)) )
       allocate( sinterp_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)) )
       allocate( sinterp_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)) )
       allocate( sinterp_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)) )
       allocate( sinterp_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)) )
       allocate( sinterp_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)) )
       allocate( sinterp_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)) )

       call bilinear_interp_input(n, input_gridDesc, &
                     sinterp_struc(n)%n111, sinterp_struc(n)%n121, &
                     sinterp_struc(n)%n211, sinterp_struc(n)%n221, &
                     sinterp_struc(n)%w111, sinterp_struc(n)%w121, &
                     sinterp_struc(n)%w211, sinterp_struc(n)%w221 )

     ! Note: Could enter "radiuspts" instead of 25 below ...
       allocate( sinterp_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )

       call conserv_interp_input(n, input_gridDesc, &
                 sinterp_struc(n)%n112, sinterp_struc(n)%n122,&
                 sinterp_struc(n)%n212, sinterp_struc(n)%n222,&
                 sinterp_struc(n)%w112, sinterp_struc(n)%w122,&
                 sinterp_struc(n)%w212, sinterp_struc(n)%w222 )

   ! Budget-bilinear (conservative) interpolation ("downscaling"):
     case( "budget-bilinear" )

     ! Note: Could enter "radiuspts" instead of 25 below ...
       allocate( sinterp_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )
       allocate( sinterp_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25) )

       call conserv_interp_input(n, input_gridDesc, &
                 sinterp_struc(n)%n112, sinterp_struc(n)%n122,&
                 sinterp_struc(n)%n212, sinterp_struc(n)%n222,&
                 sinterp_struc(n)%w112, sinterp_struc(n)%w122,&
                 sinterp_struc(n)%w212, sinterp_struc(n)%w222 )

     case ( "none" )

       write(LIS_logunit,*)"[INFO] NO spatial transform of forcing (only subset may occur)."

       call LIS_RunDomainPts( n, LIS_rc%met_proj(findex), input_gridDesc, &
                glpnc, glpnr, sinterp_struc(n)%subset_nc, sinterp_struc(n)%subset_nr, &
                sinterp_struc(n)%subset_gridDesc, &
                sinterp_struc(n)%lat_line, sinterp_struc(n)%lon_line )

     case default
       write(LIS_logunit,*) " ERR: The reprojection option selected, "&
           //trim(LIS_rc%met_interp(findex))//", is not currently supported."
       write(LIS_logunit,*) " Program stopping ... "
       call LIS_endrun()

   end select

 end subroutine genEnsFcst_interp_init


!BOP
! !ROUTINE: genEnsFcst_interp_data
!  \label{genEnsFcst_interp_data}
! 
! !INTERFACE: 
!
 subroutine genEnsFcst_interp_data( n, findex,  &
               input_nc, input_nr, tindex, input_2dvar, &
               ppt_flag, output_2dvar )
!
! !USES: 
   use LIS_coreMod, only : LIS_rc, LIS_domain
   use LIS_logMod,  only : LIS_logunit, LIS_endrun

   implicit none
!
! !ARGUMENTS: 
   integer, intent(in) :: n             ! Nest index
   integer, intent(in) :: findex        ! Forcing index
   integer, intent(in) :: input_nc      ! Input number of cols
   integer, intent(in) :: input_nr      ! Input number of rows
   logical, intent(in) :: ppt_flag      ! Precip field flag 
!   integer, intent(in) :: ntimes        ! Number of time points in a day
   integer, intent(in) :: tindex        ! Time points in a day
   real,    intent(in) :: input_2dvar(input_nc,input_nr,1)            ! Input Var
   real, intent(inout) :: output_2dvar(LIS_rc%lnc(n),LIS_rc%lnr(n),1) ! Output Var

!
! !DESCRIPTION: 
!  This routine transforms spatially the original 
!   meteorological forcing file to the target 
!   LIS output grid.
!EOP  
!
   integer   :: iret
   integer   :: numinpts
   integer   :: numoutpts
   integer   :: c,r,t
   integer   :: count1,i,j
   real      :: input_1d( input_nc*input_nr )
   real      :: output_1d( LIS_rc%lnc(n)*LIS_rc%lnr(n) )
   logical*1 :: li1d(input_nc*input_nr)
   logical*1 :: lo1d(LIS_rc%lnc(n),LIS_rc%lnr(n))

! __________________________________________________

   output_2dvar = LIS_rc%udef
   li1d = .true.  
   lo1d = .true. 

   numinpts  = input_nc * input_nr
   numoutpts = LIS_rc%lnc(n) * LIS_rc%lnr(n)

!- Loop over Number of hour increments in dataset:
!   do t = 1, ntimes
   t = 1

      ! Create input 1D array:
      count1 = 0
      do r=1, input_nr
         do c=1, input_nc
            count1 = count1 + 1
            input_1d(count1) = input_2dvar(c, r, t) 
         end do
      end do

      ! Define logical mask for input field:
      where ( input_1d == LIS_rc%udef )
         li1d = .false.
      endwhere

      ! Upscale finer scale forcing field to coarser scale run domain:
      select case( LIS_rc%met_interp(findex) )

        case( "average" )   ! Upscaling 
          call upscaleByAveraging(numinpts, numoutpts, &
                      LIS_rc%udef, sinterp_struc(n)%n111, &
                      li1d, input_1d, lo1d, output_1d )

        ! Downscale (or interpolate) or at same resolution as run domain:
        case( "neighbor" )
          call neighbor_interp( LIS_rc%gridDesc(n,:), &
                       li1d, input_1d, lo1d, output_1d,&
                       numinpts, numoutpts,&
                       LIS_domain(n)%lat, LIS_domain(n)%lon,&
                       sinterp_struc(n)%n113, LIS_rc%udef, iret )

        case( "bilinear" )
          ! If precipitation field present, use budget-bilinear:
          if( ppt_flag ) then
             call conserv_interp( LIS_rc%gridDesc(n,:), &
               li1d, input_1d, lo1d, output_1d,&
               numinpts, numoutpts, &
               LIS_domain(n)%lat,  LIS_domain(n)%lon, &
               sinterp_struc(n)%w112, sinterp_struc(n)%w122,&
               sinterp_struc(n)%w212, sinterp_struc(n)%w222,&
               sinterp_struc(n)%n112, sinterp_struc(n)%n122,&
               sinterp_struc(n)%n212, sinterp_struc(n)%n222,&
               LIS_rc%udef, iret)
          else
             call bilinear_interp( LIS_rc%gridDesc(n,:), &
               li1d, input_1d, lo1d, output_1d,&
               numinpts, numoutpts, &
               LIS_domain(n)%lat,  LIS_domain(n)%lon, &
               sinterp_struc(n)%w111, sinterp_struc(n)%w121,&
               sinterp_struc(n)%w211, sinterp_struc(n)%w221,&
               sinterp_struc(n)%n111, sinterp_struc(n)%n121,&
               sinterp_struc(n)%n211, sinterp_struc(n)%n221,&
               LIS_rc%udef, iret )
          endif

        case( "budget-bilinear" )
           call conserv_interp( LIS_rc%gridDesc(n,:), &
             li1d, input_1d, lo1d, output_1d,&
             numinpts, numoutpts, &
             LIS_domain(n)%lat,  LIS_domain(n)%lon, &
             sinterp_struc(n)%w112, sinterp_struc(n)%w122,&
             sinterp_struc(n)%w212, sinterp_struc(n)%w222,&
             sinterp_struc(n)%n112, sinterp_struc(n)%n122,&
             sinterp_struc(n)%n212, sinterp_struc(n)%n222,&
             LIS_rc%udef, iret)

        case( "none" )

        ! If input grid matches the output grid size:
          if( sinterp_struc(n)%subset_nc == input_nc .and. &
              sinterp_struc(n)%subset_nr == input_nr ) then

            output_1d = input_1d

        ! Otherwise, subset for same domain resolutions:
          else
            do r = 1, sinterp_struc(n)%subset_nr
              do c = 1, sinterp_struc(n)%subset_nc
                 output_2dvar(c,r,t) = &
                        input_2dvar(&
                            sinterp_struc(n)%lon_line(c,r),&
                            sinterp_struc(n)%lat_line(c,r),t)
              enddo ! columns
            enddo ! rows
            return
          endif

      end select

!--------------------------------------------------------------------    
! Create 2D array writing out final grid:
!--------------------------------------------------------------------    
      count1 = 0
      do r=1, LIS_rc%lnr(n)
         do c=1, LIS_rc%lnc(n)
            count1 = count1 + 1
            output_2dvar(c,r,t) = output_1d(count1)
         end do
      end do

!   end do    ! End field time loop

 end subroutine genEnsFcst_interp_data


!BOP
! !ROUTINE: genEnsFcst_interp_finalize
!  \label{genEnsFcst_interp_finalize}
! 
! !INTERFACE: 
!
 subroutine genEnsFcst_interp_finalize( findex )
!
! !USES: 
   use LIS_coreMod, only : LIS_rc, LIS_domain
   use LIS_logMod,  only : LIS_logunit, LIS_endrun

   implicit none
   integer :: n
   integer :: findex

   do n = 1, LIS_rc%nnest
      select case( LIS_rc%met_interp(findex) )

       case( "average" )   
         deallocate( sinterp_struc(n)%n111 )

       case( "neighbor" )
         deallocate( sinterp_struc(n)%n113 )

       case( "bilinear" )
         deallocate( sinterp_struc(n)%n111 )
         deallocate( sinterp_struc(n)%n121 )
         deallocate( sinterp_struc(n)%n211 )
         deallocate( sinterp_struc(n)%n221 )
         deallocate( sinterp_struc(n)%w111 )
         deallocate( sinterp_struc(n)%w121 )
         deallocate( sinterp_struc(n)%w211 )
         deallocate( sinterp_struc(n)%w221 )

         deallocate( sinterp_struc(n)%n112 )
         deallocate( sinterp_struc(n)%n122 )
         deallocate( sinterp_struc(n)%n212 )
         deallocate( sinterp_struc(n)%n222 )
         deallocate( sinterp_struc(n)%w112 )
         deallocate( sinterp_struc(n)%w122 )
         deallocate( sinterp_struc(n)%w212 )
         deallocate( sinterp_struc(n)%w222 )

       case( "budget-bilinear" )
         deallocate( sinterp_struc(n)%n112 )
         deallocate( sinterp_struc(n)%n122 )
         deallocate( sinterp_struc(n)%n212 )
         deallocate( sinterp_struc(n)%n222 )
         deallocate( sinterp_struc(n)%w112 )
         deallocate( sinterp_struc(n)%w122 )
         deallocate( sinterp_struc(n)%w212 )
         deallocate( sinterp_struc(n)%w222 )

      end select
   enddo

 end subroutine genEnsFcst_interp_finalize

end module genEnsFcst_SpatialInterpMod
