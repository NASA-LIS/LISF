!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! SUBROUTINE: read_BarlageNative_mxsnoalb
!
! REVISION HISTORY:
! 11 Sep 2017  Initial version....Eric Kemp, SSAI/NASA GSFC
!
! DESCRIPTION:
! Source code for reading and interpolating global maximum snow albedo values 
! at 0.05 deg resolution produced by Mike Barlage (NCAR) from MODIS data.  
! Requires HDF4.
!
! REFERENCE:
! Barlage, M, X Zeng, H Wei, and K E Mitchell, 2005:  A global 0.05 deg maximum
!   albedo dataset for snow-covered land based on MODIS observations.  
!   Geophysical Research Letters, 32, L17405, doi:10.1029/2005GL022881.
!
!------------------------------------------------------------------------------

#include "LDT_misc.h"

subroutine read_BarlageNative_mxsnoalb(n,albedo)

   ! Imports
   use LDT_coreMod,       only : LDT_rc
   use LDT_logMod,        only : LDT_logunit, LDT_endrun

   ! Defaults
   implicit none

   ! Arguments
   integer, intent(in) :: n
   real, intent(inout) :: albedo(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

   ! Local variables
   real, allocatable :: albedo_native(:,:)
   integer, parameter :: X_LENGTH = 7200
   integer, parameter :: Y_LENGTH = 3600

   albedo(:,:,:) = LDT_rc%udef

#if (defined USE_HDF4)

   call get_albedo(n,albedo_native) ! albedo_native allocated in subroutine
   call flip_albedo(albedo_native)  ! Flip albedo so it starts from south pole
   call interp_albedo(n,albedo_native,albedo)
   deallocate(albedo_native)
   write(LDT_logunit,*)'[INFO] Done processing Barlage max snow albedo file'

#else

   write(LDT_logunit,*) &
        '[ERR] HDF4 required to process Barlage max snow albedo!'
   write(LDT_logunit,*) &
        'Reconfigure LDT with HDF4 support, recompile, and run again!'
   call LDT_endrun

#endif

   return

contains

   ! Internal subroutine to pull albedo from HDF4 file
   subroutine get_albedo(n,albedo_native)

      ! Imports
      use LDT_albedoMod,     only : LDT_albedo_struc
      use LDT_logMod,        only : LDT_logunit, LDT_endrun
      
      ! Defaults
      implicit none

      ! Includes
#if (defined USE_HDF4)
#include "hdf.f90"
#endif

      ! Arguments
      integer, intent(in) :: n
      real, allocatable, intent(out) :: albedo_native(:,:)

#if (defined USE_HDF4)

      ! Local variables
      integer :: sd_id, sds_id, sds_index, status, dim_index, n_attrs
      integer :: start(2), edges(2), stride(2)
      integer, parameter :: MAX_VAR_DIMS = 32 ! Not in hdf.f90 for some reason
      integer :: dim_sizes(MAX_VAR_DIMS)
      integer :: rank, data_type
      integer, parameter :: MAX_NC_NAME = 256 ! Not in hdf.f90 for some reason
      character(len=MAX_NC_NAME) :: name

      ! HDF4 functions
      integer, external :: sfstart, sfn2index, sfselect, sfginfo, &
           sfrdata, sfendacc, sfend

      ! Open the file
      sd_id = sfstart(trim(LDT_albedo_struc(n)%mxsnoalbfile), DFACC_READ)
      if (sd_id == FAIL) then
         write(LDT_logunit,*)&
              '[ERR] HDF4 unable to open Barlage MaxSnowAlb map ',&
              trim(LDT_albedo_struc(n)%mxsnoalbfile)
         write(LDT_logunit,*) "[ERR] Program stopping ..."
         call LDT_endrun
      end if
      
      ! Find the index of the albedo variable in the HDF file
      sds_index = sfn2index(sd_id,"MAX_SNOW_ALBEDO")
      if (sds_index == FAIL) then
         write(LDT_logunit,*) &
              '[ERR] HDF4 unable to find MAX_SNOW_ALBEDO in file ',&
              trim(LDT_albedo_struc(n)%mxsnoalbfile)
         write(LDT_logunit,*) "[ERR] Program stopping ..."
         call LDT_endrun
      end if
      sds_id = sfselect(sd_id, sds_index)
      if (sds_id == FAIL) then
         write(LDT_logunit,*) &
              '[ERR] HDF4 unable to find MAX_SNOW_ALBEDO in file ',&
              trim(LDT_albedo_struc(n)%mxsnoalbfile)
         write(LDT_logunit,*) "[ERR] Program stopping ..."
         call LDT_endrun
      end if
      
      ! Sanity check the rank, dimensions, and type
      status = sfginfo(sds_id, name, rank, dim_sizes, data_type, n_attrs)
      if (status == FAIL) then
         write(LDT_logunit,*) &
              '[ERR] HDF4 unable to find info on MAX_SNOW_ALBEDO in file ',&
              trim(LDT_albedo_struc(n)%mxsnoalbfile)
         write(LDT_logunit,*) "[ERR] Program stopping ..."
         call LDT_endrun
      end if
      if ( rank .ne. 2 .or. &
           dim_sizes(1) .ne. X_LENGTH .or. &
           dim_sizes(2) .ne. Y_LENGTH .or. &
           (data_type .ne. DFNT_FLOAT32 .and. &
            data_type .ne. DFNT_NFLOAT32)) then
         write(LDT_logunit,*) &
              '[ERR] MAX_SNOW_ALBEDO characteristics do not match expections!'
         write(LDT_logunit,*)' Rank is ',rank,', expected 2'
         write(LDT_logunit,*)' First dimension is ',dim_sizes(1),&
              ', expected ',X_LENGTH
         write(LDT_logunit,*)' Second dimension is ',dim_sizes(2),&
              ', expected ',Y_LENGTH
         write(LDT_logunit,*)' Data type is ',data_type,&
              ', expected ',DFNT_FLOAT32,' or ',DFNT_NFLOAT32
         write(LDT_logunit,*) "[ERR] Program stopping ..."
         call LDT_endrun
      end if

      ! Now read the albedo field
      start(1) = 0
      start(2) = 0
      edges(1) = X_LENGTH
      edges(2) = Y_LENGTH
      stride(1) = 1
      stride(2) = 1
      allocate(albedo_native(X_LENGTH,Y_LENGTH))
      status = sfrdata(sds_id,start,stride,edges,albedo_native)
      if (status == FAIL) then
         write(LDT_logunit,*) &
              '[ERR] HDF4 unable to read MAX_SNOW_ALBEDO in file ',&
              trim(LDT_albedo_struc(n)%mxsnoalbfile)
         write(LDT_logunit,*) "[ERR] Program stopping ..."
         call LDT_endrun
      end if
      
      ! Terminate access to the dataset
      status = sfendacc(sds_id)
      if (status == FAIL) then
         write(LDT_logunit,*) &
              '[ERR] HDF4 unable to close dataset in file ',&
              trim(LDT_albedo_struc(n)%mxsnoalbfile)
         write(LDT_logunit,*) "[ERR] Program stopping ..."
         call LDT_endrun
      end if
      
      ! Close the file
      status = sfend(sd_id)
      if (status == FAIL) then
         write(LDT_logunit,*) &
              '[ERR] HDF4 unable to close file ',&
              trim(LDT_albedo_struc(n)%mxsnoalbfile)
         write(LDT_logunit,*) "[ERR] Program stopping ..."
         call LDT_endrun
      end if

#endif
      
   end subroutine get_albedo

   ! Internal subroutine to flip data south-to-north
   subroutine flip_albedo(albedo_native)
      
      ! Defaults
      implicit none
      
      ! Arguments
      real, intent(inout) :: albedo_native(X_LENGTH,Y_LENGTH)

      ! Local variables
      real, allocatable :: albedo_tmp(:,:)
      integer :: r,c

      allocate(albedo_tmp(X_LENGTH,Y_LENGTH))
      do r = 1, Y_LENGTH
         do c = 1, X_LENGTH
            albedo_tmp(c,r) = albedo_native(c,Y_LENGTH-r+1)
         end do ! c
      end do ! r
      albedo_native(:,:) = albedo_tmp(:,:)
      deallocate(albedo_tmp)

   end subroutine flip_albedo

   ! Internal subroutine to handle remapping albedo
   subroutine interp_albedo(n,albedo_native,albedo)

      ! Imports
      use LDT_albedoMod,     only : LDT_albedo_struc
      use LDT_coreMod,       only : LDT_rc
      use LDT_fileIOMod,     only : LDT_transform_paramgrid
      use LDT_logMod,        only : LDT_logunit, LDT_endrun

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: n
      real, intent(in) :: albedo_native(X_LENGTH,Y_LENGTH)
      real, intent(inout) :: albedo(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

      ! Local variables     
      real :: param_gridDesc(20)
      integer   :: mi
      integer   :: mo
      real,    allocatable :: gi1(:) 
      logical*1,allocatable:: li1(:) 
      real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) 
      logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) 
      integer :: r,c,i

      ! Set parameter grid array inputs.
      ! Assumes albedo has been flipped so it starts from the south pole.
      ! FIXME -- Confirm these settings with Mike Barlage
      param_gridDesc = 0
      param_gridDesc(1)  =    0.      ! Latlon
      param_gridDesc(2)  = X_LENGTH   ! ncols
      param_gridDesc(3)  = Y_LENGTH   ! nrows
      param_gridDesc(4)  =  -89.95    ! LL lat
      param_gridDesc(5)  = -180.00    ! LL lon
      param_gridDesc(6)  =  128
      param_gridDesc(7)  =   89.95    ! UR lat
      param_gridDesc(8)  =  179.95    ! UR lon
      param_gridDesc(9)  =    0.05    ! dy
      param_gridDesc(10) =    0.05    ! dx
      param_gridDesc(20) =   64

      ! Sanity check grid transform
      select case ( LDT_albedo_struc(n)%mxsnoalb_gridtransform )
      case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )
         write(LDT_logunit,*) "[INFO] Reading Barlage Max Snow Albedo "
      case default
         write(LDT_logunit,*) &
              "[ERR] The spatial transform option selected for Barlage"
         write(LDT_logunit,*) &
              "   max snow albedo file is not recognized nor recommended."
         write(LDT_logunit,*) &
              "   Please select: "
         write(LDT_logunit,*) &
              "  ==  none, neighbor, average, bilinear, budget-bilinear. "
         write(LDT_logunit,*) "Program stopping ..."
         call LDT_endrun
      end select

      !- Enter spatial downscaling/upscaling options to bring the max snow alb
      !  domain to the LIS-run domain ...
      mi = X_LENGTH*Y_LENGTH
      allocate( gi1(mi), li1(mi) )
      gi1(:) = LDT_rc%udef
      li1(:) = .false.
      mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
      lo1(:) = .false.
    
      !- Assign 2-D array to 1-D for aggregation routines:
      i = 0
      do r = 1, Y_LENGTH
         do c = 1, X_LENGTH  
            i = i + 1
            gi1(i) = albedo_native(c,r)
            if( gi1(i) > 0. )  li1(i) = .true.   ! Exclude ocean points
         enddo
      enddo

      !- Transform parameter from original grid to LIS output grid:
      call LDT_transform_paramgrid(n, &
           LDT_albedo_struc(n)%mxsnoalb_gridtransform, &
           param_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )

      ! - Convert 1D to 2D grid output arrays:
      i = 0
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            i = i + 1
            if( go1(i) <= 0.) then   ! Assign Ocean/water values as undefined
               albedo(c,r,1) = LDT_rc%udef
            else
               albedo(c,r,1) = go1(i)
            end if
         enddo
      enddo
      deallocate( li1, gi1 )

   end subroutine interp_albedo

end subroutine read_BarlageNative_mxsnoalb
