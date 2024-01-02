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
! SUBROUTINE: read_maskfile_UKMO_CAP_Netcdf
!
! REVISION HISTORY:
! 6 Oct 2017  Initial version....Eric Kemp, SSAI/NASA GSFC
!
! DESCRIPTION:
! Source code for reading and interpolating land mask data in netCDF format
! produced from UKMO_CAP_Netcdf land use data.  The logic in these
! routines borrows heavily from the read_maskfile subroutine.
!
!------------------------------------------------------------------------------

#include "LDT_misc.h"

subroutine read_maskfile_UKMO_CAP_Netcdf(n, vegtype, fgrd, localmask)

   ! Imports
   use LDT_coreMod,   only : LDT_rc
   use LDT_logMod,    only : LDT_logunit, LDT_endrun
   
   ! Defaults
   implicit none
   
   ! Arguments
   integer, intent(in)  :: n
   real,    intent(in)  :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real,    intent(in)  :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
   real,    intent(out) :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))
   
   ! Local variables
   integer :: nlon,nlat
   real, allocatable :: landmask_cap(:,:)
   
   localmask(:,:) = LDT_rc%udef
   
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   
   call get_landmask(n,nlon,nlat,landmask_cap)
   call interp_landmask(n, nlon, nlat, landmask_cap, localmask)
   deallocate(landmask_cap)
   call check_landmask(n,vegtype,fgrd,localmask)
   write(LDT_logunit,*)'[INFO] Done processing UKMO_CAP_Netcdf land mask!'
   
#else
   
   write(LDT_logunit,*) &
        '[ERR] netCDF4 required to process UKMO_CAP_Netcdf land mask!'
   write(LDT_logunit,*) &
        'Reconfigure LDT with netCDF4 support, recompile, and run again!'
   call LDT_endrun
   
#endif
   
   return
   
contains
   
   ! Internal subroutine to pull landmask from netCDF4 file
   subroutine get_landmask(n, nlon, nlat, landmask_cap)
      
      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
      use netcdf
#endif
      
      ! Defaults
      implicit none
      
      ! Arguments
      integer, intent(in) :: n
      integer, intent(out) :: nlon
      integer, intent(out) :: nlat
      real, allocatable, intent(out) :: landmask_cap(:,:)
      
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
      
      ! Local variables
      integer :: ncid, lat_dimid, lon_dimid, varid, xtype, ndims
      integer :: dimids(NF90_MAX_VAR_DIMS)
      integer :: status
      integer :: r,c
      
      ! Open file
      status = nf90_open(trim(LDT_rc%mfile(n)), NF90_NOWRITE, ncid)
      if (status .ne. NF90_NOERR) then
         write(LDT_logunit,*)'[ERR] NETCDF90 unable to open ', &
              trim(LDT_rc%mfile(n))
         write(LDT_logunit,*) trim(nf90_strerror(status))
         call LDT_endrun
      end if
      
      ! Get longitude dimension ID
      status = nf90_inq_dimid(ncid,"lon",lon_dimid)
      if (status .ne. NF90_NOERR) then
         write(LDT_logunit,*)&
              '[ERR] NETCDF90 unable to find lon dimension in ', &
              trim(LDT_rc%mfile(n))
         write(LDT_logunit,*) trim(nf90_strerror(status))
         call LDT_endrun
      end if
      
      ! Get latitude dimension ID
      status = nf90_inq_dimid(ncid,"lat",lat_dimid)
      if (status .ne. NF90_NOERR) then
         write(LDT_logunit,*)&
              '[ERR] NETCDF90 unable to find lat dimension in ', &
              trim(LDT_rc%mfile(n))
         write(LDT_logunit,*) trim(nf90_strerror(status))
         call LDT_endrun
      end if
      
      ! Get the variable ID
      status = nf90_inq_varid(ncid,"landmask",varid)
      if (status .ne. NF90_NOERR) then
         write(LDT_logunit,*)'[ERR] NETCDF90 unable to find landmask in ', &
              trim(LDT_rc%mfile(n))
         write(LDT_logunit,*) trim(nf90_strerror(status))
         call LDT_endrun
      end if
      
      ! Get longitude dimension length
      status = nf90_inquire_dimension(ncid,lon_dimid,len=nlon)
      if (status .ne. NF90_NOERR) then
         write(LDT_logunit,*)&
              '[ERR] NETCDF90 unable to find lon dimension length in ', &
              trim(LDT_rc%mfile(n))
         write(LDT_logunit,*) trim(nf90_strerror(status))
         call LDT_endrun
      end if
      
      ! Get latitude dimension length
      status = nf90_inquire_dimension(ncid,lat_dimid,len=nlat)
      if (status .ne. NF90_NOERR) then
         write(LDT_logunit,*)&
              '[ERR] NETCDF90 unable to find lat dimension length in ', &
              trim(LDT_rc%mfile(n))
         write(LDT_logunit,*) trim(nf90_strerror(status))
         call LDT_endrun
      end if
      
      ! Make sure landmask has appropriate dimensions and type
      status = nf90_inquire_variable(ncid,varid,xtype=xtype, &
           ndims=ndims, dimids=dimids)
      if (status .ne. NF90_NOERR) then
         write(LDT_logunit,*)&
              '[ERR] NETCDF90 unable to find landmask information in ', &
              trim(LDT_rc%mfile(n))
         write(LDT_logunit,*) trim(nf90_strerror(status))
         call LDT_endrun
      end if
      if (xtype .ne. NF90_FLOAT) then
         write(LDT_logunit,*)&
              '[ERR] Expected landmask to be type NF90_FLOAT, found ',xtype
         call LDT_endrun
      end if
      if (ndims .ne. 2) then
         write(LDT_logunit,*)&
              '[ERR] Expected landmask to be rank 2, found ',ndims
         call LDT_endrun
      end if
      if (dimids(1) .ne. lon_dimid) then
         write(LDT_logunit,*)&
              '[ERR] Expected landmask dimension 1 to be lon!'
         call LDT_endrun
      end if
      if (dimids(2) .ne. lat_dimid) then
         write(LDT_logunit,*)&
              '[ERR] Expected landmask dimension 2 to be lat!'
         call LDT_endrun
      end if
      
      ! Get land mask
      allocate(landmask_cap(nlon,nlat))
      status = nf90_get_var(ncid,varid,landmask_cap)
      if (status .ne. NF90_NOERR) then
         write(LDT_logunit,*)&
              '[ERR] NETCDF90 unable to find lat dimension length in ', &
              trim(LDT_rc%mfile(n))
         write(LDT_logunit,*) trim(nf90_strerror(status))
         call LDT_endrun
      end if
      
      ! Close file
      status = nf90_close(ncid)
      if (status .ne. NF90_NOERR) then
         write(LDT_logunit,*)&
              '[ERR] NETCDF90 unable to close file ', &
              trim(LDT_rc%mfile(n))
         write(LDT_logunit,*) trim(nf90_strerror(status))
         call LDT_endrun
      end if
#endif
      
   end subroutine get_landmask

   ! Internal subroutine to interpolate landmask to LIS grid
   subroutine interp_landmask(n, nlon, nlat, landmask_cap, localmask)
      
      ! Imports
      use LDT_coreMod,   only : LDT_rc
      use LDT_fileIOMod, only : LDT_transform_paramgrid
      use LDT_logMod,    only : LDT_logunit, LDT_endrun
      
      ! Defaults
      implicit none
      
      ! Arguments
      integer, intent(in)  :: n
      integer, intent(in)  :: nlon
      integer, intent(in)  :: nlat
      real,    intent(in)  :: landmask_cap(nlon,nlat)
      real,    intent(out) :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))

      ! Local variables
      integer   :: mi
      integer   :: mo
      real,    allocatable :: gi1(:) 
      logical*1,allocatable:: li1(:) 
      real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) 
      logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) 
      integer :: r,c,i

      ! Sanity check grid transform
      select case ( LDT_rc%mask_gridtransform(n) )
      case( "none", "mode", "neighbor", "tile")
         write(LDT_logunit,*) "[INFO] Reading UKMO_CAP_Netcdf land mask "
      case default
         write(LDT_logunit,*) &
              "[ERR] The spatial transform option selected for "
         write(LDT_logunit,*) &
              "   UKMO_CAP_Netcdf is not recognized nor recommended."
         write(LDT_logunit,*) &
              "   Please select: "
         write(LDT_logunit,*) &
              "  ==  none, mode, neighbor, tile "
         write(LDT_logunit,*) "Program stopping ..."
         call LDT_endrun
      end select

      ! Initialize downscaling/upscaling options
      mi = nlon*nlat
      allocate( gi1(mi), li1(mi) )
      gi1(:) = LDT_rc%udef
      li1(:) = .false.
      mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
      lo1(:) = .false.

      ! Assign 2-D array to 1-D for aggregation routines
      i = 0
      do r = 1, nlat
         do c = 1, nlon
            i = i + 1
            gi1(i) = landmask_cap(c,r)
            if (gi1(i) .ne. LDT_rc%udef) li1(i) = .true.
         enddo
      enddo

      ! Transform parameter from original grid to LIS output grid
      call LDT_transform_paramgrid(n, &
           LDT_rc%mask_gridtransform(n), &
           LDT_rc%mask_gridDesc(n,:), mi, 1, gi1, li1, mo, go1, lo1 )
      
      ! Convert 1D to 2D grid output arrays
      i = 0
      do r = 1, LDT_rc%lnr(n)
         do c = 1, LDT_rc%lnc(n)
            i = i + 1
            if (go1(i) < 0.) then
               localmask(c,r) = LDT_rc%udef
            else
               localmask(c,r) = go1(i)
            end if
         enddo
      enddo
      deallocate( li1, gi1 )
      
   end subroutine interp_landmask
   
   ! Internal subroutine for checking landmask against vegetation type and
   ! tile map
   subroutine check_landmask(n,vegtype,fgrd,localmask)
      
      ! Imports
      use LDT_coreMod,   only : LDT_rc
      use LDT_logMod,    only : LDT_logunit, LDT_endrun
      
      ! Defaults
      implicit none
      
      ! Arguments
      integer, intent(in)  :: n
      real,    intent(in)  :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
      real,    intent(in)  :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
      real,    intent(in) :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))
      
      ! Local variables
      integer :: cnt_0mask_1lc
      integer :: cnt_1mask_0lc
      integer :: c,r
      
      cnt_0mask_1lc = 0
      cnt_1mask_0lc = 0
      
      if( LDT_rc%inc_water_pts ) then
         ! Include water points
         
         !- NON-tiled option:
         select case( LDT_rc%lc_gridtransform(n) ) 
            
         case( "none", "mode" )
            
            do r = 1, LDT_rc%lnr(n)
               do c = 1, LDT_rc%lnc(n)
                  if ( vegtype(c,r) .eq. LDT_rc%waterclass .and. &
                       localmask(c,r) .ne. 0 ) then
                     cnt_1mask_0lc = cnt_1mask_0lc + 1
                  elseif( vegtype(c,r) .ne. LDT_rc%waterclass .and. &
                       localmask(c,r) .eq. 0 ) then
                     cnt_0mask_1lc = cnt_0mask_1lc + 1
                  endif
               end do ! c
            end do ! r
            
         case( "tile" )
            
            do r=1,LDT_rc%lnr(n)
               do c=1,LDT_rc%lnc(n)
                  if( fgrd(c,r,LDT_rc%waterclass) >= &
                       LDT_rc%gridcell_water_frac(n) .and. &  
                       localmask(c,r) == 1 ) then
                     cnt_1mask_0lc = cnt_1mask_0lc + 1
                  end if
                  if( fgrd(c,r,LDT_rc%waterclass) < &
                       LDT_rc%gridcell_water_frac(n) .and. & 
                       localmask(c,r) == 0 ) then
                     cnt_0mask_1lc = cnt_0mask_1lc + 1
                  endif
               end do ! c
            end do ! r
         end select
         
      else
         ! Do not include water points in tiled land cover map
         if ( LDT_rc%lc_gridtransform(n) == "tile" ) then
            do r=1,LDT_rc%lnr(n)
               do c=1,LDT_rc%lnc(n)
                  if( sum(fgrd(c,r,1:LDT_rc%nt)) <= &
                       LDT_rc%gridcell_water_frac(n)&
                       .and. localmask(c,r) == 1 ) then                     
                     cnt_1mask_0lc = cnt_1mask_0lc + 1
                  end if
               end do ! c
            end do ! r
         end if
         
      end if ! water points
      
      if( cnt_1mask_0lc > 0 .or. cnt_0mask_1lc > 0) then
         write(LDT_logunit,*)&
              " [WARN] MISMATCH between mask and landcover maps ..."
         write(LDT_logunit,*) &
              " Landmask LAND/Landcover WATER mismatch count: ", &
              cnt_1mask_0lc 
         write(LDT_logunit,*) &
              " Landmask WATER/Landcover LAND mismatch count: ", &
              cnt_0mask_1lc
         write(LDT_logunit,*) &
              " To make sure landmask and landcover fields agree, select:"
         write(LDT_logunit,*)&
              "  'neighbor' in the 'Landcover fill option:' entry ... "
      endif
      
   end subroutine check_landmask
end subroutine read_maskfile_UKMO_CAP_Netcdf
